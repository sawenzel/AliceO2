// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DetectorsBase/O2BVHSurfaceSolid.h"

#include "BoundedSurface.h"

// the third-party BVH headers plus extra kernels, shared with O2Tessellated
#include "bvh2_third_party.h"
#include "bvh2_extra_kernels.h"

#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>

using namespace o2::base;
using namespace o2::base::surface;
ClassImp(O2BVHSurfaceSolid);

namespace
{
// float BVH types following the O2Tessellated::BuildBVH pattern
using BVHScalar = float;
using BVHBBox = bvh::v2::BBox<BVHScalar, 3>;
using BVHVec3 = bvh::v2::Vec<BVHScalar, 3>;
using BVHNode = bvh::v2::Node<BVHScalar, 3>;
using BVH = bvh::v2::Bvh<BVHNode>;
using BVHRay = bvh::v2::Ray<BVHScalar, 3>;

Vec2 makeVec2(const O2BVHSurfaceSolid::Point2D& point)
{
  return {point[0], point[1]};
}

Vec3 makeVec3(const O2BVHSurfaceSolid::Point3D& point)
{
  return {point[0], point[1], point[2]};
}

Vec3 makeVec3(const Double_t* point)
{
  return {point[0], point[1], point[2]};
}

// The arbitrary skew test direction used for parity-based containment: probes all normals and
// avoids evident symmetries (same as O2Tessellated), normalized so hit distances are lengths.
const Vec3 kContainsTestDirection = normalized({1., 1.41421356237, 1.73205080757});

/// Directions the re-shoot votes over on a solid whose surface set is not a closed manifold; see
/// containsByVote. A golden-angle spiral is quasi-uniform on the sphere, so the directions are
/// mutually well separated and none aligns with a coordinate axis or a 45-degree symmetry plane.
/// Both properties matter: an earlier hand-picked triple failed on 13 of 55 known-bad points
/// precisely because two of its directions turned out correlated, and agreed on the wrong answer.
///
/// Five is where the measurement stops paying: over those 55 points the majority is right 50/55
/// with three directions, 53/55 with five, and 55/55 with thirteen. Five buys almost all of it at
/// a bounded cost, and the residue is a defect the vote is not meant to repair -- Phase 2 is.
const std::array<Vec3, 5>& reshootDirections()
{
  static const std::array<Vec3, 5> directions = [] {
    std::array<Vec3, 5> spiral{};
    for (int index = 0; index < 5; ++index) {
      const double cosTheta = 1. - 2. * (index + 0.5) / 5.;
      const double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
      const double phi = 2.399963229728653 * index; // golden angle
      spiral[index] = normalized({sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta});
    }
    return spiral;
  }();
  return directions;
}

// Ray tmax tightening in the BVH distance queries; see O2BVHSurfaceSolid::SetRayTMaxPruning.
bool gRayTMaxPruning = true;
// Per-thread diagnostic counter of leaf surface patches visited by the BVH distance queries.
thread_local long long gRayCandidateCount = 0;
// ... and by the nearest-patch queries behind Safety and ComputeNormal; see
// O2BVHSurfaceSolid::ResetSafetyCandidateCounter.
thread_local long long gSafetyCandidateCount = 0;
// Deliberately unsound node bound for the nearest-patch traversal; see
// O2BVHSurfaceSolid::SetSafetyBoundUnsoundForTest. Never true outside a test.
bool gSafetyBoundUnsound = false;

/// Squared distance from \a point to a BVH node's bounding box, evaluated in double from the
/// float box corners (which convert exactly), and then shrunk by a relative guard.
///
/// This is *the* quantity the nearest-patch traversal prunes on, so its error direction is the
/// whole safety argument. It must never come out **larger** than the true distance from the point
/// to any patch inside the node, or a nearer patch is discarded and Safety() answers too much --
/// the one failure mode that lets a navigator step through a wall. Two things could push it up
/// and both are held down here:
///
///  * the box arithmetic. Three subtractions, three squarings and two additions carry at most
///    about 4 ulps of relative error, so the result is scaled by (1 - 1e-12), three orders of
///    magnitude of headroom over that. Scaling *down* a lower bound keeps it a lower bound; it
///    costs only pruning, never correctness.
///  * the box itself, which is inflated by kBVHBoxTolerance and rounded outward in float by
///    buildBVH, i.e. it is a superset of the patch's conservative bounds. A larger box is nearer,
///    so that too errs downward.
///
/// The corresponding statement about the patch is in nearestPatchDistanceSq.
inline double boxDistanceSq(const BVHBBox& box, const Vec3& point)
{
  const double coordinates[3] = {point.xCoord, point.yCoord, point.zCoord};
  double distanceSq = 0.;
  for (int dimension = 0; dimension < 3; ++dimension) {
    const double lower = static_cast<double>(box.min[dimension]);
    const double upper = static_cast<double>(box.max[dimension]);
    const double value = coordinates[dimension];
    if (value < lower) {
      distanceSq += (lower - value) * (lower - value);
    } else if (value > upper) {
      distanceSq += (value - upper) * (value - upper);
    }
  }
  if (gSafetyBoundUnsound) {
    // The negative control: the distance to the box *centre* is not a lower bound on anything,
    // it is larger than the box distance whenever the box has any extent at all, so it prunes
    // nodes that hold the true nearest patch. Enabled only by SetSafetyBoundUnsoundForTest.
    double centreDistanceSq = 0.;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double centre =
        0.5 * (static_cast<double>(box.min[dimension]) + static_cast<double>(box.max[dimension]));
      const double gap = coordinates[dimension] - centre;
      centreDistanceSq += gap * gap;
    }
    return centreDistanceSq;
  }
  return distanceSq * (1. - 1.e-12);
}

/// Convert a double ray bound to the float the BVH traversal compares against, rounding *up* so
/// the float bound can never be smaller than the double one. Traversal in float must never
/// discard a node that the exact bound would have kept, otherwise a nearer hit is lost. Mirrors
/// the truncate_roundup lambda in O2Tessellated's distance queries.
inline BVHScalar truncateRoundUp(double bound)
{
  const double clamped = std::min(bound, static_cast<double>(std::numeric_limits<BVHScalar>::max()));
  const double biased = clamped + std::numeric_limits<BVHScalar>::epsilon() * std::abs(clamped);
  return static_cast<BVHScalar>(biased);
}

/// Lower bound on the ray parameter for the *distance* queries. Parity starts just past the
/// origin so it never re-counts the surface the point is sitting on, but a distance query
/// standing on a face has to be able to see that face's t = 0 crossing: answering "infinitely
/// far" to a navigator that is on the wall is how it tunnels through the geometry. So hits are
/// admitted from just *behind* the origin and the answer is clamped to zero, which is the
/// convention ROOT's own primitives follow. See CodeReview_Fable.md S3.
constexpr double kDistanceRayTolerance = -kRayTolerance;

/// Which side of the surface a single hit is on, as seen by a ray travelling along
/// \a rayDirection: the outward normal opposes the direction when entering. Within kTolerance of
/// tangency it is neither -- the surface is not actually crossed there.
///
/// Parity and the distance queries share this, which they did not: parity used a bare sign test
/// and the distance queries a tolerance band, so the two disagreed about exactly the hits that
/// are hardest to get right (CodeReview_Fable.md S5).
enum class CrossingSense { Entering,
                           Exiting,
                           Tangential };

/// Half-width of the window sameIntersection() treats as one intersection at ray parameter
/// \a distance, doubled for headroom. A per-surface search bound raised by this can never cut a
/// hit that would have clustered with one at \a distance.
inline double clusterMargin(double distance)
{
  return 2. * kIntersectionTolerance * std::max(1., std::abs(distance));
}

inline CrossingSense crossingSense(const RayHit& hit, const Vec3& rayDirection)
{
  const double alignment = dot(hit.normal, rayDirection);
  if (alignment < -kTolerance) {
    return CrossingSense::Entering;
  }
  if (alignment > kTolerance) {
    return CrossingSense::Exiting;
  }
  return CrossingSense::Tangential;
}

/// Sort \a hits and walk their clusters in increasing distance, calling
/// \a visitor(firstIndex, endIndex, sense) once per cluster; the visitor returns false to stop.
///
/// Near-equal hits are clustered because a shared edge or vertex is seen by several surfaces at
/// the same point. A cluster whose members all enter, or all exit, is one genuine crossing. A
/// cluster carrying both is a *graze*: the ray touched an edge and came back out on the side it
/// came from, so it crossed nothing, and the cluster reports Tangential. That distinction is what
/// keeps a rim graze from counting as a boundary -- and it is why this walk is shared with the
/// distance queries, which used to classify every hit on its own and so disagreed with parity
/// about the same ray (CodeReview_Fable.md S4).
template <typename ClusterVisitor>
void forEachCrossingCluster(std::vector<RayHit>& hits, const Vec3& rayDirection, ClusterVisitor&& visitor)
{
  std::sort(hits.begin(), hits.end(),
            [](const RayHit& firstHit, const RayHit& secondHit) { return firstHit.distance < secondHit.distance; });

  size_t hitIndex = 0;
  while (hitIndex < hits.size()) {
    bool entering = false;
    bool exiting = false;
    size_t clusterEnd = hitIndex;
    // Compared against the cluster's *first* member, not its predecessor. Chaining neighbour to
    // neighbour is transitive: N hits spaced just inside the window merge into one cluster N
    // windows wide, so a thin feature at large ray parameter -- where the window is relative and
    // therefore large -- silently collapses into a single mixed cluster and stops being a
    // boundary. Anchoring bounds every cluster to one window (CodeReview_Fable.md K8).
    while (clusterEnd < hits.size() &&
           (clusterEnd == hitIndex || sameIntersection(hits[clusterEnd].distance, hits[hitIndex].distance))) {
      switch (crossingSense(hits[clusterEnd], rayDirection)) {
        case CrossingSense::Entering:
          entering = true;
          break;
        case CrossingSense::Exiting:
          exiting = true;
          break;
        case CrossingSense::Tangential:
          break;
      }
      ++clusterEnd;
    }
    // both, or neither: nothing was crossed
    const CrossingSense sense = entering == exiting ? CrossingSense::Tangential
                                                    : (entering ? CrossingSense::Entering : CrossingSense::Exiting);
    if (!visitor(hitIndex, clusterEnd, sense)) {
      return;
    }
    hitIndex = clusterEnd;
  }
}

/// Distance to the nearest genuine entering (\a wantEntering) or exiting crossing in \a hits, or
/// Big when there is none. \a grazedFirst reports whether a cancelled (grazing) cluster was passed
/// on the way, which the BVH query needs to know: it prunes the traversal at the nearest
/// *candidate* hit, and a candidate that then cancels means the real answer may lie in the part of
/// the ray that pruning threw away.
template <bool wantEntering>
double nearestCrossingInHits(std::vector<RayHit>& hits, const Vec3& rayDirection, bool& grazedFirst)
{
  constexpr CrossingSense wanted = wantEntering ? CrossingSense::Entering : CrossingSense::Exiting;
  double distance = TGeoShape::Big();
  grazedFirst = false;
  forEachCrossingCluster(hits, rayDirection, [&](size_t firstIndex, size_t, CrossingSense sense) {
    if (sense == CrossingSense::Tangential) {
      grazedFirst = true;
      return true;
    }
    if (sense != wanted) {
      return true;
    }
    // clusters come in increasing distance, so the first match is the answer. Hits are admitted
    // from just behind the origin (kDistanceRayTolerance), so clamp: a boundary crossing is a
    // zero step, never a negative one.
    distance = std::max(0., hits[firstIndex].distance);
    return false;
  });
  return distance;
}

/// @name Persistent surface records
///
/// Translation between the public Add*Surface arguments and the flat BVHSurfaceRecord the solid
/// streams. Records are appended only on success, so replaying them reproduces exactly the
/// surface list the original solid ended up with -- including any face that the original refused,
/// which stays refused rather than reappearing.
/// @{

void fillPoint3(double (&target)[3], const O2BVHSurfaceSolid::Point3D& source)
{
  target[0] = source[0];
  target[1] = source[1];
  target[2] = source[2];
}

O2BVHSurfaceSolid::Point3D makePoint3D(const double (&source)[3])
{
  return {source[0], source[1], source[2]};
}

/// The frame-and-scalars part of a record, shared by all six surface families.
BVHSurfaceRecord makeRecord(int kind, const O2BVHSurfaceSolid::Point3D& origin,
                            const O2BVHSurfaceSolid::Point3D& axisA, const O2BVHSurfaceSolid::Point3D& axisB,
                            std::vector<double> scalars, bool innerWall, bool trimmed)
{
  BVHSurfaceRecord record;
  record.kind = kind;
  fillPoint3(record.origin, origin);
  fillPoint3(record.axisA, axisA);
  fillPoint3(record.axisB, axisB);
  record.scalars = std::move(scalars);
  record.innerWall = innerWall;
  record.trimmed = trimmed;
  return record;
}

BVHSurfaceCurveRecord makeCurveRecord(const O2BVHSurfaceSolid::PlanarBoundaryCurve& curve)
{
  BVHSurfaceCurveRecord record;
  record.kind = static_cast<int>(curve.kind);
  record.lineStart[0] = curve.lineStart[0];
  record.lineStart[1] = curve.lineStart[1];
  record.lineEnd[0] = curve.lineEnd[0];
  record.lineEnd[1] = curve.lineEnd[1];
  record.center[0] = curve.center[0];
  record.center[1] = curve.center[1];
  record.radius = curve.radius;
  record.startAngle = curve.startAngle;
  record.endAngle = curve.endAngle;
  record.degree = curve.degree;
  record.poles.reserve(2 * curve.poles.size());
  for (const auto& pole : curve.poles) {
    record.poles.push_back(pole[0]);
    record.poles.push_back(pole[1]);
  }
  record.weights = curve.weights;
  record.knots = curve.knots;
  return record;
}

O2BVHSurfaceSolid::PlanarBoundaryCurve makeBoundaryCurve(const BVHSurfaceCurveRecord& record)
{
  O2BVHSurfaceSolid::PlanarBoundaryCurve curve;
  curve.kind = static_cast<O2BVHSurfaceSolid::PlanarBoundaryCurve::Kind>(record.kind);
  curve.lineStart = {record.lineStart[0], record.lineStart[1]};
  curve.lineEnd = {record.lineEnd[0], record.lineEnd[1]};
  curve.center = {record.center[0], record.center[1]};
  curve.radius = record.radius;
  curve.startAngle = record.startAngle;
  curve.endAngle = record.endAngle;
  curve.degree = record.degree;
  curve.poles.reserve(record.poles.size() / 2);
  for (size_t index = 0; index + 1 < record.poles.size(); index += 2) {
    curve.poles.push_back({record.poles[index], record.poles[index + 1]});
  }
  curve.weights = record.weights;
  curve.knots = record.knots;
  return curve;
}

/// Store an outer wire plus its holes as one flat curve list with per-wire sizes.
void storeCurveWires(BVHSurfaceRecord& record, const std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& outerWire,
                     const std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>>& innerWires)
{
  record.wireSizes.push_back(static_cast<int>(outerWire.size()));
  for (const auto& curve : outerWire) {
    record.curves.push_back(makeCurveRecord(curve));
  }
  for (const auto& innerWire : innerWires) {
    record.wireSizes.push_back(static_cast<int>(innerWire.size()));
    for (const auto& curve : innerWire) {
      record.curves.push_back(makeCurveRecord(curve));
    }
  }
}

/// The inverse of storeCurveWires. Returns false when the per-wire sizes do not add up to the
/// stored curve count, i.e. when the record is truncated or corrupt.
bool loadCurveWires(const BVHSurfaceRecord& record, std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& outerWire,
                    std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>>& innerWires)
{
  size_t consumed = 0;
  for (size_t wireIndex = 0; wireIndex < record.wireSizes.size(); ++wireIndex) {
    const int wireSize = record.wireSizes[wireIndex];
    if (wireSize < 0 || consumed + static_cast<size_t>(wireSize) > record.curves.size()) {
      return false;
    }
    auto& wire = wireIndex == 0 ? outerWire : innerWires.emplace_back();
    for (int curveIndex = 0; curveIndex < wireSize; ++curveIndex) {
      wire.push_back(makeBoundaryCurve(record.curves[consumed + curveIndex]));
    }
    consumed += static_cast<size_t>(wireSize);
  }
  return consumed == record.curves.size();
}

/// storeCurveWires/loadCurveWires for the polygon-vertex flavour of a planar surface.
void storePolygonWires(BVHSurfaceRecord& record, const std::vector<O2BVHSurfaceSolid::Point2D>& outerWire,
                       const std::vector<std::vector<O2BVHSurfaceSolid::Point2D>>& innerWires)
{
  const auto append = [&record](const std::vector<O2BVHSurfaceSolid::Point2D>& wire) {
    record.wireSizes.push_back(static_cast<int>(wire.size()));
    for (const auto& vertex : wire) {
      record.polygonPoints.push_back(vertex[0]);
      record.polygonPoints.push_back(vertex[1]);
    }
  };
  append(outerWire);
  for (const auto& innerWire : innerWires) {
    append(innerWire);
  }
}

bool loadPolygonWires(const BVHSurfaceRecord& record, std::vector<O2BVHSurfaceSolid::Point2D>& outerWire,
                      std::vector<std::vector<O2BVHSurfaceSolid::Point2D>>& innerWires)
{
  size_t consumed = 0;
  for (size_t wireIndex = 0; wireIndex < record.wireSizes.size(); ++wireIndex) {
    const int wireSize = record.wireSizes[wireIndex];
    if (wireSize < 0 || 2 * (consumed + static_cast<size_t>(wireSize)) > record.polygonPoints.size()) {
      return false;
    }
    auto& wire = wireIndex == 0 ? outerWire : innerWires.emplace_back();
    for (int vertexIndex = 0; vertexIndex < wireSize; ++vertexIndex) {
      const size_t offset = 2 * (consumed + vertexIndex);
      wire.push_back({record.polygonPoints[offset], record.polygonPoints[offset + 1]});
    }
    consumed += static_cast<size_t>(wireSize);
  }
  return 2 * consumed == record.polygonPoints.size();
}
/// @}

// Ray-parity evaluation of a full intersection list (sorts hits in place). Near-equal
// intersections are clustered (shared edges/vertices seen by several surfaces). A cluster whose
// hits all enter or all exit is one genuine crossing; a mixed cluster means the ray grazes an
// edge or corner tangentially (e.g. the rim where a cap meets a wall) and must contribute even
// parity.
bool oddCrossingParity(std::vector<RayHit>& hits, const Vec3& rayDirection)
{
  int crossings = 0;
  forEachCrossingCluster(hits, rayDirection, [&](size_t, size_t, CrossingSense sense) {
    if (sense != CrossingSense::Tangential) {
      ++crossings;
    }
    return true;
  });
  return (crossings & 1) != 0;
}

/// A rim's state, on the scale the solid reports. The two enums carry the same four states -- a
/// rim's verdict *is* what the solid would report if that rim were its only defect -- so the solid
/// state is the worst over its rims, and there is one set of names for both.
O2BVHSurfaceSolid::NavigationReliability rimStateToReliability(RimState state)
{
  using Reliability = O2BVHSurfaceSolid::NavigationReliability;
  switch (state) {
    case RimState::Matched:
      return Reliability::Reliable;
    case RimState::Reversed:
      return Reliability::ReversedFaces;
    case RimState::Boundary:
      return Reliability::OpenSurfaceSet;
    case RimState::NonManifold:
      return Reliability::NonManifold;
  }
  return Reliability::Undetermined;
}
} // namespace

struct O2BVHSurfaceSolid::Impl {
  std::vector<std::unique_ptr<BoundedSurface>> surfaces;
  std::vector<Vec3> displayVertices;
  std::vector<std::array<int, 3>> displayTriangles;
  /// Which surface each display triangle came from, parallel to displayTriangles. The display mesh
  /// is assembled by calling appendDisplayMesh once per bounded surface, so the provenance is free
  /// at build time -- and it is the only thing that lets a point sampled inside a triangle be put
  /// back onto the *exact* patch that triangle approximates. Without it a sampled point is stuck
  /// on the polygon. See GetPointsOnSegments.
  std::vector<int> displayTriangleSurface;
  ClosureReport closure;
  /// closure.rimRecords in the public form, built once by CloseShape so the accessor can hand out
  /// a reference. The two are the same data; only the state enum and the Vec3 differ in type.
  std::vector<RimReport> rimReports;
  bool defined = false;
  std::unique_ptr<BVH> bvh; //!< acceleration structure over the surface AABBs (built in CloseShape)

  /// Build the BVH over per-surface conservative AABBs. Boxes are expanded by the documented
  /// kBVHBoxTolerance and then rounded outward during the double->float conversion so no
  /// tolerance-level boundary hit can be missed by the float traversal.
  void buildBVH()
  {
    bvh.reset();
    if (surfaces.empty()) {
      return;
    }

    std::vector<BVHBBox> primitiveBoxes;
    std::vector<BVHVec3> primitiveCenters;
    primitiveBoxes.reserve(surfaces.size());
    primitiveCenters.reserve(surfaces.size());
    for (const auto& surface : surfaces) {
      Vec3 lowerCorner{TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};
      Vec3 upperCorner{-TGeoShape::Big(), -TGeoShape::Big(), -TGeoShape::Big()};
      surface->conservativeBounds(lowerCorner, upperCorner);
      BVHBBox primitiveBox;
      for (int dimension = 0; dimension < 3; ++dimension) {
        primitiveBox.min[dimension] = std::nextafterf(
          static_cast<float>(component(lowerCorner, dimension) - kBVHBoxTolerance),
          -std::numeric_limits<float>::infinity());
        primitiveBox.max[dimension] = std::nextafterf(
          static_cast<float>(component(upperCorner, dimension) + kBVHBoxTolerance),
          std::numeric_limits<float>::infinity());
      }
      primitiveBoxes.push_back(primitiveBox);
      primitiveCenters.emplace_back(primitiveBox.get_center());
    }

    typename bvh::v2::DefaultBuilder<BVHNode>::Config config;
    config.quality = bvh::v2::DefaultBuilder<BVHNode>::Quality::High;
    // One analytic patch per leaf: patch intersections are far more expensive than node box
    // tests (unlike triangles), and the bvh2 traversal enters a leaf start node without a box
    // test, so a default-sized (up to 8 primitives) single-leaf tree would defeat all pruning
    // for small solids.
    config.max_leaf_size = 1;
    bvh = std::make_unique<BVH>(bvh::v2::DefaultBuilder<BVHNode>::build(primitiveBoxes, primitiveCenters, config));
  }

  /// Visit every surface whose BVH leaf box is traversed by the (unbounded) ray.
  template <typename SurfaceVisitor>
  void visitRayCandidates(const Vec3& rayOrigin, const Vec3& rayDirection, SurfaceVisitor&& visitor) const
  {
    BVHRay ray(BVHVec3(rayOrigin.xCoord, rayOrigin.yCoord, rayOrigin.zCoord),
               BVHVec3(rayDirection.xCoord, rayDirection.yCoord, rayDirection.zCoord), 0.f,
               std::numeric_limits<BVHScalar>::max());
    static constexpr bool useRobustTraversal = true;
    bvh::v2::GrowingStack<BVH::Index> stack;
    bvh->intersect<false, useRobustTraversal>(ray, bvh->get_root().index, stack,
                                              [&](size_t beginPrimitive, size_t endPrimitive) {
                                                for (size_t primitive = beginPrimitive; primitive < endPrimitive;
                                                     ++primitive) {
                                                  visitor(*surfaces[bvh->prim_ids[primitive]]);
                                                }
                                                return false; // keep traversing
                                              });
  }

  /// Distance to the nearest entering (\a wantEntering) or exiting crossing of the ray with the
  /// solid, bounded by \a stepmax; TGeoShape::Big() when there is none within that bound.
  ///
  /// The traversal bound shrinks to the best crossing found so far (unless pruning is switched
  /// off for benchmarking), which genuinely prunes: bvh2's node test re-reads ray.tmax on every
  /// node, so tightening it from inside the leaf callback takes effect immediately. It is safe
  /// because every hit lies inside its own leaf box, so a nearer hit's box is always entered at
  /// a ray parameter below the current best -- and the new bound is rounded up by
  /// kBVHBoxTolerance plus a float ulp so neither the double->float conversion nor the float
  /// box arithmetic can round it below a surviving candidate.
  ///
  /// The per-surface hit buffer is a function-local thread_local: template instantiation gives
  /// the entering and the exiting query a buffer each, so neither can be re-entered through the
  /// other's, and the capacity is paid once per thread rather than once per call.
  template <bool wantEntering>
  double nearestCrossing(const Vec3& rayOrigin, const Vec3& rayDirection, double stepmax) const
  {
    static thread_local std::vector<RayHit> surfaceHits;
    static thread_local std::vector<RayHit> collectedHits;
    constexpr CrossingSense wanted = wantEntering ? CrossingSense::Entering : CrossingSense::Exiting;

    // Deciding whether a hit is a real crossing needs its *neighbours*, since a cluster carrying
    // both an entering and an exiting hit is a graze and crosses nothing. So every hit is
    // collected and classified afterwards, not accepted or rejected on its own.
    //
    // That interacts with the tmax tightening: the traversal still shrinks to the nearest
    // candidate crossing, which is what makes the query cheap, but if that candidate then turns
    // out to be a graze the answer lies somewhere pruning already discarded. Rare enough to
    // handle by simply redoing the query without pruning, and doing so keeps the optimization an
    // optimization -- the two passes must return the same number.
    for (int attempt = 0; attempt < 2; ++attempt) {
      const bool pruning = gRayTMaxPruning && attempt == 0;
      collectedHits.clear();

      double bestCandidate = TGeoShape::Big();
      BVHRay ray(BVHVec3(rayOrigin.xCoord, rayOrigin.yCoord, rayOrigin.zCoord),
                 BVHVec3(rayDirection.xCoord, rayDirection.yCoord, rayDirection.zCoord), 0.f,
                 truncateRoundUp(stepmax));
      static constexpr bool useRobustTraversal = true;

      bvh::v2::GrowingStack<BVH::Index> stack;
      // ray is captured by reference on purpose: bvh2 takes it as const Ray&, but the object
      // itself is ours and mutable, and the traversal reads tmax afresh at every node test.
      bvh->intersect<false, useRobustTraversal>(
        ray, bvh->get_root().index, stack, [&](size_t beginPrimitive, size_t endPrimitive) {
          for (size_t primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
            const BoundedSurface& surface = *surfaces[bvh->prim_ids[primitive]];
            ++gRayCandidateCount;
            surfaceHits.clear();
            // the per-surface bound keeps a margin past the current candidate, so the candidate's
            // own cluster partners -- which decide whether it is a crossing at all -- are never
            // the hits that get cut
            const double bound =
              pruning ? std::min(stepmax, bestCandidate + clusterMargin(bestCandidate)) : stepmax;
            surface.appendIntersections(rayOrigin, rayDirection, kDistanceRayTolerance, bound, surfaceHits);
            for (const auto& hit : surfaceHits) {
              collectedHits.push_back(hit);
              if (crossingSense(hit, rayDirection) == wanted && hit.distance < bestCandidate) {
                bestCandidate = hit.distance;
              }
            }
          }
          if (pruning && bestCandidate < stepmax) {
            ray.tmax = std::min(ray.tmax, truncateRoundUp(bestCandidate + kBVHBoxTolerance));
          }
          return false; // keep traversing; the shrunk tmax does the pruning
        });

      bool grazedFirst = false;
      const double distance = nearestCrossingInHits<wantEntering>(collectedHits, rayDirection, grazedFirst);
      if (!pruning || !grazedFirst) {
        return distance;
      }
    }
    return TGeoShape::Big(); // unreachable: the second attempt never prunes
  }

  /// Same query without the BVH: visit every surface. Oracle and baseline for nearestCrossing.
  template <bool wantEntering>
  double nearestCrossingLoop(const Vec3& rayOrigin, const Vec3& rayDirection, double stepmax) const
  {
    static thread_local std::vector<RayHit> surfaceHits;
    static thread_local std::vector<RayHit> collectedLoopHits;

    // No pruning at all here: this is the oracle the accelerated query is checked against, so it
    // trades the shrinking upper bound for having every hit in hand and needing no retry.
    collectedLoopHits.clear();
    for (const auto& surface : surfaces) {
      surfaceHits.clear();
      surface->appendIntersections(rayOrigin, rayDirection, kDistanceRayTolerance, stepmax, surfaceHits);
      collectedLoopHits.insert(collectedLoopHits.end(), surfaceHits.begin(), surfaceHits.end());
    }
    bool grazedFirst = false;
    return nearestCrossingInHits<wantEntering>(collectedLoopHits, rayDirection, grazedFirst);
  }

  /// Parity of the crossings of the ray (\a point, \a direction) with the whole surface set: the
  /// containment answer for one shooting direction. \a useBVH selects the accelerated candidate
  /// set or the plain all-surfaces loop; both must produce the same hit multiset and hence the
  /// same answer, which is what Contains vs Contains_Loop cross-validates.
  ///
  /// Shared by Contains, Contains_Loop and the re-shoot, so that all three evaluate parity by
  /// exactly the same rule -- a re-shoot that disagreed with the primary shot for reasons other
  /// than the direction would be untestable.
  ///
  /// \a ambiguous, when given, reports whether any counted crossing carried RayHit::onTrimBoundary
  /// -- that is, whether this shot's answer rests anywhere on a tie-break the data does not
  /// support. It says nothing about whether the answer is wrong, only that this ray is one of the
  /// rays that *can* be wrong on a solid that is otherwise sound.
  bool parityAlong(const Vec3& point, const Vec3& direction, bool useBVH, bool* ambiguous = nullptr) const
  {
    // reused across calls so containment allocates nothing on the hot path; the capacity is paid
    // once per thread. Distinct from the distance queries' buffers, which are their own.
    static thread_local std::vector<RayHit> parityHits;
    parityHits.clear();
    if (useBVH) {
      visitRayCandidates(point, direction, [&](const BoundedSurface& surface) {
        surface.appendIntersections(point, direction, kRayTolerance, TGeoShape::Big(), parityHits);
      });
    } else {
      for (const auto& surface : surfaces) {
        surface->appendIntersections(point, direction, kRayTolerance, TGeoShape::Big(), parityHits);
      }
    }
    if (ambiguous != nullptr) {
      *ambiguous = std::any_of(parityHits.begin(), parityHits.end(),
                               [](const RayHit& hit) { return hit.onTrimBoundary; });
    }
    return oddCrossingParity(parityHits, direction);
  }

  /// Containment by majority vote over reshootDirections(): the re-shoot of CodeReview_Fable.md
  /// Section 4.4, for solids whose surface set is not a closed 2-manifold.
  ///
  /// A gap in the surface set costs the parity ray exactly the crossings that fall inside the gap,
  /// so a point is misclassified over the gap's whole *shadow* along the shooting direction --
  /// which is why the errors are centimetres wide and far from any surface, and equally why they
  /// are escapable: the shadow belongs to the direction, not to the point. Voting over several
  /// spread directions puts the point outside most of the shadows it is in.
  ///
  /// This mitigates a defect, it does not repair one. The answer is still undefined in the sense
  /// that IsNavigable() reports; what changes is that a single unlucky aim no longer decides it.
  /// Stops as soon as a majority is settled, so a point all directions agree on costs three shots
  /// rather than five.
  /// \a allTiedOnBoundary, when given, reports that the vote never found a direction whose parity
  /// rested on the geometry alone -- see below.
  bool containsByVote(const Vec3& point, bool useBVH, bool* allTiedOnBoundary = nullptr) const
  {
    constexpr int kMajority = 3; // of the five directions
    int inside = 0;              // shots whose parity rests on no trim-boundary tie-break
    int outside = 0;
    int insideOnBoundary = 0; // and shots that do, counted apart
    int outsideOnBoundary = 0;
    for (const auto& direction : reshootDirections()) {
      bool ambiguous = false;
      const bool answer = parityAlong(point, direction, useBVH, &ambiguous);
      if (ambiguous) {
        answer ? ++insideOnBoundary : ++outsideOnBoundary;
      } else {
        answer ? ++inside : ++outside;
      }
      if (inside >= kMajority || outside >= kMajority) {
        break;
      }
    }
    if (allTiedOnBoundary != nullptr) {
      *allTiedOnBoundary = (inside == outside);
    }
    // A shot whose parity rests on a tie-break the data does not support is not evidence of the
    // same quality as one that does not, and the vote used to count the two alike -- so a point
    // could be decided three to two by three coin flips, and report nothing. Decide among the
    // shots that stand on the geometry whenever they are not themselves tied; only a genuine tie
    // (which includes the case where every direction was ambiguous) falls back to counting all
    // five, because then there is nothing better to count.
    if (inside != outside) {
      return inside > outside;
    }
    return (inside + insideOnBoundary) > (outside + outsideOnBoundary);
  }

  /// Visit every surface whose BVH leaf box contains the point; stops early (returning true)
  /// once the visitor returns true. The leaf boxes are expanded by kBVHBoxTolerance, so a point
  /// within navigation tolerance of a patch is never missed.
  template <typename SurfaceVisitor>
  bool visitPointCandidates(const Vec3& point, SurfaceVisitor&& visitor) const
  {
    const BVHVec3 testPoint(point.xCoord, point.yCoord, point.zCoord);
    std::vector<size_t> nodeStack{0}; // start from the root node
    while (!nodeStack.empty()) {
      const auto& node = bvh->nodes[nodeStack.back()];
      nodeStack.pop_back();
      if (!bvh::v2::extra::contains(node.get_bbox(), testPoint)) {
        continue;
      }
      if (node.is_leaf()) {
        const auto beginPrimitive = node.index.first_id();
        const auto endPrimitive = beginPrimitive + node.index.prim_count();
        for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
          if (visitor(*surfaces[bvh->prim_ids[primitive]])) {
            return true;
          }
        }
      } else {
        const auto firstChild = node.index.first_id();
        for (size_t child : {firstChild, firstChild + 1}) {
          if (child < bvh->nodes.size()) {
            nodeStack.push_back(child);
          }
        }
      }
    }
    return false;
  }

  /// The brute-force nearest patch: the squared distanceSqToPatch minimum over *every* surface,
  /// and the index that realises it. The oracle for nearestPatchDistanceSq, and the definition of
  /// what Safety() and ComputeNormal() answer -- the accelerated traversal is required to return
  /// this number and this index bit for bit, not merely something close.
  ///
  /// The tie-break is part of the definition: a strict `<` means the *lowest-indexed* patch wins
  /// an exact tie, which is routine rather than exotic (the centre of a box is exactly equidistant
  /// from all six faces). ComputeNormal turns the winner into a normal, so a different winner is a
  /// different answer.
  double nearestPatchDistanceSqLoop(const Vec3& point, size_t* closestIndex) const
  {
    double bestDistanceSq = std::numeric_limits<double>::infinity();
    size_t bestIndex = surfaces.size();
    for (size_t index = 0; index < surfaces.size(); ++index) {
      const double patchDistanceSq = surfaces[index]->distanceSqToPatch(point);
      if (patchDistanceSq < bestDistanceSq) {
        bestDistanceSq = patchDistanceSq;
        bestIndex = index;
      }
    }
    if (closestIndex != nullptr) {
      *closestIndex = bestIndex;
    }
    return bestDistanceSq;
  }

  /// Same answer as nearestPatchDistanceSqLoop, reached through the BVH.
  ///
  /// Ordered depth-first traversal with a running best: a node is entered only if the distance
  /// from the point to its bounding box -- a lower bound on the distance to every patch in its
  /// subtree -- is still competitive, and children are visited nearer-box-first so the best
  /// tightens as early as possible. Each stack entry carries the bound it was pushed with and is
  /// re-tested on pop, so a node queued before a better patch was found is still dropped.
  ///
  /// **Why the bound is a lower bound, patch side.** Every distanceSqToPatch in BoundedSurface.h
  /// is realised at a point of the patch's *untrimmed* window -- the polygon or curve wire itself
  /// for the planar families, the full rim band for a cylinder or cone, the full sphere, the full
  /// torus -- and each family's conservativeBounds() encloses exactly that window. So the point
  /// realising distanceSqToPatch lies inside the surface's own leaf box, and the distance to the
  /// box cannot exceed it. Note that for a wire-trimmed quadric distanceSqToPatch is itself only a
  /// lower bound on the true distance to the trimmed patch; that is the existing, documented
  /// behaviour of the kernel and it errs the safe way. This traversal reproduces the loop's answer
  /// exactly, whatever that answer is; it does not inherit a second source of error.
  ///
  /// **Direction of every inequality here.** Both the node bound (boxDistanceSq) and the patch
  /// distance err *downward* against the true geometric distance, so Safety() can only ever come
  /// out too small. Too small costs a navigator an extra boundary computation; too large lets it
  /// walk through a wall.
  ///
  /// \a TrackIndex selects the pruning that ComputeNormal needs. Without it a node whose bound
  /// equals the current best is dropped, which cannot change the minimum. With it such a node must
  /// still be entered, because it may hold an equally-near patch of *lower index*, which the
  /// tie-break would prefer.
  template <bool TrackIndex>
  double nearestPatchDistanceSq(const Vec3& point, size_t* closestIndex) const
  {
    if (bvh == nullptr) {
      return nearestPatchDistanceSqLoop(point, closestIndex);
    }

    struct StackEntry {
      size_t node;
      double lowerBoundSq;
    };
    // reused across calls so the hot path allocates nothing; capacity is paid once per thread
    static thread_local std::vector<StackEntry> nodeStack;
    nodeStack.clear();

    double bestDistanceSq = std::numeric_limits<double>::infinity();
    size_t bestIndex = surfaces.size();

    auto pruned = [](double lowerBoundSq, double bestSoFarSq) {
      return TrackIndex ? lowerBoundSq > bestSoFarSq : lowerBoundSq >= bestSoFarSq;
    };

    nodeStack.push_back({0, boxDistanceSq(bvh->nodes[0].get_bbox(), point)});
    while (!nodeStack.empty()) {
      const StackEntry entry = nodeStack.back();
      nodeStack.pop_back();
      if (pruned(entry.lowerBoundSq, bestDistanceSq)) {
        continue;
      }
      const auto& node = bvh->nodes[entry.node];
      if (node.is_leaf()) {
        const auto beginPrimitive = node.index.first_id();
        const auto endPrimitive = beginPrimitive + node.index.prim_count();
        for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
          const size_t index = bvh->prim_ids[primitive];
          ++gSafetyCandidateCount;
          const double patchDistanceSq = surfaces[index]->distanceSqToPatch(point);
          if (patchDistanceSq < bestDistanceSq) {
            bestDistanceSq = patchDistanceSq;
            bestIndex = index;
          } else if (TrackIndex && patchDistanceSq == bestDistanceSq && index < bestIndex) {
            bestIndex = index;
          }
        }
        continue;
      }
      const size_t firstChild = node.index.first_id();
      const size_t secondChild = firstChild + 1;
      if (secondChild >= bvh->nodes.size()) {
        if (firstChild < bvh->nodes.size()) {
          nodeStack.push_back({firstChild, boxDistanceSq(bvh->nodes[firstChild].get_bbox(), point)});
        }
        continue;
      }
      double nearBound = boxDistanceSq(bvh->nodes[firstChild].get_bbox(), point);
      double farBound = boxDistanceSq(bvh->nodes[secondChild].get_bbox(), point);
      size_t nearChild = firstChild;
      size_t farChild = secondChild;
      if (farBound < nearBound) {
        std::swap(nearBound, farBound);
        std::swap(nearChild, farChild);
      }
      // farther child first: the stack is LIFO, so the nearer one is popped -- and tightens the
      // best -- before the farther one is re-tested
      if (!pruned(farBound, bestDistanceSq)) {
        nodeStack.push_back({farChild, farBound});
      }
      if (!pruned(nearBound, bestDistanceSq)) {
        nodeStack.push_back({nearChild, nearBound});
      }
    }

    if (closestIndex != nullptr) {
      *closestIndex = bestIndex;
    }
    return bestDistanceSq;
  }
};

int BVHSurfaceRecord::expectedScalarCount(int recordKind)
{
  switch (recordKind) {
    case PlanarPolygon:
    case CurvedPlanar:
      return 0;
    case Cylindrical: // radius, heightMin, heightMax, phiStart, phiSweep
    case Spherical:   // radius, thetaMin, thetaMax, phiStart, phiSweep
      return 5;
    case Conical: // radiusAtMin, radiusAtMax, heightMin, heightMax, phiStart, phiSweep
    case Toroidal: // majorRadius, minorRadius, phiStart, phiSweep, tubeStart, tubeSweep
      return 6;
    default:
      return -1;
  }
}

O2BVHSurfaceSolid::O2BVHSurfaceSolid() : TGeoBBox(), fImpl(new Impl)
{
}

O2BVHSurfaceSolid::O2BVHSurfaceSolid(const char* name) : TGeoBBox(name, 0., 0., 0.), fImpl(new Impl)
{
}

O2BVHSurfaceSolid::~O2BVHSurfaceSolid()
{
  delete fImpl;
}

bool O2BVHSurfaceSolid::AddPlanarSurface(const Point3D& origin, const Point3D& axisU, const Point3D& axisV,
                                         const std::vector<Point2D>& outerWire,
                                         const std::vector<std::vector<Point2D>>& innerWires)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddPlanarSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  std::vector<Vec2> convertedOuterWire;
  convertedOuterWire.reserve(outerWire.size());
  for (const auto& vertex : outerWire) {
    convertedOuterWire.push_back(makeVec2(vertex));
  }

  std::vector<std::vector<Vec2>> convertedInnerWires;
  convertedInnerWires.reserve(innerWires.size());
  for (const auto& innerWire : innerWires) {
    auto& convertedInnerWire = convertedInnerWires.emplace_back();
    convertedInnerWire.reserve(innerWire.size());
    for (const auto& vertex : innerWire) {
      convertedInnerWire.push_back(makeVec2(vertex));
    }
  }

  auto surface = std::make_unique<PlanarBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(origin), makeVec3(axisU), makeVec3(axisV), convertedOuterWire, convertedInnerWires,
                           errorMessage)) {
    Error("AddPlanarSurface", "%s", errorMessage.c_str());
    return false;
  }
  if (surface->wasReoriented()) {
    Warning("AddPlanarSurface", "Shape %s: planar surface %d had a wire re-oriented to match its role", GetName(),
            static_cast<int>(fImpl->surfaces.size()));
  }

  auto record = makeRecord(BVHSurfaceRecord::PlanarPolygon, origin, axisU, axisV, {}, false, false);
  storePolygonWires(record, outerWire, innerWires);
  fRecords.push_back(std::move(record));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

namespace
{
/// Translate a public PlanarBoundaryCurve wire into the internal Curve2D loop.
std::vector<Curve2D> makeCurveWire(const std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& wire)
{
  std::vector<Curve2D> curves;
  curves.reserve(wire.size());
  for (const auto& c : wire) {
    if (c.kind == O2BVHSurfaceSolid::PlanarBoundaryCurve::Arc) {
      curves.push_back(Curve2D::makeArc({c.center[0], c.center[1]}, c.radius, c.startAngle, c.endAngle));
    } else if (c.kind == O2BVHSurfaceSolid::PlanarBoundaryCurve::BSpline) {
      std::vector<Vec2> poles;
      poles.reserve(c.poles.size());
      for (const auto& pole : c.poles) {
        poles.push_back({pole[0], pole[1]});
      }
      curves.push_back(Curve2D::makeBSpline(c.degree, std::move(poles), c.weights, c.knots));
    } else {
      curves.push_back(Curve2D::makeLine({c.lineStart[0], c.lineStart[1]}, {c.lineEnd[0], c.lineEnd[1]}));
    }
  }
  return curves;
}
} // namespace

bool O2BVHSurfaceSolid::AddCurvedPlanarSurface(const Point3D& origin, const Point3D& axisU, const Point3D& axisV,
                                               const std::vector<PlanarBoundaryCurve>& outerWire,
                                               const std::vector<std::vector<PlanarBoundaryCurve>>& innerWires)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddCurvedPlanarSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  const std::vector<Curve2D> outerCurves = makeCurveWire(outerWire);
  std::vector<std::vector<Curve2D>> innerCurves;
  innerCurves.reserve(innerWires.size());
  for (const auto& innerWire : innerWires) {
    innerCurves.push_back(makeCurveWire(innerWire));
  }

  auto surface = std::make_unique<CurvedPlanarBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(origin), makeVec3(axisU), makeVec3(axisV), outerCurves, innerCurves,
                           errorMessage)) {
    Error("AddCurvedPlanarSurface", "%s", errorMessage.c_str());
    return false;
  }

  auto record = makeRecord(BVHSurfaceRecord::CurvedPlanar, origin, axisU, axisV, {}, false, false);
  storeCurveWires(record, outerWire, innerWires);
  fRecords.push_back(std::move(record));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddCylindricalSurface(const Point3D& centerPoint, const Point3D& axis,
                                              const Point3D& referenceAxisU, double radius, double heightMin,
                                              double heightMax, double phiStart, double phiSweep, bool innerWall)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddCylindricalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  auto surface = std::make_unique<CylindricalBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(centerPoint), makeVec3(axis), makeVec3(referenceAxisU), radius, heightMin,
                           heightMax, phiStart, phiSweep, innerWall, errorMessage)) {
    Error("AddCylindricalSurface", "%s", errorMessage.c_str());
    return false;
  }

  fRecords.push_back(makeRecord(BVHSurfaceRecord::Cylindrical, centerPoint, axis, referenceAxisU,
                                {radius, heightMin, heightMax, phiStart, phiSweep}, innerWall, false));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddCylindricalSurface(const Point3D& centerPoint, const Point3D& axis,
                                              const Point3D& referenceAxisU, double radius, double heightMin,
                                              double heightMax, double phiStart, double phiSweep, bool innerWall,
                                              const std::vector<PlanarBoundaryCurve>& outerTrim,
                                              const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddCylindricalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  const std::vector<Curve2D> outerCurves = makeCurveWire(outerTrim);
  std::vector<std::vector<Curve2D>> innerCurves;
  innerCurves.reserve(innerTrims.size());
  for (const auto& innerTrim : innerTrims) {
    innerCurves.push_back(makeCurveWire(innerTrim));
  }

  auto surface = std::make_unique<CylindricalBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(centerPoint), makeVec3(axis), makeVec3(referenceAxisU), radius, heightMin,
                           heightMax, phiStart, phiSweep, innerWall, outerCurves, innerCurves, errorMessage)) {
    Error("AddCylindricalSurface", "%s", errorMessage.c_str());
    return false;
  }

  auto record = makeRecord(BVHSurfaceRecord::Cylindrical, centerPoint, axis, referenceAxisU,
                           {radius, heightMin, heightMax, phiStart, phiSweep}, innerWall, true);
  storeCurveWires(record, outerTrim, innerTrims);
  fRecords.push_back(std::move(record));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddSphericalSurface(const Point3D& center, const Point3D& polarAxis,
                                            const Point3D& referenceAxisU, double radius, double thetaMin,
                                            double thetaMax, double phiStart, double phiSweep, bool innerWall)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddSphericalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  auto surface = std::make_unique<SphericalBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(center), makeVec3(polarAxis), makeVec3(referenceAxisU), radius, thetaMin,
                           thetaMax, phiStart, phiSweep, innerWall, errorMessage)) {
    Error("AddSphericalSurface", "%s", errorMessage.c_str());
    return false;
  }

  fRecords.push_back(makeRecord(BVHSurfaceRecord::Spherical, center, polarAxis, referenceAxisU,
                                {radius, thetaMin, thetaMax, phiStart, phiSweep}, innerWall, false));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddSphericalSurface(const Point3D& center, const Point3D& polarAxis,
                                            const Point3D& referenceAxisU, double radius, double thetaMin,
                                            double thetaMax, double phiStart, double phiSweep, bool innerWall,
                                            const std::vector<PlanarBoundaryCurve>& outerTrim,
                                            const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddSphericalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  const std::vector<Curve2D> outerCurves = makeCurveWire(outerTrim);
  std::vector<std::vector<Curve2D>> innerCurves;
  innerCurves.reserve(innerTrims.size());
  for (const auto& innerTrim : innerTrims) {
    innerCurves.push_back(makeCurveWire(innerTrim));
  }

  auto surface = std::make_unique<SphericalBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(center), makeVec3(polarAxis), makeVec3(referenceAxisU), radius, thetaMin,
                           thetaMax, phiStart, phiSweep, innerWall, outerCurves, innerCurves, errorMessage)) {
    Error("AddSphericalSurface", "%s", errorMessage.c_str());
    return false;
  }

  auto record = makeRecord(BVHSurfaceRecord::Spherical, center, polarAxis, referenceAxisU,
                           {radius, thetaMin, thetaMax, phiStart, phiSweep}, innerWall, true);
  storeCurveWires(record, outerTrim, innerTrims);
  fRecords.push_back(std::move(record));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddConicalSurface(const Point3D& centerPoint, const Point3D& axis,
                                          const Point3D& referenceAxisU, double radiusAtMin, double radiusAtMax,
                                          double heightMin, double heightMax, double phiStart, double phiSweep,
                                          bool innerWall)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddConicalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  auto surface = std::make_unique<ConicalBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(centerPoint), makeVec3(axis), makeVec3(referenceAxisU), radiusAtMin,
                           radiusAtMax, heightMin, heightMax, phiStart, phiSweep, innerWall, errorMessage)) {
    Error("AddConicalSurface", "%s", errorMessage.c_str());
    return false;
  }

  fRecords.push_back(makeRecord(BVHSurfaceRecord::Conical, centerPoint, axis, referenceAxisU,
                                {radiusAtMin, radiusAtMax, heightMin, heightMax, phiStart, phiSweep}, innerWall,
                                false));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddConicalSurface(const Point3D& centerPoint, const Point3D& axis,
                                          const Point3D& referenceAxisU, double radiusAtMin, double radiusAtMax,
                                          double heightMin, double heightMax, double phiStart, double phiSweep,
                                          bool innerWall, const std::vector<PlanarBoundaryCurve>& outerTrim,
                                          const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddConicalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  const std::vector<Curve2D> outerCurves = makeCurveWire(outerTrim);
  std::vector<std::vector<Curve2D>> innerCurves;
  innerCurves.reserve(innerTrims.size());
  for (const auto& innerTrim : innerTrims) {
    innerCurves.push_back(makeCurveWire(innerTrim));
  }

  auto surface = std::make_unique<ConicalBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(centerPoint), makeVec3(axis), makeVec3(referenceAxisU), radiusAtMin,
                           radiusAtMax, heightMin, heightMax, phiStart, phiSweep, innerWall, outerCurves,
                           innerCurves, errorMessage)) {
    Error("AddConicalSurface", "%s", errorMessage.c_str());
    return false;
  }

  auto record = makeRecord(BVHSurfaceRecord::Conical, centerPoint, axis, referenceAxisU,
                           {radiusAtMin, radiusAtMax, heightMin, heightMax, phiStart, phiSweep}, innerWall, true);
  storeCurveWires(record, outerTrim, innerTrims);
  fRecords.push_back(std::move(record));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddToroidalSurface(const Point3D& centerPoint, const Point3D& axis,
                                           const Point3D& referenceAxisU, double majorRadius, double minorRadius,
                                           double phiStart, double phiSweep, double tubeStart, double tubeSweep,
                                           bool innerWall)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddToroidalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  auto surface = std::make_unique<TorusBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(centerPoint), makeVec3(axis), makeVec3(referenceAxisU), majorRadius, minorRadius,
                           phiStart, phiSweep, tubeStart, tubeSweep, innerWall, errorMessage)) {
    Error("AddToroidalSurface", "%s", errorMessage.c_str());
    return false;
  }

  fRecords.push_back(makeRecord(BVHSurfaceRecord::Toroidal, centerPoint, axis, referenceAxisU,
                                {majorRadius, minorRadius, phiStart, phiSweep, tubeStart, tubeSweep}, innerWall,
                                false));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

bool O2BVHSurfaceSolid::AddToroidalSurface(const Point3D& centerPoint, const Point3D& axis,
                                           const Point3D& referenceAxisU, double majorRadius, double minorRadius,
                                           double phiStart, double phiSweep, double tubeStart, double tubeSweep,
                                           bool innerWall, const std::vector<PlanarBoundaryCurve>& outerTrim,
                                           const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (fImpl->defined) {
    Error("AddToroidalSurface", "Shape %s already fully defined. Not adding", GetName());
    return false;
  }

  const std::vector<Curve2D> outerCurves = makeCurveWire(outerTrim);
  std::vector<std::vector<Curve2D>> innerCurves;
  innerCurves.reserve(innerTrims.size());
  for (const auto& innerTrim : innerTrims) {
    innerCurves.push_back(makeCurveWire(innerTrim));
  }

  auto surface = std::make_unique<TorusBoundedSurface>();
  std::string errorMessage;
  if (!surface->initialize(makeVec3(centerPoint), makeVec3(axis), makeVec3(referenceAxisU), majorRadius, minorRadius,
                           phiStart, phiSweep, tubeStart, tubeSweep, innerWall, outerCurves, innerCurves,
                           errorMessage)) {
    Error("AddToroidalSurface", "%s", errorMessage.c_str());
    return false;
  }

  auto record = makeRecord(BVHSurfaceRecord::Toroidal, centerPoint, axis, referenceAxisU,
                           {majorRadius, minorRadius, phiStart, phiSweep, tubeStart, tubeSweep}, innerWall, true);
  storeCurveWires(record, outerTrim, innerTrims);
  fRecords.push_back(std::move(record));

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  return true;
}

void O2BVHSurfaceSolid::CloseShape(bool check)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  // An empty surface set is not a degenerate solid, it is an *unknown* one, and the difference is
  // the whole point of NavigationReliability. Every closure counter of an empty set is zero, so a
  // default ClosureReport calls it a closed, consistently oriented 2-manifold and the shape would
  // answer Reliable -- with full confidence and no surfaces to be confident about. Leave it
  // undefined (Undetermined) instead, and leave the bounding box alone: for a solid arriving from
  // a ROOT file the streamed box is better information than a zeroed one. This guard is
  // unconditional; \a check selects whether defects are *reported*, not whether they are honest.
  // See CodeReview_Fable.md S1.
  if (fImpl->surfaces.empty()) {
    Error("CloseShape", "Shape %s has no bounded surfaces; it stays undefined and reports itself not navigable",
          GetName());
    return;
  }

  ComputeBBox();
  fImpl->buildBVH();

  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  fImpl->displayTriangleSurface.clear();
  for (size_t surfaceIndex = 0; surfaceIndex < fImpl->surfaces.size(); ++surfaceIndex) {
    fImpl->surfaces[surfaceIndex]->appendDisplayMesh(fImpl->displayVertices, fImpl->displayTriangles);
    // resize() only writes the entries it adds, so each surface stamps exactly its own triangles.
    fImpl->displayTriangleSurface.resize(fImpl->displayTriangles.size(), static_cast<int>(surfaceIndex));
  }

  fImpl->closure = validateClosure(fImpl->surfaces, fModelTolerance);
  fImpl->rimReports.clear();
  fImpl->rimReports.reserve(fImpl->closure.rimRecords.size());
  for (const RimRecord& record : fImpl->closure.rimRecords) {
    RimReport report;
    report.surface = record.surfaceIndex;
    report.rimOnSurface = record.rimIndexOnSurface;
    report.closed = record.closed;
    report.chords = record.chords;
    report.unmatchedChords = record.unmatchedChords;
    report.length = record.length;
    report.unmatchedLength = record.unmatchedLength;
    report.maxIsolation = record.maxIsolation;
    report.maxIsolationPoint = {record.maxIsolationPoint.xCoord, record.maxIsolationPoint.yCoord,
                                record.maxIsolationPoint.zCoord};
    report.maxIsolationFace = record.maxIsolationFace;
    report.state = rimStateToReliability(record.state);
    fImpl->rimReports.push_back(report);
  }
  fImpl->defined = true;

  if (check) {
    const auto& closure = fImpl->closure;
    // State the *consequence*, not only the counts. A warning that reads as a quality metric gets
    // read as one: this defect went unnoticed for two sessions because "699 boundary edge(s)"
    // looked like a cosmetic mesh remark while the solid was answering Contains wrongly by
    // centimetres. See scripts/geometry/ExactTrimTopology.md.
    if (closure.edgeIdentityAvailable && closure.boundaryRims > 0) {
      // Identity, not proximity: this counts source edges with only one face on them, so it says
      // a face is missing rather than that two polylines drifted apart.
      Error("CloseShape",
            "Shape %s is NOT a closed surface: %d of its %d source edge(s) are used by only one face and %d by more "
            "than two, leaving %d of %d trim loop(s) without a neighbour. This is a count of edge identities, not a "
            "proximity test, so no tolerance can make it go away. NAVIGATION IS UNRELIABLE: parity containment is "
            "undefined in the shadow of every gap, so Contains()/DistFrom*() can be wrong arbitrarily far from any "
            "surface. Check IsNavigable()/GetNavigationReliability(); GetRimReports() names the offending loops.",
            GetName(), closure.edgeBoundaryCount, closure.edgeIncidences, closure.edgeNonManifoldCount,
            closure.boundaryRims, closure.rims);
    } else if (closure.boundaryRims > 0) {
      Error("CloseShape",
            "Shape %s is NOT a closed surface: %d of %d trim loop(s) run out of neighbouring face (declared "
            "tolerance %g cm, widened per chord by the rim sampling), leaving %g cm of its %g cm of boundary open; "
            "the most isolated open chord is %g cm from any other face. NAVIGATION IS UNRELIABLE: "
            "parity containment is undefined in the shadow of every gap, so Contains()/DistFrom*() can be wrong "
            "arbitrarily far from any surface, not just near the gap. Check IsNavigable()/"
            "GetNavigationReliability() before trusting any answer from this shape; GetRimReports() names the "
            "offending loops.",
            GetName(), closure.boundaryRims, closure.rims, closure.rimEpsilon, closure.unmatchedRimLength,
            closure.totalRimLength, closure.maxRimIsolation);
    }
    if (closure.nonManifoldRims > 0) {
      Error("CloseShape",
            "Shape %s is NOT a 2-manifold: %d of %d trim loop(s) run along two or more other faces at once "
            "(coincident or duplicated faces). NAVIGATION IS UNRELIABLE: parity containment depends on the order "
            "in which coincident hits are clustered, so Contains() is not even self-consistent between traversal "
            "orders. Check IsNavigable()/GetNavigationReliability() before trusting any answer from this shape.",
            GetName(), closure.nonManifoldRims, closure.rims);
    }
    if (!closure.orientationConsistent) {
      Error("CloseShape",
            "Shape %s has %d inconsistently oriented (reversed) trim loop(s): adjacent faces disagree about "
            "which side is out. NAVIGATION IS UNRELIABLE: the entering/exiting sign of a crossing is taken from the "
            "surface normal, so distance queries can return the wrong side. Check IsNavigable()/"
            "GetNavigationReliability() before trusting any answer from this shape.",
            GetName(), closure.reversedRims);
    }
    if (closure.closed && closure.signedVolume < 0.) {
      Warning("CloseShape",
              "Shape %s has inward-pointing surface normals (signed volume %g); navigation expects outward normals",
              GetName(), closure.signedVolume);
    }
  }
}

int O2BVHSurfaceSolid::GetNsurfaces() const
{
  return fImpl == nullptr ? 0 : static_cast<int>(fImpl->surfaces.size());
}

bool O2BVHSurfaceSolid::IsDefined() const
{
  return fImpl != nullptr && fImpl->defined;
}

void O2BVHSurfaceSolid::SetModelTolerance(double toleranceCm)
{
  if (!(toleranceCm >= 0.) || !std::isfinite(toleranceCm)) {
    Error("SetModelTolerance", "Shape %s: ignoring a non-finite or negative model tolerance %g; it stays %g",
          GetName(), toleranceCm, fModelTolerance);
    return;
  }
  fModelTolerance = toleranceCm;
}

bool O2BVHSurfaceSolid::HasBVH() const
{
  return fImpl != nullptr && fImpl->bvh != nullptr;
}

bool O2BVHSurfaceSolid::GetBVHRootBounds(Point3D& lower, Point3D& upper) const
{
  if (!HasBVH()) {
    return false;
  }
  const auto rootBox = fImpl->bvh->get_root().get_bbox();
  for (int dimension = 0; dimension < 3; ++dimension) {
    lower[dimension] = rootBox.min[dimension];
    upper[dimension] = rootBox.max[dimension];
  }
  return true;
}

int O2BVHSurfaceSolid::CountBVHRayCandidates(const Point3D& point, const Point3D& direction) const
{
  if (!HasBVH()) {
    return -1;
  }
  int candidates = 0;
  fImpl->visitRayCandidates(makeVec3(point), makeVec3(direction), [&](const BoundedSurface&) { ++candidates; });
  return candidates;
}

bool O2BVHSurfaceSolid::IsClosed() const
{
  return fImpl != nullptr && fImpl->defined && fImpl->closure.closed;
}

bool O2BVHSurfaceSolid::IsOrientationConsistent() const
{
  return fImpl != nullptr && fImpl->defined && fImpl->closure.orientationConsistent;
}

O2BVHSurfaceSolid::NavigationReliability O2BVHSurfaceSolid::GetNavigationReliability() const
{
  if (fImpl == nullptr || !fImpl->defined) {
    return NavigationReliability::Undetermined;
  }
  // worst defect wins; the enum is ordered by severity. The counts are per *rim*: a face's whole
  // trim loop, matched against the other faces as a curve. The chord counters this used to read
  // ask whether two faces emitted the same vertices along a shared edge, which they have no
  // reason to and, on real CAD, do not.
  const auto& closure = fImpl->closure;
  // When the faces state their edge identities, the question is decided by counting them and not
  // by any geometric test at all. Read those counters directly rather than the per-rim projection
  // of them: a face with no trim loop of its own (a full sphere, say) contributes no rim record,
  // so a defect on its edges would have nowhere to be reported.
  if (closure.edgeIdentityAvailable) {
    if (closure.edgeNonManifoldCount > 0) {
      return NavigationReliability::NonManifold;
    }
    if (closure.edgeBoundaryCount > 0) {
      return NavigationReliability::OpenSurfaceSet;
    }
    if (closure.edgeReversedCount > 0) {
      return NavigationReliability::ReversedFaces;
    }
    return NavigationReliability::Reliable;
  }
  if (closure.nonManifoldRims > 0) {
    return NavigationReliability::NonManifold;
  }
  if (closure.boundaryRims > 0) {
    return NavigationReliability::OpenSurfaceSet;
  }
  if (closure.reversedRims > 0) {
    return NavigationReliability::ReversedFaces;
  }
  return NavigationReliability::Reliable;
}

bool O2BVHSurfaceSolid::IsNavigable() const
{
  return GetNavigationReliability() == NavigationReliability::Reliable;
}

const char* O2BVHSurfaceSolid::GetNavigationReliabilityName(NavigationReliability reliability)
{
  switch (reliability) {
    case NavigationReliability::Undetermined:
      return "undetermined";
    case NavigationReliability::Reliable:
      return "reliable";
    case NavigationReliability::ReversedFaces:
      return "reversed-faces";
    case NavigationReliability::OpenSurfaceSet:
      return "open-surface-set";
    case NavigationReliability::NonManifold:
      return "non-manifold";
  }
  return "unknown";
}

int O2BVHSurfaceSolid::GetBoundaryEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.boundaryEdges;
}

int O2BVHSurfaceSolid::GetNonManifoldEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.nonManifoldEdges;
}

int O2BVHSurfaceSolid::GetReversedEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.reversedEdges;
}

double O2BVHSurfaceSolid::GetMaxRimIsolation() const
{
  return fImpl == nullptr ? 0. : fImpl->closure.maxRimIsolation;
}

bool O2BVHSurfaceSolid::SetSurfaceBoundaryEdges(int surfaceIndex, const std::vector<unsigned int>& edgeIds,
                                                const std::vector<unsigned char>& edgeFlags)
{
  if (fImpl == nullptr || surfaceIndex < 0 || surfaceIndex >= static_cast<int>(fImpl->surfaces.size()) ||
      surfaceIndex >= static_cast<int>(fRecords.size())) {
    Error("SetSurfaceBoundaryEdges", "Shape %s: surface index %d is out of range (%d surface(s))", GetName(),
          surfaceIndex, GetNsurfaces());
    return false;
  }
  if (edgeIds.size() != edgeFlags.size()) {
    Error("SetSurfaceBoundaryEdges", "Shape %s: surface %d was given %d edge id(s) and %d flag(s)", GetName(),
          surfaceIndex, static_cast<int>(edgeIds.size()), static_cast<int>(edgeFlags.size()));
    return false;
  }
  std::vector<BoundedSurface::BoundaryEdgeRef> refs;
  refs.reserve(edgeIds.size());
  for (size_t index = 0; index < edgeIds.size(); ++index) {
    BoundedSurface::BoundaryEdgeRef ref;
    ref.edgeId = edgeIds[index];
    ref.reversed = (edgeFlags[index] & kEdgeReversed) != 0;
    ref.degenerate = (edgeFlags[index] & kEdgeDegenerate) != 0;
    ref.anchored = (edgeFlags[index] & kEdgeAnchored) != 0;
    refs.push_back(ref);
  }
  fImpl->surfaces[static_cast<size_t>(surfaceIndex)]->setBoundaryEdges(std::move(refs));
  fRecords[static_cast<size_t>(surfaceIndex)].boundaryEdgeIds = edgeIds;
  fRecords[static_cast<size_t>(surfaceIndex)].boundaryEdgeFlags = edgeFlags;
  return true;
}

bool O2BVHSurfaceSolid::HasEdgeIdentity() const
{
  return fImpl != nullptr && fImpl->closure.edgeIdentityAvailable;
}

int O2BVHSurfaceSolid::GetSourceEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.edgeIncidences;
}

int O2BVHSurfaceSolid::GetSharedSourceEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.edgeSharedCount;
}

int O2BVHSurfaceSolid::GetBoundarySourceEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.edgeBoundaryCount;
}

int O2BVHSurfaceSolid::GetNonManifoldSourceEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.edgeNonManifoldCount;
}

int O2BVHSurfaceSolid::GetReversedSourceEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.edgeReversedCount;
}

int O2BVHSurfaceSolid::GetDegenerateSourceEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.edgeDegenerateCount;
}

double O2BVHSurfaceSolid::GetMaxSharedEdgeDeviation() const
{
  return fImpl == nullptr ? 0. : fImpl->closure.maxSharedEdgeDeviation;
}

int O2BVHSurfaceSolid::GetMeasuredSharedEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.sharedEdgesMeasured;
}

int O2BVHSurfaceSolid::GetUnmeasuredSharedEdgeCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.sharedEdgesUnmeasured;
}

double O2BVHSurfaceSolid::GetRimChordResolution() const
{
  return fImpl == nullptr ? 0. : fImpl->closure.rimChordResolution;
}

double O2BVHSurfaceSolid::GetRimMatchTolerance() const
{
  return fImpl == nullptr ? 0. : fImpl->closure.rimEpsilon;
}

double O2BVHSurfaceSolid::GetTotalRimLength() const
{
  return fImpl == nullptr ? 0. : fImpl->closure.totalRimLength;
}

double O2BVHSurfaceSolid::GetUnmatchedRimLength() const
{
  return fImpl == nullptr ? 0. : fImpl->closure.unmatchedRimLength;
}

int O2BVHSurfaceSolid::GetRimCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.rims;
}

int O2BVHSurfaceSolid::GetMatchedRimCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.matchedRims;
}

int O2BVHSurfaceSolid::GetBoundaryRimCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.boundaryRims;
}

int O2BVHSurfaceSolid::GetNonManifoldRimCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.nonManifoldRims;
}

int O2BVHSurfaceSolid::GetReversedRimCount() const
{
  return fImpl == nullptr ? 0 : fImpl->closure.reversedRims;
}

const std::vector<O2BVHSurfaceSolid::RimReport>& O2BVHSurfaceSolid::GetRimReports() const
{
  static const std::vector<RimReport> empty;
  return fImpl == nullptr ? empty : fImpl->rimReports;
}

void O2BVHSurfaceSolid::GetSurfaceCapacityContributions(std::vector<double>& contributions) const
{
  contributions.clear();
  if (fImpl == nullptr) {
    return;
  }
  contributions.reserve(fImpl->surfaces.size());
  for (const auto& surface : fImpl->surfaces) {
    contributions.push_back(surface == nullptr ? 0. : surface->capacityContribution());
  }
}

void O2BVHSurfaceSolid::ComputeBBox()
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    fDX = fDY = fDZ = 0.;
    fOrigin[0] = fOrigin[1] = fOrigin[2] = 0.;
    return;
  }

  Vec3 lowerCorner{TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};
  Vec3 upperCorner{-TGeoShape::Big(), -TGeoShape::Big(), -TGeoShape::Big()};
  for (const auto& surface : fImpl->surfaces) {
    surface->conservativeBounds(lowerCorner, upperCorner);
  }

  for (int dimension = 0; dimension < 3; ++dimension) {
    const double lowerValue = component(lowerCorner, dimension) - kTolerance;
    const double upperValue = component(upperCorner, dimension) + kTolerance;
    fOrigin[dimension] = 0.5 * (lowerValue + upperValue);
    const double halfLength = 0.5 * (upperValue - lowerValue);
    if (dimension == 0) {
      fDX = halfLength;
    } else if (dimension == 1) {
      fDY = halfLength;
    } else {
      fDZ = halfLength;
    }
  }
}

void O2BVHSurfaceSolid::GetMeshNumbers(int& nvert, int& nsegs, int& npols) const
{
  nvert = GetNmeshVertices();
  npols = fImpl == nullptr ? 0 : static_cast<int>(fImpl->displayTriangles.size());
  nsegs = 3 * npols;
}

int O2BVHSurfaceSolid::GetNmeshVertices() const
{
  return fImpl == nullptr ? 0 : static_cast<int>(fImpl->displayVertices.size());
}

namespace
{
/// A deterministic, low-discrepancy pair in [0,1)^2 -- the R2 (generalised golden-ratio) sequence.
/// Deterministic matters more here than it usually does: the point set a shape hands out decides
/// what an overlap check sees, so two runs that disagree must be a geometry difference and never a
/// seed difference. Low-discrepancy matters because a triangle receives only a handful of samples
/// and a pseudo-random pair clumps at that count.
void r2Pair(long long index, double& firstCoordinate, double& secondCoordinate)
{
  constexpr double kAlpha1 = 0.7548776662466927; // 1 / plastic number
  constexpr double kAlpha2 = 0.5698402909980532; // 1 / plastic number^2
  const double shifted = static_cast<double>(index + 1);
  firstCoordinate = std::fmod(0.5 + kAlpha1 * shifted, 1.);
  secondCoordinate = std::fmod(0.5 + kAlpha2 * shifted, 1.);
}
} // namespace

/// Newton on the patch distance along the patch normal. Two steps suffice on every family in the
/// corpus (a plane converges in one, exactly), but the loop is cheap and the *verdict* is what is
/// returned: a caller must not emit a point this could not put on the surface.
bool O2BVHSurfaceSolid::ProjectOntoPatch(int surfaceIndex, double* point) const
{
  if (fImpl == nullptr || surfaceIndex < 0 || static_cast<size_t>(surfaceIndex) >= fImpl->surfaces.size()) {
    return false;
  }
  const BoundedSurface& surface = *fImpl->surfaces[surfaceIndex];
  constexpr double kToleranceSquared = kSurfacePointTolerance * kSurfacePointTolerance;

  Vec3 current = makeVec3(point);
  double currentDistanceSq = surface.distanceSqToPatch(current);
  for (int iteration = 0; iteration < 8 && currentDistanceSq > kToleranceSquared; ++iteration) {
    const double distance = std::sqrt(currentDistanceSq);
    const Vec3 normal = surface.normalAt(current);
    const Vec3 inward{current.xCoord - distance * normal.xCoord, current.yCoord - distance * normal.yCoord,
                      current.zCoord - distance * normal.zCoord};
    const Vec3 outward{current.xCoord + distance * normal.xCoord, current.yCoord + distance * normal.yCoord,
                       current.zCoord + distance * normal.zCoord};
    const double inwardDistanceSq = surface.distanceSqToPatch(inward);
    const double outwardDistanceSq = surface.distanceSqToPatch(outward);
    const double bestDistanceSq = std::min(inwardDistanceSq, outwardDistanceSq);
    // Not converging: the nearest point of the patch is not along its own normal from here, which
    // happens when the sample falls outside the trim and the nearest point is on the trim wire.
    // Say so rather than emit whatever the last iterate happened to be.
    if (!(bestDistanceSq < currentDistanceSq)) {
      return false;
    }
    current = (inwardDistanceSq < outwardDistanceSq) ? inward : outward;
    currentDistanceSq = bestDistanceSq;
  }

  if (currentDistanceSq > kToleranceSquared) {
    return false;
  }
  point[0] = current.xCoord;
  point[1] = current.yCoord;
  point[2] = current.zCoord;
  return true;
}

Bool_t O2BVHSurfaceSolid::GetPointsOnSegments(Int_t npoints, Double_t* array) const
{
  if (array == nullptr || npoints <= 0 || fImpl == nullptr || fImpl->displayVertices.empty()) {
    return kFALSE;
  }
  const int vertexCount = static_cast<int>(fImpl->displayVertices.size());
  // Below the mesh size, hand the caller back to SetPoints(): it gets *more* points than it asked
  // for and every one of them is a patch sample rather than something this method invented. The
  // base class returns kFALSE here too; it just shouts about it first, on every shape, every run.
  if (npoints < vertexCount) {
    return kFALSE;
  }

  auto writeVertex = [&](int slot, const Vec3& vertex) {
    array[3 * slot + 0] = vertex.xCoord;
    array[3 * slot + 1] = vertex.yCoord;
    array[3 * slot + 2] = vertex.zCoord;
  };

  for (int vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex) {
    writeVertex(vertexIndex, fImpl->displayVertices[vertexIndex]);
  }

  const int extraCount = npoints - vertexCount;
  const int triangleCount = static_cast<int>(fImpl->displayTriangles.size());
  if (extraCount == 0) {
    return kTRUE;
  }
  if (triangleCount == 0 || fImpl->displayTriangleSurface.size() != fImpl->displayTriangles.size()) {
    // No triangles (or a mesh built before the provenance existed): repeat vertices rather than
    // leave the tail of the buffer uninitialised, which the caller would read as coordinates.
    for (int extraIndex = 0; extraIndex < extraCount; ++extraIndex) {
      writeVertex(vertexCount + extraIndex, fImpl->displayVertices[extraIndex % vertexCount]);
    }
    return kTRUE;
  }

  for (int extraIndex = 0; extraIndex < extraCount; ++extraIndex) {
    // Stride over the triangles rather than walking them in order, so a request that cannot cover
    // every triangle still spreads over the whole solid instead of over its first few faces.
    const int triangleIndex =
      static_cast<int>((static_cast<long long>(extraIndex) * triangleCount) / extraCount) % triangleCount;
    const auto& triangle = fImpl->displayTriangles[triangleIndex];
    const Vec3& cornerA = fImpl->displayVertices[triangle[0]];
    const Vec3& cornerB = fImpl->displayVertices[triangle[1]];
    const Vec3& cornerC = fImpl->displayVertices[triangle[2]];

    double firstCoordinate = 0.;
    double secondCoordinate = 0.;
    r2Pair(extraIndex, firstCoordinate, secondCoordinate);
    if (firstCoordinate + secondCoordinate > 1.) {
      firstCoordinate = 1. - firstCoordinate;
      secondCoordinate = 1. - secondCoordinate;
    }
    const double weightA = 1. - firstCoordinate - secondCoordinate;
    double candidate[3] = {weightA * cornerA.xCoord + firstCoordinate * cornerB.xCoord + secondCoordinate * cornerC.xCoord,
                           weightA * cornerA.yCoord + firstCoordinate * cornerB.yCoord + secondCoordinate * cornerC.yCoord,
                           weightA * cornerA.zCoord + firstCoordinate * cornerB.zCoord + secondCoordinate * cornerC.zCoord};

    if (!ProjectOntoPatch(fImpl->displayTriangleSurface[triangleIndex], candidate)) {
      // The honest fallback is a vertex of the very triangle the sample came from: it is on the
      // patch by construction, so the point set stays exactly on the solid even where the
      // projection gives up.
      candidate[0] = cornerA.xCoord;
      candidate[1] = cornerA.yCoord;
      candidate[2] = cornerA.zCoord;
    }
    array[3 * (vertexCount + extraIndex) + 0] = candidate[0];
    array[3 * (vertexCount + extraIndex) + 1] = candidate[1];
    array[3 * (vertexCount + extraIndex) + 2] = candidate[2];
  }
  return kTRUE;
}

TBuffer3D* O2BVHSurfaceSolid::MakeBuffer3D() const
{
  int nvert = 0;
  int nsegs = 0;
  int npols = 0;
  GetMeshNumbers(nvert, nsegs, npols);
  auto buff = new TBuffer3D(TBuffer3DTypes::kGeneric, nvert, 3 * nvert, nsegs, 3 * nsegs, npols, 6 * npols);
  if (buff != nullptr) {
    SetPoints(buff->fPnts);
    SetSegsAndPols(*buff);
  }
  return buff;
}

void O2BVHSurfaceSolid::Print(Option_t*) const
{
  std::cout << "=== BVH surface solid " << GetName() << " having " << GetNsurfaces() << " bounded surfaces\n";
  const auto reliability = GetNavigationReliability();
  std::cout << "    navigation: " << GetNavigationReliabilityName(reliability);
  if (reliability != NavigationReliability::Reliable && reliability != NavigationReliability::Undetermined) {
    std::cout << " (UNRELIABLE; boundary=" << GetBoundaryEdgeCount() << " non-manifold=" << GetNonManifoldEdgeCount()
              << " reversed=" << GetReversedEdgeCount() << ")";
  }
  std::cout << "\n    model tolerance: ";
  if (fModelTolerance > 0.) {
    std::cout << fModelTolerance << " cm (from the source model)";
  } else {
    std::cout << "not stated";
  }
  // The identity counts first, because when they are present they are the verdict and the rim
  // numbers below are measurement. maxSharedEdgeDeviation is the one number that does mean "how
  // far apart are these two faces", and it is printed next to the model's own tolerance so the
  // reader can see at once whether the surfaces meet as well as the CAD system claims they do.
  if (HasEdgeIdentity()) {
    std::cout << "\n    edge identity: " << GetSourceEdgeCount() << " source edge(s), shared=" << GetSharedSourceEdgeCount()
              << " boundary=" << GetBoundarySourceEdgeCount() << " non-manifold=" << GetNonManifoldSourceEdgeCount()
              << " reversed=" << GetReversedSourceEdgeCount() << " degenerate=" << GetDegenerateSourceEdgeCount()
              << "\n    shared edge deviation: max " << GetMaxSharedEdgeDeviation() << " cm over "
              << GetMeasuredSharedEdgeCount() << " measured edge(s)";
    if (GetUnmeasuredSharedEdgeCount() > 0) {
      std::cout << " (" << GetUnmeasuredSharedEdgeCount() << " not measurable: parametric-rectangle trim)";
    }
  }
  // The isolation, and the resolution that widened the band it was judged in, always together: the
  // number is how alone the loneliest chord is, not how far apart two faces are.
  if (GetRimCount() > 0) {
    std::cout << "\n    rim isolation: max " << GetMaxRimIsolation() << " cm (chord resolution "
              << GetRimChordResolution() << " cm, declared tolerance " << GetRimMatchTolerance() << " cm)"
              << "\n    rims: " << GetRimCount() << " (matched=" << GetMatchedRimCount()
              << " boundary=" << GetBoundaryRimCount() << " non-manifold=" << GetNonManifoldRimCount()
              << " reversed=" << GetReversedRimCount() << "), open " << GetUnmatchedRimLength() << " of "
              << GetTotalRimLength() << " cm";
  }
  std::cout << "\n";
}

void O2BVHSurfaceSolid::SetPoints(double* points) const
{
  if (fImpl == nullptr) {
    return;
  }
  int coordinateIndex = 0;
  for (const auto& vertex : fImpl->displayVertices) {
    points[coordinateIndex++] = vertex.xCoord;
    points[coordinateIndex++] = vertex.yCoord;
    points[coordinateIndex++] = vertex.zCoord;
  }
}

void O2BVHSurfaceSolid::SetPoints(Float_t* points) const
{
  if (fImpl == nullptr) {
    return;
  }
  int coordinateIndex = 0;
  for (const auto& vertex : fImpl->displayVertices) {
    points[coordinateIndex++] = vertex.xCoord;
    points[coordinateIndex++] = vertex.yCoord;
    points[coordinateIndex++] = vertex.zCoord;
  }
}

void O2BVHSurfaceSolid::SetSegsAndPols(TBuffer3D& buff) const
{
  if (fImpl == nullptr) {
    return;
  }

  const int color = GetBasicColor();
  int* segs = buff.fSegs;
  int* pols = buff.fPols;
  int segmentDataIndex = 0;
  int polygonDataIndex = 0;
  int segmentIndex = 0;
  for (const auto& triangle : fImpl->displayTriangles) {
    pols[polygonDataIndex++] = color;
    pols[polygonDataIndex++] = 3;
    for (int triangleEdge = 0; triangleEdge < 3; ++triangleEdge) {
      const int nextTriangleEdge = (triangleEdge + 1) % 3;
      segs[segmentDataIndex++] = color;
      segs[segmentDataIndex++] = triangle[triangleEdge];
      segs[segmentDataIndex++] = triangle[nextTriangleEdge];
      pols[polygonDataIndex + 2 - triangleEdge] = segmentIndex++;
    }
    polygonDataIndex += 3;
  }
}

const TBuffer3D& O2BVHSurfaceSolid::GetBuffer3D(int reqSections, Bool_t localFrame) const
{
  static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  FillBuffer3D(buffer, reqSections, localFrame);

  int nvert = 0;
  int nsegs = 0;
  int npols = 0;
  GetMeshNumbers(nvert, nsegs, npols);

  if (reqSections & TBuffer3D::kRawSizes) {
    if (buffer.SetRawSizes(nvert, 3 * nvert, nsegs, 3 * nsegs, npols, 6 * npols)) {
      buffer.SetSectionsValid(TBuffer3D::kRawSizes);
    }
  }
  if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
    SetPoints(buffer.fPnts);
    if (!buffer.fLocalFrame) {
      TransformPoints(buffer.fPnts, buffer.NbPnts());
    }
    SetSegsAndPols(buffer);
    buffer.SetSectionsValid(TBuffer3D::kRaw);
  }

  return buffer;
}

bool O2BVHSurfaceSolid::Contains(const Double_t* point) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return false;
  }

  if (fImpl->bvh == nullptr) {
    // Before CloseShape there is no acceleration structure *and* no bounding box: fDX/fDY/fDZ and
    // fOrigin are all still zero. So this fallback has to come first -- with the box pre-check
    // ahead of it, every point further than 1e-9 from the origin was rejected outright and the
    // documented fallback could never run, which disabled Contains() on any solid that had not
    // been closed yet. See CodeReview_Fable.md S2.
    return Contains_Loop(point);
  }

  const Vec3 testPoint = makeVec3(point);
  if (std::abs(testPoint.xCoord - fOrigin[0]) > fDX + kTolerance ||
      std::abs(testPoint.yCoord - fOrigin[1]) > fDY + kTolerance ||
      std::abs(testPoint.zCoord - fOrigin[2]) > fDZ + kTolerance) {
    return false;
  }

  // boundary policy: a point within tolerance of any surface patch counts as inside
  if (fImpl->visitPointCandidates(
        testPoint, [&](const BoundedSurface& surface) { return surface.containsPointOnSurface(testPoint); })) {
    return true;
  }

  return containsByParity(point, true);
}

bool O2BVHSurfaceSolid::ContainsAlongDirection(const Double_t* point, const Double_t* direction) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return false;
  }
  const Vec3 testPoint = makeVec3(point);
  for (const auto& surface : fImpl->surfaces) {
    if (surface->containsPointOnSurface(testPoint)) {
      return true;
    }
  }
  return fImpl->parityAlong(testPoint, normalized(makeVec3(direction)), fImpl->bvh != nullptr);
}

bool O2BVHSurfaceSolid::Contains_Loop(const Double_t* point) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return false;
  }

  const Vec3 testPoint = makeVec3(point);
  for (const auto& surface : fImpl->surfaces) {
    if (surface->containsPointOnSurface(testPoint)) {
      return true;
    }
  }

  return containsByParity(point, false);
}

bool O2BVHSurfaceSolid::containsByParity(const Double_t* point, bool useBVH) const
{
  // On a closed, consistently oriented 2-manifold parity is a topological invariant, so the answer
  // cannot depend on where the ray is aimed and one shot is the whole story. That is not a
  // convenient assumption but a measured one, and the measurement was re-run when the closure
  // criterion became rim-based, because a criterion that succeeds more often is a criterion that
  // moves solids onto this path (TolerancePolicy.md sections 3.4 and 10).
  //
  // Over both corpora -- 9 converted fixtures and 12 Bagger parts, 11k points and 13 spread
  // directions each -- the 13 Reliable parts disagree between directions at exactly **one** point
  // in 143000, on cyl_cross_cyl, where 12 of the 13 directions (including this one) agree with
  // each other and the point sits 0.1 cm from the nearest surface. That is one unlucky ray near
  // the crossing curve of two cylinders, not the shadow of a gap: a gap shadow costs the point
  // most of its directions and lies behind the gap. The parts that do disagree between directions
  // (0.55%-7% of their points) are, as before, exactly the parts the closure check rejects.
  //
  // So the fast path costs nothing where it is valid, and the vote is paid exactly where the
  // geometry has already admitted it is broken.
  //
  // That single point has since been diagnosed, and it is the reason for the re-shoot below
  // (TolerancePolicy.md section 11). It is not a root-finding defect: the ray gains a third
  // crossing on a patch that survives 6e-6 cm past the true Boolean seam, because a hit inside
  // CurveWire::boundaryBand is resolved as "inside the trim" and that tie-break is one-sided.
  // So the trigger is the tie-break itself -- RayHit::onTrimBoundary -- and not, as was tried and
  // measured before it, the coincident-hit cancellation rule. That earlier trigger cost 0.3-2% of
  // Contains on most parts and 16% on box_minus_cyl and moved none of the 143000 sweep points,
  // because it fires where the cancellation rule is already right (grazes at convex and concave
  // edges, rays through a shared edge) and stays silent here. This one fires on 1 hit in 96291,
  // which is the whole of its cost.
  const Vec3 testPoint = makeVec3(point);
  if (GetNavigationReliability() == NavigationReliability::Reliable) {
    bool ambiguous = false;
    const bool answer = fImpl->parityAlong(testPoint, kContainsTestDirection, useBVH, &ambiguous);
    if (!ambiguous) {
      return answer;
    }
    // This shot crossed a patch within its own trim accuracy, so its parity rests on a tie-break
    // rather than on the geometry. Re-aim: the sliver belongs to the ray, not to the point.
    return fImpl->containsByVote(testPoint, useBVH);
  }
  return fImpl->containsByVote(testPoint, useBVH);
}

void O2BVHSurfaceSolid::DescribeContainsCrossings(const Point3D& point,
                                                  std::vector<ContainsCrossing>& bvhCrossings,
                                                  std::vector<ContainsCrossing>& loopCrossings) const
{
  const Point3D direction{kContainsTestDirection.xCoord, kContainsTestDirection.yCoord,
                          kContainsTestDirection.zCoord};
  DescribeContainsCrossings(point, direction, bvhCrossings, loopCrossings);
}

void O2BVHSurfaceSolid::DescribeContainsCrossings(const Point3D& point, const Point3D& direction,
                                                  std::vector<ContainsCrossing>& bvhCrossings,
                                                  std::vector<ContainsCrossing>& loopCrossings) const
{
  bvhCrossings.clear();
  loopCrossings.clear();
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return;
  }
  const Vec3 testPoint = makeVec3(point.data());
  const Vec3 testDirection = normalized(makeVec3(direction.data()));

  auto collect = [&](std::vector<RayHit>& hits, std::vector<ContainsCrossing>& out) {
    std::sort(hits.begin(), hits.end(),
              [](const RayHit& first, const RayHit& second) { return first.distance < second.distance; });
    out.reserve(hits.size());
    for (const auto& hit : hits) {
      out.push_back({hit.distance, dot(hit.normal, testDirection), hit.onTrimBoundary});
    }
  };

  std::vector<RayHit> loopHits;
  for (const auto& surface : fImpl->surfaces) {
    surface->appendIntersections(testPoint, testDirection, kRayTolerance, TGeoShape::Big(), loopHits);
  }
  collect(loopHits, loopCrossings);

  if (fImpl->bvh != nullptr) {
    std::vector<RayHit> bvhHits;
    fImpl->visitRayCandidates(testPoint, testDirection, [&](const BoundedSurface& surface) {
      surface.appendIntersections(testPoint, testDirection, kRayTolerance, TGeoShape::Big(), bvhHits);
    });
    collect(bvhHits, bvhCrossings);
  }
}

Double_t O2BVHSurfaceSolid::DistFromOutside(const Double_t* point, const Double_t* dir, Int_t, Double_t stepmax,
                                           Double_t* safe) const
{
  if (safe != nullptr) {
    *safe = Safety(point, kFALSE);
  }
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }
  if (fImpl->bvh == nullptr) {
    // before CloseShape there is no acceleration structure yet; stay usable via the plain loop
    return DistFromOutside_Loop(point, dir, stepmax);
  }

  // Cheap reject before touching the BVH: along any single axis the gap between the point and
  // the bounding box is a lower bound on the distance to the solid, so a gap beyond stepmax
  // means no reachable crossing. The bounding box is the shape's own, hence the tolerance slack.
  const Double_t halfLengths[3] = {fDX, fDY, fDZ};
  for (int dimension = 0; dimension < 3; ++dimension) {
    const Double_t lower = fOrigin[dimension] - halfLengths[dimension];
    const Double_t upper = fOrigin[dimension] + halfLengths[dimension];
    if (lower - point[dimension] > stepmax + kBVHBoxTolerance ||
        point[dimension] - upper > stepmax + kBVHBoxTolerance) {
      return TGeoShape::Big();
    }
  }

  return fImpl->nearestCrossing<true>(makeVec3(point), makeVec3(dir), stepmax);
}

Double_t O2BVHSurfaceSolid::DistFromInside(const Double_t* point, const Double_t* dir, Int_t, Double_t stepmax,
                                          Double_t* safe) const
{
  if (safe != nullptr) {
    *safe = Safety(point, kTRUE);
  }
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }
  if (fImpl->bvh == nullptr) {
    return DistFromInside_Loop(point, dir, stepmax);
  }
  // no bounding-box reject here: the point is inside by contract, so the box is always reachable
  return fImpl->nearestCrossing<false>(makeVec3(point), makeVec3(dir), stepmax);
}

Double_t O2BVHSurfaceSolid::DistFromOutside_Loop(const Double_t* point, const Double_t* dir, Double_t stepmax) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }
  return fImpl->nearestCrossingLoop<true>(makeVec3(point), makeVec3(dir), stepmax);
}

Double_t O2BVHSurfaceSolid::DistFromInside_Loop(const Double_t* point, const Double_t* dir, Double_t stepmax) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }
  return fImpl->nearestCrossingLoop<false>(makeVec3(point), makeVec3(dir), stepmax);
}

void O2BVHSurfaceSolid::SetRayTMaxPruning(bool enable)
{
  gRayTMaxPruning = enable;
}

bool O2BVHSurfaceSolid::GetRayTMaxPruning()
{
  return gRayTMaxPruning;
}

void O2BVHSurfaceSolid::ResetRayCandidateCounter()
{
  gRayCandidateCount = 0;
}

long long O2BVHSurfaceSolid::GetRayCandidateCount()
{
  return gRayCandidateCount;
}

void O2BVHSurfaceSolid::ResetSafetyCandidateCounter()
{
  gSafetyCandidateCount = 0;
}

long long O2BVHSurfaceSolid::GetSafetyCandidateCount()
{
  return gSafetyCandidateCount;
}

void O2BVHSurfaceSolid::SetSafetyBoundUnsoundForTest(bool enable)
{
  gSafetyBoundUnsound = enable;
}

bool O2BVHSurfaceSolid::GetSafetyBoundUnsoundForTest()
{
  return gSafetyBoundUnsound;
}

/// The safety is the distance to the nearest patch, rounded *down* by one ulp.
///
/// The rounding is asymmetric on purpose and predates this traversal: a safety a hair too small
/// costs the navigator one boundary computation it could have skipped, a safety a hair too large
/// lets it step through a wall. Every bound inside nearestPatchDistanceSq errs the same way.
Double_t O2BVHSurfaceSolid::Safety(const Double_t* point, Bool_t) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }
  const double bestDistanceSq = fImpl->nearestPatchDistanceSq<false>(makeVec3(point), nullptr);
  return std::nextafter(std::sqrt(bestDistanceSq), 0.);
}

Double_t O2BVHSurfaceSolid::Safety_Loop(const Double_t* point, Bool_t) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }
  const double bestDistanceSq = fImpl->nearestPatchDistanceSqLoop(makeVec3(point), nullptr);
  return std::nextafter(std::sqrt(bestDistanceSq), 0.);
}

void O2BVHSurfaceSolid::ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) const
{
  computeNormalFrom(point, dir, norm, false);
}

void O2BVHSurfaceSolid::ComputeNormal_Loop(const Double_t* point, const Double_t* dir, Double_t* norm) const
{
  computeNormalFrom(point, dir, norm, true);
}

void O2BVHSurfaceSolid::computeNormalFrom(const Double_t* point, const Double_t* dir, Double_t* norm,
                                          bool useLoop) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    norm[0] = 1.;
    norm[1] = 0.;
    norm[2] = 0.;
    return;
  }

  const Vec3 testPoint = makeVec3(point);
  size_t closestIndex = fImpl->surfaces.size();
  if (useLoop) {
    fImpl->nearestPatchDistanceSqLoop(testPoint, &closestIndex);
  } else {
    fImpl->nearestPatchDistanceSq<true>(testPoint, &closestIndex);
  }

  if (closestIndex >= fImpl->surfaces.size()) {
    norm[0] = 1.;
    norm[1] = 0.;
    norm[2] = 0.;
    return;
  }

  Vec3 normal = fImpl->surfaces[closestIndex]->normalAt(testPoint);
  if (dir != nullptr) {
    const Vec3 direction = makeVec3(dir);
    if (dot(normal, direction) < 0.) {
      normal = normal * -1.;
    }
  }
  norm[0] = normal.xCoord;
  norm[1] = normal.yCoord;
  norm[2] = normal.zCoord;
}

Double_t O2BVHSurfaceSolid::Capacity() const
{
  if (fImpl == nullptr) {
    return 0.;
  }

  double capacity = 0.;
  for (const auto& surface : fImpl->surfaces) {
    capacity += surface->capacityContribution();
  }
  return std::abs(capacity);
}

bool O2BVHSurfaceSolid::RebuildFromRecords()
{
  // Add*Surface refuses to run on a defined shape and re-appends to fRecords as it replays, so
  // take the records aside and start from a fresh implementation.
  std::vector<BVHSurfaceRecord> records;
  records.swap(fRecords);
  delete fImpl;
  fImpl = new Impl;

  if (records.empty()) {
    Error("RebuildFromRecords",
          "Shape %s carries no surface records: it was either written before CloseShape() or by a version that did "
          "not persist them. The shape stays undefined and not navigable rather than pretending to be an empty "
          "closed solid; check GetNavigationReliability() before using it.",
          GetName());
    return false;
  }

  for (size_t recordIndex = 0; recordIndex < records.size(); ++recordIndex) {
    const auto& record = records[recordIndex];
    const int expectedScalars = BVHSurfaceRecord::expectedScalarCount(record.kind);
    if (expectedScalars < 0 || record.scalars.size() != static_cast<size_t>(expectedScalars)) {
      Error("RebuildFromRecords", "Shape %s: surface record %d has kind %d with %d scalar(s), expected %d", GetName(),
            static_cast<int>(recordIndex), record.kind, static_cast<int>(record.scalars.size()), expectedScalars);
      fRecords.clear();
      delete fImpl;
      fImpl = new Impl;
      return false;
    }

    const Point3D origin = makePoint3D(record.origin);
    const Point3D axisA = makePoint3D(record.axisA);
    const Point3D axisB = makePoint3D(record.axisB);
    const auto& s = record.scalars;

    std::vector<PlanarBoundaryCurve> outerWire;
    std::vector<std::vector<PlanarBoundaryCurve>> innerWires;
    std::vector<Point2D> outerPolygon;
    std::vector<std::vector<Point2D>> innerPolygons;
    const bool wiresLoaded = record.kind == BVHSurfaceRecord::PlanarPolygon
                               ? loadPolygonWires(record, outerPolygon, innerPolygons)
                               : loadCurveWires(record, outerWire, innerWires);

    bool added = false;
    if (!wiresLoaded) {
      Error("RebuildFromRecords", "Shape %s: surface record %d has inconsistent wire sizes", GetName(),
            static_cast<int>(recordIndex));
    } else {
      switch (record.kind) {
        case BVHSurfaceRecord::PlanarPolygon:
          added = AddPlanarSurface(origin, axisA, axisB, outerPolygon, innerPolygons);
          break;
        case BVHSurfaceRecord::CurvedPlanar:
          added = AddCurvedPlanarSurface(origin, axisA, axisB, outerWire, innerWires);
          break;
        case BVHSurfaceRecord::Cylindrical:
          added = record.trimmed ? AddCylindricalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4],
                                                         record.innerWall, outerWire, innerWires)
                                 : AddCylindricalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4],
                                                         record.innerWall);
          break;
        case BVHSurfaceRecord::Spherical:
          added = record.trimmed ? AddSphericalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4],
                                                       record.innerWall, outerWire, innerWires)
                                 : AddSphericalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4],
                                                       record.innerWall);
          break;
        case BVHSurfaceRecord::Conical:
          added = record.trimmed ? AddConicalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4], s[5],
                                                     record.innerWall, outerWire, innerWires)
                                 : AddConicalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4], s[5],
                                                     record.innerWall);
          break;
        case BVHSurfaceRecord::Toroidal:
          added = record.trimmed ? AddToroidalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4], s[5],
                                                      record.innerWall, outerWire, innerWires)
                                 : AddToroidalSurface(origin, axisA, axisB, s[0], s[1], s[2], s[3], s[4], s[5],
                                                      record.innerWall);
          break;
        default:
          break;
      }
    }

    if (!added) {
      // A solid that is missing a face is not the solid that was written; it is a different one
      // whose Contains() is wrong throughout that face's shadow. Discard it entirely rather than
      // hand back a plausible-looking partial shape (K7's failure mode, at the persistence layer).
      Error("RebuildFromRecords",
            "Shape %s: surface record %d (kind %d) did not rebuild. The shape is discarded and stays undefined; it "
            "would otherwise be a *different* solid than the one that was written.",
            GetName(), static_cast<int>(recordIndex), record.kind);
      fRecords.clear();
      delete fImpl;
      fImpl = new Impl;
      return false;
    }

    // The edge identities are part of the record, not derived state: replay them too, or a
    // read-back solid would silently fall back on the geometric rim verdict and could disagree
    // with the solid that was written about whether it is closed.
    if (!record.boundaryEdgeIds.empty()) {
      SetSurfaceBoundaryEdges(static_cast<int>(recordIndex), record.boundaryEdgeIds, record.boundaryEdgeFlags);
    }
  }

  // check == false: replaying a solid must not re-emit the closure diagnostics that were already
  // reported when it was first built. The report itself is recomputed, not trusted.
  CloseShape(false);
  return true;
}

void O2BVHSurfaceSolid::Streamer(TBuffer& buffer)
{
  if (buffer.IsReading()) {
    buffer.ReadClassBuffer(O2BVHSurfaceSolid::Class(), this);
    RebuildFromRecords();
  } else {
    buffer.WriteClassBuffer(O2BVHSurfaceSolid::Class(), this);
  }
}