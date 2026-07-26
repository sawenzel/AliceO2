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

// Ray tmax tightening in the BVH distance queries; see O2BVHSurfaceSolid::SetRayTMaxPruning.
bool gRayTMaxPruning = true;
// Per-thread diagnostic counter of leaf surface patches visited by the BVH distance queries.
thread_local long long gRayCandidateCount = 0;

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

/// Whether \a hit is an entering (\a wantEntering) or an exiting crossing for a ray travelling
/// along \a rayDirection. The outward normal opposes the direction when entering. Grazes within
/// kTolerance of tangency are neither and are discarded: the surface is not actually crossed
/// there, and counting one would report a spurious boundary to the navigator.
template <bool wantEntering>
inline bool isWantedCrossing(const RayHit& hit, const Vec3& rayDirection)
{
  const double alignment = dot(hit.normal, rayDirection);
  if constexpr (wantEntering) {
    return alignment < -kTolerance;
  } else {
    return alignment > kTolerance;
  }
}

// Ray-parity evaluation of a full intersection list (sorts hits in place). Near-equal
// intersections are clustered (shared edges/vertices seen by several surfaces). A cluster whose
// hits all enter or all exit is one genuine crossing; a mixed cluster means the ray grazes an
// edge or corner tangentially (e.g. the rim where a cap meets a wall) and must contribute even
// parity.
bool oddCrossingParity(std::vector<RayHit>& hits, const Vec3& rayDirection)
{
  std::sort(hits.begin(), hits.end(),
            [](const RayHit& firstHit, const RayHit& secondHit) { return firstHit.distance < secondHit.distance; });

  int crossings = 0;
  size_t hitIndex = 0;
  while (hitIndex < hits.size()) {
    bool entering = false;
    bool exiting = false;
    size_t clusterEnd = hitIndex;
    while (clusterEnd < hits.size() &&
           (clusterEnd == hitIndex || sameIntersection(hits[clusterEnd].distance, hits[clusterEnd - 1].distance))) {
      if (dot(hits[clusterEnd].normal, rayDirection) < 0.) {
        entering = true;
      } else {
        exiting = true;
      }
      ++clusterEnd;
    }
    if (entering != exiting) {
      ++crossings;
    }
    hitIndex = clusterEnd;
  }
  return (crossings & 1) != 0;
}
} // namespace

struct O2BVHSurfaceSolid::Impl {
  std::vector<std::unique_ptr<BoundedSurface>> surfaces;
  std::vector<Vec3> displayVertices;
  std::vector<std::array<int, 3>> displayTriangles;
  ClosureReport closure;
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

    double bestDistance = TGeoShape::Big();
    BVHRay ray(BVHVec3(rayOrigin.xCoord, rayOrigin.yCoord, rayOrigin.zCoord),
               BVHVec3(rayDirection.xCoord, rayDirection.yCoord, rayDirection.zCoord), 0.f,
               truncateRoundUp(stepmax));
    static constexpr bool useRobustTraversal = true;
    const bool pruning = gRayTMaxPruning;

    bvh::v2::GrowingStack<BVH::Index> stack;
    // ray is captured by reference on purpose: bvh2 takes it as const Ray&, but the object
    // itself is ours and mutable, and the traversal reads tmax afresh at every node test.
    bvh->intersect<false, useRobustTraversal>(
      ray, bvh->get_root().index, stack, [&](size_t beginPrimitive, size_t endPrimitive) {
        for (size_t primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
          const BoundedSurface& surface = *surfaces[bvh->prim_ids[primitive]];
          ++gRayCandidateCount;
          surfaceHits.clear();
          surface.appendIntersections(rayOrigin, rayDirection, kRayTolerance, std::min(stepmax, bestDistance),
                                      surfaceHits);
          for (const auto& hit : surfaceHits) {
            if (isWantedCrossing<wantEntering>(hit, rayDirection) && hit.distance < bestDistance) {
              bestDistance = hit.distance;
            }
          }
        }
        if (pruning && bestDistance < stepmax) {
          ray.tmax = std::min(ray.tmax, truncateRoundUp(bestDistance + kBVHBoxTolerance));
        }
        return false; // keep traversing; the shrunk tmax does the pruning
      });
    return bestDistance;
  }

  /// Same query without the BVH: visit every surface. Oracle and baseline for nearestCrossing.
  template <bool wantEntering>
  double nearestCrossingLoop(const Vec3& rayOrigin, const Vec3& rayDirection, double stepmax) const
  {
    static thread_local std::vector<RayHit> surfaceHits;

    double bestDistance = TGeoShape::Big();
    for (const auto& surface : surfaces) {
      surfaceHits.clear();
      // the shrinking upper bound prunes surfaces that cannot beat the current best hit
      surface->appendIntersections(rayOrigin, rayDirection, kRayTolerance, std::min(stepmax, bestDistance),
                                   surfaceHits);
      for (const auto& hit : surfaceHits) {
        if (isWantedCrossing<wantEntering>(hit, rayDirection) && hit.distance < bestDistance) {
          bestDistance = hit.distance;
        }
      }
    }
    return bestDistance;
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
};

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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
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

  fImpl->surfaces.emplace_back(std::move(surface));
  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  return true;
}

void O2BVHSurfaceSolid::CloseShape(bool check)
{
  if (fImpl == nullptr) {
    fImpl = new Impl;
  }
  if (check && fImpl->surfaces.empty()) {
    Error("CloseShape", "Shape %s has no bounded surfaces", GetName());
    return;
  }

  ComputeBBox();
  fImpl->buildBVH();

  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  for (const auto& surface : fImpl->surfaces) {
    surface->appendDisplayMesh(fImpl->displayVertices, fImpl->displayTriangles);
  }

  fImpl->closure = validateClosure(fImpl->surfaces);
  fImpl->defined = true;

  if (check && !fImpl->surfaces.empty()) {
    const auto& closure = fImpl->closure;
    // State the *consequence*, not only the counts. A warning that reads as a quality metric gets
    // read as one: this defect went unnoticed for two sessions because "699 boundary edge(s)"
    // looked like a cosmetic mesh remark while the solid was answering Contains wrongly by
    // centimetres. See scripts/geometry/ExactTrimTopology.md.
    if (closure.boundaryEdges > 0) {
      Error("CloseShape",
            "Shape %s is NOT a closed surface: %d boundary edge(s) belong to only one face, i.e. the surface set "
            "has gaps. NAVIGATION IS UNRELIABLE: parity containment is undefined in the shadow of every gap, so "
            "Contains()/DistFrom*() can be wrong arbitrarily far from any surface, not just near the gap. Check "
            "IsNavigable()/GetNavigationReliability() before trusting any answer from this shape.",
            GetName(), closure.boundaryEdges);
    }
    if (closure.nonManifoldEdges > 0) {
      Error("CloseShape",
            "Shape %s is NOT a 2-manifold: %d edge(s) are shared by more than two faces (coincident or duplicated "
            "faces). NAVIGATION IS UNRELIABLE: parity containment depends on the order in which coincident hits "
            "are clustered, so Contains() is not even self-consistent between traversal orders. Check "
            "IsNavigable()/GetNavigationReliability() before trusting any answer from this shape.",
            GetName(), closure.nonManifoldEdges);
    }
    if (!closure.orientationConsistent) {
      Error("CloseShape",
            "Shape %s has %d inconsistently oriented (reversed) boundary edge(s): adjacent faces disagree about "
            "which side is out. NAVIGATION IS UNRELIABLE: the entering/exiting sign of a crossing is taken from the "
            "surface normal, so distance queries can return the wrong side. Check IsNavigable()/"
            "GetNavigationReliability() before trusting any answer from this shape.",
            GetName(), closure.reversedEdges);
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
  // worst defect wins; the enum is ordered by severity
  const auto& closure = fImpl->closure;
  if (closure.nonManifoldEdges > 0) {
    return NavigationReliability::NonManifold;
  }
  if (closure.boundaryEdges > 0) {
    return NavigationReliability::OpenSurfaceSet;
  }
  if (closure.reversedEdges > 0) {
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

  const Vec3 testPoint = makeVec3(point);
  if (std::abs(testPoint.xCoord - fOrigin[0]) > fDX + kTolerance ||
      std::abs(testPoint.yCoord - fOrigin[1]) > fDY + kTolerance ||
      std::abs(testPoint.zCoord - fOrigin[2]) > fDZ + kTolerance) {
    return false;
  }

  if (fImpl->bvh == nullptr) {
    // before CloseShape there is no acceleration structure yet; stay usable via the plain loop
    return Contains_Loop(point);
  }

  // boundary policy: a point within tolerance of any surface patch counts as inside
  if (fImpl->visitPointCandidates(
        testPoint, [&](const BoundedSurface& surface) { return surface.containsPointOnSurface(testPoint); })) {
    return true;
  }

  // reused across calls so containment allocates nothing on the hot path; the capacity is paid
  // once per thread. Distinct from the distance queries' buffers, which are their own.
  static thread_local std::vector<RayHit> containsHits;
  containsHits.clear();
  fImpl->visitRayCandidates(testPoint, kContainsTestDirection, [&](const BoundedSurface& surface) {
    surface.appendIntersections(testPoint, kContainsTestDirection, kRayTolerance, TGeoShape::Big(), containsHits);
  });
  return oddCrossingParity(containsHits, kContainsTestDirection);
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

  static thread_local std::vector<RayHit> containsLoopHits;
  containsLoopHits.clear();
  for (const auto& surface : fImpl->surfaces) {
    surface->appendIntersections(testPoint, kContainsTestDirection, kRayTolerance, TGeoShape::Big(), containsLoopHits);
  }
  return oddCrossingParity(containsLoopHits, kContainsTestDirection);
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

Double_t O2BVHSurfaceSolid::Safety(const Double_t* point, Bool_t) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }

  const Vec3 testPoint = makeVec3(point);
  double bestDistanceSq = std::numeric_limits<double>::infinity();
  for (const auto& surface : fImpl->surfaces) {
    bestDistanceSq = std::min(bestDistanceSq, surface->distanceSqToPatch(testPoint));
  }
  return std::nextafter(std::sqrt(bestDistanceSq), 0.);
}

void O2BVHSurfaceSolid::ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    norm[0] = 1.;
    norm[1] = 0.;
    norm[2] = 0.;
    return;
  }

  const Vec3 testPoint = makeVec3(point);
  const BoundedSurface* closestSurface = nullptr;
  double bestDistanceSq = std::numeric_limits<double>::infinity();
  for (const auto& surface : fImpl->surfaces) {
    const double surfaceDistanceSq = surface->distanceSqToPatch(testPoint);
    if (surfaceDistanceSq < bestDistanceSq) {
      bestDistanceSq = surfaceDistanceSq;
      closestSurface = surface.get();
    }
  }

  if (closestSurface == nullptr) {
    norm[0] = 1.;
    norm[1] = 0.;
    norm[2] = 0.;
    return;
  }

  Vec3 normal = closestSurface->normalAt(testPoint);
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

void O2BVHSurfaceSolid::Streamer(TBuffer& buffer)
{
  if (buffer.IsReading()) {
    buffer.ReadClassBuffer(O2BVHSurfaceSolid::Class(), this);
    CloseShape(false);
  } else {
    buffer.WriteClassBuffer(O2BVHSurfaceSolid::Class(), this);
  }
}