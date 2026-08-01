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

/// \file BoundedSurface.h
/// \brief Private bounded-surface / wire / edge interface used by O2BVHSurfaceSolid.
///
/// This header is intentionally implementation-private (it lives under src/ and is only
/// meant to be included by O2BVHSurfaceSolid.cxx and its focused unit test). It defines the
/// analytic bounded-surface abstraction:
///  - SurfaceEdge: one bounded curve segment (currently a straight line) in a surface's 2D
///    parametric domain,
///  - SurfaceWire: one closed, oriented boundary loop made of edges,
///  - BoundedSurface: the abstract analytic surface patch (support surface + trim wires),
///  - PlanarBoundedSurface: the first concrete surface,
///  - DummyBoundedSurface: a trivial surface used only by unit tests.
///
/// All types live in namespace o2::base::surface and are header-inline so the header can be
/// shared between the library translation unit and the test translation unit without extra
/// linking.

#ifndef ALICEO2_BASE_BOUNDEDSURFACE_H_
#define ALICEO2_BASE_BOUNDEDSURFACE_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::base::surface
{

/// \name Numerical conventions
/// Implementation-local tolerances shared by all bounded-surface code. Keeping them in one
/// place avoids scattered magic literals (see the "numerical conventions" data-model milestone).
/// @{
inline constexpr double kTolerance = 1.e-9;              ///< generic length tolerance
inline constexpr double kToleranceSq = kTolerance * kTolerance;
inline constexpr double kAreaTolerance = 1.e-18;         ///< degenerate (zero) parametric area
inline constexpr double kRayTolerance = 1.e-9;           ///< minimum positive ray parameter t
inline constexpr double kIntersectionTolerance = 1.e-7;  ///< clustering of near-equal intersections
inline constexpr double kClosureQuantum = 1.e-7;         ///< vertex quantization for closure matching
/// Tolerance for accepting a wire loop as closed (consecutive edge/curve endpoints meeting), as a
/// **3D length in cm**, measured through the owning surface's first fundamental form (see
/// BoundedSurface::parametricMetric and ParametricMetric below).
///
/// It is deliberately looser than kTolerance: a wire is imported from a CAD extractor that
/// samples/projects each boundary curve's endpoints independently (e.g. a straight line vertex vs.
/// the end pole of a neighbouring B-spline/arc on a quadric), so the shared vertex can differ by
/// the extractor's precision even though the loop is closed. The value is that precision, ~1e-6
/// cm, and not the historical 1e-5: the old constant was compared against sqrt(du^2 + dv^2) in a
/// domain that mixes radians with centimetres, so it silently meant 1e-5 cm of arc on a 1 cm
/// cylinder and 1e-3 cm on a 100 cm one, while a 0.01 cm hole had to close to 1e-7 cm. Small holes
/// were falsely open and big cylinders falsely closed, from the one missing radius factor.
///
/// The residual gap is negligible for winding/area/navigation, which use kTolerance for
/// on-boundary classification.
inline constexpr double kWireJoinTolerance = 1.e-6;
inline constexpr double kWireJoinToleranceSq = kWireJoinTolerance * kWireJoinTolerance;
/// How flat the adaptive B-spline sampler makes each chord, in the curve's own parametric units.
/// A trim B-spline enters navigation as this polyline, so it *is* the boundary as far as winding
/// and point-to-curve distance are concerned, and nothing downstream can be more accurate than
/// this. The on-boundary band is set from it rather than from optimism (finding K5): a 1e-9 band
/// around a 1e-5 approximation is not a band, it is noise, and it made the Boundary state
/// unreachable for every B-spline trim.
inline constexpr double kBSplineFlatness = 1.e-5;
inline constexpr double kBSplineFlatnessSq = kBSplineFlatness * kBSplineFlatness;
/// How close two faces' trim boundaries must come, as a 3D length in cm, to count as the same
/// edge in the rim-based closure measurement -- used only when the source model states no
/// tolerance of its own (O2BVHSurfaceSolid::GetModelTolerance() returns 0 for "not stated").
///
/// Same value and same origin as kWireJoinTolerance: it is the CAD extractor's precision, which
/// is as close as two independently sampled boundaries of one model can be asked to come. A
/// model that declares its own tolerance gets that instead, because it knows better.
inline constexpr double kRimMatchTolerance = 1.e-6;

/// Expansion of the per-surface BVH leaf AABBs (in double, before the outward float rounding).
/// It must dominate every length tolerance used by navigation queries (kTolerance boundary
/// classification, kIntersectionTolerance clustering) so a point or ray hit within tolerance of
/// a patch is never pruned away by the BVH traversal.
inline constexpr double kBVHBoxTolerance = 1.e-3;
/// @}

/// A 2D point/vector in a surface's parametric (u, v) domain.
struct Vec2 {
  double uCoord = 0.;
  double vCoord = 0.;
};

/// A 3D point/vector in the solid's local frame.
struct Vec3 {
  double xCoord = 0.;
  double yCoord = 0.;
  double zCoord = 0.;
};

inline Vec3 operator+(const Vec3& firstVector, const Vec3& secondVector)
{
  return {firstVector.xCoord + secondVector.xCoord, firstVector.yCoord + secondVector.yCoord,
          firstVector.zCoord + secondVector.zCoord};
}

inline Vec3 operator-(const Vec3& firstVector, const Vec3& secondVector)
{
  return {firstVector.xCoord - secondVector.xCoord, firstVector.yCoord - secondVector.yCoord,
          firstVector.zCoord - secondVector.zCoord};
}

inline Vec3 operator*(const Vec3& vector, double scale)
{
  return {vector.xCoord * scale, vector.yCoord * scale, vector.zCoord * scale};
}

inline Vec3 operator*(double scale, const Vec3& vector)
{
  return vector * scale;
}

inline Vec2 operator-(const Vec2& firstPoint, const Vec2& secondPoint)
{
  return {firstPoint.uCoord - secondPoint.uCoord, firstPoint.vCoord - secondPoint.vCoord};
}

/// The 3D length squared of a parametric displacement \a delta, given a surface's first
/// fundamental form (\a gUU, \a gUV, \a gVV) at the point the displacement starts from.
///
/// This is the only honest way to compare a parametric separation against a tolerance: the (u, v)
/// domains of this file mix radians and centimetres, so `du^2 + dv^2` is a number with no physical
/// meaning, and any constant compared against it silently means a different thing on every
/// surface. See BoundedSurface::parametricMetric.
inline double parametricLengthSq(double gUU, double gUV, double gVV, const Vec2& delta)
{
  return gUU * delta.uCoord * delta.uCoord + 2. * gUV * delta.uCoord * delta.vCoord +
         gVV * delta.vCoord * delta.vCoord;
}

/// How a wire converts a parametric separation into a 3D length.
///
/// A wire is built before -- and independently of -- the surface that owns it, and both wire types
/// are deliberately usable on their own, so the surface hands its first fundamental form in as a
/// callback rather than the wire reaching back for it. This is the smaller of the two designs
/// TolerancePolicy.md section 2.2 offers, and it is the one that keeps the wires testable; the
/// alternative (moving join validation up into the six surface classes) would have needed a shared
/// helper that amounts to the same callback anyway.
///
/// A default-constructed metric is the identity, which is the only honest reading without surface
/// context: "(u, v) are already lengths in cm".
struct ParametricMetric {
  using Evaluate = void (*)(const void* context, const Vec2& uv, double& gUU, double& gUV, double& gVV);

  Evaluate evaluate = nullptr;
  const void* context = nullptr;

  /// The 3D length squared spanned by the parametric displacement \a delta starting at \a uv.
  double lengthSq(const Vec2& uv, const Vec2& delta) const
  {
    if (evaluate == nullptr) {
      return delta.uCoord * delta.uCoord + delta.vCoord * delta.vCoord;
    }
    double gUU = 1.;
    double gUV = 0.;
    double gVV = 1.;
    evaluate(context, uv, gUU, gUV, gVV);
    return parametricLengthSq(gUU, gUV, gVV, delta);
  }

  /// The 3D distance squared between two parametric points. The form varies over the domain, so
  /// it is evaluated at \a from -- for the separations this is used on (a join gap, a duplicate
  /// vertex) the two points are within tolerance of each other and the choice does not matter.
  double distanceSq(const Vec2& from, const Vec2& to) const { return lengthSq(from, to - from); }

  /// The largest 3D length a *unit* parametric displacement can span at \a uv, i.e. the square
  /// root of the larger eigenvalue of the form. Use it to convert a bound whose direction is not
  /// known -- a sampling tolerance, say -- into the length it can worst-case reach; lengthSq is
  /// the exact answer whenever the displacement itself is known.
  double maxScale(const Vec2& uv) const
  {
    if (evaluate == nullptr) {
      return 1.;
    }
    double gUU = 1.;
    double gUV = 0.;
    double gVV = 1.;
    evaluate(context, uv, gUU, gUV, gVV);
    const double trace = gUU + gVV;
    const double determinant = gUU * gVV - gUV * gUV;
    // the eigenvalues of a symmetric 2x2 form, guarded against a slightly negative discriminant
    const double discriminant = std::max(0., trace * trace - 4. * determinant);
    return std::sqrt(std::max(0., 0.5 * (trace + std::sqrt(discriminant))));
  }
};

/// A ParametricMetric that defers to \a surface, which must outlive it. Every use here is a
/// surface building its own wires inside initialize(), so that holds by construction.
template <typename Surface>
inline ParametricMetric parametricMetricOf(const Surface& surface)
{
  return {[](const void* context, const Vec2& uv, double& gUU, double& gUV, double& gVV) {
            static_cast<const Surface*>(context)->parametricMetric(uv, gUU, gUV, gVV);
          },
          &surface};
}

inline double dot(const Vec3& firstVector, const Vec3& secondVector)
{
  return firstVector.xCoord * secondVector.xCoord + firstVector.yCoord * secondVector.yCoord +
         firstVector.zCoord * secondVector.zCoord;
}

inline Vec3 cross(const Vec3& firstVector, const Vec3& secondVector)
{
  return {firstVector.yCoord * secondVector.zCoord - firstVector.zCoord * secondVector.yCoord,
          firstVector.zCoord * secondVector.xCoord - firstVector.xCoord * secondVector.zCoord,
          firstVector.xCoord * secondVector.yCoord - firstVector.yCoord * secondVector.xCoord};
}

inline double normSq(const Vec3& vector)
{
  return dot(vector, vector);
}

inline double norm(const Vec3& vector)
{
  return std::sqrt(normSq(vector));
}

inline Vec3 normalized(const Vec3& vector)
{
  const double vectorNorm = norm(vector);
  if (vectorNorm <= kTolerance) {
    return {};
  }
  return vector * (1. / vectorNorm);
}

inline double component(const Vec3& vector, int dimension)
{
  if (dimension == 0) {
    return vector.xCoord;
  }
  if (dimension == 1) {
    return vector.yCoord;
  }
  return vector.zCoord;
}

inline void assignComponent(Vec3& vector, int dimension, double value)
{
  if (dimension == 0) {
    vector.xCoord = value;
  } else if (dimension == 1) {
    vector.yCoord = value;
  } else {
    vector.zCoord = value;
  }
}

inline bool finite(const Vec2& point)
{
  return std::isfinite(point.uCoord) && std::isfinite(point.vCoord);
}

inline bool finite(const Vec3& point)
{
  return std::isfinite(point.xCoord) && std::isfinite(point.yCoord) && std::isfinite(point.zCoord);
}

inline double distanceSq(const Vec2& firstPoint, const Vec2& secondPoint)
{
  const double deltaU = firstPoint.uCoord - secondPoint.uCoord;
  const double deltaV = firstPoint.vCoord - secondPoint.vCoord;
  return deltaU * deltaU + deltaV * deltaV;
}

inline double distanceSq(const Vec3& firstPoint, const Vec3& secondPoint)
{
  return normSq(firstPoint - secondPoint);
}

inline double cross2D(const Vec2& firstVector, const Vec2& secondVector)
{
  return firstVector.uCoord * secondVector.vCoord - firstVector.vCoord * secondVector.uCoord;
}

inline double pointSegmentDistanceSq(const Vec2& point, const Vec2& segmentStart, const Vec2& segmentEnd)
{
  const Vec2 segmentVector = segmentEnd - segmentStart;
  const double segmentLengthSq = segmentVector.uCoord * segmentVector.uCoord + segmentVector.vCoord * segmentVector.vCoord;
  if (segmentLengthSq <= kToleranceSq) {
    return distanceSq(point, segmentStart);
  }
  const double pointProjection = ((point.uCoord - segmentStart.uCoord) * segmentVector.uCoord +
                                  (point.vCoord - segmentStart.vCoord) * segmentVector.vCoord) /
                                 segmentLengthSq;
  const double clampedProjection = std::max(0., std::min(1., pointProjection));
  const Vec2 closestPoint{segmentStart.uCoord + clampedProjection * segmentVector.uCoord,
                          segmentStart.vCoord + clampedProjection * segmentVector.vCoord};
  return distanceSq(point, closestPoint);
}

inline double pointSegmentDistanceSq(const Vec3& point, const Vec3& segmentStart, const Vec3& segmentEnd)
{
  const Vec3 segmentVector = segmentEnd - segmentStart;
  const double segmentLengthSq = normSq(segmentVector);
  if (segmentLengthSq <= kToleranceSq) {
    return distanceSq(point, segmentStart);
  }
  const double pointProjection = dot(point - segmentStart, segmentVector) / segmentLengthSq;
  const double clampedProjection = std::max(0., std::min(1., pointProjection));
  const Vec3 closestPoint = segmentStart + segmentVector * clampedProjection;
  return distanceSq(point, closestPoint);
}

/// \name First fundamental forms, by surface family
/// The closed forms behind BoundedSurface::parametricMetric, as free functions of the parameters
/// that define each family. The surface classes call these, and so does the sidecar reader, which
/// has a record's parameters but no surface yet and must apply the same join rule the kernel will
/// (finding S10). One definition, two callers -- the formulas are not to be written twice.
/// @{

/// Plane: the frame axes carry the domain's units and need be neither unit nor orthogonal, which
/// makes this the only family with a cross term.
inline void planeParametricMetric(const Vec3& axisU, const Vec3& axisV, double& gUU, double& gUV, double& gVV)
{
  gUU = dot(axisU, axisU);
  gUV = dot(axisU, axisV);
  gVV = dot(axisV, axisV);
}

/// Cylinder, (u, v) = (phi[rad], h[cm]).
inline void cylinderParametricMetric(double radius, double& gUU, double& gUV, double& gVV)
{
  gUU = radius * radius;
  gUV = 0.;
  gVV = 1.;
}

/// Cone, (u, v) = (phi[rad], h[cm]). \a radiusAtHeight is r(v), which reaches zero at an apex;
/// a step in h also walks along the slope, hence gVV > 1.
inline void coneParametricMetric(double radiusAtHeight, double slope, double& gUU, double& gUV, double& gVV)
{
  gUU = radiusAtHeight * radiusAtHeight;
  gUV = 0.;
  gVV = 1. + slope * slope;
}

/// Sphere, (u, v) = (phi[rad], theta[rad]). The azimuthal scale is the radius of the parallel at
/// \a theta, so it vanishes at either pole.
inline void sphereParametricMetric(double radius, double theta, double& gUU, double& gUV, double& gVV)
{
  const double parallelRadius = radius * std::sin(theta);
  gUU = parallelRadius * parallelRadius;
  gUV = 0.;
  gVV = radius * radius;
}

/// Torus, (u, v) = (phiRing[rad], phiTube[rad]). The ring scale runs from R - r to R + r.
inline void torusParametricMetric(double majorRadius, double minorRadius, double phiTube, double& gUU, double& gUV,
                                  double& gVV)
{
  const double ringRadius = majorRadius + minorRadius * std::cos(phiTube);
  gUU = ringRadius * ringRadius;
  gUV = 0.;
  gVV = minorRadius * minorRadius;
}
/// @}

inline bool sameIntersection(double firstDistance, double secondDistance)
{
  return std::abs(firstDistance - secondDistance) <=
         kIntersectionTolerance * std::max(1., std::max(std::abs(firstDistance), std::abs(secondDistance)));
}

/// One ray/surface intersection: the ray parameter and the outward-oriented surface normal at
/// the hit. Curved (quadric) patches can produce more than one hit per ray, which parity-based
/// containment and entering/exiting distance queries must see; hence intersections are reported
/// as a list rather than as a single nearest hit.
struct RayHit {
  double distance = 0.;
  Vec3 normal;
  /// The hit lies within the patch's own on-boundary band (CurveWire::boundaryBand), i.e. within
  /// the accuracy of the representation that draws the trim -- so which side of the trim it falls
  /// on is not decided by the data, it is decided by the tie-break. The band is resolved as
  /// "inside the trim", which is a choice and not a fact, and the choice is one-sided: a patch
  /// keeps a sliver of width band past its true trim curve.
  ///
  /// This is the flag TolerancePolicy.md section 2.4 asks for. It matters because on a Boolean
  /// seam that sliver lies in the solid's interior, where a crossing must not be counted, so a ray
  /// through it gains a spurious crossing and flips parity. Measured on cyl_cross_cyl: every one
  /// of 1440 sampled positions along the true crossing curve overhangs by 1.0e-5 to 1.9e-5 cm and
  /// none undercuts, and the single direction-dependent point of the section 4.2 sweep is a ray
  /// through that sliver (TolerancePolicy.md section 11).
  bool onTrimBoundary = false;
};

/// One bounded curve segment of a wire, in a surface's parametric (u, v) domain.
///
/// The first supported curve is a straight line segment. Future curve types (circle arcs,
/// etc.) can extend this by making the operations virtual; for now the concrete line segment
/// keeps the model simple and allocation-free.
struct SurfaceEdge {
  Vec2 start;
  Vec2 end;

  Vec2 direction() const { return end - start; }

  double lengthSq() const
  {
    const Vec2 delta = end - start;
    return delta.uCoord * delta.uCoord + delta.vCoord * delta.vCoord;
  }

  bool degenerate() const { return lengthSq() <= kToleranceSq; }

  /// Squared distance from a parametric point to this edge.
  double distanceSq(const Vec2& point) const { return pointSegmentDistanceSq(point, start, end); }

  /// Closest point on this edge to \a point. Returns the projected point and its clamped
  /// parameter \a parameter in [0, 1] (0 at start, 1 at end). Degenerate edges return start.
  Vec2 closestPoint(const Vec2& point, double& parameter) const
  {
    const Vec2 segmentVector = end - start;
    const double segmentLengthSq = segmentVector.uCoord * segmentVector.uCoord +
                                   segmentVector.vCoord * segmentVector.vCoord;
    if (segmentLengthSq <= kToleranceSq) {
      parameter = 0.;
      return start;
    }
    const double projection = ((point.uCoord - start.uCoord) * segmentVector.uCoord +
                               (point.vCoord - start.vCoord) * segmentVector.vCoord) /
                              segmentLengthSq;
    parameter = std::max(0., std::min(1., projection));
    return {start.uCoord + parameter * segmentVector.uCoord, start.vCoord + parameter * segmentVector.vCoord};
  }

  /// Accumulate the edge endpoints into a parametric axis-aligned bounding box.
  void extendBounds(Vec2& lower, Vec2& upper) const
  {
    lower.uCoord = std::min({lower.uCoord, start.uCoord, end.uCoord});
    lower.vCoord = std::min({lower.vCoord, start.vCoord, end.vCoord});
    upper.uCoord = std::max({upper.uCoord, start.uCoord, end.uCoord});
    upper.vCoord = std::max({upper.vCoord, start.vCoord, end.vCoord});
  }
};

/// Classification of a parametric point against a closed wire.
enum class WireClassification { Outside,
                                Boundary,
                                Inside };

/// The role a wire plays for a bounded surface. Outer wires bound the material, inner wires
/// (holes) subtract from it. The role fixes the expected winding relative to the surface normal.
enum class WireRole { Outer,
                      Inner };

/// Outcome of wire construction / validation. Valid and Reversed are both usable results;
/// Reversed additionally signals that the orientation had to be normalized (a logged repair).
enum class WireStatus {
  Valid,           ///< well-formed and already correctly oriented
  Reversed,        ///< well-formed but re-oriented to match its role (simple, logged repair)
  NonFinite,       ///< a vertex/edge contained a non-finite coordinate
  Open,            ///< an explicit edge list did not form a closed loop
  TooFewVertices,  ///< fewer than three distinct vertices after cleanup
  DegenerateVertex,///< a non-adjacent vertex coincided (self-touching / pinched loop)
  ZeroArea         ///< the loop encloses no area
};

/// Human-readable description of a wire status, for logging.
inline const char* wireStatusMessage(WireStatus status)
{
  switch (status) {
    case WireStatus::Valid:
      return "valid";
    case WireStatus::Reversed:
      return "orientation normalized to match wire role";
    case WireStatus::NonFinite:
      return "wire contains a non-finite vertex";
    case WireStatus::Open:
      return "wire edges do not form a closed loop";
    case WireStatus::TooFewVertices:
      return "wire needs at least three distinct vertices";
    case WireStatus::DegenerateVertex:
      return "wire has a coincident (pinched) vertex";
    case WireStatus::ZeroArea:
      return "wire has zero area";
  }
  return "unknown wire status";
}

/// One closed, oriented boundary loop in a surface's parametric domain.
///
/// Storage is a de-duplicated vertex ring; edges are the line segments between consecutive
/// vertices (with the last vertex connecting back to the first). Orientation is normalized so
/// that an outer wire winds counter-clockwise and an inner wire clockwise with respect to the
/// surface normal, which keeps 3D boundary edges of a closed solid consistently oriented.
struct SurfaceWire {
  std::vector<Vec2> vertices;
  WireRole role = WireRole::Outer;

  int edgeCount() const { return static_cast<int>(vertices.size()); }

  SurfaceEdge edge(int index) const
  {
    const int count = edgeCount();
    return {vertices[index % count], vertices[(index + 1) % count]};
  }

  /// Build and validate the wire from an implicitly-closed vertex ring. \a metric converts the
  /// parametric separations below into 3D lengths; a plane's frame axes need be neither unit
  /// length nor orthogonal, so without it "coincident" means a different thing per surface.
  bool initialize(const std::vector<Vec2>& inputVertices, WireRole wireRole, WireStatus& status,
                  const ParametricMetric& metric = {})
  {
    role = wireRole;
    vertices.clear();
    vertices.reserve(inputVertices.size());

    for (const auto& vertex : inputVertices) {
      if (!finite(vertex)) {
        status = WireStatus::NonFinite;
        return false;
      }
      if (vertices.empty() || metric.distanceSq(vertices.back(), vertex) > kToleranceSq) {
        vertices.push_back(vertex);
      }
    }

    // drop an explicit closing duplicate (first == last)
    if (vertices.size() > 1 && metric.distanceSq(vertices.front(), vertices.back()) <= kToleranceSq) {
      vertices.pop_back();
    }

    if (vertices.size() < 3) {
      status = WireStatus::TooFewVertices;
      return false;
    }

    // reject self-touching loops (non-adjacent coincident vertices)
    for (size_t firstIndex = 0; firstIndex < vertices.size(); ++firstIndex) {
      for (size_t secondIndex = firstIndex + 1; secondIndex < vertices.size(); ++secondIndex) {
        if (metric.distanceSq(vertices[firstIndex], vertices[secondIndex]) <= kToleranceSq) {
          status = WireStatus::DegenerateVertex;
          return false;
        }
      }
    }

    const double area = signedArea();
    if (std::abs(area) <= kAreaTolerance) {
      status = WireStatus::ZeroArea;
      return false;
    }

    // outer wires must wind CCW (positive area), inner wires CW (negative area)
    const bool wantPositiveArea = (role == WireRole::Outer);
    if ((area > 0.) != wantPositiveArea) {
      std::reverse(vertices.begin(), vertices.end());
      status = WireStatus::Reversed;
      return true;
    }

    status = WireStatus::Valid;
    return true;
  }

  /// Build and validate the wire from an explicit ordered edge list, checking loop closure.
  ///
  /// The join test is kWireJoinTolerance through \a metric, i.e. exactly the rule CurveWire uses.
  /// The two used to differ -- 1e-9 here against 1e-5 there, in incompatible units -- although the
  /// same extractor, with the same per-endpoint precision, feeds both (finding K12).
  bool initializeFromEdges(const std::vector<SurfaceEdge>& edges, WireRole wireRole, WireStatus& status,
                           const ParametricMetric& metric = {})
  {
    if (edges.size() < 3) {
      status = WireStatus::TooFewVertices;
      return false;
    }
    for (size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
      if (!finite(edges[edgeIndex].start) || !finite(edges[edgeIndex].end)) {
        status = WireStatus::NonFinite;
        return false;
      }
      const Vec2& nextStart = edges[(edgeIndex + 1) % edges.size()].start;
      if (metric.distanceSq(edges[edgeIndex].end, nextStart) > kWireJoinToleranceSq) {
        status = WireStatus::Open;
        return false;
      }
    }

    std::vector<Vec2> ringVertices;
    ringVertices.reserve(edges.size());
    for (const auto& singleEdge : edges) {
      ringVertices.push_back(singleEdge.start);
    }
    return initialize(ringVertices, wireRole, status, metric);
  }

  double signedArea() const
  {
    double area = 0.;
    for (size_t vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex) {
      const auto& currentVertex = vertices[vertexIndex];
      const auto& nextVertex = vertices[(vertexIndex + 1) % vertices.size()];
      area += currentVertex.uCoord * nextVertex.vCoord - nextVertex.uCoord * currentVertex.vCoord;
    }
    return 0.5 * area;
  }

  /// Accumulate this wire's vertices into a parametric axis-aligned bounding box. This is
  /// independent of any concrete surface so cylinders, spheres and cones can reuse it.
  void parametricBounds(Vec2& lower, Vec2& upper) const
  {
    for (const auto& vertex : vertices) {
      lower.uCoord = std::min(lower.uCoord, vertex.uCoord);
      lower.vCoord = std::min(lower.vCoord, vertex.vCoord);
      upper.uCoord = std::max(upper.uCoord, vertex.uCoord);
      upper.vCoord = std::max(upper.vCoord, vertex.vCoord);
    }
  }

  /// Ordered, closed boundary polyline of this wire in its parametric domain. For the current
  /// line-segment edges this is simply the de-duplicated vertex ring closed back to the first
  /// vertex; it gives a stable hook for visualization and for future curved edges (which would
  /// sample themselves into additional points here) without exposing the vertex storage.
  std::vector<Vec2> sampledBoundary() const
  {
    std::vector<Vec2> samples;
    if (vertices.empty()) {
      return samples;
    }
    samples.reserve(vertices.size() + 1);
    samples.insert(samples.end(), vertices.begin(), vertices.end());
    samples.push_back(vertices.front());
    return samples;
  }

  /// \a metric sizes the on-boundary band only. A polygon wire is held exactly, so its band is
  /// just kTolerance -- but kTolerance is a length in cm and this domain need not be, so it is
  /// converted through the metric's largest scale, the same rule CurveWire::boundaryBand applies.
  WireClassification classify(const Vec2& point, const ParametricMetric& metric = {}) const
  {
    const double scale = metric.maxScale(point);
    const double band = scale > kTolerance ? kTolerance / scale : 0.;
    const double bandSq = band * band;
    bool inside = false;
    for (size_t vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex) {
      const auto& segmentStart = vertices[vertexIndex];
      const auto& segmentEnd = vertices[(vertexIndex + 1) % vertices.size()];
      if (pointSegmentDistanceSq(point, segmentStart, segmentEnd) <= bandSq) {
        return WireClassification::Boundary;
      }
      const bool crossesScanline = (segmentStart.vCoord > point.vCoord) != (segmentEnd.vCoord > point.vCoord);
      if (crossesScanline) {
        const double intersectionU = segmentStart.uCoord + (point.vCoord - segmentStart.vCoord) *
                                                             (segmentEnd.uCoord - segmentStart.uCoord) /
                                                             (segmentEnd.vCoord - segmentStart.vCoord);
        if (point.uCoord < intersectionU) {
          inside = !inside;
        }
      }
    }
    return inside ? WireClassification::Inside : WireClassification::Outside;
  }
};

inline bool pointInTriangle(const Vec2& point, const Vec2& firstVertex, const Vec2& secondVertex,
                            const Vec2& thirdVertex)
{
  const double firstCross = cross2D(secondVertex - firstVertex, point - firstVertex);
  const double secondCross = cross2D(thirdVertex - secondVertex, point - secondVertex);
  const double thirdCross = cross2D(firstVertex - thirdVertex, point - thirdVertex);
  return firstCross >= -kTolerance && secondCross >= -kTolerance && thirdCross >= -kTolerance;
}

/// Ear-clipping triangulation of a simple (non-self-intersecting) parametric wire.
inline std::vector<std::array<int, 3>> triangulateSimpleWire(const SurfaceWire& wire)
{
  std::vector<int> remainingIndices;
  remainingIndices.reserve(wire.vertices.size());
  if (wire.signedArea() >= 0.) {
    for (size_t vertexIndex = 0; vertexIndex < wire.vertices.size(); ++vertexIndex) {
      remainingIndices.push_back(static_cast<int>(vertexIndex));
    }
  } else {
    for (size_t reverseIndex = wire.vertices.size(); reverseIndex > 0; --reverseIndex) {
      remainingIndices.push_back(static_cast<int>(reverseIndex - 1));
    }
  }

  std::vector<std::array<int, 3>> triangles;
  size_t guardCounter = 0;
  while (remainingIndices.size() > 3 && guardCounter++ < wire.vertices.size() * wire.vertices.size()) {
    bool clippedEar = false;
    for (size_t indexPosition = 0; indexPosition < remainingIndices.size(); ++indexPosition) {
      const int previousIndex = remainingIndices[(indexPosition + remainingIndices.size() - 1) % remainingIndices.size()];
      const int currentIndex = remainingIndices[indexPosition];
      const int nextIndex = remainingIndices[(indexPosition + 1) % remainingIndices.size()];

      const auto& previousVertex = wire.vertices[previousIndex];
      const auto& currentVertex = wire.vertices[currentIndex];
      const auto& nextVertex = wire.vertices[nextIndex];
      if (cross2D(currentVertex - previousVertex, nextVertex - currentVertex) <= kTolerance) {
        continue;
      }

      bool containsOtherVertex = false;
      for (int candidateIndex : remainingIndices) {
        if (candidateIndex == previousIndex || candidateIndex == currentIndex || candidateIndex == nextIndex) {
          continue;
        }
        if (pointInTriangle(wire.vertices[candidateIndex], previousVertex, currentVertex, nextVertex)) {
          containsOtherVertex = true;
          break;
        }
      }
      if (containsOtherVertex) {
        continue;
      }

      triangles.push_back({previousIndex, currentIndex, nextIndex});
      remainingIndices.erase(remainingIndices.begin() + indexPosition);
      clippedEar = true;
      break;
    }

    if (!clippedEar) {
      break;
    }
  }

  if (remainingIndices.size() == 3) {
    triangles.push_back({remainingIndices[0], remainingIndices[1], remainingIndices[2]});
  }
  return triangles;
}

/// \name Angular constants for parametric arc curves
/// @{
inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kTwoPi = 2. * kPi;
inline constexpr double kHalfPi = 0.5 * kPi;
/// Chord segments used when sampling a full-circle arc for visualization and for the directed
/// boundary edges of curved patches. All surfaces must use the same sampling so that shared
/// circular boundaries (e.g. a cylinder rim and its planar cap) produce matching half-edges in
/// the solid closure check. Must be divisible by 4 so quarter-turn-rotated frames still sample
/// the same point set on a shared circle.
inline constexpr int kArcSamples = 24;
/// @}

/// Angular tolerance equivalent to a kTolerance arc length at the given radius.
inline double angularTolerance(double radius)
{
  return kTolerance / std::max(radius, kTolerance);
}

/// True if \a angle lies within the angular range [start, start + sweep] (sweep in (0, 2pi]),
/// allowing \a tolerance on both ends and treating a >= 2pi sweep as the full circle.
inline bool angleInSweepRange(double angle, double start, double sweep, double tolerance)
{
  if (sweep >= kTwoPi - kTolerance) {
    return true;
  }
  double delta = angle - start;
  delta -= kTwoPi * std::floor(delta / kTwoPi); // wrap into [0, 2pi)
  return delta <= sweep + tolerance || delta >= kTwoPi - tolerance;
}

/// Compute the \a n-point Gauss-Legendre nodes \a nodes and weights \a weights on [-1, 1].
/// Used to integrate the exact signed area of polynomial (non-rational) B-spline spans and,
/// with a higher order, the numeric area of rational spans. Nodes are found by Newton iteration
/// on the Legendre polynomial P_n, which is robust for the small n used here.
inline void gaussLegendre(int n, std::vector<double>& nodes, std::vector<double>& weights)
{
  nodes.assign(std::max(n, 1), 0.);
  weights.assign(std::max(n, 1), 0.);
  if (n < 1) {
    return;
  }
  for (int i = 0; i < n; ++i) {
    double root = std::cos(kPi * (i + 0.75) / (n + 0.5)); // asymptotic initial guess
    double derivative = 1.;
    for (int iteration = 0; iteration < 100; ++iteration) {
      double previous = 1.;
      double current = root;
      for (int degreeIndex = 2; degreeIndex <= n; ++degreeIndex) {
        const double next = ((2 * degreeIndex - 1) * root * current - (degreeIndex - 1) * previous) / degreeIndex;
        previous = current;
        current = next;
      }
      derivative = n * (root * current - previous) / (root * root - 1.);
      const double delta = current / derivative;
      root -= delta;
      if (std::abs(delta) < 1.e-15) {
        break;
      }
    }
    nodes[i] = root;
    weights[i] = 2. / ((1. - root * root) * derivative * derivative);
  }
}

/// Append the real roots of the depressed cubic w^3 + P w + Q = 0. Uses Cardano's formula when
/// there is a single real root and the trigonometric (casus irreducibilis) form when there are
/// three, both numerically stable for the small, well-scaled resolvent cubics of the quartic
/// solver below.
inline void solveDepressedCubic(double coeffP, double coeffQ, std::vector<double>& roots)
{
  if (std::abs(coeffP) <= kTolerance && std::abs(coeffQ) <= kTolerance) {
    roots.push_back(0.);
    return;
  }
  const double discriminant = coeffQ * coeffQ / 4. + coeffP * coeffP * coeffP / 27.;
  if (discriminant > 0.) {
    const double sqrtDiscriminant = std::sqrt(discriminant);
    roots.push_back(std::cbrt(-0.5 * coeffQ + sqrtDiscriminant) + std::cbrt(-0.5 * coeffQ - sqrtDiscriminant));
    return;
  }
  // three real roots: coeffP < 0 here, so the trigonometric form is well defined
  const double magnitude = 2. * std::sqrt(-coeffP / 3.);
  const double cosineArgument = std::max(-1., std::min(1., 3. * coeffQ / (coeffP * magnitude)));
  const double baseAngle = std::acos(cosineArgument);
  for (int branch = 0; branch < 3; ++branch) {
    roots.push_back(magnitude * std::cos((baseAngle - kTwoPi * branch) / 3.));
  }
}

/// Real roots of a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0 (requires a4 != 0), via Ferrari's
/// method with a resolvent cubic, followed by a couple of Newton polishing steps against the
/// original quartic. A tangential (double) root is deliberately returned as a near-equal pair so
/// that crossing-parity callers can cluster and drop it (the same convention the quadric surfaces
/// use with sameIntersection); a genuinely complex pair simply contributes no real root. Used by
/// the toroidal surface, whose ray intersection is a quartic.
inline std::vector<double> solveQuarticReal(double a4, double a3, double a2, double a1, double a0)
{
  std::vector<double> roots;
  if (std::abs(a4) <= kTolerance) {
    return roots; // not a genuine quartic; the torus caller guarantees a4 = |dir|^4 > 0
  }
  // monic x^4 + b x^3 + c x^2 + d x + e
  const double coeffB = a3 / a4, coeffC = a2 / a4, coeffD = a1 / a4, coeffE = a0 / a4;
  // depress with x = y - b/4: y^4 + p y^2 + q y + r
  const double termP = coeffC - 3. * coeffB * coeffB / 8.;
  const double termQ = coeffD - coeffB * coeffC / 2. + coeffB * coeffB * coeffB / 8.;
  const double termR =
    coeffE - coeffB * coeffD / 4. + coeffB * coeffB * coeffC / 16. - 3. * coeffB * coeffB * coeffB * coeffB / 256.;
  const double shift = -coeffB / 4.;

  auto addQuadraticRoots = [&](double quadB, double quadC) {
    const double discriminant = quadB * quadB - 4. * quadC;
    if (discriminant < 0.) {
      return; // complex pair
    }
    const double sqrtDiscriminant = std::sqrt(discriminant);
    roots.push_back(shift + 0.5 * (-quadB - sqrtDiscriminant));
    roots.push_back(shift + 0.5 * (-quadB + sqrtDiscriminant));
  };

  if (std::abs(termQ) <= kTolerance) {
    // biquadratic y^4 + p y^2 + r = 0
    const double discriminant = termP * termP - 4. * termR;
    if (discriminant >= 0.) {
      const double sqrtDiscriminant = std::sqrt(discriminant);
      for (const double ySquared : {0.5 * (-termP + sqrtDiscriminant), 0.5 * (-termP - sqrtDiscriminant)}) {
        if (ySquared >= 0.) {
          const double y = std::sqrt(ySquared);
          roots.push_back(shift + y);
          roots.push_back(shift - y);
        }
      }
    }
  } else {
    // resolvent cubic m^3 + p m^2 + (p^2/4 - r) m - q^2/8 = 0; its largest real root is > 0
    const double cubicA2 = termP;
    const double cubicA1 = termP * termP / 4. - termR;
    const double cubicA0 = -termQ * termQ / 8.;
    const double cubicP = cubicA1 - cubicA2 * cubicA2 / 3.;
    const double cubicQ = 2. * cubicA2 * cubicA2 * cubicA2 / 27. - cubicA2 * cubicA1 / 3. + cubicA0;
    std::vector<double> cubicRoots;
    solveDepressedCubic(cubicP, cubicQ, cubicRoots);
    double resolvent = 0.;
    for (const double cubicRoot : cubicRoots) {
      resolvent = std::max(resolvent, cubicRoot - cubicA2 / 3.);
    }
    if (resolvent > kTolerance) {
      const double sqrtTwoResolvent = std::sqrt(2. * resolvent);
      const double linearTerm = sqrtTwoResolvent * termQ / (4. * resolvent);
      addQuadraticRoots(-sqrtTwoResolvent, termP / 2. + resolvent + linearTerm);
      addQuadraticRoots(sqrtTwoResolvent, termP / 2. + resolvent - linearTerm);
    }
  }

  // Newton polishing against the monic quartic tightens Ferrari's roots (and pulls a near-double
  // pair together so the parity clustering recognises it as a tangency)
  auto quartic = [&](double x) { return (((x + coeffB) * x + coeffC) * x + coeffD) * x + coeffE; };
  auto quarticDerivative = [&](double x) { return ((4. * x + 3. * coeffB) * x + 2. * coeffC) * x + coeffD; };
  for (double& root : roots) {
    for (int iteration = 0; iteration < 2; ++iteration) {
      const double derivative = quarticDerivative(root);
      if (std::abs(derivative) > kTolerance) {
        root -= quartic(root) / derivative;
      }
    }
  }
  return roots;
}

/// Kind of a 2D trimmed boundary curve.
enum class CurveKind { Line,   ///< straight line segment
                       Arc,    ///< circular arc
                       BSpline ///< clamped (rational) B-spline curve
};

/// One trimmed boundary curve segment in a surface's parametric (u, v) domain.
///
/// A single value type carries either a straight line segment or a circular arc, so a wire can
/// mix both kinds without heap allocation or virtual dispatch. Lines store their two endpoints;
/// arcs store a centre, radius and the start/end angles (in radians) — the signed sweep
/// (endAngle - startAngle) runs counter-clockwise when positive and clockwise when negative, and
/// a full circle is simply a +/-2pi sweep. All operations the trimming/winding machinery needs
/// are provided exactly for both kinds: endpoints, unit tangent, exact axis-aligned bounding-box
/// contribution, point projection, exact signed-area contribution (Green's theorem), and
/// horizontal-ray crossings for robust point-in-wire winding.
///
/// Keeping curves surface-independent lets planar, cylindrical, spherical and conical surfaces
/// reuse the same trimming model in their parametric domains.
struct Curve2D {
  CurveKind kind = CurveKind::Line;
  Vec2 lineStart;         ///< line: start point (unused for arcs)
  Vec2 lineEnd;           ///< line: end point (unused for arcs)
  Vec2 center;            ///< arc: circle centre (unused for lines)
  double radius = 0.;     ///< arc: circle radius
  double startAngle = 0.; ///< arc: start angle [rad]
  double endAngle = 0.;   ///< arc: end angle [rad] (sweep = endAngle - startAngle)

  /// \name B-spline data (kind == BSpline)
  /// A clamped, possibly rational B-spline curve in the (u, v) parameter domain. \a poles are
  /// the control points, \a weights the per-pole weights (empty ⇒ all one ⇒ non-rational) and
  /// \a knots the clamped flat knot vector of length poles.size() + degree + 1. A Bézier curve
  /// is the single-span special case. The curve parameter used elsewhere in Curve2D runs on
  /// [0, 1] across the whole knot domain [knots[degree], knots[nPoles]]. @{
  int degree = 0;
  std::vector<Vec2> poles;
  std::vector<double> weights;
  std::vector<double> knots;
  /// Lazily-filled flattened polyline of the B-spline (on-curve points, both ends included), so
  /// the per-point winding/distance queries used by point-in-wire classification and the numeric
  /// capacity integration do not re-run the adaptive de Boor sampling on every call. Cleared when
  /// the curve is reversed; safe to copy (the geometry it caches is unchanged by a copy).
  mutable std::vector<Vec2> bsplineCache;
  /// @}

  /// \name Loop-canonical endpoints
  /// The vertex values this curve's neighbours in a closed wire agree on. A wire's curves are
  /// imported with independently sampled endpoints, so curve i's end and curve i+1's start differ
  /// by up to kWireJoinTolerance; CurveWire::initialize therefore fixes one value per seam and
  /// hands it to both sides. The flattened polyline is built with these substituted at its ends,
  /// which is what makes winding and point-to-curve distance see the *same* boundary -- they used
  /// to disagree by the seam drift, a defect of its own hiding inside K5.
  /// @{
  Vec2 canonicalStart;
  Vec2 canonicalEnd;
  bool hasCanonicalEndpoints = false;

  void setCanonicalEndpoints(const Vec2& start, const Vec2& end)
  {
    canonicalStart = start;
    canonicalEnd = end;
    hasCanonicalEndpoints = true;
    bsplineCache.clear(); // the polyline carries them, so it has to be rebuilt
  }

  /// Where this curve begins and ends as far as the loop is concerned: the canonical seam vertex
  /// when a wire has fixed one, and the curve's own endpoint when it stands alone.
  Vec2 loopStart() const { return hasCanonicalEndpoints ? canonicalStart : startPoint(); }
  Vec2 loopEnd() const { return hasCanonicalEndpoints ? canonicalEnd : endPoint(); }
  /// @}

  static Curve2D makeLine(const Vec2& start, const Vec2& end)
  {
    Curve2D curve;
    curve.kind = CurveKind::Line;
    curve.lineStart = start;
    curve.lineEnd = end;
    return curve;
  }

  static Curve2D makeArc(const Vec2& arcCenter, double arcRadius, double arcStartAngle, double arcEndAngle)
  {
    Curve2D curve;
    curve.kind = CurveKind::Arc;
    curve.center = arcCenter;
    curve.radius = arcRadius;
    curve.startAngle = arcStartAngle;
    curve.endAngle = arcEndAngle;
    return curve;
  }

  /// Full circle as one arc curve (counter-clockwise unless \a clockwise is set).
  static Curve2D makeCircle(const Vec2& arcCenter, double arcRadius, bool clockwise = false)
  {
    return makeArc(arcCenter, arcRadius, 0., clockwise ? -kTwoPi : kTwoPi);
  }

  /// Clamped (rational) B-spline curve of degree \a splineDegree. \a splineWeights may be empty
  /// for a non-rational curve; \a splineKnots must be the clamped flat knot vector.
  static Curve2D makeBSpline(int splineDegree, std::vector<Vec2> splinePoles,
                             std::vector<double> splineWeights, std::vector<double> splineKnots)
  {
    Curve2D curve;
    curve.kind = CurveKind::BSpline;
    curve.degree = splineDegree;
    curve.poles = std::move(splinePoles);
    curve.weights = std::move(splineWeights);
    curve.knots = std::move(splineKnots);
    return curve;
  }

  bool isArc() const { return kind == CurveKind::Arc; }
  bool isBSpline() const { return kind == CurveKind::BSpline; }

  double sweep() const { return endAngle - startAngle; }

  /// \name B-spline evaluation helpers (kind == BSpline)
  /// @{
  double bsplineT0() const { return knots[degree]; }
  double bsplineT1() const { return knots[poles.size()]; }

  /// True when the knot vector is *clamped*: the first and last degree+1 knots each coincide.
  /// Only then does the curve interpolate its first and last pole, which is the shortcut
  /// startPoint()/endPoint() take. An unclamped or periodic knot vector -- exactly what OCC
  /// writes for a tube-tube intersection curve before SetNotPeriodic -- starts and ends somewhere
  /// strictly inside the control polygon, so the shortcut would return a point that is not on the
  /// curve at all: the wire then reads as Open (its edges no longer meet) and the face is
  /// rejected, or the off-curve endpoint corrupts the winding classification. See
  /// CodeReview_Fable.md K1.
  bool bsplineIsClamped() const
  {
    const size_t lastKnot = knots.size() - 1;
    for (int offset = 1; offset <= degree; ++offset) {
      if (std::abs(knots[offset] - knots[0]) > kTolerance ||
          std::abs(knots[lastKnot - offset] - knots[lastKnot]) > kTolerance) {
        return false;
      }
    }
    return true;
  }

  /// True if the curve carries non-unit weights (a rational B-spline).
  bool bsplineRational() const
  {
    for (double weight : weights) {
      if (std::abs(weight - 1.) > kTolerance) {
        return true;
      }
    }
    return false;
  }

  /// Knot span index of parameter \a knotValue for the clamped knot vector.
  int bsplineSpan(double knotValue) const
  {
    const int lastPole = static_cast<int>(poles.size()) - 1;
    if (knotValue >= knots[lastPole + 1]) {
      return lastPole;
    }
    if (knotValue <= knots[degree]) {
      return degree;
    }
    int low = degree;
    int high = lastPole + 1;
    int mid = (low + high) / 2;
    while (knotValue < knots[mid] || knotValue >= knots[mid + 1]) {
      if (knotValue < knots[mid]) {
        high = mid;
      } else {
        low = mid;
      }
      mid = (low + high) / 2;
    }
    return mid;
  }

  /// Non-zero degree-p basis functions \a basis and their first derivatives \a basisDeriv at
  /// \a knotValue in the given \a span (Piegl & Tiller, "The NURBS Book", DersBasisFuns, first
  /// derivative only).
  void bsplineBasis(int span, double knotValue, std::vector<double>& basis,
                    std::vector<double>& basisDeriv) const
  {
    const int p = degree;
    std::vector<std::vector<double>> ndu(p + 1, std::vector<double>(p + 1, 0.));
    std::vector<double> left(p + 1, 0.);
    std::vector<double> right(p + 1, 0.);
    ndu[0][0] = 1.;
    for (int j = 1; j <= p; ++j) {
      left[j] = knotValue - knots[span + 1 - j];
      right[j] = knots[span + j] - knotValue;
      double saved = 0.;
      for (int r = 0; r < j; ++r) {
        ndu[j][r] = right[r + 1] + left[j - r];
        const double temp = ndu[r][j - 1] / ndu[j][r];
        ndu[r][j] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    basis.assign(p + 1, 0.);
    basisDeriv.assign(p + 1, 0.);
    for (int j = 0; j <= p; ++j) {
      basis[j] = ndu[j][p];
    }
    // first derivative (specialization of DersBasisFuns for the k = 1 term)
    for (int r = 0; r <= p; ++r) {
      double d = 0.;
      const int pk = p - 1;
      if (r >= 1) {
        d += (1. / ndu[pk + 1][r - 1]) * ndu[r - 1][pk];
      }
      if (r <= pk) {
        d += (-1. / ndu[pk + 1][r]) * ndu[r][pk];
      }
      basisDeriv[r] = d * p;
    }
  }

  /// Evaluate the (rational) B-spline point \a pointOut and its knot-parameter derivative
  /// \a derivativeOut at knot parameter \a knotValue.
  void bsplineEval(double knotValue, Vec2& pointOut, Vec2& derivativeOut) const
  {
    const int p = degree;
    const int span = bsplineSpan(knotValue);
    std::vector<double> basis;
    std::vector<double> basisDeriv;
    bsplineBasis(span, knotValue, basis, basisDeriv);
    Vec2 weightedSum{0., 0.};
    Vec2 weightedDeriv{0., 0.};
    double weightTotal = 0.;
    double weightDeriv = 0.;
    for (int j = 0; j <= p; ++j) {
      const int idx = span - p + j;
      const double weight = weights.empty() ? 1. : weights[idx];
      weightedSum.uCoord += basis[j] * weight * poles[idx].uCoord;
      weightedSum.vCoord += basis[j] * weight * poles[idx].vCoord;
      weightTotal += basis[j] * weight;
      weightedDeriv.uCoord += basisDeriv[j] * weight * poles[idx].uCoord;
      weightedDeriv.vCoord += basisDeriv[j] * weight * poles[idx].vCoord;
      weightDeriv += basisDeriv[j] * weight;
    }
    const double invWeight = (std::abs(weightTotal) > kTolerance) ? 1. / weightTotal : 0.;
    pointOut = {weightedSum.uCoord * invWeight, weightedSum.vCoord * invWeight};
    derivativeOut = {(weightedDeriv.uCoord * weightTotal - weightedSum.uCoord * weightDeriv) * invWeight * invWeight,
                     (weightedDeriv.vCoord * weightTotal - weightedSum.vCoord * weightDeriv) * invWeight * invWeight};
  }

  /// B-spline point at curve parameter \a parameter in [0, 1].
  Vec2 bsplinePointAt(double parameter) const
  {
    const double knotValue = bsplineT0() + parameter * (bsplineT1() - bsplineT0());
    Vec2 point;
    Vec2 derivative;
    bsplineEval(knotValue, point, derivative);
    return point;
  }

  /// Adaptively sample the B-spline into an ordered on-curve polyline (including both ends),
  /// subdividing the knot domain until each chord deviates from the curve by less than
  /// sqrt(\a flatnessSq) (or the chord itself is already below the flatness scale, or the depth
  /// cap is reached). Shared by winding, distance, area-independent visualization and mesh. The
  /// bounded depth keeps the sample count small even for a near-degenerate or high-curvature
  /// span; the residual chord error stays well within the navigation length tolerances.
  void bsplineSampleInto(std::vector<Vec2>& samples, double flatnessSq = kBSplineFlatnessSq,
                         int maxDepth = 16) const
  {
    const double t0 = bsplineT0();
    const double t1 = bsplineT1();
    Vec2 startPointValue;
    Vec2 endPointValue;
    Vec2 unusedDerivative;
    bsplineEval(t0, startPointValue, unusedDerivative);
    bsplineEval(t1, endPointValue, unusedDerivative);
    samples.push_back(startPointValue);
    bsplineSampleRecursive(t0, t1, startPointValue, endPointValue, flatnessSq, maxDepth, samples);
  }

  /// Does the knot interval (\a lowT, \a highT) contain a knot strictly inside it?
  ///
  /// A B-spline is a single polynomial piece only *within* one knot span. A flatness test that
  /// straddles a knot is therefore testing a curve its own model does not describe, and two spans
  /// that mirror each other defeat it outright. So an interval that still contains an interior
  /// knot is never called flat, whatever the probes say.
  bool spansInteriorKnot(double lowT, double highT) const
  {
    // a clamped knot vector repeats its ends degree+1 times, so the interior knots are the
    // entries [degree + 1, poles.size()); a single-span (Bezier) curve has none
    const size_t firstInterior = static_cast<size_t>(degree) + 1;
    const size_t endInterior = std::min(poles.size(), knots.size());
    if (firstInterior >= endInterior) {
      return false;
    }
    const auto begin = knots.begin() + static_cast<std::ptrdiff_t>(firstInterior);
    const auto end = knots.begin() + static_cast<std::ptrdiff_t>(endInterior);
    const auto firstAbove = std::upper_bound(begin, end, lowT);
    return firstAbove != end && *firstAbove < highT;
  }

  void bsplineSampleRecursive(double t0, double t1, const Vec2& p0, const Vec2& p1, double flatnessSq,
                              int depth, std::vector<Vec2>& samples) const
  {
    const double tMid = 0.5 * (t0 + t1);
    Vec2 midPoint;
    Vec2 unusedDerivative;
    bsplineEval(tMid, midPoint, unusedDerivative);
    // A short chord alone must NOT end the recursion: a *closed* curve has p0 == p1 exactly, so
    // its chord is zero however much curve lies between them -- a full circle written as one
    // periodic B-spline (which is what a CAD kernel emits for a tube-tube intersection) would
    // otherwise flatten to two coincident points and vanish. When the chord is degenerate the only
    // meaningful test is the distance to that single point; pointSegmentDistanceSq would compare
    // against a zero-length segment.
    const bool degenerateChord = surface::distanceSq(p0, p1) <= flatnessSq;
    const auto deviationSq = [&](const Vec2& point) {
      return degenerateChord ? surface::distanceSq(point, p0) : pointSegmentDistanceSq(point, p0, p1);
    };
    // Three interior probes, not one. A single midpoint probe is blind to every curve that is
    // symmetric about its own parameter midpoint -- and that is exactly what a CAD kernel emits
    // for a tube-tube intersection seen from the *thin* tube: on Bagger/BoomCylinderInner the
    // midpoint sits 2e-7 from the chord while the curve bulges 0.16 away from it. That curve was
    // declared flat at the top level and replaced by its chord, taking the whole junction rim of
    // six Bagger parts with it (TolerancePolicy.md section 13, finding K4).
    double flatness = deviationSq(midPoint);
    for (const double fraction : {0.25, 0.75}) {
      Vec2 probePoint;
      bsplineEval(t0 + (t1 - t0) * fraction, probePoint, unusedDerivative);
      flatness = std::max(flatness, deviationSq(probePoint));
    }
    if (depth <= 0 || (flatness <= flatnessSq && !spansInteriorKnot(t0, t1))) {
      samples.push_back(p1);
      return;
    }
    bsplineSampleRecursive(t0, tMid, p0, midPoint, flatnessSq, depth - 1, samples);
    bsplineSampleRecursive(tMid, t1, midPoint, p1, flatnessSq, depth - 1, samples);
  }

  /// The cached flattened polyline (see \a bsplineCache), computed once on first use. All the
  /// per-point B-spline queries go through this so a caller that hits the same curve many times
  /// (numeric capacity integration, dense point-in-wire tests) pays the de Boor sampling only once.
  const std::vector<Vec2>& bsplineSamples() const
  {
    if (bsplineCache.empty()) {
      bsplineSampleInto(bsplineCache);
      // One polyline, and it is the canonical one: every query that goes through this -- winding,
      // closest point, the numeric capacity integrand -- then sees the identical boundary. The
      // substitution is at the ends only; the interior samples are on the curve either way.
      if (hasCanonicalEndpoints && bsplineCache.size() >= 2) {
        bsplineCache.front() = canonicalStart;
        bsplineCache.back() = canonicalEnd;
      }
    }
    return bsplineCache;
  }
  /// @}

  /// Basic structural validity (finite data, positive radius for arcs, well-formed clamped knot
  /// vector for B-splines).
  bool valid() const
  {
    if (kind == CurveKind::Line) {
      return finite(lineStart) && finite(lineEnd);
    }
    if (kind == CurveKind::Arc) {
      return finite(center) && std::isfinite(radius) && radius > kTolerance && std::isfinite(startAngle) &&
             std::isfinite(endAngle);
    }
    // B-spline
    const int nPoles = static_cast<int>(poles.size());
    if (degree < 1 || nPoles < degree + 1) {
      return false;
    }
    if (static_cast<int>(knots.size()) != nPoles + degree + 1) {
      return false;
    }
    if (!weights.empty() && static_cast<int>(weights.size()) != nPoles) {
      return false;
    }
    for (const auto& pole : poles) {
      if (!finite(pole)) {
        return false;
      }
    }
    for (double weight : weights) {
      if (!std::isfinite(weight) || weight <= kTolerance) {
        return false;
      }
    }
    for (size_t index = 1; index < knots.size(); ++index) {
      if (!std::isfinite(knots[index]) || knots[index] < knots[index - 1] - kTolerance) {
        return false;
      }
    }
    return bsplineT1() - bsplineT0() > kTolerance;
  }

  Vec2 pointAtAngle(double angle) const
  {
    return {center.uCoord + radius * std::cos(angle), center.vCoord + radius * std::sin(angle)};
  }

  /// Point at curve parameter \a parameter in [0, 1] (0 at the start, 1 at the end).
  Vec2 pointAt(double parameter) const
  {
    if (kind == CurveKind::Line) {
      return {lineStart.uCoord + parameter * (lineEnd.uCoord - lineStart.uCoord),
              lineStart.vCoord + parameter * (lineEnd.vCoord - lineStart.vCoord)};
    }
    if (kind == CurveKind::BSpline) {
      return bsplinePointAt(parameter);
    }
    return pointAtAngle(startAngle + parameter * sweep());
  }

  Vec2 startPoint() const
  {
    if (kind == CurveKind::Line) {
      return lineStart;
    }
    if (kind == CurveKind::BSpline) {
      // a clamped knot vector interpolates its first pole exactly, which is both cheaper and
      // more accurate than evaluating; anything else has to be evaluated (K1)
      return bsplineIsClamped() ? poles.front() : bsplinePointAt(0.);
    }
    return pointAtAngle(startAngle);
  }
  Vec2 endPoint() const
  {
    if (kind == CurveKind::Line) {
      return lineEnd;
    }
    if (kind == CurveKind::BSpline) {
      return bsplineIsClamped() ? poles.back() : bsplinePointAt(1.);
    }
    return pointAtAngle(endAngle);
  }

  /// dC/dt at curve parameter \a parameter in [0, 1], unnormalized. tangentAt() is this
  /// normalized; a contour integral needs the length, because dv = v'(t) dt is what it integrates
  /// against.
  Vec2 derivativeAt(double parameter) const
  {
    if (kind == CurveKind::Line) {
      return {lineEnd.uCoord - lineStart.uCoord, lineEnd.vCoord - lineStart.vCoord};
    }
    if (kind == CurveKind::BSpline) {
      const double span = bsplineT1() - bsplineT0();
      Vec2 point;
      Vec2 derivative;
      bsplineEval(bsplineT0() + parameter * span, point, derivative);
      return {derivative.uCoord * span, derivative.vCoord * span};
    }
    const double angle = startAngle + parameter * sweep();
    return {-radius * std::sin(angle) * sweep(), radius * std::cos(angle) * sweep()};
  }

  /// Unit tangent at parameter \a parameter, pointing in the direction of increasing parameter.
  Vec2 tangentAt(double parameter) const
  {
    if (kind == CurveKind::Line) {
      const Vec2 delta{lineEnd.uCoord - lineStart.uCoord, lineEnd.vCoord - lineStart.vCoord};
      const double length = std::sqrt(delta.uCoord * delta.uCoord + delta.vCoord * delta.vCoord);
      if (length <= kTolerance) {
        return {0., 0.};
      }
      return {delta.uCoord / length, delta.vCoord / length};
    }
    if (kind == CurveKind::BSpline) {
      // dC/dt scaled by the positive constant dt/ds, so the normalized direction is unchanged
      const double knotValue = bsplineT0() + parameter * (bsplineT1() - bsplineT0());
      Vec2 point;
      Vec2 derivative;
      bsplineEval(knotValue, point, derivative);
      const double length = std::sqrt(derivative.uCoord * derivative.uCoord + derivative.vCoord * derivative.vCoord);
      if (length <= kTolerance) {
        return {0., 0.};
      }
      return {derivative.uCoord / length, derivative.vCoord / length};
    }
    const double angle = startAngle + parameter * sweep();
    const double direction = sweep() >= 0. ? 1. : -1.;
    return {-direction * std::sin(angle), direction * std::cos(angle)};
  }

  /// True if \a angle lies within the arc's angular sweep (accounting for direction and wrap).
  bool angleInSweep(double angle) const
  {
    const double totalSweep = sweep();
    const double magnitude = std::abs(totalSweep);
    if (magnitude >= kTwoPi - kTolerance) {
      return true; // full circle
    }
    double delta = (totalSweep >= 0.) ? (angle - startAngle) : (startAngle - angle);
    delta -= kTwoPi * std::floor(delta / kTwoPi); // wrap into [0, 2pi)
    return delta <= magnitude + kTolerance;
  }

  /// Map an angle known to lie within the sweep to a clamped parameter in [0, 1].
  double angleParameter(double angle) const
  {
    const double totalSweep = sweep();
    if (std::abs(totalSweep) <= kTolerance) {
      return 0.;
    }
    double delta = (totalSweep >= 0.) ? (angle - startAngle) : (startAngle - angle);
    delta -= kTwoPi * std::floor(delta / kTwoPi);
    return std::max(0., std::min(1., delta / std::abs(totalSweep)));
  }

  /// Accumulate this curve's exact extent into a parametric axis-aligned bounding box.
  void extendBounds(Vec2& lower, Vec2& upper) const
  {
    auto include = [&](const Vec2& point) {
      lower.uCoord = std::min(lower.uCoord, point.uCoord);
      lower.vCoord = std::min(lower.vCoord, point.vCoord);
      upper.uCoord = std::max(upper.uCoord, point.uCoord);
      upper.vCoord = std::max(upper.vCoord, point.vCoord);
    };
    if (kind == CurveKind::BSpline) {
      // the control-point convex hull contains the curve, so its box is a conservative (exact
      // upper bound) parametric AABB — consistent with the BVH's conservative-box philosophy
      for (const auto& pole : poles) {
        include(pole);
      }
      return;
    }
    includeAnalyticExtremes(include);
  }

  /// As extendBounds, but measured on the curve rather than on the control polygon: a B-spline
  /// contributes its sampled polyline instead of its pole hull. Tight rather than conservative,
  /// for the callers that reject a wire when the box is too large (see
  /// CurveWire::tightParametricBounds).
  void extendTightBounds(Vec2& lower, Vec2& upper) const
  {
    auto include = [&](const Vec2& point) {
      lower.uCoord = std::min(lower.uCoord, point.uCoord);
      lower.vCoord = std::min(lower.vCoord, point.vCoord);
      upper.uCoord = std::max(upper.uCoord, point.uCoord);
      upper.vCoord = std::max(upper.vCoord, point.vCoord);
    };
    if (kind == CurveKind::BSpline) {
      for (const auto& sample : bsplineSamples()) {
        include(sample);
      }
      return;
    }
    includeAnalyticExtremes(include);
  }

  /// Endpoints plus, for an arc, the axis-extreme points inside the sweep: the exact extent of a
  /// line or arc. Shared by the conservative and the tight bounds, which differ only in how they
  /// treat a B-spline.
  template <typename Include>
  void includeAnalyticExtremes(const Include& include) const
  {
    include(startPoint());
    include(endPoint());
    if (kind == CurveKind::Arc) {
      // include the axis-extreme points (angles 0, pi/2, pi, 3pi/2) that fall within the sweep
      const double cardinalAngles[4] = {0., kHalfPi, kPi, 3. * kHalfPi};
      for (double cardinal : cardinalAngles) {
        if (angleInSweep(cardinal)) {
          include(pointAtAngle(cardinal));
        }
      }
    }
  }

  /// Closest point on the curve to \a point, returning the clamped parameter in \a parameter.
  Vec2 closestPoint(const Vec2& point, double& parameter) const
  {
    if (kind == CurveKind::BSpline) {
      // Distance to the cached flattened polyline: the nearest projection onto a chord. The
      // polyline is flat to the sampling tolerance, so this is accurate to ~that tolerance, which
      // is what classify's boundary band and the (deliberately non-exact) B-spline Safety need,
      // and it costs no de Boor evaluation per call (critical for dense capacity integration).
      const auto& polyline = bsplineSamples();
      if (polyline.size() < 2) {
        parameter = 0.;
        return startPoint();
      }
      const int segments = static_cast<int>(polyline.size()) - 1;
      double bestDistanceSq = std::numeric_limits<double>::infinity();
      Vec2 bestPoint = polyline.front();
      double bestParameter = 0.;
      for (int index = 0; index < segments; ++index) {
        const Vec2 segmentStart = polyline[index];
        const Vec2 segmentVector = polyline[index + 1] - segmentStart;
        const double segmentLengthSq =
          segmentVector.uCoord * segmentVector.uCoord + segmentVector.vCoord * segmentVector.vCoord;
        double projection = 0.;
        if (segmentLengthSq > kToleranceSq) {
          projection = ((point.uCoord - segmentStart.uCoord) * segmentVector.uCoord +
                        (point.vCoord - segmentStart.vCoord) * segmentVector.vCoord) /
                       segmentLengthSq;
          projection = std::max(0., std::min(1., projection));
        }
        const Vec2 candidate{segmentStart.uCoord + projection * segmentVector.uCoord,
                             segmentStart.vCoord + projection * segmentVector.vCoord};
        const double candidateDistanceSq = surface::distanceSq(point, candidate);
        if (candidateDistanceSq < bestDistanceSq) {
          bestDistanceSq = candidateDistanceSq;
          bestPoint = candidate;
          bestParameter = (index + projection) / segments;
        }
      }
      parameter = bestParameter;
      return bestPoint;
    }
    if (kind == CurveKind::Line) {
      const Vec2 segment{lineEnd.uCoord - lineStart.uCoord, lineEnd.vCoord - lineStart.vCoord};
      const double lengthSq = segment.uCoord * segment.uCoord + segment.vCoord * segment.vCoord;
      if (lengthSq <= kToleranceSq) {
        parameter = 0.;
        return lineStart;
      }
      const double projection = ((point.uCoord - lineStart.uCoord) * segment.uCoord +
                                 (point.vCoord - lineStart.vCoord) * segment.vCoord) /
                                lengthSq;
      parameter = std::max(0., std::min(1., projection));
      return {lineStart.uCoord + parameter * segment.uCoord, lineStart.vCoord + parameter * segment.vCoord};
    }
    // arc: project radially onto the circle, then clamp the angle to the sweep
    const double deltaU = point.uCoord - center.uCoord;
    const double deltaV = point.vCoord - center.vCoord;
    if (deltaU * deltaU + deltaV * deltaV <= kToleranceSq) {
      parameter = 0.; // point at the centre: every arc point is equidistant
      return startPoint();
    }
    const double angle = std::atan2(deltaV, deltaU);
    if (angleInSweep(angle)) {
      parameter = angleParameter(angle);
      return pointAtAngle(angle);
    }
    const Vec2 startCandidate = startPoint();
    const Vec2 endCandidate = endPoint();
    if (surface::distanceSq(point, startCandidate) <= surface::distanceSq(point, endCandidate)) {
      parameter = 0.;
      return startCandidate;
    }
    parameter = 1.;
    return endCandidate;
  }

  /// Squared distance from \a point to the curve.
  double distanceSq(const Vec2& point) const
  {
    double parameter = 0.;
    return surface::distanceSq(point, closestPoint(point, parameter));
  }

  /// Exact contribution of this directed curve to the enclosed signed area,
  /// i.e. (1/2) * integral of (u dv - v du) along the curve (Green's theorem).
  double signedAreaContribution() const
  {
    if (kind == CurveKind::Line) {
      return 0.5 * (lineStart.uCoord * lineEnd.vCoord - lineEnd.uCoord * lineStart.vCoord);
    }
    if (kind == CurveKind::BSpline) {
      // (1/2) integral of (u dv - v du) = (1/2) integral of (u v' - v u') dt, evaluated per knot
      // span by Gauss-Legendre. For a non-rational span the integrand is a degree 2p-1 polynomial,
      // exactly integrated with p+1 nodes; a rational span is only approximated (higher order), so
      // a bspline-trimmed face reports capacityIsExact() == false.
      const int p = degree;
      const int order = bsplineRational() ? std::max(2 * p + 2, 8) : (p + 1);
      std::vector<double> nodes;
      std::vector<double> nodeWeights;
      gaussLegendre(order, nodes, nodeWeights);
      double area = 0.;
      const int lastSpan = static_cast<int>(poles.size()) - 1;
      for (int spanIndex = p; spanIndex <= lastSpan; ++spanIndex) {
        const double spanLow = knots[spanIndex];
        const double spanHigh = knots[spanIndex + 1];
        const double halfSpan = 0.5 * (spanHigh - spanLow);
        if (halfSpan <= kTolerance) {
          continue;
        }
        const double spanMid = 0.5 * (spanLow + spanHigh);
        for (int nodeIndex = 0; nodeIndex < order; ++nodeIndex) {
          const double knotValue = spanMid + halfSpan * nodes[nodeIndex];
          Vec2 point;
          Vec2 derivative;
          bsplineEval(knotValue, point, derivative);
          area += 0.5 * (point.uCoord * derivative.vCoord - point.vCoord * derivative.uCoord) *
                  nodeWeights[nodeIndex] * halfSpan;
        }
      }
      return area;
    }
    return 0.5 * (radius * center.uCoord * (std::sin(endAngle) - std::sin(startAngle)) -
                  radius * center.vCoord * (std::cos(endAngle) - std::cos(startAngle)) +
                  radius * radius * (endAngle - startAngle));
  }

  /// How far this curve's *representation* can sit from the curve it stands for, in parametric
  /// units. Lines and arcs are held exactly and answer both winding and distance analytically, so
  /// they contribute nothing. A B-spline enters navigation as its flattened polyline, and is
  /// therefore only as good as kBSplineFlatness. Anything that classifies a point against the
  /// boundary has to widen its band by this or it is measuring noise (finding K5).
  double representationTolerance() const { return kind == CurveKind::BSpline ? kBSplineFlatness : 0.; }

  /// Number of times a horizontal ray from \a point towards +u crosses this curve, using a
  /// half-open convention so that shared endpoints and tangent extrema are not double-counted.
  int rightwardCrossings(const Vec2& point) const
  {
    return rightwardCrossings(point, startPoint(), endPoint());
  }

  /// Crossing count with caller-supplied canonical endpoints. In a closed loop the same vertex
  /// value must terminate one curve and start the next, otherwise floating-point drift at the
  /// seam (e.g. sin(2pi) != 0 for a one-arc full circle, or a within-tolerance closure gap)
  /// breaks the half-open above/below convention for scanlines passing exactly through it.
  /// CurveWire::classify therefore passes each curve's start point and the next curve's start
  /// point as the canonical endpoints.
  int rightwardCrossings(const Vec2& point, const Vec2& canonicalStart, const Vec2& canonicalEnd) const
  {
    auto segmentCrossing = [&](const Vec2& first, const Vec2& second, double exactIntersectU) {
      const bool firstAbove = first.vCoord > point.vCoord;
      const bool secondAbove = second.vCoord > point.vCoord;
      if (firstAbove == secondAbove) {
        return false;
      }
      return point.uCoord < exactIntersectU;
    };

    if (kind == CurveKind::Line) {
      const bool firstAbove = canonicalStart.vCoord > point.vCoord;
      const bool secondAbove = canonicalEnd.vCoord > point.vCoord;
      if (firstAbove == secondAbove) {
        return 0;
      }
      const double intersectU = canonicalStart.uCoord + (point.vCoord - canonicalStart.vCoord) *
                                                          (canonicalEnd.uCoord - canonicalStart.uCoord) /
                                                          (canonicalEnd.vCoord - canonicalStart.vCoord);
      return (point.uCoord < intersectU) ? 1 : 0;
    }

    if (kind == CurveKind::BSpline) {
      // Adaptively flatten the curve into an on-curve polyline, then apply the same half-open
      // segment-crossing rule as for lines. The canonical shared endpoints replace the first and
      // last polyline vertices so seams between curves stay consistent. Flattening below the chord
      // tolerance means each segment is v-monotonic to that tolerance, so a true tangency (both
      // endpoints on the same side) contributes an even count and a transversal crossing an odd
      // one, exactly as the arc's v-monotonic split does.
      const auto& polyline = bsplineSamples();
      if (polyline.size() < 2) {
        return 0;
      }
      int crossings = 0;
      for (size_t index = 0; index + 1 < polyline.size(); ++index) {
        // No substitution here: the polyline already ends on the loop-canonical vertices (see
        // setCanonicalEndpoints), so this is the same boundary closestPoint measures against.
        const Vec2 first = polyline[index];
        const Vec2 second = polyline[index + 1];
        const bool firstAbove = first.vCoord > point.vCoord;
        const bool secondAbove = second.vCoord > point.vCoord;
        if (firstAbove == secondAbove) {
          continue;
        }
        const double intersectU =
          first.uCoord + (point.vCoord - first.vCoord) * (second.uCoord - first.uCoord) /
                           (second.vCoord - first.vCoord);
        if (point.uCoord < intersectU) {
          ++crossings;
        }
      }
      return crossings;
    }

    // Split the arc into v-monotonic sub-arcs at its top/bottom extreme angles (cos(theta) = 0,
    // i.e. theta = pi/2 + k*pi). On each sub-arc cos(theta) keeps a constant sign, so the exact u
    // at v = point.v is centre.u +/- r*sqrt(1 - sin^2), with the sign taken from the sub-arc.
    const double totalSweep = sweep();
    if (std::abs(totalSweep) <= kTolerance || radius <= kTolerance) {
      return 0;
    }
    std::array<double, 8> breakParameters{};
    int breakCount = 0;
    breakParameters[breakCount++] = 0.;
    const double lowAngle = std::min(startAngle, endAngle);
    const double highAngle = std::max(startAngle, endAngle);
    const int firstK = static_cast<int>(std::floor((lowAngle - kHalfPi) / kPi)) - 1;
    const int lastK = static_cast<int>(std::ceil((highAngle - kHalfPi) / kPi)) + 1;
    for (int k = firstK; k <= lastK && breakCount < 7; ++k) {
      const double extremeAngle = kHalfPi + k * kPi;
      if (extremeAngle <= lowAngle + kTolerance || extremeAngle >= highAngle - kTolerance) {
        continue;
      }
      const double extremeParameter = (extremeAngle - startAngle) / totalSweep;
      if (extremeParameter > kTolerance && extremeParameter < 1. - kTolerance) {
        breakParameters[breakCount++] = extremeParameter;
      }
    }
    breakParameters[breakCount++] = 1.;
    std::sort(breakParameters.begin(), breakParameters.begin() + breakCount);

    double ratio = (point.vCoord - center.vCoord) / radius;
    ratio = std::max(-1., std::min(1., ratio));
    const double cosMagnitude = std::sqrt(std::max(0., 1. - ratio * ratio));

    int crossings = 0;
    for (int index = 0; index + 1 < breakCount; ++index) {
      const Vec2 subStart = (index == 0) ? canonicalStart : pointAt(breakParameters[index]);
      const Vec2 subEnd = (index + 2 == breakCount) ? canonicalEnd : pointAt(breakParameters[index + 1]);
      const double midAngle = startAngle + 0.5 * (breakParameters[index] + breakParameters[index + 1]) * totalSweep;
      const double cosSign = std::cos(midAngle) >= 0. ? 1. : -1.;
      const double intersectU = center.uCoord + cosSign * radius * cosMagnitude;
      if (segmentCrossing(subStart, subEnd, intersectU)) {
        ++crossings;
      }
    }
    return crossings;
  }

  /// Reverse the curve's direction in place (start <-> end), keeping the same geometric image.
  void reverseInPlace()
  {
    if (hasCanonicalEndpoints) {
      std::swap(canonicalStart, canonicalEnd);
      bsplineCache.clear();
    }
    if (kind == CurveKind::Line) {
      std::swap(lineStart, lineEnd);
    } else if (kind == CurveKind::Arc) {
      std::swap(startAngle, endAngle);
    } else {
      // B-spline: reverse the poles/weights and complement the knot vector about its span so the
      // parametrization runs the other way (knots stay non-decreasing and clamped).
      std::reverse(poles.begin(), poles.end());
      if (!weights.empty()) {
        std::reverse(weights.begin(), weights.end());
      }
      const double knotSum = knots.front() + knots.back();
      std::vector<double> reversedKnots(knots.size());
      for (size_t index = 0; index < knots.size(); ++index) {
        reversedKnots[index] = knotSum - knots[knots.size() - 1 - index];
      }
      knots = std::move(reversedKnots);
      bsplineCache.clear(); // geometry order changed; recompute lazily
    }
  }
};

/// One closed, oriented boundary loop built from Curve2D segments (lines and/or arcs).
///
/// This is the curved counterpart to SurfaceWire: it keeps the same role/status/orientation
/// conventions (outer loops wind counter-clockwise, inner/hole loops clockwise) but supports
/// exact circular boundaries. It is deliberately surface-independent so every analytic surface
/// can reuse it for parametric-domain trimming and winding.
struct CurveWire {
  std::vector<Curve2D> curves;
  WireRole role = WireRole::Outer;

  /// Build and validate the wire from an ordered, closed list of curves. \a metric turns the join
  /// gaps into 3D lengths so that kWireJoinTolerance means one distance on every surface; the form
  /// is evaluated at each join in turn, because it is not constant over the domain.
  bool initialize(const std::vector<Curve2D>& inputCurves, WireRole wireRole, WireStatus& status,
                  const ParametricMetric& metric = {})
  {
    role = wireRole;
    curves = inputCurves;

    if (curves.empty()) {
      status = WireStatus::TooFewVertices;
      return false;
    }
    for (size_t index = 0; index < curves.size(); ++index) {
      if (!curves[index].valid()) {
        status = WireStatus::NonFinite;
        return false;
      }
      const Vec2 currentEnd = curves[index].endPoint();
      const Vec2 nextStart = curves[(index + 1) % curves.size()].startPoint();
      if (metric.distanceSq(currentEnd, nextStart) > kWireJoinToleranceSq) {
        status = WireStatus::Open;
        return false;
      }
    }

    // Fix one vertex value per seam and give it to both curves that meet there. Every curve then
    // begins exactly where its predecessor ends, which is what lets winding and distance measure
    // against one polyline instead of two that differ by the join drift (K5).
    for (size_t index = 0; index < curves.size(); ++index) {
      curves[index].setCanonicalEndpoints(curves[index].startPoint(),
                                          curves[(index + 1) % curves.size()].startPoint());
    }

    const double area = signedArea();
    if (std::abs(area) <= kAreaTolerance) {
      status = WireStatus::ZeroArea;
      return false;
    }

    const bool wantPositiveArea = (role == WireRole::Outer);
    if ((area > 0.) != wantPositiveArea) {
      reverse();
      status = WireStatus::Reversed;
      return true;
    }
    status = WireStatus::Valid;
    return true;
  }

  /// The widest gap between this loop's representation and the boundary it stands for, in
  /// parametric units: the largest representationTolerance() over its curves. Zero for a loop of
  /// lines and arcs, which are exact.
  double representationTolerance() const
  {
    double worst = 0.;
    for (const auto& curve : curves) {
      worst = std::max(worst, curve.representationTolerance());
    }
    return worst;
  }

  /// Reverse the loop orientation in place (order and per-curve direction).
  void reverse()
  {
    std::reverse(curves.begin(), curves.end());
    for (auto& curve : curves) {
      curve.reverseInPlace();
    }
  }

  /// True if any curve of the loop is a B-spline (whose trimmed-face capacity is only numerically
  /// integrated, so the owning surface must report capacityIsExact() == false).
  bool hasBSpline() const
  {
    for (const auto& curve : curves) {
      if (curve.kind == CurveKind::BSpline) {
        return true;
      }
    }
    return false;
  }

  /// Exact signed area enclosed by the loop (positive when counter-clockwise).
  double signedArea() const
  {
    double area = 0.;
    for (const auto& curve : curves) {
      area += curve.signedAreaContribution();
    }
    return area;
  }

  /// Accumulate the loop's conservative extent into a parametric axis-aligned bounding box. For
  /// B-spline curves this is the control-point hull, which contains the curve but can be much
  /// larger than it -- the right trade for a BVH box, the wrong one for deciding whether the loop
  /// fits inside a full turn (see tightParametricBounds).
  void parametricBounds(Vec2& lower, Vec2& upper) const
  {
    for (const auto& curve : curves) {
      curve.extendBounds(lower, upper);
    }
  }

  /// The loop's extent measured on the curves themselves: exact for lines and arcs, and from the
  /// adaptively sampled polyline for B-splines. Never larger than parametricBounds() and usually
  /// smaller, so it is the bound to use for any test that *rejects* a wire for being too big --
  /// a conservative box turns "the curve fits" into "the control polygon fits", which is a
  /// different and much stronger demand (CodeReview_Fable.md K2).
  void tightParametricBounds(Vec2& lower, Vec2& upper) const
  {
    for (const auto& curve : curves) {
      curve.extendTightBounds(lower, upper);
    }
  }

  /// The half-width of the on-boundary band at \a point, in parametric units.
  ///
  /// A point classified Boundary is one whose side the representation cannot decide. That has to
  /// be set from what the representation actually is, not from optimism: a B-spline trim is a
  /// polyline flattened to kBSplineFlatness, so a 1e-9 band around it was measuring noise and the
  /// Boundary state was effectively unreachable for every B-spline trim, with in/out flipping
  /// arbitrarily within the flatness of the true curve (finding K5).
  ///
  /// The floor is kTolerance expressed as a *parametric* separation, since kTolerance is a length
  /// in cm and this domain is not: dividing by the metric's largest scale is the conservative
  /// direction (it widens the band where the surface is stretched). A degenerate scale -- a
  /// sphere pole, a cone apex -- leaves only the representation term, which is correct, because
  /// there every parametric separation spans no distance and nothing can be resolved by it.
  double boundaryBand(const Vec2& point, const ParametricMetric& metric) const
  {
    const double scale = metric.maxScale(point);
    const double lengthFloor = scale > kTolerance ? kTolerance / scale : 0.;
    return std::max(lengthFloor, representationTolerance());
  }

  /// Classify a parametric point against the closed loop (inside / outside / on-boundary).
  /// \a metric only sizes the on-boundary band (see boundaryBand); the winding count itself is a
  /// topological quantity and needs no metric at all.
  WireClassification classify(const Vec2& point, const ParametricMetric& metric = {}) const
  {
    const double band = boundaryBand(point, metric);
    const double bandSq = band * band;
    for (const auto& curve : curves) {
      if (curve.distanceSq(point) <= bandSq) {
        return WireClassification::Boundary;
      }
    }
    // Each curve's polyline already ends on the loop-canonical seam vertices, so the half-open
    // crossing convention stays consistent across seams without any substitution here.
    int crossings = 0;
    for (const auto& curve : curves) {
      crossings += curve.rightwardCrossings(point, curve.loopStart(), curve.loopEnd());
    }
    return (crossings % 2 == 1) ? WireClassification::Inside : WireClassification::Outside;
  }

  /// Ordered, closed boundary polyline; arcs are sampled into \a segmentsPerArc chords. This is a
  /// mesh-independent hook for visualization and tessellated fallback of curved boundaries.
  std::vector<Vec2> sampledBoundary(int segmentsPerArc = kArcSamples) const
  {
    std::vector<Vec2> samples;
    if (curves.empty()) {
      return samples;
    }
    for (const auto& curve : curves) {
      if (curve.kind == CurveKind::Line) {
        samples.push_back(curve.startPoint());
      } else if (curve.kind == CurveKind::BSpline) {
        // adaptively flatten and append every sample except the closing one (the next curve's
        // start reproduces it)
        std::vector<Vec2> curveSamples;
        curve.bsplineSampleInto(curveSamples);
        for (size_t index = 0; index + 1 < curveSamples.size(); ++index) {
          samples.push_back(curveSamples[index]);
        }
      } else {
        // Scale the chord count by the arc's sweep so a shared rim with a quadric wall (which
        // uses lround(kArcSamples * sweep / 2pi) segments) samples the identical vertices and
        // the closure half-edge check cancels. A full circle keeps segmentsPerArc chords.
        const int arcSteps =
          std::max(1, static_cast<int>(std::lround(segmentsPerArc * std::abs(curve.sweep()) / kTwoPi)));
        for (int step = 0; step < arcSteps; ++step) {
          samples.push_back(curve.pointAt(static_cast<double>(step) / arcSteps));
        }
      }
    }
    samples.push_back(samples.front());
    return samples;
  }
};

/// \name Curve-wire trim helpers for quadric parametric domains
/// A quadric (cylinder, cone, sphere) can carry a general line/arc CurveWire trim in its
/// periodic (u, v) parameter domain (u = phi, v = height or theta) in place of the scalar
/// parametric rectangle. These free helpers keep the outer/inner-wire classification, the
/// numerical capacity quadrature and the phi unwrapping in one place so all three quadric
/// surfaces reuse them. The trim wire is a non-wrapping loop within one phi window (<= 2pi).
/// @{

/// Shift \a angle by whole turns to lie as close as possible to the window [uMin, uMax], so a
/// raw atan2 result can be classified against a trim wire whose phi range may straddle the
/// atan2 branch cut.
inline double unwrapAngleInto(double angle, double uMin, double uMax)
{
  const double windowCenter = 0.5 * (uMin + uMax);
  return angle - kTwoPi * std::round((angle - windowCenter) / kTwoPi);
}

/// Classify a parametric point against a curve-wire trim: one outer loop plus optional inner
/// (hole) loops. Returns true when the point is inside or on the boundary of the outer loop and
/// not strictly inside any hole; \a boundary (when given) reports an on-boundary hit.
inline bool curveTrimContains(const CurveWire& outerWire, const std::vector<CurveWire>& innerWires,
                              const Vec2& point, bool* boundary = nullptr,
                              const ParametricMetric& metric = {})
{
  if (boundary != nullptr) {
    *boundary = false;
  }
  const auto outerClassification = outerWire.classify(point, metric);
  if (outerClassification == WireClassification::Outside) {
    return false;
  }
  if (outerClassification == WireClassification::Boundary) {
    if (boundary != nullptr) {
      *boundary = true;
    }
    return true;
  }
  for (const auto& innerWire : innerWires) {
    const auto innerClassification = innerWire.classify(point, metric);
    if (innerClassification == WireClassification::Boundary) {
      if (boundary != nullptr) {
        *boundary = true;
      }
      return true;
    }
    if (innerClassification == WireClassification::Inside) {
      return false;
    }
  }
  return true;
}

/// Number of Gauss-Legendre nodes per contour sub-interval in integrateOverCurveTrimByParts, and
/// the largest span in u one sub-interval may cover.
///
/// The antiderivatives the quadric surfaces supply are entire in u -- sines, cosines and u itself
/// -- so Gauss-Legendre converges geometrically rather than at any fixed order, and the only thing
/// that has to be bounded is how many oscillations sit inside one interval. Capping the span at a
/// quarter turn and spending 20 nodes on it leaves a wide margin: the error of n-node
/// Gauss-Legendre on exp(i w x) over a half-width h falls like (w h / 2)^(2n) / (2n)!, which for
/// w h = pi/8 and n = 20 is far below double precision.
inline constexpr int kContourQuadratureOrder = 20;
inline constexpr double kContourMaxSpanU = 0.25 * kPi;

/// Integrate F(u, v) dv along one directed curve of a trim wire, from \a from to \a to in the
/// curve's own [0, 1] parameter.
template <typename Antiderivative>
double contourIntegralAlongCurve(const Curve2D& curve, const Antiderivative& antiderivative, double from,
                                 double to)
{
  static thread_local std::vector<double> nodes;
  static thread_local std::vector<double> weights;
  if (static_cast<int>(nodes.size()) != kContourQuadratureOrder) {
    gaussLegendre(kContourQuadratureOrder, nodes, weights);
  }
  // split so no piece spans more than a quarter turn in u
  const double spanU = std::abs(curve.pointAt(to).uCoord - curve.pointAt(from).uCoord);
  const int pieces = std::max(1, static_cast<int>(std::ceil(spanU / kContourMaxSpanU)));
  double total = 0.;
  for (int piece = 0; piece < pieces; ++piece) {
    const double low = from + (to - from) * piece / pieces;
    const double high = from + (to - from) * (piece + 1) / pieces;
    const double half = 0.5 * (high - low);
    const double mid = 0.5 * (high + low);
    for (int nodeIndex = 0; nodeIndex < kContourQuadratureOrder; ++nodeIndex) {
      const double parameter = mid + half * nodes[nodeIndex];
      const Vec2 point = curve.pointAt(parameter);
      const Vec2 derivative = curve.derivativeAt(parameter);
      total += weights[nodeIndex] * half * antiderivative(point.uCoord, point.vCoord) * derivative.vCoord;
    }
  }
  return total;
}

/// Green's theorem for a wire-trimmed patch: given \a antiderivative F with dF/du = f,
///
///     double-integral over D of f(u, v) du dv  =  contour integral around dD of F(u, v) dv
///
/// where D is the outer loop minus the holes. The wires are already oriented for this (outer
/// counter-clockwise, inner clockwise, normalized by buildCurveTrim), so summing every loop's
/// every curve gives the signed total directly.
///
/// This replaces a 128x128 midpoint rule over a characteristic function. That rule books every
/// boundary cell whole or not at all, so it converged at O(1/N) and *oscillated* as the staircase
/// re-phased -- on Bagger/BucketLink2 it gave 16.004, 17.710, 16.927, 17.244, 17.032 cm^3 at 128
/// to 2048 samples per axis against OpenCascade's 17.079. The gate's 1e-6 relative capacity band
/// was unreachable at any practical N, and the entire capacity column of the Bagger gate was this
/// and no geometry defect at all (CodeReview_Fable.md, finding H2 in Section 13).
///
/// It is also much cheaper: 16384 point-in-wire classifications per patch, each walking every
/// curve of every loop, become a few hundred closed-form evaluations.
///
/// The contour is closed explicitly. A wire's curves are imported with independently sampled
/// endpoints that agree only to kWireJoinTolerance, and a first-order gap in the contour is a
/// first-order error in the integral -- 1e-5 parametric units against a 1e-6 target band. So each
/// seam is bridged by the straight segment between the raw endpoints, which is exactly the piece
/// that makes the loop a loop. (CurveWire::signedArea does not do this and carries the same
/// residual; it is not on this path.)
template <typename Antiderivative>
double integrateOverCurveTrimByParts(const CurveWire& outerWire, const std::vector<CurveWire>& innerWires,
                                     const Antiderivative& antiderivative)
{
  const auto loopIntegral = [&antiderivative](const CurveWire& wire) {
    double total = 0.;
    for (size_t index = 0; index < wire.curves.size(); ++index) {
      const auto& curve = wire.curves[index];
      total += contourIntegralAlongCurve(curve, antiderivative, 0., 1.);
      // seam bridge: a straight run from this curve's end to the next curve's start
      const Vec2 seamFrom = curve.endPoint();
      const Vec2 seamTo = wire.curves[(index + 1) % wire.curves.size()].startPoint();
      const double deltaV = seamTo.vCoord - seamFrom.vCoord;
      if (deltaV != 0.) {
        const Curve2D bridge = Curve2D::makeLine(seamFrom, seamTo);
        total += contourIntegralAlongCurve(bridge, antiderivative, 0., 1.);
      }
    }
    return total;
  };

  double total = loopIntegral(outerWire);
  for (const auto& innerWire : innerWires) {
    total += loopIntegral(innerWire);
  }
  return total;
}

/// Numerically integrate \a integrand(u, v) over the curve-wire-trimmed region (outer loop
/// minus holes) by midpoint quadrature on a regular grid across the outer loop's parametric
/// bounding box.
///
/// Superseded for capacity by integrateOverCurveTrimByParts, which is both exact-in-practice and
/// faster; kept because it is the only integrator that needs nothing of the integrand but its
/// value, and so is the independent check the unit tests measure the contour form against.
template <typename Integrand>
double integrateOverCurveTrim(const CurveWire& outerWire, const std::vector<CurveWire>& innerWires,
                              const Integrand& integrand, int samplesPerAxis = 128)
{
  Vec2 lower{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  Vec2 upper{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
  outerWire.parametricBounds(lower, upper);
  if (!finite(lower) || !finite(upper) || samplesPerAxis < 1) {
    return 0.;
  }
  const double stepU = (upper.uCoord - lower.uCoord) / samplesPerAxis;
  const double stepV = (upper.vCoord - lower.vCoord) / samplesPerAxis;
  const double cellArea = stepU * stepV;
  double sum = 0.;
  for (int indexU = 0; indexU < samplesPerAxis; ++indexU) {
    const double uCoord = lower.uCoord + (indexU + 0.5) * stepU;
    for (int indexV = 0; indexV < samplesPerAxis; ++indexV) {
      const double vCoord = lower.vCoord + (indexV + 0.5) * stepV;
      if (curveTrimContains(outerWire, innerWires, {uCoord, vCoord})) {
        sum += integrand(uCoord, vCoord) * cellArea;
      }
    }
  }
  return sum;
}

/// Build validated outer and inner trim wires from Curve2D loops and report the outer loop's
/// parametric bounding box. Wires are auto-reoriented (outer CCW, inner CW). Rejects a trim
/// that spans more than a full turn in phi (u), which a non-wrapping line/arc loop cannot
/// represent. Shared by the three quadric surfaces.
inline bool buildCurveTrim(const std::vector<Curve2D>& outerTrim,
                           const std::vector<std::vector<Curve2D>>& innerTrims, CurveWire& outerWire,
                           std::vector<CurveWire>& innerWires, Vec2& lower, Vec2& upper,
                           std::string& errorMessage, const ParametricMetric& metric = {})
{
  WireStatus status = WireStatus::Valid;
  if (!outerWire.initialize(outerTrim, WireRole::Outer, status, metric)) {
    errorMessage = std::string("quadric outer trim wire invalid: ") + wireStatusMessage(status);
    return false;
  }
  innerWires.clear();
  innerWires.reserve(innerTrims.size());
  for (const auto& innerLoop : innerTrims) {
    CurveWire innerWire;
    WireStatus innerStatus = WireStatus::Valid;
    if (!innerWire.initialize(innerLoop, WireRole::Inner, innerStatus, metric)) {
      errorMessage = std::string("quadric inner trim wire invalid: ") + wireStatusMessage(innerStatus);
      return false;
    }
    innerWires.push_back(std::move(innerWire));
  }
  lower = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  upper = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
  outerWire.parametricBounds(lower, upper);
  if (!finite(lower) || !finite(upper)) {
    errorMessage = "quadric trim wire has non-finite parametric bounds";
    return false;
  }
  if (upper.uCoord - lower.uCoord > kTwoPi + kTolerance) {
    // The box above is conservative -- a B-spline contributes its control-point hull, which can
    // overshoot the curve considerably. A closed intersection curve that wraps almost a full turn
    // has poles outside its own span, so the conservative box crosses 2*pi and a perfectly legal
    // through-hole host face was rejected outright. Re-measure on the curves themselves before
    // refusing (CodeReview_Fable.md K2).
    Vec2 tightLower{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Vec2 tightUpper{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    outerWire.tightParametricBounds(tightLower, tightUpper);
    if (!finite(tightLower) || !finite(tightUpper) || tightUpper.uCoord - tightLower.uCoord > kTwoPi + kTolerance) {
      errorMessage = "quadric trim wire spans more than a full turn in phi";
      return false;
    }
    // the wire is admissible; keep the tight box, since the conservative one is not a valid
    // parametric window for a periodic coordinate once it exceeds a full turn
    lower = tightLower;
    upper = tightUpper;
  }
  return true;
}

/// Sub-sample a curve-wire loop in its (u, v) parameter domain, splitting each edge so that its
/// span in u (phi) is chorded at \a segmentsPerTurn chords per full turn. A u = const edge (a
/// phi seam or a straight generator) stays a single segment, while a v = const edge (a rim
/// spanning phi) is chorded like an adjacent circular cap — so a straight parametric line, which
/// maps to a *curved* 3D edge, matches the neighbouring surfaces' rim sampling and shared rims
/// cancel in the solid-closure half-edge check. Returns the open ordered sample ring.
inline std::vector<Vec2> sampleCurveWireByU(const CurveWire& wire, int segmentsPerTurn = kArcSamples)
{
  std::vector<Vec2> samples;
  for (const auto& curve : wire.curves) {
    if (curve.kind == CurveKind::BSpline) {
      // adaptively flatten in the parameter domain; append every sample except the closing one
      std::vector<Vec2> curveSamples;
      curve.bsplineSampleInto(curveSamples);
      for (size_t index = 0; index + 1 < curveSamples.size(); ++index) {
        samples.push_back(curveSamples[index]);
      }
      continue;
    }
    Vec2 lower{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Vec2 upper{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    curve.extendBounds(lower, upper);
    const double uSpan = upper.uCoord - lower.uCoord;
    int steps = std::max(1, static_cast<int>(std::lround(segmentsPerTurn * uSpan / kTwoPi)));
    if (curve.kind == CurveKind::Arc) {
      steps = std::max(steps, static_cast<int>(std::lround(segmentsPerTurn * std::abs(curve.sweep()) / kTwoPi)));
      steps = std::max(steps, 1);
    }
    for (int step = 0; step < steps; ++step) {
      samples.push_back(curve.pointAt(static_cast<double>(step) / steps));
    }
  }
  return samples;
}

/// Append the display triangulation of a wire-trimmed quadric patch: sub-sample the outer trim
/// loop in (u, v), map each sample to 3D with \a mapUV, and ear-clip the sampled polygon (holes
/// are omitted from the display mesh, as for the planar surfaces — visualization never drives
/// navigation).
template <typename MapUV>
void appendCurveTrimMesh(const CurveWire& outerWire, const MapUV& mapUV, std::vector<Vec3>& vertices,
                         std::vector<std::array<int, 3>>& triangles)
{
  SurfaceWire sampledWire;
  sampledWire.vertices = sampleCurveWireByU(outerWire);
  if (sampledWire.vertices.size() < 3) {
    return;
  }
  const int firstVertexIndex = static_cast<int>(vertices.size());
  for (const auto& sample : sampledWire.vertices) {
    vertices.push_back(mapUV(sample.uCoord, sample.vCoord));
  }
  for (const auto& triangle : triangulateSimpleWire(sampledWire)) {
    triangles.push_back(
      {firstVertexIndex + triangle[0], firstVertexIndex + triangle[1], firstVertexIndex + triangle[2]});
  }
}

/// Append the 3D directed boundary edges of a wire-trimmed quadric patch for the solid-closure
/// half-edge check. \a orientationSign must be positive when the (u, v) -> 3D map \a mapUV is
/// orientation-consistent with the outward-of-solid normal (so a CCW parametric loop yields a
/// CCW 3D loop as seen from outside) and negative otherwise; each edge is reversed when it is
/// negative.
template <typename MapUV>
void appendCurveTrimEdges(const CurveWire& outerWire, const std::vector<CurveWire>& innerWires,
                          const MapUV& mapUV, double orientationSign,
                          std::vector<std::pair<Vec3, Vec3>>& edges)
{
  auto appendLoop = [&](const CurveWire& wire) {
    const auto samples = sampleCurveWireByU(wire);
    const size_t sampleCount = samples.size();
    for (size_t sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
      const Vec2& current = samples[sampleIndex];
      const Vec2& next = samples[(sampleIndex + 1) % sampleCount];
      const Vec3 edgeStart = mapUV(current.uCoord, current.vCoord);
      const Vec3 edgeEnd = mapUV(next.uCoord, next.vCoord);
      if (orientationSign >= 0.) {
        edges.emplace_back(edgeStart, edgeEnd);
      } else {
        edges.emplace_back(edgeEnd, edgeStart);
      }
    }
  };
  appendLoop(outerWire);
  for (const auto& innerWire : innerWires) {
    appendLoop(innerWire);
  }
}
/// @}

/// One trim loop of one face, sampled as an ordered 3D polyline.
///
/// This is the same boundary the half-edge closure check consumes as a bag of chords, kept as a
/// *curve* instead. The difference is the whole point: two faces meeting along a shared edge
/// sample that edge independently -- different chord counts, and phases that only coincide for
/// quarter-turn related frames -- so their vertices genuinely are not the same points and no
/// tolerance on vertex equality can match them. Their polylines, compared as curves, do match.
/// See scripts/geometry/TolerancePolicy.md section 1.3 (finding K9/S8).
struct SurfaceRim {
  int surfaceIndex = -1;    ///< index of the owning face in the solid's surface list
  bool closed = false;      ///< the polyline returns to its own first point
  std::vector<Vec3> points; ///< consecutive samples; a closed rim does not repeat the first point
};

/// Chain a face's directed boundary chords into rims (maximal polylines), appending them to
/// \a rims.
///
/// Chaining is by matching endpoints rather than by emission order, because the emission order is
/// not always loop-consecutive: the parametric-rectangle quadrics interleave their two rims (one
/// chord of the bottom rim, one of the top, ...) and only the wire-trimmed and planar faces emit
/// each loop in one run. Endpoint matching costs nothing extra and is right for both.
///
/// The weld tolerance is kTolerance. Within one face every shared vertex comes from the same
/// parametrisation evaluated with the same arguments, so it agrees to round-off; kTolerance is
/// many orders of magnitude above that and far below any real feature.
inline void assembleRims(const std::vector<std::pair<Vec3, Vec3>>& edges, std::vector<SurfaceRim>& rims)
{
  if (edges.empty()) {
    return;
  }
  auto quantize = [](double value) { return static_cast<int64_t>(std::llround(value / kTolerance)); };
  using VertexKey = std::tuple<int64_t, int64_t, int64_t>;
  auto keyOf = [&](const Vec3& point) {
    return VertexKey{quantize(point.xCoord), quantize(point.yCoord), quantize(point.zCoord)};
  };

  // A patch that closes on itself emits its seam twice, once in each direction: a cylinder or
  // cone whose sweep is a full turn that fullSweep() does not recognise as one (the end phi and
  // the start phi then differ only by round-off), and the same for a sphere or torus. Such a pair
  // bounds nothing -- it cancels in the half-edge check for exactly this reason -- and chaining it
  // would invent a two-point rim straddling the patch, with a gap the size of the patch. Cancel
  // reversed duplicates first. Indexing by the midpoint, which a chord shares with its reverse,
  // keeps the lookup a single 27-cell probe.
  std::vector<bool> consumed(edges.size(), false);
  std::map<VertexKey, std::vector<size_t>> edgesByMidpoint;
  for (size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
    const Vec3 midpoint = (edges[edgeIndex].first + edges[edgeIndex].second) * 0.5;
    const auto [xKey, yKey, zKey] = keyOf(midpoint);
    bool cancelled = false;
    for (int64_t dx = -1; dx <= 1 && !cancelled; ++dx) {
      for (int64_t dy = -1; dy <= 1 && !cancelled; ++dy) {
        for (int64_t dz = -1; dz <= 1 && !cancelled; ++dz) {
          const auto found = edgesByMidpoint.find(VertexKey{xKey + dx, yKey + dy, zKey + dz});
          if (found == edgesByMidpoint.end()) {
            continue;
          }
          for (const size_t candidate : found->second) {
            if (consumed[candidate] ||
                distanceSq(edges[candidate].first, edges[edgeIndex].second) > kToleranceSq ||
                distanceSq(edges[candidate].second, edges[edgeIndex].first) > kToleranceSq) {
              continue;
            }
            consumed[candidate] = true;
            consumed[edgeIndex] = true;
            cancelled = true;
            break;
          }
        }
      }
    }
    if (!cancelled) {
      edgesByMidpoint[keyOf(midpoint)].push_back(edgeIndex);
    }
  }

  std::map<VertexKey, std::vector<size_t>> edgesByStart;
  for (size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
    if (!consumed[edgeIndex]) {
      edgesByStart[keyOf(edges[edgeIndex].first)].push_back(edgeIndex);
    }
  }

  // A vertex can land either side of a lattice boundary, so probe the 27 neighbouring cells and
  // accept the first unused chord whose start really is within kTolerance.
  auto findSuccessor = [&](const Vec3& point) -> long long {
    const auto [xKey, yKey, zKey] = keyOf(point);
    for (int64_t dx = -1; dx <= 1; ++dx) {
      for (int64_t dy = -1; dy <= 1; ++dy) {
        for (int64_t dz = -1; dz <= 1; ++dz) {
          const auto found = edgesByStart.find(VertexKey{xKey + dx, yKey + dy, zKey + dz});
          if (found == edgesByStart.end()) {
            continue;
          }
          for (const size_t candidate : found->second) {
            if (!consumed[candidate] && distanceSq(edges[candidate].first, point) <= kToleranceSq) {
              return static_cast<long long>(candidate);
            }
          }
        }
      }
    }
    return -1;
  };

  for (size_t seed = 0; seed < edges.size(); ++seed) {
    if (consumed[seed]) {
      continue;
    }
    consumed[seed] = true;
    SurfaceRim rim;
    rim.points.push_back(edges[seed].first);
    rim.points.push_back(edges[seed].second);
    while (true) {
      if (distanceSq(rim.points.back(), rim.points.front()) <= kToleranceSq) {
        rim.closed = true;
        rim.points.pop_back(); // a closed rim does not repeat its first point
        break;
      }
      const long long next = findSuccessor(rim.points.back());
      if (next < 0) {
        break; // an open chain: the face's boundary is not a set of closed loops
      }
      consumed[static_cast<size_t>(next)] = true;
      rim.points.push_back(edges[static_cast<size_t>(next)].second);
    }
    if (rim.points.size() >= 2) {
      rims.push_back(std::move(rim));
    }
  }
}

/// Abstract analytic surface patch: one support surface plus its trim wires.
///
/// Concrete surfaces (planar, and later cylindrical/spherical/conical) implement the kernels
/// the O2BVHSurfaceSolid navigation needs. The interface keeps navigation independent of the
/// specific analytic surface type.
class BoundedSurface
{
 public:
  virtual ~BoundedSurface() = default;

  /// Accumulate a conservative axis-aligned bounding box of the trimmed patch.
  virtual void conservativeBounds(Vec3& lower, Vec3& upper) const = 0;

  /// True if the 3D point lies on the trimmed patch within tolerance.
  virtual bool containsPointOnSurface(const Vec3& point) const = 0;

  /// Append every intersection of the ray with the trimmed patch whose ray parameter lies in
  /// [minDistance, maxDistance], each carrying the outward-oriented surface normal at the hit.
  /// A plane contributes at most one hit, but quadric patches (cylinder, sphere, cone) can be
  /// crossed twice by the same ray — parity-based containment and entering/exiting distance
  /// queries need all of them. Tangential grazes (double roots) must not be reported, so that
  /// crossing parity stays consistent.
  virtual void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                                   double maxDistance, std::vector<RayHit>& hits) const = 0;

  /// Nearest-hit convenience wrapper around appendIntersections. Callers that only need the
  /// closest intersection (or want to prune with a shrinking maxDistance) can stop here.
  bool intersectRay(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                    double maxDistance, double& distance, Vec3& hitNormal) const
  {
    std::vector<RayHit> hits;
    appendIntersections(rayOrigin, rayDirection, minDistance, maxDistance, hits);
    if (hits.empty()) {
      return false;
    }
    const auto nearestHit = std::min_element(hits.begin(), hits.end(),
                                             [](const RayHit& firstHit, const RayHit& secondHit) {
                                               return firstHit.distance < secondHit.distance;
                                             });
    distance = nearestHit->distance;
    hitNormal = nearestHit->normal;
    return true;
  }

  /// Squared distance from a 3D point to the trimmed patch (used for Safety).
  virtual double distanceSqToPatch(const Vec3& point) const = 0;

  /// Outward-oriented normal at (or nearest to) the given point.
  virtual Vec3 normalAt(const Vec3& point) const = 0;

  /// The first fundamental form (\a gUU, \a gUV, \a gVV) of the surface at the parametric point
  /// \a uv, so that a parametric displacement (du, dv) from there spans the 3D length
  ///   sqrt(gUU*du*du + 2*gUV*du*dv + gVV*dv*dv)   (see parametricLengthSq).
  ///
  /// This is what turns the parametric coordinates -- which mix radians and centimetres, with a
  /// different mixture per surface family -- back into distances. Every tolerance applied to a
  /// parametric separation has to go through it, or it means a different physical thing on a
  /// 0.01 cm hole than on a 100 cm cylinder.
  ///
  /// Two properties callers must respect. The form is *not* constant over the domain (only the
  /// planar families have that luxury), so it must be evaluated at the point of interest rather
  /// than once per surface. And gUU vanishes at a sphere's poles and at a cone's apex, where a
  /// separation in the angular coordinate really is of zero length; that is the correct answer,
  /// but code that divides by the scale has to handle it.
  virtual void parametricMetric(const Vec2& uv, double& gUU, double& gUV, double& gVV) const = 0;

  /// The 3D length squared spanned by a parametric displacement \a delta starting at \a uv.
  double parametricLengthSqAt(const Vec2& uv, const Vec2& delta) const
  {
    double gUU = 0.;
    double gUV = 0.;
    double gVV = 0.;
    parametricMetric(uv, gUU, gUV, gVV);
    return parametricLengthSq(gUU, gUV, gVV, delta);
  }

  /// Signed divergence-theorem contribution to the enclosed volume.
  virtual double capacityContribution() const = 0;

  /// Whether capacityContribution() is analytically exact for this surface.
  virtual bool capacityIsExact() const = 0;

  /// Append this patch's visualization triangulation (navigation must never depend on it).
  virtual void appendDisplayMesh(std::vector<Vec3>& vertices,
                                 std::vector<std::array<int, 3>>& triangles) const = 0;

  /// Append the 3D directed boundary edges of the patch, for solid-closure validation.
  virtual void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const = 0;

  /// Append the patch's trim boundary as rims -- one ordered polyline per boundary loop -- for
  /// the rim-based closure measurement.
  ///
  /// This is the same data appendDirectedEdges() emits, kept as curves instead of as a bag of
  /// chords, so the default implementation just chains that output (see assembleRims). Override
  /// only for a class whose rims cannot be recovered by endpoint chaining; none currently is.
  virtual void appendRims(std::vector<SurfaceRim>& rims) const
  {
    std::vector<std::pair<Vec3, Vec3>> edges;
    appendDirectedEdges(edges);
    assembleRims(edges, rims);
  }
};

/// A bounded planar surface: an infinite plane frame trimmed by one outer wire and optional
/// inner (hole) wires expressed in the plane's local 2D coordinates.
class PlanarBoundedSurface final : public BoundedSurface
{
 public:
  bool initialize(const Vec3& surfaceOrigin, const Vec3& surfaceAxisU, const Vec3& surfaceAxisV,
                  const std::vector<Vec2>& outerWireVertices,
                  const std::vector<std::vector<Vec2>>& innerWireVertices, std::string& errorMessage)
  {
    if (!finite(surfaceOrigin) || !finite(surfaceAxisU) || !finite(surfaceAxisV)) {
      errorMessage = "surface frame contains a non-finite value";
      return false;
    }

    mOrigin = surfaceOrigin;
    mAxisU = surfaceAxisU;
    mAxisV = surfaceAxisV;
    const Vec3 normalVector = cross(mAxisU, mAxisV);
    mAreaScale = norm(normalVector);
    if (mAreaScale <= kTolerance) {
      errorMessage = "surface frame axes are degenerate";
      return false;
    }
    mNormal = normalVector * (1. / mAreaScale);

    mMetricUU = dot(mAxisU, mAxisU);
    mMetricUV = dot(mAxisU, mAxisV);
    mMetricVV = dot(mAxisV, mAxisV);
    const double metricDet = mMetricUU * mMetricVV - mMetricUV * mMetricUV;
    if (std::abs(metricDet) <= kToleranceSq) {
      errorMessage = "surface frame metric is singular";
      return false;
    }
    mInverseMetricDet = 1. / metricDet;

    WireStatus outerStatus = WireStatus::Valid;
    const ParametricMetric metric = parametricMetricOf(*this);
    if (!mOuterWire.initialize(outerWireVertices, WireRole::Outer, outerStatus, metric)) {
      errorMessage = std::string("outer wire invalid: ") + wireStatusMessage(outerStatus);
      return false;
    }
    mOuterReoriented = (outerStatus == WireStatus::Reversed);

    mInnerWires.clear();
    mInnerWires.reserve(innerWireVertices.size());
    mInnerReoriented = false;
    for (const auto& innerWireInput : innerWireVertices) {
      SurfaceWire innerWire;
      WireStatus innerStatus = WireStatus::Valid;
      if (!innerWire.initialize(innerWireInput, WireRole::Inner, innerStatus, metric)) {
        errorMessage = std::string("inner wire invalid: ") + wireStatusMessage(innerStatus);
        return false;
      }
      mInnerReoriented = mInnerReoriented || (innerStatus == WireStatus::Reversed);
      mInnerWires.emplace_back(std::move(innerWire));
    }
    return true;
  }

  /// True if either the outer or any inner wire had to be re-oriented during initialization.
  bool wasReoriented() const { return mOuterReoriented || mInnerReoriented; }

  const Vec3& normal() const { return mNormal; }
  const SurfaceWire& outerWire() const { return mOuterWire; }
  const std::vector<SurfaceWire>& innerWires() const { return mInnerWires; }

  Vec3 toGlobal(const Vec2& point) const
  {
    return mOrigin + mAxisU * point.uCoord + mAxisV * point.vCoord;
  }

  Vec2 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - mOrigin;
    const double projectionU = dot(relativePoint, mAxisU);
    const double projectionV = dot(relativePoint, mAxisV);
    return {(projectionU * mMetricVV - projectionV * mMetricUV) * mInverseMetricDet,
            (projectionV * mMetricUU - projectionU * mMetricUV) * mInverseMetricDet};
  }

  double planeDistance(const Vec3& point) const { return dot(point - mOrigin, mNormal); }

  bool containsLocal(const Vec2& point, bool* boundary = nullptr) const
  {
    if (boundary != nullptr) {
      *boundary = false;
    }

    const ParametricMetric metric = parametricMetricOf(*this);
    const auto outerClassification = mOuterWire.classify(point, metric);
    if (outerClassification == WireClassification::Outside) {
      return false;
    }
    if (outerClassification == WireClassification::Boundary) {
      if (boundary != nullptr) {
        *boundary = true;
      }
      return true;
    }

    for (const auto& innerWire : mInnerWires) {
      const auto innerClassification = innerWire.classify(point, metric);
      if (innerClassification == WireClassification::Boundary) {
        if (boundary != nullptr) {
          *boundary = true;
        }
        return true;
      }
      if (innerClassification == WireClassification::Inside) {
        return false;
      }
    }
    return true;
  }

  bool containsPointOnSurface(const Vec3& point) const override
  {
    if (std::abs(planeDistance(point)) > kTolerance) {
      return false;
    }
    return containsLocal(toLocal(point));
  }

  void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                           double maxDistance, std::vector<RayHit>& hits) const override
  {
    const double denominator = dot(mNormal, rayDirection);
    if (std::abs(denominator) <= kTolerance) {
      return;
    }
    const double candidateDistance = dot(mOrigin - rayOrigin, mNormal) / denominator;
    if (candidateDistance < minDistance || candidateDistance > maxDistance) {
      return;
    }
    const Vec3 candidatePoint = rayOrigin + rayDirection * candidateDistance;
    bool onTrimBoundary = false;
    if (!containsLocal(toLocal(candidatePoint), &onTrimBoundary)) {
      return;
    }
    hits.push_back({candidateDistance, mNormal, onTrimBoundary});
  }

  double distanceSqToEdges(const Vec3& point, const SurfaceWire& wire) const
  {
    double bestDistanceSq = std::numeric_limits<double>::infinity();
    for (size_t vertexIndex = 0; vertexIndex < wire.vertices.size(); ++vertexIndex) {
      const Vec3 segmentStart = toGlobal(wire.vertices[vertexIndex]);
      const Vec3 segmentEnd = toGlobal(wire.vertices[(vertexIndex + 1) % wire.vertices.size()]);
      bestDistanceSq = std::min(bestDistanceSq, pointSegmentDistanceSq(point, segmentStart, segmentEnd));
    }
    return bestDistanceSq;
  }

  double distanceSqToPatch(const Vec3& point) const override
  {
    const Vec2 projectedPoint = toLocal(point);
    if (containsLocal(projectedPoint)) {
      const double signedPlaneDistance = planeDistance(point);
      return signedPlaneDistance * signedPlaneDistance;
    }

    double bestDistanceSq = distanceSqToEdges(point, mOuterWire);
    for (const auto& innerWire : mInnerWires) {
      bestDistanceSq = std::min(bestDistanceSq, distanceSqToEdges(point, innerWire));
    }
    return bestDistanceSq;
  }

  Vec3 normalAt(const Vec3&) const override { return mNormal; }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    auto extendPoint = [&](const Vec2& surfacePoint) {
      const Vec3 globalPoint = toGlobal(surfacePoint);
      lower.xCoord = std::min(lower.xCoord, globalPoint.xCoord);
      lower.yCoord = std::min(lower.yCoord, globalPoint.yCoord);
      lower.zCoord = std::min(lower.zCoord, globalPoint.zCoord);
      upper.xCoord = std::max(upper.xCoord, globalPoint.xCoord);
      upper.yCoord = std::max(upper.yCoord, globalPoint.yCoord);
      upper.zCoord = std::max(upper.zCoord, globalPoint.zCoord);
    };

    for (const auto& vertex : mOuterWire.vertices) {
      extendPoint(vertex);
    }
    for (const auto& innerWire : mInnerWires) {
      for (const auto& vertex : innerWire.vertices) {
        extendPoint(vertex);
      }
    }
  }

  double area() const
  {
    double parametricArea = std::abs(mOuterWire.signedArea());
    for (const auto& innerWire : mInnerWires) {
      parametricArea -= std::abs(innerWire.signedArea());
    }
    return std::max(0., parametricArea) * mAreaScale;
  }

  /// Constant over the plane, and the only family with a cross term: the frame axes carry the
  /// domain's units and need be neither unit-length nor orthogonal. These are the same three
  /// numbers initialize() already computes to invert the frame.
  void parametricMetric(const Vec2&, double& gUU, double& gUV, double& gVV) const override
  {
    planeParametricMetric(mAxisU, mAxisV, gUU, gUV, gVV);
  }

  double capacityContribution() const override { return dot(mOrigin, mNormal) * area() / 3.; }

  bool capacityIsExact() const override { return true; }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (const auto& vertex : mOuterWire.vertices) {
      vertices.push_back(toGlobal(vertex));
    }

    const auto localTriangles = triangulateSimpleWire(mOuterWire);
    for (const auto& triangle : localTriangles) {
      triangles.push_back(
        {firstVertexIndex + triangle[0], firstVertexIndex + triangle[1], firstVertexIndex + triangle[2]});
    }
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    auto appendWire = [&](const SurfaceWire& wire) {
      for (size_t vertexIndex = 0; vertexIndex < wire.vertices.size(); ++vertexIndex) {
        const Vec3 edgeStart = toGlobal(wire.vertices[vertexIndex]);
        const Vec3 edgeEnd = toGlobal(wire.vertices[(vertexIndex + 1) % wire.vertices.size()]);
        edges.emplace_back(edgeStart, edgeEnd);
      }
    };
    appendWire(mOuterWire);
    for (const auto& innerWire : mInnerWires) {
      appendWire(innerWire);
    }
  }

 private:
  Vec3 mOrigin;
  Vec3 mAxisU;
  Vec3 mAxisV;
  Vec3 mNormal;
  double mMetricUU = 0.;
  double mMetricUV = 0.;
  double mMetricVV = 0.;
  double mInverseMetricDet = 0.;
  double mAreaScale = 0.;
  bool mOuterReoriented = false;
  bool mInnerReoriented = false;
  SurfaceWire mOuterWire;
  std::vector<SurfaceWire> mInnerWires;
};

/// A bounded planar surface trimmed by curved (line/arc) boundary loops: an infinite plane
/// trimmed by one outer CurveWire and optional inner (hole) CurveWires in the plane's local 2D
/// coordinates. This is the exact-cap counterpart to the polygonal PlanarBoundedSurface — it
/// represents disks and annuli (cylinder/cone end caps) without polygonal approximation.
///
/// The frame axes are required to be orthonormal (unit length, perpendicular). This keeps every
/// kernel exact and simple: parametric distances equal in-plane 3D distances, circles in the
/// parameter domain are circles in space, and the area scale is one.
class CurvedPlanarBoundedSurface final : public BoundedSurface
{
 public:
  bool initialize(const Vec3& surfaceOrigin, const Vec3& surfaceAxisU, const Vec3& surfaceAxisV,
                  const std::vector<Curve2D>& outerCurves,
                  const std::vector<std::vector<Curve2D>>& innerCurves, std::string& errorMessage)
  {
    if (!finite(surfaceOrigin) || !finite(surfaceAxisU) || !finite(surfaceAxisV)) {
      errorMessage = "surface frame contains a non-finite value";
      return false;
    }
    if (std::abs(norm(surfaceAxisU) - 1.) > kTolerance || std::abs(norm(surfaceAxisV) - 1.) > kTolerance ||
        std::abs(dot(surfaceAxisU, surfaceAxisV)) > kTolerance) {
      errorMessage = "curved planar surface requires orthonormal frame axes";
      return false;
    }

    mOrigin = surfaceOrigin;
    mAxisU = surfaceAxisU;
    mAxisV = surfaceAxisV;
    mNormal = cross(mAxisU, mAxisV);

    WireStatus outerStatus = WireStatus::Valid;
    const ParametricMetric metric = parametricMetricOf(*this);
    if (!mOuterWire.initialize(outerCurves, WireRole::Outer, outerStatus, metric)) {
      errorMessage = std::string("outer wire invalid: ") + wireStatusMessage(outerStatus);
      return false;
    }
    mReoriented = (outerStatus == WireStatus::Reversed);

    mInnerWires.clear();
    mInnerWires.reserve(innerCurves.size());
    for (const auto& innerCurveLoop : innerCurves) {
      CurveWire innerWire;
      WireStatus innerStatus = WireStatus::Valid;
      if (!innerWire.initialize(innerCurveLoop, WireRole::Inner, innerStatus, metric)) {
        errorMessage = std::string("inner wire invalid: ") + wireStatusMessage(innerStatus);
        return false;
      }
      mReoriented = mReoriented || (innerStatus == WireStatus::Reversed);
      mInnerWires.emplace_back(std::move(innerWire));
    }

    // A B-spline boundary makes the area (hence the capacity contribution) a numeric quadrature,
    // so flag the capacity as inexact (matching the wire-trimmed-quadric policy).
    mCapacityExact = !mOuterWire.hasBSpline();
    for (const auto& innerWire : mInnerWires) {
      mCapacityExact = mCapacityExact && !innerWire.hasBSpline();
    }
    return true;
  }

  /// True if the outer or any inner wire had to be re-oriented during initialization.
  bool wasReoriented() const { return mReoriented; }

  const Vec3& normal() const { return mNormal; }

  Vec3 toGlobal(const Vec2& point) const { return mOrigin + mAxisU * point.uCoord + mAxisV * point.vCoord; }

  Vec2 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - mOrigin;
    return {dot(relativePoint, mAxisU), dot(relativePoint, mAxisV)};
  }

  double planeDistance(const Vec3& point) const { return dot(point - mOrigin, mNormal); }

  bool containsLocal(const Vec2& point, bool* boundary = nullptr) const
  {
    if (boundary != nullptr) {
      *boundary = false;
    }

    const ParametricMetric metric = parametricMetricOf(*this);
    const auto outerClassification = mOuterWire.classify(point, metric);
    if (outerClassification == WireClassification::Outside) {
      return false;
    }
    if (outerClassification == WireClassification::Boundary) {
      if (boundary != nullptr) {
        *boundary = true;
      }
      return true;
    }

    for (const auto& innerWire : mInnerWires) {
      const auto innerClassification = innerWire.classify(point, metric);
      if (innerClassification == WireClassification::Boundary) {
        if (boundary != nullptr) {
          *boundary = true;
        }
        return true;
      }
      if (innerClassification == WireClassification::Inside) {
        return false;
      }
    }
    return true;
  }

  bool containsPointOnSurface(const Vec3& point) const override
  {
    if (std::abs(planeDistance(point)) > kTolerance) {
      return false;
    }
    return containsLocal(toLocal(point));
  }

  void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                           double maxDistance, std::vector<RayHit>& hits) const override
  {
    const double denominator = dot(mNormal, rayDirection);
    if (std::abs(denominator) <= kTolerance) {
      return;
    }
    const double candidateDistance = dot(mOrigin - rayOrigin, mNormal) / denominator;
    if (candidateDistance < minDistance || candidateDistance > maxDistance) {
      return;
    }
    bool onTrimBoundary = false;
    if (!containsLocal(toLocal(rayOrigin + rayDirection * candidateDistance), &onTrimBoundary)) {
      return;
    }
    hits.push_back({candidateDistance, mNormal, onTrimBoundary});
  }

  double distanceSqToPatch(const Vec3& point) const override
  {
    const Vec2 projectedPoint = toLocal(point);
    const double signedPlaneDistance = planeDistance(point);
    if (containsLocal(projectedPoint)) {
      return signedPlaneDistance * signedPlaneDistance;
    }

    // exact for an orthonormal frame: split into in-plane distance to the trim curves plus the
    // out-of-plane plane distance
    double bestCurveDistanceSq = std::numeric_limits<double>::infinity();
    for (const auto& curve : mOuterWire.curves) {
      bestCurveDistanceSq = std::min(bestCurveDistanceSq, curve.distanceSq(projectedPoint));
    }
    for (const auto& innerWire : mInnerWires) {
      for (const auto& curve : innerWire.curves) {
        bestCurveDistanceSq = std::min(bestCurveDistanceSq, curve.distanceSq(projectedPoint));
      }
    }
    return bestCurveDistanceSq + signedPlaneDistance * signedPlaneDistance;
  }

  Vec3 normalAt(const Vec3&) const override { return mNormal; }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    Vec2 parametricLower{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Vec2 parametricUpper{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    mOuterWire.parametricBounds(parametricLower, parametricUpper);

    // the affine image of the parametric AABB contains the patch; its corners bound the 3D AABB
    for (const double cornerU : {parametricLower.uCoord, parametricUpper.uCoord}) {
      for (const double cornerV : {parametricLower.vCoord, parametricUpper.vCoord}) {
        const Vec3 globalCorner = toGlobal({cornerU, cornerV});
        lower.xCoord = std::min(lower.xCoord, globalCorner.xCoord);
        lower.yCoord = std::min(lower.yCoord, globalCorner.yCoord);
        lower.zCoord = std::min(lower.zCoord, globalCorner.zCoord);
        upper.xCoord = std::max(upper.xCoord, globalCorner.xCoord);
        upper.yCoord = std::max(upper.yCoord, globalCorner.yCoord);
        upper.zCoord = std::max(upper.zCoord, globalCorner.zCoord);
      }
    }
  }

  double area() const
  {
    double parametricArea = std::abs(mOuterWire.signedArea());
    for (const auto& innerWire : mInnerWires) {
      parametricArea -= std::abs(innerWire.signedArea());
    }
    return std::max(0., parametricArea);
  }

  /// The identity form: initialize() rejects a frame whose axes are not orthonormal, so (u, v)
  /// here are already lengths in centimetres. (The polygon-wire plane is the general case.)
  void parametricMetric(const Vec2&, double& gUU, double& gUV, double& gVV) const override
  {
    gUU = 1.;
    gUV = 0.;
    gVV = 1.;
  }

  double capacityContribution() const override { return dot(mOrigin, mNormal) * area() / 3.; }

  bool capacityIsExact() const override { return mCapacityExact; }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    // triangulate the sampled outer boundary; holes are ignored in the display mesh (as for the
    // polygonal planar surface, visualization never influences navigation)
    auto samples = mOuterWire.sampledBoundary();
    if (samples.size() < 4) {
      return;
    }
    samples.pop_back(); // drop the closing duplicate

    SurfaceWire sampledWire;
    sampledWire.vertices = std::move(samples);

    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (const auto& vertex : sampledWire.vertices) {
      vertices.push_back(toGlobal(vertex));
    }
    for (const auto& triangle : triangulateSimpleWire(sampledWire)) {
      triangles.push_back(
        {firstVertexIndex + triangle[0], firstVertexIndex + triangle[1], firstVertexIndex + triangle[2]});
    }
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    auto appendWire = [&](const CurveWire& wire) {
      const auto samples = wire.sampledBoundary();
      for (size_t sampleIndex = 0; sampleIndex + 1 < samples.size(); ++sampleIndex) {
        edges.emplace_back(toGlobal(samples[sampleIndex]), toGlobal(samples[sampleIndex + 1]));
      }
    };
    appendWire(mOuterWire);
    for (const auto& innerWire : mInnerWires) {
      appendWire(innerWire);
    }
  }

 private:
  Vec3 mOrigin;
  Vec3 mAxisU;
  Vec3 mAxisV;
  Vec3 mNormal;
  bool mReoriented = false;
  bool mCapacityExact = true;
  CurveWire mOuterWire;
  std::vector<CurveWire> mInnerWires;
};

/// A bounded cylindrical surface: an infinite cylinder of given radius around an axis, trimmed
/// to a parametric rectangle (phi angular sweep x height range).
///
/// The parametric-rectangle trim matches TGeoTube/TGeoTubeSeg and the vast majority of
/// mechanical CAD cylinder faces; general trim wires on the periodic (phi, h) domain are left
/// for a later milestone. The frame is (axisU, axisV, axisW) orthonormal with
/// axisU x axisV = axisW (the cylinder axis); phi is measured from axisU towards axisV and the
/// height h along axisW, both relative to the reference point centerPoint (a point on the axis
/// at h = 0). An innerWall surface (the wall of a hole, e.g. the inner tube of a hollow
/// cylinder) has its outward-of-solid normal pointing towards the axis.
class CylindricalBoundedSurface final : public BoundedSurface
{
 public:
  bool initialize(const Vec3& centerPoint, const Vec3& axis, const Vec3& referenceAxisU, double radius,
                  double heightMin, double heightMax, double phiStart, double phiSweep, bool innerWall,
                  std::string& errorMessage)
  {
    if (!finite(centerPoint) || !finite(axis) || !finite(referenceAxisU) || !std::isfinite(radius) ||
        !std::isfinite(heightMin) || !std::isfinite(heightMax) || !std::isfinite(phiStart) ||
        !std::isfinite(phiSweep)) {
      errorMessage = "cylindrical surface parameter is non-finite";
      return false;
    }
    if (radius <= kTolerance) {
      errorMessage = "cylindrical surface needs a positive radius";
      return false;
    }
    if (heightMax - heightMin <= kTolerance) {
      errorMessage = "cylindrical surface needs a positive height range";
      return false;
    }
    if (phiSweep <= kTolerance || phiSweep > kTwoPi + kTolerance) {
      errorMessage = "cylindrical surface needs an angular sweep in (0, 2pi]";
      return false;
    }
    if (!makeFrame(axis, referenceAxisU, mAxisU, mAxisV, mAxisW, errorMessage)) {
      return false;
    }

    mCenter = centerPoint;
    mRadius = radius;
    mHeightMin = heightMin;
    mHeightMax = heightMax;
    mPhiStart = phiStart;
    mPhiSweep = std::min(phiSweep, kTwoPi);
    mNormalSign = innerWall ? -1. : 1.;
    return true;
  }

  /// Wire-trimmed overload: the scalar arguments define the cylinder surface (frame, radius) and
  /// a nominal window; the outer/inner Curve2D loops in the periodic (phi[rad], h[cm]) domain are
  /// authoritative for containment. The conservative window (used for the AABB, display mesh and
  /// the Safety lower bound) is tightened to the outer loop's parametric bounds.
  bool initialize(const Vec3& centerPoint, const Vec3& axis, const Vec3& referenceAxisU, double radius,
                  double heightMin, double heightMax, double phiStart, double phiSweep, bool innerWall,
                  const std::vector<Curve2D>& outerTrim, const std::vector<std::vector<Curve2D>>& innerTrims,
                  std::string& errorMessage)
  {
    if (!initialize(centerPoint, axis, referenceAxisU, radius, heightMin, heightMax, phiStart, phiSweep, innerWall,
                    errorMessage)) {
      return false;
    }
    Vec2 lower, upper;
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage,
                        parametricMetricOf(*this))) {
      return false;
    }
    mPhiStart = lower.uCoord;
    mPhiSweep = std::min(kTwoPi, upper.uCoord - lower.uCoord);
    mHeightMin = lower.vCoord;
    mHeightMax = upper.vCoord;
    mHasWireTrim = true;
    return true;
  }

  bool hasWireTrim() const { return mHasWireTrim; }

  /// True if the (phi, h) point lies in the trim wire (phi unwrapped into the wire window).
  bool pointInTrim(double phi, double height, bool* boundary = nullptr) const
  {
    const double uCoord = unwrapAngleInto(phi, mPhiStart, mPhiStart + mPhiSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, height}, boundary, parametricMetricOf(*this));
  }

  /// Build an orthonormal frame (U, V, W) with W along \a axis and U the projection of
  /// \a referenceAxisU perpendicular to W. Shared by all axis-symmetric quadric surfaces.
  static bool makeFrame(const Vec3& axis, const Vec3& referenceAxisU, Vec3& axisU, Vec3& axisV, Vec3& axisW,
                        std::string& errorMessage)
  {
    if (norm(axis) <= kTolerance) {
      errorMessage = "surface axis is degenerate";
      return false;
    }
    axisW = normalized(axis);
    const Vec3 projectedU = referenceAxisU - axisW * dot(referenceAxisU, axisW);
    if (norm(projectedU) <= kTolerance) {
      errorMessage = "surface reference axis is parallel to the main axis";
      return false;
    }
    axisU = normalized(projectedU);
    axisV = cross(axisW, axisU); // gives axisU x axisV = axisW
    return true;
  }

  bool fullSweep() const { return mPhiSweep >= kTwoPi - kTolerance; }

  Vec3 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - mCenter;
    return {dot(relativePoint, mAxisU), dot(relativePoint, mAxisV), dot(relativePoint, mAxisW)};
  }

  bool heightInRange(double height) const
  {
    return height >= mHeightMin - kTolerance && height <= mHeightMax + kTolerance;
  }

  bool phiInSweep(double phi) const
  {
    return angleInSweepRange(phi, mPhiStart, mPhiSweep, angularTolerance(mRadius));
  }

  Vec3 pointAt(double phi, double height) const
  {
    return mCenter + mAxisW * height + (mAxisU * std::cos(phi) + mAxisV * std::sin(phi)) * mRadius;
  }

  bool containsPointOnSurface(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    if (std::abs(radialDistance - mRadius) > kTolerance) {
      return false;
    }
    if (radialDistance <= kTolerance) {
      return !mHasWireTrim && heightInRange(localPoint.zCoord); // phi is undefined on the axis
    }
    const double phi = std::atan2(localPoint.yCoord, localPoint.xCoord);
    if (mHasWireTrim) {
      return pointInTrim(phi, localPoint.zCoord);
    }
    return heightInRange(localPoint.zCoord) && phiInSweep(phi);
  }

  void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                           double maxDistance, std::vector<RayHit>& hits) const override
  {
    const Vec3 localOrigin = toLocal(rayOrigin);
    const Vec3 localDirection{dot(rayDirection, mAxisU), dot(rayDirection, mAxisV), dot(rayDirection, mAxisW)};

    const double quadraticA = localDirection.xCoord * localDirection.xCoord +
                              localDirection.yCoord * localDirection.yCoord;
    if (quadraticA <= kToleranceSq) {
      return; // ray parallel to the axis: no transversal crossing of the lateral surface
    }
    const double quadraticB = 2. * (localOrigin.xCoord * localDirection.xCoord +
                                    localOrigin.yCoord * localDirection.yCoord);
    const double quadraticC = localOrigin.xCoord * localOrigin.xCoord +
                              localOrigin.yCoord * localOrigin.yCoord - mRadius * mRadius;
    const double discriminant = quadraticB * quadraticB - 4. * quadraticA * quadraticC;
    if (discriminant <= 0.) {
      return;
    }
    const double sqrtDiscriminant = std::sqrt(discriminant);
    const double firstRoot = (-quadraticB - sqrtDiscriminant) / (2. * quadraticA);
    const double secondRoot = (-quadraticB + sqrtDiscriminant) / (2. * quadraticA);
    if (sameIntersection(firstRoot, secondRoot)) {
      return; // tangential graze: report neither hit so crossing parity stays even
    }

    for (const double candidate : {firstRoot, secondRoot}) {
      if (candidate < minDistance || candidate > maxDistance) {
        continue;
      }
      const double hitU = localOrigin.xCoord + candidate * localDirection.xCoord;
      const double hitV = localOrigin.yCoord + candidate * localDirection.yCoord;
      const double hitHeight = localOrigin.zCoord + candidate * localDirection.zCoord;
      const double hitPhi = std::atan2(hitV, hitU);
      bool onTrimBoundary = false;
      if (mHasWireTrim) {
        if (!pointInTrim(hitPhi, hitHeight, &onTrimBoundary)) {
          continue;
        }
      } else if (!heightInRange(hitHeight) || !phiInSweep(hitPhi)) {
        continue;
      }
      const double radialDistance = std::hypot(hitU, hitV);
      const Vec3 hitNormal = (mAxisU * (hitU / radialDistance) + mAxisV * (hitV / radialDistance)) * mNormalSign;
      hits.push_back({candidate, hitNormal, onTrimBoundary});
    }
  }

  /// Distance to the trimmed patch. For the parametric rectangle this is exact: the (rho, h) and
  /// phi contributions separate, so the in-sweep case reduces to a 2D point/segment distance in
  /// the (rho, h) half-plane and the off-sweep case to the 3D distance to the nearest straight
  /// seam edge. For a wire trim the same value is a conservative lower bound (the wire region is
  /// contained in this rectangle window), which keeps Safety safe.
  double distanceSqToPatch(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    if (radialDistance <= kTolerance || phiInSweep(std::atan2(localPoint.yCoord, localPoint.xCoord))) {
      return pointSegmentDistanceSq(Vec2{radialDistance, localPoint.zCoord}, Vec2{mRadius, mHeightMin},
                                    Vec2{mRadius, mHeightMax});
    }
    const double distanceToStartSeam =
      pointSegmentDistanceSq(point, pointAt(mPhiStart, mHeightMin), pointAt(mPhiStart, mHeightMax));
    const double endPhi = mPhiStart + mPhiSweep;
    const double distanceToEndSeam =
      pointSegmentDistanceSq(point, pointAt(endPhi, mHeightMin), pointAt(endPhi, mHeightMax));
    return std::min(distanceToStartSeam, distanceToEndSeam);
  }

  Vec3 normalAt(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    if (radialDistance <= kTolerance) {
      return mAxisU * mNormalSign; // ill-defined on the axis; return a stable direction
    }
    return (mAxisU * (localPoint.xCoord / radialDistance) + mAxisV * (localPoint.yCoord / radialDistance)) *
           mNormalSign;
  }

  /// (u, v) = (phi[rad], h[cm]), so X_phi has length r and X_h is the unit axis. Constant over
  /// the domain, and orthogonal, but the radius factor is exactly the one whose absence made a
  /// single parametric tolerance mean 1e-5 cm on a 1 cm cylinder and 1e-3 cm on a 100 cm one.
  void parametricMetric(const Vec2&, double& gUU, double& gUV, double& gVV) const override
  {
    cylinderParametricMetric(mRadius, gUU, gUV, gVV);
  }

  /// Exact divergence-theorem contribution over the (phi, h) rectangle.
  ///
  /// For a wire trim the same integrand is reduced by Green's theorem to a contour integral around
  /// the trim wire. With X = C + r(cos phi U + sin phi V) + h W and n = s(cos phi U + sin phi V),
  ///
  ///     f(phi, h) = (1/3) X . n |X_phi x X_h| = (s r / 3) (a cos phi + b sin phi + r)
  ///
  /// with a = C.U and b = C.V -- and no dependence on h at all, because the h W term of X is
  /// orthogonal to the normal. So the antiderivative in phi is elementary:
  ///
  ///     F(phi, h) = (s r / 3) (a sin phi - b cos phi + r phi)
  ///
  /// and the double integral collapses to the contour integral of F dh around the trim.
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      const double centreU = dot(mCenter, mAxisU);
      const double centreV = dot(mCenter, mAxisV);
      const double factor = mNormalSign * mRadius / 3.;
      return integrateOverCurveTrimByParts(mTrimOuter, mTrimInner, [=](double phi, double) {
        return factor * (centreU * std::sin(phi) - centreV * std::cos(phi) + mRadius * phi);
      });
    }
    const double endPhi = mPhiStart + mPhiSweep;
    const double phiFactor = dot(mCenter, mAxisU) * (std::sin(endPhi) - std::sin(mPhiStart)) -
                             dot(mCenter, mAxisV) * (std::cos(endPhi) - std::cos(mPhiStart));
    const double height = mHeightMax - mHeightMin;
    return mNormalSign * mRadius * height * (phiFactor + mRadius * mPhiSweep) / 3.;
  }

  bool capacityIsExact() const override { return !mHasWireTrim; }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    // conservative: the AABB of the two full rim circles (partial sweeps get a larger box)
    for (const double height : {mHeightMin, mHeightMax}) {
      const Vec3 rimCenter = mCenter + mAxisW * height;
      for (int dimension = 0; dimension < 3; ++dimension) {
        const double radialExtent = mRadius * std::hypot(component(mAxisU, dimension), component(mAxisV, dimension));
        const double centerValue = component(rimCenter, dimension);
        if (dimension == 0) {
          lower.xCoord = std::min(lower.xCoord, centerValue - radialExtent);
          upper.xCoord = std::max(upper.xCoord, centerValue + radialExtent);
        } else if (dimension == 1) {
          lower.yCoord = std::min(lower.yCoord, centerValue - radialExtent);
          upper.yCoord = std::max(upper.yCoord, centerValue + radialExtent);
        } else {
          lower.zCoord = std::min(lower.zCoord, centerValue - radialExtent);
          upper.zCoord = std::max(upper.zCoord, centerValue + radialExtent);
        }
      }
    }
  }

  /// Number of chord segments used for rim sampling, consistent with CurveWire::sampledBoundary
  /// so shared circular boundaries close against curved planar caps.
  int rimSegments() const
  {
    return std::max(1, static_cast<int>(std::lround(kArcSamples * mPhiSweep / kTwoPi)));
  }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    if (mHasWireTrim) {
      appendCurveTrimMesh(mTrimOuter, [this](double phi, double height) { return pointAt(phi, height); }, vertices,
                          triangles);
      return;
    }
    const int segments = rimSegments();
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (int step = 0; step <= segments; ++step) {
      const double phi = mPhiStart + mPhiSweep * step / segments;
      vertices.push_back(pointAt(phi, mHeightMin));
      vertices.push_back(pointAt(phi, mHeightMax));
    }
    for (int step = 0; step < segments; ++step) {
      const int base = firstVertexIndex + 2 * step;
      triangles.push_back({base, base + 2, base + 3});
      triangles.push_back({base, base + 3, base + 1});
    }
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    if (mHasWireTrim) {
      // the (phi, h) -> 3D map is orientation-consistent with the outward normal, so a CCW trim
      // loop yields a CCW 3D loop for an outer wall; the sign is just mNormalSign
      appendCurveTrimEdges(mTrimOuter, mTrimInner,
                           [this](double phi, double height) { return pointAt(phi, height); }, mNormalSign, edges);
      return;
    }
    // patch boundary traversed counter-clockwise as seen along the outward normal: for an
    // outer wall the bottom rim runs with +phi and the top rim against it (reversed for an
    // inner wall), so shared rims with correctly-oriented caps cancel in the half-edge check
    const int segments = rimSegments();
    auto emitEdge = [&](const Vec3& edgeStart, const Vec3& edgeEnd) {
      if (mNormalSign > 0.) {
        edges.emplace_back(edgeStart, edgeEnd);
      } else {
        edges.emplace_back(edgeEnd, edgeStart);
      }
    };
    for (int step = 0; step < segments; ++step) {
      const double phi = mPhiStart + mPhiSweep * step / segments;
      const double nextPhi = mPhiStart + mPhiSweep * (step + 1) / segments;
      emitEdge(pointAt(phi, mHeightMin), pointAt(nextPhi, mHeightMin));
      emitEdge(pointAt(nextPhi, mHeightMax), pointAt(phi, mHeightMax));
    }
    if (!fullSweep()) {
      const double endPhi = mPhiStart + mPhiSweep;
      emitEdge(pointAt(endPhi, mHeightMin), pointAt(endPhi, mHeightMax));
      emitEdge(pointAt(mPhiStart, mHeightMax), pointAt(mPhiStart, mHeightMin));
    }
  }

 private:
  Vec3 mCenter;
  Vec3 mAxisU;
  Vec3 mAxisV;
  Vec3 mAxisW;
  double mRadius = 0.;
  double mHeightMin = 0.;
  double mHeightMax = 0.;
  double mPhiStart = 0.;
  double mPhiSweep = kTwoPi;
  double mNormalSign = 1.;
  bool mHasWireTrim = false;
  CurveWire mTrimOuter;
  std::vector<CurveWire> mTrimInner;
};

/// A bounded spherical surface: a sphere of given radius around a center, trimmed to a
/// parametric rectangle (polar angle theta from the +axisW pole x phi angular sweep), matching
/// TGeoSphere sections. A full sphere (theta in [0, pi], phi sweep 2pi) is self-closing and has
/// no boundary edges. An innerWall surface (e.g. the inner shell of a hollow ball) has its
/// outward-of-solid normal pointing towards the center.
class SphericalBoundedSurface final : public BoundedSurface
{
 public:
  bool initialize(const Vec3& center, const Vec3& polarAxis, const Vec3& referenceAxisU, double radius,
                  double thetaMin, double thetaMax, double phiStart, double phiSweep, bool innerWall,
                  std::string& errorMessage)
  {
    if (!finite(center) || !finite(polarAxis) || !finite(referenceAxisU) || !std::isfinite(radius) ||
        !std::isfinite(thetaMin) || !std::isfinite(thetaMax) || !std::isfinite(phiStart) ||
        !std::isfinite(phiSweep)) {
      errorMessage = "spherical surface parameter is non-finite";
      return false;
    }
    if (radius <= kTolerance) {
      errorMessage = "spherical surface needs a positive radius";
      return false;
    }
    if (thetaMin < -kTolerance || thetaMax > kPi + kTolerance || thetaMax - thetaMin <= kTolerance) {
      errorMessage = "spherical surface needs a polar range within [0, pi]";
      return false;
    }
    if (phiSweep <= kTolerance || phiSweep > kTwoPi + kTolerance) {
      errorMessage = "spherical surface needs an angular sweep in (0, 2pi]";
      return false;
    }
    if (!CylindricalBoundedSurface::makeFrame(polarAxis, referenceAxisU, mAxisU, mAxisV, mAxisW, errorMessage)) {
      return false;
    }

    mCenter = center;
    mRadius = radius;
    mThetaMin = std::max(0., thetaMin);
    mThetaMax = std::min(kPi, thetaMax);
    mPhiStart = phiStart;
    mPhiSweep = std::min(phiSweep, kTwoPi);
    mNormalSign = innerWall ? -1. : 1.;
    return true;
  }

  /// Wire-trimmed overload: the scalar arguments define the sphere surface (frame, radius) and a
  /// nominal window; the outer/inner Curve2D loops in the (phi[rad], theta[rad]) domain are
  /// authoritative for containment. The conservative window is tightened to the outer loop's
  /// parametric bounds (u = phi, v = theta).
  bool initialize(const Vec3& center, const Vec3& polarAxis, const Vec3& referenceAxisU, double radius,
                  double thetaMin, double thetaMax, double phiStart, double phiSweep, bool innerWall,
                  const std::vector<Curve2D>& outerTrim, const std::vector<std::vector<Curve2D>>& innerTrims,
                  std::string& errorMessage)
  {
    if (!initialize(center, polarAxis, referenceAxisU, radius, thetaMin, thetaMax, phiStart, phiSweep, innerWall,
                    errorMessage)) {
      return false;
    }
    Vec2 lower, upper;
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage,
                        parametricMetricOf(*this))) {
      return false;
    }
    mPhiStart = lower.uCoord;
    mPhiSweep = std::min(kTwoPi, upper.uCoord - lower.uCoord);
    mThetaMin = std::max(0., lower.vCoord);
    mThetaMax = std::min(kPi, upper.vCoord);
    mHasWireTrim = true;
    return true;
  }

  bool hasWireTrim() const { return mHasWireTrim; }

  /// True if the (phi, theta) point lies in the trim wire (phi unwrapped into the wire window).
  bool pointInTrim(double phi, double theta, bool* boundary = nullptr) const
  {
    const double uCoord = unwrapAngleInto(phi, mPhiStart, mPhiStart + mPhiSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, theta}, boundary, parametricMetricOf(*this));
  }

  bool fullSweep() const { return mPhiSweep >= kTwoPi - kTolerance; }

  bool fullPolarRange() const
  {
    const double thetaTolerance = angularTolerance(mRadius);
    return mThetaMin <= thetaTolerance && mThetaMax >= kPi - thetaTolerance;
  }

  Vec3 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - mCenter;
    return {dot(relativePoint, mAxisU), dot(relativePoint, mAxisV), dot(relativePoint, mAxisW)};
  }

  bool directionInTrim(const Vec3& localPoint, bool* boundary = nullptr) const
  {
    if (boundary != nullptr) {
      *boundary = false;
    }
    const double pointRadius = norm(localPoint);
    if (pointRadius <= kTolerance) {
      return true; // the center is angle-degenerate; every patch point is equidistant
    }
    const double thetaTolerance = angularTolerance(mRadius);
    const double theta = std::acos(std::max(-1., std::min(1., localPoint.zCoord / pointRadius)));
    const double transverseDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    if (mHasWireTrim) {
      if (transverseDistance <= kTolerance) {
        // on the polar axis phi is degenerate; accept by the wire's theta (v) range
        return theta >= mThetaMin - thetaTolerance && theta <= mThetaMax + thetaTolerance;
      }
      return pointInTrim(std::atan2(localPoint.yCoord, localPoint.xCoord), theta, boundary);
    }
    if (theta < mThetaMin - thetaTolerance || theta > mThetaMax + thetaTolerance) {
      return false;
    }
    if (transverseDistance <= kTolerance) {
      return true; // on the polar axis phi is degenerate
    }
    return angleInSweepRange(std::atan2(localPoint.yCoord, localPoint.xCoord), mPhiStart, mPhiSweep,
                             thetaTolerance);
  }

  Vec3 pointAt(double theta, double phi) const
  {
    const double sinTheta = std::sin(theta);
    return mCenter + (mAxisU * (sinTheta * std::cos(phi)) + mAxisV * (sinTheta * std::sin(phi)) +
                      mAxisW * std::cos(theta)) *
                       mRadius;
  }

  bool containsPointOnSurface(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    if (std::abs(norm(localPoint) - mRadius) > kTolerance) {
      return false;
    }
    return directionInTrim(localPoint);
  }

  void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                           double maxDistance, std::vector<RayHit>& hits) const override
  {
    const Vec3 relativeOrigin = rayOrigin - mCenter;
    const double quadraticA = normSq(rayDirection);
    if (quadraticA <= kToleranceSq) {
      return;
    }
    const double quadraticB = 2. * dot(relativeOrigin, rayDirection);
    const double quadraticC = normSq(relativeOrigin) - mRadius * mRadius;
    const double discriminant = quadraticB * quadraticB - 4. * quadraticA * quadraticC;
    if (discriminant <= 0.) {
      return;
    }
    const double sqrtDiscriminant = std::sqrt(discriminant);
    const double firstRoot = (-quadraticB - sqrtDiscriminant) / (2. * quadraticA);
    const double secondRoot = (-quadraticB + sqrtDiscriminant) / (2. * quadraticA);
    if (sameIntersection(firstRoot, secondRoot)) {
      return; // tangential graze
    }

    for (const double candidate : {firstRoot, secondRoot}) {
      if (candidate < minDistance || candidate > maxDistance) {
        continue;
      }
      const Vec3 localHit = toLocal(rayOrigin + rayDirection * candidate);
      bool onTrimBoundary = false;
      if (!directionInTrim(localHit, &onTrimBoundary)) {
        continue;
      }
      hits.push_back({candidate,
                      (mAxisU * localHit.xCoord + mAxisV * localHit.yCoord + mAxisW * localHit.zCoord) *
                        (mNormalSign / mRadius),
                      onTrimBoundary});
    }
  }

  /// Distance to the patch: exact when the point's direction lies inside the trim (or for a
  /// full sphere); otherwise the full-sphere distance is returned as a conservative lower
  /// bound, which keeps Safety() safe (never larger than the true distance).
  double distanceSqToPatch(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double radialOffset = norm(localPoint) - mRadius;
    return radialOffset * radialOffset;
  }

  Vec3 normalAt(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double pointRadius = norm(localPoint);
    if (pointRadius <= kTolerance) {
      return mAxisW * mNormalSign; // ill-defined at the center; return a stable direction
    }
    return (mAxisU * localPoint.xCoord + mAxisV * localPoint.yCoord + mAxisW * localPoint.zCoord) *
           (mNormalSign / pointRadius);
  }

  /// (u, v) = (phi[rad], theta[rad]) -- the trim domain's order, which is the transpose of
  /// pointAt's argument order. Both coordinates are angles, so the scale is the radius in each,
  /// but the azimuth's shrinks with the parallel it sits on: gUU vanishes at either pole, where a
  /// phi separation genuinely spans no distance at all.
  void parametricMetric(const Vec2& uv, double& gUU, double& gUV, double& gVV) const override
  {
    sphereParametricMetric(mRadius, uv.vCoord, gUU, gUV, gVV);
  }

  /// Exact divergence-theorem contribution over the (theta, phi) rectangle.
  ///
  /// For a wire trim, Green's theorem in the (u, v) = (phi, theta) trim domain. With
  /// X = C + R(sin theta cos phi U + sin theta sin phi V + cos theta W) and n = s (X - C)/R,
  ///
  ///     f(phi, theta) = (s R^2 / 3) sin theta (a sin theta cos phi + b sin theta sin phi
  ///                                            + c cos theta + R)
  ///
  /// so integrating in phi at fixed theta,
  ///
  ///     F(phi, theta) = (s R^2 / 3) sin theta (a sin theta sin phi - b sin theta cos phi
  ///                                            + (c cos theta + R) phi).
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      const double centreU = dot(mCenter, mAxisU);
      const double centreV = dot(mCenter, mAxisV);
      const double centreW = dot(mCenter, mAxisW);
      const double factor = mNormalSign * mRadius * mRadius / 3.;
      return integrateOverCurveTrimByParts(mTrimOuter, mTrimInner, [=](double phi, double theta) {
        const double sinTheta = std::sin(theta);
        return factor * sinTheta *
               (sinTheta * (centreU * std::sin(phi) - centreV * std::cos(phi)) +
                (centreW * std::cos(theta) + mRadius) * phi);
      });
    }
    const double endPhi = mPhiStart + mPhiSweep;
    const double phiFactor = dot(mCenter, mAxisU) * (std::sin(endPhi) - std::sin(mPhiStart)) -
                             dot(mCenter, mAxisV) * (std::cos(endPhi) - std::cos(mPhiStart));
    const double thetaIntegralSinSq =
      0.5 * ((mThetaMax - std::sin(mThetaMax) * std::cos(mThetaMax)) -
             (mThetaMin - std::sin(mThetaMin) * std::cos(mThetaMin)));
    const double thetaIntegralSinCos =
      0.5 * (std::sin(mThetaMax) * std::sin(mThetaMax) - std::sin(mThetaMin) * std::sin(mThetaMin));
    const double thetaIntegralSin = std::cos(mThetaMin) - std::cos(mThetaMax);
    return mNormalSign * mRadius * mRadius *
           (phiFactor * thetaIntegralSinSq + dot(mCenter, mAxisW) * mPhiSweep * thetaIntegralSinCos +
            mRadius * mPhiSweep * thetaIntegralSin) /
           3.;
  }

  bool capacityIsExact() const override { return !mHasWireTrim; }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    lower.xCoord = std::min(lower.xCoord, mCenter.xCoord - mRadius);
    lower.yCoord = std::min(lower.yCoord, mCenter.yCoord - mRadius);
    lower.zCoord = std::min(lower.zCoord, mCenter.zCoord - mRadius);
    upper.xCoord = std::max(upper.xCoord, mCenter.xCoord + mRadius);
    upper.yCoord = std::max(upper.yCoord, mCenter.yCoord + mRadius);
    upper.zCoord = std::max(upper.zCoord, mCenter.zCoord + mRadius);
  }

  int phiSegments() const
  {
    return std::max(1, static_cast<int>(std::lround(kArcSamples * mPhiSweep / kTwoPi)));
  }

  int thetaSegments() const
  {
    return std::max(1, static_cast<int>(std::lround(kArcSamples * (mThetaMax - mThetaMin) / kTwoPi)));
  }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    if (mHasWireTrim) {
      appendCurveTrimMesh(mTrimOuter, [this](double phi, double theta) { return pointAt(theta, phi); }, vertices,
                          triangles);
      return;
    }
    const int phiSteps = phiSegments();
    const int thetaSteps = thetaSegments();
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (int thetaStep = 0; thetaStep <= thetaSteps; ++thetaStep) {
      const double theta = mThetaMin + (mThetaMax - mThetaMin) * thetaStep / thetaSteps;
      for (int phiStep = 0; phiStep <= phiSteps; ++phiStep) {
        vertices.push_back(pointAt(theta, mPhiStart + mPhiSweep * phiStep / phiSteps));
      }
    }
    const int rowLength = phiSteps + 1;
    for (int thetaStep = 0; thetaStep < thetaSteps; ++thetaStep) {
      for (int phiStep = 0; phiStep < phiSteps; ++phiStep) {
        const int base = firstVertexIndex + thetaStep * rowLength + phiStep;
        triangles.push_back({base, base + 1, base + rowLength + 1});
        triangles.push_back({base, base + rowLength + 1, base + rowLength});
      }
    }
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    if (mHasWireTrim) {
      // the (phi, theta) -> 3D map is orientation-*reversed* relative to the outward normal
      // (X_phi x X_theta points inward), so the sign is -mNormalSign
      appendCurveTrimEdges(mTrimOuter, mTrimInner, [this](double phi, double theta) { return pointAt(theta, phi); },
                           -mNormalSign, edges);
      return;
    }
    // boundary of the (theta, phi) rectangle, traversed counter-clockwise for an outer wall;
    // pole rims are degenerate points and full-sweep phi seams cancel, so both are skipped
    auto emitEdge = [&](const Vec3& edgeStart, const Vec3& edgeEnd) {
      if (mNormalSign > 0.) {
        edges.emplace_back(edgeStart, edgeEnd);
      } else {
        edges.emplace_back(edgeEnd, edgeStart);
      }
    };
    const double thetaTolerance = angularTolerance(mRadius);
    const int phiSteps = phiSegments();
    const double endPhi = mPhiStart + mPhiSweep;
    if (mThetaMin > thetaTolerance) {
      for (int step = 0; step < phiSteps; ++step) {
        const double phi = mPhiStart + mPhiSweep * step / phiSteps;
        const double nextPhi = mPhiStart + mPhiSweep * (step + 1) / phiSteps;
        emitEdge(pointAt(mThetaMin, nextPhi), pointAt(mThetaMin, phi)); // -phi at the small-theta rim
      }
    }
    if (mThetaMax < kPi - thetaTolerance) {
      for (int step = 0; step < phiSteps; ++step) {
        const double phi = mPhiStart + mPhiSweep * step / phiSteps;
        const double nextPhi = mPhiStart + mPhiSweep * (step + 1) / phiSteps;
        emitEdge(pointAt(mThetaMax, phi), pointAt(mThetaMax, nextPhi)); // +phi at the large-theta rim
      }
    }
    if (!fullSweep()) {
      const int thetaSteps = thetaSegments();
      for (int step = 0; step < thetaSteps; ++step) {
        const double theta = mThetaMin + (mThetaMax - mThetaMin) * step / thetaSteps;
        const double nextTheta = mThetaMin + (mThetaMax - mThetaMin) * (step + 1) / thetaSteps;
        emitEdge(pointAt(theta, mPhiStart), pointAt(nextTheta, mPhiStart));  // +theta at phiStart
        emitEdge(pointAt(nextTheta, endPhi), pointAt(theta, endPhi));        // -theta at phiEnd
      }
    }
  }

 private:
  Vec3 mCenter;
  Vec3 mAxisU;
  Vec3 mAxisV;
  Vec3 mAxisW;
  double mRadius = 0.;
  double mThetaMin = 0.;
  double mThetaMax = kPi;
  double mPhiStart = 0.;
  double mPhiSweep = kTwoPi;
  double mNormalSign = 1.;
  bool mHasWireTrim = false;
  CurveWire mTrimOuter;
  std::vector<CurveWire> mTrimInner;
};

/// A bounded conical surface: radius varies linearly with the height along the axis,
/// r(h) = radiusAtMin + slope * (h - heightMin), trimmed to a parametric rectangle
/// (phi angular sweep x height range), matching TGeoCone/TGeoConeSeg lateral faces. The radius
/// may reach zero at one end (an apex cone); the mirror nappe of the infinite cone is always
/// rejected. Frame and innerWall conventions are the same as for CylindricalBoundedSurface;
/// slope = 0 degenerates to a cylinder and is allowed.
class ConicalBoundedSurface final : public BoundedSurface
{
 public:
  bool initialize(const Vec3& centerPoint, const Vec3& axis, const Vec3& referenceAxisU, double radiusAtMin,
                  double radiusAtMax, double heightMin, double heightMax, double phiStart, double phiSweep,
                  bool innerWall, std::string& errorMessage)
  {
    if (!finite(centerPoint) || !finite(axis) || !finite(referenceAxisU) || !std::isfinite(radiusAtMin) ||
        !std::isfinite(radiusAtMax) || !std::isfinite(heightMin) || !std::isfinite(heightMax) ||
        !std::isfinite(phiStart) || !std::isfinite(phiSweep)) {
      errorMessage = "conical surface parameter is non-finite";
      return false;
    }
    if (radiusAtMin < -kTolerance || radiusAtMax < -kTolerance ||
        std::max(radiusAtMin, radiusAtMax) <= kTolerance) {
      errorMessage = "conical surface needs non-negative radii, at least one positive";
      return false;
    }
    if (heightMax - heightMin <= kTolerance) {
      errorMessage = "conical surface needs a positive height range";
      return false;
    }
    if (phiSweep <= kTolerance || phiSweep > kTwoPi + kTolerance) {
      errorMessage = "conical surface needs an angular sweep in (0, 2pi]";
      return false;
    }
    if (!CylindricalBoundedSurface::makeFrame(axis, referenceAxisU, mAxisU, mAxisV, mAxisW, errorMessage)) {
      return false;
    }

    mCenter = centerPoint;
    mHeightMin = heightMin;
    mHeightMax = heightMax;
    mSlope = (radiusAtMax - radiusAtMin) / (heightMax - heightMin);
    mRadius0 = radiusAtMin - mSlope * heightMin; // radius at h = 0 of the linear law
    mPhiStart = phiStart;
    mPhiSweep = std::min(phiSweep, kTwoPi);
    mNormalSign = innerWall ? -1. : 1.;
    return true;
  }

  /// Wire-trimmed overload: the scalar arguments define the cone surface (frame, linear radius
  /// law r(h)) and a nominal window; the outer/inner Curve2D loops in the periodic (phi[rad],
  /// h[cm]) domain are authoritative for containment. The radius law is pinned from the scalar
  /// (radiusAtMin, radiusAtMax, heightMin, heightMax) before the window is tightened to the
  /// outer loop's parametric bounds.
  bool initialize(const Vec3& centerPoint, const Vec3& axis, const Vec3& referenceAxisU, double radiusAtMin,
                  double radiusAtMax, double heightMin, double heightMax, double phiStart, double phiSweep,
                  bool innerWall, const std::vector<Curve2D>& outerTrim,
                  const std::vector<std::vector<Curve2D>>& innerTrims, std::string& errorMessage)
  {
    if (!initialize(centerPoint, axis, referenceAxisU, radiusAtMin, radiusAtMax, heightMin, heightMax, phiStart,
                    phiSweep, innerWall, errorMessage)) {
      return false;
    }
    Vec2 lower, upper;
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage,
                        parametricMetricOf(*this))) {
      return false;
    }
    mPhiStart = lower.uCoord;
    mPhiSweep = std::min(kTwoPi, upper.uCoord - lower.uCoord);
    mHeightMin = lower.vCoord;
    mHeightMax = upper.vCoord;
    mHasWireTrim = true;
    return true;
  }

  bool hasWireTrim() const { return mHasWireTrim; }

  /// True if the (phi, h) point lies in the trim wire (phi unwrapped into the wire window).
  bool pointInTrim(double phi, double height, bool* boundary = nullptr) const
  {
    const double uCoord = unwrapAngleInto(phi, mPhiStart, mPhiStart + mPhiSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, height}, boundary, parametricMetricOf(*this));
  }

  bool fullSweep() const { return mPhiSweep >= kTwoPi - kTolerance; }

  double radiusAt(double height) const { return mRadius0 + mSlope * height; }

  double meanRadius() const { return 0.5 * (radiusAt(mHeightMin) + radiusAt(mHeightMax)); }

  Vec3 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - mCenter;
    return {dot(relativePoint, mAxisU), dot(relativePoint, mAxisV), dot(relativePoint, mAxisW)};
  }

  bool heightInRange(double height) const
  {
    return height >= mHeightMin - kTolerance && height <= mHeightMax + kTolerance;
  }

  bool phiInSweep(double phi) const
  {
    return angleInSweepRange(phi, mPhiStart, mPhiSweep, angularTolerance(meanRadius()));
  }

  Vec3 pointAt(double phi, double height) const
  {
    return mCenter + mAxisW * height + (mAxisU * std::cos(phi) + mAxisV * std::sin(phi)) * radiusAt(height);
  }

  bool containsPointOnSurface(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double surfaceRadius = radiusAt(localPoint.zCoord);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    // |rho - r(h)| overestimates the true surface distance by sqrt(1 + slope^2)
    if (std::abs(radialDistance - surfaceRadius) > kTolerance * std::sqrt(1. + mSlope * mSlope)) {
      return false;
    }
    if (radialDistance <= kTolerance) {
      return !mHasWireTrim && heightInRange(localPoint.zCoord); // phi is undefined near the apex
    }
    const double phi = std::atan2(localPoint.yCoord, localPoint.xCoord);
    if (mHasWireTrim) {
      return pointInTrim(phi, localPoint.zCoord);
    }
    return heightInRange(localPoint.zCoord) && phiInSweep(phi);
  }

  void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                           double maxDistance, std::vector<RayHit>& hits) const override
  {
    const Vec3 localOrigin = toLocal(rayOrigin);
    const Vec3 localDirection{dot(rayDirection, mAxisU), dot(rayDirection, mAxisV), dot(rayDirection, mAxisW)};

    // (ox + t dx)^2 + (oy + t dy)^2 = (radius0 + slope * (oz + t dz))^2
    const double surfaceRadiusAtOrigin = mRadius0 + mSlope * localOrigin.zCoord;
    const double quadraticA = localDirection.xCoord * localDirection.xCoord +
                              localDirection.yCoord * localDirection.yCoord -
                              mSlope * mSlope * localDirection.zCoord * localDirection.zCoord;
    const double quadraticB = 2. * (localOrigin.xCoord * localDirection.xCoord +
                                    localOrigin.yCoord * localDirection.yCoord -
                                    mSlope * localDirection.zCoord * surfaceRadiusAtOrigin);
    const double quadraticC = localOrigin.xCoord * localOrigin.xCoord +
                              localOrigin.yCoord * localOrigin.yCoord -
                              surfaceRadiusAtOrigin * surfaceRadiusAtOrigin;

    std::array<double, 2> candidates{};
    int candidateCount = 0;
    if (std::abs(quadraticA) <= kToleranceSq) {
      if (std::abs(quadraticB) <= kToleranceSq) {
        return; // ray runs along the cone surface or its asymptote: no transversal crossing
      }
      candidates[candidateCount++] = -quadraticC / quadraticB;
    } else {
      const double discriminant = quadraticB * quadraticB - 4. * quadraticA * quadraticC;
      if (discriminant <= 0.) {
        return;
      }
      const double sqrtDiscriminant = std::sqrt(discriminant);
      const double firstRoot = (-quadraticB - sqrtDiscriminant) / (2. * quadraticA);
      const double secondRoot = (-quadraticB + sqrtDiscriminant) / (2. * quadraticA);
      if (sameIntersection(firstRoot, secondRoot)) {
        return; // tangential graze (this also covers rays through the exact apex)
      }
      candidates[candidateCount++] = std::min(firstRoot, secondRoot);
      candidates[candidateCount++] = std::max(firstRoot, secondRoot);
    }

    for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex) {
      const double candidate = candidates[candidateIndex];
      if (candidate < minDistance || candidate > maxDistance) {
        continue;
      }
      const double hitHeight = localOrigin.zCoord + candidate * localDirection.zCoord;
      const double hitSurfaceRadius = radiusAt(hitHeight);
      if (hitSurfaceRadius < -kTolerance) {
        continue; // mirror nappe of the infinite cone
      }
      const double hitU = localOrigin.xCoord + candidate * localDirection.xCoord;
      const double hitV = localOrigin.yCoord + candidate * localDirection.yCoord;
      const double radialDistance = std::hypot(hitU, hitV);
      if (radialDistance <= kTolerance) {
        continue; // apex hit: the normal is undefined there
      }
      const double hitPhi = std::atan2(hitV, hitU);
      bool onTrimBoundary = false;
      if (mHasWireTrim) {
        if (!pointInTrim(hitPhi, hitHeight, &onTrimBoundary)) {
          continue;
        }
      } else if (!heightInRange(hitHeight) || !phiInSweep(hitPhi)) {
        continue;
      }
      const double normalScale = mNormalSign / std::sqrt(1. + mSlope * mSlope);
      const Vec3 hitNormal =
        (mAxisU * (hitU / radialDistance) + mAxisV * (hitV / radialDistance) - mAxisW * mSlope) * normalScale;
      hits.push_back({candidate, hitNormal, onTrimBoundary});
    }
  }

  /// Distance to the trimmed patch. For the parametric rectangle this is exact: for an in-sweep
  /// azimuth the problem reduces to the 2D distance to the straight generator segment in the
  /// (rho, h) half-plane; off-sweep the nearest point lies on one of the two straight seam
  /// generators. For a wire trim the same value is a conservative lower bound (the wire region
  /// is contained in this rectangle window), which keeps Safety safe.
  double distanceSqToPatch(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    if (radialDistance <= kTolerance || phiInSweep(std::atan2(localPoint.yCoord, localPoint.xCoord))) {
      return pointSegmentDistanceSq(Vec2{radialDistance, localPoint.zCoord},
                                    Vec2{radiusAt(mHeightMin), mHeightMin}, Vec2{radiusAt(mHeightMax), mHeightMax});
    }
    const double endPhi = mPhiStart + mPhiSweep;
    const double distanceToStartSeam =
      pointSegmentDistanceSq(point, pointAt(mPhiStart, mHeightMin), pointAt(mPhiStart, mHeightMax));
    const double distanceToEndSeam =
      pointSegmentDistanceSq(point, pointAt(endPhi, mHeightMin), pointAt(endPhi, mHeightMax));
    return std::min(distanceToStartSeam, distanceToEndSeam);
  }

  Vec3 normalAt(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    const double normalScale = mNormalSign / std::sqrt(1. + mSlope * mSlope);
    if (radialDistance <= kTolerance) {
      return (mAxisU - mAxisW * mSlope) * normalScale; // ill-defined on the axis; stable fallback
    }
    return (mAxisU * (localPoint.xCoord / radialDistance) + mAxisV * (localPoint.yCoord / radialDistance) -
            mAxisW * mSlope) *
           normalScale;
  }

  /// (u, v) = (phi[rad], h[cm]). Unlike the cylinder this varies along the axis: the azimuthal
  /// scale is the local radius, which reaches zero at an apex, and a step in h also walks along
  /// the slope, so it spans sqrt(1 + slope^2) rather than 1.
  void parametricMetric(const Vec2& uv, double& gUU, double& gUV, double& gVV) const override
  {
    coneParametricMetric(radiusAt(uv.vCoord), mSlope, gUU, gUV, gVV);
  }

  /// Exact divergence-theorem contribution over the (phi, h) rectangle.
  ///
  /// For a wire trim, the same Green's-theorem reduction as the cylinder. The sqrt(1 + slope^2) of
  /// the area element cancels against the same factor in the unit normal, leaving
  ///
  ///     f(phi, h) = (s/3) r(h) (a cos phi + b sin phi + r(h) - m (c + h))
  ///
  /// with a = C.U, b = C.V, c = C.W, m = slope and r(h) = r0 + m h. Only the cos/sin terms carry
  /// phi, so
  ///
  ///     F(phi, h) = (s/3) r(h) (a sin phi - b cos phi + (r(h) - m (c + h)) phi).
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      const double centreU = dot(mCenter, mAxisU);
      const double centreV = dot(mCenter, mAxisV);
      const double centreW = dot(mCenter, mAxisW);
      return integrateOverCurveTrimByParts(mTrimOuter, mTrimInner, [=](double phi, double height) {
        const double localRadius = radiusAt(height);
        return mNormalSign / 3. * localRadius *
               (centreU * std::sin(phi) - centreV * std::cos(phi) +
                (localRadius - mSlope * (centreW + height)) * phi);
      });
    }
    const double endPhi = mPhiStart + mPhiSweep;
    const double phiFactor = dot(mCenter, mAxisU) * (std::sin(endPhi) - std::sin(mPhiStart)) -
                             dot(mCenter, mAxisV) * (std::cos(endPhi) - std::cos(mPhiStart));
    const double radiusIntegral = mRadius0 * (mHeightMax - mHeightMin) +
                                  0.5 * mSlope * (mHeightMax * mHeightMax - mHeightMin * mHeightMin);
    return mNormalSign * radiusIntegral *
           (phiFactor + (mRadius0 - mSlope * dot(mCenter, mAxisW)) * mPhiSweep) / 3.;
  }

  bool capacityIsExact() const override { return !mHasWireTrim; }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    for (const double height : {mHeightMin, mHeightMax}) {
      const Vec3 rimCenter = mCenter + mAxisW * height;
      const double rimRadius = std::max(0., radiusAt(height));
      for (int dimension = 0; dimension < 3; ++dimension) {
        const double radialExtent =
          rimRadius * std::hypot(component(mAxisU, dimension), component(mAxisV, dimension));
        const double centerValue = component(rimCenter, dimension);
        if (dimension == 0) {
          lower.xCoord = std::min(lower.xCoord, centerValue - radialExtent);
          upper.xCoord = std::max(upper.xCoord, centerValue + radialExtent);
        } else if (dimension == 1) {
          lower.yCoord = std::min(lower.yCoord, centerValue - radialExtent);
          upper.yCoord = std::max(upper.yCoord, centerValue + radialExtent);
        } else {
          lower.zCoord = std::min(lower.zCoord, centerValue - radialExtent);
          upper.zCoord = std::max(upper.zCoord, centerValue + radialExtent);
        }
      }
    }
  }

  int rimSegments() const
  {
    return std::max(1, static_cast<int>(std::lround(kArcSamples * mPhiSweep / kTwoPi)));
  }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    if (mHasWireTrim) {
      appendCurveTrimMesh(mTrimOuter, [this](double phi, double height) { return pointAt(phi, height); }, vertices,
                          triangles);
      return;
    }
    const int segments = rimSegments();
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (int step = 0; step <= segments; ++step) {
      const double phi = mPhiStart + mPhiSweep * step / segments;
      vertices.push_back(pointAt(phi, mHeightMin));
      vertices.push_back(pointAt(phi, mHeightMax));
    }
    for (int step = 0; step < segments; ++step) {
      const int base = firstVertexIndex + 2 * step;
      // skip triangles that collapse at an apex rim
      if (radiusAt(mHeightMin) > kTolerance) {
        triangles.push_back({base, base + 2, base + 3});
      }
      if (radiusAt(mHeightMax) > kTolerance) {
        triangles.push_back({base, base + 3, base + 1});
      }
    }
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    if (mHasWireTrim) {
      // the (phi, h) -> 3D map is orientation-consistent with the outward normal (as for the
      // cylinder), so the sign is just mNormalSign
      appendCurveTrimEdges(mTrimOuter, mTrimInner,
                           [this](double phi, double height) { return pointAt(phi, height); }, mNormalSign, edges);
      return;
    }
    // same boundary orientation as the cylinder; an apex rim degenerates to a point and is
    // skipped so an apex cone closes against just one cap
    const int segments = rimSegments();
    auto emitEdge = [&](const Vec3& edgeStart, const Vec3& edgeEnd) {
      if (mNormalSign > 0.) {
        edges.emplace_back(edgeStart, edgeEnd);
      } else {
        edges.emplace_back(edgeEnd, edgeStart);
      }
    };
    for (int step = 0; step < segments; ++step) {
      const double phi = mPhiStart + mPhiSweep * step / segments;
      const double nextPhi = mPhiStart + mPhiSweep * (step + 1) / segments;
      if (radiusAt(mHeightMin) > kTolerance) {
        emitEdge(pointAt(phi, mHeightMin), pointAt(nextPhi, mHeightMin));
      }
      if (radiusAt(mHeightMax) > kTolerance) {
        emitEdge(pointAt(nextPhi, mHeightMax), pointAt(phi, mHeightMax));
      }
    }
    if (!fullSweep()) {
      const double endPhi = mPhiStart + mPhiSweep;
      emitEdge(pointAt(endPhi, mHeightMin), pointAt(endPhi, mHeightMax));
      emitEdge(pointAt(mPhiStart, mHeightMax), pointAt(mPhiStart, mHeightMin));
    }
  }

 private:
  Vec3 mCenter;
  Vec3 mAxisU;
  Vec3 mAxisV;
  Vec3 mAxisW;
  double mRadius0 = 0.;
  double mSlope = 0.;
  double mHeightMin = 0.;
  double mHeightMax = 0.;
  double mPhiStart = 0.;
  double mPhiSweep = kTwoPi;
  double mNormalSign = 1.;
  bool mHasWireTrim = false;
  CurveWire mTrimOuter;
  std::vector<CurveWire> mTrimInner;
};

/// A bounded toroidal surface: a torus of major radius R (axis to tube centre) and minor
/// (tube) radius r about a centre and axis, trimmed to a parametric rectangle in the two
/// periodic angles (phiRing around the main axis x phiTube around the tube). A full torus
/// (both sweeps 2pi) is self-closing and has no boundary edges. An innerWall surface (e.g. the
/// inner tube of a toroidal shell) has its outward-of-solid normal pointing towards the tube
/// spine. The ray/torus intersection is a quartic (solveQuarticReal); the (rho, z) meridian
/// distance to the tube circle gives the exact Safety for the full torus and a conservative
/// lower bound for a trimmed patch, matching the sphere policy.
///
/// Parametrisation (local frame U, V, W = axis): a surface point at (phiRing u, phiTube v) is
///   X(u, v) = centre + (U cos u + V sin u) (R + r cos v) + W (r sin v).
/// The map is orientation-consistent with the outward normal (X_u x X_v = r (R + r cos v) n), so
/// the directed-edge sign is +mNormalSign, as for the cylinder and cone.
class TorusBoundedSurface final : public BoundedSurface
{
 public:
  bool initialize(const Vec3& centerPoint, const Vec3& axis, const Vec3& referenceAxisU, double majorRadius,
                  double minorRadius, double phiStart, double phiSweep, double tubeStart, double tubeSweep,
                  bool innerWall, std::string& errorMessage)
  {
    if (!finite(centerPoint) || !finite(axis) || !finite(referenceAxisU) || !std::isfinite(majorRadius) ||
        !std::isfinite(minorRadius) || !std::isfinite(phiStart) || !std::isfinite(phiSweep) ||
        !std::isfinite(tubeStart) || !std::isfinite(tubeSweep)) {
      errorMessage = "toroidal surface parameter is non-finite";
      return false;
    }
    if (majorRadius <= kTolerance || minorRadius <= kTolerance) {
      errorMessage = "toroidal surface needs positive major and minor radii";
      return false;
    }
    if (phiSweep <= kTolerance || phiSweep > kTwoPi + kTolerance) {
      errorMessage = "toroidal surface needs a ring sweep in (0, 2pi]";
      return false;
    }
    if (tubeSweep <= kTolerance || tubeSweep > kTwoPi + kTolerance) {
      errorMessage = "toroidal surface needs a tube sweep in (0, 2pi]";
      return false;
    }
    if (!CylindricalBoundedSurface::makeFrame(axis, referenceAxisU, mAxisU, mAxisV, mAxisW, errorMessage)) {
      return false;
    }

    mCenter = centerPoint;
    mMajorRadius = majorRadius;
    mMinorRadius = minorRadius;
    mPhiStart = phiStart;
    mPhiSweep = std::min(phiSweep, kTwoPi);
    mTubeStart = tubeStart;
    mTubeSweep = std::min(tubeSweep, kTwoPi);
    mNormalSign = innerWall ? -1. : 1.;
    return true;
  }

  /// Wire-trimmed overload: the scalar arguments define the torus surface (frame, radii) and a
  /// nominal window; the outer/inner Curve2D loops in the periodic (phiRing[rad], phiTube[rad])
  /// domain are authoritative for containment. The conservative window is tightened to the outer
  /// loop's parametric bounds (u = phiRing, v = phiTube). As for the other quadrics the trim must
  /// be a non-wrapping loop within one turn in u (and, for the torus, in v); a trim that wraps in
  /// the tube angle is out of scope and stays on the tessellated fallback.
  bool initialize(const Vec3& centerPoint, const Vec3& axis, const Vec3& referenceAxisU, double majorRadius,
                  double minorRadius, double phiStart, double phiSweep, double tubeStart, double tubeSweep,
                  bool innerWall, const std::vector<Curve2D>& outerTrim,
                  const std::vector<std::vector<Curve2D>>& innerTrims, std::string& errorMessage)
  {
    if (!initialize(centerPoint, axis, referenceAxisU, majorRadius, minorRadius, phiStart, phiSweep, tubeStart,
                    tubeSweep, innerWall, errorMessage)) {
      return false;
    }
    Vec2 lower, upper;
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage,
                        parametricMetricOf(*this))) {
      return false;
    }
    if (upper.vCoord - lower.vCoord > kTwoPi + kTolerance) {
      errorMessage = "toroidal trim wire spans more than a full turn in the tube angle";
      return false;
    }
    mPhiStart = lower.uCoord;
    mPhiSweep = std::min(kTwoPi, upper.uCoord - lower.uCoord);
    mTubeStart = lower.vCoord;
    mTubeSweep = std::min(kTwoPi, upper.vCoord - lower.vCoord);
    mHasWireTrim = true;
    return true;
  }

  bool hasWireTrim() const { return mHasWireTrim; }

  /// True if the (phiRing, phiTube) point lies in the trim wire. Both angles are periodic for a
  /// torus, so each is unwrapped into its wire window (unlike the cylinder/cone/sphere, whose
  /// second parameter - height or theta - is not periodic). A trim that wraps across a full turn
  /// in either angle is rejected at construction, so the window is unambiguous here.
  bool pointInTrim(double phiRing, double phiTube, bool* boundary = nullptr) const
  {
    const double uCoord = unwrapAngleInto(phiRing, mPhiStart, mPhiStart + mPhiSweep);
    const double vCoord = unwrapAngleInto(phiTube, mTubeStart, mTubeStart + mTubeSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, vCoord}, boundary, parametricMetricOf(*this));
  }

  bool fullRingSweep() const { return mPhiSweep >= kTwoPi - kTolerance; }
  bool fullTubeSweep() const { return mTubeSweep >= kTwoPi - kTolerance; }

  Vec3 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - mCenter;
    return {dot(relativePoint, mAxisU), dot(relativePoint, mAxisV), dot(relativePoint, mAxisW)};
  }

  bool ringInSweep(double phiRing) const
  {
    return angleInSweepRange(phiRing, mPhiStart, mPhiSweep, angularTolerance(mMajorRadius));
  }

  bool tubeInSweep(double phiTube) const
  {
    return angleInSweepRange(phiTube, mTubeStart, mTubeSweep, angularTolerance(mMinorRadius));
  }

  Vec3 pointAt(double phiRing, double phiTube) const
  {
    const double ringRadius = mMajorRadius + mMinorRadius * std::cos(phiTube);
    return mCenter + (mAxisU * std::cos(phiRing) + mAxisV * std::sin(phiRing)) * ringRadius +
           mAxisW * (mMinorRadius * std::sin(phiTube));
  }

  /// Unit outward normal (pointing away from the tube spine) from a local surface point.
  Vec3 localNormal(const Vec3& localPoint) const
  {
    const double rho = std::hypot(localPoint.xCoord, localPoint.yCoord);
    if (rho <= kTolerance) {
      return mAxisW * (localPoint.zCoord >= 0. ? mNormalSign : -mNormalSign);
    }
    const double radialFactor = (rho - mMajorRadius) / rho;
    Vec3 normal{radialFactor * localPoint.xCoord, radialFactor * localPoint.yCoord, localPoint.zCoord};
    const double length = norm(normal);
    if (length <= kTolerance) {
      return mAxisU * mNormalSign;
    }
    return (mAxisU * normal.xCoord + mAxisV * normal.yCoord + mAxisW * normal.zCoord) * (mNormalSign / length);
  }

  bool containsPointOnSurface(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double rho = std::hypot(localPoint.xCoord, localPoint.yCoord);
    const double meridianDistance = std::hypot(rho - mMajorRadius, localPoint.zCoord) - mMinorRadius;
    if (std::abs(meridianDistance) > kTolerance) {
      return false;
    }
    const double phiTube = std::atan2(localPoint.zCoord, rho - mMajorRadius);
    if (rho <= kTolerance) {
      return false; // on the axis phiRing is undefined (only reachable on a horn/spindle torus)
    }
    const double phiRing = std::atan2(localPoint.yCoord, localPoint.xCoord);
    if (mHasWireTrim) {
      return pointInTrim(phiRing, phiTube);
    }
    return ringInSweep(phiRing) && tubeInSweep(phiTube);
  }

  void appendIntersections(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                           double maxDistance, std::vector<RayHit>& hits) const override
  {
    const Vec3 localOrigin = toLocal(rayOrigin);
    const Vec3 localDirection{dot(rayDirection, mAxisU), dot(rayDirection, mAxisV), dot(rayDirection, mAxisW)};

    // Torus implicit form (local): (|X|^2 + R^2 - r^2)^2 = 4 R^2 (x^2 + y^2). Substituting the ray
    // X = O + t D gives a quartic in t whose leading coefficient is |D|^4 > 0.
    const double dirDotDir = normSq(localDirection);
    if (dirDotDir <= kToleranceSq) {
      return; // degenerate direction
    }
    const double originDotDir = dot(localOrigin, localDirection);
    const double originDotOrigin = normSq(localOrigin);
    const double constantK = mMajorRadius * mMajorRadius - mMinorRadius * mMinorRadius;
    const double transverseE = localDirection.xCoord * localDirection.xCoord +
                               localDirection.yCoord * localDirection.yCoord;
    const double transverseF = localOrigin.xCoord * localDirection.xCoord +
                               localOrigin.yCoord * localDirection.yCoord;
    const double transverseG = localOrigin.xCoord * localOrigin.xCoord +
                               localOrigin.yCoord * localOrigin.yCoord;
    const double fourRSquared = 4. * mMajorRadius * mMajorRadius;

    const double coeff4 = dirDotDir * dirDotDir;
    const double coeff3 = 4. * dirDotDir * originDotDir;
    const double coeff2 =
      4. * originDotDir * originDotDir + 2. * dirDotDir * (originDotOrigin + constantK) - fourRSquared * transverseE;
    const double coeff1 = 4. * originDotDir * (originDotOrigin + constantK) - 2. * fourRSquared * transverseF;
    const double coeff0 = (originDotOrigin + constantK) * (originDotOrigin + constantK) - fourRSquared * transverseG;

    std::vector<double> candidates = solveQuarticReal(coeff4, coeff3, coeff2, coeff1, coeff0);
    if (candidates.empty()) {
      return;
    }
    std::sort(candidates.begin(), candidates.end());

    // Cluster near-equal roots: an even-sized cluster is a tangential (double) root that must not
    // be reported so crossing parity stays consistent; an odd-sized cluster is a genuine
    // transversal crossing reported once at its mean (the same policy the quadrics apply to their
    // pair of roots via sameIntersection).
    size_t rootIndex = 0;
    while (rootIndex < candidates.size()) {
      size_t clusterEnd = rootIndex + 1;
      double clusterSum = candidates[rootIndex];
      while (clusterEnd < candidates.size() && sameIntersection(candidates[clusterEnd], candidates[clusterEnd - 1])) {
        clusterSum += candidates[clusterEnd];
        ++clusterEnd;
      }
      const size_t clusterSize = clusterEnd - rootIndex;
      rootIndex = clusterEnd;
      if ((clusterSize & 1u) == 0u) {
        continue; // tangential graze
      }
      const double candidate = clusterSum / static_cast<double>(clusterSize);
      if (candidate < minDistance || candidate > maxDistance) {
        continue;
      }
      const Vec3 localHit = toLocal(rayOrigin + rayDirection * candidate);
      const double rho = std::hypot(localHit.xCoord, localHit.yCoord);
      if (rho <= kTolerance) {
        continue;
      }
      const double phiTube = std::atan2(localHit.zCoord, rho - mMajorRadius);
      const double phiRing = std::atan2(localHit.yCoord, localHit.xCoord);
      bool onTrimBoundary = false;
      if (mHasWireTrim) {
        if (!pointInTrim(phiRing, phiTube, &onTrimBoundary)) {
          continue;
        }
      } else if (!ringInSweep(phiRing) || !tubeInSweep(phiTube)) {
        continue;
      }
      hits.push_back({candidate, localNormal(localHit), onTrimBoundary});
    }
  }

  /// Distance to the patch: the (rho, z) meridian half-plane through the point contains the tube
  /// circle nearest to it, so |hypot(rho - R, z) - r| is the exact distance to the full torus.
  /// For a trimmed patch it is a conservative lower bound (the patch is a subset of the full
  /// torus, so the true distance is at least this), which keeps Safety safe.
  double distanceSqToPatch(const Vec3& point) const override
  {
    const Vec3 localPoint = toLocal(point);
    const double rho = std::hypot(localPoint.xCoord, localPoint.yCoord);
    const double meridianDistance = std::hypot(rho - mMajorRadius, localPoint.zCoord) - mMinorRadius;
    return meridianDistance * meridianDistance;
  }

  Vec3 normalAt(const Vec3& point) const override { return localNormal(toLocal(point)); }

  /// (u, v) = (phiRing[rad], phiTube[rad]). Both are angles; the tube's scale is the minor radius
  /// everywhere, while the ring's is the distance from the axis, which runs from R - r on the
  /// inside of the torus to R + r on the outside. The product of the two is the Jacobian the
  /// capacity integrand already uses.
  void parametricMetric(const Vec2& uv, double& gUU, double& gUV, double& gVV) const override
  {
    torusParametricMetric(mMajorRadius, mMinorRadius, uv.vCoord, gUU, gUV, gVV);
  }

  /// Exact divergence-theorem contribution (1/3) integral X . n |X_u x X_v| over the
  /// (phiRing, phiTube) rectangle, with |X_u x X_v| = r (R + r cos v).
  ///
  /// For a wire trim, Green's theorem again. Writing rho = R + r cos v and e = cos u U + sin u V,
  /// X = C + rho e + r sin v W and n = s(cos v e + sin v W), so
  ///
  ///     f(u, v) = (s r rho / 3) (cos v (a cos u + b sin u) + c sin v + rho cos v + r sin^2 v)
  ///
  /// with a = C.U, b = C.V, c = C.W. Only the (a cos u + b sin u) term carries u, so
  ///
  ///     F(u, v) = (s r rho / 3) (cos v (a sin u - b cos u)
  ///                              + (c sin v + rho cos v + r sin^2 v) u).
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      const double centreU = dot(mCenter, mAxisU);
      const double centreV = dot(mCenter, mAxisV);
      const double centreW = dot(mCenter, mAxisW);
      return integrateOverCurveTrimByParts(mTrimOuter, mTrimInner, [=](double phiRing, double phiTube) {
        const double cosTube = std::cos(phiTube);
        const double sinTube = std::sin(phiTube);
        const double rho = mMajorRadius + mMinorRadius * cosTube;
        return mNormalSign * mMinorRadius * rho / 3. *
               (cosTube * (centreU * std::sin(phiRing) - centreV * std::cos(phiRing)) +
                (centreW * sinTube + rho * cosTube + mMinorRadius * sinTube * sinTube) * phiRing);
      });
    }
    // Closed form over u in [u0, u1] (ring) and v in [v0, v1] (tube). See BVHSurfaceSolid.md.
    const double majorR = mMajorRadius;
    const double minorR = mMinorRadius;
    const double u0 = mPhiStart, u1 = mPhiStart + mPhiSweep;
    const double v0 = mTubeStart, v1 = mTubeStart + mTubeSweep;
    const double centerU = dot(mCenter, mAxisU);
    const double centerV = dot(mCenter, mAxisV);
    const double centerW = dot(mCenter, mAxisW);
    const double deltaU = u1 - u0;
    const double deltaV = v1 - v0;
    const double sinIntegralU = std::sin(u1) - std::sin(u0);       // integral cos u du
    const double cosIntegralU = std::cos(u0) - std::cos(u1);       // integral sin u du
    const double sinIntegralV = std::sin(v1) - std::sin(v0);       // integral cos v dv
    const double sinFromCosV = std::cos(v0) - std::cos(v1);        // integral sin v dv
    const double cosSquaredV = 0.5 * deltaV + 0.25 * (std::sin(2. * v1) - std::sin(2. * v0)); // integral cos^2 v dv
    const double sinCosV = 0.25 * (std::cos(2. * v0) - std::cos(2. * v1));                    // integral sin v cos v dv

    // centre-independent part, integrated over v then multiplied by the ring span
    const double centerlessV =
      minorR * ((majorR * majorR + minorR * minorR) * sinIntegralV + majorR * minorR * deltaV +
                majorR * minorR * cosSquaredV);
    // W component of the centre offset
    const double centerWpart = minorR * (majorR * sinFromCosV + minorR * sinCosV);
    // U/V components of the centre offset (ring-angle dependent)
    const double centerUVpart =
      (centerU * sinIntegralU + centerV * cosIntegralU) * minorR * (majorR * sinIntegralV + minorR * cosSquaredV);

    const double total = deltaU * centerlessV + deltaU * centerW * centerWpart + centerUVpart;
    return mNormalSign * total / 3.;
  }

  bool capacityIsExact() const override { return !mHasWireTrim; }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    // conservative: the AABB of the full torus (partial sweeps get a larger box)
    const double outerRadius = mMajorRadius + mMinorRadius;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double radialExtent = outerRadius * std::hypot(component(mAxisU, dimension), component(mAxisV, dimension)) +
                                  mMinorRadius * std::abs(component(mAxisW, dimension));
      const double centerValue = component(mCenter, dimension);
      if (dimension == 0) {
        lower.xCoord = std::min(lower.xCoord, centerValue - radialExtent);
        upper.xCoord = std::max(upper.xCoord, centerValue + radialExtent);
      } else if (dimension == 1) {
        lower.yCoord = std::min(lower.yCoord, centerValue - radialExtent);
        upper.yCoord = std::max(upper.yCoord, centerValue + radialExtent);
      } else {
        lower.zCoord = std::min(lower.zCoord, centerValue - radialExtent);
        upper.zCoord = std::max(upper.zCoord, centerValue + radialExtent);
      }
    }
  }

  int ringSegments() const
  {
    return std::max(1, static_cast<int>(std::lround(kArcSamples * mPhiSweep / kTwoPi)));
  }

  int tubeSegments() const
  {
    return std::max(1, static_cast<int>(std::lround(kArcSamples * mTubeSweep / kTwoPi)));
  }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    if (mHasWireTrim) {
      appendCurveTrimMesh(mTrimOuter, [this](double phiRing, double phiTube) { return pointAt(phiRing, phiTube); },
                          vertices, triangles);
      return;
    }
    const int ringSteps = ringSegments();
    const int tubeSteps = tubeSegments();
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (int ringStep = 0; ringStep <= ringSteps; ++ringStep) {
      const double phiRing = mPhiStart + mPhiSweep * ringStep / ringSteps;
      for (int tubeStep = 0; tubeStep <= tubeSteps; ++tubeStep) {
        vertices.push_back(pointAt(phiRing, mTubeStart + mTubeSweep * tubeStep / tubeSteps));
      }
    }
    const int rowLength = tubeSteps + 1;
    for (int ringStep = 0; ringStep < ringSteps; ++ringStep) {
      for (int tubeStep = 0; tubeStep < tubeSteps; ++tubeStep) {
        const int base = firstVertexIndex + ringStep * rowLength + tubeStep;
        triangles.push_back({base, base + rowLength, base + rowLength + 1});
        triangles.push_back({base, base + rowLength + 1, base + 1});
      }
    }
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    if (mHasWireTrim) {
      // the (phiRing, phiTube) -> 3D map is orientation-consistent with the outward normal, so
      // the sign is just mNormalSign (as for the cylinder and cone)
      appendCurveTrimEdges(mTrimOuter, mTrimInner,
                           [this](double phiRing, double phiTube) { return pointAt(phiRing, phiTube); }, mNormalSign,
                           edges);
      return;
    }
    // boundary of the (phiRing, phiTube) rectangle traversed counter-clockwise as seen along the
    // outward normal; a full sweep in either angle has no seam there, so it is skipped
    auto emitEdge = [&](const Vec3& edgeStart, const Vec3& edgeEnd) {
      if (mNormalSign > 0.) {
        edges.emplace_back(edgeStart, edgeEnd);
      } else {
        edges.emplace_back(edgeEnd, edgeStart);
      }
    };
    const int ringSteps = ringSegments();
    const int tubeSteps = tubeSegments();
    const double endRing = mPhiStart + mPhiSweep;
    const double endTube = mTubeStart + mTubeSweep;
    if (!fullTubeSweep()) {
      for (int step = 0; step < ringSteps; ++step) {
        const double phiRing = mPhiStart + mPhiSweep * step / ringSteps;
        const double nextRing = mPhiStart + mPhiSweep * (step + 1) / ringSteps;
        emitEdge(pointAt(phiRing, mTubeStart), pointAt(nextRing, mTubeStart)); // +phiRing at tubeStart
        emitEdge(pointAt(nextRing, endTube), pointAt(phiRing, endTube));       // -phiRing at tubeEnd
      }
    }
    if (!fullRingSweep()) {
      for (int step = 0; step < tubeSteps; ++step) {
        const double phiTube = mTubeStart + mTubeSweep * step / tubeSteps;
        const double nextTube = mTubeStart + mTubeSweep * (step + 1) / tubeSteps;
        emitEdge(pointAt(endRing, phiTube), pointAt(endRing, nextTube));       // +phiTube at ringEnd
        emitEdge(pointAt(mPhiStart, nextTube), pointAt(mPhiStart, phiTube));   // -phiTube at ringStart
      }
    }
  }

 private:
  Vec3 mCenter;
  Vec3 mAxisU;
  Vec3 mAxisV;
  Vec3 mAxisW;
  double mMajorRadius = 0.;
  double mMinorRadius = 0.;
  double mPhiStart = 0.;
  double mPhiSweep = kTwoPi;
  double mTubeStart = 0.;
  double mTubeSweep = kTwoPi;
  double mNormalSign = 1.;
  bool mHasWireTrim = false;
  CurveWire mTrimOuter;
  std::vector<CurveWire> mTrimInner;
};

/// A trivial bounded surface (a single 3D triangle) used only by unit tests to exercise the
/// BoundedSurface interface polymorphically without depending on a concrete analytic kernel.
class DummyBoundedSurface final : public BoundedSurface
{
 public:
  DummyBoundedSurface(const Vec3& firstVertex, const Vec3& secondVertex, const Vec3& thirdVertex)
    : mVertices{firstVertex, secondVertex, thirdVertex}
  {
    mNormal = normalized(cross(secondVertex - firstVertex, thirdVertex - firstVertex));
  }

  void conservativeBounds(Vec3& lower, Vec3& upper) const override
  {
    for (const auto& vertex : mVertices) {
      lower.xCoord = std::min(lower.xCoord, vertex.xCoord);
      lower.yCoord = std::min(lower.yCoord, vertex.yCoord);
      lower.zCoord = std::min(lower.zCoord, vertex.zCoord);
      upper.xCoord = std::max(upper.xCoord, vertex.xCoord);
      upper.yCoord = std::max(upper.yCoord, vertex.yCoord);
      upper.zCoord = std::max(upper.zCoord, vertex.zCoord);
    }
  }

  bool containsPointOnSurface(const Vec3&) const override { return false; }

  void appendIntersections(const Vec3&, const Vec3&, double, double, std::vector<RayHit>&) const override {}

  double distanceSqToPatch(const Vec3& point) const override
  {
    double bestDistanceSq = std::numeric_limits<double>::infinity();
    for (int vertexIndex = 0; vertexIndex < 3; ++vertexIndex) {
      bestDistanceSq = std::min(bestDistanceSq, pointSegmentDistanceSq(point, mVertices[vertexIndex],
                                                                       mVertices[(vertexIndex + 1) % 3]));
    }
    return bestDistanceSq;
  }

  Vec3 normalAt(const Vec3&) const override { return mNormal; }

  /// A placeholder triangle carries no parametric domain, so there is nothing to convert; the
  /// identity form keeps the interface total without inventing a scale.
  void parametricMetric(const Vec2&, double& gUU, double& gUV, double& gVV) const override
  {
    gUU = 1.;
    gUV = 0.;
    gVV = 1.;
  }

  double capacityContribution() const override { return 0.; }

  bool capacityIsExact() const override { return false; }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const override
  {
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (const auto& vertex : mVertices) {
      vertices.push_back(vertex);
    }
    triangles.push_back({firstVertexIndex, firstVertexIndex + 1, firstVertexIndex + 2});
  }

  void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const override
  {
    for (int vertexIndex = 0; vertexIndex < 3; ++vertexIndex) {
      edges.emplace_back(mVertices[vertexIndex], mVertices[(vertexIndex + 1) % 3]);
    }
  }

 private:
  std::array<Vec3, 3> mVertices;
  Vec3 mNormal;
};

/// How one rim came out of the closure measurement. The four states are exhaustive, and they are
/// exactly the four the ClosureReport rim counters tally.
enum class RimState {
  Matched = 0, ///< every chord has another face within its match band, traversed the other way
  Reversed,    ///< matched, but the partner traverses the shared curve the same way
  Boundary,    ///< some chord has no other face within its match band
  NonManifold  ///< some chord has two or more other faces within the declared tolerance
};

/// One trim loop of one face, as measureRimClosure saw it.
///
/// The aggregate counters below answer "how many rims are open, and how much length is open". They
/// cannot answer "which rim", and every diagnosis so far has had to reconstruct that from counts
/// and totals by hand (TolerancePolicy.md section 12 is entirely such a reconstruction). These
/// records name it, and say where on it the worst chord sits.
struct RimRecord {
  int surfaceIndex = -1;      ///< the owning face's index in the solid's surface list
  int rimIndexOnSurface = -1; ///< which trim loop of that face, in the order the face emits them
  bool closed = false;        ///< the polyline returns to its own first point
  int chords = 0;
  int unmatchedChords = 0;     ///< of them, how many found no other face within their match band
  double length = 0.;          ///< summed chord length, cm
  double unmatchedLength = 0.; ///< how much of it has no other face within the match band, cm
  /// The largest distance in cm from a chord midpoint of this rim to the nearest chord of a
  /// *different* face, and where that chord's midpoint is. This is the per-rim term of
  /// ClosureReport::maxRimIsolation, and it means what that field means: how alone the loneliest
  /// chord is, not how wide a seam is.
  double maxIsolation = 0.;
  Vec3 maxIsolationPoint{};
  int maxIsolationFace = -1; ///< the face owning the nearest chord there, or -1 if there was none
  RimState state = RimState::Matched;
};

/// Result of validating that a set of bounded surfaces forms a consistently-oriented, closed
/// 2-manifold. This is the solid-level counterpart to per-wire validation and detects missing
/// or reversed faces via a directed half-edge check.
struct ClosureReport {
  bool closed = true;                ///< every boundary edge is shared by exactly two faces
  bool orientationConsistent = true; ///< shared edges are traversed in opposite directions
  int boundaryEdges = 0;             ///< edges present on only one face (e.g. a missing face)
  int nonManifoldEdges = 0;          ///< edges shared by more than two faces
  int reversedEdges = 0;             ///< edges shared by two faces in the same direction
  double signedVolume = 0.;          ///< divergence-theorem volume; positive if normals point out

  /// \name Rim-based measurement
  /// The chord counters above answer "did two faces emit the same vertices"; on real CAD they
  /// answer "no" for reasons that are not gaps (see SurfaceRim), and they answer it in units of
  /// chords, which is why a 4-face solid can report 1418 boundary edges. The fields below measure
  /// the same boundary as *curves*, in centimetres, and count it per rim. **The verdict is theirs**
  /// -- validateClosure derives closed/orientationConsistent from the rim counters and the chord
  /// counters decide nothing (TolerancePolicy.md sections 9 and 10 measured what licenses that).
  ///
  /// A chord counts as matched when another face's chord is within the *match band*: the declared
  /// rimEpsilon plus the two chords' own sagittas, because two polylines sampling one shared curve
  /// disagree by that much whatever the geometry does. Only nonManifoldRims is decided on
  /// rimEpsilon alone (see measureRimClosure).
  /// @{
  /// The largest distance in cm from any rim chord to the nearest chord of a *different* face; 0
  /// when there is nothing to compare. Read it as "how alone is the loneliest chord", which is
  /// what it measures -- **not** as the width of a seam between two faces. On a solid with no open
  /// rim it is bounded by the match band and says nothing; on one with an open rim it is that
  /// rim's isolation. It was read as a face-to-face separation in four places for two sessions,
  /// and the tell that it is not one is that it has no dependence on rimEpsilon at all
  /// (TolerancePolicy.md section 12.3).
  double maxRimIsolation = 0.;
  double totalRimLength = 0.;      ///< summed length in cm of every face's trim boundary
  double unmatchedRimLength = 0.;  ///< how much of it has no other face within the match band, cm
  double rimEpsilon = 0.;          ///< the declared matching tolerance, in cm
  double rimChordResolution = 0.;  ///< the largest amount by which any rim polyline can sit off the
                                   ///< smooth rim it samples, in cm; the per-chord value of this is
                                   ///< what widens the match band
  int rims = 0;                    ///< total number of trim loops over all faces
  int matchedRims = 0;             ///< every chord has another face within the match band,
                                   ///< traversed the opposite way
  int reversedRims = 0;            ///< matched, but the partner traverses the shared curve the
                                   ///< same way (one face's outward normal points inward)
  int nonManifoldRims = 0;         ///< some chord has two or more other faces within rimEpsilon
  int boundaryRims = 0;            ///< some chord has no other face within the match band
  /// One entry per rim, in the order the faces were visited: the detail behind the counters above.
  /// The counters say how many rims are open; these say which, and where.
  std::vector<RimRecord> rimRecords;
  /// @}
};

/// Measure the face-to-face gaps of \a surfaces as curves, filling the rim fields of \a report.
///
/// \a epsilon is the declared match tolerance in cm. Each chord of each rim is probed at its
/// *midpoint* and matched against the chords of every other face: midpoints rather than vertices,
/// because a rim vertex at a box corner legitimately lies on two other faces at once and would
/// read as non-manifold, while a chord interior belongs to exactly one shared edge.
///
/// Cost note: the search is a uniform grid over chords with an expanding-shell nearest query, so
/// it is near-linear when partners are close and degrades only for rims that have no partner
/// anywhere -- which are the rims a solid should not have.
inline void measureRimClosure(const std::vector<std::unique_ptr<BoundedSurface>>& surfaces, double epsilon,
                              ClosureReport& report)
{
  report.rimEpsilon = epsilon;

  std::vector<SurfaceRim> rims;
  std::vector<int> rimIndexOnSurface;
  for (size_t surfaceIndex = 0; surfaceIndex < surfaces.size(); ++surfaceIndex) {
    if (surfaces[surfaceIndex] == nullptr) {
      continue;
    }
    const size_t firstNewRim = rims.size();
    surfaces[surfaceIndex]->appendRims(rims);
    for (size_t rimIndex = firstNewRim; rimIndex < rims.size(); ++rimIndex) {
      rims[rimIndex].surfaceIndex = static_cast<int>(surfaceIndex);
      rimIndexOnSurface.push_back(static_cast<int>(rimIndex - firstNewRim));
    }
  }
  report.rims = static_cast<int>(rims.size());
  if (rims.empty()) {
    return;
  }

  // Flatten to chords, and take each chord's own accuracy while we are here. Two faces sampling one
  // shared curve at different phases differ by the chord sagitta even when the curve is shared
  // exactly, so the sagitta is the floor of what any comparison between two rims can resolve -- and
  // it is therefore the band a chord has to be matched in. Matching below it does not measure the
  // geometry, it measures the sampler: it was a 1e-8 cm band against a 4e-9 cm sagitta that let
  // cyl_cross_cyl read as closed, and the same band against a 7e-6 cm sagitta that made it read as
  // open the moment the B-spline sampler started resolving its seam properly (section 13).
  //
  // Estimating the sagitta needs the smooth rim, which is not available here, so estimate it from
  // the polyline: a vertex whose two chords turn through a small angle is a sample of a smooth run,
  // and for a run of curvature radius r sampled every angle t the chord is 2 r sin(t/2) and the
  // sagitta r (1 - cos(t/2)) = (chord/2) tan(t/4). A vertex whose chords turn through more than
  // kMaxSmoothTurn is a genuine corner of the trim, where the polyline is exact and there is no
  // sagitta at all -- which is why the deviation of the vertex from the chord joining its
  // neighbours is the wrong measure: on a box it reports the corner offset, half a face wide.
  //
  // The bound is kept *per chord*, not just as the solid-wide maximum: a densely sampled seam has
  // no reason to inherit the uncertainty of the coarsest 24-sample arc elsewhere on the solid.
  constexpr double kMaxSmoothTurn = 0.52; // ~30 degrees; a rim sampled at kArcSamples turns by 15
  struct Chord {
    Vec3 start;
    Vec3 end;
    int surfaceIndex;
    double resolution; ///< how far this chord can sit from the smooth rim it samples, in cm
  };
  std::vector<Chord> chords;
  std::vector<std::pair<size_t, size_t>> chordRange(rims.size()); // [first, last) chord of each rim
  for (size_t rimIndex = 0; rimIndex < rims.size(); ++rimIndex) {
    const SurfaceRim& rim = rims[rimIndex];
    chordRange[rimIndex].first = chords.size();
    const size_t pointCount = rim.points.size();
    std::vector<double> vertexSagitta(pointCount, 0.);
    const size_t interiorCount = rim.closed ? pointCount : (pointCount >= 2 ? pointCount - 2 : 0);
    for (size_t offset = 0; offset < interiorCount; ++offset) {
      const size_t middle = rim.closed ? offset : offset + 1;
      const Vec3 incoming = rim.points[middle] - rim.points[(middle + pointCount - 1) % pointCount];
      const Vec3 outgoing = rim.points[(middle + 1) % pointCount] - rim.points[middle];
      const double incomingLength = norm(incoming);
      const double outgoingLength = norm(outgoing);
      if (incomingLength <= kTolerance || outgoingLength <= kTolerance) {
        continue;
      }
      const double turn = std::acos(std::clamp(dot(incoming, outgoing) / (incomingLength * outgoingLength), -1., 1.));
      if (turn > kMaxSmoothTurn) {
        continue; // a corner of the trim, not a sample of a smooth run
      }
      vertexSagitta[middle] = 0.25 * (incomingLength + outgoingLength) * std::tan(0.25 * turn);
      report.rimChordResolution = std::max(report.rimChordResolution, vertexSagitta[middle]);
    }
    const size_t chordCount = rim.closed ? pointCount : pointCount - 1;
    for (size_t pointIndex = 0; pointIndex < chordCount; ++pointIndex) {
      const size_t nextIndex = (pointIndex + 1) % pointCount;
      chords.push_back({rim.points[pointIndex], rim.points[nextIndex], rim.surfaceIndex,
                        std::max(vertexSagitta[pointIndex], vertexSagitta[nextIndex])});
    }
    chordRange[rimIndex].second = chords.size();
  }
  if (chords.empty()) {
    return;
  }

  Vec3 lower{chords.front().start};
  Vec3 upper{chords.front().start};
  auto grow = [&](const Vec3& point) {
    lower = {std::min(lower.xCoord, point.xCoord), std::min(lower.yCoord, point.yCoord),
             std::min(lower.zCoord, point.zCoord)};
    upper = {std::max(upper.xCoord, point.xCoord), std::max(upper.yCoord, point.yCoord),
             std::max(upper.zCoord, point.zCoord)};
  };
  for (const Chord& chord : chords) {
    grow(chord.start);
    grow(chord.end);
  }
  const int gridDimension =
    std::clamp(static_cast<int>(std::cbrt(static_cast<double>(chords.size()))), 1, 32);
  const Vec3 extent = upper - lower;
  const double cellSize =
    std::max({extent.xCoord, extent.yCoord, extent.zCoord, kTolerance}) / gridDimension;
  auto cellOf = [&](double coordinate, double origin) {
    return std::clamp(static_cast<int>((coordinate - origin) / cellSize), 0, gridDimension - 1);
  };
  auto cellIndex = [&](int xCell, int yCell, int zCell) {
    return (xCell * gridDimension + yCell) * gridDimension + zCell;
  };
  std::vector<std::vector<int>> cells(static_cast<size_t>(gridDimension) * gridDimension * gridDimension);
  for (size_t chordIndex = 0; chordIndex < chords.size(); ++chordIndex) {
    const Chord& chord = chords[chordIndex];
    const int xLow = cellOf(std::min(chord.start.xCoord, chord.end.xCoord), lower.xCoord);
    const int xHigh = cellOf(std::max(chord.start.xCoord, chord.end.xCoord), lower.xCoord);
    const int yLow = cellOf(std::min(chord.start.yCoord, chord.end.yCoord), lower.yCoord);
    const int yHigh = cellOf(std::max(chord.start.yCoord, chord.end.yCoord), lower.yCoord);
    const int zLow = cellOf(std::min(chord.start.zCoord, chord.end.zCoord), lower.zCoord);
    const int zHigh = cellOf(std::max(chord.start.zCoord, chord.end.zCoord), lower.zCoord);
    for (int xCell = xLow; xCell <= xHigh; ++xCell) {
      for (int yCell = yLow; yCell <= yHigh; ++yCell) {
        for (int zCell = zLow; zCell <= zHigh; ++zCell) {
          cells[cellIndex(xCell, yCell, zCell)].push_back(static_cast<int>(chordIndex));
        }
      }
    }
  }

  struct Match {
    double distance = std::numeric_limits<double>::infinity();
    int chordIndex = -1;
    /// Some other face's chord lies inside this chord's match band -- the declared tolerance plus
    /// what the two polylines can differ by while sampling the very same curve. This is what
    /// decides matched against open.
    bool withinBand = false;
    /// The distinct faces found within the declared tolerance alone. Room for three is enough:
    /// only none, one and "more than one" are distinguished, and only the last is used.
    std::array<int, 3> coincidentFaces{-1, -1, -1};
    int coincidentFaceCount = 0;
  };
  // Two questions, two bands, and they must not be the same band.
  //
  // "Is this chord's edge shared with another face" is answered in the sampling-aware band: below
  // the sagitta, two polylines of one shared curve disagree for reasons that are not gaps.
  //
  // "Do two other faces occupy this edge at once" is answered in the declared tolerance alone.
  // Widening it to the sagitta as well makes every chord next to a corner non-manifold, because a
  // corner is where a third face's rim legitimately comes within a chord's length -- probing chord
  // midpoints rather than vertices keeps a *vertex* out of that trap, but nothing keeps a whole
  // chord out of it once the band is a hundredth of a centimetre. That widening turned four
  // clean Bagger parts non-manifold and was the entire cost of the first attempt at this (section
  // 13). Neither corpus has a single coincident-face part to calibrate a wider test on, so the
  // strict one stays until one exists.
  const double maxBand = epsilon + 2. * report.rimChordResolution;
  auto nearestOtherFace = [&](const Vec3& probe, int ownSurfaceIndex, double probeResolution) {
    Match match;
    auto consider = [&](int chordIndex) {
      const Chord& chord = chords[static_cast<size_t>(chordIndex)];
      if (chord.surfaceIndex == ownSurfaceIndex) {
        return;
      }
      const double distance = std::sqrt(pointSegmentDistanceSq(probe, chord.start, chord.end));
      if (distance < match.distance) {
        match.distance = distance;
        match.chordIndex = chordIndex;
      }
      if (distance <= epsilon + probeResolution + chord.resolution) {
        match.withinBand = true;
      }
      if (distance <= epsilon && match.coincidentFaceCount < static_cast<int>(match.coincidentFaces.size())) {
        for (int seen = 0; seen < match.coincidentFaceCount; ++seen) {
          if (match.coincidentFaces[static_cast<size_t>(seen)] == chord.surfaceIndex) {
            return;
          }
        }
        match.coincidentFaces[static_cast<size_t>(match.coincidentFaceCount++)] = chord.surfaceIndex;
      }
    };
    const int xCentre = cellOf(probe.xCoord, lower.xCoord);
    const int yCentre = cellOf(probe.yCoord, lower.yCoord);
    const int zCentre = cellOf(probe.zCoord, lower.zCoord);
    for (int shell = 0; shell < gridDimension; ++shell) {
      // every cell strictly inside this shell has already been visited, so once the nearest hit is
      // closer than the shell's own inner distance nothing further out can beat it -- and the
      // shells must still reach the match band, or a partner face inside it could be missed
      const double shellReach = (shell - 1) * cellSize;
      if (shell > 0 && shellReach > std::max(match.distance, maxBand)) {
        break;
      }
      for (int xCell = xCentre - shell; xCell <= xCentre + shell; ++xCell) {
        if (xCell < 0 || xCell >= gridDimension) {
          continue;
        }
        for (int yCell = yCentre - shell; yCell <= yCentre + shell; ++yCell) {
          if (yCell < 0 || yCell >= gridDimension) {
            continue;
          }
          for (int zCell = zCentre - shell; zCell <= zCentre + shell; ++zCell) {
            if (zCell < 0 || zCell >= gridDimension) {
              continue;
            }
            const bool onShell = std::abs(xCell - xCentre) == shell || std::abs(yCell - yCentre) == shell ||
                                 std::abs(zCell - zCentre) == shell;
            if (!onShell) {
              continue; // interior of the shell: visited on an earlier pass
            }
            for (const int chordIndex : cells[cellIndex(xCell, yCell, zCell)]) {
              consider(chordIndex);
            }
          }
        }
      }
    }
    return match;
  };

  report.rimRecords.reserve(rims.size());
  for (size_t rimIndex = 0; rimIndex < rims.size(); ++rimIndex) {
    bool hasUnmatched = false;
    bool hasNonManifold = false;
    int sameDirectionVotes = 0;
    int oppositeDirectionVotes = 0;
    RimRecord record;
    record.surfaceIndex = rims[rimIndex].surfaceIndex;
    record.rimIndexOnSurface = rimIndexOnSurface[rimIndex];
    record.closed = rims[rimIndex].closed;
    record.chords = static_cast<int>(chordRange[rimIndex].second - chordRange[rimIndex].first);
    for (size_t chordIndex = chordRange[rimIndex].first; chordIndex < chordRange[rimIndex].second; ++chordIndex) {
      const Chord& chord = chords[chordIndex];
      const Vec3 along = chord.end - chord.start;
      const double chordLength = norm(along);
      report.totalRimLength += chordLength;
      record.length += chordLength;
      const Vec3 probe = chord.start + along * 0.5;
      const Match match = nearestOtherFace(probe, chord.surfaceIndex, chord.resolution);
      if (std::isfinite(match.distance)) {
        report.maxRimIsolation = std::max(report.maxRimIsolation, match.distance);
        if (match.distance > record.maxIsolation || record.maxIsolationFace < 0) {
          record.maxIsolation = match.distance;
          record.maxIsolationPoint = probe;
          record.maxIsolationFace = chords[static_cast<size_t>(match.chordIndex)].surfaceIndex;
        }
      }
      if (match.coincidentFaceCount > 1) {
        hasNonManifold = true;
      }
      if (!match.withinBand) {
        hasUnmatched = true;
        ++record.unmatchedChords;
        report.unmatchedRimLength += chordLength;
        record.unmatchedLength += chordLength;
        continue;
      }
      const Chord& partner = chords[static_cast<size_t>(match.chordIndex)];
      if (dot(along, partner.end - partner.start) < 0.) {
        ++oppositeDirectionVotes;
      } else {
        ++sameDirectionVotes;
      }
    }
    if (hasNonManifold) {
      ++report.nonManifoldRims;
      record.state = RimState::NonManifold;
    } else if (hasUnmatched) {
      ++report.boundaryRims;
      record.state = RimState::Boundary;
    } else if (sameDirectionVotes > oppositeDirectionVotes) {
      ++report.reversedRims;
      record.state = RimState::Reversed;
    } else {
      ++report.matchedRims;
      record.state = RimState::Matched;
    }
    report.rimRecords.push_back(record);
  }
}

/// Validate closure and orientation of a set of bounded surfaces using a directed half-edge map,
/// and measure the same boundary as curves (see measureRimClosure).
///
/// \a modelTolerance is the source model's own declared tolerance in cm, carried by sidecar
/// version 2; 0 means "not stated" and falls back on kRimMatchTolerance. It feeds the rim
/// measurement only -- the half-edge verdict is unchanged by it, and by everything else here.
inline ClosureReport validateClosure(const std::vector<std::unique_ptr<BoundedSurface>>& surfaces,
                                     double modelTolerance = 0.)
{
  ClosureReport report;

  auto quantize = [](double value) { return static_cast<int64_t>(std::llround(value / kClosureQuantum)); };
  using VertexKey = std::tuple<int64_t, int64_t, int64_t>;
  auto keyOf = [&](const Vec3& point) {
    return VertexKey{quantize(point.xCoord), quantize(point.yCoord), quantize(point.zCoord)};
  };

  std::vector<std::pair<Vec3, Vec3>> directedEdges;
  for (const auto& surface : surfaces) {
    if (surface != nullptr) {
      surface->appendDirectedEdges(directedEdges);
      report.signedVolume += surface->capacityContribution();
    }
  }

  // For each undirected edge, count occurrences in the forward and reverse directions.
  std::map<std::pair<VertexKey, VertexKey>, std::pair<int, int>> edgeCounts;
  for (const auto& directedEdge : directedEdges) {
    const VertexKey startKey = keyOf(directedEdge.first);
    const VertexKey endKey = keyOf(directedEdge.second);
    if (startKey == endKey) {
      continue; // degenerate edge, already flagged at wire level
    }
    const bool forward = startKey < endKey;
    const auto orderedKey = forward ? std::make_pair(startKey, endKey) : std::make_pair(endKey, startKey);
    auto& counts = edgeCounts[orderedKey];
    if (forward) {
      ++counts.first;
    } else {
      ++counts.second;
    }
  }

  for (const auto& [edgeKey, counts] : edgeCounts) {
    const int total = counts.first + counts.second;
    if (total == 1) {
      ++report.boundaryEdges; // missing neighbouring face
    } else if (total == 2) {
      if (counts.first != 1 || counts.second != 1) {
        ++report.reversedEdges; // both faces traverse the edge the same way
      }
    } else {
      ++report.nonManifoldEdges;
    }
  }

  measureRimClosure(surfaces, modelTolerance > 0. ? modelTolerance : kRimMatchTolerance, report);

  // The verdict is the rim measurement's, not the chord counters'. The counters stay, because
  // they are still the cheapest way to see *how* two faces disagree, but they answer the wrong
  // question: whether two faces emitted the same vertices along a shared edge, when each face
  // samples that edge independently and the vertices genuinely are not the same points. A check
  // that cannot say "closed" tells you nothing when it says "open".
  report.closed = (report.boundaryRims == 0) && (report.nonManifoldRims == 0);
  report.orientationConsistent = (report.reversedRims == 0);
  return report;
}

} // namespace o2::base::surface

#endif
