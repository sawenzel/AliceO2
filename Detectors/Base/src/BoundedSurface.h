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
/// Tolerance for accepting a wire loop as closed (consecutive edge/curve endpoints meeting). It is
/// deliberately looser than kTolerance: a wire is imported from a CAD extractor that samples/
/// projects each boundary curve's endpoints independently (e.g. a straight line vertex vs. the end
/// pole of a neighbouring B-spline/arc on a quadric), so the shared vertex can differ by the
/// extractor's ~1e-6 precision even though the loop is closed. The residual gap is negligible for
/// winding/area/navigation, which use kTolerance for on-boundary classification.
inline constexpr double kWireJoinTolerance = 1.e-5;
inline constexpr double kWireJoinToleranceSq = kWireJoinTolerance * kWireJoinTolerance;
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

  /// Build and validate the wire from an implicitly-closed vertex ring.
  bool initialize(const std::vector<Vec2>& inputVertices, WireRole wireRole, WireStatus& status)
  {
    role = wireRole;
    vertices.clear();
    vertices.reserve(inputVertices.size());

    for (const auto& vertex : inputVertices) {
      if (!finite(vertex)) {
        status = WireStatus::NonFinite;
        return false;
      }
      if (vertices.empty() || surface::distanceSq(vertices.back(), vertex) > kToleranceSq) {
        vertices.push_back(vertex);
      }
    }

    // drop an explicit closing duplicate (first == last)
    if (vertices.size() > 1 && surface::distanceSq(vertices.front(), vertices.back()) <= kToleranceSq) {
      vertices.pop_back();
    }

    if (vertices.size() < 3) {
      status = WireStatus::TooFewVertices;
      return false;
    }

    // reject self-touching loops (non-adjacent coincident vertices)
    for (size_t firstIndex = 0; firstIndex < vertices.size(); ++firstIndex) {
      for (size_t secondIndex = firstIndex + 1; secondIndex < vertices.size(); ++secondIndex) {
        if (surface::distanceSq(vertices[firstIndex], vertices[secondIndex]) <= kToleranceSq) {
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
  bool initializeFromEdges(const std::vector<SurfaceEdge>& edges, WireRole wireRole, WireStatus& status)
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
      if (surface::distanceSq(edges[edgeIndex].end, nextStart) > kToleranceSq) {
        status = WireStatus::Open;
        return false;
      }
    }

    std::vector<Vec2> ringVertices;
    ringVertices.reserve(edges.size());
    for (const auto& singleEdge : edges) {
      ringVertices.push_back(singleEdge.start);
    }
    return initialize(ringVertices, wireRole, status);
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

  WireClassification classify(const Vec2& point) const
  {
    bool inside = false;
    for (size_t vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex) {
      const auto& segmentStart = vertices[vertexIndex];
      const auto& segmentEnd = vertices[(vertexIndex + 1) % vertices.size()];
      if (pointSegmentDistanceSq(point, segmentStart, segmentEnd) <= kToleranceSq) {
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
  void bsplineSampleInto(std::vector<Vec2>& samples, double flatnessSq = 1.e-10, int maxDepth = 16) const
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

  void bsplineSampleRecursive(double t0, double t1, const Vec2& p0, const Vec2& p1, double flatnessSq,
                              int depth, std::vector<Vec2>& samples) const
  {
    const double tMid = 0.5 * (t0 + t1);
    Vec2 midPoint;
    Vec2 unusedDerivative;
    bsplineEval(tMid, midPoint, unusedDerivative);
    // stop when flat enough, when the chord is already shorter than the flatness scale (so no
    // finer detail is resolvable), or at the depth cap
    if (depth <= 0 || pointSegmentDistanceSq(midPoint, p0, p1) <= flatnessSq ||
        surface::distanceSq(p0, p1) <= flatnessSq) {
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
      return poles.front(); // clamped knot vector: first pole is the start point
    }
    return pointAtAngle(startAngle);
  }
  Vec2 endPoint() const
  {
    if (kind == CurveKind::Line) {
      return lineEnd;
    }
    if (kind == CurveKind::BSpline) {
      return poles.back(); // clamped knot vector: last pole is the end point
    }
    return pointAtAngle(endAngle);
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
        // substitute the caller's canonical shared endpoints at the two ends without mutating the
        // cached polyline, so loop seams stay consistent across curves
        const Vec2 first = (index == 0) ? canonicalStart : polyline[index];
        const Vec2 second = (index + 2 == polyline.size()) ? canonicalEnd : polyline[index + 1];
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

  /// Build and validate the wire from an ordered, closed list of curves.
  bool initialize(const std::vector<Curve2D>& inputCurves, WireRole wireRole, WireStatus& status)
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
      if (surface::distanceSq(currentEnd, nextStart) > kWireJoinToleranceSq) {
        status = WireStatus::Open;
        return false;
      }
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

  /// Accumulate the loop's exact extent into a parametric axis-aligned bounding box.
  void parametricBounds(Vec2& lower, Vec2& upper) const
  {
    for (const auto& curve : curves) {
      curve.extendBounds(lower, upper);
    }
  }

  /// Classify a parametric point against the closed loop (inside / outside / on-boundary).
  WireClassification classify(const Vec2& point) const
  {
    for (const auto& curve : curves) {
      if (curve.distanceSq(point) <= kToleranceSq) {
        return WireClassification::Boundary;
      }
    }
    // pass canonical shared endpoints (each curve ends exactly where the next starts) so the
    // half-open crossing convention stays consistent across seams despite floating-point drift
    int crossings = 0;
    for (size_t curveIndex = 0; curveIndex < curves.size(); ++curveIndex) {
      const Vec2 canonicalStart = curves[curveIndex].startPoint();
      const Vec2 canonicalEnd = curves[(curveIndex + 1) % curves.size()].startPoint();
      crossings += curves[curveIndex].rightwardCrossings(point, canonicalStart, canonicalEnd);
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
                              const Vec2& point, bool* boundary = nullptr)
{
  if (boundary != nullptr) {
    *boundary = false;
  }
  const auto outerClassification = outerWire.classify(point);
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
    const auto innerClassification = innerWire.classify(point);
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

/// Numerically integrate \a integrand(u, v) over the curve-wire-trimmed region (outer loop
/// minus holes) by midpoint quadrature on a regular grid across the outer loop's parametric
/// bounding box. Used for the (inexact) capacity contribution of wire-trimmed quadric patches;
/// the untrimmed rectangle path keeps its exact closed form.
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
                           std::string& errorMessage)
{
  WireStatus status = WireStatus::Valid;
  if (!outerWire.initialize(outerTrim, WireRole::Outer, status)) {
    errorMessage = std::string("quadric outer trim wire invalid: ") + wireStatusMessage(status);
    return false;
  }
  innerWires.clear();
  innerWires.reserve(innerTrims.size());
  for (const auto& innerLoop : innerTrims) {
    CurveWire innerWire;
    WireStatus innerStatus = WireStatus::Valid;
    if (!innerWire.initialize(innerLoop, WireRole::Inner, innerStatus)) {
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
    errorMessage = "quadric trim wire spans more than a full turn in phi";
    return false;
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

  /// Signed divergence-theorem contribution to the enclosed volume.
  virtual double capacityContribution() const = 0;

  /// Whether capacityContribution() is analytically exact for this surface.
  virtual bool capacityIsExact() const = 0;

  /// Append this patch's visualization triangulation (navigation must never depend on it).
  virtual void appendDisplayMesh(std::vector<Vec3>& vertices,
                                 std::vector<std::array<int, 3>>& triangles) const = 0;

  /// Append the 3D directed boundary edges of the patch, for solid-closure validation.
  virtual void appendDirectedEdges(std::vector<std::pair<Vec3, Vec3>>& edges) const = 0;
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
    if (!mOuterWire.initialize(outerWireVertices, WireRole::Outer, outerStatus)) {
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
      if (!innerWire.initialize(innerWireInput, WireRole::Inner, innerStatus)) {
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

    const auto outerClassification = mOuterWire.classify(point);
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
      const auto innerClassification = innerWire.classify(point);
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
    if (!containsLocal(toLocal(candidatePoint))) {
      return;
    }
    hits.push_back({candidateDistance, mNormal});
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
    if (!mOuterWire.initialize(outerCurves, WireRole::Outer, outerStatus)) {
      errorMessage = std::string("outer wire invalid: ") + wireStatusMessage(outerStatus);
      return false;
    }
    mReoriented = (outerStatus == WireStatus::Reversed);

    mInnerWires.clear();
    mInnerWires.reserve(innerCurves.size());
    for (const auto& innerCurveLoop : innerCurves) {
      CurveWire innerWire;
      WireStatus innerStatus = WireStatus::Valid;
      if (!innerWire.initialize(innerCurveLoop, WireRole::Inner, innerStatus)) {
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

    const auto outerClassification = mOuterWire.classify(point);
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
      const auto innerClassification = innerWire.classify(point);
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
    if (!containsLocal(toLocal(rayOrigin + rayDirection * candidateDistance))) {
      return;
    }
    hits.push_back({candidateDistance, mNormal});
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
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage)) {
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
  bool pointInTrim(double phi, double height) const
  {
    const double uCoord = unwrapAngleInto(phi, mPhiStart, mPhiStart + mPhiSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, height});
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
      if (mHasWireTrim) {
        if (!pointInTrim(hitPhi, hitHeight)) {
          continue;
        }
      } else if (!heightInRange(hitHeight) || !phiInSweep(hitPhi)) {
        continue;
      }
      const double radialDistance = std::hypot(hitU, hitV);
      const Vec3 hitNormal = (mAxisU * (hitU / radialDistance) + mAxisV * (hitV / radialDistance)) * mNormalSign;
      hits.push_back({candidate, hitNormal});
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

  /// Exact divergence-theorem contribution over the (phi, h) rectangle; for a wire trim the same
  /// integrand (1/3) X . n |X_phi x X_h| (with |X_phi x X_h| = radius) is integrated numerically.
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      return integrateOverCurveTrim(mTrimOuter, mTrimInner, [this](double phi, double height) {
        const Vec3 surfacePoint = pointAt(phi, height);
        return dot(surfacePoint, normalAt(surfacePoint)) * mRadius / 3.;
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
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage)) {
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
  bool pointInTrim(double phi, double theta) const
  {
    const double uCoord = unwrapAngleInto(phi, mPhiStart, mPhiStart + mPhiSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, theta});
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

  bool directionInTrim(const Vec3& localPoint) const
  {
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
      return pointInTrim(std::atan2(localPoint.yCoord, localPoint.xCoord), theta);
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
      if (!directionInTrim(localHit)) {
        continue;
      }
      hits.push_back({candidate, (mAxisU * localHit.xCoord + mAxisV * localHit.yCoord + mAxisW * localHit.zCoord) *
                                   (mNormalSign / mRadius)});
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

  /// Exact divergence-theorem contribution over the (theta, phi) rectangle; for a wire trim the
  /// same integrand (1/3) X . n |X_phi x X_theta| (with |X_phi x X_theta| = R^2 sin theta) is
  /// integrated numerically over the (phi, theta) domain.
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      return integrateOverCurveTrim(mTrimOuter, mTrimInner, [this](double phi, double theta) {
        const Vec3 surfacePoint = pointAt(theta, phi);
        return dot(surfacePoint, normalAt(surfacePoint)) * mRadius * mRadius * std::sin(theta) / 3.;
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
    if (!buildCurveTrim(outerTrim, innerTrims, mTrimOuter, mTrimInner, lower, upper, errorMessage)) {
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
  bool pointInTrim(double phi, double height) const
  {
    const double uCoord = unwrapAngleInto(phi, mPhiStart, mPhiStart + mPhiSweep);
    return curveTrimContains(mTrimOuter, mTrimInner, {uCoord, height});
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
      if (mHasWireTrim) {
        if (!pointInTrim(hitPhi, hitHeight)) {
          continue;
        }
      } else if (!heightInRange(hitHeight) || !phiInSweep(hitPhi)) {
        continue;
      }
      const double normalScale = mNormalSign / std::sqrt(1. + mSlope * mSlope);
      const Vec3 hitNormal =
        (mAxisU * (hitU / radialDistance) + mAxisV * (hitV / radialDistance) - mAxisW * mSlope) * normalScale;
      hits.push_back({candidate, hitNormal});
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

  /// Exact divergence-theorem contribution over the (phi, h) rectangle; for a wire trim the same
  /// integrand (1/3) X . n |X_phi x X_h| (with |X_phi x X_h| = r(h) sqrt(1 + slope^2)) is
  /// integrated numerically.
  double capacityContribution() const override
  {
    if (mHasWireTrim) {
      const double slopeFactor = std::sqrt(1. + mSlope * mSlope);
      return integrateOverCurveTrim(mTrimOuter, mTrimInner, [this, slopeFactor](double phi, double height) {
        const Vec3 surfacePoint = pointAt(phi, height);
        return dot(surfacePoint, normalAt(surfacePoint)) * radiusAt(height) * slopeFactor / 3.;
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
};

/// Validate closure and orientation of a set of bounded surfaces using a directed half-edge map.
inline ClosureReport validateClosure(const std::vector<std::unique_ptr<BoundedSurface>>& surfaces)
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

  report.closed = (report.boundaryEdges == 0) && (report.nonManifoldEdges == 0);
  report.orientationConsistent = (report.reversedEdges == 0);
  return report;
}

} // namespace o2::base::surface

#endif
