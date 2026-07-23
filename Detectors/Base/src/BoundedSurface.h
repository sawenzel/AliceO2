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

/// Kind of a 2D trimmed boundary curve.
enum class CurveKind { Line, ///< straight line segment
                       Arc   ///< circular arc
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

  bool isArc() const { return kind == CurveKind::Arc; }

  double sweep() const { return endAngle - startAngle; }

  /// Basic structural validity (finite data, positive radius for arcs).
  bool valid() const
  {
    if (kind == CurveKind::Line) {
      return finite(lineStart) && finite(lineEnd);
    }
    return finite(center) && std::isfinite(radius) && radius > kTolerance && std::isfinite(startAngle) &&
           std::isfinite(endAngle);
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
    return pointAtAngle(startAngle + parameter * sweep());
  }

  Vec2 startPoint() const { return kind == CurveKind::Line ? lineStart : pointAtAngle(startAngle); }
  Vec2 endPoint() const { return kind == CurveKind::Line ? lineEnd : pointAtAngle(endAngle); }

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
      if (surface::distanceSq(currentEnd, nextStart) > kToleranceSq) {
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
      if (curve.kind == CurveKind::Line) {
        std::swap(curve.lineStart, curve.lineEnd);
      } else {
        std::swap(curve.startAngle, curve.endAngle);
      }
    }
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
    const int arcSteps = std::max(1, segmentsPerArc);
    for (const auto& curve : curves) {
      if (curve.kind == CurveKind::Line) {
        samples.push_back(curve.startPoint());
      } else {
        for (int step = 0; step < arcSteps; ++step) {
          samples.push_back(curve.pointAt(static_cast<double>(step) / arcSteps));
        }
      }
    }
    samples.push_back(samples.front());
    return samples;
  }
};

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

  bool capacityIsExact() const override { return true; }

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
    return heightInRange(localPoint.zCoord) &&
           (radialDistance <= kTolerance || phiInSweep(std::atan2(localPoint.yCoord, localPoint.xCoord)));
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
      if (!heightInRange(hitHeight) || !phiInSweep(std::atan2(hitV, hitU))) {
        continue;
      }
      const double radialDistance = std::hypot(hitU, hitV);
      const Vec3 hitNormal = (mAxisU * (hitU / radialDistance) + mAxisV * (hitV / radialDistance)) * mNormalSign;
      hits.push_back({candidate, hitNormal});
    }
  }

  /// Exact distance to the trimmed patch: the (rho, h) and phi contributions separate, so the
  /// in-sweep case reduces to a 2D point/segment distance in the (rho, h) half-plane and the
  /// off-sweep case to the 3D distance to the nearest straight seam edge.
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

  /// Exact divergence-theorem contribution: (1/3) Int (X . n) dA over the (phi, h) rectangle.
  double capacityContribution() const override
  {
    const double endPhi = mPhiStart + mPhiSweep;
    const double phiFactor = dot(mCenter, mAxisU) * (std::sin(endPhi) - std::sin(mPhiStart)) -
                             dot(mCenter, mAxisV) * (std::cos(endPhi) - std::cos(mPhiStart));
    const double height = mHeightMax - mHeightMin;
    return mNormalSign * mRadius * height * (phiFactor + mRadius * mPhiSweep) / 3.;
  }

  bool capacityIsExact() const override { return true; }

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
    if (theta < mThetaMin - thetaTolerance || theta > mThetaMax + thetaTolerance) {
      return false;
    }
    const double transverseDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
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

  /// Exact divergence-theorem contribution over the (theta, phi) rectangle.
  double capacityContribution() const override
  {
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

  bool capacityIsExact() const override { return true; }

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
    if (!heightInRange(localPoint.zCoord)) {
      return false;
    }
    const double surfaceRadius = radiusAt(localPoint.zCoord);
    const double radialDistance = std::hypot(localPoint.xCoord, localPoint.yCoord);
    // |rho - r(h)| overestimates the true surface distance by sqrt(1 + slope^2)
    if (std::abs(radialDistance - surfaceRadius) > kTolerance * std::sqrt(1. + mSlope * mSlope)) {
      return false;
    }
    return radialDistance <= kTolerance || phiInSweep(std::atan2(localPoint.yCoord, localPoint.xCoord));
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
      if (!heightInRange(hitHeight)) {
        continue;
      }
      const double hitU = localOrigin.xCoord + candidate * localDirection.xCoord;
      const double hitV = localOrigin.yCoord + candidate * localDirection.yCoord;
      const double radialDistance = std::hypot(hitU, hitV);
      if (radialDistance <= kTolerance) {
        continue; // apex hit: the normal is undefined there
      }
      if (!phiInSweep(std::atan2(hitV, hitU))) {
        continue;
      }
      const double normalScale = mNormalSign / std::sqrt(1. + mSlope * mSlope);
      const Vec3 hitNormal =
        (mAxisU * (hitU / radialDistance) + mAxisV * (hitV / radialDistance) - mAxisW * mSlope) * normalScale;
      hits.push_back({candidate, hitNormal});
    }
  }

  /// Exact distance to the trimmed patch: for an in-sweep azimuth the problem reduces to the
  /// 2D distance to the straight generator segment in the (rho, h) half-plane; off-sweep the
  /// nearest point lies on one of the two straight seam generators.
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

  /// Exact divergence-theorem contribution: (1/3) Int (X . n) dA over the (phi, h) rectangle;
  /// the area element and the normal's sqrt(1 + slope^2) factors cancel.
  double capacityContribution() const override
  {
    const double endPhi = mPhiStart + mPhiSweep;
    const double phiFactor = dot(mCenter, mAxisU) * (std::sin(endPhi) - std::sin(mPhiStart)) -
                             dot(mCenter, mAxisV) * (std::cos(endPhi) - std::cos(mPhiStart));
    const double radiusIntegral = mRadius0 * (mHeightMax - mHeightMin) +
                                  0.5 * mSlope * (mHeightMax * mHeightMax - mHeightMin * mHeightMin);
    return mNormalSign * radiusIntegral *
           (phiFactor + (mRadius0 - mSlope * dot(mCenter, mAxisW)) * mPhiSweep) / 3.;
  }

  bool capacityIsExact() const override { return true; }

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
