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

  /// Intersect a ray with the trimmed patch. On success returns the ray parameter \a distance
  /// (within [minDistance, maxDistance]) and the outward-oriented surface normal at the hit.
  virtual bool intersectRay(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance,
                            double maxDistance, double& distance, Vec3& hitNormal) const = 0;

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

  bool intersectRay(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance, double maxDistance,
                    double& distance, Vec3& hitNormal) const override
  {
    const double denominator = dot(mNormal, rayDirection);
    if (std::abs(denominator) <= kTolerance) {
      return false;
    }
    const double candidateDistance = dot(mOrigin - rayOrigin, mNormal) / denominator;
    if (candidateDistance < minDistance || candidateDistance > maxDistance) {
      return false;
    }
    const Vec3 candidatePoint = rayOrigin + rayDirection * candidateDistance;
    if (!containsLocal(toLocal(candidatePoint))) {
      return false;
    }
    distance = candidateDistance;
    hitNormal = mNormal;
    return true;
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

  bool intersectRay(const Vec3&, const Vec3&, double, double, double&, Vec3&) const override { return false; }

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
