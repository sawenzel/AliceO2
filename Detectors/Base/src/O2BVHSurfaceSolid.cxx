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

#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

using namespace o2::base;
ClassImp(O2BVHSurfaceSolid);

namespace
{
constexpr double kTolerance = 1.e-9;
constexpr double kToleranceSq = kTolerance * kTolerance;
constexpr double kAreaTolerance = 1.e-18;
constexpr double kRayTolerance = 1.e-9;
constexpr double kIntersectionTolerance = 1.e-7;

struct Vec2 {
  double uCoord = 0.;
  double vCoord = 0.;
};

struct Vec3 {
  double xCoord = 0.;
  double yCoord = 0.;
  double zCoord = 0.;
};

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

Vec3 operator+(const Vec3& firstVector, const Vec3& secondVector)
{
  return {firstVector.xCoord + secondVector.xCoord, firstVector.yCoord + secondVector.yCoord,
          firstVector.zCoord + secondVector.zCoord};
}

Vec3 operator-(const Vec3& firstVector, const Vec3& secondVector)
{
  return {firstVector.xCoord - secondVector.xCoord, firstVector.yCoord - secondVector.yCoord,
          firstVector.zCoord - secondVector.zCoord};
}

Vec3 operator*(const Vec3& vector, double scale)
{
  return {vector.xCoord * scale, vector.yCoord * scale, vector.zCoord * scale};
}

Vec3 operator*(double scale, const Vec3& vector)
{
  return vector * scale;
}

double dot(const Vec3& firstVector, const Vec3& secondVector)
{
  return firstVector.xCoord * secondVector.xCoord + firstVector.yCoord * secondVector.yCoord +
         firstVector.zCoord * secondVector.zCoord;
}

Vec3 cross(const Vec3& firstVector, const Vec3& secondVector)
{
  return {firstVector.yCoord * secondVector.zCoord - firstVector.zCoord * secondVector.yCoord,
          firstVector.zCoord * secondVector.xCoord - firstVector.xCoord * secondVector.zCoord,
          firstVector.xCoord * secondVector.yCoord - firstVector.yCoord * secondVector.xCoord};
}

double normSq(const Vec3& vector)
{
  return dot(vector, vector);
}

double norm(const Vec3& vector)
{
  return std::sqrt(normSq(vector));
}

Vec3 normalized(const Vec3& vector)
{
  const double vectorNorm = norm(vector);
  if (vectorNorm <= kTolerance) {
    return {};
  }
  return vector * (1. / vectorNorm);
}

double component(const Vec3& vector, int dimension)
{
  if (dimension == 0) {
    return vector.xCoord;
  }
  if (dimension == 1) {
    return vector.yCoord;
  }
  return vector.zCoord;
}

void assignComponent(Vec3& vector, int dimension, double value)
{
  if (dimension == 0) {
    vector.xCoord = value;
  } else if (dimension == 1) {
    vector.yCoord = value;
  } else {
    vector.zCoord = value;
  }
}

bool finite(const Vec2& point)
{
  return std::isfinite(point.uCoord) && std::isfinite(point.vCoord);
}

bool finite(const Vec3& point)
{
  return std::isfinite(point.xCoord) && std::isfinite(point.yCoord) && std::isfinite(point.zCoord);
}

double distanceSq(const Vec2& firstPoint, const Vec2& secondPoint)
{
  const double deltaU = firstPoint.uCoord - secondPoint.uCoord;
  const double deltaV = firstPoint.vCoord - secondPoint.vCoord;
  return deltaU * deltaU + deltaV * deltaV;
}

double distanceSq(const Vec3& firstPoint, const Vec3& secondPoint)
{
  return normSq(firstPoint - secondPoint);
}

double cross2D(const Vec2& firstVector, const Vec2& secondVector)
{
  return firstVector.uCoord * secondVector.vCoord - firstVector.vCoord * secondVector.uCoord;
}

Vec2 operator-(const Vec2& firstPoint, const Vec2& secondPoint)
{
  return {firstPoint.uCoord - secondPoint.uCoord, firstPoint.vCoord - secondPoint.vCoord};
}

double pointSegmentDistanceSq(const Vec2& point, const Vec2& segmentStart, const Vec2& segmentEnd)
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

double pointSegmentDistanceSq(const Vec3& point, const Vec3& segmentStart, const Vec3& segmentEnd)
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

enum class WireClassification { Outside, Boundary, Inside };

struct SurfaceWire {
  std::vector<Vec2> vertices;

  bool initialize(const std::vector<Vec2>& inputVertices, std::string& errorMessage)
  {
    vertices.clear();
    vertices.reserve(inputVertices.size());

    for (const auto& vertex : inputVertices) {
      if (!finite(vertex)) {
        errorMessage = "wire contains a non-finite vertex";
        return false;
      }
      if (vertices.empty() || distanceSq(vertices.back(), vertex) > kToleranceSq) {
        vertices.push_back(vertex);
      }
    }

    if (vertices.size() > 1 && distanceSq(vertices.front(), vertices.back()) <= kToleranceSq) {
      vertices.pop_back();
    }

    if (vertices.size() < 3) {
      errorMessage = "wire needs at least three distinct vertices";
      return false;
    }
    if (std::abs(signedArea()) <= kAreaTolerance) {
      errorMessage = "wire has zero area";
      return false;
    }
    return true;
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

bool pointInTriangle(const Vec2& point, const Vec2& firstVertex, const Vec2& secondVertex, const Vec2& thirdVertex)
{
  const double firstCross = cross2D(secondVertex - firstVertex, point - firstVertex);
  const double secondCross = cross2D(thirdVertex - secondVertex, point - secondVertex);
  const double thirdCross = cross2D(firstVertex - thirdVertex, point - thirdVertex);
  return firstCross >= -kTolerance && secondCross >= -kTolerance && thirdCross >= -kTolerance;
}

std::vector<std::array<int, 3>> triangulateSimpleWire(const SurfaceWire& wire)
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

struct PlanarBoundedSurface {
  Vec3 origin;
  Vec3 axisU;
  Vec3 axisV;
  Vec3 normal;
  double metricUU = 0.;
  double metricUV = 0.;
  double metricVV = 0.;
  double inverseMetricDet = 0.;
  double areaScale = 0.;
  SurfaceWire outerWire;
  std::vector<SurfaceWire> innerWires;

  bool initialize(const Vec3& surfaceOrigin, const Vec3& surfaceAxisU, const Vec3& surfaceAxisV,
                  const std::vector<Vec2>& outerWireVertices, const std::vector<std::vector<Vec2>>& innerWireVertices,
                  std::string& errorMessage)
  {
    if (!finite(surfaceOrigin) || !finite(surfaceAxisU) || !finite(surfaceAxisV)) {
      errorMessage = "surface frame contains a non-finite value";
      return false;
    }

    origin = surfaceOrigin;
    axisU = surfaceAxisU;
    axisV = surfaceAxisV;
    const Vec3 normalVector = cross(axisU, axisV);
    areaScale = norm(normalVector);
    if (areaScale <= kTolerance) {
      errorMessage = "surface frame axes are degenerate";
      return false;
    }
    normal = normalVector * (1. / areaScale);

    metricUU = dot(axisU, axisU);
    metricUV = dot(axisU, axisV);
    metricVV = dot(axisV, axisV);
    const double metricDet = metricUU * metricVV - metricUV * metricUV;
    if (std::abs(metricDet) <= kToleranceSq) {
      errorMessage = "surface frame metric is singular";
      return false;
    }
    inverseMetricDet = 1. / metricDet;

    if (!outerWire.initialize(outerWireVertices, errorMessage)) {
      errorMessage = "outer " + errorMessage;
      return false;
    }

    innerWires.clear();
    innerWires.reserve(innerWireVertices.size());
    for (const auto& innerWireInput : innerWireVertices) {
      SurfaceWire innerWire;
      if (!innerWire.initialize(innerWireInput, errorMessage)) {
        errorMessage = "inner " + errorMessage;
        return false;
      }
      innerWires.emplace_back(std::move(innerWire));
    }
    return true;
  }

  Vec3 toGlobal(const Vec2& point) const
  {
    return origin + axisU * point.uCoord + axisV * point.vCoord;
  }

  Vec2 toLocal(const Vec3& point) const
  {
    const Vec3 relativePoint = point - origin;
    const double projectionU = dot(relativePoint, axisU);
    const double projectionV = dot(relativePoint, axisV);
    return {(projectionU * metricVV - projectionV * metricUV) * inverseMetricDet,
            (projectionV * metricUU - projectionU * metricUV) * inverseMetricDet};
  }

  double planeDistance(const Vec3& point) const
  {
    return dot(point - origin, normal);
  }

  bool containsLocal(const Vec2& point, bool* boundary = nullptr) const
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

  bool containsPointOnSurface(const Vec3& point) const
  {
    if (std::abs(planeDistance(point)) > kTolerance) {
      return false;
    }
    return containsLocal(toLocal(point));
  }

  bool intersectRay(const Vec3& rayOrigin, const Vec3& rayDirection, double minDistance, double maxDistance,
                    double& distance) const
  {
    const double denominator = dot(normal, rayDirection);
    if (std::abs(denominator) <= kTolerance) {
      return false;
    }
    const double candidateDistance = dot(origin - rayOrigin, normal) / denominator;
    if (candidateDistance < minDistance || candidateDistance > maxDistance) {
      return false;
    }
    const Vec3 candidatePoint = rayOrigin + rayDirection * candidateDistance;
    if (!containsLocal(toLocal(candidatePoint))) {
      return false;
    }
    distance = candidateDistance;
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

  double distanceSqToPatch(const Vec3& point) const
  {
    const Vec2 projectedPoint = toLocal(point);
    if (containsLocal(projectedPoint)) {
      const double signedPlaneDistance = planeDistance(point);
      return signedPlaneDistance * signedPlaneDistance;
    }

    double bestDistanceSq = distanceSqToEdges(point, outerWire);
    for (const auto& innerWire : innerWires) {
      bestDistanceSq = std::min(bestDistanceSq, distanceSqToEdges(point, innerWire));
    }
    return bestDistanceSq;
  }

  void extendBounds(Vec3& lowerCorner, Vec3& upperCorner) const
  {
    auto extendPoint = [&](const Vec2& surfacePoint) {
      const Vec3 globalPoint = toGlobal(surfacePoint);
      lowerCorner.xCoord = std::min(lowerCorner.xCoord, globalPoint.xCoord);
      lowerCorner.yCoord = std::min(lowerCorner.yCoord, globalPoint.yCoord);
      lowerCorner.zCoord = std::min(lowerCorner.zCoord, globalPoint.zCoord);
      upperCorner.xCoord = std::max(upperCorner.xCoord, globalPoint.xCoord);
      upperCorner.yCoord = std::max(upperCorner.yCoord, globalPoint.yCoord);
      upperCorner.zCoord = std::max(upperCorner.zCoord, globalPoint.zCoord);
    };

    for (const auto& vertex : outerWire.vertices) {
      extendPoint(vertex);
    }
    for (const auto& innerWire : innerWires) {
      for (const auto& vertex : innerWire.vertices) {
        extendPoint(vertex);
      }
    }
  }

  double area() const
  {
    double parametricArea = std::abs(outerWire.signedArea());
    for (const auto& innerWire : innerWires) {
      parametricArea -= std::abs(innerWire.signedArea());
    }
    return std::max(0., parametricArea) * areaScale;
  }

  double capacityContribution() const
  {
    return dot(origin, normal) * area() / 3.;
  }

  void appendDisplayMesh(std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& triangles) const
  {
    const int firstVertexIndex = static_cast<int>(vertices.size());
    for (const auto& vertex : outerWire.vertices) {
      vertices.push_back(toGlobal(vertex));
    }

    const auto localTriangles = triangulateSimpleWire(outerWire);
    for (const auto& triangle : localTriangles) {
      triangles.push_back({firstVertexIndex + triangle[0], firstVertexIndex + triangle[1], firstVertexIndex + triangle[2]});
    }
  }
};

bool sameIntersection(double firstDistance, double secondDistance)
{
  return std::abs(firstDistance - secondDistance) <=
         kIntersectionTolerance * std::max(1., std::max(std::abs(firstDistance), std::abs(secondDistance)));
}

} // namespace

struct O2BVHSurfaceSolid::Impl {
  std::vector<PlanarBoundedSurface> surfaces;
  std::vector<Vec3> displayVertices;
  std::vector<std::array<int, 3>> displayTriangles;
  bool defined = false;
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

  PlanarBoundedSurface surface;
  std::string errorMessage;
  if (!surface.initialize(makeVec3(origin), makeVec3(axisU), makeVec3(axisV), convertedOuterWire, convertedInnerWires,
                          errorMessage)) {
    Error("AddPlanarSurface", "%s", errorMessage.c_str());
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

  fImpl->displayVertices.clear();
  fImpl->displayTriangles.clear();
  for (const auto& surface : fImpl->surfaces) {
    surface.appendDisplayMesh(fImpl->displayVertices, fImpl->displayTriangles);
  }

  fImpl->defined = true;
}

int O2BVHSurfaceSolid::GetNsurfaces() const
{
  return fImpl == nullptr ? 0 : static_cast<int>(fImpl->surfaces.size());
}

bool O2BVHSurfaceSolid::IsDefined() const
{
  return fImpl != nullptr && fImpl->defined;
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
    surface.extendBounds(lowerCorner, upperCorner);
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

  for (const auto& surface : fImpl->surfaces) {
    if (surface.containsPointOnSurface(testPoint)) {
      return true;
    }
  }

  const Vec3 testDirection = normalized({1., 1.41421356237, 1.73205080757});
  std::vector<double> intersections;
  intersections.reserve(fImpl->surfaces.size());
  for (const auto& surface : fImpl->surfaces) {
    double intersectionDistance = 0.;
    if (surface.intersectRay(testPoint, testDirection, kRayTolerance, TGeoShape::Big(), intersectionDistance)) {
      intersections.push_back(intersectionDistance);
    }
  }

  std::sort(intersections.begin(), intersections.end());
  int crossings = 0;
  double previousIntersection = 0.;
  bool hasPreviousIntersection = false;
  for (const double intersection : intersections) {
    if (!hasPreviousIntersection || !sameIntersection(intersection, previousIntersection)) {
      ++crossings;
      previousIntersection = intersection;
      hasPreviousIntersection = true;
    }
  }
  return (crossings & 1) != 0;
}

Double_t O2BVHSurfaceSolid::DistFromOutside(const Double_t* point, const Double_t* dir, Int_t, Double_t stepmax,
                                           Double_t* safe) const
{
  if (safe != nullptr) {
    *safe = Safety(point, kFALSE);
  }
  if (fImpl == nullptr) {
    return TGeoShape::Big();
  }

  const Vec3 rayOrigin = makeVec3(point);
  const Vec3 rayDirection = makeVec3(dir);
  double bestDistance = TGeoShape::Big();
  const double maxDistance = std::min(stepmax, bestDistance);
  for (const auto& surface : fImpl->surfaces) {
    if (dot(surface.normal, rayDirection) >= -kTolerance) {
      continue;
    }
    double intersectionDistance = 0.;
    if (surface.intersectRay(rayOrigin, rayDirection, kRayTolerance, maxDistance, intersectionDistance)) {
      bestDistance = std::min(bestDistance, intersectionDistance);
    }
  }
  return bestDistance;
}

Double_t O2BVHSurfaceSolid::DistFromInside(const Double_t* point, const Double_t* dir, Int_t, Double_t stepmax,
                                          Double_t* safe) const
{
  if (safe != nullptr) {
    *safe = Safety(point, kTRUE);
  }
  if (fImpl == nullptr) {
    return TGeoShape::Big();
  }

  const Vec3 rayOrigin = makeVec3(point);
  const Vec3 rayDirection = makeVec3(dir);
  double bestDistance = TGeoShape::Big();
  const double maxDistance = std::min(stepmax, bestDistance);
  for (const auto& surface : fImpl->surfaces) {
    if (dot(surface.normal, rayDirection) <= kTolerance) {
      continue;
    }
    double intersectionDistance = 0.;
    if (surface.intersectRay(rayOrigin, rayDirection, kRayTolerance, maxDistance, intersectionDistance)) {
      bestDistance = std::min(bestDistance, intersectionDistance);
    }
  }
  return bestDistance;
}

Double_t O2BVHSurfaceSolid::Safety(const Double_t* point, Bool_t) const
{
  if (fImpl == nullptr || fImpl->surfaces.empty()) {
    return TGeoShape::Big();
  }

  const Vec3 testPoint = makeVec3(point);
  double bestDistanceSq = std::numeric_limits<double>::infinity();
  for (const auto& surface : fImpl->surfaces) {
    bestDistanceSq = std::min(bestDistanceSq, surface.distanceSqToPatch(testPoint));
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
  const PlanarBoundedSurface* closestSurface = nullptr;
  double bestDistanceSq = std::numeric_limits<double>::infinity();
  for (const auto& surface : fImpl->surfaces) {
    const double surfaceDistanceSq = surface.distanceSqToPatch(testPoint);
    if (surfaceDistanceSq < bestDistanceSq) {
      bestDistanceSq = surfaceDistanceSq;
      closestSurface = &surface;
    }
  }

  if (closestSurface == nullptr) {
    norm[0] = 1.;
    norm[1] = 0.;
    norm[2] = 0.;
    return;
  }

  Vec3 normal = closestSurface->normal;
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
    capacity += surface.capacityContribution();
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