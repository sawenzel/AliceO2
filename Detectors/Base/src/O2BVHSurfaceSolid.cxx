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
} // namespace

struct O2BVHSurfaceSolid::Impl {
  std::vector<std::unique_ptr<BoundedSurface>> surfaces;
  std::vector<Vec3> displayVertices;
  std::vector<std::array<int, 3>> displayTriangles;
  ClosureReport closure;
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
    surface->appendDisplayMesh(fImpl->displayVertices, fImpl->displayTriangles);
  }

  fImpl->closure = validateClosure(fImpl->surfaces);
  if (check && !fImpl->surfaces.empty()) {
    const auto& closure = fImpl->closure;
    if (!closure.closed) {
      Warning("CloseShape", "Shape %s is not a closed manifold: %d boundary edge(s), %d non-manifold edge(s)",
              GetName(), closure.boundaryEdges, closure.nonManifoldEdges);
    }
    if (!closure.orientationConsistent) {
      Warning("CloseShape", "Shape %s has %d inconsistently oriented (reversed) boundary edge(s)", GetName(),
              closure.reversedEdges);
    }
    if (closure.closed && closure.signedVolume < 0.) {
      Warning("CloseShape",
              "Shape %s has inward-pointing surface normals (signed volume %g); navigation expects outward normals",
              GetName(), closure.signedVolume);
    }
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

bool O2BVHSurfaceSolid::IsClosed() const
{
  return fImpl != nullptr && fImpl->defined && fImpl->closure.closed;
}

bool O2BVHSurfaceSolid::IsOrientationConsistent() const
{
  return fImpl != nullptr && fImpl->defined && fImpl->closure.orientationConsistent;
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
    if (surface->containsPointOnSurface(testPoint)) {
      return true;
    }
  }

  const Vec3 testDirection = normalized({1., 1.41421356237, 1.73205080757});
  std::vector<double> intersections;
  intersections.reserve(fImpl->surfaces.size());
  for (const auto& surface : fImpl->surfaces) {
    double intersectionDistance = 0.;
    Vec3 hitNormal;
    if (surface->intersectRay(testPoint, testDirection, kRayTolerance, TGeoShape::Big(), intersectionDistance,
                              hitNormal)) {
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
    double intersectionDistance = 0.;
    Vec3 hitNormal;
    if (surface->intersectRay(rayOrigin, rayDirection, kRayTolerance, maxDistance, intersectionDistance, hitNormal)) {
      // entering: the outward normal opposes the ray direction
      if (dot(hitNormal, rayDirection) >= -kTolerance) {
        continue;
      }
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
    double intersectionDistance = 0.;
    Vec3 hitNormal;
    if (surface->intersectRay(rayOrigin, rayDirection, kRayTolerance, maxDistance, intersectionDistance, hitNormal)) {
      // exiting: the outward normal is aligned with the ray direction
      if (dot(hitNormal, rayDirection) <= kTolerance) {
        continue;
      }
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