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

#define BOOST_TEST_MODULE Test O2BVHSurfaceSolid class
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "DetectorsBase/O2BVHSurfaceSolid.h"

#include "../src/BoundedSurface.h"

#include "TGeoBBox.h"
#include "TGeoShape.h"

#include <array>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace
{
using SurfaceSolid = o2::base::O2BVHSurfaceSolid;
using Point2D = SurfaceSolid::Point2D;
using Point3D = SurfaceSolid::Point3D;
namespace surf = o2::base::surface;

std::vector<Point2D> rectangleWire(double extentU, double extentV)
{
  return {{0., 0.}, {extentU, 0.}, {extentU, extentV}, {0., extentV}};
}

// Add a single box face by index (0:+x 1:-x 2:+y 3:-y 4:+z 5:-z). When "reversed" is set the
// face's parametric axes are swapped, which flips the outward normal inward without changing the
// covered rectangle - used to build an orientation-inconsistent fixture.
bool addBoxFace(SurfaceSolid& solid, int faceIndex, double halfX, double halfY, double halfZ, bool reversed = false)
{
  Point3D origin{};
  Point3D axisU{};
  Point3D axisV{};
  double extentU = 0.;
  double extentV = 0.;
  switch (faceIndex) {
    case 0:
      origin = {halfX, -halfY, -halfZ};
      axisU = {0., 1., 0.};
      axisV = {0., 0., 1.};
      extentU = 2. * halfY;
      extentV = 2. * halfZ;
      break;
    case 1:
      origin = {-halfX, -halfY, -halfZ};
      axisU = {0., 0., 1.};
      axisV = {0., 1., 0.};
      extentU = 2. * halfZ;
      extentV = 2. * halfY;
      break;
    case 2:
      origin = {-halfX, halfY, -halfZ};
      axisU = {0., 0., 1.};
      axisV = {1., 0., 0.};
      extentU = 2. * halfZ;
      extentV = 2. * halfX;
      break;
    case 3:
      origin = {-halfX, -halfY, -halfZ};
      axisU = {1., 0., 0.};
      axisV = {0., 0., 1.};
      extentU = 2. * halfX;
      extentV = 2. * halfZ;
      break;
    case 4:
      origin = {-halfX, -halfY, halfZ};
      axisU = {1., 0., 0.};
      axisV = {0., 1., 0.};
      extentU = 2. * halfX;
      extentV = 2. * halfY;
      break;
    default:
      origin = {-halfX, -halfY, -halfZ};
      axisU = {0., 1., 0.};
      axisV = {1., 0., 0.};
      extentU = 2. * halfY;
      extentV = 2. * halfX;
      break;
  }
  if (reversed) {
    std::swap(axisU, axisV);
    std::swap(extentU, extentV);
  }
  return solid.AddPlanarSurface(origin, axisU, axisV, rectangleWire(extentU, extentV));
}

void addBoxSurfaces(SurfaceSolid& solid, double halfX, double halfY, double halfZ)
{
  for (int faceIndex = 0; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(solid, faceIndex, halfX, halfY, halfZ));
  }
}

surf::SurfaceWire makeWire(const std::vector<surf::Vec2>& vertices, surf::WireRole role, surf::WireStatus& status)
{
  surf::SurfaceWire wire;
  wire.initialize(vertices, role, status);
  return wire;
}

void checkClose(double value, double reference, double tolerance = 1.e-9)
{
  BOOST_CHECK_SMALL(value - reference, tolerance);
}

} // namespace

BOOST_AUTO_TEST_CASE(PlanarBoxNavigationMatchesTGeoBBox)
{
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;

  SurfaceSolid solid("planarBox");
  addBoxSurfaces(solid, halfX, halfY, halfZ);
  solid.CloseShape();

  TGeoBBox reference("referenceBox", halfX, halfY, halfZ);

  BOOST_CHECK(solid.IsDefined());
  BOOST_CHECK_EQUAL(solid.GetNsurfaces(), 6);

  int meshVertices = 0;
  int meshSegments = 0;
  int meshPolygons = 0;
  solid.GetMeshNumbers(meshVertices, meshSegments, meshPolygons);
  BOOST_CHECK_EQUAL(meshVertices, 24);
  BOOST_CHECK_EQUAL(meshSegments, 36);
  BOOST_CHECK_EQUAL(meshPolygons, 12);

  const std::array<std::array<double, 3>, 5> insidePoints{{{0., 0., 0.}, {0.9, 0., 0.}, {1., 0., 0.},
                                                          {1., 2., 3.}, {-1., -2., -3.}}};
  for (const auto& point : insidePoints) {
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      BOOST_CHECK(solid.Contains(point.data()));
      BOOST_CHECK(reference.Contains(point.data()));
    }
  }

  const std::array<std::array<double, 3>, 4> outsidePoints{{{1.1, 0., 0.}, {0., 2.1, 0.}, {0., 0., -3.1},
                                                           {2., 3., 4.}}};
  for (const auto& point : outsidePoints) {
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      BOOST_CHECK(!solid.Contains(point.data()));
      BOOST_CHECK(!reference.Contains(point.data()));
    }
  }

  const double fromLeft[3] = {-3., 0., 0.};
  const double toRight[3] = {1., 0., 0.};
  checkClose(solid.DistFromOutside(fromLeft, toRight, 3), reference.DistFromOutside(fromLeft, toRight, 3));

  const double fromFront[3] = {0., -5., 0.};
  const double toBack[3] = {0., 1., 0.};
  checkClose(solid.DistFromOutside(fromFront, toBack, 3), reference.DistFromOutside(fromFront, toBack, 3));

  const double fromCenter[3] = {0., 0., 0.};
  const double alongX[3] = {1., 0., 0.};
  const double alongZ[3] = {0., 0., -1.};
  checkClose(solid.DistFromInside(fromCenter, alongX, 3), reference.DistFromInside(fromCenter, alongX, 3));
  checkClose(solid.DistFromInside(fromCenter, alongZ, 3), reference.DistFromInside(fromCenter, alongZ, 3));

  const double outsideSafetyPoint[3] = {2.5, 0., 0.};
  checkClose(solid.Safety(fromCenter, kTRUE), reference.Safety(fromCenter, kTRUE));
  checkClose(solid.Safety(outsideSafetyPoint, kFALSE), reference.Safety(outsideSafetyPoint, kFALSE));

  const double normalPoint[3] = {halfX, 0., 0.};
  double normal[3] = {0., 0., 0.};
  solid.ComputeNormal(normalPoint, alongX, normal);
  checkClose(normal[0], 1.);
  checkClose(normal[1], 0.);
  checkClose(normal[2], 0.);

  checkClose(solid.Capacity(), 8. * halfX * halfY * halfZ);

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());
}

BOOST_AUTO_TEST_CASE(WireValidationAndOrientation)
{
  using surf::WireRole;
  using surf::WireStatus;

  // outer square given clockwise (negative area) must be re-oriented to CCW
  WireStatus reversedStatus = WireStatus::Valid;
  auto reversedOuter = makeWire({{0., 0.}, {0., 1.}, {1., 1.}, {1., 0.}}, WireRole::Outer, reversedStatus);
  BOOST_CHECK(reversedStatus == WireStatus::Reversed);
  BOOST_CHECK_GT(reversedOuter.signedArea(), 0.);

  // outer square already CCW stays valid
  WireStatus outerStatus = WireStatus::Valid;
  auto outerWire = makeWire({{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}, WireRole::Outer, outerStatus);
  BOOST_CHECK(outerStatus == WireStatus::Valid);
  BOOST_CHECK_GT(outerWire.signedArea(), 0.);

  // inner (hole) wire must end up clockwise (negative area)
  WireStatus innerStatus = WireStatus::Valid;
  auto innerWire = makeWire({{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}, WireRole::Inner, innerStatus);
  BOOST_CHECK(innerStatus == WireStatus::Reversed);
  BOOST_CHECK_LT(innerWire.signedArea(), 0.);

  // degenerate / invalid inputs are rejected with a specific status
  surf::SurfaceWire scratch;
  WireStatus status = WireStatus::Valid;
  BOOST_CHECK(!scratch.initialize({{0., 0.}, {1., 0.}}, WireRole::Outer, status));
  BOOST_CHECK(status == WireStatus::TooFewVertices);

  BOOST_CHECK(!scratch.initialize({{0., 0.}, {1., 0.}, {2., 0.}}, WireRole::Outer, status));
  BOOST_CHECK(status == WireStatus::ZeroArea);

  // self-touching (pinched) loop: a non-adjacent vertex repeats
  BOOST_CHECK(!scratch.initialize({{0., 0.}, {1., 0.}, {0., 0.}, {1., 1.}}, WireRole::Outer, status));
  BOOST_CHECK(status == WireStatus::DegenerateVertex);

  // explicit edge list that does not close is flagged as open
  const std::vector<surf::SurfaceEdge> openEdges{{{0., 0.}, {1., 0.}}, {{1., 0.}, {1., 1.}}, {{1., 1.}, {0.5, 0.5}}};
  BOOST_CHECK(!scratch.initializeFromEdges(openEdges, WireRole::Outer, status));
  BOOST_CHECK(status == WireStatus::Open);

  // a closed edge list is accepted
  const std::vector<surf::SurfaceEdge> closedEdges{
    {{0., 0.}, {1., 0.}}, {{1., 0.}, {1., 1.}}, {{1., 1.}, {0., 1.}}, {{0., 1.}, {0., 0.}}};
  BOOST_CHECK(scratch.initializeFromEdges(closedEdges, WireRole::Outer, status));

  // point classification: inside, outside, and on-edge
  BOOST_CHECK(outerWire.classify({0.5, 0.5}) == surf::WireClassification::Inside);
  BOOST_CHECK(outerWire.classify({1.5, 0.5}) == surf::WireClassification::Outside);
  BOOST_CHECK(outerWire.classify({0.5, 0.}) == surf::WireClassification::Boundary);
}

BOOST_AUTO_TEST_CASE(DummyBoundedSurfaceInterface)
{
  auto dummy = std::make_unique<surf::DummyBoundedSurface>(surf::Vec3{0., 0., 0.}, surf::Vec3{1., 0., 0.},
                                                           surf::Vec3{0., 1., 0.});

  surf::Vec3 lower{surf::Vec3{1.e30, 1.e30, 1.e30}};
  surf::Vec3 upper{surf::Vec3{-1.e30, -1.e30, -1.e30}};
  dummy->conservativeBounds(lower, upper);
  checkClose(lower.xCoord, 0.);
  checkClose(upper.xCoord, 1.);
  checkClose(upper.yCoord, 1.);

  const surf::Vec3 normal = dummy->normalAt({0., 0., 0.});
  checkClose(std::abs(normal.zCoord), 1.);
  BOOST_CHECK(!dummy->capacityIsExact());

  std::vector<surf::Vec3> vertices;
  std::vector<std::array<int, 3>> triangles;
  dummy->appendDisplayMesh(vertices, triangles);
  BOOST_CHECK_EQUAL(vertices.size(), 3u);
  BOOST_CHECK_EQUAL(triangles.size(), 1u);

  // a single open triangle is not a closed manifold
  std::vector<std::unique_ptr<surf::BoundedSurface>> surfaces;
  surfaces.emplace_back(std::move(dummy));
  const surf::ClosureReport report = surf::validateClosure(surfaces);
  BOOST_CHECK(!report.closed);
  BOOST_CHECK_EQUAL(report.boundaryEdges, 3);
}

BOOST_AUTO_TEST_CASE(SolidClosureDetectsMissingAndReversedFaces)
{
  constexpr double halfX = 1.;
  constexpr double halfY = 1.5;
  constexpr double halfZ = 2.;

  // missing face: only five of the six box faces are added
  SurfaceSolid missing("missingFaceBox");
  for (int faceIndex = 1; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(missing, faceIndex, halfX, halfY, halfZ));
  }
  missing.CloseShape(false);
  BOOST_CHECK(!missing.IsClosed());

  // reversed face: the +x face keeps its geometry but has an inward normal
  SurfaceSolid reversed("reversedFaceBox");
  BOOST_REQUIRE(addBoxFace(reversed, 0, halfX, halfY, halfZ, true));
  for (int faceIndex = 1; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(reversed, faceIndex, halfX, halfY, halfZ));
  }
  reversed.CloseShape(false);
  BOOST_CHECK(reversed.IsClosed());
  BOOST_CHECK(!reversed.IsOrientationConsistent());
}