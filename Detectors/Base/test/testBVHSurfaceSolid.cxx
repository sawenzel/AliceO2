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

BOOST_AUTO_TEST_CASE(NumericalConventions)
{
  using surf::WireRole;
  using surf::WireStatus;
  using surf::WireClassification;

  // near-boundary point classification: a unit square wire, points offset from the bottom edge.
  WireStatus status = WireStatus::Valid;
  auto square = makeWire({{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}, WireRole::Outer, status);
  BOOST_REQUIRE(status == WireStatus::Valid);

  // within tolerance of an edge -> Boundary, on both sides
  BOOST_CHECK(square.classify({0.5, 0.5 * surf::kTolerance}) == WireClassification::Boundary);
  BOOST_CHECK(square.classify({0.5, -0.5 * surf::kTolerance}) == WireClassification::Boundary);
  // clearly beyond tolerance -> Inside / Outside
  BOOST_CHECK(square.classify({0.5, 1.e3 * surf::kTolerance}) == WireClassification::Inside);
  BOOST_CHECK(square.classify({0.5, -1.e3 * surf::kTolerance}) == WireClassification::Outside);

  // near-tangent rays against a planar surface in the z = 0 plane.
  surf::PlanarBoundedSurface plane;
  std::string planeError;
  BOOST_REQUIRE(plane.initialize({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.},
                                 {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}, {}, planeError));

  const surf::Vec3 origin{0.5, 0.5, 1.};
  double distance = 0.;
  surf::Vec3 hitNormal{};

  // direction almost parallel to the plane (tiny z component) -> grazing miss
  const surf::Vec3 grazing = surf::normalized({1., 0., 0.1 * surf::kTolerance});
  BOOST_CHECK(!plane.intersectRay(origin, grazing, 0., 1.e30, distance, hitNormal));

  // steeper direction -> real intersection at the plane
  const surf::Vec3 steep = surf::normalized({0., 0., -1.});
  BOOST_REQUIRE(plane.intersectRay(origin, steep, 0., 1.e30, distance, hitNormal));
  checkClose(distance, 1.);

  // a hit rejected when it falls below the minimum ray parameter
  BOOST_CHECK(!plane.intersectRay(origin, steep, 2., 1.e30, distance, hitNormal));

  // duplicate-intersection clustering respects kIntersectionTolerance.
  BOOST_CHECK(surf::sameIntersection(1., 1. + 0.1 * surf::kIntersectionTolerance));
  BOOST_CHECK(!surf::sameIntersection(1., 1. + 1.e3 * surf::kIntersectionTolerance));
}

BOOST_AUTO_TEST_CASE(WireDataModel)
{
  using surf::WireRole;
  using surf::WireStatus;
  using surf::WireClassification;

  // outer square wire: area, orientation, parametric AABB, and boundary sampling.
  WireStatus status = WireStatus::Valid;
  auto square = makeWire({{0., 0.}, {2., 0.}, {2., 3.}, {0., 3.}}, WireRole::Outer, status);
  BOOST_REQUIRE(status == WireStatus::Valid);
  checkClose(square.signedArea(), 6.);

  surf::Vec2 lower{1.e30, 1.e30};
  surf::Vec2 upper{-1.e30, -1.e30};
  square.parametricBounds(lower, upper);
  checkClose(lower.uCoord, 0.);
  checkClose(lower.vCoord, 0.);
  checkClose(upper.uCoord, 2.);
  checkClose(upper.vCoord, 3.);

  // the sampled boundary is the closed vertex ring (first vertex repeated at the end).
  const auto samples = square.sampledBoundary();
  BOOST_CHECK_EQUAL(samples.size(), square.vertices.size() + 1);
  checkClose(samples.front().uCoord, samples.back().uCoord);
  checkClose(samples.front().vCoord, samples.back().vCoord);

  // edge distance / projection (closest point) on the bottom edge.
  const surf::SurfaceEdge bottom{{0., 0.}, {2., 0.}};
  double parameter = -1.;
  const surf::Vec2 projected = bottom.closestPoint({1., 5.}, parameter);
  checkClose(projected.uCoord, 1.);
  checkClose(projected.vCoord, 0.);
  checkClose(parameter, 0.5);
  // projection is clamped to the segment endpoints.
  bottom.closestPoint({-5., 1.}, parameter);
  checkClose(parameter, 0.);
  bottom.closestPoint({5., 1.}, parameter);
  checkClose(parameter, 1.);
  checkClose(std::sqrt(bottom.distanceSq({1., 4.})), 4.);

  // reversed wire: same shape, opposite winding sign, identical parametric AABB.
  WireStatus reversedStatus = WireStatus::Valid;
  auto reversed = makeWire({{0., 0.}, {0., 3.}, {2., 3.}, {2., 0.}}, WireRole::Outer, reversedStatus);
  BOOST_CHECK(reversedStatus == WireStatus::Reversed);
  checkClose(reversed.signedArea(), 6.); // normalized back to CCW (positive)
  surf::Vec2 reversedLower{1.e30, 1.e30};
  surf::Vec2 reversedUpper{-1.e30, -1.e30};
  reversed.parametricBounds(reversedLower, reversedUpper);
  checkClose(reversedUpper.uCoord, 2.);
  checkClose(reversedUpper.vCoord, 3.);

  // open wire via an explicit non-closing edge list is rejected.
  surf::SurfaceWire scratch;
  const std::vector<surf::SurfaceEdge> openEdges{{{0., 0.}, {2., 0.}}, {{2., 0.}, {2., 3.}}, {{2., 3.}, {1., 1.}}};
  BOOST_CHECK(!scratch.initializeFromEdges(openEdges, WireRole::Outer, status));
  BOOST_CHECK(status == WireStatus::Open);

  // point-on-edge classification.
  BOOST_CHECK(square.classify({1., 0.}) == WireClassification::Boundary);
  BOOST_CHECK(square.classify({1., 1.5}) == WireClassification::Inside);
  BOOST_CHECK(square.classify({3., 1.5}) == WireClassification::Outside);

  // square-with-hole: a planar surface with one inner (hole) wire.
  surf::PlanarBoundedSurface holedFace;
  std::string faceError;
  const std::vector<surf::Vec2> outer{{0., 0.}, {4., 0.}, {4., 4.}, {0., 4.}};
  const std::vector<std::vector<surf::Vec2>> holes{{{1., 1.}, {3., 1.}, {3., 3.}, {1., 3.}}};
  BOOST_REQUIRE(holedFace.initialize({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, outer, holes, faceError));

  bool boundary = false;
  BOOST_CHECK(holedFace.containsLocal({0.5, 0.5}, &boundary)); // in material, outside the hole
  BOOST_CHECK(!boundary);
  BOOST_CHECK(!holedFace.containsLocal({2., 2.})); // inside the hole -> not on the patch
  BOOST_CHECK(holedFace.containsLocal({2., 1.}, &boundary)); // on the hole boundary -> on the patch
  BOOST_CHECK(boundary);

  // the trimmed area accounts for the hole (16 - 4).
  checkClose(holedFace.area(), 12.);
}

BOOST_AUTO_TEST_CASE(TrimmedCurveBoundaries)
{
  using surf::Curve2D;
  using surf::CurveWire;
  using surf::WireClassification;
  using surf::WireRole;
  using surf::WireStatus;

  // --- line curve: endpoint, tangent, bounds, projection -------------------------------------
  const Curve2D line = Curve2D::makeLine({0., 0.}, {4., 0.});
  checkClose(line.startPoint().uCoord, 0.);
  checkClose(line.endPoint().uCoord, 4.);
  const surf::Vec2 lineTangent = line.tangentAt(0.5);
  checkClose(lineTangent.uCoord, 1.);
  checkClose(lineTangent.vCoord, 0.);

  double lineParameter = -1.;
  const surf::Vec2 lineProjection = line.closestPoint({1., 5.}, lineParameter);
  checkClose(lineProjection.uCoord, 1.);
  checkClose(lineProjection.vCoord, 0.);
  checkClose(lineParameter, 0.25);
  checkClose(std::sqrt(line.distanceSq({1., 5.})), 5.);

  // --- arc curve: endpoint, tangent, exact bounds, projection --------------------------------
  // quarter circle of radius 2 centred at the origin, from angle 0 to pi/2.
  const Curve2D quarter = Curve2D::makeArc({0., 0.}, 2., 0., surf::kHalfPi);
  checkClose(quarter.startPoint().uCoord, 2.);
  checkClose(quarter.startPoint().vCoord, 0.);
  checkClose(quarter.endPoint().uCoord, 0.);
  checkClose(quarter.endPoint().vCoord, 2.);
  // tangent at the start of a CCW arc points in +v.
  const surf::Vec2 arcTangent = quarter.tangentAt(0.);
  checkClose(arcTangent.uCoord, 0.);
  checkClose(arcTangent.vCoord, 1.);

  // the quarter arc's exact bounding box is [0, 2] x [0, 2] (no cardinal extreme inside).
  surf::Vec2 arcLower{1.e30, 1.e30};
  surf::Vec2 arcUpper{-1.e30, -1.e30};
  quarter.extendBounds(arcLower, arcUpper);
  checkClose(arcLower.uCoord, 0.);
  checkClose(arcLower.vCoord, 0.);
  checkClose(arcUpper.uCoord, 2.);
  checkClose(arcUpper.vCoord, 2.);

  // projection of a far radial point lands on the circle (distance = |d - r|).
  double arcParameter = -1.;
  const surf::Vec2 arcProjection = quarter.closestPoint({5., 5.}, arcParameter);
  checkClose(std::hypot(arcProjection.uCoord, arcProjection.vCoord), 2.);
  checkClose(arcParameter, 0.5);

  // a full circle's exact bounding box spans the whole diameter in both axes.
  const Curve2D circle = Curve2D::makeCircle({1., -1.}, 3.);
  surf::Vec2 circleLower{1.e30, 1.e30};
  surf::Vec2 circleUpper{-1.e30, -1.e30};
  circle.extendBounds(circleLower, circleUpper);
  checkClose(circleLower.uCoord, -2.);
  checkClose(circleUpper.uCoord, 4.);
  checkClose(circleLower.vCoord, -4.);
  checkClose(circleUpper.vCoord, 2.);

  // --- disk: one full-circle outer wire ------------------------------------------------------
  WireStatus status = WireStatus::Valid;
  CurveWire disk;
  BOOST_REQUIRE(disk.initialize({Curve2D::makeCircle({0., 0.}, 2.)}, WireRole::Outer, status));
  BOOST_CHECK(status == WireStatus::Valid);
  // exact area of the disk is pi * r^2.
  checkClose(disk.signedArea(), surf::kPi * 4., 1.e-9);
  BOOST_CHECK(disk.classify({0., 0.}) == WireClassification::Inside);
  BOOST_CHECK(disk.classify({1.5, 0.}) == WireClassification::Inside);
  BOOST_CHECK(disk.classify({3., 0.}) == WireClassification::Outside);
  BOOST_CHECK(disk.classify({0., 3.}) == WireClassification::Outside);
  BOOST_CHECK(disk.classify({2., 0.}) == WireClassification::Boundary);
  BOOST_CHECK(disk.classify({0., -2.}) == WireClassification::Boundary);

  // a clockwise circle used as an outer wire is re-oriented to counter-clockwise.
  WireStatus reversedStatus = WireStatus::Valid;
  CurveWire reversedDisk;
  BOOST_REQUIRE(reversedDisk.initialize({Curve2D::makeCircle({0., 0.}, 2., true)}, WireRole::Outer, reversedStatus));
  BOOST_CHECK(reversedStatus == WireStatus::Reversed);
  checkClose(reversedDisk.signedArea(), surf::kPi * 4., 1.e-9);

  // --- annulus: outer disk (CCW) minus an inner hole wire (CW) --------------------------------
  WireStatus outerStatus = WireStatus::Valid;
  WireStatus holeStatus = WireStatus::Valid;
  CurveWire outerRing;
  CurveWire innerRing;
  BOOST_REQUIRE(outerRing.initialize({Curve2D::makeCircle({0., 0.}, 3.)}, WireRole::Outer, outerStatus));
  BOOST_REQUIRE(innerRing.initialize({Curve2D::makeCircle({0., 0.}, 1.)}, WireRole::Inner, holeStatus));
  BOOST_CHECK(holeStatus == WireStatus::Reversed); // CCW circle normalized to CW for a hole
  BOOST_CHECK_LT(innerRing.signedArea(), 0.);
  // net annulus area = pi * (R^2 - r^2).
  checkClose(outerRing.signedArea() + innerRing.signedArea(), surf::kPi * (9. - 1.), 1.e-9);

  // a point in the material (between radii) is inside the outer ring and outside the inner hole.
  const surf::Vec2 materialPoint{2., 0.};
  BOOST_CHECK(outerRing.classify(materialPoint) == WireClassification::Inside);
  BOOST_CHECK(innerRing.classify(materialPoint) == WireClassification::Outside);
  // a point inside the hole is inside both rings (so subtracted from the material).
  const surf::Vec2 holePoint{0.2, 0.};
  BOOST_CHECK(outerRing.classify(holePoint) == WireClassification::Inside);
  BOOST_CHECK(innerRing.classify(holePoint) == WireClassification::Inside);

  // --- mixed line + arc loop: a stadium / half-disk closed by a diameter ---------------------
  // upper half-disk: diameter along v = 0 from (-2,0) to (2,0), closed by a CCW semicircle.
  WireStatus halfStatus = WireStatus::Valid;
  CurveWire halfDisk;
  const std::vector<Curve2D> halfDiskCurves{Curve2D::makeLine({-2., 0.}, {2., 0.}),
                                            Curve2D::makeArc({0., 0.}, 2., 0., surf::kPi)};
  BOOST_REQUIRE(halfDisk.initialize(halfDiskCurves, WireRole::Outer, halfStatus));
  BOOST_CHECK(halfStatus == WireStatus::Valid);
  checkClose(halfDisk.signedArea(), 0.5 * surf::kPi * 4., 1.e-9); // half of pi*r^2
  BOOST_CHECK(halfDisk.classify({0., 1.}) == WireClassification::Inside);
  BOOST_CHECK(halfDisk.classify({0., -1.}) == WireClassification::Outside);
  BOOST_CHECK(halfDisk.classify({0., 0.}) == WireClassification::Boundary);

  // an open curve loop is rejected.
  WireStatus openStatus = WireStatus::Valid;
  CurveWire openWire;
  const std::vector<Curve2D> openCurves{Curve2D::makeLine({0., 0.}, {2., 0.}),
                                        Curve2D::makeLine({2., 0.}, {2., 2.})};
  BOOST_CHECK(!openWire.initialize(openCurves, WireRole::Outer, openStatus));
  BOOST_CHECK(openStatus == WireStatus::Open);
}