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
#include "DetectorsBase/O2SurfaceSolidIO.h"

#include "../src/BoundedSurface.h"

#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoCone.h"
#include "TGeoShape.h"
#include "TGeoSphere.h"
#include "TGeoTorus.h"
#include "TGeoTube.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
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

using BoundaryCurve = SurfaceSolid::PlanarBoundaryCurve;

// A full-circle boundary wire centred at (0,0) as a single +/-2pi arc (clockwise for holes).
std::vector<BoundaryCurve> circleWire(double radius, bool clockwise = false)
{
  return {BoundaryCurve::makeArc({0., 0.}, radius, 0., clockwise ? -surf::kTwoPi : surf::kTwoPi)};
}

// A rectangular trim loop in a quadric's (u, v) parametric domain, as four line boundary curves
// (u = phi, v = height or theta). Wound counter-clockwise; the kernel reorients as needed.
std::vector<BoundaryCurve> paramRectWire(double uMin, double uMax, double vMin, double vMax)
{
  return {BoundaryCurve::makeLine({uMin, vMin}, {uMax, vMin}), BoundaryCurve::makeLine({uMax, vMin}, {uMax, vMax}),
          BoundaryCurve::makeLine({uMax, vMax}, {uMin, vMax}), BoundaryCurve::makeLine({uMin, vMax}, {uMin, vMin})};
}

// Same rectangle as paramRectWire but as internal Curve2D segments, for direct kernel-level tests.
std::vector<surf::Curve2D> paramRectWireCurves(double uMin, double uMax, double vMin, double vMax)
{
  return {surf::Curve2D::makeLine({uMin, vMin}, {uMax, vMin}), surf::Curve2D::makeLine({uMax, vMin}, {uMax, vMax}),
          surf::Curve2D::makeLine({uMax, vMax}, {uMin, vMax}), surf::Curve2D::makeLine({uMin, vMax}, {uMin, vMin})};
}

// Add a planar disk (or annulus when holeRadius > 0) via the general curved-planar API,
// replacing the retired AddPlanarDiskSurface convenience.
bool addDiskSurface(SurfaceSolid& solid, const Point3D& center, const Point3D& axisU, const Point3D& axisV,
                    double radius, double holeRadius = 0.)
{
  std::vector<std::vector<BoundaryCurve>> inners;
  if (holeRadius > 0.) {
    inners.push_back(circleWire(holeRadius, true)); // clockwise hole: no reorientation needed
  }
  return solid.AddCurvedPlanarSurface(center, axisU, axisV, circleWire(radius), inners);
}

// Local frame (origin + parametric axes + rectangle extents) of a box face by index
// (0:+x 1:-x 2:+y 3:-y 4:+z 5:-z), for a box centred at the origin.
struct FaceFrame {
  Point3D origin;
  Point3D axisU;
  Point3D axisV;
  double extentU;
  double extentV;
};

FaceFrame boxFaceFrame(int faceIndex, double halfX, double halfY, double halfZ)
{
  switch (faceIndex) {
    case 0:
      return {{halfX, -halfY, -halfZ}, {0., 1., 0.}, {0., 0., 1.}, 2. * halfY, 2. * halfZ};
    case 1:
      return {{-halfX, -halfY, -halfZ}, {0., 0., 1.}, {0., 1., 0.}, 2. * halfZ, 2. * halfY};
    case 2:
      return {{-halfX, halfY, -halfZ}, {0., 0., 1.}, {1., 0., 0.}, 2. * halfZ, 2. * halfX};
    case 3:
      return {{-halfX, -halfY, -halfZ}, {1., 0., 0.}, {0., 0., 1.}, 2. * halfX, 2. * halfZ};
    case 4:
      return {{-halfX, -halfY, halfZ}, {1., 0., 0.}, {0., 1., 0.}, 2. * halfX, 2. * halfY};
    default:
      return {{-halfX, -halfY, -halfZ}, {0., 1., 0.}, {1., 0., 0.}, 2. * halfY, 2. * halfX};
  }
}

// Add a single box face by index of a box centred at "center". When "reversed" is set the
// face's parametric axes are swapped, which flips the outward normal inward without changing
// the covered rectangle - used to build an orientation-inconsistent fixture.
bool addBoxFace(SurfaceSolid& solid, int faceIndex, double halfX, double halfY, double halfZ, bool reversed = false,
                const Point3D& center = {0., 0., 0.})
{
  FaceFrame frame = boxFaceFrame(faceIndex, halfX, halfY, halfZ);
  if (reversed) {
    std::swap(frame.axisU, frame.axisV);
    std::swap(frame.extentU, frame.extentV);
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    frame.origin[dimension] += center[dimension];
  }
  return solid.AddPlanarSurface(frame.origin, frame.axisU, frame.axisV, rectangleWire(frame.extentU, frame.extentV));
}

void addBoxSurfaces(SurfaceSolid& solid, double halfX, double halfY, double halfZ,
                    const Point3D& center = {0., 0., 0.})
{
  for (int faceIndex = 0; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(solid, faceIndex, halfX, halfY, halfZ, false, center));
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

std::array<double, 3> unitDirection(double x, double y, double z)
{
  const double length = std::sqrt(x * x + y * y + z * z);
  return {x / length, y / length, z / length};
}

// Compare Contains against a reference ROOT shape on a regular grid. The fractional offsets keep
// grid points away from exact shape boundaries, where inside/outside conventions may differ.
void compareContainsGrid(const SurfaceSolid& solid, const TGeoShape& reference, double extent, int samples)
{
  for (int stepX = 0; stepX < samples; ++stepX) {
    for (int stepY = 0; stepY < samples; ++stepY) {
      for (int stepZ = 0; stepZ < samples; ++stepZ) {
        const double point[3] = {-extent + 2. * extent * (stepX + 0.517) / samples,
                                 -extent + 2. * extent * (stepY + 0.263) / samples,
                                 -extent + 2. * extent * (stepZ + 0.741) / samples};
        BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
        {
          BOOST_CHECK_EQUAL(solid.Contains(point), reference.Contains(point));
        }
      }
    }
  }
}

// Compare the direction-appropriate distance function against a reference ROOT shape.
void compareDistance(const SurfaceSolid& solid, const TGeoShape& reference, const std::array<double, 3>& point,
                     const std::array<double, 3>& direction, double tolerance = 1.e-9)
{
  BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ") direction = ("
                                 << direction[0] << ", " << direction[1] << ", " << direction[2] << ")")
  {
    const bool inside = reference.Contains(point.data());
    BOOST_CHECK_EQUAL(solid.Contains(point.data()), inside);
    if (inside) {
      checkClose(solid.DistFromInside(point.data(), direction.data(), 3),
                 reference.DistFromInside(point.data(), direction.data(), 3), tolerance);
    } else {
      checkClose(solid.DistFromOutside(point.data(), direction.data(), 3),
                 reference.DistFromOutside(point.data(), direction.data(), 3), tolerance);
    }
  }
}

/// @name Closed fixtures for the navigation sweeps
///
/// The distance tests exercise the same solids from several directions rather than one shape per
/// case, so the fixtures are built once here. Each is closed and therefore BVH-backed. The solid
/// is neither copyable nor movable, hence the unique_ptr.
/// @{

std::unique_ptr<SurfaceSolid> makeBoxSolid(const char* name, double halfX, double halfY, double halfZ)
{
  auto solid = std::make_unique<SurfaceSolid>(name);
  addBoxSurfaces(*solid, halfX, halfY, halfZ);
  solid->CloseShape();
  return solid;
}

// innerRadius > 0 gives a hollow tube (an inner wall plus annular caps).
std::unique_ptr<SurfaceSolid> makeTubeSolid(const char* name, double innerRadius, double outerRadius,
                                            double halfHeight)
{
  auto solid = std::make_unique<SurfaceSolid>(name);
  BOOST_REQUIRE(solid->AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, outerRadius, -halfHeight,
                                             halfHeight));
  if (innerRadius > 0.) {
    BOOST_REQUIRE(solid->AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, innerRadius, -halfHeight,
                                               halfHeight, 0., surf::kTwoPi, true));
  }
  BOOST_REQUIRE(addDiskSurface(*solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, outerRadius, innerRadius));
  BOOST_REQUIRE(addDiskSurface(*solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, outerRadius, innerRadius));
  solid->CloseShape();
  return solid;
}

std::unique_ptr<SurfaceSolid> makeConeSolid(const char* name, double radiusAtBottom, double radiusAtTop,
                                            double halfHeight)
{
  auto solid = std::make_unique<SurfaceSolid>(name);
  BOOST_REQUIRE(solid->AddConicalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radiusAtBottom, radiusAtTop,
                                         -halfHeight, halfHeight));
  BOOST_REQUIRE(addDiskSurface(*solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radiusAtTop));
  BOOST_REQUIRE(addDiskSurface(*solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radiusAtBottom));
  solid->CloseShape();
  return solid;
}

std::unique_ptr<SurfaceSolid> makeSphereSolid(const char* name, double radius)
{
  auto solid = std::make_unique<SurfaceSolid>(name);
  BOOST_REQUIRE(solid->AddSphericalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius));
  solid->CloseShape();
  return solid;
}

std::unique_ptr<SurfaceSolid> makeTorusSolid(const char* name, double majorRadius, double minorRadius)
{
  auto solid = std::make_unique<SurfaceSolid>(name);
  BOOST_REQUIRE(solid->AddToroidalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, majorRadius, minorRadius));
  solid->CloseShape();
  return solid;
}

// Cylinder barrel closed by two hemispherical endcaps; a mixed-quadric solid with no ROOT
// primitive equivalent, so the loop oracle is the only reference it has.
std::unique_ptr<SurfaceSolid> makeCapsuleSolid(const char* name, double radius, double halfHeight)
{
  auto solid = std::make_unique<SurfaceSolid>(name);
  BOOST_REQUIRE(solid->AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                             halfHeight));
  BOOST_REQUIRE(solid->AddSphericalSurface({0., 0., halfHeight}, {0., 0., 1.}, {1., 0., 0.}, radius, 0.,
                                           surf::kPi / 2.));
  BOOST_REQUIRE(solid->AddSphericalSurface({0., 0., -halfHeight}, {0., 0., -1.}, {1., 0., 0.}, radius, 0.,
                                           surf::kPi / 2.));
  solid->CloseShape();
  return solid;
}

/// @}

// Directions probed by the navigation sweeps: the three axes both ways, face and body diagonals,
// and a few skew directions that align with no symmetry of any fixture.
const std::vector<std::array<double, 3>>& probeDirections()
{
  static const std::vector<std::array<double, 3>> directions{
    {1., 0., 0.},   {-1., 0., 0.},   {0., 1., 0.},    {0., -1., 0.},  {0., 0., 1.},
    {0., 0., -1.},  unitDirection(1., 1., 0.),        unitDirection(1., 0., 1.),
    unitDirection(0., 1., 1.),       unitDirection(1., 1., 1.),       unitDirection(-1., 1., -1.),
    unitDirection(0.37, -0.82, 0.44), unitDirection(-0.91, 0.13, 0.39), unitDirection(0.21, 0.55, -0.81)};
  return directions;
}

// A deterministic point grid over the cube of half-side "extent". The fractional offsets are the
// same irrational-looking shifts the Contains sweeps use, which keeps samples off exact symmetry
// planes and shape boundaries.
std::vector<std::array<double, 3>> probeGrid(double extent, int samples)
{
  std::vector<std::array<double, 3>> points;
  points.reserve(static_cast<size_t>(samples) * samples * samples);
  for (int stepX = 0; stepX < samples; ++stepX) {
    for (int stepY = 0; stepY < samples; ++stepY) {
      for (int stepZ = 0; stepZ < samples; ++stepZ) {
        points.push_back({-extent + 2. * extent * (stepX + 0.517) / samples,
                          -extent + 2. * extent * (stepY + 0.263) / samples,
                          -extent + 2. * extent * (stepZ + 0.741) / samples});
      }
    }
  }
  return points;
}

/// The BVH distance queries must return *exactly* what the all-surfaces loop returns. Both run
/// the same analytic kernels on the same patches and take a minimum over the same hit set; they
/// differ only in the order surfaces are visited and in which of them the BVH lets them skip. So
/// any difference at all -- not merely one above a tolerance -- is a traversal or pruning bug,
/// and exact comparison is the sharpest available oracle. Independent of any mesh reference.
///
/// Ray tmax tightening is an optimization and nothing else, so both settings are checked and must
/// agree with the same loop value.
void checkDistanceAgainstLoop(const SurfaceSolid& solid, const std::array<double, 3>& point,
                              const std::array<double, 3>& direction, double stepmax = TGeoShape::Big())
{
  const double loopOutside = solid.DistFromOutside_Loop(point.data(), direction.data(), stepmax);
  const double loopInside = solid.DistFromInside_Loop(point.data(), direction.data(), stepmax);
  for (const bool pruning : {true, false}) {
    SurfaceSolid::SetRayTMaxPruning(pruning);
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ") direction = ("
                                   << direction[0] << ", " << direction[1] << ", " << direction[2]
                                   << ") stepmax = " << stepmax << " pruning = " << pruning)
    {
      BOOST_CHECK_EQUAL(solid.DistFromOutside(point.data(), direction.data(), 3, stepmax), loopOutside);
      BOOST_CHECK_EQUAL(solid.DistFromInside(point.data(), direction.data(), 3, stepmax), loopInside);
    }
  }
  SurfaceSolid::SetRayTMaxPruning(true);
}

// Sweep every grid point against every probe direction, cross-checking BVH against the loop.
void sweepDistanceAgainstLoop(const SurfaceSolid& solid, double extent, int samples)
{
  for (const auto& point : probeGrid(extent, samples)) {
    for (const auto& direction : probeDirections()) {
      checkDistanceAgainstLoop(solid, point, direction);
    }
  }
}

// Sweep both distance functions against a reference ROOT primitive, using each point in the role
// (inside/outside) the reference itself assigns it. Points closer than "skin" to the reference
// boundary are skipped: there the two shapes may legitimately disagree on which side the point is
// on, and the resulting distances are then answers to different questions.
void sweepDistanceAgainstReference(const SurfaceSolid& solid, const TGeoShape& reference, double extent, int samples,
                                   double tolerance = 1.e-9, double skin = 1.e-6)
{
  for (const auto& point : probeGrid(extent, samples)) {
    const bool inside = reference.Contains(point.data());
    if (reference.Safety(point.data(), inside) < skin) {
      continue;
    }
    for (const auto& direction : probeDirections()) {
      BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ") direction = ("
                                     << direction[0] << ", " << direction[1] << ", " << direction[2] << ")")
      {
        if (inside) {
          checkClose(solid.DistFromInside(point.data(), direction.data(), 3),
                     reference.DistFromInside(point.data(), direction.data(), 3), tolerance);
        } else {
          checkClose(solid.DistFromOutside(point.data(), direction.data(), 3),
                     reference.DistFromOutside(point.data(), direction.data(), 3), tolerance);
        }
      }
    }
  }
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

  // safeties against analytic distances (TGeo safeties may be weaker underestimates)
  const double outsideSafetyPoint[3] = {2.5, 0., 0.};
  checkClose(solid.Safety(fromCenter, kTRUE), halfX);
  checkClose(solid.Safety(outsideSafetyPoint, kFALSE), 2.5 - halfX);

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

// The queryable navigation-reliability state (scripts/geometry/ExactTrimTopology.md, item 4).
// A caller must be able to ask "can I trust this solid's navigation answers" and get a single
// answer, rather than having to notice a printed warning; the state must also survive being
// closed with check==false, since diagnostics and reporting are separate concerns.
BOOST_AUTO_TEST_CASE(NavigationReliabilityIsQueryable)
{
  using Reliability = SurfaceSolid::NavigationReliability;
  constexpr double halfX = 1.;
  constexpr double halfY = 1.5;
  constexpr double halfZ = 2.;

  // before CloseShape there are no diagnostics at all
  SurfaceSolid fresh("freshBox");
  for (int faceIndex = 0; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(fresh, faceIndex, halfX, halfY, halfZ));
  }
  BOOST_CHECK(fresh.GetNavigationReliability() == Reliability::Undetermined);
  BOOST_CHECK(!fresh.IsNavigable());

  fresh.CloseShape(false);
  BOOST_CHECK(fresh.GetNavigationReliability() == Reliability::Reliable);
  BOOST_CHECK(fresh.IsNavigable());
  BOOST_CHECK_EQUAL(fresh.GetBoundaryEdgeCount(), 0);
  BOOST_CHECK_EQUAL(fresh.GetNonManifoldEdgeCount(), 0);
  BOOST_CHECK_EQUAL(fresh.GetReversedEdgeCount(), 0);

  // a missing face leaves boundary edges: the gap case that motivates the whole state
  SurfaceSolid missing("missingFaceBoxState");
  for (int faceIndex = 1; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(missing, faceIndex, halfX, halfY, halfZ));
  }
  missing.CloseShape(false);
  BOOST_CHECK(missing.GetNavigationReliability() == Reliability::OpenSurfaceSet);
  BOOST_CHECK(!missing.IsNavigable());
  BOOST_CHECK(missing.GetBoundaryEdgeCount() > 0);

  // a reversed face is closed but inconsistently oriented
  SurfaceSolid reversed("reversedFaceBoxState");
  BOOST_REQUIRE(addBoxFace(reversed, 0, halfX, halfY, halfZ, true));
  for (int faceIndex = 1; faceIndex < 6; ++faceIndex) {
    BOOST_REQUIRE(addBoxFace(reversed, faceIndex, halfX, halfY, halfZ));
  }
  reversed.CloseShape(false);
  BOOST_CHECK(reversed.GetNavigationReliability() == Reliability::ReversedFaces);
  BOOST_CHECK(!reversed.IsNavigable());
  BOOST_CHECK(reversed.GetReversedEdgeCount() > 0);

  // duplicated faces: every edge is now shared by four faces. Non-manifold outranks the boundary
  // and orientation cases because parity is not even order-independent on such input.
  SurfaceSolid duplicated("duplicatedFaceBox");
  for (int pass = 0; pass < 2; ++pass) {
    for (int faceIndex = 0; faceIndex < 6; ++faceIndex) {
      BOOST_REQUIRE(addBoxFace(duplicated, faceIndex, halfX, halfY, halfZ));
    }
  }
  duplicated.CloseShape(false);
  BOOST_CHECK(duplicated.GetNavigationReliability() == Reliability::NonManifold);
  BOOST_CHECK(!duplicated.IsNavigable());
  BOOST_CHECK(duplicated.GetNonManifoldEdgeCount() > 0);

  BOOST_CHECK_EQUAL(std::string(SurfaceSolid::GetNavigationReliabilityName(Reliability::Reliable)), "reliable");
  BOOST_CHECK_EQUAL(std::string(SurfaceSolid::GetNavigationReliabilityName(Reliability::OpenSurfaceSet)),
                    "open-surface-set");
  BOOST_CHECK_EQUAL(std::string(SurfaceSolid::GetNavigationReliabilityName(Reliability::NonManifold)), "non-manifold");
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

BOOST_AUTO_TEST_CASE(BSplineTrimCurveKernels)
{
  using surf::Curve2D;
  using surf::CurveWire;
  using surf::Vec2;
  using surf::WireClassification;
  using surf::WireRole;
  using surf::WireStatus;

  // --- non-rational cubic B-spline: validity, clamped endpoints, convex-hull bounds ------------
  const std::vector<Vec2> poles{{0., 0.}, {1., 2.}, {2., -1.}, {3., 1.}, {4., 0.}};
  const std::vector<double> knots{0., 0., 0., 0., 0.5, 1., 1., 1., 1.};
  const Curve2D spline = Curve2D::makeBSpline(3, poles, {}, knots);
  BOOST_CHECK(spline.valid());
  checkClose(spline.startPoint().uCoord, 0.);
  checkClose(spline.startPoint().vCoord, 0.);
  checkClose(spline.endPoint().uCoord, 4.);
  checkClose(spline.endPoint().vCoord, 0.);
  // extendBounds returns the (conservative) control-point convex-hull box
  Vec2 lower{1.e30, 1.e30};
  Vec2 upper{-1.e30, -1.e30};
  spline.extendBounds(lower, upper);
  checkClose(lower.uCoord, 0.);
  checkClose(upper.uCoord, 4.);
  checkClose(lower.vCoord, -1.);
  checkClose(upper.vCoord, 2.);

  // --- rational quadratic B-spline: an exact NURBS quarter circle -----------------------------
  const std::vector<Vec2> circlePoles{{1., 0.}, {1., 1.}, {0., 1.}};
  const std::vector<double> circleWeights{1., std::sqrt(0.5), 1.};
  const std::vector<double> circleKnots{0., 0., 0., 1., 1., 1.};
  const Curve2D quarter = Curve2D::makeBSpline(2, circlePoles, circleWeights, circleKnots);
  BOOST_CHECK(quarter.valid());
  for (int index = 0; index <= 8; ++index) {
    const Vec2 point = quarter.pointAt(static_cast<double>(index) / 8);
    checkClose(std::hypot(point.uCoord, point.vCoord), 1., 1.e-9);
  }

  // --- closed loop (B-spline top + three closing lines): area vs a fine-polygon reference ------
  const std::vector<Curve2D> loop{spline, Curve2D::makeLine({4., 0.}, {4., -3.}),
                                  Curve2D::makeLine({4., -3.}, {0., -3.}), Curve2D::makeLine({0., -3.}, {0., 0.})};
  WireStatus status = WireStatus::Valid;
  CurveWire wire;
  BOOST_REQUIRE(wire.initialize(loop, WireRole::Outer, status));
  double referenceArea = 0.;
  const auto boundarySamples = wire.sampledBoundary();
  for (size_t k = 0; k + 1 < boundarySamples.size(); ++k) {
    referenceArea += 0.5 * (boundarySamples[k].uCoord * boundarySamples[k + 1].vCoord -
                            boundarySamples[k + 1].uCoord * boundarySamples[k].vCoord);
  }
  // wire.signedArea() is the exact Gauss-Legendre value; referenceArea is a chord-polyline
  // approximation of it, so compare at the sampling-accuracy level rather than machine precision
  checkClose(wire.signedArea(), std::abs(referenceArea), 1.e-4);

  // classify inside / outside / boundary
  BOOST_CHECK(wire.classify({2., -1.5}) == WireClassification::Inside);
  BOOST_CHECK(wire.classify({2., -2.9}) == WireClassification::Inside);
  BOOST_CHECK(wire.classify({-1., -1.}) == WireClassification::Outside);
  BOOST_CHECK(wire.classify({2., 5.}) == WireClassification::Outside);
  BOOST_CHECK(wire.classify({0., 0.}) == WireClassification::Boundary);  // on the B-spline start
  BOOST_CHECK(wire.classify({2., -3.}) == WireClassification::Boundary); // on the bottom line

  // --- horizontal-tangent case: a scanline tangent to a smooth apex must not flip parity -------
  // downward arch: quadratic B-spline (0,0) -> apex (1,1) -> (2,0), closed by the baseline.
  const Curve2D arch = Curve2D::makeBSpline(2, {{0., 0.}, {1., 2.}, {2., 0.}}, {}, {0., 0., 0., 1., 1., 1.});
  checkClose(arch.pointAt(0.5).vCoord, 1.); // apex height
  WireStatus archStatus = WireStatus::Valid;
  CurveWire archRegion;
  BOOST_REQUIRE(archRegion.initialize({arch, Curve2D::makeLine({2., 0.}, {0., 0.})}, WireRole::Outer, archStatus));
  BOOST_CHECK(archRegion.classify({1., 0.5}) == WireClassification::Inside);
  BOOST_CHECK(archRegion.classify({1., 1.5}) == WireClassification::Outside);
  // scanline v = 1 is tangent to the apex to the right of these points: a robust kernel counts an
  // even number of crossings so both points classify Outside.
  BOOST_CHECK(archRegion.classify({-1., 1.}) == WireClassification::Outside);
  BOOST_CHECK(archRegion.classify({3., 1.}) == WireClassification::Outside);

  // --- reversal keeps the same geometric image (poles/knots complemented) ----------------------
  Curve2D reversed = spline;
  reversed.reverseInPlace();
  checkClose(reversed.startPoint().uCoord, 4.);
  checkClose(reversed.endPoint().uCoord, 0.);
  checkClose(reversed.pointAt(0.25).uCoord, spline.pointAt(0.75).uCoord, 1.e-9);
  checkClose(reversed.pointAt(0.25).vCoord, spline.pointAt(0.75).vCoord, 1.e-9);
}

BOOST_AUTO_TEST_CASE(CurvedPlanarDiskKernels)
{
  using surf::Curve2D;

  // annulus in the z = 0 plane: outer radius 2, hole radius 1
  surf::CurvedPlanarBoundedSurface annulus;
  std::string error;
  BOOST_REQUIRE(annulus.initialize({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.},
                                   {Curve2D::makeCircle({0., 0.}, 2.)},
                                   {{Curve2D::makeCircle({0., 0.}, 1., true)}}, error));
  BOOST_CHECK(!annulus.wasReoriented()); // outer CCW, hole CW: both already correctly oriented

  // a skewed (non-orthonormal) frame is rejected
  surf::CurvedPlanarBoundedSurface skewed;
  BOOST_CHECK(!skewed.initialize({0., 0., 0.}, {1., 0., 0.}, {0.5, 1., 0.},
                                 {Curve2D::makeCircle({0., 0.}, 2.)}, {}, error));

  // on-surface classification: material, hole, outside, off-plane
  BOOST_CHECK(annulus.containsPointOnSurface({1.5, 0., 0.}));
  BOOST_CHECK(!annulus.containsPointOnSurface({0.5, 0., 0.}));
  BOOST_CHECK(!annulus.containsPointOnSurface({3., 0., 0.}));
  BOOST_CHECK(!annulus.containsPointOnSurface({1.5, 0., 0.5}));

  // ray intersections: one hit through the material, none through the hole
  std::vector<surf::RayHit> hits;
  annulus.appendIntersections({1.5, 0., 1.}, {0., 0., -1.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 1u);
  checkClose(hits.front().distance, 1.);
  checkClose(hits.front().normal.zCoord, 1.);
  hits.clear();
  annulus.appendIntersections({0.5, 0., 1.}, {0., 0., -1.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());

  // exact patch distances: in-hole (to the hole rim), in-material off-plane, outside the rim,
  // and the combined in-plane plus out-of-plane case
  checkClose(annulus.distanceSqToPatch({0., 0., 0.}), 1.);
  checkClose(annulus.distanceSqToPatch({1.5, 0., 2.}), 4.);
  checkClose(annulus.distanceSqToPatch({4., 0., 0.}), 4.);
  checkClose(annulus.distanceSqToPatch({0.5, 0., 1.}), 1.25);

  // divergence-theorem contribution of an offset disk: (origin . normal) * area / 3
  surf::CurvedPlanarBoundedSurface offsetDisk;
  BOOST_REQUIRE(offsetDisk.initialize({0., 0., 2.}, {1., 0., 0.}, {0., 1., 0.},
                                      {Curve2D::makeCircle({0., 0.}, 1.5)}, {}, error));
  checkClose(offsetDisk.capacityContribution(), 2. * surf::kPi * 1.5 * 1.5 / 3., 1.e-9);
  BOOST_CHECK(offsetDisk.capacityIsExact());
}

BOOST_AUTO_TEST_CASE(CylindricalSurfaceKernels)
{
  // full lateral cylinder, radius 2, height [-3, 3], axis z
  surf::CylindricalBoundedSurface cylinder;
  std::string error;
  BOOST_REQUIRE(cylinder.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -3., 3., 0., surf::kTwoPi,
                                    false, error));

  // a transversal ray crosses the lateral surface twice: both hits must be reported
  std::vector<surf::RayHit> hits;
  cylinder.appendIntersections({-5., 0.5, 0.}, {1., 0., 0.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 2u);
  const double chordHalf = std::sqrt(4. - 0.25);
  checkClose(hits[0].distance, 5. - chordHalf);
  checkClose(hits[1].distance, 5. + chordHalf);
  // entering hit: outward normal opposes the ray direction; exiting hit: aligned
  BOOST_CHECK_LT(hits[0].normal.xCoord, 0.);
  BOOST_CHECK_GT(hits[1].normal.xCoord, 0.);

  // tangential graze reports no hits (keeps crossing parity even)
  hits.clear();
  cylinder.appendIntersections({-5., 2., 0.}, {1., 0., 0.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());

  // axis-parallel ray never crosses the lateral surface
  hits.clear();
  cylinder.appendIntersections({0., 0., -5.}, {0., 0., 1.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());

  // exact patch distances: radial, above the rim, and the diagonal rim case
  checkClose(cylinder.distanceSqToPatch({4., 0., 0.}), 4.);
  checkClose(cylinder.distanceSqToPatch({0., 0., 5.}), 8.);
  checkClose(cylinder.distanceSqToPatch({3., 0., 4.}), 2.);

  // half cylinder (phi in [0, pi]): the phi trim filters hits and surface points
  surf::CylindricalBoundedSurface halfCylinder;
  BOOST_REQUIRE(halfCylinder.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -3., 3., 0., surf::kPi,
                                        false, error));
  BOOST_CHECK(halfCylinder.containsPointOnSurface({0., 2., 0.}));
  BOOST_CHECK(!halfCylinder.containsPointOnSurface({0., -2., 0.}));
  hits.clear();
  halfCylinder.appendIntersections({0., -5., 0.}, {0., 1., 0.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 1u);
  checkClose(hits.front().distance, 7.);
}

BOOST_AUTO_TEST_CASE(ClosedCylinderMatchesTGeoTube)
{
  constexpr double radius = 2.;
  constexpr double halfHeight = 3.;

  SurfaceSolid solid("closedCylinder");
  BOOST_REQUIRE(solid.AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                            halfHeight));
  // cap frames: outward normal is axisU x axisV, so the bottom cap flips axisV
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radius));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radius));
  solid.CloseShape();

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoTube reference("referenceTube", 0., radius, halfHeight);
  compareContainsGrid(solid, reference, 4., 9);

  compareDistance(solid, reference, {5., 0., 0.}, {-1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 5.}, {0., 0., -1.});
  compareDistance(solid, reference, {5., 0.5, 1.}, {-1., 0., 0.});
  compareDistance(solid, reference, {-4., -1., -2.}, unitDirection(1., 0.3, 0.5));
  compareDistance(solid, reference, {0., 0., 0.}, {1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 0.}, {0., 0., 1.});
  compareDistance(solid, reference, {0., 0., 0.}, unitDirection(1., 1., 1.));
  compareDistance(solid, reference, {1., 0.5, -2.}, unitDirection(0.3, -0.4, 0.5));
  compareDistance(solid, reference, {5., 2.5, 0.}, {-1., 0., 0.}); // grazing miss

  // safeties against analytic distances (TGeo safeties may be weaker underestimates, so they
  // are not compared directly)
  const double center[3] = {0., 0., 0.};
  const double insidePoint[3] = {1., 0.5, 1.};
  const double radialOutside[3] = {4., 0., 0.};
  const double axialOutside[3] = {0., 0., 5.};
  const double cornerOutside[3] = {4., 0., 5.};
  checkClose(solid.Safety(center, kTRUE), radius);
  checkClose(solid.Safety(insidePoint, kTRUE), radius - std::sqrt(1.25));
  checkClose(solid.Safety(radialOutside, kFALSE), 2.);
  checkClose(solid.Safety(axialOutside, kFALSE), 2.);
  checkClose(solid.Safety(cornerOutside, kFALSE), std::sqrt(8.)); // exact corner distance

  // normals on the lateral surface and the caps
  double normal[3] = {0., 0., 0.};
  const double sidePoint[3] = {radius, 0., 1.};
  const double alongX[3] = {1., 0., 0.};
  solid.ComputeNormal(sidePoint, alongX, normal);
  checkClose(normal[0], 1.);
  checkClose(normal[1], 0.);
  checkClose(normal[2], 0.);
  const double capPoint[3] = {0.5, 0.5, halfHeight};
  const double alongZ[3] = {0., 0., 1.};
  solid.ComputeNormal(capPoint, alongZ, normal);
  checkClose(normal[2], 1.);

  checkClose(solid.Capacity(), reference.Capacity(), 1.e-9);

  int meshVertices = 0;
  int meshSegments = 0;
  int meshPolygons = 0;
  solid.GetMeshNumbers(meshVertices, meshSegments, meshPolygons);
  BOOST_CHECK_GT(meshVertices, 0);
  BOOST_CHECK_GT(meshPolygons, 0);
}

BOOST_AUTO_TEST_CASE(HollowCylinderMatchesTGeoTube)
{
  constexpr double innerRadius = 1.;
  constexpr double outerRadius = 2.;
  constexpr double halfHeight = 3.;

  SurfaceSolid solid("hollowCylinder");
  BOOST_REQUIRE(solid.AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, outerRadius, -halfHeight,
                                            halfHeight));
  BOOST_REQUIRE(solid.AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, innerRadius, -halfHeight,
                                            halfHeight, 0., surf::kTwoPi, true));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, outerRadius,
                                           innerRadius));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, outerRadius,
                                           innerRadius));
  solid.CloseShape();

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoTube reference("referenceHollowTube", innerRadius, outerRadius, halfHeight);
  compareContainsGrid(solid, reference, 4., 9);

  // from the hole the solid is entered through the inner wall
  compareDistance(solid, reference, {0., 0., 0.}, {1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 2.}, unitDirection(0.4, 0.2, -1.));
  // inside the material both walls are exit candidates
  compareDistance(solid, reference, {1.5, 0., 0.}, {1., 0., 0.});
  compareDistance(solid, reference, {1.5, 0., 0.}, {-1., 0., 0.});
  compareDistance(solid, reference, {-1.2, 0.8, 1.}, unitDirection(-0.2, 0.9, 0.4));
  compareDistance(solid, reference, {5., 0., 0.}, {-1., 0., 0.});

  // analytic safety in the middle of the material: 0.5 to either wall
  const double materialPoint[3] = {1.5, 0., 0.};
  checkClose(solid.Safety(materialPoint, kTRUE), 0.5);

  checkClose(solid.Capacity(), reference.Capacity(), 1.e-9);
}

BOOST_AUTO_TEST_CASE(SphereMatchesTGeoSphere)
{
  constexpr double radius = 2.5;

  SurfaceSolid solid("fullSphere");
  BOOST_REQUIRE(solid.AddSphericalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius));
  solid.CloseShape();

  // a full sphere is self-closing: no boundary edges at all
  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoSphere reference("referenceSphere", 0., radius);
  compareContainsGrid(solid, reference, 3.5, 9);

  compareDistance(solid, reference, {5., 0., 0.}, {-1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 0.}, {1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 0.}, unitDirection(1., 1., 1.));
  compareDistance(solid, reference, {1., 1., 1.}, unitDirection(-0.3, 0.5, 0.8));
  compareDistance(solid, reference, {-4., 0.5, 0.5}, {1., 0., 0.});
  compareDistance(solid, reference, {-4., 2.6, 0.}, {1., 0., 0.}); // clean miss

  // analytic safeties: |distance to center - radius|
  const double insidePoint[3] = {1., 0., 0.};
  const double outsidePoint[3] = {4., 0., 0.};
  checkClose(solid.Safety(insidePoint, kTRUE), radius - 1.);
  checkClose(solid.Safety(outsidePoint, kFALSE), 4. - radius);

  double normal[3] = {0., 0., 0.};
  const double surfacePoint[3] = {radius, 0., 0.};
  const double alongX[3] = {1., 0., 0.};
  solid.ComputeNormal(surfacePoint, alongX, normal);
  checkClose(normal[0], 1.);

  checkClose(solid.Capacity(), reference.Capacity(), 1.e-9);

  int meshVertices = 0;
  int meshSegments = 0;
  int meshPolygons = 0;
  solid.GetMeshNumbers(meshVertices, meshSegments, meshPolygons);
  BOOST_CHECK_GT(meshVertices, 0);
  BOOST_CHECK_GT(meshPolygons, 0);
}

BOOST_AUTO_TEST_CASE(SphericalSectionKernels)
{
  // upper hemisphere shell of radius 2 (theta in [0, pi/2], full phi)
  surf::SphericalBoundedSurface hemisphere;
  std::string error;
  BOOST_REQUIRE(hemisphere.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., 0., surf::kHalfPi, 0.,
                                      surf::kTwoPi, false, error));

  // divergence contribution of a centred hemisphere shell: 2 pi R^3 / 3
  checkClose(hemisphere.capacityContribution(), 2. * surf::kPi * 8. / 3., 1.e-9);

  // the polar-axis ray meets the sphere twice but only the upper hit is on the patch
  std::vector<surf::RayHit> hits;
  hemisphere.appendIntersections({0., 0., 5.}, {0., 0., -1.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 1u);
  checkClose(hits.front().distance, 3.);
  checkClose(hits.front().normal.zCoord, 1.);

  // a transversal ray at z = 1 stays in the upper hemisphere: both hits reported
  hits.clear();
  hemisphere.appendIntersections({-5., 0., 1.}, {1., 0., 0.}, 0., 1.e30, hits);
  BOOST_CHECK_EQUAL(hits.size(), 2u);

  // the mirrored ray at z = -1 misses the trimmed patch entirely
  hits.clear();
  hemisphere.appendIntersections({-5., 0., -1.}, {1., 0., 0.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());

  // trim-aware surface point classification (equator lies on the trim boundary)
  BOOST_CHECK(hemisphere.containsPointOnSurface({0., 0., 2.}));
  BOOST_CHECK(hemisphere.containsPointOnSurface({2., 0., 0.}));
  BOOST_CHECK(!hemisphere.containsPointOnSurface({0., 0., -2.}));

  // patch distance: exact radially above the pole, conservative lower bound below the equator
  checkClose(hemisphere.distanceSqToPatch({0., 0., 5.}), 9.);
  BOOST_CHECK_LE(hemisphere.distanceSqToPatch({0., 0., -4.}), 4. + 1.e-9);
}

BOOST_AUTO_TEST_CASE(TruncatedConeMatchesTGeoCone)
{
  constexpr double halfHeight = 3.;
  constexpr double radiusAtBottom = 2.;
  constexpr double radiusAtTop = 1.;

  SurfaceSolid solid("truncatedCone");
  BOOST_REQUIRE(solid.AddConicalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radiusAtBottom, radiusAtTop,
                                        -halfHeight, halfHeight));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radiusAtTop));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radiusAtBottom));
  solid.CloseShape();

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoCone reference("referenceCone", halfHeight, 0., radiusAtBottom, 0., radiusAtTop);
  compareContainsGrid(solid, reference, 3.5, 9);

  compareDistance(solid, reference, {5., 0., 0.}, {-1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 5.}, {0., 0., -1.});
  compareDistance(solid, reference, {0., 0., 0.}, {1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 0.}, {0., 0., 1.});
  compareDistance(solid, reference, {0., 0., 0.}, {0., 0., -1.});
  compareDistance(solid, reference, {0.5, -0.3, 1.}, unitDirection(0.6, 0.4, 0.2));
  compareDistance(solid, reference, {-4., 0.2, -2.}, unitDirection(1., 0.05, 0.3));

  // central safety: exact distance to the lateral generator segment (2,-3)-(1,3) in (rho, z);
  // TGeoCone's safety degenerates to 0 on the axis of an rmin = 0 cone, so no direct comparison
  const double center[3] = {0., 0., 0.};
  checkClose(solid.Safety(center, kTRUE), 9. / std::sqrt(37.));

  // lateral-surface normal against the ROOT cone
  double normal[3] = {0., 0., 0.};
  double referenceNormal[3] = {0., 0., 0.};
  const double sidePoint[3] = {1.5, 0., 0.};
  const double alongX[3] = {1., 0., 0.};
  solid.ComputeNormal(sidePoint, alongX, normal);
  reference.ComputeNormal(sidePoint, alongX, referenceNormal);
  checkClose(normal[0], referenceNormal[0]);
  checkClose(normal[1], referenceNormal[1]);
  checkClose(normal[2], referenceNormal[2]);

  checkClose(solid.Capacity(), reference.Capacity(), 1.e-9);
}

BOOST_AUTO_TEST_CASE(ApexConeClosesWithSingleCap)
{
  // full cone: radius 3 at z = -1.5 shrinking to the apex at z = +1.5, closed by one cap
  constexpr double halfHeight = 1.5;
  constexpr double baseRadius = 3.;

  SurfaceSolid solid("apexCone");
  BOOST_REQUIRE(solid.AddConicalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, baseRadius, 0., -halfHeight,
                                        halfHeight));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, baseRadius));
  solid.CloseShape();

  // the apex rim degenerates to a point, so one cap suffices for a closed manifold
  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  // analytic containment: inside iff |z| < halfHeight and rho < r(z) = halfHeight - z
  const auto analyticInside = [&](double x, double y, double z) {
    return std::abs(z) < halfHeight && std::hypot(x, y) < halfHeight - z;
  };
  const std::array<std::array<double, 3>, 7> probePoints{{{0., 0., 0.},
                                                          {1., 0., 0.},
                                                          {1.4, 0., 0.5},
                                                          {0., 0., 1.4},
                                                          {0., 0., 1.6},
                                                          {2., 2., -1.},
                                                          {2., 0., -1.}}};
  for (const auto& probe : probePoints) {
    BOOST_TEST_CONTEXT("point = (" << probe[0] << ", " << probe[1] << ", " << probe[2] << ")")
    {
      BOOST_CHECK_EQUAL(solid.Contains(probe.data()), analyticInside(probe[0], probe[1], probe[2]));
    }
  }

  // radial exit through the slanted surface
  const double insidePoint[3] = {0., 0., -1.};
  const double alongX[3] = {1., 0., 0.};
  checkClose(solid.DistFromInside(insidePoint, alongX, 3), 2.5);

  // central safety is the exact distance to the slanted line rho + z = halfHeight
  const double center[3] = {0., 0., 0.};
  checkClose(solid.Safety(center, kTRUE), halfHeight / std::sqrt(2.), 1.e-9);

  // exact capacity of a full cone: pi R^2 H / 3
  checkClose(solid.Capacity(), surf::kPi * baseRadius * baseRadius * 2. * halfHeight / 3., 1.e-9);
}

BOOST_AUTO_TEST_CASE(ToroidalSurfaceKernels)
{
  // full torus, major radius 3, minor (tube) radius 1, axis z
  constexpr double majorR = 3.;
  constexpr double minorR = 1.;
  surf::TorusBoundedSurface torus;
  std::string error;
  BOOST_REQUIRE(torus.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, majorR, minorR, 0., surf::kTwoPi, 0.,
                                 surf::kTwoPi, false, error));

  // a ray along +x through the centre crosses the donut four times: at rho = -(R+r), -(R-r),
  // (R-r), (R+r), i.e. distances 6, 8, 12, 14 from the origin at x = -10
  std::vector<surf::RayHit> hits;
  torus.appendIntersections({-10., 0., 0.}, {1., 0., 0.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 4u);
  std::sort(hits.begin(), hits.end(),
            [](const surf::RayHit& a, const surf::RayHit& b) { return a.distance < b.distance; });
  checkClose(hits[0].distance, 6., 1.e-7);
  checkClose(hits[1].distance, 8., 1.e-7);
  checkClose(hits[2].distance, 12., 1.e-7);
  checkClose(hits[3].distance, 14., 1.e-7);
  // crossings alternate enter/exit/enter/exit along the ray
  BOOST_CHECK_LT(hits[0].normal.xCoord, 0.);
  BOOST_CHECK_GT(hits[1].normal.xCoord, 0.);
  BOOST_CHECK_LT(hits[2].normal.xCoord, 0.);
  BOOST_CHECK_GT(hits[3].normal.xCoord, 0.);

  // a z-ray tangent to the outer equator (rho = R + r) touches at a single double root: no hit
  hits.clear();
  torus.appendIntersections({majorR + minorR, 0., -10.}, {0., 0., 1.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());

  // a ray passing above the whole torus (z = 2 r) misses entirely
  hits.clear();
  torus.appendIntersections({-10., 0., 2. * minorR}, {1., 0., 0.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());

  // outward normals: +x at the outer equator, -x (towards the axis) at the inner equator
  const surf::Vec3 outerNormal = torus.normalAt({majorR + minorR, 0., 0.});
  checkClose(outerNormal.xCoord, 1.);
  const surf::Vec3 innerNormal = torus.normalAt({majorR - minorR, 0., 0.});
  checkClose(innerNormal.xCoord, -1.);
  // top of the tube: normal points along +z
  const surf::Vec3 topNormal = torus.normalAt({majorR, 0., minorR});
  checkClose(topNormal.zCoord, 1.);

  // exact meridian distances: radially outside the outer equator and inside the hole
  checkClose(torus.distanceSqToPatch({majorR + minorR + 2., 0., 0.}), 4.);
  checkClose(torus.distanceSqToPatch({0., 0., 0.}), (majorR - minorR) * (majorR - minorR));

  // surface-point classification
  BOOST_CHECK(torus.containsPointOnSurface({majorR + minorR, 0., 0.}));
  BOOST_CHECK(torus.containsPointOnSurface({majorR, 0., minorR}));
  BOOST_CHECK(!torus.containsPointOnSurface({majorR, 0., 0.}));       // tube spine (interior)
  BOOST_CHECK(!torus.containsPointOnSurface({majorR + 5., 0., 0.}));  // off the surface

  // exact divergence-theorem capacity of a full torus: 2 pi^2 R r^2
  checkClose(torus.capacityContribution(), 2. * surf::kPi * surf::kPi * majorR * minorR * minorR, 1.e-9);
  BOOST_CHECK(torus.capacityIsExact());

  // partial tube section (a quarter-tube fillet-like patch, phiTube in [0, pi/2], full ring):
  // the trim filters intersections and surface points
  surf::TorusBoundedSurface quarterTube;
  BOOST_REQUIRE(quarterTube.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, majorR, minorR, 0., surf::kTwoPi, 0.,
                                       surf::kHalfPi, false, error));
  BOOST_CHECK(quarterTube.containsPointOnSurface({majorR + minorR, 0., 0.})); // phiTube = 0 boundary
  BOOST_CHECK(quarterTube.containsPointOnSurface({majorR, 0., minorR}));      // phiTube = pi/2 boundary
  BOOST_CHECK(!quarterTube.containsPointOnSurface({majorR, 0., -minorR}));    // phiTube = -pi/2, off patch
  hits.clear();
  quarterTube.appendIntersections({majorR, 0., -10.}, {0., 0., 1.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 1u); // only the top (+z) tube point is on the quarter patch
  checkClose(hits.front().distance, 10. + minorR, 1.e-7);
}

BOOST_AUTO_TEST_CASE(FullTorusMatchesTGeoTorus)
{
  constexpr double majorR = 3.;
  constexpr double minorR = 1.;

  SurfaceSolid solid("fullTorus");
  BOOST_REQUIRE(solid.AddToroidalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, majorR, minorR));
  solid.CloseShape();

  // a full torus is self-closing: no boundary edges
  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  // TGeoTorus(R, Rmin, Rmax): a solid torus has Rmin = 0, Rmax = tube radius
  TGeoTorus reference("referenceTorus", majorR, 0., minorR);
  compareContainsGrid(solid, reference, 4.5, 11);

  // analytic x-axis crossings (see the kernel test): from outside and from inside the material
  const double outsidePoint[3] = {-10., 0., 0.};
  const double alongX[3] = {1., 0., 0.};
  checkClose(solid.DistFromOutside(outsidePoint, alongX, 3), 6., 1.e-7);
  const double materialPoint[3] = {majorR + minorR - 0.25, 0., 0.}; // inside the tube on the +x side
  BOOST_CHECK(solid.Contains(materialPoint));
  checkClose(solid.DistFromInside(materialPoint, alongX, 3), 0.25, 1.e-7);

  // a couple of oblique rays cross-checked against the ROOT torus
  compareDistance(solid, reference, {-10., 0.3, 0.2}, {1., 0., 0.}, 1.e-6);
  compareDistance(solid, reference, {0., 0., 5.}, {0., 0., -1.}, 1.e-6); // clean axial miss (through the hole)

  // exact capacity: 2 pi^2 R r^2
  checkClose(solid.Capacity(), reference.Capacity(), 1.e-7);
  checkClose(solid.Capacity(), 2. * surf::kPi * surf::kPi * majorR * minorR * minorR, 1.e-9);

  int meshVertices = 0;
  int meshSegments = 0;
  int meshPolygons = 0;
  solid.GetMeshNumbers(meshVertices, meshSegments, meshPolygons);
  BOOST_CHECK_GT(meshVertices, 0);
  BOOST_CHECK_GT(meshPolygons, 0);
}

BOOST_AUTO_TEST_CASE(WireTrimmedTorusMatchesSection)
{
  // A partial toroidal patch (ring [0, pi/2], tube [0.2, 2.0] - a non-wrapping fillet-like arc)
  // built two ways must classify points identically: with the scalar parametric rectangle and
  // with an equivalent (phiRing, phiTube) line-wire trim. This exercises the wire-trim path
  // (numeric capacity, conservative Safety) and the periodic-in-both-angles unwrapping.
  constexpr double majorR = 4.;
  constexpr double minorR = 1.5;
  constexpr double tubeLow = 0.2;
  constexpr double tubeHigh = 2.0;
  std::string error;

  surf::TorusBoundedSurface scalarSection;
  BOOST_REQUIRE(scalarSection.initialize({0.2, -0.1, 0.3}, {0., 0., 1.}, {1., 0., 0.}, majorR, minorR, 0.,
                                         surf::kHalfPi, tubeLow, tubeHigh - tubeLow, false, error));

  surf::TorusBoundedSurface wireSection;
  const auto wire = paramRectWireCurves(0., surf::kHalfPi, tubeLow, tubeHigh);
  BOOST_REQUIRE(wireSection.initialize({0.2, -0.1, 0.3}, {0., 0., 1.}, {1., 0., 0.}, majorR, minorR, 0., surf::kHalfPi,
                                       tubeLow, tubeHigh - tubeLow, false, wire, {}, error));
  BOOST_CHECK(wireSection.hasWireTrim());

  // classification agrees across a set of on-surface probes at several ring/tube angles
  for (double ring : {0.1, 0.7, 1.2, 1.7, 2.5}) {
    for (double tube : {0.3, 0.8, 1.5, 1.9, 2.6}) {
      const surf::Vec3 probe = scalarSection.pointAt(ring, tube);
      BOOST_TEST_CONTEXT("ring = " << ring << " tube = " << tube)
      {
        BOOST_CHECK_EQUAL(scalarSection.containsPointOnSurface(probe), wireSection.containsPointOnSurface(probe));
      }
    }
  }

  // wire-trim capacity is numeric (flagged inexact) but must approximate the exact scalar value
  BOOST_CHECK(scalarSection.capacityIsExact());
  BOOST_CHECK(!wireSection.capacityIsExact());
  BOOST_CHECK_SMALL(wireSection.capacityContribution() - scalarSection.capacityContribution(),
                    1.e-2 * std::abs(scalarSection.capacityContribution()));
}

BOOST_AUTO_TEST_CASE(BVHConstructionAndTraversal)
{
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;

  SurfaceSolid solid("bvhBox");
  addBoxSurfaces(solid, halfX, halfY, halfZ);
  BOOST_CHECK(!solid.HasBVH()); // built only in CloseShape
  solid.CloseShape();
  BOOST_REQUIRE(solid.HasBVH());

  // the BVH root box must enclose the exact solid bounds and stay conservative: not tighter
  // than the exact bounds, not looser than the documented expansion (plus float rounding)
  Point3D lower{};
  Point3D upper{};
  BOOST_REQUIRE(solid.GetBVHRootBounds(lower, upper));
  const Point3D exactLower{-halfX, -halfY, -halfZ};
  const Point3D exactUpper{halfX, halfY, halfZ};
  constexpr double boxSlack = 2. * surf::kBVHBoxTolerance;
  for (int dimension = 0; dimension < 3; ++dimension) {
    BOOST_TEST_CONTEXT("dimension = " << dimension)
    {
      BOOST_CHECK(lower[dimension] <= exactLower[dimension]);
      BOOST_CHECK(lower[dimension] >= exactLower[dimension] - boxSlack);
      BOOST_CHECK(upper[dimension] >= exactUpper[dimension]);
      BOOST_CHECK(upper[dimension] <= exactUpper[dimension] + boxSlack);
    }
  }

  // a ray through the box must traverse (at least) the entry and exit face leaves ...
  BOOST_CHECK_GE(solid.CountBVHRayCandidates({-2., 0., 0.}, {1., 0., 0.}), 2);
  // ... while a ray pointing away from the solid reaches no leaf at all
  BOOST_CHECK_EQUAL(solid.CountBVHRayCandidates({0., 5., 0.}, {0., 1., 0.}), 0);

  // two disjoint boxes: BVH pruning with well-separated primitive clusters. The union of two
  // closed manifolds is still a closed manifold, and parity containment handles it naturally.
  constexpr double half = 1.;
  constexpr double centerX = 3.;
  SurfaceSolid twoBoxes("twoBoxes");
  addBoxSurfaces(twoBoxes, half, half, half, {-centerX, 0., 0.});
  addBoxSurfaces(twoBoxes, half, half, half, {centerX, 0., 0.});
  twoBoxes.CloseShape();
  BOOST_REQUIRE(twoBoxes.HasBVH());
  BOOST_CHECK_EQUAL(twoBoxes.GetNsurfaces(), 12);
  BOOST_CHECK(twoBoxes.IsClosed());
  BOOST_CHECK(twoBoxes.IsOrientationConsistent());

  const auto analyticInside = [&](const double* point) {
    return (std::abs(std::abs(point[0]) - centerX) < half) && std::abs(point[1]) < half && std::abs(point[2]) < half;
  };
  constexpr int samples = 9;
  constexpr double extent = 5.;
  for (int stepX = 0; stepX < samples; ++stepX) {
    for (int stepY = 0; stepY < samples; ++stepY) {
      for (int stepZ = 0; stepZ < samples; ++stepZ) {
        const double point[3] = {-extent + 2. * extent * (stepX + 0.517) / samples,
                                 -extent + 2. * extent * (stepY + 0.263) / samples,
                                 -extent + 2. * extent * (stepZ + 0.741) / samples};
        BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
        {
          const bool bvhInside = twoBoxes.Contains(point);
          BOOST_CHECK_EQUAL(bvhInside, twoBoxes.Contains_Loop(point));
          BOOST_CHECK_EQUAL(bvhInside, analyticInside(point));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ContainsBoundaryPointsAndCapsule)
{
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;

  SurfaceSolid box("boundaryBox");
  addBoxSurfaces(box, halfX, halfY, halfZ);
  box.CloseShape();

  // boundary policy: points exactly on faces, edges and vertices count as inside,
  // in the BVH-accelerated path and in the trivial loop alike
  const std::array<std::array<double, 3>, 6> boundaryPoints{{{halfX, 0., 0.},          // face
                                                             {0., -halfY, 0.},        // face
                                                             {halfX, halfY, 0.},      // edge
                                                             {-halfX, 0., halfZ},     // edge
                                                             {halfX, halfY, halfZ},   // vertex
                                                             {-halfX, -halfY, -halfZ} // vertex
                                                            }};
  for (const auto& point : boundaryPoints) {
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      BOOST_CHECK(box.Contains(point.data()));
      BOOST_CHECK(box.Contains_Loop(point.data()));
    }
  }

  // capsule: cylinder barrel closed by two spherical endcaps - a mixed quadric fixture with no
  // ROOT primitive equivalent, cross-validated against the trivial loop and the analytic shape
  constexpr double radius = 1.;
  constexpr double halfHeight = 1.5;
  SurfaceSolid capsule("capsule");
  BOOST_REQUIRE(capsule.AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                              halfHeight));
  BOOST_REQUIRE(capsule.AddSphericalSurface({0., 0., halfHeight}, {0., 0., 1.}, {1., 0., 0.}, radius, 0.,
                                            surf::kPi / 2.));
  BOOST_REQUIRE(capsule.AddSphericalSurface({0., 0., -halfHeight}, {0., 0., -1.}, {1., 0., 0.}, radius, 0.,
                                            surf::kPi / 2.));
  capsule.CloseShape();
  BOOST_REQUIRE(capsule.HasBVH());
  BOOST_CHECK(capsule.IsClosed());
  BOOST_CHECK(capsule.IsOrientationConsistent());

  const auto capsuleInside = [&](const double* point) {
    const double axialDistance = std::max(0., std::abs(point[2]) - halfHeight);
    return std::hypot(point[0], point[1], axialDistance) < radius;
  };
  constexpr int samples = 9;
  constexpr double extent = 3.;
  for (int stepX = 0; stepX < samples; ++stepX) {
    for (int stepY = 0; stepY < samples; ++stepY) {
      for (int stepZ = 0; stepZ < samples; ++stepZ) {
        const double point[3] = {-extent + 2. * extent * (stepX + 0.517) / samples,
                                 -extent + 2. * extent * (stepY + 0.263) / samples,
                                 -extent + 2. * extent * (stepZ + 0.741) / samples};
        BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
        {
          const bool bvhInside = capsule.Contains(point);
          BOOST_CHECK_EQUAL(bvhInside, capsule.Contains_Loop(point));
          BOOST_CHECK_EQUAL(bvhInside, capsuleInside(point));
        }
      }
    }
  }

  // a few characteristic capsule points, incl. points exactly on the barrel and cap surfaces
  const double onBarrel[3] = {radius, 0., 0.5};
  const double onCapApex[3] = {0., 0., halfHeight + radius};
  const double aboveApex[3] = {0., 0., halfHeight + radius + 0.01};
  const double onRim[3] = {radius, 0., halfHeight}; // shared cylinder/sphere rim
  BOOST_CHECK(capsule.Contains(onBarrel));
  BOOST_CHECK(capsule.Contains(onCapApex));
  BOOST_CHECK(!capsule.Contains(aboveApex));
  BOOST_CHECK(capsule.Contains(onRim));

  // exact capacity: cylinder plus a full sphere from the two hemispheres
  checkClose(capsule.Capacity(),
             surf::kPi * radius * radius * 2. * halfHeight + 4. * surf::kPi * radius * radius * radius / 3., 1.e-9);
}

BOOST_AUTO_TEST_CASE(DistanceBVHMatchesLoopOnAllFixtures)
{
  // The BVH distance queries against their all-surfaces oracle, over every fixture family and a
  // dense point x direction sweep. This is the correctness guard that does not depend on any
  // reference shape: it isolates traversal and pruning from the analytic kernels, which the
  // per-shape cases above already validate against ROOT.
  const std::array<std::pair<std::unique_ptr<SurfaceSolid>, double>, 7> fixtures{{
    {makeBoxSolid("loopBox", 1., 2., 3.), 4.},
    {makeTubeSolid("loopTube", 0., 2., 3.), 4.},
    {makeTubeSolid("loopHollowTube", 1., 2., 3.), 4.},
    {makeConeSolid("loopCone", 2., 1., 3.), 4.},
    {makeSphereSolid("loopSphere", 2.5), 3.5},
    {makeTorusSolid("loopTorus", 3., 1.), 4.5},
    {makeCapsuleSolid("loopCapsule", 1., 1.5), 3.},
  }};

  for (const auto& [solid, extent] : fixtures) {
    BOOST_TEST_CONTEXT("fixture = " << solid->GetName())
    {
      BOOST_REQUIRE(solid->HasBVH());
      sweepDistanceAgainstLoop(*solid, extent, 5);
    }
  }
}

BOOST_AUTO_TEST_CASE(DistanceSweepsMatchRootPrimitives)
{
  // Systematic point x direction sweeps against the ROOT primitives, for both the entering and
  // the exiting query. The per-shape cases above check a handful of hand-picked rays; this walks
  // a grid, so it also covers rays that miss, that graze, and that cross a hole.
  constexpr int samples = 5;

  const auto box = makeBoxSolid("sweepBox", 1., 2., 3.);
  TGeoBBox boxReference("sweepBoxReference", 1., 2., 3.);
  sweepDistanceAgainstReference(*box, boxReference, 4., samples);

  const auto tube = makeTubeSolid("sweepTube", 0., 2., 3.);
  TGeoTube tubeReference("sweepTubeReference", 0., 2., 3.);
  sweepDistanceAgainstReference(*tube, tubeReference, 4., samples);

  const auto hollowTube = makeTubeSolid("sweepHollowTube", 1., 2., 3.);
  TGeoTube hollowTubeReference("sweepHollowTubeReference", 1., 2., 3.);
  sweepDistanceAgainstReference(*hollowTube, hollowTubeReference, 4., samples);

  const auto cone = makeConeSolid("sweepCone", 2., 1., 3.);
  TGeoCone coneReference("sweepConeReference", 3., 0., 2., 0., 1.);
  sweepDistanceAgainstReference(*cone, coneReference, 4., samples);

  const auto sphere = makeSphereSolid("sweepSphere", 2.5);
  TGeoSphere sphereReference("sweepSphereReference", 0., 2.5);
  sweepDistanceAgainstReference(*sphere, sphereReference, 3.5, samples);

  // the torus kernel solves a quartic, so it carries more rounding than the quadric shapes
  const auto torus = makeTorusSolid("sweepTorus", 3., 1.);
  TGeoTorus torusReference("sweepTorusReference", 3., 0., 1.);
  sweepDistanceAgainstReference(*torus, torusReference, 4.5, samples, 1.e-6);
}

BOOST_AUTO_TEST_CASE(DistanceHardCases)
{
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;
  const auto box = makeBoxSolid("hardCaseBox", halfX, halfY, halfZ);
  TGeoBBox reference("hardCaseBoxReference", halfX, halfY, halfZ);

  // --- rays through a shared edge and a shared vertex -----------------------------------------
  // Both are seen by more than one patch, so the same crossing is reported several times. Taking
  // the minimum over entering hits is insensitive to that, but the BVH and the loop must still
  // see the same set, which is what the loop cross-check asserts.
  const std::array<std::array<double, 3>, 4> throughFeature{{
    {-5., halfY, 0.},        // straight at the x = -1 / y = +2 edge
    {-5., halfY, halfZ},     // straight at the (-1, +2, +3) vertex
    {0., 0., 0.},            // from the centre out through the +y/+z edge
    {-5., -5., -5.},         // body diagonal through the (-1,-2,-3) vertex
  }};
  const std::array<std::array<double, 3>, 4> throughFeatureDirection{{
    {1., 0., 0.},
    {1., 0., 0.},
    unitDirection(0., 1., 1.5),
    unitDirection(1., 2., 3.),
  }};
  for (size_t index = 0; index < throughFeature.size(); ++index) {
    checkDistanceAgainstLoop(*box, throughFeature[index], throughFeatureDirection[index]);
  }

  // --- grazing / tangent rays ------------------------------------------------------------------
  // A ray in the plane of a face never enters: every hit it can report is tangential, and a
  // tangential hit is not a crossing. Both queries must agree with the loop and find nothing.
  const std::array<std::array<double, 3>, 3> grazing{{
    {-5., halfY, 0.},   // in the plane of the y = +2 face
    {halfX, -5., 0.},   // in the plane of the x = +1 face
    {0., 0., halfZ},    // in the plane of the z = +3 face
  }};
  const std::array<std::array<double, 3>, 3> grazingDirection{{
    {1., 0., 0.},
    {0., 1., 0.},
    unitDirection(1., 1., 0.),
  }};
  for (size_t index = 0; index < grazing.size(); ++index) {
    checkDistanceAgainstLoop(*box, grazing[index], grazingDirection[index]);
  }
  // a cylinder tangent ray: the double root must not be reported as two crossings
  const auto tube = makeTubeSolid("hardCaseTube", 0., 2., 3.);
  checkDistanceAgainstLoop(*tube, {-5., 2., 0.}, {1., 0., 0.});
  checkDistanceAgainstLoop(*tube, {-5., 2. - 1.e-7, 0.}, {1., 0., 0.}); // just inside tangency
  checkDistanceAgainstLoop(*tube, {-5., 2. + 1.e-7, 0.}, {1., 0., 0.}); // just outside tangency

  // --- rays starting exactly on a surface -------------------------------------------------------
  // The on-surface convention (a crossing at t = 0 is below kRayTolerance and is not reported) is
  // inherited from the analytic kernels; what matters here is that the BVH reproduces it exactly.
  const std::array<std::array<double, 3>, 4> onSurface{{
    {halfX, 0., 0.},         // on a face
    {-halfX, 0.5, -1.},      // on the opposite face
    {halfX, halfY, 0.},      // on an edge
    {halfX, halfY, halfZ},   // on a vertex
  }};
  for (const auto& point : onSurface) {
    for (const auto& direction : probeDirections()) {
      checkDistanceAgainstLoop(*box, point, direction);
    }
  }
  // just off the surface the answers must be the ordinary ones: entering after ~1e-6 from
  // outside, exiting after the full traversal from inside
  const std::array<double, 3> justOutside{halfX + 1.e-6, 0., 0.};
  const std::array<double, 3> justInside{halfX - 1.e-6, 0., 0.};
  const std::array<double, 3> inward{-1., 0., 0.};
  checkClose(box->DistFromOutside(justOutside.data(), inward.data(), 3), 1.e-6, 1.e-12);
  checkClose(box->DistFromInside(justInside.data(), inward.data(), 3), 2. * halfX - 1.e-6, 1.e-12);
  checkClose(box->DistFromOutside(justOutside.data(), inward.data(), 3),
             reference.DistFromOutside(justOutside.data(), inward.data(), 3), 1.e-12);

  // --- stepmax ----------------------------------------------------------------------------------
  const std::array<double, 3> farOutside{-5., 0., 0.};
  const std::array<double, 3> alongX{1., 0., 0.};
  const double entryDistance = box->DistFromOutside(farOutside.data(), alongX.data(), 3);
  checkClose(entryDistance, 4.);
  // a hit beyond stepmax must not be reported ...
  BOOST_CHECK_EQUAL(box->DistFromOutside(farOutside.data(), alongX.data(), 3, entryDistance * 0.5),
                    TGeoShape::Big());
  // ... including when it lies only just beyond, and the cheap bounding-box reject must agree
  BOOST_CHECK_EQUAL(box->DistFromOutside(farOutside.data(), alongX.data(), 3, entryDistance - 1.e-3),
                    TGeoShape::Big());
  // ... while a stepmax past the hit changes nothing
  checkClose(box->DistFromOutside(farOutside.data(), alongX.data(), 3, entryDistance + 1.e-3), entryDistance);
  checkClose(box->DistFromOutside(farOutside.data(), alongX.data(), 3, 100.), entryDistance);
  // the same for the exiting query
  const std::array<double, 3> center{0., 0., 0.};
  const double exitDistance = box->DistFromInside(center.data(), alongX.data(), 3);
  checkClose(exitDistance, halfX);
  BOOST_CHECK_EQUAL(box->DistFromInside(center.data(), alongX.data(), 3, exitDistance * 0.5), TGeoShape::Big());
  checkClose(box->DistFromInside(center.data(), alongX.data(), 3, exitDistance * 2.), exitDistance);
  // and the loop must honour stepmax identically, at and around the hit
  for (const double stepmax : {entryDistance * 0.5, entryDistance - 1.e-9, entryDistance, entryDistance + 1.e-9,
                               entryDistance * 2.}) {
    checkDistanceAgainstLoop(*box, farOutside, alongX, stepmax);
    checkDistanceAgainstLoop(*box, center, alongX, stepmax);
  }

  // --- a ray that cannot reach the solid at all --------------------------------------------------
  const std::array<double, 3> wayOff{-1000., 0., 0.};
  BOOST_CHECK_EQUAL(box->DistFromOutside(wayOff.data(), alongX.data(), 3, 10.), TGeoShape::Big());
  checkClose(box->DistFromOutside(wayOff.data(), alongX.data(), 3), 999.);
}

BOOST_AUTO_TEST_CASE(RayTMaxPruningIsOptimizationOnly)
{
  // A row of well-separated boxes: a ray along the row enters the first one, after which every
  // node behind it is beyond the tightened bound and must not be visited. Turning the tightening
  // off must cost candidates without changing a single answer.
  constexpr int boxCount = 8;
  constexpr double half = 0.5;
  constexpr double spacing = 4.;

  SurfaceSolid row("prunedRow");
  for (int boxIndex = 0; boxIndex < boxCount; ++boxIndex) {
    addBoxSurfaces(row, half, half, half, {boxIndex * spacing, 0., 0.});
  }
  row.CloseShape();
  BOOST_REQUIRE(row.HasBVH());
  BOOST_CHECK(row.IsClosed());
  BOOST_CHECK_EQUAL(row.GetNsurfaces(), 6 * boxCount);

  const std::array<double, 3> beforeRow{-5., 0., 0.};
  const std::array<double, 3> alongRow{1., 0., 0.};

  BOOST_CHECK(SurfaceSolid::GetRayTMaxPruning()); // on by default

  SurfaceSolid::ResetRayCandidateCounter();
  const double prunedDistance = row.DistFromOutside(beforeRow.data(), alongRow.data(), 3);
  const long long prunedCandidates = SurfaceSolid::GetRayCandidateCount();

  SurfaceSolid::SetRayTMaxPruning(false);
  SurfaceSolid::ResetRayCandidateCounter();
  const double unprunedDistance = row.DistFromOutside(beforeRow.data(), alongRow.data(), 3);
  const long long unprunedCandidates = SurfaceSolid::GetRayCandidateCount();
  SurfaceSolid::SetRayTMaxPruning(true);

  // same answer, and it is the entry face of the first box
  BOOST_CHECK_EQUAL(prunedDistance, unprunedDistance);
  checkClose(prunedDistance, 5. - half);
  // ... reached after strictly less work
  BOOST_CHECK_GT(prunedCandidates, 0);
  BOOST_CHECK_LT(prunedCandidates, unprunedCandidates);

  // the answers stay identical over a full sweep, which is the property that lets the benchmark
  // treat the switch as a pure cost knob
  sweepDistanceAgainstLoop(row, 1.2 * boxCount * spacing / 2., 4);

  // the counter is not touched by the _Loop variants, which visit everything by construction
  SurfaceSolid::ResetRayCandidateCounter();
  row.DistFromOutside_Loop(beforeRow.data(), alongRow.data());
  BOOST_CHECK_EQUAL(SurfaceSolid::GetRayCandidateCount(), 0);
}

BOOST_AUTO_TEST_CASE(RayTMaxPruningKeepsNearTies)
{
  // Two entering candidates a controlled hair apart, one of them behind a very loose bounding
  // box: the geometry in which a mis-set tmax would do its damage.
  //
  // Why this shape of test. A node is culled when the ray *enters its box* beyond tmax, and a box
  // is always entered no later than the patch inside it is hit. A candidate nearer than the
  // current best therefore has a box entered earlier than the current best's hit, and survives
  // any bound at or above that hit -- which is why the implementation's bound (the best hit,
  // rounded up, plus the box inflation) can be argued safe rather than merely measured safe. The
  // narrow window that is left needs a loose box visited first and a tight one entered between
  // the loose patch's box and its hit, so that is what this fixture builds.
  //
  // Fixture: a sphere hit by a near-limb ray far behind where its bounding box starts, plus a
  // small flat patch just in front of that hit, swept over several decades of separation. It is
  // deliberately not a closed manifold (the patch clips into the sphere) and is closed with the
  // diagnostics off: it exists to place the two candidates, not to model a solid. That is
  // legitimate because the oracle is DistFromOutside_Loop, which minimises over the same hits.
  //
  // Scope, honestly: mutation-testing this suite showed that a bound scaled by 0.5 is caught
  // loudly by the sweeps above, while one scaled by 0.999 is caught by neither them nor this
  // case -- with so few primitives, both leaves are box-tested in the same inner-node visit,
  // before any leaf callback has run and tightened anything. So this pins the near-tie geometry
  // and the pruning-on == pruning-off == loop identity; the guarantee against a subtly tight
  // bound rests on the argument above, not on this test.
  constexpr double radius = 2.;
  constexpr double rayOffsetY = 1.9; // near the limb: box entered at x = -2, surface at x = -0.62
  const double sphereHitX = -std::sqrt(radius * radius - rayOffsetY * rayOffsetY);
  const std::array<double, 3> rayOrigin{-10., rayOffsetY, 0.};
  const std::array<double, 3> alongX{1., 0., 0.};
  const double sphereDistance = sphereHitX - rayOrigin[0];

  // relative offsets spanning several decades below the sphere hit, so any tmax that is too
  // tight by anything in that range is caught by at least one of them regardless of how the
  // builder happens to lay out the tree
  for (const double relativeOffset : {1.e-5, 3.e-5, 1.e-4, 3.e-4, 1.e-3, 3.e-3, 1.e-2}) {
    const double patchX = sphereHitX - relativeOffset * sphereDistance;
    BOOST_TEST_CONTEXT("relativeOffset = " << relativeOffset << " patchX = " << patchX)
    {
      SurfaceSolid solid("nearTie");
      BOOST_REQUIRE(solid.AddSphericalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius));
      // axisU x axisV = z x y = -x: the patch faces the incoming ray, so crossing it enters
      BOOST_REQUIRE(solid.AddPlanarSurface({patchX, rayOffsetY - 0.1, -0.1}, {0., 0., 1.}, {0., 1., 0.},
                                           rectangleWire(0.2, 0.2)));
      solid.CloseShape(false);
      BOOST_REQUIRE(solid.HasBVH());

      // the flat patch is the nearer entering crossing, by construction
      const double expected = patchX - rayOrigin[0];
      for (const bool pruning : {true, false}) {
        SurfaceSolid::SetRayTMaxPruning(pruning);
        BOOST_TEST_CONTEXT("pruning = " << pruning)
        {
          const double distance = solid.DistFromOutside(rayOrigin.data(), alongX.data(), 3);
          checkClose(distance, expected, 1.e-9);
          BOOST_CHECK_EQUAL(distance, solid.DistFromOutside_Loop(rayOrigin.data(), alongX.data()));
        }
      }
      SurfaceSolid::SetRayTMaxPruning(true);
    }
  }
}

BOOST_AUTO_TEST_CASE(CurvedPlanarStadiumPrism)
{
  // A stadium (rectangle with two semicircular ends) extruded along z: the two end caps are
  // planar faces with mixed line+arc wires - the general curved-planar case a disk cannot
  // express. Straight sides are flat rectangles; the round ends are half-cylinders.
  constexpr double halfLen = 3.;    // straight half-length along x
  constexpr double radius = 2.;     // corner radius and half-width along y
  constexpr double halfHeight = 4.; // half-height along z

  // Stadium cross-section boundary in the cap's local (u=x, v=y) frame, CCW: bottom line,
  // right semicircle, top line, left semicircle.
  const std::vector<BoundaryCurve> stadiumWire{
    BoundaryCurve::makeLine({-halfLen, -radius}, {halfLen, -radius}),
    BoundaryCurve::makeArc({halfLen, 0.}, radius, -surf::kHalfPi, surf::kHalfPi),
    BoundaryCurve::makeLine({halfLen, radius}, {-halfLen, radius}),
    BoundaryCurve::makeArc({-halfLen, 0.}, radius, surf::kHalfPi, 3. * surf::kHalfPi)};

  SurfaceSolid solid("stadiumPrism");
  // caps (outward +z / -z: the bottom cap flips axisV)
  BOOST_REQUIRE(solid.AddCurvedPlanarSurface({0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, stadiumWire));
  BOOST_REQUIRE(solid.AddCurvedPlanarSurface({0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, stadiumWire));
  // flat side walls at y = +/- radius (outward +/- y)
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfLen, radius, -halfHeight}, {0., 0., 1.}, {1., 0., 0.},
                                       rectangleWire(2. * halfHeight, 2. * halfLen)));
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfLen, -radius, -halfHeight}, {1., 0., 0.}, {0., 0., 1.},
                                       rectangleWire(2. * halfLen, 2. * halfHeight)));
  // round ends as half-cylinders (outer walls)
  BOOST_REQUIRE(solid.AddCylindricalSurface({halfLen, 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                            halfHeight, -surf::kHalfPi, surf::kPi));
  BOOST_REQUIRE(solid.AddCylindricalSurface({-halfLen, 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                            halfHeight, surf::kHalfPi, surf::kPi));
  solid.CloseShape();

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  // exact capacity: (rectangle 2L x 2R + full circle pi R^2) x height 2H
  checkClose(solid.Capacity(), (4. * halfLen * radius + surf::kPi * radius * radius) * 2. * halfHeight, 1.e-6);

  const auto stadiumInside = [&](const double* point) {
    if (std::abs(point[2]) > halfHeight) {
      return false;
    }
    const double ax = std::abs(point[0]);
    const double dx = ax > halfLen ? ax - halfLen : 0.;
    return dx * dx + point[1] * point[1] <= radius * radius;
  };
  // deterministic grid spanning well beyond the solid on every axis
  constexpr int samples = 21;
  const double extentX = 6., extentY = 3.5, extentZ = 5.;
  for (int stepX = 0; stepX < samples; ++stepX) {
    for (int stepY = 0; stepY < samples; ++stepY) {
      for (int stepZ = 0; stepZ < samples; ++stepZ) {
        const double point[3] = {-extentX + 2. * extentX * (stepX + 0.517) / samples,
                                 -extentY + 2. * extentY * (stepY + 0.263) / samples,
                                 -extentZ + 2. * extentZ * (stepZ + 0.741) / samples};
        BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
        {
          BOOST_CHECK_EQUAL(solid.Contains(point), stadiumInside(point));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(WireTrimmedCylinderMatchesTube)
{
  constexpr double radius = 2.;
  constexpr double halfHeight = 3.;

  // lateral wall via the wire-trim overload: the trim is the full parametric rectangle
  // phi in [0, 2pi] x h in [-hh, hh] expressed as four line edges, which must behave exactly like
  // the scalar rectangle path (equivalence check).
  SurfaceSolid solid("wireTrimmedCylinder");
  BOOST_REQUIRE(solid.AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                            halfHeight, 0., surf::kTwoPi, false,
                                            paramRectWire(0., surf::kTwoPi, -halfHeight, halfHeight)));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radius));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radius));
  solid.CloseShape();

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoTube reference("wireTrimTube", 0., radius, halfHeight);
  compareContainsGrid(solid, reference, 4., 9);
  compareDistance(solid, reference, {5., 0., 0.}, {-1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 5.}, {0., 0., -1.});
  compareDistance(solid, reference, {-4., -1., -2.}, unitDirection(1., 0.3, 0.5));
  compareDistance(solid, reference, {0., 0., 0.}, unitDirection(1., 1., 1.));
  compareDistance(solid, reference, {5., 2.5, 0.}, {-1., 0., 0.}); // grazing miss

  // capacity is numerically integrated for a wire trim; the wall integrand is constant here so it
  // stays accurate, but compare with a relaxed tolerance to reflect the quadrature
  checkClose(solid.Capacity(), reference.Capacity(), 1.e-6);

  double normal[3] = {0., 0., 0.};
  const double sidePoint[3] = {radius, 0., 1.};
  const double alongX[3] = {1., 0., 0.};
  solid.ComputeNormal(sidePoint, alongX, normal);
  checkClose(normal[0], 1.);
  checkClose(normal[1], 0.);
  checkClose(normal[2], 0.);
}

BOOST_AUTO_TEST_CASE(WireTrimmedConeMatchesCone)
{
  constexpr double halfHeight = 3.;
  constexpr double radiusAtBottom = 2.;
  constexpr double radiusAtTop = 1.;

  SurfaceSolid solid("wireTrimmedCone");
  BOOST_REQUIRE(solid.AddConicalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radiusAtBottom, radiusAtTop,
                                        -halfHeight, halfHeight, 0., surf::kTwoPi, false,
                                        paramRectWire(0., surf::kTwoPi, -halfHeight, halfHeight)));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radiusAtTop));
  BOOST_REQUIRE(addDiskSurface(solid, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radiusAtBottom));
  solid.CloseShape();

  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoCone reference("wireTrimCone", halfHeight, 0., radiusAtBottom, 0., radiusAtTop);
  compareContainsGrid(solid, reference, 3.5, 9);
  compareDistance(solid, reference, {5., 0., 0.}, {-1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 0.}, {1., 0., 0.});
  compareDistance(solid, reference, {-4., 0.2, -2.}, unitDirection(1., 0.05, 0.3));
  // varying integrand -> looser quadrature tolerance
  checkClose(solid.Capacity(), reference.Capacity(), 1.e-3);
}

BOOST_AUTO_TEST_CASE(WireTrimmedQuadricKernels)
{
  using surf::Curve2D;
  using surf::Vec3;
  std::string error;

  const auto onCylinder = [](double phi, double height) {
    return Vec3{2. * std::cos(phi), 2. * std::sin(phi), height};
  };

  // (1) cylinder wall with a rectangular window (hole) in (phi, h): phi in [2.0, 2.5], h in [-1, 1]
  surf::CylindricalBoundedSurface windowed;
  const std::vector<Curve2D> outer{Curve2D::makeLine({0., -3.}, {surf::kTwoPi, -3.}),
                                   Curve2D::makeLine({surf::kTwoPi, -3.}, {surf::kTwoPi, 3.}),
                                   Curve2D::makeLine({surf::kTwoPi, 3.}, {0., 3.}),
                                   Curve2D::makeLine({0., 3.}, {0., -3.})};
  const std::vector<Curve2D> hole{Curve2D::makeLine({2.0, -1.}, {2.5, -1.}), Curve2D::makeLine({2.5, -1.}, {2.5, 1.}),
                                  Curve2D::makeLine({2.5, 1.}, {2.0, 1.}), Curve2D::makeLine({2.0, 1.}, {2.0, -1.})};
  BOOST_REQUIRE(windowed.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -3., 3., 0., surf::kTwoPi, false,
                                    outer, {hole}, error));
  BOOST_CHECK(windowed.containsPointOnSurface(onCylinder(0.5, 0.)));    // material
  BOOST_CHECK(windowed.containsPointOnSurface(onCylinder(2.25, 2.5)));  // material above the window
  BOOST_CHECK(!windowed.containsPointOnSurface(onCylinder(2.25, 0.)));  // inside the window
  BOOST_CHECK(windowed.containsPointOnSurface(onCylinder(2.25, 1.)));   // on the window edge (boundary)

  // a radial ray into the window is filtered out; a radial ray into material registers one hit
  std::vector<surf::RayHit> hits;
  windowed.appendIntersections({0., 0., 0.}, {std::cos(2.25), std::sin(2.25), 0.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());
  hits.clear();
  windowed.appendIntersections({0., 0., 0.}, {std::cos(0.5), std::sin(0.5), 0.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 1u);
  checkClose(hits.front().distance, 2.);

  // (2) arc trim: a parametric circle (disk in (phi, h)) centred at (pi, 0), radius 0.5
  surf::CylindricalBoundedSurface arcTrim;
  const std::vector<Curve2D> arcOuter{Curve2D::makeCircle({surf::kPi, 0.}, 0.5)};
  BOOST_REQUIRE(arcTrim.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -1., 1., 0., surf::kTwoPi, false,
                                   arcOuter, {}, error));
  BOOST_CHECK(arcTrim.containsPointOnSurface(onCylinder(surf::kPi, 0.)));         // centre of the disk
  BOOST_CHECK(!arcTrim.containsPointOnSurface(onCylinder(surf::kPi, 0.6)));       // outside in h
  BOOST_CHECK(!arcTrim.containsPointOnSurface(onCylinder(surf::kPi + 0.6, 0.)));  // outside in phi
  BOOST_CHECK_GT(std::abs(arcTrim.capacityContribution()), 0.);

  // (3) sphere section reproduced as a (phi, theta) rectangle wire must match the scalar section
  const auto onSphere = [](double theta, double phi) {
    return Vec3{2. * std::sin(theta) * std::cos(phi), 2. * std::sin(theta) * std::sin(phi), 2. * std::cos(theta)};
  };
  surf::SphericalBoundedSurface sphereWire;
  BOOST_REQUIRE(sphereWire.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., surf::kHalfPi / 2., surf::kHalfPi,
                                      0., surf::kHalfPi, false,
                                      paramRectWireCurves(0., surf::kHalfPi, surf::kHalfPi / 2., surf::kHalfPi), {},
                                      error));
  surf::SphericalBoundedSurface sphereScalar;
  BOOST_REQUIRE(sphereScalar.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., surf::kHalfPi / 2.,
                                        surf::kHalfPi, 0., surf::kHalfPi, false, error));
  BOOST_CHECK(sphereWire.containsPointOnSurface(onSphere(surf::kPi / 3., surf::kPi / 4.)));      // inside the section
  BOOST_CHECK(!sphereWire.containsPointOnSurface(onSphere(surf::kPi / 6., surf::kPi / 4.)));     // theta too small
  BOOST_CHECK(!sphereWire.containsPointOnSurface(onSphere(surf::kPi / 3., 3. * surf::kPi / 4.))); // phi outside
  checkClose(sphereWire.capacityContribution(), sphereScalar.capacityContribution(), 1.e-3);

  // (4) a trim spanning more than a full turn in phi is rejected
  surf::CylindricalBoundedSurface tooWide;
  const std::vector<Curve2D> wideOuter{Curve2D::makeLine({0., -1.}, {7., -1.}), Curve2D::makeLine({7., -1.}, {7., 1.}),
                                       Curve2D::makeLine({7., 1.}, {0., 1.}), Curve2D::makeLine({0., 1.}, {0., -1.})};
  BOOST_CHECK(!tooWide.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -1., 1., 0., surf::kTwoPi, false,
                                  wideOuter, {}, error));
}

namespace
{
// Helpers writing the surface sidecar binary format (version 1) documented in
// scripts/geometry/BVHSurfaceSolid.md. Kept independent of the loader implementation so the
// test is a true round-trip through the documented byte layout.
void appendU32(std::vector<char>& bytes, uint32_t value)
{
  const char* raw = reinterpret_cast<const char*>(&value);
  bytes.insert(bytes.end(), raw, raw + sizeof(value));
}

void appendDoubles(std::vector<char>& bytes, std::initializer_list<double> values)
{
  for (const double value : values) {
    const char* raw = reinterpret_cast<const char*>(&value);
    bytes.insert(bytes.end(), raw, raw + sizeof(value));
  }
}

void appendSidecarHeader(std::vector<char>& bytes, uint32_t nSurfaces)
{
  bytes.insert(bytes.end(), {'O', '2', 'S', 'S'});
  appendU32(bytes, 1); // version
  appendU32(bytes, nSurfaces);
  appendU32(bytes, 0); // reserved
}

// plane record (type 1) with a single rectangular outer wire of four line-segment edges
void appendPlaneRecord(std::vector<char>& bytes, const FaceFrame& frame)
{
  appendU32(bytes, 1); // surfaceType plane
  appendU32(bytes, 0); // flags
  appendU32(bytes, 9); // nParams
  appendDoubles(bytes, {frame.origin[0], frame.origin[1], frame.origin[2], frame.axisU[0], frame.axisU[1],
                        frame.axisU[2], frame.axisV[0], frame.axisV[1], frame.axisV[2]});
  appendU32(bytes, 1); // nWires
  appendU32(bytes, 0); // wireRole outer
  appendU32(bytes, 4); // nEdges
  const double extentU = frame.extentU;
  const double extentV = frame.extentV;
  const std::array<std::array<double, 4>, 4> edges{{{0., 0., extentU, 0.},
                                                    {extentU, 0., extentU, extentV},
                                                    {extentU, extentV, 0., extentV},
                                                    {0., extentV, 0., 0.}}};
  for (const auto& edge : edges) {
    appendU32(bytes, 0); // curveType line
    appendU32(bytes, 4); // nCurveParams
    appendDoubles(bytes, {edge[0], edge[1], edge[2], edge[3]});
  }
}

// plane record (type 1) for a disk/annulus cap: one full-circle outer arc wire, plus a
// clockwise inner arc wire when holeRadius > 0. Exercises the arc-wire reader path.
void appendDiskPlaneRecord(std::vector<char>& bytes, const Point3D& center, const Point3D& axisU,
                           const Point3D& axisV, double radius, double holeRadius = 0.)
{
  appendU32(bytes, 1); // surfaceType plane
  appendU32(bytes, 0); // flags
  appendU32(bytes, 9); // nParams
  appendDoubles(bytes, {center[0], center[1], center[2], axisU[0], axisU[1], axisU[2], axisV[0], axisV[1], axisV[2]});
  const uint32_t nWires = holeRadius > 0. ? 2u : 1u;
  appendU32(bytes, nWires);
  appendU32(bytes, 0); // outer wire role
  appendU32(bytes, 1); // one edge
  appendU32(bytes, 1); // curveType arc
  appendU32(bytes, 5); // nCurveParams
  appendDoubles(bytes, {0., 0., radius, 0., 2. * surf::kPi}); // cu cv radius phiStart phiSweep (CCW full circle)
  if (holeRadius > 0.) {
    appendU32(bytes, 1); // inner wire role
    appendU32(bytes, 1);
    appendU32(bytes, 1); // arc
    appendU32(bytes, 5);
    appendDoubles(bytes, {0., 0., holeRadius, 0., -2. * surf::kPi}); // clockwise hole
  }
}

std::filesystem::path writeSidecarFile(const std::string& name, const std::vector<char>& bytes)
{
  const auto path = std::filesystem::temp_directory_path() / name;
  std::ofstream out(path, std::ios::binary);
  out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  BOOST_REQUIRE(out.good());
  return path;
}
} // namespace

BOOST_AUTO_TEST_CASE(SurfaceSidecarRoundTrip)
{
  // planar box: six plane records with polygon wires, loaded and compared against TGeoBBox
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;

  std::vector<char> boxBytes;
  appendSidecarHeader(boxBytes, 6);
  for (int faceIndex = 0; faceIndex < 6; ++faceIndex) {
    appendPlaneRecord(boxBytes, boxFaceFrame(faceIndex, halfX, halfY, halfZ));
  }
  const auto boxPath = writeSidecarFile("o2_sidecar_roundtrip_box.bin", boxBytes);

  SurfaceSolid box("sidecarBox");
  BOOST_REQUIRE(o2::base::LoadSurfaceSolid(boxPath.string(), box));
  std::filesystem::remove(boxPath);
  BOOST_CHECK_EQUAL(box.GetNsurfaces(), 6);
  box.CloseShape();
  BOOST_CHECK(box.IsClosed());
  BOOST_CHECK(box.IsOrientationConsistent());

  TGeoBBox referenceBox("referenceBox", halfX, halfY, halfZ);
  compareContainsGrid(box, referenceBox, 4., 7);
  compareDistance(box, referenceBox, {5., 0.5, 0.5}, {-1., 0., 0.});
  compareDistance(box, referenceBox, {0., 0., 0.}, unitDirection(1., 1., 1.));
  checkClose(box.Capacity(), referenceBox.Capacity(), 1.e-9);

  // quadric + arc-wire caps: closed cylinder (lateral wall + two disk caps) against TGeoTube
  constexpr double radius = 2.;
  constexpr double halfHeight = 3.;

  std::vector<char> tubeBytes;
  appendSidecarHeader(tubeBytes, 3);
  appendU32(tubeBytes, 2);  // surfaceType cylinder
  appendU32(tubeBytes, 0);  // flags (outer wall)
  appendU32(tubeBytes, 14); // nParams
  appendDoubles(tubeBytes, {0., 0., 0., 0., 0., 1., 1., 0., 0., radius, -halfHeight, halfHeight, 0., 2. * surf::kPi});
  appendU32(tubeBytes, 0); // nWires
  // caps as arc-wire plane records: outward normal is axisU x axisV, so the bottom cap flips axisV
  appendDiskPlaneRecord(tubeBytes, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radius);
  appendDiskPlaneRecord(tubeBytes, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radius);
  const auto tubePath = writeSidecarFile("o2_sidecar_roundtrip_tube.bin", tubeBytes);

  SurfaceSolid tube("sidecarTube");
  BOOST_REQUIRE(o2::base::LoadSurfaceSolid(tubePath.string(), tube));
  std::filesystem::remove(tubePath);
  BOOST_CHECK_EQUAL(tube.GetNsurfaces(), 3);
  tube.CloseShape();
  BOOST_CHECK(tube.IsClosed());
  BOOST_CHECK(tube.IsOrientationConsistent());

  TGeoTube referenceTube("referenceTube", 0., radius, halfHeight);
  compareContainsGrid(tube, referenceTube, 4., 7);
  compareDistance(tube, referenceTube, {5., 0.5, 1.}, {-1., 0., 0.});
  compareDistance(tube, referenceTube, {0., 0., 0.}, unitDirection(1., 1., 1.));
  checkClose(tube.Capacity(), referenceTube.Capacity(), 1.e-9);

  // malformed input must be rejected without loading surfaces
  const auto badPath = writeSidecarFile("o2_sidecar_bad_magic.bin", {'X', 'X', 'X', 'X', 0, 0, 0, 0});
  SurfaceSolid bad("sidecarBad");
  BOOST_CHECK(!o2::base::LoadSurfaceSolid(badPath.string(), bad));
  std::filesystem::remove(badPath);
  BOOST_CHECK_EQUAL(bad.GetNsurfaces(), 0);
  BOOST_CHECK(!o2::base::LoadSurfaceSolid("/nonexistent/o2_sidecar_missing.bin", bad));

  // truncated file: valid header announcing a surface that never follows
  std::vector<char> truncatedBytes;
  appendSidecarHeader(truncatedBytes, 1);
  const auto truncatedPath = writeSidecarFile("o2_sidecar_truncated.bin", truncatedBytes);
  SurfaceSolid truncated("sidecarTruncated");
  BOOST_CHECK(!o2::base::LoadSurfaceSolid(truncatedPath.string(), truncated));
  std::filesystem::remove(truncatedPath);
}

BOOST_AUTO_TEST_CASE(WireTrimmedSidecarRoundTrip)
{
  // a cylinder record carrying a (line) trim wire block in its (phi, h) domain must load through
  // the wire-taking Add* overload and navigate like the equivalent scalar cylinder
  constexpr double radius = 2.;
  constexpr double halfHeight = 3.;

  std::vector<char> bytes;
  appendSidecarHeader(bytes, 3);
  appendU32(bytes, 2);  // surfaceType cylinder
  appendU32(bytes, 0);  // flags (outer wall)
  appendU32(bytes, 14); // nParams
  appendDoubles(bytes, {0., 0., 0., 0., 0., 1., 1., 0., 0., radius, -halfHeight, halfHeight, 0., 2. * surf::kPi});
  appendU32(bytes, 1); // nWires
  appendU32(bytes, 0); // outer wire role
  appendU32(bytes, 4); // nEdges
  const std::array<std::array<double, 4>, 4> edges{{{0., -halfHeight, 2. * surf::kPi, -halfHeight},
                                                    {2. * surf::kPi, -halfHeight, 2. * surf::kPi, halfHeight},
                                                    {2. * surf::kPi, halfHeight, 0., halfHeight},
                                                    {0., halfHeight, 0., -halfHeight}}};
  for (const auto& edge : edges) {
    appendU32(bytes, 0); // curveType line
    appendU32(bytes, 4); // nCurveParams
    appendDoubles(bytes, {edge[0], edge[1], edge[2], edge[3]});
  }
  appendDiskPlaneRecord(bytes, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radius);
  appendDiskPlaneRecord(bytes, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radius);
  const auto path = writeSidecarFile("o2_sidecar_wiretrim_cylinder.bin", bytes);

  SurfaceSolid solid("sidecarWireCylinder");
  BOOST_REQUIRE(o2::base::LoadSurfaceSolid(path.string(), solid));
  std::filesystem::remove(path);
  BOOST_CHECK_EQUAL(solid.GetNsurfaces(), 3);
  solid.CloseShape();
  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoTube reference("wireTrimSidecarTube", 0., radius, halfHeight);
  compareContainsGrid(solid, reference, 4., 7);
  compareDistance(solid, reference, {5., 0.5, 1.}, {-1., 0., 0.});
  compareDistance(solid, reference, {0., 0., 0.}, unitDirection(1., 1., 1.));
  checkClose(solid.Capacity(), reference.Capacity(), 1.e-6);
}

BOOST_AUTO_TEST_CASE(TorusSidecarRoundTrip)
{
  // a full-torus record (surfaceType 5, 15 params, empty wire block) must load through the
  // scalar AddToroidalSurface path and navigate like TGeoTorus
  constexpr double majorR = 3.;
  constexpr double minorR = 1.;

  std::vector<char> bytes;
  appendSidecarHeader(bytes, 1);
  appendU32(bytes, 5);  // surfaceType torus
  appendU32(bytes, 0);  // flags (outer wall)
  appendU32(bytes, 15); // nParams
  appendDoubles(bytes, {0., 0., 0., 0., 0., 1., 1., 0., 0., majorR, minorR, 0., 2. * surf::kPi, 0., 2. * surf::kPi});
  appendU32(bytes, 0); // nWires (full torus: scalar path)
  const auto path = writeSidecarFile("o2_sidecar_roundtrip_torus.bin", bytes);

  SurfaceSolid solid("sidecarTorus");
  BOOST_REQUIRE(o2::base::LoadSurfaceSolid(path.string(), solid));
  std::filesystem::remove(path);
  BOOST_CHECK_EQUAL(solid.GetNsurfaces(), 1);
  solid.CloseShape();
  BOOST_CHECK(solid.IsClosed());
  BOOST_CHECK(solid.IsOrientationConsistent());

  TGeoTorus reference("sidecarTorusRef", majorR, 0., minorR);
  compareContainsGrid(solid, reference, 4.5, 9);
  checkClose(solid.Capacity(), reference.Capacity(), 1.e-7);
}

namespace
{
// Append a plane record whose rectangular outer wire has its bottom edge as a degree-3 B-spline
// with collinear poles — geometrically identical to the straight edge, so the box still closes,
// but it exercises the whole B-spline sidecar pipeline (curveType 2 reader -> kernel).
void appendBSplineEdgePlaneRecord(std::vector<char>& bytes, const FaceFrame& frame)
{
  appendU32(bytes, 1); // surfaceType plane
  appendU32(bytes, 0); // flags
  appendU32(bytes, 9); // nParams
  appendDoubles(bytes, {frame.origin[0], frame.origin[1], frame.origin[2], frame.axisU[0], frame.axisU[1],
                        frame.axisU[2], frame.axisV[0], frame.axisV[1], frame.axisV[2]});
  appendU32(bytes, 1); // nWires
  appendU32(bytes, 0); // wireRole outer
  appendU32(bytes, 4); // nEdges
  const double extentU = frame.extentU;
  const double extentV = frame.extentV;
  // edge 0: collinear cubic B-spline from (0, 0) to (extentU, 0)
  appendU32(bytes, 2);  // curveType bspline
  appendU32(bytes, 22); // nCurveParams = 2 + 2*4 + 4 + 8
  appendDoubles(bytes, {3., 4.,                                                          // degree, nPoles
                        0., 0., extentU / 3., 0., 2. * extentU / 3., 0., extentU, 0.,    // poles
                        1., 1., 1., 1.,                                                  // weights
                        0., 0., 0., 0., 1., 1., 1., 1.});                                // clamped knots
  const std::array<std::array<double, 4>, 3> lines{
    {{extentU, 0., extentU, extentV}, {extentU, extentV, 0., extentV}, {0., extentV, 0., 0.}}};
  for (const auto& edge : lines) {
    appendU32(bytes, 0); // curveType line
    appendU32(bytes, 4); // nCurveParams
    appendDoubles(bytes, {edge[0], edge[1], edge[2], edge[3]});
  }
}

// A rational quadratic B-spline (NURBS) quarter circle from angle a0 to a0 + pi/2, in the (u, v)
// domain, centred at (cu, cv) with radius r. Four of these form an exact circle.
/// A full circle as ONE closed rational B-spline (the standard 9-pole degree-2 NURBS circle),
/// i.e. a wire whose single edge starts and ends at the same point. This is the shape a CAD
/// kernel writes for a tube-tube intersection curve, and it is structurally different from the
/// same circle spelled as four separate quarter arcs.
surf::Curve2D fullCircleBSpline(double cu, double cv, double r)
{
  const double w = std::sqrt(0.5);
  const std::vector<surf::Vec2> poles{{cu + r, cv},     {cu + r, cv + r}, {cu, cv + r},
                                      {cu - r, cv + r}, {cu - r, cv},     {cu - r, cv - r},
                                      {cu, cv - r},     {cu + r, cv - r}, {cu + r, cv}};
  return surf::Curve2D::makeBSpline(2, poles, {1., w, 1., w, 1., w, 1., w, 1.},
                                    {0., 0., 0., 1., 1., 2., 2., 3., 3., 4., 4., 4.});
}

surf::Curve2D quarterCircleBSpline(double cu, double cv, double r, double a0)
{
  const double a1 = a0 + surf::kHalfPi;
  const double aMid = 0.5 * (a0 + a1);
  const std::vector<surf::Vec2> poles{{cu + r * std::cos(a0), cv + r * std::sin(a0)},
                                      {cu + r * std::sqrt(2.) * std::cos(aMid), cv + r * std::sqrt(2.) * std::sin(aMid)},
                                      {cu + r * std::cos(a1), cv + r * std::sin(a1)}};
  return surf::Curve2D::makeBSpline(2, poles, {1., std::sqrt(0.5), 1.}, {0., 0., 0., 1., 1., 1.});
}
} // namespace

BOOST_AUTO_TEST_CASE(BSplineSidecarRoundTrip)
{
  // a closed box whose first face carries a (collinear) B-spline boundary edge round-trips through
  // the sidecar reader and navigates identically to TGeoBBox
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;

  std::vector<char> bytes;
  appendSidecarHeader(bytes, 6);
  appendBSplineEdgePlaneRecord(bytes, boxFaceFrame(0, halfX, halfY, halfZ));
  for (int faceIndex = 1; faceIndex < 6; ++faceIndex) {
    appendPlaneRecord(bytes, boxFaceFrame(faceIndex, halfX, halfY, halfZ));
  }
  const auto path = writeSidecarFile("o2_sidecar_bspline_box.bin", bytes);

  SurfaceSolid box("sidecarBSplineBox");
  BOOST_REQUIRE(o2::base::LoadSurfaceSolid(path.string(), box));
  std::filesystem::remove(path);
  BOOST_CHECK_EQUAL(box.GetNsurfaces(), 6);
  box.CloseShape();
  BOOST_CHECK(box.IsClosed());
  BOOST_CHECK(box.IsOrientationConsistent());

  TGeoBBox reference("bsplineBoxRef", halfX, halfY, halfZ);
  compareContainsGrid(box, reference, 4., 7);
  compareDistance(box, reference, {5., 0.5, 0.5}, {-1., 0., 0.});
  compareDistance(box, reference, {0., 0., 0.}, unitDirection(1., 1., 1.));
  checkClose(box.Capacity(), reference.Capacity(), 1.e-6);
}

BOOST_AUTO_TEST_CASE(BSplineWindowInCylinderWall)
{
  using surf::Curve2D;
  using surf::Vec3;
  std::string error;

  const auto onCylinder = [](double phi, double height) {
    return Vec3{2. * std::cos(phi), 2. * std::sin(phi), height};
  };

  // an exact circular trim in (phi, h) built from four NURBS quarter arcs must classify identically
  // to the same circle expressed as one exact arc Curve2D — validating the B-spline trim path on a
  // quadric against the closed-form arc path.
  const double centrePhi = surf::kPi;
  const double trimRadius = 0.5;
  surf::CylindricalBoundedSurface bsplineDisk;
  const std::vector<Curve2D> bsplineOuter{quarterCircleBSpline(centrePhi, 0., trimRadius, 0.),
                                          quarterCircleBSpline(centrePhi, 0., trimRadius, surf::kHalfPi),
                                          quarterCircleBSpline(centrePhi, 0., trimRadius, surf::kPi),
                                          quarterCircleBSpline(centrePhi, 0., trimRadius, 3. * surf::kHalfPi)};
  BOOST_REQUIRE(bsplineDisk.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -1., 1., 0., surf::kTwoPi, false,
                                       bsplineOuter, {}, error));
  BOOST_CHECK(!bsplineDisk.capacityIsExact()); // B-spline (wire) trim -> numeric capacity

  surf::CylindricalBoundedSurface arcDisk;
  BOOST_REQUIRE(arcDisk.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -1., 1., 0., surf::kTwoPi, false,
                                   {Curve2D::makeCircle({centrePhi, 0.}, trimRadius)}, {}, error));

  // classification agrees across a grid of the (phi, h) neighbourhood of the trim (skip a thin band
  // around the boundary, where the exact-arc and sampled-B-spline classifications can legitimately
  // differ by the sampling tolerance)
  int compared = 0;
  for (int phiStep = -12; phiStep <= 12; ++phiStep) {
    const double phi = centrePhi + 0.09 * phiStep;
    for (int hStep = -12; hStep <= 12; ++hStep) {
      const double height = 0.09 * hStep;
      // radial distance in (phi, h) from the trim centre; skip the boundary band
      const double distToCentre = std::hypot(phi - centrePhi, height);
      if (std::abs(distToCentre - trimRadius) < 5.e-3) {
        continue;
      }
      const Vec3 point = onCylinder(phi, height);
      BOOST_CHECK_EQUAL(bsplineDisk.containsPointOnSurface(point), arcDisk.containsPointOnSurface(point));
      ++compared;
    }
  }
  BOOST_CHECK_GT(compared, 100);

  // a radial ray into the B-spline window (its centre) registers exactly one wall hit; a ray well
  // outside the window (opposite side of the cylinder) misses the trimmed patch
  std::vector<surf::RayHit> hits;
  bsplineDisk.appendIntersections({0., 0., 0.}, {std::cos(centrePhi), std::sin(centrePhi), 0.}, 0., 1.e30, hits);
  BOOST_REQUIRE_EQUAL(hits.size(), 1u);
  checkClose(hits.front().distance, 2.);
  hits.clear();
  bsplineDisk.appendIntersections({0., 0., 0.}, {std::cos(0.), std::sin(0.), 0.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());
}

// A B-spline wire used as an *inner hole*, and in particular one spelled as a single CLOSED
// B-spline edge. This is what a tube-tube intersection produces and what the converter emits
// constantly: where a boom tube is planted on a fat tube, the fat tube's wall keeps a
// full-rectangle outer wire and carries the intersection curve as one closed B-spline hole.
//
// It regresses a bug that silently deleted such a wire outright. bsplineSampleRecursive used to
// end the recursion when the chord p0->p1 was shorter than the flatness scale; a closed curve has
// p0 == p1 exactly, so a full circle flattened to two coincident points and every polyline-based
// query (winding, closest point, boundary band, display mesh) saw an empty curve. The wire still
// validated and still reported the correct enclosed area, because signedAreaContribution
// integrates the curve by Gauss-Legendre rather than from the polyline -- which is exactly why
// this survived: every check that could have caught it used the analytic path.
//
// Impact: on Bagger/BoomCylinderOuter_0_1_1_9 a point 0.026 cm inside such a hole was reported as
// lying on the face, and a whole face whose outer wire was one closed B-spline did not exist at
// all. `WireTrimmedQuadricKernels` covers a *line* hole and `BSplineWindowInCylinderWall` covers a
// B-spline outer wire built from four *open* quarter arcs, so neither could see it.
// See scripts/geometry/ExactTrimTopology.md.
BOOST_AUTO_TEST_CASE(BSplineHoleInCylinderWall)
{
  using surf::Curve2D;
  using surf::Vec3;
  std::string error;

  const auto onCylinder = [](double phi, double height) {
    return Vec3{2. * std::cos(phi), 2. * std::sin(phi), height};
  };

  // full-sweep outer wire (what the converter writes for an untrimmed cylinder wall) ...
  const std::vector<Curve2D> outer{Curve2D::makeLine({0., -3.}, {surf::kTwoPi, -3.}),
                                   Curve2D::makeLine({surf::kTwoPi, -3.}, {surf::kTwoPi, 3.}),
                                   Curve2D::makeLine({surf::kTwoPi, 3.}, {0., 3.}),
                                   Curve2D::makeLine({0., 3.}, {0., -3.})};
  // ... with a circular hole punched in it, expressed once as four NURBS quarter arcs and once as
  // the equivalent exact arc. The arc form is the oracle: it is the already-trusted path.
  const double centrePhi = surf::kPi;
  const double trimRadius = 0.5;
  const std::vector<Curve2D> bsplineHole{quarterCircleBSpline(centrePhi, 0., trimRadius, 0.),
                                         quarterCircleBSpline(centrePhi, 0., trimRadius, surf::kHalfPi),
                                         quarterCircleBSpline(centrePhi, 0., trimRadius, surf::kPi),
                                         quarterCircleBSpline(centrePhi, 0., trimRadius, 3. * surf::kHalfPi)};

  surf::CylindricalBoundedSurface bsplineHoled;
  BOOST_REQUIRE(bsplineHoled.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -3., 3., 0., surf::kTwoPi, false,
                                        outer, {bsplineHole}, error));
  surf::CylindricalBoundedSurface arcHoled;
  BOOST_REQUIRE(arcHoled.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -3., 3., 0., surf::kTwoPi, false,
                                    outer, {{Curve2D::makeCircle({centrePhi, 0.}, trimRadius)}}, error));

  // the defining property of a hole: its interior is NOT part of the face
  BOOST_CHECK(!arcHoled.containsPointOnSurface(onCylinder(centrePhi, 0.)));     // oracle
  BOOST_CHECK(!bsplineHoled.containsPointOnSurface(onCylinder(centrePhi, 0.))); // the case under test
  // and material well away from the hole still is
  BOOST_CHECK(bsplineHoled.containsPointOnSurface(onCylinder(0.5, 0.)));
  BOOST_CHECK(bsplineHoled.containsPointOnSurface(onCylinder(centrePhi, 2.5)));

  // the two spellings of the same hole must classify identically away from the boundary band
  int compared = 0;
  for (int phiStep = -12; phiStep <= 12; ++phiStep) {
    const double phi = centrePhi + 0.09 * phiStep;
    for (int hStep = -12; hStep <= 12; ++hStep) {
      const double height = 0.09 * hStep;
      if (std::abs(std::hypot(phi - centrePhi, height) - trimRadius) < 5.e-3) {
        continue;
      }
      const Vec3 point = onCylinder(phi, height);
      BOOST_CHECK_EQUAL(bsplineHoled.containsPointOnSurface(point), arcHoled.containsPointOnSurface(point));
      ++compared;
    }
  }
  BOOST_CHECK_GT(compared, 100);

  // a radial ray aimed through the hole must not register a wall hit; one aimed at material must
  std::vector<surf::RayHit> hits;
  bsplineHoled.appendIntersections({0., 0., 0.}, {std::cos(centrePhi), std::sin(centrePhi), 0.}, 0., 1.e30, hits);
  BOOST_CHECK(hits.empty());
  hits.clear();
  bsplineHoled.appendIntersections({0., 0., 0.}, {std::cos(0.5), std::sin(0.5), 0.}, 0., 1.e30, hits);
  BOOST_CHECK_EQUAL(hits.size(), 1u);

  // The same hole as ONE closed B-spline edge rather than four arc segments. This is what the
  // converter actually emits for a tube-tube seam (`_quadric_trim_wire` writes one B-spline per
  // BREP edge, and the intersection curve is a single closed edge).
  surf::CylindricalBoundedSurface singleEdgeHoled;
  BOOST_REQUIRE(singleEdgeHoled.initialize({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, 2., -3., 3., 0., surf::kTwoPi,
                                           false, outer, {{fullCircleBSpline(centrePhi, 0., trimRadius)}}, error));
  BOOST_CHECK(!singleEdgeHoled.containsPointOnSurface(onCylinder(centrePhi, 0.)));
  BOOST_CHECK(singleEdgeHoled.containsPointOnSurface(onCylinder(0.5, 0.)));
  BOOST_CHECK(singleEdgeHoled.containsPointOnSurface(onCylinder(centrePhi, 2.5)));
  for (int phiStep = -12; phiStep <= 12; ++phiStep) {
    const double phi = centrePhi + 0.09 * phiStep;
    for (int hStep = -12; hStep <= 12; ++hStep) {
      const double height = 0.09 * hStep;
      if (std::abs(std::hypot(phi - centrePhi, height) - trimRadius) < 5.e-3) {
        continue;
      }
      const Vec3 point = onCylinder(phi, height);
      BOOST_CHECK_EQUAL(singleEdgeHoled.containsPointOnSurface(point), arcHoled.containsPointOnSurface(point));
    }
  }
}

namespace
{
// The public-API mirror of quarterCircleBSpline, for building a NURBS trim through Add*Surface.
BoundaryCurve quarterCircleBoundaryCurve(double cu, double cv, double r, double a0)
{
  const double a1 = a0 + surf::kHalfPi;
  const double aMid = 0.5 * (a0 + a1);
  const std::vector<Point2D> poles{{cu + r * std::cos(a0), cv + r * std::sin(a0)},
                                   {cu + r * std::sqrt(2.) * std::cos(aMid), cv + r * std::sqrt(2.) * std::sin(aMid)},
                                   {cu + r * std::cos(a1), cv + r * std::sin(a1)}};
  return BoundaryCurve::makeBSpline(2, poles, {1., std::sqrt(0.5), 1.}, {0., 0., 0., 1., 1., 1.});
}

// Assert that two solids are the *same* solid, not merely similar ones: identical closure
// diagnostics and reliability, identical bounding box and capacity, and bit-identical answers
// from all four navigation kernels over the standard probe grid and direction set. This is the
// acceptance criterion for persistence -- a solid that survives a write/read cycle must be
// indistinguishable through the public interface.
void checkSolidsIdentical(const SurfaceSolid& solid, const SurfaceSolid& other, double extent, int samples)
{
  BOOST_CHECK_EQUAL(other.GetNsurfaces(), solid.GetNsurfaces());
  BOOST_CHECK_EQUAL(other.IsDefined(), solid.IsDefined());
  BOOST_CHECK_EQUAL(other.HasBVH(), solid.HasBVH());
  BOOST_CHECK_EQUAL(other.IsClosed(), solid.IsClosed());
  BOOST_CHECK_EQUAL(other.IsOrientationConsistent(), solid.IsOrientationConsistent());
  BOOST_CHECK_EQUAL(static_cast<int>(other.GetNavigationReliability()),
                    static_cast<int>(solid.GetNavigationReliability()));
  BOOST_CHECK_EQUAL(other.GetBoundaryEdgeCount(), solid.GetBoundaryEdgeCount());
  BOOST_CHECK_EQUAL(other.GetNonManifoldEdgeCount(), solid.GetNonManifoldEdgeCount());
  BOOST_CHECK_EQUAL(other.GetReversedEdgeCount(), solid.GetReversedEdgeCount());

  BOOST_CHECK_EQUAL(other.GetDX(), solid.GetDX());
  BOOST_CHECK_EQUAL(other.GetDY(), solid.GetDY());
  BOOST_CHECK_EQUAL(other.GetDZ(), solid.GetDZ());
  for (int dimension = 0; dimension < 3; ++dimension) {
    BOOST_CHECK_EQUAL(other.GetOrigin()[dimension], solid.GetOrigin()[dimension]);
  }
  BOOST_CHECK_EQUAL(other.Capacity(), solid.Capacity());

  for (const auto& point : probeGrid(extent, samples)) {
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      BOOST_CHECK_EQUAL(other.Contains(point.data()), solid.Contains(point.data()));
      BOOST_CHECK_EQUAL(other.Safety(point.data(), solid.Contains(point.data())),
                        solid.Safety(point.data(), solid.Contains(point.data())));
      for (const auto& direction : probeDirections()) {
        BOOST_CHECK_EQUAL(other.DistFromOutside(point.data(), direction.data(), 3),
                          solid.DistFromOutside(point.data(), direction.data(), 3));
        BOOST_CHECK_EQUAL(other.DistFromInside(point.data(), direction.data(), 3),
                          solid.DistFromInside(point.data(), direction.data(), 3));
      }
    }
  }
}

// Write "solid" to a ROOT file and read it back as an independent object.
std::unique_ptr<SurfaceSolid> writeAndReadBack(const SurfaceSolid& solid)
{
  const auto path = std::filesystem::temp_directory_path() /
                    (std::string("o2_bvhsurfacesolid_persist_") + solid.GetName() + ".root");
  {
    TFile file(path.string().c_str(), "RECREATE");
    BOOST_REQUIRE(!file.IsZombie());
    // WriteObject takes a non-const pointer; the call does not modify the solid.
    file.WriteObject(const_cast<SurfaceSolid*>(&solid), "solid");
  }
  std::unique_ptr<SurfaceSolid> restored;
  {
    TFile file(path.string().c_str(), "READ");
    BOOST_REQUIRE(!file.IsZombie());
    restored.reset(file.Get<SurfaceSolid>("solid"));
  }
  std::filesystem::remove(path);
  return restored;
}
} // namespace

// ROOT persistence round trip. The kernel objects behind the solid (BoundedSurface, the BVH, the
// display mesh) are all *derived* state; what has to survive a write/read cycle is the sequence of
// Add*Surface calls the solid was built from, after which CloseShape() reconstructs the rest.
//
// It regresses CodeReview_Fable.md S1, where nothing at all was streamed: fImpl was transient, so
// a read-back solid came back with zero surfaces, CloseShape(false) then zeroed the streamed
// bounding box, and an *empty* ClosureReport defaults to closed/consistent -- so the husk reported
// NavigationReliability::Reliable and answered "outside" everywhere with full confidence. Any
// TGeoManager::Export/Import of a geometry containing one of these solids silently replaced it by
// an authoritatively-reliable empty point.
BOOST_AUTO_TEST_CASE(PersistenceRoundTrip)
{
  // every surface family and both trim flavours (scalar range and wire trim, the latter with
  // line, arc and B-spline curves) must survive, so each record field is exercised
  const auto box = makeBoxSolid("persistBox", 1., 2., 3.);
  const auto tube = makeTubeSolid("persistTube", 1., 2., 3.); // inner wall + annular arc-wire caps
  const auto cone = makeConeSolid("persistCone", 2., 1., 3.);
  const auto sphere = makeSphereSolid("persistSphere", 2.);
  const auto torus = makeTorusSolid("persistTorus", 3., 1.);
  const auto capsule = makeCapsuleSolid("persistCapsule", 2., 3.);

  for (const auto* solid : {box.get(), tube.get(), cone.get(), sphere.get(), torus.get(), capsule.get()}) {
    BOOST_TEST_CONTEXT("solid = " << solid->GetName())
    {
      const auto restored = writeAndReadBack(*solid);
      BOOST_REQUIRE(restored != nullptr);
      checkSolidsIdentical(*solid, *restored, 4.5, 5);
    }
  }

  // a wire-trimmed cylinder whose window is a NURBS loop: the B-spline degree, poles, weights and
  // knots all have to make the round trip, and the trimmed overload has to be the one replayed
  SurfaceSolid trimmed("persistTrimmed");
  constexpr double radius = 2.;
  constexpr double halfHeight = 3.;
  const std::vector<BoundaryCurve> window{quarterCircleBoundaryCurve(surf::kPi, 0., 0.5, 0.),
                                          quarterCircleBoundaryCurve(surf::kPi, 0., 0.5, surf::kHalfPi),
                                          quarterCircleBoundaryCurve(surf::kPi, 0., 0.5, surf::kPi),
                                          quarterCircleBoundaryCurve(surf::kPi, 0., 0.5, 3. * surf::kHalfPi)};
  BOOST_REQUIRE(trimmed.AddCylindricalSurface({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}, radius, -halfHeight,
                                              halfHeight, 0., surf::kTwoPi, false, window));
  BOOST_REQUIRE(addDiskSurface(trimmed, {0., 0., halfHeight}, {1., 0., 0.}, {0., 1., 0.}, radius));
  BOOST_REQUIRE(addDiskSurface(trimmed, {0., 0., -halfHeight}, {1., 0., 0.}, {0., -1., 0.}, radius));
  trimmed.CloseShape(false);
  BOOST_CHECK_EQUAL(trimmed.GetNsurfaces(), 3);

  const auto restoredTrimmed = writeAndReadBack(trimmed);
  BOOST_REQUIRE(restoredTrimmed != nullptr);
  checkSolidsIdentical(trimmed, *restoredTrimmed, 4.5, 5);

  // an unnavigable solid must come back unnavigable: the failure mode S1 describes is precisely a
  // defective solid that acquires a clean bill of health by losing its surfaces on the way
  SurfaceSolid openBox("persistOpenBox");
  for (int faceIndex = 0; faceIndex < 5; ++faceIndex) { // deliberately missing the sixth face
    BOOST_REQUIRE(addBoxFace(openBox, faceIndex, 1., 2., 3.));
  }
  openBox.CloseShape(false);
  BOOST_REQUIRE(!openBox.IsNavigable());
  BOOST_CHECK_EQUAL(static_cast<int>(openBox.GetNavigationReliability()),
                    static_cast<int>(SurfaceSolid::NavigationReliability::OpenSurfaceSet));

  const auto restoredOpenBox = writeAndReadBack(openBox);
  BOOST_REQUIRE(restoredOpenBox != nullptr);
  BOOST_CHECK(!restoredOpenBox->IsNavigable());
  checkSolidsIdentical(openBox, *restoredOpenBox, 4.5, 5);
}

// A solid that reaches the reader with no surface records -- a file written by an older version,
// or a solid streamed before CloseShape() -- must report Undetermined rather than manufacture a
// clean ClosureReport out of an empty surface set. "I do not know" is the only honest answer, and
// the difference matters: NavigationReliability is the flag callers are told to check.
BOOST_AUTO_TEST_CASE(EmptySolidIsNotReliable)
{
  SurfaceSolid empty("emptySolid");
  BOOST_CHECK_EQUAL(static_cast<int>(empty.GetNavigationReliability()),
                    static_cast<int>(SurfaceSolid::NavigationReliability::Undetermined));
  BOOST_CHECK(!empty.IsNavigable());

  // CloseShape on an empty surface set must not define the shape, with or without checking
  empty.CloseShape(false);
  BOOST_CHECK(!empty.IsDefined());
  BOOST_CHECK_EQUAL(static_cast<int>(empty.GetNavigationReliability()),
                    static_cast<int>(SurfaceSolid::NavigationReliability::Undetermined));
  BOOST_CHECK(!empty.IsNavigable());
  BOOST_CHECK(!empty.IsClosed());
}


namespace
{
// A golden-angle spiral of unit directions: quasi-uniform on the sphere, so no two are
// near-parallel and none aligns with a coordinate axis or a 45-degree symmetry plane. Used to
// test the invariant that containment does not depend on where the parity ray is aimed.
std::vector<std::array<double, 3>> spiralDirections(int count)
{
  std::vector<std::array<double, 3>> directions;
  directions.reserve(count);
  for (int index = 0; index < count; ++index) {
    const double cosTheta = 1. - 2. * (index + 0.5) / count;
    const double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
    const double phi = 2.399963229728653 * index;
    directions.push_back({sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta});
  }
  return directions;
}
} // namespace

// Parity containment answers a topological question, so on a closed, consistently oriented
// 2-manifold it cannot depend on where the ray is aimed. That invariant is what licenses the
// single-shot fast path: Contains() casts one fixed direction and stops.
//
// It is also the sharpest available oracle for the surface set itself -- no reference shape is
// involved, only the solid disagreeing with itself. Measured over the Phase 0 corpus, every part
// the closure check calls Reliable has *zero* direction disagreements in 11k points, and every
// part with disagreements is one the closure check already rejects; see
// scripts/geometry/CodeReview_Fable.md Section 4.4.
BOOST_AUTO_TEST_CASE(ContainsIsDirectionIndependentOnClosedSolids)
{
  const auto box = makeBoxSolid("dirBox", 1., 2., 3.);
  const auto tube = makeTubeSolid("dirTube", 1., 2., 3.);
  const auto cone = makeConeSolid("dirCone", 2., 1., 3.);
  const auto sphere = makeSphereSolid("dirSphere", 2.);
  const auto torus = makeTorusSolid("dirTorus", 3., 1.);
  const auto capsule = makeCapsuleSolid("dirCapsule", 2., 3.);

  const auto directions = spiralDirections(13);
  for (const auto* solid : {box.get(), tube.get(), cone.get(), sphere.get(), torus.get(), capsule.get()}) {
    BOOST_TEST_CONTEXT("solid = " << solid->GetName())
    {
      BOOST_REQUIRE(solid->IsNavigable());
      for (const auto& point : probeGrid(4.5, 7)) {
        const bool reference = solid->Contains(point.data());
        for (const auto& direction : directions) {
          BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ") direction = ("
                                         << direction[0] << ", " << direction[1] << ", " << direction[2] << ")")
          {
            BOOST_CHECK_EQUAL(solid->ContainsAlongDirection(point.data(), direction.data()), reference);
          }
        }
      }
    }
  }
}

// Section 4.4's re-shoot, on the defect it exists for. A gap in the surface set costs the parity
// ray exactly the crossings that fall inside the gap, so a point is misclassified over the whole
// *shadow* of the gap along the shooting direction -- centimetres of wrong answers arbitrarily far
// from any surface. Aiming the ray somewhere else escapes that shadow, which is why a majority
// over several directions recovers the right answer: measured over the 55 points where the single
// fixed direction disagrees with the OpenCascade oracle on the Phase 0 corpus, not one point is
// wrong in every direction.
//
// The fixture makes the mechanism explicit rather than statistical: the +x face of a box is split
// into two rectangles with a thin strip left out, so a ray leaving along +x from inside sees no
// crossing at all and reports "outside".
BOOST_AUTO_TEST_CASE(ContainsReshootsThroughSurfaceGaps)
{
  constexpr double halfX = 1.;
  constexpr double halfY = 2.;
  constexpr double halfZ = 3.;
  constexpr double gap = 0.05; // half-width in z of the missing strip on the +x face

  SurfaceSolid gapped("gappedBox");
  for (int faceIndex = 1; faceIndex < 6; ++faceIndex) { // every face but +x
    BOOST_REQUIRE(addBoxFace(gapped, faceIndex, halfX, halfY, halfZ));
  }
  // the +x face as two rectangles, leaving z in (-gap, +gap) uncovered. Frame of face 0:
  // origin (halfX, -halfY, -halfZ), axisU = +y, axisV = +z.
  BOOST_REQUIRE(gapped.AddPlanarSurface({halfX, -halfY, -halfZ}, {0., 1., 0.}, {0., 0., 1.},
                                        rectangleWire(2. * halfY, halfZ - gap)));
  BOOST_REQUIRE(gapped.AddPlanarSurface({halfX, -halfY, gap}, {0., 1., 0.}, {0., 0., 1.},
                                        rectangleWire(2. * halfY, halfZ - gap)));
  gapped.CloseShape(false);

  // the gap is what makes the solid unnavigable, and only an unnavigable solid re-shoots
  BOOST_REQUIRE(!gapped.IsNavigable());
  BOOST_CHECK_EQUAL(static_cast<int>(gapped.GetNavigationReliability()),
                    static_cast<int>(SurfaceSolid::NavigationReliability::OpenSurfaceSet));

  // a point deep inside whose +x ray leaves straight through the gap
  const std::array<double, 3> insidePoint{0., 0.3, 0.};
  const std::array<double, 3> throughGap{1., 0., 0.};
  BOOST_CHECK(!gapped.ContainsAlongDirection(insidePoint.data(), throughGap.data())); // the defect itself
  BOOST_CHECK(gapped.Contains(insidePoint.data()));                                   // the re-shoot recovers it
  BOOST_CHECK(gapped.Contains_Loop(insidePoint.data()));                              // ... on both paths

  // the same point in the same box *without* the gap is inside from every direction, so the
  // fixture isolates the gap and not some accident of the point
  const auto intact = makeBoxSolid("intactBox", halfX, halfY, halfZ);
  BOOST_CHECK(intact->Contains(insidePoint.data()));
  BOOST_CHECK(intact->ContainsAlongDirection(insidePoint.data(), throughGap.data()));

  // and the BVH and loop parities still agree everywhere on the defective solid: the re-shoot is
  // applied by one shared helper, so it can never make the two paths differ
  for (const auto& point : probeGrid(4.5, 7)) {
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      BOOST_CHECK_EQUAL(gapped.Contains(point.data()), gapped.Contains_Loop(point.data()));
    }
  }
}


namespace
{
// An L-shaped prism: the adversarial fixture the review found missing (CodeReview_Fable.md S4,
// "no concave fixture exists"). Its footprint is ([0,3]x[0,1]) union ([0,1]x[1,2]) extruded over
// z in [0, height], so it has a *reflex* (concave) vertical edge at x = 1, y = 1 -- the one place
// where a ray can touch the boundary from inside and stay inside, which no convex fixture can
// reproduce. Built from eight planar faces with outward normals, so it is closed and consistently
// oriented.
std::unique_ptr<SurfaceSolid> makeLPrismSolid(const char* name, double height = 1.)
{
  // footprint, counter-clockwise; (1,1) is the reflex vertex
  const std::vector<Point2D> footprint{{0., 0.}, {3., 0.}, {3., 1.}, {1., 1.}, {1., 2.}, {0., 2.}};

  auto solid = std::make_unique<SurfaceSolid>(name);

  // bottom (outward normal -z) and top (+z); axisU x axisV fixes the normal
  std::vector<Point2D> bottomWire;
  bottomWire.reserve(footprint.size());
  for (const auto& vertex : footprint) {
    bottomWire.push_back({vertex[1], vertex[0]}); // (u, v) = (y, x) so that axisU x axisV = -z
  }
  BOOST_REQUIRE(solid->AddPlanarSurface({0., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}, bottomWire));
  BOOST_REQUIRE(solid->AddPlanarSurface({0., 0., height}, {1., 0., 0.}, {0., 1., 0.}, footprint));

  // one vertical wall per footprint edge; axisU along the edge and axisV = +z put the normal at
  // (dy, -dx, 0), which points out of a counter-clockwise footprint
  for (size_t index = 0; index < footprint.size(); ++index) {
    const auto& start = footprint[index];
    const auto& end = footprint[(index + 1) % footprint.size()];
    const double deltaU = end[0] - start[0];
    const double deltaV = end[1] - start[1];
    const double length = std::hypot(deltaU, deltaV);
    BOOST_REQUIRE(solid->AddPlanarSurface({start[0], start[1], 0.}, {deltaU / length, deltaV / length, 0.},
                                          {0., 0., 1.}, rectangleWire(length, height)));
  }
  solid->CloseShape();
  return solid;
}
} // namespace

// S2: the bounding-box pre-check ran *before* the documented "no BVH yet, fall back to the plain
// loop" branch. Before CloseShape() the box is still all zeros, so the pre-check rejected every
// point outside a 1e-9 cube at the origin and the fallback was unreachable -- Contains() was
// effectively disabled on any solid that had not been closed yet.
BOOST_AUTO_TEST_CASE(ContainsWorksBeforeCloseShape)
{
  SurfaceSolid box("preCloseBox");
  addBoxSurfaces(box, 1., 2., 3.);
  BOOST_REQUIRE(!box.HasBVH()); // the premise: no acceleration structure and no bounding box yet

  for (const auto& point : probeGrid(4.5, 5)) {
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      const bool inside = std::abs(point[0]) < 1. && std::abs(point[1]) < 2. && std::abs(point[2]) < 3.;
      BOOST_CHECK_EQUAL(box.Contains(point.data()), inside);
      BOOST_CHECK_EQUAL(box.Contains(point.data()), box.Contains_Loop(point.data()));
    }
  }
}

// S3: a point *on* a face is inside for Contains, but its t = 0 exit was below the minimum ray
// parameter and therefore invisible to DistFromInside, which then returned Big. A navigator that
// asks "how far to the wall" while standing on the wall and is told "never" tunnels straight
// through the geometry. ROOT's own primitives answer 0 here, and so must this.
BOOST_AUTO_TEST_CASE(BoundaryPointsAgreeBetweenContainsAndDistances)
{
  const auto box = makeBoxSolid("boundaryPolicyBox", 1., 2., 3.);
  const std::array<double, 3> onFace{1., 0.5, 0.5}; // exactly on the +x face
  const std::array<double, 3> outward{1., 0., 0.};
  const std::array<double, 3> inward{-1., 0., 0.};

  BOOST_CHECK(box->Contains(onFace.data())); // documented policy: on a face counts as inside

  TGeoBBox reference("boundaryPolicyReference", 1., 2., 3.);
  BOOST_CHECK_EQUAL(box->DistFromInside(onFace.data(), outward.data(), 3), 0.);
  BOOST_CHECK_EQUAL(box->DistFromOutside(onFace.data(), inward.data(), 3), 0.);
  checkClose(box->DistFromInside(onFace.data(), outward.data(), 3),
             reference.DistFromInside(onFace.data(), outward.data(), 3));
  checkClose(box->DistFromOutside(onFace.data(), inward.data(), 3),
             reference.DistFromOutside(onFace.data(), inward.data(), 3));

  // going the other way the far wall is still the answer, so the fix is not "always return 0"
  checkClose(box->DistFromInside(onFace.data(), inward.data(), 3), 2.);

  // the BVH and loop paths must agree on all of it
  for (const auto& direction : {outward, inward}) {
    checkDistanceAgainstLoop(*box, onFace, direction);
  }
}

// S4 / S5: a ray that only *touches* the boundary has not crossed it. Contains() knows this --
// near-equal hits are clustered and a cluster carrying both an entering and an exiting hit
// contributes even parity -- but the distance queries classified every hit on its own, so they
// reported the touch as a crossing. The two then disagree about the same ray: DistFromOutside
// hands the navigator a step to the touch point, Contains says it is still outside once it gets
// there, and the navigator takes zero-length steps forever.
//
// Both flavours are covered. A convex edge graze (box) is the outside-facing case, and the L-prism
// reflex edge is the inside-facing one, which no convex solid can produce: there the ray leaves
// and re-enters the material at a single point and must be reported as never having left.
BOOST_AUTO_TEST_CASE(EdgeGrazesAreNotCrossings)
{
  const double invSqrt2 = 1. / std::sqrt(2.);

  // --- convex: touch the box edge x = +1, y = +2 and stay outside on both sides of the touch
  const auto box = makeBoxSolid("grazeBox", 1., 2., 3.);
  const std::array<double, 3> grazeDirection{invSqrt2, -invSqrt2, 0.};
  const std::array<double, 3> grazeOrigin{1. - 5. * invSqrt2, 2. + 5. * invSqrt2, 0.};
  BOOST_REQUIRE(!box->Contains(grazeOrigin.data()));

  // a point just past the touch is still outside, so nothing was entered ...
  const std::array<double, 3> pastTouch{1. + 1.e-3 * invSqrt2, 2. - 1.e-3 * invSqrt2, 0.};
  BOOST_REQUIRE(!box->Contains(pastTouch.data()));
  // ... and the distance query must say so too
  BOOST_CHECK_EQUAL(box->DistFromOutside(grazeOrigin.data(), grazeDirection.data(), 3), TGeoShape::Big());
  BOOST_CHECK_EQUAL(box->DistFromOutside_Loop(grazeOrigin.data(), grazeDirection.data()), TGeoShape::Big());

  // --- concave: the L-prism's reflex edge at x = 1, y = 1
  const auto prism = makeLPrismSolid("grazePrism");
  BOOST_REQUIRE(prism->IsNavigable());
  checkClose(prism->Capacity(), 4.); // 3x1 plus 1x1, extruded over unit height

  // a ray through the reflex edge along (1,-1): inside before the touch, inside after it, so the
  // touch is not an exit. The real exit is where it leaves the long arm at y = 0.
  const std::array<double, 3> reflexDirection{invSqrt2, -invSqrt2, 0.};
  const std::array<double, 3> reflexOrigin{1. - 0.5 * invSqrt2, 1. + 0.5 * invSqrt2, 0.5};
  BOOST_REQUIRE(prism->Contains(reflexOrigin.data()));
  const std::array<double, 3> pastReflex{1. + 1.e-3 * invSqrt2, 1. - 1.e-3 * invSqrt2, 0.5};
  BOOST_REQUIRE(prism->Contains(pastReflex.data())); // still inside: the touch was not an exit

  const double touchDistance = 0.5;
  const double exitDistance = 0.5 + std::sqrt(2.); // on to (2, 0, 0.5), in the middle of the y = 0 wall
  const double reported = prism->DistFromInside(reflexOrigin.data(), reflexDirection.data(), 3);
  BOOST_CHECK_GT(reported, touchDistance + 1.e-6); // the touch is not the answer ...
  checkClose(reported, exitDistance);              // ... the far wall is
  BOOST_CHECK_EQUAL(prism->DistFromInside_Loop(reflexOrigin.data(), reflexDirection.data()), reported);
}

// The concave fixture earns its keep beyond the single grazing ray: the whole sweep battery is
// run on it, since every invariant the convex fixtures pin (BVH == loop, direction-independent
// parity, Contains consistent with the distance answers) is weaker on shapes with no reflex edge.
BOOST_AUTO_TEST_CASE(LPrismSweeps)
{
  const auto prism = makeLPrismSolid("sweepPrism");
  BOOST_REQUIRE(prism->IsNavigable());

  sweepDistanceAgainstLoop(*prism, 3.5, 5);

  const auto directions = spiralDirections(13);
  for (const auto& point : probeGrid(3.5, 7)) {
    const bool inside = prism->Contains(point.data());
    BOOST_TEST_CONTEXT("point = (" << point[0] << ", " << point[1] << ", " << point[2] << ")")
    {
      BOOST_CHECK_EQUAL(prism->Contains_Loop(point.data()), inside);
      for (const auto& direction : directions) {
        BOOST_CHECK_EQUAL(prism->ContainsAlongDirection(point.data(), direction.data()), inside);
      }
      // the closed-form answer for the extruded L footprint
      const bool expected = point[2] > 0. && point[2] < 1. &&
                            ((point[0] > 0. && point[0] < 3. && point[1] > 0. && point[1] < 1.) ||
                             (point[0] > 0. && point[0] < 1. && point[1] >= 1. && point[1] < 2.));
      BOOST_CHECK_EQUAL(inside, expected);
    }
  }
}
