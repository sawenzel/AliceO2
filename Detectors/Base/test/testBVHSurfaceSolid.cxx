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

#include "TGeoBBox.h"
#include "TGeoShape.h"

#include <array>
#include <cmath>
#include <vector>

namespace
{
using SurfaceSolid = o2::base::O2BVHSurfaceSolid;
using Point2D = SurfaceSolid::Point2D;
using Point3D = SurfaceSolid::Point3D;

std::vector<Point2D> rectangleWire(double extentU, double extentV)
{
  return {{0., 0.}, {extentU, 0.}, {extentU, extentV}, {0., extentV}};
}

void addBoxSurfaces(SurfaceSolid& solid, double halfX, double halfY, double halfZ)
{
  BOOST_REQUIRE(solid.AddPlanarSurface({halfX, -halfY, -halfZ}, {0., 1., 0.}, {0., 0., 1.},
                                       rectangleWire(2. * halfY, 2. * halfZ)));
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfX, -halfY, -halfZ}, {0., 0., 1.}, {0., 1., 0.},
                                       rectangleWire(2. * halfZ, 2. * halfY)));
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfX, halfY, -halfZ}, {0., 0., 1.}, {1., 0., 0.},
                                       rectangleWire(2. * halfZ, 2. * halfX)));
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfX, -halfY, -halfZ}, {1., 0., 0.}, {0., 0., 1.},
                                       rectangleWire(2. * halfX, 2. * halfZ)));
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfX, -halfY, halfZ}, {1., 0., 0.}, {0., 1., 0.},
                                       rectangleWire(2. * halfX, 2. * halfY)));
  BOOST_REQUIRE(solid.AddPlanarSurface({-halfX, -halfY, -halfZ}, {0., 1., 0.}, {1., 0., 0.},
                                       rectangleWire(2. * halfY, 2. * halfX)));
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
}