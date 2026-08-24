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

#define BOOST_TEST_MODULE Test O2FlatCSG class
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "DetectorsBase/O2FlatCSG.h"

#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TMath.h"

#include <cmath>
#include <vector>

namespace
{
using o2::base::O2FlatCSG;

/// A small deterministic generator, so a failing case is reproducible from its seed alone.
class Rng
{
 public:
  explicit Rng(unsigned long long seed) : mState(seed) {}
  double uniform(double low, double high)
  {
    mState = mState * 6364136223846793005ULL + 1442695040888963407ULL;
    const double unit = static_cast<double>((mState >> 11) & ((1ULL << 53) - 1)) / static_cast<double>(1ULL << 53);
    return low + unit * (high - low);
  }

 private:
  unsigned long long mState;
};

/// The quadric of the plane with outward unit normal \a n through \a p: Q(x) = n.(x - p).
void planeQuadric(const double n[3], const double p[3], double coeff[10])
{
  for (int index = 0; index < 6; ++index) {
    coeff[index] = 0.;
  }
  coeff[6] = 0.5 * n[0];
  coeff[7] = 0.5 * n[1];
  coeff[8] = 0.5 * n[2];
  coeff[9] = -(n[0] * p[0] + n[1] * p[1] + n[2] * p[2]);
}

/// The quadric of the cylinder of radius \a r about the z axis: Q(x) = x^2 + y^2 - r^2.
void zCylinderQuadric(double r, double coeff[10])
{
  const double values[10] = {1., 0., 0., 1., 0., 0., 0., 0., 0., -r * r};
  for (int index = 0; index < 10; ++index) {
    coeff[index] = values[index];
  }
}

/// A box of half-extents (dx, dy, dz) centred on the origin, as one cell of six planes.
void addBoxCell(O2FlatCSG& solid, double dx, double dy, double dz)
{
  const double half[3] = {dx, dy, dz};
  const int first = solid.GetNhalfspaces();
  for (int axis = 0; axis < 3; ++axis) {
    for (int sense = -1; sense <= 1; sense += 2) {
      double normal[3] = {0., 0., 0.};
      double through[3] = {0., 0., 0.};
      normal[axis] = static_cast<double>(sense);
      through[axis] = sense * half[axis];
      double coeff[10];
      planeQuadric(normal, through, coeff);
      solid.AddQuadric(1., coeff);
    }
  }
  solid.AddCell(first, 6, 8. * dx * dy * dz);
}
} // namespace

BOOST_AUTO_TEST_CASE(box_from_six_planes_contains_like_TGeoBBox)
{
  O2FlatCSG solid("box");
  addBoxCell(solid, 3., 4., 5.);
  BOOST_CHECK_EQUAL(solid.GetNcells(), 1);
  BOOST_CHECK_EQUAL(solid.GetNhalfspaces(), 6);

  TGeoBBox reference(3., 4., 5.);
  Rng rng(20260824ULL);
  int scored = 0;
  for (int trial = 0; trial < 20000; ++trial) {
    const double point[3] = {rng.uniform(-6., 6.), rng.uniform(-7., 7.), rng.uniform(-8., 8.)};
    // skip the boundary shell, where the two shapes are allowed to disagree by tolerance
    if (std::abs(std::abs(point[0]) - 3.) < 1.e-9 || std::abs(std::abs(point[1]) - 4.) < 1.e-9 ||
        std::abs(std::abs(point[2]) - 5.) < 1.e-9) {
      continue;
    }
    BOOST_REQUIRE_EQUAL(solid.Contains_Loop(point), reference.Contains(point));
    BOOST_REQUIRE_EQUAL(solid.Contains(point), solid.Contains_Loop(point));
    ++scored;
  }
  BOOST_CHECK_GT(scored, 19000);
}

BOOST_AUTO_TEST_CASE(tube_from_two_cylinders_and_two_planes_contains_like_TGeoTube)
{
  // rmin = 2, rmax = 5, dz = 7: the inner cylinder is a COMPLEMENTED halfspace, which is what
  // makes this cell non-convex and is the case the whole class exists for.
  O2FlatCSG solid("tube");
  double coeff[10];
  zCylinderQuadric(5., coeff);
  solid.AddQuadric(1., coeff);
  zCylinderQuadric(2., coeff);
  solid.AddQuadric(-1., coeff);
  const double up[3] = {0., 0., 1.};
  const double down[3] = {0., 0., -1.};
  const double top[3] = {0., 0., 7.};
  const double bottom[3] = {0., 0., -7.};
  planeQuadric(up, top, coeff);
  solid.AddQuadric(1., coeff);
  planeQuadric(down, bottom, coeff);
  solid.AddQuadric(1., coeff);
  solid.AddCell(0, 4, TMath::Pi() * (25. - 4.) * 14.);

  TGeoTube reference(2., 5., 7.);
  Rng rng(777ULL);
  for (int trial = 0; trial < 20000; ++trial) {
    const double point[3] = {rng.uniform(-6., 6.), rng.uniform(-6., 6.), rng.uniform(-8., 8.)};
    const double radius = std::hypot(point[0], point[1]);
    if (std::abs(radius - 2.) < 1.e-9 || std::abs(radius - 5.) < 1.e-9 ||
        std::abs(std::abs(point[2]) - 7.) < 1.e-9) {
      continue;
    }
    BOOST_REQUIRE_EQUAL(solid.Contains_Loop(point), reference.Contains(point));
  }
}

BOOST_AUTO_TEST_CASE(two_disjoint_cells_are_a_union)
{
  O2FlatCSG solid("two_boxes");
  addBoxCell(solid, 1., 1., 1.);
  // a second box, centred at x = +10, as six planes of its own
  const int first = solid.GetNhalfspaces();
  const double centre = 10.;
  for (int axis = 0; axis < 3; ++axis) {
    for (int sense = -1; sense <= 1; sense += 2) {
      double normal[3] = {0., 0., 0.};
      double through[3] = {centre, 0., 0.};
      normal[axis] = static_cast<double>(sense);
      through[axis] += (axis == 0 ? sense * 1. : 0.);
      if (axis != 0) {
        through[axis] = sense * 1.;
      }
      double coeff[10];
      planeQuadric(normal, through, coeff);
      solid.AddQuadric(1., coeff);
    }
  }
  solid.AddCell(first, 6, 8.);

  const double inFirst[3] = {0., 0., 0.};
  const double inSecond[3] = {10., 0., 0.};
  const double between[3] = {5., 0., 0.};
  BOOST_CHECK(solid.Contains_Loop(inFirst));
  BOOST_CHECK(solid.Contains_Loop(inSecond));
  BOOST_CHECK(!solid.Contains_Loop(between));
}
