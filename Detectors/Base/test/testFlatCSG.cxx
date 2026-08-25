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
#include "TGeoShape.h"
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

BOOST_AUTO_TEST_CASE(box_distances_match_TGeoBBox)
{
  O2FlatCSG solid("box_dist");
  addBoxCell(solid, 3., 4., 5.);
  TGeoBBox reference(3., 4., 5.);

  Rng rng(4242ULL);
  for (int trial = 0; trial < 20000; ++trial) {
    double point[3] = {rng.uniform(-12., 12.), rng.uniform(-12., 12.), rng.uniform(-12., 12.)};
    double dir[3];
    double norm = 0.;
    do {
      for (int index = 0; index < 3; ++index) {
        dir[index] = rng.uniform(-1., 1.);
      }
      norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    } while (norm < 1.e-3);
    for (int index = 0; index < 3; ++index) {
      dir[index] /= norm;
    }
    const bool inside = reference.Contains(point);
    if (inside != static_cast<bool>(solid.Contains_Loop(point))) {
      continue; // a boundary point; task 1 already covers classification
    }
    const double mine = inside ? solid.DistFromInside_Loop(point, dir, TGeoShape::Big())
                               : solid.DistFromOutside_Loop(point, dir, TGeoShape::Big());
    const double theirs = inside ? reference.DistFromInside(point, dir, 3, TGeoShape::Big(), nullptr)
                                 : reference.DistFromOutside(point, dir, 3, TGeoShape::Big(), nullptr);
    if (theirs >= TGeoShape::Big()) {
      BOOST_REQUIRE_GE(mine, TGeoShape::Big());
    } else {
      BOOST_REQUIRE_SMALL(mine - theirs, 1.e-9);
    }
  }
}

BOOST_AUTO_TEST_CASE(tube_distances_match_TGeoTube_through_the_bore)
{
  // the complemented inner cylinder makes the occupancy along a ray TWO intervals for a ray that
  // crosses the bore, which is the case a convexity assumption would get wrong
  O2FlatCSG solid("tube_dist");
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
  solid.AddCell(0, 4, 0.);

  TGeoTube reference(2., 5., 7.);
  // a ray straight along +x at z = 0 enters the wall at x = -5, leaves it at x = -2, re-enters at
  // x = +2 and leaves at x = +5
  const double origin[3] = {-9., 0., 0.};
  const double dir[3] = {1., 0., 0.};
  BOOST_CHECK_SMALL(solid.DistFromOutside_Loop(origin, dir, TGeoShape::Big()) - 4., 1.e-12);

  const double inWall[3] = {-4., 0., 0.};
  BOOST_CHECK_SMALL(solid.DistFromInside_Loop(inWall, dir, TGeoShape::Big()) - 2., 1.e-12);

  const double inBore[3] = {0., 0., 0.};
  BOOST_CHECK(!solid.Contains_Loop(inBore));
  BOOST_CHECK_SMALL(solid.DistFromOutside_Loop(inBore, dir, TGeoShape::Big()) - 2., 1.e-12);

  Rng rng(99ULL);
  for (int trial = 0; trial < 20000; ++trial) {
    double point[3] = {rng.uniform(-9., 9.), rng.uniform(-9., 9.), rng.uniform(-10., 10.)};
    double direction[3];
    double norm = 0.;
    do {
      for (int index = 0; index < 3; ++index) {
        direction[index] = rng.uniform(-1., 1.);
      }
      norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
    } while (norm < 1.e-3);
    for (int index = 0; index < 3; ++index) {
      direction[index] /= norm;
    }
    const bool inside = reference.Contains(point);
    if (inside != static_cast<bool>(solid.Contains_Loop(point))) {
      continue;
    }
    const double mine = inside ? solid.DistFromInside_Loop(point, direction, TGeoShape::Big())
                               : solid.DistFromOutside_Loop(point, direction, TGeoShape::Big());
    const double theirs = inside ? reference.DistFromInside(point, direction, 3, TGeoShape::Big(), nullptr)
                                 : reference.DistFromOutside(point, direction, 3, TGeoShape::Big(), nullptr);
    if (theirs >= TGeoShape::Big()) {
      BOOST_REQUIRE_GE(mine, TGeoShape::Big());
    } else {
      BOOST_REQUIRE_SMALL(mine - theirs, 1.e-8);
    }
  }
}

BOOST_AUTO_TEST_CASE(a_ray_leaving_one_cell_into_a_touching_one_does_not_stop_between_them)
{
  // two unit boxes sharing the face at x = 1: the union's DistFromInside from the origin along +x
  // is 3, not 1. This is why DistFromInside needs the union across cells and not one cell's exit.
  O2FlatCSG solid("touching");
  addBoxCell(solid, 1., 1., 1.);
  const int first = solid.GetNhalfspaces();
  const double planes[6][2][3] = {{{1., 0., 0.}, {3., 0., 0.}},
                                  {{-1., 0., 0.}, {1., 0., 0.}},
                                  {{0., 1., 0.}, {0., 1., 0.}},
                                  {{0., -1., 0.}, {0., -1., 0.}},
                                  {{0., 0., 1.}, {0., 0., 1.}},
                                  {{0., 0., -1.}, {0., 0., -1.}}};
  for (const auto& plane : planes) {
    double coeff[10];
    planeQuadric(plane[0], plane[1], coeff);
    solid.AddQuadric(1., coeff);
  }
  solid.AddCell(first, 6, 8.);

  const double origin[3] = {0., 0., 0.};
  const double dir[3] = {1., 0., 0.};
  BOOST_CHECK_SMALL(solid.DistFromInside_Loop(origin, dir, TGeoShape::Big()) - 3., 1.e-12);
}

BOOST_AUTO_TEST_CASE(tangential_ray_on_a_cylinder_from_a_point_on_its_surface_has_no_nan_root)
{
  // a ray tangential to a cylinder, starting exactly on its surface, has beta == 0 and gamma == 0
  // together in HalfspaceRoots' quadratic -- the q == 0 case that used to divide 0./0. into a
  // NaN second root instead of recognising the single double root at t = 0
  O2FlatCSG solid("tangent_ray");
  double coeff[10];
  zCylinderQuadric(5., coeff);
  solid.AddQuadric(1., coeff);
  const auto& cylinder = solid.GetHalfspace(0);

  const double origin[3] = {5., 0., 0.};
  const double dir[3] = {0., 1., 0.};
  double roots[4];
  const int found = O2FlatCSG::HalfspaceRoots(cylinder, origin, dir, roots);

  BOOST_REQUIRE_EQUAL(found, 1);
  BOOST_CHECK(std::isfinite(roots[0]));
  BOOST_CHECK_SMALL(roots[0], 1.e-12);

  // the twin: an independent check that the reported root really is one, by plugging it back
  // into the surface equation directly rather than trusting the root-finder's own algebra
  const double hit[3] = {origin[0] + roots[0] * dir[0], origin[1] + roots[0] * dir[1],
                         origin[2] + roots[0] * dir[2]};
  BOOST_CHECK_SMALL(O2FlatCSG::EvalHalfspace(cylinder, hit), 1.e-9);
}
