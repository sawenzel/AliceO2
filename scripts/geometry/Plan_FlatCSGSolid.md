# O2FlatCSG Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `o2::base::O2FlatCSG`, a `TGeoShape` that stores a solid as a union of
intersection-cells over signed implicit halfspaces with a BVH over *sub-boxes of cells*, and ship
through it the 21 parts that rung R4 measured as over the boolean-tree budget.

**Architecture:** A part is a flat two-level DNF: one array of halfspaces (a ten-double general
quadric, or a torus in canonical form), one array of cells indexing into it. `CloseShape()`
recursively subdivides each cell's AABB, using a rigorous range bound of each halfspace over a box
to drop boxes that are wholly outside the cell and to drop halfspaces that hold everywhere inside
it; the surviving boxes go into a `bvh::v2::Bvh`. Queries traverse the BVH and, **inside each
box's slab only**, evaluate that box's short active list. Every accelerated query has a `_Loop`
twin that walks all cells and all halfspaces with no BVH and no pruning; the tests demand bit
identity.

**Tech Stack:** C++17, ROOT 6.36 (`TGeoBBox`, `TGeoShape`), the bvh2 v2 headers already vendored
at `Detectors/Base/src/bvh2_third_party.h`, `solveQuarticReal` from
`Detectors/Base/src/BoundedSurface.h`, Boost.Test via `o2_add_test`, and on the converter side
Python 3.10 with pythonOCC.

**Spec:** [`Design_FlatCSGSolid.md`](Design_FlatCSGSolid.md) — read it alongside this plan; every
task below argues from one of its sections.

## Global Constraints

- Branch `swenzel/bvhsurfacesolid`. Every task ends in a commit; commit messages follow
  `~/.claude/CLAUDE.md` (plain title, one lead sentence, short item list, no headings in the body).
- **Environment, eval and command in ONE shell invocation** (`Handoff_FlatCSG.md` §4):
  ```
  export ALIBUILD_WORK_DIR=$HOME/alisw/sw
  B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
  cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"
  export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
  ```
  The converter env adds OCC on top; the **sim env must stay separate** — the pythonOCC
  `PYTHONPATH` makes `o2-sim` segfault at startup.
- **One ninja at a time.** Reconfiguring CMake needs Clang on the prefix path:
  `export CMAKE_PREFIX_PATH=$HOME/alisw/sw/ubuntu2404_aarch64/Clang/v20.1.7-local1:$CMAKE_PREFIX_PATH`.
- **`rm` is blocked by a repo hook** (`.claude/hooks/deny-deletions.py`). Move files aside.
- **Conversions run strictly serially** — parallel `--csg auto` runs race and lose shapes.
- Floors that may not shrink, checked at every task that could touch them (spec §10): corpora
  PIPE 166/9/1 · ITS 250/14/0 · TPC 163/9/0 · ABSO 29/0/0 · TRD 1167/3/0; known-source 1775/1775;
  fixtures gate 10/10 csg exit 0; Bagger 13/13 + 7/7 exit 0 with CSG 7 bit-identical;
  `csg/emit.py --self-test` 301 with six digest tables; `checkKnownSource.py` 17;
  `O2_TGeoToCAD.py` 107; `O2_CADtoTGeo.py` 54; ctest `BVHSurfaceSolid` 113, `BVHAssembly` 22.
- Material side convention, used verbatim everywhere: a halfspace is inside when
  `sign * f(x) <= 0`.
- Never write into `STEP_examples/` or `ALICE_3_example/`.

---

### Task 1: The halfspace block and `Contains`, with no BVH

Builds the data model and the point query as a plain loop over every cell and every halfspace.
There is no acceleration yet on purpose: this loop becomes `Contains_Loop`, the twin that defines
the answer for every later task.

**Files:**
- Create: `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`
- Create: `Detectors/Base/src/O2FlatCSG.cxx`
- Create: `Detectors/Base/test/testFlatCSG.cxx`
- Modify: `Detectors/Base/CMakeLists.txt` (add the source at line 45 after `src/CADGeometryUtils.cxx`; add the header to the install list after line 90; add an `o2_add_test` block after the `BVHAssembly` one)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `struct o2::base::FlatCSGHalfspace { int kind; double sign; double c[11]; }` with
    `enum { kQuadric = 0, kTorus = 1 }`. Quadric uses `c[0..9]` as
    `(a00,a01,a02,a11,a12,a22,b0,b1,b2,cc)` for `Q(x)=xᵀAx+2bᵀx+cc`; torus uses `c[0..7]` as
    `(px,py,pz,dx,dy,dz,R,r)`.
  - `struct o2::base::FlatCSGCell { int first; int count; double volume; }`
  - `int O2FlatCSG::AddQuadric(double sign, const double coeff[10])` → halfspace index
  - `int O2FlatCSG::AddCell(int first, int count, double volume)` → cell index
  - `static double O2FlatCSG::EvalHalfspace(const FlatCSGHalfspace&, const double p[3])` →
    `sign * f(p)`; inside is `<= 0`
  - `Bool_t O2FlatCSG::Contains_Loop(const Double_t* p) const`
  - `Bool_t O2FlatCSG::Contains(const Double_t* p) const` (identical to the twin in this task)
  - `int O2FlatCSG::GetNcells() const`, `int O2FlatCSG::GetNhalfspaces() const`

- [ ] **Step 1: Write the failing test**

Create `Detectors/Base/test/testFlatCSG.cxx`:

```cpp
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
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)" \
  && export CMAKE_PREFIX_PATH=$HOME/alisw/sw/ubuntu2404_aarch64/Clang/v20.1.7-local1:$CMAKE_PREFIX_PATH \
  && cmake . >/dev/null && ninja O2test-detectorsbase-FlatCSG
```

Expected: FAIL — `DetectorsBase/O2FlatCSG.h: No such file or directory`.

- [ ] **Step 3: Write the header**

Create `Detectors/Base/include/DetectorsBase/O2FlatCSG.h` with the O2 copyright block (copy the
ten-line header verbatim from `O2BVHAssembly.h`), then:

```cpp
#ifndef ALICEO2_BASE_O2FLATCSG_
#define ALICEO2_BASE_O2FLATCSG_

#include "TGeoBBox.h"

#include <vector>

namespace o2
{
namespace base
{

/// One signed implicit halfspace: the region where `sign * f(x) <= 0`.
///
/// `kQuadric` stores `f(x) = x^T A x + 2 b^T x + c` in ten doubles,
/// `(a00, a01, a02, a11, a12, a22, b0, b1, b2, c)`. Plane, sphere, cylinder, cone and elliptic
/// cylinder are all this one block; see the carrier table in
/// scripts/geometry/Design_FlatCSGSolid.md section 3.1.
///
/// `kTorus` stores the canonical `(px, py, pz, dx, dy, dz, R, r)` in the first eight, because the
/// class wants the torus's exact signed distance for the range bound and its quartic for the ray,
/// and neither is nicer read off expanded quartic coefficients. No reference direction: a carrier
/// torus is a full surface of revolution.
struct FlatCSGHalfspace {
  enum Kind : int { kQuadric = 0,
                    kTorus = 1 };
  int kind = kQuadric;
  double sign = 1.;
  double c[11] = {};
};

/// One cell of the DNF: `[first, first + count)` of the shape's halfspace array, intersected.
///
/// `volume` is the cell's own volume, written by the converter from OCCT's `GProp` on the piece
/// the cell came from. The cells of a decomposition are disjoint, so `Capacity()` is their sum.
struct FlatCSGCell {
  int first = 0;
  int count = 0;
  double volume = 0.;
};

/// A solid stored as a union of intersection-cells over signed implicit halfspaces.
///
/// This is the flat two-level DNF of scripts/geometry/Stream_AA_FlatCSG.md section 5 -- depth is
/// two, always -- and it is what a part decomposed by `csg/decompose.py` actually is. Two things
/// motivate it over the `TGeoCompositeShape` the converter emits today. It is cheaper: a 66-leaf
/// boolean tree costs 66 virtual calls per `Contains` whatever the point is, whereas this class
/// costs the halfspaces that are still undecided in the box the query lands in. And it is more
/// faithful: ROOT has no usable halfspace, so the composite emission has to bound each one into a
/// padded native primitive sized to the part's neighbourhood, and this one does not.
///
/// Every accelerated query has a `_Loop` twin that walks all cells and all halfspaces with no BVH
/// and no active-list pruning. The twin is the definition of the answer and the tests require bit
/// identity; see scripts/geometry/Design_FlatCSGSolid.md section 6.
class O2FlatCSG : public TGeoBBox
{
 public:
  O2FlatCSG();
  explicit O2FlatCSG(const char* name);
  ~O2FlatCSG() override;

  // ---- building -------------------------------------------------------------------------
  /// Append a quadric halfspace; returns its index. `sign` is +1 or -1, inside is `sign*Q <= 0`.
  int AddQuadric(double sign, const double coeff[10]);
  /// Append a cell over `[first, first + count)` of the halfspace array; returns its index.
  int AddCell(int first, int count, double volume);

  int GetNhalfspaces() const { return static_cast<int>(fHalfspaces.size()); }
  int GetNcells() const { return static_cast<int>(fCells.size()); }
  const FlatCSGHalfspace& GetHalfspace(int index) const { return fHalfspaces[index]; }
  const FlatCSGCell& GetCell(int index) const { return fCells[index]; }

  /// `sign * f(point)`; the halfspace contains the point when this is `<= 0`.
  static double EvalHalfspace(const FlatCSGHalfspace& halfspace, const double* point);

  // ---- the TGeoShape contract ------------------------------------------------------------
  Bool_t Contains(const Double_t* point) const override;

  // ---- the reference twins ---------------------------------------------------------------
  Bool_t Contains_Loop(const Double_t* point) const;

 protected:
  /// True when every halfspace of cell `index` contains `point`.
  bool CellContains(int index, const double* point) const;

  std::vector<FlatCSGHalfspace> fHalfspaces; ///< the flat halfspace array
  std::vector<FlatCSGCell> fCells;           ///< the DNF's cells, indexing into it

  ClassDefOverride(O2FlatCSG, 1) // flat-DNF halfspace shape class
};

} // namespace base
} // namespace o2

#endif
```

- [ ] **Step 4: Write the minimal implementation**

Create `Detectors/Base/src/O2FlatCSG.cxx` with the same copyright block, then:

```cpp
#include "DetectorsBase/O2FlatCSG.h"

#include <cmath>

ClassImp(o2::base::O2FlatCSG);

namespace o2
{
namespace base
{

O2FlatCSG::O2FlatCSG() : TGeoBBox(0., 0., 0.) {}

O2FlatCSG::O2FlatCSG(const char* name) : TGeoBBox(name, 0., 0., 0.) {}

O2FlatCSG::~O2FlatCSG() = default;

int O2FlatCSG::AddQuadric(double sign, const double coeff[10])
{
  FlatCSGHalfspace halfspace;
  halfspace.kind = FlatCSGHalfspace::kQuadric;
  halfspace.sign = sign < 0. ? -1. : 1.;
  for (int index = 0; index < 10; ++index) {
    halfspace.c[index] = coeff[index];
  }
  fHalfspaces.push_back(halfspace);
  return static_cast<int>(fHalfspaces.size()) - 1;
}

int O2FlatCSG::AddCell(int first, int count, double volume)
{
  FlatCSGCell cell;
  cell.first = first;
  cell.count = count;
  cell.volume = volume;
  fCells.push_back(cell);
  return static_cast<int>(fCells.size()) - 1;
}

double O2FlatCSG::EvalHalfspace(const FlatCSGHalfspace& halfspace, const double* point)
{
  // the torus branch arrives in task 3; until then only quadrics are constructible
  const double* c = halfspace.c;
  const double x = point[0];
  const double y = point[1];
  const double z = point[2];
  const double quadratic = c[0] * x * x + c[3] * y * y + c[5] * z * z +
                           2. * (c[1] * x * y + c[2] * x * z + c[4] * y * z);
  const double linear = 2. * (c[6] * x + c[7] * y + c[8] * z);
  return halfspace.sign * (quadratic + linear + c[9]);
}

bool O2FlatCSG::CellContains(int index, const double* point) const
{
  const FlatCSGCell& cell = fCells[index];
  for (int offset = 0; offset < cell.count; ++offset) {
    if (EvalHalfspace(fHalfspaces[cell.first + offset], point) > 0.) {
      return false;
    }
  }
  return true;
}

Bool_t O2FlatCSG::Contains_Loop(const Double_t* point) const
{
  for (int index = 0; index < GetNcells(); ++index) {
    if (CellContains(index, point)) {
      return kTRUE;
    }
  }
  return kFALSE;
}

Bool_t O2FlatCSG::Contains(const Double_t* point) const
{
  // no acceleration until task 5; the twin is the answer
  return Contains_Loop(point);
}

} // namespace base
} // namespace o2
```

- [ ] **Step 5: Wire the build**

In `Detectors/Base/CMakeLists.txt`, add `src/O2FlatCSG.cxx` to the `SOURCES` list immediately
after `src/CADGeometryUtils.cxx`; add
`include/DetectorsBase/O2FlatCSG.h` to the install list immediately after
`include/DetectorsBase/O2BVHAssembly.h`; and add this block immediately after the `BVHAssembly`
`o2_add_test`:

```cmake
o2_add_test(
  FlatCSG
  SOURCES test/testFlatCSG.cxx
  COMPONENT_NAME DetectorsBase
  PUBLIC_LINK_LIBRARIES O2::DetectorsBase ROOT::Geom ROOT::RIO
  LABELS detectorsbase)
```

- [ ] **Step 6: Run the tests to verify they pass**

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)" \
  && export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH \
  && ninja O2test-detectorsbase-FlatCSG && ctest -R 'FlatCSG|BVHSurfaceSolid|BVHAssembly' --output-on-failure
```

Expected: PASS, 3 FlatCSG cases; `BVHSurfaceSolid` and `BVHAssembly` unchanged and green.

- [ ] **Step 7: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2FlatCSG.h \
  Detectors/Base/src/O2FlatCSG.cxx Detectors/Base/test/testFlatCSG.cxx Detectors/Base/CMakeLists.txt
git commit -m "$(cat <<'EOF'
Add the flat-DNF halfspace shape and its point query

This adds O2FlatCSG, which stores a solid as a union of intersection-cells
over signed implicit halfspaces, and its Contains.

- A halfspace is a ten-double general quadric with a sign; inside is sign*Q <= 0.
- A cell is a range of the flat halfspace array, and a part is its cells.
- Contains walks every cell and every halfspace; it is the twin later tasks are
  measured against, so it is deliberately not accelerated.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Ray roots and interval clipping — the distances, still with no BVH

The generic distance machinery, and the reason the class needs no convexity assumption: a cell's
occupancy along a ray is read off the sorted roots of its halfspaces by classifying midpoints.

**Files:**
- Modify: `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`
- Modify: `Detectors/Base/src/O2FlatCSG.cxx`
- Modify: `Detectors/Base/test/testFlatCSG.cxx`

**Interfaces:**
- Consumes: `FlatCSGHalfspace`, `FlatCSGCell`, `EvalHalfspace`, `CellContains` (Task 1).
- Produces:
  - `static int O2FlatCSG::HalfspaceRoots(const FlatCSGHalfspace&, const double o[3], const double d[3], double* t)` — writes up to 4 roots, unsorted, returns the count
  - `int O2FlatCSG::CellIntervals(int cell, const int* active, int nActive, const double o[3], const double d[3], double tlo, double thi, double* out, int maxOut) const` — writes `[enter, exit]` pairs into `out`, returns the pair count
  - `Double_t O2FlatCSG::DistFromOutside_Loop(const Double_t* p, const Double_t* d, Double_t step) const`
  - `Double_t O2FlatCSG::DistFromInside_Loop(const Double_t* p, const Double_t* d, Double_t step) const`
  - the `DistFromOutside` / `DistFromInside` overrides, delegating to the twins in this task

- [ ] **Step 1: Write the failing test**

Append to `Detectors/Base/test/testFlatCSG.cxx`:

```cpp
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
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)" \
  && ninja O2test-detectorsbase-FlatCSG
```

Expected: FAIL — `no member named 'DistFromInside_Loop' in 'o2::base::O2FlatCSG'`.

- [ ] **Step 3: Declare the new members in the header**

Add to the public section of `O2FlatCSG`, after `EvalHalfspace`:

```cpp
  /// Real roots of `sign * f(origin + t*dir) = 0`, unsorted, at most four; returns the count.
  static int HalfspaceRoots(const FlatCSGHalfspace& halfspace, const double* origin,
                            const double* dir, double* roots);

  /// The occupancy of cell \a cell along the ray, restricted to `[tlo, thi]`.
  ///
  /// \a active lists which of the cell's halfspaces are still undecided; pass `nullptr` with
  /// \a nActive `< 0` to use all of them, which is what the twins do. Writes `[enter, exit]`
  /// pairs into \a out and returns the pair count.
  ///
  /// No convexity is assumed anywhere: a complemented halfspace makes a cell non-convex and the
  /// occupancy several intervals, which is exactly what every part with a hole in it looks like.
  int CellIntervals(int cell, const int* active, int nActive, const double* origin,
                    const double* dir, double tlo, double thi, double* out, int maxOut) const;

  Double_t DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;
  Double_t DistFromInside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                          Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;

  Double_t DistFromOutside_Loop(const Double_t* point, const Double_t* dir,
                                Double_t step = TGeoShape::Big()) const;
  Double_t DistFromInside_Loop(const Double_t* point, const Double_t* dir,
                               Double_t step = TGeoShape::Big()) const;
```

and add `#include "TGeoShape.h"` next to the `TGeoBBox.h` include, plus a file-scope constant in
the `.cxx`:

```cpp
/// The most roots one cell can contribute to one ray: four per torus halfspace.
constexpr int kMaxRootsPerHalfspace = 4;
```

- [ ] **Step 4: Implement**

Add to `Detectors/Base/src/O2FlatCSG.cxx`, inside `namespace o2::base`:

```cpp
int O2FlatCSG::HalfspaceRoots(const FlatCSGHalfspace& halfspace, const double* origin,
                              const double* dir, double* roots)
{
  // the torus branch arrives in task 3
  const double* c = halfspace.c;
  // A d
  const double ad[3] = {c[0] * dir[0] + c[1] * dir[1] + c[2] * dir[2],
                        c[1] * dir[0] + c[3] * dir[1] + c[4] * dir[2],
                        c[2] * dir[0] + c[4] * dir[1] + c[5] * dir[2]};
  // A o + b
  const double aob[3] = {c[0] * origin[0] + c[1] * origin[1] + c[2] * origin[2] + c[6],
                         c[1] * origin[0] + c[3] * origin[1] + c[4] * origin[2] + c[7],
                         c[2] * origin[0] + c[4] * origin[1] + c[5] * origin[2] + c[8]};
  const double alpha = dir[0] * ad[0] + dir[1] * ad[1] + dir[2] * ad[2];
  const double beta = dir[0] * aob[0] + dir[1] * aob[1] + dir[2] * aob[2];
  const double gamma = EvalHalfspace(halfspace, origin) * halfspace.sign; // the unsigned Q(o)

  // scale-aware degeneracy: a plane has alpha exactly 0, and a ray parallel to a cylinder axis
  // has it near 0 -- both are the linear equation, not a badly conditioned quadratic
  const double reference = std::abs(beta) + std::abs(gamma) + 1.e-300;
  if (std::abs(alpha) <= 1.e-14 * reference) {
    if (std::abs(beta) <= 1.e-300) {
      return 0;
    }
    roots[0] = -0.5 * gamma / beta;
    return 1;
  }
  const double disc = beta * beta - alpha * gamma;
  if (disc < 0.) {
    return 0;
  }
  const double root = std::sqrt(disc);
  // the numerically stable pair, so a grazing ray does not lose the near root to cancellation
  const double q = -(beta + (beta >= 0. ? root : -root));
  roots[0] = q / alpha;
  roots[1] = gamma / q;
  return 2;
}

int O2FlatCSG::CellIntervals(int cell, const int* active, int nActive, const double* origin,
                             const double* dir, double tlo, double thi, double* out,
                             int maxOut) const
{
  const FlatCSGCell& description = fCells[cell];
  const int count = nActive < 0 ? description.count : nActive;
  if (thi <= tlo) {
    return 0;
  }

  // every root of every active halfspace, clipped to the window
  double breaks[2 + 64 * kMaxRootsPerHalfspace];
  int nBreaks = 0;
  breaks[nBreaks++] = tlo;
  breaks[nBreaks++] = thi;
  for (int slot = 0; slot < count; ++slot) {
    const int index = active != nullptr ? active[slot] : description.first + slot;
    double roots[kMaxRootsPerHalfspace];
    const int found = HalfspaceRoots(fHalfspaces[index], origin, dir, roots);
    for (int root = 0; root < found; ++root) {
      if (roots[root] > tlo && roots[root] < thi &&
          nBreaks < static_cast<int>(sizeof(breaks) / sizeof(breaks[0]))) {
        breaks[nBreaks++] = roots[root];
      }
    }
  }
  std::sort(breaks, breaks + nBreaks);

  // classify the midpoint of each sub-interval and merge the runs that are inside
  int pairs = 0;
  bool open = false;
  for (int index = 0; index + 1 < nBreaks; ++index) {
    const double lo = breaks[index];
    const double hi = breaks[index + 1];
    if (hi <= lo) {
      continue;
    }
    const double middle = 0.5 * (lo + hi);
    double probe[3] = {origin[0] + middle * dir[0], origin[1] + middle * dir[1],
                       origin[2] + middle * dir[2]};
    bool inside = true;
    for (int slot = 0; slot < count && inside; ++slot) {
      const int halfspace = active != nullptr ? active[slot] : description.first + slot;
      inside = EvalHalfspace(fHalfspaces[halfspace], probe) <= 0.;
    }
    if (inside) {
      if (open) {
        out[2 * (pairs - 1) + 1] = hi;
      } else if (pairs < maxOut) {
        out[2 * pairs] = lo;
        out[2 * pairs + 1] = hi;
        ++pairs;
        open = true;
      }
    } else {
      open = false;
    }
  }
  return pairs;
}
```

and the four distance entry points. Both twins gather the occupancy of **every** cell and merge:

```cpp
namespace
{
/// Merge `[enter, exit]` pairs in place, joining ones that touch within \a glue.
int mergeIntervals(double* pairs, int count, double glue)
{
  if (count < 2) {
    return count;
  }
  // sort by entry
  for (int outer = 1; outer < count; ++outer) {
    const double lo = pairs[2 * outer];
    const double hi = pairs[2 * outer + 1];
    int inner = outer - 1;
    while (inner >= 0 && pairs[2 * inner] > lo) {
      pairs[2 * (inner + 1)] = pairs[2 * inner];
      pairs[2 * (inner + 1) + 1] = pairs[2 * inner + 1];
      --inner;
    }
    pairs[2 * (inner + 1)] = lo;
    pairs[2 * (inner + 1) + 1] = hi;
  }
  int kept = 1;
  for (int index = 1; index < count; ++index) {
    if (pairs[2 * index] <= pairs[2 * (kept - 1) + 1] + glue) {
      pairs[2 * (kept - 1) + 1] = std::max(pairs[2 * (kept - 1) + 1], pairs[2 * index + 1]);
    } else {
      pairs[2 * kept] = pairs[2 * index];
      pairs[2 * kept + 1] = pairs[2 * index + 1];
      ++kept;
    }
  }
  return kept;
}
} // namespace

Double_t O2FlatCSG::DistFromOutside_Loop(const Double_t* point, const Double_t* dir,
                                         Double_t step) const
{
  double best = TGeoShape::Big();
  double pairs[2 * 64];
  for (int cell = 0; cell < GetNcells(); ++cell) {
    const int found = CellIntervals(cell, nullptr, -1, point, dir, 0., step, pairs, 64);
    for (int pair = 0; pair < found; ++pair) {
      // a point exactly on the boundary is already inside; only a real entry counts
      if (pairs[2 * pair + 1] > TGeoShape::Tolerance() && pairs[2 * pair] < best) {
        best = std::max(pairs[2 * pair], 0.);
      }
    }
  }
  return best;
}

Double_t O2FlatCSG::DistFromInside_Loop(const Double_t* point, const Double_t* dir,
                                        Double_t step) const
{
  // the union's occupancy, so a ray that leaves one cell into a touching one keeps going
  double pairs[2 * 256];
  int count = 0;
  for (int cell = 0; cell < GetNcells() && count < 128; ++cell) {
    count += CellIntervals(cell, nullptr, -1, point, dir, 0., step, pairs + 2 * count,
                           128 - count);
  }
  count = mergeIntervals(pairs, count, TGeoShape::Tolerance());
  for (int pair = 0; pair < count; ++pair) {
    if (pairs[2 * pair] <= TGeoShape::Tolerance()) {
      return pairs[2 * pair + 1];
    }
  }
  return 0.;
}

Double_t O2FlatCSG::DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact,
                                    Double_t step, Double_t* safe) const
{
  if (iact < 3 && safe != nullptr) {
    *safe = Safety(point, kFALSE);
    if (iact == 0) {
      return TGeoShape::Big();
    }
    if (iact == 1 && step < *safe) {
      return TGeoShape::Big();
    }
  }
  return DistFromOutside_Loop(point, dir, step); // accelerated in task 5
}

Double_t O2FlatCSG::DistFromInside(const Double_t* point, const Double_t* dir, Int_t iact,
                                   Double_t step, Double_t* safe) const
{
  if (iact < 3 && safe != nullptr) {
    *safe = Safety(point, kTRUE);
    if (iact == 0) {
      return TGeoShape::Big();
    }
    if (iact == 1 && step < *safe) {
      return TGeoShape::Big();
    }
  }
  return DistFromInside_Loop(point, dir, step); // accelerated in task 5
}
```

Add `#include <algorithm>` to the `.cxx`. `Safety` does not exist until Task 6, so for this task
add a placeholder override in the header and `.cxx` returning `0.` — always a legal safety — with
the comment `// task 6 replaces this; 0 is always sound`.

- [ ] **Step 5: Run the tests to verify they pass**

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)" \
  && export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH \
  && ninja O2test-detectorsbase-FlatCSG && ctest -R FlatCSG --output-on-failure
```

Expected: PASS, 6 cases.

- [ ] **Step 6: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2FlatCSG.h \
  Detectors/Base/src/O2FlatCSG.cxx Detectors/Base/test/testFlatCSG.cxx
git commit -m "$(cat <<'EOF'
Give the flat solid its distances by interval clipping

This adds the ray query to O2FlatCSG: a cell's occupancy along a ray is read
off the sorted roots of its halfspaces.

- Each quadric contributes at most two roots, solved in the stable pair form.
- A midpoint classification between consecutive roots gives the occupancy, so
  no convexity is assumed; a complemented halfspace makes cells non-convex.
- DistFromInside merges the occupancy across cells, so a ray leaving one cell
  into a touching one is not stopped between them.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: The torus halfspace

The one carrier that is not a quadric. Its ray intersection is a quartic and its range bound (Task
4) wants its exact signed distance, so it is stored in canonical form and gets its own branch in
the two static evaluators.

**Files:**
- Modify: `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`
- Modify: `Detectors/Base/src/O2FlatCSG.cxx`
- Modify: `Detectors/Base/test/testFlatCSG.cxx`

**Interfaces:**
- Consumes: `EvalHalfspace`, `HalfspaceRoots` (Tasks 1–2), `o2::base::solveQuarticReal` from `Detectors/Base/src/BoundedSurface.h`.
- Produces: `int O2FlatCSG::AddTorus(double sign, const double centre[3], const double axis[3], double major, double minor)` → halfspace index.

- [ ] **Step 1: Write the failing test**

Append to `Detectors/Base/test/testFlatCSG.cxx` (and add `#include "TGeoTorus.h"` at the top):

```cpp
BOOST_AUTO_TEST_CASE(torus_contains_and_distances_match_TGeoTorus)
{
  // a full torus, R = 10, r = 3, about z -- one cell of one halfspace
  O2FlatCSG solid("torus");
  const double centre[3] = {0., 0., 0.};
  const double axis[3] = {0., 0., 1.};
  solid.AddTorus(1., centre, axis, 10., 3.);
  solid.AddCell(0, 1, 2. * TMath::Pi() * TMath::Pi() * 10. * 9.);

  TGeoTorus reference(10., 0., 3.);
  Rng rng(31415ULL);
  int scoredPoints = 0;
  for (int trial = 0; trial < 20000; ++trial) {
    const double point[3] = {rng.uniform(-15., 15.), rng.uniform(-15., 15.), rng.uniform(-5., 5.)};
    const double radial = std::hypot(point[0], point[1]);
    const double distance = std::hypot(radial - 10., point[2]) - 3.;
    if (std::abs(distance) < 1.e-9) {
      continue;
    }
    BOOST_REQUIRE_EQUAL(solid.Contains_Loop(point), reference.Contains(point));
    ++scoredPoints;
  }
  BOOST_CHECK_GT(scoredPoints, 19000);

  for (int trial = 0; trial < 20000; ++trial) {
    double point[3] = {rng.uniform(-20., 20.), rng.uniform(-20., 20.), rng.uniform(-8., 8.)};
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
      continue;
    }
    const double mine = inside ? solid.DistFromInside_Loop(point, dir, TGeoShape::Big())
                               : solid.DistFromOutside_Loop(point, dir, TGeoShape::Big());
    const double theirs = inside ? reference.DistFromInside(point, dir, 3, TGeoShape::Big(), nullptr)
                                 : reference.DistFromOutside(point, dir, 3, TGeoShape::Big(), nullptr);
    if (theirs >= TGeoShape::Big()) {
      BOOST_REQUIRE_GE(mine, TGeoShape::Big());
    } else {
      // the quartic is the looser of the two solvers; 1e-6 cm on a 10 cm torus
      BOOST_REQUIRE_SMALL(mine - theirs, 1.e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(a_tilted_torus_is_the_same_solid_as_an_upright_one_rotated)
{
  // the frame handling is where a torus block goes wrong silently, so it gets its own case
  const double axis[3] = {0., 1. / std::sqrt(2.), 1. / std::sqrt(2.)};
  const double centre[3] = {1., 2., 3.};
  O2FlatCSG solid("tilted_torus");
  solid.AddTorus(1., centre, axis, 8., 2.);
  solid.AddCell(0, 1, 0.);

  Rng rng(2718ULL);
  for (int trial = 0; trial < 20000; ++trial) {
    const double point[3] = {rng.uniform(-14., 16.), rng.uniform(-13., 17.), rng.uniform(-12., 18.)};
    // the closed-form signed distance is the reference: sqrt((rho - R)^2 + z^2) - r
    const double offset[3] = {point[0] - centre[0], point[1] - centre[1], point[2] - centre[2]};
    const double along = offset[0] * axis[0] + offset[1] * axis[1] + offset[2] * axis[2];
    double radialVec[3];
    for (int index = 0; index < 3; ++index) {
      radialVec[index] = offset[index] - along * axis[index];
    }
    const double rho = std::sqrt(radialVec[0] * radialVec[0] + radialVec[1] * radialVec[1] +
                                 radialVec[2] * radialVec[2]);
    const double signedDistance = std::hypot(rho - 8., along) - 2.;
    if (std::abs(signedDistance) < 1.e-9) {
      continue;
    }
    BOOST_REQUIRE_EQUAL(static_cast<bool>(solid.Contains_Loop(point)), signedDistance < 0.);
  }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run the Task 2 build command. Expected: FAIL — `no member named 'AddTorus'`.

- [ ] **Step 3: Implement**

Declare in the header, next to `AddQuadric`:

```cpp
  /// Append a torus halfspace: inside is `sign * (sqrt((rho - major)^2 + z^2) - minor) <= 0`,
  /// with `rho` the distance from the axis through \a centre along \a axis and `z` the coordinate
  /// along it. \a axis must be a unit vector. Returns the halfspace index.
  int AddTorus(double sign, const double* centre, const double* axis, double major, double minor);
```

In the `.cxx`, add `#include "BoundedSurface.h"` (it is a header in `src/`, so the plain include
works from there) and:

```cpp
int O2FlatCSG::AddTorus(double sign, const double* centre, const double* axis, double major,
                        double minor)
{
  FlatCSGHalfspace halfspace;
  halfspace.kind = FlatCSGHalfspace::kTorus;
  halfspace.sign = sign < 0. ? -1. : 1.;
  for (int index = 0; index < 3; ++index) {
    halfspace.c[index] = centre[index];
    halfspace.c[3 + index] = axis[index];
  }
  halfspace.c[6] = major;
  halfspace.c[7] = minor;
  fHalfspaces.push_back(halfspace);
  return static_cast<int>(fHalfspaces.size()) - 1;
}
```

Add a torus branch at the top of `EvalHalfspace`, before the quadric code:

```cpp
  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    const double* c = halfspace.c;
    const double offset[3] = {point[0] - c[0], point[1] - c[1], point[2] - c[2]};
    const double along = offset[0] * c[3] + offset[1] * c[4] + offset[2] * c[5];
    const double radial[3] = {offset[0] - along * c[3], offset[1] - along * c[4],
                              offset[2] - along * c[5]};
    const double rho = std::sqrt(radial[0] * radial[0] + radial[1] * radial[1] +
                                 radial[2] * radial[2]);
    // the exact signed distance, which is 1-Lipschitz -- task 4's range bound needs that
    return halfspace.sign * (std::hypot(rho - c[6], along) - c[7]);
  }
```

Add a torus branch at the top of `HalfspaceRoots`. Working in the torus frame with `P` the origin
offset and `D` the direction, both resolved into (axis-perpendicular, axis) components:

```cpp
  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    const double* c = halfspace.c;
    const double axis[3] = {c[3], c[4], c[5]};
    const double major = c[6];
    const double minor = c[7];
    const double offset[3] = {origin[0] - c[0], origin[1] - c[1], origin[2] - c[2]};
    // components along the axis, and the perpendicular parts
    const double pz = offset[0] * axis[0] + offset[1] * axis[1] + offset[2] * axis[2];
    const double dz = dir[0] * axis[0] + dir[1] * axis[1] + dir[2] * axis[2];
    double pPerp[3];
    double dPerp[3];
    for (int index = 0; index < 3; ++index) {
      pPerp[index] = offset[index] - pz * axis[index];
      dPerp[index] = dir[index] - dz * axis[index];
    }
    const double pp = pPerp[0] * pPerp[0] + pPerp[1] * pPerp[1] + pPerp[2] * pPerp[2];
    const double dd = dPerp[0] * dPerp[0] + dPerp[1] * dPerp[1] + dPerp[2] * dPerp[2];
    const double pd = pPerp[0] * dPerp[0] + pPerp[1] * dPerp[1] + pPerp[2] * dPerp[2];
    // (|X|^2 + R^2 - r^2)^2 - 4 R^2 (X_perp . X_perp) = 0 with X = P + tD, |D| = 1
    const double e = pp + pz * pz + major * major - minor * minor;
    const double f = pd + pz * dz;
    const double a4 = 1.;
    const double a3 = 4. * f;
    const double a2 = 2. * e + 4. * f * f - 4. * major * major * dd;
    const double a1 = 4. * e * f - 8. * major * major * pd;
    const double a0 = e * e - 4. * major * major * pp;
    const std::vector<double> found = solveQuarticReal(a4, a3, a2, a1, a0);
    int count = 0;
    for (double value : found) {
      if (count < kMaxRootsPerHalfspace) {
        roots[count++] = value;
      }
    }
    return count;
  }
```

Note for the implementer: `solveQuarticReal` lives in `o2::base` in `src/BoundedSurface.h` — check
the enclosing namespace when you open it and qualify accordingly. It requires `a4 != 0`, which
holds here because the leading coefficient is exactly 1 for a unit direction; if the caller ever
passes a non-unit direction the whole class breaks, so assert `|d| = 1` in a debug build.

- [ ] **Step 4: Run the tests to verify they pass**

Run the Task 2 test command. Expected: PASS, 8 cases.

- [ ] **Step 5: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2FlatCSG.h \
  Detectors/Base/src/O2FlatCSG.cxx Detectors/Base/test/testFlatCSG.cxx
git commit -m "$(cat <<'EOF'
Add the torus halfspace to the flat solid

This adds the one carrier that is not a quadric, in canonical form rather than
as expanded quartic coefficients.

- Containment reads the torus's exact signed distance, which is 1-Lipschitz and
  is what the sub-cell range bound will need.
- The ray intersection is the standard quartic, solved by solveQuarticReal from
  BoundedSurface.h rather than by a second solver.
- A tilted torus gets its own case, because the frame is where a torus block
  goes wrong silently.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: The sub-cell box builder

Spec §4.2. Recursive subdivision of each cell's AABB with a rigorous range bound, dropping boxes
that are wholly outside the cell and halfspaces that hold everywhere inside a box. No BVH yet —
this task produces the boxes and proves they are sound.

**Files:**
- Modify: `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`
- Modify: `Detectors/Base/src/O2FlatCSG.cxx`
- Modify: `Detectors/Base/test/testFlatCSG.cxx`

**Interfaces:**
- Consumes: `EvalHalfspace`, `CellContains` (Tasks 1–3).
- Produces:
  - `struct o2::base::FlatCSGBox { double min[3]; double max[3]; int cell; int firstActive; int nActive; }` — `nActive == 0` means the box is wholly inside its cell
  - `static void O2FlatCSG::HalfspaceRange(const FlatCSGHalfspace&, const double lo[3], const double hi[3], double& rangeLo, double& rangeHi)` — a rigorous enclosure of `sign*f` over the box
  - `void O2FlatCSG::SetCellBBox(int cell, const double lo[3], const double hi[3])` — the converter supplies each cell's AABB, since the halfspaces alone do not bound it
  - `void O2FlatCSG::CloseShape()` — builds the boxes (and, from Task 5, the BVH)
  - `int O2FlatCSG::GetNboxes() const`, `const FlatCSGBox& O2FlatCSG::GetBox(int) const`
  - `void O2FlatCSG::SetSplitDepth(int)`, `void O2FlatCSG::SetMinBoxFraction(double)`

- [ ] **Step 1: Write the failing test**

Append to `Detectors/Base/test/testFlatCSG.cxx`:

```cpp
BOOST_AUTO_TEST_CASE(the_range_bound_encloses_the_sampled_range)
{
  // the bound must be an ENCLOSURE: over-wide is safe, under-wide is a wrong solid
  O2FlatCSG solid("range");
  double coeff[10];
  zCylinderQuadric(5., coeff);
  const int cylinder = solid.AddQuadric(1., coeff);
  const double normal[3] = {0., 0., 1.};
  const double through[3] = {0., 0., 2.};
  planeQuadric(normal, through, coeff);
  const int plane = solid.AddQuadric(-1., coeff);
  const double centre[3] = {1., 0., 0.};
  const double axis[3] = {0., 0., 1.};
  const int torus = solid.AddTorus(1., centre, axis, 7., 2.);

  Rng rng(555ULL);
  for (int trial = 0; trial < 3000; ++trial) {
    double lo[3];
    double hi[3];
    for (int index = 0; index < 3; ++index) {
      const double a = rng.uniform(-12., 12.);
      const double b = a + rng.uniform(0.01, 6.);
      lo[index] = a;
      hi[index] = b;
    }
    for (int which : {cylinder, plane, torus}) {
      double rangeLo = 0.;
      double rangeHi = 0.;
      O2FlatCSG::HalfspaceRange(solid.GetHalfspace(which), lo, hi, rangeLo, rangeHi);
      BOOST_REQUIRE_LE(rangeLo, rangeHi);
      for (int sample = 0; sample < 200; ++sample) {
        const double point[3] = {rng.uniform(lo[0], hi[0]), rng.uniform(lo[1], hi[1]),
                                 rng.uniform(lo[2], hi[2])};
        const double value = O2FlatCSG::EvalHalfspace(solid.GetHalfspace(which), point);
        BOOST_REQUIRE_GE(value, rangeLo - 1.e-9);
        BOOST_REQUIRE_LE(value, rangeHi + 1.e-9);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(the_boxes_cover_the_solid_and_their_active_lists_are_sound)
{
  // rmin = 2, rmax = 5, dz = 7 again, so there is a bore for the boxes to carve around
  O2FlatCSG solid("boxes");
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
  const double lo[3] = {-5., -5., -7.};
  const double hi[3] = {5., 5., 7.};
  solid.SetCellBBox(0, lo, hi);
  solid.CloseShape();

  BOOST_CHECK_GT(solid.GetNboxes(), 1);

  Rng rng(8080ULL);
  int insideSamples = 0;
  for (int trial = 0; trial < 50000; ++trial) {
    const double point[3] = {rng.uniform(-6., 6.), rng.uniform(-6., 6.), rng.uniform(-8., 8.)};
    if (!solid.Contains_Loop(point)) {
      continue;
    }
    ++insideSamples;
    // COVERAGE: every point of the solid is in some box
    bool covered = false;
    for (int index = 0; index < solid.GetNboxes() && !covered; ++index) {
      const auto& box = solid.GetBox(index);
      covered = point[0] >= box.min[0] && point[0] <= box.max[0] && point[1] >= box.min[1] &&
                point[1] <= box.max[1] && point[2] >= box.min[2] && point[2] <= box.max[2];
    }
    BOOST_REQUIRE(covered);
  }
  BOOST_CHECK_GT(insideSamples, 5000);

  // SOUNDNESS of the active lists: in every box, the active list alone decides membership
  for (int index = 0; index < solid.GetNboxes(); ++index) {
    const auto& box = solid.GetBox(index);
    for (int sample = 0; sample < 200; ++sample) {
      const double point[3] = {rng.uniform(box.min[0], box.max[0]),
                               rng.uniform(box.min[1], box.max[1]),
                               rng.uniform(box.min[2], box.max[2])};
      bool byActive = true;
      for (int slot = 0; slot < box.nActive && byActive; ++slot) {
        byActive = O2FlatCSG::EvalHalfspace(
                     solid.GetHalfspace(solid.GetActive(box.firstActive + slot)), point) <= 0.;
      }
      BOOST_REQUIRE_EQUAL(byActive, solid.CellContainsPublic(box.cell, point));
    }
  }
}

BOOST_AUTO_TEST_CASE(a_box_wholly_inside_a_cell_carries_no_active_halfspaces)
{
  O2FlatCSG solid("solid_boxes");
  addBoxCell(solid, 4., 4., 4.);
  const double lo[3] = {-4., -4., -4.};
  const double hi[3] = {4., 4., 4.};
  solid.SetCellBBox(0, lo, hi);
  solid.CloseShape();
  // a box has six planes and is convex, so subdivision must find interior boxes with an empty list
  int solidBoxes = 0;
  for (int index = 0; index < solid.GetNboxes(); ++index) {
    if (solid.GetBox(index).nActive == 0) {
      ++solidBoxes;
    }
  }
  BOOST_CHECK_GT(solidBoxes, 0);
}
```

The test uses two accessors that exist only to let it look inside: add
`int GetActive(int index) const { return fActive[index]; }` and
`bool CellContainsPublic(int cell, const double* point) const { return CellContains(cell, point); }`
to the public section, both documented as `/// For the tests: the box structure is the thing being
proved sound, so it has to be readable.`

- [ ] **Step 2: Run the test to verify it fails**

Run the Task 2 build command. Expected: FAIL — `no member named 'HalfspaceRange'`.

- [ ] **Step 3: Implement the range bound**

In the header, next to `EvalHalfspace`:

```cpp
  /// A rigorous enclosure `[rangeLo, rangeHi]` of `sign * f` over the box `[lo, hi]`.
  ///
  /// Conservative in the only safe direction: an over-wide enclosure loses pruning, never
  /// correctness. For a quadric it is the centred form -- with `m` the centre, `h` the
  /// half-extents and `g = A m + b`, `|Q(x) - Q(m)| <= 2 sum|g_i| h_i + sum |A_ij| h_i h_j` --
  /// which tightens quadratically as the boxes shrink. For a torus it is the 1-Lipschitz signed
  /// distance, so the enclosure is `f(m) +/- |h|`.
  static void HalfspaceRange(const FlatCSGHalfspace& halfspace, const double* lo, const double* hi,
                             double& rangeLo, double& rangeHi);
```

In the `.cxx`:

```cpp
void O2FlatCSG::HalfspaceRange(const FlatCSGHalfspace& halfspace, const double* lo,
                               const double* hi, double& rangeLo, double& rangeHi)
{
  double centre[3];
  double half[3];
  for (int index = 0; index < 3; ++index) {
    centre[index] = 0.5 * (lo[index] + hi[index]);
    half[index] = 0.5 * (hi[index] - lo[index]);
  }
  const double middle = EvalHalfspace(halfspace, centre);

  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    const double reach = std::sqrt(half[0] * half[0] + half[1] * half[1] + half[2] * half[2]);
    rangeLo = middle - reach;
    rangeHi = middle + reach;
    return;
  }

  const double* c = halfspace.c;
  const double a[3][3] = {{c[0], c[1], c[2]}, {c[1], c[3], c[4]}, {c[2], c[4], c[5]}};
  const double b[3] = {c[6], c[7], c[8]};
  double slack = 0.;
  for (int row = 0; row < 3; ++row) {
    double gradient = b[row];
    for (int column = 0; column < 3; ++column) {
      gradient += a[row][column] * centre[column];
      slack += std::abs(a[row][column]) * half[row] * half[column];
    }
    slack += 2. * std::abs(gradient) * half[row];
  }
  // the sign multiplies f, so it multiplies the deviation too -- but the slack is symmetric
  rangeLo = middle - slack;
  rangeHi = middle + slack;
}
```

- [ ] **Step 4: Implement the box builder**

Add to the header:

```cpp
/// One axis-aligned box of the sub-cell subdivision.
///
/// `nActive == 0` means every halfspace of `cell` holds everywhere in the box, so the box is
/// wholly inside the solid. A box is only ever created when the cell might intersect it: a box
/// the range bound proves wholly outside is dropped.
///
/// An active list is a correct description of the cell ONLY INSIDE ITS BOX. Every ray query must
/// clip to the box's slab before using it; see scripts/geometry/Design_FlatCSGSolid.md 4.4.
struct FlatCSGBox {
  double min[3] = {};
  double max[3] = {};
  int cell = -1;
  int firstActive = 0;
  int nActive = 0;
};
```

and, in the class:

```cpp
  /// The AABB of cell \a cell. The halfspaces alone do not bound a cell -- an intersection of
  /// halfspaces can be unbounded -- so the converter supplies the box the decomposition measured.
  void SetCellBBox(int cell, const double* lo, const double* hi);

  /// Build the sub-cell boxes (and, from task 5, the BVH). Call once, after the last AddCell.
  void CloseShape();
  bool IsClosed() const { return fClosed; }

  int GetNboxes() const { return static_cast<int>(fBoxes.size()); }
  const FlatCSGBox& GetBox(int index) const { return fBoxes[index]; }
  int GetActive(int index) const { return fActive[index]; }
  /// For the tests: the box structure is the thing being proved sound, so it has to be readable.
  bool CellContainsPublic(int cell, const double* point) const { return CellContains(cell, point); }

  /// Subdivision depth cap. Default 6; task 10 measures where it belongs.
  void SetSplitDepth(int depth) { fSplitDepth = depth; }
  /// Stop splitting a box narrower than this fraction of the part's bounding-box diagonal.
  /// Default 0.01; task 10 measures where it belongs.
  void SetMinBoxFraction(double fraction) { fMinBoxFraction = fraction; }
```

with the private members `std::vector<FlatCSGBox> fBoxes;`, `std::vector<int> fActive;`,
`std::vector<double> fCellLo, fCellHi;` (three doubles per cell), `bool fClosed = false;`,
`int fSplitDepth = 6;`, `double fMinBoxFraction = 0.01;`, and a private
`void SplitBox(int cell, const double* lo, const double* hi, const std::vector<int>& active, int depth, double minSize);`.

In the `.cxx`:

```cpp
void O2FlatCSG::SetCellBBox(int cell, const double* lo, const double* hi)
{
  if (static_cast<int>(fCellLo.size()) < 3 * GetNcells()) {
    fCellLo.resize(3 * GetNcells(), 0.);
    fCellHi.resize(3 * GetNcells(), 0.);
  }
  for (int index = 0; index < 3; ++index) {
    fCellLo[3 * cell + index] = lo[index];
    fCellHi[3 * cell + index] = hi[index];
  }
}

void O2FlatCSG::SplitBox(int cell, const double* lo, const double* hi,
                         const std::vector<int>& active, int depth, double minSize)
{
  std::vector<int> stillActive;
  stillActive.reserve(active.size());
  for (int halfspace : active) {
    double rangeLo = 0.;
    double rangeHi = 0.;
    HalfspaceRange(fHalfspaces[halfspace], lo, hi, rangeLo, rangeHi);
    if (rangeLo > 0.) {
      return; // the box is wholly outside this halfspace, hence wholly outside the cell
    }
    if (rangeHi > 0.) {
      stillActive.push_back(halfspace); // undecided; it stays
    }
    // rangeHi <= 0: the halfspace holds everywhere in the box, so it is dropped
  }

  double longest = 0.;
  int axis = 0;
  for (int index = 0; index < 3; ++index) {
    if (hi[index] - lo[index] > longest) {
      longest = hi[index] - lo[index];
      axis = index;
    }
  }
  const bool keep = stillActive.empty() || depth <= 0 || longest <= minSize;
  if (keep) {
    FlatCSGBox box;
    for (int index = 0; index < 3; ++index) {
      box.min[index] = lo[index];
      box.max[index] = hi[index];
    }
    box.cell = cell;
    box.firstActive = static_cast<int>(fActive.size());
    box.nActive = static_cast<int>(stillActive.size());
    fActive.insert(fActive.end(), stillActive.begin(), stillActive.end());
    fBoxes.push_back(box);
    return;
  }

  const double middle = 0.5 * (lo[axis] + hi[axis]);
  double childLo[3] = {lo[0], lo[1], lo[2]};
  double childHi[3] = {hi[0], hi[1], hi[2]};
  childHi[axis] = middle;
  SplitBox(cell, childLo, childHi, stillActive, depth - 1, minSize);
  childHi[axis] = hi[axis];
  childLo[axis] = middle;
  SplitBox(cell, childLo, childHi, stillActive, depth - 1, minSize);
}

void O2FlatCSG::CloseShape()
{
  fBoxes.clear();
  fActive.clear();
  fCellLo.resize(3 * GetNcells(), 0.);
  fCellHi.resize(3 * GetNcells(), 0.);

  double partLo[3] = {TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};
  double partHi[3] = {-TGeoShape::Big(), -TGeoShape::Big(), -TGeoShape::Big()};
  for (int cell = 0; cell < GetNcells(); ++cell) {
    for (int index = 0; index < 3; ++index) {
      partLo[index] = std::min(partLo[index], fCellLo[3 * cell + index]);
      partHi[index] = std::max(partHi[index], fCellHi[3 * cell + index]);
    }
  }
  const double diagonal = std::sqrt((partHi[0] - partLo[0]) * (partHi[0] - partLo[0]) +
                                    (partHi[1] - partLo[1]) * (partHi[1] - partLo[1]) +
                                    (partHi[2] - partLo[2]) * (partHi[2] - partLo[2]));
  const double minSize = fMinBoxFraction * diagonal;

  for (int cell = 0; cell < GetNcells(); ++cell) {
    std::vector<int> active;
    active.reserve(fCells[cell].count);
    for (int offset = 0; offset < fCells[cell].count; ++offset) {
      active.push_back(fCells[cell].first + offset);
    }
    SplitBox(cell, &fCellLo[3 * cell], &fCellHi[3 * cell], active, fSplitDepth, minSize);
  }
  fClosed = true;
  ComputeBBox();
}
```

`ComputeBBox` arrives properly in Task 6; for this task add a minimal override setting `fDX/fDY/
fDZ/fOrigin` from the union of `fBoxes`, so `CloseShape` compiles.

- [ ] **Step 5: Run the tests to verify they pass**

Run the Task 2 test command. Expected: PASS, 11 cases.

- [ ] **Step 6: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2FlatCSG.h \
  Detectors/Base/src/O2FlatCSG.cxx Detectors/Base/test/testFlatCSG.cxx
git commit -m "$(cat <<'EOF'
Subdivide each cell into boxes with a rigorous range bound

This adds the sub-cell subdivision the flat solid is accelerated on, one level
below the cells.

- A quadric's range over a box comes from the centred form and tightens
  quadratically as boxes shrink; a torus uses its 1-Lipschitz signed distance.
- A box the bound proves wholly outside a halfspace is dropped, and a halfspace
  that holds everywhere in a box is dropped from that box's list.
- A box whose list empties is wholly inside the cell.
- The tests prove coverage and that a box's active list alone decides membership
  inside that box.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: The BVH, and the accelerated queries

Spec §4.4 is the whole risk of this task: an active list is valid only inside its box, so every
ray query must clip to the box's slab before using it. The `_Loop` twins from Tasks 1–3 are what
prove it.

**Files:**
- Modify: `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`
- Modify: `Detectors/Base/src/O2FlatCSG.cxx`
- Modify: `Detectors/Base/test/testFlatCSG.cxx`

**Interfaces:**
- Consumes: `FlatCSGBox`, `CloseShape`, `CellIntervals`, all three `_Loop` twins.
- Produces: accelerated `Contains`, `DistFromOutside_Loop`-identical accelerated `DistFromOutside`
  and `DistFromInside`; `size_t O2FlatCSG::GetBVHMemory() const`.

- [ ] **Step 1: Write the failing test**

Append to `Detectors/Base/test/testFlatCSG.cxx` a helper that builds a deliberately awkward solid
and then demands bit identity:

```cpp
namespace
{
/// An L-shaped bracket with a bore: three cells, a complemented cylinder, a long diagonal extent.
/// Deliberately the shape a cell-level BVH would handle badly.
void buildBracket(O2FlatCSG& solid)
{
  double coeff[10];
  // cell 0: the long arm, x in [-10, 10], y in [-1, 1], z in [-1, 1]
  const double arm[6][2][3] = {{{1., 0., 0.}, {10., 0., 0.}},   {{-1., 0., 0.}, {-10., 0., 0.}},
                               {{0., 1., 0.}, {0., 1., 0.}},    {{0., -1., 0.}, {0., -1., 0.}},
                               {{0., 0., 1.}, {0., 0., 1.}},    {{0., 0., -1.}, {0., 0., -1.}}};
  int first = solid.GetNhalfspaces();
  for (const auto& plane : arm) {
    planeQuadric(plane[0], plane[1], coeff);
    solid.AddQuadric(1., coeff);
  }
  solid.AddCell(first, 6, 8. * 10. * 1. * 1.);
  const double armLo[3] = {-10., -1., -1.};
  const double armHi[3] = {10., 1., 1.};
  solid.SetCellBBox(0, armLo, armHi);

  // cell 1: the upright, x in [8, 10], y in [-1, 1], z in [1, 12]
  const double upright[6][2][3] = {{{1., 0., 0.}, {10., 0., 0.}}, {{-1., 0., 0.}, {8., 0., 0.}},
                                   {{0., 1., 0.}, {0., 1., 0.}},  {{0., -1., 0.}, {0., -1., 0.}},
                                   {{0., 0., 1.}, {0., 0., 12.}}, {{0., 0., -1.}, {0., 0., 1.}}};
  first = solid.GetNhalfspaces();
  for (const auto& plane : upright) {
    planeQuadric(plane[0], plane[1], coeff);
    solid.AddQuadric(1., coeff);
  }
  solid.AddCell(first, 6, 2. * 2. * 11.);
  const double uprightLo[3] = {8., -1., 1.};
  const double uprightHi[3] = {10., 1., 12.};
  solid.SetCellBBox(1, uprightLo, uprightHi);

  // cell 2: a washer around z at x = -8, with a bore -- a complemented cylinder, so non-convex
  first = solid.GetNhalfspaces();
  const double centreShift = -8.;
  // outer cylinder about the axis through (-8, 0, *): translate by completing the square
  const double outer[10] = {1., 0., 0., 1., 0., 0., -centreShift, 0., 0.,
                            centreShift * centreShift - 9.};
  solid.AddQuadric(1., outer);
  const double inner[10] = {1., 0., 0., 1., 0., 0., -centreShift, 0., 0.,
                            centreShift * centreShift - 1.};
  solid.AddQuadric(-1., inner);
  const double washer[2][2][3] = {{{0., 0., 1.}, {0., 0., 1.}}, {{0., 0., -1.}, {0., 0., -1.}}};
  for (const auto& plane : washer) {
    planeQuadric(plane[0], plane[1], coeff);
    solid.AddQuadric(1., coeff);
  }
  solid.AddCell(first, 4, TMath::Pi() * (9. - 1.) * 2.);
  const double washerLo[3] = {-11., -3., -1.};
  const double washerHi[3] = {-5., 3., 1.};
  solid.SetCellBBox(2, washerLo, washerHi);
}
} // namespace

BOOST_AUTO_TEST_CASE(the_accelerated_contains_is_bit_identical_to_its_twin)
{
  O2FlatCSG solid("bracket");
  buildBracket(solid);
  solid.CloseShape();
  BOOST_CHECK_GT(solid.GetNboxes(), 3);
  BOOST_CHECK_GT(solid.GetBVHMemory(), 0u);

  Rng rng(123456ULL);
  for (int trial = 0; trial < 200000; ++trial) {
    const double point[3] = {rng.uniform(-13., 13.), rng.uniform(-5., 5.), rng.uniform(-3., 14.)};
    BOOST_REQUIRE_EQUAL(solid.Contains(point), solid.Contains_Loop(point));
  }
}

BOOST_AUTO_TEST_CASE(the_accelerated_distances_are_bit_identical_to_their_twins)
{
  O2FlatCSG solid("bracket_dist");
  buildBracket(solid);
  solid.CloseShape();

  Rng rng(654321ULL);
  for (int trial = 0; trial < 200000; ++trial) {
    double point[3] = {rng.uniform(-16., 16.), rng.uniform(-8., 8.), rng.uniform(-6., 17.)};
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
    if (solid.Contains_Loop(point)) {
      BOOST_REQUIRE_EQUAL(solid.DistFromInside(point, dir, 3, TGeoShape::Big(), nullptr),
                          solid.DistFromInside_Loop(point, dir, TGeoShape::Big()));
    } else {
      BOOST_REQUIRE_EQUAL(solid.DistFromOutside(point, dir, 3, TGeoShape::Big(), nullptr),
                          solid.DistFromOutside_Loop(point, dir, TGeoShape::Big()));
    }
  }
}

BOOST_AUTO_TEST_CASE(a_ray_along_the_long_arm_crosses_every_cell_it_should)
{
  // the case a per-box clip gets wrong if it forgets to clip: a ray running the length of the
  // bracket passes through many boxes of the same cell, and must see ONE interval, not many
  O2FlatCSG solid("bracket_long");
  buildBracket(solid);
  solid.CloseShape();
  const double origin[3] = {-20., 0., 0.};
  const double dir[3] = {1., 0., 0.};
  BOOST_CHECK_SMALL(solid.DistFromOutside(origin, dir, 3, TGeoShape::Big(), nullptr) - 9., 1.e-12);
  const double inArm[3] = {0., 0., 0.};
  // inside the arm at the origin, the exit is x = 10 (the arm and the upright touch at x = 8..10
  // only for z > 1, so along z = 0 the arm alone decides)
  BOOST_CHECK_SMALL(solid.DistFromInside(inArm, dir, 3, TGeoShape::Big(), nullptr) - 10., 1.e-12);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run the Task 2 build command. Expected: FAIL — `no member named 'GetBVHMemory'`.

- [ ] **Step 3: Implement**

In the header, add to the class:

```cpp
  /// Bytes held by the BVH nodes and the primitive-index permutation.
  size_t GetBVHMemory() const;
```

with the private member `void* fBVH = nullptr; //! bvh::v2::Bvh over the sub-cell boxes` — the
`//!` matters, the BVH is rebuilt on load and never streamed (spec §7). Declare a copy constructor
and assignment operator as `= delete`, exactly as `O2BVHAssembly` does, because of the raw pointer.

In the `.cxx`, follow `O2BVHAssembly.cxx`'s pattern verbatim: include `bvh2_third_party.h` and
`bvh2_extra_kernels.h`, define the `BVH`, `BVHNode`, `BVHBBox`, `BVHVec3` aliases the same way,
and at the end of `CloseShape()`:

```cpp
  delete static_cast<BVH*>(fBVH);
  fBVH = nullptr;
  if (!fBoxes.empty()) {
    std::vector<BVHBBox> boxes;
    std::vector<BVHVec3> centers;
    boxes.reserve(fBoxes.size());
    centers.reserve(fBoxes.size());
    for (const auto& box : fBoxes) {
      BVHBBox bounds;
      for (int index = 0; index < 3; ++index) {
        bounds.min[index] = box.min[index];
        bounds.max[index] = box.max[index];
      }
      boxes.push_back(bounds);
      centers.emplace_back(bounds.get_center());
    }
    typename bvh::v2::DefaultBuilder<BVHNode>::Config config;
    config.quality = bvh::v2::DefaultBuilder<BVHNode>::Quality::High;
    // one box per leaf: resolving a box costs a slab clip plus an interval scan, far more than a
    // node box test, and the bvh2 traversal enters a leaf's start node without testing its box
    config.max_leaf_size = 1;
    fBVH = static_cast<void*>(new BVH(bvh::v2::DefaultBuilder<BVHNode>::build(boxes, centers, config)));
  }
```

`Contains` becomes a BVH point query: descend the stack testing `bvh::v2::extra::contains` on each
node box, and at a leaf take the box, answer `kTRUE` immediately if `box.nActive == 0`, otherwise
evaluate the box's active list and answer `kTRUE` on the first box that accepts.

The distances become a traversal that visits boxes and, **for each box**, clips the ray to that
box's slab and calls `CellIntervals(box.cell, &fActive[box.firstActive], box.nActive, ...)` with
`tlo`/`thi` set to the slab entry and exit. Write the slab clip as a small helper:

```cpp
/// The ray's parameter window inside \a box; returns false when it misses.
bool slab(const FlatCSGBox& box, const double* origin, const double* dir, double& tlo, double& thi)
{
  tlo = 0.;
  thi = TGeoShape::Big();
  for (int index = 0; index < 3; ++index) {
    if (std::abs(dir[index]) < 1.e-300) {
      if (origin[index] < box.min[index] || origin[index] > box.max[index]) {
        return false;
      }
      continue;
    }
    const double inverse = 1. / dir[index];
    double near = (box.min[index] - origin[index]) * inverse;
    double far = (box.max[index] - origin[index]) * inverse;
    if (near > far) {
      std::swap(near, far);
    }
    tlo = std::max(tlo, near);
    thi = std::min(thi, far);
    if (tlo > thi) {
      return false;
    }
  }
  return true;
}
```

For `DistFromOutside`, gather every box the ray meets, sort the intervals it yields, and take the
first entry — the twin's answer. For `DistFromInside`, gather the intervals from every box the ray
meets, merge them with `mergeIntervals`, and return the far end of the one containing `t = 0`.

**A warning for the implementer, and it is the point of the task:** merging must happen across
boxes *after* each box has produced its intervals from its own clipped window, never by
concatenating active lists across boxes. Boxes of the same cell tile a region, so one cell's
occupancy arrives in several adjacent pieces that `mergeIntervals` rejoins — that is why the glue
tolerance exists and why `a_ray_along_the_long_arm_crosses_every_cell_it_should` is in the test
set.

If bit identity fails, the twin is right and the traversal is wrong; do not relax the test.

- [ ] **Step 4: Run the tests to verify they pass**

Run the Task 2 test command. Expected: PASS, 14 cases.

- [ ] **Step 5: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2FlatCSG.h \
  Detectors/Base/src/O2FlatCSG.cxx Detectors/Base/test/testFlatCSG.cxx
git commit -m "$(cat <<'EOF'
Accelerate the flat solid with a BVH over its sub-cell boxes

This puts the sub-cell boxes into a BVH and routes the queries through it.

- A ray is clipped to a box's slab before that box's active list is used, since
  an active list describes the cell only inside its own box.
- One box per leaf, as in O2BVHAssembly: resolving a box costs more than a node
  box test.
- The tests require the accelerated queries to be bit-identical to the loop
  twins over 200 000 points and 200 000 rays on an L-bracket with a bore.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: `Safety`, `Capacity`, `ComputeBBox`, `ComputeNormal`

The rest of the `TGeoShape` contract. `Safety` comes entirely from the box structure (spec §5.4),
which is why it waited for Task 5.

**Files:**
- Modify: `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`
- Modify: `Detectors/Base/src/O2FlatCSG.cxx`
- Modify: `Detectors/Base/test/testFlatCSG.cxx`

**Interfaces:**
- Consumes: `fBoxes`, `fBVH`, `fCells[].volume`.
- Produces: `Safety`, `Safety_Loop`, `Capacity`, `ComputeBBox`, `ComputeNormal` overrides.

- [ ] **Step 1: Write the failing test**

```cpp
BOOST_AUTO_TEST_CASE(safety_is_sound_and_matches_its_twin)
{
  O2FlatCSG solid("bracket_safety");
  buildBracket(solid);
  solid.CloseShape();

  Rng rng(24680ULL);
  for (int trial = 0; trial < 50000; ++trial) {
    double point[3] = {rng.uniform(-16., 16.), rng.uniform(-8., 8.), rng.uniform(-6., 17.)};
    const bool inside = solid.Contains_Loop(point);
    const double safety = solid.Safety(point, inside);
    BOOST_REQUIRE_GE(safety, 0.);
    BOOST_REQUIRE_EQUAL(safety, solid.Safety_Loop(point, inside));

    // SOUNDNESS: no point within `safety` of `point` may have the opposite classification
    for (int probe = 0; probe < 40; ++probe) {
      double dir[3];
      double norm = 0.;
      do {
        for (int index = 0; index < 3; ++index) {
          dir[index] = rng.uniform(-1., 1.);
        }
        norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
      } while (norm < 1.e-3);
      const double reach = safety * rng.uniform(0., 0.999) / norm;
      const double near[3] = {point[0] + reach * dir[0], point[1] + reach * dir[1],
                              point[2] + reach * dir[2]};
      BOOST_REQUIRE_EQUAL(static_cast<bool>(solid.Contains_Loop(near)), inside);
    }
  }
}

BOOST_AUTO_TEST_CASE(capacity_is_the_sum_of_the_cell_volumes)
{
  O2FlatCSG solid("bracket_capacity");
  buildBracket(solid);
  solid.CloseShape();
  const double expected = 8. * 10. * 1. * 1. + 2. * 2. * 11. + TMath::Pi() * (9. - 1.) * 2.;
  BOOST_CHECK_SMALL(solid.Capacity() - expected, 1.e-12);
}

BOOST_AUTO_TEST_CASE(the_bounding_box_is_tight_around_the_retained_boxes)
{
  O2FlatCSG solid("bracket_bbox");
  buildBracket(solid);
  solid.CloseShape();
  // the bracket spans x in [-11, 10], y in [-3, 3], z in [-1, 12]
  BOOST_CHECK_SMALL(solid.GetDX() - 10.5, 1.e-9);
  BOOST_CHECK_SMALL(solid.GetDY() - 3., 1.e-9);
  BOOST_CHECK_SMALL(solid.GetDZ() - 6.5, 1.e-9);
  BOOST_CHECK_SMALL(solid.GetOrigin()[0] + 0.5, 1.e-9);
  BOOST_CHECK_SMALL(solid.GetOrigin()[2] - 5.5, 1.e-9);

  // and every point of the solid is inside it
  Rng rng(13579ULL);
  for (int trial = 0; trial < 50000; ++trial) {
    const double point[3] = {rng.uniform(-13., 13.), rng.uniform(-5., 5.), rng.uniform(-3., 14.)};
    if (solid.Contains_Loop(point)) {
      BOOST_REQUIRE(solid.TGeoBBox::Contains(point));
    }
  }
}

BOOST_AUTO_TEST_CASE(the_normal_on_a_face_is_the_face_normal)
{
  O2FlatCSG solid("box_normal");
  addBoxCell(solid, 3., 4., 5.);
  const double lo[3] = {-3., -4., -5.};
  const double hi[3] = {3., 4., 5.};
  solid.SetCellBBox(0, lo, hi);
  solid.CloseShape();

  const double onFace[3] = {3., 1., 1.};
  const double dir[3] = {1., 0., 0.};
  double normal[3] = {0., 0., 0.};
  solid.ComputeNormal(onFace, dir, normal);
  BOOST_CHECK_SMALL(normal[0] - 1., 1.e-12);
  BOOST_CHECK_SMALL(normal[1], 1.e-12);
  BOOST_CHECK_SMALL(normal[2], 1.e-12);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run the Task 2 build command. Expected: FAIL — `no member named 'Safety_Loop'`, and the
`Capacity`/`ComputeNormal` cases fail on the placeholders.

- [ ] **Step 3: Implement**

Replace the Task 2 `Safety` placeholder. Both `Safety` and `Safety_Loop` compute the same number;
the twin walks all boxes and the accelerated one uses the BVH with `bvh::v2::extra::SafetySqToNode`
to prune, so bit identity means the pruning is not losing a nearer box.

```cpp
Double_t O2FlatCSG::Safety_Loop(const Double_t* point, Bool_t in) const
{
  if (!in) {
    // outside: every point of the solid is in some box, so the distance to the nearest box is a
    // lower bound on the distance to the solid
    double best = TGeoShape::Big();
    for (const auto& box : fBoxes) {
      double squared = 0.;
      for (int index = 0; index < 3; ++index) {
        const double value = point[index];
        if (value < box.min[index]) {
          squared += (box.min[index] - value) * (box.min[index] - value);
        } else if (value > box.max[index]) {
          squared += (value - box.max[index]) * (value - box.max[index]);
        }
      }
      best = std::min(best, squared);
    }
    return best >= TGeoShape::Big() ? 0. : std::sqrt(best);
  }

  // inside: the box holding the point bounds how far the answer can be trusted. A box with an
  // empty active list is wholly inside the cell, so the distance to its faces is a sound bound;
  // an undecided box gives the same bound, which is still sound because the box's own faces are
  // reached no later than the solid's boundary within it.
  double best = 0.;
  for (const auto& box : fBoxes) {
    bool holds = true;
    for (int index = 0; index < 3 && holds; ++index) {
      holds = point[index] >= box.min[index] && point[index] <= box.max[index];
    }
    if (!holds || box.nActive != 0) {
      continue;
    }
    double toFace = TGeoShape::Big();
    for (int index = 0; index < 3; ++index) {
      toFace = std::min(toFace, std::min(point[index] - box.min[index], box.max[index] - point[index]));
    }
    best = std::max(best, toFace);
  }
  return std::max(best, 0.);
}
```

A note the implementer must respect: the inside branch deliberately considers **only** boxes with
an empty active list. An undecided box may contain boundary anywhere, so its face distance is not
a bound at all; returning `0.` for a point in one is the sound answer, and Task 10 measures what
that costs.

`Capacity` is `std::accumulate` over `fCells[].volume`. `ComputeBBox` takes the union of
`fBoxes` — tighter than the union of cell AABBs, which is the point of §5.5 — and sets
`fDX/fDY/fDZ/fOrigin`. `ComputeNormal` finds the cell whose halfspace is closest to zero at the
point and returns the normalised `sign * 2(Ax + b)` for a quadric, or the normalised gradient of
the torus's signed distance, flipped so it points along `dir`, matching the `TGeoShape` contract
(the normal is oriented outward with respect to the direction of travel).

- [ ] **Step 4: Run the tests to verify they pass**

Run the Task 2 test command. Expected: PASS, 18 cases.

- [ ] **Step 5: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2FlatCSG.h \
  Detectors/Base/src/O2FlatCSG.cxx Detectors/Base/test/testFlatCSG.cxx
git commit -m "$(cat <<'EOF'
Complete the flat solid's TGeoShape contract

This adds Safety, Capacity, ComputeBBox and ComputeNormal to O2FlatCSG.

- Safety comes from the box structure alone, so no point-to-quadric distance
  formula is needed; a point in an undecided box gets the sound answer of zero.
- Capacity is the sum of the cells' own volumes, which the converter measures on
  the source piece; the cells are disjoint.
- ComputeBBox takes the union of the retained boxes, which is tighter than the
  union of the cell bounding boxes.
- The safety test checks soundness directly: no point within the returned radius
  may be classified the other way.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: The sidecar format, the loader, the streamer and the dictionary

Spec §7. The shape reaches `geom.C` through a binary sidecar, as the surface path already does.

**Files:**
- Modify: `Detectors/Base/include/DetectorsBase/O2SurfaceSolidIO.h`
- Modify: `Detectors/Base/src/O2SurfaceSolidIO.cxx`
- Modify: `Detectors/Base/src/DetectorsBaseLinkDef.h`
- Modify: `Detectors/Base/test/testFlatCSG.cxx`
- Modify: `scripts/geometry/BVHSurfaceSolid.md` (document the format next to the surface one)

**Interfaces:**
- Consumes: `AddQuadric`, `AddTorus`, `AddCell`, `SetCellBBox`, `CloseShape`.
- Produces: `bool o2::base::LoadFlatCSG(const std::string& file, O2FlatCSG& solid)`.

**The format, version 1, little-endian throughout:**

```
magic          char[8]   "O2FLTCSG"
version        uint32    1
nHalfspaces    uint32
nCells         uint32
halfspaces     nHalfspaces * { int32 kind; float64 sign; float64 c[11] }
cells          nCells     * { int32 first; int32 count; float64 volume;
                              float64 lo[3]; float64 hi[3] }
```

- [ ] **Step 1: Write the failing test**

```cpp
BOOST_AUTO_TEST_CASE(a_sidecar_round_trip_reproduces_the_solid)
{
  O2FlatCSG original("bracket_io");
  buildBracket(original);
  original.CloseShape();

  const std::string path = "testFlatCSG_roundtrip.bin";
  BOOST_REQUIRE(o2::base::WriteFlatCSG(path, original)); // test-only writer, see step 3

  O2FlatCSG loaded("bracket_io_loaded");
  BOOST_REQUIRE(o2::base::LoadFlatCSG(path, loaded));
  loaded.CloseShape();

  BOOST_CHECK_EQUAL(loaded.GetNhalfspaces(), original.GetNhalfspaces());
  BOOST_CHECK_EQUAL(loaded.GetNcells(), original.GetNcells());
  BOOST_CHECK_EQUAL(loaded.GetNboxes(), original.GetNboxes());
  BOOST_CHECK_EQUAL(loaded.Capacity(), original.Capacity());

  Rng rng(97531ULL);
  for (int trial = 0; trial < 100000; ++trial) {
    const double point[3] = {rng.uniform(-13., 13.), rng.uniform(-5., 5.), rng.uniform(-3., 14.)};
    BOOST_REQUIRE_EQUAL(loaded.Contains(point), original.Contains(point));
  }
  std::filesystem::remove(path);
}

BOOST_AUTO_TEST_CASE(a_truncated_sidecar_is_refused_rather_than_half_loaded)
{
  O2FlatCSG original("bracket_trunc");
  buildBracket(original);
  original.CloseShape();
  const std::string path = "testFlatCSG_truncated.bin";
  BOOST_REQUIRE(o2::base::WriteFlatCSG(path, original));
  std::filesystem::resize_file(path, std::filesystem::file_size(path) - 17);

  O2FlatCSG loaded("bracket_trunc_loaded");
  BOOST_CHECK(!o2::base::LoadFlatCSG(path, loaded));
  std::filesystem::remove(path);
}

BOOST_AUTO_TEST_CASE(the_shape_survives_a_ROOT_file_without_its_sidecar)
{
  O2FlatCSG original("bracket_root");
  buildBracket(original);
  original.CloseShape();

  const std::string path = "testFlatCSG_shape.root";
  {
    TFile file(path.c_str(), "RECREATE");
    file.WriteObject(&original, "shape");
  }
  O2FlatCSG* restored = nullptr;
  {
    TFile file(path.c_str(), "READ");
    file.GetObject("shape", restored);
  }
  BOOST_REQUIRE(restored != nullptr);
  restored->CloseShape(); // the BVH is not streamed; it is rebuilt

  Rng rng(11223ULL);
  for (int trial = 0; trial < 100000; ++trial) {
    const double point[3] = {rng.uniform(-13., 13.), rng.uniform(-5., 5.), rng.uniform(-3., 14.)};
    BOOST_REQUIRE_EQUAL(restored->Contains(point), original.Contains(point));
  }
  std::filesystem::remove(path);
}
```

Add `#include "TFile.h"`, `#include <filesystem>` and `#include "DetectorsBase/O2SurfaceSolidIO.h"`
at the top of the test.

- [ ] **Step 2: Run the test to verify it fails**

Run the Task 2 build command. Expected: FAIL — `no member named 'LoadFlatCSG' in namespace 'o2::base'`.

- [ ] **Step 3: Implement the loader and the writer**

Declare both in `O2SurfaceSolidIO.h`, next to `LoadFacetSolid`:

```cpp
class O2FlatCSG;

/// Load a flat-CSG sidecar (flatcsg_*.bin, version-1 format documented in
/// scripts/geometry/BVHSurfaceSolid.md, "Flat-CSG sidecar format") into \a solid by dispatching
/// to its AddQuadric / AddTorus / AddCell / SetCellBBox methods. The caller is expected to call
/// CloseShape() afterwards. Returns false (with an error message) on I/O or format problems; the
/// solid may then be partly filled and should be discarded.
bool LoadFlatCSG(const std::string& file, O2FlatCSG& solid);

/// Write \a solid in the same format. Used by the converter's tests and by the round-trip case;
/// the production writer is scripts/geometry/csg/flat.py, and the two must agree byte for byte.
bool WriteFlatCSG(const std::string& file, const O2FlatCSG& solid);
```

Implement both in `O2SurfaceSolidIO.cxx` following `LoadFacetSolid`'s error handling: check the
magic, check the version, and check that the remaining file length is exactly what the counts
imply **before** reading any record, so a truncated file is refused rather than half-loaded.
Validate `first >= 0`, `count > 0` and `first + count <= nHalfspaces` for every cell, and return
false naming the offending cell if not.

Add to `DetectorsBaseLinkDef.h`, after the `O2BVHAssembly` line:

```cpp
#pragma link C++ class o2::base::FlatCSGHalfspace + ;
#pragma link C++ class o2::base::FlatCSGCell + ;
#pragma link C++ class o2::base::FlatCSGBox + ;
#pragma link C++ class std::vector < o2::base::FlatCSGHalfspace> + ;
#pragma link C++ class std::vector < o2::base::FlatCSGCell> + ;
#pragma link C++ class std::vector < o2::base::FlatCSGBox> + ;
#pragma link C++ class o2::base::O2FlatCSG + ;
```

The `+` (rather than `O2BVHSurfaceSolid`'s `-`) is correct here because everything streamed is a
POD vector; `fBVH` is excluded by its `//!` comment.

- [ ] **Step 4: Document the format**

Add a "Flat-CSG sidecar format" section to `scripts/geometry/BVHSurfaceSolid.md` beside the
existing "Surface sidecar format" section, with the table from this task's header verbatim.

- [ ] **Step 5: Run the tests to verify they pass**

Run the Task 2 test command, then also `ctest -R 'FlatCSG|BVHSurfaceSolid|BVHAssembly'`.
Expected: PASS, 21 FlatCSG cases; `BVHSurfaceSolid` 113 and `BVHAssembly` 22 unchanged.

- [ ] **Step 6: Commit**

```bash
cd $HOME/alisw/O2 && git add Detectors/Base/include/DetectorsBase/O2SurfaceSolidIO.h \
  Detectors/Base/src/O2SurfaceSolidIO.cxx Detectors/Base/src/DetectorsBaseLinkDef.h \
  Detectors/Base/test/testFlatCSG.cxx scripts/geometry/BVHSurfaceSolid.md
git commit -m "$(cat <<'EOF'
Give the flat solid a sidecar format and a dictionary

This adds LoadFlatCSG and WriteFlatCSG next to the surface-solid loaders, and
the ROOT dictionary entries for the shape.

- The format is documented in BVHSurfaceSolid.md beside the surface one.
- A truncated file is refused before any record is read, rather than half-loaded.
- The BVH is not streamed; it is rebuilt by CloseShape after a read.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: The Python flat emitter, and the `_cell_leaf` bridge test

Spec §3.3 and §8. The riskiest thing in the rung: an inverted halfspace still produces a solid, so
the sign convention is measured against the corpus-validated `_cell_leaf` rather than trusted.

**Files:**
- Create: `scripts/geometry/csg/flat.py`
- Modify: `scripts/geometry/csg/emit.py` (self-test cases only in this task)

**Interfaces:**
- Consumes: `recognise._halfspace_carriers`, `recognise._cell_leaf`, `recognise._CellBox`, `primitives.build_occ`, `decompose.split_into_cells`.
- Produces:
  - `csg.flat.quadric_from_carrier(carrier) -> (sign, [10 floats])` — raises `recognise.Declined` for a carrier kind it cannot express
  - `csg.flat.torus_from_carrier(carrier) -> (sign, centre, axis, major, minor)`
  - `csg.flat.blocks_from_carriers(carriers) -> list[dict]` with each dict `{"kind": "quadric"|"torus", "sign": float, "c": [...]}`
  - `csg.flat.flat_contains(blocks, point) -> bool`
  - `csg.flat.write_sidecar(path, blocks, cells) -> None`, `cells` a list of `{"first", "count", "volume", "lo", "hi"}`
  - `csg.flat.SIDECAR_MAGIC = b"O2FLTCSG"`, `csg.flat.SIDECAR_VERSION = 1`

- [ ] **Step 1: Write the failing self-test cases**

Add to `csg/emit.py`'s `--self-test`, in the same `check(...)` style the file already uses:

```python
    # --- the flat emitter's sign convention, measured against the shipped cell emitter --------
    from csg import flat as flatmod

    def _flat_agrees_with_cell_leaf(solid, seed=20260824, samples=4000):
        """Do the flat halfspaces and _cell_leaf's padded primitives classify the same points?

        `_cell_leaf` has shipped 1775 known-source-clean parts, so it is the oracle for the one
        thing the flat emitter can get silently wrong: an inverted halfspace is still a solid.
        """
        import random
        tol = recognise._tolerance_for(solid)
        diag = decompose.bbox_diagonal(solid)
        carriers = recognise._halfspace_carriers(solid, tol)
        box = recognise._CellBox(solid, diag)
        blocks = flatmod.blocks_from_carriers(carriers)
        leaves = [recognise._cell_leaf(c, box) for c in carriers]
        cand = prim.cell("intersection" if len(leaves) > 1 else "primitive", leaves)
        padded = prim.build_occ(cand)
        rng = random.Random(seed)
        (xlo, ylo, zlo, xhi, yhi, zhi) = recognise._bbox_of(solid)
        disagreements = 0
        for _ in range(samples):
            point = (rng.uniform(xlo, xhi), rng.uniform(ylo, yhi), rng.uniform(zlo, zhi))
            if flatmod.flat_contains(blocks, point) != accept.occ_contains(padded, point):
                disagreements += 1
        return disagreements, samples

    for name, maker in (("a box", make_boolean_fixtures.box_solid),
                        ("a blind bore", make_boolean_fixtures.blind_bore_solid),
                        ("a grooved block", make_boolean_fixtures.grooved_block_solid),
                        ("a cone frustum", make_boolean_fixtures.cone_solid),
                        ("a torus ply", make_boolean_fixtures.torus_ply_solid)):
        bad, scored = _flat_agrees_with_cell_leaf(maker())
        check(f"the flat halfspaces of {name} classify exactly as _cell_leaf's primitives",
              bad == 0, f"{bad} of {scored} points disagree")
```

If a named fixture maker does not exist under that name in `make_boolean_fixtures.py`, use the
closest one that produces a single cell of that carrier kind and rename the check accordingly —
what matters is that all five carrier kinds (plane, cylinder, cone, sphere, torus) are covered by
at least one case, and that at least one of them uses an `exterior` carrier so a complemented sign
is exercised.

- [ ] **Step 2: Run the self-test to verify it fails**

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)" \
  && export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH \
  && SW=$HOME/alisw/sw/ubuntu2404_aarch64 \
  && export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH \
  && export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH \
  && cd $HOME/alisw/O2 && python3 scripts/geometry/csg/emit.py --self-test
```

Expected: FAIL — `ModuleNotFoundError: No module named 'csg.flat'`.

- [ ] **Step 3: Write `csg/flat.py`**

```python
"""The flat-DNF emitter: a cell's carriers as signed implicit halfspaces.

`recognise._halfspace_carriers` already produces exactly what this module needs -- an oriented
carrier per face with the material side decided -- and `recognise._cell_leaf` then *bounds* each
one into a padded native primitive, because ROOT has no usable halfspace and a
`TGeoCompositeShape` is the only thing the converter could emit. `O2FlatCSG` removes that
constraint, so this module maps a carrier to the halfspace it always was.

The sign is the whole risk
--------------------------
An inverted halfspace still produces a solid, so a sign error here is silent. Two things compose
into the single sign: the carrier's own orientation (a plane's normal is already flipped for
`TopAbs_REVERSED` upstream) and the `side` field. `csg/emit.py --self-test` measures the result
against `_cell_leaf`'s padded conjunction on sampled points, because `_cell_leaf` has shipped 1775
known-source-clean parts and is the strongest oracle available for it.
"""

import math
import struct

SIDECAR_MAGIC = b"O2FLTCSG"
SIDECAR_VERSION = 1


def _outer(u, v):
    return [[u[i] * v[j] for j in range(3)] for i in range(3)]


def _quadric(a, b, c):
    """Pack A (3x3 symmetric), b (3) and c into the ten-double block the shape stores."""
    return [a[0][0], a[0][1], a[0][2], a[1][1], a[1][2], a[2][2], b[0], b[1], b[2], c]


def quadric_from_carrier(carrier):
    """`(sign, block)` for a plane, sphere, cylinder, cone or elliptic-cylinder carrier.

    The material side is `sign * Q(x) <= 0`. `side == "exterior"` means the material is on the
    far side of the carrier from its own inside, which is exactly a flipped sign.
    """
    from csg import recognise
    sign = -1.0 if carrier["side"] == "exterior" else 1.0
    kind = carrier["kind"]

    if kind == "plane":
        n = carrier["n"]
        p = carrier["p"]
        # Q(x) = n.(x - p); the material side of an outward normal is Q <= 0
        return sign, _quadric([[0.0] * 3 for _ in range(3)],
                              [0.5 * n[0], 0.5 * n[1], 0.5 * n[2]],
                              -(n[0] * p[0] + n[1] * p[1] + n[2] * p[2]))

    if kind == "sphere":
        p = carrier["p"]
        r = carrier["r"]
        identity = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        return sign, _quadric(identity, [-p[0], -p[1], -p[2]],
                              p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - r * r)

    if kind in ("cylinder", "cone"):
        d = carrier["d"]
        p = carrier["p"]
        r = carrier["r"]
        k = 0.0 if kind == "cylinder" else math.tan(carrier["a"])
        scale = 1.0 + k * k
        dd = _outer(d, d)
        a = [[(1.0 if i == j else 0.0) - scale * dd[i][j] for j in range(3)] for i in range(3)]
        ap = [sum(a[i][j] * p[j] for j in range(3)) for i in range(3)]
        pd = sum(p[i] * d[i] for i in range(3))
        b = [-ap[i] - r * k * d[i] for i in range(3)]
        c = sum(p[i] * ap[i] for i in range(3)) + 2.0 * r * k * pd - r * r
        return sign, _quadric(a, b, c)

    if kind == "ellipticCylinder":
        p = carrier["p"]
        xhat = carrier["x"]
        yhat = carrier["y"]
        semi_a = carrier["a"]
        semi_b = carrier["b"]
        xx = _outer(xhat, xhat)
        yy = _outer(yhat, yhat)
        a = [[xx[i][j] / (semi_a * semi_a) + yy[i][j] / (semi_b * semi_b) for j in range(3)]
             for i in range(3)]
        ap = [sum(a[i][j] * p[j] for j in range(3)) for i in range(3)]
        return sign, _quadric(a, [-ap[0], -ap[1], -ap[2]],
                              sum(p[i] * ap[i] for i in range(3)) - 1.0)

    raise recognise.Declined(f"a {kind} carrier has no quadric form")


def torus_from_carrier(carrier):
    """`(sign, centre, axis, major, minor)` for a torus carrier."""
    sign = -1.0 if carrier["side"] == "exterior" else 1.0
    return sign, list(carrier["p"]), list(carrier["d"]), carrier["r"], carrier["rt"]


def blocks_from_carriers(carriers):
    """One halfspace block per carrier, in the carriers' own order."""
    blocks = []
    for carrier in carriers:
        if carrier["kind"] == "torus":
            sign, centre, axis, major, minor = torus_from_carrier(carrier)
            blocks.append({"kind": "torus", "sign": sign,
                           "c": centre + axis + [major, minor, 0.0, 0.0, 0.0]})
        else:
            sign, block = quadric_from_carrier(carrier)
            blocks.append({"kind": "quadric", "sign": sign, "c": block + [0.0]})
    return blocks


def eval_block(block, point):
    """`sign * f(point)`; the halfspace contains the point when this is `<= 0`.

    The same arithmetic as `O2FlatCSG::EvalHalfspace`, and it exists so the self-test can
    classify a point without a built shape. The C++ side is the authority; if the two ever
    disagree the test that finds it is `emit.py --self-test`.
    """
    c = block["c"]
    x, y, z = point
    if block["kind"] == "torus":
        offset = (x - c[0], y - c[1], z - c[2])
        along = offset[0] * c[3] + offset[1] * c[4] + offset[2] * c[5]
        radial = tuple(offset[i] - along * c[3 + i] for i in range(3))
        rho = math.sqrt(sum(v * v for v in radial))
        return block["sign"] * (math.hypot(rho - c[6], along) - c[7])
    quadratic = (c[0] * x * x + c[3] * y * y + c[5] * z * z +
                 2.0 * (c[1] * x * y + c[2] * x * z + c[4] * y * z))
    return block["sign"] * (quadratic + 2.0 * (c[6] * x + c[7] * y + c[8] * z) + c[9])


def flat_contains(blocks, point):
    """True when every block contains the point: one cell's membership test."""
    return all(eval_block(block, point) <= 0.0 for block in blocks)


def write_sidecar(path, blocks, cells):
    """Write the version-1 flat-CSG sidecar. Must stay byte-compatible with `WriteFlatCSG`."""
    with open(path, "wb") as handle:
        handle.write(SIDECAR_MAGIC)
        handle.write(struct.pack("<III", SIDECAR_VERSION, len(blocks), len(cells)))
        for block in blocks:
            handle.write(struct.pack("<i", 1 if block["kind"] == "torus" else 0))
            handle.write(struct.pack("<d", block["sign"]))
            coefficients = list(block["c"]) + [0.0] * (11 - len(block["c"]))
            handle.write(struct.pack("<11d", *coefficients))
        for cell in cells:
            handle.write(struct.pack("<ii", cell["first"], cell["count"]))
            handle.write(struct.pack("<d", cell["volume"]))
            handle.write(struct.pack("<3d", *cell["lo"]))
            handle.write(struct.pack("<3d", *cell["hi"]))
```

Note for the implementer: `_halfspace_carriers` does not currently emit an `ellipticCylinder`
kind — R2 taught `_match_eltu` to read an extruded ellipse, but the carrier reader canonicalises
to the five kinds. Check `csg/tier0.py`'s `canonicalise` return before writing that branch: if no
elliptic kind reaches here, delete the branch and let the `Declined` cover it, and say so in the
task report. Do not invent a carrier kind that nothing produces.

- [ ] **Step 4: Run the self-test to verify it passes**

Run the Step 2 command. Expected: PASS, 301 + the new checks, with all six digest tables still
green.

- [ ] **Step 5: Verify the two writers agree byte for byte**

```bash
cd $HOME/alisw/O2 && python3 - <<'PY'
# the Python writer and WriteFlatCSG must produce identical bytes for the same solid
import subprocess, sys
sys.path.insert(0, "scripts/geometry")
from csg import flat
blocks = [{"kind": "quadric", "sign": 1.0, "c": [1., 0., 0., 1., 0., 0., 0., 0., 0., -25., 0.]}]
cells = [{"first": 0, "count": 1, "volume": 1.5, "lo": [-5., -5., -5.], "hi": [5., 5., 5.]}]
flat.write_sidecar("/tmp/claude-1000/flat_py.bin", blocks, cells)
print("python sidecar written")
PY
```

Then write the same content from a two-line C++ snippet using `WriteFlatCSG` and `cmp` the files.
If they differ, the Python writer is wrong — the C++ format is the one the loader reads.

- [ ] **Step 6: Commit**

```bash
cd $HOME/alisw/O2 && git add scripts/geometry/csg/flat.py scripts/geometry/csg/emit.py
git commit -m "$(cat <<'EOF'
Emit a cell's carriers as signed implicit halfspaces

This adds csg/flat.py, which maps the carriers the cell emitter already reads
to the halfspaces they always were, and writes the flat-CSG sidecar.

- A plane, sphere, cylinder, cone or elliptic cylinder becomes one ten-double
  quadric block; a torus keeps its canonical parameters.
- The sign composes the carrier's orientation with its material side, which is
  the one thing here that can be silently wrong.
- The self-test measures the flat halfspaces against _cell_leaf's padded
  primitives on sampled points, over all five carrier kinds.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Route the over-budget parts through the flat path

The rung's purpose: the 21 parts stop declining. Spec §8 and §10.

**Files:**
- Modify: `scripts/geometry/csg/recognise.py` (a flat branch and its budget)
- Modify: `scripts/geometry/csg/primitives.py` (a `flatCells` description and its OCCT builder)
- Modify: `scripts/geometry/O2_CADtoTGeo.py` (macro emission + sidecar write)
- Modify: `scripts/geometry/csg/emit.py` (routing self-tests)
- Modify: `scripts/geometry/csg/decline_catalogue.py` (the new decline reasons)

**Interfaces:**
- Consumes: everything from Tasks 1–8.
- Produces:
  - `primitives.flat_cells(cells, recogniser, notes=None)` → `{"op": "flatCells", "cells": [...], "recogniser", "notes"}`, where each cell is `{"blocks": [...], "volume": float, "lo": [3], "hi": [3]}`
  - `recognise._match_flat_cells(solid, records, tol, diag)` and `recognise.recognise_flat_cells(solid)`
  - `recognise._PART_MAX_FLAT_CELLS = 256`, `recognise._PART_MAX_FLAT_HALFSPACES = 1024`

- [ ] **Step 1: Write the failing self-test cases**

Add to `csg/emit.py`'s `--self-test`:

```python
    # --- the flat path takes what the tree budget refuses, and nothing else ------------------
    wide = make_boolean_fixtures.many_celled_solid()   # more cells than _PART_MAX_LEAVES allows
    tree_declined, tree_why = recognise.recognise_union_of_cells(wide, max_leaves=8)
    check("a part over the tree's leaf budget still declines on the tree path",
          tree_declined is None, tree_why)
    flat_record, flat_why = recognise.recognise_flat_cells(wide)
    check("the same part is accepted on the flat path",
          flat_record is not None and flat_record["recogniser"] == "flat-cells", flat_why)
    check("the flat record carries a bounding box per cell",
          flat_record is not None and
          all(len(c["lo"]) == 3 and len(c["hi"]) == 3 for c in flat_record["cells"]),
          "a cell is missing its box")

    over_cells, over_why = recognise.recognise_flat_cells(wide, max_cells=2)
    check("a part over the FLAT cell budget declines naming the bound",
          over_cells is None and "flat" in over_why and "2" in over_why, over_why)
    over_halfspaces, over_hs_why = recognise.recognise_flat_cells(wide, max_halfspaces=3)
    check("a part over the FLAT halfspace budget declines naming the bound",
          over_halfspaces is None and "3" in over_hs_why, over_hs_why)

    # the false-accept guard must run here too: OCCT Cut can report IsDone with zero solids both
    # ways, which is how ST0923290_01#b19 got through the symmetric difference in rung 4
    check("the flat path runs the containment corroboration",
          "containsScored" in (flat_record or {}).get("notes", {}),
          "no containment corroboration on the flat record")
```

- [ ] **Step 2: Run the self-test to verify it fails**

Run Task 8's Step 2 command. Expected: FAIL — `module 'csg.recognise' has no attribute 'recognise_flat_cells'`.

- [ ] **Step 3: Implement the description and its OCCT builder**

In `primitives.py`, beside `union_of_cells`:

```python
def flat_cells(cells, recogniser, notes=None):
    """A union of halfspace cells, for `O2FlatCSG`: the DNF with nothing bounded.

    `{op: "flatCells", cells: [...], recogniser, notes}`, and deliberately no `leaves` and no
    `op` on a cell, so nothing that reads a one- or two-level *primitive* description can read
    half of this one and get a wrong solid rather than a `KeyError`. A cell here is
    `{blocks, volume, lo, hi}` -- the halfspace blocks `csg/flat.py` produced, the cell's own
    volume, and the bounding box the decomposition measured, because an intersection of
    halfspaces does not bound itself.
    """
```

with validation mirroring `union_of_cells`'s strictness (every cell a dict with exactly those four
keys, at least one block per cell, `lo[i] <= hi[i]`), and **one cell is legal here** — unlike
`union_of_cells`, whose two-cell minimum stays as it is.

Its OCCT builder realises each cell by intersecting the *bounded* forms — reuse `_cell_leaf` via
the carriers the record kept — so the acceptance test still has a `TopoDS_Shape` to difference
against. Record in the docstring that the OCCT realisation is the padded one and the shipped
solid is not, and that this is sound because the padding covers the part's bounding box, which is
where the symmetric difference is measured.

- [ ] **Step 4: Implement the recogniser branch**

In `recognise.py`, add beside `_match_union_of_cells`:

```python
# The flat path's budgets. Not the tree's: `_PART_MAX_LEAVES` bounds how wide a
# `TGeoCompositeShape` may be, and `O2FlatCSG` is not one. These bound the sidecar and the box
# build instead, and are set an order of magnitude above the measured demand -- the worst part in
# `Handoff_FlatCSG.md`'s R5 table is `IBCYSSFlangeA` at about 59 cells.
_PART_MAX_FLAT_CELLS = 256
_PART_MAX_FLAT_HALFSPACES = 1024
```

`_match_flat_cells` reuses `decompose.split_into_cells` unchanged, then for each piece calls
`_halfspace_carriers` and `flat.blocks_from_carriers`, records the piece's OCCT volume and bounding
box, and runs the **same** acceptance as the union path: `_measured_gap`, then
`accept.contains_disagreements`. Do not skip the corroboration; §10 of the spec says why.

Routing: the flat path is tried **only after** the tree path declines, so no part that converts
today changes representation. Write that as a comment at the call site, not just in a report.

- [ ] **Step 5: Emit the shape into `geom.C`**

In `O2_CADtoTGeo.py`, add a `flatCells` branch to the shape emission that writes
`flatcsg_<part>.bin` next to the other sidecars and emits

```cpp
  auto* shape_<part> = new o2::base::O2FlatCSG("<part>");
  if (!o2::base::LoadFlatCSG(sidecarPath("flatcsg_<part>.bin"), *shape_<part>)) {
    ::Fatal("geom", "flat-CSG sidecar for <part> failed to load");
  }
  shape_<part>->CloseShape();
```

The `::Fatal` is deliberate and is the lesson of `NEXT.md` item 2: a module that fails to load must
be loud, never silently absent.

- [ ] **Step 6: Run the self-tests and the gates**

```bash
cd $HOME/alisw/O2 && python3 scripts/geometry/csg/emit.py --self-test
python3 scripts/geometry/O2_CADtoTGeo.py --self-test
python3 scripts/geometry/checkKnownSource.py --self-test
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/claude-1000/gate_r5 --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/claude-1000/gate_r5b \
    --model scripts/geometry/STEP_examples/Bagger.step
```

Expected: emit self-test 301 + the new checks with six digest tables green; converter 54;
checker 17; fixtures gate exit 0 with shape 10/10; Bagger exit 0, 13/13 + 7/7, CSG 7
bit-identical.

- [ ] **Step 7: Re-convert the five corpora and score them**

Strictly serially, one module at a time, in the sim env for the geometry and the converter env for
the conversion (see Global Constraints). For each of PIPE, ITS, TPC, ABSO, TRD:

```bash
# sim env
export O2_ROOT=$B/stage && $B/stage/bin/o2-sim-serial -n 0 -g boxgen -m <MOD>
# converter env
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root <MOD>.step --report <MOD>_writer_report.json
python3 scripts/geometry/O2_CADtoTGeo.py <MOD>.step -o geom.C --exact-surfaces auto --csg auto
python3 scripts/geometry/checkKnownSource.py --original o2sim_geometry.root \
    --writer-report <MOD>_writer_report.json --converted .
```

Expected: every floor of Global Constraints holds or improves, known-source rises above
1775/1775 by exactly the parts that newly convert, and every newly converted part is at
`dV_sym = 0`. **A floor that drops is a stop, not a note.** Diff the per-part reports against the
pre-task run and confirm that no previously converting part changed tier or evidence.

- [ ] **Step 8: Commit**

```bash
cd $HOME/alisw/O2 && git add scripts/geometry/csg/recognise.py scripts/geometry/csg/primitives.py \
  scripts/geometry/csg/emit.py scripts/geometry/csg/decline_catalogue.py scripts/geometry/O2_CADtoTGeo.py
git commit -m "$(cat <<'EOF'
Ship a part over the tree budget as a flat halfspace solid

This routes the parts the boolean-tree budget refuses through O2FlatCSG, which
is what rung R5 exists to do.

- The flat path is tried only after the tree path declines, so no part that
  converts today changes representation.
- Its budgets bound the sidecar and the box build, not the width of a
  TGeoCompositeShape, and are named in the decline messages.
- The containment corroboration runs here too, because OCCT's Cut can report
  IsDone with zero solids in both directions.
- A sidecar that fails to load is fatal rather than a silently absent module.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Measure the crossover and the split knobs, and write the record

Spec §9. The rung is not done when the parts convert; it is done when the emission policy rests on
a measurement.

**Files:**
- Modify: `Detectors/Base/test/runXRayBenchmark.cxx` (accept an `O2FlatCSG` sidecar as a subject)
- Create: `scripts/geometry/Stream_AK_FlatCSG.md`
- Modify: `scripts/geometry/NEXT.md`, `scripts/geometry/Handoff_FlatCSG.md`, `scripts/geometry/INDEX.md`
- Modify: `scripts/geometry/website_data/decline_reasons.json` (regenerate)

- [ ] **Step 1: Measure the crossover**

For every part in `Handoff_FlatCSG.md`'s R5 table that the tree budget *can* also express (raise
`_PART_MAX_LEAVES` temporarily in a scratch run — do **not** commit the raise), emit it both ways
and X-ray-bench the two. Report cost per `Contains`, per `DistFromOutside`, per `Safety`, and the
transport time, against both the cell count and the halfspace count.

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)" \
  && export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH \
  && $B/stage/bin/o2-detectors-base-xray --help
```

Record the crossover in cells and in halfspaces, and say which is the better predictor.

- [ ] **Step 2: Measure the split knobs**

Sweep `SetSplitDepth` over 3..10 and `SetMinBoxFraction` over 0.002..0.05 on the three largest
parts, recording box count, mean and maximum leaf active-list length, `GetBVHMemory()`, build
time, and query cost. Pick the defaults from the data and change them in the header in this task,
with a one-line comment naming the measurement.

- [ ] **Step 3: Measure the safety quality**

On the fixtures, compare `Safety` against the true distance (available from the fixtures' analytic
forms) and report the ratio distribution, plus the fraction of inside queries that land in an
undecided box and therefore get `0.`. If that fraction is large enough to hurt transport, say so
and record it as the first item of the successor stream — do not fix it in this task.

- [ ] **Step 4: Write `Stream_AK_FlatCSG.md`**

The successor stream document `Handoff_FlatCSG.md` §3 asks for. It records what was built, every
number from Steps 1–3, the emission policy that follows from the crossover, the parts that now
convert with their cell and halfspace counts, and an honest concerns section in the style of the
R4 report — including the `G4MultiUnion` correspondence from spec §11.5 and the fact that the 24
boundary-gap declines are still declined and now affordable to attack.

- [ ] **Step 5: Rewrite `NEXT.md` and update the handoff and index**

`NEXT.md` is the standing entry point and whoever finishes a session rewrites it: new corpus rows,
the new ctest floor including `FlatCSG`, R5 marked done with its measured crossover, and R6 named
as next. Add `Design_FlatCSGSolid.md`, `Plan_FlatCSGSolid.md` and `Stream_AK_FlatCSG.md` to
`INDEX.md`. Mark R5 done in `Handoff_FlatCSG.md` §1.1 in the same style as R1–R4.

- [ ] **Step 6: Regenerate the decline catalogue and the website data**

```bash
cd $HOME/alisw/O2 && python3 scripts/geometry/csg/decline_catalogue.py --regenerate
```

- [ ] **Step 7: Run every floor one last time**

```bash
ctest -R 'FlatCSG|BVHSurfaceSolid|BVHAssembly' --output-on-failure
python3 scripts/geometry/csg/emit.py --self-test
python3 scripts/geometry/checkKnownSource.py --self-test
python3 scripts/geometry/O2_CADtoTGeo.py --self-test
python3 scripts/geometry/O2_TGeoToCAD.py --self-test
python3 scripts/geometry/runOracleGate.py --self-test
```

Expected: every count at or above its Global Constraints floor.

- [ ] **Step 8: Commit**

```bash
cd $HOME/alisw/O2 && git add scripts/geometry/Stream_AK_FlatCSG.md scripts/geometry/NEXT.md \
  scripts/geometry/Handoff_FlatCSG.md scripts/geometry/INDEX.md \
  scripts/geometry/website_data/decline_reasons.json \
  Detectors/Base/test/runXRayBenchmark.cxx Detectors/Base/include/DetectorsBase/O2FlatCSG.h
git commit -m "$(cat <<'EOF'
Record the flat solid's measurements and the emission policy

This adds Stream_AK_FlatCSG.md, the record of rung R5, and sets the subdivision
defaults from the measurement rather than from a guess.

- The crossover against plain-composite emission is measured in cells and in
  halfspaces, and the better predictor is named.
- The split depth and minimum box size are swept on the three largest parts.
- Safety is compared against the true distance on the fixtures, and the share of
  inside queries that fall in an undecided box is reported.
- NEXT.md is rewritten and R5 is marked done in the flat-CSG hand-off.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

## Self-review

**Spec coverage.** §1 → Tasks 1–5, 9. §3.1 the quadric block → Task 1; the torus block → Task 3.
§3.2 → Task 8 (the emitter that drops the padding) and Task 9's OCCT-realisation note. §3.3 the
sign convention and its `_cell_leaf` oracle → Task 8 Step 1. §4.1–4.3 → Task 4. §4.4 the
clip-inside-the-box invariant → Task 5 Step 3, with the long-arm case as its falsifier. §5.1 →
Task 1 and Task 5. §5.2 → Task 2. §5.3 `Capacity` → Task 6. §5.4 `Safety` → Task 6. §5.5
`ComputeBBox`/`ComputeNormal` → Task 6. §6 the twins → Tasks 1, 2, 5, 6. §7 persistence → Task 7.
§8 the converter → Tasks 8 and 9. §9 the measurements → Task 10. §10 the floors → Global
Constraints, re-checked in Tasks 1, 7, 9 and 10. §11's risks map to Task 8 Step 1 (risk 1), Task 5
Step 3 (risk 2), Task 10 Step 2 (risk 3), Task 6 (risk 4), Task 10 Step 4 (risks 5 and 6).

**Type consistency.** `FlatCSGHalfspace`/`FlatCSGCell`/`FlatCSGBox` field names are used
identically in Tasks 1, 4, 5, 6 and 7. `EvalHalfspace` returns `sign * f` and "inside" is `<= 0`
in every task, in both C++ and `csg/flat.py`'s `eval_block`. `CellIntervals`'s `active`/`nActive`
contract (`nullptr` with `-1` means "all of them") is used by the twins in Task 2 and by the
per-box calls in Task 5. The sidecar field order in Task 7's table matches `write_sidecar`'s
`struct.pack` calls in Task 8, and Task 8 Step 5 is the check that they do.

**Two things the implementer must not paper over**, both flagged in place: Task 3's note that
`solveQuarticReal`'s enclosing namespace has to be read rather than assumed, and Task 8's note
that `ellipticCylinder` may not be a carrier kind anything actually produces — in which case the
branch is deleted, not invented around.
