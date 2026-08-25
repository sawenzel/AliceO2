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

#ifndef ALICEO2_BASE_O2FLATCSG_
#define ALICEO2_BASE_O2FLATCSG_

#include "TGeoBBox.h"
#include "TGeoShape.h"

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

  /// Real roots of `sign * f(origin + t*dir) = 0`, unsorted, at most four; returns the count.
  static int HalfspaceRoots(const FlatCSGHalfspace& halfspace, const double* origin,
                            const double* dir, double* roots);

  /// The occupancy of cell \a cell along the ray, restricted to `[tlo, thi]`.
  ///
  /// \a active lists which of the cell's halfspaces are still undecided; pass `nullptr` with
  /// \a nActive `< 0` to use all of them, which is what the twins do. Writes `[enter, exit]`
  /// pairs into \a out and returns the pair count -- or a negative value if \a maxOut was too
  /// small to hold every pair this cell produced along this ray. A mis-sized caller must not
  /// silently receive a truncated list, so overflow fails loudly instead of returning a short
  /// count that looks like a valid answer.
  ///
  /// No convexity is assumed anywhere: a complemented halfspace makes a cell non-convex and the
  /// occupancy several intervals, which is exactly what every part with a hole in it looks like.
  int CellIntervals(int cell, const int* active, int nActive, const double* origin,
                    const double* dir, double tlo, double thi, double* out, int maxOut) const;

  // ---- the TGeoShape contract ------------------------------------------------------------
  Bool_t Contains(const Double_t* point) const override;

  Double_t DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;
  Double_t DistFromInside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                          Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;

  /// task 6 replaces this; 0 is always sound
  Double_t Safety(const Double_t* point, Bool_t in = kTRUE) const override;

  // ---- the reference twins ---------------------------------------------------------------
  Bool_t Contains_Loop(const Double_t* point) const;

  Double_t DistFromOutside_Loop(const Double_t* point, const Double_t* dir,
                                Double_t step = TGeoShape::Big()) const;
  Double_t DistFromInside_Loop(const Double_t* point, const Double_t* dir,
                               Double_t step = TGeoShape::Big()) const;

 protected:
  /// True when every halfspace of cell `index` contains `point`.
  bool CellContains(int index, const double* point) const;

  std::vector<FlatCSGHalfspace> fHalfspaces; ///< the flat halfspace array
  std::vector<FlatCSGCell> fCells;           ///< the DNF's cells, indexing into it

  // `CellIntervals`, `DistFromOutside_Loop` and `DistFromInside_Loop` each grow a scratch buffer
  // on demand; those live as `thread_local` function-local statics in the .cxx, not as members --
  // TGeo shares one shape object across every navigator under `TGeoManager::SetMaxThreads`, and a
  // `std::vector` resized on one thread while another holds its `data()` is undefined behaviour,
  // not merely a stale read.

  ClassDefOverride(O2FlatCSG, 1) // flat-DNF halfspace shape class
};

} // namespace base
} // namespace o2

#endif
