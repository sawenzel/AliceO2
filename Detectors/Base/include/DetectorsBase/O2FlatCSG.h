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

  // The shape owns a raw `bvh::v2::Bvh` behind `fBVH`, so a compiler-written copy would hand two
  // shapes the same BVH and then free it twice; same treatment as O2BVHAssembly.
  O2FlatCSG(const O2FlatCSG&) = delete;
  O2FlatCSG& operator=(const O2FlatCSG&) = delete;

  // ---- building -------------------------------------------------------------------------
  /// Append a quadric halfspace; returns its index. `sign` is +1 or -1, inside is `sign*Q <= 0`.
  int AddQuadric(double sign, const double coeff[10]);
  /// Append a torus halfspace: inside is `sign * (sqrt((rho - major)^2 + z^2) - minor) <= 0`,
  /// with `rho` the distance from the axis through \a centre along \a axis and `z` the coordinate
  /// along it. \a axis must be a unit vector. Returns the halfspace index.
  int AddTorus(double sign, const double* centre, const double* axis, double major, double minor);
  /// Append a cell over `[first, first + count)` of the halfspace array; returns its index.
  int AddCell(int first, int count, double volume);

  int GetNhalfspaces() const { return static_cast<int>(fHalfspaces.size()); }
  int GetNcells() const { return static_cast<int>(fCells.size()); }
  const FlatCSGHalfspace& GetHalfspace(int index) const { return fHalfspaces[index]; }
  const FlatCSGCell& GetCell(int index) const { return fCells[index]; }

  /// The AABB of cell \a cell. The halfspaces alone do not bound a cell -- an intersection of
  /// halfspaces can be unbounded -- so the converter supplies the box the decomposition measured.
  void SetCellBBox(int cell, const double* lo, const double* hi);
  /// The AABB `SetCellBBox` recorded for cell \a cell, for the sidecar writer. Reads back zeros
  /// for a cell whose box was never set.
  void GetCellBBox(int cell, double* lo, double* hi) const;

  /// Build the sub-cell boxes (and, from task 5, the BVH). Call once, after the last AddCell.
  void CloseShape();
  bool IsClosed() const { return fClosed; }

  /// Bytes held by the BVH nodes and the primitive-index permutation.
  size_t GetBVHMemory() const;

  int GetNboxes() const { return static_cast<int>(fBoxes.size()); }
  const FlatCSGBox& GetBox(int index) const { return fBoxes[index]; }
  /// For the tests: the box structure is the thing being proved sound, so it has to be readable.
  int GetActive(int index) const { return fActive[index]; }
  /// For the tests: the box structure is the thing being proved sound, so it has to be readable.
  bool CellContainsPublic(int cell, const double* point) const { return CellContains(cell, point); }

  /// Subdivision depth cap. See `fSplitDepth` for where the default comes from.
  void SetSplitDepth(int depth) { fSplitDepth = depth; }
  /// Stop splitting a box narrower than this fraction of the part's bounding-box diagonal.
  /// See `fMinBoxFraction` for where the default comes from.
  void SetMinBoxFraction(double fraction) { fMinBoxFraction = fraction; }

  /// `sign * f(point)`; the halfspace contains the point when this is `<= 0`.
  static double EvalHalfspace(const FlatCSGHalfspace& halfspace, const double* point);

  /// A rigorous enclosure `[rangeLo, rangeHi]` of `sign * f` over the box `[lo, hi]`.
  ///
  /// Conservative in the only safe direction: an over-wide enclosure loses pruning, never
  /// correctness. For a quadric it is the centred form -- with `m` the centre, `h` the
  /// half-extents and `g = A m + b`, `|Q(x) - Q(m)| <= 2 sum|g_i| h_i + sum |A_ij| h_i h_j` --
  /// which tightens quadratically as the boxes shrink. For a torus it is the 1-Lipschitz signed
  /// distance, so the enclosure is `f(m) +/- |h|`. Both are padded by a small multiple of the
  /// magnitude actually accumulated computing `f(m)`, not of `f(m)` itself, so the margin does
  /// not collapse to nothing on a box straddling the surface, where `f(m)` is heavily cancelled.
  ///
  /// REQUIRES `lo[i] <= hi[i]` on every axis (so every half-extent is `>= 0`) and every component
  /// of `lo`/`hi` finite. `CloseShape` enforces both on the boxes it feeds this function, because
  /// the quadric branch's cross-term bound `sum |A_ij| h_i h_j` is only an over-estimate of the
  /// true deviation when every `h_i, h_j >= 0` -- see the AM-GM identity at the cross-term
  /// accumulation in the .cxx -- and because a NaN half-extent or centre would produce a NaN
  /// range that neither drop test in `SplitBox` can act on, silently keeping a box that should
  /// have been rejected. A caller that builds boxes without going through `CloseShape` must keep
  /// both invariants itself; this function does not check them.
  static void HalfspaceRange(const FlatCSGHalfspace& halfspace, const double* lo, const double* hi,
                             double& rangeLo, double& rangeHi);

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
  //
  // All three accelerated queries fall back to their `_Loop` twin when `IsClosed()` is false.
  // `CloseShape` refuses, and builds nothing, when a cell's bounding box is missing, inverted or
  // non-finite; an accelerated query that then walked an empty box array would answer "no material
  // anywhere", which is the silent vanishing that refusal exists to prevent.
  Bool_t Contains(const Double_t* point) const override;

  Double_t DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;
  Double_t DistFromInside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                          Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;

  /// A rigorous lower bound on the distance to the boundary, from the box structure alone -- no
  /// point-to-quadric distance formula exists in general, so this class does not try to have one.
  /// Outside, the distance to the nearest box is a bound because every point of the solid is in
  /// some box. Inside, only a box with an empty active list (wholly inside its cell, a hard
  /// guarantee from Task 4) contributes: the distance to its faces is a bound. A point in an
  /// undecided box gets the sound answer `0.`, because the box's active list says nothing about
  /// where the boundary sits within it. Falls back to `Safety_Loop` when `!IsClosed()`, like the
  /// other accelerated queries.
  Double_t Safety(const Double_t* point, Bool_t in = kTRUE) const override;

  /// The union of the retained sub-cell boxes -- tighter than the union of the cell AABBs, which
  /// is the point of design section 5.5.
  void ComputeBBox() override;

  /// The sum of the cells' own volumes. The cells of a decomposition are disjoint by construction
  /// (`decompose`'s volume guard checks it), so there is no inclusion-exclusion to do.
  Double_t Capacity() const override;

  /// The outward normal (oriented with respect to `dir`, per the `TGeoShape` contract) at the
  /// halfspace closest to being satisfied with equality at `point`: `sign * 2(Ax + b)` for a
  /// quadric, or the gradient of the torus's signed distance, normalised and flipped along `dir`.
  void ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) const override;

  // ---- the reference twins ---------------------------------------------------------------
  Bool_t Contains_Loop(const Double_t* point) const;

  Double_t DistFromOutside_Loop(const Double_t* point, const Double_t* dir,
                                Double_t step = TGeoShape::Big()) const;
  Double_t DistFromInside_Loop(const Double_t* point, const Double_t* dir,
                               Double_t step = TGeoShape::Big()) const;
  /// `Safety`'s twin: walks all boxes with no BVH. Required to equal `Safety` (proving the BVH
  /// pruning loses no nearer box) and, independently, to be a SOUND bound -- the soundness leg is
  /// the one that matters, since a safety merely self-consistent with its twin could still be
  /// wrong.
  Double_t Safety_Loop(const Double_t* point, Bool_t in = kTRUE) const;

 protected:
  /// True when every halfspace of cell `index` contains `point`.
  bool CellContains(int index, const double* point) const;

  /// The accelerated bodies of `DistFromOutside`/`DistFromInside`, called once the BVH is there.
  /// Both walk the boxes and, per box, clip the ray to that box's slab BEFORE running the interval
  /// clipping over that box's active list; see scripts/geometry/Design_FlatCSGSolid.md 4.4.
  Double_t DistFromOutsideBVH(const Double_t* point, const Double_t* dir, Double_t step) const;
  Double_t DistFromInsideBVH(const Double_t* point, const Double_t* dir, Double_t step) const;

  /// For every sub-cell box the ray meets within `[0, step]`, that box's own occupancy over that
  /// box's own clipped window: \a pairs gets `[enter, exit]` per piece and \a cells the cell each
  /// piece came from, one entry per pair. Both are cleared first.
  ///
  /// The pieces are deliberately NOT merged here. Boxes of one cell tile a region, so a cell's
  /// occupancy arrives in several adjacent pieces, and how they are rejoined differs between the
  /// two callers: `DistFromOutside` needs the per-cell intervals its twin computes, while
  /// `DistFromInside` needs the union across cells.
  ///
  /// Returns false when a `CellIntervals` call overflowed its output buffer. That is a bug in the
  /// buffer sizing rather than live data, so it must not be accumulated into an answer.
  bool GatherRayPieces(const Double_t* point, const Double_t* dir, Double_t step,
                       std::vector<double>& pairs, std::vector<int>& cells) const;

  /// Recursively split `[lo, hi]` for `cell`, dropping halfspaces the range bound proves hold
  /// everywhere in the box and returning (dropping the box) when the bound proves it wholly
  /// outside one. Recursion stops -- and the box is kept -- when `active` empties, `depth`
  /// reaches zero, the box's longest side is no longer than `minSize`, or the box is STILL far
  /// from cubic once `cubifyBudget` also reaches zero (a near-cubic box never checks
  /// `cubifyBudget` at all -- see below).
  ///
  /// `depth` and `cubifyBudget` are two separate purses, because a level cap denominated purely
  /// in tree depth charges full price for a split that only equalises aspect ratio: a 20x2x2 cell
  /// spends four of six levels just becoming roughly cubic, leaving two to actually detach a leaf
  /// from six faces. So a split that leaves the box still far from cubic (`longest > 2 *
  /// shortest`, checked on the box BEFORE that split) draws from `cubifyBudget` instead of
  /// `depth` -- `depth` is untouched until the box is within a factor of two on every axis, at
  /// which point it behaves exactly as before. `cubifyBudget` is the hard ceiling design section
  /// 4.2 describes, on the aspect-ratio-equalising splits one root-to-leaf PATH may spend: it only
  /// ever decreases going down the recursion, so it bounds one branch, not the whole cell -- a
  /// near-cubic cell never draws on it at all (design section 4.2's invariance argument: halving
  /// the longest of three extents with `c <= 2a` keeps the new ratio `<= 2`, so a box that starts
  /// near-cubic can never become far-from-cubic), and is therefore unaffected by its value.
  void SplitBox(int cell, const double* lo, const double* hi, const std::vector<int>& active,
               int depth, double minSize, int cubifyBudget);

  std::vector<FlatCSGHalfspace> fHalfspaces; ///< the flat halfspace array
  std::vector<FlatCSGCell> fCells;           ///< the DNF's cells, indexing into it

  /// The sub-cell boxes produced by `CloseShape`. Not streamed -- like `fBVH`, this is rebuilt
  /// from `fHalfspaces`, `fCells`, `fCellLo` and `fCellHi` by `CloseShape`, and streaming it would
  /// only be a second thing that can disagree with the data it describes (design section 7).
  std::vector<FlatCSGBox> fBoxes; //!
  /// The boxes' active-halfspace lists, concatenated. Derived alongside `fBoxes`; not streamed
  /// for the same reason.
  std::vector<int> fActive; //!
  std::vector<double> fCellLo; ///< each cell's AABB low corner, 3 doubles per cell
  std::vector<double> fCellHi; ///< each cell's AABB high corner, 3 doubles per cell
  /// Whether `SetCellBBox` was ever called for a given cell; `CloseShape` refuses to build a
  /// solid missing one rather than silently drop that cell -- see `CloseShape`'s implementation.
  std::vector<bool> fCellBBoxSet;
  /// Set only by a successful `CloseShape`; not streamed, since it describes `fBoxes`/`fBVH`
  /// rather than the shape's own source data, and the caller is always expected to call
  /// `CloseShape` again after a read (design section 7, and the loader's own contract).
  bool fClosed = false; //!
  /// Subdivision depth cap, and the minimum box size as a fraction of the part's bounding-box
  /// diagonal. Both set from the sweep in scripts/geometry/Stream_AK_FlatCSG.md section 5
  /// (depth 3..10 x fraction 0.002..0.2 on six shipped parts): every kernel gets cheaper as
  /// the boxes get coarser, because a cell of this rung carries 5 to 16 halfspaces and
  /// evaluating the ones a coarse box leaves undecided is cheaper than the BVH descent that
  /// would have removed them. 0.05 is the coarsest setting that still leaves a real sub-cell
  /// structure (2.5 to 19 boxes per cell); the depth cap only binds on one of the six parts
  /// there, and 4 costs it 15 % of one kernel to keep a level of headroom for a part whose
  /// cells are large relative to the part, which is the case the cap exists for.
  int fSplitDepth = 4;
  double fMinBoxFraction = 0.05;

  /// The BVH over `fBoxes`, built by `CloseShape`. Never streamed and never persisted -- it is
  /// rebuilt from the boxes, which are themselves rebuilt from the halfspaces (design section 7),
  /// so a stored BVH would only be a second thing that can disagree with the data it describes.
  void* fBVH = nullptr; //! bvh::v2::Bvh over the sub-cell boxes

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
