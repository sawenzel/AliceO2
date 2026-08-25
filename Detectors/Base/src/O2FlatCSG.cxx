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

#include "DetectorsBase/O2FlatCSG.h"

#include "BoundedSurface.h"

// the same third-party BVH2 entry point O2Tessellated, O2BVHSurfaceSolid and O2BVHAssembly use
#include "bvh2_third_party.h"
#include "bvh2_extra_kernels.h"

#include "TGeoShape.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

ClassImp(o2::base::O2FlatCSG);

/// The most roots one cell can contribute to one ray: four per torus halfspace.
constexpr int kMaxRootsPerHalfspace = 4;

/// A hard ceiling on the aspect-ratio-equalising splits `SplitBox` may spend along one
/// root-to-leaf PATH before it gives up and keeps that box as-is, however far from cubic --
/// design section 4.2's two-purse split rule. It is a per-path bound, not a per-cell one:
/// `cubifyBudget` only ever decreases going down the recursion, so the worst case for one cell is
/// `2^kMaxCubifySplits` leaves along its widest branch, not `kMaxCubifySplits` splits in total.
///
/// The ceiling bites only once a cell's aspect ratio survives that many halvings of its longest
/// axis, and how far that reaches depends on the SHAPE of the cell, not just its worst ratio: a
/// rod `(r, 1, 1)` keeps the same axis longest every split, so N splits buy a full `2^N` reduction
/// and this value alone covers a `2^11 = 2048:1` rod; a plate `(r, r, 1)` alternates between its
/// two long axes, so only every other split reduces either one, covering `2^6 = 64:1` -- four
/// orders of magnitude short of the rod case for the same budget. 10 is chosen over the fix round
/// 1 value of 16 because 16 covers rods well past anything this class ships (`2^17 ~= 131000:1`)
/// at a plate cost the reviewer judged not worth carrying; 10 keeps a generous margin over the
/// 1000:1 case design section 4.2 cites (as a rod) while still leaving real cover for a plate.
/// Either way this is a worst-case cap, not a target: an ordinary cell resolves in a handful of
/// splits (see `SplitBox`'s header doc comment for the near-cubic invariance argument), and it is
/// this constant, not `fSplitDepth`, that Task 10 should revisit if a shipped part's cells turn
/// out to need more.
constexpr int kMaxCubifySplits = 10;

namespace
{
/// A safe upper bound on the number of `[enter, exit]` pairs one cell can produce along a ray:
/// each of its halfspaces contributes at most `kMaxRootsPerHalfspace` roots to the break list,
/// and a run of merged "inside" sub-intervals cannot outnumber the gaps between breaks.
int maxPairsForCell(int halfspaceCount)
{
  return 2 + kMaxRootsPerHalfspace * halfspaceCount;
}
} // namespace

namespace o2
{
namespace base
{

namespace
{
// float BVH types, following the O2BVHAssembly / O2Tessellated pattern. Float is enough here
// because the BVH only ever *nominates* candidates: every box a leaf hands over is then clipped
// and evaluated against its own double `FlatCSGBox` bounds, so a node box has to be a superset of
// the geometry and nothing more. roundOutward is what makes it one.
using BVHScalar = float;
using BVHBBox = bvh::v2::BBox<BVHScalar, 3>;
using BVHVec3 = bvh::v2::Vec<BVHScalar, 3>;
using BVHNode = bvh::v2::Node<BVHScalar, 3>;
using BVH = bvh::v2::Bvh<BVHNode>;

/// Round a double outward into float, away from the interval the box encloses.
inline float roundOutward(double value, bool up)
{
  return std::nextafterf(static_cast<float>(value), up ? std::numeric_limits<float>::infinity()
                                                       : -std::numeric_limits<float>::infinity());
}

/// Narrow `[tlo, thi]` to the ray's parameter window inside the axis-aligned box
/// `[boxMin, boxMax]`; returns false when nothing survives.
///
/// This is the clip design section 4.4 is about. An active list describes its cell only inside its
/// own box, so this window -- computed per box, and never pooled across boxes -- is what bounds
/// the `CellIntervals` call that uses that list.
///
/// It divides by `dir[i]` rather than multiplying by a precomputed `1/dir[i]`, and that is not a
/// stylistic choice. A cell's bounding box routinely has a face lying exactly on one of that
/// cell's own axis-aligned plane halfspaces -- the arm of an L-bracket is six such faces -- and
/// there the box face and the halfspace surface are the same plane, crossed at the same parameter.
/// `HalfspaceRoots` reaches it as `-0.5*gamma/beta`, whose exact powers of two cancel, so it is
/// `(v - o_k) / d_k`; `(v - o_k) * (1/d_k)` rounds twice and can land an ulp away, which is enough
/// to make an accelerated exit distance differ from the twin's in the last bit. One division is
/// also simply the more accurate of the two.
bool slabWindow(const double* boxMin, const double* boxMax, const double* origin, const double* dir,
                double& tlo, double& thi)
{
  for (int index = 0; index < 3; ++index) {
    if (std::abs(dir[index]) < 1.e-300) {
      // parallel to this pair of faces: the ray is either inside the slab for every t or outside
      // it for every t
      if (origin[index] < boxMin[index] || origin[index] > boxMax[index]) {
        return false;
      }
      continue;
    }
    double low = (boxMin[index] - origin[index]) / dir[index];
    double high = (boxMax[index] - origin[index]) / dir[index];
    if (low > high) {
      std::swap(low, high);
    }
    tlo = std::max(tlo, low);
    thi = std::min(thi, high);
    if (tlo > thi) {
      return false;
    }
  }
  return true;
}

/// The same clip against a BVH node's (float, outward-rounded) box.
inline bool nodeWindow(const BVHBBox& box, const double* origin, const double* dir, double& tlo,
                       double& thi)
{
  const double lo[3] = {box.min[0], box.min[1], box.min[2]};
  const double hi[3] = {box.max[0], box.max[1], box.max[2]};
  return slabWindow(lo, hi, origin, dir, tlo, thi);
}

/// Whether \a point is in the box's own double bounds, closed on every face.
inline bool boxHoldsPoint(const FlatCSGBox& box, const double* point)
{
  return point[0] >= box.min[0] && point[0] <= box.max[0] && point[1] >= box.min[1] &&
         point[1] <= box.max[1] && point[2] >= box.min[2] && point[2] <= box.max[2];
}

/// Gradient (unnormalised) of `sign * f` at \a point, for `ComputeNormal`: `2(Ax + b)` for a
/// quadric, scaled by `sign`; the gradient of the torus's signed distance, scaled by `sign`,
/// otherwise. The torus derivation: with `offset = point - centre`, `along = offset . axis`,
/// `radial = offset - along * axis`, `rho = |radial|`, `u = rho - major`, `s = hypot(u, along)`,
/// the unsigned distance is `h = s - minor`, and the chain rule through `d(rho)/dx = radial / rho`
/// (radial is already perpendicular to axis) and `d(along)/dx = axis` gives
/// `grad h = (u / s) * (radial / rho) + (along / s) * axis`.
void halfspaceGradient(const FlatCSGHalfspace& halfspace, const double* point, double grad[3])
{
  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    const double* c = halfspace.c;
    const double axis[3] = {c[3], c[4], c[5]};
    const double major = c[6];
    const double offset[3] = {point[0] - c[0], point[1] - c[1], point[2] - c[2]};
    const double along = offset[0] * axis[0] + offset[1] * axis[1] + offset[2] * axis[2];
    double radial[3];
    for (int index = 0; index < 3; ++index) {
      radial[index] = offset[index] - along * axis[index];
    }
    const double rho = std::sqrt(radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]);
    const double u = rho - major;
    const double s = std::hypot(u, along);
    if (s < 1.e-300 || rho < 1.e-300) {
      // degenerate: exactly on the revolution axis, or at the surface's own kissing point where
      // the signed distance is not differentiable -- leave the gradient at zero and let the
      // caller fall back
      grad[0] = grad[1] = grad[2] = 0.;
      return;
    }
    const double du = u / s;
    const double dv = along / s;
    for (int index = 0; index < 3; ++index) {
      grad[index] = halfspace.sign * (du * (radial[index] / rho) + dv * axis[index]);
    }
    return;
  }
  const double* c = halfspace.c;
  const double a[3][3] = {{c[0], c[1], c[2]}, {c[1], c[3], c[4]}, {c[2], c[4], c[5]}};
  const double b[3] = {c[6], c[7], c[8]};
  for (int row = 0; row < 3; ++row) {
    double value = b[row];
    for (int column = 0; column < 3; ++column) {
      value += a[row][column] * point[column];
    }
    grad[row] = halfspace.sign * 2. * value;
  }
}

/// Hand every leaf primitive whose node box the ray meets within `[0, tmax]` to \a visit.
///
/// The node boxes are outward-rounded supersets of the sub-cell boxes, so this can only ever
/// nominate too many; the caller decides with the box's own double bounds.
template <typename Visit>
void traverseRay(const BVH& bvh, const double* origin, const double* dir, double tmax, Visit&& visit)
{
  // thread_local rather than a member or a fresh vector per call: TGeo shares one shape object
  // across every navigator under TGeoManager::SetMaxThreads, and this is not re-entered
  thread_local std::vector<size_t> stack;
  stack.clear();
  stack.push_back(0); // the bvh2 root node
  while (!stack.empty()) {
    const size_t current = stack.back();
    stack.pop_back();
    const auto& node = bvh.nodes[current];
    double tlo = 0.;
    double thi = tmax;
    if (!nodeWindow(node.get_bbox(), origin, dir, tlo, thi)) {
      continue;
    }
    if (node.is_leaf()) {
      const auto beginPrimitive = node.index.first_id();
      const auto endPrimitive = beginPrimitive + node.index.prim_count();
      for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
        visit(static_cast<int>(bvh.prim_ids[primitive]));
      }
    } else {
      const auto firstChild = node.index.first_id();
      for (size_t child : {firstChild, firstChild + 1}) {
        if (child < bvh.nodes.size()) {
          stack.push_back(child);
        }
      }
    }
  }
}
} // namespace

O2FlatCSG::O2FlatCSG() : TGeoBBox(0., 0., 0.) {}

O2FlatCSG::O2FlatCSG(const char* name) : TGeoBBox(name, 0., 0., 0.) {}

O2FlatCSG::~O2FlatCSG()
{
  delete static_cast<BVH*>(fBVH);
  fBVH = nullptr;
}

size_t O2FlatCSG::GetBVHMemory() const
{
  const auto* bvh = static_cast<const BVH*>(fBVH);
  if (bvh == nullptr) {
    return 0;
  }
  return bvh->nodes.size() * sizeof(BVHNode) + bvh->prim_ids.size() * sizeof(size_t);
}

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

int O2FlatCSG::AddTorus(double sign, const double* centre, const double* axis, double major,
                        double minor)
{
  FlatCSGHalfspace halfspace;
  halfspace.kind = FlatCSGHalfspace::kTorus;
  halfspace.sign = sign < 0. ? -1. : 1.;
  // both evaluators assume |axis| == 1; normalising here (rather than asserting it and trusting
  // the caller) means a converter that hands in an OCCT axis with ordinary floating-point noise
  // is corrected once, silently, instead of being wrong in every EvalHalfspace/HalfspaceRoots
  // call for the life of the shape -- a genuinely degenerate (zero) axis is still a caller bug,
  // so that case still asserts rather than dividing by zero
  const double axisNorm = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  assert(axisNorm > 0. && "O2FlatCSG::AddTorus: axis must not be the zero vector");
  for (int index = 0; index < 3; ++index) {
    halfspace.c[index] = centre[index];
    halfspace.c[3 + index] = axis[index] / axisNorm;
  }
  halfspace.c[6] = major;
  halfspace.c[7] = minor;
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

void O2FlatCSG::SetCellBBox(int cell, const double* lo, const double* hi)
{
  if (cell < 0 || cell >= GetNcells()) {
    // "call after the last AddCell" is only a doc comment; without this check a cell index that
    // predates its AddCell, or is simply wrong, would write past the end of fCellLo/fCellHi --
    // silent heap corruption in release, not a bounds error anyone would see at the call site
    Error("SetCellBBox", "Shape %s: cell %d is out of range (%d cell(s) so far); ignoring",
          GetName(), cell, GetNcells());
    return;
  }
  if (static_cast<int>(fCellBBoxSet.size()) < GetNcells()) {
    fCellLo.resize(3 * GetNcells(), 0.);
    fCellHi.resize(3 * GetNcells(), 0.);
    fCellBBoxSet.resize(GetNcells(), false);
  }
  for (int index = 0; index < 3; ++index) {
    fCellLo[3 * cell + index] = lo[index];
    fCellHi[3 * cell + index] = hi[index];
  }
  fCellBBoxSet[cell] = true;
}

double O2FlatCSG::EvalHalfspace(const FlatCSGHalfspace& halfspace, const double* point)
{
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
  const double* c = halfspace.c;
  const double x = point[0];
  const double y = point[1];
  const double z = point[2];
  const double quadratic = c[0] * x * x + c[3] * y * y + c[5] * z * z +
                           2. * (c[1] * x * y + c[2] * x * z + c[4] * y * z);
  const double linear = 2. * (c[6] * x + c[7] * y + c[8] * z);
  return halfspace.sign * (quadratic + linear + c[9]);
}

void O2FlatCSG::HalfspaceRange(const FlatCSGHalfspace& halfspace, const double* lo,
                               const double* hi, double& rangeLo, double& rangeHi)
{
  // Preconditions this function does not check at runtime (see the header doc comment): every
  // half-extent nonnegative (lo[i] <= hi[i]) and every bound finite. CloseShape enforces both on
  // the boxes it feeds SplitBox/HalfspaceRange; this catches only a caller going around it.
  assert(std::isfinite(lo[0]) && std::isfinite(lo[1]) && std::isfinite(lo[2]) &&
        std::isfinite(hi[0]) && std::isfinite(hi[1]) && std::isfinite(hi[2]) &&
        lo[0] <= hi[0] && lo[1] <= hi[1] && lo[2] <= hi[2] &&
        "O2FlatCSG::HalfspaceRange: lo/hi must be finite and lo[i] <= hi[i] on every axis");

  double centre[3];
  double half[3];
  for (int index = 0; index < 3; ++index) {
    centre[index] = 0.5 * (lo[index] + hi[index]);
    half[index] = 0.5 * (hi[index] - lo[index]);
  }
  const double middle = EvalHalfspace(halfspace, centre);

  // The pad below tracks mag, the sum of the MAGNITUDES of the terms EvalHalfspace(centre) adds,
  // not |middle| (the cancelled result). Deep in a subdivision, on a box straddling the surface,
  // middle and halfWidth both go to zero, but the terms that summed to middle do not -- so a pad
  // built from |middle| would collapse to nothing exactly on the boxes whose nActive == 0 Task 6
  // trusts as a hard guarantee. The quadric chain is roughly 6 diagonal and 6 cross-term
  // multiplications plus about 10 additions -- nearer 16 to 26 operations than a clean dozen, so
  // a rigorous Higham bound sits close to 16u; kPadFactor is set well past that edge rather than
  // riding it.
  constexpr double kPadFactor = 32. * std::numeric_limits<double>::epsilon();

  double halfWidth;
  double mag;
  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    // f is the torus's exact signed distance, proved 1-Lipschitz in task 3: |f(x) - f(m)| <=
    // |x - m|, and the farthest point of the box from its own centre is the corner, at distance
    // |h|. sign only flips the sign of the deviation, never its magnitude (|sign| == 1), so the
    // same reach bounds sign*f - sign*f(m) too.
    halfWidth = std::sqrt(half[0] * half[0] + half[1] * half[1] + half[2] * half[2]);
    const double* c = halfspace.c;
    const double offset[3] = {centre[0] - c[0], centre[1] - c[1], centre[2] - c[2]};
    const double along = offset[0] * c[3] + offset[1] * c[4] + offset[2] * c[5];
    const double radial[3] = {offset[0] - along * c[3], offset[1] - along * c[4],
                              offset[2] - along * c[5]};
    const double rho = std::sqrt(radial[0] * radial[0] + radial[1] * radial[1] +
                                 radial[2] * radial[2]);
    // mag does not need a term for centre/c[0..2]'s own scale, however far the torus sits from
    // the origin: point and centre are always within O(major + minor) of each other near the
    // core circle, so by Sterbenz's lemma (subtracting two doubles within a factor of two of one
    // another is exact) offset carries no rounding error regardless of |c[0..2]|. The only
    // cancellation left is local, at rho - major, which rho + |major| already covers -- the
    // magnitudes feeding hypot(rho - major, along) - minor, the terms EvalHalfspace's torus
    // branch actually forms.
    mag = rho + std::abs(c[6]) + std::abs(along) + std::abs(c[7]);
  } else {
    const double* c = halfspace.c;
    const double a[3][3] = {{c[0], c[1], c[2]}, {c[1], c[3], c[4]}, {c[2], c[4], c[5]}};
    const double b[3] = {c[6], c[7], c[8]};
    double slack = 0.;
    mag = std::abs(c[9]);
    for (int row = 0; row < 3; ++row) {
      double gradient = b[row];
      mag += 2. * std::abs(b[row] * centre[row]);
      for (int column = 0; column < 3; ++column) {
        gradient += a[row][column] * centre[column];
        // sum |A_ij| h_i h_j over-estimates the true |Q(x) - Q(m)| cross-term deviation ONLY
        // when every half-extent is >= 0 -- summed over the (i,j)/(j,i) pair, the mis-indexed
        // h_j^2-for-h_i*h_j substitution a reviewer once proposed here exceeds the correct term
        // by |A_ij| (h_i - h_j)^2 >= 0, which needs h_i, h_j real and the inequality direction
        // needs them nonnegative to mean "over", not merely "different". CloseShape's refusal of
        // an inverted cell bbox is what guarantees that for every box reaching this function from
        // SplitBox; a caller going around CloseShape must keep the invariant itself.
        slack += std::abs(a[row][column]) * half[row] * half[column];
        mag += std::abs(a[row][column] * centre[row] * centre[column]);
      }
      slack += 2. * std::abs(gradient) * half[row];
    }
    // slack bounds |Q(x) - Q(m)| (the unsigned Q, built from the unsigned A/b above). middle is
    // sign*Q(m), and sign*Q(x) - sign*Q(m) = sign*(Q(x) - Q(m)), whose magnitude is |Q(x) - Q(m)|
    // because |sign| == 1 -- so the same slack bounds the signed deviation too, on either sign.
    halfWidth = slack;
  }
  // SplitBox's drop tests (rangeLo > 0., rangeHi <= 0.) treat this bound as exact, and Task 6
  // leans on nActive == 0 as a HARD guarantee -- but slack/reach and middle are both floating-
  // point sums, so their true values can round a few ulps off. Only ever widens the enclosure,
  // which is the safe direction.
  halfWidth += kPadFactor * mag;
  rangeLo = middle - halfWidth;
  rangeHi = middle + halfWidth;
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

void O2FlatCSG::SplitBox(int cell, const double* lo, const double* hi,
                         const std::vector<int>& active, int depth, double minSize,
                         int cubifyBudget)
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
  double shortest = TGeoShape::Big();
  int axis = 0;
  for (int index = 0; index < 3; ++index) {
    const double extent = hi[index] - lo[index];
    if (extent > longest) {
      longest = extent;
      axis = index;
    }
    shortest = std::min(shortest, extent);
  }
  // whether THIS box, before any further split, is still far from cubic -- checked before the
  // split it gates, per the header doc comment: a split out of a far-from-cubic box draws from
  // `cubifyBudget` rather than `depth`, so `depth` stays untouched until aspect ratio is within a
  // factor of two on every axis, which is exactly what keeps a near-cubic cell's boxes identical
  // to what a depth-only budget would have produced (design section 4.2 states the invariance
  // this relies on: halving the longest of three extents with `c <= 2a` keeps the new ratio
  // `<= 2`, so a cell that starts near-cubic never becomes far-from-cubic and never touches
  // `cubifyBudget` at all).
  //
  // `shortest` is floored at `minSize` for this test only, not for `keep`'s own `longest <=
  // minSize` check below: a cell bbox with a zero (or merely tiny) extent on one axis -- valid
  // input as far as `CloseShape`'s bbox validation is concerned, which rejects only unset,
  // inverted or non-finite boxes, not degenerate-but-flat ones -- would otherwise leave `shortest`
  // pinned at (near) zero on an axis that is never the longest and so never gets split, making
  // `farFromCubic` permanently true and burning the entire `cubifyBudget` ceiling on a cell that a
  // depth-only budget would have resolved in a handful of splits. Once the other axes have shrunk
  // to `minSize` there is nothing left worth cubifying towards, so the floor lets the box fall
  // through to the ordinary `longest <= minSize` stop instead of exhausting the ceiling first.
  const bool farFromCubic = longest > 2. * std::max(shortest, minSize);
  const bool keep = stillActive.empty() || depth <= 0 || longest <= minSize ||
                    (farFromCubic && cubifyBudget <= 0);
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

  const int childDepth = farFromCubic ? depth : depth - 1;
  const int childCubifyBudget = farFromCubic ? cubifyBudget - 1 : cubifyBudget;
  const double middle = 0.5 * (lo[axis] + hi[axis]);
  double childLo[3] = {lo[0], lo[1], lo[2]};
  double childHi[3] = {hi[0], hi[1], hi[2]};
  childHi[axis] = middle;
  SplitBox(cell, childLo, childHi, stillActive, childDepth, minSize, childCubifyBudget);
  childHi[axis] = hi[axis];
  childLo[axis] = middle;
  SplitBox(cell, childLo, childHi, stillActive, childDepth, minSize, childCubifyBudget);
}

void O2FlatCSG::CloseShape()
{
  fBoxes.clear();
  fActive.clear();
  fClosed = false;
  // dropped before the validation below can return: a BVH left over from an earlier CloseShape
  // would describe boxes that no longer exist, and the queries key off `fBVH != nullptr`
  delete static_cast<BVH*>(fBVH);
  fBVH = nullptr;

  if (static_cast<int>(fCellBBoxSet.size()) < GetNcells()) {
    fCellLo.resize(3 * GetNcells(), 0.);
    fCellHi.resize(3 * GetNcells(), 0.);
    fCellBBoxSet.resize(GetNcells(), false);
  }
  // A cell an intersection of halfspaces does not bound itself, so a cell whose bbox was never
  // handed in by SetCellBBox has no box built for it: it would silently vanish from the solid,
  // Contains() returning false inside real material with nothing to signal it. Fail loudly and
  // build nothing rather than emit a partial, silently-wrong solid. An inverted box (e.g. from a
  // converter that swapped lo/hi) is caught here too: an all-axes-inverted box would otherwise
  // never grow longest past its 0. initialiser in SplitBox, so it would be kept immediately with
  // whatever active list the (invalid, negative-half-extent) range bound happened to compute --
  // possibly nActive == 0, which reads as solid material -- with IsClosed() true and nothing
  // printed. A non-finite (NaN/Inf) bound gets its own check rather than folding into the
  // inverted-box test: every comparison against NaN is false, so `hi < lo` alone would let a NaN
  // box sail through silently and go on to produce a NaN range in HalfspaceRange, which fails
  // both of SplitBox's drop tests and is kept -- possibly as another spuriously "solid"
  // nActive == 0 box. Same treatment as the missing-bbox case: refuse the whole shape.
  bool anyProblem = false;
  for (int cell = 0; cell < GetNcells(); ++cell) {
    if (!fCellBBoxSet[cell]) {
      Error("CloseShape",
            "Shape %s cell %d has no bounding box (SetCellBBox was never called for it); it would "
            "silently vanish from the solid. Not building any boxes -- IsClosed() stays false.",
            GetName(), cell);
      anyProblem = true;
      continue;
    }
    for (int index = 0; index < 3; ++index) {
      const double loValue = fCellLo[3 * cell + index];
      const double hiValue = fCellHi[3 * cell + index];
      if (!std::isfinite(loValue) || !std::isfinite(hiValue)) {
        Error("CloseShape",
              "Shape %s cell %d has a non-finite bounding box on axis %d (lo %g, hi %g). Not "
              "building any boxes -- IsClosed() stays false.",
              GetName(), cell, index, loValue, hiValue);
        anyProblem = true;
        continue;
      }
      if (hiValue < loValue) {
        Error("CloseShape",
              "Shape %s cell %d has an inverted bounding box on axis %d (lo %g > hi %g); "
              "SetCellBBox's arguments look swapped. Not building any boxes -- IsClosed() stays "
              "false.",
              GetName(), cell, index, loValue, hiValue);
        anyProblem = true;
      }
    }
  }
  if (anyProblem) {
    return;
  }

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

#ifndef NDEBUG
  // A cell's bounding box must CONTAIN the cell, and since task 5 that is a correctness
  // precondition rather than a cost one: boxes are only ever built inside it, so a cell that
  // spills past its own bbox is material the accelerated Contains cannot find and the accelerated
  // distances truncate, while the twins -- which know nothing of boxes -- still see it. The
  // converter owes this (design section 4.2); the check samples each face from just outside and
  // requires the cell to have ended.
  {
    const double reach = 1.e-6 * (diagonal > 0. ? diagonal : 1.);
    for (int cell = 0; cell < GetNcells(); ++cell) {
      const double* cellLo = &fCellLo[3 * cell];
      const double* cellHi = &fCellHi[3 * cell];
      for (int axis = 0; axis < 3; ++axis) {
        const int first = (axis + 1) % 3;
        const int second = (axis + 2) % 3;
        for (int side = 0; side < 2; ++side) {
          for (int step1 = 0; step1 <= 4; ++step1) {
            for (int step2 = 0; step2 <= 4; ++step2) {
              double probe[3];
              probe[axis] = side == 0 ? cellLo[axis] - reach : cellHi[axis] + reach;
              probe[first] = cellLo[first] + 0.25 * step1 * (cellHi[first] - cellLo[first]);
              probe[second] = cellLo[second] + 0.25 * step2 * (cellHi[second] - cellLo[second]);
              assert(!CellContains(cell, probe) &&
                     "O2FlatCSG::CloseShape: a cell reaches past the bounding box SetCellBBox was "
                     "given; the accelerated queries build boxes only inside it and would lose that "
                     "material");
            }
          }
        }
      }
    }
  }
#endif

  for (int cell = 0; cell < GetNcells(); ++cell) {
    std::vector<int> active;
    active.reserve(fCells[cell].count);
    for (int offset = 0; offset < fCells[cell].count; ++offset) {
      active.push_back(fCells[cell].first + offset);
    }
    SplitBox(cell, &fCellLo[3 * cell], &fCellHi[3 * cell], active, fSplitDepth, minSize,
             kMaxCubifySplits);
  }

  if (!fBoxes.empty()) {
    std::vector<BVHBBox> boxes;
    std::vector<BVHVec3> centers;
    boxes.reserve(fBoxes.size());
    centers.reserve(fBoxes.size());
    for (const auto& box : fBoxes) {
      BVHBBox bounds;
      for (int index = 0; index < 3; ++index) {
        // outward, so a float node box is a superset of the double box it stands for and the
        // traversal can only ever nominate too many candidates -- never drop one
        bounds.min[index] = roundOutward(box.min[index], false);
        bounds.max[index] = roundOutward(box.max[index], true);
      }
      boxes.push_back(bounds);
      centers.emplace_back(bounds.get_center());
    }
    typename bvh::v2::DefaultBuilder<BVHNode>::Config config;
    config.quality = bvh::v2::DefaultBuilder<BVHNode>::Quality::High;
    // One box per leaf, as in O2BVHAssembly: resolving a box costs a slab clip plus an interval
    // scan over its active list, far more than a node box test, and the bvh2 traversal enters a
    // leaf's start node without testing its box -- so a multi-box leaf would pay several clips
    // for one box test. It also means a box lives in exactly one leaf and is visited at most once
    // per traversal, which is what lets the per-cell merge below assume no duplicate pieces.
    config.max_leaf_size = 1;
    fBVH = static_cast<void*>(
      new BVH(bvh::v2::DefaultBuilder<BVHNode>::build(boxes, centers, config)));
  }

  fClosed = true;
  ComputeBBox();
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

////////////////////////////////////////////////////////////////////////////////
/// Contains -- a BVH point query.
///
/// Inside a box, the box's active list IS the cell: every halfspace the split dropped holds
/// everywhere in that box, so testing the survivors is testing all of them. That is why the point
/// has to be in the box's own double bounds before its list is consulted -- the float node box the
/// traversal descends is a deliberate superset and proves nothing.
///
/// Cells are disjoint by construction, so the first box that accepts is the answer and the visit
/// order does not matter.

Bool_t O2FlatCSG::Contains(const Double_t* point) const
{
  if (!fClosed || fBVH == nullptr) {
    // CloseShape refused (or was never called), so there are no boxes to walk. An empty box array
    // in an accelerated path reads as "no material anywhere", which is exactly the silent
    // vanishing CloseShape's validation exists to prevent -- so answer from the twin instead.
    return Contains_Loop(point);
  }
  const BVH& bvh = *static_cast<const BVH*>(fBVH);
  const BVHVec3 query(static_cast<float>(point[0]), static_cast<float>(point[1]),
                      static_cast<float>(point[2]));

  thread_local std::vector<size_t> stack;
  stack.clear();
  stack.push_back(0); // the bvh2 root node
  while (!stack.empty()) {
    const size_t current = stack.back();
    stack.pop_back();
    const auto& node = bvh.nodes[current];
    if (!bvh::v2::extra::contains(node.get_bbox(), query)) {
      continue;
    }
    if (node.is_leaf()) {
      const auto beginPrimitive = node.index.first_id();
      const auto endPrimitive = beginPrimitive + node.index.prim_count();
      for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
        const FlatCSGBox& box = fBoxes[bvh.prim_ids[primitive]];
        if (!boxHoldsPoint(box, point)) {
          continue;
        }
        if (box.nActive == 0) {
          return kTRUE; // wholly inside its cell: nothing left to test
        }
        bool inside = true;
        for (int slot = 0; slot < box.nActive && inside; ++slot) {
          inside = EvalHalfspace(fHalfspaces[fActive[box.firstActive + slot]], point) <= 0.;
        }
        if (inside) {
          return kTRUE;
        }
      }
    } else {
      const auto firstChild = node.index.first_id();
      for (size_t child : {firstChild, firstChild + 1}) {
        if (child < bvh.nodes.size()) {
          stack.push_back(child);
        }
      }
    }
  }
  return kFALSE;
}

int O2FlatCSG::HalfspaceRoots(const FlatCSGHalfspace& halfspace, const double* origin,
                              const double* dir, double* roots)
{
  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    // the quartic derivation below takes the leading coefficient a4 = |dir|^4 to be exactly 1;
    // a non-unit direction silently returns wrong roots instead of failing, so catch it here
    assert(std::abs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] - 1.) < 1.e-9 &&
          "O2FlatCSG::HalfspaceRoots: torus branch requires a unit direction");
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
    // solveQuarticReal is scale-normalised internally, so the torus branch needs no degeneracy
    // guard of its own -- unlike the quadric path below, whose alpha threshold is a substitute
    // for exactly that normalisation
    const std::vector<double> found = o2::base::surface::solveQuarticReal(a4, a3, a2, a1, a0);
    int count = 0;
    for (double value : found) {
      if (count < kMaxRootsPerHalfspace) {
        roots[count++] = value;
      }
    }
    return count;
  }
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
  const double gamma = EvalHalfspace(halfspace, origin) * halfspace.sign; // sign*sign==1: the unsigned Q(o)

  // scale-aware degeneracy: a plane has alpha exactly 0, and a ray parallel to a cylinder axis
  // has it near 0 -- both are the linear equation, not a badly conditioned quadratic
  //
  // NOTE: 1.e-14 is dimensionally an accident, not a derivation: alpha has units [A], beta
  // [A*L], gamma [A*L^2], so this comparison depends on the length unit. It assumes coordinates
  // in cm (ALICE's native unit) -- a root this branch discards sits at |t| >= ~1e6 cm, outside
  // any ALICE geometry. Task 3's torus branch must not copy this constant blind.
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
  if (q == 0.) {
    // q == 0 happens only when beta == 0 and gamma == 0 together (origin exactly on the
    // surface, direction exactly tangential): disc == 0 too, so the quadratic has one double
    // root at t = -beta/alpha = 0. The general second-root formula would divide 0./0. here.
    roots[0] = 0.;
    return 1;
  }
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

  // every root of every active halfspace, clipped to the window; the buffer is sized from this
  // cell's own halfspace count, worst case, so no root is ever silently dropped. thread_local so
  // that a std::vector resize on one thread can never race a data() held by another -- TGeo
  // shares one shape object across every navigator under TGeoManager::SetMaxThreads.
  thread_local std::vector<double> breakBuffer;
  const std::size_t needed = 2 + static_cast<std::size_t>(kMaxRootsPerHalfspace) * static_cast<std::size_t>(count);
  if (breakBuffer.size() < needed) {
    breakBuffer.resize(needed);
  }
  double* breaks = breakBuffer.data();
  int nBreaks = 0;
  breaks[nBreaks++] = tlo;
  breaks[nBreaks++] = thi;
  for (int slot = 0; slot < count; ++slot) {
    const int index = active != nullptr ? active[slot] : description.first + slot;
    double roots[kMaxRootsPerHalfspace];
    const int found = HalfspaceRoots(fHalfspaces[index], origin, dir, roots);
    for (int root = 0; root < found; ++root) {
      if (roots[root] > tlo && roots[root] < thi) {
        breaks[nBreaks++] = roots[root];
      }
    }
  }
  std::sort(breaks, breaks + nBreaks);

  // classify the midpoint of each sub-interval and merge the runs that are inside
  int pairs = 0;
  bool open = false;
  bool overflow = false;
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
      } else {
        // maxOut was too small for this cell along this ray: fail loudly (a negative count)
        // rather than hand the caller a silently truncated list that reads as a valid answer
        overflow = true;
        open = false;
      }
    } else {
      open = false;
    }
  }
  return overflow ? -1 : pairs;
}

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
  // thread_local: see the comment on the scratch-buffer members it replaced in the header
  thread_local std::vector<double> pairBuffer;
  double best = TGeoShape::Big();
  for (int cell = 0; cell < GetNcells(); ++cell) {
    // sized from this cell's own halfspace count, so a busy cell's intervals are never truncated
    const int capacity = maxPairsForCell(fCells[cell].count);
    if (static_cast<int>(pairBuffer.size()) < 2 * capacity) {
      pairBuffer.resize(2 * capacity);
    }
    const int found = CellIntervals(cell, nullptr, -1, point, dir, 0., step,
                                    pairBuffer.data(), capacity);
    // capacity is provably sufficient (maxPairsForCell), so CellIntervals cannot overflow here;
    // a negative found would mean that bound itself is wrong, which is a bug, not live data
    for (int pair = 0; pair < found; ++pair) {
      // a point exactly on the boundary is already inside; only a real entry counts
      if (pairBuffer[2 * pair + 1] > TGeoShape::Tolerance() && pairBuffer[2 * pair] < best) {
        best = std::max(pairBuffer[2 * pair], 0.);
      }
    }
  }
  return best;
}

Double_t O2FlatCSG::DistFromInside_Loop(const Double_t* point, const Double_t* dir,
                                        Double_t step) const
{
  // the union's occupancy, so a ray that leaves one cell into a touching one keeps going; the
  // buffer is sized to fit every cell's worst case at once, so nothing gathered here is ever
  // lost. thread_local: see the comment on the scratch-buffer members it replaced in the header
  thread_local std::vector<double> pairBuffer;
  int totalCapacity = 0;
  for (int cell = 0; cell < GetNcells(); ++cell) {
    totalCapacity += maxPairsForCell(fCells[cell].count);
  }
  if (static_cast<int>(pairBuffer.size()) < 2 * totalCapacity) {
    pairBuffer.resize(2 * totalCapacity);
  }
  int count = 0;
  for (int cell = 0; cell < GetNcells(); ++cell) {
    // each call's capacity (totalCapacity - count) is at least that cell's own maxPairsForCell,
    // so this is not expected to overflow -- but a negative count must never reach the pointer
    // arithmetic below, so check it explicitly rather than trust the bound silently
    const int found = CellIntervals(cell, nullptr, -1, point, dir, 0., step,
                                    pairBuffer.data() + 2 * count, totalCapacity - count);
    if (found < 0) {
      Error("DistFromInside_Loop",
            "CellIntervals overflowed for cell %d: the maxPairsForCell bound no longer holds",
            cell);
      return TGeoShape::Big();
    }
    count += found;
  }
  count = mergeIntervals(pairBuffer.data(), count, TGeoShape::Tolerance());
  for (int pair = 0; pair < count; ++pair) {
    if (pairBuffer[2 * pair] <= TGeoShape::Tolerance()) {
      return pairBuffer[2 * pair + 1];
    }
  }
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// GatherRayPieces -- the clip-inside-the-box invariant, in one place.
///
/// Every path into the accelerated distances comes through here, and here the window handed to
/// `CellIntervals` is always this box's own slab intersected with `[0, step]`. No active list is
/// ever used outside the box it was computed for, and no window is ever pooled across boxes.

bool O2FlatCSG::GatherRayPieces(const Double_t* point, const Double_t* dir, Double_t step,
                                std::vector<double>& pairs, std::vector<int>& cells) const
{
  pairs.clear();
  cells.clear();
  const BVH& bvh = *static_cast<const BVH*>(fBVH);

  // one box's intervals; thread_local for the reason the header's scratch-buffer comment gives
  thread_local std::vector<double> boxPairs;
  bool overflowed = false;
  traverseRay(bvh, point, dir, step, [&](int index) {
    const FlatCSGBox& box = fBoxes[index];
    double tlo = 0.;
    double thi = step;
    if (!slabWindow(box.min, box.max, point, dir, tlo, thi) || thi <= tlo) {
      return;
    }
    // sized from THIS box's active-list length, which is the count CellIntervals will walk, so
    // the bound it is asked to respect is the one it was given
    const int capacity = maxPairsForCell(box.nActive);
    if (static_cast<int>(boxPairs.size()) < 2 * capacity) {
      boxPairs.resize(2 * capacity);
    }
    // nActive == 0 means the box is wholly inside its cell; CellIntervals then has no halfspace
    // to break on and returns the whole window, which is exactly the right answer
    const int* active = box.nActive > 0 ? fActive.data() + box.firstActive : nullptr;
    const int found = CellIntervals(box.cell, active, box.nActive, point, dir, tlo, thi,
                                    boxPairs.data(), capacity);
    if (found < 0) {
      overflowed = true;
      return;
    }
    for (int pair = 0; pair < found; ++pair) {
      pairs.push_back(boxPairs[2 * pair]);
      pairs.push_back(boxPairs[2 * pair + 1]);
      cells.push_back(box.cell);
    }
  });
  return !overflowed;
}

////////////////////////////////////////////////////////////////////////////////
/// DistFromOutsideBVH
///
/// The twin's rule is "the smallest entry over the per-CELL occupancy intervals whose exit clears
/// the tolerance", so the pieces are rejoined per cell before it is applied. Boxes of one cell
/// tile it, so a piece that a box boundary cut short abuts the next piece exactly -- the slab
/// bound of one box and of its neighbour are the same expression on the same shared face
/// coordinate -- and the run merge below closes them back into the interval the twin computes.
///
/// Merging across cells would NOT be the twin: two cells that touch within the tolerance are one
/// interval after a merge but two intervals to the twin, and the exit-clears-tolerance test can
/// then pick a different entry.

Double_t O2FlatCSG::DistFromOutsideBVH(const Double_t* point, const Double_t* dir,
                                       Double_t step) const
{
  thread_local std::vector<double> pairs;
  thread_local std::vector<int> cells;
  if (!GatherRayPieces(point, dir, step, pairs, cells)) {
    Error("DistFromOutside",
          "Shape %s: CellIntervals overflowed a per-box buffer sized from that box's own active "
          "list; the maxPairsForCell bound no longer holds. Answering from the loop twin.",
          GetName());
    return DistFromOutside_Loop(point, dir, step);
  }

  // sort the pieces by (cell, entry) through a permutation, so the run merge below sees each
  // cell's pieces contiguously and in order
  const int count = static_cast<int>(cells.size());
  thread_local std::vector<int> order;
  order.resize(count);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int left, int right) {
    if (cells[left] != cells[right]) {
      return cells[left] < cells[right];
    }
    return pairs[2 * left] < pairs[2 * right];
  });

  double best = TGeoShape::Big();
  int index = 0;
  while (index < count) {
    const int cell = cells[order[index]];
    const double enter = pairs[2 * order[index]];
    double exit = pairs[2 * order[index] + 1];
    ++index;
    // join what is only one interval of this cell, cut into pieces by the boxes that tile it
    while (index < count && cells[order[index]] == cell && pairs[2 * order[index]] <= exit) {
      exit = std::max(exit, pairs[2 * order[index] + 1]);
      ++index;
    }
    // DistFromOutside_Loop's rule, unchanged: a point exactly on the boundary is already inside,
    // so only an interval that really extends past the tolerance counts as an entry
    if (exit > TGeoShape::Tolerance() && enter < best) {
      best = std::max(enter, 0.);
    }
  }
  return best;
}

////////////////////////////////////////////////////////////////////////////////
/// DistFromInsideBVH
///
/// The far end of the UNION's occupancy interval containing t = 0: a point can leave one cell into
/// an adjacent one, and the part ends where the union does. So unlike DistFromOutside this one
/// merges everything the ray met, cells included, with the same glue tolerance the twin uses.

Double_t O2FlatCSG::DistFromInsideBVH(const Double_t* point, const Double_t* dir,
                                      Double_t step) const
{
  thread_local std::vector<double> pairs;
  thread_local std::vector<int> cells;
  if (!GatherRayPieces(point, dir, step, pairs, cells)) {
    Error("DistFromInside",
          "Shape %s: CellIntervals overflowed a per-box buffer sized from that box's own active "
          "list; the maxPairsForCell bound no longer holds. Answering from the loop twin.",
          GetName());
    return DistFromInside_Loop(point, dir, step);
  }
  const int count = mergeIntervals(pairs.data(), static_cast<int>(cells.size()),
                                   TGeoShape::Tolerance());
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
  if (!fClosed || fBVH == nullptr) {
    // no boxes to walk: see the note on Contains. The twin is the definition of the answer, and
    // an empty box array in the accelerated path would silently report empty space.
    return DistFromOutside_Loop(point, dir, step);
  }
  return DistFromOutsideBVH(point, dir, step);
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
  if (!fClosed || fBVH == nullptr) {
    return DistFromInside_Loop(point, dir, step);
  }
  return DistFromInsideBVH(point, dir, step);
}

////////////////////////////////////////////////////////////////////////////////
/// Safety_Loop -- the definition of the answer, from the box structure alone.
///
/// Outside: every point of the solid lies in some retained box, so the distance to the nearest
/// box (over ALL of them, decided or not) is a lower bound on the distance to the solid.
///
/// Inside: only a box with an empty active list is a sound bound, because it is wholly inside its
/// cell (a hard guarantee established by Task 4). An undecided box may carry boundary anywhere
/// within it, so it contributes nothing; a point that is not covered by any solid-marked box gets
/// the sound answer `0.`.

Double_t O2FlatCSG::Safety_Loop(const Double_t* point, Bool_t in) const
{
  if (!in) {
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

  double best = 0.;
  for (const auto& box : fBoxes) {
    if (!boxHoldsPoint(box, point) || box.nActive != 0) {
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

////////////////////////////////////////////////////////////////////////////////
/// Safety -- the same box-structure computation as `Safety_Loop`, sped up with the BVH.
///
/// Outside, this is a branch-and-bound nearest-box search: `bvh::v2::extra::SafetySqToNode`
/// prunes a subtree once its own box-to-point squared distance cannot beat the standing best, and
/// a leaf is scored with the exact same double-precision box-to-point formula the twin uses, so
/// the two compute the identical `min` over a set the pruning never shrinks past the true nearest
/// box -- bit identity is the check that no nearer box was pruned away.
///
/// Inside, this is the same `Contains`-style containment descent: `bvh::v2::extra::contains` on
/// the (outward-rounded, float) node box only ever nominates too many candidates, and each
/// nominated box is then re-tested against its own double bounds and `nActive == 0` exactly as
/// the twin does, so the `max` runs over the same qualifying set either way.

Double_t O2FlatCSG::Safety(const Double_t* point, Bool_t in) const
{
  if (!fClosed || fBVH == nullptr) {
    // no boxes to walk: see the note on Contains -- the twin is the definition of the answer.
    return Safety_Loop(point, in);
  }
  const BVH& bvh = *static_cast<const BVH*>(fBVH);
  const BVHVec3 query(static_cast<float>(point[0]), static_cast<float>(point[1]),
                      static_cast<float>(point[2]));

  thread_local std::vector<size_t> stack;
  stack.clear();
  stack.push_back(0); // the bvh2 root node

  if (!in) {
    // The prune test needs a genuine LOWER bound on the true (double) distance to whatever real
    // box sits under a node, or it can throw away the nearest one -- which is exactly what
    // happens if both the node box and the query point are rounded to float first: the query
    // rounds by up to half a float ulp, and that alone can push a float distance above the exact
    // double distance to the box that holds the true minimum, no matter how the box itself is
    // padded. So this reads the node's (float, outward-rounded) box back as double -- an exact
    // upcast, still a superset of every double box under it -- and measures it against the
    // point's own double coordinates, with the same SafetySqToNode kernel bvh2_extra_kernels.h
    // gives the leaf test, just instantiated at double instead of float.
    using DVec3 = bvh::v2::Vec<double, 3>;
    using DBBox = bvh::v2::BBox<double, 3>;
    const DVec3 dpoint(point[0], point[1], point[2]);
    double best = TGeoShape::Big();
    while (!stack.empty()) {
      const size_t current = stack.back();
      stack.pop_back();
      const auto& node = bvh.nodes[current];
      const auto& fbox = node.get_bbox();
      const DBBox dbox(DVec3(static_cast<double>(fbox.min[0]), static_cast<double>(fbox.min[1]),
                             static_cast<double>(fbox.min[2])),
                       DVec3(static_cast<double>(fbox.max[0]), static_cast<double>(fbox.max[1]),
                             static_cast<double>(fbox.max[2])));
      const double nodeSquared = bvh::v2::extra::SafetySqToNode(dbox, dpoint);
      if (nodeSquared >= best) {
        continue; // this subtree cannot hold anything nearer than what is already found
      }
      if (node.is_leaf()) {
        const auto beginPrimitive = node.index.first_id();
        const auto endPrimitive = beginPrimitive + node.index.prim_count();
        for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
          const FlatCSGBox& box = fBoxes[bvh.prim_ids[primitive]];
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
      } else {
        const auto firstChild = node.index.first_id();
        for (size_t child : {firstChild, firstChild + 1}) {
          if (child < bvh.nodes.size()) {
            stack.push_back(child);
          }
        }
      }
    }
    return best >= TGeoShape::Big() ? 0. : std::sqrt(best);
  }

  double best = 0.;
  while (!stack.empty()) {
    const size_t current = stack.back();
    stack.pop_back();
    const auto& node = bvh.nodes[current];
    if (!bvh::v2::extra::contains(node.get_bbox(), query)) {
      continue;
    }
    if (node.is_leaf()) {
      const auto beginPrimitive = node.index.first_id();
      const auto endPrimitive = beginPrimitive + node.index.prim_count();
      for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
        const FlatCSGBox& box = fBoxes[bvh.prim_ids[primitive]];
        if (!boxHoldsPoint(box, point) || box.nActive != 0) {
          continue;
        }
        double toFace = TGeoShape::Big();
        for (int index = 0; index < 3; ++index) {
          toFace = std::min(toFace, std::min(point[index] - box.min[index], box.max[index] - point[index]));
        }
        best = std::max(best, toFace);
      }
    } else {
      const auto firstChild = node.index.first_id();
      for (size_t child : {firstChild, firstChild + 1}) {
        if (child < bvh.nodes.size()) {
          stack.push_back(child);
        }
      }
    }
  }
  return std::max(best, 0.);
}

Double_t O2FlatCSG::Capacity() const
{
  // the cells of a decomposition are disjoint by construction, so their own volumes just sum
  return std::accumulate(fCells.begin(), fCells.end(), 0.,
                         [](double sum, const FlatCSGCell& cell) { return sum + cell.volume; });
}

////////////////////////////////////////////////////////////////////////////////
/// ComputeNormal -- the gradient of whichever halfspace is closest to being satisfied with
/// equality at \a point, oriented per the `TGeoShape` convention (with respect to \a dir).
///
/// Two things a naive global argmin over `|EvalHalfspace|` gets wrong, both fixed here per design
/// section 5.5:
///
/// 1. `|EvalHalfspace|` is an ALGEBRAIC value, not a distance, and its gain per unit distance
///    differs per halfspace -- 1 for a unit plane, `2R` for a cylinder of radius `R`, and whatever
///    a rescaled plane's coefficients happen to carry (design section 3.1). A point a hair off an
///    R = 100 cylinder's wall can have `|f|` two orders of magnitude larger than a plane it is
///    genuinely much farther from, so an unscaled argmin picks the far plane. The fix is to select
///    on the first-order distance `|f| / |grad f|` instead, which both halves of this function
///    already compute.
/// 2. Comparing across EVERY halfspace of EVERY cell lets a distant cell's halfspace win at all.
///    Restricting the scan to the active list of the box containing `point` -- design section
///    5.5's own words -- removes that false winner as well as most of the cost. This is sound for
///    a stronger reason than "a dropped halfspace cannot be the boundary": `HalfspaceRange` pads
///    its enclosure outward by `32 * eps * mag` (see its own doc comment), so the `rangeHi <= 0`
///    test that drops a halfspace from a box's active list proves the true maximum of `sign * f`
///    over the box is STRICTLY negative, not merely non-positive -- no point of the box lies on a
///    dropped halfspace's surface at all, so the active list is not just the halfspaces still
///    undecided there (section 4.2), it is provably every halfspace that could possibly be at
///    equality anywhere in the box.

void O2FlatCSG::ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) const
{
  norm[0] = norm[1] = norm[2] = 0.;
  if (fHalfspaces.empty()) {
    return;
  }

  const int* candidates = nullptr;
  int nCandidates = 0;
  // a full cell's halfspace run, built only when the containing box turns out to be wholly
  // inside its cell (nActive == 0): such a box's OWN active list is empty by construction, but a
  // point genuinely on a boundary cannot really be in one (see the reasoning at the call site
  // below), so this is a defensive fallback rather than a path real navigation should reach
  thread_local std::vector<int> wholeCell;

  if (fClosed && fBVH != nullptr) {
    const BVH& bvh = *static_cast<const BVH*>(fBVH);
    const BVHVec3 query(static_cast<float>(point[0]), static_cast<float>(point[1]),
                        static_cast<float>(point[2]));
    thread_local std::vector<size_t> stack;
    stack.clear();
    stack.push_back(0); // the bvh2 root node
    while (!stack.empty() && candidates == nullptr) {
      const size_t current = stack.back();
      stack.pop_back();
      const auto& node = bvh.nodes[current];
      if (!bvh::v2::extra::contains(node.get_bbox(), query)) {
        continue;
      }
      if (node.is_leaf()) {
        const auto beginPrimitive = node.index.first_id();
        const auto endPrimitive = beginPrimitive + node.index.prim_count();
        for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
          const FlatCSGBox& box = fBoxes[bvh.prim_ids[primitive]];
          if (!boxHoldsPoint(box, point)) {
            continue;
          }
          if (box.nActive > 0) {
            candidates = fActive.data() + box.firstActive;
            nCandidates = box.nActive;
          } else {
            const FlatCSGCell& description = fCells[box.cell];
            wholeCell.resize(description.count);
            for (int slot = 0; slot < description.count; ++slot) {
              wholeCell[slot] = description.first + slot;
            }
            candidates = wholeCell.data();
            nCandidates = static_cast<int>(wholeCell.size());
          }
          break; // cells are disjoint; the first box that holds the point is the answer
        }
      } else {
        const auto firstChild = node.index.first_id();
        for (size_t child : {firstChild, firstChild + 1}) {
          if (child < bvh.nodes.size()) {
            stack.push_back(child);
          }
        }
      }
    }
  }

  // no box claimed the point -- the shape never closed, or the point sits outside every retained
  // box -- so fall back to every halfspace of every cell. `ComputeNormal` has no `_Loop` twin of
  // its own to hand this off to; this scan IS that fallback.
  thread_local std::vector<int> allHalfspaces;
  if (candidates == nullptr) {
    allHalfspaces.resize(GetNhalfspaces());
    for (int index = 0; index < GetNhalfspaces(); ++index) {
      allHalfspaces[index] = index;
    }
    candidates = allHalfspaces.data();
    nCandidates = GetNhalfspaces();
  }

  int best = -1;
  double bestValue = std::numeric_limits<double>::infinity();
  double bestGrad[3] = {0., 0., 0.};
  for (int slot = 0; slot < nCandidates; ++slot) {
    const FlatCSGHalfspace& halfspace = fHalfspaces[candidates[slot]];
    const double f = EvalHalfspace(halfspace, point);
    double grad[3];
    halfspaceGradient(halfspace, point, grad);
    const double gradLength = std::sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
    if (gradLength < 1.e-300) {
      continue; // degenerate gradient (see halfspaceGradient); this halfspace cannot win
    }
    const double value = std::abs(f) / gradLength; // the first-order distance to this surface
    if (value < bestValue) {
      bestValue = value;
      best = candidates[slot];
      bestGrad[0] = grad[0] / gradLength;
      bestGrad[1] = grad[1] / gradLength;
      bestGrad[2] = grad[2] / gradLength;
    }
  }

  if (best < 0) {
    // every candidate had a degenerate gradient -- only possible on a torus's own core-circle
    // singularity or its revolution axis; fall back to the travel direction, at least a unit
    // vector oriented the way the contract requires
    const double dirLength = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    if (dirLength > 1.e-300) {
      for (int index = 0; index < 3; ++index) {
        norm[index] = dir[index] / dirLength;
      }
    }
    return;
  }

  for (int index = 0; index < 3; ++index) {
    norm[index] = bestGrad[index];
  }
  const double dot = norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2];
  if (dot < 0.) {
    for (int index = 0; index < 3; ++index) {
      norm[index] = -norm[index];
    }
  }
}

void O2FlatCSG::ComputeBBox()
{
  // the union of the retained sub-cell boxes -- tighter than the union of the cell AABBs, which
  // is the point of design section 5.5
  if (fBoxes.empty()) {
    return;
  }
  double lo[3] = {TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};
  double hi[3] = {-TGeoShape::Big(), -TGeoShape::Big(), -TGeoShape::Big()};
  for (const FlatCSGBox& box : fBoxes) {
    for (int index = 0; index < 3; ++index) {
      lo[index] = std::min(lo[index], box.min[index]);
      hi[index] = std::max(hi[index], box.max[index]);
    }
  }
  for (int index = 0; index < 3; ++index) {
    fOrigin[index] = 0.5 * (lo[index] + hi[index]);
  }
  fDX = 0.5 * (hi[0] - lo[0]);
  fDY = 0.5 * (hi[1] - lo[1]);
  fDZ = 0.5 * (hi[2] - lo[2]);
}

} // namespace base
} // namespace o2
