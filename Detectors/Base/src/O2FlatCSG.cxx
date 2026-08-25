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

#include "TGeoShape.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>

ClassImp(o2::base::O2FlatCSG);

/// The most roots one cell can contribute to one ray: four per torus halfspace.
constexpr int kMaxRootsPerHalfspace = 4;

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
  double centre[3];
  double half[3];
  for (int index = 0; index < 3; ++index) {
    centre[index] = 0.5 * (lo[index] + hi[index]);
    half[index] = 0.5 * (hi[index] - lo[index]);
  }
  const double middle = EvalHalfspace(halfspace, centre);

  if (halfspace.kind == FlatCSGHalfspace::kTorus) {
    // f is the torus's exact signed distance, proved 1-Lipschitz in task 3: |f(x) - f(m)| <=
    // |x - m|, and the farthest point of the box from its own centre is the corner, at distance
    // |h|. sign only flips the sign of the deviation, never its magnitude (|sign| == 1), so the
    // same reach bounds sign*f - sign*f(m) too.
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
  // slack bounds |Q(x) - Q(m)| (the unsigned Q, built from the unsigned A/b above). middle is
  // sign*Q(m), and sign*Q(x) - sign*Q(m) = sign*(Q(x) - Q(m)), whose magnitude is |Q(x) - Q(m)|
  // because |sign| == 1 -- so the same slack bounds the signed deviation too, on either sign.
  rangeLo = middle - slack;
  rangeHi = middle + slack;
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
  fClosed = false;

  if (static_cast<int>(fCellBBoxSet.size()) < GetNcells()) {
    fCellLo.resize(3 * GetNcells(), 0.);
    fCellHi.resize(3 * GetNcells(), 0.);
    fCellBBoxSet.resize(GetNcells(), false);
  }
  // A cell an intersection of halfspaces does not bound itself, so a cell whose bbox was never
  // handed in by SetCellBBox has no box built for it: it would silently vanish from the solid,
  // Contains() returning false inside real material with nothing to signal it. Fail loudly and
  // build nothing rather than emit a partial, silently-wrong solid.
  bool anyMissing = false;
  for (int cell = 0; cell < GetNcells(); ++cell) {
    if (!fCellBBoxSet[cell]) {
      Error("CloseShape",
            "Shape %s cell %d has no bounding box (SetCellBBox was never called for it); it would "
            "silently vanish from the solid. Not building any boxes -- IsClosed() stays false.",
            GetName(), cell);
      anyMissing = true;
    }
  }
  if (anyMissing) {
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

Double_t O2FlatCSG::Safety(const Double_t* /*point*/, Bool_t /*in*/) const
{
  // task 6 replaces this; 0 is always sound
  return 0.;
}

void O2FlatCSG::ComputeBBox()
{
  // Minimal override, just enough for CloseShape to compile: the union of the retained boxes.
  // Task 6 writes the real one.
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
