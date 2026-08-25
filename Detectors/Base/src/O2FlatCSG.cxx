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

#include "TGeoShape.h"

#include <algorithm>
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

} // namespace base
} // namespace o2
