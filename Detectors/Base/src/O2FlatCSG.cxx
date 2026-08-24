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
