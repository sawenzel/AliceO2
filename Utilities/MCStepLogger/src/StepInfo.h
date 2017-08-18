// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//****************************************************************************
//* This file is free software: you can redistribute it and/or modify        *
//* it under the terms of the GNU General Public License as published by     *
//* the Free Software Foundation, either version 3 of the License, or        *
//* (at your option) any later version.                                      *
//*                                                                          *
//* Primary Authors: Sandro Wenzel <sandro.wenzel@cern.ch>                   *
//*                                                                          *
//* The authors make no claims about the suitability of this software for    *
//* any purpose. It is provided "as is" without express or implied warranty. *
//****************************************************************************

#ifndef O2_STEPINFO
#define O2_STEPINFO

#include <Rtypes.h>
#include <chrono>
 
class TVirtualMC;

namespace o2
{
// class collecting info about one MC step done
class StepInfo {
 public:

  StepInfo() = default;
  // construct directly using virtual mc
  StepInfo(TVirtualMC *mc);  
//  ~StepInfo() { if(secondaryprocesses) { delete[] secondaryprocesses; } }

  int   stepid;
  int   volId;
  int   trackID; // reduntant
  int   pdg; // reduntant
  float x;
  float y;
  float z;
  float E;
  float step;
  float maxstep;
  float cputimestamp;
  int   nsecondaries;
  int*  secondaryprocesses; //[nsecondaries]

  static int stepcounter;

  ClassDefNV(StepInfo, 1);
};

}

#endif
