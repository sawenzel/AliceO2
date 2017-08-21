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
class TGeoVolume;
class TGeoMedium;

namespace o2
{
// class collecting info about one MC step done

struct VolInfoContainer {
  VolInfoContainer() = default;

  // keeps info about volumes (might exist somewhere else already)
  // essentially a mapping from volumeID x copyNo to TGeoVolumes
  std::vector<std::vector<TGeoVolume const *> *> volumes; // sparse container

  void insert(int id, int copyNo, TGeoVolume const *vol) {
    if(volumes.size() <= id) {
      volumes.resize(id + 1, nullptr);
    }
    if (volumes[id]==nullptr) {
      volumes[id] = new std::vector<TGeoVolume const *>;
    }
    if (volumes[id]->size() <= copyNo) {
      volumes[id]->resize(copyNo + 1, nullptr);
    }
    (*volumes[id])[copyNo] = vol;
  }

  TGeoVolume const * get(int id, int copy) const {return (*volumes[id])[copy];}

  ClassDefNV(VolInfoContainer, 1);
};

struct StepInfo {
  StepInfo() = default;
  // construct directly using virtual mc
  StepInfo(TVirtualMC *mc);  

  long  cputimestamp;
  int   stepid;
  int   volId; // keep another branch somewhere mapping this to name, medium, etc.
  int   copyNo;
  int   trackID;
  int   pdg;
  float x;
  float y;
  float z;
  float E;
  float step;
  float maxstep;
  int   nsecondaries;
  int*  secondaryprocesses; //[nsecondaries]
  int   nprocessesactive; // number of active processes
  bool  stopped; //

  // somehow I can't serialize this:
  // TGeoVolume* geovolume; //->
  std::string volname;
//  TGeoMedium* medium; //->
  std::string mediumname;

  char const* getVolName();
  char const* getMediumName();

  // for cuts?
  bool isVolume(const char *);

  static int stepcounter;
  static std::chrono::time_point<std::chrono::high_resolution_clock> starttime;
//  static VolInfoContainer volinfos;

  ClassDefNV(StepInfo, 2);
};

struct MagCallInfo {
  MagCallInfo() = default;
  MagCallInfo(TVirtualMC *mc, float x, float y, float z, float Bx, float By, float Bz);
  
  long id;
  long stepid; // cross-reference to current MC stepid (if any??)
  float x;
  float y;
  float z;
  float Bx;
  float By;
  float Bz;

  static int stepcounter;
  ClassDefNV(MagCallInfo, 1);
};
 
}
#endif
