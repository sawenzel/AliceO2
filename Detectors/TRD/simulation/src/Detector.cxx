// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TClonesArray.h"
#include "FairVolume.h"
#include "FairRootManager.h"
#include "TRDSimulation/Detector.h"
#include "TRDBase/TRDGeometry.h"

using namespace o2::trd;

Detector::Detector(const char* Name, Bool_t Active):
o2::Base::Detector(Name, Active),
  mHitCollection(new TClonesArray("o2::BasicXYZEHit<float>"))
{
}

void Detector::Initialize(){
}

Bool_t Detector::ProcessHits(FairVolume* v) {
  return true;
}

void Detector::Register(){
  FairRootManager::Instance()->Register("TRDHit", "TRD", mHitCollection, kTRUE);
}

TClonesArray* Detector::GetCollection(Int_t iColl) const {
  if(iColl == 0) return mHitCollection;
  return nullptr;
}

void Detector::Reset() {
}

// setting up the geometry
void Detector::ConstructGeometry() {
  std::cerr << "TRD geom called\n";
  TRDGeometry geom;
  geom.CreateGeometry(nullptr);
}


ClassImp(Detector);
