// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TOF_DETECTOR_H_
#define ALICEO2_TOF_DETECTOR_H_

#include "DetectorsBase/Detector.h"
#include "TOFBase/Geo.h"

#include "SimulationDataFormat/BaseHits.h"
#include "DetectorsBase/VMCUtilities.h"
#include "FairVolume.h"

namespace o2
{
namespace tof
{
using HitType = o2::BasicXYZEHit<float>;

class Detector : public o2::Base::DetImpl<Detector>
{
 public:
  enum TOFMaterial {
    kAir = 1,
    kNomex = 2,
    kG10 = 3,
    kFiberGlass = 4,
    kAlFrame = 5,
    kHoneycomb = 6,
    kFre = 7,
    kCuS = 8,
    kGlass = 9,
    kWater = 10,
    kCable = 11,
    kCableTubes = 12,
    kCopper = 13,
    kPlastic = 14,
    kCrates = 15,
    kHoneyHoles = 16
  };

  Detector() = default;

  Detector(Bool_t active);

  ~Detector() override = default;

  void Initialize() final;

  // Bool_t ProcessHits(FairVolume* v = nullptr) final;

  void Register() override;

  std::vector<HitType>* getHits(int iColl) const
  {
    if (iColl == 0) {
      return mHits;
    }
    return nullptr;
  }

  void Reset() final;
  void EndOfEvent() final;

  void CreateMaterials();
  void ConstructGeometry() final;
  void addAlignableVolumes() const override;

  void setTOFholes(Bool_t flag = kTRUE) { mTOFHoles = flag; }
 protected:
  virtual void DefineGeometry(Float_t xtof, Float_t ytof, Float_t zlenA) final;
  virtual void MaterialMixer(Float_t* p, const Float_t* const a, const Float_t* const m, Int_t n) const final;

 private:
  /// copy constructor (used in MT)
  Detector(const Detector& rhs);

  void createModules(Float_t xtof, Float_t ytof, Float_t zlenA, Float_t xFLT, Float_t yFLT, Float_t zFLTA) const;
  void makeStripsInModules(Float_t ytof, Float_t zlenA) const;
  void createModuleCovers(Float_t xtof, Float_t zlenA) const;
  void createBackZone(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void makeFrontEndElectronics(Float_t xtof) const;
  void makeFEACooling(Float_t xtof) const;
  void makeNinoMask(Float_t xtof) const;
  void makeSuperModuleCooling(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void makeSuperModuleServices(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void makeReadoutCrates(Float_t ytof) const;

  void makeModulesInBTOFvolumes(Float_t ytof, Float_t zlenA) const;
  void makeCoversInBTOFvolumes() const;
  void makeBackInBTOFvolumes(Float_t ytof) const;

  bool isMergable(HitType hit1, HitType hit2)
  {
    if (hit1.GetTrackID() != hit2.GetTrackID()) {
      return false;
    }

    if (std::abs(hit1.GetTime() - hit2.GetTime()) > 1.0 /*1 ns*/) {
      return false;
    }

    return true;
  }

  template<typename E>
  Bool_t ProcessHitsKernel(FairVolume* v)
  {
    using namespace vmc_helpers;
    // This method is called from the MC stepping for the sensitive volume only
    if (static_cast<int>(TrackCharge<E>(fMC)) == 0) {
      // set a very large step size for neutral particles
      return kFALSE; // take only charged particles
    }

//    float pos2x, pos2y, pos2z;
 //   fMC->TrackPosition(pos2x, pos2y, pos2z);
  //  Float_t radius = std::sqrt(pos2x * pos2x + pos2y * pos2y);
  //  LOG(DEBUG) << "Process hit in TOF volume ar R=" << radius << " - Z=" << pos2z;

    Float_t enDep = Edep<E>(fMC);
    if (enDep < 1E-8) {
      return kFALSE; // wo se need a threshold?
    }

    // ADD HIT
    float posx, posy, posz;
    TrackPosition<E>(fMC, posx, posy, posz);
    float time = TrackTime<E>(fMC) * 1.0e09;
    int trackID = mO2Stack->GetCurrentTrackNumber();
    int sensID = v->getMCid();
    Int_t det[5];
    Float_t pos[3] = { posx, posy, posz };
    Float_t delta[3];
    Geo::getPadDxDyDz(pos, det, delta);
    auto channel = Geo::getIndex(det);
    HitType newhit(posx, posy, posz, time, enDep, trackID, sensID);
    if (channel != mLastChannelID || !isMergable(newhit, mHits->back())) {
      mHits->push_back(newhit);
    } else {
      mHits->back().SetEnergyLoss(mHits->back().GetEnergyLoss() + newhit.GetEnergyLoss());
      // LOG(INFO)<<"Merging hit "<<"\n";
      //  <<mHits->back().GetId()<<"with new hit "<<newhit.GetId()<<"\n";
    }
    mLastChannelID = channel;
    return kTRUE;
  }


  Int_t mEventNr; // event count
  Int_t mTOFSectors[o2::tof::Geo::NSECTORS];
  Bool_t mTOFHoles; // flag to allow for holes in front of the PHOS
  Int_t mLastChannelID = -1;///< Last channel seen by the hit

  /// container for data points
  std::vector<HitType>* mHits; //!

  template <typename Det>
  friend class o2::Base::DetImpl;
  ClassDefOverride(Detector, 1);
};
}
}
#endif
