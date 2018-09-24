// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_ZDC_DETECTOR_H_
#define ALICEO2_ZDC_DETECTOR_H_

#include <vector>                             // for vector
#include "Rtypes.h"                           // for Int_t, Double_t, Float_t, Bool_t, etc
#include "TGeoManager.h"                      // for gGeoManager, TGeoManager (ptr only)
#include "DetectorsBase/GeometryManager.h"    // for getSensID
#include "DetectorsBase/Detector.h"           // for Detector
#include "DetectorsCommonDataFormats/DetID.h" // for Detector
#include "ZDCBase/Geometry.h"

#include "SimulationDataFormat/BaseHits.h"

class FairVolume;
class FairModule;

class TParticle;

namespace o2
{
namespace zdc
{
using Hit = o2::BasicXYZEHit<float>;

class Detector : public o2::Base::DetImpl<Detector>
{
 public:
  enum ZDCMaterial {
    kWalloy = 1,
    kCuZn = 2,
    kSiO2pmc = 3,
    kSiO2pmq = 4,
    kPb = 5,
    kCu = 6,
    kFe = 7,
    kAl = 8,
    kGraphite = 9,
    kVoidNoField = 10,
    kVoidwField = 11,
    kAir = 12
  };

  Detector() = default;

  Detector(Bool_t active);

  ~Detector() override = default;

  FairModule* CloneModule() const override;

  void Initialize() final;

  Bool_t ProcessHits(FairVolume* v = nullptr) final;

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


  void ConstructGeometry() final;
  void CreateMaterials();
  void addAlignableVolumes() const override;

  Hit* AddHit(Int_t trackID, Int_t trackPDG, Int_t parentID, Int_t sFlag, Double_t primaryEnergy, Int_t& detID,
              Double_t& pos, Double_t& mom, Double_t tof, Double_t& xImpact, Double_t energyloss, Int_t nphe);

   private:
    /// copy constructor
    Detector(const Detector& rhs);
    void CreateAsideBeamLine();
    void CreateCsideBeamLine();
    void CreateSupports();
    void CreateMagnets();
    void CreateDetectors();
    /// Define sensitive volumes
    void defineSensitiveVolumes();

    bool isMergeable(Hit hit1, Hit hit2)
    {
      if (hit1.GetTrackID() != hit2.GetTrackID()) {
        return false;
      }

      return true;
    }

    Int_t mZDCdetID[2]; //detector|tower in ZDC
    Int_t mPcMother; // track mother 0
    Int_t mSecondaryFlag;
    Float_t mPrimaryEnergy;
    Float_t mXImpact[3];
    Float_t mTrackTOF;
    Float_t mTotDepEnergy;
    Float_t mTotLight[2]; //[0]PMC [1]sumPMQi

    /// container for data points
    std::vector<Hit>* mHits; //!

  ClassDefOverride(Detector, 1);
};
}
}
#endif
