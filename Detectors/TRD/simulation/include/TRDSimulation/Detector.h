// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TRD_DETECTOR_H_
#define ALICEO2_TRD_DETECTOR_H_

#include <vector>
#include "DetectorsBase/Detector.h"
#include "SimulationDataFormat/BaseHits.h"
#include "CommonUtils/ShmAllocator.h"

class FairVolume;
class TClonesArray;

#ifdef USESHM
namespace std
{
template<> class allocator<o2::BasicXYZEHit<float>> : public o2::utils::ShmAllocator<o2::BasicXYZEHit<float>>
{
};
}
#endif

namespace o2
{
namespace trd
{
class TRDGeometry;

// define TRD hit type
using HitType = o2::BasicXYZEHit<float>;

class Detector : public o2::Base::DetImpl<Detector>
{
  friend class o2::Base::DetImpl<Detector>;

 public:
  Detector(Bool_t active=true);

  ~Detector() override;

  FairModule* CloneModule() const override;

  void Initialize() override;

  bool ProcessHits(FairVolume* v = nullptr) override;

  void Register() override;

  std::vector<HitType>* getHits(int iColl) const
  {
    if (iColl == 0) {
      return mHits;
    }
    return nullptr;
  }

 private:
   bool setHits(int i, std::vector<HitType>* ptr)
   {
     if (i == 0) {
       mHits = ptr;
       return false;
     }
     return false;
   }

   void createHitBuffers()
   {
     for (int buffer = 0; buffer < NHITBUFFERS; ++buffer) {
       int probe = 0;
       bool more{ false };
       do {
         auto ptr = o2::utils::createSimVector<HitType>();
         more = setHits(probe, ptr);
         mCachedPtr[buffer].emplace_back(ptr);
         probe++;
       } while (more);
     }
   }

  public:
   void Reset() override;
   void EndOfEvent() override;

   void createMaterials();
   void ConstructGeometry() override;

  private:
   /// copy constructor (used in MT)
   Detector(const Detector& rhs);

   // defines/sets-up the sensitive volumes
   void defineSensitiveVolumes();

   // addHit
   template <typename T>
   void addHit(T x, T y, T z, T time, T energy, int trackId, int detId);

   std::vector<HitType>* mHits = nullptr; ///!< Collection of TRD hits

   float mFoilDensity;
   float mGasNobleFraction;
   float mGasDensity;

   TRDGeometry* mGeom = nullptr;

   ClassDefOverride(Detector, 1)
};

template <typename T>
void Detector::addHit(T x, T y, T z, T time, T energy, int trackId, int detId)
{
  mHits->emplace_back(x, y, z, time, energy, trackId, detId);
}

} // end namespace trd
} // end global namespace

#ifdef USESHM
namespace o2
{
namespace Base
{
template <>
struct UseShm<o2::trd::Detector> {
  static constexpr bool value = true;
};
}
}
#endif
#endif
