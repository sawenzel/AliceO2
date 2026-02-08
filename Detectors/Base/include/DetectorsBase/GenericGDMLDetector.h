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

/// \file GenericGDMLDetector.h
/// \brief Definition of the a generic GDML Detector class

#ifndef ALICEO2_BASE_GDMLDETECTOR_H_
#define ALICEO2_BASE_GDMLDETECTOR_H_

#include <DetectorsBase/Detector.h>
#include <TGDMLParse.h> // from ROOT
#include <SimulationDataFormat/BaseHits.h>
#include <iostream>
#include <TGeoManager.h>

namespace o2::base {

// options used to configure a generic plug and play detector from GDML
struct GDMLDetectorOptions {
  std::string gdml_file;
  std::string top_volume; // the volume to be added
  std::string anchor_volume; // the volume into which we want to add
  std::string anchor_transform;
  std::string sensitiveActionMacro;
};

// just a Helper function which can be compiled in a different translation unit
class GDMLHelper {

};

// template <typename HitType>
class GenericGDMLDetector : public DetImpl<GenericGDMLDetector> {
  public:
    using DetImpl<GenericGDMLDetector>::DetImpl;
    // GenericGDMLDetector();
    virtual ~GenericGDMLDetector() = default;
    using Hit_t = o2::BasicXYZEHit<float>;

    void setOptions(GDMLDetectorOptions const& opts) {
      mOptions = opts;
    }

    void EndOfEvent() final;
    void Register() final;
    void Reset() final;

    bool ProcessHits(FairVolume *) final;
    
    std::vector<Hit_t>* getHits(Int_t iColl) const
    {
      if (iColl == 0) {
        return mHits;
      }
      return nullptr;
    }

    void addSensitiveVolumes();

    // a special detector init hook called from the O2 framework -- not FairRoot
    void InitializeO2Detector() final {
      addSensitiveVolumes();
      mHits = o2::utils::createSimVector<Hit_t>();
    }

    // an interface from FairROOT; needs to be overwritten
    // should define the materials and media and geometry and register with all services
    void ConstructGeometry() final {
       std::cout << "Executing GDML construction\n";
       if (!gGeoManager) {
         LOG(error) << " No GeoManager present";
         return;
       }

       if (mOptions.gdml_file.size() == 0) {
         LOG(error) << " No GDML file given";
         return;
       }

       // now obtain the wanted volume from GDML
       auto gdml_world = TGDMLParse::StartGDML(mOptions.gdml_file.c_str());
       
       if (!gdml_world) {
         LOG(error) << " Importing GDML volumes unsuccessful";
       }

       // now we need to extract specific gdml_vol that we are interested in
       auto volume_to_add = gGeoManager->FindVolumeFast(mOptions.top_volume.c_str());
       if (!volume_to_add) {
         LOG(error) << "GDML volume " << mOptions.top_volume << " not found ";
       }

       // now we need to find the volume to anchor into; it should already be here
       auto anchor_volume = gGeoManager->FindVolumeFast(mOptions.anchor_volume.c_str());
       if (!anchor_volume) {
         LOG(error) << "Anchor volume " << mOptions.anchor_volume << " not found ";
       }

       //
       if (volume_to_add) {
         // TODO: treat
         anchor_volume->AddNode(volume_to_add, 1); 
       }
    
    
       // register the materials and media with the MaterialManager
       // apply scaling etc.
    }

  // protected:
    std::vector<Hit_t>* mHits;
    GDMLDetectorOptions mOptions;

    ClassDef(GenericGDMLDetector, 1);
};

} // end namespace

#ifdef USESHM
namespace o2
{
namespace base
{
template <>
struct UseShm<o2::base::GenericGDMLDetector> {
  static constexpr bool value = false;
};
} // namespace base
} // namespace o2
#endif


#endif