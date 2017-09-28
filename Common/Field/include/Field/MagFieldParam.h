// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MagFieldParam.h
/// \brief Definition of the MagFieldParam: container for ALICE mag. field parameters
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_FIELD_MAGFIELDPARAM_H_
#define ALICEO2_FIELD_MAGFIELDPARAM_H_

#include "FairParGenericSet.h"
#include <TString.h>

class FairParamList;


namespace o2 {
  namespace field {

class MagneticField;
    
class MagFieldParam : public FairParGenericSet
{
  public:
    enum BMap_t
    {
      k2kG, k5kG, k5kGUniform, kNFieldTypes
    };
    enum BeamType_t
    {
        kNoBeamField, kBeamTypepp, kBeamTypeAA, kBeamTypepA, kBeamTypeAp
    };

    MagFieldParam(const char* name="", const char* title="", const char* context="");

    void setParam(const MagneticField* field);
    
    BMap_t     getMapType()                   const {return mMapType;}
    BeamType_t getBeamType()                  const {return mBeamType;}
    Int_t                     getDefInt()     const {return mDefaultIntegration;}
    Double_t                  getFactorSol()  const {return mFactorSol;}
    Double_t                  getFactorDip()  const {return mFactorDip;}
    Double_t                  getBeamEnergy() const {return mBeamEnergy;}
    Double_t                  getMaxField()   const {return mMaxField;}
    const char*               getMapPath()    const {return mMapPath.Data();}

    void   putParams(FairParamList* list) override;
    Bool_t getParams(FairParamList* list) override;
    
  protected:
    BMap_t     mMapType;  ///< map type ID
    BeamType_t mBeamType; ///< beam type ID
    Int_t    mDefaultIntegration;        ///< field integration type for MC
    Double_t mFactorSol;                 ///< solenoid current factor
    Double_t mFactorDip;                 ///< dipole current factor
    Double_t mBeamEnergy;                ///< beam energy
    Double_t mMaxField;                  ///< max field for geant
    TString  mMapPath;                   ///< path to map file
    
    ClassDefOverride(MagFieldParam,1)
};
}
}

#endif
