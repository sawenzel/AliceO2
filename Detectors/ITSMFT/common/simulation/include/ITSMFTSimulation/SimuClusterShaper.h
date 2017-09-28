// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SimuClusterShaper.h
/// \brief Cluster shaper for the ALPIDE response simulation

#ifndef ALICEO2_ITSMFT_SIMUCLUSTERSHAPER_H_
#define ALICEO2_ITSMFT_SIMUCLUSTERSHAPER_H_

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to generate the cluster shape in the ITSU simulation    //
// Author: Davide Pagano                                         //
///////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <sstream>

#include "ITSMFTSimulation/ClusterShape.h"

namespace o2 {
  namespace ITSMFT {

    class SimuClusterShaper : public TObject {

    public:
      SimuClusterShaper();
      SimuClusterShaper(const UInt_t &cs);
      ~SimuClusterShaper() override;
      void fillClusterRandomly();
      void fillClusterSorted();
      inline void setFireCenter(Bool_t v) {
        mFireCenter = v;
      }
      void addNoisePixel();

      inline void    setHit(Int_t ix, Int_t iz, Float_t x, Float_t z, const SegmentationPixel* seg) {
        mHitC = ix;
        mHitR = iz;
        mHitX = x;
        mHitZ = z;
        mSeg  = seg;
      }
      inline UInt_t  getNRows() {return mCShape->getNRows();}
      inline UInt_t  getNCols() {return mCShape->getNCols();}
      inline void    getShape(std::vector<UInt_t>& v) {mCShape->getShape(v);}
      inline UInt_t  getCenterR() {return mCShape->getCenterR();}
      inline UInt_t  getCenterC() {return mCShape->getCenterC();}
      inline size_t  getCS() {return mCShape->getNFiredPixels();}

      inline std::string shapeSting(UInt_t cs, UInt_t *cshape) const {
        std::stringstream out;
        for (Int_t i = 0; i < cs; ++i) {
          out << cshape[i];
          if (i < cs-1) out << " ";
        }
        return out.str();
      }

    private:
      void reComputeCenters();

      Float_t mHitX;
      Float_t mHitZ;
      Int_t   mHitC;
      Int_t   mHitR;
      Bool_t  mFireCenter;
      UInt_t  mNpixOn;
      const SegmentationPixel* mSeg;
      ClusterShape *mCShape;

      ClassDefOverride(SimuClusterShaper,1);

    };
  }
}
#endif
