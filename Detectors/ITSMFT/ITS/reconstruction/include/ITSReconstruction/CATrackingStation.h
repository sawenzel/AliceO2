// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE ITSU Project       *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIITSUCATRACKINGSTATION_H
#define ALIITSUCATRACKINGSTATION_H

#include <algorithm>
#include <vector>
#include <cassert>
#include "ITSReconstruction/CAaux.h"

class TClonesArray;

namespace o2 {
  namespace ITS {
    class GeometryTGeo;
    namespace CA {
      class TrackingStation {
        public:
          //
          struct ClBinInfo { // info on bin clusters start, number of clusters
            int ncl;    // number of clusters
            int first;  // entry of 1st cluster in sorted vector of ClsInfo
            int index;  // index in the vector containing cells with non-0 occupancy
          };
          typedef struct ClBinInfo ClBinInfo_t;
          //
          TrackingStation();
          TrackingStation(int id,float zMin, float zMax, int nzbins,int nphibins);

          virtual ~TrackingStation();
          const ClsInfo_t& operator[](int i)    const {return mSortedClInfo[i];}
          //
          int     getVIDOffset()                const {return mVIDOffset;}
          int     getNClusters()                const {return mNClusters;}
          int     getNZBins()                   const {return mNZBins;}
          int     getNPhiBins()                 const {return mNPhiBins;}
          float   getZMin()                     const {return mZMin;}
          float   getZMax()                     const {return mZMax;}
          //
          void    setNZBins(int v)                    {mNZBins = v;}
          void    setNPhiBins(int v)                  {mNPhiBins = v;}
          void    setZMin(float v)                    {mZMin = v;}
          void    setZMax(float v)                    {mZMax = v;}
          //
          void init(TClonesArray* points, o2::ITS::GeometryTGeo* geom);
          //
          void sortClusters(const float vertex[3]);
          int  getPhiBin(float phi)             const {return phi * mDPhiInv;}
          int  getZBin  (float z)               const {return (z - mZMin) * mDZInv;}
          int  getBinIndex(int iz, int iphi)    const {return iphi * mNZBins + iz;}
          int  getBinZ(int ipz)                 const {return ipz % mNZBins;}
          int  getBinPhi(int ipz)               const {return ipz / mNZBins;}
          void getBinZPhi(int ipz,int &iz,int &iphi) const {iz = getBinZ(ipz); iphi=getBinPhi(ipz);}
          //
          int  selectClusters(float zmin,float zmax,float phimin,float phimax);
          int  getNFoundBins()                  const {return mFoundBins.size();}
          int  getFoundBin(int i)               const {return mFoundBins[i];}
          int  getFoundBinClusters(int i, int &first)  const;
          void resetFoundIterator();
          const ClsInfo_t& getClusterInfo(int i) const {return mSortedClInfo[i];}
          ClsInfo_t* getNextClusterInfo();
          int                     getNextClusterInfoID();
          //
          ITSDetInfo_t& getDetInfo(int id)     const
          {assert(mIndex[id] > -1 && "Empty sensor");return (ITSDetInfo_t&)mDetectors[mIndex[id]];}
          int                   getNDetectors()           const
          {return mDetectors.size();}
          //
          void         clearSortedInfo();
          void clear();
          //virtual void Print(Option_t *opt="")  const;

        protected:
          TrackingStation(const TrackingStation &);
          TrackingStation &operator=(const TrackingStation &);
          int   mID;                  // id of the layer
          int   mVIDOffset;           // offset of VID for detectors of this layer
          int   mNClusters;           // N clusters
          //
          float mZMin;                // Zmin
          float mZMax;                // Zmax
          float mDZInv;               // inverse size of Z bin
          float mDPhiInv;             // inverse size of Phi bin
          int   mNZBins;             // N cells in Z
          int   mNPhiBins;           // N cells in Phi
          //
          int   mQueryZBmin;         // min bin in Z of the query
          int   mQueryZBmax;         // max bin in Z of the query
          int   mQueryPhiBmin;       // min bin in phi of the query
          int   mQueryPhiBmax;       // max bin in phi of the query
          ClBinInfo_t* mBins;           // 2D (z,phi) grid of clusters binned in z,phi
          int* mOccBins;              // id's of bins with non-0 occupancy
          int  mNOccBins;             // number of occupied bins
          int  mNFoundClusters;       // number of found clusters in the query zone
          int  mFoundClusterIterator; // at which cluster within the bin we are?
          int  mFoundBinIterator;     // at which foune bin we are?
          std::vector<int>     mIndex;
          std::vector<int>     mFoundBins;    // occupied bins satisfying to query
          std::vector<ClsInfo_t> mSortedClInfo; // processed cluster info
          std::vector<ITSDetInfo_t> mDetectors; // detector params
          //
      };

      inline int TrackingStation::getFoundBinClusters(int i, int &first)  const {
        // set the entry of the first cl.info in the mSortedClInfo
        // and return n clusters in the bin
        ClBinInfo_t& bin = mBins[getFoundBin(i)];
        first = bin.first;
        return bin.ncl;
      }

      inline ClsInfo_t* TrackingStation::getNextClusterInfo() {
        // return cluster info for next matching cluster
        int id = getNextClusterInfoID();
        return id < 0 ? nullptr : (ClsInfo_t*)&mSortedClInfo[id];
      }
    } // namespace CA
  } // namespace ITS
} // namespace AliceO2
#endif

