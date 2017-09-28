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

#ifndef ALIITSUCATRACKER_H
#define ALIITSUCATRACKER_H

#include <vector>
#include <array>

#include "ITSReconstruction/CAaux.h"
#include "ITSReconstruction/CATrackingStation.h"
#include "DetectorsBase/Track.h"

namespace o2 {
  namespace ITS {
    namespace CA {
      typedef o2::Base::Track::TrackParCov TrackPC;

      class Tracker {
        public:
          Tracker(TrackingStation *stations[7]);
          // These functions must be implemented
          int clusters2Tracks();
          int  propagateBack();
          int  refitInward();
          int  loadClusters();
          void unloadClusters();
          // Possibly, other public functions
          float    getMaterialBudget(const double* p0, const double* p1, double& x2x0, double& rhol) const;
          bool     GetSAonly() const { return mSAonly; }
          void     setChi2Cut(float cut) { mChi2Cut = cut; }
          void     setPhiCut(float cut) { mPhiCut = cut; }
          void     setSAonly(bool sa = true) { mSAonly = sa; }
          void     setZCut(float cut) { mZCut = cut; }
          //
          float    getX() const { return mVertex[0]; }
          float    getY() const { return mVertex[1]; }
          float    getZ() const { return mVertex[2]; }
          template<typename F> void setVertex(F v[3]) { for(int i=0;i<3;++i) mVertex[i]=v[i]; }
        private:
          Tracker(const Tracker&);
          Tracker &operator=(const Tracker &tr);
          //
          bool   CellParams(int l, const Cluster& c1, const Cluster& c2, const Cluster& c3, float &curv, std::array<float,3> &np);
          void   cellsTreeTraversal(std::vector<Road> &roads, const int &iD, const int &doubl);
          void   findTracksCA(int iteration);
          void   makeCells(int iteration);
          bool   RefitAt(float xx, Track* t);
          void   setCuts(int it);
          void   setLabel(Track &t, float wrong);
          //
          TrackingStation**     mLayer;
          std::vector<bool>          mUsedClusters[7];
          float                 mChi2Cut;
          float                 mPhiCut;
          float                 mZCut;
          std::vector<Doublets>      mDoublets[6];
          std::vector<Cell>          mCells[5];
          std::vector<Track>         mCandidates[4];
          bool                  mSAonly;             // true if the standalone tracking only
          // Cuts
          float mCPhi;
          float mCDTanL;
          float mCDPhi;
          float mCZ;
          float mCDCAz[5];
          float mCDCAxy[5];
          float mCDN[4];
          float mCDP[4];
          float mCDZ[6];
          //
          float mVertex[3];
          float mBz;
          //
          static const float              mkChi2Cut;      // chi2 cut during track merging
          static const int                mkNumberOfIterations;
          static const float              mkR[7];
          //
      };
    } // namespace CA
  } // namespace ITS
} // namespace AliceO2

#endif // ALIITSUCATRACKER_H
