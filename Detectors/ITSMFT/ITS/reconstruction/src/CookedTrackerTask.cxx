// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CookedTrackerTask.cxx
/// \brief Implementation of the ITS "Cooked Matrix" tracker task
/// \author iouri.belikov@cern.ch

#include "ITSReconstruction/CookedTrackerTask.h"
#include "ITSMFTReconstruction/Cluster.h"
#include "DetectorsBase/Utils.h"
#include "MathUtils/Cartesian3D.h"

#include "FairLogger.h"      // for LOG
#include "FairRootManager.h" // for FairRootManager
#include "TClonesArray.h"    // for TClonesArray

ClassImp(o2::ITS::CookedTrackerTask)

using namespace o2::ITS;
using namespace o2::Base;
using namespace o2::Base::Utils;

//_____________________________________________________________________
CookedTrackerTask::CookedTrackerTask(Int_t n) : FairTask("ITSCookedTrackerTask"), mNumOfThreads(n), mClustersArray(nullptr), mTracksArray(nullptr) {}

//_____________________________________________________________________
CookedTrackerTask::~CookedTrackerTask()
{
  if (mTracksArray) {
    mTracksArray->Delete();
    delete mTracksArray;
  }
}

//_____________________________________________________________________
/// \brief Init function
/// Inititializes the tracker and connects input and output container
InitStatus CookedTrackerTask::init()
{
  FairRootManager* mgr = FairRootManager::Instance();
  if (!mgr) {
    LOG(ERROR) << "Could not instantiate FairRootManager. Exiting ..." << FairLogger::endl;
    return kERROR;
  }

  mClustersArray = dynamic_cast<const TClonesArray*>(mgr->GetObject("ITSCluster"));
  if (!mClustersArray) {
    LOG(ERROR) << "ITS clusters not registered in the FairRootManager. Exiting ..." << FairLogger::endl;
    return kERROR;
  }

  // Register output container
  mTracksArray = new TClonesArray("o2::ITS::CookedTrack");
  mgr->Register("ITSTrack", "ITS", mTracksArray, kTRUE);

  GeometryTGeo* geom = GeometryTGeo::Instance(instance;
  geom->fillMatrixCache( bit2Mask(TransformType::T2GRot) ); // make sure T2GRot matrices are loaded
  mTracker.setGeometry(geom);

  mTracker.setNumberOfThreads(mNumOfThreads);
  
  return kSUCCESS;
}

//_____________________________________________________________________
void CookedTrackerTask::exec(Option_t* option)
{
  mTracksArray->Clear();
  LOG(DEBUG) << "Running digitization on new event" << FairLogger::endl;

  mTracker.process(*mClustersArray, *mTracksArray);
}
