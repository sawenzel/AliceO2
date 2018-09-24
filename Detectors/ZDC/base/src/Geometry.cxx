// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ZDCBase/Geometry.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "FairLogger.h"

ClassImp(o2::zdc::Geometry);

using namespace o2::zdc;


void Geometry::Init()
{
  Info("zdc::Geometry::Init", "Initialization of ZDC rotation parameters");
  
  Double_t rotAnglesPipe1[6] =  {90., -1.0027, 0., 90., 90., 1.0027, 180.};
  Double_t rotAnglesPipe2[6] =  {90.,  1.0027, 0., 90., 90., 1.0027, 0.};

}
