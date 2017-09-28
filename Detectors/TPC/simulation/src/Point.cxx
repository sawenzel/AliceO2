// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#include "TPCSimulation/Point.h"
#include <iostream>

using std::cout;
using std::endl;

using namespace o2::TPC;

void Point::print(const Option_t* opt) const
{
  cout << "-I- Point: O2tpc point for track " << getTrackID()
       << " in detector " << getDetectorID() << endl;
  cout << "    Position (" << getX() << ", " << getY() << ", " << getZ()
       << ") cm" << endl;
  cout << "    Time " << getTime() << " ns, n electrons " << getEnergyLoss() << endl;
}

ClassImp(Point)
ClassImp(LinkableHitGroup)
ClassImp(ElementalHit)
