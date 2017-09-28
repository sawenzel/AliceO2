// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "EMCALBase/Hit.h"

ClassImp(o2::EMCAL::Hit)

using namespace o2::EMCAL;

void Hit::printStream(std::ostream &stream) const {
  stream  << "EMCAL point: Track " << getTrackID() << " in detector segment " << getDetectorID()
          << " at position (" << getX() << "|" << getY() << "|" << getZ() << "), energy loss " << getEnergyLoss()
          << ", parent " << mParent << " with energy " << mInitialEnergy;
}

Bool_t Hit::operator<(const Hit &rhs) const {
  if(mParent != rhs.mParent) return mParent < rhs.mParent;
  return getDetectorID() < rhs.getDetectorID();
}

Bool_t Hit::operator==(const Hit &rhs) const {
  return (getDetectorID() == getDetectorID()) && (mParent == rhs.mParent);
}

Hit &Hit::operator+=(const Hit &rhs) {
  setEnergyLoss(getEnergyLoss() + rhs.getEnergyLoss());
  return *this;
}

Hit Hit::operator+(const Hit &rhs) const {
  Hit result(*this);
  result.setEnergyLoss(result.getEnergyLoss() + rhs.getEnergyLoss());
  return *this;
}

std::ostream &operator<<(std::ostream &stream, const Hit &p) {
  p.printStream(stream);
  return stream;
}
