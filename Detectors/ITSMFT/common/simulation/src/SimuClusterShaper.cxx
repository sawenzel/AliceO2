// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SimuClusterShaper.cxx
/// \brief Cluster shaper for the ALPIDE response simulation

#include <iostream>
#include <map>
#include <TBits.h>
#include <TRandom.h>

#include "ITSMFTSimulation/SimuClusterShaper.h"

using namespace o2::ITSMFT;

ClassImp(o2::ITSMFT::SimuClusterShaper);

//______________________________________________________________________
SimuClusterShaper::SimuClusterShaper() :
mHitX(0.f),
mHitZ(0.f),
mHitC(0),
mHitR(0),
mFireCenter(false),
mNpixOn(0),
mSeg(nullptr),
mCShape(nullptr) {}


//______________________________________________________________________
SimuClusterShaper::SimuClusterShaper(const UInt_t &cs) {
  mHitX = 0.f;
  mHitZ = 0.f;
  mHitC = 0;
  mHitR = 0;
  mFireCenter = false;
  mNpixOn = cs;
  UInt_t nRows = cs;
  UInt_t nCols = cs;

  mSeg = nullptr;
  mCShape = new ClusterShape(nRows, nCols);
}


//______________________________________________________________________
SimuClusterShaper::~SimuClusterShaper() {
  delete mCShape;
}


//______________________________________________________________________
void SimuClusterShaper::fillClusterRandomly() {
  Int_t matrixSize = mCShape->getNRows()*mCShape->getNCols();

  // generate UNIQUE random numbers
  UInt_t i = 0, j = 0;
  auto *bits = new TBits(mNpixOn);

  if (mFireCenter) {
    bits->SetBitNumber(mCShape->getCenterIndex());
    i++;
  }
  while (i < mNpixOn) {
    j = gRandom->Integer(matrixSize); // [0, matrixSize-1]
    if (bits->TestBitNumber(j)) continue;
    bits->SetBitNumber(j);
    i++;
  }

  Int_t bit = 0;
  for (i = 0; i < mNpixOn; ++i) {
    j = bits->FirstSetBit(bit);
    mCShape->addShapeValue(j);
    bit = j+1;
  }
  delete bits;
}


//______________________________________________________________________
void SimuClusterShaper::fillClusterSorted() {
  UInt_t matrixSize = mCShape->getNRows()*mCShape->getNCols();
  if (matrixSize == 1) {
    mCShape->addShapeValue(mCShape->getCenterIndex());
    return;
  }

  reComputeCenters();

  std::map<Double_t, UInt_t> sortedpix;
  Float_t pX = 0.f, pZ = 0.f;

  for (UInt_t i = 0; i < matrixSize; ++i) {
    UInt_t r = i / mCShape->getNRows();
    UInt_t c = i % mCShape->getNRows();
    UInt_t nx = mHitC - mCShape->getCenterC() + c;
    UInt_t nz = mHitR - mCShape->getCenterR() + r;
    mSeg->detectorToLocal(nx, nz, pX, pZ);
    Double_t d = sqrt(pow(mHitX-pX,2)+pow(mHitZ-pZ,2));

    // what to do when you reached the border?
    if (d > 1.) continue;
    sortedpix[d] = i;
  }

  // border case
  if (sortedpix.size() < mNpixOn) mNpixOn = sortedpix.size();
  for (std::map<Double_t, UInt_t>::iterator it = sortedpix.begin(); it != std::next(sortedpix.begin(),mNpixOn); ++it) {
    // std::cout << "  " << it->second << std::endl;
    mCShape->addShapeValue(it->second);
  }
}


//______________________________________________________________________
void SimuClusterShaper::addNoisePixel() {
  Int_t matrixSize = mCShape->getNRows()*mCShape->getNCols();
  UInt_t j = gRandom->Integer(matrixSize); // [0, matrixSize-1]
  while (mCShape->hasElement(j)) {
    j = gRandom->Integer(matrixSize);
  }
  //fCShape->SetShapeValue(i, j);
}


//______________________________________________________________________
void SimuClusterShaper::reComputeCenters() {
  UInt_t  r  = 0,   c = 0;
  Float_t pX = 0.f, pZ = 0.f;
  mSeg->detectorToLocal(mHitC, mHitR, pX, pZ);

  // c is even
  if (mCShape->getNCols() % 2 == 0) {
    if (mHitX > pX) { // n/2 - 1
      c = mCShape->getNCols()/2 - 1;
    } else { // n/2
      c = mCShape->getNCols()/2;
    }
  } else { // c is odd
    c = (mCShape->getNCols()-1)/2;
  }

  // r is even
  if (mCShape->getNRows() % 2 == 0) {
    if (mHitZ > pZ) { // n/2 - 1
      r = mCShape->getNRows()/2 - 1;
    } else { // n/2
      r = mCShape->getNRows()/2;
    }
  } else { // r is odd
    r = (mCShape->getNRows()-1)/2;
  }

  mCShape->setCenter(r, c);
}
