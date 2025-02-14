// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define _USE_MATH_DEFINES

#include <cmath>
#include <memory>

// root includes
#include "TFile.h"
#include <TH1.h>

// o2 includes
#include "TPCQC/GPUErrorQA.h"
#include "GPUErrors.h"

ClassImp(o2::tpc::qc::GPUErrorQA);

using namespace o2::tpc::qc;

//______________________________________________________________________________
void GPUErrorQA::initializeHistograms()
{
  TH1::AddDirectory(false);
  mHist = std::make_unique<TH1F>("ErrorCounter", "ErrorCounter", o2::gpu::GPUErrors::getMaxErrors(), 0, o2::gpu::GPUErrors::getMaxErrors());
}
//______________________________________________________________________________
void GPUErrorQA::resetHistograms()
{
  mHist->Reset();
}
//______________________________________________________________________________
void GPUErrorQA::processErrors(gsl::span<const std::array<uint32_t, 4>> errors)
{
  for (const auto& error : errors) {
    uint32_t errorCode = error[0];
    mHist->Fill(static_cast<float>(errorCode));
  }
}

//______________________________________________________________________________
void GPUErrorQA::dumpToFile(const std::string filename)
{
  auto f = std::unique_ptr<TFile>(TFile::Open(filename.c_str(), "recreate"));
  mHist->Write();
  f->Close();
}
