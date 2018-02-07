// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Digitizer.cxx
/// \brief Implementation of the ALICE TPC digitizer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "TPCSimulation/Digitizer.h"
#include "TPCSimulation/ElectronTransport.h"
#include "TPCSimulation/GEMAmplification.h"
#include "TPCSimulation/PadResponse.h"
#include "TPCSimulation/Point.h"
#include "TPCSimulation/SAMPAProcessing.h"
#include "TPCBase/ParameterDetector.h"
#include "TPCBase/ParameterElectronics.h"
#include "TPCBase/ParameterGas.h"

#include "TPCBase/Mapper.h"

#include "FairLogger.h"

ClassImp (o2::TPC::Digitizer)

using namespace o2::TPC;
bool o2::TPC::Digitizer::mIsContinuous = true;

Digitizer::Digitizer()
    : mDigitContainer(nullptr)
{
}

Digitizer::~Digitizer()
{
  delete mDigitContainer;
}

void Digitizer::init()
{
  /// Initialize the task and the output container
  /// \todo get rid of new? check with Mohammad
  mDigitContainer = new DigitContainer();

//  mDebugTreePRF = std::unique_ptr<TTree> (new TTree("PRFdebug", "PRFdebug"));
//  mDebugTreePRF->Branch("GEMresponse", &GEMresponse, "CRU:timeBin:row:pad:nElectrons");
}

template<typename T> int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

DigitContainer* Digitizer::Process(const Sector &sector, const std::vector<o2::TPC::HitGroup>& hits, int eventID,
    float eventTime)
{
//  mDigitContainer->reset();

  if (!mIsContinuous)
    eventTime = 0.f; /// transform in us

  /// \todo static_thread for thread savety?

  static size_t hitCounter = 0;
  for (auto& inputgroup : hits) {
    ProcessHitGroup(inputgroup, sector, eventTime, eventID);
  }
  /// end of loop over points

  return mDigitContainer;
}

DigitContainer* Digitizer::ProcessNEW(const Sector &sector, const std::vector<std::vector<o2::TPC::HitGroup>*> &hits,
    const std::vector<o2::TPC::TPCHitGroupID>& hitids, const o2::steer::RunContext& context)
{
  const auto& interactRecords = context.getEventRecords();

  for (auto& id : hitids) {
    auto entryvector = hits[id.entry];
    auto& group = (*entryvector)[id.groupID];
    auto& MCrecord = interactRecords[id.entry];
    ProcessHitGroup(group, sector, MCrecord.timeNS * 0.001f, id.entry);
  }

  return mDigitContainer;
}

void Digitizer::ProcessHitGroup(const HitGroup &inputgroup, const Sector &sector, const float eventTime, const int eventID)
{

  const static Mapper& mapper = Mapper::instance();
  const static ParameterDetector &detParam = ParameterDetector::defaultInstance();
  const static ParameterElectronics &eleParam = ParameterElectronics::defaultInstance();
  static GEMAmplification gemAmplification;
  static ElectronTransport electronTransport;
  static PadResponse padResponse;

  const int nShapedPoints = eleParam.getNShapedPoints();
  static std::vector<float> signalArray;
  signalArray.resize(nShapedPoints);

  const int MCTrackID = inputgroup.GetTrackID();
  for (size_t hitindex = 0; hitindex < inputgroup.getSize(); ++hitindex) {
    const auto& eh = inputgroup.getHit(hitindex);

    const GlobalPosition3D posEle(eh.GetX(), eh.GetY(), eh.GetZ());
    std::cout << Sector::ToSector(eh.GetX(), eh.GetY(), eh.GetZ()) << "\n";
    if (electronTransport.isCompletelyOutOfSectorCourseElectronDrift(posEle, sector)) {
      std::cout << "not processing hit ... " << Sector::ToSector(eh.GetX(), eh.GetY(), eh.GetZ()) << "\n";
      continue;
    }
    std::cout << "PROCESSING hit ... " << Sector::ToSector(eh.GetX(), eh.GetY(), eh.GetZ()) << "\n";


    // The energy loss stored is really nElectrons
    const int nPrimaryElectrons = static_cast<int>(eh.GetEnergyLoss());

    /// Loop over electrons
    /// \todo can be vectorized?
    /// \todo split transport and signal formation in two separate loops?
    for (int iEle = 0; iEle < nPrimaryElectrons; ++iEle) {

      /// Drift and Diffusion
      GlobalPosition3D posEleDiff = electronTransport.getElectronDrift(posEle);
      if (sgn(posEleDiff.Z()) != sgn(posEle.Z())) {
        posEleDiff.SetZ(posEle.Z());
      }

      /// \todo Time management in continuous mode (adding the time of the event?)
      const float driftTime = SAMPAProcessing::getDriftTime(posEleDiff.Z()) + eh.GetTime() * 0.001; /// in us
      const float absoluteTime = driftTime + eventTime;

      /// Attachment
      if (electronTransport.isElectronAttachment(driftTime))
        continue;

      /// Remove electrons that end up outside the active volume
      /// \todo should go to mapper?
      if (std::abs(posEleDiff.Z()) > detParam.getTPClength())
        continue;

      const DigitPos digiPadPos = mapper.findDigitPosFromGlobalPosition(posEleDiff);
      if (!digiPadPos.isValid())
        continue;

      // check that the Digit is in the sector of interest
      if(digiPadPos.getCRU().sector() != sector) continue;

      const int nElectronsGEM = gemAmplification.getStackAmplification();
      if (nElectronsGEM == 0)
        continue;

      /// Loop over all individual pads with signal due to pad response function
      /// Currently the PRF is not applied yet due to some problems with the mapper
      /// which results in most of the cases in a normalized pad response = 0
      /// \todo Problems of the mapper to be fixed
      /// \todo Mapper should provide a functionality which finds the adjacent pads of a given pad
      // for(int ipad = -2; ipad<3; ++ipad) {
      //   for(int irow = -2; irow<3; ++irow) {
      //     PadPos padPos(digiPadPos.getPadPos().getRow() + irow, digiPadPos.getPadPos().getPad() + ipad);
      //     DigitPos digiPos(digiPadPos.getCRU(), padPos);

      DigitPos digiPos = digiPadPos;
      if (!digiPos.isValid())
        continue;
      // const float normalizedPadResponse = padResponse.getPadResponse(posEleDiff, digiPos);

      const float normalizedPadResponse = 1.f;
      if (normalizedPadResponse <= 0)
        continue;
      const int pad = digiPos.getPadPos().getPad();
      const int row = digiPos.getPadPos().getRow();
      const GlobalPadNumber globalPad = mapper.getPadNumberInROC(
          PadROCPos(digiPadPos.getCRU().roc(), PadPos(row, pad)));

      const float ADCsignal = SAMPAProcessing::getADCvalue(nElectronsGEM * normalizedPadResponse);
      SAMPAProcessing::getShapedSignal(ADCsignal, absoluteTime, signalArray);
      for (float i = 0; i < nShapedPoints; ++i) {
        const float signal = signalArray[i];
        if (signal == 0)
          continue;
        const float time = absoluteTime + i * eleParam.getZBinWidth();
        const auto digisector = digiPos.getCRU().sector();
        if (digisector == sector || digisector == Sector::getLeft(sector) || digisector == Sector::getRight(sector)) {
          mDigitContainer->addDigit(eventID, MCTrackID, digiPos.getCRU(), SAMPAProcessing::getTimeBinFromTime(time),
              globalPad, signal);
        } else {
          LOG(INFO) << "WARNING: Digit in unexpected sector " << digisector << " " << posEle.Z() << "\t"
              << posEleDiff.Z() << "\n";
        }
      }

      // }
      // }
      /// end of loop over prf
    }
    /// end of loop over electrons
  }
}
