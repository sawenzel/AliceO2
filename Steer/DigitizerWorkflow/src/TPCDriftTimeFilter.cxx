// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCDriftTimeFilter.h"
#include <FairMQLogger.h>
#include <TMessage.h> // object serialization
#include <cassert>
#include <cstring> // memcpy
#include <memory>  // std::unique_ptr
#include <string>  // std::string
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Headers/DataHeader.h"
#include "Steer/HitProcessingManager.h"
#include "TPCSimulation/Digitizer.h"
#include <TPCSimulation/Digitizer.h>
#include <TPCSimulation/DigitizerTask.h>
#include <cassert>
#include <functional>
#include "ITSMFTSimulation/Hit.h"
#include "TPCSimulation/Point.h"
#include "TPCSimulation/SAMPAProcessing.h"
#include <sstream>
#include <algorithm>

using DataProcessorSpec = o2::framework::DataProcessorSpec;
using Inputs = o2::framework::Inputs;
using Outputs = o2::framework::Outputs;
using Options = o2::framework::Options;
using InputSpec = o2::framework::InputSpec;
using OutputSpec = o2::framework::OutputSpec;
using AlgorithmSpec = o2::framework::AlgorithmSpec;
using InitContext = o2::framework::InitContext;
using ProcessingContext = o2::framework::ProcessingContext;
using VariantType = o2::framework::VariantType;
using ControlService = o2::framework::ControlService;
using SubSpecificationType = o2::framework::DataAllocator::SubSpecificationType;
using DataRefUtils = o2::framework::DataRefUtils;

namespace o2
{
namespace steer
{

template <typename Collection>
void getHits(TChain& chain, const Collection& eventrecords, std::vector<std::vector<o2::TPC::HitGroup>*>& hitvectors,
             std::vector<o2::TPC::TPCHitGroupID>& hitids, const char* branchname, float tmin /*NS*/, float tmax /*NS*/,
             std::function<float(float, float, float)>&& f)
{
  // f is some function taking event time + z of hit and returns final "digit" time
  LOG(INFO) << "BR NAME " << branchname;
  auto br = chain.GetBranch(branchname);
  if (!br) {
    std::cerr << "No branch found\n";
    return;
  }

  auto nentries = br->GetEntries();
  hitvectors.resize(nentries, nullptr);

  // do the filtering
  for (int entry = 0; entry < nentries; ++entry) {
    if (tmin > f(eventrecords[entry].timeNS, 0, 0)) {
      continue;
    }
    if (tmax < f(eventrecords[entry].timeNS, 0, 250)) {
      break;
    }

    br->SetAddress(&hitvectors[entry]);
    br->GetEntry(entry);

    int groupid = -1;
    auto groups = hitvectors[entry];
    for (auto& singlegroup : *groups) {
      std::cout << "This Group is in sector " << o2::TPC::Sector::ToSector(singlegroup.getHit(0).getPos()) << "\n";
      groupid++;
      auto zmax = singlegroup.mZAbsMax;
      auto zmin = singlegroup.mZAbsMin;
      // in case of secondaries, the time ordering may be reversed
      if (zmax < zmin)
        std::swap(zmax, zmin);
      assert(zmin <= zmax);
      // auto tof = singlegroup.
      float tmaxtrack = f(eventrecords[entry].timeNS, 0., zmin);
      float tmintrack = f(eventrecords[entry].timeNS, 0., zmax);
      std::cout << tmintrack << " & " << tmaxtrack << "\n";
      assert(tmaxtrack >= tmintrack);
      if (tmin > tmaxtrack || tmax < tmintrack) {
        std::cout << "DISCARDING " << groupid << " OF ENTRY " << entry << "\n";
        continue;
      }
      // need to record index of the group
      hitids.emplace_back(entry, groupid);
    }
  }
}

// TPC hit selection lambda
auto fTPC = [](float tNS, float tof, float z) {
  // returns time in NS
  return tNS + o2::TPC::SAMPAProcessing::getDriftTime(z) * 1000 + tof;
};

DataProcessorSpec getTPCDriftTimeDigitizer(int sector, int channel, bool cachehits)
{
  TChain* simChain = new TChain("o2sim");
  std::stringstream branchnamestream;
  branchnamestream << "TPCHitsShiftedSector" << sector;
  std::string branchname = branchnamestream.str();

  auto doit = [simChain, branchname, sector](ProcessingContext& pc) {
    // have to do a loop over drift times
    std::vector<std::vector<o2::TPC::HitGroup>*> hitvectors; // "TPCHitVector"
    std::vector<o2::TPC::TPCHitGroupID> hitids;              // "TPCHitIDs"

    // obtain collision times
    // access data
    auto dataref = pc.inputs().get("timeinput");
    auto header = o2::header::get<const o2::header::DataHeader>(dataref.header);
    LOG(INFO) << "PAYLOAD SIZE " << header->payloadSize;

    auto timesview = DataRefUtils::as<o2::MCInteractionRecord>(dataref);

    // detect possible drift times
    double maxtime = 0;
    for (auto e : timesview) {
      maxtime = std::max(maxtime, e.timeNS);
    }
    auto ndrifts = 1 + (int) maxtime / 100000;
    LOG(INFO) << "NDRIFTS " << ndrifts << "\n";

    for(int drift = 1; drift <= ndrifts; ++drift) {

    // load filtered hits
    auto starttime = (drift-1) * 100000;
	auto endtime = drift * 100000;
    getHits(*simChain, timesview, hitvectors, hitids, branchname.c_str(), starttime, endtime, fTPC);

    LOG(INFO) << "DRIFTTIME " << drift << " SECTOR " << sector << " : SELECTED " << hitids.size() << " IDs\n ";

    // perform digitization
    // digitizer.Process(hitvectors, hitids);
    }
  };

  // init function return a lambda taking a ProcessingContext
  auto initIt = [simChain, doit](InitContext& ctx) {
    // setup the input chain
    simChain->AddFile("o2sim.root");

    //digitizer.Init();
    //digitizer.SetOutput("o2sim_digit_sector.root");
    return doit;
  };

  std::stringstream id;
  id << "TPCDigitizer" << sector;
  return DataProcessorSpec{
    id.str().c_str(), Inputs{ InputSpec{ "timeinput", "SIM", "EVENTTIMES", static_cast<SubSpecificationType>(channel),
                                         InputSpec::Timeframe } },
    Outputs{
      // define channel by triple of (origin, type id of data to be sent on this channel, subspecification)
    },
    AlgorithmSpec{ initIt }, Options{ /*{ "simFile", VariantType::String, "o2sim.root", { "Sim input filename" } }*/ }
  };
}
}
}
