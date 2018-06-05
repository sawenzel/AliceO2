// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_DEVICES_HITMERGER_H_
#define ALICEO2_DEVICES_HITMERGER_H_

#include <memory>
#include "FairMQMessage.h"
#include <FairMQDevice.h>
#include <FairLogger.h>
#include "../macro/o2sim.C"
#include "TVirtualMC.h"
#include <SimulationDataFormat/Stack.h>
#include <SimulationDataFormat/PrimaryChunk.h>
#include <DetectorsCommonDataFormats/DetID.h>
#include <gsl/gsl>
#include "TFile.h"
#include "TTree.h"
#include <memory>
#include <TMessage.h>
#include <FairMQParts.h>
#include <ctime>
#include <TStopwatch.h>
#include <sstream>
#include <cassert>
#include "FairSystemInfo.h"
#include "Steer/InteractionSampler.h"

#include "O2HitMerger.h"
#include <DetectorsCommonDataFormats/DetID.h>
#include <TPCSimulation/Detector.h>
#include <ITSSimulation/Detector.h>
#include <MFTSimulation/Detector.h>
#include <EMCALSimulation/Detector.h>
#include <TOFSimulation/Detector.h>
#include <TRDSimulation/Detector.h>
#include <FITSimulation/Detector.h>
#include <HMPIDSimulation/Detector.h>
#include <PHOSSimulation/Detector.h>

namespace o2 {
namespace devices {

class O2HitMerger : public FairMQDevice
{

  class TMessageWrapper : public TMessage
  {
   public:
    TMessageWrapper(void* buf, Int_t len) : TMessage(buf, len) { ResetBit(kIsOwner); }
    ~TMessageWrapper() override = default;
  };

 public:
  /// Default constructor
  O2HitMerger()
  {
    initDetInstances();

    // ideally we have a specialized handler function per data channel
    OnData("mctracks", &O2HitMerger::HandleMCTrackData);
    OnData("itshits", &O2HitMerger::HandleITSHits);
    OnData("tpchits", &O2HitMerger::HandleTPCHits);
    OnData("simdata", &O2HitMerger::handleSimData);
    mTimer.Start();
  }

  /// Default destructor
  ~O2HitMerger()
  {
    FairSystemInfo sysinfo;
    mOutTree->SetEntries(mEntries);
    mOutTree->Write();
    mOutFile->Close();
    LOG(INFO) << "TIME-STAMP " << mTimer.RealTime() << "\t";
    mTimer.Continue();
    LOG(INFO) << "MEM-STAMP " << sysinfo.GetCurrentMemory() / (1024. * 1024) << " "
              << sysinfo.GetMaxMemory() << " MB\n";
  }

 private:
  /// Overloads the InitTask() method of FairMQDevice
  void InitTask() final
  {
    mOutFile = new TFile("o2sim_merged_hits.root", "RECREATE");
    mOutTree = new TTree("o2sim", "o2sim");

    // make some branches
    //mOutTree->Branch("MCTracks", &mTracksPtr);
    //mOutTree->Branch("ITSHits", &mITSHits);

    //for (int s = 0; s < 36; ++s) {
    //  std::stringstream sec;
    //  sec << "TPCHitsShiftedSector" << s;
    //  mOutTree->Branch(sec.str().c_str(), &mTPCHits);
    //}
  }

  template <typename T, typename I>
  static gsl::span<T> getViewAndMetaData(FairMQMessagePtr& data, I& meta)
  {
    // unpack the message to retreive the actual vector + the meta information
    // that was pickybacked with the message
    const auto view = gsl::span<T>(static_cast<T*>(data->GetData()), data->GetSize() / sizeof(T));
    auto ntrailingelements = o2::Data::getTrailingObjectFromVector(view, meta);
    auto actualview = gsl::span<T>(static_cast<T*>(data->GetData()), data->GetSize() / sizeof(T) - ntrailingelements);
    return actualview;
  }

  // general flush function, merging and flushing to ROOT file
  // works for std::vector data
  // cleans up maps
  template <typename T>
  static void flush(int eventid,
                    TTree& outtree,
                    const char* brname,
                    std::vector<T>& mergeVector,
                    std::map<uint32_t, std::vector<FairMQMessagePtr>>& messages,
                    size_t totalsize)
  {
    LOG(INFO) << "Flushing event " << eventid << " for branch " << brname << "\n";
    auto br = outtree.GetBranch(brname);
    if (!br) {
      LOG(WARNING) << "Did not find brname " << brname << " ... not filling \n";
      return;
    }

    // for the moment copying all buffers into one vector
    // and then flushing as FairRootManager would do
    auto iter = messages.find(eventid);
    auto& messagelist = iter->second;

    mergeVector.reserve(totalsize);
    mergeVector.clear();
    for (auto& message : messagelist) {
      FairMQMessagePtr localmessage;
      localmessage = std::move(message);
      o2::Data::SubEventInfo tmp;
      auto view = getViewAndMetaData<T>(localmessage, tmp);
      // interpret as view and copy data to merge container
      std::copy(view.begin(), view.end(), std::back_inserter(mergeVector));

      // NEED TO FIX THE TRACK IDS UNLESS WE WRITE 1ENTRY = 1PART

      // here the original message will
      // be destroyed/freed since it was moved to
      // localmessage
    }
    LOG(INFO) << "Filled branch " << brname << " with " << mergeVector.size() << " items \n";
    br->Fill();
    outtree.SetEntries(eventid);
    // finally we can cleanup the map by removing this event
    messages.erase(iter);
  }

  template <typename T, typename V>
  V insertAdd(std::map<T, V>& m, T const& key, V value)
  {
    const auto iter = m.find(key);
    V accum{ 0 };
    if (iter != m.end()) {
      iter->second += value;
      accum = iter->second;
    } else {
      m.insert(std::make_pair(key, value));
      accum = value;
    }
    return accum;
  }

  template <typename T>
  bool isDataComplete(T checksum, T nparts)
  {
    return checksum == nparts * (nparts + 1) / 2;
  }

  bool HandleMCTrackData(FairMQMessagePtr& data, int /*index*/)
  {
    o2::Data::SubEventInfo info;
    auto mctracks = getViewAndMetaData<o2::MCTrack>(data, info);
    auto accum = insertAdd<uint32_t, uint32_t>(mTracksPartsCheckSum, info.eventID, (uint32_t)info.part);

    mTrackMessages[info.eventID].push_back(std::move(data));
    auto totalsize = insertAdd<uint32_t, size_t>(mTracksTotalSize, info.eventID, mctracks.size());

    // are tracks complete for one event?
    if (isDataComplete<uint32_t>(accum, info.nparts)) {
      flush(info.eventID, *mOutTree, "MCTracks", mMergedTracks, mTrackMessages, totalsize);

      // last event seen -> exit
      if (info.eventID == info.maxEvents) {
        return false;
      }
    }
    FairSystemInfo sysinfo;
    LOG(INFO) << "TIME-STAMP " << mTimer.RealTime() << "\t";
    mTimer.Continue();
    LOG(INFO) << "MEM-STAMP " << sysinfo.GetCurrentMemory() / (1024. * 1024) << " " << sysinfo.GetMaxMemory() << " MB\n";
    return true;
  }

  bool HandleITSHits(FairMQMessagePtr& data, int /*index*/)
  {
    // FIXME: this code is almost same as for MCTracks -> generalize
    o2::Data::SubEventInfo info;
    auto itshits = getViewAndMetaData<o2::ITSMFT::Hit>(data, info);
    auto accum = insertAdd<uint32_t, uint32_t>(mITSPartsCheckSum, info.eventID, (uint32_t)info.part);
    mITSMessages[info.eventID].push_back(std::move(data));
    auto totalsize = insertAdd<uint32_t, size_t>(mITSTotalSize, info.eventID, itshits.size());

    // TOY FOR THE MOMENT; MAKE A TIMESTAMP AND FORWARD TO DIGITIZER
    // TODO: NEED TO SEND ALL DATA FOR EVENT i BEFORE ANY DATA FOR NEXT EVENT j
#ifdef DIRECTFORWARD
    FairMQParts parts;

    auto t = mInteractionSampler.generateCollisionTime();
    LOG(INFO) << "SAMPLED TIMED " << t.timeNS << "\n";

    // STAMP TIME (same event for all)
    info.eventtime = t.timeNS;
    parts.AddPart(NewSimpleMessage(info));
    itshits[itshits.size() - 1].Print(nullptr);
    auto buffer = (char*)&itshits[0];
    auto buffersize = itshits.size() * sizeof(o2::ITSMFT::Hit);
    FairMQMessagePtr message(NewMessage(buffer, buffersize,
                                        [](void* x, void* hint) {}, buffer));
    parts.AddPart(std::move(message));
    Send(parts, "itsdigitizerchannel");
#endif

    // are tracks complete for one event?
    if (isDataComplete<uint32_t>(accum, info.nparts)) {
      flush(info.eventID, *mOutTree, "ITSHits", mITSHits, mITSMessages, totalsize);

      if (info.eventID == info.maxEvents) {
        // here we are done -- send an empty message to digitizer to tell him that he can stop
        FairMQParts final;
        final.AddPart(NewSimpleMessage(info));
        Send(final, "itsdigitizerchannel");
        return false;
      }
    }
    return true;
  }

  void flushTPC(int eventid)
  {
    LOG(INFO) << "Flushing TPC\n";
    for (int s = 0; s < 36; ++s) {
      std::stringstream sec;
      sec << "TPCHitsShiftedSector" << s;
      auto br = mOutTree->GetBranch(sec.str().c_str());
      assert(br);
      auto iter = mTPCMessages[s].find(eventid);
      if (iter == mTPCMessages[s].end()) {
        continue;
      }
      auto& messagelist = iter->second;
      for (auto& message : messagelist) {
        // message should be a pointer to a TPC
        br->SetAddress(&message);
        br->Fill();
        delete message;
        break;
      }
    }

    mOutTree->SetEntries(eventid);
    // finally we can cleanup the map by removing this event
  }

  bool consumeHits(FairMQParts& data, int& index) {
	auto detIDmessage = std::move(data.At(index++));
    // this should be a detector ID
    LOG(INFO) << detIDmessage->GetSize();
    if (detIDmessage->GetSize() == 4) {
      auto ptr = (int*)detIDmessage->GetData();
      o2::detectors::DetID id(ptr[0]);
      LOG(INFO) << "I1 " << ptr[0] << " NAME " << id.getName();

      // get the detector than can interpret it
      auto detector = mDetectorInstances[id].get();
      detector->fillHitBranch(*mOutTree, data, index);
    }
    return true;
  }

  bool handleSimData(FairMQParts& data, int /*index*/)
  {
    LOG(INFO) << "SIMDATA channel got " << data.Size() << " parts\n";

    // extract the event info
    auto infomessage = std::move(data.At(0));
    o2::Data::SubEventInfo info;
    // could actually avoid the copy
    memcpy((void*)&info, infomessage->GetData(), infomessage->GetSize());

    auto accum = insertAdd<uint32_t, uint32_t>(mITSPartsCheckSum, info.eventID, (uint32_t)info.part);
    // auto totalsize = insertAdd<uint32_t, size_t>(mITSTotalSize, info.eventID, itshits.size());


    // consumeMCTracks(data);
    // consumeTrackReferences();
    // indices 0 and 1 are MCData and TrackRef
    int index = 2;
    while (index < data.Size()) {
      LOG(INFO) << "now consuming at " << index;
      consumeHits(data, index);
    }
    mEntries++;

    if (isDataComplete<uint32_t>(accum, info.nparts)) {
      LOG(INFO) << "EVERYTHING IS HERE FOR EVENT " << info.eventID << "\n";
      if (info.eventID == info.maxEvents) {
        return false;
      }
    }

    return true;
  }

  bool HandleTPCHits(FairMQParts& data, int /*index*/)
  {
    // retrieve the info object
    auto infomessage = std::move(data.At(0));
    o2::Data::SubEventInfo info;
    // could actually avoid the copy
    memcpy((void*)&info, infomessage->GetData(), infomessage->GetSize());
    auto accum = insertAdd<uint32_t, uint32_t>(mTPCPartsCheckSum, info.eventID, (uint32_t)info.part);

    // loop over sectors
    for (int s = 0; s < 36; ++s) {
      // in this case we need to deserialize
      auto tpcmessage = std::move(data.At(1 + s));

      auto message = std::make_unique<TMessageWrapper>(tpcmessage->GetData(), tpcmessage->GetSize());
      auto tpchits = static_cast<std::vector<o2::TPC::HitGroup>*>(message.get()->ReadObjectAny(message.get()->GetClass()));
      mTPCMessages[s][info.eventID].push_back(tpchits);
    }
    if (isDataComplete<uint32_t>(accum, info.nparts)) {
      flushTPC(info.eventID);
      if (info.eventID == info.maxEvents) {
        return false;
      }
    }
    return true;
  }

  std::vector<o2::MCTrack> mMergedTracks; //! accumulator for actual tracks
  std::vector<o2::MCTrack>* mTracksPtr = &mMergedTracks;

  // I think we can get rid of the maps at some stage
  std::map<uint32_t, std::vector<FairMQMessagePtr>> mTrackMessages; //! accumulator for track messages

  std::map<uint32_t, uint32_t> mTracksPartsCheckSum; //! mapping event id -> part checksum used to detect when all info

  // -- ITS --
  std::vector<o2::ITSMFT::Hit> mITSHits;
  std::map<uint32_t, std::vector<FairMQMessagePtr>> mITSMessages;
  std::map<uint32_t, uint32_t> mITSPartsCheckSum; //! mapping event id -> part checksum used to detect when all info

  // -- TPC --
  std::vector<o2::TPC::HitGroup> mTPCHits;

  // for each of the sectors I keep a map to
  // collect messages for each sector
  std::map<uint32_t, std::vector<std::vector<o2::TPC::HitGroup>*>> mTPCMessages[36];
  std::map<uint32_t, uint32_t> mTPCPartsCheckSum; //! mapping event id -> part checksum used to detect when all info

  //  for a given event is here
  std::map<uint32_t, size_t> mTracksTotalSize; //! total number of tracks for a given event

  std::map<uint32_t, size_t> mITSTotalSize; //! total number of tracks for a given event

  std::map<uint32_t, size_t> mTPCTotalSize; //! total number of tracks for a given event

  TFile* mOutFile; //!
  TTree* mOutTree; //!
  int mEntries = 0; //! counts the number of entries in the branches
  TStopwatch mTimer;
  steer::InteractionSampler mInteractionSampler;

  //
  std::vector<std::unique_ptr<o2::Base::Detector>> mDetectorInstances;

  // init detector instances
  void initDetInstances();

};

// init detector instances
void O2HitMerger::initDetInstances()
{
  using o2::detectors::DetID;

  mDetectorInstances.resize(DetID::Last);
  // like a factory of detector objects

  int counter = 0;
  for (int i = DetID::First; i <= DetID::Last; ++i) {
    if (i == DetID::TPC) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::TPC::Detector>(true));
      counter++;
    }
    if (i == DetID::ITS) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::ITS::Detector>(true));
      counter++;
    }
    if (i == DetID::MFT) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::MFT::Detector>());
      counter++;
    }
    if (i == DetID::TRD) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::trd::Detector>(true));
      counter++;
    }
    if (i == DetID::PHS) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::phos::Detector>(true));
      counter++;
    }
    if (i == DetID::EMC) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::EMCAL::Detector>(true));
      counter++;
    }
    if (i == DetID::HMP) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::hmpid::Detector>(true));
      counter++;
    }
    if (i == DetID::TOF) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::tof::Detector>(true));
      counter++;
    }
    if (i == DetID::FIT) {
      mDetectorInstances[i] = std::move(std::make_unique<o2::fit::Detector>(true));
      counter++;
    }
    if ( counter != DetID::Last ) {
      LOG(WARNING) << " O2HitMerger: Some Detectors are potentially missing in this initialization ";
    }
  }
}


} // namespace devices
} // namespace o2

#endif
