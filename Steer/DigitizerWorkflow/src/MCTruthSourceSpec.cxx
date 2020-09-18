// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Framework/Lifetime.h"
#include <SimulationDataFormat/MCTruthContainer.h>
#include <SimulationDataFormat/ConstMCTruthContainer.h>
#include <SimulationDataFormat/MCCompLabel.h>
#include "Framework/Task.h"
#include "Framework/Logger.h"

using namespace o2::framework;

namespace o2
{

class MCTruthSourceTask : public o2::framework::Task
{
 public:
  MCTruthSourceTask(bool newmctruth) : mNew{newmctruth} {}

  using TruthElement = o2::MCCompLabel;
  using Container = o2::dataformats::MCTruthContainer<TruthElement>;

  void init(framework::InitContext& ic) override
  {
    LOG(INFO) << "Initializing MCTruth source";
    mSize = ic.options().get<int>("size");
  }

  void run(framework::ProcessingContext& pc) override
  {
    if (mFinished) {
      // modify the buffer a bit
      sleep(2);
      LOG(INFO) << "MODIFYING BUFFER";
      mBuffer[0] = -112;
      if (mBuffer2Addr) {
        // trying to modify the labels in shared memory
        TruthElement* element = (TruthElement*)mBuffer2Addr;
        TruthElement e(-111, -111, -111);
        LOG(INFO) << "WRITING NEW LABEL " << e;
        std::memcpy(element, &e, sizeof(TruthElement));
      }
      sleep(5);
      return;
    }
    LOG(INFO) << "Creating MCTruth container";

    // create a very large container and stream it to TTree
    for (int i = 0; i < mSize; ++i) {
      mLabels1.addElement(i, TruthElement(i, i, i));
      mLabels1.addElement(i, TruthElement(i + 1, i, i));
      mLabels2.addElement(i, TruthElement(-i, -i, -i));
      mLabels2.addElement(i, TruthElement(-i - 1, -i, -i));
    }

    if (mNew) {
      LOG(INFO) << "New serialization";
      // we need to flatten it and write to managed shared memory container
      auto& sharedlabels = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<TruthElement>>(Output{"TST", "LABELS", 0, Lifetime::Timeframe});
      mLabels1.flatten_to(sharedlabels);
      auto s = sharedlabels.size();
      // check last truth element
      TruthElement checkelement; // = (TruthElement)sharedlabels[s-sizeof(TruthElement)]; // last byte
      std::memcpy(&checkelement, &sharedlabels[s - sizeof(TruthElement)], sizeof(TruthElement));

      // remember the address of the last label
      mBuffer2Addr = (void*)(&(sharedlabels[s - sizeof(TruthElement)]));

      LOG(INFO) << "CHECK " << s << " SIZE AND ELEMENT " << checkelement << " vs " << mLabels1.getLabels(mSize - 1)[1];
      sleep(1);
    } else {
      LOG(INFO) << "Old serialization";
      pc.outputs().snapshot({"TST", "LABELS", 0, Lifetime::Timeframe}, mLabels1);
      sleep(1);
    }

    int size = 111;
    auto buffer = pc.outputs().make<char>(Output{"TST", "SHM", 0, Lifetime::Timeframe}, size);
    for (int i = 0; i < size; ++i)
      buffer[i] = size - i;
    LOG(INFO) << "SENDING BUFFER ADDRESS " << (void*)&(buffer[0]) << " AND SIZE " << buffer.size();
    mBuffer = &(buffer[0]);

    // we should be only called once; tell DPL that this process is ready to exit
    // pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    mFinished = true;
  }

 private:
  bool mFinished = false;
  int mSize = 0;
  bool mNew = false;
  char* mBuffer = nullptr;
  void* mBuffer2Addr = nullptr;

  o2::dataformats::MCTruthContainer<TruthElement> mLabels1; // labels which get filled
  o2::dataformats::MCTruthContainer<TruthElement> mLabels2; // labels which get filled
};

o2::framework::DataProcessorSpec getMCTruthSourceSpec(bool newmctruth)
{
  // create the full data processor spec using
  //  a name identifier
  //  input description
  //  algorithmic description (here a lambda getting called once to setup the actual processing function)
  //  options that can be used for this processor (here: input file names where to take the hits)
  std::vector<OutputSpec> outputs;
  outputs.emplace_back("TST", "LABELS", 0, Lifetime::Timeframe);
  outputs.emplace_back("TST", "SHM", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    "MCTruthSource",
    Inputs{},
    outputs,
    AlgorithmSpec{adaptFromTask<MCTruthSourceTask>(newmctruth)},
    Options{
      {"size", VariantType::Int, 100000, {"Sample size"}}}};
}

} // end namespace o2
