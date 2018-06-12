// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/ExternalFairMQDeviceProxy.h"
#include "FairMQLogger.h"
#include "Headers/HeartbeatFrame.h"

#include "FairMQMessage.h"
#include "FairMQChannel.h"
#include "FairMQTransportFactory.h"

using namespace o2::framework;

#include "Framework/runDataProcessing.h"

using DataHeader = o2::header::DataHeader;
using DataOrigin = o2::header::DataOrigin;

// A simple workflow which takes heartbeats from
// a DPL device as invocation of its own algorithm
// and fetching tasks/data from an external task server
// whenever it wants
WorkflowSpec defineDataProcessing(ConfigContext const&specs) {
  auto outspecHB = OutputSpec{ {"hbout"}, o2::header::DataOrigin("SMPL"), o2::header::gDataDescriptionHeartbeatFrame, 0, Lifetime::Timeframe};

  auto inspecHB = InputSpec{"heartbeat",
                          o2::header::DataOrigin("SMPL"),
                          o2::header::gDataDescriptionHeartbeatFrame, 0, Lifetime::Timeframe};

  auto outspecTask = OutputSpec{ "TSK", "TASKS", 0, Lifetime::Timeframe };
  auto inspecTask = InputSpec{ "task", "TSK", "TASKS", 0, Lifetime::Timeframe };

  return WorkflowSpec{

    // TASKSERVER
    specifyExternalFairMQDeviceProxy("TaskServer",
                                     { outspecTask },
                                     "type=req,method=connect,address=tcp://localhost:5450,rateLogging=0",
                                     o2DataModelAdaptor(outspecTask, 0, 1)),

    // WORKER
    DataProcessorSpec{
      "Worker",
      Inputs{ inspecHB, inspecTask }, // we get the heartbeat as well as task data
      {},
      AlgorithmSpec{
        [](ProcessingContext& ctx) {

          auto factory = FairMQTransportFactory::CreateTransportFactory("zeromq");
          auto channel = FairMQChannel("", "req", factory);
          channel.Connect("tcp://localhost:5450");
          channel.ValidateChannel();

        } } },

    // HB
    DataProcessorSpec{
      "HeartBeat",
      Inputs{},
      Outputs{ outspecHB },
      AlgorithmSpec{
        [](ProcessingContext& ctx) {
          static int hbcounter = 0;
          hbcounter++;
          ctx.outputs().snapshot(OutputRef{ "hbout" }, hbcounter);
        } } }
  };
}
