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
#include "Framework/CompletionPolicy.h"
#include "Framework/DeviceSpec.h"

using namespace o2::framework;

void customize(std::vector<CompletionPolicy> &policies) {
  auto matcher = [](DeviceSpec const &device) -> bool {
    return device.name == "Worker";
  };
  auto policy = [](gsl::span<PartRef const> const &inputs) -> CompletionPolicy::CompletionOp {
    return CompletionPolicy::CompletionOp::Consume;
  };
  policies.push_back({CompletionPolicy{"process-any", matcher, policy}});
}

#include "Framework/runDataProcessing.h"

using DataHeader = o2::header::DataHeader;
using DataOrigin = o2::header::DataOrigin;

// A simple workflow which takes heartbeats from
// a DPL device as invocation of its own algorithm
// and fetching tasks/data from an external task server
// whenever it wants
WorkflowSpec defineDataProcessing(ConfigContext const&specs) {

  auto outspecTask = OutputSpec{ "TSK", "TASKS", 0, Lifetime::Timeframe };
  auto inspecTask = InputSpec{ "task", "TSK", "TASKS", 0, Lifetime::Timeframe };

  return WorkflowSpec{
    specifyExternalFairMQDeviceProxy("TaskServer",
                                     { outspecTask },
                                     "type=req,method=connect,address=tcp://localhost:5450,rateLogging=0",
                                     o2DataModelAdaptor(outspecTask, 0, 1)),

    // WORKER
    DataProcessorSpec{
      "Worker",
      Inputs{ InputSpec{ "heartbeat", "BAR", "FOO" }, inspecTask }, // we get the heartbeat as well as task data
      {},
      AlgorithmSpec{

        [](InitContext& ictx) {
          auto factory = FairMQTransportFactory::CreateTransportFactory("zeromq");
          auto channel = std::make_shared<FairMQChannel>("requestchannel", "req", factory);
          channel->Connect("tcp://localhost:5450");
          channel->ValidateChannel();

          return [channel](ProcessingContext& ctx) {
            if (ctx.inputs().isValid("heartbeat")) {
        	  LOG(INFO) << "INVOKED FROM HEARTBEAT";
            }
            if (ctx.inputs().isValid("task")) {
        	  LOG(INFO) << "INVOKED FROM TASK";
            }

            static int counter = 0;
            counter++;
            if (counter % 10 == 0) {

              int secret = 231;
              FairMQMessagePtr request(channel->NewSimpleMessage(secret));

              FairMQMessagePtr reply(channel->NewMessage());

              if (channel->Send(request, 2000) > 0) {
                if (channel->Receive(reply, 2000) > 0) {
                  LOG(INFO) << "Received answer with " << reply->GetSize();
                } else {
                  LOG(INFO) << "No answer";
                }
              } else {
                LOG(INFO) << "DID not send";
              }
            }
          };
        }

      } },

    // HB
    DataProcessorSpec{
      "HeartBeat",
      Inputs{},
      Outputs{ OutputSpec{ { "label" }, "BAR", "FOO" } },
      AlgorithmSpec{
        [](ProcessingContext& ctx) {
          static int hbcounter = 0;
          hbcounter++;
          ctx.outputs().snapshot( OutputRef{"label"}, hbcounter);

        } } }
  };
}
