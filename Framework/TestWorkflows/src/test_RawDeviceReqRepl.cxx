// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/ExternalFairMQDeviceProxy.h"
#include "FairMQLogger.h"
#include "Headers/HeartbeatFrame.h"

#include "FairMQMessage.h"
#include "FairMQChannel.h"
#include "FairMQTransportFactory.h"

using namespace o2::framework;


using DataHeader = o2::header::DataHeader;
using DataOrigin = o2::header::DataOrigin;

// A simple workflow which takes heartbeats from
// a DPL device as invocation of its own algorithm
// and fetching tasks/data from an external task server
// whenever it wants
WorkflowSpec defineDataProcessing(ConfigContext const&specs) {

  return WorkflowSpec{

    // WORKER
    DataProcessorSpec{
      "Worker",
      Inputs{ InputSpec{ "heartbeat", "BAR", "FOO" } }, // we get the heartbeat as well as task data
      {},
      AlgorithmSpec{

        [](InitContext& ictx) {
          auto factory = FairMQTransportFactory::CreateTransportFactory("zeromq");
          auto channel = std::make_shared<FairMQChannel>("requestchannel", "req", factory);
          channel->Connect("tcp://localhost:5450");
          channel->ValidateChannel();

          return [channel](ProcessingContext& ctx) {
            LOG(INFO) << "INVOKED";
            static int counter = 0;
            counter++;
            if (counter % 10 == 0) {

              int secret = 231;
              FairMQMessagePtr request(channel->NewSimpleMessage(secret));

              FairMQMessagePtr reply(channel->NewMessage());

              if (channel->Send(request, 2000)) {
                if (channel->Receive(reply, 2000)) {
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
