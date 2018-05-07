// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/InputSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataSampling.h"
#include "Framework/ParallelContext.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ControlService.h"

using namespace o2::framework;

// This test is testing the outputs().flush() API
// This demonstrates that the algorithm of dataProducer is only executed once
// .. while that of dataSink multiple times
void defineDataProcessing(std::vector<DataProcessorSpec>& specs)
{
  DataProcessorSpec dataProducer{ "dataProducer",
                                  Inputs{},
                                  { OutputSpec{ "FOO", "NUMBER", 0, Lifetime::Condition } },
                                  AlgorithmSpec{ [](ProcessingContext& ctx) {
                                    // this is testing the flush method
                                    // send some data in a loop
                                    static bool finished = false;
                                    if (finished) {
                                      return;
                                    }

                                    int i = -1;
                                    for (int j = 10; j >= -1; --j) {
                                      i = j;
                                      LOG(INFO) << "SENDING MESSAGE " << i;
                                      ctx.outputs().snapshot(Output{ "FOO", "NUMBER", 0 }, i);
                                      ctx.outputs().flush();
                                    }
                                    ctx.services().get<ControlService>().readyToQuit(false);
                                    finished = true;
                                    return;
                                  } } };

  specs.push_back(dataProducer);

  DataProcessorSpec dataSink{ "dataSink", Inputs{ InputSpec{ "input", "FOO", "NUMBER", 0, Lifetime::Condition } },
                              Outputs{}, AlgorithmSpec{ [](ProcessingContext& ctx) {
                                static bool finished = false;
                                if (finished) {
                                  return;
                                }
                                auto iptr = ctx.inputs().get<int>("input");
                                const int i = *(iptr.get());
                                LOG(INFO) << "GOT MESSAGE " << i;

                                // take value of -1 as signal to stop
                                if (i == -1) {
                                  finished = true;
                                  ctx.services().get<ControlService>().readyToQuit(false);
                                }
                                return;
                              } } };
  specs.push_back(dataSink);
}
