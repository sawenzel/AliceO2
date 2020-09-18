// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/DataRefUtils.h"
#include "Framework/ServiceRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/runDataProcessing.h"
#include "Headers/DataHeader.h"
#include "Framework/Logger.h"
#include <gsl/span>

using namespace o2::framework;
using DataHeader = o2::header::DataHeader;
using DataOrigin = o2::header::DataOrigin;
using DataDescription = o2::header::DataDescription;

constexpr size_t SIZE = 10;

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    {"producer",
     {},
     {OutputSpec{"TST", "VECCHAR"},
      OutputSpec{"TST", "CHARBUFF"}},
     AlgorithmSpec{
       [](ProcessingContext& ctx) {
         static int counter = 0;
         static char* buff1ptr = nullptr;
         static char* buff2ptr = nullptr;
         if (counter > 0) {
           sleep(1);
           // second time coming in
           LOG(INFO) << "SECOND TIME";

           // change the content
           for (int i = 0; i < SIZE; ++i) {
             buff1ptr[i] = -i;
             buff2ptr[i] = -i;
           }
           sleep(10);
           ctx.services().get<ControlService>().readyToQuit(QuitRequest::Me);
           return;
         }

         auto& charbuff = ctx.outputs().make<char>(Output{"TST", "CHARBUFF", 0}, SIZE);
         auto& aVec = ctx.outputs().make<std::vector<char>>(Output{"TST", "VECCHAR", 0}, SIZE);

         // fill the buffer
         for (int i = 0; i < SIZE; ++i) {
           charbuff[i] = i;
           aVec[i] = i;
         }
         // remember the buffers for next time
         buff1ptr = &charbuff[0];
         buff2ptr = &aVec[0];
         counter++;
       }}},
    {"consumer",
     {
       InputSpec{"vecchar", "TST", "VECCHAR"},
       InputSpec{"charbuff", "TST", "CHARBUFF"},
     },
     {},
     AlgorithmSpec{
       [](ProcessingContext& ctx) {
         LOG(INFO) << "CONSUMER CALLED";
         auto c = ctx.inputs().get<const char*>("charbuff"); // OK. No copy

         // OBSERVATION .get<const gsl::span<char>> --> wrong result
         // OBSERVATION .get<const gsl::span<const char>> --> wrong type

         auto v = ctx.inputs().get<gsl::span<char>>("vecchar");          // OK. No copy
         auto v2 = ctx.inputs().get<const std::vector<char>>("vecchar"); // COPY MADE
         auto c2 = ctx.inputs().get<gsl::span<char>>("charbuff");        // OK

         LOG(INFO) << "GOT " << v.size() << " BYTES";
         LOG(INFO) << "GOT " << v2.size() << " BYTES";
         LOG(INFO) << "GOT " << c2.size() << " BYTES";

         LOG(INFO) << "FIRST C " << (int)c[SIZE - 1];
         LOG(INFO) << "FIRST C2 " << (int)c2[SIZE - 1];
         LOG(INFO) << "FIRST V " << (int)v[SIZE - 1];
         LOG(INFO) << "FIRST V2 " << (int)v2[SIZE - 1];
         sleep(4);
         LOG(INFO) << "SECOND C " << (int)c[SIZE - 1];
         LOG(INFO) << "SECOND C2 " << (int)c2[SIZE - 1];
         LOG(INFO) << "SECOND V " << (int)v[SIZE - 1];
         LOG(INFO) << "SECOND V2 " << (int)v2[SIZE - 1];

         // SLEEP FOR A BIT TO ALLOW THE SENDER TO MODIFY THE DATA
         ctx.services().get<ControlService>().readyToQuit(QuitRequest::Me);
       }}}};
}
