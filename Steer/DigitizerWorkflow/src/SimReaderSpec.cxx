// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "SimReaderSpec.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Framework/ControlService.h"
#include "Headers/DataHeader.h"
#include "Steer/HitProcessingManager.h"
#include <FairMQLogger.h>
#include <TMessage.h> // object serialization
#include <memory>  // std::unique_ptr
#include <cstring> // memcpy
#include <string>  // std::string
#include <cassert>

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

namespace o2 {
namespace steer {
DataProcessorSpec getSimReaderSpec() {
  // set up the processing function
	 // return processing function
	     // capture by reference ok since global instance

  auto doit = [](ProcessingContext& pc) {
    auto& mgr = steer::HitProcessingManager::instance();
    auto eventrecords = mgr.getRunContext().getEventRecords();

    // send data via data allocator
    // for times a flat buffer
   // pc.allocator().snapshot(OutputSpec{ "SIM", "EVENTTIMES", 0, OutputSpec::Timeframe }, eventrecords);

    auto msg = new TMessage();
    auto cl = TClass::GetClass(typeid(decltype(eventrecords)));
    assert(cl);
    msg->WriteObjectAny(&eventrecords, cl);

    pc.allocator().adopt(OutputSpec{ "SIM", "EVENTTIMES", 0, OutputSpec::Timeframe }, msg);

    //static int counter = 0;
    //pc.allocator().snapshot(OutputSpec{ "SIM", "EVENTTIMES", 0, OutputSpec::Timeframe }, counter++);

    // do this only one
    // pc.services().get<ControlService>().readyToQuit(true);
  };

  // init function return a lambda taking a ProcessingContext
  auto initIt = [doit](InitContext& ctx) {
    // initialize fundamental objects
    auto& mgr = steer::HitProcessingManager::instance();
    mgr.addInputFile(ctx.options().get<std::string>("simFile").c_str());
    mgr.setupRun();

    LOG(INFO) << "Initializing Spec ... have " << mgr.getRunContext().getEventRecords().size() << " times ";
    return doit;
  };

  return DataProcessorSpec{
	  /*ID*/ "SimReader",
	  /*INPUT CHANNELS*/ Inputs{},
	  /*OUTPUT CHANNELS*/ Outputs{
		 // define channel by triple of (origin, type id of data to be sent on this channel, subspecification)
		 OutputSpec{"SIM", "EVENTTIMES", 0, OutputSpec::Timeframe},

		 // OutputSpec{"SIM", "TPCHITS", 0, OutputSpec::Timeframe},

		 // channel for kinematics information
		 // OutputSpec{"SIM", "MCTRACKS", 0, OutputSpec::Timeframe}
	  },
     /* ALGORITHM */
      AlgorithmSpec{ initIt },
     /* OPTIONS */
     Options {
       { "simFile", VariantType::String, "o2sim.root", { "Sim input filename" }}
     }
  };
}

}
}



