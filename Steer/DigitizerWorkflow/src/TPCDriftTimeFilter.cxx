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
#include <sstream>

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
namespace o2
{
namespace steer
{
DataProcessorSpec getTPCDriftTimeDigitizer(int sector, int channel, bool cachehits)
{
  auto doit = [](ProcessingContext& pc) {
    // load hits

    // perform filtering

    // perform digitization
  };

  // init function return a lambda taking a ProcessingContext
  auto initIt = [doit](InitContext& ctx) {
    // initialize fundamental objects; such as digitizer
	//std::cout << "DIGITIZER " << ctx.options().get<std::string>("simFile2").c_str() << "\n";
 	return doit;
  };


  std::stringstream id;
  id << "TPCDigitizer" << sector;
  return DataProcessorSpec{
    id.str().c_str(), Inputs{ /*InputSpec{ "hitinput", "SIM", "TPCHITFILES", 0, InputSpec::Timeframe },*/
                            InputSpec{ "timeinput", "SIM", "EVENTTIMES",  static_cast<SubSpecificationType>(channel), InputSpec::Timeframe } },
    Outputs{
      // define channel by triple of (origin, type id of data to be sent on this channel, subspecification)
    },
    AlgorithmSpec{ initIt }, Options{ /*{ "simFile", VariantType::String, "o2sim.root", { "Sim input filename" } }*/ }
  };
}
}
}
