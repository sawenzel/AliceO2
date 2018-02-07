// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "CollisionTimePrinter.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Headers/DataHeader.h"
#include "Steer/HitProcessingManager.h"
#include <FairMQLogger.h>
#include <TMessage.h> // object serialization
#include <memory>  // std::unique_ptr
#include <cstring> // memcpy
#include <string>  // std::string

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

namespace o2 {
namespace steer {
DataProcessorSpec getCollisionTimePrinter() {
  // set up the processing function

  // init function return a lambda taking a ProcessingContext
  auto doIt = [](ProcessingContext& pc) {
    LOG(INFO) << "doit in print";
  };

  return DataProcessorSpec{
	  /*ID*/ "CollTimePrinter",
	  /*INPUT CHANNELS*/ Inputs{
		InputSpec{"input", "SIM", "EVENTTIMES", 0}
	  },
	  /*OUTPUT CHANNELS*/ Outputs{},
     /* ALGORITHM */
      AlgorithmSpec{ doIt },
     /* OPTIONS */
     Options {}
  };
}

}
}



