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
#include <Steer/InteractionSampler.h>
#include "Headers/DataHeader.h"

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
using DataRefUtils = o2::framework::DataRefUtils;

namespace o2 {
namespace steer {
DataProcessorSpec getCollisionTimePrinter() {
  // set up the processing function

  // init function return a lambda taking a ProcessingContext
  auto doIt = [](ProcessingContext& pc) {
    std::cout << " ######### ";

    // access data
    auto dataref = pc.inputs().get("input");
    auto header = o2::header::get<const o2::header::DataHeader>(dataref.header);
	LOG(INFO) << "PAYLOAD SIZE " << header->payloadSize;

    //auto view = DataRefUtils::as<int>(dataref);
    //LOG(INFO) << "## " << view[0] << "\n";

	auto msg = DataRefUtils::as<TMessage>(dataref);
	msg->Print();
	return;
    using T = std::vector<o2::MCInteractionRecord>;
    auto cl = TClass::GetClass(typeid(T));
    assert(cl);
    auto records = static_cast<T*>(msg->ReadObjectAny(cl));
    assert(records);

    LOG(INFO) << "GOT " << records->size() << "times";
    int counter=0;
    for (auto& collrecord : *records) {
      LOG(INFO) << "TIME " << counter++ << " : " << collrecord.timeNS;
    }
  };

  return DataProcessorSpec{
	  /*ID*/ "CollTimePrinter",
	  /*INPUT CHANNELS*/ Inputs{
		InputSpec{"input", "SIM", "EVENTTIMES", 0, InputSpec::Timeframe}
	  },
	  /*OUTPUT CHANNELS*/ Outputs{},
     /* ALGORITHM */
      AlgorithmSpec(doIt),
     /* OPTIONS */
     Options {}
  };
}

}
}



