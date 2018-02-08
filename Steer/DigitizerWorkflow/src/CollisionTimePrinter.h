#ifndef STEER_DIGITIZERWORKFLOW_SRC_COLLISIONTIMEPRINTER_H_
#define STEER_DIGITIZERWORKFLOW_SRC_COLLISIONTIMEPRINTER_H_

#include "Framework/DataProcessorSpec.h"

namespace o2 {
namespace steer {
  o2::framework::DataProcessorSpec getCollisionTimePrinter(int subchannel);
}
}

#endif /* STEER_DIGITIZERWORKFLOW_SRC_COLLISIONTIMEPRINTER_H_ */
