#ifndef STEER_DIGITIZERWORKFLOW_SRC_HMPIDDIGITIZERSPEC_H_
#define STEER_DIGITIZERWORKFLOW_SRC_HMPIDDIGITIZERSPEC_H_

#include "Framework/DataProcessorSpec.h"

namespace o2
{
namespace hmpid
{

o2::framework::DataProcessorSpec getHMPIDDigitizerSpec(int channel);

} // end namespace hmpid
} // end namespace o2

#endif /* STEER_DIGITIZERWORKFLOW_SRC_HMPIDDIGITIZERSPEC_H_ */
