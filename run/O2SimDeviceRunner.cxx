// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "runFairMQDevice.h"
#include "O2SimDevice.h"

#include <string>

namespace bpo = boost::program_options;

void addCustomOptions(bpo::options_description& options)
{

}

FairMQDevice* getDevice(const FairMQProgOptions& config)
{
  return new o2::devices::O2SimDevice();
}
