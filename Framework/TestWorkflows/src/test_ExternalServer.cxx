/*
 * test_ExternalServer.cxx
 *
 *  Created on: Jun 12, 2018
 *      Author: sandro
 */


#include <FairMQDevice.h>
#include <FairMQMessage.h>
#include <TMessage.h>
#include <TClass.h>
#include <typeinfo>

namespace o2
{
namespace devices
{

class ExtServer : public FairMQDevice
{
 public:
  /// Default constructor
  ExtServer()
  {
    OnData("requestchannel", &ExtServer::HandleRequest);
  }

  /// Default destructor
  ~ExtServer() final = default;

 protected:
  void InitTask() final
  {
    LOG(INFO) << "Init Server device ";
  }

  /// Overloads the ConditionalRun() method of FairMQDevice
  bool HandleRequest(FairMQMessagePtr& request, int /*index*/)
  {
    LOG(INFO) << "GOT A REQUEST WITH SIZE " << request->GetSize();
    std::string requeststring(static_cast<char*>(request->GetData()), request->GetSize());


    int replysecret = -11;
    FairMQMessagePtr reply(NewSimpleMessage(replysecret));
    Send(reply, "requestchannel");

    return true;
  }
};

} // namespace devices
} // namespace o2

#include "runFairMQDevice.h"
namespace bpo = boost::program_options;

void addCustomOptions(bpo::options_description& options)
{
}

FairMQDevice* getDevice(const FairMQProgOptions& config)
{
  return new o2::devices::ExtServer();
}
