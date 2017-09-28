// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Backend.cxx
/// \brief Implementation of the Backend class
/// \author Charis Kouzinopoulos <charalampos.kouzinopoulos@cern.ch>

#include "CCDB/Backend.h"
#include "request.pb.h"

using namespace o2::CDB;
using namespace std;

void Backend::serialize(std::string*& messageString, const std::string& key, const std::string& operationType,
                        const std::string& dataSource, const std::string& object /*= std::vector<char>()*/)
{
  messaging::RequestMessage* requestMessage = new messaging::RequestMessage;
  requestMessage->setCommand(operationType);
  requestMessage->setDatasource(dataSource);
  requestMessage->setKey(key);

  if (object.length() > 0) {
    requestMessage->setValue(object);
  }

  requestMessage->SerializeToString(messageString);

  delete requestMessage;
}
