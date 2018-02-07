// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/WorkflowSpec.h"
#include "Framework/runDataProcessing.h"

#include "SimReaderSpec.h"
#include "CollisionTimePrinter.h"

using namespace o2::framework;

/// This function is required to be implemented to define the workflow
/// specifications
void defineDataProcessing(WorkflowSpec &specs) {
  specs.clear();

  //
  specs.emplace_back(o2::steer::getSimReaderSpec());
  specs.emplace_back(o2::steer::getCollisionTimePrinter());

  // specs.emplace_back(o2::steer::getTPCHitFilterSpec());
  // specs.emplace_back(o2::steer::getTPCDigitizerSpec());
}
