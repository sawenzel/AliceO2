// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TOFClusterWriterSpec2.h"
#include <SimulationDataFormat/MCCompLabel.h>
#include <SimulationDataFormat/MCTruthContainer.h>
#include "Utils/MakeRootTreeWriterSpec.h"
#include "DataFormatsTOF/Cluster.h"
#include <vector>

using namespace o2::framework;
using namespace o2::tof;

template<typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;

DataProcessorSpec o2::tof::getTOFClusterWriterSpec2()
{
  return MakeRootTreeWriterSpec(
    "TOFClusterWriter2", // process name
    "tofclusters2.root", // default file name
    "o2sim",             // default tree name
    1,                   // default number of events (how often does it receive data?)
    BranchDefinition<std::vector<o2::tof::Cluster>>{ InputSpec{ "tofclusters", "TOF", "CLUSTERS" }, "cluster" },
    BranchDefinition<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>{
      InputSpec{ "toflabels", "TOF", "CLUSTERSMCTR" }, "labels" } // branch config
    )();
};
