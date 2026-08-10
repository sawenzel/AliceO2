// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test MCUtils
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "SimulationDataFormat/MCUtils.h"

using namespace o2::mcutils;

// the navigator has to reject indices that do not address a track in the
// container, in particular the one just past its end
BOOST_AUTO_TEST_CASE(MCTrackNavigator_outOfRange_test)
{
  std::vector<o2::MCTrack> container(3);

  o2::MCTrack probe;

  probe.SetMotherTrackId(static_cast<int>(container.size()));
  BOOST_CHECK(MCTrackNavigator::getMother(probe, container) == nullptr);

  probe.SetFirstDaughterTrackId(static_cast<int>(container.size()));
  BOOST_CHECK(MCTrackNavigator::getDaughter(probe, container) == nullptr);
  BOOST_CHECK(MCTrackNavigator::getDaughter0(probe, container) == nullptr);

  probe.SetLastDaughterTrackId(static_cast<int>(container.size()));
  BOOST_CHECK(MCTrackNavigator::getDaughter1(probe, container) == nullptr);

  // a negative index means "no such track" and must be rejected as well
  probe.SetMotherTrackId(-1);
  BOOST_CHECK(MCTrackNavigator::getMother(probe, container) == nullptr);

  // an index inside the container still resolves
  probe.SetMotherTrackId(2);
  BOOST_CHECK(MCTrackNavigator::getMother(probe, container) == &container[2]);
}
