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

#define BOOST_TEST_MODULE Test MCStack class
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "DetectorsBase/Stack.h"
#include "TFile.h"
#include "TMCProcess.h"

using namespace o2;

// unit tests on MC stack
BOOST_AUTO_TEST_CASE(Stack_test)
{
  o2::data::Stack st;
  int a;
  TMCProcess proc{kPPrimary};
  // add a 2 primary particles
  st.PushTrack(1, -1, 0, 0, 0., 0., 10., 5., 5., 5., 0.1, 0., 0., 0., proc, a, 1., 1);
  st.PushTrack(1, -1, 0, 0, 0., 0., 10., 5., 5., 5., 0.1, 0., 0., 0., proc, a, 1., 1);
  BOOST_CHECK(st.getPrimaries().size() == 2);

  {
    // serialize it
    TFile f("StackOut.root", "RECREATE");
    f.WriteObject(&st, "Stack");
    f.Close();
  }

  {
    o2::data::Stack* inst = nullptr;
    TFile f("StackOut.root", "OPEN");
    f.GetObject("Stack", inst);
    BOOST_CHECK(inst->getPrimaries().size() == 2);
  }
}

// convenience wrapper to push a track and return the assigned trackID
static int pushTrack(o2::data::Stack& st, int parentId, TMCProcess proc)
{
  int trackId;
  st.PushTrack(1, parentId, 0, 0., 0., 0., 10., 5., 5., 5., 0.1, 0., 0., 0., proc, trackId, 1., 1);
  return trackId;
}

// unit test for the mother-chain queries
BOOST_AUTO_TEST_CASE(Stack_motherChain_test)
{
  o2::data::Stack st;

  // two primaries; note that primaries do not enter mParticles, only secondaries do
  const auto prim0 = pushTrack(st, -1, kPPrimary);
  const auto prim1 = pushTrack(st, -1, kPPrimary);

  // a secondary of the second primary, and a secondary of that secondary
  const auto secondary = pushTrack(st, prim1, kPHadronic);
  const auto grandSecondary = pushTrack(st, secondary, kPHadronic);

  // a primary has no mother in the transported track numbering. this is required
  // rather than checked, because an ancestry walk that does not recognise a
  // primary never terminates and would hang the checks below
  BOOST_REQUIRE_EQUAL(st.getMotherTrackId(prim0), -1);
  BOOST_REQUIRE_EQUAL(st.getMotherTrackId(prim1), -1);

  // secondaries report the trackID of their direct mother
  BOOST_CHECK_EQUAL(st.getMotherTrackId(secondary), prim1);
  BOOST_CHECK_EQUAL(st.getMotherTrackId(grandSecondary), secondary);

  // direct and indirect descendants of the second primary
  BOOST_CHECK(st.isTrackDaughterOf(secondary, prim1));
  BOOST_CHECK(st.isTrackDaughterOf(grandSecondary, prim1));
  BOOST_CHECK(st.isTrackDaughterOf(grandSecondary, secondary));

  // nothing here descends from the first primary; the search has to terminate
  BOOST_CHECK(!st.isTrackDaughterOf(prim1, prim0));
  BOOST_CHECK(!st.isTrackDaughterOf(secondary, prim0));
  BOOST_CHECK(!st.isTrackDaughterOf(grandSecondary, prim0));
}
