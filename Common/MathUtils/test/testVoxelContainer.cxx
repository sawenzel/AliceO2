// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/*
 *  Created on: Feb 9, 2019
 *      Author: swenzel
 */

#define BOOST_TEST_MODULE Test VoxelContainer
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "MathUtils/VoxelContainer.h"

using namespace o2;

BOOST_AUTO_TEST_CASE(VoxelContainer_test1)
{
  VoxelContainer<int> cont;
  int NX = 10;
  int NY = 10;
  int NZ = 10;
  double Lz = 100;
  double Lx = 20;
  double Ly = 20;
  cont.init(NX, NY, NZ, Lx, Ly, Lz);
  BOOST_TEST(cont.getBinX(-Lx/2.) == 0);
  BOOST_TEST(cont.getBinY(-Ly/2.) == 0);
  BOOST_TEST(cont.getBinZ(-Lz/2.) == 0);
  BOOST_TEST(cont.getBinX(Lx/2.) == NX);
  BOOST_TEST(cont.getBinY(Ly/2.) == NY);
  BOOST_TEST(cont.getBinZ(Lz/2.) == NZ);

  BOOST_TEST(cont.getVoxelContent(0.,0.,0.) == 0);
  cont.getVoxelContent(0.,0.,0.) = 1;
  BOOST_TEST(cont.getVoxelContent(0.,0.,0.) == 1);

  BOOST_TEST(cont.getVoxelContentByBin(0, 0, 0) == 0);
}
