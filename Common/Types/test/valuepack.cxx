// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   valuepack.cxx
/// @author Matthias Richter
/// @since  2017-09-28
/// @brief  Unit test for ValuePack type

#define BOOST_TEST_MODULE Test Type Value Pack
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <bitset>
#include "../include/Types/ValuePack.h"

using TwoFourNinePack = o2::types::ValuePack<unsigned int, 2, 4, 9>;

union Foo {
  unsigned int value;
private:
  struct {
    unsigned int x:2;
    unsigned int y:4;
    unsigned int z:9;
  };
public:
  unsigned int getX() const { return x; }
  unsigned int getY() const { return y; }
  unsigned int getZ() const { return z; }
};

__attribute__((noinline))
unsigned int Get2(Foo const &f) {
  return f.getY();
}

__attribute__((noinline))
unsigned int Get2(TwoFourNinePack const &f) {
  return f.get<1,unsigned int>();
}

BOOST_AUTO_TEST_CASE(test_valuepack)
{
  std::cout << TwoFourNinePack::nbits << " "
            << TwoFourNinePack::nfields << " "
            << TwoFourNinePack::size << " "
            << std::endl;

  BOOST_CHECK(TwoFourNinePack::nbits == 15);
  BOOST_CHECK(TwoFourNinePack::nfields == 3);

  TwoFourNinePack pack(3, 6, 1022);
  std::cout << "pack from values 3, 6, 1022:" << std::endl;

  std::cout << "0x" << std::hex << pack
            << ": " << std::bitset<TwoFourNinePack::nbits>(pack) << std::endl;

  BOOST_CHECK((pack.get<0, unsigned int>()) == 3);
  BOOST_CHECK((pack.get<1, unsigned int>()) == 6);
  BOOST_CHECK((pack.get<2, unsigned int>()) == 510);

  pack.set<1>(0xa);
  std::cout << "0x" << std::hex << pack
            << ": " << std::bitset<TwoFourNinePack::nbits>(pack) << std::endl;
  BOOST_CHECK((pack.get<1, unsigned int>()) == 10);

  // check consistency with union
  Foo f; f.value = pack;
  BOOST_CHECK(Get2(f) == Get2(pack));
  
  TwoFourNinePack copy = pack;
}

