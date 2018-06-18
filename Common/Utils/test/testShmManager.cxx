// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test TreeShmManager
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "CommonUtils/ShmManager.h"
#include "CommonUtils/ShmAllocator.h"
#include <FairLogger.h>

using namespace o2::utils;

// testing ShmManager
BOOST_AUTO_TEST_CASE(TreeShmManager_test1)
{
  auto& instance = ShmManager::Instance();
  instance.printFreeBlocks();

  // VERY EASY SCENARIO: malloc an array and remove it
  {
    double* array = (double*)instance.getmemblock(sizeof(double) * 100);
    instance.printSummary();

    BOOST_CHECK(instance.getNumAllocedBlocks() == 1);
    BOOST_CHECK(instance.getNumFreeBlocks() == 1);
    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());


    // free this allocation
    instance.freememblock(array);
    instance.printSummary();
    BOOST_CHECK(instance.getNumAllocedBlocks() == 0);
    BOOST_CHECK(instance.getNumFreeBlocks() == 1);
    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());
  }

  // SECOND SCENARIO: malloc 2 arrays and remove first
  {
    double* array1 = (double*)instance.getmemblock(sizeof(double) * 100);
    double* array2 = (double*)instance.getmemblock(sizeof(double) * 50);
    instance.printSummary();

    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());

    BOOST_CHECK(instance.getNumAllocedBlocks() == 2);
    BOOST_CHECK(instance.getNumFreeBlocks() == 1);

    // free one allocation
    instance.freememblock(array1);
    instance.printSummary();
    BOOST_CHECK(instance.getNumAllocedBlocks() == 1);
    // this should be now split
    BOOST_CHECK(instance.getNumFreeBlocks() == 2);

    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());

    instance.freememblock(array2);
    BOOST_CHECK(instance.getNumAllocedBlocks() == 0);
    BOOST_CHECK(instance.getNumFreeBlocks() == 1);

    instance.printSummary();
  }


  // THIRD SCENARIO: malloc 3 arrays and remove first 2
  {
    double* array1 = (double*)instance.getmemblock(sizeof(double) * 100);
    double* array2 = (double*)instance.getmemblock(sizeof(double) * 50);
    double* array3 = (double*)instance.getmemblock(sizeof(double) * 75);
    instance.printSummary();

    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());

    BOOST_CHECK(instance.getNumAllocedBlocks() == 3);
    BOOST_CHECK(instance.getNumFreeBlocks() == 1);

    // free one allocation
    instance.freememblock(array1);
    instance.freememblock(array2);
    instance.printSummary();
    BOOST_CHECK(instance.getNumAllocedBlocks() == 1);
    // this should be now split
    BOOST_CHECK(instance.getNumFreeBlocks() == 2);

    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());

    // remove remaining one
    instance.freememblock(array3);
    BOOST_CHECK(instance.getNumAllocedBlocks() == 0);
    // this should be now split
    BOOST_CHECK(instance.getNumFreeBlocks() == 1);
    BOOST_CHECK(instance.getManagedSize() == instance.getFreeSize() + instance.getAllocedSize());
  }
}


// testing ShmAllocator
BOOST_AUTO_TEST_CASE(TreeShmManager_test2)
{
  auto& instance = ShmManager::Instance();
  {
    std::vector<double, o2::utils::ShmAllocator<double>> v;
    v.emplace_back(10);
    instance.printSummary();
    BOOST_CHECK(instance.getNumAllocedBlocks() == 1);
  }
  // here the vector is out of scope and should be freed
  BOOST_CHECK(instance.getNumAllocedBlocks() == 0);
}

// testing ShmAllocator: place the whole vector in shared mem
BOOST_AUTO_TEST_CASE(TreeShmManager_test3)
{
  auto& instance = ShmManager::Instance();
  {
    using vector_t = std::vector<double, o2::utils::ShmAllocator<double>>;
    auto placement = instance.getmemblock(sizeof(vector_t));
    vector_t* vptr = new (placement) vector_t;
    vptr->emplace_back(10);
    instance.printSummary();

    // we should have 2 alloced blocks
    BOOST_CHECK(instance.getNumAllocedBlocks() == 2);

    // COOL: So this thing is entirely in shared memory

    // complicated removal; how to do it better?
    vptr->clear();
    vptr->shrink_to_fit();
    instance.freememblock(placement);
  }
  // here the vector is out of scope and should be freed
  BOOST_CHECK(instance.getNumAllocedBlocks() == 0);
}


class Foo {
public:
	  Foo(double x) : mX{x} {}
	  double mX = 0;
};

// this tells the system to put all Foo objects in shared mem (if allocator is used)
namespace std
{
template <>
class allocator<Foo> : public o2::utils::ShmAllocator<Foo>
{
};
}

// check std::allocator specialization
BOOST_AUTO_TEST_CASE(TreeShmManager_test4)
{
  LOG(INFO) << "Starting test4 with specialized std::allocator";
  auto& instance = ShmManager::Instance();
  {
    using vector_t = std::vector<Foo>;
    auto placement = instance.getmemblock(sizeof(vector_t));
    vector_t* vptr = new (placement) vector_t;
    vptr->emplace_back(10);
    instance.printSummary();

    // we should have 2 alloced blocks
    BOOST_CHECK(instance.getNumAllocedBlocks() == 2);

    // COOL: So this thing is entirely in shared memory

    // complicated removal; how to do it better?
    vptr->clear();
    vptr->shrink_to_fit();
    instance.freememblock(placement);
  }
  // here the vector is out of scope and should be freed
  BOOST_CHECK(instance.getNumAllocedBlocks() == 0);
}


