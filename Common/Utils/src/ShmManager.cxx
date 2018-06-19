/*
 * ShmManager.cxx
 *
 *  Created on: Jun 17, 2018
 *      Author: swenzel
 */

#include "CommonUtils/ShmManager.h"
#include <FairLogger.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include <algorithm>

namespace o2 {
namespace utils {

ShmManager::ShmManager()
{
  if ((mShmID = shmget(IPC_PRIVATE, SHMPOOLSIZE, IPC_CREAT | 0666)) == -1) {
    perror("shmget: shmget failed");
  } else {
    // we will put a secret message in form of an array
    auto addr = shmat(mShmID, nullptr, 0);

    // insert the initial free block
    mFreeBlocks.push_back({ addr, SHMPOOLSIZE });
    mMappedPtr = addr;
  }
  LOG(INFO) << "SHARED MEM INITIALIZED AT ID " << mShmID;
}

ShmManager::~ShmManager()
{
  LOG(INFO) << "REMOVING SHARED MEM SEGMENT ID" << mShmID;
  shmctl(mShmID, IPC_RMID, nullptr);
}

void ShmManager::printBlocks(std::list<MemBlock> const& blocks) const
{
  for (auto& block : blocks) {
    LOG(INFO) << " BLOCKSTART " << block.startptr << " SIZE " << block.bytes;
  }
}

void ShmManager::printAllocedBlocks() const
{
  LOG(INFO) << "=== ALLOCATED BLOCKS ====";
  printBlocks(mAllocedBlocks);
  LOG(INFO) << "-------------------------";
}

void ShmManager::printFreeBlocks() const
{
  LOG(INFO) << "=== FREE BLOCKS ======";
  printBlocks(mFreeBlocks);
  LOG(INFO) << "----------------------";
}

void ShmManager::printSummary() const
{
  printAllocedBlocks();
  printFreeBlocks();
}

// This implements a very rudimentary malloc/free mechanism

void* ShmManager::getmemblock(size_t size)
{
  if (size == 0) {
    LOG(WARNING) << "someone asked for zero size block";
  }
  // iterate through free blocks and see if we find one that is large enough
  // for the requested size ... then split it
  auto iter = std::find_if(mFreeBlocks.begin(), mFreeBlocks.end(),
                           [size](MemBlock const& mem) { return size <= mem.bytes; });
  if (iter == mAllocedBlocks.end()) {
    LOG(INFO) << "DID NOT FIND LARGE ENOUGH MEMORY BLOCK";
    return nullptr;
  }
  auto addr = iter->startptr;
  // add taken block to allocated blocks
  mAllocedBlocks.push_back({ addr, size });

  // change free block info to new information
  iter->startptr = (char*)iter->startptr + size;
  iter->bytes = iter->bytes - size;
  // remove zero free blocks
  if (iter->bytes == 0) {
    mFreeBlocks.erase(iter);
  }
  return addr;
}

// brute force merge sweep
// TODO: be more clever by directly using iterator input
// returns true of merge was done;
bool ShmManager::mergeFreeBlocks()
{
  auto iter = mFreeBlocks.begin();
  auto previous = iter;
  for (; iter != mFreeBlocks.end(); ++iter) {
    if (iter == mFreeBlocks.begin()) {
      continue;
    }
    if (iter->startptr == (char*)previous->startptr + previous->bytes) {
      previous->bytes += iter->bytes;
      mFreeBlocks.erase(iter);
      return true;
    }
    previous = iter;
  }
  return false;
}

void ShmManager::freememblock(void* ptr)
{
  // look for the ptr and add this block to the list of free blocks
  auto iter = std::find_if(mAllocedBlocks.begin(), mAllocedBlocks.end(), [ptr](MemBlock const& mem) { return mem.startptr == ptr; });
  if (iter == mAllocedBlocks.end()) {
    LOG(INFO) << "FREE CORRUPTED (did not find address: " << ptr << " was pointer ok? " << isPointerOk(ptr);
  } else {
    MemBlock memblock{ iter->startptr, iter->bytes };
    if (iter->bytes > 0) {
      auto cmp = [](MemBlock const& block1, MemBlock const& block2) { return block1.startptr < block2.startptr; };
      auto insertpos = std::upper_bound(mFreeBlocks.begin(), mFreeBlocks.end(), memblock, cmp);
      mFreeBlocks.insert(insertpos, memblock);
    }
    mAllocedBlocks.erase(iter);

    while (mergeFreeBlocks()) {
    }
  }
}

void ShmManager::release()
{
  LOG(INFO) << "REMOVING SHARED MEM SEGMENT ID" << mShmID;
  shmctl(mShmID, IPC_RMID, nullptr);
  mShmID = -1;
  mAllocedBlocks.clear();
  mFreeBlocks.clear();
}

} // end namespace utils
} // end namespace o2

