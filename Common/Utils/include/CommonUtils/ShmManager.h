/*
 * ShmManager.h
 *
 *  Created on: Jun 17, 2018
 *      Author: swenzel
 */

#ifndef COMMON_UTILS_INCLUDE_COMMONUTILS_SHMMANAGER_H_
#define COMMON_UTILS_INCLUDE_COMMONUTILS_SHMMANAGER_H_

#include <list>

namespace o2 {
namespace utils {
namespace {
struct MemBlock {
  void* startptr;
  size_t bytes;
};
}

constexpr size_t SHMPOOLSIZE = 1024*1024; // 1MB of shared memory

// class creating a shared memory pool
// and manges allocations within the ppol
class ShmManager {
public:
 static ShmManager& Instance()
 {
   static ShmManager instance;
   return instance;
 }

 // the equivalent of malloc
 void* getmemblock(size_t size);
 // the equivalent of free
 void freememblock(void*);

 void printAllocedBlocks() const;
 void printFreeBlocks() const;
 void printSummary() const;
 void release();
 int getShmID() const { return mShmID; }
 size_t getPointerOffset(void*ptr) const { return (size_t)((char*)ptr - (char*)mMappedPtr); }

 size_t getNumAllocedBlocks() const { return mAllocedBlocks.size(); }
 size_t getNumFreeBlocks() const { return mFreeBlocks.size(); }
 size_t getManagedSize() const { return SHMPOOLSIZE; }
 size_t getAllocedSize() const
 {
   size_t accum{ 0 };
   for (auto& b : mAllocedBlocks) {
     accum += b.bytes;
   }
   return accum;
 }
 size_t getFreeSize() const
 {
   size_t accum{ 0 };
   for (auto& b : mFreeBlocks) {
     accum += b.bytes;
   }
   return accum;
 }

private:
 ShmManager();
 ~ShmManager();
 int mShmID = -1; // id of shared mem created
 void* mMappedPtr = nullptr; // the mapped ptr of the segment
 // sweep over memory and merge free blocks if they "touch"
 bool mergeFreeBlocks();
 void printBlocks(std::list<MemBlock> const& block) const;

 std::list<MemBlock> mAllocedBlocks; // allocedblocks
 std::list<MemBlock> mFreeBlocks;    // free blocks
};

}
} // end namespace

#endif /* COMMON_UTILS_INCLUDE_COMMONUTILS_SHMMANAGER_H_ */
