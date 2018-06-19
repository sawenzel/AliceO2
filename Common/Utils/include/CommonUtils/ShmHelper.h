/*
 * ShmHelper.h
 *
 *  Created on: Jun 19, 2018
 *      Author: sandro
 */

#ifndef COMMON_UTILS_INCLUDE_COMMONUTILS_SHMHELPER_H_
#define COMMON_UTILS_INCLUDE_COMMONUTILS_SHMHELPER_H_

#include <sys/shm.h>
#include <FairLogger.h>

namespace o2 { namespace utils {

// class creating a shared memory pool
// and manges allocations within the ppol
class ShmHelper {
public:
 static ShmHelper& Instance()
 {
   static ShmHelper instance;
   return instance;
 }

 using address_offset_type = std::pair<void*, int>;

 address_offset_type attachOrGiveBack(int id, void* base_ptr)
 {
   auto iter = mShmTable.find(id);
   if (iter != mShmTable.end()) {
     return iter->second;
   }
   auto addr = shmat(id, base_ptr, 0);
   LOG(INFO) << " SHM ADDRESS " << addr << " VS WANTED " << base_ptr;
   if (addr != base_ptr) {
     LOG(WARNING) << " Trying a second time without constraint ";
     addr = shmat(id, nullptr, 0);
     LOG(INFO) << " SHM ADDRESS " << addr;
   }
   int offset = (int)((char*)addr - (char*)base_ptr);
   auto pair = address_offset_type(addr, offset);
   mShmTable[id] = pair;
   return pair;
 }

private:
  std::map<int, address_offset_type> mShmTable; // map of attached segments
};

}
} // end namespace


#endif /* COMMON_UTILS_INCLUDE_COMMONUTILS_SHMHELPER_H_ */
