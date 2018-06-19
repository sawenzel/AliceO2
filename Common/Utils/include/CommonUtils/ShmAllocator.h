/*
 * ShmAllocator.h
 *
 *  Created on: Jun 18, 2018
 *      Author: swenzel
 */

#ifndef COMMON_UTILS_INCLUDE_COMMONUTILS_SHMALLOCATOR_H_
#define COMMON_UTILS_INCLUDE_COMMONUTILS_SHMALLOCATOR_H_

#include "CommonUtils/ShmManager.h"
#include <FairLogger.h>

namespace o2 {
namespace utils {

// an allocator placing objects in shared memory managed
// by ShmManager
template <typename T>
class ShmAllocator
{
 public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  typedef T* pointer;
  typedef const T* const_pointer;

  typedef T& reference;
  typedef const T& const_reference;

 public:
  inline ShmAllocator() throw() {}

  template <typename T2>
  inline ShmAllocator(const ShmAllocator<T2>&) throw()
  {
  }

  inline ~ShmAllocator() throw() {}

  inline pointer adress(reference r) { return &r; }

  inline const_pointer adress(const_reference r) const { return &r; }

  // the actually important functions:
  inline pointer allocate(size_type n) { return (pointer)ShmManager::Instance().getmemblock(sizeof(value_type) * n); }
  inline void deallocate(pointer p, size_type) { ShmManager::Instance().freememblock(p); }

  inline void construct(pointer p, const value_type& value) {
	//LOG(FATAL) << " oka="
	  new (p) value_type(value);
  }

  template <class U, class... Args>
  void construct(U* p, Args&&... args)
  {
	// LOG(INFO) << "IS POINTER OK ? : " << ShmManager::Instance().isPointerOk((void*)p);
    ::new ((void*)p) U(std::forward<Args>(args)...);
  }

  inline void destroy(pointer p) { p->~value_type(); }

  inline size_type max_size() const throw() { return size_type(-1) / sizeof(value_type); }

  template <typename T2>
  struct rebind {
    typedef ShmAllocator<T2> other;
  };

  bool operator!=(const ShmAllocator<T>& other) const { return !(*this == other); }

  // Returns true if and only if storage allocated from *this
  // can be deallocated from other, and vice versa.
  // Always returns true for stateless allocators.
  bool operator==(const ShmAllocator<T>& /*other*/) const { return true; }
};

//  template <class T, std::size_t = sizeof(T)>
//  std::true_type is_complete_impl(T*);
//  std::false_type is_complete_impl(...);
//  template <class T>
//  using is_complete = decltype(is_complete_impl(std::declval<T*>()));

template <typename T>
std::vector<T>* createSimVector()
{
  using vector_t = std::vector<T>;
#ifdef USESHM
  auto& instance = o2::utils::ShmManager::Instance();
  auto placement = instance.getmemblock(sizeof(vector_t));

  // at this moment we have to trust that std::
  return new (placement) vector_t;
#else
  return new vector_t;
#endif
}

template <typename T>
void freeSimVector(std::vector<T>* ptr) {
  using vector_t = std::vector<T>;
#ifdef USESHM
  auto& instance = o2::utils::ShmManager::Instance();
  ptr->clear();
  ptr->shrink_to_fit();
  instance.freememblock(ptr);
  // at this moment we have to trust that std::
#else
  delete ptr;
#endif
}
}
}

#endif /* COMMON_UTILS_INCLUDE_COMMONUTILS_SHMALLOCATOR_H_ */
