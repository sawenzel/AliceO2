#pragma once
#include "TClonesArray.h"
#include "TObject.h"
#include <type_traits>
#include <cassert>

namespace AliceO2 {

// A smart and type-safe wrapper around TClonesArray
template <typename T> class TypedTClonesArray {
  // need to impose that T is a pointer type (as is the case for TClonesArray)
  static_assert(std::is_pointer<T>::value,
                "template argument of TypedTClonesArray needs to be pointer");

  // imposing that T is derived from ROOT TObject class
  static_assert(std::is_base_of<TObject, typename std::remove_pointer<T>::type>::value,
                "template argument of TypedTClonesArray needs to inherit from TObject");

  using nopointerT = typename std::remove_pointer<T>::type;

  // internal iterator class
  // http://stackoverflow.com/questions/7758580/writing-your-own-stl-container/7759622#7759622
  class iterator {
  public:
    typedef T value_type;
    typedef T &reference;
    typedef T pointer; /* this is because T is itself a pointer */
    iterator(size_t index, TClonesArray *array)
        : mIndex(index), mArray(array) {}
    iterator(const iterator &other)
        : mIndex(other.mIndex), mArray(other.mArray) {}
    ~iterator(){/* nothing to be done */};

    iterator &operator=(const iterator &other) {
      mIndex = other.mIndex;
      mArray = other.mArray;
    }
    bool operator==(const iterator &other) const {
      return mArray == other.mArray && mIndex == other.mIndex;
    }
    bool operator!=(const iterator &other) const {
      return mArray != other.mArray || mIndex != other.mIndex;
    }

    // neglecting empty slots
    iterator &operator++() {
      do {
        mIndex++;
      } while ( mArray->UncheckedAt(mIndex) == nullptr && mIndex < mArray->GetLast()-1 );
      return *this;
    }

    // dereferencing the iterator:
    pointer operator*() const {
      return static_cast<T>(mArray->UncheckedAt(mIndex));
    }
    // reference operator->() const;

  private:
    // the state of the iterator
    size_t mIndex;
    TClonesArray *mArray;
  };

public:
  // implicit constructor
  TypedTClonesArray(TClonesArray &p) : mBare(&p) {
    assert(nopointerT::Class() == p.GetClass());
  }
  TypedTClonesArray(TClonesArray *ptr) : mBare(ptr) {
    assert(nopointerT::Class() == ptr->GetClass());
  }

  // explicit construction (so that user never has to say TClonesArray)
  TypedTClonesArray(size_t size) : mBare(new TClonesArray(nopointerT::Class())) {

  }

  // assignment operator from a TClonesArray pointer
  TypedTClonesArray operator=(TClonesArray *ptr) {
    assert(nopointerT::Class() == ptr->GetClass());
    mBare = ptr;
    return *this;
  }

  TypedTClonesArray(TypedTClonesArray const &other) : mBare(other.mBare) {}

  iterator begin() const { return iterator(0, mBare); }
  iterator end() const { return iterator(mBare->GetLast() + 1, mBare); }

  // A type safe operator[] to read elements
  // We could make a fast version of this (without checks realocation etc.)
  // by directly accessing the internal data of mBare
  T operator[](size_t i) const {
    return static_cast<T>(mBare->UncheckedAt(i));
  }

  // we delete this function: It is evil !!
  // This also means we will dissallow placemenent news in favour of ConstructedAt
  T& operator[](size_t) = delete;

  // construct object in place at slot i or return existing object
  T ConstructedAt(size_t i) {
    return static_cast<T>(mBare->ConstructedAt(i));
  }

  // returns if slot i in the array has been constructed with an object
  bool IsSlotConstructed(size_t i) const {
    return this->operator [](i) != nullptr;
  }

  // we can provide here the better and fast sorting:
  // void Sort();

  // access to the raw container
  TClonesArray *raw() { return mBare; }
  TClonesArray const *raw() const { return mBare; }

  // should we allow implicit castings to TClonesArray (or prefer raw())??
  operator TClonesArray*() const { return mBare; }

  // Ideas: offer reasonable push_back, emplace_back

  void Delete() {
    mBare->Delete();
  }

private:
  TClonesArray *mBare;
};


} // end namespace
