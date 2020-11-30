// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_ROOTSERKEYVALUESTORE_H
#define ALICEO2_ROOTSERKEYVALYESTORE_H

#include <map>
#include <string>
#include <TClass.h>
#include <Rtypes.h>
#include <typeinfo>
#include <iostream>

namespace o2
{
namespace utils
{

/// Put doc here
///
class RootSerializableKeyValueStore
{
 public:
  RootSerializableKeyValueStore() = default;

  template <typename T>
  void put(std::string const& key, T value);

  /// a pointer would be better for indicating errors
  template <typename T>
  T get(std::string const& key) const;

 private:
  std::map<std::string, std::pair<char*, std::type_info const*>> mStore;

  ClassDefNV(RootSerializableKeyValueStore, 1);
};

template <typename T>
inline void RootSerializableKeyValueStore::put(std::string const& key, T value)
{
  auto ptr = new T(value);
  mStore.insert(std::pair<std::string, std::pair<char*, std::type_info const*>>(key, std::pair<char*, std::type_info const*>((char*)ptr, &typeid(value))));
}

template <typename T>
inline T RootSerializableKeyValueStore::get(std::string const& key) const
{
  auto iter = mStore.find(key);
  if (iter != mStore.end()) {
    if (typeid(T) == *(iter->second.second)) {
      T rv(*(T*)(iter->second.first));
      return rv;
    } else {
      std::cerr << "typeid does not match \n";
      return T(0);
    }
  }
  std::cerr << "no such key " << key << "\n";
  return T(0);
}

} // namespace utils
} // namespace o2

#endif
