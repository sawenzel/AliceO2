// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_DATA_PRIMARYCHUNK_H_
#define ALICEO2_DATA_PRIMARYCHUNK_H_

#include <cstring>

namespace o2 {
namespace Data {

struct SubEventInfo {
  float eventtime = 0.;
  uint32_t eventID = 0; // which event ID
  int32_t  maxEvents = -1; // the number of events in this run (if known otherwise set to -1)
  int32_t  runID = 0; // the runID of this run
  uint16_t part = 0;    // which part of the eventID
  uint16_t nparts = 0; // out of how many parts
  uint32_t seed = 0; // seed for RNG
  uint32_t index = 0;
  // might add more fields (such as which process treated this chunk etc)
  ClassDefNV(SubEventInfo, 1);
};

 inline
 bool operator<(SubEventInfo const &a, SubEventInfo const &b) {
   return a.eventID <= b.eventID && (a.part < b.part);
 }
 
 template <typename V, typename T>
 size_t getNecessaryElements() {
  using VT = typename V::value_type;  
  const auto ints = sizeof(T) / sizeof(VT);
  const auto rest = sizeof(T) % sizeof(VT);
  return (rest == 0) ? ints : ints + 1;
 }
 
// utility function that goes somewhere else? 
template <typename V, typename T>
void addTrailingObjectToVector(V &v, T const &t) {
  // function is attaching a T object at the end of a vector V holding types V::value_type
  // this should work by expanding the vector by a few elements of enough size
  // and then using the space of those elements to hold one T
  const auto oldsize = v.size();
  for(auto i = 0; i< getNecessaryElements<V,T>(); i++) {
    v.emplace_back(); // emplace_back default constructed objects
  }
  // now copy t to vector
  // for the moment only works in cases the alignment matches
  static_assert(alignof(typename V::value_type) % alignof(T) == 0, "alignement has to match");
  memcpy((void*)&v[oldsize],(const void*)&t, sizeof(t));
}

// inverse function to unpack
// utility function that goes somewhere else? 
template <typename V, typename T>
void removeTrailingObjectFromVector(V &v, T &t) {
  const auto oldsize = v.size();
  const auto tpos = oldsize - getNecessaryElements<V,T>();
  memcpy(&t, (const void*)&v[tpos], sizeof(t));
  v.resize(tpos); 
}

// gets the trailing object and the number of elements that
// can be removed at the end of the vector
template <typename V, typename T>
size_t getTrailingObjectFromVector(const V &v, T &t) {
  const auto oldsize = v.size();
  const auto tpos = oldsize - getNecessaryElements<V,T>();
  memcpy(&t, (const void*)&v[tpos], sizeof(t));
  return (oldsize - tpos);
}


// Encapsulating primaries/tracks
// to be processed by the simulation processors.
// May be part of an event or multiple events. Contains
// seed information
struct PrimaryChunk {
  std::vector<SubEventInfo> mEventIDs;
  std::vector<TParticle> mParticles;
  ClassDefNV(PrimaryChunk, 1);
};

}
}

#endif

