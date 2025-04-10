// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUDefParametersConstants.h
/// \author David Rohr

// This file contains compile-time constants, independent from the backend

#ifndef GPUDEFPARAMETERSCONSTANTS_H
#define GPUDEFPARAMETERSCONSTANTS_H
// clang-format off

#define GPUCA_THREAD_COUNT_SCAN 512 // TODO: WARNING!!! Must not be GPUTYPE-dependent right now! // TODO: Fix!

#if defined(__CUDACC__) || defined(__HIPCC__)
  #define GPUCA_SPECIALIZE_THRUST_SORTS
#endif

#define GPUCA_MAX_THREADS 1024
#define GPUCA_MAX_STREAMS 36

#if defined(GPUCA_GPUCODE)
  #define GPUCA_SORT_STARTHITS                                         // Sort the start hits when running on GPU
#endif

#define GPUCA_ROWALIGNMENT 16                                          // Align of Row Hits and Grid
#define GPUCA_BUFFER_ALIGNMENT 64                                      // Alignment of buffers obtained from SetPointers
#define GPUCA_MEMALIGN (64 * 1024)                                     // Alignment of allocated memory blocks

// Default maximum numbers
#define GPUCA_MAX_CLUSTERS           ((size_t)     1024 * 1024 * 1024) // Maximum number of TPC clusters
#define GPUCA_MAX_TRD_TRACKLETS      ((size_t)             128 * 1024) // Maximum number of TRD tracklets
#define GPUCA_MAX_ITS_FIT_TRACKS     ((size_t)              96 * 1024) // Max number of tracks for ITS track fit
#define GPUCA_MEMORY_SIZE            ((size_t) 6 * 1024 * 1024 * 1024) // Size of memory allocated on Device
#define GPUCA_HOST_MEMORY_SIZE       ((size_t) 1 * 1024 * 1024 * 1024) // Size of memory allocated on Host
#define GPUCA_GPU_STACK_SIZE         ((size_t)               8 * 1024) // Stack size per GPU thread
#define GPUCA_GPU_HEAP_SIZE          ((size_t)       16 * 1025 * 1024) // Stack size per GPU thread

#ifdef GPUCA_GPUCODE
  #ifndef GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP
     #define GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP 6
  #endif
  #ifndef GPUCA_TRACKLET_SELECTOR_HITS_REG_SIZE
     #define GPUCA_TRACKLET_SELECTOR_HITS_REG_SIZE 12
  #endif
  #ifndef GPUCA_ALTERNATE_BORDER_SORT
     #define GPUCA_ALTERNATE_BORDER_SORT 0
  #endif
  #ifndef GPUCA_SORT_BEFORE_FIT
     #define GPUCA_SORT_BEFORE_FIT 0
  #endif
  #ifndef GPUCA_MERGER_SPLIT_LOOP_INTERPOLATION
     #define GPUCA_MERGER_SPLIT_LOOP_INTERPOLATION 0
  #endif
  #ifndef GPUCA_COMP_GATHER_KERNEL
     #define GPUCA_COMP_GATHER_KERNEL 0
  #endif
  #ifndef GPUCA_COMP_GATHER_MODE
     #define GPUCA_COMP_GATHER_MODE 2
  #endif
#else
  #define GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP 0
  #define GPUCA_TRACKLET_SELECTOR_HITS_REG_SIZE 0
  #define GPUCA_ALTERNATE_BORDER_SORT 0
  #define GPUCA_SORT_BEFORE_FIT 0
  #define GPUCA_MERGER_SPLIT_LOOP_INTERPOLATION 0
  #define GPUCA_THREAD_COUNT_FINDER 1
  #define GPUCA_COMP_GATHER_KERNEL 0
  #define GPUCA_COMP_GATHER_MODE 0
#endif
#ifndef GPUCA_DEDX_STORAGE_TYPE
  #define GPUCA_DEDX_STORAGE_TYPE float
#endif
#ifndef GPUCA_MERGER_INTERPOLATION_ERROR_TYPE
  #define GPUCA_MERGER_INTERPOLATION_ERROR_TYPE float
#endif

// clang-format on
#endif // GPUDEFPARAMETERSCONSTANTS_H
