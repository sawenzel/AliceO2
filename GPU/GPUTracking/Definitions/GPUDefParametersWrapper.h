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

/// \file GPUDefParametersWrapper.h
/// \author David Rohr

// Wrapper file to load all compile-time parameters (architecture / rtc - dependent ones, and constant ones)
// Compile-time constants affecting the tracking algorithms / results are located in GPUDefConstantsAndSettings.h

#ifndef GPUDEFPARAMETERSWRAPPER_H
#define GPUDEFPARAMETERSWRAPPER_H
// clang-format off

#include "GPUCommonDef.h"
#include "GPUDefMacros.h"

#ifndef GPUCA_GPUCODE_GENRTC
#include "GPUDefParametersDefaults.h"
#endif
#include "GPUDefParametersConstants.h"

#ifdef GPUCA_GPUCODE
  #define GPUCA_GET_THREAD_COUNT(...) GPUCA_M_FIRST(__VA_ARGS__)
#else
  #define GPUCA_GET_THREAD_COUNT(...) 1 // On the host, a thread is a block, and we run 1 "device thread" per block.
#endif

#define GPUCA_GET_WARP_COUNT(...) (GPUCA_GET_THREAD_COUNT(__VA_ARGS__) / GPUCA_WARP_SIZE)

#define GPUCA_MERGER_INTERPOLATION_ERROR_TYPE_A GPUCA_DETERMINISTIC_CODE(float, GPUCA_MERGER_INTERPOLATION_ERROR_TYPE)
#define GPUCA_DEDX_STORAGE_TYPE_A GPUCA_DETERMINISTIC_CODE(float, GPUCA_DEDX_STORAGE_TYPE)

// #define GPUCA_TRACKLET_CONSTRUCTOR_DO_PROFILE                       // Output Profiling Data for Tracklet Constructor Tracklet Scheduling

// #define GPUCA_KERNEL_DEBUGGER_OUTPUT

// Some assertions to make sure the parameters are not invalid
#if defined(GPUCA_GPUCODE)
  static_assert(GPUCA_MAXN >= GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP, "Invalid GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP");
  static_assert(GPUCA_ROW_COUNT >= GPUCA_TRACKLET_SELECTOR_HITS_REG_SIZE, "Invalid GPUCA_TRACKLET_SELECTOR_HITS_REG_SIZE");
  static_assert(GPUCA_M_FIRST(GPUCA_LB_GPUTPCCompressionKernels_step1unattached) * 2 <= GPUCA_TPC_COMP_CHUNK_SIZE, "Invalid GPUCA_TPC_COMP_CHUNK_SIZE");
#endif

// Derived parameters
#ifdef GPUCA_USE_TEXTURES
  #define GPUCA_TEXTURE_FETCH_CONSTRUCTOR                              // Fetch data through texture cache
#endif

// clang-format on
#endif // GPUDEFPARAMETERSWRAPPER_H
