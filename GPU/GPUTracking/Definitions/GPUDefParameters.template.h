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

/// \file GPUDefParameters.h
/// \author David Rohr

#ifndef GPUDEFPARAMETERS_H
#define GPUDEFPARAMETERS_H

namespace o2::gpu
{
struct GPUDefParameters {  // clang-format off
  int32_t par_LB_maxThreads[$<LIST:LENGTH,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>>] = {};
  int32_t par_LB_minBlocks[$<LIST:LENGTH,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>>] = {};
  int32_t par_LB_forceBlocks[$<LIST:LENGTH,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>>] = {};
};  // clang-format on
}  // namespace o2::gpu

#endif
