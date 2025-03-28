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

/// \file GPUDefParametersLoad.inc
/// \author David Rohr

#include "GPUDefParameters.h"
#include "GPUDefMacros.h"
#include <sstream>
#include <string>

namespace o2::gpu::internal
{

static GPUDefParameters GPUDefParametersLoad()
{
  return GPUDefParameters{
    // clang-format off
    {$<JOIN:$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>,REPLACE,[^A-Za-z0-9]+,_>,PREPEND,GPUCA_M_FIRST(GPUCA_LB_>,APPEND,)>,$<COMMA>>},
    {$<JOIN:$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>,REPLACE,[^A-Za-z0-9]+,_>,PREPEND,GPUCA_M_FIRST(GPUCA_M_SHIFT(GPUCA_LB_>,APPEND,$<COMMA>0))>,$<COMMA>>},
    {$<JOIN:$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>,REPLACE,[^A-Za-z0-9]+,_>,PREPEND,GPUCA_M_FIRST(GPUCA_M_SHIFT(GPUCA_M_SHIFT(GPUCA_LB_>,APPEND,$<COMMA>0$<COMMA>0)))>,$<COMMA>>}
    // clang-format on
  };
}

#define GPUCA_EXPORT_KERNEL(name)                                                                                  \
  if (par.par_LB_maxThreads[i] > 0) {                                                                              \
    o << "#define GPUCA_" << (forRTC ? "RTC_" : "") << "LB_" << GPUCA_M_STR(name) " " << par.par_LB_maxThreads[i]; \
    if (par.par_LB_minBlocks[i] > 0) {                                                                             \
      o << ", " << par.par_LB_minBlocks[i];                                                                        \
    }                                                                                                              \
    if (!forRTC && par.par_LB_forceBlocks[i] > 0) {                                                                \
      o << ", " << par.par_LB_forceBlocks[i];                                                                      \
    }                                                                                                              \
    o << "\n";                                                                                                     \
  }                                                                                                                \
  i++;

static std::string GPUDefParametersExport(const GPUDefParameters& par, bool forRTC)
{
  std::stringstream o; // clang-format off
  int32_t i = 0;
  $<JOIN:$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<LIST:TRANSFORM,$<TARGET_PROPERTY:O2_GPU_KERNELS,O2_GPU_KERNEL_NAMES>,REPLACE,[^A-Za-z0-9]+,_>,PREPEND,GPUCA_EXPORT_KERNEL(>,APPEND,)>,
  >
  return o.str(); // clang-format on
}

} // namespace o2::gpu::internal
