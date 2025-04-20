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

/// \file GPUTPCNNClusterizerHost.cxx
/// \author Christian Sonnabend

#include <CommonUtils/StringUtils.h>

#include "GPUTPCNNClusterizerHost.h"
#include "GPUTPCNNClusterizer.h"
#include "GPUSettings.h"
#include "ML/3rdparty/GPUORTFloat16.h"
#include "GPUReconstruction.h"

#ifdef GPUCA_HAS_ONNX
#include <onnxruntime_cxx_api.h>
#endif

using namespace o2::gpu;

void GPUTPCNNClusterizerHost::init(const GPUSettingsProcessingNNclusterizer& settings)
{
  std::string class_model_path = settings.nnClassificationPath, reg_model_path = settings.nnRegressionPath;
  std::vector<std::string> reg_model_paths;
  std::vector<std::string> evalMode = o2::utils::Str::tokenize(settings.nnEvalMode, ':');

  if (settings.nnLoadFromCCDB) {
    reg_model_path = settings.nnLocalFolder + "/net_regression_c1.onnx"; // Needs to be set identical to NeuralNetworkClusterizer.cxx, otherwise the networks might be loaded from the wrong place
    if (evalMode[0] == "c1") {
      class_model_path = settings.nnLocalFolder + "/net_classification_c1.onnx";
    } else if (evalMode[0] == "c2") {
      class_model_path = settings.nnLocalFolder + "/net_classification_c2.onnx";
    }

    if (evalMode[1] == "r2") {
      reg_model_path += ":" + settings.nnLocalFolder + "/net_regression_c2.onnx";
    }
  }

  OrtOptions = {
    {"model-path", class_model_path},
    {"device-type", settings.nnInferenceDevice},
    {"allocate-device-memory", std::to_string(settings.nnInferenceAllocateDevMem)},
    {"intra-op-num-threads", std::to_string(settings.nnInferenceIntraOpNumThreads)},
    {"inter-op-num-threads", std::to_string(settings.nnInferenceInterOpNumThreads)},
    {"enable-optimizations", std::to_string(settings.nnInferenceEnableOrtOptimization)},
    {"enable-profiling", std::to_string(settings.nnInferenceOrtProfiling)},
    {"profiling-output-path", settings.nnInferenceOrtProfilingPath},
    {"logging-level", std::to_string(settings.nnInferenceVerbosity)},
    {"onnx-environment-name", "c1"}};

  model_class.initOptions(OrtOptions);
  modelsUsed[0] = true;

  reg_model_paths = o2::utils::Str::tokenize(reg_model_path, ':');

  if (!settings.nnClusterizerUseCfRegression) {
    if (reg_model_paths.size() == 1) {
      OrtOptions["model-path"] = reg_model_paths[0];
      OrtOptions["onnx-environment-name"] = "r1";
      model_reg_1.initOptions(OrtOptions);
      modelsUsed[1] = true;
    } else {
      OrtOptions["model-path"] = reg_model_paths[0];
      OrtOptions["onnx-environment-name"] = "r1";
      model_reg_1.initOptions(OrtOptions);
      modelsUsed[1] = true;
      OrtOptions["model-path"] = reg_model_paths[1];
      OrtOptions["onnx-environment-name"] = "r2";
      model_reg_2.initOptions(OrtOptions);
      modelsUsed[2] = true;
    }
  }
}

void GPUTPCNNClusterizerHost::initClusterizer(const GPUSettingsProcessingNNclusterizer& settings, GPUTPCNNClusterizer& clustererNN)
{
  clustererNN.nnClusterizerUseCfRegression = settings.nnClusterizerUseCfRegression;
  clustererNN.nnClusterizerSizeInputRow = settings.nnClusterizerSizeInputRow;
  clustererNN.nnClusterizerSizeInputPad = settings.nnClusterizerSizeInputPad;
  clustererNN.nnClusterizerSizeInputTime = settings.nnClusterizerSizeInputTime;
  clustererNN.nnClusterizerAddIndexData = settings.nnClusterizerAddIndexData;
  clustererNN.nnClusterizerElementSize = ((2 * settings.nnClusterizerSizeInputRow + 1) * (2 * settings.nnClusterizerSizeInputPad + 1) * (2 * settings.nnClusterizerSizeInputTime + 1)) + (settings.nnClusterizerAddIndexData ? 3 : 0);
  clustererNN.nnClusterizerBatchedMode = settings.nnClusterizerBatchedMode;
  clustererNN.nnClusterizerBoundaryFillValue = settings.nnClusterizerBoundaryFillValue;
  clustererNN.nnSigmoidTrafoClassThreshold = settings.nnSigmoidTrafoClassThreshold;
  if (clustererNN.nnSigmoidTrafoClassThreshold) {
    clustererNN.nnClassThreshold = (float)std::log(settings.nnClassThreshold / (1.f - settings.nnClassThreshold));
  } else {
    clustererNN.nnClassThreshold = settings.nnClassThreshold;
  }
  if (settings.nnClusterizerVerbosity < 0) {
    clustererNN.nnClusterizerVerbosity = settings.nnInferenceVerbosity;
  } else {
    clustererNN.nnClusterizerVerbosity = settings.nnClusterizerVerbosity;
  }
  clustererNN.nnInferenceInputDType = settings.nnInferenceInputDType.find("32") != std::string::npos;
  clustererNN.nnInferenceOutputDType = settings.nnInferenceOutputDType.find("32") != std::string::npos;
  clustererNN.nnClusterizerModelClassNumOutputNodes = model_class.getNumOutputNodes()[0][1];
  if (!settings.nnClusterizerUseCfRegression) {
    if (model_class.getNumOutputNodes()[0][1] == 1 || !model_reg_2.isInitialized()) {
      clustererNN.nnClusterizerModelReg1NumOutputNodes = model_reg_1.getNumOutputNodes()[0][1];
    } else {
      clustererNN.nnClusterizerModelReg1NumOutputNodes = model_reg_1.getNumOutputNodes()[0][1];
      clustererNN.nnClusterizerModelReg2NumOutputNodes = model_reg_2.getNumOutputNodes()[0][1];
    }
  }
}

// MockedOrtAllocator implementation to be able to use volatile assignment
struct MockedOrtAllocator : OrtAllocator {
  MockedOrtAllocator(GPUReconstruction* = nullptr, OrtMemoryInfo* = nullptr);
  ~MockedOrtAllocator();

  void* Alloc(size_t size);
  void Free(void* p);
  const OrtMemoryInfo* Info() const;
  void* Reserve(size_t size);
  size_t NumAllocations() const;
  size_t NumReserveAllocations() const;

  void LeakCheck();

 private:
  MockedOrtAllocator(const MockedOrtAllocator&) = delete;
  MockedOrtAllocator& operator=(const MockedOrtAllocator&) = delete;

  std::atomic<size_t> memory_inuse{0};
  std::atomic<size_t> num_allocations{0};
  std::atomic<size_t> num_reserve_allocations{0};
  OrtMemoryInfo* memory_info;
  GPUReconstruction* rec;
};

MockedOrtAllocator::MockedOrtAllocator(GPUReconstruction* r, OrtMemoryInfo* info)
{
  OrtAllocator::version = ORT_API_VERSION;
  OrtAllocator::Alloc = [](OrtAllocator* this_, size_t size) { return static_cast<MockedOrtAllocator*>(this_)->Alloc(size); };
  OrtAllocator::Free = [](OrtAllocator* this_, void* p) { static_cast<MockedOrtAllocator*>(this_)->Free(p); };
  OrtAllocator::Info = [](const OrtAllocator* this_) { return static_cast<const MockedOrtAllocator*>(this_)->Info(); };
  OrtAllocator::Reserve = [](OrtAllocator* this_, size_t size) { return static_cast<MockedOrtAllocator*>(this_)->Reserve(size); };
  rec = r;
  memory_info = info;
}

MockedOrtAllocator::~MockedOrtAllocator()
{
  // Ort::GetApi().ReleaseMemoryInfo(memory_info);
}

void* MockedOrtAllocator::Alloc(size_t size)
{
  // LOG(info) << "(ORT) Allocating volatile memory of size " << size << " bytes";
  return rec->AllocateVolatileDeviceMemory(size);
}

void* MockedOrtAllocator::Reserve(size_t size)
{
  // LOG(info) << "(ORT) Reserving volatile memory of size " << size << " bytes";
  return rec->AllocateVolatileDeviceMemory(size);
}

void MockedOrtAllocator::Free(void* p)
{
  // LOG(info) << "(ORT) Freeing volatile memory " << p;
  rec->ReturnVolatileDeviceMemory();
}

const OrtMemoryInfo* MockedOrtAllocator::Info() const
{
  return memory_info;
}

size_t MockedOrtAllocator::NumAllocations() const
{
  return num_allocations.load();
}

size_t MockedOrtAllocator::NumReserveAllocations() const
{
  return num_reserve_allocations.load();
}

void MockedOrtAllocator::LeakCheck()
{
  if (memory_inuse.load())
    LOG(warning) << "memory leak!!!";
}

void GPUTPCNNClusterizerHost::volatileOrtAllocator(Ort::Env* env, Ort::MemoryInfo* memInfo, GPUReconstruction* rec, bool recreate)
{
  mockedAlloc = std::make_shared<MockedOrtAllocator>(rec, (OrtMemoryInfo*)(*memInfo));
  if (recreate) {
    Ort::ThrowOnError(Ort::GetApi().UnregisterAllocator((OrtEnv*)(*env), (OrtMemoryInfo*)(*memInfo)));
  }
  Ort::ThrowOnError(Ort::GetApi().RegisterAllocator((OrtEnv*)(*env), mockedAlloc.get()));
  memInfo = (Ort::MemoryInfo*)mockedAlloc->Info();
}

const OrtMemoryInfo* GPUTPCNNClusterizerHost::getMockedMemoryInfo()
{
  return mockedAlloc->Info();
}

MockedOrtAllocator* GPUTPCNNClusterizerHost::getMockedAllocator()
{
  return mockedAlloc.get();
}
