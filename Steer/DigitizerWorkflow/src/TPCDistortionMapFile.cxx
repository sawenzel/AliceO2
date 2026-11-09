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

#include "TPCDistortionMapFile.h"
#include "Framework/InputRecord.h"
#include "Framework/Logger.h"
#include "Framework/ProcessingContext.h"
#include "TPCCalibration/CorrMapParam.h"
#include "TPCSimulation/Digitizer.h"
#include "TPCSpaceCharge/SpaceCharge.h"
#include <Compression.h>
#include <TFile.h>
#include <TNamed.h>
#include <TStopwatch.h>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#include <filesystem>
#include <memory>
#include <fmt/format.h>

using namespace o2::tpc;
using SC = o2::tpc::SpaceCharge<float>;

namespace
{
constexpr const char* ObjectName = "ccdb_object";
constexpr const char* KeyName = "tpc_distortion_key";

/// exclusive lock on a file, released at destruction
class FileLock
{
 public:
  explicit FileLock(const std::string& path)
  {
    mFD = ::open(path.c_str(), O_CREAT | O_RDWR, 0644);
    if (mFD < 0 || ::flock(mFD, LOCK_EX) != 0) {
      LOGP(fatal, "TPC: cannot lock {}", path);
    }
  }
  ~FileLock()
  {
    ::flock(mFD, LOCK_UN);
    ::close(mFD);
  }
  FileLock(const FileLock&) = delete;
  FileLock& operator=(const FileLock&) = delete;

 private:
  int mFD = -1;
};

std::string getETag(o2::framework::ProcessingContext& pc, const char* binding)
{
  const auto& headers = pc.inputs().get<o2::framework::CCDBMetadataExtractor>(binding);
  const auto it = headers.find("ETag");
  return it == headers.end() ? std::string{} : it->second;
}

std::unique_ptr<SC> readMap(const std::string& path, const std::string& key)
{
  if (!std::filesystem::exists(path)) {
    return nullptr;
  }
  TFile file(path.c_str(), "READ");
  std::unique_ptr<TNamed> storedKey(file.IsZombie() ? nullptr : file.Get<TNamed>(KeyName));
  if (!storedKey || key != storedKey->GetTitle()) {
    LOGP(warning, "TPC: {} was written for other distortion inputs and will be replaced", path);
    return nullptr;
  }
  return std::unique_ptr<SC>(file.Get<SC>(ObjectName));
}

void writeMap(const std::string& path, const std::string& key, const SC& sc)
{
  // written under another name and renamed, so that a crash never leaves a partial file
  const auto tmpPath = fmt::format("{}.{}.tmp", path, ::getpid());
  {
    TFile file(tmpPath.c_str(), "RECREATE", "", ROOT::RCompressionSetting::EDefaults::kUseAnalysis);
    file.WriteObjectAny(&sc, "o2::tpc::SpaceCharge<float>", ObjectName);
    TNamed(KeyName, key.c_str()).Write();
  }
  std::filesystem::rename(tmpPath, path);
}
} // namespace

void TPCDistortionMapFile::prepare(o2::framework::ProcessingContext& pc, Digitizer& digitizer, float meanLumi, float meanLumiDerivative,
                                   const std::function<void()>& compute)
{
  const auto etagDist = getETag(pc, "tpcdistortions");
  const auto etagDeriv = getETag(pc, "tpcdistortionsderiv");
  if (etagDist.empty() || etagDeriv.empty()) {
    LOGP(warning, "TPC: no ETag for the distortion maps, the scaled distortions are computed and not shared through {}", mPath);
    compute();
    return;
  }
  const auto key = fmt::format("dist {} deriv {} lumiInst {} meanLumi {} meanLumiDerivative {}", etagDist, etagDeriv,
                               CorrMapParam::Instance().lumiInst, meanLumi, meanLumiDerivative);
  if (key == mKey) {
    return; // the digitizer holds them already
  }

  FileLock lock(mPath + ".lock");
  TStopwatch timer;
  if (auto sc = readMap(mPath, key)) {
    digitizer.setScaledSCDistortions(std::move(sc));
    LOGP(info, "TPC: scaled distortions read from {} in {:.1f} s", mPath, timer.RealTime());
  } else {
    compute();
    timer.Start();
    writeMap(mPath, key, *digitizer.getSCDistortions());
    LOGP(info, "TPC: scaled distortions written to {} in {:.1f} s", mPath, timer.RealTime());
  }
  mKey = key;
}
