// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//****************************************************************************
//* This file is free software: you can redistribute it and/or modify        *
//* it under the terms of the GNU General Public License as published by     *
//* the Free Software Foundation, either version 3 of the License, or        *
//* (at your option) any later version.                                      *
//*                                                                          *
//* Primary Authors: Sandro Wenzel <sandro.wenzel@cern.ch>                   *
//*                                                                          *
//* The authors make no claims about the suitability of this software for    *
//* any purpose. It is provided "as is" without express or implied warranty. *
//****************************************************************************

//  @file   MCStepLoggerImpl.cxx
//  @author Sandro Wenzel
//  @since  2017-06-29
//  @brief  A logging service for MCSteps (hooking into Stepping of TVirtualMCApplication's)

#include <TTree.h>
#include <TVirtualMC.h>
#include <TVirtualMCApplication.h>
#include <TVirtualMagField.h>
#include <TClonesArray.h>
#include <sstream>
#include <StepInfo.h>
#include <TFile.h>
#include <TBranch.h>

#include <dlfcn.h>
#include <iostream>
#include <map>
#include <set>
#include <cstdlib>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

namespace o2
{
namespace {
  struct LogOutFile {
    TFile *outfile;
    TTree *outtree;
  
    static LogOutFile & Instance() {
      static LogOutFile instance;
      return instance;
    }

    private:
      LogOutFile() {
        outfile=new TFile("MCStepLogResult.root", "RECREATE");
        outtree=new TTree("StepLogTree", "Tree produced from MC Step Logger");
      }
  };
}


// a class collecting field access per volume
class FieldLogger
{
  int counter = 0;
  std::map<int, int> volumetosteps;
  std::map<int, std::string> idtovolname;

  bool mTTreeIO = false;

  std::vector<MagCallInfo> callcontainer;
  TBranch* outbranch = nullptr;
 public:

  FieldLogger() {
    // check if streaming or interactive
    // configuration done via env variable
    if (std::getenv("MCSTEPLOG_TTREE")) {
      mTTreeIO = true;
    }
  }
  
  void addStep(TVirtualMC* mc, const double *x, const double *b)
  {
    callcontainer.emplace_back(mc, x[0], x[1], x[2] ,b[0], b[1], b[2]);
    return;
    
    std::cerr << "field called at " << x[0] << " " << x[1] << " " << x[2] << " with " << b[0] << " " << b[1] << " " << b[2] << "\n";
    counter++;
    int copyNo;
    auto id = mc->CurrentVolID(copyNo);
    if (volumetosteps.find(id) == volumetosteps.end()) {
      volumetosteps.insert(std::pair<int, int>(id, 0));
    } else {
      volumetosteps[id]++;
    }
    if (idtovolname.find(id) == idtovolname.end()) {
      idtovolname.insert(std::pair<int, std::string>(id, std::string(mc->CurrentVolName())));
    }
  }

  void clear()
  {
    counter = 0;
    volumetosteps.clear();
    idtovolname.clear();
    callcontainer.clear();
  }

  void flush()
  {
    if (mTTreeIO) {
      auto logfile = LogOutFile::Instance();  
      if (logfile.outtree->GetBranch("Calls") == nullptr)
      {
        LogOutFile::Instance().outtree->Branch("Calls", &callcontainer);    
      }
    }
    // std::cerr << "[FIELDLOGGER]: did " << counter << " steps \n";
    // summarize steps per volume
    //for (auto& p : volumetosteps) {
    //  std::cerr << "[FIELDLOGGER]: VolName " << idtovolname[p.first] << " COUNT " << p.second;
    // std::cerr << "\n";
    //}
    //clear();
    //std::cerr << "[FIELDLOGGER]: ----- END OF EVENT ------\n";
  }
};

class StepLogger
{
  int stepcounter = 0;

  std::set<int> trackset;
  std::set<int> pdgset;
  std::map<int, int> volumetosteps;
  std::map<int, std::string> idtovolname;
  std::map<int, int> volumetoNSecondaries;            // number of secondaries created in this volume
  std::map<std::pair<int, int>, int> volumetoProcess; // mapping of volumeid x processID to secondaries produced

  std::vector<StepInfo> container;
  bool mTTreeIO = false;

 public:
  StepLogger() {
    // check if streaming or interactive
    // configuration done via env variable
    if (std::getenv("MCSTEPLOG_TTREE")) {
      mTTreeIO = true;
    }
  }

  void addStep(TVirtualMC* mc)
  {
    if (mTTreeIO) {
      container.emplace_back(mc);
    }
    else {
    assert(mc);
    stepcounter++;

    auto stack = mc->GetStack();
    assert(stack);
    trackset.insert(stack->GetCurrentTrackNumber());
    pdgset.insert(mc->TrackPid());
    int copyNo;
    auto id = mc->CurrentVolID(copyNo);

TArrayI procs;
    mc->StepProcesses(procs);
    for (int i=0; i<procs.GetSize(); ++i){
//      std::cerr << TMCProcessName[procs.At(i)] << "\t";
    }
    
    if (volumetosteps.find(id) == volumetosteps.end()) {
      volumetosteps.insert(std::pair<int, int>(id, 0));
    } else {
      volumetosteps[id]++;
    }
    if (idtovolname.find(id) == idtovolname.end()) {
      idtovolname.insert(std::pair<int, std::string>(id, std::string(mc->CurrentVolName())));
    }

    // for the secondaries
    if (volumetoNSecondaries.find(id) == volumetoNSecondaries.end()) {
      volumetoNSecondaries.insert(std::pair<int, int>(id, mc->NSecondaries()));
    } else {
      volumetoNSecondaries[id] += mc->NSecondaries();
    }

    // for the processes
    for (int i = 0; i < mc->NSecondaries(); ++i) {
      auto process = mc->ProdProcess(i);
      auto p = std::pair<int, int>(id, process);
      if (volumetoProcess.find(p) == volumetoProcess.end()) {
        volumetoProcess.insert(std::pair<std::pair<int, int>, int>(p, 1));
      } else {
        volumetoProcess[p]++;
      }
    }
    }
  }

  void clear()
  {
    stepcounter = 0;
    trackset.clear();
    pdgset.clear();
    volumetosteps.clear();
    idtovolname.clear();
    volumetoNSecondaries.clear();
    volumetoProcess.clear();
    container.clear();
  }

  // prints list of processes for volumeID
  void printProcesses(int volid)
  {
    for (auto& p : volumetoProcess) {
      if (p.first.first == volid) {
        std::cerr << "P[" << TMCProcessName[p.first.second] << "]:" << p.second << "\t";
      }
    }
  }

  void flush()
  {
    if (!mTTreeIO) {
     std::cerr << "[STEPLOGGER]: did " << stepcounter << " steps \n";
     std::cerr << "[STEPLOGGER]: transported " << trackset.size() << " different tracks \n";
     std::cerr << "[STEPLOGGER]: transported " << pdgset.size() << " different types \n";
     // summarize steps per volume
     for (auto& p : volumetosteps) {
       std::cerr << "[STEPLOGGER]: VolName " << idtovolname[p.first] << " COUNT " << p.second << " SECONDARIES "
                << volumetoNSecondaries[p.first] << " ";
       // loop over processes
       printProcesses(p.first);
       std::cerr << "\n";
     }
     clear();
     std::cerr << "[STEPLOGGER]: ----- END OF EVENT ------\n";
    }
    else {
      auto logfile = LogOutFile::Instance();  
      if (logfile.outtree->GetBranch("Steps") == nullptr)
      {
         LogOutFile::Instance().outtree->Branch("Steps", &container);    
      }
    }
  }
};

StepLogger logger;
FieldLogger fieldlogger;

} // end namespace

// a helper template kernel describing generically the redispatching prodecure
template <typename Object /* the original object type */, typename MethodType /* member function type */,
          typename... Args /* original arguments to function */>
void dispatchOriginalKernel(Object* obj, char const* libname, char const* origFunctionName, Args... args)
{
  // Object, MethodType, and Args are of course related so we could do some static_assert checks or automatic deduction

  // static map to avoid having to lookup the right symbols in the shared lib at each call
  // (We could do this outside of course)
  static std::map<const char*, MethodType> functionNameToSymbolMap;
  MethodType origMethod = nullptr;

  auto iter = functionNameToSymbolMap.find(origFunctionName);
  if (iter == functionNameToSymbolMap.end()) {
    auto libHandle = dlopen(libname, RTLD_NOW);
    // try to make the library loading a bit more portable:
    if (!libHandle) {
      // try appending *.so
      std::stringstream stream;
      stream << libname << ".so";
      libHandle = dlopen(stream.str().c_str(), RTLD_NOW);
    }
    if (!libHandle) {
      // try appending *.dylib
      std::stringstream stream;
      stream << libname << ".dylib";
      libHandle = dlopen(stream.str().c_str(), RTLD_NOW);
    }
    assert(libHandle);
    void* symbolAddress = dlsym(libHandle, origFunctionName);
    assert(symbolAddress);
// Purposely ignore compiler warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsizeof-pointer-memaccess"
    // hack since C++ does not allow casting to C++ member function pointers
    // thanks to gist.github.com/mooware/1174572
    memcpy(&origMethod, &symbolAddress, sizeof(&symbolAddress));
#pragma GCC diagnostic pop
    functionNameToSymbolMap[origFunctionName] = origMethod;
  } else {
    origMethod = iter->second;
  }
  // the final C++ member function call redispatch
  (obj->*origMethod)(args...);
}

// a generic function that can dispatch to the original method of a TVirtualMCApplication
extern "C" void dispatchOriginal(TVirtualMCApplication* app, char const* libname, char const* origFunctionName)
{
  typedef void (TVirtualMCApplication::*StepMethodType)();
  dispatchOriginalKernel<TVirtualMCApplication, StepMethodType>(app, libname, origFunctionName);
}

// a generic function that can dispatch to the original method of a TVirtualMagField
extern "C" void dispatchOriginalField(TVirtualMagField* field, char const* libname, char const* origFunctionName,
                                      const double x[3], double* B)
{
  typedef void (TVirtualMagField::*MethodType)(const double[3], double*);
  dispatchOriginalKernel<TVirtualMagField, MethodType>(field, libname, origFunctionName, x, B);
}

extern "C" void performLogging(TVirtualMCApplication* app)
{
  static TVirtualMC* mc = TVirtualMC::GetMC();
  o2::logger.addStep(mc);
}

extern "C" void logField(double const *p, double const *b)
{
  static TVirtualMC* mc = TVirtualMC::GetMC();
  o2::fieldlogger.addStep(mc, p, b);
}

extern "C" void flushLog()
{
  o2::logger.flush();
  o2::fieldlogger.flush();
  o2::LogOutFile::Instance().outtree->Fill();
  o2::LogOutFile::Instance().outtree->Write();
}
