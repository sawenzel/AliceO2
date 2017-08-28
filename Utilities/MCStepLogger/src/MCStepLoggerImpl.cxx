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
#include <TGeoManager.h>
#include <TClonesArray.h>
#include <sstream>
#include <StepInfo.h>
#include <TFile.h>
#include <TBranch.h>
#include <FairModule.h>
#include <TGeoVolume.h>

#include <dlfcn.h>
#include <iostream>
#include <map>
#include <set>
#include <cstdlib>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include <execinfo.h> // for backtrace
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
// This function produces a stack backtrace with demangled function & method names.
std::string Backtrace(int skip = 1)
{
    void *callstack[128];
    const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
    char buf[1024];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    std::ostringstream trace_buf;
    for (int i = skip; i < nFrames; i++) {
        printf("%s\n", symbols[i]);

        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = NULL;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                     info.dli_sname == 0 ? symbols[i] : info.dli_sname,
                     (char *)callstack[i] - (char *)info.dli_saddr);
            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
        }
        trace_buf << buf;
    }
    free(symbols);
    if (nFrames == nMaxFrames)
        trace_buf << "[truncated]\n";
    return trace_buf.str();
}

namespace o2
{

const char* getLogFileName() 
{
  if(const char* f = std::getenv("MCSTEPLOG_OUTFILE")) 
  {
    return f;
  }
  else {
    return "MCStepLoggerOutput.root";
  }
}

template <typename T>
void flushToTTree(const char* branchname, T* address) {
  TFile *f = new TFile(getLogFileName(), "UPDATE");
  const char* treename = "StepLoggerTree";
  auto tree = (TTree*)f->Get(treename);
  if (!tree) {
    // create tree
    tree = new TTree(treename, "Tree container information from MC step logger");
  }
  auto branch = tree->GetBranch(branchname);
  if (!branch) {
   branch = tree->Branch(branchname, &address);
  }
  branch->SetAddress(&address);
  branch->Fill();
  tree->SetEntries(branch->GetEntries());
  f->Write();
  f->Close();
  delete f;
}


// a class collecting field access per volume
class FieldLogger
{
  int counter = 0;
  std::map<int, int> volumetosteps;
  std::map<int, std::string> idtovolname;
  bool mTTreeIO = false;
  std::vector<MagCallInfo> callcontainer;
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
    if (mTTreeIO){
      callcontainer.emplace_back(mc, x[0], x[1], x[2] ,b[0], b[1], b[2]);
      return;
    }
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
    if (mTTreeIO) {
      callcontainer.clear();
    }
  }

  void flush()
  {
    if (mTTreeIO) {
     flushToTTree("Calls", &callcontainer);
    }
    std::cerr << "[FIELDLOGGER]: did " << counter << " steps \n";
    // summarize steps per volume
    for (auto& p : volumetosteps) {
      std::cerr << "[FIELDLOGGER]: VolName " << idtovolname[p.first] << " COUNT " << p.second;
      std::cerr << "\n";
    }
    std::cerr << "[FIELDLOGGER]: ----- END OF EVENT ------\n";
    clear();
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
    if (mTTreeIO){
      container.clear();
    }
    StepInfo::resetCounter();
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
     std::cerr << "[STEPLOGGER]: ----- END OF EVENT ------\n";
    }
    else {
      flushToTTree("Steps", &container);
    }
    clear();
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

extern "C" void dispatchConstructGeom(FairModule* module, char const* libname, char const* origFunctionName) {
  auto list = (TObjArray*)gGeoManager->GetListOfVolumes();
  auto sizebefore = list->GetEntries();
    
  typedef void (FairModule::*MethodType)();
  dispatchOriginalKernel<FairModule, MethodType>(module, libname, origFunctionName);
  // record changes

  auto sizeafter = list->GetEntries();
  std::cerr << module->GetName() << " added " << sizeafter-sizebefore << " volumes \n";
}

extern "C" void dispatchAddVolume(TGeoManager* geom, char const* libname, char const* origFunctionName, TGeoVolume *v) {
  typedef int (TGeoManager::*MethodType)(TGeoVolume*);
  Backtrace();
  dispatchOriginalKernel<TGeoManager, MethodType>(geom, libname, origFunctionName, v);
}

extern "C" void dispatchVolume(TGeoVolume* v, char const* libname, char const* origFunctionName,
			       const char *name, const TGeoShape* s, const TGeoMedium* m)
{
  std::cerr << "Adding volume " << name << "\n";  
  typedef void (TGeoVolume::*MethodType)(const char *, const TGeoShape*, const TGeoMedium*);
  Backtrace();
  dispatchOriginalKernel<TGeoVolume, MethodType>(v, libname, origFunctionName, name, s, m);
}

extern "C" void dispatchBuilder(TGeoManager* m, char const* libname, char const* origFunctionName,
				const char *name, const char *shape, int nmed, double *p, int npar)
{
  std::cerr << "Adding volume " << name << "\n";  
  typedef void (TGeoManager::*MethodType)(const char *, const char *, int, double *, int);
  Backtrace();
  dispatchOriginalKernel<TGeoManager, MethodType>(m, libname, origFunctionName, name, shape, nmed, p, npar);
}

extern "C" void dispatchBuilderF(TGeoManager* m, char const* libname, char const* origFunctionName,
				const char *name, const char *shape, int nmed, float *p, int npar)
{
  std::cerr << "Adding volume COOL " << name << "\n";  
  typedef void (TGeoManager::*MethodType)(const char *, const char *, int, float *, int);
  Backtrace();
  dispatchOriginalKernel<TGeoManager, MethodType>(m, libname, origFunctionName, name, shape, nmed, p, npar);
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
  std::cerr << "--- START FLUSHING ----\n";
  o2::logger.flush();
  o2::fieldlogger.flush();
  std::cerr << "--- END FLUSHING ----\n";
}
