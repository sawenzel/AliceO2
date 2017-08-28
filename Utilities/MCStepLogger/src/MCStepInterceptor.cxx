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

//  @file   MCStepInterceptor.cxx
//  @author Sandro Wenzel
//  @since  2017-06-29
//  @brief  A LD_PRELOAD logger hooking into Stepping of TVirtualMCApplication's

class TVirtualMCApplication;
class TVirtualMagField;
#include <iostream>

// (re)declare symbols to be able to hook into them
#define DECLARE_INTERCEPT_SYMBOLS(APP) \
  class APP                            \
  {                                    \
   public:                             \
    void Stepping();                   \
    void FinishEvent();                \
  };

DECLARE_INTERCEPT_SYMBOLS(FairMCApplication)
DECLARE_INTERCEPT_SYMBOLS(AliMC)

// same for field
#define DECLARE_INTERCEPT_FIELD_SYMBOLS(FIELD)                 \
  class FIELD                                                  \
  {                                                            \
   public:                                                     \
    void Field(const double* point, double* bField);           \
  };

namespace o2
{
namespace field
{
DECLARE_INTERCEPT_FIELD_SYMBOLS(MagneticField);
}
}
DECLARE_INTERCEPT_FIELD_SYMBOLS(AliMagF);

extern "C" void performLogging(TVirtualMCApplication*);
extern "C" void logField(const double*, const double*);
extern "C" void dispatchOriginal(TVirtualMCApplication*, char const* libname, char const*);
extern "C" void dispatchOriginalField(TVirtualMagField*, char const* libname, char const*, const double x[3],
                                      double* B);
extern "C" void flushLog();


#define INTERCEPT_STEPPING(APP, LIB, SYMBOL)                       \
  void APP::Stepping()                                             \
  {                                                                \
    auto baseptr = reinterpret_cast<TVirtualMCApplication*>(this); \
    performLogging(baseptr);                                       \
    dispatchOriginal(baseptr, LIB, SYMBOL);                        \
  }

#define INTERCEPT_FINISHEVENT(APP, LIB, SYMBOL)                    \
  void APP::FinishEvent()                                          \
  {                                                                \
    auto baseptr = reinterpret_cast<TVirtualMCApplication*>(this); \
    flushLog();                                                    \
    dispatchOriginal(baseptr, LIB, SYMBOL);                        \
  }

// the runtime will now dispatch to these functions due to LD_PRELOAD
INTERCEPT_STEPPING(FairMCApplication, "libBase", "_ZN17FairMCApplication8SteppingEv")
INTERCEPT_STEPPING(AliMC, "libSTEER", "_ZN5AliMC8SteppingEv")

INTERCEPT_FINISHEVENT(FairMCApplication, "libBase", "_ZN17FairMCApplication11FinishEventEv")
INTERCEPT_FINISHEVENT(AliMC, "libSTEER", "_ZN5AliMC11FinishEventEv")

#define INTERCEPT_FIELD(FIELD, LIB, SYMBOL)                        \
  void FIELD::Field(const double* point, double* bField)           \
  {                                                                \
    auto baseptr = reinterpret_cast<TVirtualMagField*>(this);      \
    dispatchOriginalField(baseptr, LIB, SYMBOL, point, bField);    \
    logField(point, bField);           				   \
  }

namespace o2
{
namespace field
{
INTERCEPT_FIELD(MagneticField, "libField", "_ZN2o25field13MagneticField5FieldEPKdPd");
}
}
// for AliRoot
INTERCEPT_FIELD(AliMagF, "libSTEERBase", "_ZN7AliMagF5FieldEPKdPd")



// intercept geometry construction since
// we can then create a mapping of volumename -> module useful for histogramming
//class FairModule {
//  void ConstructGeometry();
//};
//extern "C" void dispatchConstructGeom(FairModule*, char const* libname, char const*);
//void FairModule::ConstructGeometry(){
//  auto baseptr = reinterpret_cast<FairModule*>(this);
//  dispatchConstructGeom(baseptr, "libBase", "_ZN10FairModule17ConstructGeometryEv");
//};


// intercept geometry construction since
// we can then create a mapping of volumename -> module useful for histogramming
class TGeoShape;
class TGeoMedium;
class TGeoVolume {
  TGeoVolume(const char *name, const TGeoShape*, const TGeoMedium*);
};

class TGeoManager {
  int AddVolume(TGeoVolume *);
  TGeoVolume* Volume(const char *, const char *, int,
		     double *, int);
  TGeoVolume* Volume(const char *, const char *, int,
		     float *, int);
};

//Int_t TGeoManager::AddVolume(TGeoVolume *volume)
extern "C" void dispatchAddVolume(TGeoManager*, char const* libname, char const*, TGeoVolume *);

int TGeoManager::AddVolume(TGeoVolume *v){
  std::cerr << "HUHU\n";
  //auto baseptr = reinterpret_cast<TGeoManager*>(this);
  dispatchAddVolume(this, "libGeom", "_ZN11TGeoManager9AddVolumeEP10TGeoVolume", v);
};

extern "C" void dispatchBuilder(TGeoManager*, char const* libname, char const*, char const *, const char *, int, double *, int);
TGeoVolume* TGeoManager::Volume(const char *name, const char *shape, int nmed,
				double *par, int npar){
  dispatchBuilder(this,  "libGeom", "_ZN11TGeoManager6VolumeEPKcS1_iPdi", name, shape, nmed, par, npar);
}

extern "C" void dispatchBuilderF(TGeoManager*, char const* libname, char const*, char const *, const char *, int, float *, int);
TGeoVolume* TGeoManager::Volume(const char *name, const char *shape, int nmed,
				float *par, int npar){
  dispatchBuilderF(this,  "libGeom", "_ZN11TGeoManager6VolumeEPKcS1_iPfi", name, shape, nmed, par, npar);
}



extern "C" void dispatchVolume(TGeoVolume*, char const* libname, char const*, char const *, const TGeoShape *sp, const TGeoMedium *m);
TGeoVolume::TGeoVolume(const char *name, const TGeoShape* shape, const TGeoMedium* medium) {
  dispatchVolume(this, "libGeom", "_ZN10TGeoVolumeC1EPKcPK9TGeoShapePK10TGeoMedium", name, shape, medium);
};
