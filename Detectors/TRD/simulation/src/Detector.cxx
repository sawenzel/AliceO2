// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TClonesArray.h"
#include "FairVolume.h"
#include "FairRootManager.h"
#include "TRDSimulation/Detector.h"
#include "TRDBase/TRDGeometry.h"
#include "TRDBase/TRDCommonParam.h"
#include <vector>

using namespace o2::trd;

Detector::Detector(const char* Name, Bool_t Active):
o2::Base::Detector(Name, Active),
  mHitCollection(new TClonesArray("o2::BasicXYZEHit<float>"))
{
}

void Detector::Initialize(){
}

Bool_t Detector::ProcessHits(FairVolume* v) {
  return true;
}

void Detector::Register(){
  FairRootManager::Instance()->Register("TRDHit", "TRD", mHitCollection, kTRUE);
}

TClonesArray* Detector::GetCollection(Int_t iColl) const {
  if(iColl == 0) return mHitCollection;
  return nullptr;
}

void Detector::Reset() {
}

//_____________________________________________________________________________
void Detector::createMaterials()
{
  //
  // Create the materials for the TRD
  //

  Int_t   isxfld = 2;//((AliMagF *) TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = 10.;//((AliMagF *) TGeoGlobalMagField::Instance()->GetField())->Max();
  
  //////////////////////////////////////////////////////////////////////////
  //     Define Materials 
  //////////////////////////////////////////////////////////////////////////

  // Aluminum
  Material( 1,"Al",  26.98, 13.0, 2.7,    8.9,  37.2);
  // Copper
  Material( 2,"Cu",  63.54, 29.0, 8.96,   1.43, 14.8);
  // Carbon
  Material( 3,"C" ,  12.01,  6.0, 2.265, 18.8,  74.4);
  // Carbon for fiber mats
  Material( 4,"C2",  12.01,  6.0, 1.75,  18.8,  74.4);
  // Zinc
  Material( 5,"Sn", 118.71, 50.0, 7.31,   1.21, 14.8);
  // Silicon
  Material( 6,"Si",  28.09, 14.0, 2.33,   9.36, 37.2);
  // Iron
  Material( 7,"Fe",  55.85, 26.0, 7.87,   1.76, 14.8);

  // Air  
  Float_t aAir[4]   = { 12.011   , 14.0     , 15.9994  , 36.0      };
  Float_t zAir[4]   = {  6.0     ,  7.0     ,  8.0     , 18.0      };
  Float_t wAir[4]   = {  0.000124,  0.755267,  0.231781,  0.012827 };
  Float_t dAir      = 1.20479e-03;
  Mixture(51,"Air",          aAir,   zAir,   dAir,    4, wAir  );
  // Polyethilene (CH2) 
  Float_t ape[2]    = { 12.011 ,  1.0079 };
  Float_t zpe[2]    = {  6.0   ,  1.0    };
  Float_t wpe[2]    = {  1.0   ,  2.0    };
  Float_t dpe       = 0.95;
  Mixture(52,"Polyethilene", ape,    zpe,    dpe,    -2, wpe   );
  // Gas mixtures
  // Xe/CO2-gas-mixture (85% / 15%) 
  Float_t aXeCO2[3] = { 131.29   ,  12.0107 ,  15.9994  };
  Float_t zXeCO2[3] = {  54.0    ,   6.0    ,   8.0     };
  Float_t wXeCO2[3] = {   8.5    ,   1.5    ,   3.0     }; 
  Float_t fxc       = 0.85;
  Float_t dxe       = 0.00549; // at 20C
  Float_t dco       = 0.00186; // at 20C
  Float_t dgmXe     = fxc * dxe + (1.0 - fxc) * dco;
  // Ar/CO2-gas-mixture
  Float_t aArCO2[3] = {  39.948  ,  12.0107 ,  15.9994  };
  Float_t zArCO2[3] = {  18.0    ,   6.0    ,   8.0     };
 Float_t wArCO2[3] = {   8.2    ,   1.8    ,   3.6     }; 
  Float_t fac       = 0.82;
  Float_t dar       = 0.00166; // at 20C
  Float_t dgmAr     = fac * dar + (1.0 - fac) * dco;
  if      (TRDCommonParam::Instance()->IsXenon()) {
    Mixture(53,"XeCO2",        aXeCO2, zXeCO2, dgmXe,  -3, wXeCO2);
  }
  else if (TRDCommonParam::Instance()->IsArgon()) {
    //AliInfo("Gas mixture: Ar C02 (80/20)");
    LOG(INFO) << "Gas mixture: Ar C02 (80/20)\n";
    Mixture(53,"ArCO2",        aArCO2, zArCO2, dgmAr,  -3, wArCO2);
  }
  else {
    //AliFatal("Wrong gas mixture");
    exit(1);
  }
  // G10
  Float_t aG10[4]   = {  1.0079, 12.011 , 15.9994, 28.086  };
  Float_t zG10[4]   = {  1.0   ,  6.0   ,  8.0   , 14.0    };
  Float_t wG10[4]   = {  0.023 ,  0.194 ,  0.443 ,  0.340  };
  Float_t dG10      = 2.0;
  Mixture(54,"G10",          aG10,  zG10,  dG10,   4,wG10  );
  // Water
  Float_t awa[2]    = {  1.0079, 15.9994 };
  Float_t zwa[2]    = {  1.0   ,  8.0    };
  Float_t wwa[2]    = {  2.0   ,  1.0    };
  Float_t dwa       = 1.0;
  Mixture(55,"Water",        awa,   zwa,   dwa,   -2,wwa   );
  // Rohacell (C5H8O2), X0 = 535.005cm
  Float_t arh[3]    = { 12.011 ,  1.0079, 15.9994 };
  Float_t zrh[3]    = {  6.0   ,  1.0   ,  8.0    };
  Float_t wrh[3]    = {  5.0   ,  8.0   ,  2.0    };
  Float_t drh       = 0.075;   
  Mixture(56,"Rohacell",     arh,   zrh,   drh,   -3,wrh   );
  // Epoxy (C18H19O3)
  Float_t aEpoxy[3] = { 15.9994,  1.0079, 12.011  }; 
  Float_t zEpoxy[3] = {  8.0   ,  1.0   ,  6.0    }; 
  Float_t wEpoxy[3] = {  3.0   , 19.0   , 18.0    }; 
  Float_t dEpoxy    = 1.8 ; 
  Mixture(57,"Epoxy",        aEpoxy,zEpoxy,dEpoxy,-3,wEpoxy);
  // Araldite, low density epoxy (C18H19O3)
  Float_t aAral[3]  = { 15.9994,  1.0079, 12.011  }; 
  Float_t zAral[3]  = {  8.0   ,  1.0   ,  6.0    }; 
  Float_t wAral[3]  = {  3.0   , 19.0   , 18.0    }; 
  Float_t dAral     = 1.12; // Hardener: 1.15, epoxy: 1.1, mixture: 1/2
  Mixture(58,"Araldite",     aAral, zAral, dAral, -3,wAral );
  // Mylar
  Float_t aMy[3]    = { 12.011 ,   1.0  , 15.9994 };
  Float_t zMy[3]    = {  6.0   ,   1.0  ,  8.0    };
  Float_t wMy[3]    = {  5.0   ,   4.0  ,  2.0    };
  Float_t dMy       = 1.39;
  Mixture(59,"Mylar",        aMy,   zMy,   dMy,   -3,wMy   );
  // Polypropylene (C3H6) for radiator fibers
  Float_t app[2]    = { 12.011 ,  1.0079 };
  Float_t zpp[2]    = {  6.0   ,  1.0    };
  Float_t wpp[2]    = {  3.0   ,  6.0    };
  Float_t dpp       = 0.068;
  Mixture(60,"Polypropylene",app,   zpp,   dpp,   -2,wpp   );
  // Aramide for honeycomb
  Float_t aAra[4]   = {  1.0079, 12.011 , 15.9994, 14.0067 };
  Float_t zAra[4]   = {  1.0   ,  6.0   ,  8.0   ,  7.0    };
  Float_t wAra[4]   = {  3.0   ,  1.0   ,  1.0   ,  1.0    };
  Float_t dAra      = 0.032;
  Mixture(61,"Aramide",      aAra,  zAra,  dAra,  -4,wAra  );
  // GFK for Wacosit (Epoxy + Si)
  Float_t aGFK[4]   = {  1.0079, 12.011 , 15.9994, 28.086  };
  Float_t zGFK[4]   = {  1.0   ,  6.0   ,  8.0   , 14.0    };
  Float_t wGFK[4]   = {  0.0445,  0.5031,  0.1118,  0.340  };
  Float_t dGFK      = 2.0;
  Mixture(62,"GFK",          aGFK,  zGFK,  dGFK,   4,wGFK  );

  //////////////////////////////////////////////////////////////////////////
  //     Tracking Media Parameters 
  //////////////////////////////////////////////////////////////////////////

  // General tracking parameter
  Float_t tmaxfd    = -10.0;
  Float_t stemax    = -1.0e10;
  Float_t deemax    = -0.1;
  Float_t epsil     =  1.0e-4;
  Float_t stmin     = -0.001;

  // Al Frame 
  Medium( 1,"Al Frame"   , 1,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Air 
  Medium( 2,"Air"        ,51,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Wires
  Medium( 3,"Wires"      , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // All other ROB materials (caps, etc.)
  Medium( 4,"ROB Other"  , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu pads 
  Medium( 5,"Padplane"   , 2,1,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Fee + cables 
  Medium( 6,"Readout"    , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // C frame (Wacosit) 
  Medium( 7,"Wacosit"    ,62,0,isxfld,sxmgmx
             ,tmaxfd,stemax,deemax,epsil,stmin);
  // INOX of cooling bus bars
  Medium( 8,"Cooling bus", 7,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Gas-mixture (Xe/CO2) 
  Medium( 9,"Gas-mix"    ,53,1,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Honeycomb
  Medium(10,"Honeycomb"  ,61,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Araldite glue
  Medium(11,"Glue"       ,58,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // G10-plates
  Medium(13,"G10-plates" ,54,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cooling water
  Medium(14,"Water"      ,55,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Rohacell for the radiator
  Medium(15,"Rohacell"   ,56,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Al layer in MCMs
  Medium(16,"MCM-Al"     , 1,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Sn layer in MCMs
  Medium(17,"MCM-Sn"     , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu layer in MCMs
  Medium(18,"MCM-Cu"     , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // G10 layer in MCMs
  Medium(19,"MCM-G10"    ,54,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Si in readout chips
  Medium(20,"Chip-Si"    , 6,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Epoxy in readout chips
  Medium(21,"Chip-Ep"    ,57,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // PE in connectors
  Medium(22,"Conn-PE"    ,52,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu in connectors
  Medium(23,"Chip-Cu"    , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Al of cooling pipes
  Medium(24,"Cooling"    , 1,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu in services
  Medium(25,"Serv-Cu"    , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Carbon fiber mat
  Medium(26,"Carbon"     , 4,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Mylar foil
  Medium(27,"Mylar"      ,59,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Polypropylene fibers
  Medium(28,"Fiber"      ,60,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);

  // Save the density values for the TRD absorbtion
  Float_t dmy  = 1.39;
  fFoilDensity = dmy;
  if      (TRDCommonParam::Instance()->IsXenon()) {
    fGasDensity       = dgmXe;
    fGasNobleFraction = fxc;
  }
  else if (TRDCommonParam::Instance()->IsArgon()) {
    fGasDensity       = dgmAr;
    fGasNobleFraction = fac;
  }
}

// setting up the geometry
void Detector::ConstructGeometry() {
  createMaterials();

  std::vector<int> medmapping;
  // now query the medium mapping and fill a vector to be passed along
  // to the geometry creation
  getMediumIDMappingAsVector(medmapping);

  TRDGeometry geom;
  geom.CreateGeometry(medmapping);
}


ClassImp(Detector);
