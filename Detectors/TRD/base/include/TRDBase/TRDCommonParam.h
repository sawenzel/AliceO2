#ifndef O2_TRDCOMMONPARAM_H
#define O2_TRDCOMMONPARAM_H

#include "TRDBase/TRDSimParam.h"

class TRDPadPlane;

namespace o2 {
namespace trd {

class TRDCommonParam
{
  public:
    enum { kNlayer  = 6
         , kNstack  = 5
         , kNsector = 18
         , kNdet    = 540 };

    enum { kXenon =   0
	 , kArgon =   1   };
    
    TRDCommonParam(const TRDCommonParam &p);   
    TRDCommonParam &operator=(const TRDCommonParam &p); 
    ~TRDCommonParam();

    static TRDCommonParam *Instance();
    static  void    Terminate();
    
    void            SetExB(Int_t exbOn = 1)                        { fExBOn             = exbOn;    }
    void            SetSamplingFrequency(Float_t freq)             { fSamplingFrequency = freq;     }
    void            SetXenon()                                     { fGasMixture        = kXenon; 
                                                                     TRDSimParam::Instance()->ReInit(); }
    void            SetArgon()                                     { fGasMixture        = kArgon; 
                                                                     TRDSimParam::Instance()->ReInit(); }

    Bool_t          ExBOn() const                                  { return fExBOn;                 }
    Bool_t          IsXenon() const                                { return (fGasMixture == kXenon) 
                                                                     ? kTRUE : kFALSE;              }
    Bool_t          IsArgon() const                                { return (fGasMixture == kArgon) 
                                                                     ? kTRUE : kFALSE;              }

    Int_t           GetGasMixture() const                          { return fGasMixture;            }
    Float_t         GetSamplingFrequency() const                   { return fSamplingFrequency;     }

    Float_t         GetOmegaTau(Float_t vdrift);
    Bool_t          GetDiffCoeff(Float_t &dl, Float_t &dt, Float_t vdrift);

    Double_t        TimeStruct(Float_t vdrift, Double_t xd, Double_t z);

  protected:

    void            SampleTimeStruct(Float_t vdrift);

    static TRDCommonParam *fgInstance;          //  Instance of this class (singleton implementation)
    static Bool_t             fgTerminated;        //  Defines if this class has already been terminated    

    Int_t                     fExBOn;              //  Switch for the ExB effects

    Float_t                   fDiffusionT;         //  Transverse drift coefficient
    Float_t                   fDiffusionL;         //  Longitudinal drift coefficient
    Float_t                   fDiffLastVdrift;     //  The structures are valid for fLastVdrift (caching)

    Float_t                  *fTimeStruct1;        //! Time Structure of Drift Cells
    Float_t                  *fTimeStruct2;        //! Time Structure of Drift Cells
    Float_t                   fVDlo;               //  Lower drift velocity, for interpolation
    Float_t                   fVDhi;               //  Higher drift velocity, for interpolation
    Float_t                   fTimeLastVdrift;     //  The structures are valid for fLastVdrift (caching)

    Float_t                   fSamplingFrequency;  //  Sampling Frequency in MHz

    Int_t                     fGasMixture;         //  Gas mixture: 0-Xe/C02 1-Ar/CO2. 
  
  private:

    // This is a singleton, constructor is private!  
    TRDCommonParam();
  
    ClassDef(TRDCommonParam,1)                  // The constant parameters common to simulation and reconstruction

};
}
}
#endif
