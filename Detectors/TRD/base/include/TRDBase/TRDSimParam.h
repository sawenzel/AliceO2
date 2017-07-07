#ifndef O2_TRDSIMPARAM_H
#define O2_TRDSIMPARAM_H

namespace o2 {
namespace trd {

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Class containing constant simulation parameters                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
class TRDSimParam
{
  
 public:
  
          enum { kNplan =   6
               , kNcham =   5
               , kNsect =  18
               , kNdet  = 540 };

  static  TRDSimParam *Instance();
  static  void     Terminate();

          void     SetGasGain(float gasgain)               { fGasGain           = gasgain;          }
          void     SetNoise(float noise)                   { fNoise             = noise;            }
          void     SetChipGain(float chipgain)             { fChipGain          = chipgain;         }
          void     SetADCoutRange(float range)             { fADCoutRange       = range;            }
          void     SetADCinRange(float range)              { fADCinRange        = range;            }
          void     SetADCbaseline(int basel)               { fADCbaseline       = basel;            }   
          void     SetDiffusion(int diffOn = 1)            { fDiffusionOn       = diffOn;           }
          void     SetElAttach(int elOn = 1)               { fElAttachOn        = elOn;             }
          void     SetElAttachProp(float prop)             { fElAttachProp      = prop;             }
          void     SetTimeResponse(int trfOn = 1)          { fTRFOn             = trfOn; ReInit();  }  
          void     SetCrossTalk(int ctOn = 1)              { fCTOn              = ctOn; ReInit();   }
          void     SetPadCoupling(float v)                 { fPadCoupling       = v;                }
          void     SetTimeCoupling(float v)                { fTimeCoupling      = v;                }
          void     SetTimeStruct(bool tsOn = 1)            { fTimeStructOn      = tsOn;             }
          void     SetPadResponse(int prfOn = 1)           { fPRFOn             = prfOn;            }
          void     SetNTimeBins(int ntb)                   { fNTimeBins         = ntb;              }
          void     SetNTBoverwriteOCDB(bool over = kTRUE)  { fNTBoverwriteOCDB  = over;             }

          float  GetGasGain() const                        { return fGasGain;                       }
          float  GetNoise() const                          { return fNoise;                         }
          float  GetChipGain() const                       { return fChipGain;                      }
          float  GetADCoutRange() const                    { return fADCoutRange;                   }
          float  GetADCinRange() const                     { return fADCinRange;                    }
          int    GetADCbaseline() const                    { return fADCbaseline;                   }
          float  GetTRFlo() const                          { return fTRFlo;                         }
          float  GetTRFhi() const                          { return fTRFhi;                         }
          float  GetPadCoupling() const                    { return fPadCoupling;                   }
          float  GetTimeCoupling() const                   { return fTimeCoupling;                  }
          int    GetNTimeBins() const                      { return fNTimeBins;                     }
          bool   GetNTBoverwriteOCDB() const               { return fNTBoverwriteOCDB;              }

          bool   DiffusionOn() const                       { return fDiffusionOn;                   }
          bool   ElAttachOn() const                        { return fElAttachOn;                    } 
          float  GetElAttachProp() const                   { return fElAttachProp;                  }
          bool   TRFOn() const                             { return fTRFOn;                         }
          bool   CTOn() const                              { return fCTOn;                          }
          bool   TimeStructOn() const                      { return fTimeStructOn;                  }
          bool   PRFOn() const                             { return fPRFOn;                         }

  Double_t TimeResponse(Double_t time) const;  
  Double_t CrossTalk(Double_t time) const; 

  void     ReInit();
  
 protected:

  static TRDSimParam* fgInstance;   //  Instance of this class (singleton implementation)
  static bool          fgTerminated; //  Defines if this class has already been terminated and
                                       //  therefore does not return instances in GetInstance anymore
  
          // Digitization parameter
          float  fGasGain;           //  Gas gain
          float  fNoise;             //  Electronics noise
          float  fChipGain;          //  Electronics gain
  
          float  fADCoutRange;       //  ADC output range (number of channels)
          float  fADCinRange;        //  ADC input range (input charge)
          int    fADCbaseline;       //  ADC intrinsic baseline in ADC channel
  
          int    fDiffusionOn;       //  Switch for the diffusion
  
          int    fElAttachOn;        //  Switch for the electron attachment
          float  fElAttachProp;      //  Propability for electron attachment (for 1m)
  
          int    fTRFOn;             //  Switch for the time response
          float *fTRFsmp;            //! Integrated time response
          int    fTRFbin;            //  Number of bins for the TRF
          float  fTRFlo;             //  Lower boundary of the TRF
          float  fTRFhi;             //  Higher boundary of the TRF
          float  fTRFwid;            //  Bin width of the integrated TRF
  
          int    fCTOn;              //  Switch for cross talk
          float *fCTsmp;             //! Integrated cross talk
  
          float  fPadCoupling;       //  Pad coupling factor
          float  fTimeCoupling;      //  Time coupling factor (image charge of moving ions)
          int    fTimeStructOn;      //  Switch for cell time structure
  
          int    fPRFOn;             //  Switch for the pad response

          int    fNTimeBins;         //  Number of time bins (only used it fNTBoverwriteOCDB = true)
          bool   fNTBoverwriteOCDB;  //  Switch to overwrite number of time bins from PCDB

 private:

  // This is a singleton, constructor is private!  
  TRDSimParam();
  ~TRDSimParam();

  void Init();
  void SampleTRF();
  
  ClassDefNV(TRDSimParam,1)          // The TRD simulation parameters

};
}}
#endif
