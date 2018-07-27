// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#ifndef AliceO2_TPC_Detector_H_
#define AliceO2_TPC_Detector_H_

#include "DetectorsBase/Detector.h"   // for Detector
#include "DetectorsBase/VMCUtilities.h"
#include "Rtypes.h"          // for Int_t, Double32_t, Double_t, Bool_t, etc
#include "TLorentzVector.h"  // for TLorentzVector
#include "TString.h"

#include "TPCSimulation/Point.h"
#include "TPCBase/Sector.h"
#include "TPCBase/ParameterGas.h"
#include "FairVolume.h"

//class FairVolume;  // lines 10-10

namespace o2 {
namespace TPC {

class Detector: public o2::Base::DetImpl<Detector> {

  public:
   /** Local material/media IDs for TPC */
   enum EMedium {
     kAir = 0,
     kDriftGas1 = 1,
     kDriftGas2 = 2,
     kCO2 = 3,
     kDriftGas3 = 20,
     kAl = 4,
     kKevlar = 5,
     kNomex = 6,
     kMakrolon = 7,
     kMylar = 8,
     kTedlar = 9,
     kPrepreg1 = 10,
     kPrepreg2 = 11,
     kPrepreg3 = 12,
     kEpoxy = 13,
     kCu = 14,
     kSi = 15,
     kG10 = 16,
     kPlexiglas = 17,
     kSteel = 18,
     kPeek = 19,
     kAlumina = 21,
     kWater = 22,
     kBrass = 23,
     kEpoxyfm = 24,
     kEpoxy1 = 25,
     kAlumina1 = 26
   };
   /**      Name :  Detector Name
    *       Active: kTRUE for active detectors (ProcessHits() will be called)
    *               kFALSE for inactive detectors
    */
   Detector(Bool_t Active);

   /**      default constructor    */
   Detector();

   /**       destructor     */
   ~Detector() override;

   /**      Initialization of the detector is done here    */
   void Initialize() override;

   /**       this method is called for each step during simulation
    *       (see FairMCApplication::Stepping())
    */
   //     virtual Bool_t ProcessHitsOrig( FairVolume* v=0);
   // Bool_t ProcessHits(FairVolume* v = nullptr) override;

   template <typename VMCEngine = TVirtualMC>
   bool ProcessHitsKernel(FairVolume* v = nullptr);

   /**       Registers the produced collections in FAIRRootManager.     */
   void Register() override;

   /** Get the produced hits */
   std::vector<HitGroup>* getHits(Int_t iColl) const
   {
     if (iColl >= 0 && iColl < Sector::MAXSECTOR) {
       return mHitsPerSectorCollection[iColl];
     }
     return nullptr;
    }

    /** tell the branch names corresponding to hits **/
    std::string getHitBranchNames(int coll) const override;

    /**      has to be called after each event to reset the containers      */
    void   Reset() override;

    /**      Create the detector geometry        */
    void ConstructGeometry() override;

    /**      This method is an example of how to add your own point
     *       of type DetectorPoint to the clones array
    */
    Point* addHit(float x, float y, float z, float time, float nElectrons, float trackID, float detID);
    

    /// Copied from AliRoot - should go to someplace else
    
    /// Empirical ALEPH parameterization of the Bethe-Bloch formula, normalized to 1 at the minimum.
    /// @param bg Beta*Gamma of the incident particle
    /// @param kp* Parameters for the ALICE TPC
    /// @return Bethe-Bloch value in MIP units
    template <typename T>
    T BetheBlochAleph(T bg, T kp1, T kp2, T kp3, T kp4, T kp5);

    /// Copied from AliRoot - should go to someplace else
    /// Function to generate random numbers according to Gamma function 
    /// From Hisashi Tanizaki:
    /// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.158.3866&rep=rep1&type=pdf
    /// Implemented by A. Morsch 14/01/2014    
    /// @k is the mean and variance
    Double_t Gamma(Double_t k);

    
    /** The following methods can be implemented if you need to make
     *  any optional action in your detector during the transport.
    */

    /// Special Geant3? limits and definitions
    /// \todo Check how to deal with this in O2 compared to AliRoot
    /// \todo Discuss in a wider scope
    /// \todo Check correctness of the implementation
    void   SetSpecialPhysicsCuts() override;

    void   EndOfEvent() override;
    void   FinishPrimary() override {;}
    void   FinishRun() override {;}
    void   BeginPrimary() override {;}
    void   PostTrack() override {;}
    void   PreTrack() override {;}
    void   BeginEvent() override {;}

    void SetGeoFileName(const TString file) { mGeoFileName=file;   }
    const TString& GetGeoFileName() const   { return mGeoFileName; }

  private:
    int mHitCounter = 0;
    int mElectronCounter = 0;
    int mStepCounter = 0;

    /// Create the detector materials
    virtual void CreateMaterials();

    /// Geant settings hack
    /// \todo Check if still needed see comment in \ref SetSpecialPhysicsCuts
    void GeantHack();

    /// Construct the detector geometry
    void LoadGeometryFromFile();
    /// Construct the detector geometry
    void ConstructTPCGeometry();

    /** Define the sensitive volumes of the geometry */
    void DefineSensitiveVolumes() override;

    /** container for produced hits */
    std::vector<HitGroup>*  mHitsPerSectorCollection[Sector::MAXSECTOR]; //! container that keeps track-grouped hits per sector

    TString mGeoFileName;                  ///< Name of the file containing the TPC geometry
    // size_t mEventNr;                       //!< current event number
    // Events are not successive in MT mode

    /// copy constructor (used in MT)
    Detector(const Detector& rhs);
    Detector& operator=(const Detector&);

    template <typename Det>
    friend class o2::Base::DetImpl;
    ClassDefOverride(Detector, 1)
};

template<typename T>
inline
T Detector::BetheBlochAleph(T bg, T kp1, T kp2, T kp3, T kp4, T kp5){
  T beta = bg/std::sqrt(static_cast<T>(1.)+ bg*bg);

  T aa = std::pow(beta,kp4);
  T bb = std::pow(static_cast<T>(1.)/bg,kp5);
  bb=std::log(kp3+bb);

  return (kp2-aa-bb)*kp1/aa;
}

template<typename VMCEngine>
inline
bool Detector::ProcessHitsKernel(FairVolume *vol) {
  // This method is called from the MC stepping for the sensitive volume only
  using namespace vmc_helpers;
  mStepCounter++;
  const static ParameterGas& gasParam = ParameterGas::defaultInstance();

  const double trackCharge = TrackCharge<VMCEngine>(fMC);
  if (static_cast<int>(trackCharge) == 0) {
    // set a very large step size for neutral particles
    fMC->SetMaxStep(1.e10);
    return false; // take only charged particles
  }

  // ===| SET THE LENGTH OF THE NEXT ENERGY LOSS STEP |=========================
  // (We do this first so we can skip out of the method in the following
  // Eloss->Ionization part.)
  //
  // In all cases we have multiple collisions and we use 2mm (+ some
  // random shift to avoid binning effects), which was tuned for GEANT4, see
  // https://indico.cern.ch/event/316891/contributions/732168/

  const double rnd = fMC->GetRandom()->Rndm();
  fMC->SetMaxStep(0.2 + (2. * rnd - 1.) * 0.05); // 2 mm +- rndm*0.5mm step


  // ===| check active sector |=================================================
  //
  // Get the sector ID and check if the sector is active
  static thread_local TLorentzVector position;
  TrackPosition<VMCEngine>(fMC, position);
  // for processing reasons in the digitizer, the sectors are shifted by -10deg, so sector 0 will be
  //   from -10 - +10 deg instead of 0-20 and so on
  const int sectorID = static_cast<int>(Sector::ToShiftedSector(position.X(), position.Y(), position.Z()));
  // const int sectorID = static_cast<int>(Sector::ToSector(position.X(), position.Y(), position.Z()));
  // TODO: Temporary hack to process only one sector
  // if (sectorID != 0) return kFALSE;

  // ---| momentum and beta gamma |---
  static TLorentzVector momentum; // static to make avoid creation/deletion of this expensive object
  TrackMomentum<VMCEngine>(fMC, momentum);

  const float time = TrackTime<VMCEngine>(fMC) * 1.0e9;
  const int trackID = mO2Stack->GetCurrentTrackNumber();
  const int detID = vol->getMCid();
  if (IsTrackEntering<VMCEngine>(fMC) || IsTrackExiting<VMCEngine>(fMC)) {
     mO2Stack->addTrackReference(
	 o2::TrackReference(position.X(), position.Y(), position.Z(), momentum.X(), momentum.Y(),
	 momentum.Z(), TrackLength<VMCEngine>(fMC), time, trackID, GetDetId()));
  }

  // ===| CONVERT THE ENERGY LOSS TO IONIZATION ELECTRONS |=====================
  //
  // The energy loss is implemented directly below and taken GEANT3,
  // ILOSS model 5 (in gfluct.F), which gives
  // the energy loss in a single collision (NA49 model).
  // TODO: Add discussion about drawback

  Int_t numberOfElectrons = 0;
  // ---| Stepsize in cm |---
  const double stepSize = TrackStep<VMCEngine>(fMC);
  double betaGamma = momentum.P() / TrackMass<VMCEngine>(fMC);
  betaGamma = TMath::Max(betaGamma, 7.e-3); // protection against too small bg
  // ---| number of primary ionisations per cm |---
  const double primaryElectronsPerCM =
      gasParam.getNprim() * BetheBlochAleph(static_cast<float>(betaGamma), gasParam.getBetheBlochParam(0),
                                            gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2),
                                            gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));

  // ---| mean number of collisions and random for this event |---
  const double meanNcoll = stepSize * trackCharge * trackCharge * primaryElectronsPerCM;
  const int nColl = static_cast<int>(fMC->GetRandom()->Poisson(meanNcoll));

  // Variables needed to generate random powerlaw distributed energy loss
  const double alpha_p1 = 1. - gasParam.getExp(); // NA49/G3 value
  const double oneOverAlpha_p1 = 1. / alpha_p1;
  const double eMin = gasParam.getIpot();
  const double eMax = gasParam.getEend();
  const double kMin = TMath::Power(eMin, alpha_p1);
  const double kMax = TMath::Power(eMax, alpha_p1);
  const double wIon = gasParam.getWion();

  for (Int_t n = 0; n < nColl; n++) {
     // Use GEANT3 / NA49 expression:
     // P(eDep) ~ k * edep^-gasParam.getExp()
     // eMin(~I) < eDep < eMax(300 electrons)
     // k fixed so that Int_Emin^EMax P(Edep) = 1.
    const double rndm = fMC->GetRandom()->Rndm();
    const double eDep = TMath::Power((kMax - kMin) * rndm + kMin, oneOverAlpha_p1);
    int nel_step = static_cast<int>(((eDep - eMin) / wIon) + 1);
    nel_step = TMath::Min(nel_step, 300); // 300 electrons corresponds to 10 keV
    numberOfElectrons += nel_step;
  }

  if (numberOfElectrons <= 0) {
	// Could maybe be smaller than 0 due to the Gamma function
    return false;
  }

  // ADD HIT
  static thread_local int oldTrackId = trackID;
  static thread_local int oldDetId = detID;
  static thread_local int groupCounter = 0;
  static thread_local int oldSectorId = sectorID;

  //  a new group is starting -> put it into the container
  static thread_local HitGroup* currentgroup = nullptr;
  if (groupCounter == 0) {
    mHitsPerSectorCollection[sectorID]->emplace_back(trackID);
    currentgroup = &(mHitsPerSectorCollection[sectorID]->back());
  }
  if (trackID == oldTrackId && oldSectorId == sectorID) {
    groupCounter++;
    mHitCounter++;
    mElectronCounter += numberOfElectrons;
    currentgroup->addHit(position.X(), position.Y(), position.Z(), time, numberOfElectrons);
  }
  // finish group
  else {
   oldTrackId = trackID;
   oldSectorId = sectorID;
   groupCounter = 0;
  }
  return true;
}

}
}

#endif // AliceO2_TPC_Detector_H_
