#include "Generators/BoxGenerator.h"
#include "TRandom.h"
#include "TDatabasePDG.h"

using namespace o2::eventgen;

/*
if (fMult <= 0)
   59         return kFALSE;
   60     TDatabasePDG *pid = TDatabasePDG::Instance();
   61     TParticlePDG *p = pid->GetParticle(fPDGType);
   62     if (p != nullptr) {
   63         LOG(info) << this->ClassName() << ": particle with PDG =" << GetPDGType() << " Found";
   64         fPDGMass = p->Mass();
   65     }
   66     return kTRUE;
*/

double GetPDGMass(int pdg) {
  static TDatabasePDG *pid = TDatabasePDG::Instance();
  TParticlePDG *p = pid->GetParticle(pdg);
  if (p != nullptr) {
    // LOG(info) << this->ClassName() << ": particle with PDG =" << GetPDGType() << " Found";
    return p->Mass(); // fPDGMass = p->Mass();
  }
  // LOG(warn) << "pdg not known";
  return 0.;
}

TParticle o2::eventgen::BoxGenerator::sampleParticle() const {
    // Primary particles are distributed uniformly along
    // those kinematics variables which were limitted by setters.
    // if SetCosTheta() function is used, the distribution will be uniform in
    // cos(theta)

    static double mass = GetPDGMass(mPDG);
   
    double pabs = 0, phi, pt = 0, theta = 0, eta, y, mt, px, py, pz = 0;
    phi = gRandom->Uniform(mPhiMin, mPhiMax) * TMath::DegToRad();
    if (mPRangeIsSet) {
        pabs = gRandom->Uniform(mPMin, mPMax);
    } else if (mPtRangeIsSet) {
        pt = gRandom->Uniform(mPtMin, mPtMax);
    }
    if (mThetaRangeIsSet) {
        if (mCosThetaIsSet)
            theta = acos(gRandom->Uniform(cos(mThetaMin * TMath::DegToRad()), cos(mThetaMax * TMath::DegToRad())));
        else {
            theta = gRandom->Uniform(mThetaMin, mThetaMax) * TMath::DegToRad();
        }
    } else if (mEtaRangeIsSet) {
        eta = gRandom->Uniform(mEtaMin, mEtaMax);
        theta = 2 * TMath::ATan(TMath::Exp(-eta));
    } else if (mYRangeIsSet) {
        y = gRandom->Uniform(mYMin, mYMax);
        mt = TMath::Sqrt(mass * mass + pt * pt);
        pz = mt * TMath::SinH(y);
    }
   
    if (mThetaRangeIsSet || mEtaRangeIsSet) {
        if (mPRangeIsSet) {
            pz = pabs * TMath::Cos(theta);
            pt = pabs * TMath::Sin(theta);
        } else if (mPtRangeIsSet) {
            pz = pt / TMath::Tan(theta);
        }
    }
    px = pt * TMath::Cos(phi);
    py = pt * TMath::Sin(phi);
   
    double vx = 0., vy = 0., vz = 0.;
    double etot = 1.;
    return TParticle (mPDG, 0 /*status*/, -1 /* mother1 */, -1 /* mother2 */, 
                      -1 /* daughter1 */, -1 /* daughter2 */, px, py, pz, etot, vx, vy, vz, 0. /*time*/);

}