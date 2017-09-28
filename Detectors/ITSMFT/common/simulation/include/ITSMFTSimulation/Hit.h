// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Hit.h
/// \brief Definition of the ITSMFT Hit class

#ifndef ALICEO2_ITSMFT_POINT_H_
#define ALICEO2_ITSMFT_POINT_H_

#include "SimulationDataFormat/BaseHits.h"     // for BasicXYZEHit
#include "SimulationDataFormat/MCCompLabel.h"
#include "Rtypes.h"       // for Bool_t, Double_t, Int_t, Double32_t, etc
#include "TVector3.h"     // for TVector3
#include <iostream>

namespace o2 {
namespace ITSMFT {

class Hit : public o2::BasicXYZEHit<Float_t,Float_t>
{

  public:
    using Label = o2::MCCompLabel;
    enum HitStatus_t
    {
        kTrackEntering = 0x1,
        kTrackInside   = 0x1<<1,
        kTrackExiting  = 0x1<<2,
        kTrackOut      = 0x1<<3,
        kTrackStopped  = 0x1<<4,
        kTrackAlive    = 0x1<<5
    };

    /// Default constructor
    Hit() = default;

    /// Class Constructor
    /// \param trackID Index of MCTrack
    /// \param detID Detector ID
    /// \param startPos Coordinates at entrance to active volume [cm]
    /// \param pos Coordinates to active volume [cm]
    /// \param mom Momentum of track at entrance [GeV]
    /// \param endTime Time at entrance [ns]
    /// \param time Time since event start [ns]
    /// \param eLoss Energy deposit [GeV]
    /// \param startStatus: status at entrance
    /// \param endStatus: status at exit
    inline Hit(int trackID, unsigned short detID, TVector3 startPos, TVector3 pos, TVector3 mom, double startE,
		 double endTime, double eLoss,unsigned char statusStart, unsigned char status);


    // Entrance position getters
    Point3D<Float_t> getPosStart() const { return mPosStart; }
    Float_t getStartX() const { return mPosStart.X(); }
    Float_t getStartY() const { return mPosStart.Y(); }
    Float_t getStartZ() const { return mPosStart.Z(); }  
    template<typename F> void getStartPosition(F &x, F &y, F &z) const
    {
      x = getStartX();
      y = getStartY();
      z = getStartZ();
    }
    // momentum getters
    Vector3D<Float_t> getMomentum() const { return mMomentum; }
    Vector3D<Float_t>& getMomentum()      { return mMomentum; }
    Float_t getPx() const { return mMomentum.X(); }
    Float_t getPy() const { return mMomentum.Y(); }
    Float_t getPz() const { return mMomentum.Z(); }
    Float_t getE()  const { return mE; }
    Float_t getTotalEnergy() const { return getE(); }
    
    UChar_t getStatusEnd()   const  { return mTrackStatusEnd; }
    UChar_t getStatusStart() const  { return mTrackStatusStart; }

    Bool_t isEntering()      const  { return mTrackStatusEnd & kTrackEntering; }
    Bool_t isInside()        const  { return mTrackStatusEnd & kTrackInside; }
    Bool_t isExiting()       const  { return mTrackStatusEnd & kTrackExiting; }
    Bool_t isOut()           const  { return mTrackStatusEnd & kTrackOut; }
    Bool_t isStopped()       const  { return mTrackStatusEnd & kTrackStopped; }
    Bool_t isAlive()         const  { return mTrackStatusEnd & kTrackAlive; }

    Bool_t isEnteringStart() const  { return mTrackStatusStart & kTrackEntering; }
    Bool_t isInsideStart()   const  { return mTrackStatusStart & kTrackInside; }
    Bool_t isExitingStart()  const  { return mTrackStatusStart & kTrackExiting; }
    Bool_t isOutStart()      const  { return mTrackStatusStart & kTrackOut; }
    Bool_t isStoppedStart()  const  { return mTrackStatusStart & kTrackStopped; }
    Bool_t isAliveStart()    const  { return mTrackStatusStart & kTrackAlive; }

    /// Output to screen
    void print(const Option_t *opt) const override;
    friend std::ostream &operator<<(std::ostream &of, const Hit &point)
    {
      of << "-I- Hit: O2its point for track " << point.getTrackID() << " in detector " << point.getDetectorID() << std::endl;
      /*
      of << "    Position (" << point.fX << ", " << point.fY << ", " << point.fZ << ") cm" << std::endl;
      of << "    Momentum (" << point.fPx << ", " << point.fPy << ", " << point.fPz << ") GeV" << std::endl;
      of << "    Time " << point.fTime << " ns,  Length " << point.fLength << " cm,  Energy loss "
      << point.fELoss * 1.0e06 << " keV" << std::endl;
      */
      return of;
    }

    void setSrcEvID(int srcID, int evID) {
      /// RS: ATTENTION! this is just a trick until we clarify how the hits from different source are
      // provided and identified. At the moment we just create a combined identifier from eventID
      // and sourceID and store it TEMPORARILY in the cached Point's TObject UniqueID
      SetUniqueID( ( srcID << Label::nbitsEvID ) | evID);
    }

    Label getCombLabel() const {
      /// RS: ATTENTION! this is just a trick until we clarify how the hits from different source are
      // provided and identified. At the moment we just create on the fly the label from the track ID
      // and SrcEv id stored as UniqueID
      int srcID = ( GetUniqueID()>>Label::nbitsEvID ) & Label::maskSrcID;
      int evID = GetUniqueID() & Label::maskEvID;
      return Label( getTrackID(), evID, srcID );
    }
    
  private:
    /// Copy constructor
    Hit(const Hit &point);
    Hit operator=(const Hit &point);
    Vector3D<Float_t> mMomentum;              ///< momentum at entrance
    Point3D<Float_t> mPosStart;               ///< position at entrance (base mPos give position on exit)
    Float_t mE;                               ///< total energy at entrance
    UChar_t mTrackStatusEnd;                  ///< MC status flag at exit
    UChar_t mTrackStatusStart;                ///< MC status at starting point

  ClassDefOverride(Hit, 3)
};

Hit::Hit(int trackID, unsigned short detID, TVector3 startPos, TVector3 endPos, TVector3 startMom,
             double startE,double endTime, double eLoss, unsigned char startStatus, unsigned char endStatus)
  : BasicXYZEHit(endPos.X(),endPos.Y(),endPos.Z(),endTime,eLoss,trackID,detID),
    mMomentum(startMom.Px(),startMom.Py(),startMom.Pz()),
    mPosStart(startPos.X(),startPos.Y(),startPos.Z()),
    mE(startE),
    mTrackStatusEnd(endStatus),
    mTrackStatusStart(startStatus)
{}


}
}

#endif
