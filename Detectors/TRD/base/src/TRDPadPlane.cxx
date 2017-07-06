///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Describes a pad plane of a TRD ROC                                       //
//                                                                           //
//  Contains the information on pad postions, pad dimensions,                //
//  tilting angle, etc.                                                      //
//  It also provides methods to identify the current pad number from         //
//  global coordinates.                                                      //
//  The numbering and coordinates should follow the official convention      //
//  (see David Emschermanns note on TRD convention                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include "TRDBase/TRDPadPlane.h"

using namespace o2::trd;

//_____________________________________________________________________________
TRDPadPlane::TRDPadPlane() :
   fLayer(0)
  ,fStack(0)
  ,fLength(0)
  ,fWidth(0)
  ,fLengthRim(0)
  ,fWidthRim(0)
  ,fLengthOPad(0)
  ,fWidthOPad(0)
  ,fLengthIPad(0)
  ,fWidthIPad(0)
  ,fRowSpacing(0)
  ,fColSpacing(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fTiltingAngle(0)
  ,fTiltingTan(0)
  ,fPadRow(0)
  ,fPadCol(0)
  ,fPadRowSMOffset(0)
  ,fAnodeWireOffset(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
TRDPadPlane::TRDPadPlane(Int_t layer, Int_t stack)
  : fLayer(layer)
  ,fStack(stack)
  ,fLength(0)
  ,fWidth(0)
  ,fLengthRim(0)
  ,fWidthRim(0)
  ,fLengthOPad(0)
  ,fWidthOPad(0)
  ,fLengthIPad(0)
  ,fWidthIPad(0)
  ,fRowSpacing(0)
  ,fColSpacing(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fTiltingAngle(0)
  ,fTiltingTan(0)
  ,fPadRow(0)
  ,fPadCol(0)
  ,fPadRowSMOffset(0)
  ,fAnodeWireOffset(0)
{
  //
  // Constructor
  //

}

//_____________________________________________________________________________
TRDPadPlane::~TRDPadPlane()
{
  //
  // TRDPadPlane destructor
  //

  if (fPadRow) {
    delete [] fPadRow;
    fPadRow = 0;
  }

  if (fPadCol) {
    delete [] fPadCol;
    fPadCol = 0;
  }

}

/*
//_____________________________________________________________________________
void TRDPadPlane::Copy(TObject &p) const
{
  //
  // Copy function
  //

  Int_t iBin = 0;

  ((TRDPadPlane &) p).fLayer           = fLayer;
  ((TRDPadPlane &) p).fStack           = fStack;

  ((TRDPadPlane &) p).fLength          = fLength;
  ((TRDPadPlane &) p).fWidth           = fWidth;
  ((TRDPadPlane &) p).fLengthRim       = fLengthRim;
  ((TRDPadPlane &) p).fWidthRim        = fWidthRim;
  ((TRDPadPlane &) p).fLengthOPad      = fLengthOPad;
  ((TRDPadPlane &) p).fWidthOPad       = fWidthOPad;
  ((TRDPadPlane &) p).fLengthIPad      = fLengthIPad;
  ((TRDPadPlane &) p).fWidthIPad       = fWidthIPad;

  ((TRDPadPlane &) p).fRowSpacing      = fRowSpacing;
  ((TRDPadPlane &) p).fColSpacing      = fColSpacing;

  ((TRDPadPlane &) p).fNrows           = fNrows;
  ((TRDPadPlane &) p).fNcols           = fNcols;

  ((TRDPadPlane &) p).fTiltingAngle    = fTiltingAngle;
  ((TRDPadPlane &) p).fTiltingTan      = fTiltingTan;

  ((TRDPadPlane &) p).fPadRowSMOffset  = fPadRowSMOffset;
  ((TRDPadPlane &) p).fAnodeWireOffset = fAnodeWireOffset;

  if (((TRDPadPlane &) p).fPadRow) {
    delete [] ((TRDPadPlane &) p).fPadRow;
  }
  ((TRDPadPlane &) p).fPadRow = new Double_t[fNrows];
  for (iBin = 0; iBin < fNrows; iBin++) {
    ((TRDPadPlane &) p).fPadRow[iBin] = fPadRow[iBin];
  }                                                                             

  if (((TRDPadPlane &) p).fPadCol) {
    delete [] ((TRDPadPlane &) p).fPadCol;
  }
  ((TRDPadPlane &) p).fPadCol = new Double_t[fNrows];
  for (iBin = 0; iBin < fNrows; iBin++) {
    ((TRDPadPlane &) p).fPadCol[iBin] = fPadCol[iBin];
  }                                                                             

  TObject::Copy(p);

}
*/

//_____________________________________________________________________________
void TRDPadPlane::SetTiltingAngle(Double_t t)
{
  //
  // Set the tilting angle of the pads
  //
 
  fTiltingAngle = t; 
  fTiltingTan   = TMath::Tan(TMath::Pi()/180.0 * fTiltingAngle); 

}

//_____________________________________________________________________________
Int_t TRDPadPlane::GetPadRowNumber(Double_t z) const
{
  //
  // Finds the pad row number for a given z-position in local supermodule system
  //

  Int_t row    = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t middle = 0;

  if ((z > GetRow0()  ) || 
      (z < GetRowEnd())) {

    row = -1;

  }
  else {

    nabove = fNrows + 1;
    nbelow = 0;
    while (nabove - nbelow > 1) {
      middle = (nabove + nbelow) / 2;
      if (z == (fPadRow[middle-1] + fPadRowSMOffset)) {
        row    = middle;
      }
      if (z  > (fPadRow[middle-1] + fPadRowSMOffset)) {
        nabove = middle;
      }
      else {
        nbelow = middle;
      }
    }
    row = nbelow - 1;

  }

  return row;

}

//_____________________________________________________________________________
Int_t TRDPadPlane::GetPadRowNumberROC(Double_t z) const
{
  //
  // Finds the pad row number for a given z-position in local ROC system
  //

  Int_t row    = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t middle = 0;

  if ((z > GetRow0ROC()  ) || 
      (z < GetRowEndROC())) {

    row = -1;

  }
  else {

    nabove = fNrows + 1;
    nbelow = 0;
    while (nabove - nbelow > 1) {
      middle = (nabove + nbelow) / 2;
      if (z == fPadRow[middle-1]) {
        row    = middle;
      }
      if (z  > fPadRow[middle-1]) {
        nabove = middle;
      }
      else {
        nbelow = middle;
      }
    }
    row = nbelow - 1;

  }

  return row;

}

//_____________________________________________________________________________
Int_t TRDPadPlane::GetPadColNumber(Double_t rphi) const
{
  //
  // Finds the pad column number for a given rphi-position
  //

  Int_t col    = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t middle = 0;

  if ((rphi < GetCol0()  ) || 
      (rphi > GetColEnd())) {

    col = -1;

  }
  else {

    nabove = fNcols;
    nbelow = 0;
    while (nabove - nbelow > 1) {
      middle = (nabove + nbelow) / 2;
      if (rphi == fPadCol[middle]) {
        col    = middle;
      }
      if (rphi  > fPadCol[middle]) {
        nbelow = middle;
      }
      else {
        nabove = middle;
      }
    }
    col = nbelow;

  }

  return col;

}

ClassImp(TRDPadPlane)
