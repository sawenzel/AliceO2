// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_ZDC_GEOMETRY_H
#define ALICEO2_ZDC_GEOMETRY_H

namespace o2
{
namespace zdc
{

class Geo
{
 public:
 
  static constexpr ZNDIMENSION[3] = {3.52, 3.52, 50.};  
  static constexpr ZNFIBRE[3] = {0., 0.01825, 50.};
  static constexpr ZNSECTORS[2] = {2, 2};
  static constexpr ZNDIVISION[2] = {11, 11};
  static constexpr ZNGROOVES[3] = {0.03, 0.03, 50.}; 
  static constexpr ZNCPOSITION[3] = {0., 0., -11261.3}; 
  static constexpr ZNAPOSITION[3] = {0., 0., 11259.8};  
  static constexpr ZNSUPPORTBASE[3] = {6.3, 4.57, 71.2}; //Basement of ZN table (thick one)
  static constexpr ZNSUPPORTBASEPOS[3] = {0., -14., 21.2};
  static constexpr ZNSUPPORTSCINTILH[3] = {4.32-0.8, 0.8, 50.}; //Scintillator container: top&bottom
  static constexpr ZNSUPPORTSCINTILV[3] = {0.8, 1.955, 50.}; //Scintillator container: sides
  static constexpr ZNSUPPORTWALLSUD[3] = {3.52,1.,50.}; //Top and bottom walls
  static constexpr ZNSUPPORTWALLSLAT[3] = {0.4, 5.52, 50.}; //Side walls
  //
  static constexpr ZPDIMENSION[3] = {11.2, 6., 75.};    
  static constexpr ZPFIBRE[3] = {0., 0.0275, 75.};
  static constexpr ZPSECTORS[2] = {4, 1};
  static constexpr ZPDIVISION[2] = {7, 15};
  static constexpr ZPGROOVES[3] = {0.04, 0.04, 75.}; 
  static constexpr ZPCPOSITION[3] = {24.35, 0., -11253.3}; 
  static constexpr ZPAPOSITION[3] = {24.35, 0., 11251.8}; 
  static constexpr ZPSUPPORTBASE1[3] = {12.5, 1.4, 75.};//Bottom basement of ZP table (thinner one)
  static constexpr ZPSUPPORTBASE1POS[3] = {0., -17., 0.};
  static constexpr ZPSUPPORTBASE2[3] = {12.5, 2.5, 75.};//Upper basement of ZP table (thicker one)
  static constexpr ZPSUPPORTBASE2POS[3] = {0., -9., 0.};
  static constexpr ZPSUPPORTBASE3[3] = {1.5, 2.05, 75.};//support table heels (piedini)
  static constexpr ZPSUPPORTWALLBOTTOM[3] = {11.2, 0.25, 75.}; //Bottom wall
  static constexpr ZPSUPPORTWALLUP[3] = {11.2, 1., 75.}; //Top wall
  //static constexpr ZPSUPPORTWALLSLAT[3] = {0.5, 7.25, 75.}; //Side walls (original)
  static constexpr ZPSUPPORTWALLSLAT[3] = {0.5, 6., 75.}; //Side walls (modified)
 
  //
  static constexpr ZEMDIMENSION[6] = {20.62/2., 3.5, 3.5, 45., 0., 0.};    
  static constexpr ZEMFIBRE[3] = {0., 0.0275, 7.};
  static constexpr ZEMDIVISION[3] = {92, 0, 20};
  static constexpr ZEMPBSLICE[6] = {0.15*(TMath::Sqrt(2.)), 3.5, 3.5, 45., 0., 0.};
  static constexpr ZEMVOIDLAYER[6] = {(20.62/20.)/2., 3.5, 3.5, 45., 0., 0.};
  static constexpr ZEMPOSITION[3] = {9.69, 0., 760.};
  static constexpr ZEMSUPPORTTABLE[3] =  {55./2., 1.5/2., 110./2.};
  static constexpr ZEMSUPPORTBOX[6] = {10.5/2., 100./2., 95./2., 0.25/2., 2./2., 2./2.};
  static constexpr ZEMSUPPORT1[3] = {15./2, 3./2., 95./2.}; //support table 
  static constexpr ZEMSUPPORT2[3] = {2./2, 5./2., 95./2.};  //support table heels (piedini)
  static constexpr ZEMSUPPORT3[3] = {3.5, 2./2., 20./2.};   //screens around ZEM
  static constexpr ZEMSUPPORT4[6] = {20./2., 3.5, 1.5/2., 45., 0., 0.}; //detector box walls (side)
  static constexpr ZEMWALLH[3] = {10.5/2., bthickness[1], 95./2.}; //box walls
  static constexpr ZEMWALLVFWD[3] = {10.5/2., (100.-2.)/2., 0.2};
  static constexpr ZEMWALLVBKW[3] = {10.5/2., (100.-2.)/2., 2./2.};
  static constexpr ZEMWALLVSIDE[3] = {0.25/2., (100.-2.)/2., (95.-2.)/2.};
  //
  static constexpr PIPEROTANGLE1[6] = {90.-1.0027, 0., 90., 90., 1.0027, 180.};
  static constexpr PIPEROTANGLE2[6] = {90.+1.0027, 0., 90., 90., 1.0027, 0.};
  static constexpr PIPEROTANGLE3[6] = {90.-1.8934,0.,90.,90.,1.8934,180.};
  static constexpr PIPEROTANGLE4[6] = {90.-3.8,0.,90.,90.,3.8,180.};
  static constexpr PIPEROTANGLE5[6] = {90.+9.8,0.,90.,90.,9.8,0.};
  static constexpr ZNCROTANGLE[6] = {90., 180., 90., 90., 180., 0.};
  static constexpr ZEMROTANGLE1[6] = {0.,0.,90.,90.,-90.,0.};
  static constexpr ZEMROTANGLE2[6] = {180.,0.,90.,45.+90.,90.,45.};
  //
  static constexpr ZBEGINCSIDE = -1947.2; 
  static constexpr ZBEGINASIDE = 1910.2;  
  static constexpr ZENDMOTHERVOLUME = 13500.;
  static constexpr ZCOMPENSATOR = -1974.;
  static constexpr ZCOMPMAGFIELD = 1972.5;
  static constexpr ZITMAGFIELD = 2296.5;
  static constexpr ZD1MAGFIELD = 5838.3;
  static constexpr ZD2MAGFIELD = 12167.8;

 private:
  static void Init(); 

  static Int_t getDetector(const Float_t* pos);
  static Int_t getSector(const Float_t* pos);

  static void fromGlobalToSector(Float_t* pos, Int_t isector); // change coords to Sector reference

  static Float_t mBeamPipeRotationMatrix[2]; // rotation matrixes

  ClassDefNV(Geometry, 1);
};
}
}

#endif
