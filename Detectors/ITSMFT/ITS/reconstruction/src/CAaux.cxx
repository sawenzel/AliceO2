// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITSReconstruction/CAaux.h"
#include "DetectorsBase/Constants.h"

using namespace o2::ITS::CA;
using o2::Base::Constants::kPI;
using std::array;

Cell::Cell(int xx,int yy, int zz, int dd0, int dd1, float curv, array<float,3> n)
  : m1OverR{curv},
  md0{dd0},
  md1{dd1},
  mN{n[0],n[1],n[2]},
  mVector{xx,yy,zz,1u} {
}

bool Cell::Combine(Cell &neigh, int idd) {
  // From outside inward
  if (this->y() == neigh.z() && this->x() == neigh.y()) { // Cells sharing two points
    mVector.push_back(idd);
    if (neigh.getLevel() + 1 > getLevel()) {
      setLevel(neigh.getLevel() + 1u);
    }
    return true;
  }
  return false;
}

Track::Track(float x, float a, array<float,Base::Track::kNParams> p, array<float,Base::Track::kCovMatSize> c, int *cl) :
  mT{x,a,p,c},
  mCl{},
  mLabel{-1},
  mChi2{0.f} {
    for (int i = 0; i < 7; ++i) mCl[i] = cl[i];
  }

bool Track::Update(const Cluster &cl) {
  array<float,2> p{cl.y,cl.z};
  const float dChi2 = mT.getPredictedChi2(p,cl.cov);
  if (!mT.Update(p,cl.cov)) return false;
  else mChi2 += dChi2;
  return true;
}

bool Track::GetPhiZat(float r, float bfield,float &phi, float &z) const {
  float rp4=mT.getCurvature(bfield);

  float xt=mT.getY(), yt=mT.getY();
  float x = 0.f, y = 0.f;
  float sn=sin(mT.getAlpha()), cs=cos(mT.getAlpha());
  float a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  sn=rp4*xt - mT.getSnp(); cs=rp4*yt + sqrt((1.- mT.getSnp())*(1.+mT.getSnp()));
  a=2*(xt*mT.getSnp() - yt*sqrt((1.-mT.getSnp())*(1.+mT.getSnp())))-rp4*(xt*xt + yt*yt);
  float d =  -a/(1 + sqrt(sn*sn + cs*cs));

  if (fabs(d) > r) {
    if (r>1e-1) return false;
    r = fabs(d);
  }

  float rcurr=sqrt(mT.getX()*mT.getX() + mT.getY()*mT.getY());
  float phicurr=mT.getPhi();

  if (mT.getX()>=0.) {
    phi=phicurr+asin(d/r)-asin(d/rcurr);
  } else {
    phi=phicurr+asin(d/r)+asin(d/rcurr)-kPI;
  }

  //return a phi in [0,2pi
  if (phi<0.) phi+=2.*kPI;
  else if (phi>=2.*kPI) phi-=2.*kPI;
  z=mT.getZ()+mT.getTgl()*(sqrt((r-d)*(r+d))-sqrt((rcurr-d)*(rcurr+d)));

  return true;
}

