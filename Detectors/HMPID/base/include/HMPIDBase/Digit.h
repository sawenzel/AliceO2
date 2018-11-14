#ifndef DETECTORS_HMPID_BASE_INCLUDE_HMPIDBASE_DIGIT_H_
#define DETECTORS_HMPID_BASE_INCLUDE_HMPIDBASE_DIGIT_H_

#include "CommonDataFormat/TimeStamp.h"
#include "HMPIDBase/Hit.h" // for hit
#include "HMPIDBase/Param.h" // for param
#include "TMath.h"

namespace o2
{
namespace hmpid
{

/// \class Digit
/// \brief HMPID digit implementation
using DigitBase = o2::dataformats::TimeStamp<double>;
class Digit : public DigitBase
{
 public:
  Digit() = default;
  Digit(float charge) : mQ(charge) {}
  float getCharge() const { return mQ; }

  Digit(HitType const& hit, int pad)
  {
    int pc, px, py;
    float localX;
    float localY;

    // chamber number is in detID
    const auto chamber = hit.GetDetectorID();

    double tmp[3] = { hit.GetX(), hit.GetY(), hit.GetZ() };
    o2::hmpid::Param::Instance()->Mars2Lors(chamber, tmp, localX, localY);

    o2::hmpid::Param::Lors2Pad( localX, localY, pc, px, py);

    // calculate pad id
    mPad = o2::hmpid::Param::Abs(chamber, pc, px, py);

    // calculate charge
    mQ = /*fQHit **/ InMathieson(localX, localY);
  }

 private:
  float mQ = 0.;
  int mPad = 0.; // -1 indicates invalid digit

  float LorsX() const { return Param::LorsX(Param::A2P(mPad), Param::A2X(mPad)); } //center of the pad x, [cm]
  float LorsY() const { return Param::LorsY(Param::A2P(mPad), Param::A2Y(mPad)); } //center of the pad y, [cm]

  float IntPartMathiX(float x) const
  {
    // Integration of Mathieson.
    // This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
    // Arguments: x,y- position of the center of Mathieson distribution
    //  Returns: a charge fraction [0-1] imposed into the pad
    auto shift1 = -LorsX() + 0.5 * o2::hmpid::Param::SizePadX();
    auto shift2 = -LorsX() - 0.5 * o2::hmpid::Param::SizePadX();

    auto ux1 = o2::hmpid::Param::SqrtK3x() * TMath::TanH(o2::hmpid::Param::K2x() * (x + shift1) / o2::hmpid::Param::PitchAnodeCathode());
    auto ux2 = o2::hmpid::Param::SqrtK3x() * TMath::TanH(o2::hmpid::Param::K2x() * (x + shift2) / o2::hmpid::Param::PitchAnodeCathode());

    return o2::hmpid::Param::K4x() * (TMath::ATan(ux2) - TMath::ATan(ux1));
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Double_t IntPartMathiY(Double_t y) const
  {
    // Integration of Mathieson.
    // This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
    // Arguments: x,y- position of the center of Mathieson distribution
    //  Returns: a charge fraction [0-1] imposed into the pad
    Double_t shift1 = -LorsY() + 0.5 * o2::hmpid::Param::SizePadY();
    Double_t shift2 = -LorsY() - 0.5 * o2::hmpid::Param::SizePadY();

    Double_t uy1 = o2::hmpid::Param::SqrtK3y() * TMath::TanH(o2::hmpid::Param::K2y() * (y + shift1) / o2::hmpid::Param::PitchAnodeCathode());
    Double_t uy2 = o2::hmpid::Param::SqrtK3y() * TMath::TanH(o2::hmpid::Param::K2y() * (y + shift2) / o2::hmpid::Param::PitchAnodeCathode());

    return o2::hmpid::Param::K4y() * (TMath::ATan(uy2) - TMath::ATan(uy1));
  }

  float InMathieson(float localX, float localY) const
  {
    return 4. * IntPartMathiX(localX) * IntPartMathiY(localY);
  }

  ClassDefNV(Digit, 1);
};

}
}

#endif /* DETECTORS_HMPID_BASE_INCLUDE_HMPIDBASE_DIGIT_H_ */
