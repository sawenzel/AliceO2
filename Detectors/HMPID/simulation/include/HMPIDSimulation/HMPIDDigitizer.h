#ifndef DETECTORS_HMPID_SIMULATION_INCLUDE_HMPIDSIMULATION_HMPIDDIGITIZER_H_
#define DETECTORS_HMPID_SIMULATION_INCLUDE_HMPIDSIMULATION_HMPIDDIGITIZER_H_

#include "HMPIDBase/Digit.h"
#include "HMPIDSimulation/Detector.h" // for the hit
#include <vector>

namespace o2
{
namespace hmpid
{

class HMPIDDigitizer
{
 public:
  void setEventTime(double timeNS) { mTime = timeNS; }
  void setEventID(int eventID) { mEventID = eventID; }
  void setSrcID(int sID) { mSrcID = sID; }

  // this will process hits and fill the digit vector with digits which are finalized
  void process(std::vector<o2::hmpid::HitType> const&, std::vector<o2::hmpid::Digit>& digit);

 private:
  double mTime = 0.;
  int mEventID = 0;
  int mSrcID = 0;

  // internal buffers for digits
  std::vector<o2::hmpid::Digit> mSummable;
  std::vector<o2::hmpid::Digit> mFinal;

  // other stuff needed for digitizaton

  ClassDefNV(HMPIDDigitizer, 1);
};
}
}

#endif /* DETECTORS_HMPID_SIMULATION_INCLUDE_HMPIDSIMULATION_HMPIDDIGITIZER_H_ */
