#include <DetectorsBase/GenericGDMLDetector.h>
#include <regex>
#include <TVirtualMC.h>
#include <DetectorsBase/Stack.h>

using namespace o2::base;

bool GenericGDMLDetector::ProcessHits(FairVolume*) {
  // we create this function from a ROOT macro
  // mOptions
  //
  std::cout << "GDML: Process Hits function called\n";
  float x, y, z, e;
  fMC->TrackPosition(x, y, z);
  e = fMC->Edep();
  auto t = fMC->TrackTime();
  mHits->emplace_back(Hit_t(x, y, z, t, e, 0, 0));
  auto stack = (o2::data::Stack*)fMC->GetStack();
  stack->addHit(GetDetId()); // tells that current track caused a hit
  return true;
}

void GenericGDMLDetector::addSensitiveVolumes() {
  LOG(info) << "Adding sensitive volumes for GenericGDML detector";
  
  std::regex pattern("^ITSUSensor.*");

  // loop over all logical volumes from TGeoManager
  auto volume_list = gGeoManager->GetListOfVolumes();

  for (int i = 0; i < volume_list->GetEntries(); ++i) {
    auto vol = static_cast<TGeoVolume*>(volume_list->At(i));
    std::string vol_name(vol->GetName());
    if (std::regex_match(vol_name, pattern)) {
      std::cout << "Match " << vol_name << "\n";
      AddSensitiveVolume(vol); // from FairRoot
    }
  }
}

void GenericGDMLDetector::EndOfEvent() { Reset(); }

void GenericGDMLDetector::Register()
{
  // This will create a branch in the output tree called Hit, setting the last
  // parameter to kFALSE means that this collection will not be written to the file,
  // it will exist only during the simulation

  if (FairRootManager::Instance()) {
    FairRootManager::Instance()->RegisterAny(addNameTo("Hit").data(), mHits, kTRUE);
  }
}

void GenericGDMLDetector::Reset()
{
  // if (!o2::utils::ShmManager::Instance().isOperational()) {
  mHits->clear();
  //}
}

ClassImp(GenericGDMLDector);
