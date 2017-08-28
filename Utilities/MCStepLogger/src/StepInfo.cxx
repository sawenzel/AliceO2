#include <StepInfo.h>
#include <TVirtualMC.h>
#include <TParticle.h>
#include <chrono>
#include <TArrayI.h>

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMedium.h>
#include <TDatabasePDG.h>
#include <FairRunSim.h>
#include <FairModule.h>
#include <iostream>

ClassImp(o2::StepInfo);
ClassImp(o2::MagCallInfo);

namespace o2 {

// construct directly using virtual mc
StepInfo::StepInfo(TVirtualMC *mc) {
  // init base time point
  if (stepcounter==-1) {
    starttime = std::chrono::high_resolution_clock::now();
  }
  stepcounter++;
  stepid = stepcounter;
  currentinstance=this;

  eventid = mc->CurrentEvent();
  auto stack = mc->GetStack();
  trackID = stack->GetCurrentTrackNumber();
  pdg = mc->TrackPid();
  pname = TDatabasePDG::Instance()->GetParticle(pdg)->GetName();
  auto id = mc->CurrentVolID(copyNo);
  volId=id;
  primary = stack->GetCurrentTrack()->IsPrimary();
  
  auto geovolume = gGeoManager->GetCurrentVolume();
  volname = geovolume ? geovolume->GetName() : "null";
  mediumname = geovolume ? geovolume->GetMedium()->GetName() : nullptr;

  auto* run = FairRunSim::Instance();
  moduleid = run->GetVolToModulesMap().at(volId);
  auto modulelist = run->GetListOfModules();
  auto module=(FairModule*)modulelist->At(moduleid);
  modulename = module->GetName();
    
  //auto v2 = gGeoManager->GetCurrentNavigator()->GetCurrentVolume();
  //if (strcmp(mc->CurrentVolName(), v2->GetName())!=0){
  //  std::cerr << "inconsistent state\n";
  //}
  
  double xd,yd,zd;
  mc->TrackPosition(xd,yd,zd);
  x=xd;
  y=yd;
  z=zd;
  E=mc->GetStack()->GetCurrentTrack()->Energy();    
  auto now = std::chrono::high_resolution_clock::now();
  cputimestamp = std::chrono::duration_cast<std::chrono::nanoseconds>
                             (now-starttime).count();
  nsecondaries = mc->NSecondaries();
  secondaryprocesses = new int[nsecondaries];
  // for the processes
  for (int i = 0; i < nsecondaries; ++i) {
    secondaryprocesses[i] = mc->ProdProcess(i);
  }

  TArrayI procs;
  mc->StepProcesses(procs);
  nprocessesactive = procs.GetSize();

  // is track entering volume

  // is track exiting volume

  // was track stoped due to energy limit ??
  stopped = mc->IsTrackStop();
}

std::chrono::time_point<std::chrono::high_resolution_clock> StepInfo::starttime;
int StepInfo::stepcounter = -1;
StepInfo* StepInfo::currentinstance = nullptr; 

MagCallInfo::MagCallInfo(TVirtualMC *mc, float ax, float ay, float az, float aBx, float aBy, float aBz) :
  x{ax}
 ,y{ay}
 ,z{az}
 ,Bx{aBx}
 ,By{aBy}
 ,Bz{aBz}
 {
   stepcounter++;
   id = stepcounter;
   stepid = StepInfo::stepcounter;
   // copy the stepinfo
   if (StepInfo::currentinstance) {
     // stepinfo = *StepInfo::currentinstance;
   }
 }

int MagCallInfo::stepcounter = -1;

}
