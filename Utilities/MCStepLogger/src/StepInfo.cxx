#include <StepInfo.h>
#include <TVirtualMC.h>
#include <TParticle.h>
#include <chrono>
#include <TArrayI.h>

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <iostream>

namespace o2 {

// construct directly using virtual mc
StepInfo::StepInfo(TVirtualMC *mc) {
  // init base time point
  if (stepcounter==-1) {
    starttime = std::chrono::high_resolution_clock::now();
  }
  stepcounter++;
  stepid = stepcounter;

  auto stack = mc->GetStack();
  trackID = stack->GetCurrentTrackNumber();
  pdg = mc->TrackPid();
  auto id = mc->CurrentVolID(copyNo);
  volId=id;

  //auto v1 = gGeoManager->GetVolume(mc->CurrentVolName());
  //if (copyNo != 1) {
    //  std::cerr << "SCHEISSE " << copyNo << "\n";
  //}
  auto v2 = gGeoManager->GetCurrentNavigator()->GetCurrentVolume();
  if (std::strcmp(mc->CurrentVolName(), v2->GetName())!=0){
    std::cerr << "inconsistent state\n";
  }
  
//  std::cerr << "id " << id << " MaxStep " << mc->MaxStep() << " " << mc->TrackStep() << "\n"; 
//  TArrayI procs;
//  mc->StepProcesses(procs);
//  std::cerr << "##id " << id << " ";
//  for (int i=0; i<procs.GetSize(); ++i){
//    std::cerr << TMCProcessName[procs.At(i)] << "\t";
//  }
//  std::cerr << "\n";
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
  //ClassImp(StepInfo);


  MagCallInfo::MagCallInfo(TVirtualMC *mc, float x, float y, float z, float Bx, float By, float Bz) :
   x{x}
  ,y{y}
  ,z{z}
  ,Bx{Bx}
  ,By{By}
  ,Bz{Bz}
  {
    stepcounter++;
    id = stepcounter;
    stepid = StepInfo::stepcounter; // or -1 ??
}

int MagCallInfo::stepcounter = -1;

}
