#include <StepInfo.h>
#include <TVirtualMC.h>
#include <TParticle.h>

namespace o2 {

  // construct directly using virtual mc
  StepInfo::StepInfo(TVirtualMC *mc) {
    stepid = stepcounter++;

    auto stack = mc->GetStack();
    trackID = stack->GetCurrentTrackNumber();
    pdg = mc->TrackPid();
    int copyNo;
    auto id = mc->CurrentVolID(copyNo);
    volId=id;

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

    nsecondaries = mc->NSecondaries();
    secondaryprocesses = new int[nsecondaries];
    // for the processes
    for (int i = 0; i < nsecondaries; ++i) {
      secondaryprocesses[i] = mc->ProdProcess(i);
    }
  }

  int StepInfo::stepcounter = 0;
//  ClassImp(StepInfo);
}
