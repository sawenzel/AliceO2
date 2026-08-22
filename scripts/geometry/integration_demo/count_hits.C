// External-detector hits from an o2-sim-serial run.
//
//   root -l -b -q 'count_hits.C("o2sim.root")'
//
// o2-sim-serial leaves external-detector hits in the monolithic o2sim.root under branches
// named <NAME>Hit. macro/migrateSimFiles.C only splits off detectors that the GRP marks as
// read out and whose branch names SimTraits knows, and it knows nothing about external
// detectors -- so no o2sim_Hits<DETID>.root appears in serial mode. In parallel mode
// (-j >= 2) O2HitMerger writes them like any other detector's.
R__LOAD_LIBRARY(libO2ExternalDetectors)
#include "ExternalDetectors/Hit.h"

void count_hits(const char* file)
{
  TFile f(file);
  auto* t = (TTree*)f.Get("o2sim");
  if (!t) { printf("NOTREE\n"); return; }
  for (auto* o : *t->GetListOfBranches()) {
    TString bn = o->GetName();
    if (!bn.EndsWith("Hit")) continue;
    std::vector<o2::ext::Hit>* v = nullptr;
    t->SetBranchAddress(bn, &v);
    long n = 0;
    double rmin = 1e30, rmax = -1e30, zmin = 1e30, zmax = -1e30, edep = 0;
    for (long i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      if (!v) continue;
      n += v->size();
      for (auto& h : *v) {
        double r = std::hypot(h.GetX(), h.GetY());
        rmin = std::min(rmin, r); rmax = std::max(rmax, r);
        zmin = std::min(zmin, (double)h.GetZ()); zmax = std::max(zmax, (double)h.GetZ());
        edep += h.GetEnergyLoss();
      }
    }
    printf("HITS %-8s n=%6ld  r[%8.3f,%8.3f] z[%9.3f,%9.3f] sumEdep=%.6g GeV\n",
           bn.Data(), n, rmin, rmax, zmin, zmax, edep);
    t->ResetBranchAddresses();
  }
}
