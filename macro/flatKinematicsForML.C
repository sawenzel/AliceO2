// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "SimulationDataFormat/MCTrack.h"
#include <iostream>
#endif


// the purpose is to prepare completely flat branches (associated to each feature) for ML starting
// with the more AOS layout std::vector<Tracks> that we write in simulation
void flatKinematicsForML(std::string filename = "o2sim.root")
{
  // prepare outfile / outtree
  auto outfile = new TFile("flat_kinematics.root", "RECREATE");
  auto outTree = new TTree("o2simsecondaries", "o2simsecondaries");

  // define variables to put
  float x; outTree->Branch("x", &x);
  float y; outTree->Branch("y", &y);
  float z; outTree->Branch("z", &z);
  float t; outTree->Branch("t", &t);
  float px; outTree->Branch("px", &px);
  float py; outTree->Branch("py", &py);
  float pz; outTree->Branch("pz", &pz);

  float pdg; outTree->Branch("pdg", &pdg);
  float pt; outTree->Branch("pt", &pt);
  float e; outTree->Branch("e", &e);
  float mass; outTree->Branch("mass", &e);
  float vol; outTree->Branch("vol", &vol); // volume where created
  float eta; outTree->Branch("eta", &eta);
  float classification; outTree->Branch("classification", &classification);

  // read the data
  TFile infile(filename.c_str());
  auto intree = static_cast<TTree*>(infile.Get("o2sim"));
  if (!intree) {
    std::cerr << "no input tree found\n"; return;
  }

  std::vector<o2::MCTrack>* mctracks = nullptr;
  auto br = intree->GetBranch("MCTrack");
  br->SetAddress(&mctracks);

  for (int i = 0; i < br->GetEntries(); ++i) {
    br->GetEntry(i);
    std::cout << "adding " << mctracks->size() << " tracks \n";
    for (auto& track : *mctracks) {
      x = track.X();
      y = track.Y();
      z = track.Z();
      t = track.T();
      px = track.PX();
      py = track.PY();
      pz = track.PZ();
      pdg = track.GetPdgCode();
      pt = track.GetPt();
      e = track.GetEnergy();
      mass = track.GetMass();
      vol = track.getVolumeID(); // volume where created
      eta = track.GetRapidity();
      classification = track.getStore();

      outTree->Fill();
    }
  }
  infile.Close();
  // outTree->Write();
  outfile->WriteTObject(outTree);
  outfile->Close();
}
