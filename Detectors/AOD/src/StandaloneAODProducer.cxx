// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief Produce AOD from reconstructed objects

/// This is a standalone demonstrator producer of the AOD format
/// from MC-RECO. For the moment not a DPL device but in future this
/// promotion is targeted. The executable produces the aod.root file which
/// can be fed into analysis workflows.

#include "Framework/AnalysisDataModel.h"
#include "Framework/TableBuilder.h"
#include "Framework/TableTreeHelpers.h"
#include "Framework/Logger.h"
#include "TFile.h"

#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

#include <vector>

using namespace o2;
using namespace o2::framework;

using GIndex = o2::dataformats::VtxTrackIndex;

template <typename TrackT>
std::vector<TrackT>* fetchTracks(const char* filename, const char* treename, const char* branchname)
{
  TFile file(filename, "OPEN");
  auto tree = (TTree*)file.Get(treename);
  auto br = tree->GetBranch(branchname);
  std::vector<TrackT>* tracks = nullptr;
  br->SetAddress(&tracks);
  br->GetEntry(0);
  return tracks;
}

// add vertices/collisions/tracks
void fillCollisionAndTrackTable()
{
  // open the file for vertices
  TFile f("o2_primary_vertex.root", "OPEN");
  auto t = (TTree*)f.Get("o2sim");

  // fetch the tracks (these names are not following any convention!!)
  auto tpctracks = fetchTracks<o2::tpc::TrackTPC>("tpctracks.root", "events", "Tracks");
  auto itstracks = fetchTracks<o2::its::TrackITS>("o2trac_its.root", "o2sim", "ITSTrack");
  auto itstpctracks = fetchTracks<o2::dataformats::TrackTPCITS>("o2match_itstpc.root", "matchTPCITS", "TPCITS");
  LOG(INFO) << "FOUND " << tpctracks->size() << " TPC tracks";
  LOG(INFO) << "FOUND " << itstracks->size() << " ITS tracks";
  LOG(INFO) << "FOUND " << itstpctracks->size() << " ITCTPC tracks";

  if (t) {
    auto br = t->GetBranch("PrimaryVertex");
    std::vector<o2::dataformats::PrimaryVertex>* vertices = nullptr;
    br->SetAddress(&vertices);
    br->GetEntry(0);

    // this referes to actual tracks
    auto indexbr = t->GetBranch("PVTrackIndices");
    std::vector<GIndex>* vertexTrackIDs = nullptr;
    indexbr->SetAddress(&vertexTrackIDs);
    indexbr->GetEntry(0);

    // this makes the connection of vertex to track indices
    auto v2totrackrefbr = t->GetBranch("PV2TrackRefs");
    std::vector<o2::dataformats::VtxTrackRef>* v2trackref = nullptr;
    v2totrackrefbr->SetAddress(&v2trackref);
    v2totrackrefbr->GetEntry(0);

    if (vertices && vertexTrackIDs) {
      TableBuilder collBuilder;
      auto collCursor = collBuilder.cursor<o2::aod::Collisions>();

      TableBuilder trackBuilder;
      auto trackCursor = trackBuilder.cursor<o2::aod::Tracks>();

      int index = 0;
      for (auto& v : *vertices) {
        //DECLARE_SOA_TABLE(Collisions, "AOD", "COLLISION", o2::soa::Index<>,
        // collision::BCId, collision::PosX, collision::PosY, collision::PosZ,
        // collision::CovXX, collision::CovXY, collision::CovXZ, collision::CovYY, collision::CovYZ, collision::CovZZ,
        // collision::Chi2, collision::NumContrib, collision::CollisionTime, collision::CollisionTimeRes, collision::CollisionTimeMask);
        int BCid = 0;
        auto& cov = v.getCov();
        auto& ts = v.getTimeStamp();

        // TODO: figure out BC + CollisionTimeMask
        collCursor(0, BCid, v.getX(), v.getY(), v.getZ(),
                   cov[0], cov[1], cov[2], cov[3], cov[4], cov[6],
                   v.getChi2(), v.getNContributors(), ts.getTimeStamp(), ts.getTimeStampError(), 1);

        // get the track for each vertex and fill the tracks table
        // now go over tracks via the indices
        auto& trackref = (*v2trackref)[index];
        int start = trackref.getFirstEntryOfSource(0);
        int ntracks = trackref.getEntriesOfSource(0);
        for (int ti = 0; ti < ntracks; ++ti) {
          auto trackindex = (*vertexTrackIDs)[start + ti];

          // now we need to fetch the actual track and fill the table
          const auto source = trackindex.getSource();
          o2::track::TrackParCov* track = nullptr;
          if (source == o2::dataformats::VtxTrackIndex::Source::TPC) {
            track = &((*tpctracks)[trackindex.getIndex()]);
          } else if (source == o2::dataformats::VtxTrackIndex::Source::ITS) {
            track = &((*itstracks)[trackindex.getIndex()]);
          } else if (source == o2::dataformats::VtxTrackIndex::Source::TPCITS) {
            track = &((*itstpctracks)[trackindex.getIndex()]);
          } else {
            LOG(WARNING) << "Unsupported track source";
          }

          //DECLARE_SOA_TABLE_FULL(StoredTracks, "Tracks", "AOD", "TRACK:PAR",
          //                       o2::soa::Index<>, track::CollisionId, track::TrackType,
          //                       track::X, track::Alpha,
          //                       track::Y, track::Z, track::Snp, track::Tgl,
          //                       track::Signed1Pt,
          //                       track::NormalizedPhi<track::RawPhi>,
          //                       track::Px<track::Signed1Pt, track::Snp, track::Alpha>,
          //                       track::Py<track::Signed1Pt, track::Snp, track::Alpha>,
          //                       track::Pz<track::Signed1Pt, track::Tgl>,
          //                      track::Charge<track::Signed1Pt>);

          std::array<float, 3> pxpypz;
          track->getPxPyPzGlo(pxpypz);
          trackCursor(0, index, 0 /* CORRECT THIS */, track->getX(), track->getAlpha(), track->getY(), track->getZ(), track->getSnp(), track->getTgl(),
                      track->getPt() /*CHECK!!*/, track->getPhi(), pxpypz[0], pxpypz[1], pxpypz[2]);
        }
        index++;
      }
      auto colltable = collBuilder.finalize();
      auto tracktable = trackBuilder.finalize();

      f.Close();
      TFile outfile("aod.root", "RECREATE");
      {
        TableToTree t2t(colltable, &outfile, aod::MetadataTrait<o2::aod::Collisions>::metadata::tableLabel());
        t2t.addAllBranches();
        t2t.process();
      }
      {
        TableToTree t2t(tracktable, &outfile, "Tracks" /* aod::MetadataTrait<o2::aod::Tracks>::metadata::tableLabel() */);
        t2t.addAllBranches();
        t2t.process();
      }
    }
  }
}

// TODO: add MCparticles

int main()
{
  fillCollisionAndTrackTable();
  return 0;
}
