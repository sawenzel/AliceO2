// Stage-3 probe for the CAD->TGeo Geant integration test (Handoff_IntegrationTest.md).
//
// Loads one converted geom.C exactly the way o2::passive::ExternalModule does
// (o2::base::buildCADVolumeFromMacro), places it into a probe world with the same
// rotation/translation convention as the ExternalModule JSON (RotateX,Y,Z in degrees, then
// translation in cm), closes the geometry, and reports:
//   - the module's world-frame bounding box (the §4 rotation/shift question, answered by
//     measurement rather than taken on faith),
//   - TGeo volume/node counts and the CloseGeometry wall time,
//   - the media used, with volume counts per medium (the Stage-2 material check),
//   - optionally TGeoManager::CheckOverlaps (a HINT only: that checker is measured-unreliable
//     on O2BVHSurfaceSolid, see Handoff §8.3).
//
// Build (NEXT.md "quick kernel probe" pattern):
//   g++ -O2 -o cadModuleProbe cadModuleProbe.cxx $(root-config --cflags --libs) \
//       -I<srcdir>/Detectors/Base/include -L$B/stage/lib -lO2DetectorsBase -lGeom
//
// Usage:
//   cadModuleProbe <geom.C> <tag> [--rot rx ry rz] [--trans x y z] [--overlaps <tol_cm>]

#include <DetectorsBase/CADGeometryUtils.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TStopwatch.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <sys/resource.h>

int main(int argc, char** argv)
{
  if (argc < 3) {
    fprintf(stderr, "usage: %s <geom.C> <tag> [--rot rx ry rz] [--trans x y z] [--overlaps tol]\n", argv[0]);
    return 2;
  }
  const std::string macro = argv[1];
  const std::string tag = argv[2];
  double rot[3] = {0, 0, 0}, trans[3] = {0, 0, 0};
  double overlapTol = -1.0;
  for (int i = 3; i < argc; ++i) {
    if (!strcmp(argv[i], "--rot") && i + 3 < argc) {
      for (int k = 0; k < 3; ++k) {
        rot[k] = atof(argv[++i]);
      }
    } else if (!strcmp(argv[i], "--trans") && i + 3 < argc) {
      for (int k = 0; k < 3; ++k) {
        trans[k] = atof(argv[++i]);
      }
    } else if (!strcmp(argv[i], "--overlaps") && i + 1 < argc) {
      overlapTol = atof(argv[++i]);
    }
  }

  auto geom = new TGeoManager("probe", "cadModuleProbe world");
  auto vacMat = new TGeoMaterial("Vacuum", 0, 0, 0);
  auto vacMed = new TGeoMedium("Vacuum", 1, vacMat);
  auto world = geom->MakeBox("world", vacMed, 4000., 4000., 4000.);
  geom->SetTopVolume(world);

  TStopwatch swBuild;
  swBuild.Start();
  auto top = o2::base::buildCADVolumeFromMacro(macro, tag);
  swBuild.Stop();
  if (!top) {
    fprintf(stderr, "FAILED to build module from %s\n", macro.c_str());
    return 1;
  }

  // the ExternalModule placement convention, exactly (makePlacementFromJSON)
  auto combi = new TGeoCombiTrans();
  combi->RotateX(rot[0]);
  combi->RotateY(rot[1]);
  combi->RotateZ(rot[2]);
  combi->SetDx(trans[0]);
  combi->SetDy(trans[1]);
  combi->SetDz(trans[2]);
  world->AddNode(top, 1, combi);

  TStopwatch swClose;
  swClose.Start();
  geom->CloseGeometry();
  swClose.Stop();

  // world-frame bounding box of the placed module: transform the 8 corners of the module's
  // own (assembly-computed) bbox
  auto bb = dynamic_cast<TGeoBBox*>(top->GetShape());
  double lo[3] = {1e30, 1e30, 1e30}, hi[3] = {-1e30, -1e30, -1e30};
  if (bb) {
    const double* o = bb->GetOrigin();
    double dx = bb->GetDX(), dy = bb->GetDY(), dz = bb->GetDZ();
    for (int cx = -1; cx <= 1; cx += 2) {
      for (int cy = -1; cy <= 1; cy += 2) {
        for (int cz = -1; cz <= 1; cz += 2) {
          double local[3] = {o[0] + cx * dx, o[1] + cy * dy, o[2] + cz * dz};
          double master[3];
          combi->LocalToMaster(local, master);
          for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], master[k]);
            hi[k] = std::max(hi[k], master[k]);
          }
        }
      }
    }
    printf("module_local_bbox_cm: origin=(%.3f %.3f %.3f) half=(%.3f %.3f %.3f)\n",
           o[0], o[1], o[2], dx, dy, dz);
    printf("module_world_bbox_cm: x=[%.3f, %.3f] y=[%.3f, %.3f] z=[%.3f, %.3f]\n",
           lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
  } else {
    printf("module top volume has no TGeoBBox-derived shape (%s)\n",
           top->GetShape() ? top->GetShape()->ClassName() : "null");
  }

  printf("volumes: %d\n", geom->GetListOfVolumes()->GetEntries());
  printf("nodes_total: %d\n", geom->CountNodes(world, 100));
  printf("build_module_s: %.3f  close_geometry_s: %.3f\n", swBuild.RealTime(), swClose.RealTime());

  // media census: which media are used, by how many volumes (Stage-2 material verification)
  std::map<std::string, int> mediaCount;
  TIter next(geom->GetListOfVolumes());
  while (auto* obj = next()) {
    auto* vol = static_cast<TGeoVolume*>(obj);
    if (vol->IsAssembly()) {
      continue;
    }
    auto* med = vol->GetMedium();
    mediaCount[med ? med->GetName() : "<none>"]++;
  }
  printf("media_census (volumes per medium):\n");
  for (auto& kv : mediaCount) {
    printf("  %-40s %d\n", kv.first.c_str(), kv.second);
  }

  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  printf("max_rss_mb: %.1f\n", ru.ru_maxrss / 1024.0);

  if (overlapTol >= 0) {
    TStopwatch swOv;
    swOv.Start();
    geom->CheckOverlaps(overlapTol);
    swOv.Stop();
    auto* ov = geom->GetListOfOverlaps();
    printf("check_overlaps: tol_cm=%.4f n=%d time_s=%.1f  (HINT ONLY on O2BVHSurfaceSolid, see Handoff §8.3)\n",
           overlapTol, ov ? ov->GetEntries() : -1, swOv.RealTime());
    if (ov) {
      TIter onext(ov);
      int shown = 0;
      while (auto* o = onext()) {
        if (shown++ >= 10) {
          printf("  ... (%d more)\n", ov->GetEntries() - 10);
          break;
        }
        printf("  %s\n", o->GetName());
      }
    }
  }
  // skip TGeo teardown; process exit is cheaper and avoids atexit issues (same reason as o2-sim)
  fflush(stdout);
  _exit(0);
}
