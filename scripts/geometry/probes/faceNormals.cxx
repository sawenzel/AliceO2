// Stream L, stage 2: at every interior sample point of every source face, compare the kernel's
// outward normal (O2BVHSurfaceSolid::ComputeNormal with dir = nullptr, i.e. unflipped) with
// OpenCascade's ORIENTED outward normal there.
//
// A face whose two normals are antiparallel has its outward direction inverted in the sidecar.
// Closure and edge identity are both blind to that (a global sign cancels in every count they
// make); DistFromOutside / DistFromInside are not, because both select hits by exactly that sign.

#include "DetectorsBase/O2BVHSurfaceSolid.h"
#include "DetectorsBase/O2SurfaceSolidIO.h"

#include "TGeoManager.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;
using namespace o2::base;

int main(int argc, char** argv)
{
  std::string surfaces, samplesFile;
  bool quiet = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--surfaces") {
      surfaces = argv[++i];
    } else if (a == "--samples") {
      samplesFile = argv[++i];
    } else if (a == "--quiet") {
      quiet = true;
    }
  }

  json doc;
  {
    std::ifstream in(samplesFile);
    if (!in) {
      std::cerr << "cannot open " << samplesFile << "\n";
      return 2;
    }
    in >> doc;
  }

  auto* manager = new TGeoManager("xL3", "face normal probe");
  auto* solid = new O2BVHSurfaceSolid("part");
  if (!LoadSurfaceSolid(surfaces, *solid)) {
    std::cerr << "LoadSurfaceSolid failed for " << surfaces << "\n";
    return 2;
  }
  solid->CloseShape(true);
  (void)manager;

  const int nFaces = doc.at("nFaces").get<int>();
  std::printf("# %s : %d OCCT faces, %d sidecar surfaces %s | %s\n",
              samplesFile.c_str(), nFaces, solid->GetNsurfaces(),
              nFaces == solid->GetNsurfaces() ? "(MATCH)" : "(*** MISMATCH ***)",
              O2BVHSurfaceSolid::GetNavigationReliabilityName(solid->GetNavigationReliability()));

  int flipped = 0, ok = 0, undecided = 0;
  std::vector<int> flippedFaces;
  for (const auto& f : doc.at("faces")) {
    const int idx = f.at("index").get<int>();
    int nFlip = 0, nOk = 0, nFar = 0;
    double worstAlign = 1.;
    for (const auto& s : f.at("samples")) {
      double p[3] = {s.at("p")[0].get<double>(), s.at("p")[1].get<double>(),
                     s.at("p")[2].get<double>()};
      const double nOcct[3] = {s.at("n")[0].get<double>(), s.at("n")[1].get<double>(),
                               s.at("n")[2].get<double>()};
      // The point must actually be on the kernel's boundary, or the comparison is meaningless.
      const bool in = solid->Contains(p);
      if (solid->Safety(p, in ? kTRUE : kFALSE) > 1.e-6) {
        ++nFar;
        continue;
      }
      double nKernel[3];
      solid->ComputeNormal(p, nullptr, nKernel);
      const double align = nKernel[0] * nOcct[0] + nKernel[1] * nOcct[1] + nKernel[2] * nOcct[2];
      if (align < 0.) {
        ++nFlip;
      } else {
        ++nOk;
      }
      worstAlign = std::min(worstAlign, align);
    }
    if (nFlip > 0) {
      ++flipped;
      flippedFaces.push_back(idx);
    } else if (nOk > 0) {
      ++ok;
    } else {
      ++undecided;
    }
    if (!quiet && (nFlip > 0 || nOk == 0)) {
      std::string kinds;
      for (const auto& e : f.at("edges")) {
        kinds += e.at("kind").get<std::string>() + " ";
      }
      std::printf("  face %4d %-10s rev=%d faceTol=%.2e edges=%2d [%s] "
                  "-> flipped=%d ok=%d offSurface=%d worstAlign=%+.4f\n",
                  idx, f.at("type").get<std::string>().c_str(), (int)f.at("reversed").get<bool>(),
                  f.at("faceTol").get<double>(), f.at("nEdges").get<int>(), kinds.c_str(),
                  nFlip, nOk, nFar, worstAlign);
    }
  }
  std::printf("  => faces with an INVERTED outward normal: %d / %d   (consistent %d, "
              "not on the kernel boundary %d)\n", flipped, nFaces, ok, undecided);
  if (flipped) {
    std::printf("     flipped face indices:");
    for (int i : flippedFaces) {
      std::printf(" %d", i);
    }
    std::printf("\n");
  }
  return 0;
}
