// Diagnostic probe (Stream L): per-ray transport diff against the OpenCascade crossing list, with
// the kernel's own hit list along the same ray attached to every reference crossing.
//
// Uses the SAME stepping loop the benchmark runs (XRayTransport.h), so this is not a second
// implementation of the idea; the summary line must reproduce the benchmark's per-part counts.
//
// For every crossing OpenCascade found we record what the kernel's ray/surface intersector says
// at that distance: whether it produced a hit there at all, the sign of dot(outward normal, ray
// direction) it attached to it, and whether the hit sat in the patch's on-trim-boundary band.
// That is what separates "the intersection is never computed" (a solver or trim defect) from "it
// is computed and then filtered by its normal" (an orientation defect).
//
// The python side (attributeToFaces.py) joins this against OCCT's own face index.

#include "XRayTransport.h"
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
using namespace o2::base::xray;

int main(int argc, char** argv)
{
  std::string surfaces, crossingsFile, outFile;
  bool alsoClean = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--surfaces") {
      surfaces = argv[++i];
    } else if (a == "--crossings") {
      crossingsFile = argv[++i];
    } else if (a == "--out") {
      outFile = argv[++i];
    } else if (a == "--all-rays") {
      alsoClean = true;
    }
  }
  if (surfaces.empty() || crossingsFile.empty() || outFile.empty()) {
    std::cerr << "usage: faceAttrib --surfaces <bin> --crossings <crossings.json> --out <json> "
                 "[--all-rays]\n";
    return 2;
  }

  json doc;
  {
    std::ifstream in(crossingsFile);
    if (!in) {
      std::cerr << "cannot open " << crossingsFile << "\n";
      return 2;
    }
    in >> doc;
  }
  const double tolerance = doc.value("tolerance", 1.e-7);

  auto* manager = new TGeoManager("xL", "attribution probe");
  auto* solid = new O2BVHSurfaceSolid("part");
  if (!LoadSurfaceSolid(surfaces, *solid)) {
    std::cerr << "LoadSurfaceSolid failed for " << surfaces << "\n";
    return 2;
  }
  solid->CloseShape(true);
  (void)manager;

  StepConfig cfg;
  cfg.matchTolerance = std::max(tolerance, 1.e-6);

  json out;
  out["surfaces"] = surfaces;
  out["nSurfaces"] = solid->GetNsurfaces();
  out["modelTolerance"] = solid->GetModelTolerance();
  out["matchTolerance"] = cfg.matchTolerance;
  out["oracleTolerance"] = tolerance;
  out["reliability"] =
    O2BVHSurfaceSolid::GetNavigationReliabilityName(solid->GetNavigationReliability());
  out["navigable"] = solid->IsNavigable();
  json rays = json::array();

  long long nLost = 0, nExtra = 0, nKind = 0, nDisp = 0, nUnterm = 0, nRays = 0, nIdent = 0;

  const auto& rayDocs = doc.at("rays");
  for (size_t ri = 0; ri < rayDocs.size(); ++ri) {
    const auto& r = rayDocs[ri];
    if (r.value("amb", false)) {
      continue;
    }
    Point3D origin{r.at("o")[0].get<double>(), r.at("o")[1].get<double>(), r.at("o")[2].get<double>()};
    Point3D dir{r.at("d")[0].get<double>(), r.at("d")[1].get<double>(), r.at("d")[2].get<double>()};
    const double tMax = r.at("tmax").get<double>();

    std::vector<Crossing> ref;
    const auto& ts = r.at("t");
    const auto& ks = r.at("k");
    for (size_t i = 0; i < ts.size(); ++i) {
      ref.push_back({ts[i].get<double>(), ks[i].get<int>()});
    }

    Robustness stats;
    auto cand = stepWithShapeApi(solid, origin, dir, tMax, cfg, stats);
    ++nRays;
    nUnterm += stats.unterminated;

    // Pair the two lists exactly as XRayTransport.h::compareLists does, but keep WHICH reference
    // crossing was unmatched instead of only how many.
    std::vector<std::string> refVerdict(ref.size(), "matched");
    json extras = json::array();
    bool sameShape = cand.size() == ref.size();
    for (size_t i = 0; sameShape && i < cand.size(); ++i) {
      sameShape = cand[i].kind == ref[i].kind;
    }
    bool identical = true;
    if (sameShape) {
      for (size_t i = 0; i < cand.size(); ++i) {
        if (std::fabs(cand[i].t - ref[i].t) > cfg.matchTolerance) {
          ++nDisp;
          identical = false;
          refVerdict[i] = "displaced";
        }
      }
    } else {
      identical = false;
      size_t i = 0, j = 0;
      while (i < cand.size() && j < ref.size()) {
        const double delta = cand[i].t - ref[j].t;
        if (std::fabs(delta) <= cfg.matchTolerance) {
          if (cand[i].kind != ref[j].kind) {
            ++nKind;
            refVerdict[j] = "kind";
          }
          ++i;
          ++j;
        } else if (delta < 0.) {
          ++nExtra;
          extras.push_back(cand[i].t);
          ++i;
        } else {
          ++nLost;
          refVerdict[j] = "lost";
          ++j;
        }
      }
      for (; i < cand.size(); ++i) {
        ++nExtra;
        extras.push_back(cand[i].t);
      }
      for (; j < ref.size(); ++j) {
        ++nLost;
        refVerdict[j] = "lost";
      }
    }
    nIdent += identical;
    if (identical && !alsoClean) {
      continue;
    }

    // The kernel's own ray/surface hit list along this very ray, from the ray origin. This is the
    // multiset DistFromOutside/DistFromInside select from, so a reference crossing that appears
    // here and is still missing from the stepped list was FILTERED, not missed.
    std::vector<O2BVHSurfaceSolid::ContainsCrossing> bvh, loop;
    solid->DescribeContainsCrossings(origin, dir, bvh, loop);

    json refs = json::array();
    for (size_t j = 0; j < ref.size(); ++j) {
      json e;
      e["t"] = ref[j].t;
      e["kind"] = ref[j].kind;
      e["verdict"] = refVerdict[j];
      double p[3];
      for (int k = 0; k < 3; ++k) {
        p[k] = origin[k] + ref[j].t * dir[k];
      }
      e["point"] = {p[0], p[1], p[2]};
      const bool in = solid->Contains(p);
      e["safety"] = solid->Safety(p, in ? kTRUE : kFALSE);
      // The kernel hit nearest this reference distance, if any.
      int best = -1;
      double bestD = 1.e300;
      for (size_t h = 0; h < bvh.size(); ++h) {
        const double d = std::fabs(bvh[h].distance - ref[j].t);
        if (d < bestD) {
          bestD = d;
          best = static_cast<int>(h);
        }
      }
      if (best >= 0 && bestD <= cfg.matchTolerance) {
        e["hit"] = true;
        e["hitDelta"] = bestD;
        e["normalDotDir"] = bvh[best].normalAlignment;
        e["onTrimBoundary"] = bvh[best].onTrimBoundary;
        // dot(n, d) < 0 means the ray goes against the outward normal, i.e. ENTERS. The reference
        // kind is +1 for enter. These must agree; when they do not, the face's outward normal
        // points the wrong way.
        const int impliedKind = bvh[best].normalAlignment < 0. ? +1 : -1;
        e["impliedKind"] = impliedKind;
        e["normalFlipped"] = (impliedKind != ref[j].kind);
      } else {
        e["hit"] = false;
        e["hitDelta"] = bestD;
      }
      refs.push_back(std::move(e));
    }

    json rj;
    rj["index"] = ri;
    rj["o"] = {origin[0], origin[1], origin[2]};
    rj["d"] = {dir[0], dir[1], dir[2]};
    rj["tmax"] = tMax;
    rj["identical"] = identical;
    rj["unterminated"] = stats.unterminated > 0;
    json ct = json::array(), ck = json::array();
    for (const auto& c : cand) {
      ct.push_back(c.t);
      ck.push_back(c.kind);
    }
    rj["candT"] = ct;
    rj["candK"] = ck;
    rj["extraT"] = extras;
    rj["ref"] = refs;
    rj["nKernelHits"] = bvh.size();
    rj["nKernelHitsLoop"] = loop.size();
    rays.push_back(std::move(rj));
  }

  out["rays"] = rays;
  out["summary"] = {{"rays", nRays}, {"identical", nIdent}, {"lost", nLost}, {"extra", nExtra},
                    {"kind", nKind}, {"displaced", nDisp}, {"unterminated", nUnterm}};
  std::ofstream o(outFile);
  o << out.dump();
  std::printf("%s: rays=%lld identical=%lld LOST=%lld extra=%lld kind=%lld displaced=%lld unterm=%lld"
              " (matchTol=%.3g, %d surfaces, %s)\n",
              outFile.c_str(), nRays, nIdent, nLost, nExtra, nKind, nDisp, nUnterm,
              cfg.matchTolerance, solid->GetNsurfaces(),
              O2BVHSurfaceSolid::GetNavigationReliabilityName(solid->GetNavigationReliability()));
  return 0;
}
