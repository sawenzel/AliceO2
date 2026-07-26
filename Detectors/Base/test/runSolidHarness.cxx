// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file runSolidHarness.cxx
/// \brief Front-end for the O2SolidHarness validation / performance comparison harness.
///
/// Built as o2-bench-detectorsbase-solid-harness (see Detectors/Base/CMakeLists.txt). Loads
/// paired surfaces_*.bin / facets_*.bin parts from a test-part database (see
/// scripts/geometry/makeTestPartDB.py) and, for each, validates and times
/// O2BVHSurfaceSolid (candidate) against O2Tessellated (reference). Design, ground rules and the
/// `--only` single-kernel `perf record` use case are documented in
/// scripts/geometry/SolidNavigationHarness.md, Step 3.

#include "DetectorsBase/O2SolidHarness.h"
#include "DetectorsBase/O2BVHSurfaceSolid.h"
#include "DetectorsBase/O2Tessellated.h"
#include "DetectorsBase/O2SurfaceSolidIO.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using json = nlohmann::json;
using namespace o2::base;
using namespace o2::base::harness;

namespace
{

struct Options {
  std::string db;
  std::string explicitSurfaces;
  std::string explicitFacets;
  std::string partsPattern;
  int points = 5000;
  int rays = 5000;
  uint64_t seed = 1;
  std::set<std::string> only = {"contains", "distout", "distin", "safety"};
  bool loopCrosscheck = false;
  bool pruningAb = false;
  std::string jsonOut;
  int warmup = 1;
  int repeat = 3;
};

struct Part {
  std::string id;
  std::string model;
  std::string surfaces;
  std::string facets;
};

void printUsage(const char* argv0)
{
  std::cout <<
    "Usage: " << argv0 << " --db <dir> [--parts <substring>] [--points N] [--rays N] [--seed N]\n"
    "                              [--only contains,distout,distin,safety] [--loop-crosscheck]\n"
    "                              [--pruning-ab] [--json <out.json>] [--warmup N] [--repeat N]\n"
    "   or: " << argv0 << " --surfaces <file> --facets <file> [options as above]\n\n"
    "  --loop-crosscheck  also run the surface solid's non-BVH _Loop twins and require exact\n"
    "                     agreement; this is the correctness guard that does not involve the mesh\n"
    "  --pruning-ab       re-run the distance kernels with ray tmax pruning disabled, reporting\n"
    "                     the BVH candidate counts and ns/call both ways (prices the optimization)\n\n"
    "perf record entry point (single kernel, one part):\n"
    "  perf record -g " << argv0 << " --db <db> --parts oTOF --only distout --rays 200000\n";
}

std::set<std::string> splitCsv(const std::string& s)
{
  std::set<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (!tok.empty()) {
      out.insert(tok);
    }
  }
  return out;
}

bool parseArgs(int argc, char** argv, Options& opt)
{
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](const char* flag) -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + flag);
      }
      return argv[++i];
    };
    if (a == "--db") {
      opt.db = next("--db");
    } else if (a == "--surfaces") {
      opt.explicitSurfaces = next("--surfaces");
    } else if (a == "--facets") {
      opt.explicitFacets = next("--facets");
    } else if (a == "--parts") {
      opt.partsPattern = next("--parts");
    } else if (a == "--points") {
      opt.points = std::stoi(next("--points"));
    } else if (a == "--rays") {
      opt.rays = std::stoi(next("--rays"));
    } else if (a == "--seed") {
      opt.seed = std::stoull(next("--seed"));
    } else if (a == "--only") {
      opt.only = splitCsv(next("--only"));
    } else if (a == "--loop-crosscheck") {
      opt.loopCrosscheck = true;
    } else if (a == "--pruning-ab") {
      opt.pruningAb = true;
    } else if (a == "--json") {
      opt.jsonOut = next("--json");
    } else if (a == "--warmup") {
      opt.warmup = std::stoi(next("--warmup"));
    } else if (a == "--repeat") {
      opt.repeat = std::stoi(next("--repeat"));
    } else if (a == "-h" || a == "--help") {
      printUsage(argv[0]);
      return false;
    } else {
      throw std::runtime_error("unrecognized option: " + a);
    }
  }
  if (opt.db.empty() && (opt.explicitSurfaces.empty() || opt.explicitFacets.empty())) {
    throw std::runtime_error("either --db <dir> or both --surfaces/--facets are required");
  }
  return true;
}

std::vector<Part> collectParts(const Options& opt)
{
  std::vector<Part> parts;
  if (!opt.explicitSurfaces.empty()) {
    parts.push_back({"adhoc", "adhoc", opt.explicitSurfaces, opt.explicitFacets});
    return parts;
  }
  const std::string manifestPath = opt.db + "/manifest.json";
  std::ifstream in(manifestPath);
  if (!in) {
    throw std::runtime_error("cannot open " + manifestPath);
  }
  json manifest;
  in >> manifest;
  for (const auto& p : manifest.at("parts")) {
    Part part;
    part.id = p.at("id").get<std::string>();
    part.model = p.at("model").get<std::string>();
    part.surfaces = p.at("surfaces").get<std::string>();
    part.facets = p.at("facets").get<std::string>();
    if (!opt.partsPattern.empty()) {
      const bool idMatch = part.id.find(opt.partsPattern) != std::string::npos;
      const bool modelMatch = part.model.find(opt.partsPattern) != std::string::npos;
      if (!idMatch && !modelMatch) {
        continue;
      }
    }
    parts.push_back(std::move(part));
  }
  return parts;
}

json validationToJson(const ValidationResult& r)
{
  json j;
  j["nSamples"] = r.nSamples;
  j["nAgree"] = r.nAgree;
  j["nMismatchWithinBand"] = r.nMismatchWithinBand;
  j["nMismatchUnexplained"] = r.nMismatchUnexplained;
  j["worstDeviation"] = r.worstDeviation;
  json offenders = json::array();
  for (const auto& o : r.worstOffenders) {
    offenders.push_back({{"point", {o.point[0], o.point[1], o.point[2]}},
                         {"dir", {o.dir[0], o.dir[1], o.dir[2]}},
                         {"candidateValue", o.candidateValue},
                         {"referenceValue", o.referenceValue},
                         {"deviation", o.deviation},
                         {"referenceSafety", o.referenceSafety}});
  }
  j["worstOffenders"] = offenders;
  return j;
}

json timingToJson(const TimingResult& t)
{
  return {{"nCalls", t.nCalls}, {"nsPerCall", t.nsPerCall}, {"checksum", t.checksum}};
}

void printValidation(const std::string& name, const ValidationResult& r)
{
  const double agreePct = r.nSamples ? 100. * static_cast<double>(r.nAgree) / static_cast<double>(r.nSamples) : 0.;
  std::printf("  %-10s samples=%-7zu agree=%6.2f%%  mismatch(band=%zu,unexplained=%zu)  worstDev=%.6g\n",
             name.c_str(), r.nSamples, agreePct, r.nMismatchWithinBand, r.nMismatchUnexplained, r.worstDeviation);
  if (r.nMismatchUnexplained > 0) {
    const size_t nShow = std::min<size_t>(3, r.worstOffenders.size());
    for (size_t i = 0; i < nShow; ++i) {
      const auto& o = r.worstOffenders[i];
      std::printf("      offender[%zu]: point=(%.6g,%.6g,%.6g) dir=(%.6g,%.6g,%.6g) cand=%.6g ref=%.6g dev=%.6g refSafety=%.6g\n",
                 i, o.point[0], o.point[1], o.point[2], o.dir[0], o.dir[1], o.dir[2], o.candidateValue,
                 o.referenceValue, o.deviation, o.referenceSafety);
    }
  }
}

void printTiming(const std::string& name, const TimingResult& candidate, const TimingResult& reference)
{
  const double ratio = reference.nsPerCall > 0. ? candidate.nsPerCall / reference.nsPerCall : 0.;
  std::printf("  %-10s candidate=%9.1f ns/call   reference=%9.1f ns/call   ratio(cand/ref)=%.2fx\n", name.c_str(),
             candidate.nsPerCall, reference.nsPerCall, ratio);
}

// What the BVH traversal buys over the all-surfaces loop on the *same* shape: unlike the
// candidate/reference ratio this compares like with like, so it prices the acceleration structure
// alone rather than analytic patches against triangles.
void printLoopSpeedup(const std::string& name, const TimingResult& bvh, const TimingResult& loop)
{
  const double speedup = bvh.nsPerCall > 0. ? loop.nsPerCall / bvh.nsPerCall : 0.;
  std::printf("  %-10s   BVH=%9.1f ns/call   _Loop=%9.1f ns/call   speedup(loop/bvh)=%.2fx\n", name.c_str(),
             bvh.nsPerCall, loop.nsPerCall, speedup);
}

double toSeconds(std::chrono::steady_clock::time_point t0, std::chrono::steady_clock::time_point t1)
{
  return std::chrono::duration<double>(t1 - t0).count();
}

} // namespace

int main(int argc, char** argv)
{
  Options opt;
  try {
    if (!parseArgs(argc, argv, opt)) {
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    printUsage(argv[0]);
    return 1;
  }

  std::vector<Part> parts;
  try {
    parts = collectParts(opt);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  if (parts.empty()) {
    std::cerr << "no parts matched (pattern='" << opt.partsPattern << "')\n";
    return 1;
  }

  json jsonReport = json::array();

  for (const auto& part : parts) {
    std::printf("=== %s (%s) ===\n", part.id.c_str(), part.model.c_str());

    O2BVHSurfaceSolid surf(part.id.c_str());
    if (!LoadSurfaceSolid(part.surfaces, surf)) {
      std::cerr << "  skip: LoadSurfaceSolid failed for " << part.surfaces << "\n";
      continue;
    }
    auto t0 = std::chrono::steady_clock::now();
    surf.CloseShape(true);
    auto t1 = std::chrono::steady_clock::now();
    const double surfCloseSeconds = toSeconds(t0, t1);

    O2Tessellated mesh(part.id.c_str());
    if (!LoadFacetSolid(part.facets, mesh)) {
      std::cerr << "  skip: LoadFacetSolid failed for " << part.facets << "\n";
      continue;
    }
    t0 = std::chrono::steady_clock::now();
    mesh.CloseShape();
    t1 = std::chrono::steady_clock::now();
    const double meshCloseSeconds = toSeconds(t0, t1);

    const TGeoShape* candidate = &surf;
    const TGeoShape* reference = &mesh;

    const Point3D bboxMin{mesh.GetOrigin()[0] - mesh.GetDX(), mesh.GetOrigin()[1] - mesh.GetDY(),
                          mesh.GetOrigin()[2] - mesh.GetDZ()};
    const Point3D bboxMax{mesh.GetOrigin()[0] + mesh.GetDX(), mesh.GetOrigin()[1] + mesh.GetDY(),
                          mesh.GetOrigin()[2] + mesh.GetDZ()};

    std::printf("  surfaces=%d triangles=%d  closeShape: surface=%.4fs mesh=%.4fs\n", surf.GetNsurfaces(),
               mesh.GetNfacets(), surfCloseSeconds, meshCloseSeconds);

    SampleConfig cfg;
    cfg.nBulk = opt.points;
    cfg.nBoundary = opt.points;
    cfg.nInside = std::max(1, opt.points / 2);
    cfg.nOutsideRays = opt.rays;
    cfg.nInsideRays = std::max(1, opt.rays / 2);
    cfg.seed = opt.seed;
    const SampleSet samples = generateSamples(reference, bboxMin, bboxMax, cfg);

    long long candidatesSampled = 0;
    const size_t nProbe = std::min<size_t>(200, samples.outsideRays.size());
    for (size_t i = 0; i < nProbe; ++i) {
      const auto& r = samples.outsideRays[i];
      const int n = surf.CountBVHRayCandidates(r.origin, r.dir);
      if (n > 0) {
        candidatesSampled += n;
      }
    }
    std::printf("  BVH ray candidates: sum=%lld over %zu probe rays\n", candidatesSampled, nProbe);

    json partJson;
    partJson["id"] = part.id;
    partJson["model"] = part.model;
    partJson["nSurfaces"] = surf.GetNsurfaces();
    partJson["nTriangles"] = mesh.GetNfacets();
    partJson["closeShapeSecondsSurface"] = surfCloseSeconds;
    partJson["closeShapeSecondsMesh"] = meshCloseSeconds;
    partJson["bvhRayCandidatesSampled"] = candidatesSampled;
    partJson["bvhRayCandidatesProbeRays"] = nProbe;

    std::vector<Point3D> allPoints = samples.bulkPoints;
    allPoints.insert(allPoints.end(), samples.boundaryPoints.begin(), samples.boundaryPoints.end());
    allPoints.insert(allPoints.end(), samples.insidePoints.begin(), samples.insidePoints.end());

    if (opt.only.count("contains")) {
      auto v = validateContains(candidate, reference, allPoints);
      printValidation("contains", v);
      auto tc = timeContains(candidate, allPoints, opt.warmup, opt.repeat);
      auto tr = timeContains(reference, allPoints, opt.warmup, opt.repeat);
      printTiming("contains", tc, tr);
      partJson["contains"] = {{"validation", validationToJson(v)},
                              {"timingCandidate", timingToJson(tc)},
                              {"timingReference", timingToJson(tr)}};
    }
    if (opt.only.count("distout")) {
      auto v = validateDistFromOutside(candidate, reference, samples.outsideRays);
      printValidation("distout", v);
      auto tc = timeDistFromOutside(candidate, samples.outsideRays, opt.warmup, opt.repeat);
      auto tr = timeDistFromOutside(reference, samples.outsideRays, opt.warmup, opt.repeat);
      printTiming("distout", tc, tr);
      // the all-surfaces baseline: what the BVH traversal buys over visiting every patch
      auto tl = timeRayKernel(samples.outsideRays, opt.warmup, opt.repeat,
                              [&surf](const Point3D& o, const Point3D& d) {
                                return surf.DistFromOutside_Loop(o.data(), d.data());
                              });
      printLoopSpeedup("distout", tc, tl);
      partJson["distout"] = {{"validation", validationToJson(v)},
                             {"timingCandidate", timingToJson(tc)},
                             {"timingReference", timingToJson(tr)},
                             {"timingCandidateLoop", timingToJson(tl)}};
    }
    if (opt.only.count("distin")) {
      auto v = validateDistFromInside(candidate, reference, samples.insideRays);
      printValidation("distin", v);
      auto tc = timeDistFromInside(candidate, samples.insideRays, opt.warmup, opt.repeat);
      auto tr = timeDistFromInside(reference, samples.insideRays, opt.warmup, opt.repeat);
      printTiming("distin", tc, tr);
      auto tl = timeRayKernel(samples.insideRays, opt.warmup, opt.repeat,
                              [&surf](const Point3D& o, const Point3D& d) {
                                return surf.DistFromInside_Loop(o.data(), d.data());
                              });
      printLoopSpeedup("distin", tc, tl);
      partJson["distin"] = {{"validation", validationToJson(v)},
                            {"timingCandidate", timingToJson(tc)},
                            {"timingReference", timingToJson(tr)},
                            {"timingCandidateLoop", timingToJson(tl)}};
    }
    if (opt.only.count("safety")) {
      // Never compared against each other (see ground rules): each shape's Safety() is checked
      // against its own DistFrom{Inside,Outside} contract independently.
      auto vc = validateSafety(candidate, allPoints);
      auto vr = validateSafety(reference, allPoints);
      printValidation("safety(cand)", vc);
      printValidation("safety(ref)", vr);
      auto tc = timeSafety(candidate, allPoints, opt.warmup, opt.repeat);
      auto tr = timeSafety(reference, allPoints, opt.warmup, opt.repeat);
      printTiming("safety", tc, tr);
      partJson["safety"] = {{"validationCandidate", validationToJson(vc)},
                            {"validationReference", validationToJson(vr)},
                            {"timingCandidate", timingToJson(tc)},
                            {"timingReference", timingToJson(tr)}};
    }

    if (opt.loopCrosscheck) {
      // Independent of the tessellated reference entirely: separates BVH/traversal bugs from
      // surface-kernel bugs (see SolidNavigationHarness.md, "Why these three are one piece of
      // work"). The distance twins must agree *exactly*, not within a tolerance: both take a
      // minimum over the same hits from the same kernels and differ only in which surfaces the
      // BVH lets them skip, so any difference at all is a traversal or pruning bug.
      size_t containsAgree = 0;
      for (const auto& p : allPoints) {
        if (surf.Contains(p.data()) == surf.Contains_Loop(p.data())) {
          ++containsAgree;
        }
      }
      std::printf("  loop-crosscheck contains: BVH==Loop for %zu/%zu points\n", containsAgree, allPoints.size());
      partJson["loopCrosscheckContains"] = {{"agree", containsAgree}, {"total", allPoints.size()}};

      size_t outAgree = 0;
      double worstOutDeviation = 0.;
      for (const auto& r : samples.outsideRays) {
        const double bvh = surf.DistFromOutside(r.origin.data(), r.dir.data(), 3);
        const double loop = surf.DistFromOutside_Loop(r.origin.data(), r.dir.data());
        if (bvh == loop) {
          ++outAgree;
        } else {
          worstOutDeviation = std::max(worstOutDeviation, std::fabs(bvh - loop));
        }
      }
      std::printf("  loop-crosscheck distout : BVH==Loop for %zu/%zu rays (worstDev=%.6g)\n", outAgree,
                 samples.outsideRays.size(), worstOutDeviation);
      partJson["loopCrosscheckDistOutside"] = {
        {"agree", outAgree}, {"total", samples.outsideRays.size()}, {"worstDeviation", worstOutDeviation}};

      size_t inAgree = 0;
      double worstInDeviation = 0.;
      for (const auto& r : samples.insideRays) {
        const double bvh = surf.DistFromInside(r.origin.data(), r.dir.data(), 3);
        const double loop = surf.DistFromInside_Loop(r.origin.data(), r.dir.data());
        if (bvh == loop) {
          ++inAgree;
        } else {
          worstInDeviation = std::max(worstInDeviation, std::fabs(bvh - loop));
        }
      }
      std::printf("  loop-crosscheck distin  : BVH==Loop for %zu/%zu rays (worstDev=%.6g)\n", inAgree,
                 samples.insideRays.size(), worstInDeviation);
      partJson["loopCrosscheckDistInside"] = {
        {"agree", inAgree}, {"total", samples.insideRays.size()}, {"worstDeviation", worstInDeviation}};
    }

    if (opt.pruningAb) {
      // Prices the ray tmax tightening: the same rays run with it on and off, reporting both the
      // surface patches the traversal actually handed to the leaf callback and the wall time. The
      // answers must be bit-identical -- the switch is a cost knob, never a semantic one, and a
      // mismatch here is a bug in the tightening rather than a measurement.
      json pruningJson;
      size_t identical = 0;
      std::vector<double> prunedValues;
      prunedValues.reserve(samples.outsideRays.size());

      O2BVHSurfaceSolid::SetRayTMaxPruning(true);
      O2BVHSurfaceSolid::ResetRayCandidateCounter();
      for (const auto& r : samples.outsideRays) {
        prunedValues.push_back(surf.DistFromOutside(r.origin.data(), r.dir.data(), 3));
      }
      const long long prunedCandidates = O2BVHSurfaceSolid::GetRayCandidateCount();
      auto tPruned = timeDistFromOutside(candidate, samples.outsideRays, opt.warmup, opt.repeat);

      O2BVHSurfaceSolid::SetRayTMaxPruning(false);
      O2BVHSurfaceSolid::ResetRayCandidateCounter();
      for (size_t i = 0; i < samples.outsideRays.size(); ++i) {
        const auto& r = samples.outsideRays[i];
        if (surf.DistFromOutside(r.origin.data(), r.dir.data(), 3) == prunedValues[i]) {
          ++identical;
        }
      }
      const long long unprunedCandidates = O2BVHSurfaceSolid::GetRayCandidateCount();
      auto tUnpruned = timeDistFromOutside(candidate, samples.outsideRays, opt.warmup, opt.repeat);
      O2BVHSurfaceSolid::SetRayTMaxPruning(true);

      const double candidateRatio =
        unprunedCandidates > 0 ? static_cast<double>(prunedCandidates) / static_cast<double>(unprunedCandidates) : 0.;
      const double speedup = tPruned.nsPerCall > 0. ? tUnpruned.nsPerCall / tPruned.nsPerCall : 0.;
      std::printf("  tmax-pruning A/B (distout, %zu rays): identical=%zu/%zu\n", samples.outsideRays.size(), identical,
                 samples.outsideRays.size());
      std::printf("      candidates: pruned=%lld  unpruned=%lld  (%.1f%% of the work)\n", prunedCandidates,
                 unprunedCandidates, 100. * candidateRatio);
      std::printf("      time      : pruned=%9.1f ns/call  unpruned=%9.1f ns/call  speedup=%.2fx\n",
                 tPruned.nsPerCall, tUnpruned.nsPerCall, speedup);

      pruningJson["identical"] = identical;
      pruningJson["total"] = samples.outsideRays.size();
      pruningJson["candidatesPruned"] = prunedCandidates;
      pruningJson["candidatesUnpruned"] = unprunedCandidates;
      pruningJson["timingPruned"] = timingToJson(tPruned);
      pruningJson["timingUnpruned"] = timingToJson(tUnpruned);
      partJson["tmaxPruningAB"] = std::move(pruningJson);
    }

    jsonReport.push_back(std::move(partJson));
  }

  if (!opt.jsonOut.empty()) {
    std::ofstream out(opt.jsonOut);
    out << jsonReport.dump(1);
    std::printf("\nWrote %s\n", opt.jsonOut.c_str());
  }

  return 0;
}
