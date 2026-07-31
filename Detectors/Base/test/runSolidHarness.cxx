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
#include <array>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <map>
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
  std::string dumpSamples; ///< directory to write per-part sample sets into, for the OCCT oracle
  std::string refAnswers;  ///< directory holding the oracle's answers for those sample sets
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
    "                     the BVH candidate counts and ns/call both ways (prices the optimization)\n"
    "  --dump-samples D   write each part's sample set to D/samples_<part>.json\n"
    "  --ref-answers D    validate against D/answers_<part>.json instead of the mesh; those are\n"
    "                     produced by scripts/geometry/occtOracle.py from the part's .brep, so a\n"
    "                     disagreement outside the model tolerance is a defect, not chording\n\n"
    "OCCT oracle round trip:\n"
    "  " << argv0 << " --db <db> --dump-samples /tmp/o\n"
    "  occtOracle.py --brep <part>.brep --samples /tmp/o/samples_<part>.json \\\n"
    "                --out /tmp/o/answers_<part>.json\n"
    "  " << argv0 << " --db <db> --ref-answers /tmp/o\n\n"
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
    } else if (a == "--dump-samples") {
      opt.dumpSamples = next("--dump-samples");
    } else if (a == "--ref-answers") {
      opt.refAnswers = next("--ref-answers");
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
  j["nMismatchMissedSurface"] = r.nMismatchMissedSurface;
  j["nMismatchUnexplained"] = r.nMismatchUnexplained;
  j["nNoVerdict"] = r.nNoVerdict;
  j["nRelabelled"] = r.nRelabelled;
  j["worstDeviation"] = r.worstDeviation;
  json offenders = json::array();
  for (const auto& o : r.worstOffenders) {
    offenders.push_back({{"point", {o.point[0], o.point[1], o.point[2]}},
                         {"dir", {o.dir[0], o.dir[1], o.dir[2]}},
                         {"candidateValue", o.candidateValue},
                         {"referenceValue", o.referenceValue},
                         {"deviation", o.deviation},
                         {"referenceSafety", o.referenceSafety},
                         {"incidenceCosine", o.incidenceCosine}});
  }
  j["worstOffenders"] = offenders;
  return j;
}

// The sample/answer JSON contract shared with scripts/geometry/occtOracle.py. Bump on both sides
// together; the oracle refuses a version it does not speak rather than guessing.
constexpr int kOracleFormatVersion = 1;

/// Part ids carry '/' and other path-hostile characters; the oracle round trip pairs files by
/// this sanitized form on both sides.
std::string sanitizePartId(const std::string& id)
{
  std::string out;
  out.reserve(id.size());
  for (const char c : id) {
    out.push_back((std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '.') ? c : '_');
  }
  return out;
}

json pointsToJson(const std::vector<Point3D>& points)
{
  json array = json::array();
  for (const auto& p : points) {
    array.push_back({p[0], p[1], p[2]});
  }
  return array;
}

json raysToJson(const std::vector<Ray>& rays)
{
  json array = json::array();
  for (const auto& r : rays) {
    array.push_back({{"o", {r.origin[0], r.origin[1], r.origin[2]}},
                     {"d", {r.dir[0], r.dir[1], r.dir[2]}}});
  }
  return array;
}

/// Serialize a sample set so an external oracle can answer exactly the same queries. The samples
/// come from a seeded mt19937_64 inside the harness, so nothing outside can regenerate them;
/// dumping is the only way to ask another implementation about the same points.
void writeSamples(const std::string& dir, const std::string& partId, const SampleSet& samples)
{
  json doc;
  doc["version"] = kOracleFormatVersion;
  doc["part"] = partId;
  doc["bboxMin"] = {samples.bboxMin[0], samples.bboxMin[1], samples.bboxMin[2]};
  doc["bboxMax"] = {samples.bboxMax[0], samples.bboxMax[1], samples.bboxMax[2]};
  doc["points"] = {{"bulk", pointsToJson(samples.bulkPoints)},
                   {"boundary", pointsToJson(samples.boundaryPoints)},
                   {"inside", pointsToJson(samples.insidePoints)}};
  doc["rays"] = {{"outside", raysToJson(samples.outsideRays)},
                 {"inside", raysToJson(samples.insideRays)}};
  const std::string path = dir + "/samples_" + sanitizePartId(partId) + ".json";
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  out << doc.dump(1);
  std::printf("  wrote samples: %s\n", path.c_str());
}

/// Oracle answers for one part, or `has == false` when no answer file exists for it.
struct OracleAnswers {
  bool has = false;
  double tolerance = 0.;
  double capacity = 0.;
  bool valid = false;
  std::map<std::string, std::vector<int>> containsState;
  /// Per ray category, the oracle's classification of each ray *origin* (1/0/-1). This is what
  /// makes the distance columns soundly categorised rather than categorised by the reference mesh.
  std::map<std::string, std::vector<int>> originContains;
  std::map<std::string, std::vector<double>> boundaryDistance;
  std::map<std::string, std::vector<double>> distOutside;
  std::map<std::string, std::vector<double>> distInside;
};

template <typename T>
std::map<std::string, std::vector<T>> readColumns(const json& parent, const char* key)
{
  std::map<std::string, std::vector<T>> columns;
  if (!parent.contains(key)) {
    return columns;
  }
  for (const auto& [category, values] : parent.at(key).items()) {
    columns[category] = values.template get<std::vector<T>>();
  }
  return columns;
}

OracleAnswers loadOracleAnswers(const std::string& dir, const std::string& partId)
{
  OracleAnswers answers;
  const std::string path = dir + "/answers_" + sanitizePartId(partId) + ".json";
  std::ifstream in(path);
  if (!in) {
    std::printf("  oracle: no answers file %s, skipping oracle validation\n", path.c_str());
    return answers;
  }
  json doc;
  in >> doc;
  const int version = doc.value("version", -1);
  if (version != kOracleFormatVersion) {
    throw std::runtime_error(path + ": answer format version " + std::to_string(version) +
                             ", this harness speaks " + std::to_string(kOracleFormatVersion));
  }
  answers.has = true;
  answers.tolerance = doc.value("tolerance", 0.);
  answers.capacity = doc.value("capacity", 0.);
  answers.valid = doc.value("valid", false);
  answers.containsState = readColumns<int>(doc, "contains");
  answers.originContains = readColumns<int>(doc, "originContains");
  answers.boundaryDistance = readColumns<double>(doc, "safetyUpperBound");
  answers.distOutside = readColumns<double>(doc, "distFromOutside");
  answers.distInside = readColumns<double>(doc, "distFromInside");
  return answers;
}

/// Concatenate the oracle's per-category columns in the same order the harness concatenates its
/// point categories, so index i of the merged column belongs to point i of `allPoints`.
///
/// Each column is padded to its category's *point count* before the next one is appended. That
/// padding is not cosmetic: the oracle answers `contains` for every point but caps the expensive
/// boundary-distance query, so its columns have different lengths per category. Concatenating
/// them raw would shift every later category's answers onto the wrong points -- a silent,
/// systematic mis-scoring rather than an error.
template <typename T>
std::vector<T> mergeCategories(const std::map<std::string, std::vector<T>>& columns,
                               const std::array<size_t, 3>& categorySizes, T missing)
{
  static constexpr std::array<const char*, 3> kOrder = {"bulk", "boundary", "inside"};
  std::vector<T> merged;
  for (size_t categoryIndex = 0; categoryIndex < kOrder.size(); ++categoryIndex) {
    const size_t expected = categorySizes[categoryIndex];
    const auto it = columns.find(kOrder[categoryIndex]);
    const size_t available = it == columns.end() ? 0 : std::min(expected, it->second.size());
    for (size_t i = 0; i < available; ++i) {
      merged.push_back(it->second[i]);
    }
    merged.insert(merged.end(), expected - available, missing);
  }
  return merged;
}

json timingToJson(const TimingResult& t)
{
  return {{"nCalls", t.nCalls}, {"nsPerCall", t.nsPerCall}, {"checksum", t.checksum}};
}

void printValidation(const std::string& name, const ValidationResult& r)
{
  // Scored = everything the reference was willing to answer. Reporting the percentage against
  // nSamples would let a reference that abstains on half the points look like agreement.
  const size_t scored = r.nSamples - r.nNoVerdict;
  const double agreePct = scored ? 100. * static_cast<double>(r.nAgree) / static_cast<double>(scored) : 0.;
  std::printf("  %-10s scored=%-7zu agree=%6.2f%%  mismatch(band=%zu,missed=%zu,unexplained=%zu)"
             "  noVerdict=%zu  worstDev=%.6g\n",
             name.c_str(), scored, agreePct, r.nMismatchWithinBand, r.nMismatchMissedSurface,
             r.nMismatchUnexplained, r.nNoVerdict, r.worstDeviation);
  if (r.nRelabelled > 0) {
    // Not a candidate result: it says how many rays the sample generator had put in the wrong
    // category, which is a statement about the reference mesh. Printed so an improvement in these
    // columns is never mistaken for a kernel improvement.
    std::printf("  %-10s relabelled=%zu ray(s) by the oracle's own origin classification\n",
               name.c_str(), r.nRelabelled);
  }
  if (r.nMismatchUnexplained > 0 || r.nMismatchMissedSurface > 0) {
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
  std::vector<std::string> unreliableParts;

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

    // Label every measurement with whether its subject is a closed manifold at all. Without this
    // line an "unexplained" mismatch count below is unattributable: on an open surface set the
    // exact solid's own answers are undefined in the shadow of each gap, so the number says
    // nothing about mesh chording (see scripts/geometry/ExactTrimTopology.md, item 4).
    const auto reliability = surf.GetNavigationReliability();
    const char* reliabilityName = O2BVHSurfaceSolid::GetNavigationReliabilityName(reliability);
    const bool navigable = surf.IsNavigable();
    std::printf("  navigation: %s%s  (boundary=%d non-manifold=%d reversed=%d)\n", reliabilityName,
               navigable ? "" : "  *** UNRELIABLE: results below are not a measurement of accuracy ***",
               surf.GetBoundaryEdgeCount(), surf.GetNonManifoldEdgeCount(), surf.GetReversedEdgeCount());
    // The same boundary measured as curves, in cm. Nothing above derives from it yet; it is here
    // to answer "how far apart are the faces", which a chord count cannot. Always with the chord
    // resolution next to it -- a gap below that is how the rims were sampled, not a real gap.
    std::printf("  rim gap: max %.3g cm (chord resolution %.3g cm, matched at %.3g cm); rims %d "
                "(matched=%d boundary=%d non-manifold=%d reversed=%d), open %.3g of %.3g cm\n",
                surf.GetMaxRimGap(), surf.GetRimChordResolution(), surf.GetRimMatchTolerance(), surf.GetRimCount(),
                surf.GetMatchedRimCount(), surf.GetBoundaryRimCount(), surf.GetNonManifoldRimCount(),
                surf.GetReversedRimCount(), surf.GetUnmatchedRimLength(), surf.GetTotalRimLength());
    if (!navigable) {
      unreliableParts.push_back(part.id + " (" + reliabilityName + ")");
    }

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
    partJson["navigation"] = {{"reliability", reliabilityName},
                              {"navigable", navigable},
                              {"boundaryEdges", surf.GetBoundaryEdgeCount()},
                              {"nonManifoldEdges", surf.GetNonManifoldEdgeCount()},
                              {"reversedEdges", surf.GetReversedEdgeCount()},
                              {"maxRimGap", surf.GetMaxRimGap()},
                              {"rimChordResolution", surf.GetRimChordResolution()},
                              {"rimMatchTolerance", surf.GetRimMatchTolerance()},
                              {"totalRimLength", surf.GetTotalRimLength()},
                              {"unmatchedRimLength", surf.GetUnmatchedRimLength()},
                              {"rims", surf.GetRimCount()},
                              {"matchedRims", surf.GetMatchedRimCount()},
                              {"boundaryRims", surf.GetBoundaryRimCount()},
                              {"nonManifoldRims", surf.GetNonManifoldRimCount()},
                              {"reversedRims", surf.GetReversedRimCount()}};

    std::vector<Point3D> allPoints = samples.bulkPoints;
    allPoints.insert(allPoints.end(), samples.boundaryPoints.begin(), samples.boundaryPoints.end());
    allPoints.insert(allPoints.end(), samples.insidePoints.begin(), samples.insidePoints.end());
    const std::array<size_t, 3> categorySizes{samples.bulkPoints.size(), samples.boundaryPoints.size(),
                                              samples.insidePoints.size()};

    if (!opt.dumpSamples.empty()) {
      writeSamples(opt.dumpSamples, part.id, samples);
    }

    // Ground-truth validation, when the oracle has answered this part. Kept separate from the
    // mesh comparison below rather than replacing it: the mesh columns stay comparable with every
    // measurement recorded so far, while these columns are the ones a gate can be written against.
    if (!opt.refAnswers.empty()) {
      const OracleAnswers oracle = loadOracleAnswers(opt.refAnswers, part.id);
      if (oracle.has) {
        ValidationOptions oracleOpt;
        // The band is now the model's own declared tolerance instead of a guessed mesh sagitta.
        // A floor keeps a perfectly-toleranced synthetic fixture from demanding bit equality.
        oracleOpt.meshBand = std::max(oracle.tolerance, oracleOpt.distanceTolerance);
        const auto boundaryDistance =
          mergeCategories<double>(oracle.boundaryDistance, categorySizes, -1.);
        const auto containsState = mergeCategories<int>(oracle.containsState, categorySizes, -1);

        std::printf("  oracle: %s tolerance=%.3g capacity=%.6g cm^3 (band=%.3g)\n",
                   oracle.valid ? "valid" : "*** NOT BRepCheck-VALID ***", oracle.tolerance,
                   oracle.capacity, oracleOpt.meshBand);
        json oracleJson;
        oracleJson["tolerance"] = oracle.tolerance;
        oracleJson["capacity"] = oracle.capacity;
        oracleJson["valid"] = oracle.valid;
        const double capacity = surf.Capacity();
        oracleJson["capacityCandidate"] = capacity;
        oracleJson["capacityRelativeDeviation"] =
          oracle.capacity != 0. ? (capacity - oracle.capacity) / oracle.capacity : 0.;
        std::printf("  oracle: capacity candidate=%.6g reference=%.6g relDev=%.3g\n", capacity,
                   oracle.capacity, oracleJson["capacityRelativeDeviation"].get<double>());

        if (opt.only.count("contains")) {
          auto v = validateContainsAgainstOracle(candidate, allPoints, containsState,
                                                 boundaryDistance, oracleOpt);
          printValidation("O:contains", v);
          oracleJson["contains"] = validationToJson(v);
        }
        const auto originStateFor = [&oracle](const char* category) {
          const auto it = oracle.originContains.find(category);
          return it == oracle.originContains.end() ? std::vector<int>{} : it->second;
        };
        if (opt.only.count("distout")) {
          const auto it = oracle.distOutside.find("outside");
          if (it != oracle.distOutside.end()) {
            auto v = validateDistanceAgainstOracle(candidate, samples.outsideRays, it->second,
                                                   /*wantInside=*/false, oracleOpt,
                                                   originStateFor("outside"));
            printValidation("O:distout", v);
            oracleJson["distout"] = validationToJson(v);
          }
        }
        if (opt.only.count("distin")) {
          const auto it = oracle.distInside.find("inside");
          if (it != oracle.distInside.end()) {
            auto v = validateDistanceAgainstOracle(candidate, samples.insideRays, it->second,
                                                   /*wantInside=*/true, oracleOpt,
                                                   originStateFor("inside"));
            printValidation("O:distin", v);
            oracleJson["distin"] = validationToJson(v);
          }
        }
        if (opt.only.count("safety")) {
          auto v = validateSafetyAgainstOracle(candidate, allPoints, boundaryDistance, oracleOpt);
          printValidation("O:safety", v);
          oracleJson["safety"] = validationToJson(v);
        }
        partJson["oracle"] = oracleJson;
      }
    }

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
      size_t crossingDumps = 0;
      constexpr size_t kMaxCrossingDumps = 3;
      std::vector<O2BVHSurfaceSolid::ContainsCrossing> bvhCrossings;
      std::vector<O2BVHSurfaceSolid::ContainsCrossing> loopCrossings;
      for (const auto& p : allPoints) {
        if (surf.Contains(p.data()) == surf.Contains_Loop(p.data())) {
          ++containsAgree;
          continue;
        }
        // A parity disagreement between two paths over the same kernels means the two hit lists
        // differ. Print them: the difference is the diagnosis, and guessing at it has already
        // cost this project one three-item plan built on a wrong premise.
        if (crossingDumps++ >= kMaxCrossingDumps) {
          continue;
        }
        surf.DescribeContainsCrossings(p, bvhCrossings, loopCrossings);
        std::printf("      BVH!=Loop at (%.9g,%.9g,%.9g): BVH=%d (%zu crossings) Loop=%d (%zu crossings)\n",
                   p[0], p[1], p[2], static_cast<int>(surf.Contains(p.data())), bvhCrossings.size(),
                   static_cast<int>(surf.Contains_Loop(p.data())), loopCrossings.size());
        const size_t nShow = std::max(bvhCrossings.size(), loopCrossings.size());
        for (size_t i = 0; i < nShow; ++i) {
          const char* bvhKind = i < bvhCrossings.size()
                                  ? (bvhCrossings[i].normalAlignment < 0. ? "ENTER" : "EXIT ") : "-----";
          const char* loopKind = i < loopCrossings.size()
                                   ? (loopCrossings[i].normalAlignment < 0. ? "ENTER" : "EXIT ") : "-----";
          const double bvhT = i < bvhCrossings.size() ? bvhCrossings[i].distance : -1.;
          const double loopT = i < loopCrossings.size() ? loopCrossings[i].distance : -1.;
          std::printf("        [%2zu] BVH %s t=%-18.12g   Loop %s t=%-18.12g%s\n", i, bvhKind, bvhT,
                     loopKind, loopT,
                     (i < bvhCrossings.size() && i < loopCrossings.size() &&
                      std::fabs(bvhT - loopT) > 1.e-12) ? "   <-- differs" : "");
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

  // Repeated at the end because per-part lines scroll away in a 19-part run, and because the whole
  // point of item 4 is that no future reader can attribute an "unexplained" column to mesh
  // chording without first seeing whether the subject was a closed manifold at all.
  if (!unreliableParts.empty()) {
    std::printf("\n*** %zu of %zu part(s) are NOT navigable; their accuracy columns above measure an\n"
               "*** undefined answer, not the exact solid's error. See scripts/geometry/ExactTrimTopology.md.\n",
               unreliableParts.size(), parts.size());
    for (const auto& id : unreliableParts) {
      std::printf("***   %s\n", id.c_str());
    }
  } else {
    std::printf("\nAll %zu part(s) closed consistently oriented manifolds: navigation results are meaningful.\n",
               parts.size());
  }

  if (!opt.jsonOut.empty()) {
    std::ofstream out(opt.jsonOut);
    out << jsonReport.dump(1);
    std::printf("\nWrote %s\n", opt.jsonOut.c_str());
  }

  return 0;
}
