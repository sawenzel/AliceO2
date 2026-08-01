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

/// \file O2SolidHarness.h
/// \brief Reusable validation / performance-comparison harness for TGeoShape navigation.
///
/// Design and usage are documented in scripts/geometry/SolidNavigationHarness.md; the ground
/// rules there (tessellated mesh is a reference not the truth; never compare Safety() between two
/// shapes; classify distance mismatches via the reference's own Safety band) are load-bearing and
/// this header follows them directly. Everything here is typed on plain `TGeoShape*` so the same
/// code compares O2BVHSurfaceSolid against O2Tessellated *and* against ROOT primitives
/// (TGeoBBox/TGeoTube/...) in the unit test.

#ifndef ALICEO2_BASE_O2SOLIDHARNESS_
#define ALICEO2_BASE_O2SOLIDHARNESS_

#include "TGeoShape.h"

#include <array>
#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

namespace o2
{
namespace base
{
namespace harness
{

using Point3D = std::array<double, 3>;

struct Ray {
  Point3D origin{};
  Point3D dir{}; // unit vector by convention (TGeo contract); not renormalized by the harness
};

/// Parameters controlling `generateSamples`. All counts are targets: rejection-sampled
/// categories (boundary/inside/rays) may return fewer if `maxRejectionAttempts` is exhausted,
/// which is reported as a smaller vector rather than a silent hang or a thrown error.
struct SampleConfig {
  int nBulk = 2000;              ///< uniform points over the inflated bbox
  int nBoundary = 2000;          ///< points within `boundaryBand` of the reference surface
  int nInside = 1000;            ///< points accepted by the reference Contains()
  int nOutsideRays = 4000;       ///< rays from outside origins, for DistFromOutside
  int nInsideRays = 2000;        ///< rays from inside origins, for DistFromInside
  double bboxInflate = 0.15;     ///< fractional bbox half-extent padding for bulk/outside sampling
  double boundaryBand = -1.;     ///< absolute distance (cm); <0 auto-picks 1e-3 * bbox diagonal
  double aimedRayFraction = 0.5; ///< fraction of rays aimed at a random interior bbox point rather
                                  ///< than an isotropic direction (keeps DistFromOutside hit rates
                                  ///< non-degenerate -- see SolidNavigationHarness.md, Step 2)
  int maxRejectionAttempts = 200;///< attempts per accepted sample before giving up on that category
  uint64_t seed = 1;             ///< every SampleSet is fully determined by this and the bbox
};

struct SampleSet {
  Point3D bboxMin{};
  Point3D bboxMax{};
  std::vector<Point3D> bulkPoints;
  std::vector<Point3D> boundaryPoints;
  std::vector<Point3D> insidePoints;
  std::vector<Ray> outsideRays;
  std::vector<Ray> insideRays;
};

/// Build a deterministic sample set from `cfg.seed` and the given bbox. `reference` classifies
/// bulk points into the boundary-band / inside-biased subsets via its own Contains()/Safety(), and
/// picks outside-ray origins via Contains(); it should be the trusted reference shape (the
/// tessellated solid), never the shape under test -- see the ground rules in
/// SolidNavigationHarness.md.
SampleSet generateSamples(const TGeoShape* reference, const Point3D& bboxMin, const Point3D& bboxMax,
                          const SampleConfig& cfg = {});

// ------------------------------------------------------------------------------------------
// Validation
// ------------------------------------------------------------------------------------------

/// One worst-case disagreement, with enough state (point/direction/values) to reproduce it
/// directly outside the harness.
struct Offender {
  Point3D point{};
  Point3D dir{};               // zero for point-only queries (Contains, Safety)
  double candidateValue = 0.;
  double referenceValue = 0.;
  double deviation = 0.;
  double referenceSafety = 0.;   // point queries: reference distance to its own surface
  double incidenceCosine = 1.;   // ray queries: |cos| between ray and surface normal at the hit,
                                  // i.e. how much surface uncertainty this ray amplifies
};

struct ValidationResult {
  size_t nSamples = 0;
  size_t nAgree = 0;
  size_t nMismatchWithinBand = 0; // explainable by the reference's own imprecision (see below)
  size_t nMismatchMissedSurface = 0; // one side found no crossing where the other did
  size_t nMismatchUnexplained = 0;
  size_t nNoVerdict = 0;          // oracle mode only: the reference declined to answer
  size_t nRelabelled = 0;         // ray queries, oracle mode: origins whose category the oracle
                                  // contradicted, so the other TGeo entry point was asked
  double worstDeviation = 0.;
  std::vector<Offender> worstOffenders; // bounded by opt.maxOffenders, worst-first
};

/// `nMismatchMissedSurface` exists because the previous two-way split silently absorbed the worst
/// failure mode there is. When the candidate returns Big and the reference hits, the old
/// classifier probed at the *reference's own* hit point, where the reference Safety is ~0 by
/// construction, so "candidate missed an entire wall" was always counted as "explained by mesh
/// chording". The same held for tunneling: a candidate that skips the near wall and stops on a
/// far one also stops *on the reference mesh*. Both are now their own category and can never be
/// explained away. See CodeReview_Fable.md, finding S9.
struct ValidationOptions {
  double distanceTolerance = 1.e-6; ///< absolute agreement tolerance for distances (cm)
  double meshBand = 1.e-2;          ///< the reference surface's own positional uncertainty (cm).
                                     ///< For the tessellated reference this is the chord sagitta;
                                     ///< in oracle mode it is the model's declared tolerance.
  /// A surface displaced by `meshBand` moves a ray's crossing by `meshBand / |cos(incidence)|`,
  /// so a grazing ray legitimately shows a much larger distance difference than a perpendicular
  /// one. The allowed difference is scaled by the measured incidence, floored here so a tangent
  /// ray cannot excuse an unbounded error.
  double minIncidenceCosine = 1.e-2;
  double stepmax = TGeoShape::Big();
  size_t maxOffenders = 10;
};

ValidationResult validateContains(const TGeoShape* candidate, const TGeoShape* reference,
                                  const std::vector<Point3D>& points, const ValidationOptions& opt = {});

ValidationResult validateDistFromOutside(const TGeoShape* candidate, const TGeoShape* reference,
                                         const std::vector<Ray>& rays, const ValidationOptions& opt = {});

ValidationResult validateDistFromInside(const TGeoShape* candidate, const TGeoShape* reference,
                                        const std::vector<Ray>& rays, const ValidationOptions& opt = {});

/// Checks the Safety() lower-bound contract of a *single* shape against its own
/// DistFrom{Inside,Outside} along a handful of deterministic probe directions. Never compares
/// Safety() between two shapes -- per the ground rules, two correct implementations may disagree
/// on Safety arbitrarily; only `safety <= true distance` and `safety >= 0` are contracted. Here
/// `candidateValue` is the reported safety and `referenceValue` is the smallest probed distance
/// (a violation is candidateValue > referenceValue + tolerance); `nMismatchWithinBand` is unused
/// (always 0) since there is no mesh-band concept for a single shape's own contract.
ValidationResult validateSafety(const TGeoShape* shape, const std::vector<Point3D>& points,
                                const ValidationOptions& opt = {});

// ------------------------------------------------------------------------------------------
// Validation against an external oracle (OpenCascade)
// ------------------------------------------------------------------------------------------
//
// The tessellated reference above can only ever say "you differ from a chorded approximation of
// the truth". These entry points take answers computed by scripts/geometry/occtOracle.py on the
// *same* BREP the sidecar was extracted from, so a disagreement outside the model's declared
// tolerance is a defect rather than a discussion. The oracle's per-point boundary distance is
// what makes the band principled: a point closer to the boundary than the model tolerance has no
// defined answer and is counted as `nNoVerdict` instead of being scored either way.

/// `oracleState`: 1 inside, 0 outside, -1 the oracle declined (point ON the boundary).
/// `oracleBoundaryDistance` may be shorter than `points` (the oracle caps that expensive query);
/// points without a distance are scored, but only their -1 state can abstain.
ValidationResult validateContainsAgainstOracle(const TGeoShape* candidate,
                                               const std::vector<Point3D>& points,
                                               const std::vector<int>& oracleState,
                                               const std::vector<double>& oracleBoundaryDistance,
                                               const ValidationOptions& opt = {});

/// `oracleDistance` holds the nearest positive ray/boundary crossing, or >= TGeoShape::Big() for
/// a miss. The oracle uses one computation for both directions (see occtOracle.py), so which TGeo
/// entry point the candidate must be asked is decided by where the *origin* is -- and that is the
/// one thing the sample generator cannot be trusted about.
///
/// `oracleOriginState` is the oracle's own classification of each ray origin: 1 inside, 0 outside,
/// -1 ON the boundary. When it is present it decides the entry point, per ray, overriding
/// `wantInside`; a -1 origin abstains (`nNoVerdict`) because neither entry point is defined there.
/// When it is empty -- an older answer file -- `wantInside` decides as before.
///
/// This exists because the generator labels a ray from the *tessellated* reference, and that mesh
/// can be wrong: on Bagger/BucketLink2 it is not watertight and its left plate sits 0.2 cm from
/// the BREP's, so all 20 recorded offenders had the opposite category, the harness asked the
/// candidate for an entry distance from a point that was inside, compared it against an exit
/// distance, and booked a correct answer as a missed surface. Re-labelling is sound precisely
/// because it never consults the candidate: the oracle is ground truth for its own origins.
/// See CodeReview_Fable.md, finding H1 in Section 13.
ValidationResult validateDistanceAgainstOracle(const TGeoShape* candidate,
                                               const std::vector<Ray>& rays,
                                               const std::vector<double>& oracleDistance,
                                               bool wantInside, const ValidationOptions& opt = {},
                                               const std::vector<int>& oracleOriginState = {});

/// Safety's contract against ground truth: `0 <= safety <= trueDistance`. Unlike
/// `validateSafety`, which can only probe six directions of the shape's own distance functions,
/// the oracle supplies the exact distance to the boundary, so an over-estimate cannot hide
/// between the probes.
ValidationResult validateSafetyAgainstOracle(const TGeoShape* candidate,
                                             const std::vector<Point3D>& points,
                                             const std::vector<double>& oracleBoundaryDistance,
                                             const ValidationOptions& opt = {});

// ------------------------------------------------------------------------------------------
// Timing
// ------------------------------------------------------------------------------------------

struct TimingResult {
  size_t nCalls = 0;
  double nsPerCall = 0.;
  uint64_t checksum = 0; ///< accumulated from the results so the optimizer cannot elide the calls
};

namespace detail
{
/// Checksum mixer the timing loops accumulate results through, so the optimizer cannot elide the
/// measured calls. Exposed only because `timeRayKernel` below is a template.
uint64_t mixDouble(uint64_t acc, double value);
} // namespace detail

/// Time an arbitrary per-ray distance kernel with exactly the methodology of the `timeDistFrom*`
/// functions below (same warmup, same timed repeats, same checksum accumulation), so the two are
/// directly comparable. `kernel(origin, dir)` returns a distance in cm.
///
/// This exists so a caller can time entry points this TGeoShape-typed header knows nothing about
/// -- O2BVHSurfaceSolid's `_Loop` baselines, or the same query with ray tmax pruning switched off
/// -- against the ordinary virtual ones without reimplementing the measurement.
template <typename RayKernel>
TimingResult timeRayKernel(const std::vector<Ray>& rays, int warmupRepeats, int timedRepeats, RayKernel&& kernel)
{
  for (int warmup = 0; warmup < warmupRepeats; ++warmup) {
    for (const auto& ray : rays) {
      volatile double sink = kernel(ray.origin, ray.dir);
      (void)sink;
    }
  }
  uint64_t checksum = 0;
  const auto start = std::chrono::steady_clock::now();
  for (int repeat = 0; repeat < timedRepeats; ++repeat) {
    for (const auto& ray : rays) {
      checksum = detail::mixDouble(checksum, kernel(ray.origin, ray.dir));
    }
  }
  const auto stop = std::chrono::steady_clock::now();
  TimingResult result;
  result.nCalls = rays.size() * static_cast<size_t>(timedRepeats);
  const double nanoseconds = std::chrono::duration<double, std::nano>(stop - start).count();
  result.nsPerCall = result.nCalls > 0 ? nanoseconds / static_cast<double>(result.nCalls) : 0.;
  result.checksum = checksum;
  return result;
}

TimingResult timeContains(const TGeoShape* shape, const std::vector<Point3D>& points, int warmupRepeats,
                          int timedRepeats);
TimingResult timeDistFromOutside(const TGeoShape* shape, const std::vector<Ray>& rays, int warmupRepeats,
                                 int timedRepeats, double stepmax = TGeoShape::Big());
TimingResult timeDistFromInside(const TGeoShape* shape, const std::vector<Ray>& rays, int warmupRepeats,
                                int timedRepeats);
TimingResult timeSafety(const TGeoShape* shape, const std::vector<Point3D>& points, int warmupRepeats,
                        int timedRepeats);

/// Per-part structural stats. The harness does not compute these itself -- they come from each
/// shape's own introspection API (O2BVHSurfaceSolid::GetNsurfaces/CountBVHRayCandidates,
/// O2Tessellated::GetNfacets, ...) and externally measured build timings -- but the type lives
/// here so callers report a consistent shape. See SolidNavigationHarness.md, Step 2: primitive
/// count and build cost are half the comparison.
struct PartStats {
  std::string partId;
  int nSurfaces = 0;
  int nTriangles = 0;
  double closeShapeSeconds = 0.;
  long long bvhRayCandidatesSampled = 0; // summed over a probe ray set, caller-defined
};

// ------------------------------------------------------------------------------------------
// The `shape_<part>.root` sidecar: handing an arbitrary TGeoShape to the harness and the gate
// ------------------------------------------------------------------------------------------
//
// The gate could previously score exactly one thing: an O2BVHSurfaceSolid loaded from
// `surfaces_<part>.bin`. The four scored queries are TGeoShape virtuals, so nothing but the
// loading was actually specific to that class. These two functions are the third way to obtain a
// candidate -- a plain ROOT-serialised TGeoShape -- which is how the planned CSG emitter will
// hand over TGeoTube / TGeoBBox / TGeoCompositeShape trees.
//
// The convention, stated in full so an emitter can be written against it without reading this
// code (and again in scripts/geometry/Stream_G_AnyShape.md):
//
//   * one file per part, named `shape_<VOL>_<LID>.root`, sitting in the same converter output
//     directory as `surfaces_<VOL>_<LID>.bin`, `facets_<VOL>_<LID>.bin` and
//     `brep_<VOL>_<LID>.brep`;
//   * it holds exactly one object inheriting from TGeoShape, stored under the key name "shape".
//     A file whose "shape" key is missing is searched for the first TGeoShape-derived key, so a
//     hand-written file with a different key name still loads, but emitters must write "shape";
//   * lengths are in CENTIMETRES and the shape is expressed in the part's own LOCAL frame -- the
//     same frame as the sidecar, the mesh and the .brep, i.e. the leaf solid as the converter
//     emits it, with no placement matrix applied. This is the whole correctness precondition:
//     the oracle answers questions about the .brep in that frame, and a TGeoShape answers in its
//     own, so a shape written in the assembly frame produces a full table of plausible nonsense.
//     The harness measures the deviation between the shape's bounding box and the oracle's
//     (`bboxDeviationFromOracle`) and reports it per representation so the mistake is visible
//     rather than silent;
//   * a TGeoCompositeShape is written whole -- its TGeoBoolNode and component shapes stream with
//     it -- and needs no TGeoManager on either side.

/// Read the single TGeoShape from a `shape_<part>.root` sidecar. Returns nullptr on any failure,
/// with a human-readable reason in `*error` when that is non-null. The caller owns the result;
/// it is detached from the file, which is closed before returning.
TGeoShape* loadShapeFromRootFile(const std::string& path, std::string* error = nullptr);

/// Write one, so that producer and consumer of the convention cannot drift apart: the unit test
/// and any fixture generator go through the same function the harness reads with.
bool saveShapeToRootFile(const std::string& path, const TGeoShape& shape, std::string* error = nullptr);

} // namespace harness
} // namespace base
} // namespace o2

#endif
