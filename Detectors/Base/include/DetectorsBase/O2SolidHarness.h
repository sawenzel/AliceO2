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
  double referenceSafety = 0.; // the value used to classify this offender as mesh-band or not
};

struct ValidationResult {
  size_t nSamples = 0;
  size_t nAgree = 0;
  size_t nMismatchWithinBand = 0; // |referenceSafety| < opt.meshBand: explained by mesh chording
  size_t nMismatchUnexplained = 0;
  double worstDeviation = 0.;
  std::vector<Offender> worstOffenders; // bounded by opt.maxOffenders, worst-first
};

struct ValidationOptions {
  double distanceTolerance = 1.e-6; ///< absolute agreement tolerance for distances (cm)
  double meshBand = 1.e-2;          ///< |reference Safety| below this classifies a mismatch as
                                     ///< mesh-band rather than unexplained (cm)
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

} // namespace harness
} // namespace base
} // namespace o2

#endif
