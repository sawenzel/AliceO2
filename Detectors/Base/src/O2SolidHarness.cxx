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

#include "DetectorsBase/O2SolidHarness.h"

#include "TClass.h"
#include "TFile.h"
#include "TGeoMatrix.h"
#include "TKey.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <memory>
#include <random>

namespace o2
{
namespace base
{
namespace harness
{

namespace
{

// Both DistFrom{Outside,Inside} implementations in this project (O2Tessellated, and the BVH
// O2BVHSurfaceSolid) ignore `iact` beyond the "compute distance" case; iact=3 is the convention
// O2Tessellated documents as the one it implements. Used uniformly here for both shapes.
constexpr Int_t kIact = 3;

Point3D add(const Point3D& a, const Point3D& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Point3D sub(const Point3D& a, const Point3D& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Point3D scale(const Point3D& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double normSq(const Point3D& a) { return a[0] * a[0] + a[1] * a[1] + a[2] * a[2]; }

Point3D sampleUniform(std::mt19937_64& rng, const Point3D& lo, const Point3D& hi)
{
  std::uniform_real_distribution<double> ux(lo[0], hi[0]);
  std::uniform_real_distribution<double> uy(lo[1], hi[1]);
  std::uniform_real_distribution<double> uz(lo[2], hi[2]);
  return {ux(rng), uy(rng), uz(rng)};
}

Point3D isotropicDir(std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> uCos(-1., 1.);
  std::uniform_real_distribution<double> uPhi(0., 2. * M_PI);
  const double cosTheta = uCos(rng);
  const double sinTheta = std::sqrt(std::max(0., 1. - cosTheta * cosTheta));
  const double phi = uPhi(rng);
  return {sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
}

bool isBig(double d) { return d >= 0.9 * TGeoShape::Big(); }

} // namespace

namespace detail
{
uint64_t mixDouble(uint64_t acc, double value)
{
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  bits += 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
  return acc ^ bits;
}
} // namespace detail

SampleSet generateSamples(const TGeoShape* reference, const Point3D& bboxMin, const Point3D& bboxMax,
                          const SampleConfig& cfg)
{
  SampleSet out;
  out.bboxMin = bboxMin;
  out.bboxMax = bboxMax;

  const Point3D center = scale(add(bboxMin, bboxMax), 0.5);
  const Point3D halfExtent = scale(sub(bboxMax, bboxMin), 0.5);
  const Point3D inflatedLo = sub(center, scale(halfExtent, 1. + cfg.bboxInflate));
  const Point3D inflatedHi = add(center, scale(halfExtent, 1. + cfg.bboxInflate));

  double band = cfg.boundaryBand;
  if (band < 0.) {
    const double diag = std::sqrt(normSq(sub(bboxMax, bboxMin)));
    band = 1.e-3 * diag;
  }

  std::mt19937_64 rng(cfg.seed);

  out.bulkPoints.reserve(cfg.nBulk);
  for (int i = 0; i < cfg.nBulk; ++i) {
    out.bulkPoints.push_back(sampleUniform(rng, inflatedLo, inflatedHi));
  }

  out.boundaryPoints.reserve(cfg.nBoundary);
  {
    const long long budget = static_cast<long long>(cfg.nBoundary) * cfg.maxRejectionAttempts;
    long long attempts = 0;
    while (static_cast<int>(out.boundaryPoints.size()) < cfg.nBoundary && attempts < budget) {
      ++attempts;
      const Point3D p = sampleUniform(rng, bboxMin, bboxMax);
      const bool in = reference->Contains(p.data());
      const double s = reference->Safety(p.data(), in);
      if (s < band) {
        out.boundaryPoints.push_back(p);
      }
    }
  }

  out.insidePoints.reserve(cfg.nInside);
  {
    const long long budget = static_cast<long long>(cfg.nInside) * cfg.maxRejectionAttempts;
    long long attempts = 0;
    while (static_cast<int>(out.insidePoints.size()) < cfg.nInside && attempts < budget) {
      ++attempts;
      const Point3D p = sampleUniform(rng, bboxMin, bboxMax);
      if (reference->Contains(p.data())) {
        out.insidePoints.push_back(p);
      }
    }
  }

  out.outsideRays.reserve(cfg.nOutsideRays);
  {
    std::uniform_real_distribution<double> u01(0., 1.);
    const long long budget = static_cast<long long>(cfg.nOutsideRays) * cfg.maxRejectionAttempts;
    long long attempts = 0;
    while (static_cast<int>(out.outsideRays.size()) < cfg.nOutsideRays && attempts < budget) {
      ++attempts;
      const Point3D origin = sampleUniform(rng, inflatedLo, inflatedHi);
      if (reference->Contains(origin.data())) {
        continue;
      }
      Point3D dir;
      if (u01(rng) < cfg.aimedRayFraction) {
        Point3D target = sampleUniform(rng, bboxMin, bboxMax);
        Point3D delta = sub(target, origin);
        double len = std::sqrt(normSq(delta));
        if (len < 1.e-12) {
          dir = isotropicDir(rng);
        } else {
          dir = scale(delta, 1. / len);
        }
      } else {
        dir = isotropicDir(rng);
      }
      out.outsideRays.push_back({origin, dir});
    }
  }

  out.insideRays.reserve(cfg.nInsideRays);
  {
    const long long budget = static_cast<long long>(cfg.nInsideRays) * cfg.maxRejectionAttempts;
    long long attempts = 0;
    while (static_cast<int>(out.insideRays.size()) < cfg.nInsideRays && attempts < budget) {
      ++attempts;
      const Point3D origin = sampleUniform(rng, bboxMin, bboxMax);
      if (!reference->Contains(origin.data())) {
        continue;
      }
      out.insideRays.push_back({origin, isotropicDir(rng)});
    }
  }

  return out;
}

// ------------------------------------------------------------------------------------------
// Validation
// ------------------------------------------------------------------------------------------

namespace
{

enum class MismatchClass { WithinBand, MissedSurface, Unexplained };

void recordOffender(ValidationResult& result, const ValidationOptions& opt, Offender&& off,
                    MismatchClass mismatchClass)
{
  switch (mismatchClass) {
    case MismatchClass::WithinBand:
      ++result.nMismatchWithinBand;
      break;
    case MismatchClass::MissedSurface:
      ++result.nMismatchMissedSurface;
      break;
    case MismatchClass::Unexplained:
      ++result.nMismatchUnexplained;
      break;
  }
  result.worstDeviation = std::max(result.worstDeviation, std::fabs(off.deviation));
  result.worstOffenders.push_back(std::move(off));
  std::sort(result.worstOffenders.begin(), result.worstOffenders.end(),
           [](const Offender& a, const Offender& b) { return std::fabs(a.deviation) > std::fabs(b.deviation); });
  if (result.worstOffenders.size() > opt.maxOffenders) {
    result.worstOffenders.resize(opt.maxOffenders);
  }
}

/// How far a ray's crossing may legitimately move when the reference surface itself is uncertain
/// by `opt.meshBand`. A surface displaced by d shifts the crossing of a ray meeting it at
/// incidence angle theta by d/|cos theta|, so a grazing hit is allowed much more slack than a
/// perpendicular one. `normalSource` supplies the surface normal at the probe point; when it
/// cannot (no shape, or a degenerate normal) the perpendicular case is assumed, which is the
/// strict end of the range.
double allowedCrossingShift(const TGeoShape* normalSource, const Point3D& probePoint,
                            const Point3D& dir, const ValidationOptions& opt, double& cosIncidence)
{
  cosIncidence = 1.;
  if (normalSource != nullptr) {
    double normal[3] = {0., 0., 0.};
    // ComputeNormal takes non-const pointers in the TGeo API but does not modify the inputs.
    auto* mutableShape = const_cast<TGeoShape*>(normalSource);
    mutableShape->ComputeNormal(probePoint.data(), dir.data(), normal);
    const double normalNorm =
      std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (normalNorm > 0.) {
      const double dotProduct =
        (normal[0] * dir[0] + normal[1] * dir[1] + normal[2] * dir[2]) / normalNorm;
      cosIncidence = std::fabs(dotProduct);
    }
  }
  const double effectiveCosine = std::max(cosIncidence, opt.minIncidenceCosine);
  return std::max(opt.distanceTolerance, opt.meshBand / effectiveCosine);
}

/// Shared classification for both distance queries. `dc`/`dr` are the candidate and reference
/// distances; `reference` (may be null) is used only to measure the incidence angle.
MismatchClass classifyDistanceMismatch(const TGeoShape* reference, const Ray& ray, double dc,
                                       double dr, bool dcBig, bool drBig,
                                       const ValidationOptions& opt, double& cosIncidence)
{
  cosIncidence = 1.;
  // One side found a crossing where the other found none. No amount of surface uncertainty
  // explains a missing wall, so this can never be counted as "within band".
  if (dcBig != drBig) {
    return MismatchClass::MissedSurface;
  }
  const Point3D probePoint = add(ray.origin, scale(ray.dir, dr));
  const double allowed = allowedCrossingShift(reference, probePoint, ray.dir, opt, cosIncidence);
  return std::fabs(dc - dr) <= allowed ? MismatchClass::WithinBand : MismatchClass::Unexplained;
}

} // namespace

ValidationResult validateContains(const TGeoShape* candidate, const TGeoShape* reference,
                                  const std::vector<Point3D>& points, const ValidationOptions& opt)
{
  ValidationResult result;
  result.nSamples = points.size();
  for (const auto& p : points) {
    const bool bc = candidate->Contains(p.data());
    const bool br = reference->Contains(p.data());
    if (bc == br) {
      ++result.nAgree;
      continue;
    }
    const double refSafety = reference->Safety(p.data(), br);
    Offender off;
    off.point = p;
    off.candidateValue = bc ? 1. : 0.;
    off.referenceValue = br ? 1. : 0.;
    off.deviation = refSafety; // rank Contains mismatches by how deep into the "unambiguous" region they are
    off.referenceSafety = refSafety;
    // A point closer to the reference surface than the reference's own positional uncertainty
    // genuinely has no defined reference answer; further out, the reference is authoritative.
    recordOffender(result, opt, std::move(off),
                   refSafety < opt.meshBand ? MismatchClass::WithinBand : MismatchClass::Unexplained);
  }
  return result;
}

ValidationResult validateDistFromOutside(const TGeoShape* candidate, const TGeoShape* reference,
                                         const std::vector<Ray>& rays, const ValidationOptions& opt)
{
  ValidationResult result;
  result.nSamples = rays.size();
  for (const auto& r : rays) {
    const double dc = candidate->DistFromOutside(r.origin.data(), r.dir.data(), kIact, opt.stepmax);
    const double dr = reference->DistFromOutside(r.origin.data(), r.dir.data(), kIact, opt.stepmax);
    const bool dcBig = isBig(dc);
    const bool drBig = isBig(dr);
    if (dcBig && drBig) {
      ++result.nAgree;
      continue;
    }
    if (!dcBig && !drBig && std::fabs(dc - dr) <= opt.distanceTolerance) {
      ++result.nAgree;
      continue;
    }
    double cosIncidence = 1.;
    const MismatchClass mismatchClass =
      classifyDistanceMismatch(reference, r, dc, dr, dcBig, drBig, opt, cosIncidence);
    Offender off;
    off.point = r.origin;
    off.dir = r.dir;
    off.candidateValue = dcBig ? opt.stepmax : dc;
    off.referenceValue = drBig ? opt.stepmax : dr;
    off.deviation = off.candidateValue - off.referenceValue;
    off.incidenceCosine = cosIncidence;
    recordOffender(result, opt, std::move(off), mismatchClass);
  }
  return result;
}

ValidationResult validateDistFromInside(const TGeoShape* candidate, const TGeoShape* reference,
                                        const std::vector<Ray>& rays, const ValidationOptions& opt)
{
  ValidationResult result;
  result.nSamples = rays.size();
  for (const auto& r : rays) {
    const double dc = candidate->DistFromInside(r.origin.data(), r.dir.data(), kIact, opt.stepmax);
    const double dr = reference->DistFromInside(r.origin.data(), r.dir.data(), kIact, opt.stepmax);
    const bool dcBig = isBig(dc);
    const bool drBig = isBig(dr);
    if (dcBig && drBig) {
      ++result.nAgree;
      continue;
    }
    if (!dcBig && !drBig && std::fabs(dc - dr) <= opt.distanceTolerance) {
      ++result.nAgree;
      continue;
    }
    double cosIncidence = 1.;
    const MismatchClass mismatchClass =
      classifyDistanceMismatch(reference, r, dc, dr, dcBig, drBig, opt, cosIncidence);
    Offender off;
    off.point = r.origin;
    off.dir = r.dir;
    off.candidateValue = dcBig ? opt.stepmax : dc;
    off.referenceValue = drBig ? opt.stepmax : dr;
    off.deviation = off.candidateValue - off.referenceValue;
    off.incidenceCosine = cosIncidence;
    recordOffender(result, opt, std::move(off), mismatchClass);
  }
  return result;
}

ValidationResult validateSafety(const TGeoShape* shape, const std::vector<Point3D>& points,
                                const ValidationOptions& opt)
{
  static const std::array<Point3D, 6> kProbeDirs = {
    Point3D{1., 0., 0.}, Point3D{-1., 0., 0.}, Point3D{0., 1., 0.},
    Point3D{0., -1., 0.}, Point3D{0., 0., 1.}, Point3D{0., 0., -1.}};

  ValidationResult result;
  result.nSamples = points.size();
  for (const auto& p : points) {
    const bool in = shape->Contains(p.data());
    const double s = shape->Safety(p.data(), in);

    double minProbed = TGeoShape::Big();
    for (const auto& d : kProbeDirs) {
      const double dist = in ? shape->DistFromInside(p.data(), d.data(), kIact, opt.stepmax)
                             : shape->DistFromOutside(p.data(), d.data(), kIact, opt.stepmax);
      minProbed = std::min(minProbed, isBig(dist) ? opt.stepmax : dist);
    }

    const bool violatesLowerBound = s < -opt.distanceTolerance;
    const bool violatesUpperBound = s > minProbed + opt.distanceTolerance;
    if (!violatesLowerBound && !violatesUpperBound) {
      ++result.nAgree;
      continue;
    }
    Offender off;
    off.point = p;
    off.candidateValue = s;
    off.referenceValue = minProbed;
    off.deviation = s - minProbed;
    off.referenceSafety = s;
    recordOffender(result, opt, std::move(off), MismatchClass::Unexplained);
  }
  return result;
}

// ------------------------------------------------------------------------------------------
// Validation against an external oracle
// ------------------------------------------------------------------------------------------

namespace
{
/// The oracle computes the exact boundary distance only for a capped prefix of each category
/// (it is the most expensive query it has). Points beyond the cap return a negative value,
/// meaning "unknown", and are then scored without a distance-based abstention.
constexpr double kUnknownDistance = -1.;

double oracleDistanceAt(const std::vector<double>& distances, size_t index)
{
  return index < distances.size() ? distances[index] : kUnknownDistance;
}
} // namespace

ValidationResult validateContainsAgainstOracle(const TGeoShape* candidate,
                                               const std::vector<Point3D>& points,
                                               const std::vector<int>& oracleState,
                                               const std::vector<double>& oracleBoundaryDistance,
                                               const ValidationOptions& opt)
{
  ValidationResult result;
  result.nSamples = points.size();
  for (size_t index = 0; index < points.size(); ++index) {
    const int state = index < oracleState.size() ? oracleState[index] : -1;
    const double boundaryDistance = oracleDistanceAt(oracleBoundaryDistance, index);
    // Two ways the oracle abstains: it classified the point as ON the boundary, or the point is
    // closer to the boundary than the model's own tolerance. Neither is a defect in the
    // candidate -- the truth itself is undefined there.
    if (state < 0 || (boundaryDistance >= 0. && boundaryDistance < opt.meshBand)) {
      ++result.nNoVerdict;
      continue;
    }
    const bool candidateInside = candidate->Contains(points[index].data());
    if (candidateInside == (state == 1)) {
      ++result.nAgree;
      continue;
    }
    Offender off;
    off.point = points[index];
    off.candidateValue = candidateInside ? 1. : 0.;
    off.referenceValue = state == 1 ? 1. : 0.;
    // Rank by how far into unambiguous territory the disagreement sits: a wrong answer 1 cm from
    // any surface is a different animal from one 1 um away.
    off.deviation = boundaryDistance >= 0. ? boundaryDistance : 0.;
    off.referenceSafety = boundaryDistance;
    recordOffender(result, opt, std::move(off), MismatchClass::Unexplained);
  }
  return result;
}

ValidationResult validateDistanceAgainstOracle(const TGeoShape* candidate,
                                               const std::vector<Ray>& rays,
                                               const std::vector<double>& oracleDistance,
                                               bool wantInside, const ValidationOptions& opt,
                                               const std::vector<int>& oracleOriginState)
{
  ValidationResult result;
  result.nSamples = rays.size();
  for (size_t index = 0; index < rays.size(); ++index) {
    if (index >= oracleDistance.size()) {
      ++result.nNoVerdict;
      continue;
    }
    // The origin's own classification, from the oracle, decides which entry point is defined
    // here. The category the sample generator assigned is only a hint, and on a non-watertight
    // reference mesh it is a wrong one.
    bool askInside = wantInside;
    if (index < oracleOriginState.size()) {
      const int state = oracleOriginState[index];
      if (state < 0) {
        ++result.nNoVerdict; // origin ON the boundary: neither entry point is defined
        continue;
      }
      askInside = state == 1;
      if (askInside != wantInside) {
        ++result.nRelabelled;
      }
    }
    const auto& ray = rays[index];
    const double dc = askInside
                        ? candidate->DistFromInside(ray.origin.data(), ray.dir.data(), kIact, opt.stepmax)
                        : candidate->DistFromOutside(ray.origin.data(), ray.dir.data(), kIact, opt.stepmax);
    const double dr = oracleDistance[index];
    const bool dcBig = isBig(dc);
    const bool drBig = isBig(dr);
    if (dcBig && drBig) {
      ++result.nAgree;
      continue;
    }
    if (!dcBig && !drBig && std::fabs(dc - dr) <= opt.distanceTolerance) {
      ++result.nAgree;
      continue;
    }
    // No reference shape to take a normal from, so the incidence cannot be measured here and the
    // strict (perpendicular) allowance applies. That is the right default against ground truth:
    // the oracle's crossing is exact, only the model tolerance is in play.
    double cosIncidence = 1.;
    const MismatchClass mismatchClass =
      classifyDistanceMismatch(nullptr, ray, dc, dr, dcBig, drBig, opt, cosIncidence);
    Offender off;
    off.point = ray.origin;
    off.dir = ray.dir;
    off.candidateValue = dcBig ? opt.stepmax : dc;
    off.referenceValue = drBig ? opt.stepmax : dr;
    off.deviation = off.candidateValue - off.referenceValue;
    off.incidenceCosine = cosIncidence;
    recordOffender(result, opt, std::move(off), mismatchClass);
  }
  return result;
}

ValidationResult validateSafetyAgainstOracle(const TGeoShape* candidate,
                                             const std::vector<Point3D>& points,
                                             const std::vector<double>& oracleBoundaryDistance,
                                             const ValidationOptions& opt)
{
  ValidationResult result;
  result.nSamples = points.size();
  for (size_t index = 0; index < points.size(); ++index) {
    const double trueDistance = oracleDistanceAt(oracleBoundaryDistance, index);
    if (trueDistance < 0.) {
      ++result.nNoVerdict;
      continue;
    }
    const bool inside = candidate->Contains(points[index].data());
    const double safety = candidate->Safety(points[index].data(), inside);
    // Safety is contractually a lower bound on the true distance, and never negative. Being too
    // small is legal (if pessimistic); being too large is a navigation bug that lets a track step
    // through a surface.
    const bool violatesLowerBound = safety < -opt.distanceTolerance;
    const bool violatesUpperBound = safety > trueDistance + opt.distanceTolerance;
    if (!violatesLowerBound && !violatesUpperBound) {
      ++result.nAgree;
      continue;
    }
    Offender off;
    off.point = points[index];
    off.candidateValue = safety;
    off.referenceValue = trueDistance;
    off.deviation = safety - trueDistance;
    off.referenceSafety = trueDistance;
    recordOffender(result, opt, std::move(off), MismatchClass::Unexplained);
  }
  return result;
}

// ------------------------------------------------------------------------------------------
// Timing
// ------------------------------------------------------------------------------------------

TimingResult timeContains(const TGeoShape* shape, const std::vector<Point3D>& points, int warmupRepeats,
                          int timedRepeats)
{
  for (int w = 0; w < warmupRepeats; ++w) {
    for (const auto& p : points) {
      volatile bool sink = shape->Contains(p.data());
      (void)sink;
    }
  }
  uint64_t checksum = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < timedRepeats; ++r) {
    for (const auto& p : points) {
      checksum = detail::mixDouble(checksum, shape->Contains(p.data()) ? 1. : 0.);
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  TimingResult result;
  result.nCalls = points.size() * static_cast<size_t>(timedRepeats);
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  result.nsPerCall = result.nCalls > 0 ? ns / static_cast<double>(result.nCalls) : 0.;
  result.checksum = checksum;
  return result;
}

TimingResult timeDistFromOutside(const TGeoShape* shape, const std::vector<Ray>& rays, int warmupRepeats,
                                 int timedRepeats, double stepmax)
{
  for (int w = 0; w < warmupRepeats; ++w) {
    for (const auto& r : rays) {
      volatile double sink = shape->DistFromOutside(r.origin.data(), r.dir.data(), kIact, stepmax);
      (void)sink;
    }
  }
  uint64_t checksum = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < timedRepeats; ++r) {
    for (const auto& ray : rays) {
      checksum = detail::mixDouble(checksum, shape->DistFromOutside(ray.origin.data(), ray.dir.data(), kIact, stepmax));
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  TimingResult result;
  result.nCalls = rays.size() * static_cast<size_t>(timedRepeats);
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  result.nsPerCall = result.nCalls > 0 ? ns / static_cast<double>(result.nCalls) : 0.;
  result.checksum = checksum;
  return result;
}

TimingResult timeDistFromInside(const TGeoShape* shape, const std::vector<Ray>& rays, int warmupRepeats,
                                int timedRepeats)
{
  for (int w = 0; w < warmupRepeats; ++w) {
    for (const auto& r : rays) {
      volatile double sink = shape->DistFromInside(r.origin.data(), r.dir.data(), kIact, TGeoShape::Big());
      (void)sink;
    }
  }
  uint64_t checksum = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < timedRepeats; ++r) {
    for (const auto& ray : rays) {
      checksum = detail::mixDouble(checksum, shape->DistFromInside(ray.origin.data(), ray.dir.data(), kIact, TGeoShape::Big()));
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  TimingResult result;
  result.nCalls = rays.size() * static_cast<size_t>(timedRepeats);
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  result.nsPerCall = result.nCalls > 0 ? ns / static_cast<double>(result.nCalls) : 0.;
  result.checksum = checksum;
  return result;
}

TimingResult timeSafety(const TGeoShape* shape, const std::vector<Point3D>& points, int warmupRepeats,
                        int timedRepeats)
{
  for (int w = 0; w < warmupRepeats; ++w) {
    for (const auto& p : points) {
      volatile double sink = shape->Safety(p.data(), shape->Contains(p.data()));
      (void)sink;
    }
  }
  uint64_t checksum = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < timedRepeats; ++r) {
    for (const auto& p : points) {
      checksum = detail::mixDouble(checksum, shape->Safety(p.data(), shape->Contains(p.data())));
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  TimingResult result;
  result.nCalls = points.size() * static_cast<size_t>(timedRepeats);
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  result.nsPerCall = result.nCalls > 0 ? ns / static_cast<double>(result.nCalls) : 0.;
  result.checksum = checksum;
  return result;
}

// ------------------------------------------------------------------------------------------
// The `shape_<part>.root` sidecar
// ------------------------------------------------------------------------------------------

namespace
{
/// The key an emitter is required to write. Kept here rather than duplicated at both call sites
/// so reader and writer cannot disagree about it.
constexpr const char* kShapeKeyName = "shape";
/// The optional companion key: the shape's rigid placement, `local -> part`. Absent means
/// identity (scripts/geometry/Stream_N_PlacedPrimitives.md).
constexpr const char* kPlacementKeyName = "placement";
} // namespace

TGeoShape* loadShapeFromRootFile(const std::string& path, std::string* error)
{
  const auto fail = [error](const std::string& why) -> TGeoShape* {
    if (error != nullptr) {
      *error = why;
    }
    return nullptr;
  };
  // TFile is chatty and *fatal-looking* on a missing file; check first so a part that simply has
  // no shape sidecar (the overwhelmingly common case today) is silent rather than alarming.
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    return fail(path + ": cannot be opened as a ROOT file");
  }
  TObject* object = file->Get<TObject>(kShapeKeyName);
  if (object == nullptr) {
    // Fall back to the first TGeoShape-derived key. Emitters must write "shape"; this exists so a
    // file produced by hand (`root -e '...'`) with a different key name is still usable, which is
    // exactly how the first fixtures for this path were made.
    TIter next(file->GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
      TClass* cl = TClass::GetClass(key->GetClassName());
      if (cl != nullptr && cl->InheritsFrom(TGeoShape::Class())) {
        object = key->ReadObj();
        break;
      }
    }
  }
  if (object == nullptr) {
    return fail(path + ": holds no object inheriting from TGeoShape (expected key \"" +
                kShapeKeyName + "\")");
  }
  auto* shape = dynamic_cast<TGeoShape*>(object);
  if (shape == nullptr) {
    const std::string className = object->ClassName();
    delete object;
    return fail(path + ": key \"" + kShapeKeyName + "\" holds a " + className +
                ", which does not inherit from TGeoShape");
  }
  // The object was read out of a TDirectory but is not a TDirectory-owned type (TGeoShape is not
  // a histogram/tree), so we own it and it stays valid past the file's destruction.
  return shape;
}

TGeoHMatrix* loadShapePlacementFromRootFile(const std::string& path)
{
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    return nullptr;
  }
  auto* stored = file->Get<TGeoHMatrix>(kPlacementKeyName);
  if (stored == nullptr) {
    return nullptr;
  }
  // Copied out rather than detached: TGeoHMatrix registers itself with gGeoManager when one
  // exists, and handing the file's own object to a caller that may outlive the file has bitten
  // this project before. A copy is 16 doubles.
  auto* placement = new TGeoHMatrix(*stored);
  return placement;
}

bool saveShapeToRootFile(const std::string& path, const TGeoShape& shape, std::string* error)
{
  return saveShapeToRootFile(path, shape, nullptr, error);
}

bool saveShapeToRootFile(const std::string& path, const TGeoShape& shape,
                         const TGeoMatrix* placement, std::string* error)
{
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) {
    if (error != nullptr) {
      *error = path + ": cannot be opened for writing";
    }
    return false;
  }
  const int written = file->WriteTObject(&shape, kShapeKeyName);
  // An identity placement is deliberately NOT written: "no key" must keep meaning "identity", and
  // an artefact that records the identity explicitly would be a second way of saying the same
  // thing for a reader to get wrong.
  if (placement != nullptr && !placement->IsIdentity()) {
    TGeoHMatrix stored(*placement);
    stored.SetName(kPlacementKeyName);
    file->WriteTObject(&stored, kPlacementKeyName);
  }
  file->Close();
  if (written <= 0) {
    if (error != nullptr) {
      *error = path + ": WriteTObject wrote 0 bytes";
    }
    return false;
  }
  return true;
}

} // namespace harness
} // namespace base
} // namespace o2
