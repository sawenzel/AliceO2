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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
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

void recordOffender(ValidationResult& result, const ValidationOptions& opt, Offender&& off, bool withinBand)
{
  if (withinBand) {
    ++result.nMismatchWithinBand;
  } else {
    ++result.nMismatchUnexplained;
  }
  result.worstDeviation = std::max(result.worstDeviation, std::fabs(off.deviation));
  result.worstOffenders.push_back(std::move(off));
  std::sort(result.worstOffenders.begin(), result.worstOffenders.end(),
           [](const Offender& a, const Offender& b) { return std::fabs(a.deviation) > std::fabs(b.deviation); });
  if (result.worstOffenders.size() > opt.maxOffenders) {
    result.worstOffenders.resize(opt.maxOffenders);
  }
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
    recordOffender(result, opt, std::move(off), refSafety < opt.meshBand);
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
    const double probeDist = dcBig ? dr : dc;
    const Point3D probePoint = add(r.origin, scale(r.dir, probeDist));
    const double refSafety = reference->Safety(probePoint.data(), reference->Contains(probePoint.data()));
    Offender off;
    off.point = r.origin;
    off.dir = r.dir;
    off.candidateValue = dcBig ? opt.stepmax : dc;
    off.referenceValue = drBig ? opt.stepmax : dr;
    off.deviation = off.candidateValue - off.referenceValue;
    off.referenceSafety = refSafety;
    recordOffender(result, opt, std::move(off), refSafety < opt.meshBand);
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
    const double probeDist = dcBig ? dr : dc;
    const Point3D probePoint = add(r.origin, scale(r.dir, probeDist));
    const double refSafety = reference->Safety(probePoint.data(), reference->Contains(probePoint.data()));
    Offender off;
    off.point = r.origin;
    off.dir = r.dir;
    off.candidateValue = dcBig ? opt.stepmax : dc;
    off.referenceValue = drBig ? opt.stepmax : dr;
    off.deviation = off.candidateValue - off.referenceValue;
    off.referenceSafety = refSafety;
    recordOffender(result, opt, std::move(off), refSafety < opt.meshBand);
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
    recordOffender(result, opt, std::move(off), /*withinBand=*/false);
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

} // namespace harness
} // namespace base
} // namespace o2
