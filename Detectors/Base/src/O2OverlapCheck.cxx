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

/// \file O2OverlapCheck.cxx
/// \brief An overlap census that asks the shapes rather than their pictures.
///
/// `TGeoManager::CheckOverlaps` is the only overlap check most users will ever run, and its
/// algorithm is: take the display-mesh points of one solid, ask whether the other Contains() them,
/// and call the deepest such point the overlap. Both halves of that were broken on the shapes this
/// branch ships -- the point set was not on the shapes (fixed in
/// O2BVHSurfaceSolid::GetPointsOnSegments) and nothing checked that it was.
///
/// This is the same idea done so that the answer is evidence:
///  * every sampled point is verified to be on the boundary of the solid it came from, and
///    discarded if not, so a shape with a poor display mesh loses coverage instead of producing a
///    phantom;
///  * the discriminator is the *depth*, not containment, which is what separates the 11 Bagger
///    pairs that share a face exactly (legal, and the normal case) from the 3 that interpenetrate;
///  * a positive depth is a proof of positive shared volume, so there are no false positives
///    beyond the exactness of Contains()/Safety() themselves -- what sampling bounds is the false
///    *negatives*.
///
/// See scripts/geometry/Stream_V_OverlapCheck.md.

#include "DetectorsBase/O2OverlapCheck.h"

#include "TGeoShape.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>

namespace o2
{
namespace base
{

const char* OverlapVerdictName(OverlapVerdict verdict)
{
  switch (verdict) {
    case OverlapVerdict::Disjoint:
      return "disjoint";
    case OverlapVerdict::Touching:
      return "touching";
    case OverlapVerdict::Interpenetrating:
      return "INTERPENETRATING";
    case OverlapVerdict::Contained:
      return "CONTAINED";
  }
  return "unknown";
}

namespace
{

/// The master-frame axis-aligned box of a shape's local bounding box under \a matrix, inflated by
/// \a pad. Conservative for a rotation because it takes the box of the eight transformed corners.
struct MasterBox {
  double lower[3] = {0., 0., 0.};
  double upper[3] = {0., 0., 0.};
  bool valid = false;
};

MasterBox masterBox(const TGeoShape* shape, const TGeoMatrix* matrix, double pad)
{
  MasterBox box;
  const auto* boundingBox = dynamic_cast<const TGeoBBox*>(shape);
  if (boundingBox == nullptr) {
    return box;
  }
  const double* origin = boundingBox->GetOrigin();
  const double halfLengths[3] = {boundingBox->GetDX(), boundingBox->GetDY(), boundingBox->GetDZ()};
  for (int dimension = 0; dimension < 3; ++dimension) {
    box.lower[dimension] = std::numeric_limits<double>::max();
    box.upper[dimension] = -std::numeric_limits<double>::max();
  }
  for (int corner = 0; corner < 8; ++corner) {
    const double local[3] = {origin[0] + ((corner & 1) ? halfLengths[0] : -halfLengths[0]),
                             origin[1] + ((corner & 2) ? halfLengths[1] : -halfLengths[1]),
                             origin[2] + ((corner & 4) ? halfLengths[2] : -halfLengths[2])};
    double master[3] = {0., 0., 0.};
    matrix->LocalToMaster(local, master);
    for (int dimension = 0; dimension < 3; ++dimension) {
      box.lower[dimension] = std::min(box.lower[dimension], master[dimension]);
      box.upper[dimension] = std::max(box.upper[dimension], master[dimension]);
    }
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    box.lower[dimension] -= pad;
    box.upper[dimension] += pad;
  }
  box.valid = true;
  return box;
}

bool boxesOverlap(const MasterBox& first, const MasterBox& second)
{
  if (!first.valid || !second.valid) {
    return true; // no box means no rejection; test the pair
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    if (first.upper[dimension] < second.lower[dimension] || second.upper[dimension] < first.lower[dimension]) {
      return false;
    }
  }
  return true;
}

/// A deterministic 64-bit hash used as the volume estimator's point source. Deterministic because
/// two runs that disagree must be a geometry difference and never a seed difference.
inline double halton(unsigned int index, unsigned int base)
{
  double result = 0.;
  double fraction = 1.;
  while (index > 0) {
    fraction /= base;
    result += fraction * (index % base);
    index /= base;
  }
  return result;
}

} // namespace

int SampleBoundaryPoints(const TGeoShape* shape, int npoints, double residualTolerance,
                         std::vector<double>& points, int& rejected, double& worstResidual,
                         bool* usedPointsOnSegments)
{
  points.clear();
  rejected = 0;
  worstResidual = 0.;
  if (usedPointsOnSegments != nullptr) {
    *usedPointsOnSegments = false;
  }
  if (shape == nullptr || npoints <= 0) {
    return 0;
  }
  auto* mutableShape = const_cast<TGeoShape*>(shape);

  int meshVertices = 0;
  int meshSegments = 0;
  int meshPolygons = 0;
  mutableShape->GetMeshNumbers(meshVertices, meshSegments, meshPolygons);

  // Exactly the choice TGeoChecker::MakeCheckOverlap makes, and for the same reason: a shape that
  // declines to generate points on request still has a display mesh, and its vertices are usually
  // the better sample anyway.
  const int capacity = std::max(npoints, meshVertices);
  std::vector<double> raw(3 * static_cast<size_t>(std::max(capacity, 1)), 0.);
  int rawCount = 0;
  if (mutableShape->GetPointsOnSegments(npoints, raw.data())) {
    rawCount = npoints;
    if (usedPointsOnSegments != nullptr) {
      *usedPointsOnSegments = true;
    }
  } else {
    if (meshVertices <= 0) {
      return 0;
    }
    mutableShape->SetPoints(raw.data());
    rawCount = meshVertices;
  }

  points.reserve(3 * static_cast<size_t>(rawCount));
  for (int index = 0; index < rawCount; ++index) {
    const double* candidate = &raw[3 * static_cast<size_t>(index)];
    // Safety() is a lower bound on the distance to the boundary, so a large value is a proof that
    // the point is *not* on it. That is the direction this filter needs.
    const double residual = mutableShape->Safety(candidate, mutableShape->Contains(candidate));
    if (!(residual <= residualTolerance)) {
      rejected++;
      continue;
    }
    worstResidual = std::max(worstResidual, residual);
    points.push_back(candidate[0]);
    points.push_back(candidate[1]);
    points.push_back(candidate[2]);
  }
  return static_cast<int>(points.size() / 3);
}

namespace
{

/// One direction of the pair test: every accepted boundary point of \a points (in \a matFrom's
/// local frame) against \a target.
struct DirectionResult {
  int contained = 0;
  int deep = 0;
  double maxDepth = 0.;
  double deepestMaster[3] = {0., 0., 0.};
  double minSeparation = std::numeric_limits<double>::max();
};

DirectionResult probeDirection(const std::vector<double>& points, const TGeoMatrix* matFrom,
                               const TGeoShape* target, const TGeoMatrix* matTo, double depthTolerance)
{
  DirectionResult result;
  auto* mutableTarget = const_cast<TGeoShape*>(target);
  const size_t count = points.size() / 3;
  for (size_t index = 0; index < count; ++index) {
    double master[3] = {0., 0., 0.};
    double local[3] = {0., 0., 0.};
    matFrom->LocalToMaster(&points[3 * index], master);
    matTo->MasterToLocal(master, local);
    if (mutableTarget->Contains(local)) {
      result.contained++;
      const double depth = mutableTarget->Safety(local, kTRUE);
      if (depth > depthTolerance) {
        result.deep++;
      }
      if (depth > result.maxDepth) {
        result.maxDepth = depth;
        std::memcpy(result.deepestMaster, master, 3 * sizeof(double));
      }
    } else {
      result.minSeparation = std::min(result.minSeparation, mutableTarget->Safety(local, kFALSE));
    }
  }
  return result;
}

} // namespace

OverlapPair CheckPairOverlap(const TGeoShape* shapeA, const TGeoMatrix* matA, const std::string& nameA,
                             const TGeoShape* shapeB, const TGeoMatrix* matB, const std::string& nameB,
                             const OverlapOptions& options)
{
  OverlapPair pair;
  pair.nameA = nameA;
  pair.nameB = nameB;
  if (shapeA == nullptr || shapeB == nullptr || matA == nullptr || matB == nullptr) {
    return pair;
  }

  int rejectedA = 0;
  int rejectedB = 0;
  double residualA = 0.;
  double residualB = 0.;
  std::vector<double> pointsA;
  std::vector<double> pointsB;
  pair.sampledA = SampleBoundaryPoints(shapeA, options.pointsPerSolid, options.residualTolerance, pointsA, rejectedA,
                                       residualA);
  pair.sampledB = SampleBoundaryPoints(shapeB, options.pointsPerSolid, options.residualTolerance, pointsB, rejectedB,
                                       residualB);

  const DirectionResult aInB = probeDirection(pointsA, matA, shapeB, matB, options.depthTolerance);
  const DirectionResult bInA = probeDirection(pointsB, matB, shapeA, matA, options.depthTolerance);

  pair.pointsAInsideB = aInB.contained;
  pair.pointsBInsideA = bInA.contained;
  pair.deepPointsAInsideB = aInB.deep;
  pair.deepPointsBInsideA = bInA.deep;

  if (aInB.maxDepth >= bInA.maxDepth) {
    pair.depthCm = aInB.maxDepth;
    std::copy(aInB.deepestMaster, aInB.deepestMaster + 3, pair.deepestPoint.begin());
    pair.deepestPointFrom = nameA;
  } else {
    pair.depthCm = bInA.maxDepth;
    std::copy(bInA.deepestMaster, bInA.deepestMaster + 3, pair.deepestPoint.begin());
    pair.deepestPointFrom = nameB;
  }

  const bool anyContained = (aInB.contained > 0) || (bInA.contained > 0);
  const bool anyDeep = (aInB.deep > 0) || (bInA.deep > 0);
  // Containment: every boundary point of one solid is inside the other, and none of them is merely
  // on its boundary. Legal only as a declared mother/daughter, which a flat conversion never emits.
  const bool allAInside = pair.sampledA > 0 && aInB.contained == pair.sampledA && aInB.deep == pair.sampledA;
  const bool allBInside = pair.sampledB > 0 && bInA.contained == pair.sampledB && bInA.deep == pair.sampledB;

  if (allAInside || allBInside) {
    pair.verdict = OverlapVerdict::Contained;
  } else if (anyDeep) {
    pair.verdict = OverlapVerdict::Interpenetrating;
  } else if (anyContained) {
    pair.verdict = OverlapVerdict::Touching;
  } else {
    pair.verdict = OverlapVerdict::Disjoint;
    const double separation = std::min(aInB.minSeparation, bInA.minSeparation);
    if (separation < std::numeric_limits<double>::max()) {
      pair.separationCm = separation;
    }
  }

  if (options.volumeSamples > 0 &&
      (pair.verdict == OverlapVerdict::Interpenetrating || pair.verdict == OverlapVerdict::Contained)) {
    const MasterBox boxA = masterBox(shapeA, matA, 0.);
    const MasterBox boxB = masterBox(shapeB, matB, 0.);
    if (boxA.valid && boxB.valid) {
      double lower[3];
      double upper[3];
      double boxVolume = 1.;
      for (int dimension = 0; dimension < 3; ++dimension) {
        lower[dimension] = std::max(boxA.lower[dimension], boxB.lower[dimension]);
        upper[dimension] = std::min(boxA.upper[dimension], boxB.upper[dimension]);
        boxVolume *= std::max(0., upper[dimension] - lower[dimension]);
      }
      if (boxVolume > 0.) {
        auto* mutableA = const_cast<TGeoShape*>(shapeA);
        auto* mutableB = const_cast<TGeoShape*>(shapeB);
        int hits = 0;
        for (int sample = 0; sample < options.volumeSamples; ++sample) {
          const double master[3] = {lower[0] + (upper[0] - lower[0]) * halton(sample + 1, 2),
                                    lower[1] + (upper[1] - lower[1]) * halton(sample + 1, 3),
                                    lower[2] + (upper[2] - lower[2]) * halton(sample + 1, 5)};
          double local[3];
          matA->MasterToLocal(master, local);
          if (!mutableA->Contains(local)) {
            continue;
          }
          matB->MasterToLocal(master, local);
          if (mutableB->Contains(local)) {
            hits++;
          }
        }
        const double fraction = double(hits) / options.volumeSamples;
        pair.sharedVolumeHits = hits;
        pair.sharedVolumeCm3 = fraction * boxVolume;
        pair.sharedVolumeErrCm3 = std::sqrt(std::max(1., double(hits))) / options.volumeSamples * boxVolume;
      }
    }
  }

  return pair;
}

OverlapCensus CheckWorldOverlaps(const TGeoVolume* volume, const OverlapOptions& options)
{
  const auto startTime = std::chrono::steady_clock::now();
  OverlapCensus census;
  if (volume == nullptr) {
    return census;
  }
  auto* mutableVolume = const_cast<TGeoVolume*>(volume);
  const int daughters = mutableVolume->GetNdaughters();
  census.nSolids = daughters;
  census.nPairsTotal = daughters * (daughters - 1) / 2;

  std::vector<const TGeoShape*> shapes(daughters, nullptr);
  std::vector<const TGeoMatrix*> matrices(daughters, nullptr);
  std::vector<std::string> names(daughters);
  std::vector<MasterBox> boxes(daughters);
  std::vector<std::vector<double>> points(daughters);

  for (int index = 0; index < daughters; ++index) {
    TGeoNode* node = mutableVolume->GetNode(index);
    shapes[index] = node->GetVolume()->GetShape();
    matrices[index] = node->GetMatrix();
    names[index] = node->GetVolume()->GetName();
    boxes[index] = masterBox(shapes[index], matrices[index], options.padCm);

    OverlapSolidReport report;
    report.name = names[index];
    report.shapeClass = shapes[index] != nullptr ? shapes[index]->ClassName() : "none";
    report.requested = options.pointsPerSolid;
    bool usedSegments = false;
    report.accepted = SampleBoundaryPoints(shapes[index], options.pointsPerSolid, options.residualTolerance,
                                           points[index], report.rejected, report.worstResidualCm, &usedSegments);
    report.usedPointsOnSegments = usedSegments;
    census.nPointsRejected += report.rejected;
    census.worstResidualCm = std::max(census.worstResidualCm, report.worstResidualCm);
    census.solids.push_back(report);
  }

  for (int first = 0; first < daughters; ++first) {
    for (int second = first + 1; second < daughters; ++second) {
      if (!boxesOverlap(boxes[first], boxes[second])) {
        continue;
      }
      census.nPairsTested++;
      // Reuse the point sets: sampling is the expensive part and it does not depend on the partner.
      OverlapPair pair;
      pair.nameA = names[first];
      pair.nameB = names[second];
      pair.sampledA = static_cast<int>(points[first].size() / 3);
      pair.sampledB = static_cast<int>(points[second].size() / 3);
      const DirectionResult aInB =
        probeDirection(points[first], matrices[first], shapes[second], matrices[second], options.depthTolerance);
      const DirectionResult bInA =
        probeDirection(points[second], matrices[second], shapes[first], matrices[first], options.depthTolerance);
      pair.pointsAInsideB = aInB.contained;
      pair.pointsBInsideA = bInA.contained;
      pair.deepPointsAInsideB = aInB.deep;
      pair.deepPointsBInsideA = bInA.deep;
      if (aInB.maxDepth >= bInA.maxDepth) {
        pair.depthCm = aInB.maxDepth;
        std::copy(aInB.deepestMaster, aInB.deepestMaster + 3, pair.deepestPoint.begin());
        pair.deepestPointFrom = names[first];
      } else {
        pair.depthCm = bInA.maxDepth;
        std::copy(bInA.deepestMaster, bInA.deepestMaster + 3, pair.deepestPoint.begin());
        pair.deepestPointFrom = names[second];
      }
      const bool allAInside = pair.sampledA > 0 && aInB.contained == pair.sampledA && aInB.deep == pair.sampledA;
      const bool allBInside = pair.sampledB > 0 && bInA.contained == pair.sampledB && bInA.deep == pair.sampledB;
      if (allAInside || allBInside) {
        pair.verdict = OverlapVerdict::Contained;
        census.nContained++;
      } else if (aInB.deep > 0 || bInA.deep > 0) {
        pair.verdict = OverlapVerdict::Interpenetrating;
        census.nInterpenetrating++;
      } else if (aInB.contained > 0 || bInA.contained > 0) {
        pair.verdict = OverlapVerdict::Touching;
        census.nTouching++;
      } else {
        pair.verdict = OverlapVerdict::Disjoint;
        const double separation = std::min(aInB.minSeparation, bInA.minSeparation);
        if (separation < std::numeric_limits<double>::max()) {
          pair.separationCm = separation;
        }
        census.nDisjoint++;
      }

      if (options.volumeSamples > 0 && (pair.verdict == OverlapVerdict::Interpenetrating ||
                                        pair.verdict == OverlapVerdict::Contained)) {
        const OverlapPair measured = CheckPairOverlap(shapes[first], matrices[first], names[first], shapes[second],
                                                      matrices[second], names[second], options);
        pair.sharedVolumeCm3 = measured.sharedVolumeCm3;
        pair.sharedVolumeErrCm3 = measured.sharedVolumeErrCm3;
        pair.sharedVolumeHits = measured.sharedVolumeHits;
      }
      census.pairs.push_back(pair);
    }
  }

  // Extrusion: a daughter's boundary point outside its mother. The mother's own boundary points are
  // irrelevant here -- they are outside every daughter by construction on a legal geometry, and
  // inside one on an illegal one, which is the sibling test's job.
  if (options.checkExtrusion && mutableVolume->GetShape() != nullptr && !mutableVolume->IsAssembly()) {
    TGeoIdentity identity;
    for (int index = 0; index < daughters; ++index) {
      OverlapPair pair;
      pair.nameA = names[index];
      pair.nameB = mutableVolume->GetName();
      pair.sampledA = static_cast<int>(points[index].size() / 3);
      auto* mother = const_cast<TGeoShape*>(mutableVolume->GetShape());
      double worst = 0.;
      int outside = 0;
      double worstMaster[3] = {0., 0., 0.};
      for (size_t point = 0; point < points[index].size() / 3; ++point) {
        double master[3] = {0., 0., 0.};
        matrices[index]->LocalToMaster(&points[index][3 * point], master);
        if (!mother->Contains(master)) {
          const double depth = mother->Safety(master, kFALSE);
          if (depth > options.depthTolerance) {
            outside++;
            if (depth > worst) {
              worst = depth;
              std::memcpy(worstMaster, master, 3 * sizeof(double));
            }
          }
        }
      }
      if (outside > 0) {
        pair.verdict = OverlapVerdict::Interpenetrating;
        pair.depthCm = worst;
        pair.deepPointsAInsideB = outside;
        pair.deepestPointFrom = names[index];
        std::copy(worstMaster, worstMaster + 3, pair.deepestPoint.begin());
        census.extrusions.push_back(pair);
        census.nExtruding++;
      }
    }
  }

  census.elapsedSeconds =
    std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
  return census;
}

} // namespace base
} // namespace o2
