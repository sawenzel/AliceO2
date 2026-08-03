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

#ifndef ALICEO2_BASE_O2OVERLAPCHECK_
#define ALICEO2_BASE_O2OVERLAPCHECK_

#include <array>
#include <string>
#include <vector>

class TGeoShape;
class TGeoMatrix;
class TGeoVolume;

namespace o2
{
namespace base
{

/// Whether two placed solids may legally coexist in a TGeo/Geant4 world.
///
/// The three states are not degrees of one quantity, they are different situations, and an
/// assembly's normal state is the middle one: 11 of `Bagger.step`'s 78 pairs share a face exactly
/// and are perfectly legal, while 3 interpenetrate and are not. A checker that cannot tell those
/// apart would condemn every real assembly ever authored.
enum class OverlapVerdict {
  Disjoint,         ///< no sampled boundary point of either solid lies inside the other
  Touching,         ///< boundary points coincide, but none is deeper than the depth tolerance
  Interpenetrating, ///< a boundary point of one solid lies strictly inside the other: illegal
  Contained         ///< every sampled boundary point of the smaller solid is inside the other
};

const char* OverlapVerdictName(OverlapVerdict verdict);

struct OverlapOptions {
  /// Boundary points sampled per solid. Coverage, not accuracy: every individual answer is exact,
  /// so this bounds the *false negatives* and nothing else.
  int pointsPerSolid = 20000;
  /// A containment shallower than this is a shared boundary, not an overlap. In cm.
  double depthTolerance = 1.e-6;
  /// A sampled point further than this from the boundary of the solid it was sampled from is not
  /// evidence about anything and is discarded (and counted). In cm.
  double residualTolerance = 1.e-6;
  /// Bounding boxes are inflated by this before the pairwise rejection. It does not change which
  /// pairs are found to overlap -- an overlapping pair's boxes overlap at pad 0 -- it decides which
  /// *disjoint* pairs get their separation measured. In cm.
  double padCm = 0.1;
  /// Monte-Carlo samples in the bounding-box intersection used to estimate the shared volume of an
  /// interpenetrating pair. 0 disables it, which is the default: the depth already proves the
  /// shared volume is positive, and the estimate costs 2 x N Contains() calls per pair.
  int volumeSamples = 0;
  /// Also test every daughter against the mother it sits in (ROOT's "extrusion" case). Silently a
  /// no-op when the mother is an assembly, which has no shape to be extruded from.
  bool checkExtrusion = true;
};

/// One pair of placed solids, and everything measured about it.
struct OverlapPair {
  std::string nameA;
  std::string nameB;
  OverlapVerdict verdict = OverlapVerdict::Disjoint;

  /// The largest distance, over every sampled boundary point of one solid that lies inside the
  /// other, from that point to the other solid's boundary. **This is the verdict's evidence and
  /// the quantity that separates touching from interpenetrating.**
  ///
  /// It is a *lower bound* on the penetration depth (the sample may miss the deepest point) and it
  /// is a proof of illegality when positive: a point of dA lying at distance d > 0 from dB and
  /// inside B means a ball of radius d around it is inside B, and that ball meets the interior of
  /// A, so the two interiors share positive volume.
  double depthCm = 0.;
  std::array<double, 3> deepestPoint{{0., 0., 0.}}; ///< in the master frame
  std::string deepestPointFrom;                     ///< which solid's boundary the deepest point came from

  int pointsAInsideB = 0; ///< sampled points of A found inside B at any depth
  int pointsBInsideA = 0;
  int deepPointsAInsideB = 0; ///< ... of which deeper than depthTolerance
  int deepPointsBInsideA = 0;
  int sampledA = 0; ///< accepted (on-boundary) sample counts actually used
  int sampledB = 0;

  /// Smallest sampled distance from a boundary point of one solid to the other solid, in cm. Only
  /// meaningful when the verdict is Disjoint, and it is a sampled estimate rather than a bound:
  /// the sample may miss the closest point (pushing it up) and Safety() is a lower bound on the
  /// true distance (pushing it down).
  double separationCm = -1.;

  double sharedVolumeCm3 = -1.;   ///< Monte-Carlo estimate; < 0 when not measured
  double sharedVolumeErrCm3 = 0.; ///< its 1-sigma statistical error
  int sharedVolumeHits = 0;
};

/// One solid's own sampling report. A point that its own shape does not consider to be on its
/// boundary is discarded before any pair sees it, so this is where a shape with a bad display mesh
/// -- ROOT's composites, above all -- shows up as reduced coverage instead of as a wrong verdict.
struct OverlapSolidReport {
  std::string name;
  std::string shapeClass;
  int requested = 0;
  int accepted = 0;
  int rejected = 0;
  double worstResidualCm = 0.; ///< the largest own-boundary distance among the *accepted* points
  bool usedPointsOnSegments = false;
};

struct OverlapCensus {
  std::vector<OverlapSolidReport> solids;
  std::vector<OverlapPair> pairs; ///< only the pairs that survived the bounding-box rejection
  std::vector<OverlapPair> extrusions;

  int nSolids = 0;
  int nPairsTotal = 0;  ///< N (N - 1) / 2
  int nPairsTested = 0; ///< after the bounding-box rejection
  int nDisjoint = 0;
  int nTouching = 0;
  int nInterpenetrating = 0;
  int nContained = 0;
  int nExtruding = 0;
  int nPointsRejected = 0;
  double worstResidualCm = 0.;
  double elapsedSeconds = 0.;

  /// The one-line answer: nInterpenetrating + nContained + nExtruding.
  int illegalCount() const { return nInterpenetrating + nContained + nExtruding; }
};

/// Sample \a npoints points on \a shape's own boundary into \a points (3 doubles each), keeping
/// only those the shape's own Safety() puts within \a residualTolerance of its boundary. Returns
/// the number kept; \a rejected and \a worstResidual report what the filter saw.
///
/// This is the one place the whole check depends on the display mesh, and it depends on it only
/// for *where* to look. Any shape whose mesh is not on its surface loses coverage here rather than
/// producing a wrong answer downstream.
int SampleBoundaryPoints(const TGeoShape* shape, int npoints, double residualTolerance,
                         std::vector<double>& points, int& rejected, double& worstResidual,
                         bool* usedPointsOnSegments = nullptr);

/// Test one placed pair. \a matA / \a matB take each shape's local frame to the common frame.
OverlapPair CheckPairOverlap(const TGeoShape* shapeA, const TGeoMatrix* matA, const std::string& nameA,
                             const TGeoShape* shapeB, const TGeoMatrix* matB, const std::string& nameB,
                             const OverlapOptions& options = OverlapOptions());

/// Census every pair of \a volume's immediate daughters, plus (optionally) each daughter against
/// \a volume itself. Nested levels are out of scope by design: this answers the question a flat
/// CAD -> TGeo conversion asks, which is whether the parts it placed side by side are legal.
OverlapCensus CheckWorldOverlaps(const TGeoVolume* volume, const OverlapOptions& options = OverlapOptions());

} // namespace base
} // namespace o2

#endif
