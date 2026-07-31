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

#ifndef ALICEO2_BASE_O2BVHSURFACESOLID_
#define ALICEO2_BASE_O2BVHSURFACESOLID_

#include "TGeoBBox.h"

#include <array>
#include <iosfwd>
#include <vector>

class TBuffer3D;

namespace o2
{
namespace base
{

/// One boundary curve of a BVHSurfaceRecord, in the flat form ROOT streams: a straight segment,
/// a circular arc or a (rational) B-spline. This is the persistent mirror of
/// O2BVHSurfaceSolid::PlanarBoundaryCurve and is deliberately a separate type, so that the
/// on-disk schema is a decision of its own rather than an accident of the in-memory API.
struct BVHSurfaceCurveRecord {
  int kind = 0;                   ///< PlanarBoundaryCurve::Kind: 0 = Line, 1 = Arc, 2 = BSpline
  double lineStart[2] = {0., 0.};
  double lineEnd[2] = {0., 0.};
  double center[2] = {0., 0.};
  double radius = 0.;
  double startAngle = 0.;
  double endAngle = 0.;
  int degree = 0;                 ///< B-spline degree
  std::vector<double> poles;      ///< B-spline control points, flattened (u, v) pairs
  std::vector<double> weights;    ///< B-spline weights (empty => non-rational)
  std::vector<double> knots;      ///< B-spline clamped flat knot vector
};

/// The persistent record of one successful O2BVHSurfaceSolid::Add*Surface call.
///
/// The kernel objects a solid is made of (the bounded surfaces, the BVH, the display mesh) are
/// all *derived* state, cheaply reconstructed by replaying the calls that produced them. So the
/// records are what the solid persists, and reading one back means replaying them and calling
/// CloseShape() -- which also means the closure diagnostics of a read-back solid are recomputed
/// rather than trusted.
struct BVHSurfaceRecord {
  enum Kind { PlanarPolygon = 0, CurvedPlanar = 1, Cylindrical = 2, Spherical = 3, Conical = 4, Toroidal = 5 };

  int kind = PlanarPolygon;
  double origin[3] = {0., 0., 0.}; ///< origin / centerPoint / center
  double axisA[3] = {0., 0., 0.};  ///< axisU / axis / polarAxis
  double axisB[3] = {0., 0., 0.};  ///< axisV / referenceAxisU

  /// The remaining scalar arguments, in the declaration order of the Add*Surface signature:
  ///   PlanarPolygon, CurvedPlanar : (none)
  ///   Cylindrical : radius, heightMin, heightMax, phiStart, phiSweep
  ///   Spherical   : radius, thetaMin, thetaMax, phiStart, phiSweep
  ///   Conical     : radiusAtMin, radiusAtMax, heightMin, heightMax, phiStart, phiSweep
  ///   Toroidal    : majorRadius, minorRadius, phiStart, phiSweep, tubeStart, tubeSweep
  /// The count is checked against the kind on replay, so a truncated or foreign record is
  /// rejected instead of being silently reinterpreted.
  std::vector<double> scalars;

  bool innerWall = false;
  bool trimmed = false;                      ///< the wire-trim overload was used (quadrics only)

  /// The wires, outer wire first and then the inner wires (holes). PlanarPolygon surfaces store
  /// their vertices in \a polygonPoints as flattened (u, v) pairs; every other kind stores its
  /// boundary curves in \a curves. \a wireSizes says how many points/curves each wire owns.
  std::vector<double> polygonPoints;
  std::vector<BVHSurfaceCurveRecord> curves;
  std::vector<int> wireSizes;

  /// How many entries \a scalars must hold for \a kind, or -1 for an unknown kind.
  static int expectedScalarCount(int recordKind);
};

class O2BVHSurfaceSolid : public TGeoBBox
{
 public:
  using Point2D = std::array<double, 2>;
  using Point3D = std::array<double, 3>;

  O2BVHSurfaceSolid();
  explicit O2BVHSurfaceSolid(const char* name);
  ~O2BVHSurfaceSolid() override;

  O2BVHSurfaceSolid(const O2BVHSurfaceSolid&) = delete;
  O2BVHSurfaceSolid& operator=(const O2BVHSurfaceSolid&) = delete;

  bool AddPlanarSurface(const Point3D& origin, const Point3D& axisU, const Point3D& axisV,
                        const std::vector<Point2D>& outerWire,
                        const std::vector<std::vector<Point2D>>& innerWires = {});

  /// One boundary curve of a curved planar surface, in the surface's local (u, v) frame. It is a
  /// straight line segment (\a lineStart -> \a lineEnd), a circular arc (\a center, \a radius,
  /// from \a startAngle to \a endAngle in radians; a full circle is a +/-2pi sweep) or a clamped
  /// (rational) B-spline (\a degree, \a poles, optional \a weights, clamped flat \a knots). This
  /// is the public, kernel-independent mirror of the internal Curve2D type.
  struct PlanarBoundaryCurve {
    enum Kind { Line, Arc, BSpline };
    Kind kind = Line;
    Point2D lineStart{0., 0.};
    Point2D lineEnd{0., 0.};
    Point2D center{0., 0.};
    double radius = 0.;
    double startAngle = 0.;
    double endAngle = 0.;
    int degree = 0;                    ///< B-spline degree
    std::vector<Point2D> poles;        ///< B-spline control points
    std::vector<double> weights;       ///< B-spline weights (empty ⇒ non-rational)
    std::vector<double> knots;         ///< B-spline clamped flat knot vector

    static PlanarBoundaryCurve makeLine(const Point2D& start, const Point2D& end)
    {
      PlanarBoundaryCurve curve;
      curve.kind = Line;
      curve.lineStart = start;
      curve.lineEnd = end;
      return curve;
    }
    static PlanarBoundaryCurve makeArc(const Point2D& c, double r, double start, double end)
    {
      PlanarBoundaryCurve curve;
      curve.kind = Arc;
      curve.center = c;
      curve.radius = r;
      curve.startAngle = start;
      curve.endAngle = end;
      return curve;
    }
    static PlanarBoundaryCurve makeBSpline(int splineDegree, std::vector<Point2D> splinePoles,
                                           std::vector<double> splineWeights, std::vector<double> splineKnots)
    {
      PlanarBoundaryCurve curve;
      curve.kind = BSpline;
      curve.degree = splineDegree;
      curve.poles = std::move(splinePoles);
      curve.weights = std::move(splineWeights);
      curve.knots = std::move(splineKnots);
      return curve;
    }
  };

  /// Add an exact planar surface bounded by general line/arc wires (e.g. a rounded rectangle,
  /// a slot, a disk or an annulus). \a axisU and \a axisV must be orthonormal; the outward
  /// normal is axisU x axisV. Wires are auto-reoriented (outer CCW, inner CW) as needed.
  bool AddCurvedPlanarSurface(const Point3D& origin, const Point3D& axisU, const Point3D& axisV,
                              const std::vector<PlanarBoundaryCurve>& outerWire,
                              const std::vector<std::vector<PlanarBoundaryCurve>>& innerWires = {});

  /// Add a cylindrical lateral surface of given \a radius around \a axis. \a centerPoint is the
  /// height reference (h = 0) on the axis, the height range runs along the (normalized) axis and
  /// phi is measured from \a referenceAxisU (projected perpendicular to the axis). With
  /// \a innerWall set, the surface bounds a hole and its outward normal points towards the axis.
  bool AddCylindricalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                             double radius, double heightMin, double heightMax, double phiStart = 0.,
                             double phiSweep = 6.283185307179586, bool innerWall = false);

  /// Wire-trimmed cylindrical surface: as AddCylindricalSurface, but the domain is the general
  /// line/arc trim \a outerTrim (plus optional \a innerTrims holes) in the periodic parametric
  /// (phi[rad], h[cm]) domain instead of the scalar phi/height rectangle. The scalar arguments
  /// still pin the frame/radius and a nominal window; the wire is authoritative for containment.
  bool AddCylindricalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                             double radius, double heightMin, double heightMax, double phiStart, double phiSweep,
                             bool innerWall, const std::vector<PlanarBoundaryCurve>& outerTrim,
                             const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims = {});

  /// Add a spherical surface of given \a radius around \a center, trimmed to a polar range
  /// (theta from the +polarAxis pole) and a phi sweep from \a referenceAxisU. The default
  /// arguments give a full, self-closing sphere.
  bool AddSphericalSurface(const Point3D& center, const Point3D& polarAxis, const Point3D& referenceAxisU,
                           double radius, double thetaMin = 0., double thetaMax = 3.141592653589793,
                           double phiStart = 0., double phiSweep = 6.283185307179586, bool innerWall = false);

  /// Wire-trimmed spherical surface: as AddSphericalSurface, but the domain is the general
  /// line/arc trim \a outerTrim (plus optional \a innerTrims holes) in the parametric
  /// (phi[rad], theta[rad]) domain instead of the scalar theta/phi rectangle.
  bool AddSphericalSurface(const Point3D& center, const Point3D& polarAxis, const Point3D& referenceAxisU,
                           double radius, double thetaMin, double thetaMax, double phiStart, double phiSweep,
                           bool innerWall, const std::vector<PlanarBoundaryCurve>& outerTrim,
                           const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims = {});

  /// Add a conical lateral surface whose radius runs linearly from \a radiusAtMin at
  /// \a heightMin to \a radiusAtMax at \a heightMax along \a axis; one radius may be zero (apex
  /// cone). Frame and innerWall conventions as for AddCylindricalSurface.
  bool AddConicalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                         double radiusAtMin, double radiusAtMax, double heightMin, double heightMax,
                         double phiStart = 0., double phiSweep = 6.283185307179586, bool innerWall = false);

  /// Wire-trimmed conical surface: as AddConicalSurface, but the domain is the general line/arc
  /// trim \a outerTrim (plus optional \a innerTrims holes) in the periodic parametric
  /// (phi[rad], h[cm]) domain instead of the scalar phi/height rectangle. The scalar radii pin
  /// the linear radius law r(h); the wire is authoritative for containment.
  bool AddConicalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                         double radiusAtMin, double radiusAtMax, double heightMin, double heightMax, double phiStart,
                         double phiSweep, bool innerWall, const std::vector<PlanarBoundaryCurve>& outerTrim,
                         const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims = {});

  /// Add a toroidal surface of major radius \a majorRadius (axis to tube centre) and minor
  /// (tube) radius \a minorRadius about \a centerPoint / \a axis. \a referenceAxisU fixes the
  /// phiRing = 0 direction (projected perpendicular to the axis). The patch is trimmed to a
  /// parametric rectangle in the two periodic angles: phiRing (around the axis) from
  /// \a phiStart over \a phiSweep and phiTube (around the tube, measured from the outer equator
  /// towards the +axis pole) from \a tubeStart over \a tubeSweep. The default arguments give a
  /// full, self-closing torus. With \a innerWall set the outward normal points towards the tube
  /// spine (e.g. the inner tube of a toroidal shell).
  bool AddToroidalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                          double majorRadius, double minorRadius, double phiStart = 0.,
                          double phiSweep = 6.283185307179586, double tubeStart = 0.,
                          double tubeSweep = 6.283185307179586, bool innerWall = false);

  /// Wire-trimmed toroidal surface: as AddToroidalSurface, but the domain is the general
  /// line/arc/bspline trim \a outerTrim (plus optional \a innerTrims holes) in the periodic
  /// parametric (phiRing[rad], phiTube[rad]) domain instead of the scalar phiRing/phiTube
  /// rectangle. The scalar arguments still pin the frame/radii and a nominal window; the wire is
  /// authoritative for containment. The trim must not wrap more than a full turn in either angle.
  bool AddToroidalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                          double majorRadius, double minorRadius, double phiStart, double phiSweep, double tubeStart,
                          double tubeSweep, bool innerWall, const std::vector<PlanarBoundaryCurve>& outerTrim,
                          const std::vector<std::vector<PlanarBoundaryCurve>>& innerTrims = {});

  /// Finalize the shape: compute the bounding box, the display mesh, the closure/orientation
  /// diagnostics (when \a check is set) and build the BVH acceleration structure over the
  /// bounded-surface AABBs used by the navigation queries.
  void CloseShape(bool check = true);

  int GetNsurfaces() const;
  bool IsDefined() const;

  /// \name The model's own tolerance
  /// How well the CAD model that produced this solid says its boundary is defined, in cm. It is
  /// the source model's declared tolerance (the largest BRep sub-shape tolerance), carried through
  /// the sidecar rather than guessed: without it the kernel has no way to know what epsilon two
  /// faces of the same imported solid should be expected to agree to, and every closure or
  /// adjacency decision has to fall back on a constant that nobody chose for this geometry.
  ///
  /// Zero means "not stated": a solid built through the Add*Surface API directly, or read from a
  /// sidecar older than version 2. Consumers must treat zero as unknown and fall back on their own
  /// documented constant rather than on a tolerance of nothing.
  /// @{
  void SetModelTolerance(double toleranceCm);
  double GetModelTolerance() const { return fModelTolerance; }
  /// @}

  /// Whether the BVH acceleration structure has been built (after CloseShape).
  bool HasBVH() const;
  /// Fill the BVH root-node bounding box; returns false when no BVH has been built.
  bool GetBVHRootBounds(Point3D& lower, Point3D& upper) const;
  /// Diagnostic/test hook: number of candidate surface patches whose BVH leaf boxes the given
  /// ray traverses (counted with multiplicity per leaf primitive). Returns -1 without a BVH.
  int CountBVHRayCandidates(const Point3D& point, const Point3D& direction) const;

  /// Ray tmax tightening in the BVH-accelerated distance queries: as a hit is found the
  /// traversal bound shrinks to it, so nodes lying entirely beyond the current best hit are
  /// never visited. Enabled by default. Switching it off must not change any returned distance;
  /// the switch exists so a benchmark can price the optimization by running the same rays both
  /// ways. Process-wide and not thread safe: flip it between measurement passes, never during
  /// one.
  static void SetRayTMaxPruning(bool enable);
  static bool GetRayTMaxPruning();

  /// Per-thread counter of surface patches handed to the BVH leaf callback by DistFromOutside
  /// and DistFromInside since the last reset. Together with SetRayTMaxPruning this measures how
  /// much work the pruning actually avoids. Diagnostic only; not incremented by the _Loop
  /// variants, which visit every surface by construction.
  static void ResetRayCandidateCounter();
  static long long GetRayCandidateCount();

  /// One crossing of the containment parity ray, as seen by Contains().
  struct ContainsCrossing {
    double distance = 0.;      ///< ray parameter of the hit
    double normalAlignment = 0.; ///< dot(hit normal, test direction): < 0 enters, > 0 exits
  };

  /// Diagnostic hook: the parity ray's crossing list at \a point, from the BVH traversal into
  /// \a bvhCrossings and from the all-surfaces loop into \a loopCrossings, both sorted by
  /// distance. Contains() and Contains_Loop() are the parities of these two lists, so when they
  /// disagree this shows *why* -- a crossing present in one list and not the other, or a pair
  /// close enough that one path clustered it and the other did not.
  ///
  /// It exists because the disagreement on non-manifold parts was attributed to order-dependent
  /// clustering, which reading the code does not support: the clusterer sorts first, so it is a
  /// function of the sorted distances alone. That leaves differing hit *multisets* as the
  /// explanation, which is what this dumps. Measure before theorizing (CodeReview_Fable.md, S6).
  void DescribeContainsCrossings(const Point3D& point, std::vector<ContainsCrossing>& bvhCrossings,
                                 std::vector<ContainsCrossing>& loopCrossings) const;

  /// Whether the closed shape forms a closed 2-manifold (every boundary edge shared by two faces).
  /// Meaningful only after CloseShape(); detects e.g. missing faces.
  bool IsClosed() const;
  /// Whether all shared boundary edges are traversed in opposite directions after CloseShape();
  /// detects e.g. reversed faces (inconsistent outward normals).
  bool IsOrientationConsistent() const;

  /// How far the navigation queries can be trusted on this solid, as determined by CloseShape's
  /// half-edge closure check. Parity-based containment answers a *topological* question ("does a
  /// ray from this point cross the boundary an odd number of times"), so it is only defined on a
  /// closed, consistently oriented 2-manifold. On anything else the answers are not merely
  /// imprecise near the defect: a single sliver gap flips Contains over the gap's whole shadow
  /// along the parity test direction, which can be centimetres of wrong answers far from any
  /// surface (see scripts/geometry/ExactTrimTopology.md). Callers that care about correctness
  /// must check this; the solid still answers queries either way (see CloseShape).
  ///
  /// The values are ordered by increasing severity, so `reliability > Reliable` is "do not trust
  /// this" and the worst defect present is the one reported.
  enum class NavigationReliability {
    Undetermined = 0,       ///< CloseShape() has not run yet: no diagnostics exist
    Reliable,               ///< closed, consistently oriented 2-manifold: parity is well defined
    ReversedFaces,          ///< closed, but some shared edge is traversed the same way by both
                            ///< faces: at least one face's outward normal points inward
    OpenSurfaceSet,         ///< boundary edges (missing faces / trim gaps): parity is undefined in
                            ///< the shadow of every gap along the parity test direction
    NonManifold             ///< edges shared by more than two faces (coincident/duplicated faces):
                            ///< parity depends on the order hits are clustered in
  };

  /// The reliability state derived from the last CloseShape(); Undetermined before it has run.
  NavigationReliability GetNavigationReliability() const;
  /// Shorthand for GetNavigationReliability() == NavigationReliability::Reliable. False means the
  /// navigation answers of this solid are not to be trusted anywhere, not just near the defect.
  bool IsNavigable() const;
  /// Short stable identifier of a reliability state ("reliable", "open-surface-set", ...), for
  /// logs and machine-readable reports.
  static const char* GetNavigationReliabilityName(NavigationReliability reliability);

  /// Closure-check counts behind GetNavigationReliability(); all zero on a navigable solid.
  int GetBoundaryEdgeCount() const;
  int GetNonManifoldEdgeCount() const;
  int GetReversedEdgeCount() const;

  void ComputeBBox() override;

  int DistancetoPrimitive(int, int) override { return 99999; }
  const TBuffer3D& GetBuffer3D(int reqSections, Bool_t localFrame) const override;
  void GetMeshNumbers(int& nvert, int& nsegs, int& npols) const override;
  int GetNmeshVertices() const override;
  void InspectShape() const override {}
  TBuffer3D* MakeBuffer3D() const override;
  void Print(Option_t* option = "") const override;
  void SavePrimitive(std::ostream&, Option_t*) override {}
  void SetPoints(double* points) const override;
  void SetPoints(Float_t* points) const override;
  void SetSegsAndPols(TBuffer3D& buff) const override;
  void Sizeof3D() const override {}

  Double_t DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;
  Double_t DistFromInside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                          Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;
  bool Contains(const Double_t* point) const override;
  /// Trivial non-BVH Contains looping over all surfaces; kept for debugging and
  /// cross-validation of the BVH-accelerated path (see O2Tessellated::Contains_Loop).
  bool Contains_Loop(const Double_t* point) const;
  /// Diagnostic hook: the parity answer for one explicit shooting \a direction (normalized
  /// internally), bypassing the re-shoot policy that Contains() applies. On a closed, consistently
  /// oriented 2-manifold the answer is independent of the direction, so a point where it *does*
  /// depend on the direction is one where the surface set has a defect the ray happens to hit.
  /// That distinction is what separates "a re-shoot would rescue this" from "the geometry has a
  /// hole"; measure with it before attributing either.
  bool ContainsAlongDirection(const Double_t* point, const Double_t* direction) const;
  /// Trivial non-BVH DistFromOutside/DistFromInside looping over all surfaces. They are both the
  /// oracle the BVH paths are cross-validated against (they must agree exactly: same hits, same
  /// arithmetic, only a different visiting order) and the baseline the BVH speedup is measured
  /// against. \a stepmax is honoured exactly as by the overrides.
  Double_t DistFromOutside_Loop(const Double_t* point, const Double_t* dir,
                                Double_t stepmax = TGeoShape::Big()) const;
  Double_t DistFromInside_Loop(const Double_t* point, const Double_t* dir,
                               Double_t stepmax = TGeoShape::Big()) const;
  Double_t Safety(const Double_t* point, Bool_t in = kTRUE) const override;
  void ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) const override;
  Double_t Capacity() const override;

  /// The Add*Surface calls this solid was built from, in order. Public and const so that the
  /// converter tooling and the tests can inspect what a solid actually carries; the only way to
  /// add to it is to call Add*Surface.
  const std::vector<BVHSurfaceRecord>& GetSurfaceRecords() const { return fRecords; }

 private:
  /// The containment answer shared by Contains() and Contains_Loop(): one parity shot on a solid
  /// the closure check calls Reliable, and a majority vote over several spread directions on one
  /// it does not (CodeReview_Fable.md Section 4.4). \a useBVH selects the accelerated candidate
  /// set. Both public entry points route through here so that the re-shoot policy can never make
  /// them differ -- their bit-identical agreement is the branch's strongest self-check.
  bool containsByParity(const Double_t* point, bool useBVH) const;

  /// Discard the derived state and replay fRecords through the Add*Surface methods, then
  /// CloseShape(). Used by the Streamer when reading. Returns false (having left the solid
  /// undefined, i.e. NavigationReliability::Undetermined) when there is nothing to replay or when
  /// a record does not rebuild -- an unusable solid must never present itself as a usable one.
  bool RebuildFromRecords();

  struct Impl;
  Impl* fImpl = nullptr; //! private bounded-surface implementation

  /// The persistent state: everything else is rebuilt from it. See BVHSurfaceRecord.
  std::vector<BVHSurfaceRecord> fRecords;

  /// The source model's declared tolerance in cm; 0 when unknown. See SetModelTolerance.
  double fModelTolerance = 0.;

  ClassDefOverride(O2BVHSurfaceSolid, 3) // BVH surface-bounded shape class
};

} // namespace base
} // namespace o2

#endif