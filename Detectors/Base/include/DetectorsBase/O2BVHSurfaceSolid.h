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

  /// Whether the BVH acceleration structure has been built (after CloseShape).
  bool HasBVH() const;
  /// Fill the BVH root-node bounding box; returns false when no BVH has been built.
  bool GetBVHRootBounds(Point3D& lower, Point3D& upper) const;
  /// Diagnostic/test hook: number of candidate surface patches whose BVH leaf boxes the given
  /// ray traverses (counted with multiplicity per leaf primitive). Returns -1 without a BVH.
  int CountBVHRayCandidates(const Point3D& point, const Point3D& direction) const;

  /// Whether the closed shape forms a closed 2-manifold (every boundary edge shared by two faces).
  /// Meaningful only after CloseShape(); detects e.g. missing faces.
  bool IsClosed() const;
  /// Whether all shared boundary edges are traversed in opposite directions after CloseShape();
  /// detects e.g. reversed faces (inconsistent outward normals).
  bool IsOrientationConsistent() const;

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
  Double_t Safety(const Double_t* point, Bool_t in = kTRUE) const override;
  void ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) const override;
  Double_t Capacity() const override;

 private:
  struct Impl;
  Impl* fImpl = nullptr; //! private bounded-surface implementation

  ClassDefOverride(O2BVHSurfaceSolid, 1) // BVH surface-bounded shape class
};

} // namespace base
} // namespace o2

#endif