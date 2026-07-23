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

  /// Add an exact planar disk (or annulus when holeRadius > 0) centred at \a center, e.g. a
  /// cylinder or cone end cap. \a axisU and \a axisV must be orthonormal; the outward normal is
  /// axisU x axisV.
  bool AddPlanarDiskSurface(const Point3D& center, const Point3D& axisU, const Point3D& axisV, double radius,
                            double holeRadius = 0.);

  /// Add a cylindrical lateral surface of given \a radius around \a axis. \a centerPoint is the
  /// height reference (h = 0) on the axis, the height range runs along the (normalized) axis and
  /// phi is measured from \a referenceAxisU (projected perpendicular to the axis). With
  /// \a innerWall set, the surface bounds a hole and its outward normal points towards the axis.
  bool AddCylindricalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                             double radius, double heightMin, double heightMax, double phiStart = 0.,
                             double phiSweep = 6.283185307179586, bool innerWall = false);

  /// Add a spherical surface of given \a radius around \a center, trimmed to a polar range
  /// (theta from the +polarAxis pole) and a phi sweep from \a referenceAxisU. The default
  /// arguments give a full, self-closing sphere.
  bool AddSphericalSurface(const Point3D& center, const Point3D& polarAxis, const Point3D& referenceAxisU,
                           double radius, double thetaMin = 0., double thetaMax = 3.141592653589793,
                           double phiStart = 0., double phiSweep = 6.283185307179586, bool innerWall = false);

  /// Add a conical lateral surface whose radius runs linearly from \a radiusAtMin at
  /// \a heightMin to \a radiusAtMax at \a heightMax along \a axis; one radius may be zero (apex
  /// cone). Frame and innerWall conventions as for AddCylindricalSurface.
  bool AddConicalSurface(const Point3D& centerPoint, const Point3D& axis, const Point3D& referenceAxisU,
                         double radiusAtMin, double radiusAtMax, double heightMin, double heightMax,
                         double phiStart = 0., double phiSweep = 6.283185307179586, bool innerWall = false);

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