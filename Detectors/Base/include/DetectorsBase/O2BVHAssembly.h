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

#ifndef ALICEO2_BASE_O2BVHASSEMBLY_
#define ALICEO2_BASE_O2BVHASSEMBLY_

#include "TGeoShapeAssembly.h"

class TGeoVolumeAssembly;

namespace o2
{
namespace base
{

/// A BVH-accelerated stand-in for ROOT's `TGeoShapeAssembly`.
///
/// An assembly volume answers three navigation questions about its daughters -- am I in one of
/// them, which one do I reach next along this ray, and how far may I move without touching any --
/// and `TGeoShapeAssembly` answers all three by walking the daughters, either linearly or through
/// the volume's `TGeoVoxelFinder`. This class answers the same three through a BVH over the
/// daughter placements' axis-aligned bounding boxes, so the cost grows with log N rather than N.
///
/// It *derives* from `TGeoShapeAssembly` rather than replacing it, which is what makes it a
/// drop-in: `TGeoVolumeAssembly::AddNode` casts the shape to `TGeoShapeAssembly*` to invalidate
/// the bounding box, the navigator asks `IsAssembly()`, and the daughter identity keeps flowing
/// out the way ROOT's does -- through `fVolume->SetCurrentNodeIndex` / `SetNextNodeIndex`, which
/// is how `TGeoNavigator` learns which daughter to descend into and therefore how a hit keeps its
/// sensitive volume. Everything not on the hot path (`ComputeBBox`, `ComputeNormal`,
/// `DistFromInside`, `InspectShape`, the drawing stubs) is inherited unchanged.
///
/// Use it either by construction on an existing assembly volume,
/// \code
///   auto* shape = o2::base::O2BVHAssembly::MakeBVHAssembly(assemblyVolume);
/// \endcode
/// which swaps the volume's shape and drops the now-redundant voxel finder, or by handing the
/// volume to the constructor and calling `TGeoVolume::SetShape` yourself. Call it *after*
/// `TGeoManager::CloseGeometry()`: closing re-voxelizes every volume, which would re-create the
/// voxel finder this class makes unnecessary (it does not restore ROOT's shape, so the
/// acceleration itself survives).
///
/// Every accelerated query has a `_Loop` twin that walks all daughters in index order. The twins
/// are the definition of the answer, including the tie-break: the *lowest-indexed* daughter wins,
/// which is what ROOT's linear loop also yields. The tests require bit identity between each
/// query and its twin.
///
/// See scripts/geometry/Stream_AE_BVHAssembly.md for the measurements and for the one place this
/// class deliberately does *not* reproduce ROOT: `TGeoShapeAssembly::DistFromOutside` returns
/// `Big()` for a point outside the assembly bounding box whenever the volume is voxelized, and
/// this class returns the crossing the daughters actually have.
class O2BVHAssembly : public TGeoShapeAssembly
{
 public:
  O2BVHAssembly();
  /// Build over the daughters \a volume has *now*; later daughters trigger a lazy rebuild.
  explicit O2BVHAssembly(TGeoVolumeAssembly* volume);
  ~O2BVHAssembly() override;

  O2BVHAssembly(const O2BVHAssembly&) = delete;
  O2BVHAssembly& operator=(const O2BVHAssembly&) = delete;

  /// (Re)build the acceleration structure from the volume's current daughter list.
  void BuildBVH();
  /// Number of daughter placements the current BVH covers, -1 if it was never built.
  int GetNbuilt() const { return fNbuilt; }
  /// Bytes held by the BVH nodes and the primitive-index permutation.
  size_t GetBVHMemory() const;
  /// The daughter index the last accelerated query resolved to, -1 if none did.
  int GetLastNodeIndex() const { return fLastNode; }

  // ---- the accelerated part of the TGeoShapeAssembly contract ----------------------------
  Bool_t Contains(const Double_t* point) const override;
  Double_t DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t* safe = nullptr) const override;
  Double_t Safety(const Double_t* point, Bool_t in = kTRUE) const override;

  // ---- the reference twins: same answer, all daughters, index order ----------------------
  Bool_t Contains_Loop(const Double_t* point) const;
  Double_t DistFromOutside_Loop(const Double_t* point, const Double_t* dir,
                                Double_t step = TGeoShape::Big()) const;
  Double_t Safety_Loop(const Double_t* point, Bool_t in = kTRUE) const;

  /// Replace \a volume's shape by an O2BVHAssembly and return it. With \a dropVoxels the
  /// volume's `TGeoVoxelFinder` is deleted: nothing reads it once this shape is in place, and on
  /// a 62 628-daughter assembly building it costs half a second.
  static O2BVHAssembly* MakeBVHAssembly(TGeoVolumeAssembly* volume, bool dropVoxels = true);

 private:
  /// Rebuild if the volume gained or lost daughters since the last build (or was never built).
  /// Same lazy contract as the base class's `fBBoxOK`, and equally not thread-safe while the
  /// geometry is still being assembled.
  void EnsureBuilt() const;

  void* fBVH = nullptr;       //! bvh::v2::Bvh over the daughter placement boxes
  int fNbuilt = -1;           //! daughter count the BVH was built for; -1 = never built
  mutable int fLastNode = -1; //! daughter resolved by the last accelerated query

  ClassDefOverride(O2BVHAssembly, 1) // BVH-accelerated assembly shape
};

} // namespace base
} // namespace o2

#endif
