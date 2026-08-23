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

#include "DetectorsBase/O2BVHAssembly.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoVoxelFinder.h"

// the same third-party BVH2 entry point O2Tessellated and O2BVHSurfaceSolid use
#include "bvh2_third_party.h"
#include "bvh2_extra_kernels.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace o2::base;
ClassImp(O2BVHAssembly);

namespace
{
// float BVH types, following the O2Tessellated::BuildBVH pattern
using BVHScalar = float;
using BVHBBox = bvh::v2::BBox<BVHScalar, 3>;
using BVHVec3 = bvh::v2::Vec<BVHScalar, 3>;
using BVHNode = bvh::v2::Node<BVHScalar, 3>;
using BVH = bvh::v2::Bvh<BVHNode>;
using BVHRay = bvh::v2::Ray<BVHScalar, 3>;

/// Every daughter box is widened by this before it is rounded outward into float, the same
/// discipline (and the same value) O2BVHSurfaceSolid applies to its cover boxes: navigation
/// answers are only ever asked to a tolerance of this order, so a leaf box must never be the
/// thing that decides a boundary case.
constexpr double kBoxTolerance = 1.e-3;

/// Round a double outward into float, away from the interval the box encloses.
inline float roundOutward(double value, bool up)
{
  return std::nextafterf(static_cast<float>(value),
                         up ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity());
}

/// A float ray bound that is never *shorter* than the double distance it stands for.
inline float truncateRoundUp(double value)
{
  const float rounded = static_cast<float>(value);
  return rounded < value ? std::nextafterf(rounded, std::numeric_limits<float>::infinity()) : rounded;
}

/// Squared Euclidean distance from \a point to a BVH node box, evaluated in double from the float
/// corners (which convert exactly) and shrunk by a relative guard. Multiply by kSafetyBoundShare
/// before comparing it against a safety; see below.
inline double boxDistanceSq(const BVHBBox& box, const double* point)
{
  double distanceSq = 0.;
  for (int dimension = 0; dimension < 3; ++dimension) {
    const double lower = static_cast<double>(box.min[dimension]);
    const double upper = static_cast<double>(box.max[dimension]);
    const double coordinate = point[dimension];
    if (coordinate < lower) {
      const double gap = lower - coordinate;
      distanceSq += gap * gap;
    } else if (coordinate > upper) {
      const double gap = coordinate - upper;
      distanceSq += gap * gap;
    }
  }
  return distanceSq * (1. - 1.e-12);
}

/// The share of a node's squared box distance that is a sound lower bound on the *Safety* any
/// daughter inside that node can return. One third, and the third is not conservatism -- it is
/// what ROOT's safety convention costs.
///
/// The nearest-daughter traversal must never prune on a value larger than a daughter's answer, or
/// it discards the winner and Safety() returns too much. The natural bound, the Euclidean distance
/// to the node box, is *not* sound, because `TGeoBBox::Safety(point, kFALSE)` does not return the
/// Euclidean distance: it returns the largest per-axis gap, in the daughter's *own* frame. A point
/// off a box corner by (3, 4, 0) is 5 away and ROOT says 4; rotate the daughter and the two
/// numbers stop being related by anything but an inequality.
///
/// The inequality that does hold: for gaps d_i, max_i d_i >= sqrt(sum_i max(d_i, 0)^2) / sqrt(3),
/// so a daughter's axis-max answer is at least its Euclidean distance over sqrt(3); a rigid node
/// matrix preserves that distance; and the daughter's world bounding box is inside the node box,
/// so the node's box distance is no larger still. Squared, that is one third. It holds through
/// nested assemblies too, because each level only re-expresses the same Euclidean distance.
///
/// This is exactly the step `TGeoShapeAssembly::Safety` omits -- it prunes on the unscaled
/// Euclidean gap and therefore returns more than the minimum over its own daughters; see
/// scripts/geometry/Stream_AE_BVHAssembly.md.
constexpr double kSafetyBoundShare = 1. / 3.;

inline bool boxContains(const BVHBBox& box, const double* point)
{
  return point[0] >= static_cast<double>(box.min[0]) && point[0] <= static_cast<double>(box.max[0]) &&
         point[1] >= static_cast<double>(box.min[1]) && point[1] <= static_cast<double>(box.max[1]) &&
         point[2] >= static_cast<double>(box.min[2]) && point[2] <= static_cast<double>(box.max[2]);
}
} // namespace

O2BVHAssembly::O2BVHAssembly() : TGeoShapeAssembly() {}

O2BVHAssembly::O2BVHAssembly(TGeoVolumeAssembly* volume) : TGeoShapeAssembly(volume)
{
  if (volume != nullptr) {
    BuildBVH();
  }
}

O2BVHAssembly::~O2BVHAssembly()
{
  delete static_cast<BVH*>(fBVH);
  fBVH = nullptr;
}

size_t O2BVHAssembly::GetBVHMemory() const
{
  const auto* bvh = static_cast<const BVH*>(fBVH);
  if (bvh == nullptr) {
    return 0;
  }
  return bvh->nodes.size() * sizeof(BVHNode) + bvh->prim_ids.size() * sizeof(size_t);
}

////////////////////////////////////////////////////////////////////////////////
/// BuildBVH -- one primitive per daughter placement.
///
/// The primitive is the daughter shape's own bounding box carried into the assembly frame by the
/// node matrix: its eight corners are transformed and re-bounded, exactly the way
/// TGeoShapeAssembly::ComputeBBox builds the assembly box out of the same corners. It is then
/// widened by kBoxTolerance and rounded outward in float, so the float traversal is a superset of
/// the double geometry and can only ever hand on *too many* candidates.
///
/// One primitive per daughter is what makes the per-query dedup marker of O2BVHSurfaceSolid
/// unnecessary here: a primitive lives in exactly one leaf, so each daughter is handed to the
/// callbacks exactly once per traversal.

void O2BVHAssembly::BuildBVH()
{
  delete static_cast<BVH*>(fBVH);
  fBVH = nullptr;
  fNbuilt = -1;
  fLastNode = -1;
  if (fVolume == nullptr) {
    return;
  }
  ComputeBBox();
  const int nDaughters = fVolume->GetNdaughters();
  fNbuilt = nDaughters;
  if (nDaughters == 0) {
    return;
  }

  std::vector<BVHBBox> boxes;
  std::vector<BVHVec3> centers;
  boxes.reserve(nDaughters);
  centers.reserve(nDaughters);

  double corners[24];
  double master[3];
  for (int index = 0; index < nDaughters; ++index) {
    TGeoNode* node = fVolume->GetNode(index);
    TGeoShape* shape = node->GetVolume()->GetShape();
    // an assembly daughter, or one whose box was never computed, has to produce it first --
    // the same guard TGeoShapeAssembly::RecomputeBoxLast uses
    if (node->GetVolume()->IsAssembly() || TGeoShape::IsSameWithinTolerance(((TGeoBBox*)shape)->GetDX(), 0.)) {
      shape->ComputeBBox();
    }
    ((TGeoBBox*)shape)->SetBoxPoints(corners);
    double lower[3] = {TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};
    double upper[3] = {-TGeoShape::Big(), -TGeoShape::Big(), -TGeoShape::Big()};
    for (int corner = 0; corner < 8; ++corner) {
      node->LocalToMaster(&corners[3 * corner], master);
      for (int dimension = 0; dimension < 3; ++dimension) {
        lower[dimension] = std::min(lower[dimension], master[dimension]);
        upper[dimension] = std::max(upper[dimension], master[dimension]);
      }
    }
    BVHBBox box;
    for (int dimension = 0; dimension < 3; ++dimension) {
      box.min[dimension] = roundOutward(lower[dimension] - kBoxTolerance, false);
      box.max[dimension] = roundOutward(upper[dimension] + kBoxTolerance, true);
    }
    boxes.push_back(box);
    centers.emplace_back(box.get_center());
  }

  typename bvh::v2::DefaultBuilder<BVHNode>::Config config;
  config.quality = bvh::v2::DefaultBuilder<BVHNode>::Quality::High;
  // One daughter per leaf. Resolving a daughter costs a matrix transform plus a full shape query,
  // far more than a node box test (unlike a triangle), and the bvh2 traversal enters a leaf's
  // start node without testing its box -- so a multi-primitive leaf would pay several shape
  // queries for one box test.
  config.max_leaf_size = 1;
  auto* built = new BVH(bvh::v2::DefaultBuilder<BVHNode>::build(boxes, centers, config));
  fBVH = static_cast<void*>(built);
}

void O2BVHAssembly::EnsureBuilt() const
{
  const int nDaughters = fVolume != nullptr ? fVolume->GetNdaughters() : 0;
  if (fNbuilt == nDaughters && (fBVH != nullptr || nDaughters == 0)) {
    return;
  }
  const_cast<O2BVHAssembly*>(this)->BuildBVH();
}

////////////////////////////////////////////////////////////////////////////////
/// Contains

Bool_t O2BVHAssembly::Contains(const Double_t* point) const
{
  EnsureBuilt();
  if (!fBBoxOK) {
    const_cast<O2BVHAssembly*>(this)->ComputeBBox();
  }
  if (!TGeoBBox::Contains(point)) {
    fLastNode = -1;
    return kFALSE;
  }
  const auto* bvh = static_cast<const BVH*>(fBVH);
  if (bvh == nullptr) {
    fLastNode = -1;
    return kFALSE;
  }

  int best = -1;
  double local[3];
  std::vector<size_t> stack;
  stack.reserve(64);
  stack.push_back(0); // the bvh2 root node
  while (!stack.empty()) {
    const auto& node = bvh->nodes[stack.back()];
    stack.pop_back();
    if (!boxContains(node.get_bbox(), point)) {
      continue;
    }
    if (node.is_leaf()) {
      const auto beginPrimitive = node.index.first_id();
      const auto endPrimitive = beginPrimitive + node.index.prim_count();
      for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
        const int daughter = static_cast<int>(bvh->prim_ids[primitive]);
        // the loop twin takes the lowest-indexed daughter that contains the point, so a candidate
        // that cannot beat the standing answer need not be resolved at all
        if (best >= 0 && daughter > best) {
          continue;
        }
        TGeoNode* geoNode = fVolume->GetNode(daughter);
        geoNode->MasterToLocal(point, local);
        if (geoNode->GetVolume()->GetShape()->Contains(local)) {
          best = daughter;
        }
      }
    } else {
      const auto firstChild = node.index.first_id();
      for (size_t child : {firstChild, firstChild + 1}) {
        if (child < bvh->nodes.size()) {
          stack.push_back(child);
        }
      }
    }
  }

  fLastNode = best;
  if (best < 0) {
    return kFALSE;
  }
  // this is how the daughter identity reaches TGeoNavigator, and through it the hit
  fVolume->SetCurrentNodeIndex(best);
  fVolume->SetNextNodeIndex(best);
  return kTRUE;
}

Bool_t O2BVHAssembly::Contains_Loop(const Double_t* point) const
{
  if (!fBBoxOK) {
    const_cast<O2BVHAssembly*>(this)->ComputeBBox();
  }
  if (!TGeoBBox::Contains(point)) {
    fLastNode = -1;
    return kFALSE;
  }
  double local[3];
  const int nDaughters = fVolume != nullptr ? fVolume->GetNdaughters() : 0;
  for (int index = 0; index < nDaughters; ++index) {
    TGeoNode* geoNode = fVolume->GetNode(index);
    geoNode->MasterToLocal(point, local);
    if (geoNode->GetVolume()->GetShape()->Contains(local)) {
      fLastNode = index;
      fVolume->SetCurrentNodeIndex(index);
      fVolume->SetNextNodeIndex(index);
      return kTRUE;
    }
  }
  fLastNode = -1;
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// DistFromOutside
///
/// The distance from \a point to the nearest daughter crossing within \a step, and the daughter
/// that realises it. Unlike ROOT's, this answer does not depend on whether the point is inside
/// the assembly bounding box: the BVH is entered with the ray as given.
///
/// Every daughter is asked with the *query* bound \a step, never with a bound that shrinks as the
/// traversal proceeds. That costs the daughter shapes their own early exit, but it is what makes
/// the answer independent of visit order and therefore bit-identical to DistFromOutside_Loop --
/// including the tie-break, where the lowest-indexed daughter wins. The pruning that matters is
/// done one level up, on the ray bound, where it is order-independent.

Double_t O2BVHAssembly::DistFromOutside(const Double_t* point, const Double_t* dir, Int_t iact, Double_t step,
                                        Double_t* safe) const
{
  EnsureBuilt();
  if (!fBBoxOK) {
    const_cast<O2BVHAssembly*>(this)->ComputeBBox();
  }
  if (iact < 3 && safe != nullptr) {
    *safe = Safety(point, kFALSE);
    if (iact == 0) {
      return TGeoShape::Big();
    }
    if (iact == 1 && step <= *safe) {
      return TGeoShape::Big();
    }
  }
  const auto* bvh = static_cast<const BVH*>(fBVH);
  if (bvh == nullptr) {
    fLastNode = -1;
    return TGeoShape::Big();
  }

  double best = TGeoShape::Big();
  int bestIndex = -1;
  BVHRay ray(BVHVec3(static_cast<float>(point[0]), static_cast<float>(point[1]), static_cast<float>(point[2])),
             BVHVec3(static_cast<float>(dir[0]), static_cast<float>(dir[1]), static_cast<float>(dir[2])), 0.f,
             truncateRoundUp(step + kBoxTolerance));
  static constexpr bool useRobustTraversal = true;
  bvh::v2::GrowingStack<BVH::Index> stack;
  auto* volume = fVolume;
  bvh->intersect<false, useRobustTraversal>(
    ray, bvh->get_root().index, stack, [&](size_t beginPrimitive, size_t endPrimitive) {
      double local[3];
      double localDir[3];
      for (size_t primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
        const int daughter = static_cast<int>(bvh->prim_ids[primitive]);
        TGeoNode* geoNode = volume->GetNode(daughter);
        geoNode->MasterToLocal(point, local);
        geoNode->MasterToLocalVect(dir, localDir);
        const double distance = geoNode->GetVolume()->GetShape()->DistFromOutside(local, localDir, 3, step);
        if (distance < best) {
          best = distance;
          bestIndex = daughter;
        } else if (distance == best && daughter < bestIndex) {
          bestIndex = daughter;
        }
      }
      // A daughter whose box the ray only meets beyond best + kBoxTolerance cannot cross nearer,
      // and cannot tie either: its true crossing is at least its box entry distance.
      if (bestIndex >= 0) {
        ray.tmax = std::min(ray.tmax, truncateRoundUp(best + kBoxTolerance));
      }
      return false; // keep traversing
    });

  fLastNode = bestIndex;
  if (bestIndex < 0 || best >= step) {
    return TGeoShape::Big();
  }
  volume->SetNextNodeIndex(bestIndex);
  return best;
}

Double_t O2BVHAssembly::DistFromOutside_Loop(const Double_t* point, const Double_t* dir, Double_t step) const
{
  if (!fBBoxOK) {
    const_cast<O2BVHAssembly*>(this)->ComputeBBox();
  }
  double best = TGeoShape::Big();
  int bestIndex = -1;
  double local[3];
  double localDir[3];
  const int nDaughters = fVolume != nullptr ? fVolume->GetNdaughters() : 0;
  for (int index = 0; index < nDaughters; ++index) {
    TGeoNode* geoNode = fVolume->GetNode(index);
    geoNode->MasterToLocal(point, local);
    geoNode->MasterToLocalVect(dir, localDir);
    const double distance = geoNode->GetVolume()->GetShape()->DistFromOutside(local, localDir, 3, step);
    if (distance < best) {
      best = distance;
      bestIndex = index;
    }
  }
  fLastNode = bestIndex;
  if (bestIndex < 0 || best >= step) {
    return TGeoShape::Big();
  }
  fVolume->SetNextNodeIndex(bestIndex);
  return best;
}

////////////////////////////////////////////////////////////////////////////////
/// Safety
///
/// From inside there is nothing to accelerate: ROOT resolves the safety of the daughter the
/// navigator is currently in, by walking the current-node chain, and that is inherited unchanged.
/// From outside the answer is the smallest daughter safety, and that is a nearest-neighbour query
/// the BVH answers by ordered descent with a running best.

Double_t O2BVHAssembly::Safety(const Double_t* point, Bool_t in) const
{
  if (in) {
    return TGeoShapeAssembly::Safety(point, in);
  }
  EnsureBuilt();
  if (!fBBoxOK) {
    const_cast<O2BVHAssembly*>(this)->ComputeBBox();
  }
  const auto* bvh = static_cast<const BVH*>(fBVH);
  if (bvh == nullptr) {
    return TGeoShape::Big();
  }

  double best = TGeoShape::Big();
  std::vector<std::pair<double, size_t>> stack;
  stack.reserve(64);
  stack.emplace_back(boxDistanceSq(bvh->nodes[0].get_bbox(), point), size_t(0));
  while (!stack.empty()) {
    const auto entry = stack.back();
    stack.pop_back();
    if (entry.first * kSafetyBoundShare >= best * best) {
      continue;
    }
    const auto& node = bvh->nodes[entry.second];
    if (node.is_leaf()) {
      const auto beginPrimitive = node.index.first_id();
      const auto endPrimitive = beginPrimitive + node.index.prim_count();
      for (auto primitive = beginPrimitive; primitive < endPrimitive; ++primitive) {
        const int daughter = static_cast<int>(bvh->prim_ids[primitive]);
        const double safety = fVolume->GetNode(daughter)->Safety(point, kFALSE);
        if (safety <= 0.) {
          return 0.;
        }
        best = std::min(best, safety);
      }
    } else {
      const auto firstChild = node.index.first_id();
      const size_t children[2] = {firstChild, firstChild + 1};
      double distancesSq[2] = {TGeoShape::Big(), TGeoShape::Big()};
      for (int side = 0; side < 2; ++side) {
        if (children[side] < bvh->nodes.size()) {
          distancesSq[side] = boxDistanceSq(bvh->nodes[children[side]].get_bbox(), point);
        }
      }
      // push the farther child first so the nearer one is popped, and prunes, first
      const bool leftIsFarther = distancesSq[0] >= distancesSq[1];
      const int order[2] = {leftIsFarther ? 0 : 1, leftIsFarther ? 1 : 0};
      for (int side = 0; side < 2; ++side) {
        const int which = order[side];
        if (children[which] < bvh->nodes.size() && distancesSq[which] * kSafetyBoundShare < best * best) {
          stack.emplace_back(distancesSq[which], children[which]);
        }
      }
    }
  }
  return best;
}

Double_t O2BVHAssembly::Safety_Loop(const Double_t* point, Bool_t in) const
{
  if (in) {
    return TGeoShapeAssembly::Safety(point, in);
  }
  double best = TGeoShape::Big();
  const int nDaughters = fVolume != nullptr ? fVolume->GetNdaughters() : 0;
  for (int index = 0; index < nDaughters; ++index) {
    const double safety = fVolume->GetNode(index)->Safety(point, kFALSE);
    if (safety <= 0.) {
      return 0.;
    }
    best = std::min(best, safety);
  }
  return best;
}

////////////////////////////////////////////////////////////////////////////////
/// MakeBVHAssembly

O2BVHAssembly* O2BVHAssembly::MakeBVHAssembly(TGeoVolumeAssembly* volume, bool dropVoxels)
{
  if (volume == nullptr) {
    return nullptr;
  }
  auto* shape = new O2BVHAssembly(volume);
  volume->SetShape(shape);
  if (dropVoxels && volume->GetVoxels() != nullptr) {
    // measured to be a bad idea outside a benchmark: TGeoNavigator::SearchNode reads the finder
    // itself once it is inside the assembly. The same clone guard TGeoVolume::Voxelize uses keeps
    // a cloned volume's shared finder alive.
    if (!volume->TestBit(TGeoVolume::kVolumeClone)) {
      delete volume->GetVoxels();
    }
    volume->SetVoxelFinder(nullptr);
  }
  return shape;
}
