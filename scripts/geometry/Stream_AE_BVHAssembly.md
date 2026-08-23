# Stream AE — O2BVHAssembly, a BVH behind the assembly-shape contract

Roadmap item (j). Sandro, 2026-08-22:

> *"Assemblies in TGeo are not good CPU + mem wise for many parts (same as TGeoNavigator). I
> would suggest to make a modern high performance O2BVHAssembly class that we can use to wrap
> complicated assemblies such as the oTOF one. Should be straightforward to implement (just
> follow O2Tessellated)."*

This document records the measurement that came before the class, the design the measurement
argued for, and what the class then measured.

## 1. What the navigator actually asks an assembly

Read from the ROOT sources (`v6-36-10-alice2`, `geom/geom/src/TGeoShapeAssembly.cxx` and
`TGeoNavigator.cxx`), because the answer decides whether a shape-level class can help at all.

`TGeoNavigator` never *rests* inside an assembly: `FindNextBoundaryAndStep` opens with
`while (fCurrentNode->GetVolume()->IsAssembly() && fLevel) CdUp();`, and `SearchNode` /
`CrossBoundaryAndLocate` descend through one whenever they reach it. So from the navigator's point
of view an assembly volume is a *daughter node like any other*, and the three questions it is
asked all go through its shape:

| navigator site | shape call |
| --- | --- |
| `SearchNode` locating a point | `TGeoShapeAssembly::Contains` |
| `FindNextDaughterBoundary` stepping towards it | `TGeoShapeAssembly::DistFromOutside` |
| `TGeoNavigator::Safety` | `TGeoShapeAssembly::Safety(point, kFALSE)` |

The daughter identity comes back not as a return value but through the volume:
`Contains` calls `fVolume->SetCurrentNodeIndex(id)` / `SetNextNodeIndex(id)` and `DistFromOutside`
calls `SetNextNodeIndex(i)`; the navigator then reads `GetNextNodeIndex()` and `CdDown`s. That is
the whole mechanism by which a hit in a nested assembly keeps its sensitive volume, and it is a
mechanism a *shape* owns.

**Conclusion: the shape level is the right hook point.** No navigator change is needed, and a
class deriving from `TGeoShapeAssembly` is a drop-in.

## 2. The corpus

The converted oTOF (`scripts/geometry/O2_CADtoTGeo.py`, `Stream_AC_OTOFTraversal.md`): 20
prototype solids, **62 628 leaf placements**, of which 62 560 are `TGeoBBox` and 68 are
`O2BVHSurfaceSolid`. As converted the placements sit in a **nested** tree — 207 assembly volumes,
maximum depth 4, at most 68 daughters each, 66 574 nodes in the expanded tree. Because assemblies
are transparent containers, the same 62 628 placements can equally be written as one **flat**
assembly, which is the shape roadmap item (j) describes and the one that stresses N.

Both were built and measured. World: a box 1.2× the assembly extent, assembly placed at the
identity, geometry closed (so ROOT voxelizes every assembly volume).

## 3. Baseline

Single-threaded, machine load average 0.5, ROOT v6.36.10, `-O2`, aarch64. Query points are
sampled uniformly inside randomly chosen leaf boxes and carried to the assembly frame, so they lie
*on the content* (99.9 % are inside a daughter); rays run from inside the world through such a
point. Every ROOT voxel finder is warmed with 20 000 queries before timing, because `Import`
streams the optimisations back in a state that rebuilds on first use.

"Leaf loop" is the brute-force reference: the minimum over all 62 628 placements, each asked in
its own frame. It is slow on purpose — it is the definition of the right answer, not a candidate
implementation.

| query | nested, ROOT | flat, ROOT | leaf loop (reference) |
| --- | --- | --- | --- |
| `Contains` | 0.71 µs | 4.02 µs | 145 µs |
| `DistFromOutside` | 0.11 µs **wrong** | 0.17 µs **wrong** | 903 µs |
| `Safety(out)` | 2.12 µs | 43.6 µs | 165 µs |
| `TGeoManager::FindNode` | 0.79 µs | 4.29 µs | — |
| transport, `FindNextBoundaryAndStep` | 3.38 µs/step | 27.2 µs/step | — |
| `Import` of the closed geometry | 331 ms | 359 ms | — |
| RSS after import | 281 MB | 337 MB | — |
| `TGeoVolume::Voxelize` of the assembly alone | — | 0.48 s | — |

Flattening costs a factor 5.7 on `Contains`, 21 on `Safety` and 8 on transport. That is the linear
term the BVH is meant to remove, and it is why the nested tree — which is really 207 small
assemblies, each within reach of ROOT's voxel finder — looks so much healthier than the roadmap
implied.

## 4. A correctness defect found while measuring

`TGeoShapeAssembly::DistFromOutside` **returns `Big()` for a point outside the assembly's bounding
box whenever the volume is voxelized**. Mechanism, from the source:

```cpp
   if (!TGeoBBox::Contains(point)) {
      snext = TGeoBBox::DistFromOutside(point, dir, 3, stepmax);
      // Approach bounding box to minimize errors
      snext = TMath::Min(0.01 * snext, 1.E-6);      // <- does not move the point to the box
      ...
      for (i = 0; i < 3; i++) pt[i] += snext * dir[i];
```

The comment says "approach bounding box"; the statement assigns a *distance of about 1e-6* instead
of subtracting it from `snext`, so `pt` stays where `point` was, far outside. The voxel branch
below (`nd >= 5 && voxels`) then calls `TGeoVoxelFinder::SortCrossedVoxels(pt, dir, td)`, whose
`GetIndices` returns false for a point outside the slices, so no candidate is ever produced and
the function falls through to `Big()`. Under the linear branch (`nd < 5 || !voxels`) the same
query is answered correctly.

Measured, not inferred:

* minimal repro (`scratchpad/bvhassembly/minrepro.cxx`): an assembly of *n* unit boxes in a row,
  queried from 40 cm away. `n = 3` (linear branch) → 47, correct. `n = 10` (voxel branch) → `Big()`,
  and also `Big()` from 0.5 cm outside the box.
* oTOF, nested: the shape returns `Big()` on **300/300** rays where the leaf loop finds a crossing
  at ~390 cm.
* oTOF, nested, voxel finders removed from all 207 assembly volumes: the shape agrees with the
  leaf loop on **0/300 disagreements**, at 956 µs/query.

Why it is not catastrophic today: the assembly bounding box is large, so a track that is already
inside it takes the working path. With ray origins forced inside the world, all 500 transport rays
still entered the assembly in both configurations. The exposure is a track approaching from
outside the bounding box with a step long enough to cross it — which is exactly the situation in a
full ALICE 3 world where oTOF is one detector among many.

**O2BVHAssembly deliberately does not reproduce this.** Its `DistFromOutside` enters the BVH with
the ray as given and answers the crossing the daughters actually have; the `_Loop` twin, not ROOT,
is its reference. This is worth a ROOT bug report.

## 5. Design

`O2BVHAssembly : public TGeoShapeAssembly`, in `Detectors/Base`.

* **One BVH primitive per daughter placement.** The daughter shape's bounding box, eight corners
  carried through the node matrix and re-bounded — the same construction
  `TGeoShapeAssembly::ComputeBBox` uses — widened by `kBoxTolerance` = 1e-3 cm and rounded outward
  into float with `nextafterf`. The float tree is therefore a superset of the double geometry and
  can only ever hand on too many candidates. Same discipline as `O2BVHSurfaceSolid::buildBVH`.
* **No per-query dedup marker.** `O2BVHSurfaceSolid` needs one because a surface owns several
  cover boxes; here a daughter owns exactly one primitive, which lives in exactly one leaf, so
  each daughter is handed on exactly once per traversal by construction.
* **`max_leaf_size = 1`.** Resolving a daughter costs a matrix transform plus a full shape query,
  far more than a node box test, and bvh2 enters a leaf's start node without testing its box.
* **Derive, don't replace.** `TGeoVolumeAssembly::AddNode` casts the shape to
  `TGeoShapeAssembly*`; the navigator asks `IsAssembly()`; `ComputeBBox`, `ComputeNormal`,
  `DistFromInside`, `Safety(point, kTRUE)` and the drawing stubs are all inherited unchanged.
  Only `Contains`, `DistFromOutside` and `Safety(point, kFALSE)` are overridden.
* **Identity flows the ROOT way**, through `fVolume->SetCurrentNodeIndex` / `SetNextNodeIndex`.
* **Every accelerated query has a `_Loop` twin** over all daughters in index order, and the twin
  is the definition of the answer — including the tie-break, where the lowest-indexed daughter
  wins. To make that identity hold *by construction* rather than by luck, `DistFromOutside` asks
  every daughter with the query bound `step`, never with a bound that shrinks as the traversal
  proceeds: a shrinking bound would make a daughter's answer depend on visit order. The pruning is
  done one level up, on the ray's `tmax`, where it is order-independent (a box met beyond
  `best + kBoxTolerance` cannot hold a nearer or an equal crossing).
* **Lazy rebuild.** The BVH is transient and the daughter count it was built for is remembered; a
  changed count (or a geometry read back from a file) rebuilds on the next query, the same lazy
  contract as the base class's `fBBoxOK`.

* **`Safety(point, kFALSE)` prunes on one third of the squared box distance, not on the box
  distance itself.** ROOT's `TGeoBBox::Safety` from outside returns the largest per-axis gap, in
  the daughter's own frame, not the Euclidean distance — a point off a corner by (3, 4, 0) is 5
  away and ROOT says 4. Pruning on the Euclidean distance therefore discards daughters that would
  have answered less. The sound bound is Euclidean/sqrt(3), because max_i d_i is at least the norm
  over sqrt(3), a rigid node matrix preserves that norm, and the daughter's box is inside the node
  box. Squared, one third. `TGeoShapeAssembly::Safety` omits exactly this step; section 7 measures
  what it costs ROOT.

Usage: `O2BVHAssembly::MakeBVHAssembly(assemblyVolume)` swaps the shape, and by default leaves the
volume's `TGeoVoxelFinder` in place — section 6 says why, and it is not the reason one would
guess. Call it **after** `TGeoManager::CloseGeometry()`; closing re-voxelizes every volume, which
does not restore ROOT's shape, so an earlier call still leaves the acceleration in place.

## 6. What the shape level does and does not reach

Measured, and it corrects section 1 in one important place.

`TGeoNavigator::SearchNode`, once it has *descended into* the assembly node, does not go back
through the shape: it reads `vol->GetVoxels()` and walks the check list itself, falling back to a
linear loop over all daughters when there is no voxel finder. The shape's `Contains` is called on
the way *in* — when the assembly is a daughter being tested from its mother, and from
`IsSameLocation` / `FindInCluster` — but the point-location descent inside the assembly is
navigator-level code.

The first version of `MakeBVHAssembly` deleted the voxel finder, on the reasoning that nothing
reads it once the shape answers the queries. That reasoning was wrong, and the measurement says so
plainly (oTOF, flat):

| | voxel finder kept | voxel finder dropped |
| --- | --- | --- |
| `FindNode` | 4.29 µs | 814 µs |
| transport | 4.23 µs/step | 746 µs/step |

So `dropVoxels` now defaults to false and is documented as a benchmarking switch. What the shape
level *does* reach is `DistFromOutside`, `Safety(point, kFALSE)`, and the mother-side `Contains`.

## 7. Acceptance

Same corpus, same conditions as section 3 (single-threaded, load average 0.3–0.9), 100 000 query
points, `MakeBVHAssembly` applied to every assembly volume in the geometry.

**Flat — one assembly, 62 628 daughters.**

| query | ROOT | O2BVHAssembly | |
| --- | --- | --- | --- |
| `Contains` | 4.05 µs | **0.95 µs** | 4.3× |
| `DistFromOutside` | 0.17 µs, 0/100 000 crossings found | 4.05 µs, 99 995/100 000 found | correct where ROOT is not; 224× faster than the 905 µs leaf loop it agrees with (0/300 disagreements) |
| `Safety(out)` | 43.4 µs | **1.51 µs** | 28.7× |
| `FindNode` | 4.21 µs | 4.29 µs | 1.0 — section 6 |
| transport | 27.3 µs/step | **4.23 µs/step** | 6.5× |
| build | `Voxelize` 0.48 s | 51.5 ms, 66 592 primitives over 208 volumes | 9× cheaper |
| structure size | — | 4.06 MB | |

**Nested — 207 assemblies, at most 68 daughters each (oTOF as converted).**

| query | ROOT | O2BVHAssembly | |
| --- | --- | --- | --- |
| `Contains` | 0.49 µs | 0.68 µs | 0.72× — *slower* |
| `DistFromOutside` | 0.11 µs, 57/100 000 found | 4.05 µs, 99 995/100 000 found | correct where ROOT is not |
| `Safety(out)` | 2.07 µs | **1.00 µs** | 2.1× |
| `FindNode` | 0.77 µs | 0.89 µs | 0.87× — *slower* |
| transport | 3.32 µs/step | **1.73 µs/step** | 1.9× |
| build | — | 1.3 ms, 3 964 primitives | |
| structure size | — | 0.24 MB | |

Read plainly: **the class earns its keep where N is large.** On the flat 62 628-daughter assembly
it is 4–29× faster on the queries it owns and 6.5× faster in transport. On the nested tree, where
no single assembly has more than 68 daughters, ROOT's voxel finder is *better* at `Contains` than
a BVH traversal — the tree walk costs more than 68 box tests do — and only `Safety` and transport
improve. Transport improves on both, and by more than the per-query numbers suggest, because the
navigator's `Safety` is a large share of a step and because the crossings ROOT loses are found.

`DistFromOutside` looks 20–35× *slower* than ROOT in these tables. It is not: ROOT is returning
`Big()` without doing the work (section 4). The honest comparison is against the leaf loop, which
computes the same answer: 905 µs against 4.05 µs, a factor 224.

## 8. Limitations, and what to decide next

* **No Geant4 transport test.** The class is exercised by ROOT navigation (`FindNode`,
  `FindNextBoundaryAndStep`) and by the 22-case unit suite, not by an `o2-sim` run with hits. The
  identity path is the one ROOT uses, so hits should follow, but that is an argument, not a
  measurement.
* **Single-threaded only.** The BVH is read-only after construction and holds no per-query state,
  so it should be thread-safe where ROOT's shape is; not measured. The lazy rebuild is *not*
  thread-safe, exactly as the base class's `fBBoxOK` is not.
* **The nested corpus does not want this class.** Applying it to every assembly of the converted
  oTOF makes `Contains` and `FindNode` slightly worse. If it is used there at all it should be
  applied selectively, to the assemblies above some daughter count.
* **`Safety(point, kTRUE)` is inherited unchanged** — it resolves the daughter the navigator is
  already in and has nothing to accelerate.
* **The pruning bound costs a factor sqrt(3)** because ROOT's box safety is an axis maximum rather
  than a distance. A tighter bound would need per-shape knowledge of how loose each `Safety` is.

**For Sandro, the open question:** the shape level cannot reach point location inside the
assembly, because `TGeoNavigator::SearchNode` enumerates the assembly's daughters through the
volume's `TGeoVoxelFinder` rather than through the shape. Three ways out, in increasing order of
intrusiveness:

1. Leave it. The voxel finder stays, costs 0.48 s to build on the flat oTOF, and does that job
   about as well as a BVH would; the BVH takes the other three queries.
2. Teach `TGeoVolume` to expose a pluggable daughter-search interface that `SearchNode` asks
   before the voxel finder — a small, upstreamable ROOT change, and the natural place for an
   Embree back end too (roadmap item d).
3. Replace `TGeoNavigator` for these branches, which is the larger conversation the roadmap
   already has open.

And separately: **two ROOT defects fell out of this** (sections 4 and the `Safety` over-pruning
in section 7's test messages, 69 of 2000 grid points). Both deserve a report upstream.
