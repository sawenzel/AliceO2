# Stream AM — carving mothers out, and what the back conversion then gets

Roadmap items (a) and (c), raised by Sandro on 2026-08-27: a STEP meant for real CAD editing
wants **disjoint** solids, TGeo wants **nested** ones, and the question is what the second
direction gets if the export is carved. Measured on two modules — MAG (19 leaf solids, shallow)
and ITS (261, ten levels deep) — with the writer's `--carve-mothers`, the shipped back-conversion
cascade, and `ClosureTest/matbudget_diff.py` against the native `o2-sim` run of the same module in
the same hall. Everything below was run.

**The one-line verdict:** carving and nesting are not two independent switches — **each mother
needs exactly one of them**, the writer is the only side that knows which, and it now says so.
With that, a mostly-flat back conversion is materially correct; with either mechanism applied
blindly to every mother, ITS loses 54 % or 35 % of its material.

---

## 1. The design

Two binary choices, so four cells, plus the one that came out of them:

| | mothers carved? | TGeo nesting restored? | |
| --- | --- | --- | --- |
| **A** | no | yes | the shipped configuration: TGeo carves implicitly, the sidecar puts it back |
| **B** | yes | yes | both mechanisms at once |
| **C** | yes | no | the CAD-native shape — disjoint solids as peers |
| **D** | no | no | negative control |
| **E** | **where it can** | **where it must** | per mother, from the writer's own verdict |

C and D reuse B's and A's STEP, so a variant comparison is a comparison of representations and
never of two exports. The scripts are in `agent-workspace/carve-study/`
(`run_variants.sh`, `run_variant_E.sh`, `score_variants.sh`, `timing_sweep.sh`).

## 2. The result

2000 Fibonacci rays from the origin to r = 800 cm, against the native module in the same hall.
`us/cr` is transport (`FindNextBoundary` + `Step`) per boundary crossing; `us/FN` is point
location, sampled inside the detector (r < 45 cm for ITS). Quiet box, best of three.

**MAG**

| | CSG | leaves | placements | max depth | mean leaf depth | x/X₀ | rays > 1 % | crossings | us/cr | us/FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | 19/19 | 19 | 193 | 10 | 9.19 | −0.000 % | 0 | 13.9 | 1.55 | 0.74 |
| B | 16/20 | 20 | 195 | 9 | 8.13 | −0.055 % | 43 | 13.9 | 2.98 | 0.79 |
| C | 16/20 | 20 | 195 | 7 | 6.82 | −0.000 % | 0 | 13.9 | 14.89 | 1.79 |
| D | 19/19 | 19 | 193 | 7 | 6.63 | **−19.014 %** | 1273 | 7.2 | 11.29 | 0.96 |
| **E** | 16/20 | 20 | 195 | 8 | 7.81 | **−0.000 %** | **0** | 13.9 | **2.83** | 0.83 |

**ITS**

| | CSG | leaves | placements | max depth | mean leaf depth | x/X₀ | rays > 1 % | crossings | us/cr | us/FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | 252/261 | 261 | 1997 | 14 | 12.25 | −0.018 % | 2 | 198.9 | 1.83 | 1.38 |
| B | 279/316 | 316 | 2070 | 12 | 10.60 | **−53.681 %** | 1978 | 64.4 | 5.27 | — |
| C | 279/316 | 316 | 2070 | 10 | 9.22 | **−34.973 %** | 1984 | 14.3 | 8.14 | 4.03 |
| D | 252/261 | 261 | 1997 | 10 | 9.20 | **−30.142 %** | 1984 | 13.2 | 6.72 | — |
| **E** | 273/308 | 308 | 2060 | 12 | 10.46 | **−0.137 %** | 122 | 198.7 | **2.46** | 2.19 |

Native reference: MAG mean x/X₀ 27.356672 at 13.9 crossings, ITS 0.202923 at 198.9.

**The structure is never lost, only the navigation.** ITS sensor placements are
108 / 144 / 180 / 2688 / 3360 / 8232 / 9408 for layers 0–6 in the native geometry and in *every*
variant, C and D included. A geometry that has every sensor and still shows 14.3 crossings against
198.9 is not missing parts; it is being navigated wrongly.

## 3. Why B is the worst cell, and it is not intuition

A daughter placed inside the cavity its mother was carved for is **outside its mother**, and the
navigator never enters it. Minimal reproducer (`carve-study/cavity_probe.py`): a 10 cm cube minus a
2 cm cube at the centre, with that same 2 cm cube of iron nested at the centre, and a ray along +x:

```
W@-19.00 -> MOTHER@-5.00 -> W@-1.00 -> MOTHER@1.00 -> W@5.00     KID entered: False
```

The cavity reads back as *world*. On ITS, where 675 placements were still nested after carving,
that costs 53.7 % of the material and takes the crossings from 198.9 to 64.4.

**Carving also breaks the nesting table on its own.** `bodyOfAssembly` is keyed by name, and
carving splits 18 ITS mothers into several disjoint bodies — `IBConnectorASide` into 8,
`IBConnectorCSide` into 6, `ServicesContainer0/1/2` into 4 each — after which the name no longer
identifies one body and the converter refuses to nest rather than guess. Nesting therefore
collapsed from 1497 of 1996 placements (54 mothers) to 675 of 2069 (36) without anyone asking for
it. B is not a clean cell; it is a hybrid nobody chose.

## 4. Why C fails on ITS and not on MAG

`_carve` subtracts daughter **solids**. An assembly daughter has none — `placed_children` carries
`None` for it — so it was silently skipped, and the mother kept that daughter's whole volume.
Flattened, such a mother becomes a sibling that overlaps everything inside it and wins.

Every ITS mother in that position is an **air envelope**: `ITSUWrapVol0/1/2__body`,
`SpaceFrameVolumeLay3–6__body`, `ServicesContainer0/1/2__body`, `IBConnectorASide/CSide__body` are
all `ITS_AIR$`. `ITSUWrapVol1` and `ITSUWrapVol2` are 38.8 % and 36.2 % occupied by assembly
daughters. Four more — `SpaceFrameVolumeLay3–6` — fail the boolean outright ("OCCT result contains
no solid"), which at least logs.

**Completing the carve is not the way out.** `carve-study/carve_cost.py` counts the leaves a
complete carve would have to fuse: `ITSUWrapVol2` **211 422**, `ITSUWrapVol1` 75 076,
`ITSUWrapVol0` 10 059. That is not an OCCT workload.

MAG escapes only by luck. Its one incomplete mother, `L3CM`, fails the boolean because its 168
daughters *fully consume* it, so the uncarved body it keeps has nothing left to swallow — and TGeo
happened to resolve the overlap in the daughters' favour. C's `−0.000 %` on MAG is a coincidence,
not a property.

## 5. What was changed

**`O2_TGeoToCAD.py`.** `_carve` now returns whether it subtracted *every* daughter, and the media
sidecar carries `carvedComplete` per mother whenever `--carve-mothers` ran. A carve that cannot
finish is **discarded**, not shipped half done: a partially carved mother still has to be nested for
the daughters that were not subtracted, which would put the ones that *were* into cavities — B's
failure, reintroduced through the back door. Either every daughter or none.

**`O2_CADtoTGeo.py`.** A mother the sidecar marks complete is left flat; one marked incomplete keeps
its nesting. No new flag — the sidecar carries the decision, so the two sides cannot disagree.

On ITS this splits 54 mother bodies into **47 complete and 7 incomplete**, and the 7 are exactly
`ITSUWrapVol0/1/2` and `SpaceFrameVolumeLay3–6`. Nesting drops from 1497 placements to 364.

Self-tests: `O2_TGeoToCAD.py --self-test` **129 checks, 0 failures** (four new: complete verdict,
incomplete verdict, a complete carve returns a smaller body, an incomplete carve returns the mother
whole); `O2_CADtoTGeo.py --self-test` **54 checks, 0 failures**, unchanged.

## 6. What carving costs the export and the recognition

| | MAG uncarved | MAG carved | ITS uncarved | ITS carved |
| --- | --- | --- | --- | --- |
| STEP size | 1.7 MB | 2.4 MB | 15.4 MB | 40.3 MB |
| `ADVANCED_FACE` | 486 | 704 | 3332 | 10080 |
| `NEXT_ASSEMBLY_USAGE_OCCURRENCE` | 192 | 192 | 1996 | 1996 |
| writer time | 1.1 s | 26.7 s | 27.6 s | 99.4 s |

(The ITS figures reproduce the Roadmap's 2026-08-27 measurement exactly.) Carving does not touch
instance count, which is what a healing-and-tree-building CAD importer chokes on, so it is still
not the fix for the FreeCAD stall.

**Recognition falls, exactness does not.** ITS goes from 252/261 CSG (96.6 %) to 273/308 (88.6 %),
MAG from 19/19 to 16/20 — and **not one part falls to the mesh tier** in either. What a carved
mother becomes is a solid with many more planar faces: MAG's `L3CCI__body`, a `TGeoPgon` before, is
100 planar faces afterwards and blows the 64-cell decomposition budget, so it ships as an exact
surface solid instead. The leaf count rises (261 → 308 on ITS) because carving splits mothers into
several bodies.

## 7. Depth, and what flattening costs to navigate

E is genuinely flatter than A — ITS mean leaf depth 12.25 → 10.46, MAG 9.19 → 7.81 — and full
flattening (C) goes further, to 9.22 and 6.82. But **in TGeo a flatter tree is slower, not faster**:
MAG transport goes 1.55 → 14.89 µs per crossing from A to C on *identical* solids, and ITS point
location inside the detector goes 1.38 → 4.03 µs. Nesting in a hand-built tree does double duty —
material precedence *and* spatial grouping — and flattening throws the second away while
`TGeoNavigator` still walks siblings.

E keeps the grouping where the fan-out is large and pays only 1.7× (MAG) and 1.34× (ITS) against A.

**`O2BVHAssembly` does not recover it, and Stream AE §6 says why.** Wrapping every assembly of the
flat worlds (verified: `L3CM`, 169 daughters, `TGeoShapeAssembly` → `o2::base::O2BVHAssembly`):

| | transport, plain | transport, BVH | FindNode, plain | FindNode, BVH |
| --- | --- | --- | --- | --- |
| MAG C | 14.89 | 15.23 | 1.79 | 2.07 |
| ITS C | 8.14 | 9.46 | 4.03 | 4.24 |
| ITS E | 2.46 | 2.42 | 2.19 | 2.44 |

Nothing, in both directions. Once `TGeoNavigator` has descended *into* an assembly it enumerates
that assembly's daughters through the volume's `TGeoVoxelFinder`, never through the shape, and
that is where this geometry spends its time — the oTOF corpus gained 6.5× because rays approached
one 62 628-daughter assembly *from outside*. The class still belongs in the pipeline for the
queries it does own; it is not the answer to flattening. Stream AE §8's option 2 — a pluggable
daughter-search interface `SearchNode` asks before the voxel finder — is.

## 7b. How much of the cost is the representation, and does a bigger `O2FlatCSG` budget help?

Asked directly, because §7 leaves two candidates: the flatter tree, and the carved mothers landing
in the exact-surface tier instead of CSG. Measured on MAG's E, whose tree is fixed while only the
solids change. The carved mothers decline the flat path on two separate budgets, and both were
raised in turn:

| budgets | tiers | conversion | transport |
| --- | --- | --- | --- |
| shipped: `PART_MAX_CELLS` 64, `_PART_MAX_FLAT_HALFSPACES` 1024 | CSG 16, surface 4 | — | 3.05 µs/crossing |
| cells 256 | CSG 18, surface 2 | 34 s | — |
| cells 256, halfspaces 4096 | **CSG 19, surface 1** | 68 s | **2.79 µs/crossing** |
| (uncarved + nested, for scale) | CSG 19, surface 0 | — | 1.69 µs/crossing |

So **the representation is the smaller term**: recovering three of the four surface parts into
`O2FlatCSG` — `L3CD__body` as 72 cells / 378 halfspaces, `L3CCO__body` as 168 / 1028 — buys 8.7 %,
about a fifth of the gap to the nested tree. The other four fifths are the flattening.

The one part that still declines does so on **accuracy, not budget**: `L3CCI__body`'s flat cells sit
0.005 cm from the part's own boundary, 2.98e-06 of its 1677 cm diagonal against a 1e-06 bound. That
is exactly the failure `NEXT.md`'s flat-CSG item 3 predicted for a raised budget, and it is the
pipeline refusing correctly.

The budgets are left at their shipped values. The measurement says what raising them is worth, and
it is not worth 2x the conversion time on its own.

## 7c. It is not flatness, and it is not fan-out — it is a large sibling

Sandro's suggestion was to buy the grouping back with **nested assemblies**, which are free of the
objections a solid envelope has: an assembly carries no material and no containment obligation, so
two of them may overlap and only their contents may not. `carve-study/regroup.py` does it — a
median split on daughter box centres, the same construction a BVH build uses, applied before
`CloseGeometry` because that is the only point a TGeo tree can still be restructured.

It makes things **worse**, at both a deep and a shallow setting:

| MAG, flat variant C | transport | point location |
| --- | --- | --- |
| as converted, one assembly of 169 | 15.0 µs/crossing | 1.84 µs |
| regrouped, leaf 16 | 22.1 | 2.15 |
| regrouped, leaf 64 | 21.4 | — |
| (uncarved + nested, for scale) | 1.64 | 0.63 |

ROOT already voxelizes a volume's daughters, so a flat assembly of 169 is not searched linearly;
replacing that one voxel finder with a hierarchy of small ones only adds per-level assembly descent.
This agrees with Stream AE §7, where ROOT's voxel finder beats a BVH traversal on `Contains` below
about 68 daughters.

So the cost is **not** fan-out, and not flatness as such. The decisive measurement is MAG E against
MAG C, which differ in exactly one mother — `L3CM`, nested in E and flat in C:

| | transport |
| --- | --- |
| A, every mother nested | 1.64 µs/crossing |
| E, three mothers flat, `L3CM` nested | 2.83 |
| C, `L3CM` flat too | 15.0 |

One mother of 168 daughters, demoted from mother to sibling, costs **5.3×**. The mechanism is its
bounding box: a mother's solid spans its whole content by construction, so as a sibling it is a
primitive whose AABB is hit by every ray that goes anywhere near the group, and no acceleration —
ROOT's voxels, `O2BVHAssembly`'s BVH, an Embree back end — can prune a box that always tests true.
**Carving does not fix this**: a carved shell has the same AABB as the mother it came from. That is
why the BVH numbers in §7 moved nothing.

The practical consequence is that flattening should be judged per mother by how much of the group
its body spans, not by tier or by depth — and that a big thin shell is the one case where nesting
is worth keeping even though carving succeeded.

## 8. Open

1. **E's residual on ITS: −0.137 %, 122 rays of 2000 above 1 %, and 198.7 crossings against
   198.9.** The median is 4.8e-11, so almost every ray is exact; the tail is at carved boundaries.
   Not diagnosed. A is still the better geometry for simulation, and E is for CAD.
2. **A third option, unmeasured**: separate the mother's two roles — a cheap envelope volume of the
   surrounding medium keeps the navigator's descent, and the mother's own material ships as a
   carved sibling. No material precedence, no depth penalty, no ROOT change.
3. **`--carve-mothers` and assembly daughters remain incompatible in principle**, not just in this
   implementation. Anything that wants a genuinely fully-disjoint STEP has to flatten those
   assemblies in the *writer*, and the leaf counts in §4 say what that costs.
4. **`DistFromInside` on `O2BVHAssembly`** would be new capability, not a fix: ROOT's
   `TGeoShapeAssembly::DistFromInside` is a stub that prints and returns `Big()`, but the navigator
   never reaches it — the two hot call sites in `TGeoNavigator::FindNextBoundary` are guarded by
   `if (!mother->IsAssembly())` and the third is reached only after
   `while (fCurrentNode->GetVolume()->IsAssembly() && fLevel) CdUp();`. Worth implementing for a
   navigator that would ask, which is the VecGeom case, not for TGeo today.
5. **Flatten by span, not by tier**: §7c says a mother's body is affordable as a sibling only when
   its bounding box does not span the group. Unmeasured as a policy.
6. **`O2BVHAssembly` emission from the converter** is drafted but not landed
   (`--bvh-assemblies [MINDAUGHTERS]`, default 68, emitted at the end of `build()`); on this
   evidence it should wait for a corpus where it measurably helps.
