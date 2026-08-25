# Stream AK — `O2FlatCSG`, the flat-DNF halfspace solid, and what it is measured to be worth

Rung **R5** of [`Handoff_FlatCSG.md`](Handoff_FlatCSG.md), executed 2026-08-25 in ten tasks against
[`Design_FlatCSGSolid.md`](Design_FlatCSGSolid.md) and [`Plan_FlatCSGSolid.md`](Plan_FlatCSGSolid.md).
This is the programme record: what was built, every number the design's §9 asked for, the emission
policy those numbers imply, and the honest remainder. The per-task engineering reports live in the
git-ignored `.superpowers/sdd/Plan_FlatCSGSolid/` workspace and are not the record — this is.

## 1. The headline

`o2::base::O2FlatCSG` is a `TGeoShape` that stores a solid as a **union of intersection-cells over
signed implicit halfspaces** — a flat two-level DNF, depth two always — accelerated by a BVH over
sub-boxes of those cells. It is what a part decomposed by `csg/decompose.py` actually *is*, with
nothing padded into a native primitive to make ROOT able to hold it.

**Ten detector parts that declined against the tree budget now convert, every one at `dV_sym = 0`
exactly**, and every previously accepted candidate is bit-identical. Corpora went

| corpus | before R5 | after R5 |
| --- | --- | --- |
| PIPE | 166 / 9 / 1 | **166 / 9 / 1** |
| ITS | 250 / 14 / 0 | **255 / 9 / 0** |
| TPC | 163 / 9 / 0 | **165 / 7 / 0** |
| ABSO | 29 / 0 / 0 | **29 / 0 / 0** |
| TRD | 1167 / 3 / 0 | **1170 / 0 / 0** (100 %) |
| known-source | 1775 / 1775 | **1785 / 1785** |

and the composite-sourced recognition rate moved 153/188 → **163/188 (87 %)**.

The measurement the rung existed for: against the *same ten parts emitted as plain
`TGeoCompositeShape` unions of the same cells*, the flat solid is **3.2× to 34× faster on
`Contains`**, **4.5× to 30× faster on `Safety`**, **0.9× to 4.5× on `DistFromInside`**, **0.4× to
2.3× on `DistFromOutside`** — the one kernel where it is not uniformly better — and **0.83× to
7.1× on full geantino transport**. The crossover expressed in the quantity the router actually
tests, the composite's own leaf count, is at **about 45 leaves**; the router's threshold is 64, so
every part the cascade currently sends this way is on the winning side of the measured crossover
with margin, and the policy of §7 is therefore *unchanged from what shipped*, now on evidence
rather than on `Stream_AA`'s prediction.

## 2. What was built

| piece | where |
| --- | --- |
| the shape class: halfspace/cell/box tables, `Contains`, both distances by interval clipping, `Safety`, `Capacity`, `ComputeBBox`, `ComputeNormal`, and a `_Loop` twin for every accelerated query | `Detectors/Base/include/DetectorsBase/O2FlatCSG.h`, `Detectors/Base/src/O2FlatCSG.cxx` |
| the quadric block (plane, sphere, cylinder, cone, elliptic cylinder — ten doubles) and the torus block (canonical `(p, d, R, r)`, quartic roots, 1-Lipschitz range bound) | same |
| the sub-cell subdivision with a rigorous centred-form range bound and the two-purse depth/cubify split rule | `SplitBox`, `HalfspaceRange` |
| the `flatcsg_*.bin` sidecar, v1, written from Python and read in C++ | `csg/flat.py::write_sidecar`, `O2SurfaceSolidIO.cxx::LoadFlatCSG`/`WriteFlatCSG` |
| the emitter: carriers → signed halfspaces, the cone apex plane, the exterior-cone refusal, the per-cell box and its outward probe | `csg/flat.py`, `csg/recognise.py::_match_flat_cells` |
| the description, the two builders, the ROOT realisation *through the sidecar* | `csg/primitives.py::flat_cells`, `_root_flat_csg` |
| the routing (last in the cascade, after the tree path has declined) and the `geom.C` branch with its two `::Fatal`s | `csg/recognise.py::_cascade`, `csg/hook.py`, `O2_CADtoTGeo.py` |
| the twin-parity acceptance gate on **both** emission paths | `csg/emit.py::twin_parity`, `csg/hook.py`, `csg/emit.py::from_json` |
| `--flatcsg` as an X-ray benchmark subject, and the two split knobs as bench options | `Detectors/Base/test/runXRayBenchmark.cxx` |
| 31 ctest cases | `Detectors/Base/test/testFlatCSG.cxx` |

Two structural decisions are worth restating because everything else rests on them.

**The twins are the definition of the answer.** Every accelerated query has a `_Loop` twin that
walks all cells and all halfspaces with no BVH and no active-list pruning, and the tests require
**bit identity**, not agreement to a tolerance. That is what makes an acceleration bug — the kind
that only bites on some rays — findable at all.

**The cell bounding box is a correctness obligation, not a cost one.** An intersection of
halfspaces does not bound itself, so the converter is the only thing that can say where a cell
ends, and `CloseShape` builds boxes strictly inside the box it is given. A cell reaching past its
box makes the accelerated queries and the twins describe different solids. The converter therefore
probes 216 points outside each declared box and **refuses the part** when one of them is still
inside the cell — never widens the box, which would ship the phantom material (§9.7).

## 3. The ten parts

`dV_sym = 0` exactly on every one, zero ROOT-vs-CAD containment disagreements over 4000
corroborated points each, zero known-source disagreements against the original `TGeoShape`.

| corpus | part | cells | halfspaces | tree leaves | splits | boundary gap | sidecar |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ITS | `IBCYSSFlangeC` | 12 | 96 | 72 | 8 | 4.2e-13 cm | 10 388 B |
| ITS | `SpaceFrameVolumeLay3__body` | 17 | 111 | 91 | 16 | 1.7e-12 cm | 12 208 B |
| ITS | `SpaceFrameVolumeLay4__body` | 17 | 111 | 91 | 16 | 1.7e-12 cm | 12 208 B |
| ITS | `SpaceFrameVolumeLay5__body` | 17 | 111 | 91 | 16 | 1.7e-12 cm | 12 208 B |
| ITS | `SpaceFrameVolumeLay6__body` | 17 | 111 | 91 | 16 | 1.7e-12 cm | 12 208 B |
| TPC | `TPC_OCGEM` | 12 | 72 | 72 | 11 | 1.8e-11 cm | 7 988 B |
| TPC | `TPC_OCGEM__mirrored` | 12 | 72 | 72 | 11 | 1.8e-11 cm | 7 988 B |
| TRD | `BREF1` | 47 | 282 | 282 | 35 | 1.7e-10 cm | 31 228 B |
| TRD | `VolTOFrail` | 37 | 227 | 92 | 19 | 2.3e-12 cm | 25 088 B |
| TRD | `B077__body` | 20 | 122 | 122 | 13 | 1.7e-10 cm | 13 500 B |

"Tree leaves" is what the same cells cost as a `TGeoCompositeShape`, and it is **not** the
halfspace count: `_cell_leaf` folds six planes into one `TGeoBBox`, so `VolTOFrail`'s 227
halfspaces are 92 leaves while `BREF1`'s 282 are 282. That gap is itself a finding — the two
numbers the design asked the crossover to be reported against are related by a per-cell folding
factor between 1.0 and 2.5 on this corpus, so they are not interchangeable predictors (§4).

Per-cell halfspace counts run 5 to 16; the largest is `IBCYSSFlangeC`'s three 16-halfspace cells.

### Ten, not the twenty-one the demand table names

`Design_FlatCSGSolid.md` §9 item 1 asks for "the 21 parts". That table spans the five detector
corpora **plus ALICE3 and Bagger**, and this rung converted the five corpora. All eleven
non-shipping parts, in one place:

| part(s) | why it did not ship |
| --- | --- |
| ALICE3 `ST1782525_01`, `ST1A38495_01`, `ST1A38526_01`, `ST0923290_01#b11`, `ST1829909_01` ×4 (**8**) | not in the five detector corpora. ALICE3 was converted separately for the decline catalogue and **none of the eight ships there either** — all eight still decline, and none produces a `flatCells` candidate. |
| Bagger `Bucket` (**1**) | not in the five corpora; Bagger is the gate model. It is also in the boundary-gap class §11 puts out of scope — 4.11 cm, 0.211 of its diagonal. |
| ITS `BPSupportLowerCollar` (**1**) | refused by the cell-box containment probe: cell 6 of 21 still holds 0.67 cm past the CAD piece's own box. §9.7. |
| ITS `IBCYSSFlangeA` (**1**) | refused by `decompose.PART_MAX_CELLS = 64`, which design §8 froze for this rung. §9.12 prices the raise, and finds it lands in the boundary-gap class too. |

So: **10 shipped, 9 out of corpus scope, 2 refused by acceptance working as designed.** Nothing was
waived, and design §2 and §9 now carry the same reconciliation.

## 4. Measurement 1 — the crossover

Every one of the ten was emitted **both ways** and X-ray-benched against itself. The composite
emission is the `unionOfCells` description the union path would have produced with
`_PART_MAX_LEAVES` raised; it is reconstructed from the `notes["occCells"]` the flat candidate
already carries, which are literally the same `_cell_leaves(piece, tol, diag, whole_part=False)`
results the union matcher uses, and the reconstruction was **verified against a real raised-budget
scratch conversion** of both `TPC_OCGEM` parts — the twelve cells came back identical, digest for
digest. The raise was never committed.

Both subjects answer one sample set per part, drawn once and partitioned by the composite;
`--raster 32` (3072 rays), 4096 points and 4096 rays per kernel, median of 9 timed passes after 2
warmups, warm cache, single-threaded. **Both subjects returned exactly the same crossing and step
counts on every part**, which is the check that the two are the same solid.

Ratios are composite ÷ flat, so **> 1 means the flat solid is faster**.

| part | cells | halfspaces | leaves | `Contains` | `DistFromOutside` | `DistFromInside` | `Safety` | transport |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `TPC_OCGEM` | 12 | 72 | 72 | **6.8** | 0.87 | 2.4 | **7.7** | **1.70** |
| `TPC_OCGEM__mirrored` | 12 | 72 | 72 | **6.8** | 0.87 | 2.5 | **7.9** | **1.62** |
| `IBCYSSFlangeC` | 12 | 96 | 72 | **3.2** | 0.53 | 0.90 | **4.5** | 0.83 |
| `SpaceFrameVolumeLay3__body` | 17 | 111 | 91 | **23.8** | 1.39 | 2.3 | **9.3** | **3.29** |
| `SpaceFrameVolumeLay4__body` | 17 | 111 | 91 | **23.0** | 1.41 | 2.4 | **9.4** | **3.25** |
| `SpaceFrameVolumeLay5__body` | 17 | 111 | 91 | **33.3** | 2.29 | 2.5 | **8.6** | **4.10** |
| `SpaceFrameVolumeLay6__body` | 17 | 111 | 91 | **33.6** | 2.28 | 2.5 | **8.6** | **3.77** |
| `B077__body` | 20 | 122 | 122 | **6.9** | 1.09 | 1.3 | **11.9** | **2.27** |
| `VolTOFrail` | 37 | 227 | 92 | **5.3** | 0.39 | 1.8 | **13.8** | 0.89 |
| `BREF1` | 47 | 282 | 282 | **30.8** | 2.33 | 4.5 | **30.4** | **7.12** |

Absolute costs for the two extremes, ns per call, composite → flat:

```
TPC_OCGEM   Contains  413 ->   61   DistOut 1297 -> 1483   DistIn 1250 ->  530   Safety  2961 ->  384   transport  4552 -> 2685 ns/ray
BREF1       Contains 2063 ->   67   DistOut 2868 -> 1233   DistIn 4897 -> 1088   Safety 17162 ->  565   transport 12661 -> 1778 ns/ray
```

### Which count predicts it

Regressing `log(ratio)` on `log(count)` over the ten parts:

| kernel | vs `log cells` | vs `log halfspaces` | vs `log tree leaves` |
| --- | --- | --- | --- |
| `Contains` | r +0.28, R² 0.08 | r +0.25, R² 0.06 | r +0.48, R² 0.23 |
| `DistFromOutside` | r +0.13, R² 0.02 | r +0.09, R² 0.01 | r +0.48, R² 0.23 |
| `DistFromInside` | r +0.43, R² 0.18 | r +0.32, R² 0.11 | r +0.54, R² 0.29 |
| `Safety` | **r +0.91, R² 0.82** | r +0.84, R² 0.70 | r +0.91, R² 0.82 |
| transport | r +0.33, R² 0.11 | r +0.29, R² 0.08 | **r +0.65, R² 0.43** |

**Read the R² column before the crossover column.** For transport it is 0.11 against cells and
0.08 against halfspaces: the fitted crossovers below summarise ten points, they do not establish a
law, and R6's census is what would turn them into one.

**Of the two counts the design named, the cell count is the better predictor** — it beats the
halfspace count on every one of the five kernels, because a cell is what the BVH indexes and a
halfspace is only what a leaf then evaluates, and the per-cell halfspace count barely varies
(5 to 16, mean about 6). But neither explains the transport spread. **The best predictor is the
composite's own leaf count**, which is not a property of the flat solid at all — it is a property
of the thing it is being compared against, and that is the honest reading: what varies across
these ten parts is mostly how expensive the *composite* is, not how cheap the flat solid is.

Crossovers, from the fitted lines (ratio = 1):

| kernel | in cells | in halfspaces | in tree leaves |
| --- | ---: | ---: | ---: |
| `Contains` | below 1 | below 1 | 10 |
| `DistFromOutside` | 8 | 41 | 81 |
| `DistFromInside` | 3 | 11 | 26 |
| `Safety` | 2 | 10 | 12 |
| transport | 3 | 18 | **45** |

`Stream_AA_FlatCSG.md` predicted a crossover at ≥ 10 cells. The measurement says **3 cells on
transport** — the prediction was conservative by a factor of three. In leaves, the number the
router tests, it is 45 against a routing threshold of 64.

## 5. Measurement 2 — the split knobs, and the defaults they set

`SetSplitDepth` swept over 3..10 and `SetMinBoxFraction` over 0.002..0.2, on the three largest
shipped parts (`BREF1` 47 cells, `VolTOFrail` 37, `B077__body` 20) and three more for control.
Every knob setting answers one sample set drawn once per part, so the rows are comparable by
construction.

`BREF1`, 47 cells / 282 halfspaces, four representative rows:

| depth | minfrac | boxes | mean active | max active | empty | build | BVH | `Contains` | `Safety`(out) | `DistFromOutside` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 0.002 | 120 968 | 2.08 | 4 | 0 | 96 ms | 7 561 kB | 91 ns | 1795 ns | 2169 ns |
| 6 | 0.010 | 18 911 | 2.24 | 6 | 0 | 13.9 ms | 1 182 kB | 97 ns | 1028 ns | 1836 ns |
| **4** | **0.050** | **772** | **3.32** | **6** | **0** | **0.5 ms** | **48 kB** | **98 ns** | **535 ns** | **1293 ns** |
| 3 | 0.200 | 101 | 4.93 | 6 | 0 | 0.06 ms | 6 kB | 85 ns | 337 ns | 809 ns |

`VolTOFrail`, 37 cells / 227 halfspaces:

| depth | minfrac | boxes | mean active | max active | build | BVH | `Contains` | `Safety`(out) | `DistFromOutside` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 0.002 | 47 872 | 1.55 | 4 | 36 ms | 2 992 kB | 267 ns | 785 ns | 3910 ns |
| 6 | 0.010 | 816 | 4.19 | 6 | 0.5 ms | 51 kB | 185 ns | 446 ns | 1358 ns |
| **4** | **0.050** | **204** | **4.46** | **6** | **0.13 ms** | **13 kB** | **151 ns** | **333 ns** | **1029 ns** |
| 3 | 0.200 | 51 | 5.55 | 7 | 0.03 ms | 3 kB | 116 ns | 234 ns | 945 ns |

`B077__body`, 20 cells / 122 halfspaces, the only part in the set where the depth cap still binds
at `minfrac = 0.05`:

| depth | minfrac | boxes | mean active | empty | build | BVH | `Contains` | `Safety`(out) | `DistFromOutside` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 0.002 | 79 836 | 1.13 | 12 734 | 62 ms | 4 990 kB | 200 ns | 1160 ns | 6813 ns |
| 6 | 0.010 | 7 000 | 1.44 | 536 | 4.9 ms | 438 kB | 180 ns | 753 ns | 3975 ns |
| **4** | **0.050** | **1 320** | **2.22** | **0** | **0.8 ms** | **83 kB** | **162 ns** | **526 ns** | **2973 ns** |
| 3 | 0.050 | 884 | 2.36 | 0 | 0.6 ms | 55 kB | 147 ns | 443 ns | 2576 ns |
| 3 | 0.200 | 124 | 4.05 | 0 | 0.10 ms | 8 kB | 123 ns | 319 ns | 1730 ns |

**The result is monotone and it is not what the design expected: on these parts, every extra box
costs and none pays.** Query cost rises with the box count on all four kernels, memory and build
time rise with it linearly, and the leaf active-list length — the thing subdivision exists to
shrink — is already down at 1 to 5 halfspaces before the subdivision does anything, because a cell
of this rung carries only 5 to 16 halfspaces in the first place. Evaluating the two or three a
coarse box leaves undecided is cheaper than the extra BVH descent that would have removed them.
This is design §11 risk 3 arriving from the opposite direction: not "the bound stays undecided so
the box count grows", but "the box count grows and buys nothing".

**The chosen defaults are `fSplitDepth = 4` (was 6) and `fMinBoxFraction = 0.05` (was 0.01).**
Two reservations belong with that choice and are recorded in full at **§9.10a** (the size
floor is a fraction of the *part's* diagonal, not the cell's) and **§9.10b** (the cost curve
was monotone to the edge of the swept range, so I stopped on an argument, not on the data).

`0.05` rather than `0.1` or `0.2`, which measured better still: at 0.2 the subdivision produces
2 to 3 boxes per cell, which is a BVH over the cell AABBs with rounding — the arrangement design
§4.1 deliberately rejected, and one with no headroom at all for a part whose cells are large
relative to the part. `0.05` is the coarsest setting that still leaves a real sub-cell structure
(2.5 to 19 boxes per cell here) and it sits on the flat part of the cost curve: 0.01 → 0.05 wins
on every kernel on all six parts and costs nothing anywhere.

Depth `4` rather than `6` because at `minfrac = 0.05` the cap binds on exactly one of the six
parts, `B077__body`, where 6 buys 2.7× the boxes and 37 % more `DistFromOutside` for nothing; and
rather than `3`, which is better again on that one part by 15 %, because the cap exists for the
part whose cells are large relative to the part and none of these ten is that part.

Measured effect of the change, on the ten-part crossover with everything else held fixed:
transport ratios moved from **0.29–4.92 (median 1.31)** to **0.83–7.12 (median 2.76)**, and the
flat solid's structural memory fell by 4–19× — `BREF1` 2.6 MB → 141 kB, `IBCYSSFlangeC` 2.4 MB →
132 kB, `TPC_OCGEM` 112 kB → 27 kB, which is *below* the composite's 29 kB. Two parts
(`VolTOFrail` 0.29 → 0.89, `IBCYSSFlangeC` 0.42 → 0.83) moved from clearly losing on transport to
within 20 % of parity.

**Two ctest cases had to be pinned rather than defaulted**, and this is the right change rather
than a workaround: `a_box_wholly_inside_a_cell_carries_no_active_halfspaces` and
`safety_is_sound_when_the_inside_bound_is_actually_nonzero` are about what the subdivision *can*
produce, which is a different question from what the shipped defaults are tuned for. Both now set
both knobs explicitly.

## 6. Measurement 3 — `Safety` quality, and why the depth knob is the wrong lever

`Safety` is a **structural lower bound**, deliberately: there is no point-to-quadric distance
formula in general and this class does not pretend to one. Outside, it is the distance to the
nearest sub-cell box — a real bound, since every point of the solid is in some box. Inside, only a
box with an *empty* active list contributes; a point in an undecided box gets the sound answer
`0.`, because that box's active list says nothing about where the boundary sits within it.

Measured against the exact closed-form distance on five analytic fixtures, 400 000 uniform trials
each, at the previous defaults (depth 6, minfrac 0.01):

| fixture | boxes | empty boxes | inside queries getting `0.` | mean `Safety`/true | median where positive |
| --- | ---: | ---: | ---: | ---: | ---: |
| box 20×2×2 | 512 | 120 | **76.6 %** | 0.023 | 0.085 |
| cube 4×4×4 | 64 | 8 | **87.6 %** | 0.012 | 0.084 |
| cylinder r3 h10 | 64 | 8 | **84.2 %** | 0.019 | 0.099 |
| sphere r5 | 64 | 8 | **76.2 %** | 0.029 | 0.104 |
| plate 10×10×4 with an r1.5 hole | 256 | 40 | **83.2 %** | 0.017 | 0.087 |

On the **shipped parts** it is worse, because the size floor stops the split before any leaf
detaches — their cells are small relative to the *part*, and the floor is a fraction of the part's
diagonal, not the cell's (§9.10a):

| part | boxes | empty | inside queries getting `0.` |
| --- | ---: | ---: | ---: |
| `B077__body` | 7 000 | 536 | 77.6 % |
| `SpaceFrameVolumeLay3__body` | 6 540 | 252 | 92.0 % |
| `IBCYSSFlangeC` | 17 964 | 164 | 95.2 % |
| `SpaceFrameVolumeLay5__body`, `TPC_OCGEM`, `BREF1`, `VolTOFrail` | 3 134 / 730 / 18 911 / 816 | **0** | **100 %** |

### This is a tessellation-resolution property, not a depth limit

Because each cell's bounding box is its exact extent, **every leaf that touches a face keeps that
face's plane undecided**. With `n_i` splits on axis `i`, the fraction of the cell's volume in a box
that touches at least one face is

```
1 - Π_i (1 - 2^(1 - n_i))
```

and the measurement matches the closed form: the 20×2×2 slab splits (5, 2, 2) at depth 6, which
predicts 76.6 % against a measured 76.6 %; the cube splits (2, 2, 2), predicting 87.5 % against a
measured 87.6 %.

The consequence is arithmetic, and the depth sweep confirms it exactly:

| depth (minfrac 1e-7) | undecided, cube | boxes, cube | undecided, slab | boxes, slab |
| ---: | ---: | ---: | ---: | ---: |
| ≤ 5 | 100 % | 32 | 100 % | 256 |
| 6 | 87.6 % | 64 | 76.6 % | 512 |
| 8 | 72.0 % | 224 | 63.7 % | 1 680 |
| 10 | 50.7 % | 704 | 44.6 % | 4 848 |
| 12 | 33.0 % | 2 096 | 24.7 % | 14 304 |
| 14 | 23.2 % | 5 760 | 18.7 % | 38 128 |
| 16 | **14.9 %** | **15 312** | **12.5 %** | **96 064** |

Getting below 15 % needs about five splits per axis — depth ≈ 15, and **two orders of magnitude
more boxes**, which §5 has just measured to be the single most expensive thing this class can do.
And it would not even buy a *good* safety: even at depth 16 the mean `Safety`/true is 0.09 to 0.11.

**So: tuning `fSplitDepth` to chase this number is the wrong conclusion, and the defaults chosen in
§5 deliberately do not chase it.** The real fix is design §11's per-halfspace 1-Lipschitz safety —
for each halfspace, a bound on the distance to its own surface (exact for a plane and a sphere,
available from the normalised gradient for the rest, already exact for the torus since the class
stores its signed distance), and the minimum over the cell's halfspaces. That gives a positive
bound at every interior point regardless of the tessellation, and it is the first item of the
successor stream.

### The one thing this measurement did change immediately

`checkKnownSource.py` skipped any point within `--skin` of *either* shape's boundary, and it asked
the emitted shape for that distance. For an `O2FlatCSG` that question has no answer, so 78–100 %
of interior points were being skipped: `SpaceFrameVolumeLay5__body` scored 3 383 of 20 000 points
and `VolTOFrail` 6 172, against a median of 0 skipped over the other 255 ITS parts. The band is now
taken on the **source alone** whenever the emitted shape's `Safety` is a structural lower bound
(`safety_is_a_true_distance`), and the result records `skinnedBoth` so the number can never be read
as the symmetric one when it was not. **This makes the test strictly stricter**, not looser: every
point the source is confident about is now scored, and a disagreement there is a failure rather
than a skip. Measured on the re-run corpora: all ten flat parts now score **20 000 of 20 000 with
0 skipped and 0 disagreements** — `SpaceFrameVolumeLay5__body` went 3 383 → 20 000 — while the
other 1 775 parts keep `skinnedBoth = true` and are untouched.

### And one cheap acceleration the measurement paid for

`Safety`'s outside branch pushed both BVH children unconditionally in tree order, so the standing
best improved in tree order rather than in distance order and subtrees a nearer box would have cut
off were expanded anyway. It now computes both children's node-box lower bounds — which the pop
would have computed regardless — visits the nearer first, and drops a child whose bound already
meets the standing best. The answer is untouched (it is the same `min`, and bit identity with
`Safety_Loop` is what the ctest asserts). **Measured 3.7× to 7.3× on `Safety`, median 5.6×**, over
the ten shipped parts.

## 7. The emission policy

**A part ships as `o2::base::O2FlatCSG` only when the cascade's tree path has already declined
it**, which today means its cells would make a `TGeoCompositeShape` wider than
`_PART_MAX_LEAVES = 64` leaves. Below that it ships as the plain `unionOfCells` composite, exactly
as it did before this rung. Nothing that converted before R5 changed representation.

That policy is now a measurement rather than a guess, and the measurement supports it:

- the transport crossover is at **45 composite leaves**; the threshold is 64, so every part routed
  this way is measured to be on the winning side, with margin;
- below the crossover the composite really is faster on transport, so the flat path must stay
  **last** in the cascade — it is not a strict improvement and must not be treated as one;
- `Contains` and `Safety` are faster flat everywhere in the range measured, including well below
  the crossover, so a workload dominated by those two would justify a lower threshold. Transport is
  not that workload, and transport is what a simulation pays.

Two honest qualifications the policy carries:

1. **`DistFromOutside` is the kernel where the flat solid is not uniformly better** — 0.39× to
   2.33×, losing on 3 of 10 parts. Its BVH walk gathers pieces from every box the ray meets within
   `[0, step]` with no front-to-back ordering and no early termination once the nearest entry is
   already found. That is the same defect §6 just fixed in `Safety`, in a harder setting, and it is
   the second item of the successor stream.
2. **Memory is not free.** The flat solid's structural bytes at the new defaults run 27 kB to
   195 kB against the composite's 29 kB to 113 kB — comparable on the small parts, up to 4× on
   `B077__body`. The sidecar on disk is 8 to 31 kB. At the old defaults it was 4–19× worse, which
   is most of why the defaults moved.

## 8. Acceptance, floors, and byte-identity

The three-test discipline of `Handoff_Recognition.md` §3 binds unchanged: the OCCT symmetric
difference, the oracle gate, and `checkKnownSource.py` against the original `TGeoShape`. The false-
accept guard `accept.contains_disagreements` runs on this path too — `BRepAlgoAPI_Cut` can report
`IsDone()` with zero solids in **both** directions (`ST0923290_01#b19`), so a symmetric difference
of 0.0 is not on its own evidence — and every one of the ten parts reports 4000 corroborated
points.

On top of those, this rung adds a fourth gate that is the flat path's own: **twin parity**.
`Contains` is compared against `Contains_Loop` over the union of the shape's declared cell boxes
doubled about its centre, so 7/8 of the points fall outside every box, which is where a cell that
reaches past its own box lives. It runs on **both** emission paths — the live one in `csg/hook.py`
and the deferred completion in `emit.from_json`, which is the path `runOracleGate.py` forces — and
a part it refuses gets no `shape_*.root` and no `flatcsg_*.bin`, so `geom.C` ships it one tier
down. The sample count now scales as `max(4000, 500 × cells)` instead of a flat 20 000: the defect
is per-cell and the points are spread over the union of every cell's box, so a fixed count thinned
out as a part gained cells (`BREF1`'s 47 cells got 425 points each) while a one-cell part got
20 000 it did not need.

**Nothing that converted before this rung changed.** Every corpus was run twice — once on this
branch and once against a pristine `git archive` of `scripts/geometry` at the pre-rung commit — and
the `sha256` of each part's candidate description compared: **1775 of 1775 bit-identical**, plus
the ten new ones. Bagger separately: 7 of 7. The split-default change of §5 was re-checked the same
way against the pre-change conversion, since it alters `CloseShape` and therefore the `Contains`
the twin-parity gate samples: **1785 of 1785 bit-identical**, and every corpus row unmoved.

### The floors, at handover

| gate | floor | measured |
| --- | --- | --- |
| `csg/emit.py --self-test` | 353, six digest tables green | **353**, all six green and unchanged |
| `checkKnownSource.py --self-test` | 17 | **17** |
| `O2_CADtoTGeo.py --self-test` | 54 (18+8+10+12+6) | **54** |
| `O2_TGeoToCAD.py --self-test` | 107 | **107** |
| `runOracleGate.py --self-test` | 17 | **17/17** |
| ctest `testFlatCSG` | 31 | **31** (10 908 407 assertions) |
| ctest `testBVHSurfaceSolid` | 113 | **113** (323 406 assertions) |
| ctest `testBVHAssembly` | 22 | **22** (115 638 assertions) |
| fixtures gate | exit 0, shape 10/10 | **exit 0, 10/10** |
| Bagger gate | exit 0, 13/13 + 7/7, CSG 7 bit-identical | **exit 0, 13/13 + 7/7**, bit-identical |

## 9. The honest remainder

Everything below was found during the rung, checked against the code at the close of it, and left
in place deliberately. Nothing here is a suspicion; each one is a measurement or a reading of the
source.

**9.1 A shape read off a file IS closed, and this was measured rather than assumed.** `fBoxes`,
`fActive`, `fBVH` and `fClosed` are transient by design (§7 of the spec: a stored BVH is a second
thing that can disagree with the data it describes), so the object ROOT reconstructs has no boxes
of its own. The gap is closed by a **`#pragma read` rule** in `DetectorsBaseLinkDef.h` — not a
hand-written `Streamer`, which would have meant `#pragma link C++ class … -` and giving up
automatic schema evolution over the two member vectors — whose code calls `CloseShape()` on every
object read and shouts if it refuses.

**It fires.** Measured at the close of the rung on a current build, raw ROOT with no harness
re-close, on both paths and from both languages:

```
PyROOT  TFile::Get("shape")      class=o2::base::O2FlatCSG closed=1 nboxes=129 ncells=12 capacity=87.6233
PyROOT  TGeoManager::Import      class=o2::base::O2FlatCSG closed=1 nboxes=129 ncells=12 capacity=87.6233
C++     TGeoManager::Import      class=o2::base::O2FlatCSG closed=1 nboxes=129 ncells=12
C++     TFile::Get<O2FlatCSG>    class=o2::base::O2FlatCSG closed=1 nboxes=129 ncells=12
```

(`TPC_OCGEM`, written by the converter, put into a `TGeoManager`, exported, and read back in a
fresh process.) **The `closed=0` reading this section used to quote is stale** — it predates the
read rule, and it is withdrawn. The three explicit re-closes in `checkKnownSource.py`,
`harness::loadShapeFromRootFile` and `geom.C` are now belt-and-braces rather than the mechanism,
and they are kept: each is cheap, idempotent (`CloseShape` on a closed shape rebuilds the same
boxes) and each guards a reader that must not be silently slow. **Nothing here is open work.**

**9.2 `Safety`'s inside branch is a lower bound and mostly returns `0.`** — 78 % to 100 % on the
shipped parts (§6). Sound, and a transport-speed defect no correctness test can see. The fix is the
1-Lipschitz per-halfspace safety, not a depth bump.

**9.3 `DistFromOutside` has no front-to-back box ordering.** It gathers pieces from every sub-cell
box the ray meets in `[0, step]` and merges afterwards. That is why it is the one kernel the flat
solid can lose (§7 item 1). `Safety`'s outside branch had the same shape and gained 5.6× from
fixing it (§6).

**9.4 The cubify ceiling bounds a root-to-leaf PATH, not a cell**, and its reach depends on the
cell's shape rather than on its worst aspect ratio: at `kMaxCubifySplits = 10` it covers a 2048:1
rod but only a 64:1 plate. That gap is a factor of **32**, about one and a half orders of magnitude
— not the "four orders of magnitude" the source comment and design §4.2 claimed until this rung.
Both texts are corrected; at 16 the gap would be 256, about two and a half orders.

**9.5 A torus halfspace is stored in absolute coordinates, so grazing chords are lost in proportion
to the lever arm.** Task 3 measured it: a 0.5 cm torus with the ray starting 375 cm away loses
chords up to ~0.009 cm, and the quartic's root residuals grow to 4e-5 cm there against 3e-11 close
in. No shipped part is affected (no flat part carries a torus yet), but **any part placed far from
the origin is**, and the fix — storing the torus in a local frame and carrying the ray into it — is
a contained change nobody has needed yet.

**9.6 Two debug-only guards are weaker than they read.** `HalfspaceRoots`'s unit-direction check is
a blunt `|d·d − 1| < 1e-9` **absolute** tolerance and is `assert`-only, so it is absent from a
release build; and the `CloseShape` cell-AABB face sampler has only **negative** evidence — Task 5's
`#undef NDEBUG` run shows it does not false-fire, and nothing exercises a deliberately
out-of-bounds cell through it. The release-path equivalent of the second one does exist and is
tested both ways (`recognise._flat_box_holds_cell`, §9.7), which is why this is a note and not an
item.

**9.7 Two real parts are refused by the cell-box containment probe, and that refusal is right.**
ITS `BPSupportLowerCollar` (cell 6 of 21, the cell still holds 0.67 cm past its own box) and TPC
`TPC_ORH` (cell 1 of 5, 8.44 cm). The declared box is the CAD piece's **own** bounding box, so a
point outside it is outside the part: the accelerated `Contains`, which says there is no material
there, is right, and `Contains_Loop` is inventing material because the cell's halfspaces do not
close the cell up. **Widening the box would ship that phantom material**, and it is the obvious
wrong fix. The probe is a sampler and does not claim to be a proof — 216 points on the six faces
will find a cell that runs off through a face, not one that pokes out through a corner smaller than
the grid spacing.

**9.8 The two emission paths do not run byte-identical gates.** Both run `twin_parity`. The live
path in `csg/hook.py` additionally runs `crosscheck_contains` against the CAD solid, which the
deferred `emit.from_json` path cannot, because it has PyROOT and no OCCT. The shared gate is the one
that catches the defect `SetCellBBox` cannot detect for itself; the extra one is a second,
independent sampling of the same property. Worth knowing when reading two runs of the same part.

**9.9 `Capacity` is inherited from OCCT.** It is the sum of the cells' `GProp` volumes on the
pieces, per spec §11 item 4 — OCCT's measurement of the CAD body, not an independent measurement of
the shipped solid. Every capacity quoted anywhere in this document carries that caveat.

**9.10 `crosscheck_bbox` is systematically the cell-box margin for a flat part**, because
`ComputeBBox` is the union of the retained sub-cell boxes and those reach the declared box. It is
reported and never gated, and it is conservative in the safe direction, but a reader expecting
~1e-7 cm for a CSG part should know why it is ~1e-2 to ~1 cm here.

**9.10a `fMinBoxFraction` is a fraction of the PART's bounding-box diagonal, not the cell's.**
`minSize = fMinBoxFraction × diagonal`, and `diagonal` is measured over the union of every cell's
box. For a part whose cells are scattered over a large volume — `BREF1`'s union diagonal is
1101 cm, `B077__body`'s 1428 cm, against cells of a few cm — the same fraction is a wildly
different constraint per cell, and it is a blunt instrument. **This directly qualifies §6's
"their cells are small relative to the part": that is a statement about this knob's denominator as
much as about the parts**, and it is part of why coarse settings win so uniformly in §5. Making it
per-cell is a behaviour change beyond "pick the defaults from the data" and I did not take it.

**9.10b The split defaults are the number I would most want overruled or confirmed.** §5 gives the
structural argument for stopping at `0.05`, and I stand behind it, but I want to be plain that the
data did not say stop: the cost curve is monotone all the way to `0.2`, the coarsest setting swept,
and I stopped one notch earlier because beyond it the subdivision degenerates into a BVH over the
cell AABBs with no headroom for a part whose cells are large relative to the part — the case the
whole mechanism exists for, and one that none of these ten parts happens to be. Six parts from
three detectors chose this default. If R6's census turns up a large-celled part, **that** part
should choose it, not these.

**9.11 The rescale test's two bounds bind in different regimes.** The relative assertion
(`1e-13 × max(1, |twin|)`) binds below `|twin| = 10` and the aggregate absolute bound (`1e-12`)
above it; both keep 8×–80× headroom over the measured `1.24e-14`. Deliberate, and recorded so
nobody reads the looser one as the whole statement.

**9.12 `decompose.PART_MAX_CELLS = 64` is now the binding constraint, not the flat budgets — and
raising it alone buys nothing.** The largest part shipped is 47 of 256 cells and 282 of 1024
halfspaces — the flat budgets are a quarter used at worst — while `IBCYSSFlangeA` is refused by the
*decomposition's* cell budget at 59 terminal cells plus a pending queue. Design §8 froze
`csg/decompose.py` for this rung, so it stayed frozen.

Priced with a scratch run at `PART_MAX_CELLS = 256`: `IBCYSSFlangeA` **does** decompose, into **157
cells**, well inside the flat budget of 256, in 26.6 s of wall clock for the two `IBCYSSFlange`
parts at 869 MB peak. And it then declines for a *different* reason — `the flat cells's boundary is
1.23 cm from the part's (0.0356 of the part's 34.5445 cm diagonal, over 1e-06)` — which is the
**boundary-gap class of §11**, not a budget. So the two next-rung items are coupled: the raise is a
prerequisite for splitting at every carrier crossing, and splitting at every carrier crossing is
what makes the raise worth anything. Neither is worth doing alone.

## 10. Geant4: what maps and what does not

For whoever picks up the `TVirtualGeoConverter` seam (spec §11.5):

A cell whose halfspaces are **all** `side == "interior"` and whose part is a pure union maps onto
`G4MultiUnion` directly — that is exactly what `G4MultiUnion` is, a union of placed solids with a
voxelised accelerator, and it is the shape this class generalises. The general case does **not**
map. A complemented halfspace makes a cell non-convex, and there is no `G4` counterpart: the
choices are a `G4BooleanSolid` tree per cell — which reintroduces precisely the cost this class
removes, since a boolean tree charges per node whatever the point is — or a `G4VSolid` subclass
mirroring this one. **The second is the honest option and it is its own rung.**

Measured on the ten parts shipped here, the correspondence covers more than expected: **nine of
the ten are pure unions of all-interior cells** — 111 of `SpaceFrameVolumeLay3__body`'s 111
halfspaces are interior, and so are all of `BREF1`'s 282, `VolTOFrail`'s 227, `B077__body`'s 122
and `TPC_OCGEM`'s 72. Only `IBCYSSFlangeC` complements anything, and it complements a lot: 45 of
its 96 halfspaces, in every one of its 12 cells. So on this corpus the union-of-convex-cells case
is the *common* one and the general case is the exception — which makes the `G4MultiUnion` route
worth taking first and the `G4VSolid` subclass the thing that finishes the job rather than the
thing that starts it.

One caveat before anyone counts on that. An all-interior cell is convex, but `G4MultiUnion` unions
*placed `G4VSolid`s*, and a convex cell over six planes plus a cylinder is a `G4VSolid` only if
some native class happens to be it. `_cell_leaf` already answers that question for the ROOT side —
it is what produces the "tree leaves" column of §3 — and its answer there is an *intersection* of
native primitives, not one primitive. So the correspondence is between a cell and a
`G4IntersectionSolid`, and `G4MultiUnion` buys the union level, not the cell level.

## 11. The 24 boundary-gap declines, and why they are now affordable

Twenty-four parts still decline because a dihedrally-convex splitter piece is **not** a cell of the
carrier arrangement — its boundary sits 2–23 % off, median 2 %, and Bagger's rams are the type
specimen. `Design_FlatCSGSolid.md` §2 put them explicitly out of this rung's scope and they are
still out.

What changed is the price of the fix. Splitting at **every carrier crossing** rather than only at
trusted concave edges would make every piece a genuine cell, and it multiplies the cell count —
which, before this rung, was the thing that refused a part, because a wide union of cells is a wide
boolean tree and the tree budget is 64 leaves. It is no longer: a flat solid does not care how many
cells it has, the flat budgets are a quarter used at worst, and §5 has just measured that the
per-cell cost of this class is dominated by the BVH descent, which grows as `log(cells)`. **The
cell-count multiplication the fix causes is no longer what refuses a part.** That is the whole
reason spec §11 item 6 says the fix has become affordable, and it is now a measured statement
rather than an expectation.

The one thing the fix will need that this rung does not provide: `decompose.PART_MAX_CELLS`
(§9.12) has to move with it — and §9.12's price probe shows the dependency runs the other way too.
`IBCYSSFlangeA` with the budget raised decomposes into 157 cells and then declines *here*, on the
boundary gap, at 3.56 % of its diagonal. The budget and the splitter are one change, not two.

## 12. What the successor stream should do, in order

1. **The 1-Lipschitz per-halfspace `Safety`** (spec §11, §6 above). A positive bound at every
   interior point, independent of the tessellation. It retires §9.2, it retires the
   `checkKnownSource` skin special case, and it is the single largest transport-speed item this
   class has.
2. **Front-to-back ordering and early exit in `DistFromOutside`** (§9.3). The same fix `Safety`
   just took, in a harder setting; it is the one kernel the flat solid can lose.
3. **Split at every carrier crossing** (§11), with `decompose.PART_MAX_CELLS` raised to match
   (§9.12). This is where the remaining 24 declines are.
4. **R6, the breadth census** — MAG, TOF, MFT, MCH, MID, FT0, FV0, ZDC, EMC, PHS, CPV, HMP. The
   crossover of §4 rests on ten parts from three detectors, at transport R² 0.11 against the cell
   count; R6 is what makes it a statement about the geometry rather than about ten parts. It is
   also what would settle the split defaults (§9.10b) on a part the current six do not represent.

**Not on this list, and it was on an earlier draft of it:** deciding the streamer. §9.1 measures
the `#pragma read` rule firing on every read path, so there is nothing to build.
