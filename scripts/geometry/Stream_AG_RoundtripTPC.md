# Stream AG — the TGeo → STEP → TGeo round trip on the TPC

`Stream_AD_TGeoToStep.md` built the bench and ran it on the beam pipe. This runs the same bench,
unchanged, on the **TPC** — the largest single Run 3 detector by placement count — to find out what
a second, differently-built detector does to it. The source TGeo is again the exact oracle: the
STEP was generated from it, so the right answer is known by construction.

Measured 2026-08-23 on `swenzel/bvhsurfacesolid`, aarch64, OCCT/pythonOCC 7.9.0, ROOT 6.36.10.
Every number below was run on this machine. No source file was changed for this stream; the bugs
it found are its results.

**The one-line verdict.** Every part that reaches the STEP comes back exact — 1 289 990 `Contains`
samples over 129 parts, **zero disagreements**, and 73/73 recognised CSG shapes at symmetric
difference exactly 0 cm³ — but **only 56 % of the TPC's leaf placements reach the STEP at all**,
because the TPC builds its A/C-side symmetry out of reflecting placements of *assemblies*, and the
exporter's bake-the-solid workaround cannot bake a subtree. The parts it *does* bake are wrong by
up to 1.06 % in volume, and that is an OCCT call in `O2_TGeoToCAD.py`, measured and reproduced in
nine lines.

---

## 1. The source geometry

```bash
# sim env only -- NEVER with the pythonOCC prepends
VMCWORKDIR=$B/stage/share \
$B/stage/bin/o2-sim -n 0 -g boxgen -m TPC --configKeyValues "align-geom.mDetectors=none"
```

1.3 s, 533 MB. Everything below uses the **ideal** `o2sim_geometry.root`.

### Census

| | |
| --- | --- |
| `TGeoVolume` objects in the file | 168 |
| distinct volume *names* | 165 |
| reachable from `cave` | **165** |
| DAG placement edges (`AddNode`) | 2 385 |
| fully expanded placements | **38 983** (ROOT's own "38983 nodes") |
| max depth | 9, max placements of one volume 735 |
| leaf-solid placements / mother-body placements | 36 851 / 464 = **37 315** |

| shape class | reachable volumes |
| --- | ---: |
| `TGeoTubeSeg` | 35 |
| `TGeoShapeAssembly` | 30 |
| `TGeoTube` | 21 |
| `TGeoCompositeShape` | 20 |
| `TGeoTrd1` | 16 |
| `TGeoPcon` | 14 |
| `TGeoBBox` | 13 |
| `TGeoPgon` | 7 |
| `TGeoArb8` | 5 |
| `TGeoCone` | 4 |

Ten classes, all ten mapped by `O2_TGeoToCAD.py`. The 20 composites are 48 subtractions and
29 unions (no intersections) over operands that are `TGeoBBox` (25), **`TGeoHalfSpace` (15)**,
`TGeoTube` (14), `TGeoTrd1` (10), `TGeoPcon` (4), `TGeoPgon` (2), `TGeoTubeSeg` (1).

**The half-spaces are the TPC-specific coverage hole.** Five composites cut themselves with
unbounded planes — `TPC_RR_CU` (2), `TPC_MMHC` (6), `TPC_OMH5` (5), `TPC_IHSS` (1), `TPC_OHS5` (1)
— and the mapper declines all five, because a half-space has no B-rep body. PIPE has none of these;
the full Run 3 census (Stream AD §5) declines 5 composites for this reason in the whole
experiment, and they are exactly these five — `FULL_report.json` names `TPC_IHSS`, `TPC_MMHC`,
`TPC_OHS5`, `TPC_OMH5`, `TPC_RR_CU` and nothing else.

**Three names are carried by two `TGeoVolume` objects each**: `TPC_ENDCAP`, `TPC_SECT` (the A- and
C-side twins, `Detector.cxx:1923-1924`) and `TPC_Strip`. In TPC the twins are geometrically
identical (`TPC_Strip`: two `TGeoPgon`s, both 1350.387288727189 cm³), so nothing is lost — but
`O2_TGeoToCAD.py` keys its solid cache on the volume *name*, so two same-named volumes with
different shapes would silently share one solid. Benign here, latent in general.

### The drift volume

`TPC_Drift` is the detector's **only sensitive volume** (`Detector.cxx:3165-3166`). It is a `TGeoPcon`
with 735 daughters, capacity 99 283 715.31615923 cm³, and the OCCT solid reproduces that to
**relative deviation exactly 0.0** — the same double bit for bit. Its enclosure `TPC_M` is a
`TGeoPcon` at 3.4e-16. The inner-field-cage sandwich (`TPC_IFC1`…`TPC_IFC13`, `TPC_IFEPOX1-4`,
`TPC_PRSTR1-3`) is 35 `TGeoTubeSeg`s and 21 `TGeoTube`s, and every one of them comes back through
the whole round trip as a **recognised CSG `TGeoTubeSeg`/`TGeoTube` at `dV_sym = 0`**. The
precision-critical parts survive exactly.

## 2. The coincidence walk — one duplicate, and it is a typo

The PIPE walk found 81 coincident (volume, world transform) pairs. The same walk on TPC:

```
leaf-solid placements: 36851   mother-body placements: 464   total 37315
distinct (volume, world transform) pairs: 37314
coincident duplicates: 1 extra placement over 1 pair, multiplicity {1: 37313, 2: 1}
   TPC_PRSTR1  x2
```

```
/cave/barrel_1/TPC_M_1/TPC_Drift_1/TPC_IFC_1/TPC_PRSTR1_1
/cave/barrel_1/TPC_M_1/TPC_Drift_1/TPC_IFC_1/TPC_PRSTR1_2
```

both at world translation (0, 0, −177.925) cm, a `TGeoTubeSeg` of 16.335755 cm³. Unlike PIPE's
plies, these are **direct siblings under one assembly**, and the cause is a copy-paste in
`Detectors/TPC/simulation/src/Detector.cxx:1388-1389`:

```cpp
tv100->AddNode(tvpr1, 1, new TGeoTranslation(0., 0., -177.925)); // prepreg strip
tv100->AddNode(tvpr1, 2, new TGeoTranslation(0., 0., -177.925));
```

The two lines directly above it (`tv2`, `tvep1`) and the *same block repeated for the 120° and 240°
segments* (lines 1412-1413 and 1436-1437) all pair −177.925 with +177.925. So the inner field
cage's first segment has **two prepreg strips at z = −177.925 and none at z = +177.925**. It is one
16.3 cm³ `TGeoTubeSeg` of Tedlar in the wrong place, in a 99 million cm³ drift volume.

**`TGeoManager::CheckOverlaps(1e-4)` reports 0 overlaps in the whole of TPC** — including this pair.
The PIPE case was blind because the coincident discs sat in *sibling assemblies* the checker never
cross-compares; this one is blind for a different and more basic reason: two exactly coincident
copies of one shape have zero *penetration depth*, and a checker that measures penetration cannot
see them. Two mechanisms, one conclusion: the placement-coincidence invariant that
`O2_CADtoTGeo.py` enforces finds a class of defect `CheckOverlaps` structurally cannot.

## 3. The write

```bash
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root TPC_dedup.step \
        --dedup-world --report TPC_dedup_report.json
```

| | DAG (`--no-step`) | expanded (`--dedup-world`) |
| --- | ---: | ---: |
| solids built | 130 | 130 |
| volumes with daughters / pure assemblies | 69 / 30 | 69 / 30 |
| components | 2 217 | **38 676** |
| declined | 5 | 5 |
| coincident placements dropped | 0 | **1** (`TPC_PRSTR1_2`, §2) |
| reflecting placements / of which baked | 40 / 3 | **74 / 37** |
| capacity check, median / max | 2.00e-15 / 1.99e-02 | 2.00e-15 / 1.93e-02 |
| wall (converter-reported) | 1.2 s | **6.39 s** |
| wall (process), peak RSS | 2.3 s, 727 MB | **7.5 s, 1.16 GB** |
| file | — | **37.06 MB, 505 426 STEP entities** |

**The STEP writer did not segfault.** Stream AD left the full-geometry
`STEPCAFControl_Writer` SIGSEGV undiagnosed at 21 150 definitions / 74 601 components. TPC writes
130 definitions / 38 676 components / 505 426 entities in 6.4 seconds. That brackets the failure:
it is above TPC's size and at or below the full model's, and it is not simply "many components".

### Declines: 5, all one reason

| n | class | reason | placements affected |
| ---: | --- | --- | ---: |
| 5 | `TGeoCompositeShape` | *unbounded solid: a half-space has no B-rep body of its own* | 732 |

`TPC_RR_CU` alone is 660 of those 732. Mapping `TGeoHalfSpace` needs a bounding box to clip
against, which the enclosing composite always supplies — the operand is never used unbounded in
practice. This is the cheapest coverage win the bench has surfaced so far.

### The real coverage loss: reflected placements of assemblies

The TPC builds its A/C-side symmetry with `diag(1, 1, −1)` placements. `build_world` recurses into
such a child first, then discards the subtree it just built and tries to bake a *single mirrored
solid* instead; when the child is an assembly there is no such solid, so the whole subtree is
dropped (`O2_TGeoToCAD.py:1169-1177`). The discarded assembly label is still an XCAF `NewShape()`,
so it survives in the file as an **orphaned free shape at the un-reflected, un-placed position**.

Reading the written STEP back:

```
free shapes (roots): 38          # cave, plus 37 orphans
leaf definitions:   167          # 130 solids + 37 baked mirrored copies
expanded leaf occurrences: 36 582
   root cave:        20 872
   root TPC_ENDCAP:  13 500      <- the whole readout endcap, one side
   root TPC_rrod:       663  x2
   root TPC_mrod:        26  x34
```

An independent TGeo walk that replays the same rule accounts for every placement, and lands on the
STEP's number exactly:

| | placements |
| --- | ---: |
| kept — leaf solids | 20 483 |
| kept — mother bodies | 352 |
| kept — baked mirror leaves | 37 |
| **kept, total** | **20 872** ← matches the STEP's `cave` root exactly |
| lost — reflected subtree dropped | 16 040 |
| lost — declined half-space solid | 402 |
| lost — coincident duplicate (`--dedup-world`) | 1 |
| **source total** | **37 315** |

**44 % of the TPC's leaf placements do not reach the assembly tree of the STEP.** One reflected
placement of `TPC_ENDCAP` costs 13 500 of them. PIPE has zero reflected placements, so Stream AD
could not see this; the full Run 3 run recorded "98 placements dropped" as a footnote — TPC alone
contributes 37 of them, and they are the expensive ones.

### The baked copies are wrong by up to 1.06 %, and an analytic oracle says so

For the 37 reflecting placements of *leaf* volumes the exporter does bake a mirrored solid, through
`apply_tgeo_matrix` → `BRepBuilderAPI_GTransform` with a `gp_GTrsf` (`O2_TGeoToCAD.py:530-534`).
Comparing each `__mirrored` leaf against its un-mirrored twin **in the same file**:

| baked copy | n | mirrored vol | un-mirrored vol | rel | faces before | faces after |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `TPC_CDCE__mirrored` | 1 | 1085.736636 | 1076.191048 | **8.87e-03** | 6 cylinder, 4 plane | **10 bspline** |
| `TPC_INPLL__mirrored` | 18 | 14.491711 | 14.339313 | **1.06e-02** | 12 cyl, 10 plane, 2 cone | **24 bspline** |
| `TPC_OPLL__mirrored` | 18 | 14.643517 | 14.572017 | **4.91e-03** | 14 cyl, 10 plane, 2 cone | **26 bspline** |

`TPC_CDCE` is a `TGeoPcon`, so this needs no Monte Carlo at all: its analytic
`TGeoShape::Capacity()` is 1076.191048 cm³, the un-mirrored OCCT solid reproduces it to 3e-15, and
**the mirrored copy of the same solid is 0.887 % too large**. A reflection is an isometry; it cannot
change a volume. The nine-line reproducer, on a plain hollow cylinder:

```
original           vol=1005309.649148734  faces={'cylinder': 2, 'plane': 2}
GTransform(z->-z)  vol=1013638.003727279  rel=8.284e-03  faces={'bspline': 4}
Trsf mirror        vol=1005309.649148734  rel=0.000e+00  faces={'cylinder': 2, 'plane': 2}
```

`BRepBuilderAPI_GTransform` re-approximates every analytic carrier as a B-spline and moves the
volume by 0.83 % doing it. `gp_Trsf::SetMirror` + `BRepBuilderAPI_Transform` is exact and keeps the
carriers, and `gp_Trsf` carries a mirror perfectly well — the reason `tgeo_matrix_to_trsf` refuses a
reflection is that a STEP *component location* must be a proper rigid motion, which is a constraint
on placements, not on baking. This is the same defect class as the ThruSections B-splines of
Stream AD §2, on the one path PIPE never exercised, and unlike `_boolean` this path carries **no
volume invariant** — `_check` only asks whether the shape is valid.

## 4. What came back

```bash
python3 scripts/geometry/O2_CADtoTGeo.py TPC_dedup.step --output-folder rt_tpc_out -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report ... --csg-report ...
```

**2 min 40.5 s, peak RSS 866 MB**, exit 0. 167 leaf definitions in, 168 parts out —
`TPC_ORH`, a `TGeoCompositeShape` of a box with five box cuts and two tube cuts, is a compound of
**two disjoint bodies** and the multibody-leaf rule (Stream AC) correctly splits it.

```
Exact-surface extraction (auto): 167/168 leaf solids represented exactly, 1 falls back
tiers: CSG 73, exact surfaces 94, tessellated 1  (of 168 leaf solids)
```

| | as shipped | excluding the 37 mirrored artefacts |
| --- | ---: | ---: |
| parts | 168 | 131 |
| CSG | **73 (43.5 %)** | **73 (55.7 %)** |
| exact surface | 94 | 57 |
| tessellated | 1 | 1 |

**99.4 % exact surfaces, one mesh.** The one mesh is `TPC_WSEG`, a `TGeoTubeSeg − TGeoPgon`, and
the surface extractor's reason is specific: *"2 face(s): planar polygon wire has fewer than
3 edges"* — the cut leaves two degenerate two-edge planar wires. That is a converter-side item with
a named, reproducible cause, not a free-form surface.

The reverse converter also warned about the export's own damage, correctly:

```
WARNING: this CAD model DECLARES 1535 leaf solid placement(s) that coincide exactly with another
placement of the same solid.
  Leaf placements: 36618 declared (multiplicity {1: 34408, 2: 649, 34: 12, 36: 14}) -> 35083 distinct.
  Suppressed 1535 placement edge(s) ... (0 by root-containment, 1535 by coincident-occurrence)
```

The multiplicities 34 and 36 are the 34 orphaned `TPC_mrod` roots piling on top of each other at the
identity. **That is an artefact of §3's dropped-subtree bug, not of the TPC** — TPC's own
coincidence count is 1.

### CSG acceptance is exactly the Tier-1 primitive set

| source class | parts | csg | surface | mesh |
| --- | ---: | ---: | ---: | ---: |
| `TGeoTubeSeg` | 35 | **35** | 0 | 0 |
| `TGeoTube` | 21 | **21** | 0 | 0 |
| `TGeoBBox` | 13 | **13** | 0 | 0 |
| `TGeoCone` | 4 | **4** | 0 | 0 |
| `TGeoTrd1` | 16 | 0 | 16 | 0 |
| `TGeoPcon` | 14 | 0 | 14 | 0 |
| `TGeoPgon` | 7 | 0 | 7 | 0 |
| `TGeoArb8` | 5 | 0 | 5 | 0 |
| `TGeoCompositeShape` | 16 | 0 | 15 | 1 |
| mirrored copies | 37 | 0 | 37 | 0 |

100 % of the four Tier-1 classes, 0 % of everything else. There is no partial credit anywhere: the
recogniser's boundary is exactly the shape vocabulary it was built for.

### The 95 declines, by family and by the source class behind them

| n | reason family | source classes |
| ---: | --- | --- |
| **37** | free-form faces (surface kind outside plane/cylinder/cone/sphere) | 37 mirrored copies — **all of them**, and only them |
| **27** | *a box face has no opposite partner* | `TGeoTrd1` 16, `TGeoPgon` 6, `TGeoArb8` 5 |
| **14** | *N cap plane(s) perpendicular to the axis, expected 2* | `TGeoPcon` 13, composite 1 |
| 8 | *a planar face is neither a cap nor a wedge of any axis cluster* | composite 8 |
| 6 | *N planar faces: not a six-plane box* | composite 5, `TGeoPgon` 1 |
| 2 | *N axis clusters: beyond the recogniser's scope (Tier 3)* | composite 2 |
| 1 | *mixed lateral surface kinds* | `TGeoPcon` 1 |

Three things this says that PIPE could not.

**(a) The free-form class is 100 % self-inflicted.** All 37 free-form declines are the mirrored
copies of §3, and every mirrored copy declines free-form. On PIPE the free-form class was 23 genuine
elliptical cylinders and a twisted `Arb8`. Fix the mirror and this whole column disappears — and the
37 duplicate definitions of 3 shapes with it.

**(b) The largest genuine gap is `TGeoTrd1`, and it is a small job.** 16 tapered boxes decline
because the box rule wants three pairs of parallel faces. A `Trd1`/`Trd2` is six planes with one
parallel pair and four tapered sides; it is a strictly easier detector than the revolved-profile one,
and with `TGeoPgon` (6 more on the same rule) and `TGeoArb8` (5) it is worth **27 of 95 declines**.

**(c) Roadmap (h), "are we able to recognise Pcons, Polyhedra?", answered a second time and the same
way.** No. 13 of 14 `TGeoPcon`s decline on *"more than two cap planes perpendicular to the axis"* —
the identical reason that took 46 of PIPE's 106 declines. Two independent detectors, one missing
piece: the revolved-profile detector (`Workstreams.md` Stream D, `O2RevolvedSolid`), *all carriers
coaxial ⇒ revolution*. TPC's 14 Pcons include `TPC_Drift` and `TPC_M` themselves.

## 5. Is the geometry the same? The three instruments

10 000 uniform points per part in the part's own bbox and frame; **1.68 M points drawn** (168 parts),
points within 1e-6 cm of a surface skipped as unfair (`TGeoShape::Safety`; only 10 points in the
whole run were skipped). The 37 mirrored copies are excluded from the `Contains` instruments — a
reflected solid is deliberately *not* the source shape as a point set — and the 2 multibody
`TPC_ORH` pieces have no single source shape.

**A — the STEP read back, against the source `TGeoShape`.** `BRepClass3d_SolidClassifier` on the
shape OCCT reads out of our own file, against `TGeoShape::Contains`.

> **129 parts, 1 289 990 points scored, 129 with zero disagreements, 0 mismatching points.**

PIPE had one disagreeing point in 1.69 M. TPC has none in 1.29 M.

**B — the recognised CSG, against the source `TGeoShape`.** The `TGeoShape` the reverse converter
emitted into `shape_*.root`, put back in the part frame with the converter's own `shapePlacement`.

> **73 parts, 729 990 points scored, 0 disagreements.** (B re-uses A's points on the CSG subset.)

And by the converter's own acceptance measurement, **all 73 recognised CSG shapes have symmetric
difference exactly 0 cm³** against the CAD part — 73/73, max 0.000e+00.

**C — capacity**, `TGeoShape::Capacity()` of the source against `BRepGProp` on the STEP read-back,
all 167 leaf definitions:

| source shape class | n | median | max |
| --- | ---: | ---: | ---: |
| `TGeoTubeSeg` | 35 | 3.5e-13 | 1.5e-10 |
| `TGeoTube` | 21 | 1.8e-15 | 8.7e-13 |
| `TGeoTrd1` | 16 | 1.9e-14 | 7.9e-13 |
| `TGeoPcon` | 14 | 1.7e-14 | 8.3e-13 |
| `TGeoBBox` | 13 | 1.6e-16 | 4.5e-16 |
| `TGeoPgon` | 7 | 3.2e-13 | 4.8e-11 |
| `TGeoArb8` | 5 | 1.9e-14 | 1.6e-13 |
| `TGeoCone` | 4 | 2.4e-12 | 4.4e-12 |
| **all non-composite, non-mirrored** | **115** | per-class 1.6e-16 – 2.4e-12 | **1.5e-10** |
| `TGeoCompositeShape` | 15 | 5.5e-03 | 1.5e-02 |
| mirrored copies | 37 | 7.6e-03 | 2.6e-02 |

The mirrored row is §3's defect. The composite row is **the oracle's noise, not ours**, and TPC
proves it twice over.

First, an independent 2 000 000-point Monte Carlo of the *source* shape, on the eight worst
composites:

| volume | `Capacity()` | our OCCT | MC ± σ (N = 2e6) | our σ | `Capacity()` σ |
| --- | ---: | ---: | ---: | ---: | ---: |
| `TPC_OPLL` | 14.296651 | 14.572017 | 14.579269 ± 0.015475 | **0.5** | 18.3 |
| `TPC_OCGEM` | 86.347091 | 87.623316 | 87.674599 ± 0.288640 | **0.2** | 4.6 |
| `TPC_rods` | 83.404395 | 84.424880 | 84.488456 ± 0.105060 | **0.6** | 10.3 |
| `TPC_ORH` | 4.699719 | 4.751730 | 4.752123 ± 0.006918 | **0.1** | 7.6 |
| `TPC_cliholdv` | 1.256247 | 1.269469 | 1.270731 ± 0.000982 | **1.3** | 14.7 |
| `TPC_INPLL` | 14.205732 | 14.339313 | 14.330567 ± 0.016277 | **0.5** | 7.7 |
| `TPC_ORCAP` | 2721.023371 | 2746.105975 | 2740.350045 ± 5.353748 | **1.1** | 3.6 |
| `TPC_rodl` | 99.886127 | 99.276088 | 99.255530 ± 0.122948 | **0.2** | 5.1 |

Our volumes sit at 0.1–1.3 σ of the truth; `TGeoCompositeShape::Capacity()` sits at 3.6–18.3 σ.

Second, and more directly: the capacity pass calls `Capacity()` once per part, and `TPC_INPLL` has
18 mirrored copies, so it was called 18 times **for the same shape, in one process**. The six
largest of those calls came back 14.1497, 14.1540, 14.1544, 14.2077, 14.2246, 14.2503 — a different
answer every time, spread over 0.7 %. Our OCCT volume for that shape is one number, 14.339313, to
the last digit on every call.

## 6. Wall time and memory

| stage | wall | peak RSS | output |
| --- | ---: | ---: | --- |
| `o2-sim -n 0 -m TPC` | 1.3 s | 533 MB | 208 kB geometry |
| `O2_TGeoToCAD.py --no-step` | 2.3 s | 727 MB | coverage table |
| `O2_TGeoToCAD.py --dedup-world` | 7.5 s | 1.16 GB | **37.06 MB STEP**, 505 426 entities |
| `O2_CADtoTGeo.py --exact-surfaces auto --csg auto` | **2 min 40.5 s** | 866 MB | 29 MB, 579 files |
| scoring: 2 M-point MC of 8 composites | < 30 s | — | §5 table |

## 7. Findings, in the order I would act on them

1. **`BRepBuilderAPI_GTransform` loses 0.5–1.1 % of the volume and turns every analytic face into a
   B-spline.** Analytic oracle (`TPC_CDCE`, a `TGeoPcon`): 0.887 % too large. Nine-line reproducer
   in §3. The fix is `gp_Trsf::SetMirror` + `BRepBuilderAPI_Transform`, measured exact and
   carrier-preserving. This path carries no volume invariant; it should carry the same one
   `_boolean` does.
2. **A reflecting placement of a volume with daughters drops the subtree and leaves it in the file
   as an orphaned root.** 37 of them on TPC, 16 040 leaf placements — 44 % of the detector, one
   whole readout endcap. It is silent apart from a decline record, and the orphans then trip the
   reverse converter's coincidence warning 1 535 times. Mirroring a subtree means pushing the
   reflection down to the leaves, which the per-occurrence `--dedup-world` walk is already the right
   place for.
3. **`TGeoHalfSpace` costs 5 composites and 732 placements on TPC** and is the whole of the mapper's
   decline list here. The enclosing composite always bounds it.
4. **`TGeoTrd1` in the reverse recogniser: 16 declines, and `TGeoPgon`+`TGeoArb8` 11 more on the same
   rule.** Six planes, one parallel pair, four tapered sides.
5. **The revolved-profile detector**, again: 13 of 14 TPC `TGeoPcon`s, 46 of 106 PIPE declines. Two
   detectors, one missing piece.
6. **`TPC_PRSTR1` is placed twice at z = −177.925 and never at +177.925**
   (`Detector.cxx:1388-1389`). One misplaced 16.3 cm³ Tedlar strip in the inner field cage. Worth
   reporting upstream; `CheckOverlaps` cannot see it.
7. **`TPC_WSEG`'s two degenerate two-edge planar wires** are the only tessellated part in the
   detector and have a named cause.
8. **The STEP writer's size limit is between 38 676 and 74 601 components.** TPC writes; the full
   model segfaults. A bisection now has both ends.

## 8. Dead premises, recorded

- **"TPC will be `Pcon`/`Pgon`-heavy like the beam pipe."** It is `TubeSeg`-heavy (35) with only 14
  `Pcon`s, and by *placement* count it is boxes: `TPC_RCCONO` and `TPC_RCCONOB` at 8 352 each.
- **"The full-geometry writer segfault will bite here too."** It did not, at 505 426 entities in
  6.4 s. Assuming it would have cost nothing; assuming it would *not* have would have cost the
  stream.
- **"The reverse converter will refuse TPC the way it refused PIPE."** It refused PIPE on a genuine
  TGeo defect. On TPC it warned about 1 535 coincidences that our own exporter had manufactured.
  Stream AD's rule — *arbitrate on the TGeo side first* — worked in both directions: the independent
  walk said TPC's own count is 1, so the other 1 534 had to be ours.
- **"A reflection is an isometry, so baking one is safe."** Refuted by an analytic oracle at
  0.887 %. Isometry is a property of the *transform*; whether the CAD kernel's implementation of it
  is an isometry is a separate, measurable question.
- **"A composite disagreeing with `Capacity()` at 1e-2 means we are wrong."** Refuted for the second
  time, and this time without any Monte Carlo of my own: 18 calls to `Capacity()` on one unchanged
  shape returned 18 different numbers.
- **"Volume-only checking would have caught the mirror bug."** It nearly did not. The 37 mirrored
  copies sit in the report at 7.6e-03 median — *inside* the composite row's ordinary 5.5e-03 noise.
  What actually separated them was the face-carrier table (10 planes+cylinders becoming 10
  B-splines) and the one mirrored `TGeoPcon` with an analytic `Capacity()`. Stream AD's suite 6
  exists for exactly this, and it was right again.

## 9. Reproducing it

```bash
# 1. the geometry (sim env; NEVER with the pythonOCC prepends)
VMCWORKDIR=$B/stage/share $B/stage/bin/o2-sim -n 0 -g boxgen -m TPC \
        --configKeyValues "align-geom.mDetectors=none"

# 2. out (converter env = O2 env PLUS the OCC prepends, in one process)
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root TPC_dedup.step \
        --dedup-world --report TPC_dedup_report.json          # 7.5 s, 37.06 MB

# 3. back
python3 scripts/geometry/O2_CADtoTGeo.py TPC_dedup.step --output-folder rt_tpc_out -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report rt_tpc_out/surface_report.json --csg-report rt_tpc_out/csg_report.json
```

The scoring scripts are Stream AD's (`rtcompare.py`, `rtcap.py`, `dupcheck.py`, `census.py`), with
two one-line changes for TPC: PyROOT's `GetVolume` returns a *null-pointer proxy*, not `None`, for
the `__mirrored` names, so the guard must be `if not vol`; and the mirrored copies must be excluded
from the `Contains` instruments. The scripts written for this stream (`lost.py`, `account.py`,
`stepstat.py`, `mccomp.py`, `mirror.py`, `mirrorfix.py`, `dupaths.py`, `comps.py`) live in this
session's scratch and are not committed; every number they produced is quoted above with the sample
size it came from.
