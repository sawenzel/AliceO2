# Stream AF — the TGeo → STEP → TGeo round trip on ITS

`Stream_AD_TGeoToStep.md` ran the round trip on the beam pipe and left a question open: how much
of what it measured was PIPE, and how much was the method. This document runs the same
methodology, with the same scripts, on **ITS** — a detector with two hundred times PIPE's
placements, no torus and no elliptical cylinder anywhere, and the first one whose shape mix is
dominated by `TGeoXtru` and by deep boolean chains.

The point of the exercise is unchanged: the STEP is generated from a known `TGeo`, so the
**original `TGeoShape` is an exact oracle** and there is no tolerance argument to have about the
model. This is Roadmap item (i) applied a second time, and it is the Pcon/Pgon evidence Roadmap
item (h) asked for, now on a corpus that PIPE could not supply.

Measured 2026-08-23 on `swenzel/bvhsurfacesolid`, aarch64, OCCT/pythonOCC 7.9.0, ROOT 6.36.10.
Every number below was run on this machine.

---

## 1. The source geometry and its census

```bash
# the sim env only -- NEVER with the pythonOCC prepends, they segfault o2-sim
VMCWORKDIR=$B/stage/share \
  $B/stage/bin/o2-sim -n 0 -g boxgen -m ITS --configKeyValues "align-geom.mDetectors=none"
```

2.1 s, and the ideal `o2sim_geometry.root` is 802 kB. `VMCWORKDIR=$B/stage/share` is needed here
even though PIPE did not need it.

| | |
| --- | --- |
| logical volumes in the file | 409 |
| reachable from `cave` | **270** |
| DAG placement edges (`AddNode`) | 1 936 |
| fully expanded placements | **296 943** (ROOT's own "296943 nodes") |
| max depth | 10, max placements of one volume 126 |

| shape class | reachable volumes |
| --- | ---: |
| `TGeoBBox` | 75 |
| `TGeoCompositeShape` | 70 |
| `TGeoTube` | 39 |
| `TGeoShapeAssembly` | 29 |
| `TGeoXtru` | 23 |
| `TGeoPcon` | 19 |
| `TGeoTubeSeg` | 13 |
| `TGeoTrd1` / `TGeoCone` | 1 each |

**ITS is a different corpus from PIPE, and deliberately so.** Nine shape classes against PIPE's
ten, but the overlap is only partial: there is no `TGeoTorus`, no `TGeoEltu`, no `TGeoSphere` and
no free-form anything, and there are 23 `TGeoXtru` where PIPE had none. The 70 composites are
**442 subtractions and 77 unions** over operands that are `TGeoTube` (70), `TGeoBBox` (58),
`TGeoXtru` (31), `TGeoPcon` (20), `TGeoTubeSeg` (13), `TGeoArb8` (2) and `TGeoPgon` (1).

### The boolean trees are chains, not trees

Measured over all 70 composites — depth and operator count of each boolean tree:

| | |
| --- | --- |
| median depth | 2 |
| max depth | **60** (`EndWheelCBasis2`) |
| composites deeper than 32 | **5** |

For every one of the deep ones the operator count *equals* the depth — 60 ops at depth 60, 52 at
52, 50, 44, 35. They are left-linear spines of successive `TGeoSubtraction`s, one hole at a time,
not balanced trees. This matters in §3.

## 2. The coincidence walk: ITS has none

PIPE's round trip was refused by `O2_CADtoTGeo.py` because the beam pipe places 81 volumes exactly
on identical copies of themselves. The same walk on ITS, composing every world matrix over the
fully expanded tree:

```
leaf-solid placements: 261164   mother-body placements: 35557   total 296721
distinct (volume, world transform) pairs: 296721
coincident duplicates: 0 extra placements over 0 pairs
multiplicity: {1: 296721}
```

**Zero.** All 296 721 (volume, world transform) pairs are distinct, so `--dedup-world` is not
needed and was not used. Run twice with two independently written transform composers — the
`TGeoHMatrix` walk from Stream AD (`dupcheck.py`) and a numpy 4×4 rewrite — which agree to the
placement. That is the useful control: the 81 plies were a property of the beam pipe, not an
artefact of the instrument.

It is also the reason ITS could be exported at all as one file. `--dedup-world` expands per
occurrence; on ITS that would have turned 1 982 STEP components into 296 712, which is four times
the component count that already segfaults `STEPCAFControl_Writer` on the full geometry
(Stream AD §5).

## 3. TGeo → STEP: 235 of 241 shaped volumes, and a policy limit

```bash
python3 scripts/geometry/O2_TGeoToCAD.py sim/o2sim_geometry.root ITS.step --report ITS_report.json
```

| | |
| --- | --- |
| solids built | **235** |
| volumes with daughters | 84 |
| pure assemblies | 29 |
| components | **1 982** |
| declined | **6** |
| baked mirror copies / reflected placements | **0 / 0** |
| wall | **19.5 s** |
| peak RSS | **858 MB** |
| the STEP | **12.84 MB, 223 652 entities, written without incident** |

235 of the 241 shaped reachable volumes convert — **97.5 %**. The write that dies on the full Run 3
geometry completes here in seconds, which localises Stream AD's open item 2: it is not
`STEPCAFControl_Writer` per se, it is scale, and 1 982 components is comfortably inside it.

### The six declines, and what they say

| n | volume(s) | reason as reported |
| ---: | --- | --- |
| 5 | `EndWheelCBasis0/1/2`, `IBCYSSFlangeA`, `IBCYSSFlangeC` | *TGeoCompositeShape: boolean tree deeper than 32* |
| 1 | `IBGammaConvWireOuterSupport` | *TGeoArb8: sections have [2, 4] distinct vertices; a prism needs the same count in every section* |

The five are exactly the five composites §1 measured at depth 35–60, and the guard they hit is a
**bare constant in the mapper** — `if depth > 32: raise ShapeDeclined(...)` at
`O2_TGeoToCAD.py:829`, with no CLI knob and no OCCT limit behind it. It was never exercised before
because PIPE's deepest composite is nowhere near it. Since the trees are linear chains of cuts
rather than nested trees, raising the constant is the whole fix, and the volume invariant already
in `_boolean` is what would have to price the result. This is recorded as a finding, not fixed:
no source-code edits were made for this study.

The sixth is a **degenerate `TGeoArb8`** — one of its two z sections has collapsed from four
distinct vertices to two, i.e. it is a wedge that ends on an edge. The prism sewer refuses it
because it wants the same vertex count in both sections. That is a real geometry idiom (a knife
edge), and it is a second cheap win.

**Nine placements are lost**, out of 296 721: the six volumes are placed 1, 1, 1, 2, 2, 2 times.
The reverse converter independently counts 296 712 leaf placements, and 296 721 − 9 = 296 712, so
the bookkeeping closes exactly across the two programs.

### Capacity of what was built, against `TGeoShape::Capacity()`

| source shape class | n | median | max |
| --- | ---: | ---: | ---: |
| `TGeoBBox` | 75 | 0.0 | 4.3e-16 |
| `TGeoTube` | 39 | 4.6e-16 | 2.2e-14 |
| `TGeoXtru` | 23 | 2.0e-16 | 3.8e-15 |
| `TGeoPcon` | 19 | 2.7e-15 | 4.0e-14 |
| `TGeoTubeSeg` | 13 | 2.3e-15 | 4.7e-14 |
| `TGeoCone` / `TGeoTrd1` | 1 each | 0.0 | 0.0 |
| **all non-composite** | **171** | **1.9e-16** | **4.7e-14** |
| `TGeoCompositeShape` | 64 | 5e-03 | 1.7e-02 |

The non-composite column is the double-precision floor, and the 23 `TGeoXtru` sit in it — the
face-by-face prism sewer (Stream AD §2) reproduces an ITS extruded profile exactly.

**The composite column is the oracle's noise, and ITS prices it harder than PIPE did.**
`TGeoCompositeShape::Capacity()` is itself a Monte Carlo. Scoring the six worst against an
independent 2 000 000-point MC of the *source* shape:

| volume | `Capacity()` | our OCCT volume | MC ± σ (N = 2e6) | our σ | `Capacity()` σ |
| --- | ---: | ---: | ---: | ---: | ---: |
| `CageCoverRib` | 249.276 | 245.097 | 244.683 ± 0.940 | **0.4** | 4.9 |
| `ITSServicesWaterIB` | 943.849 | 957.860 | 986.283 ± 28.27 | **1.0** | 1.5 |
| `CageSidePanelGuide` | 7862.82 | 7747.92 | 7747.30 ± 10.10 | **0.1** | 11.4 |
| `ITSServicesCopperOB` | 196099 | 193749 | 193515 ± 349.3 | **0.7** | 7.4 |
| `SpaceFrameVolumeLay3` | 997.111 | 985.438 | 984.784 ± 0.617 | **1.1** | **20.0** |
| `BPSupportClamp` | 5.53704 | 5.60000 | 5.60092 ± 0.0055 | **0.2** | 11.6 |

Our volumes sit at 0.1–1.1 σ of the truth; `TGeo`'s own `Capacity()` sits at up to **20 σ**. The
1e-2 column is the estimator, not the geometry, and this is the second detector to say so.

## 4. STEP → TGeo: the tier table, and the first zero

```bash
python3 scripts/geometry/O2_CADtoTGeo.py ITS.step --output-folder rt -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report rt/surface_report.json --csg-report rt/csg_report.json
```

92.6 s, peak RSS 961 MB.

```
Placement check: 296712 leaf placement(s), all at distinct world transforms.
Surface report: 235/235 logical volumes eligible for exact O2BVHSurfaceSolid conversion
Exact-surface extraction (auto): 235/235 leaf solids represented exactly, 0 fall back
  tiers: CSG 128, exact surfaces 107, tessellated 0  (of 235 leaf solids)
```

| tier | parts | share |
| --- | ---: | ---: |
| CSG (native `TGeoShape`) | **128** | 54.5 % |
| exact surfaces (`O2BVHSurfaceSolid`) | **107** | 45.5 % |
| tessellated | **0** | 0 % |

**Zero tessellated parts.** PIPE had 23. The reason is visible in the face census of the whole
file, as the reverse converter classifies it:

| plane | cylinder | cone | torus | bspline | extrusion |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 292 | 459 | 27 | 0 | 0 | 0 |

2 778 faces, **every single one a quadric**. PIPE's 23 tessellated parts were 22 `TGeoEltu` and one
twisted `TGeoArb8`; ITS contains neither, so the surface tier covers the file completely and this
is the first detector to round-trip with no mesh anywhere. Worth separating carefully: the
*carriers* are all quadrics, but 360 of the 12 650 trim **curves** are B-splines (plus 32
ellipses) — the exact-trim machinery handles those, and they cost nothing at the tier level.

## 5. The CSG study — what is recognised, and why not more

**128 of 235 parts accepted as CSG, 54.5 %.** Broken down by the *source* `TGeo` class each part
came from — which is the whole value of a round trip against a known oracle:

| source class | parts | CSG | rate |
| --- | ---: | ---: | ---: |
| `TGeoBBox` | 75 | 75 | **100 %** |
| `TGeoTube` | 39 | 39 | **100 %** |
| `TGeoTubeSeg` | 13 | 13 | **100 %** |
| `TGeoCone` | 1 | 1 | **100 %** |
| `TGeoCompositeShape` | 64 | 0 | 0 % |
| `TGeoXtru` | 23 | 0 | **0 %** |
| `TGeoPcon` | 19 | 0 | **0 %** |
| `TGeoTrd1` | 1 | 0 | 0 % |

Every quadric primitive comes back, all 128 of them, and every one at symmetric difference 0 by
the converter's own acceptance measurement. Nothing else comes back at all.

### The decline histogram, 107 parts

| n | reason |
| ---: | --- |
| **47** | *N planar faces: not a six-plane box* |
| 19 | *a planar face is neither a cap nor a wedge of any axis cluster* |
| **18** | *N cap plane(s) perpendicular to the axis, expected 2* |
| 12 | *N axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built)* |
| 5 | *mixed lateral surface kinds `['cone','cylinder']` on one axis* |
| 4 | *symmetric difference exceeds the band* |
| 2 | *a box face has no opposite partner* |

Crossed with the source class, the histogram stops being a list of recogniser messages and becomes
a statement about ALICE geometry:

| n | source class | reason |
| ---: | --- | --- |
| 25 | `TGeoCompositeShape` | not a six-plane box |
| **22** | `TGeoXtru` | **not a six-plane box** |
| 19 | `TGeoCompositeShape` | planar face is neither cap nor wedge |
| **18** | `TGeoPcon` | **N cap planes, expected 2** |
| 12 | `TGeoCompositeShape` | 3–14 axis clusters, Tier 3 |
| 4 | `TGeoCompositeShape` | mixed lateral kinds on one axis |
| 4 | `TGeoCompositeShape` | symmetric difference exceeds the band |
| 1 | `TGeoXtru` | a box face has no opposite partner |
| 1 | `TGeoPcon` | mixed lateral kinds on one axis |
| 1 | `TGeoTrd1` | a box face has no opposite partner |

### This is the answer to Roadmap (h), and ITS answers the half PIPE could not

**Pcons: no, and for exactly the reason PIPE measured.** All 19 `TGeoPcon` parts decline, 18 of
them as *"N cap planes perpendicular to the axis, expected 2"* with N distributed
{1: 1, 3: 6, 4: 6, 5: 3, 6: 1, 8: 1} — broadly one cap per z plane of the profile, which is what a
Pcon is: `OBEndWheelCLowerRing5` reaches 8 caps over 22 faces, and the one part with a *single*
cap, `IBCYSSConeFoam`, is an 11-face cone/cylinder stack where only one of its three planes is
perpendicular to the axis at all.
Every individual face is recognised (they are in the cylinder and cone columns of §4's table); it
is the whole-part rule that expects one axis cluster bounded by two caps. PIPE said this with 58
Pcons and one detector; ITS says it with 19 Pcons and a second, independent one. The missing piece
is the same: the revolved-profile detector, *all carriers coaxial ⇒ revolution*
(`Workstreams.md` Stream D, `O2RevolvedSolid`).

The nineteenth Pcon, `IBCYSSCone`, declines differently and usefully: *mixed lateral surface kinds
`['cone','cylinder']` on one axis, 14 faces: 4 cone, 6 cylinder, 4 plane, 1 axis cluster*. It has
both conical and cylindrical bands, so even a two-cap Pcon would not have been recognised. A
revolved-profile detector must accept a **mixed** cone/cylinder lateral set on one axis, not just a
homogeneous one — PIPE could not have shown that, and it is a concrete requirement on Stream D.

**Polyhedra: no, and this is the half of (h) that only ITS could measure.** PIPE had no `TGeoXtru`
and one `TGeoPgon` operand. ITS has 23 `TGeoXtru` volumes, and **22 of them decline as
*"N planar faces: not a six-plane box"***, N distributed over
{5, 7, 8, 10, 11, 16, 18, 20, 22, 24, 25, 26, 30, 33, 50, 54, 58, 74} — 10 planar faces is the
mode (9 parts), 74 the largest (4 parts). These are all-planar prisms and the recogniser's only
all-planar rule is the six-plane box. The gap is therefore precisely a **prism recogniser**: a
polygonal profile extruded along one axis, which is Tier-3 prism territory in `Stream_AA_FlatCSG.md`
and would also take the 25 all-planar composites in the same class. Together that is **47 of the
107 declines, the single largest class**.

The 23rd Xtru and the sole `TGeoTrd1` decline as *"a box face has no opposite partner"*: they are
six-plane solids whose tapered sides are not parallel. A `TGeoTrd1` is a stock TGeo primitive and
comes back as six planes with the right topology — the box rule just insists on three pairs of
parallel faces. That is the cheapest gap on the list.

### The four the recogniser refused itself

`ITSServicesWaterOB`, `ITSServicesPolymerOB`, `ITSServicesCopperOB` and `ITSServicesCarbonOB`
produced a *candidate* that the acceptance measurement then rejected, e.g.

```
symmetric difference 2345.92 cm^3 exceeds the band 0.00954243 cm^3
(= 1e-07 cm x area 95424.3 cm^2); extra 2345.92, missing 0
```

`extra 2345.92, missing 0` says the candidate strictly contains the part: the recogniser proposed
the enclosing barrel of what is really a barrel with material removed. On a 196 099 cm³ part that
is a 1.2 % overstatement of the volume, and it would have been shipped silently by a recogniser
that trusted its own pattern match. It was not, because acceptance is measured against the CAD
solid and not asserted. This is the self-checking-parity principle earning its keep on a detector
it was not designed against, and one more entry on this branch's list of things that
`IsDone()`-style trust would have gotten wrong.

## 6. Is the geometry the same? The scoring

Two instruments from Stream AD, reused unchanged (`rtcompare.py`), sampling **8 000 uniform points
per part** in the part's own bbox and frame — 235 parts, so a **1.88 M point budget**, inside the
2 M cap set for this study. Points within 1e-6 cm of a surface are skipped as unfair
(`TGeoShape::Safety`).

### A — the STEP read back, against the source `TGeoShape`

`BRepClass3d_SolidClassifier` on the shape OCCT reads out of our own STEP file, against
`TGeoShape::Contains` on the original volume. 235 parts, 1 880 000 points drawn, 270 skipped as
within 1e-6 cm of a surface, **1 879 730 scored**:

> **235 parts with zero disagreements. Not one mismatching point in 1.88 million.**

PIPE's number was 168 of 169 parts clean, with one point in 10 000 on a composite of tori. ITS has
no tori and no free-form carriers, and the instrument returns a clean sweep.

### B — the recognised CSG, against the source `TGeoShape`

For the 128 parts accepted as CSG, the `TGeoShape` the converter emitted into `shape_*.root`, put
back into the part frame with the converter's own `shapePlacement`, against the original shape:
**1 023 760 points over 128 parts — the same sample set, less the on-surface skips —
and 0 disagreements.** That is the complete
loop — C++ `TGeo` → STEP → recognised `TGeo` — closed in `TGeo`'s own terms.

Separately, all 128 have **symmetric difference exactly 0 cm³** against the CAD part, by the
converter's own acceptance measurement.

### C — capacity through the write and the read

`TGeoShape::Capacity()` of the source against `BRepGProp` on the STEP *read-back* — so this one
prices the file format as well as the mapping:

| source shape class | n | median | max |
| --- | ---: | ---: | ---: |
| `TGeoBBox` | 75 | 1.3e-16 | 3.2e-13 |
| `TGeoCompositeShape` | 64 | 4.9e-03 | 1.6e-02 |
| `TGeoTube` | 39 | 2.7e-15 | 2.6e-14 |
| `TGeoXtru` | 23 | 2.6e-14 | 5.6e-12 |
| `TGeoPcon` | 19 | 3.2e-14 | 5.5e-12 |
| `TGeoTubeSeg` | 13 | 4.4e-14 | 5.3e-13 |
| `TGeoTrd1` / `TGeoCone` | 1 each | 1.3e-14 / 3.1e-14 | — |
| **all non-composite** | **171** | **1.2e-15** | **5.6e-12** |
| all 235 | | 1.7e-14 | 1.6e-02 |

Compare §3: in memory the non-composite median was 1.9e-16, and after a STEP write and read it is
1.2e-15, with the tail moving from 4.7e-14 to 5.6e-12 (a `TGeoXtru`). **About an order of magnitude
in the median, and two in the tail, are paid to the file rather than to the geometry** — AP214
writes decimal text, and that is what survives it. Far below any physical tolerance either way, but
worth knowing which stage costs it.

The composite row is unchanged from §3 and is the Monte-Carlo oracle, as the 2 M-point cross-check
there establishes.

### The verdict

The round trip is exact on ITS. Every one of the 235 parts is the shape the TGeo author wrote, to
1.88 million classifier points and to the double-precision floor in volume; the 107 parts that are
not CSG are not *wrong*, they are carried by `O2BVHSurfaceSolid` exactly, and nothing at all was
approximated by a mesh. The only geometry lost anywhere in the pipeline is the 6 volumes / 9
placements the mapper declined in §3.

## 7. Wall, memory, and where the scale actually bites

| step | wall | peak RSS | output |
| --- | ---: | ---: | --- |
| `o2-sim -n 0 -m ITS` | 2.1 s | — | 802 kB |
| census (pyROOT) | 2.5 s | — | — |
| coincidence walk, 296 721 placements | 5.5 s | — | — |
| `O2_TGeoToCAD.py` | **19.5 s** | **858 MB** | 12.84 MB STEP |
| `O2_CADtoTGeo.py --exact-surfaces auto --csg auto` | **92.6 s** | **961 MB** | 8.0 MB: `geom.C` 1.37 MB, 235 surface sidecars, 128 `shape_*.root` |
| scoring, 1.88 M classifier points | 6 min 28 s | 757 MB | — |

**The 296 943-placement scale is not where the cost is.** It costs 1 982 STEP components and one
`Placement check` line; everything expensive scales with the 235 *definitions* and their face
counts. That is the DAG economy both converters were built around, and ITS is the first detector
large enough to demonstrate it rather than assert it.

## 8. Dead premises, recorded

- **"ITS is large, so the round trip will break on scale."** It is large in *placements* —
  296 943, 200× PIPE — and that turned out to be free. What actually broke was a **depth-32
  constant** meeting a 60-deep chain of cuts, which has nothing to do with size.
- **"A detector with many staves will be Arb8-heavy."** ITS has two `TGeoArb8` in the whole
  geometry, both as composite operands, and one of them is degenerate. The polyhedral load is
  `TGeoXtru`, 23 of them, which the census said in 2.5 seconds and the expectation did not.
- **"The composites' 1e-2 capacity spread is a converter problem."** Same refutation as PIPE, with
  a sharper number: on `SpaceFrameVolumeLay3` our volume is 1.1 σ from an independent 2 M-point MC
  and `TGeoCompositeShape::Capacity()` is 20 σ from it.
- **"The reverse converter's tessellated tier is unavoidable at detector scale."** On ITS it is
  empty. Free-form output is a property of specific TGeo classes — `TGeoEltu`, twisted `TGeoArb8` —
  not of the pipeline.
- **`--dedup-world` was *not* used**, and should not be reached for by reflex: on a geometry with
  no coincidences it would have multiplied the component count by 150 for nothing.

## 9. What is open

1. **The depth-32 guard** (`O2_TGeoToCAD.py:829`). Five ITS volumes, chains not trees. Not fixed
   here.
2. **Degenerate `TGeoArb8` sections** — a wedge ending on an edge, one volume.
3. **A prism recogniser.** 47 of 107 declines, of which 22 are `TGeoXtru` with a measured
   planar-face distribution to build against.
4. **The revolved-profile detector**, with a new requirement from `IBCYSSCone`: mixed cone and
   cylinder lateral faces on one axis must be accepted.
5. **`TGeoTrd1` and the box rule** — six planes, no parallel partners. One line of recogniser.
6. The `geom.C` this study produced was **not loaded into ROOT or transported through Geant**;
   scoring is against the shapes, as in Stream AD. NEXT item 1's JIT-namespace bug is unaddressed
   and would bite anyone who tried.

## 10. Reproducing it

```bash
# 1. the geometry (sim env; NEVER with the pythonOCC prepends)
VMCWORKDIR=$B/stage/share \
  $B/stage/bin/o2-sim -n 0 -g boxgen -m ITS --configKeyValues "align-geom.mDetectors=none"

# 2/3. converter env = the O2 env PLUS the OCC prepends, in one process
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root ITS.step --report ITS_report.json

# 4. back again
python3 scripts/geometry/O2_CADtoTGeo.py ITS.step --output-folder rt -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report rt/surface_report.json --csg-report rt/csg_report.json
```

The census, coincidence-walk and scoring scripts are Stream AD's (`census.py`, `dupcheck.py`,
`rtcompare.py`), unchanged, plus a numpy rewrite of the coincidence walk and a boolean-depth
counter written for this stream; they live in this session's scratch and are not committed. Every
number above is quoted with the sample size it came from.
