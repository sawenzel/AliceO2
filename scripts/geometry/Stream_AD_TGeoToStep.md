# Stream AD — TGeo → STEP, and the round trip back

The ALICE Run 3 geometry is written in C++ and only ever exists as `TGeo`.
`O2_TGeoToCAD.py` writes it out as a STEP file with its assembly tree intact, which makes it
CAD-editable and — the reason it lives on this branch — makes the **original TGeo an exact
oracle** for `O2_CADtoTGeo.py`: the STEP was generated from a known geometry, so the right answer
is known by construction and there is no tolerance argument to have about the model.

This is Roadmap item (i), and it delivers the test bench item (h) asked for.

Measured 2026-08-22/23 on `swenzel/bvhsurfacesolid`, aarch64, OCCT/pythonOCC 7.9.0, ROOT 6.36.10.
Every number below was run on this machine.

---

## 1. The source geometry

```bash
# the sim env only -- NEVER with the pythonOCC prepends, they segfault o2-sim
$B/stage/bin/o2-sim -n 0 -g boxgen -m PIPE --configKeyValues "align-geom.mDetectors=none"
```

`-n 0` writes the geometry without transporting anything (1.7 s for PIPE). o2-sim writes **both**
`o2sim_geometry.root` (ideal) and `o2sim_geometry-aligned.root`; everything here uses the **ideal**
file.

**Environment trap, new.** For the *full* geometry `VMCWORKDIR=$HOME/alisw/O2` is not enough:
`Detectors/Geometry/MFT/data/Geometry.xml` does not exist in the source tree, only in the install,
and o2-sim aborts with `Could not parse Geometry XML File`. Use
`VMCWORKDIR=$B/stage/share`. PIPE alone does not need it.

### PIPE census

| | |
| --- | --- |
| logical volumes in the file | 208 |
| reachable from `cave` | **198** |
| DAG placement edges (`AddNode`) | 333 |
| fully expanded placements | **1 475** (ROOT's own "1475 nodes") |
| max depth | 6 |

| shape class | reachable volumes |
| --- | --- |
| `TGeoPcon` | 58 |
| `TGeoTube` | 55 |
| `TGeoShapeAssembly` | 28 |
| `TGeoCompositeShape` | 25 |
| `TGeoEltu` | 22 |
| `TGeoBBox` | 6 |
| `TGeoTrd1` / `TGeoCone` / `TGeoTorus` / `TGeoTubeSeg` | 1 each |

The 25 composites are 21 subtractions, 14 intersections and 6 unions over operands that are
`TGeoTube` (38), `TGeoTorus` (20), `TGeoBBox` (9), `TGeoTrd1` (2), `TGeoArb8` (2), `TGeoPgon` (1).

### Full Run 3 census

| | |
| --- | --- |
| logical volumes in the file | 46 132 |
| reachable from `cave` | **27 351** |
| DAG placement edges | 73 858 |
| fully expanded placements | **5 465 306** |
| max depth | 13, max placements of one volume 735 |

| shape class | reachable | mapped? |
| --- | ---: | --- |
| `TGeoBBox` | 17 823 | yes |
| `TGeoShapeAssembly` | 6 134 | assembly, no solid |
| `TGeoCompositeShape` | 1 924 | yes, recursively |
| `TGeoTube` | 623 | yes |
| `TGeoTorus` | 298 | yes |
| `TGeoPcon` | 124 | yes |
| `TGeoTubeSeg` | 96 | yes |
| `TGeoTrd1` | 86 | yes |
| `TGeoTrap` | 78 | yes |
| `TGeoCone` | 44 | yes |
| `TGeoXtru` | 37 | yes |
| `TGeoPgon` | 29 | yes |
| `TGeoEltu` | 26 | yes |
| `TGeoPara` | 13 | **no** |
| `TGeoCtub` | 8 | yes |
| `TGeoArb8` | 7 | yes |
| `TGeoConeSeg` | 1 | yes |

Seventeen shape classes in the whole Run 3 geometry, and after this stream **sixteen of them are
mapped**. `TGeoPara` — 13 volumes, 0.05 % — is the only one that is not, and it is eight vertices
like `TGeoArb8`, so it is a short job rather than a research one.

## 2. The mapper

`scripts/geometry/O2_TGeoToCAD.py`. pyROOT walks the `TGeoVolume` DAG, pythonOCC builds an OCCT
solid per volume, and one XCAF document carries the tree:

```
TGeoVolume, no daughters      -> one XCAF simple shape (a definition)
TGeoVolume, with daughters    -> one XCAF assembly label,
                                 one component per TGeoNode with its TGeoMatrix
   ... and its own solid      -> plus one component `<name>__body` at the identity
TGeoVolumeAssembly            -> a pure assembly, no solid
```

A volume is converted **once** and referenced from every node that places it, so a shared volume
costs one solid — the same economy TGeo gets from its DAG. Lengths and translations are scaled
cm → mm (×10); the file is AP214 with `write.step.unit = MM`.

### Three construction decisions that turned out to matter

**Solids of revolution are revolved, not fused.** `TGeoTube`/`TubeSeg`/`Cone`/`ConeSeg`/`Pcon` are
built by revolving their r–z profile in one `BRepPrimAPI_MakeRevol`. `rmin > 0` then comes out as a
real inner face and an annular cap, rather than a boolean cut of two coaxial cylinders. The parked
`attic/O2_TGeoToCAD_carved*.py` fused per-band cones and had the inner cut commented out; this
route has no bands and no booleans.

**Polyhedral shapes are sewn from explicit faces.** `Trd1`, `Trd2`, `Arb8`, `Trap`, `Xtru`, `Pgon`
were first built with `BRepOffsetAPI_ThruSections`, which is correct in volume and **wrong in
carrier**: it writes a B-spline surface even for a flat section. Measured on the beam pipe, the
reverse converter then said *"free-form faces: 2 of 10"* for the `TRD` volume, *"8 of 10"* for
`ARB8` and *"2 of 6"* for `TDRPlate`. They are now sewn from faces built one at a time — planar
where the four corners are coplanar, ruled (`brepfill.Face`) only where they are not. That moved
**8 faces from bspline to plane** on PIPE, and the only 4 free-form faces left in the whole file
are the genuinely twisted sides of one `TGeoArb8`.

**Every boolean is priced against a volume invariant.** On the MFT part `MF21`,
`BRepAlgoAPI_Fuse` returned `IsDone() == True` and a valid solid that was **just the second
operand** — a 4.38 cm³ operand vanished, with nothing reporting an error, and the enclosing volume
came out 25 % wrong. `_boolean` therefore checks what the result must satisfy
(`fuse: max(vA,vB) ≤ v ≤ vA+vB`, `cut: vA−vB ≤ v ≤ vA`, `common: 0 ≤ v ≤ min(vA,vB)`, at 1e-4
relative slack), retries once with the operands run through `ShapeUpgrade_UnifySameDomain` — the
cause is seam faces accumulated by a chain of fuses — and declines with the numbers if it still
fails. Measured effect on the eight affected MFT volumes, against a 1 000 000-point Monte Carlo of
the source shape:

| volume | before the fix | after | MC ± σ |
| --- | ---: | ---: | ---: |
| `MF21` | 19.85882 (443 σ) | **16.11139 (1.7 σ)** | 16.13114 ± 0.01194 |
| `MF31` | — | **23.47602 (0.1 σ)** | 23.47476 ± 0.01555 |

and in the aggregate, the worst composite deviation over the whole Run 3 geometry fell from
**2.49e-1 to 2.82e-2** — the latter being the Monte-Carlo oracle's own noise (§4.3).

**Two OCCT primitives do not mean what TGeo means.**

- `BRepPrimAPI_MakeSphere(ax, R, lat1, lat2, dphi)` cuts theta with **planes** — a spherical zone.
  TGeo cuts it with **cones through the centre** — a spherical cone. For `TGeoSphere(0, 2, 30, 120)`
  that is 27.843367 cm³ against the right answer 22.887975 cm³, a 21.7 % error that no volume band
  would have forgiven. A theta-sectioned sphere is therefore revolved from an arc profile.
- `TGeoCtub`'s cut planes **replace** the ±dz end faces rather than being intersected with them, so
  its z extent follows the planes. Clipping the tube at ±dz first and then cutting gave 5.283185 cm³
  where the answer is 6.283185 cm³. The base tube is now built long enough to reach past both planes.

### Coverage

Mapped: `TGeoBBox`, `TGeoTube`, `TGeoTubeSeg`, `TGeoCtub`, `TGeoCone`, `TGeoConeSeg`, `TGeoPcon`,
`TGeoPgon`, `TGeoSphere`, `TGeoTorus`, `TGeoEltu`, `TGeoTrd1`, `TGeoTrd2`, `TGeoArb8`, `TGeoTrap`,
`TGeoXtru`, `TGeoScaledShape`, and `TGeoCompositeShape` recursively (`BRepAlgoAPI_Fuse` / `Cut` /
`Common`, with the boolean node's left and right matrices applied).

Declined, with the reason that lands in the JSON report: `TGeoHalfSpace` (unbounded),
`TGeoGtra` (the lateral twist is not a ruled loft of the eight Arb8 vertices), `TGeoParaboloid`,
`TGeoHype`, `TGeoPara`, `TGeoTessellated` (already a mesh).

A reflecting or non-uniformly scaling `TGeoMatrix` cannot be a STEP component location, because
those are rigid by construction. Such a placement is **baked** into a private mirrored copy of the
child solid, counted in the report. PIPE has none; the full Run 3 geometry has **9 309**, of which
**9 211 are baked** and 98 are not, because a reflected placement of a volume that *has daughters*
has no single solid to bake — those 98 placements are dropped and reported. Fixing that means
mirroring a whole subtree, and it is not done.

### Options that exist because the beam pipe needed them

- `--dedup-world` — expand per occurrence and drop repeated (volume, world transform) pairs; see §4.
- `--no-step` — build every solid and write the report without writing the STEP; see §5.
- `--carve-mothers` — subtract the placed daughters from each mother, for a CAD-facing export.
  **Off by default**: TGeo's implicit-subtraction semantics are deliberately not applied, so every
  part in the file is the shape the TGeo author wrote and per-part capacity stays comparable
  against `TGeoShape::Capacity()`.
- `--no-mother-bodies`, `--skip-top-body`, `--top`, `--max-depth`, `--include-name`.

## 3. The self-test — 71 checks, 0 failures

```
$ python3 scripts/geometry/O2_TGeoToCAD.py --self-test
--- primitives vs TGeoShape::Capacity(), band 1e-9: 28 checks, 0 failures
--- composites vs an independent MC of TGeo (N=200k, 4 sigma): 7 checks, 0 failures
--- Contains agreement, TGeo vs OCCT classifier (3000 pts each): 3 checks, 0 failures
--- placement transforms: 5 checks, 0 failures
--- negative controls (a wrong parameter must be rejected): 4 checks, 0 failures
--- analytic surface types of every face: 16 checks, 0 failures
--- XCAF assembly document, written and read back: 8 checks, 0 failures

71 checks, 0 failures
```

Never quote the last line alone; the six suites are 28 + 7 + 3 + 5 + 4 + 16 + 8.

**Suite 1** compares `BRepGProp` on our solid against the analytic `TGeoShape::Capacity()` at
1e-9 relative, over every mapped primitive including `rmin > 0` tubes and cones, phi wedges, theta
and phi sectioned spheres, hollow tori, a hollow `Pgon` stack, and a `Pcon` with a zero-thickness
radius jump. Typical achieved deviation is **1e-16**, i.e. the double-precision floor.

**Suite 2** cannot use `Capacity()` — `TGeoCompositeShape` estimates it by Monte Carlo, and so does
`TGeoCtub`'s (inherited) one. It runs its own MC with a stated N = 200 000, quotes the sigma, and
requires the OCCT volume within 4 sigma; the achieved separations are 0.02–1.30 sigma. The control
on the control prices a deliberate +2 % error at 8.1 sigma, so the test could have failed.

**Suite 3** classifies 3 000 random points per shape with `BRepClass3d_SolidClassifier` against
`TGeoShape::Contains`, skipping points within 1e-6 cm of a surface (`TGeoShape::Safety`). Zero
disagreements on a subtraction, an intersection, and a `Pcon` with `rmin > 0`.

**Suite 4** checks `TGeoMatrix::LocalToMaster` against the `gp_Trsf` we build from it on 200 random
points each, worst |Δ| = 2.1e-14 mm, and that a reflecting rotation is *refused* as a rigid
placement and survives as a baked general transform.

**Suite 5** is the negative-control suite: a wrong `rmin`, a wrong `Pcon` inner radius (0.8 → 0.9)
and a wrong `Pgon` edge count (6 → 7) must each be **rejected** by the 1e-9 band, at 2.5e-1, 1.3e-2
and 2.8e-2 relative. Plus one control on the control — the right answer is accepted.

**Suite 6** asserts the analytic surface type of every face of every mapped primitive, because a
volume that is right on a carrier that is wrong is exactly the defect the ThruSections lofts had.
It is the suite that would catch a regression back to B-splines.

**Suite 7** writes a small in-memory geometry to STEP and reads it back: one free shape, names
preserved, a shared volume emitted once and placed twice, the mother body present as
`world__body`, leaf occurrences equal to placements.

`--self-test` leaves via `os._exit` after printing. PyROOT double-frees the loose `TGeoShape`s at
interpreter teardown (`gGeoManager` owns them too) and aborts the process *after* every check has
run; leaving early keeps the exit status meaningful.

## 4. The PIPE round trip — the acceptance

### 4.1 The write, and a defect in the source geometry

The straightforward DAG export is 170 solids / 78 volumes with daughters / 28 pure assemblies /
383 components, 3.47 MB, 3.5 s. `O2_CADtoTGeo.py` **refused** it:

```
RuntimeError: placement invariant violated: 1170 leaf placements hold only 1146 distinct
(definition, world transform) pairs (multiplicity {1: 1122, 2: 24})
```

That is a finding, and the cheap arbiter puts it on the **TGeo** side, before any STEP exists. An
independent walk of the source geometry composing every world matrix:

```
leaf-solid placements: 1152   mother-body placements: 71   total 1223
distinct (volume, world transform) pairs: 1142
coincident duplicates: 81 extra placements, multiplicity {1: 1061, 2: 81}
```

**The ALICE Run 3 beam pipe places 81 volumes exactly on top of an identical copy of themselves.**
The mechanism is visible in the paths:

```
/barrel_1/BeamPipeCsideSection_1/bellows1BellowUS_1/bellows1Wiggle_1/bellows1LowerPlie_1
/barrel_1/BeamPipeCsideSection_1/bellows1BellowUS_1/bellows1Wiggle_2/bellows1LowerPlie_1
```

adjacent bellows convolutions each carry a "plie" disc at their shared edge, so consecutive
`Wiggle` placements put two identical discs in the same place. By volume:
`RB26s2PlieConn1` ×28, `RB26s2LowerPlie` ×15, `RB26s2UpperPlie` ×14, and 6 each of
`bellows1LowerPlie`, `bellows1vacLowerPlie`, `bellows2LowerPlie`, `bellows2vacLowerPlie`.
`O2_CADtoTGeo.py` is right to refuse: its own dedup is keyed by definition and cannot drop these,
because the same DAG edge is needed on one path and redundant on another.

`--dedup-world` was added for it: expand per occurrence (leaf definitions stay shared, so there is
still one prototype per `TGeoVolume`), drop every repeated (volume, world transform), and count
them. On PIPE that is exactly **81 dropped**, matching the independent walk, and gives
1 450 components / **4.36 MB** / 3.5 s.

### 4.2 What came back

`O2_CADtoTGeo.py --exact-surfaces auto --csg auto`, no `--mesh`, 40 s:

```
Placement check: 1146 leaf placement(s), all at distinct world transforms.
Exact-surface extraction (auto): 149/172 leaf solids represented exactly, 23 fall back
```

| | source | round trip |
| --- | ---: | ---: |
| leaf solid definitions | 170 | **172** |
| leaf placements | 1 223 − 81 = 1 142 | **1 146** |
| tiers | — | **CSG 62 / exact surface 87 / tessellated 23** |

The 170 → 172 is one part: `BPS_supportBarCarbon` is a `TGeoCompositeShape` whose OCCT result is a
compound of **three** disjoint solid bodies, and the reverse converter's multibody-leaf rule
(Stream AC) correctly splits it into three volumes. 170 − 1 + 3 = 172, and 1 142 + 4 = 1 146.

**The faces in the file, as the reverse converter classifies them:**

| plane | cylinder | cone | torus | bspline | extrusion |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 504 | 445 | 80 | 38 | 4 | 22 |

1 093 faces, of which **26 are not quadrics — and all 26 are faces TGeo itself does not define as
quadrics**: 22 elliptical-cylinder sides of the `TGeoEltu` volumes, and the 4 twisted sides of one
`TGeoArb8`. Those same 26 are the entire reason for the 23 tessellated parts
(`why_not_surface`: `extrusion face extraction not implemented yet` ×22,
`bspline face extraction not implemented yet` ×1).

### 4.3 Is the geometry the same? Three instruments

**A — the STEP read back, against the source `TGeoShape`.** 10 000 uniform points per part in the
part's own bbox and frame, `BRepClass3d_SolidClassifier` on the shape OCCT reads out of our file
against `TGeoShape::Contains`, points within 1e-6 cm of a surface skipped as unfair
(`TGeoShape::Safety`). 169 parts (the 3 multibody ones are excluded — they have no single source
shape), **1 689 997 points scored**:

> **168 parts with zero disagreements; one part with one disagreement.**
> `bellows2vacLowerPlie`, 1 point in 10 000, a `TGeoCompositeShape` of tori.

**B — the recognised CSG, against the source `TGeoShape`.** For the 62 parts `O2_CADtoTGeo.py`
accepted as CSG, the `TGeoShape` it emitted into `shape_*.root`, put back into the part frame with
the converter's own `shapePlacement`, against the original shape. 62 parts × 10 000 points =
**620 000 points, 0 disagreements**. This is the complete loop scored in TGeo terms.

**C — capacity.** `TGeoShape::Capacity()` of the source against `BRepGProp` on the STEP read-back,
all 170 parts:

| source shape class | n | median | max |
| --- | ---: | ---: | ---: |
| `TGeoPcon` | 58 | 1.0e-13 | 1.8e-08 |
| `TGeoTube` | 55 | 1.6e-13 | 4.2e-12 |
| `TGeoEltu` | 22 | 6.6e-15 | 1.1e-14 |
| `TGeoBBox` | 6 | 1.4e-13 | 1.6e-13 |
| `TGeoCone` / `Torus` / `TubeSeg` / `Trd1` | 1 each | — | 3.4e-12 |
| **all non-composite** | **145** | **1.5e-13** | **1.8e-08** |
| `TGeoCompositeShape` | 25 | 5.5e-03 | 2.8e-02 |

The single 1.8e-08 is `aluBeforeBellowsVacuum`, a `TGeoPcon` whose consecutive radii differ in the
eighth decimal (2.070000171661377 against 2.0699999332427979 — float32 round-off in the O2 source).
We reproduce every one of those micro-steps; 3.6e-06 cm³ on 198.86 cm³ is `BRepGProp`'s integration
floor over the resulting 7-band revolve, not a geometry difference.

**The composites' 1e-2 is the oracle, not us, and here is the proof.** `TGeoCompositeShape::Capacity()`
is itself a Monte Carlo — it returns a different number on each call (50.4996 and 49.6531 on two
calls for the same shape). Scoring the five worst against an independent 2 000 000-point MC of the
*source* shape:

| volume | `Capacity()` | our OCCT volume | MC ± σ (N = 2e6) | our σ | `Capacity()` σ |
| --- | ---: | ---: | ---: | ---: | ---: |
| `RB24VMABCRBI` | 49.653051 | 49.082427 | 48.812382 ± 0.169071 | 1.6 | 0.6 |
| `bellows1UpperPlie` | 0.246133 | 0.250486 | 0.250049 ± 0.001141 | 0.4 | 0.6 |
| `bellows2vacLowerPlie` | 2.379669 | 2.370584 | 2.370499 ± 0.001111 | **0.1** | **10.6** |
| `BPS_supportBarCarbon` | 37.136273 | 37.151921 | 37.169124 ± 0.024436 | **0.7** | **7.4** |
| `ARB8` | 26.975265 | 26.737681 | 26.702017 ± 0.106453 | **0.3** | **4.5** |

Our volumes sit at 0.1–1.6 σ of the truth; `Capacity()` sits at up to 10.6 σ. Read together with
instrument A's 1 point in 1.69 million, the composites are right and the 1e-2 column is the
estimator's noise.

Separately, all **62 recognised CSG shapes have symmetric difference exactly 0 cm³** against the CAD
part, by the converter's own acceptance measurement.

### 4.4 What came back as what, and why not more

62 of 172 parts were accepted as CSG. The 106 declines, by reason:

| n | reason |
| ---: | --- |
| 46 | *N cap plane(s) perpendicular to the axis, expected 2* |
| 23 | *free-form faces (surface kind outside plane/cylinder/cone/sphere)* |
| 15 | *toroidal faces (out of the recogniser's scope)* |
| 11 | *mixed lateral surface kinds* |
| 6 | *a planar face is neither a cap nor a wedge of any axis cluster* |
| 2 | *N planar faces: not a six-plane box* |
| 1 | *a box face has no opposite partner* |
| 1 | *N coaxial cone faces, expected 1 or 2* |
| 1 | *3 axis clusters: beyond the recogniser's scope (Tier 3 territory)* |

**This is the direct answer to Roadmap item (h), "are we able to recognise Pcons, Polyhedra?"**
No, and the reason is now measured and single: **the largest decline class, 46 of 106, is
"more than two cap planes perpendicular to the axis"** — which is precisely what a `TGeoPcon`
becomes. Every individual face of a Pcon is recognised (they show up as cones, cylinders and
planes in the face table above), but the whole-part rule expects one axis cluster bounded by two
caps, and a Pcon stack has one cap per z plane. The missing piece is exactly the revolved-profile
detector (`Workstreams.md` Stream D, `O2RevolvedSolid`): *all carriers coaxial ⇒ revolution*.
With PIPE's 58 Pcons this bench now has a large, exactly-known corpus to build it against.

The 23 free-form declines are the 22 `TGeoEltu` and the twisted `TGeoArb8` — a recogniser gap for
elliptical cylinders, not a converter defect on either side. The 15 toroidal are the bellows plies,
already known to be out of scope.

## 5. The full Run 3 geometry — stretch, and where it stops

`O2_TGeoToCAD.py full/o2sim_geometry.root FULL.step` on the whole thing:

| | |
| --- | --- |
| logical volumes walked | **27 351** |
| solids built | **21 150** |
| pure assemblies | **6 134** |
| declined | **67** |
| components | **74 601** |
| build + verify wall time | **3 min 08 s** |
| peak RSS | **1.60 GB** |
| the STEP write | **`STEPCAFControl_Writer` dies with SIGSEGV** |

**21 150 of 21 217 shaped volumes convert — 99.7 %.** The 67 declines, every one with its reason
in `FULL_report.json`:

| n | class | reason |
| ---: | --- | --- |
| 43 | `TGeoCompositeShape` | an operand that is itself a composite failed |
| 13 | `TGeoPara` | parallelepiped not mapped |
| 5 | `TGeoCompositeShape` | operand is a `TGeoHalfSpace` (unbounded) |
| 2 | `TGeoTorus` | the hollow-torus cut failed |
| 1 each | `TGeoCompositeShape` | a `Pgon`, a `Trd1`, an `Arb8` operand failed; one union hit the volume invariant |

Volume agreement over the whole geometry, converter's own check: **19 276 non-composite volumes at
median 4.3e-16, max 7.1e-06** (one `TGeoConeSeg`, `DCS021`, which an independent MC puts at 1.3 σ —
so the 7e-06 is the integrator, not the geometry); 1 872 composites at median 7.6e-04, max 2.8e-02
against TGeo's Monte-Carlo `Capacity()`.

So the mapper scales to the whole experiment — every solid is built and volume-checked in three
minutes — and it is **OCCT's STEP writer** that gives out on a 21 150-definition,
74 601-component XCAF document, not the mapping. `--no-step` exists to get the coverage table out
of a model whose write will not complete, and nothing was diagnosed further tonight: the write was
not bisected, no smaller subtree was tried, and no OCCT-side cause was established. Writing
per-detector STEPs with `--top` is the obvious next thing to try, and is untested.

## 6. Dead premises, recorded

- **`BRepPrimAPI_MakeSphere`'s theta arguments.** Assumed to be TGeo's theta. They are latitudes of
  a *plane*-bounded zone, and the resulting solid is 21.7 % too large. Refuted by suite 1 before it
  ever reached a geometry.
- **`TGeoCtub` = tube ∩ half-spaces.** Reasonable, and wrong by 16 %: the planes replace the end
  faces. What made this hard to see is that `TGeoCtub::Capacity()` *is* analytically right
  (π R² · 2dz, independent of the tilt), so it agreed with the wrong construction's *bounding*
  intuition and disagreed with the wrong construction's volume. It took the MC and the Contains
  instrument together to say which side was wrong.
- **"A ruled loft of flat sections gives planar faces."** It gives B-splines. Volume-only testing
  would never have found it — suite 6 exists because of this.
- **"The reverse converter choking means our STEP is wrong."** It meant the ALICE beam pipe has 81
  coincident duplicate placements. Always arbitrate on the TGeo side first; it costs one script.
- **"`IsDone()` means the boolean worked."** It does not. `BRepAlgoAPI_Fuse` returned a valid,
  `IsDone()` solid that had silently lost a whole operand. Nothing but an independent invariant
  would have caught it, and the first version of that invariant was itself wrong — at 1e-7 relative
  it fired on `BRepGProp`'s own ~1e-6 integration noise for a near-disjoint union (0.01 mm³ in
  5 751). A safety net for catastrophic failure has to be priced for catastrophic failure.
- **Carving mother volumes** (the approach the parked `attic/O2_TGeoToCAD_carved*.py` took) is
  *not* used. It makes the STEP a better CAD model and a worse oracle, because no part's capacity
  is comparable against its `TGeoShape::Capacity()` any more. It is available behind
  `--carve-mothers` and was not exercised tonight.

## 7. Reproducing it

```bash
# 1. the geometry (sim env; NEVER with the pythonOCC prepends)
$B/stage/bin/o2-sim -n 0 -g boxgen -m PIPE --configKeyValues "align-geom.mDetectors=none"
#    full geometry: drop -m PIPE and set VMCWORKDIR=$B/stage/share

# 2/3. the converter env is the O2 env PLUS the OCC prepends, in one process
python3 scripts/geometry/O2_TGeoToCAD.py --self-test                     # 71 checks
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root PIPE.step \
        --dedup-world --report PIPE_report.json                          # 3.5 s, 4.36 MB

# 4. back again, and score it
python3 scripts/geometry/O2_CADtoTGeo.py PIPE.step --output-folder rt_pipe -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report rt_pipe/surface_report.json --csg-report rt_pipe/csg_report.json

# 5. the full geometry's coverage table without the write that segfaults
python3 scripts/geometry/O2_TGeoToCAD.py full/o2sim_geometry.root FULL.step \
        --no-step --report FULL_report.json
```

The scoring scripts used for §4.3 (`rtcompare.py`, `rtcap.py`, `dupcheck.py`, `census.py`) live in
this session's scratch and are not committed; every number they produced is quoted above with the
sample size it came from.

## 8. What is open

1. **`TGeoPara`** — the last unmapped class in the Run 3 geometry, 13 volumes. Eight vertices, so
   the `TGeoArb8` path is most of it. With it, the 43 composites that decline because an operand
   declined would also shrink.
1b. **98 reflected placements of volumes with daughters** are dropped in the full geometry (§2).
2. **The full-geometry STEP write.** Bisect it, or write per-detector files with `--top`. Neither
   was attempted.
3. **The revolved-profile detector** on the reverse side. PIPE's 58 Pcons are now an
   exactly-known corpus for it, and it is worth 46 of the 106 CSG declines.
4. **Elliptical cylinders** in the reverse recogniser: 22 of PIPE's 23 tessellated parts.
5. **The 81 coincident placements in the ALICE beam pipe** are a defect in the source geometry, not
   in either converter. Worth reporting upstream.
6. **The one disagreeing point** in 1.69 million, on `bellows2vacLowerPlie`. Not localised.
