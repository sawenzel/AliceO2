# Stream N — implicit (co-surface) trims: is a fitted curve needed at all?

Date: 2026-08-02. Branch `swenzel/bvhsurfacesolid`. **Measurement only — no production code was
changed and nothing was built.** Companions: [`Stream_K_Tier0.md`](Stream_K_Tier0.md) §2/§10.1 (the
plan this measurement was asked to test), [`NEXT.md`](NEXT.md) item 3, [`Tutorial.md`](Tutorial.md).

Two housekeeping notes, recorded because they affect how this document should be read.
**(i)** `Stream_N_PlacedPrimitives.md` was written in parallel and also calls itself "Stream N";
the two are unrelated and neither is a continuation of the other. **(ii)** `O2_CADtoTGeo.py` was
being edited on disk by that parallel work while these measurements ran. It never touched the
recognition or trim path — `_recognize_analytic_surface`, `_recognized_quadric_wire_block`,
`recognize_and_extract_face`, `_quadric_trim_wire` and `extract_planar_face` do not appear in its
diff at all — and the census was re-run against the changed file at the end and returned the
identical 4533 / 3770 / 763 / B 691 / C 56 / A 16. Every number below is from that final run.

---

## 0. The answer first

**Of the 16 ALICE3 leaf solids that are surface-eligible but emit no sidecar, 15 become fully
exactly-trimmable if a boundary edge may be expressed as the intersection of two analytic
surfaces. Not one of their 443 rejected edges needs a fitted curve.**

The 16th, `ST0923290_021`, has **zero** rejected edges: it is declined by a different line
(`recognized quadric trim wire has fewer than 3 edges`) and no trim representation touches it.

Two smaller corpora move as well, through a *second* code path — `extract_planar_face`'s trim-curve
vocabulary, which declines an **ellipse**:

| corpus | today | with implicit co-surface trims | what changes |
| --- | --- | --- | --- |
| ALICE3 `CAD_noETA.stp` | 20 / 55 sidecars, 36 eligible | **35 / 55** | 15 solids, 443 edges, all bucket B |
| `Bagger.step` | 12 / 13 | **13 / 13** | `Bucket`, 4 ellipse edges, plane∩cylinder, dev ≤ 3e-11 cm |
| ladder fixtures | 9 / 10 | **10 / 10** | `oblique_cut_cyl`, 1 ellipse edge, plane∩cylinder, dev 4.6e-13 cm |
| `as1-oc-214.stp` | 5 / 5 | 5 / 5 | nothing to gain; 112 edges, 0 rejected |
| IRIS `ST2487728_01` | 11 / 21 | **11 / 21** | *nothing* — see §7, this is the honest negative |

`Stream_K_Tier0.md` §2 concluded that these edges are "genuinely free-form" and that "there is no
exact representation available", so carrying them required a **fitted** B-spline with a recorded 3D
deviation. That conclusion does not survive the measurement. The edges are not arbitrary curves;
they are the intersections of surfaces we already recognise, and 691 of the 763 currently-rejected
edges lie on **both** of their two analytic surfaces to within the source B-rep's own declared
tolerance.

---

## 1. The question, precisely

`_recognized_quadric_wire_block` accepts a boundary edge only if it is *iso* in the recognised
(φ, other) domain — constant φ (a generator/meridian) or constant `other` (a rim/cap) — and emits
it as a straight segment between the endpoint vertices. Everything else is declined, and one
declined edge declines the whole face, and one declined face declines the whole solid.

A machined part's boundary edge is usually not an arbitrary curve. It is where two surfaces meet.
If both of those surfaces are analytic — and on this corpus they overwhelmingly are — then the edge
has an exact description that needs no parametric curve at all: *the set where this cylinder meets
that plane*. This document measures how much of the corpus that covers.

Every rejected edge is put in exactly one bucket:

| bucket | meaning |
| --- | --- |
| **A** | iso after all — in the shipped parametrisation at the level of the edge's own declared tolerance (**A1**), or in an alternative but equally valid parametrisation of the *same* recognised surface (**A2**) |
| **B** | exactly the intersection of two analytic surfaces (the neighbouring face's surface is analytic, stored *or* recognised, and the edge lies on both) |
| **C** | the neighbour is free-form — no exact co-surface trim exists |
| **D** | anything else, named individually |

---

## 2. The probes, and how to re-run every number here

Both live in [`probes/`](probes/) and are standalone. Nothing here is part of the build.

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
SW=$HOME/alisw/sw/ubuntu2404_aarch64
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
cd $HOME/alisw/O2/scripts/geometry

# the census: every model in this document, ~25 s for ALICE3, 581 MB
python3 probes/trimEdgeCensus.py \
    --model ALICE_3_example/CAD_noETA.stp \
    --model IRIS:IRIS/ST2487728_01-03032026.stp \
    --model STEP_examples/as1-oc-214.stp \
    --model STEP_examples/Bagger.step \
    --fixtures --json /tmp/n/all.json

# every table below
python3 probes/trimEdgeCensusReport.py /tmp/n/all.json

# the reconciliation with Stream_K_Tier0.md §2 (see §6)
python3 probes/trimEdgeCensus.py --model ALICE_3_example/CAD_noETA.stp \
    --legacy-recognizer --json /tmp/n/alice3_legacy.json
python3 probes/trimEdgeCensusReport.py --compare /tmp/n/alice3_legacy.json /tmp/n/all.json
```

| probe | what it measures |
| --- | --- |
| `probes/trimEdgeCensus.py` | loads each model through the converter's own `extract_graph` (so the leaf solids are the converter's `def_shapes`), re-implements the shipped per-edge iso test so it does not stop at the first bad edge, and classifies every rejected edge into A/B/C/D with the measured deviation from both implicit surfaces |
| `probes/trimEdgeCensusReport.py` | turns that JSON into the tables here; no OCC needed |

### 2.1 The method, in one paragraph

For every face the converter cannot extract directly, run the shipping
`_recognize_analytic_surface`. If it returns a cylinder/cone/sphere, build the *same* `project`
closure `recognize_and_extract_face` builds, and run the shipped iso test on all of the face's
boundary edges. For each rejected edge, find the adjacent faces through
`TopExp::MapShapesAndAncestors(EDGE, FACE)` — identity, never proximity — and take the neighbour's
analytic surface, from the STEP file where it is stored analytically and from the same recogniser
where it is not. Then sample the edge's own 3D curve at **129 points** and measure the largest
distance to *both* implicit surfaces, in cm and normalised three ways.

---

## 3. Two self-checks, run before anything else was believed

**(1) The probe's face verdict must be the shipped one.** The probe re-implements the iso test;
if its re-implementation drifted, the census would be about a criterion nothing ships. Over all
998 ALICE3 faces that reach the wire block, and over IRIS, Bagger, as1-oc-214 and the ten fixtures:
**0 verdict mismatches**. (The one deliberate exception is `--legacy-recognizer`, where the probe's
recogniser and the shipped extractor differ by construction; 35 mismatches there, all expected.)

**(2) The deviation instrument must be calibrated on edges we already believe.** The edges the
shipped test *accepts* are rims and generators — a rim circle is exactly where the recognised
quadric meets its neighbouring cap plane, so it is a known-exact co-surface intersection. Measured
with the identical instrument, over **3425** accepted ALICE3 edges that have an analytic neighbour:

```
   accepted-edge control, dev / declared BRep tolerance
   min 6.3e-09   median 1.1e-04   p90 1.02   p99 1.16   max 1.39
```

That is the number that fixes the threshold. **A measurement that only ever sees small numbers on
the population it is arguing about has not been controlled**, and this one does not: the control
itself reaches 1.39 declared tolerances, so any bucket-B threshold below that would reject edges
this project already treats as exact. The worst control edges are named: `ST0923290_016` f#56/f#57
and `ST0923290_019` f#9/f#10/f#27/f#30, all `cylinder|plane`, dev 2.24e-06 cm against a declared
1.60e-06 cm.

---

## 4. The census

### 4.1 Per model

| model | leaf solids | emitted | faces reaching the wire block | boundary edges there | iso | non-iso | A | B | C | D |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ALICE3 `CAD_noETA.stp` | 55 | 20 | 998 (786 cyl / 176 cone / 36 sph) | 4533 | 3770 | **763** | 16 | **691** | 56 | **0** |
| IRIS `ST2487728_01` | 21 | 11 | 998 (same parts) | 4533 | 3770 | 763 | 16 | 691 | 56 | 0 |
| `as1-oc-214.stp` | 5 | 5 | 28 | 112 | 112 | **0** | — | — | — | — |
| `Bagger.step` | 13 | 12 | **0** | 0 | 0 | 0 | — | — | — | — |
| ladder fixtures (10) | 10 | 9 | **0** | 0 | 0 | 0 | — | — | — | — |

**Bagger and the fixtures reach the recognition path zero times**, which reproduces
`Stream_K_Tier0.md` §9.2's "0 of 244 gated faces reach the recognition pre-pass". They are not
silent about implicit trims, though — see §5.

IRIS's per-edge census is cell-for-cell identical to ALICE3's because the two models share the
parts (`Stream_A_CSG.md` §1.6). Its *solid* rollup is not; see §7.

### 4.2 Bucket detail, ALICE3

| bucket | edges | detail |
| ---: | ---: | --- |
| A | **16** | A1 — iso within the edge's own declared tolerance, rejected by the shipped fixed constant |
| B | **691** | edge lies on both implicit surfaces within 2× the declared BRep tolerance |
| C | **56** | the adjacent face is free-form |
| D | **0** | — |

**A2 is empty, and that is a refuted hypothesis, not an absent one.** See §8.1.

### 4.3 Per solid, ALICE3 (only solids with faces reaching the wire block)

`badF` = faces carrying at least one rejected edge; `unsupF` = faces with no analytic surface at
all, which no trim work can rescue.

| solid | nFaces | emitted | recF | badF | edges | non-iso | A | B | C | D | unsupF |
| --- | ---: | :-- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `SOLID` | 50 | no | 2 | 0 | 8 | 0 | 0 | 0 | 0 | 0 | 2 |
| `ST0923290_01` | 352 | no | 246 | 76 | 1058 | 152 | 0 | 144 | 8 | 0 | 28 |
| `ST0923290_002` | 24 | no | 9 | 3 | 48 | 18 | 0 | **18** | 0 | 0 | 0 |
| `ST0923290_003` | 24 | no | 9 | 3 | 48 | 16 | 0 | **16** | 0 | 0 | 0 |
| `ST0923290_004` | 24 | no | 9 | 3 | 48 | 18 | 0 | **18** | 0 | 0 | 0 |
| `ST0923290_005` | 24 | no | 9 | 3 | 48 | 16 | 0 | **16** | 0 | 0 | 0 |
| `ST0923290_006` | 45 | no | 24 | 9 | 132 | 30 | 4 | 26 | 0 | 0 | 5 |
| `ST0923290_007` | 45 | no | 24 | 9 | 132 | 30 | 4 | 26 | 0 | 0 | 5 |
| `ST0923290_008` | 45 | no | 24 | 9 | 132 | 30 | 4 | 26 | 0 | 0 | 5 |
| `ST0923290_009` | 45 | no | 24 | 9 | 132 | 30 | 4 | 26 | 0 | 0 | 5 |
| `ST0923290_010` | 86 | no | 23 | 0 | 92 | 0 | 0 | 0 | 0 | 0 | 6 |
| `ST0923290_011` | 193 | no | 161 | 54 | 779 | 200 | 0 | **200** | 0 | 0 | 0 |
| `ST0923290_012` | 10 | **yes** | 4 | 0 | 16 | 0 | 0 | 0 | 0 | 0 | 0 |
| `ST0923290_013` | 20 | **yes** | 9 | 0 | 36 | 0 | 0 | 0 | 0 | 0 | 0 |
| `ST0923290_014` | 138 | no | 101 | 6 | 408 | 12 | 0 | **12** | 0 | 0 | 0 |
| `ST0923290_015` | 246 | no | 117 | 12 | 502 | 48 | 0 | 0 | **48** | 0 | 68 |
| `ST0923290_016` | 59 | no | 41 | 5 | 172 | 6 | 0 | **6** | 0 | 0 | 0 |
| `ST0923290_017` | 24 | no | 9 | 3 | 48 | 15 | 0 | **15** | 0 | 0 | 0 |
| `ST0923290_018` | 48 | **yes** | 27 | 0 | 108 | 0 | 0 | 0 | 0 | 0 | 0 |
| `ST0923290_019` | 44 | **yes** | 30 | 0 | 120 | 0 | 0 | 0 | 0 | 0 | 0 |
| `ST0923290_020` | 37 | **yes** | 17 | 0 | 68 | 0 | 0 | 0 | 0 | 0 | 0 |
| `ST0923290_021` | 33 | no | 16 | 0 | 62 | **0** | 0 | 0 | 0 | 0 | 0 |
| `ST0923290_022` … `_026` | 24 each | no | 9 | 3 | 48 | 22 | 0 | **22** | 0 | 0 | 0 |
| `ST0923290_027`, `_028` | 24 each | no | 9 | 3 | 48 | 16 | 0 | **16** | 0 | 0 | 0 |

The 26 leaf solids not listed have no face that reaches the wire block at all.

### 4.4 Surface-pair frequency, bucket B (ALICE3)

| pair | edges | max dev (cm) | max dev / declared tol |
| --- | ---: | ---: | ---: |
| `cone ∩ cone` | **196** | 5.6e-11 | 0.0056 |
| `plane ∩ sphere` | **148** | 4.7e-05 | 1.010 |
| `cylinder ∩ cylinder` | **122** | 8.1e-05 | 1.024 |
| `cylinder ∩ plane` | 80 | 3.9e-05 | 1.013 |
| `cone ∩ cylinder` | 62 | 3.3e-05 | 1.021 |
| `sphere ∩ sphere` | 44 | 5.4e-11 | 0.0054 |
| `cylinder ∩ sphere` | 33 | 8.1e-12 | 0.00081 |
| `cone ∩ plane` | 6 | 1.4e-05 | 1.006 |

**No single pair dominates.** The five common quadric pairs plus `plane ∩ *` cover everything; a
first useful implementation cannot be one special case. Restricted to the **15 solids that the
rollup actually turns on** (443 edges): `plane∩sphere` 148, `cylinder∩cylinder` 122,
`cone∩cone` 64, `sphere∩sphere` 44, `cylinder∩sphere` 33, `cone∩cylinder` 26, `cone∩plane` 6 —
seven pairs, no `cylinder∩plane` at all. 289 of those 443 neighbours are themselves **recognised**
NURBS quadrics and 154 are stored analytic surfaces, so the co-surface trim has to compose with the
recognition path, not merely with the STEP file.

---

## 5. The second mechanism: an ellipse on a plane

`extract_planar_face` declines a planar face whose boundary carries a curve outside
{line, circle, bspline, bezier}. On ALICE3 that never fires. On the two *other* corpora it is the
**only** thing standing between them and full coverage, and the curve it declines is an **ellipse**
— which is exactly the brief's motivating example: a cylinder cut by an oblique plane.

| solid | face | curve | neighbour | dev from both surfaces | declared tol | dev / tol | dev / edge length |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| `oblique_cut_cyl` | f#1 | ellipse | cylinder (stored) | **4.6e-13 cm** | 1e-08 cm | 4.6e-05 | 3.9e-14 |
| `Bagger/Bucket` | f#4 | ellipse | cylinder (stored) | 3.0e-11 cm | 1e-08 cm | 0.0030 | 1.2e-10 |
| `Bagger/Bucket` | f#4 | ellipse | cylinder (stored) | 2.2e-12 cm | 1e-08 cm | 0.00022 | 8.5e-12 |
| `Bagger/Bucket` | f#31 | ellipse | cylinder (stored) | 2.2e-12 cm | 1e-08 cm | 0.00022 | 8.5e-12 |
| `Bagger/Bucket` | f#31 | ellipse | cylinder (stored) | 3.0e-11 cm | 1e-08 cm | 0.0030 | 1.2e-10 |

These are exact at machine precision, not "within tolerance". `oblique_cut_cyl` is the ladder
fixture [`NEXT.md`](NEXT.md) records as *"has no sidecar and never has"*; **one ellipse edge is the
entire reason**, and it is a plane∩cylinder intersection.

`Bagger/Bucket` is more surprising. [`Tutorial.md`](Tutorial.md) §6 attributes its tessellated
fallback to "97 faces, 4 spherical + 2 toroidal". Measured today, `Bucket` has **zero** faces
without an analytic surface: its four rejected edges are the whole story, and Bagger's exact-surface
coverage is one trim feature away from **13/13**.

---

## 6. Re-deriving §2's 1891 / 1053: it reproduces, but only against the recogniser of the day

Against the **shipping** converter the figures are different, and the difference is not small:

| | Stream K §2 | this probe, shipping recogniser |
| --- | ---: | ---: |
| faces carrying ≥1 non-iso edge | 373 | **225** |
| boundary edges on those faces (non-degenerate) | 1891 | **1303** |
| iso | 834 | 540 |
| non-iso | 1057 (= 1053 free-form + 4 arcs) | **763** |

Run against the **pre-fix** recogniser — the code as of `git 237be7f81a^`, copied verbatim into the
probe and reachable with `--legacy-recognizer` — every cell comes back:

| | Stream K §2 | `--legacy-recognizer` |
| --- | ---: | ---: |
| recognised faces reaching the wire block | 1180 (786 cyl / 358 cone / 36 sph) | **1180 (786 / 358 / 36)** |
| faces carrying ≥1 non-iso edge | 373 | 372 |
| boundary edges on those faces, **excluding degenerate** | 1891 | **1891** |
| iso | 834 | **834** |
| non-iso | 1057 | **1057** |

So the 1891 / 834 / 1057 triple is exact (the face count is off by one, and the denominator
excludes the 99 degenerate edges on those faces — 1990 including them). What changed the numbers is
`Stream_K_Tier0.md` §5's own acceptance-criterion fix, which landed *after* §2 was measured and
removed 182 phantom cones on `ST2487458_01`. **287 of §2's 1053 "genuinely free-form" edges were on
those phantom cones** and are no longer recognised at all:

| | pre-fix | shipping |
| --- | ---: | ---: |
| bucket B | 691 | **691** |
| bucket C | 343 (`ST2487458_01` 287, `ST0923290_015` 48, `ST0923290_01` 8) | **56** |
| bucket A | 16 | 16 |
| bucket D | 7 (all `ST2487458_01`) | **0** |

Bucket B is **identical** in both runs, which is the point: the co-surface finding is not an
artefact of the recogniser change. What the recogniser change removed was noise in bucket C.

**The honest reading of `Stream_K_Tier0.md` §2's headline sentence** — *"1053 are genuinely
free-form there … there is no exact representation available"* — is that 691 of them are exact
co-surface intersections, 343 have a free-form neighbour (287 of those on faces that are no longer
claimed to be quadrics), 16 are iso within the model's own tolerance, and 7 are neither.

---

## 7. The solid-level rollup — the number that decides

Population: leaf solids that emit no sidecar and whose **every** declined face was declined for a
*trim* reason. A solid with even one genuinely unsupported surface is out of scope, because no trim
work rescues it.

### ALICE3 — **15 of 16**

| solid | rejected edges | buckets | verdict |
| --- | ---: | --- | --- |
| `ST0923290_002` | 18 | B 18 | **fully covered** |
| `ST0923290_003` | 16 | B 16 | **fully covered** |
| `ST0923290_004` | 18 | B 18 | **fully covered** |
| `ST0923290_005` | 16 | B 16 | **fully covered** |
| `ST0923290_011` | 200 | B 200 | **fully covered** |
| `ST0923290_014` | 12 | B 12 | **fully covered** |
| `ST0923290_016` | 6 | B 6 | **fully covered** |
| `ST0923290_017` | 15 | B 15 | **fully covered** |
| `ST0923290_022` … `_026` | 22 each | B 22 | **fully covered** |
| `ST0923290_027`, `_028` | 16 each | B 16 | **fully covered** |
| `ST0923290_021` | **0** | — | **still fails** — `recognized quadric trim wire has fewer than 3 edges` |

443 edges, every one of them bucket B. Deviation over that whole set: **max 8.1e-05 cm**,
max 1.024 declared tolerances, max 1.2e-04 patch diagonals, max 7.7e-03 edge lengths.

`ST0923290_021` is the one solid the implicit trim does **not** reach, and it is not a geometry
problem: one recognised cone's trim wire has fewer than three edges (a full-circle cap written as
one or two edges), and the wire block requires three. That is an implementation limit worth its own
line in the ledger, not a representation gap.

### Bagger — **1 of 1**, and the fixtures — **1 of 1**

| model | solid | rejected | buckets | verdict |
| --- | --- | ---: | --- | --- |
| `Bagger.step` | `Bucket` | 4 | B 4 | **fully covered** → Bagger 13/13 |
| fixtures | `oblique_cut_cyl` | 1 | B 1 | **fully covered** → fixtures 10/10 |

### IRIS — **0 of 0**, and this is the load-bearing negative

IRIS carries the same parts and produces a **bit-identical edge census** (4533 / 763 / 691 / 56 /
16), yet its rollup gains nothing, because its XCAF definition labels group the geometry into 21
leaf solids instead of 55: every recognised face lives inside one 1734-face solid, `ST0923290_01`,
which also carries **122 faces with no analytic surface** and would still fall back to tessellation.

**Coverage is a property of how the assembly decomposes into leaf solids, not only of the trim.**
The same trim feature is worth 15 solids on ALICE3 and 0 on IRIS. Any coverage claim made from this
work must name its model *and* its decomposition.

---

## 8. Being adversarial about this result

### 8.1 A2 was a plausible free win and it is refuted

The converter hard-codes a recognised sphere's polar axis to (0, 0, 1)
(`recognize_and_extract_face`, the `sphere` branch). A circle of constant θ about an axis **n** is
exactly the sphere cut by a plane perpendicular to **n**, so *every planar edge on a sphere is iso
for some polar axis*. On ALICE3, **225 of the 763 rejected edges are planar to within their own
declared tolerance** and sit on recognised spheres. That looked like free coverage for the existing
code path with no new machinery at all.

It is not, and the face-level test says so. `_recognized_quadric_wire_block` needs **one** frame for
the whole face, and these faces are bounded by circles in *different* planes. Of the **36** sphere
faces carrying rejected edges, **0** are rescued by any single polar axis — the probe tries the
shipped z, every rejected edge's own plane normal, and every pairwise intersection of two edge
planes (36–37 candidates per face) and none serves all edges.

Named, so it can be re-run: `ST0923290_002` f#8 has 8 boundary edges, 6 rejected, whose neighbours
are one sphere, one cylinder and four planes. No axis makes a sphere-circle, a cylinder-circle and
four plane-circles simultaneously iso.

**A2 is therefore reported as 0.** All 225 edges are in bucket B instead — which is the stronger
statement anyway: their deviation from both implicit surfaces is ≤ 1.010 declared tolerances.

### 8.2 The threshold is calibrated, and the answer does not depend on it

`--b-factor` sweeps the only free constant. The answer is flat:

| b-factor | buckets | solids fully covered (of 16) |
| ---: | --- | ---: |
| 1.0 | B 398, D 293, C 56, A 16 | **1** |
| 1.5 | B 691, C 56, A 16 | **15** |
| 2.0 (default) | B 691, C 56, A 16 | **15** |
| 3.0 | B 691, C 56, A 16 | **15** |

The transition sits at **1.025** — the largest dev/tol among all rejected edges with an analytic
neighbour is 1.024 — and there is nothing between 1.03 and ∞. The accepted-edge control (§3)
independently *forces* the threshold to be at least **1.395**. So the interval that the control
demands and the interval that produces the answer are disjoint from the interval that would change
it. The constant is over-determined, not tuned.

And the threshold cannot manufacture a B: the 56 bucket-C edges have **no** analytic neighbour, so
no value of `--b-factor` reaches them.

### 8.3 The tail, not the mean

| quantity | n | min | p50 | p90 | p99 | **max** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| B, dev (cm) | 691 | 1.3e-14 | 3.0e-06 | 3.1e-05 | 6.5e-05 | **8.1e-05** |
| B, dev / edge length | 691 | 5.4e-14 | 2.6e-05 | 3.2e-04 | 1.0e-03 | **7.7e-03** |
| B, dev / patch diagonal | 691 | 2.6e-14 | 1.2e-06 | 2.6e-05 | 6.7e-05 | **1.2e-04** |
| B, dev / declared tolerance | 691 | 1.7e-09 | 1.00 | 1.00 | 1.01 | **1.02** |
| control (accepted edges), dev / declared tol | 3425 | 6.3e-09 | 1.1e-04 | 1.02 | 1.16 | **1.39** |

The ten worst bucket-B edges by dev/tol, so any one can be re-run alone:

```
ST0923290_011 f#102  cylinder|cylinder  dev 1.33e-05 cm  tol 1.30e-05 cm  ratio 1.024
ST0923290_011 f#129  cylinder|cylinder  dev 1.33e-05 cm  tol 1.30e-05 cm  ratio 1.024
ST0923290_011 f#72   cone|cylinder      dev 2.44e-05 cm  tol 2.39e-05 cm  ratio 1.021
ST0923290_011 f#109  cone|cylinder      dev 2.44e-05 cm  tol 2.39e-05 cm  ratio 1.021
ST0923290_006 f#19   cylinder|plane     dev 5.24e-06 cm  tol 5.17e-06 cm  ratio 1.013
ST0923290_007 f#19   cylinder|plane     dev 5.24e-06 cm  tol 5.17e-06 cm  ratio 1.013
ST0923290_008 f#19   cylinder|plane     dev 5.24e-06 cm  tol 5.17e-06 cm  ratio 1.013
ST0923290_009 f#19   cylinder|plane     dev 5.24e-06 cm  tol 5.17e-06 cm  ratio 1.013
ST0923290_006 f#21   cylinder|plane     dev 5.24e-06 cm  tol 5.17e-06 cm  ratio 1.013
ST0923290_008 f#21   cylinder|plane     dev 5.24e-06 cm  tol 5.17e-06 cm  ratio 1.013
```

The distribution is **bimodal and the reason is structural.** Half of bucket B sits at machine
precision (dev 1e-14 … 1e-11 cm — the `cone∩cone`, `sphere∩sphere` and `cylinder∩sphere` families,
where the declared edge tolerance is OCCT's floor of 1e-08 cm). The other half sits at
dev/tol ≈ 1.000, and that is not a coincidence: **OCCT enlarges an edge's tolerance to cover the
discrepancy it measured**, so on those edges the declared tolerance *is* the deviation. Reading the
ratio as "just inside a band" is the wrong reading; the correct one is that the deviation and the
tolerance are the same measurement made twice.

Which term dominates is also measured: on 470 of 691 edges it is the distance from the edge to its
**own** recognised surface, not to the neighbour's. That says the residual is the source B-rep's own
edge-versus-surface slack — the same slack that shows up on the edges the converter already accepts
— and not a failure of the co-surface hypothesis. The largest `dev / edge length` values
(7.7e-03) are all on 15–35 µm edges on `ST0923290_011` where the absolute deviation is the declared
tolerance; the ratio is small edges, not large errors.

### 8.4 Is anything rejected for a reason other than non-isoness? Yes — 16 edges, a tolerance defect

The shipped test uses `tol_other = 1e-6 × max(1.0, span of *other* over the whole face)` and
`tol_phi = 1e-7` rad. Both are **fixed constants unrelated to the model's own declared tolerance**,
which ALICE3 carries per edge and which the sidecar already stores as `modelTolerance`.

16 edges — 4 each on `ST0923290_006`, `_007`, `_008`, `_009`, faces f#6, f#8 (×2) and f#10, all
recognised cones — are cap circles that are iso *within their own declared tolerance* and are
rejected anyway:

| solid | face | edge length | span in *h* | edge's declared tolerance | shipped `tol_other` |
| --- | --- | ---: | ---: | ---: | ---: |
| `ST0923290_006` | f#6 | 0.478 cm | 2.6e-05 cm | **5.2e-05 cm** | 1e-06 cm |
| `ST0923290_006` | f#8 | 0.237 cm | 9.8e-06 cm | **1.6e-05 cm** | 1e-06 cm |
| `ST0923290_006` | f#8 | 0.241 cm | 9.2e-06 cm | **1.2e-05 cm** | 1e-06 cm |
| `ST0923290_006` | f#10 | 0.478 cm | 8.3e-06 cm | **1.2e-05 cm** | 1e-06 cm |

(identical on `_007`, `_008`, `_009` — the four are the same part placed four times).

The constant is 12× to 52× tighter than what the model says about itself. Whether to widen it is a
**policy** decision of exactly the kind `Stream_K_Tier0.md` §10.2 describes for the 8 planes on
`ST2487195_01`, and it should be decided by measuring what the resulting trim error does to
transport — not by picking a number. It is recorded here because it is a rejection that has nothing
to do with the edge being non-iso.

**The other 747 rejected edges are nowhere near the constant.** The 10th percentile of
`min(span_other/tol_other, span_phi/tol_phi)` is **4.6e+03** and the median **4.5e+04**. Only 12
edges are within 10× of the shipped constant and only 16 within 100×. Relaxing the constant is
therefore worth 16 edges and no solids; it is a correctness/consistency item, not a coverage lever.

### 8.5 Things that are *not* wrong

- **No periodicity or seam artefacts.** Zero rejected edges (and zero accepted edges) have an
  unwrapped φ span exceeding 2π; the largest is exactly π. The continuous per-sample unwrap in
  `_recognized_quadric_wire_block` is doing its job.
- **No bucket D at all on the shipping recogniser.** Every rejected edge has exactly one
  neighbouring face; there is no seam edge, no non-manifold edge and no edge without a 3D curve in
  the whole rejected population.
- **Bucket C is confined to two solids that fail anyway.** All 56 are cylinders on
  `ST0923290_015` (48 edges, faces f#7/34/39/41/139/144/226/229/231/232/237/240) and
  `ST0923290_01` (8 edges, faces f#116/117/120/121/122/123/125/126). Both solids carry 68 and 28
  faces with no analytic surface at all, so neither is in the rollup population and neither would
  emit even if every one of those 56 edges were carried exactly.

---

## 9. What this does and does not license us to build

### It licenses

1. **Dropping the fitted-curve plan** (`Stream_K_Tier0.md` §2/§10.1, `NEXT.md` item 3) as the route
   to the 16 ALICE3 solids. It is not needed for any of them. The project's standing bargain —
   *exact, or tessellated, never quietly approximate* — does not have to be spent here.
2. **Investigating an implicit / co-surface trim representation** as the replacement, with the
   measured scope in hand: 8 surface-pair kinds, none dominant, 289 of 443 neighbours themselves
   recognised NURBS quadrics rather than stored analytic surfaces.
3. **A cheap, separate win on the planar path**: `plane ∩ cylinder` alone (an ellipse boundary on a
   planar face) takes Bagger to 13/13 and the ladder to 10/10, at machine precision, and it is a
   *stored* analytic neighbour in every case measured. That is the smallest useful first version
   and it is testable against two already-gated corpora.

### It does not license

1. **"Implicit trims give exact geometry."** They give geometry as exact as the *source B-rep* —
   which on 293 of the 691 ALICE3 edges means a deviation equal to the edge's declared tolerance,
   up to **8.1e-05 cm**. That is not machine precision. It is the same standing as sidecar v3's
   edge identity and as `modelTolerance`: exact *relative to what the model asserts about itself*.
   Any acceptance criterion must be written against the declared tolerance, and the achieved
   deviation must be recorded per edge exactly as §10.1 required of the fitted version.
2. **Any statement about the kernel.** Nothing here has been through `O2BVHSurfaceSolid`, the
   oracle gate or the X-ray benchmark. The 15 solids would produce sidecars nothing has ever
   loaded; ALICE3 already has 2 emitted sidecars that do not load at all
   (`Stream_L_ALICE3Defect.md` mechanism 3, the fixed 1e-06 cm wire-join tolerance) and 15 new ones
   would be 15 new chances to hit it. **The rollup is a converter-side count, not a coverage
   claim.**
3. **Any statement about correctness of the resulting solids.** `Stream_L_ALICE3Defect.md` is
   explicit that `reliable` and `navigable` were not evidence, and that four ALICE3 parts lost 418
   crossings while reporting both. A solid that newly *emits* has not been shown to *transport*.
4. **The IRIS number, or any number, without its decomposition.** §7: identical edges, zero solids.
5. **Retiring `ST0923290_021`.** It needs a wire block that accepts a trim wire of fewer than three
   edges. Different problem, different fix, one solid.
6. **Anything about the 19 ALICE3 solids with genuinely free-form *surfaces*.** They are untouched
   by this and remain `NEXT.md` item 4.

### The order I would take it in

1. **The planar ellipse** (§5). Smallest, exact at machine precision, moves two *gated* corpora
   (Bagger 12→13, fixtures 9→10) so the oracle gate and the X-ray benchmark can actually score it.
   That makes it the only piece of this work that can be verified with the instruments that exist.
2. **Design the implicit trim's acceptance before its representation.** The natural bound is the
   edge's own declared tolerance; the achieved deviation per edge is what this probe already
   measures, and `probes/trimEdgeCensus.py` is the instrument that would score it.
3. **`plane ∩ sphere` and `cylinder ∩ cylinder`** (148 + 122 edges) — the two largest families among
   the 15 solids, and between them they carry `ST0923290_002/003/004/005/017/022…028` and
   `ST0923290_011`.
4. **The `tol_other` / `tol_phi` constants** (§8.4). 16 edges, 0 solids, but it is a guard compared
   against a quantity the model has a better answer for — the same defect class as
   `Stream_E_Scale.md` §5 and `Stream_K_Tier0.md` §3.

---

## 10. What could not be measured

- **Whether the resulting solids are navigable, closed or transport-correct.** No sidecar was
  written; nothing was gated. See §9's second and third disclaimers.
- **Whether an implicit trim can be *evaluated* cheaply in the kernel.** This document says the
  information exists and is exact; it says nothing about the cost of a point-in-trim test expressed
  as "inside every neighbouring half-space", nor about how such a trim orders its boundary, nor
  about how `CloseShape` would count its edges. Those are design questions and none of them was
  touched.
- **`oTOF System V3-R92cm.step`**, the fourth corpus in `STEP_examples/`, was not censused.
- **The 4 arcs** `Stream_K_Tier0.md` §2 counts separately (non-iso but an exact arc in (φ, other)).
  This probe classifies by neighbour and deviation, not by the curve's shape in the recognised
  domain, so it neither confirms nor contradicts that cell; the 1057 total it belongs to does
  reproduce (§6).
