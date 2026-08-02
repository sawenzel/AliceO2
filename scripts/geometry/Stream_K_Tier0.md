# Stream K — Tier 0: canonical recognition of NURBS-encoded quadrics

Date: 2026-08-02. Branch `swenzel/bvhsurfacesolid`. Converter-side work only; nothing under
`Detectors/Base/**` is touched. Companions: [`Stream_A_CSG.md`](Stream_A_CSG.md) §2.3 (the census
that scoped this), [`Stream_L_ALICE3Defect.md`](Stream_L_ALICE3Defect.md) (the inverted-normal trap
this path already carries), [`Stream_J_XRay.md`](Stream_J_XRay.md) (the verification instrument).

**Read §1 first. The brief this stream was given rests on a premise that is no longer true, and the
value has moved somewhere else.**

---

## 1. The dead premise, stated before anything else

> "CAD exporters routinely write exact analytic surfaces as rational NURBS patches. […] the
> converter dispatches on the *stored* type, so the analytic geometry is discarded and the solid
> falls out of the exact path."

**The converter has not dispatched on the stored type since before this stream started.**
`O2_CADtoTGeo.py` already carries a full canonical-recognition pre-pass — `_recognize_analytic_surface`
(a normal-field differential recognizer ported from `analyze_surface_geometry.py`),
`recognize_and_extract_face`, `_recognized_quadric_wire_block` and the `--recognize-surfaces
exact|off` switch — and it **fires**. Measured, on this branch, before any change:

| model | stored `bspline` faces | recognized by the converter today |
| --- | ---: | ---: |
| `as1-oc-214.stp` | 28 | **28** (all cylinders) |
| `CAD_noETA.stp` (ALICE3) | 2377 | **1180** (786 cylinder, 358 cone, 36 sphere) |

Two of the brief's acceptance targets are therefore already met, and were met before this session:

- **`as1-oc-214.stp` is `5/5`, not `0/5`.** All five leaf solids emit an exact sidecar today.
  The `0/5 → 5/5` figure in the brief is `Stream_A_CSG.md`'s **CSG quadric-only census column**
  (`quadricOnly` 0 → `quadricOnlyAfterTier0` 5), which is a statement about what a *CSG*
  recogniser could reach, not about sidecar emission. The two were conflated.
- ALICE3's `n_eligible` in `surface_report.json` is **36/55** — exactly the census's
  `quadricOnlyAfterTier0`. The surface half of Tier 0 is done and its number is already banked.

What is *not* done, and where the whole remaining coverage lever sits, is §2. What is *wrong*, and
where the remaining correctness risk sits, is §3.

---

## 2. Where the coverage actually is: the trim, not the surface

> **Superseded in its conclusion (2026-08-02) — see [`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md).**
> The diagnosis below is right: the 16 solids fail on the trim, not on recognition, and the edges
> are not iso and not straight in (φ, other). The *conclusion* — "no exact representation is
> available, therefore a fitted curve" — is not. A per-edge census with the neighbouring face's
> surface in hand finds that **691 of the currently-rejected 763 edges are exactly the intersection
> of two analytic surfaces we already recognise**, and that **15 of the 16 solids** are covered
> entirely by that, on 443 edges, none of which needs a fit. The 1891 / 834 / 4 / 1053 table below
> reproduces cell-for-cell **only against the pre-fix recogniser** (§5 landed after it was
> measured, and removed 287 of the "free-form" edges along with `ST2487458_01`'s phantom cones);
> against the shipping converter it is 1303 / 540 / — / 763.

ALICE3 emits **20** sidecars while **36** solids are surface-eligible. The 16 missing solids do not
fail on recognition. They fail *after* it, in `_recognized_quadric_wire_block`.

Per leaf solid, running the converter's own `extract_surfaces_for_shape` over the converter's own
`def_shapes` (this matters — a `STEPControl_Reader` load of the same file heals differently and
gives 17, not 20):

| outcome | solids |
| --- | ---: |
| every face extracts → sidecar emitted | **20** |
| fails **only** on recognized-quadric trim | **16** |
| fails on trim *and* on genuinely free-form faces | 7 |
| fails only on genuinely unsupported surfaces (bspline / extrusion / revolution) | 12 |

The 16 are `ST0923290_002/_003/_004/_005/_011/_014/_016/_017/_021/_022/_023/_024/_025/_026/_027/_028`.
Fifteen of them are declined by exactly one line:

```
recognized as {cone,cylinder,sphere} but recognized quadric boundary edge is not axis-aligned in (phi, h/theta)
```

`_recognized_quadric_wire_block` accepts a boundary edge only if it is iso in the recognized
(phi, other) domain — constant `phi` (a generator/meridian) or constant `other` (a rim/cap) — and
then emits it as a straight segment between the two endpoint vertices. 373 recognized faces on
ALICE3 carry at least one edge that is neither.

**Those edges are not near-iso, and they are not straight either.** Sampled at 33 points and
classified in the recognized domain:

| class | edges | max residual |
| --- | ---: | ---: |
| iso (already accepted) | 834 | — |
| straight in (phi, other) but not iso | **0** | — |
| an arc in (phi, other) | 4 | 8.1e-07 |
| **genuinely free-form in (phi, other)** | **1053** | 6.7e-01 relative |

So relaxing the iso test to a straightness test recovers nothing, and this needs saying because it
is the obvious first move.

**And there is no exact representation available for them in the current sidecar format.** This is
structural, not a missing feature. On an OCC *canonical* cylinder the surface parameter `u` **is**
the angle, so a pcurve lives in (phi, h) already and `_quadric_trim_wire` maps it under an affine
map — and a B-spline is closed under affine maps, so the trim is carried exactly. On a
NURBS-*encoded* cylinder the angle is a rational (tangent-half-angle) function of `u`, so
`phi(u)` is transcendental and the image of a polynomial pcurve under it is not a line, an arc, or
a B-spline in (phi, h). The sidecar's trim vocabulary is line / arc / B-spline in that domain.

Carrying these 1053 edges therefore means **fitting** a B-spline in (phi, other) to a declared
tolerance and recording the achieved 3D deviation of the reconstructed boundary against the source
edge — a real, bounded approximation of the *trim* (not of the surface), of the same standing as
the pcurve approximations OCC itself stores within edge tolerance. That is a defensible design and
it is the single largest remaining converter-side coverage item (20 → 36 ALICE3 sidecars, and it is
the only thing between `n_eligible = 36` and `emitted = 20`). **It is scoped here and deliberately
not built in this stream**, because the brief's standing instruction is that a declined face costs
coverage and a wrongly-accepted face costs correctness, and this project trades the first for the
second every time.

---

## 3. What *is* wrong: the recognizer's acceptance criterion is not a criterion

`_recognize_analytic_surface` scores each candidate model with whatever expression falls out of its
own linear solve, and then compares those numbers against one tolerance (`_RECOGNIZE_TOL_EXACT
= 1e-9`) **and against each other**, to select a model. They are not the same quantity:

| candidate | what its "residual" was | dimension |
| --- | --- | --- |
| plane | `max ‖N·N₀‖ − 1‖` | an angle |
| sphere | `max ‖‖P−C‖ − r‖ / diag` | a relative distance |
| cylinder | `max ‖‖P⊥‖ − R‖ / diag` | a relative distance |
| **cone** | `max ‖u·N‖`, then `max(that, std ‖u·axis‖)` | **an angle** |

For the cone this is not conservative, and the consequence is measured, not argued.

### 3.1 The measurement

An independent instrument (`/tmp` probe, self-checked below) computes **one** number for every face
the converter recognizes: the largest distance from the recognition samples to the recognized
surface, in cm and relative to the patch's own bounding-box diagonal. Alongside it, OCCT's
`ShapeAnalysis_CanonicalRecognition` — a completely separate recogniser — is asked the same
question. On ALICE3, all 55 leaf solids, 1188 rows:

| converter says | OCCT says | faces | max gap (cm) | max gap / diagonal | median gap / diagonal |
| --- | --- | ---: | ---: | ---: | ---: |
| cylinder | cylinder | 786 | 3.8e-11 | 2.4e-10 | 4.3e-13 |
| cone | cone | 174 | 2.3e-12 | 4.0e-12 | 1.8e-13 |
| sphere | sphere | 36 | 1.2e-11 | 1.1e-11 | 1.0e-11 |
| **cone** | **nothing** | **184** | **7.9e+01** | **3.7e+01** | **1.2e+00** |
| nothing | plane | 8 | — | — | — |

786 + 174 + 36 + 8 = **1004**, which is `Stream_A_CSG.md`'s census number **cell for cell**. Three
instruments — the STEP-text parse, the OCCT census, and this one — agree on which faces are
quadrics in disguise. The converter agrees with all three on the cylinders and the spheres, misses
8 planes, and **accepts 184 faces that none of the others do, whose recognized cone misses the real
surface by up to 79 cm — thirty-seven patch diagonals.**

The instrument's own positive control, run first: the same recognizer and the same gap measurement
applied to faces whose *stored* surface already is the analytic thing must return that same kind at
a gap of zero. Measured over a 1-in-7 sample of ALICE3's natively-analytic faces —
plane 465 faces max gap **4.4e-15 cm**, cylinder 395 max **1.9e-11 cm**, cone 62 max **1.9e-14 cm**.

### 3.2 The mechanism, named

The 184 are all on `ST2487458_01`. Anatomy of one (`f#500`) against a good cone (`ST0923290_011 f#1`):

| | bad | good |
| --- | ---: | ---: |
| patch diagonal | 2.19e+01 | 2.34e+00 |
| fitted apex distance / diagonal | **3.3e+11** | 3.5e-01 |
| recognized half-angle | **0.0000** | 1.0472 |
| `max ‖u·N‖` (the cone test) | 1.5e-12 ✓ | 1.9e-12 ✓ |
| `std ‖u·axis‖` (the half-angle test) | 1.0e-16 ✓ | 3.6e-12 ✓ |
| cond(N) | **2.5e+12** | 1.1e+01 |
| smallest singular value of the normal field | 4.0e-13 (relative) | — |
| **measured gap** | **78 cm** | 9.4e-13 cm |

The face is an **extruded free-form profile** — a "generalized cylinder". Its normals are coplanar
to 7.8e-13 (it passes the cylinder branch's *coplanarity* gate) but the best-fit circle misses the
profile by 7.0e-04 mm on a 9.8 mm radius, so the cylinder branch correctly declines at a relative
residual of 3e-05. The normal field is then nearly rank-2, the least-squares apex runs off along
the extrusion direction to 3e11 patch diagonals, the half-angle collapses to zero — and **both cone
tests pass vacuously**: with every `u` collapsed onto one direction, `u·N ≈ 0` merely restates that
the normals are perpendicular to *that* direction, and the spread of `u·axis` is the spread of a
constant.

This is `Stream_E_Scale.md` §5's defect class in a new place: a guard whose value is compared
against a quantity of the wrong dimension. Here it is worse than a wrong constant — the quantity
being guarded stops being informative exactly when the fit stops being valid.

### 3.3 It is latent, not shipped — and that qualifier belongs next to the number

All 182 faces whose relative gap exceeds 1e-9 are on **`ST2487458_01`, which emits no sidecar**
(889 of its 2034 faces fail, 742 of them genuinely free-form). Both gated corpora reach the
recognition pre-pass **zero** times (`Stream_L_ALICE3Defect.md` §6 measured this: 0 of 244 faces).
So nothing that ships today carries a bogus cone.

What makes it worth fixing now rather than later is that it is **exactly** the trap
`Stream_L_ALICE3Defect.md` describes, one level up: it is invisible to every existing check, it
fires only on NURBS-encoded surfaces, and **the §2 trim work would multiply the number of faces
that reach it** — the moment a recognized face can carry a general trim, `ST2487458_01`'s 182
phantom cones become emittable.

---

## 4. The self-test, and the red it started from

`python3 scripts/geometry/O2_CADtoTGeo.py --self-test` builds every control in-process from OCC
primitives — no CAD file, no fixture — and asserts the shipping `_recognize_analytic_surface`
against answers known in closed form. The model is `csg/census.py --self-test`, and so is its
lesson: **a recogniser's positive control is worthless without a negative one.**

| control | expectation |
| --- | --- |
| `BRepBuilderAPI_NurbsConvert` of a `gp_Cylinder` / `gp_Cone` / `gp_Sphere` / `gp_Pln` patch | recognized as that kind — this *is* the exporter artefact, rebuilt |
| a free-form saddle; a narrow free-form ridge | declined |
| a NURBS-encoded torus | declined (no torus model — a stated limitation, asserted so it cannot drift) |
| **a swept non-circular profile** (`r = R(1 + b·cos 3a)(1 + ε·t)`), 6 combinations | **declined** |
| the invariant: every accepted face's **measured** gap < the acceptance tolerance | holds |

The swept-profile control is `ST2487458_01` reduced to a closed form: `b` makes it genuinely
non-circular, `ε` tunes how nearly rank-2 the normal field is — the knob that sends the apex to
infinity. It reproduces the defect at every combination tried.

**Before the fix — 18 checks, 7 failures:**

```
  [FAIL] swept non-circular profile (bulge 1e-03, taper 1e-04) is declined -- got cone
  [FAIL] swept non-circular profile (bulge 1e-03, taper 1e-06) is declined -- got cone
  [FAIL] swept non-circular profile (bulge 1e-03, taper 1e-08) is declined -- got cone
  [FAIL] swept non-circular profile (bulge 1e-02, taper 1e-04) is declined -- got cone
  [FAIL] swept non-circular profile (bulge 1e-02, taper 1e-06) is declined -- got cone
  [FAIL] swept non-circular profile (bulge 1e-02, taper 1e-08) is declined -- got cone
  [FAIL] every accepted face's MEASURED gap is below the acceptance tolerance
         -- worst cone at 2.11e-01 against 1e-09
```

Every positive control passes, at gaps of 8e-17 to 6e-16 of the patch diagonal, which is what makes
the failures mean something: the recogniser is not broken, its *acceptance* is.

---

## 5. The fix, and what it moves

One change, in both the shipping recognizer (`O2_CADtoTGeo.py::_recognize_analytic_surface`) and
the reference implementation it was ported from (`analyze_surface_geometry.py::classify_surface`):

> **Every candidate is scored by exactly one quantity** — `_analytic_surface_gap` (the largest
> distance from a recognition sample to the candidate surface) divided by the sample bounding-box
> diagonal — and that same number is the acceptance test against `_RECOGNIZE_TOL_EXACT = 1e-9`.

The linear solves still *propose* the models; none of them decides any more. A degenerate proposal
(the apex at 1e11 diagonals) is left in and simply scored like any other — the gap sees what the
solve's own residual cannot.

`--self-test`: **18 checks, 0 failures** (was 18 checks, 7 failures). The invariant line now reads
`worst cylinder at 5.60e-16 against 1e-09`.

### 5.1 Effect on the three corpora — totals and disagreements together

| | before | after |
| --- | --- | --- |
| `as1-oc-214.stp` sidecars emitted / leaf solids | 5 / 5 | **5 / 5** |
| `as1-oc-214.stp` faces recognized | 28 cylinder | **28 cylinder** |
| `as1-oc-214.stp` sidecar bytes changed | — | **0 of 5 files** |
| `Bagger.step` sidecars emitted / leaf solids | 12 / 13 | **12 / 13** |
| `Bagger.step` faces reaching recognition | 0 of 191 | **0 of 191** |
| `Bagger.step` sidecar bytes changed | — | **0 of 12 files** |
| ALICE3 `n_eligible` (surface report) | 36 / 55 | **36 / 55** |
| ALICE3 sidecars emitted | 20 / 55 | **20 / 55** |
| ALICE3 sidecar bytes changed | — | **0 of 20 files** |
| ALICE3 faces recognized | 1180 | **998** |
| — cylinder | 786 | **786** |
| — cone | **358** | **176** |
| — sphere | 36 | **36** |
| ALICE3 recognized faces with measured relative gap > 1e-9 | **182** | **0** |

Every emitted sidecar on every corpus is **byte-identical**. The only thing that moved is that 182
faces which no longer claim to be cones — and every one of them is on `ST2487458_01`, which emits
no sidecar either way, so no representation changed anywhere. That is the intended shape of this
change: it removes a wrong answer that was not yet being used, before the §2 work makes it usable.

*A comparison that cannot fail is not evidence*, so: the byte-identity is backed by the reason it
holds (Bagger and the fixtures reach recognition 0 times; ALICE3's changed faces are all on a solid
that fails for 742 other reasons), and the self-test's positive controls all still pass at
8e-17…6e-16 relative, which is what shows the recogniser itself was not weakened.

### 5.2 Converter vs OCCT after the fix — and why they now differ by exactly ten faces

| converter | OCCT `ShapeAnalysis_CanonicalRecognition` | faces |
| --- | --- | ---: |
| cylinder | cylinder | 786 |
| cone | cone | 174 |
| sphere | sphere | 36 |
| cone | *declines* | 2 |
| *declines* | plane | 8 |
| declines | declines | 1456 |

**Zero disagreements of the dangerous kind remain** — there is no face left where the converter
claims a quadric that OCCT rejects *and* the measured gap is large. The two residual classes are
both explained, and both are the converter being **stricter**, which is the direction this project
wants:

- the **8 planes** are all on `ST2487195_01` and are flat to only **1.5e-07 of the patch diagonal**
  (3.4e-08 mm absolute — they are, exactly, the source of `Stream_A_CSG.md`'s headline
  `maxCanonicalGap = 3.4e-8 mm`). OCCT accepts them against its own *absolute* 1e-07 mm tolerance;
  this converter declines them at 150× its *relative* bound. A declined face costs coverage —
  `ST2487195_01` has 138 genuinely free-form faces of 182 and never emits — and a wrongly-accepted
  one costs correctness.
- the **2 cones** are geometrically sound (measured relative gap ≤ 1e-09) and are declined by OCCT
  under its absolute-tolerance criterion.

So the honest statement of the census's "1004 quadrics in disguise" for *this* converter is
**996 + 2 = 998 at a relative gap below 1e-9**, with the remaining 8 sitting between 1e-9 and
1.6e-7 and reported rather than taken.

---

## 6. `surface_report.json`: what it now records, and the two numbers that must be quoted together

The report used to answer one question — *are this solid's surfaces individually representable?* —
and ALICE3's answer to it (`n_eligible = 36`) has been quoted as coverage. The number that
describes what is actually written is `emitted = 20`. **Reporting one without the other is how §2's
16-solid gap stayed invisible**, so the report now carries both, per solid and in the summary,
along with the evidence each recognition was accepted on.

| field | per | meaning |
| --- | --- | --- |
| `recognized_gap_cm`, `recognized_gap_relative` | face | the **achieved** gap the acceptance was made on |
| `recognized_counts`, `recognized_max_gap_cm` | solid | which target types, and the worst gap in it |
| `eligible_without_recognition` | solid | the same solid with the pre-pass off |
| `emitted`, `extraction_reasons` | solid | whether a sidecar was written, and if not, why |
| `recognized_max_gap_cm` (by type), `recognized_acceptance_tolerance_relative` | model | the acceptance criterion and what it achieved |
| `n_eligible_without_recognition`, `n_rescued_by_recognition` | model | the **coverage delta** recognition is responsible for |
| `n_emitted`, `n_eligible_but_not_emitted`, `n_emitted_carrying_recognized_faces` | model | what was written, and the §2 gap as a number |

### 6.1 The coverage table

| model | leaf solids | eligible **without** recognition | rescued by recognition | eligible | **emitted** | eligible but declined at extraction | faces recognized | max achieved gap (cm) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| `as1-oc-214.stp` | 5 | **0** | **5** | 5 | **5** | 0 | 28 cylinder | 5.5e-12 |
| `Bagger.step` | 13 | 12 | **0** | 12 | **12** | 0 | — | — |
| `CAD_noETA.stp` | 55 | **15** | **21** | 36 | **20** | **16** | 786 cyl / 176 cone / 36 sph | 3.8e-11 |
| `ST2487728_01` (IRIS) | 21 | 11 | 0 | 11 | **11** | 0 | 786 cyl / 176 cone / 36 sph | 3.8e-11 |

`15 + 21 = 36` reproduces `Stream_A_CSG.md`'s `quadricOnly = 15` and `tier0Rescues = 21`
**exactly**, from the converter's own extractor rather than from the census — a fourth independent
instrument agreeing on the same cells. as1-oc-214 is `0 → 5` **on the eligibility column**, which is
the sense in which the brief's `0/5 → 5/5` is true; on the sidecar column it was 5/5 already (§1).

*(IRIS's leaf-solid count is 21 here against the census's 62 prototypes: the converter dedups by
XCAF definition label and the census by `TShape` identity, so the denominators are not the same
quantity. Its per-face recognition counts are identical to ALICE3's because the two models share
parts, which `Stream_A_CSG.md` §1.6 already records.)*

---

## 7. The outward-normal audit

`Stream_L_ALICE3Defect.md`'s criterion — *no face's outward normal may be antiparallel to the source
face's* — applied to **every face of every emitted sidecar**, before and after this work, with the
existing probes (`probes/faceNormalSamples.py` + `probes/faceNormals.cxx`; no third instrument was
written).

| corpus | sidecars audited | faces checked | **antiparallel** | of which recognized faces |
| --- | ---: | ---: | ---: | ---: |
| `CAD_noETA.stp`, before | 18 of 20 | 1936 | **0** | 87 |
| `CAD_noETA.stp`, after | 18 of 20 | 1936 | **0** | 87 |
| `as1-oc-214.stp`, after | 5 of 5 | 53 | **0** | 28 |
| `Bagger.step` | 12 of 12 | 191 | **0** | 0 |

**115 recognized faces ship inside sidecars across the two models, and none of them is inverted.**

Two things this audit is not:

- It is **not** evidence that recognition "went well" — `Stream_L_ALICE3Defect.md` is explicit that
  `reliable` was not evidence before and is not now. It is one specific check, the only one in the
  project that is sign-sensitive.
- The 2 ALICE3 sidecars that do not load are excluded, not clean: `ST1829909_004` and
  `ST1829909_01` still fail at load on the fixed 1e-06 cm wire-join tolerance (mechanism 3, gaps
  4.00e-06 and 5.41e-06 cm), unchanged by this work and unchanged in count (20 emitted / 18 load,
  before and after). Nothing here touches that.

---

## 8. The control that says how much of this the recognition path is responsible for

`--recognize-surfaces off` is the positive control for the whole coverage claim, and it is decisive
on `as1-oc-214.stp`:

```
$ O2_CADtoTGeo.py as1-oc-214.stp --exact-surfaces auto --recognize-surfaces off
Surface report: 0/5 logical volumes eligible for exact O2BVHSurfaceSolid conversion
Exact-surface extraction (auto): 0/5 leaf solids represented exactly, 5 fall back to tessellation
  emitted 0/5
```

**0/5 with the pre-pass off, 5/5 with it on.** That model has no analytic curved surface stored at
all; every one of its five sidecars exists only because of canonical recognition. It is the
sharpest available demonstration that the mechanism works, and it is the true reading of the
brief's `0/5 → 5/5` — a statement about the *pre-pass*, not about this session's changes.

---

## 9. Acceptance: the whole board, totals and disagreement counts together

### 9.1 The invariant

| | baseline | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 91 cases, green | **92 cases, green** (one added, §9.5) |
| fixtures gate, shipped | 9/9 | **9/9** (9 of 10 leaf solids scored) |
| fixtures tiers | CSG 2 / surface 7 | **CSG 2 / surface 7** |
| Bagger gate, shipped | 12/12 | **12/12** (12 of 13 leaf solids scored) |
| Bagger tiers | CSG 7 / surface 5 | **CSG 7 / surface 5** (+ tessellated 1, unscored) |
| unexplained disagreements, **surface**, fixtures | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| unexplained disagreements, **surface**, Bagger | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| unexplained disagreements, **shape**, fixtures / Bagger | 0/0/0/0 (2 parts) / 0/0/0/0 (7 parts) | **unchanged** |
| unexplained disagreements, **mesh**, fixtures / Bagger | 283/6936/5504/5561 / 418/6721/7973/10299 | **unchanged** |
| `runOracleGate.py --self-test` | 17/17 | **17/17** |
| `o2-bench-detectorsbase-xray --self-test` | 17 checks, 0 failures | **17 checks, 0 failures** |
| `compareGateRuns.py --self-test` | 4/4 injected defects caught | **4/4** |
| `O2_CADtoTGeo.py --self-test` | *(did not exist)* | **18 checks, 0 failures** |

### 9.2 Inertness on the two gated corpora

`compareGateRuns.py` on before/after `gate.json`, full field-by-field:

| corpus | parts | differing fields | of which are absolute workdir `source` paths | **real changes** |
| --- | ---: | ---: | ---: | ---: |
| fixtures | 9 | 38 | 38 | **0** |
| Bagger | 12 | 55 | 55 | **0** |

and the stronger statement underneath it: **every `surfaces_*.bin` is byte-identical** — 9 of 9
fixtures, 12 of 12 Bagger, 5 of 5 as1-oc-214, 20 of 20 ALICE3 (both the plain conversion and the
X-ray DB's `--mesh --csg auto` build). *A comparison that cannot fail is not evidence*, so the
reason it holds is stated as a measurement: **0 of 244 gated faces reach the recognition pre-pass
at all** (`Stream_L_ALICE3Defect.md` §6), and the differ was given a positive control — it catches
4/4 injected defects, including a 1e-8 relative nudge of a single numeric leaf.

### 9.3 X-ray transport, `--beams 96` — `as1-oc-214.stp`

This model is where recognition actually ships (28 of 53 faces in the emitted sidecars), so it is
the sharpest transport test of the recognition path in the corpus. Fibonacci fan, **96 directions**,
raster 12 → 69120 rays, both stepping modes:

| | before | after |
| --- | --- | --- |
| rays identical to OpenCascade, mode (a) shape API | 69104 / 69120 | **69104 / 69120** |
| rays identical to OpenCascade, mode (b) navigator | 69104 / 69120 | **69104 / 69120** |
| **LOST** / displaced / sense (`kind`) | 0 / 0 / 0 | **0 / 0 / 0** |
| extra | 24 | **24** |
| zero / non-advancing / unstick / cap / unterminated / odd / duplicate / parity | all 0 | **all 0** |
| mode (a) vs mode (b) | 0 of 69120 rays disagree | **0 of 69120** |
| parts fully clean | 2 / 5 | **2 / 5** |

The two full reports are **identical after stripping timings and workdir paths** — not "agree
within", identical.

**The 24 `extra` crossings are pre-existing and are a sampling artefact, and that is measured
rather than assumed.** They scale linearly with the ray count at a fixed rate, which a defect
localised to a face would not:

| raster | rays | extra | rate |
| ---: | ---: | ---: | ---: |
| 11 | 58080 | 15 | 2.6e-04 |
| 12 | 69120 | 24 | 3.5e-04 |
| 13 | 81120 | 32 | 3.9e-04 |

`LOST = 0` at every density, worst |Δt| 9e-12…2e-11 cm, and they are on the three parts with
curved walls — consistent with near-tangential grazes where the kernel admits a hit OCCT drops.
Not investigated further; flagged because it is an unexplained disagreement with the oracle on a
shipped model, and it should not be quietly inherited.

### 9.4 X-ray transport, `--beams 96` — ALICE3 `CAD_noETA.stp`

Fibonacci fan, **96 directions**, raster 6 → 3456 rays per part, 18 loadable sidecars of the 20
emitted, 62207 rays, both stepping modes:

| | before | after |
| --- | ---: | ---: |
| rays identical to OpenCascade, mode (a) / mode (b) | 62202 / 62207 each | **62202 / 62207** each |
| **LOST** / displaced / sense (`kind`) | 0 / 0 / 0 | **0 / 0 / 0** |
| extra | 10 | **10** (3 parts: `ST0923290_019` 6, `ST1829909_002` 2, `ST1829909_003` 2) |
| parts fully clean | 15 / 18 | **15 / 18** |
| zero / non-advancing / unstick / cap / unterminated / odd / parity / parityNB | all 0 | **all 0** |
| duplicate crossings | 2 | **2** |
| mode (a) vs mode (b) | 0 of 62208 rays disagree | **0 of 62208** |
| sidecars emitted / loading | 20 / 18 | **20 / 18** |

The two full reports are **identical after stripping timings and workdir paths**, and underneath
that the two X-ray DBs' 20 `surfaces_*.bin` are byte-identical, so they could not have differed.

**A number that moved against the standing baseline, and it is not this work.**
`Stream_L_ALICE3Defect.md` §2.5 records 13814/13822 rays, 18/18 parts clean and every robustness
counter zero, measured with a **three-axis raster at N=16**. This fan sees 15/18 parts clean, 10
`extra` crossings and 2 duplicate crossings on parts that raster calls clean — **before** the change
as well as after. That is `Stream_J_XRay.md` §6.2's warning landing again: a parallel beam is
direction-poor, and "clean" from it means "clean at this sampling". The 18/18 figure should be
quoted with its raster.

**`LOST = 0` here is not a claim that the torus quartic defect (`Stream_M_Quartic.md`) is gone.**
`Stream_L_ALICE3Defect.md` §2.5 measured 14 LOST on `ST2487462_01` with a *three-axis* raster at
N=16; this is a different and much more direction-diverse sampling at a lower transverse density,
and it does not hit that configuration. A benchmark that does not reproduce a known defect has not
retired it — it has failed to sample it, which is exactly the standing warning in
`Tutorial.md` §4.5.

The two sidecars that do not load are unchanged and are `Stream_L_ALICE3Defect.md`'s **mechanism 3**
— a fixed 1e-06 cm wire-join tolerance against declared edge tolerances of 8.6e-05 and 3.1e-04 cm:

```
surfaces_ST1829909_004_...bin: surface 371:  wire edge 0 end does not join the next edge start (gap 4.00e-06 cm, tolerance 1e-06 cm)
surfaces_ST1829909_01_...bin:  surface 1006: wire edge 1 end does not join the next edge start (gap 5.41e-06 cm, tolerance 1e-06 cm)
```

Neither part contains a single B-spline *surface*, so nothing in this stream can have caused or
cured them. **No newly-emitted sidecar fails to load, because no sidecar is newly emitted**: the
count is 20 before and 20 after, and the 20 files are byte-identical.

### 9.5 The one C++ test

The recognition work is entirely converter-side and its controls are `--self-test`. What was added
to `Detectors/Base/test/testBVHSurfaceSolid.cxx` is the **kernel-side contract the converter
measures against**: `_recognized_inner_wall()` decides a recognized quadric's `inner_wall` by
comparing the face's own outward normal against "away from the axis", and that measurement is only
correct if the kernel's convention is the one it assumes. This work multiplies the number of faces
going through that path, so the convention is pinned for the cylinder, the cone **and the sphere** —
`Stream_L_ALICE3Defect.md` §9 records that recognized planes and spheres are exercised by nothing,
and this is now the one thing that exercises the sphere and cone branches' sign.

---

## 10. What this leaves open, in the order it is worth doing

1. **The trim for recognized quadrics — the whole remaining coverage lever (§2).**
   *(The "B-spline fit" prescription in this item is superseded — see
   [`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md) and the box at the head of §2. The
   acceptance-criterion advice in the last two sentences still stands, and its own measurement
   confirms them: the natural bound really is the edge's own declared tolerance, and landing §5
   first really did matter — 287 of §2's "free-form" edges were on `ST2487458_01`'s phantom
   cones.)* 16 ALICE3 solids
   recognize completely and are then declined by one line, and they are the entire difference
   between `n_eligible = 36` and `emitted = 20`. It needs a B-spline **fit** in the recognized
   (phi, other) domain with the achieved 3D deviation of the reconstructed boundary measured
   against the source edge and recorded per edge, because the exact composition is transcendental
   and no exact representation exists in the sidecar's trim vocabulary. Design the acceptance
   before the fit: the natural bound is the edge's *own declared tolerance*, which ALICE3 carries
   (2e-06 … 4e-04 cm) and which the sidecar already stores as `modelTolerance`. **Land the
   acceptance-criterion fix (§5) first** — it is already landed — because without it
   `ST2487458_01`'s 182 phantom cones become emittable the moment the trim can carry them.
   1-2 days with verification.

2. **The 8 planes on `ST2487195_01` (§5.2).** They are flat to 1.5e-07 of the patch diagonal and
   are declined at 1e-09. Whether that is right is a *policy* question about a relative bound, not
   a bug, and it should be decided by measuring what a 1.5e-07 relative surface error does to
   transport rather than by choosing a constant. It is the only place where this converter is
   measurably less permissive than OCCT.

3. **A recognition fixture.** `Stream_L_ALICE3Defect.md` §7 asks for a ladder fixture containing a
   NURBS-encoded cylinder, and it is still missing: the gated corpora reach the recognition
   pre-pass **zero** times, so every statement in this document about the shipped gates is a
   statement about a path they do not exercise. `make_boolean_fixtures.py` builds canonical quadrics
   only; `--self-test`'s controls (built with `BRepBuilderAPI_NurbsConvert`) are exactly the
   ingredients a fixture would need, so this is now a small job rather than a half-day one.

4. **The face-normal criterion as a gate column.** §7 is a hand-run audit over 2180 faces with two
   probes. `Stream_L_ALICE3Defect.md` §9 already asks for `--face-normals` in the harness and a
   gate column; this stream is the second consumer of it, and the second time it had to be run by
   hand.

5. **The `extra` crossings (§9.3, §9.4).** 24 on as1-oc-214 and 10 on ALICE3, `LOST = 0`
   everywhere, scaling linearly with the ray count. Almost certainly near-tangential grazes, but
   "almost certainly" is not this project's standard and they are unexplained disagreements with
   the oracle on shipped models. The per-face attribution probe (`probes/faceAttrib.cxx`) would
   name them in an hour.

6. **Not touched, by instruction and by scope:** the loader's fixed 1e-06 cm wire-join tolerance
   (`Stream_L_ALICE3Defect.md` mechanism 3, still 20 emitted / 18 loading on ALICE3), the torus
   quartic (`Stream_M_Quartic.md`), and anything under `Detectors/Base/src/**`.
