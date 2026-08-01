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
