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
