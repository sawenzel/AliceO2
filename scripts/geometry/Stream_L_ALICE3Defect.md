# Stream L — the ALICE3 transport defect, named per face

Date: 2026-08-01. Branch `swenzel/bvhsurfacesolid`. Companion to
[`Stream_J_XRay.md`](Stream_J_XRay.md) §7 (which found it and said so),
[`Stream_F_EdgeIdentity.md`](Stream_F_EdgeIdentity.md) §8 (whose stated caveat this is the first
instance of) and [`Stream_E_Scale.md`](Stream_E_Scale.md) §5 (the quartic defect, here measured at
production scale for the first time).

**In one paragraph.** The X-ray benchmark reported that 4 of 18 loaded ALICE3 parts lose 418
boundary crossings, invert the sense of 236, leave 70 transports unterminated and contradict
`Contains()` on 336 interval midpoints — while all 18 report `reliable` and `navigable` with zero
non-manifold edges. This stream attributes every one of those events to the **source face** that
produced it and finds **three distinct mechanisms, not one**. 404 of the 418 lost crossings, all
236 sense inversions and 68 of the 70 unterminated transports come from **21 faces whose outward
normal is exactly inverted** — NURBS-encoded cylinders whose `inner_wall` flag the converter's
canonical-recognition pre-pass infers from `TopoDS` face orientation, which on a NURBS-encoded
quadric means nothing. The remaining 14 lost crossings and 2 unterminated transports are the
**known, unfixed torus quartic defect**, reproduced here down to the coefficients of a single
`solveQuarticReal` call, at ×1 production scale on real geometry. The 2 sidecars that do not load
at all are a **third**, unrelated defect: a hard-coded 1e-06 cm wire-join tolerance against a model
whose edges are defined to 3e-04 cm. Mechanism 1 is fixed here, in one commit, in the converter;
it is provably inert on both gated corpora, which contain **zero** recognized faces. Mechanisms 2
and 3 are diagnosed, reproduced and scoped, and deliberately not fixed.

---

## 1. The instrument, and its self-checks

The project's rule is to build the thing that *names* the object before theorising. `LOST = 418`
is an aggregate; four instruments were built to turn it into faces.

| probe | what it answers | self-check |
| --- | --- | --- |
| **per-ray transport diff** (`XRayTransport.h`'s own `stepWithShapeApi`, plus the same list pairing `compareLists` uses) | which *reference* crossing was lost / sense-inverted, and what the kernel's own hit list says at that distance | must reproduce the benchmark's per-part counts **exactly**, and report zero on a clean part |
| **per-face attribution** (OCCT `IntCurvesFace_ShapeIntersector::Face(i)`) | which source face produced the crossing the kernel lost | face count must equal the sidecar's surface count (`TopologyExplorer(shape).faces()` is the converter's own iteration, so face *i* is record *i*) |
| **per-face outward normal** (`ComputeNormal(p, nullptr, n)` vs OCCT's oriented normal at interior samples of every face) | which faces have an inverted outward direction — **independently of any ray** | must read 0 on all 14 transport-clean ALICE3 parts and on both gated corpora |
| **surface + trim residual** (the sidecar's own `GetSurfaceRecords()`) | whether a crossing the kernel never produced lies on the trimmed patch at all | the implicit residual must be ~0 on a face that is present |

Every one passed its self-check. In particular the transport diff reproduces
Stream J §7's table part for part — 768/768, 700/768, 696/768, 672/768, 768/768, 760/768 with
LOST 0/68/144/192/0/14 — before anything was concluded from it.

**The normal probe is strictly sharper than the ray attribution**, which is the point of building
it: on `ST0923290_013` the rays name 7 inverted faces and the normal probe names **9**. The other
two are simply not hit by a 768-ray, three-direction raster. An instrument that depends on the
sampling cannot say "no defect", only "no defect here".

---

## 2. Mechanism 1 — an inverted outward normal on a recognized NURBS quadric

### 2.1 What the transport sees

`ST0923290_013`, ray 146, along +x:

```
    cand: 8.41905897(+1)
    ref : 1.72494103(+1) 8.41905897(-1)
    LOST t=1.724941032 (+1)  from t=0   Contains=0
        need dt=1.724941032 | DistOut bvh=8.41905897 loop=8.41905897 | DistIn bvh=1.724941032
        parity-ray hits from there (bvh 2 / loop 2):
           bvh  d=1.724941032  n.d=+0.726  onTrimBoundary=0
           bvh  d=8.41905897   n.d=-0.726  onTrimBoundary=0
```

The kernel **finds both crossings**, at the right distances, on the BVH and the `_Loop` twin
alike. It attaches the wrong sign to both: `dot(n, d) = +0.726` at the entry and `−0.726` at the
exit. `DistFromOutside` accepts only hits with `dot(n, d) < 0`, so it returns the *exit* as the
entry; `DistFromInside` accepts only the others. The transport enters at the far wall, calls it an
entry (`kind`), and never leaves (`unterm`).

`Safety()` at every one of the 418 lost crossings is between **8.9e-16 and 5.6e-12 cm**. The
boundary is exactly where OpenCascade says it is. **Nothing is mis-trimmed and nothing is missing.**

`Contains()` is a *parity* count and is blind to the sign, which is why it stays right and why the
`parity` counter fires 336 times: `Contains()` and the crossing list contradict each other because
the list is wrong and `Contains()` is not.

### 2.2 The per-face table

Kernel outward normal vs OpenCascade's oriented outward normal, at 4 interior samples of every
face of every ALICE3 part. `align` is `dot(n_kernel, n_OCCT)`.

| part | faces | inverted | which | stored type | `TopoDS` orientation | face tol | trim curves | edge tol max | degenerate | worst align |
| --- | ---: | ---: | --- | --- | --- | ---: | --- | ---: | ---: | ---: |
| `ST0923290_013` | 20 | **9** | 1, 2, 3, 5, 7, 8, 10, 11, 13 | `bspline` → recognized **cylinder** | **REVERSED** (9/9) | 1.0e-07 | 4 per face, b-spline (3 faces carry one line) | 3.30e-06 | 0 | **−1.0000** |
| `ST0923290_018` | 48 | **6** | 7, 8, 9, 39, 40, 41 | `bspline` → recognized **cylinder** | **FORWARD** (6/6) | 1.0e-07 | 4 per face, b-spline | 2.18e-06 | 0 | **−1.0000** |
| `ST0923290_019` | 44 | **6** | 13, 14, 15, 16, 17, 18 | `bspline` → recognized **cylinder** | **FORWARD** (6/6) | 1.0e-07 | 4 per face, b-spline | 2.42e-06 | 0 | **−1.0000** |
| `ST0923290_012` | 10 | 0 | — | 4 `bspline` → cylinder, FORWARD | — | | | 3.29e-07 | 0 | +1 |
| `ST0923290_020` | 37 | 0 | — | 17 `bspline` → cylinder, REVERSED | — | | | 3.26e-06 | 0 | +1 |
| `ST2487462_01` | 80 | **0** | — | — | — | | | 5.28e-05 | 0 | +1 |
| the other 12 loaded parts | 1–965 | **0** | — | — | — | | | 1e-07 … 4.3e-04 | 0–20 | +1 |

Read it:

- **`align` is exactly −1.0000.** Not a geometry error with a sign consequence: a pure sign
  inversion. The surface, its frame and its trim are all right.
- **Every inverted face is a `bspline` face recognized as a cylinder.** No natively-analytic face
  anywhere is inverted, on either corpus.
- **The `TopoDS` orientation does not predict it.** `_013`'s nine are all REVERSED and `_018`/`_019`'s
  twelve are all FORWARD, and `_012`'s four FORWARD and `_020`'s seventeen REVERSED recognized
  cylinders are all *correct*. Whichever way the flag points, it is right about half the time.
- **The trim-curve kinds and the edge tolerances are not the discriminator either.** These parts
  carry the *tightest* edge tolerances on ALICE3 (2.2e-06 to 3.3e-06 cm); the parts with the
  loosest (4.3e-04 cm) are clean.
- Ray attribution agrees face for face: the crossings lost on faces 7/8/9/39/40/41 of `_018` are
  exactly the ones flagged `FLIPPED`, and the losses attributed to *planes* and to REVERSED
  b-spline faces on the same rays carry no `FLIPPED` flag — they are **consequential**, lost after
  the stepping state was already inverted.

### 2.3 The cause, in the converter

`recognize_and_extract_face` (`scripts/geometry/O2_CADtoTGeo.py`) set

```python
inner_wall = face.Orientation() == TopAbs_REVERSED  # same convention as the native extractors:
# the recognized surface's canonical normal is always "away from the axis/center", exactly
# like OCC's own canonical quadric normal, regardless of the recognized frame's axis/refu sign.
```

The premise in that comment is true for a face whose **stored** surface is an OCC canonical
quadric: OCC's canonical quadric normal does point away from the axis, so FORWARD means "away" and
REVERSED means "towards". It is false for a **NURBS-encoded** quadric, which is exactly what the
recognition pre-pass exists to handle. There the stored surface is a B-spline; its normal is
`∂S/∂u × ∂S/∂v`, whose direction is a choice the exporter made when it ordered the parametric
directions, and FORWARD/REVERSED is relative to *that*. The flag therefore carries no information
about which side of the recognized cylinder is outside, and reading it as if it did inverts an
arbitrary subset of the recognized faces.

Reproduced directly, by running the real extractors over the real faces:

```
  face    1 bspline    rev=1 recognized->cylinder   -> record type=cylinder inner_wall=1 align=-1.0000  <== INVERTED
  face    2 bspline    rev=1 recognized->cylinder   -> record type=cylinder inner_wall=1 align=-1.0000  <== INVERTED
  ...
  face    4 plane      rev=1 native                 -> record type=plane                 (consistent)
```

**Why closure cannot see it.** Edge identity counts each source edge's occurrences and their
relative sense; a *global* sign flip on a face cancels in every count it makes. `maxRimIsolation`,
`unmatchedRimLength` and `maxSharedEdgeDeviation` are all distances and are sign-blind too. So
`ST0923290_013` reports `reliable`, `navigable`, 0 boundary / 0 non-manifold / 0 reversed edges —
and a transported particle walks through nine of its faces. This is
`Stream_F_EdgeIdentity.md` §8's caveat, *"identity certifies that the source B-rep's topology
survived conversion; it does not certify each face's geometry"*, measured for the first time.

### 2.4 The fix, and the four independent numbers that confirm it

The fix measures what the flag cannot say: take the face's own outward normal — OCC's rule,
`(∂S/∂u × ∂S/∂v)` negated for a REVERSED face — at every recognition sample, and compare its sign
against the recognized quadric's canonical outward direction there. `_recognized_inner_wall()`, one
function, ~25 lines, called from the one place that had the wrong inference. It returns `None`
when the samples do not decide, and the orientation flag is then still all there is.

Before and after, on the three offending parts, reconverted from the same `.brep`:

| part | | faces with inverted normal | rays identical | LOST | kind | unterm | `Capacity()` vs OCCT |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ST0923290_013` | before | **9 / 20** | 700 / 768 | 68 | 68 | 68 | **−3.263e-01** |
| | after | **0 / 20** | **768 / 768** | **0** | **0** | **0** | **+3.7e-07** |
| `ST0923290_018` | before | **6 / 48** | 696 / 768 | 144 | 72 | 0 | **+1.979e-02** |
| | after | **0 / 48** | **768 / 768** | **0** | **0** | **0** | **−1.9e-07** |
| `ST0923290_019` | before | **6 / 44** | 672 / 768 | 192 | 96 | 0 | **+2.733e-02** |
| | after | **0 / 44** | **768 / 768** | **0** | **0** | **0** | **+4.3e-07** |
| `ST0923290_012` | before / after | 0 / 0 | 768 / 768 both | 0 | 0 | 0 | unchanged |
| `ST0923290_020` | before / after | 0 / 0 | 768 / 768 both | 0 | 0 | 0 | unchanged |

Four instruments, three of them not used to find the defect, all agree: the per-face normal probe,
the transport crossing lists, the divergence-theorem `Capacity()` (a **−32.6 %** volume error on
`ST0923290_013` becomes 3.7e-07), and the chord-integral volume. The two clean members of the same
family are bit-for-bit unchanged, which is the negative control.

*(The residual ~4e-07 on `Capacity()` is the rounding of the reconversion path used for this
measurement — it reads the dumped `.brep`, already in cm, where the converter scales from mm — not
of the fix. The gate reconverts from the STEP and is unaffected.)*

### 2.5 The whole model, reconverted from the STEP and re-transported

Not the three parts in isolation: `runXRayBench.py --model CAD_noETA.stp --raster 16
--representations surface`, end to end, before and after.

| | before | after |
| --- | --- | --- |
| `surface` (a) shape | 13578 / 13822 rays identical | **13814 / 13822** |
| `surface` (b) nav | 13580 / 13822 | **13816 / 13822** |
| parts fully clean | 14 / 18 | **17 / 18** |
| **LOST** | 418 | **14** |
| extra / displaced | 0 / 0 | **0 / 0** |
| **kind** (sense inverted) | 236 | **0** |
| **unterminated** | 70 | **2** |
| **parity** (`Contains()` vs the list) | 336 (`parityNB` 0) | **0** (`parityNB` 0) |
| zero / non-advancing / stalled / capped steps | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| mode (a) vs mode (b) | 6 of 13824 rays | 6 of 13824 rays |
| sidecars emitted / loaded | 20 / 18 | 20 / 18 |

Everything mechanism 1 explains is gone, including all 236 sense inversions and every one of the
336 parity contradictions. What remains — 14 LOST, 2 unterminated, 6 mode-(a)-vs-(b) disagreements,
all on `ST2487462_01` — is mechanism 2, untouched by instruction.

---

## 3. Mechanism 2 — the torus quartic, at ×1, on real geometry

`ST2487462_01` has **zero** inverted normals. Its 14 lost crossings are a different thing, and the
instrument says so without being asked: 8 of them are positions where the kernel produced **no hit
at all**, where mechanism 1 always produces one.

### 3.1 The faces

| face | type | orientation | face tol | trim curves | edge tol max | events |
| ---: | --- | --- | ---: | --- | ---: | --- |
| 32, 38, 40, 47, 48, 54, 55, 58 | **torus** | 3 FORWARD / 5 REVERSED | 1.0e-07 | 2–3 circles + 1–2 b-splines | 1.9e-05 … 3.9e-05 | **no kernel hit** (8 primary losses) |
| 31, 57 | torus | REVERSED | 1.0e-07 | circles + b-spline | 2.0e-05 | lost, hit present — consequential |
| 45, 46 | cone | FORWARD | 1.0e-07 | 2 circles + 2 b-splines | 3.3e-07 | lost, hit present — consequential |
| 6, 76, 79 | plane | | 1.0e-07 | circles | 1.0e-07 | matched |

All eight primary losses are **toroidal**. The lead is confirmed by measurement, not assumed; had
they been planes or cylinders it would have been excluded.

### 3.2 It is not the trim

At each missing crossing, evaluated against the sidecar's own record:

```
missing crossing at (-0.792770358, -4.3883125, -3.0380625)  t=375.33923 kind=+1
  torus surface  47: |implicit| = 1.679e-14 cm   R=5.3 r=0.1
      phi=+0.605545 in [+0.560644, +0.114477] -> INSIDE | tube=-1.188214 in [+4.712389, +1.570796] -> INSIDE
```

The point is on the untrimmed torus to **1.7e-14 cm** and inside **both** parameter windows. The
patch is there, the trim admits it, and the ray/surface solver returns nothing.

### 3.3 Reduced to one `solveQuarticReal` call

Walking the ray's *start point* along the ray — the same physical intersection, a differently
conditioned quartic — gives an erratic, knife-edge pattern, not a smooth degradation:

| start offset (cm) | remaining (cm) | crossing found? |
| ---: | ---: | --- |
| 0 | 375.34 | **missing** |
| 1 | 374.34 | **missing** |
| 10 | 365.34 | found |
| 100 | 275.34 | **missing** |
| 200 | 175.34 | **missing** |
| 300, 350, 370, 374, 375 | 75.3 … 0.34 | found |

The real coefficients, at offset 0:

```
a4=1  a3=-1501.7280000044018  a2=845808.25396968238  a1=-211752288.545858  a0=19882619385.616932
   ->  p=113.1342207   q=-5.960464478e-08   r=-0.9737548828
       resolvent = 7.105427358e-15   <== guard `resolvent > kTolerance` FAILS: no roots returned
```

The true roots are 375.3392298 and 375.5247703. `solveQuarticReal` returns **none**.

**The mechanism, exactly.** This quartic is *biquadratic*: its true `termQ` is 0 (the ray is
perpendicular to the torus axis). `termQ` is computed as `d − b·c/2 + b³/8` from terms of magnitude
~1e8, so its floating-point value is not 0 but **−6.0e-08** — and the branch test is
`std::abs(termQ) <= kTolerance` with `kTolerance = 1e-9 cm`. The symmetric quartic is therefore
misrouted into the resolvent-cubic branch, where the resolvent evaluates to 7.1e-15, the
`resolvent > kTolerance` guard rejects it, and the function returns an empty list. Both guards are
**absolute lengths applied to quantities of dimension L³ and L²**, which is precisely the defect
`Stream_E_Scale.md` §5 records; whether the cancellation residual lands above or below 1e-09
depends on the coefficient magnitudes, i.e. on the ray's lever arm — hence the erratic table.

**This is the same defect, at ×1.** Stream E measured it at ×0.1 on a fixture and concluded "fix it
before ALICE3". The trigger is not the model's scale: it is the ratio of the ray-origin distance to
the tube radius, here **375 cm / 0.1 cm = 3750**. On a real detector model, where the raster window
is the part's bounding box and the toroidal features are fillets and knuckles of ≤1 cm, that ratio
is normal.

Not fixed here, by instruction. `Stream_E_Scale.md` §5.3 scopes it; §5 below adds what this
measurement contributes to that scope.

---

## 4. Mechanism 3 — the two sidecars that do not load

Distinct from both, and confirmed rather than assumed:

```
surfaces_ST1829909_004_...bin: surface 371:  wire edge 0 end does not join the next edge start (gap 4.00e-06 cm, tolerance 1e-06 cm)
surfaces_ST1829909_01_...bin:  surface 1006: wire edge 1 end does not join the next edge start (gap 5.41e-06 cm, tolerance 1e-06 cm)
```

| part | face | type | trim curves | that face's edge tolerances | gap | loader tolerance |
| --- | ---: | --- | --- | ---: | ---: | ---: |
| `ST1829909_004` | 371 | cylinder | 5 b-spline pcurves | 1e-07 … **8.63e-05** | 4.00e-06 | 1e-06 |
| `ST1829909_01` | 1006 | cylinder | 5 b-spline pcurves | 1e-07 … **3.06e-04** | 5.41e-06 | 1e-06 |

Neither part contains a single B-spline *surface*, so **mechanism 1 cannot apply to them**; neither
gap is on a torus. The gaps are 14× to 57× *smaller* than the tolerance those very edges declare —
the source data is self-consistent, and the loader's fixed absolute 1e-06 cm join tolerance is what
rejects it. Same family as the argument that retired the geometric closure band
(`Stream_F_EdgeIdentity.md` §1.2: ALICE3 edge tolerances are 4700×–43000× Bagger's), different
defect, different symptom.

---

## 5. Leads that died, and they should be recorded as such

Every one of these was a plausible, documented candidate. All were tested.

1. **"ALICE3's loose edge tolerances (4.7e-04 … 4.3e-03 cm) are what breaks it."** *Refuted for
   mechanisms 1 and 2.* The four offending parts carry edge tolerances of **2.2e-06 to 5.3e-05 cm**
   — among the tightest on the model — while `ST1829909_002` and `_003`, at 4.3e-04 and 4.2e-04,
   are clean to 768/768 rays. It is confirmed only for mechanism 3, where it is the whole story.
2. **"`ST2487462_01`'s max shared-edge deviation (6.8e-05) exceeds its declared tolerance
   (5.3e-05), and it still reports reliable."** True, and worth acting on separately — but *not*
   the cause of its 14 lost crossings, which are a torus root solve that never returns.
3. **"ALICE3's degenerate edges exercise an untested path."** *Refuted, and the premise it rested
   on is now out of date.* `Stream_F_EdgeIdentity.md` §8.3 says the degenerate-edge path has never
   been exercised on a converted part. It has been, since this conversion:
   **`ST1829909_002` (12), `ST1829909_003` (20) and `ST2487459_01` (6)** carry degenerate edges,
   load, close as `reliable`, and are **768/768 transport-clean**. The four offenders have **zero**
   degenerate edges. The path works; it is not the defect and it is no longer untested.
4. **"A mis-trimmed face."** *Refuted by the sharpest available measurement.* `Safety()` at every
   one of the 418 lost crossings is ≤ 5.6e-12 cm, and the eight torus positions sit on their
   untrimmed surface to 1e-14 cm and inside both parameter windows. No trim rejects a point it
   should keep anywhere in this data set.
5. **"It is a scale effect — 9266 faces, metre coordinates."** *Refuted.* The 965-patch
   `ST1829909_002` and the 236-patch `_003` are clean to 2.5e-13 and 6.9e-13 cm. The defect count
   does not correlate with part size, face count or coordinate magnitude; it correlates with
   *recognized NURBS cylinders* and with *toroidal faces at a long lever arm*.
6. **A documentation discrepancy, unrelated but worth correcting.**
   `o2-bench-detectorsbase-xray --self-test` prints and passes **17** checks, not the 19 that
   `Stream_J_XRay.md` §1/§6.4 and `Tutorial.md` §5.7 state. Measured before any file on this
   branch was touched, so it is pre-existing. `xrayOracle.py --self-test` still reports its 11.

---

## 6. Why the shipped corpora are unaffected, and when they would stop being

**Mechanism 1 cannot fire on either gated corpus, and this is a measured statement, not an
argument.** `recognize_and_extract_face` runs only for a face whose *stored* surface type has no
direct extractor. Over all 21 gated parts — 9 fixtures and 12 Bagger, **244 faces** — the
recognition pre-pass is invoked **zero** times:

| corpus | parts | faces | faces reaching recognition |
| --- | ---: | ---: | ---: |
| fixture ladder | 9 | 53 | **0** |
| Bagger | 12 | 191 | **0** |
| ALICE3, 20 emitted sidecars | 20 | 3708 | **87** (all `bspline`, recognized as cylinders, all on the `ST0923290` family) |

So the fix is inert on both gates by construction, and the gates cannot regress on it — which is
also the reason the defect survived to ALICE3: **the shipped corpora contain no NURBS-encoded
quadric at all**, and it is exactly and only the NURBS-encoded quadric that the flag lies about.
It switches on the moment a model is exported by a system that writes quadrics as NURBS, which is
what ALICE3's exporter does for 1004 of its 2377 "B-spline" surfaces (`Stream_A_CSG.md`), and it
would have hit **Tier 0** — the highest-value item on the board — head-on, since Tier 0 is precisely
"decode NURBS-encoded quadrics" and would have multiplied the affected face count.

**Mechanism 2 does not switch on at a scale; it switches on at a ratio.** Stream E measured it at
model scale ×0.1; the trigger measured here is (ray-origin distance) / (tube minor radius) ≈ 3750
at ×1. Bagger has 2 toroidal faces in 288 and a 30 cm bounding box; the ladder's
`torus_union_cyl` has a tube radius of order 1 cm in a 5 cm box — ratios of order 10, where the
`termQ` cancellation residual stays under 1e-09 and the branch test happens to be right. **It will
bite the shipped corpora as soon as anything places them in a world large enough to fire rays from
far away**, which an assembly-level transport (`Stream_J_XRay.md` §9) does by definition. Say this
loudly: the guard is not conservative, it is silently non-conservative, and the failure is a
missing wall, not a wrong distance.

---

## 7. Reproducers

**Mechanism 1 — full, and it is now a regression the gate can carry.** No fixture needed: the
`.brep` of `ST0923290_013` plus the per-face normal probe is a two-line check, and the criterion
is corpus-independent —

> for every face of a converted solid, the kernel's outward normal at an interior point must not be
> antiparallel to OpenCascade's oriented outward normal there.

That is the **checkable criterion for face geometry** that `Stream_F_EdgeIdentity.md` §8 says edge
identity does not provide. It is cheap (one `ComputeNormal` and one `Safety` per sample), it needs
no rays, it is sign-sensitive where every existing closure measurement is sign-blind, and it caught
2 faces that a 768-ray raster did not. **Recommended as a harness flag (`--face-normals`) and as a
gate column.** A ladder fixture would need a STEP file containing a NURBS-encoded cylinder;
`make_boolean_fixtures.py` builds canonical quadrics only, so producing one means converting a
cylinder to a B-spline surface (`GeomConvert`/`BRepBuilderAPI_NurbsConvert`) — half a day, and it
would also be the first fixture that exercises the recognition path at all, which nothing does
today.

**Mechanism 2 — minimal, and it needs no CAD.** A single call:

```
solveQuarticReal(1.0, -1501.7280000044018, 845808.25396968238,
                 -211752288.545858, 19882619385.616932)
```

must return roots at 375.3392298 and 375.5247703 and today returns **an empty vector**. That is a
five-line, pure-arithmetic regression test for `testBVHSurfaceSolid.cxx` and it belongs in the same
commit as the quartic fix — **not before it**, because a test that fails on landing is not a test.
It is recorded here so that whoever takes `Stream_E_Scale.md` §5.3 has a production-scale case to
land it against, alongside Stream E's ×0.1 fixture ray.

**Mechanism 3 — the two `.brep` files themselves**, with the loader's error message naming the
offending surface index (371 and 1006).

---

## 8. What was changed, and what was not

**Changed:** `scripts/geometry/O2_CADtoTGeo.py` — `_recognized_inner_wall()` (new), the two
recogniser returns now also carry the sampled normal field `N`, and `recognize_and_extract_face`
consults the measurement instead of the orientation flag. Nothing else.

**Not changed, deliberately:** the quartic guards (scheduled separately, `Stream_E_Scale.md` §5.3);
the loader's wire-join tolerance; anything under `Detectors/Base/src/**`; the harness; the gate;
`XRayTransport.h`; the benchmark. No new surface type, no new recognition, no Tier 0.

### Acceptance, before and after — totals and disagreement counts together

| | before | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 86 cases, green | **86 cases, green** |
| fixtures gate, shipped | 9/9 | **9/9** |
| fixtures gate, surface (historical) | 9/9 | **9/9** |
| fixtures scored | 9 of 10 leaf solids | **9 of 10** |
| Bagger gate, shipped | 12/12 | **12/12** |
| Bagger gate, surface (historical) | 9/12 | **9/12** |
| Bagger scored | 12 of 13 leaf solids | **12 of 13** |
| unexplained oracle disagreements, **surface**, fixtures | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| unexplained oracle disagreements, **surface**, Bagger | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| disagreements, **shape**, fixtures / Bagger | 0/0/0/0 | **0/0/0/0** |
| disagreements, **mesh**, Bagger | 19 / 1011 / 922 / 573 | **19 / 1011 / 922 / 573** |
| `runOracleGate.py --self-test` | 17/17 | **17/17** |
| `o2-bench-detectorsbase-xray --self-test` | 17 checks, 0 failures | **17 checks, 0 failures** |

`gate.json` was compared field by field, with `timing*`, `*Seconds`, `nsPerCall` and the absolute
`source` paths stripped (the two runs use different workdirs), on both corpora, before and after:

| corpus | leaf fields | removed | **changed** | added |
| --- | ---: | ---: | ---: | ---: |
| fixtures (9 parts) | 8487 | 0 | **0** | 0 |
| Bagger (12 parts) | 15140 | 0 | **0** | 0 |

and, the stronger statement, **all 21 `surfaces_*.bin` are byte-identical** between the before and
after conversions — 9 of 9 fixtures, 12 of 12 Bagger. The converter emits the same bytes it did.

That is the expected outcome, because the recognition pre-pass is never reached on either corpus
(§6) — and *a comparison that cannot fail is not evidence*, so two things back it up: the reason it
cannot fail is a measurement (0 of 244 faces reach recognition), and the differ was given a
positive control (one numeric leaf nudged by 1) and reported exactly that one change on both
corpora.

---

## 9. What this leaves open

* **The quartic.** Diagnosed to the guard and reproduced to five coefficients; not fixed. The fix
  is to make both branch tests **relative** to the coefficient scale — e.g.
  `|termQ| <= eps * max(|b|³, |b·c|, |d|)` and the resolvent likewise — which is what
  `Stream_E_Scale.md` §5.3 already scopes. Cost: the arithmetic is an hour; the risk is that it
  moves torus results on the shipped corpora, so it needs both gates and the ladder X-ray, and
  `torus_union_cyl`'s 2 parity mismatches at 262144 rays are the natural acceptance target. Half a
  day with verification.
* **The loader's wire-join tolerance** (mechanism 3). It should scale with the model's declared
  tolerance, which the sidecar already carries (`modelTolerance`), instead of being a fixed 1e-06
  cm. That would take the ALICE3 sidecar count from 18 to 20. Small, but it changes an acceptance
  criterion, so it needs its own before/after on both corpora. Half a day.
* **The face-normal criterion should become a gate column** (§7). It is the missing half of
  `Stream_F_EdgeIdentity.md` §8's honest two-number summary of a part: *closed by topology*, *the
  faces meet to X cm* — and now *every face faces outward*.
* **`GetNavigationReliability` still cannot express this.** A part with nine inverted faces reports
  `reliable` and `navigable`. Whatever the fallback policy ends up being, "reliable" must stop
  meaning "the topology parsed".
* **Recognized planes and spheres are untested.** All 87 recognized faces on ALICE3 are cylinders.
  The plane branch derives its frame from the sampled normal itself and is correct by construction;
  the sphere and cone branches now go through the same measurement as the cylinder, but nothing in
  reach exercises them.
* **The raster is thin.** 768 rays in three directions per part. `Stream_J_XRay.md` §6.2 is the
  standing warning: `--beams` is what finds numerical defects, and it was not used for the ALICE3
  scan. The clean parts should be read as "no defect at this density".
