# Stream Q — the ellipse boundary on a planar face

Date: 2026-08-02. Branch `swenzel/bvhsurfacesolid`. **Converter-only change; nothing was built, no
C++ was touched.** Companion: [`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md) §5, whose
measurement this acts on, and [`NEXT.md`](NEXT.md) item 3's "smallest first step".

---

## 0. The answer first

`extract_planar_face` declined a planar face whose boundary carried an **ellipse**. That single line
was the entire reason two solids fell back to a mesh. The fix is **one guard**: the sidecar format
already carries the curve exactly, because *a conic **is** a rational quadratic B-spline* and the
sidecar's `bspline` record is rational. This is `Stream_O`'s outcome (1) — the format expresses it
and only the Python guard was wrong. No kernel change, no new curve kind, no fit.

| | before | after |
| --- | --- | --- |
| ladder fixtures | 9/9 scored of 10 leaf solids | **10/10 scored, 10/10 pass**, `surface` 0/0/0/0 |
| `Bagger.step` | 12/12 shipped, 1 unscored (`Bucket`, mesh) | **13/13 shipped, 13/13 pass**, `surface` 0/0/0/0 |
| Bagger cascade | CSG 7, exact surfaces 5, **tessellated 1** | CSG 7, exact surfaces 6, **tessellated 0** |
| `O2_CADtoTGeo.py --self-test` | 26 checks, 0 failures | **36 checks, 0 failures** (18 + 8 + **10**) |

Achieved trim deviation on the newly-accepted faces, measured against the CAD 3D boundary curve
both ways (§4): **7.9e-15 cm** on `oblique_cut_cyl` f#1 and **2.6e-14 / 2.2e-12 cm** on the four
`Bagger/Bucket` ellipse arcs. The ellipse is never the limiting term on either face.

---

## 1. The reproduction, before anything was changed

`Stream_O_ImplicitTrims.md` §5 was run by another agent. Reproduced here from the shipped
converter's own `--surface-report`, on the unmodified file:

```
oblique_cut_cyl   3 faces: 1 cylinder + 2 plane
  f#0 cylinder  wire {ellipse 1, line 2, circle 1}   supported, trim_kind "general"
  f#1 plane     wire {ellipse 1}                     NOT supported:
                "plane with unsupported boundary curves: ['ellipse']"
  f#2 plane     wire {circle 1}                      supported
  -> emitted: false
     "planar boundary edge is a ellipse curve (only line/circle/bspline supported)"

Bagger/Bucket   97 faces: 69 plane + 22 cylinder + 4 sphere + 2 torus
  unsupported faces: exactly two -- f#4 and f#31, both planes, both
                     "plane with unsupported boundary curves: ['ellipse']"
  -> emitted: false, same single reason
```

Both claims hold, including the surprising one. **`Bucket` has zero faces without an analytic
surface**: its 4 spherical and 2 toroidal faces all extract. Two ellipse-bounded planes were the
whole story.

**This invalidates [`Tutorial.md`](Tutorial.md) §6**, which attributes `Bucket`'s tessellated
fallback to "97 faces, 4 spherical + 2 toroidal". That explanation was wrong; the spheres and the
tori were never the problem. `NEXT.md`'s "`oblique_cut_cyl` has no sidecar and never has" and
`CodeReview_Fable.md`'s row for it are now out of date in the same way.

---

## 2. Which of the three outcomes this is, and why it is not an approximation

**Outcome (1): the format already expresses it.**

The planar sidecar wire takes `line`, `arc` and `bspline` segments, and the `bspline` record carries
**weights** — it is a rational B-spline. Every conic has an exact rational quadratic B-spline form,
and OCCT's `GeomConvert::CurveToBSplineCurve(..., Convert_TgtThetaOver2)` writes precisely that
form. The converter already called it, from `_planar_bspline_edge_params`, for B-spline and Bezier
boundary edges; that helper reads the edge's 3D curve and is indifferent to which conic it is. The
plane projection is an **isometry**, so the control poles land in the plane's (u, v) frame unchanged.

The quadric path had reached the same conclusion long ago and its comment says so: `ellipse` is in
`_SUPPORTED_QUADRIC_CURVES` and `_quadric_trim_wire` converts circle/ellipse/Bezier/B-spline
pcurves through the identical route. The planar vocabulary was simply never widened to match.

So the change is:

```python
_SUPPORTED_PLANAR_CURVES = {"line", "circle", "ellipse", "bspline", "bezier"}
...
if gt not in (GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse,
              GeomAbs_BSplineCurve, GeomAbs_BezierCurve):
```

and nothing else. Ellipses fall into the existing B-spline branch, which then runs the existing
canonical pre-pass — so an "ellipse" whose two semi-axes happen to be equal still comes out as the
cheap exact `arc`, not as a heavier B-spline.

**Hyperbola and parabola stay declined.** Both are also exactly writable as (rational) B-splines,
so this is not a representational limit — it is the standing bargain: no corpus in this project
contains one, so nothing has measured one, and the project does not ship a representation it has
not measured. Both are negative controls in the self-test (§3).

---

## 3. The checks, and the order they were written in

`O2_CADtoTGeo.py --self-test` grew a third block, `run_planar_trim_self_test()`. It was written
**before** the guard was changed and watched to fail for the stated reason:

```
Planar trim vocabulary: the ellipse boundary (Stream_Q_EllipseTrim.md)
  [ok ] the oblique cut of a cylinder really does produce an ellipse-bounded planar face -- found
  [FAIL] an oblique planar cut of a cylinder is accepted -- declined: planar boundary edge is a
         ellipse curve (only line/circle/bspline supported)
  [FAIL] an ellipse ARC boundary (the Bucket case) is accepted -- declined: ...
  [ok ] a hyperbola boundary edge is declined
  [ok ] a parabola boundary edge is declined

5 checks, 2 failure(s)
```

and after:

```
  [ok ] the oblique cut of a cylinder really does produce an ellipse-bounded planar face -- found
  [ok ] an oblique planar cut of a cylinder is accepted -- accepted
  [ok ] the ellipse is stored as a B-spline segment -- segments: ['bspline']
  [ok ] it is a RATIONAL B-spline -- the exact conic form, not a polynomial fit -- weight spread 0.500
  [ok ] the stored trim reproduces the CAD boundary at machine precision
        -- max deviation 7.94e-15 cm = 1.48e-15 patch diagonals
  [ok ] an ellipse ARC boundary (the Bucket case) is accepted -- accepted
  [ok ] the ellipse arc's stored trim reproduces the CAD boundary at machine precision
        -- max deviation 1.13e-15 cm = 2.59e-16 patch diagonals
 the deviation instrument must be able to return a large number
  [ok ] a deliberately wrong conic is caught by the same measurement
        -- max deviation 2.48e-01 cm = 5.67e-02 patch diagonals
 negative controls: a boundary that is not an exactly-writable conic is still declined
  [ok ] a hyperbola boundary edge is declined -- planar boundary edge is a hyperbola curve
        (only line/circle/ellipse/bspline supported)
  [ok ] a parabola boundary edge is declined -- ...

10 checks, 0 failure(s)
```

Total: **36 checks, 0 failures** (18 recognition + 8 placement + 10 planar trim).

Four things in there are deliberate, and each exists because a cheaper version would have been
worthless:

1. **The positive control is the real fixture**, built in-process: a cylinder cut by a box rotated
   60° about x. The face is a genuine B-rep face from `BRepAlgoAPI_Cut`, not a hand-built ellipse.
2. **The deviation is measured both ways** — every CAD sample projected onto the stored curve, and
   every stored sample projected back onto the CAD 3D curve. One-way agreement is not agreement: a
   stored curve covering half the CAD edge is close to it at every one of *its own* points.
3. **The instrument is shown capable of a large answer.** The same measurement applied to a
   deliberately squashed conic returns 2.48e-01 cm. A measurement that only ever sees small numbers
   on the population it is arguing about has not been controlled.
4. **The rationality check.** A *polynomial* B-spline fit of an ellipse would also pass a loose
   deviation test; asserting the weights are not all equal proves the stored object is the conic's
   exact form.

One implementation note that is easy to get wrong: `GeomAPI_ProjectPointOnCurve` returns only the
*perpendicular* extrema and finds **nothing** when the nearest point is an end of the curve, which
is exactly what happens on the deliberately-wrong control. The helper falls back to the endpoint
distance, which is an **upper** bound on the true distance — it can over-report but never
under-report, so it cannot turn an approximation into a pass.

**No C++ test was added and nothing was built.** The kernel is unchanged; the rational-B-spline
planar trim path it exercises is not new to the kernel, and it is now covered end to end by two
gated corpora and 5.1M X-ray rays (§5). Adding an unbuilt, unrun C++ case would have been weaker
evidence than what is below, not stronger.

---

## 4. Achieved deviation, per newly-accepted face, on the real corpora

Measured with the same two-way instrument, 513 samples per edge, on the converter's own leaf solids:

| solid | face | boundary | max deviation | / patch diagonal | patch diagonal |
| --- | --- | --- | ---: | ---: | ---: |
| `oblique_cut_cyl` | f#1 | 1 edge, closed ellipse | **7.944e-15 cm** | 1.480e-15 | 5.3666 cm |
| `Bagger/Bucket` | f#4 | 10 edges, 8 line + 2 ellipse arc | 2.887e-11 cm | 2.900e-12 | 9.9565 cm |
| `Bagger/Bucket` | f#31 | 10 edges, 8 line + 2 ellipse arc | 2.887e-11 cm | 2.900e-12 | 9.9565 cm |

`oblique_cut_cyl` is at machine precision: 1.5e-15 patch diagonals is one double ulp.

`Bucket` is **not**, and the localisation says the ellipse is not why. Per segment:

```
Bucket f#4                                   Bucket f#31
  CAD line    -> line     0.000e+00 cm         CAD line    -> line     7.105e-16 cm
  CAD line    -> line     2.887e-11 cm  <--     CAD line    -> line     2.543e-14 cm
  CAD line    -> line     8.360e-12 cm         CAD line    -> line     2.155e-12 cm
  CAD ellipse -> bspline  2.605e-14 cm         CAD line    -> line     5.217e-12 cm
  CAD line    -> line     2.543e-14 cm         CAD line    -> line     1.237e-11 cm
  CAD ellipse -> bspline  2.172e-12 cm         CAD ellipse -> bspline  2.172e-12 cm
  CAD line    -> line     1.237e-11 cm         CAD line    -> line     2.542e-14 cm
  CAD line    -> line     5.217e-12 cm         CAD ellipse -> bspline  2.605e-14 cm
  CAD line    -> line     2.155e-12 cm         CAD line    -> line     8.360e-12 cm
  CAD line    -> line     2.542e-14 cm         CAD line    -> line     2.887e-11 cm  <--
```

The face maximum, 2.887e-11 cm, is on a **line** segment. That is the converter's long-standing
vertex-chaining for straight edges — each `line` runs from its own start vertex to the *next*
edge's start vertex — and it carries the source B-rep's own vertex-versus-curve slack. The two
ellipse arcs sit at 2.6e-14 and 2.2e-12 cm. Against `Bucket`'s declared BRep tolerance of 1e-08 cm,
the whole face is 2.9e-03 declared tolerances.

`Stream_O` §5 quotes different numbers for the same faces (4.6e-13 cm and ≤3.0e-11 cm) because it
measures a different thing: the CAD edge's distance from *two implicit surfaces*. Both are correct
and they agree on the order of magnitude for `Bucket`.

---

## 5. Acceptance

Environment and commands as in `NEXT.md`. Pre-change runs used the converter at `b7fcf9afd2`,
copied into a scratch directory; every run is a full reconvert (`--skip-convert` would have been
invalid, the converter changed).

### 5.1 Self-tests

| | |
| --- | --- |
| `O2_CADtoTGeo.py --self-test` | **36 checks, 0 failures** (18 + 8 + 10) |
| `runOracleGate.py --self-test` | **17/17** |
| `csg/emit.py --self-test` | **33/33** |

### 5.2 Fixtures gate — 10/10, and the disagreement counts with them

```
10/10 scored part(s) pass on the representation they ship in
10/10 scored part(s) pass on the surface representation
oracle disagreements outside tolerance (surface representation):
    contains=0  distout=0  distin=0  safety=0
  totals per representation:
    surface  contains=0 distout=0 distin=0 safety=0   (10/10 parts with zero disagreements)
    shape    contains=0 distout=0 distin=0 safety=0   (2/2 parts with zero disagreements)
 *oblique_cut_cyl  surface  O2BVHSurfaceSolid  0  0  0  0   capacity 4.55e-10  reliable
```

Gate exit code **0**, where the pre-change run exited 1. (The gate exits non-zero whenever any part
fails *or* is unscored, so this is a consequence of `oblique_cut_cyl` now being scored, not an
independent result. Read the counts, not the exit code.)

### 5.3 Bagger gate — 13/13, and Bagger now ships no mesh at all

```
tiers: CSG 7, exact surfaces 6, tessellated 0  (of 13 leaf solids)

13/13 scored part(s) pass on the representation they ship in
10/13 scored part(s) pass on the surface representation
oracle disagreements outside tolerance (surface representation):
    contains=0  distout=0  distin=0  safety=0
  totals per representation:
    surface  contains=0 distout=0 distin=0 safety=0   (13/13 parts with zero disagreements)
    shape    contains=0 distout=0 distin=0 safety=0   (7/7 parts with zero disagreements)
 *Bagger/Bucket_0_1_1_6  surface  O2BVHSurfaceSolid  0  0  0  0  capacity 3.97e-10  reliable
```

The `10/13 on the surface representation` is the three rams whose *surface* capacity is 1.3–1.4e-06
against a 1e-06 band; they ship as exact CSG (`dV_sym = 0`) and the historical column is unchanged
in definition. It was 9/12 before, i.e. the same three parts.

### 5.4 `compareGateRuns.py`, pre vs post — every moved field accounted for

| corpus | parts | differences | what they are |
| --- | ---: | ---: | --- |
| fixtures | 10 vs 9 | **39** | 1 `[EXTRA] oblique_cut_cyl` + **38 path strings** (`representations.*.source`, `verdict.shipped.source`, `verdict.shippedVerdict.source` — the two runs have different workdirs) |
| Bagger | 13 vs 12 | **56** | 1 `[EXTRA] Bagger/Bucket_0_1_1_6` + **55 path strings**, same three field names |

**Not one numeric or integer field moved on any pre-existing part.** Independently confirmed at the
strongest possible level: every `surfaces_*.bin` and `facets_*.bin` of the 9 unchanged fixtures and
the 12 unchanged Bagger parts is **byte-identical** between the two runs. The only new files are
`surfaces_oblique_cut_cyl_0_1_1_1.bin` and `surfaces_Bucket_0_1_1_6.bin`.

That also settles what the X-ray counters below mean: any counter on an unchanged part is
pre-existing by construction, not by argument.

### 5.5 X-ray, `--beams 96` (96 Fibonacci directions × 48² rays = 221184 rays per part)

| corpus | repr | rays identical to OpenCascade | LOST | extra | displaced | worst dt | parts clean |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| fixtures | surface | **2211835 / 2211840** | 0 | 2 | 3 | 8.8e-06 cm | 9/10 |
| fixtures | shape | 442368 / 442368 | 0 | 0 | 0 | 2.4e-13 cm | 2/2 |
| Bagger | surface | **2875391 / 2875392** | 0 | 1 | 0 | 5.4e-10 cm | 12/13 |
| Bagger | shape | 1548288 / 1548288 | 0 | 0 | 0 | 1.9e-10 cm | 7/7 |

Robustness counters, surface representation, both stepping modes, both corpora:
`zeroStep=0 noAdv=0 unstick=0 capHit=0 nonAlt=0 dupXing=0 parityNB=0 noTrans=0 outWorld=0 orgIn=0`,
and **mode (a) shape-API vs mode (b) `TGeoNavigator`: 0 of 2211840 and 0 of 2875392 rays disagree**.
The non-zero ones:

- **fixtures: `unterm=2 oddList=2 parity=2`, all on `cyl_inter_cyl`** (221179/221184, extra=2,
  displaced=3). Pre-existing — that part's sidecar is byte-identical across the change and it has
  no planar face at all. It is the ×1 ladder residual `Tutorial.md` §5.7 already records.
- **`oblique_cut_cyl`: 221184/221184 rays identical, every counter zero in both modes**, worst
  dt 4.46e-14 cm. Chord-integral volume 11.311442 cm³ against OCCT's chord integral over the same
  rays: **3.3e-14 relative**.
- **Bagger: `unterm=1 oddList=1 parity=1`, all on `Bucket`** — one extra crossing in 221184 rays
  (4.5e-06), worst dt 5.41e-10 cm, LOST=0. This is **new**, and it is the one number in this
  document that a reviewer should push on (§7). Chord-integral volume 58.454077 cm³ vs OCCT's
  58.454077 cm³ over the same rays.

### 5.6 Harness

`oblique_cut_cyl` and `Bucket` are both `reliable`, closed, consistently oriented; edge identity
3/3 and 247/247 source edges shared, 0 boundary, 0 non-manifold, 0 reversed, 0 degenerate; rims
4/4 and 109/109 matched. `--loop-crosscheck` at 50000 points / 20000 rays: **BVH == _Loop
bit-for-bit on every query of both parts**.

---

## 6. The one number that looks bad, and what it actually is

`GetMaxSharedEdgeDeviation()` on `oblique_cut_cyl` reads **7.635e-03 cm** — six orders larger than
any other fixture (next worst 4.7e-09). `Bucket` reads 4.567e-05 cm against a Bagger next-worst of
2.1e-07. Both parts are new, so both numbers are new, and on a project whose central claim is
edge-sharing this needs an answer rather than a footnote.

**It is a polyline-flattening artefact and it is on the cap circle, not on the ellipse.** Three
measurements, in the order they were made:

1. **The stored curves agree.** Both faces' realisations of every shared edge of `oblique_cut_cyl`
   were measured against the same CAD 3D curve. Plane side: ellipse 7.9e-15 cm. Cylinder side:
   ellipse 3.9e-09 cm (OCCT's pcurve for the sinusoid h(φ), which is genuinely not a polynomial in
   (φ, h) — this is the *cylinder's* pre-existing standing, not the plane's), seam lines 3.4e-16
   and 5.0e-13 cm, cap circle 5.0e-13 cm. So the two faces' stored boundaries cannot disagree by
   more than ~4e-09 cm. Whatever 7.6e-03 is, it is produced **downstream** of the stored geometry.

2. **The kernel's flattening was reimplemented and it reproduces exactly.** Running the kernel's
   `bsplineSampleRecursive` rule (three interior probes, `flatnessSq = 1e-10`, `maxDepth 16`, never
   across an interior knot) on both stored curves gives **842 chords on the cylinder side and 1274
   on the plane side** — precisely the counts the harness reports — and their symmetric Hausdorff
   is **3.61e-04 cm**, which is precisely the harness's reported *rim isolation* (3.59e-04 cm). So
   the ellipse rim is not the 7.6e-03 either.

3. **The cap circle is, and the mechanism is a quarter-chord phase offset.** The bottom cap edge is
   an `arc` on the plane and an iso `line` in (φ, h) on the cylinder; the harness flattens both into
   24 chords of a circle of radius 1.196 cm, whose sagitta is 1.028e-02 cm. Sampling the two
   realisations from the stored records and offsetting their phase:

   | phase offset | symmetric Hausdorff |
   | ---: | ---: |
   | 0 chord | 6.49e-14 cm |
   | **0.25 chord** | **7.697e-03 cm** |
   | 0.5 chord | 1.027e-02 cm |

   against the kernel's 7.635e-03 cm. Two equal-resolution polyline flattenings of the *same exact
   circle*, started at different azimuths because the cap plane's OCC frame and the cylinder's
   `XDirection` do not align on this fixture.

This is therefore a **resolution artefact of the reported measurement**, not a geometry
disagreement, and not something the ellipse work introduced — it is `kBSplineFlatness` /
rim-sampling, which `Tutorial.md` §5.5 already lists as open. It decides nothing: closure decides on
**identity** (§5.4 of the Tutorial), the number is reported, and the parts pass the oracle gate
0/0/0/0 and the X-ray 221184/221184.

**What it does show is that the number is not comparable across parts** — it scales with the
coarsest rim's chord sagitta, which for `oblique_cut_cyl` is 1.03e-02 cm and for `box` is zero
because a box's rims are straight. Quoting it without `GetRimChordResolution()` beside it is the
same error as quoting a gate total without its disagreement counts.

---

## 7. Being adversarial about this result

**The one thing to challenge: `Bucket`'s single extra crossing.** 1 ray in 221184 produces an odd
crossing list, hence one unterminated transport and one parity mismatch, at a worst `dt` of
5.41e-10 cm — a grazing hit, and both stepping modes agree with each other about it, so it is a
property of the solid and not of the stepping. It is not localised to a face here. It could be the
new ellipse trim, or it could be any of `Bucket`'s 97 patches finally being exercised for the first
time — nothing in the corpus separates those two hypotheses, because before this change *no* ray
had ever been stepped through `Bucket`'s exact representation. Honest reading: **221183 of 221184
rays are identical to OpenCascade and the one that is not is unexplained.** A per-face localiser for X-ray crossings is what
would settle it, and it does not exist.

**The second thing: the deviation is machine precision on one part and 1e-11 on the other, and the
1e-11 is not the ellipse's fault** (§4). The line-chaining residual it comes from is pre-existing
and applies to every planar face this converter has ever emitted; `Bucket` is simply the first part
where anyone measured it per segment. If that residual matters, it is its own item.

**The third: coverage claims.** This moves two corpora that the gate can score. It moves **nothing**
on ALICE3 — `Stream_O` §4.1 measured `plane`-face ellipse rejections there as **zero**. The 15
ALICE3 solids in `Stream_O`'s rollup need the *other* mechanism (co-surface trims on quadric faces),
which is untouched here.

**What was not done:** hyperbola and parabola boundaries (§2, deliberate); any C++ test (§3, no
kernel change); ALICE3 (nothing to gain); `oTOF` (still unreachable from the converter,
`NEXT.md` item 2).

---

## 8. Claims in other documents that this invalidates

- **[`Tutorial.md`](Tutorial.md) §6**: `Bucket` ships as a mesh because of "97 faces, 4 spherical +
  2 toroidal". Wrong — `Bucket` has **zero** faces without an analytic surface, and its two
  ellipse-bounded planes were the only blockers. The table row and the "tessellated 1" tier count
  are both now out of date; Bagger is CSG 7 / exact surfaces 6 / **tessellated 0**.
- **[`NEXT.md`](NEXT.md)** status table: "ladder fixtures 9/9 scored — of 10 leaf solids;
  `oblique_cut_cyl` has no sidecar and never has" and "`Bagger.step` 12/12 shipped … 1 unscored
  (`Bucket`, ships as mesh)". Both superseded. `--self-test` is 36 checks in three blocks, not 26 in
  two.
- **[`CodeReview_Fable.md`](CodeReview_Fable.md)** ladder table row for `oblique_cut_cyl`
  ("does not convert: `planar boundary edge is a ellipse curve`") and
  **[`CodeReview_Fable_v2.md`](CodeReview_Fable_v2.md)** §"G1 synthetic Boolean ladder"
  ("8/9 converted"), and **[`BVHSurfaceSolid.md`](BVHSurfaceSolid.md)** line 1177 and
  **[`TolerancePolicy.md`](TolerancePolicy.md)** line 313 ("`oblique_cut_cyl` still does not
  convert"). All four now read 10/10.
- **[`Stream_I_Verdict.md`](Stream_I_Verdict.md)** §"`Bucket` and `oblique_cut_cyl` are still
  unscored". Both are scored, and both pass.

I have not edited any of those; they are for central reconciliation.
