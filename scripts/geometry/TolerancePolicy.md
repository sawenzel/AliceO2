# Tolerance policy and the closure criterion — plan for Phase 1 item 4

Status: **done** (2026-07-31/08-01, commits `1bc5c4fbc9`, `9f45887ef7`, `cacd64e4a5`,
`f612f895a9`, `213b8aca8c` and 5b). §5 lists the five steps; **§9 and §10 are the measurements**,
and they are the part to read — 5a's numbers contradict this document's own prediction by four
orders of magnitude, and 5b's sweep leaves one offending point that is a live lead for K6.

This is the last outstanding Phase 1 item; Phase 2 should not start before it, because Phase 2's
whole premise is that "closed" becomes a statement you can act on.

This document covers findings **K3, K5, K9, K12, S8 and S10** of the review, which are all one
problem seen from six places: *the code compares distances that are not distances*, and then
judges closure by exact equality of numbers it never had reason to expect to be equal. K3, K5,
K12 and S10 are closed. **K9 and S8 are what step 5 is for**: 5a measures the gap and reports it
(done, §9), 5b lets it decide.

Read [`CodeReview_Fable.md`](CodeReview_Fable.md) Sections 5, 6 and 12 first.

---

## 0. `BucketLink2` — done, and it was not the kernel

This section used to ask for that diagnosis. It was carried out on 2026-07-31 with three throwaway
probes; the full account, with every measurement, is in
[`CodeReview_Fable.md`](CodeReview_Fable.md) **Section 13**. The short form:

- The **24 + 48 distance disagreements are an artifact of the gate.** OpenCascade classifies all
  20 recorded offender ray origins *opposite* to the category the sample generator assigned, the
  exact solid agrees with OpenCascade on 20 of 20, and where the category is right the exact
  solid's distances match the oracle to every printed digit. The categoriser is the tessellated
  reference, and this mesh is not watertight (its own `Check` reports hundreds of two-neighbour
  facets) and puts the part's left plate half a centimetre off.
- The **6.3% capacity drift is quadrature, not geometry.** A 4M-point Monte Carlo of the exact
  solid gives 17.061 ± 0.052 cm³ against OpenCascade's 17.079, while `Capacity()` returns 16.004.
  `integrateOverCurveTrim` is a fixed-128 midpoint rule over a characteristic function, so it
  converges at O(1/N); the whole capacity column of the Bagger gate is this and nothing else.

So `BucketLink2` is **not** evidence for K4 or K6, and it is no longer a lead. Two precisely-scoped
items replace it, both recorded in Section 13: gate ray-category soundness, and Green's-theorem
capacity for wire-trimmed quadrics. Neither blocks the work below.

---

## 1. What is actually wrong

> **Historical, as of the state before the five steps.** All three are fixed and their file:line
> references are pre-fix, so they will not resolve against the current tree. Kept as written
> because the diagnosis is the argument for the design — but see §9.2 for the one claim in §1.3
> that the measurement did *not* bear out: on both corpora the two faces of a shared edge do in
> fact sample it at matching phases, so vertex matching happens to work here. The structural
> defect is real (a unit test exhibits it); this corpus does not exercise it.

### 1.1 Distances in mixed units (K3, K12, S10) — fixed in step 2

Quadric trims live in parametric domains whose two coordinates have different units and different
physical meaning:

| surface | (u, v) | u unit | v unit |
| --- | --- | --- | --- |
| cylinder | (phi, h) | rad | cm |
| cone | (phi, h) | rad | cm |
| sphere | (phi, theta) | rad | rad |
| torus | (phiRing, phiTube) | rad | rad |
| plane / curved plane | (u, v) along axisU, axisV | cm | cm (axes need not be orthonormal) |

Every tolerance below is applied to `sqrt(du^2 + dv^2)` in those coordinates, i.e. to a quantity
with no physical meaning at all:

- `kWireJoinTolerance = 1e-5` (`BoundedSurface.h:63`), used at `:1505` — wire closure acceptance.
- `kToleranceSq` in `SurfaceWire`'s polygon join (`BoundedSurface.h:417`) — a *different* constant
  for the same job on the same extractor output, which is K12.
- `kToleranceSq` in `CurveWire::classify`'s on-boundary band (`BoundedSurface.h:1585`) — K5.
- `kJoinTolerance = 1e-5` in `O2SurfaceSolidIO.cxx:162`, checked at `:175` — S10, the IO twin.

The consequence is not academic and has a measured reproducer. On the ALICE3 assembly,
`ST1829909_01` (1052 surfaces) extracts but does **not** load: `LoadSurfaceSolid` rejects it with
*"surface 1006: wire edge 1 end does not join the next edge start"*. Over all 15 sidecars (17212
edge joins) exactly this one volume exceeds the tolerance, at **6 joins, worst 2.96e-5**, all of
them on cylinder wire trims. In arc length on those cylinders that drift is far below any
tolerance worth having; in radians it is 3x over. Meanwhile a large-radius cylinder passing the
same test at 1e-5 rad is accepting a real 3D gap of `r * 1e-5`, which on a 92 cm oTOF cylinder is
about a *micron per radian-unit of radius* — small, but the point is that nobody chose it.

So: small holes are falsely open, big cylinders falsely closed, and both follow from the same
missing radius factor.

Note the model is **not** on this machine — `scripts/geometry/STEP_examples/` holds only
`as1-oc-214.stp`, `Bagger.step` and `oTOF System V3-R92cm.step`. Reproduce the *class*
synthetically instead (see the first test in §4.3); the measured numbers above are the
acceptance target.

### 1.2 A boundary band tighter than the thing it measures (K5) — fixed in step 4

`CurveWire::classify` returns `WireClassification::Boundary` when a point is within `kTolerance`
(1e-9) of a trim curve. For a B-spline trim the curve it measures against is the *flattened
polyline*, which `bsplineSampleInto` produces to a flatness of ~1e-5. A 1e-9 band around a 1e-5
approximation is not a band, it is noise: the `Boundary` state is effectively unreachable for
B-spline trims, and in/out flips arbitrarily within +-1e-5 of the true curve.

Worse, winding and distance use two slightly different polylines (canonical endpoints vs. the raw
cache), so the two can disagree about the same point.

This is the same defect as K3 wearing a different hat, and it needs the same metric: the band has
to be expressed as a *length* and set from the representation's actual accuracy.

### 1.3 A closure criterion that cannot succeed (K9, S8) — fixed in step 5

`validateClosure` (`BoundedSurface.h:3919`) quantizes 3D vertices onto a 1e-7 cm lattice
(`kClosureQuantum`) and matches *chord endpoints* by exact equality of the lattice key. Two
independent failure modes, both structural:

1. **Exact bin matching, no neighbour bins.** Two vertices 1.0e-7 apart can land in different
   cells and never match, however close they are.
2. **The chords are not the same chords.** Each face samples its own rim independently —
   `sampleCurveWireByU` with `kArcSamples = 24` per turn, and B-spline edges resampled per
   pcurve. Chord *counts* are synchronized but sample *phases* only coincide for quarter-turn
   related frames. So the two faces of a shared edge emit different point sets for the same curve,
   and no tolerance on vertex equality can fix that, because the vertices genuinely are not the
   same points.

That is why "no real CAD part closes" — it is a property of the check, not only of the geometry.
The counts are also chord-inflated: one unmatched rim reports ~`kArcSamples` boundary edges, which
is why `tube_window` reports **1418** of them for a 4-face solid.

A closure check that cannot say "closed" is a check that tells you nothing when it says "open".

---

## 2. The fix, part A — one metric, applied everywhere

### 2.1 The metric

Give `BoundedSurface` a first fundamental form at a parametric point:

```cpp
/// The first fundamental form (guu, guv, gvv) of the surface at parametric point \a uv, so that
/// a parametric displacement (du, dv) spans the 3D length
///   sqrt(guu*du*du + 2*guv*du*dv + gvv*dv*dv).
/// This is what turns the parametric coordinates -- which mix radians and centimetres -- back
/// into distances, and it is the only honest way to compare a trim tolerance against anything.
virtual void parametricMetric(const Vec2& uv, double& gUU, double& gUV, double& gVV) const = 0;
```

with the closed forms below. All four quadric families are orthogonal (`gUV = 0`); the planar
families are the only ones that need the cross term, and they already compute exactly this metric
at `BoundedSurface.h:1947-1950` — it is simply not used for any tolerance today:

| surface | gUU | gUV | gVV |
| --- | --- | --- | --- |
| plane / curved plane | `dot(axisU, axisU)` | `dot(axisU, axisV)` | `dot(axisV, axisV)` |
| cylinder | `r^2` | 0 | 1 |
| cone | `r(h)^2` | 0 | `1 + k^2`, `k = dr/dh` |
| sphere | `(R sin theta)^2` | 0 | `R^2` |
| torus | `(Rmajor + r cos phiTube)^2` | 0 | `r^2` |

Add the free helper

```cpp
inline double parametricLengthSq(double gUU, double gUV, double gVV, const Vec2& delta);
```

Two traps in the table. The sphere's `gUU` vanishes at the poles and the cone's at an apex — a
*zero-length* separation in phi there is geometrically zero, which is correct, but any code that
divides by the scale must handle it. And the metric varies over the domain, so a join test must
evaluate it at the join point, not once per surface.

### 2.2 Where to apply it

Each of these becomes "convert the parametric separation to a length with the surface's metric,
then compare against one declared length tolerance":

1. `CurveWire::initialize` join check (`BoundedSurface.h:1505`) — currently `kWireJoinToleranceSq`.
2. `SurfaceWire::initialize` join check (`BoundedSurface.h:417`) — currently `kToleranceSq`. Same
   extractor feeds both, so after this they use the same rule and **K12 is closed by
   construction**.
3. `CurveWire::classify` on-boundary band (`BoundedSurface.h:1585`) — K5, see §2.4.
4. `O2SurfaceSolidIO.cxx:162/175` — S10. The loader knows the surface type and its radii from the
   record parameters, so it can build the same metric before validating the wire. Do not duplicate
   the formulas: expose them from the kernel and call them from the reader.

**Design decision to make explicitly:** `CurveWire` does not know its surface. Either
(a) pass a metric callable into `initialize`, defaulting to the identity metric for wires with no
surface context, or (b) move join validation up to the surface. (a) is smaller and keeps the wire
testable on its own; (b) is tidier. Pick one and say why in the commit — do not do both.

### 2.3 One tolerance family, and where the number comes from

Replace the four scattered constants with a documented family, all in **cm**:

- `kLengthTolerance` — the existing 1e-9 generic tolerance, unchanged in meaning.
- `kJoinTolerance` — wire closure acceptance. Today's 1e-5 was chosen for a mixed-unit quantity;
  once it is a real length, state what it means and pick it from the *extractor's* precision
  (~1e-6 cm per the 2026-07-24 note), not from the old value.
- `kTrimBandTolerance` — the on-boundary band of §2.4.

Then carry the **model's own tolerance** in the sidecar so the kernel stops guessing:

- Bump the sidecar to **version 2**: keep the `version, nSurfaces, reserved` uint32 header and
  append a `float64 modelTolerance` (STEP's declared tolerance; the oracle already extracts the
  max sub-shape tolerance — see `occtOracle.py`, which reports it as `tolerance` and uses it as
  the gate band).
- The reader accepts v1 (with a documented fallback value and a warning) and v2.
- `O2_CADtoTGeo.py` writes it.
- Store it on the solid so `CloseShape` can use it as the declared epsilon of §3.

This is the item the review flags as the reason the C++ side "literally cannot know what epsilon
to glue with".

### 2.4 The on-boundary band (K5)

Set the band from the representation's real accuracy rather than from optimism:

- line and arc trims are exact: band = `kLengthTolerance`.
- B-spline trims are polylines flattened to `sqrt(flatnessSq)`: band = that flatness, converted
  through the metric, plus the model tolerance from §2.3.

Have `Curve2D` report its own `representationTolerance()` so the wire can take the max over its
curves rather than assuming. And make winding and distance use **one** polyline — the discrepancy
between the canonical-endpoint and raw-cache versions is a bug of its own hiding inside K5.

This is also the "ambiguity band" of `CodeReview_Fable.md` Section 4.4. Note that Phase 1 item 1
deliberately did **not** build it: the measurement showed a reliability-gated re-shoot was both
sound and free on clean geometry, and a band in mixed units would have been neither. With a real
metric the band becomes worth having as the *finer* trigger, so that a `Reliable` solid can still
notice a locally ambiguous hit. Add it as a refinement of `containsByParity`, and re-measure
(§4) — do not assume it improves anything.

---

## 3. The fix, part B — a closure criterion that can succeed

Replace vertex-equality matching with **rim matching**, and report a number.

### 3.1 Shape of it

- Each surface emits its boundary as **rims** — one ordered 3D polyline per wire — instead of a
  flat list of chords. `appendDirectedEdges` becomes `appendRims`, or gains a sibling; the seven
  existing overrides (`BoundedSurface.h:2139, 2381, 2697, 3010, 3382, 3784, 3894`) already build
  exactly this data before flattening it.
- Build a spatial hash over all rim sample points with cell size ~ the declared epsilon.
- For each rim, and each of its sample points, find the nearest point lying on a rim **of a
  different face**. The rim's `gap` is the maximum of those distances along it.
- Classify per rim: matched (gap <= epsilon, exactly one partner face, opposite traversal),
  reversed (matched but same traversal), non-manifold (>= 2 distinct partner faces), boundary
  (unmatched).

Because this compares *curves* rather than vertices, two faces that sample the same edge with
different counts and phases match — which is precisely the structural blocker in §1.3.

### 3.2 What it must report

`ClosureReport` keeps its counts (callers and `NavigationReliability` depend on them) and gains:

- `maxGap` — the largest rim gap in cm. **This is the number the whole item exists to produce.**
- `unmatchedRimLength` and its fraction of total rim length — an honest replacement for
  chord-inflated counts. `tube_window` should stop saying "1418 boundary edges" about a 4-face
  solid and start saying how much boundary is open and by how much.
- counts expressed per *rim*, not per chord.

Report `maxGap` in `Print()`, in the harness's per-part line and in its `--json` under
`navigation`, next to the existing reliability fields.

### 3.3 Two things not to do

- **Do not tune epsilon until parts pass.** The oracle gate is the arbiter of correctness; the
  closure check is a *diagnostic*. If a Bagger cylinder part reports a gap, the right outcome is
  that it stays unnavigable and now says how badly — not that epsilon grows until it passes.
  > **Superseded in its numbers, not in its rule.** This paragraph used to predict, from the
  > 2026-07-26 shared-edge pcurve measurement, that the cyl-cyl parts would stay open at
  > ~1.3e-5 cm. Step 5a measured them: they are open by **0.25 to 0.75 cm**. See §9.1. The
  > prohibition stands and matters more, not less — but nobody should carry the 1.3e-5 forward.
- **Do not let the counts drift from the semantics `NavigationReliability` promises.** Its doc
  comment states what each state means; if rim matching changes what "boundary edge" counts, update
  the enum documentation in the same commit.

### 3.4 The coupling nobody will expect — read this

Phase 1 item 1 made `Contains` take a **single parity shot when the solid is `Reliable`** and vote
over five directions otherwise (`O2BVHSurfaceSolid.cxx`, `containsByParity`). That fast path is
justified by a measurement: on the Phase 0 corpus every `Reliable` part had zero disagreements
between shooting directions, and every part that disagreed was already rejected by the closure
check. The argument for safety is that the closure check **under-reports** `Reliable`.

A closure criterion that succeeds more often moves solids onto the fast path. So:

> **Any change to the closure criterion invalidates the measurement that licenses the single-shot
> fast path, and must re-run it.**

That is §4's direction-independence sweep. If a newly-`Reliable` part turns out to have direction
disagreements, the correct response is to fix the closure criterion (it is over-reporting) or to
drop the gating and always vote — not to leave it.

---

## 4. Acceptance

### 4.1 The gate, before and after every sub-item

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
ninja o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH ctest -R BVHSurfaceSolid

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
```

**Starting line (2026-07-31, after Phase 1 items 1-3):** `ctest` 48 cases green; gate **G1 6/10**
(6 of 9 converting fixtures pass; `oblique_cut_cyl` still does not convert), **G3 4/12**;
`contains` disagreements **0** on fixtures and **2** on Bagger; BVH == Loop bit-identical
everywhere. Any of those going backwards is a regression, whatever the closure numbers do.

### 4.2 The direction-independence sweep (mandatory after §3)

For every part in both models, over the full sample set, `ContainsAlongDirection` must agree with
`Contains` for a spread of directions on every part now reported `Reliable`. See §6.

### 4.3 New tests this item owes

- **Join metric on a large and a small radius.** Two cylinder trims with the *same* parametric
  join drift (say 2e-5 rad) on radius 0.01 cm and radius 100 cm: the small one must be accepted
  (sub-nanometre in arc length), the large one rejected (2e-3 cm is a real gap). Today's rule
  cannot tell them apart, which is the whole point. This is the synthetic form of the
  `ST1829909_01` reproducer of §1.1.
- **Sphere pole and cone apex**, where `gUU` vanishes.
- **Non-orthonormal planar axes**, where `gUV != 0`.
- **Sidecar v1 and v2 both load**, and v2's tolerance reaches `CloseShape`.
- **Closure reports a gap size.** Build a box with one face shifted by a known delta; `maxGap`
  must equal delta. Build the same box with one face resampled at a different chord count; it must
  still report **closed** — that is the structural fix of §1.3, and it fails today.
- **`tube_window` reports a rim-based count**, not ~1418 chords.

---

## 5. The five steps, and where they stand

Each step is a commit, gated before and after.

1. **[done, `1bc5c4fbc9`] `parametricMetric` on all six surface classes, plus tests.** Pure
   addition. Gate bit-identical, as required. The test checks each family against central
   differences of the surface's own parametrisation rather than restating the closed form.
2. **[done, `9f45887ef7`] Apply it to the wire join checks** (kernel and IO). `kWireJoinTolerance`
   is now a length in cm and moved 1e-5 → 1e-6; K3, K12 and S10 are closed. The value was measured,
   not assumed: 690 curve-wire joins across the Bagger and fixture sidecars have a worst residual
   of **4.06e-11 cm**, 471 of them exactly zero. Gate bit-identical, which is the expected outcome
   — no face on this corpus was ever near either threshold. Both tests were verified to fail
   against the pre-fix rule.
   - Design decision recorded: `CurveWire`/`SurfaceWire` take an optional `ParametricMetric`
     callback (option (a)), not moved-up validation (option (b)), because it keeps both wire types
     usable on their own and (b) needs a shared helper that amounts to the same callback.
3. **[done, `cacd64e4a5`] Sidecar v2 + model tolerance**, reader and writer, both versions loading.
   `Bagger.step` declares **1e-8 cm**. Nothing consults it yet — step 5 is where `CloseShape` does.
4. **[done, `f612f895a9`] The on-boundary band (K5)**, including collapsing the two polylines into
   one. Gate bit-identical; that is *expected rather than reassuring*, because a 1e-5-wide shell
   around a trim boundary is not something ~5000 bbox-spread samples per part will land in. The
   change is pinned by tests instead.
5. **Rim-based closure + `maxGap` (K9/S8)**, split by §8 into 5a and 5b.
   - **[done] 5a — measure and report, decide nothing.** `SurfaceRim`, a default `appendRims`,
     `measureRimClosure`, and the numbers on `ClosureReport`, `Print()`, the harness line and
     `--json`. Gate bit-identical on both corpora, which was the required outcome and is a real
     check here. **The measurement is in §9 and it did not come out as this document predicted.**
   - **[done] 5b — let it decide.** `NavigationReliability`, `IsClosed()` and
     `IsOrientationConsistent()` read the rim counts; the §4.2 sweep ran in the same commit and is
     in §10.1. Gate bit-identical, because the two criteria agree part for part here. The
     §2.4/Section 4.4 ambiguity band is the one piece **not** delivered: the free form of it was
     built, measured and dropped (§10.2), and the form §2.4 asks for is deferred behind diagnosing
     the single offending point of §10.1.

---

## 6. How to write a probe in a minute

The fastest diagnostic loop found in this project. `ContainsAlongDirection` is public, so a
standalone `.cxx` against the built library beats a ROOT macro — cling trips over unrelated O2
headers on this machine and will waste ten minutes before it fails.

```bash
source <your env>            # ALIBUILD_WORK_DIR, alienv, B=.../BUILD/O2-latest/O2
g++ -std=c++20 -O2 -o /tmp/probe /tmp/probe.cxx \
    -I$B/stage/include -I$HOME/alisw/O2/Detectors/Base/include \
    $(root-config --cflags --libs) -L$B/stage/lib -lO2DetectorsBase
```

The probe body: `LoadSurfaceSolid(sidecar, solid)`, `solid.CloseShape(false)`, then whatever you
want to ask. Sidecars and `.brep` references are left in the gate's workdir under
`<workdir>/db/<model>/`, and the oracle's own answers and sample sets under `<workdir>/oracle/`;
the offending points of any failing column are in `<workdir>/gate.json` under
`oracle.<column>.worstOffenders`.

If the same question gets asked twice, promote it to a harness flag
(`Detectors/Base/test/runSolidHarness.cxx`) instead of keeping a scratch binary — that is where
`--loop-crosscheck` and `--dump-samples` came from.

---

## 7. Environment traps

- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or `ctest` and the harness
  resolve `libO2DetectorsBase` from the installed prefix and you get an undefined symbol — or
  worse, a silently stale measurement.
- pythonOCC needs the alibuild **python3.10**; the system `python3` is 3.12 and cannot import
  `OCC`. `runOracleGate.py` sets this itself; a hand-run `occtOracle.py` does not.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`. Every session so far
  has used a scratch folder; keep it that way.
- The branch's targets appear only after a CMake re-run (`cmake .` in `$B`).

---

## 8. Step 5, after steps 1-4 — split it in two

§3 still describes what to build. What the first four steps changed is the *risk*, and there is a
decomposition that removes most of it. Do not do step 5 in one commit.

**The risk.** §3.4 is the whole problem: a closure criterion that succeeds more often moves solids
onto `Contains`'s single-shot fast path, which is licensed by a measurement taken when the closure
check under-reported `Reliable`. So "build rim matching" and "let rim matching decide navigability"
are two changes, and only the second one invalidates that measurement. Landing them together means
that if the §4.2 sweep cannot be finished, the solid is left with a fast path whose justification
no longer applies — which is exactly the half-finished state to avoid.

**5a — measure, report, decide nothing.** Add rims and the gap number, and change no verdict.

- `ClosureReport` gains `maxGap`, `unmatchedRimLength`, `totalRimLength` and per-rim counts.
  `closed`, `boundaryEdges`, `nonManifoldEdges`, `reversedEdges` and everything
  `NavigationReliability` derives from them stay exactly as they are.
- Report `maxGap` in `Print()`, in the harness's per-part line and in its `--json` under
  `navigation`, next to the existing reliability fields.
- The gate **must be bit-identical**, which is a real check here: if it is not, rim matching has
  leaked into a verdict.
- Cheap way in: `appendDirectedEdges` already emits each loop's edges consecutively, so a default
  `appendRims` on `BoundedSurface` can recover the rims by chaining maximal runs where one edge's
  end equals the next one's start, instead of overriding it in all seven classes. Same data, no
  per-class churn. Override later only if a class turns out to need it.
- Use the sidecar's model tolerance (step 3) as the declared epsilon, falling back on a documented
  constant when it is zero — `GetModelTolerance()` returns zero for "not stated" on purpose.

That commit alone answers the question the item exists to answer: **how far apart are the faces,
in cm?** It is done, and the answer is in §9: the cyl-cyl parts are open by a quarter to three
quarters of a centimetre, not by the ~1.3e-5 cm this document expected, and `tube_window` reports
seven rims of which three are open instead of 1418 boundary edges.

**5b — let it decide.** Only then switch `NavigationReliability` over to the rim verdict, run the
§4.2 direction-independence sweep in the same commit, and reconcile `containsByParity`'s gating
with what the sweep says. §3.3's two prohibitions apply here and nowhere else: do not tune epsilon
until parts pass, and do not let the counts drift from what the enum documentation promises.

The Section 4.4 ambiguity-band refinement of `containsByParity` (deferred out of step 4) belongs
in 5b for the same reason: it changes which points get re-shot, so it needs the same sweep.

---

## 9. Step 5a executed — what the gap actually is

5a is committed and both gates are bit-identical (fixtures 6/9, Bagger 4/12, `contains` 0 and 2;
`gate.json` matches the baseline field for field once the timing keys and the new `navigation`
keys are removed). `ctest -R BVHSurfaceSolid` is green at **57 cases**, from 53.

### 9.1 The headline: §3.3's prediction was wrong by four orders of magnitude

§3.3 says "the 2026-07-26 measurement puts the shared-edge pcurve disagreement at ~1.3e-5 model
units, so expect the cyl-cyl parts to remain open" — and treats that as the *reason* they stay
open. They do stay open, but not for that reason and not by that much. Matching at the model's own
declared 1e-8 cm:

| Bagger part | rims | matched | boundary | `maxGap` cm | chord res. cm | open / total cm |
| --- | --- | --- | --- | --- | --- | --- |
| `BasePin` | 4 | 4 | 0 | 9.5e-16 | 8.6e-3 | 0 / 25.1 |
| `Base` | 62 | 62 | 0 | 1.2e-12 | 0.262 | 0 / 660 |
| `Boom` | 52 | 52 | 0 | 9.7e-12 | 0.704 | 0 / 953 |
| `Stick` | 45 | 45 | 0 | 1.3e-11 | 0.691 | 0 / 631 |
| `BucketLink2` | 40 | 40 | 0 | 4.1e-11 | 0.057 | 0 / 270 |
| `BucketLink1` | 38 | 32 | 6 | **1.6e-3** | 0.041 | 20 / 215 |
| `BucketCylinderOuter` | 17 | 14 | 3 | **0.252** | 8.6e-3 | 5.9 / 60.1 |
| `BoomCylinderOuter` | 15 | 12 | 3 | **0.475** | 0.013 | 13.9 / 92 |
| `StickCylinderOuter` | 15 | 12 | 3 | **0.475** | 0.013 | 13.9 / 92 |
| `BoomCylinderInner` | 11 | 10 | 1 | **0.747** | 0.010 | 3.8 / 59 |
| `StickCylinderInner` | 11 | 10 | 1 | **0.747** | 0.010 | 3.8 / 59 |
| `BucketCylinderInner` | 11 | 10 | 1 | **0.748** | 7.7e-3 | 2.2 / 41.7 |

The six failing cylinder parts are open by **a quarter to three quarters of a centimetre**, over
4-15% of their rim length. That is not two pcurves disagreeing in the fifth decimal; it is a face
that is missing or trimmed to the wrong curve, and it is visible at the scale of the part. Whatever
Phase 2 does for these parts, "reconcile the shared-edge pcurves" is not it — the geometry is not
nearly closed and then spoiled by tolerance.

`BucketLink1` is the one part in a different class: 1.6e-3 cm, which is *below* its own chord
resolution of 4.1e-2 cm. That is the honest reading — for that part this measurement does not
resolve the gap, and it is the only part where the sagitta ceiling of §9.3 actually bites.

The fixtures agree, and add the tidiest result of the item:

| fixture | rims | matched | boundary | `maxGap` cm | chord res. cm | open / total cm |
| --- | --- | --- | --- | --- | --- | --- |
| `box`, `box_union_box` | 6, 10 | all | 0 | 0 | 0 | 0 |
| `box_minus_cyl` | 10 | 10 | 0 | 2.1e-14 | 6.8e-3 | 0 / 116 |
| `cyl_cross_cyl` | 12 | 12 | 0 | 3.9e-9 | 8.6e-3 | 0 / 80.7 |
| `cyl_inter_cyl` | 6 | 6 | 0 | 3.9e-9 | 1.7e-4 | 0 / 38.6 |
| `cyl_plus_cone` | 6 | 6 | 0 | 2.8e-13 | 8.6e-3 | 0 / 31.3 |
| `sphere_minus_cyl` | 4 | 4 | 0 | 8.3e-13 | 5.1e-3 | 0 / 15 |
| `torus_union_cyl` | 10 | 10 | 0 | 3.3e-13 | 0.028 | 0 / 142 |
| `tube_window` | 7 | 4 | **3** | **2.61** | 0.013 | **9.94 / 53** |

`tube_window` was the motivating case for "the counts are chord-inflated": **1418 boundary edges**
for a solid with seven boundary loops, of which three are open, accounting for 9.94 cm of 53 cm.
§3.2 asked for exactly that replacement and it is now what the harness prints.

### 9.2 What this means for 5b — it is a no-op on this corpus

Across both models, **every part the half-edge check calls closed has `boundaryRims == 0`, and
every part it calls open has `boundaryRims > 0`.** The two criteria agree part for part here, so
switching `NavigationReliability` to the rim verdict changes no verdict on either corpus, no solid
moves onto `Contains`'s single-shot fast path, and the §4.2 sweep should reproduce its previous
result exactly. §3.4's warning is real and still governs — it just does not fire here.

The reason the criteria agree is worth stating, because §1.3 predicted they would not: on both
corpora the converter emits **matching chord counts and phases** for the two faces of a shared
edge, so vertex matching happens to work. The structural blocker §1.3 describes is real — a unit
test builds a box whose last face is sampled at twice the chord count and the half-edge check calls
that closed box open — but no part of either corpus exercises it. Take that as "this corpus does
not test the thing", not as "the thing was not a problem".

### 9.3 Four things the plan did not anticipate

- **§8's "cheap way in" does not work as written.** It says `appendDirectedEdges` "already emits
  each loop's edges consecutively", so runs of consecutive edges can be chained. The wire-trimmed
  and planar faces do; the parametric-rectangle quadrics **interleave** their two rims, one bottom
  chord then one top chord. Chaining by matching endpoints instead costs nothing more, works for
  both, and still needs no change in any of the seven classes — so the goal of §8's suggestion is
  met, by a different mechanism.
- **A full-turn patch emits its seam twice, once each way**, when `fullSweep()` does not recognise
  the sweep as full. That pair bounds nothing and cancels in the half-edge check, which is why
  nobody had noticed it; chained naively it becomes a **two-point rim straddling the patch**, and
  it alone produced a spurious boundary rim and a gap the size of the patch (2 cm on
  `box_minus_cyl`) on five of the nine fixtures. Reversed duplicates within a face are cancelled
  before chaining, for the same reason the half-edge check cancels them.
- **The sampling noise floor has to be estimated from the turn angle, not from the polyline's
  deviation.** The obvious estimator — how far a rim vertex sits from the chord joining its
  neighbours — reports the *corner offset* on a box: 2.4 cm on a 4 cm box, where the true sampling
  error is zero. What separates a corner from a sample of a smooth run is the turn angle, so
  vertices turning by more than ~30° are skipped and the rest contribute `(chord/2)·tan(turn/4)`.
  A unit test checks that against the closed form `r(1 - cos(π/kArcSamples))` on a bare cylinder.
- **Probe chord midpoints, not rim vertices.** A box corner vertex legitimately lies on the rims
  of two other faces at once and would read as non-manifold; a chord interior belongs to exactly
  one shared edge.

### 9.4 The ceiling this measurement has, stated once

A rim is a polyline sampled at `kArcSamples = 24` per turn, radius-independent. Two faces sampling
one shared *curved* edge at different phases differ by up to the chord sagitta, `r(1-cos(π/24))` ≈
`8.6e-3 · r` cm, whatever the true gap is. So `maxGap` on a curved rim cannot resolve anything
below that, and `rimChordResolution` is reported next to it so nobody reads it as more than it is.

On this corpus it does not bite — the phases coincide, and every closed part measures below
4.1e-11 cm — but it is the reason `maxGap` must never be quoted alone, and it is the thing to fix
(rim sampling driven by a target sagitta in cm, not a fixed count per turn) if a future model
lands between 1e-8 and 1e-2 cm and needs an answer.

---

## 10. Step 5b executed — the rim verdict decides, and the sweep that licenses it

`IsClosed()`, `IsOrientationConsistent()` and `GetNavigationReliability()` now read the rim counts.
The chord counters stay as a diagnostic — they are still the cheapest way to see *how* two faces
disagree along a shared edge — but nothing derives a verdict from them, and the enum documentation
says what each state now means, per §3.3's second prohibition. `CloseShape`'s error messages report
the gap in centimetres and the open fraction of the boundary rather than a chord count.

Both gates bit-identical again (fixtures 6/9, Bagger 4/12, `contains` 0 and 2), which §9.2
predicted: the two criteria agree part for part on both corpora. `ctest` green at 57 cases.

### 10.1 The §4.2 direction-independence sweep

Run over both corpora, before and after the switch, from a standalone probe (§6): 21 parts, 11000
bbox-spread points each, **13** golden-spiral directions per point, `Contains` against
`ContainsAlongDirection`.

- **13 parts are `Reliable`; 143000 points; one disagreement.** Identical before and after the
  switch — the same 13 parts, the same single offending point. That is the measurement §3.4
  demands, and it says the fast path's licence is untouched, because no solid moved onto it.
- The parts that disagree between directions are exactly the parts the closure check rejects,
  at **0.55%-7.0%** of their points. The correlation §3.4 relies on still holds.
- **The one offender is worth carrying forward.** `cyl_cross_cyl`, point
  `(0.65334264649649720, 0.88394684996026007, 0.97463122696308724)`, `Contains` = outside, **1 of
  13** directions disagrees, `Safety` = 0.0992 cm. It is not a gap shadow: a gap costs a point most
  of its directions and puts it behind the gap, and this point is a millimetre from the nearest
  surface with 12 of 13 directions — including the fixed one — agreeing. One unlucky ray near the
  curve where the two cylinders cross. The earlier Phase 1 measurement reported "zero" over a
  smaller, differently seeded sample; this is a finer sample of the same thing, not a
  contradiction. **It is a live lead for K6** (cancellation in the naive quadratic formula), and
  the cheapest known reproducer of whatever it is.

### 10.2 What the ambiguity band was, and why nothing shipped

§2.4 defers "the Section 4.4 ambiguity band" into 5b: a `Reliable` solid should still notice a
*locally* ambiguous hit, and re-shoot for that point alone. §2.4 also says, in its own words, to
re-measure and **not to assume it improves anything**. It does not, in the form that was tried.

The form tried is the one the shot can supply for free: re-shoot whenever a parity cluster held
more than one coincident hit — that is, whenever the answer came out of the coincident-hit
cancellation rule rather than out of the geometry. It is self-checking in exactly Section 4.4's
sense, needs no change to any surface class, and the clusters are walked anyway. Measured:

| | box | box_minus_cyl | box_union_box | cyl_cross_cyl | cyl_inter_cyl | cyl_plus_cone | sphere_minus_cyl | torus_union_cyl |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| off, ns/call | 219.8 | 343.8 | 235.5 | 3010.9 | 8564.8 | 202.4 | 179.1 | 524.1 |
| on, ns/call | 226.9 | 398.1 | 236.2 | 3031.6 | 8673.5 | 202.9 | 180.5 | 525.6 |

0.2-1.3% on seven parts and **+15.8%** on `box_minus_cyl`, and it moved **not one** of the 143000
sweep points — including the single offender of §10.1, which it was aimed at. That is the whole
case against it: it fires where the cancellation rule is already right (grazes at convex and
concave edges, a ray straight through a shared edge — all of which the entering/exiting
cancellation handles correctly) and stays silent where the residual defect is. It is not kept, and
`containsByParity` carries the measurement in a comment so nobody rebuilds it.

**The instrument §2.4 actually asks for is a different one and is still unbuilt**: a band around
the trim curve in the *parametric* domain, sized through the surface's first fundamental form, so
that a hit landing within the representation's own accuracy of a trim boundary is flagged by the
patch that produced it. That needs `RayHit` to carry the flag and each of the seven
`appendIntersections` to set it — the hottest code in the solid. §10.1 sets its priority: the
measured headroom on the fast path is **one point in 143000**, and that point is a K6 suspect
rather than a trim-boundary one. Diagnose the offender first; build the band only if the offender
turns out to be what it is for.

---

## 11. The offender diagnosed — it is the trim boundary, and the band is now built

§10.2 ended by saying: diagnose the single direction-dependent point before building §2.4's
parametric ambiguity band, because the band is the fix for a *trim-boundary* ambiguity and that
point might not be one. It is one. The band is built.

### 11.1 What the offending ray actually does

`cyl_cross_cyl` is two unit cylinders on the z and x axes, fused, so the truth is analytic. At the
offending point the truth is **outside**, and the golden-spiral direction #10 of 13 is the one that
says inside. Its crossing list, against the two crossings the truth has:

| | t | n·d | hit | ρ_z − 1 | ρ_x − 1 |
| --- | --- | --- | --- | --- | --- |
| truth | 0.33836336623124036 | — | — | — | — |
| solid | 0.33836336623124036 | −0.203 | (0.76639, 0.64238, 0.76641) | −1.1e-16 | **+1.70e-05** |
| solid | 0.33838167219971466 | −0.930 | (0.76639, 0.64237, 0.76640) | **−3.71e-06** | +0.0 |
| truth / solid | 2.4325853070001155 | +0.930 | (1.46604, −0.85273, −0.52234) | +0.696 | −3.3e-16 |

The first hit is the z-cylinder and is genuine — it sits 1.7e-5 cm *outside* the x-cylinder, so it
is on the union's boundary. The second is the x-cylinder, 1.83e-5 further along, and it sits
3.71e-6 cm **inside** the z-cylinder. A point of the x-cylinder that is inside the z-cylinder is
interior to the fused solid; that patch should have been trimmed away there. It was not. Three
crossings instead of two, odd parity, "inside".

**The two hits are not clustered and never could be.** They are 5.4e-5 apart in relative terms and
`sameIntersection` merges at 1e-7·|t|. That is why §10.2's coincident-hit trigger stayed silent on
the one point it was aimed at — the mechanism was never coincident hits.

### 11.2 The trim overhang, measured directly rather than waited for

The two cylinders meet on `P(θ) = (cos θ, sin θ, ±cos θ)`. Walking along the x-cylinder's surface
away from that curve by a signed arc length and asking the solid whether its patch is still there
(shoot a short ray along the normal; a hit at the expected distance means yes), bisected at each of
720 θ on both branches:

| | value |
| --- | --- |
| sampled positions | 1440 |
| positions that **overhang** the true seam | **1440 (all of them)** |
| positions that undercut it | **0** |
| worst overhang | **1.95e-05 cm** |
| overhang floor | **1.00e-05 cm** |
| mean | 1.33e-05 cm |

So the error is **systematic and one-sided**: every patch keeps a sliver past the true seam, none
stops short. That is not what an approximation error looks like — an approximation would err both
ways — and the floor says why.

### 11.3 The cause is the on-boundary band's tie-break, i.e. step 4's own choice

`CurveWire::classify` returns `Boundary` within `boundaryBand`, whose floor for a B-spline trim is
`representationTolerance()` = `kBSplineFlatness` = **1e-5** — exactly the measured floor. And
`curveTrimContains` resolves `Boundary` as **inside the trim**:

```
  if (outerClassification == WireClassification::Boundary) { *boundary = true; return true; }
```

That is a tie-break, not a fact, and it is one-sided by construction. The excess above the floor,
0 to 1e-5 and oscillating with θ, is the polyline flattening of the B-spline seam on top of it.

Step 4 was right to widen the band — a 1e-9 band around a 1e-5 polyline was measuring noise (K5).
What step 4 did not do is say what a hit *in* the band means to a ray, and "accept it silently" is
the one answer that turns an undecidable point into a wrong answer with no trace.

**This is K5's other half, not K6.** K6 (cancellation in the naive quadratic) remains untouched and
has lost its only reproducer: the roots here are accurate to the last printed digit, and both
cylinders are hit exactly where they should be. The defect is entirely in which hits are *kept*.

### 11.4 What shipped

`RayHit` carries `onTrimBoundary`. The five `appendIntersections` that can produce a wire-trimmed
hit (planar, curved-planar, cylinder, cone, sphere, torus) set it from the `bool*` out-param
`curveTrimContains` already had and every caller was passing `nullptr` to. `parityAlong` reports
whether any counted crossing carried it, and `containsByParity` — on a `Reliable` solid, where it
would otherwise return a single shot — re-shoots through `containsByVote` when one did.

This is §2.4's instrument, and it is much cheaper than §2.4 feared: no new geometry, no per-class
trim mathematics, one bool through five call sites, because the classification already knew.

Measured on `cyl_cross_cyl`, against the analytic truth:

| | uniform, 5e6 points | seam-aimed, 2e6 points |
| --- | --- | --- |
| re-shoot fires | 48 (9.6e-06 of shots) | — |
| flagged hits / all hits | 58 / 3145428 (1.8e-05) | — |
| fires **and changes** the answer | 24 | 453008 (22.7%) |
| — of which corrected | **24** | **453008** |
| — of which broken | **0** | **0** |
| `Contains` wrong before | 24 (4.8e-06) | 453008 (22.7%) |
| `Contains` wrong after | **0** | **0** |

"Seam-aimed" starts each ray on the true crossing curve and backs off along the fixed direction, so
the ray passes through the sliver by construction — the adversarial case, where the old single shot
is wrong 22.7% of the time and the new answer is never wrong.

**Cost.** The re-shoot fires on one shot in 104000 and costs at most five extra parity shots when
it does, so ~5e-5 of `Contains` amortised. The gate's own timing cannot resolve that: its
`contains` column moved by −29% to +46% between two runs of the *same* code path, which is the
noise floor of this machine, not an effect. Do not quote the gate's timings for anything at this
scale.

### 11.5 Acceptance

- `ctest -R BVHSurfaceSolid` green at **58 cases**, from 57 (`TrimBoundaryHitsAreFlaggedAsAmbiguous`).
- Both gates **bit-identical** to the baseline: fixtures 6/9, Bagger 4/12, `contains` 0 and 2, and
  `gate.json` matching field for field once the `timing*`/`*Seconds` keys are removed. The gate
  cannot see this change — its `contains` column already agreed 5000/5000 on every `Reliable` part,
  and the defect is at the 1e-5 rate.
- §4.2 direction-independence sweep re-run over both corpora: 13 `Reliable` parts, 143000 points,
  13 directions — **0** direction-split points, and `Contains` differs from the raw fixed shot at
  **0** of them.

  Read that carefully. The sweep's split metric is `ContainsAlongDirection` against itself, and
  that entry point documents that it bypasses the re-shoot policy, so **the metric is blind to this
  change by construction**. The previous session's 1-in-143000 and this session's 0-in-143000 are
  different samples of the same rate, not a before/after. The improvement is the §11.4 table; the
  sweep's job here is only to confirm that no solid moved and nothing regressed.

### 11.6 What this does not fix

The sliver is still there. `Contains` now declines to be decided by it; `DistFromOutside`,
`DistFromInside` and `ComputeNormal` still consume the flagged hit silently, and `nearestCrossing`
has no equivalent policy (S4 is the same gap seen from the distance side). The honest scope is:
**parity is now self-checking about its own tie-breaks; the distance queries are not.**

Removing the sliver rather than labelling it needs the seam to be known better than either face's
own chart knows it — which is Phase 2, adjacency-based exact trims, and is what `onTrimBoundary`
would become the acceptance test for: on a model with exact adjacency, this flag should never fire.

---
