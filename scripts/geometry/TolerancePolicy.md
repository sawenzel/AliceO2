# Tolerance policy and the closure criterion — plan for Phase 1 item 4

Status: **steps 1-4 of five are done** (2026-07-31, commits `1bc5c4fbc9`, `9f45887ef7`,
`cacd64e4a5`, `f612f895a9`). **Step 5 — rim-based closure and `maxGap` — is not started**, and it
is the one with design content. See §5 for the state of each step and §8 for what step 5 now looks
like after the first four; §8 is the part to read before starting it.

This is the last outstanding Phase 1 item; Phase 2 should not start before it, because Phase 2's
whole premise is that "closed" becomes a statement you can act on.

This document covers findings **K3, K5, K9, K12, S8 and S10** of the review, which are all one
problem seen from six places: *the code compares distances that are not distances*, and then
judges closure by exact equality of numbers it never had reason to expect to be equal. K3, K5,
K12 and S10 are closed; **K9 and S8 are what step 5 is for.**

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

> **Historical, as of the state before steps 1-4.** §1.1 and §1.2 are fixed; their file:line
> references are pre-fix and will not resolve against the current tree. §1.3 is still true and is
> what step 5 addresses. Kept as written because the diagnosis is the argument for the design.

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

### 1.3 A closure criterion that cannot succeed (K9, S8) — still open, this is step 5

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
  closure check is a *diagnostic*. If a Bagger cylinder part reports a 1.3e-5 cm gap, the right
  outcome is that it stays unnavigable and now says how badly — not that epsilon becomes 1e-4. The
  2026-07-26 measurement puts the shared-edge pcurve disagreement at ~1.3e-5 model units, so
  expect the cyl-cyl parts to remain open. That is the correct answer until Phase 2.
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
5. **[not started] Rim-based closure + `maxGap` (K9/S8)**, then the §4.2 sweep, then reconcile
   `containsByParity`'s gating with what the sweep says. **Read §8 before starting.**

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
in cm?** On this corpus expect the cyl-cyl parts to stay open at ~1.3e-5 cm (the 2026-07-26
shared-edge pcurve measurement) and `tube_window` to stop claiming 1418 boundary edges about a
4-face solid.

**5b — let it decide.** Only then switch `NavigationReliability` over to the rim verdict, run the
§4.2 direction-independence sweep in the same commit, and reconcile `containsByParity`'s gating
with what the sweep says. §3.3's two prohibitions apply here and nowhere else: do not tune epsilon
until parts pass, and do not let the counts drift from what the enum documentation promises.

The Section 4.4 ambiguity-band refinement of `containsByParity` (deferred out of step 4) belongs
in 5b for the same reason: it changes which points get re-shot, so it needs the same sweep.

---
