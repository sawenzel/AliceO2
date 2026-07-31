# NEXT — session-start instruction for the BVHSurfaceSolid work

This file is the current hand-over. Whoever finishes a session should **rewrite it** for the next
one. Last updated 2026-08-01, after Phase 1 item 4 step 5 (commits `213b8aca8c`, `d32196029e`).

Branch `swenzel/bvhsurfacesolid`. The tree is clean; everything below is committed.

**Phase 1 is finished.** Every item of it is done and K9/S8 are closed. What follows is a choice
between four independent items, not a queue.

## Read first, in this order

1. The last handoff note at the end of [`BVHSurfaceSolid.md`](BVHSurfaceSolid.md).
2. [`TolerancePolicy.md`](TolerancePolicy.md) **§9 and §10** — the measurements step 5 produced.
   Two of them contradict what the rest of that document predicted, so read them before trusting
   any number stated elsewhere in it.
3. [`CodeReview_Fable.md`](CodeReview_Fable.md) **Section 14** (Phase 1 item 4 in full) and
   **Section 13** (why `BucketLink2` is closed and what replaced it).

## Do — pick one, in whichever order the evidence favours

None of these blocks the others. The first is new and is the smallest.

**1. Diagnose the one direction-dependent point.** The §4.2 sweep leaves exactly one disagreement
in 143000 points on `Reliable` parts: `cyl_cross_cyl` at
`(0.65334264649649720, 0.88394684996026007, 0.97463122696308724)`, `Contains` = outside, **1 of 13**
directions disagreeing, `Safety` = 0.0992 cm. Twelve directions including the fixed one agree, and
the point is a millimetre from any surface, so it is not a gap shadow. It is one unlucky ray near
the curve where the two cylinders cross — the cheapest known reproducer for **K6** (cancellation in
the naive quadratic formula), and probably a half-hour with the §6 probe recipe. Do this before
building §2.4's parametric ambiguity band: that band is the fix for a *trim-boundary* ambiguity,
and this point may not be one. The measured headroom on the fast path is this single point, so the
band's whole justification rests on what this turns out to be.

**2. Green's-theorem capacity for wire-trimmed quadrics.** `integrateOverCurveTrim` is a fixed-128
midpoint rule over a characteristic function, so it converges at O(1/N) and the gate's 1e-6
relative capacity band is unreachable at any practical N. The whole capacity column of the Bagger
gate is this and nothing else — measured, see Section 13. The integrand does not depend on `h` for
a cylinder, so an antiderivative in phi exists in closed form and the area integral collapses to a
contour integral around the trim wire, evaluated by the same Gauss-Legendre the planar
`signedAreaContribution` already runs. Same reduction for cone, sphere and torus. This would also
make `capacityIsExact()` mean something.

**3. Gate ray-category soundness.** The sample generator assigns a ray to `distout` or `distin`
using the tessellated reference, which on `BucketLink2` is not watertight and is half a centimetre
out of place. The oracle should report its own containment for each ray origin, and the harness
should re-label — or at minimum decline to score — a ray whose category the oracle contradicts. The
rule is sound because it never consults the candidate. Expect the Bagger gate to improve when this
lands; that is a *gate* change and must be recorded as one, not read as a kernel improvement.

**4. Phase 2 — adjacency-based exact trims.** Its premise was that "closed" becomes a statement you
can act on, and it now is. But **read `TolerancePolicy.md` §9.1 before scoping it**: the six failing
Bagger cylinder parts are open by 0.25-0.75 cm over 4-15% of their rim length, not by the ~1.3e-5 cm
of shared-edge pcurve disagreement the plan assumed. Reconciling two pcurves will not close them —
a face is missing or trimmed to the wrong curve, and finding out which is the real first step.

**K4 and K6** remain untouched. Item 1 above is the live lead for K6; K4
(`bsplineSampleRecursive`'s degenerate-chord recursion probes only the parametric midpoint) still
has no reproducer.

Run the oracle gate **before and after every item** so each effect is attributable, and update the
docs and commit as you go.

## Baseline to beat

- `ctest -R BVHSurfaceSolid` green, **57 cases**.
- Oracle gate: fixtures **6/9** converting parts (`oblique_cut_cyl` still does not convert),
  Bagger **4/12**.
- `contains` disagreements outside tolerance: **0** on fixtures, **2** on Bagger.
- BVH == `_Loop` bit-identical everywhere.
- Direction-independence sweep: **1** disagreement in 143000 points over the 13 `Reliable` parts.

Gate totals and `contains` counts are separate numbers. Do not quote one without the other.

## Three things that will not be obvious

**`maxGap` has a resolution, and it is printed next to it.** A rim is a polyline sampled at
`kArcSamples = 24` per turn, radius-independent, so two faces sampling one shared *curved* edge at
different phases differ by up to the chord sagitta `r(1-cos(π/24)) ≈ 8.6e-3·r` cm whatever the true
gap is. `rimChordResolution` reports that floor. On this corpus it does not bite — the phases
coincide and every closed part measures below 4.1e-11 cm — but never quote `maxGap` alone. The fix,
if a model ever lands between 1e-8 and 1e-2 cm, is rim sampling driven by a target sagitta in cm
rather than a fixed count per turn.

**The fast path is licensed by a measurement, not by an argument.** `Contains` takes a single
parity shot when the solid is `Reliable` and votes over five directions otherwise. **Any change to
the closure criterion must re-run the §4.2 direction-independence sweep**, because a criterion that
succeeds more often moves solids onto that path. The probe is thirty lines against the built library
(§6 recipe); a copy of the one used this session is worth rewriting rather than hunting for.

**Bit-identical is sometimes the answer and sometimes a red flag.** Steps 1, 3 and 5a *had* to be
bit-identical. Steps 2, 4 and 5b came out bit-identical because nothing on this corpus was near the
thresholds they changed — worth stating, but not evidence the change was inert, and each is pinned
by tests instead. 5b's case is the sharpest: the rim and chord criteria agree part for part here, so
the gate cannot see the switch at all, and only a unit test can.

## Traps, each of which cost time before

- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or `ctest` and the harness
  silently use stale installed libraries.
- The test binary is at `$B/stage/tests/`, not `$B/stage/bin/`; the harness is at `$B/stage/bin/`.
- pythonOCC needs the alibuild **python3.10**, not the system 3.12. `runOracleGate.py` sets this
  itself; a hand-run `occtOracle.py` does not — export `PYTHONPATH` to
  `$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/...` and
  `LD_LIBRARY_PATH` to `$SW/OCCT/latest/lib:$SW/Python/latest/lib`.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`; use a scratch folder.
- For a quick kernel probe, compile a standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase`
  (recipe in `TolerancePolicy.md` §6, and add `-lGeom` — `root-config --libs` alone does not pull
  it in and the link fails on `TGeoBBox`). A ROOT macro trips over unrelated O2 headers here.
- `--skip-convert` reuses `<workdir>/db`, which saves minutes — but only when the *converter* did
  not change. It is also the fast way to re-measure a kernel change: copy a finished workdir and
  re-run with `--skip-convert`.
- The gate exits non-zero whenever any part fails, which is the normal state. Read the counts, not
  the exit code.

## Commands

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

To compare two gate runs, diff `<workdir>/gate.json` with the timing fields (`timing*`,
`*Seconds`) stripped — everything else, including the offender lists and the query checksums, is
deterministic and must match when a change is meant to be inert.

To see the rim numbers for one part without a full gate run:

```bash
$B/stage/bin/o2-bench-detectorsbase-solid-harness --db <workdir>/db --parts <part> \
    --points 1 --rays 1 | grep -E 'navigation|rim gap'
```
