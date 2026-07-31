# NEXT — session-start instruction for the BVHSurfaceSolid work

This file is the current hand-over. Whoever finishes a session should **rewrite it** for the next
one. Last updated 2026-07-31, after the `BucketLink2` diagnosis and Phase 1 item 4 steps 1-4
(commits `5716d0070d`, `1bc5c4fbc9`, `9f45887ef7`, `cacd64e4a5`, `f612f895a9`).

Branch `swenzel/bvhsurfacesolid`. The tree is clean; everything below is committed.

## Read first, in this order

1. The last handoff note at the end of [`BVHSurfaceSolid.md`](BVHSurfaceSolid.md).
2. [`TolerancePolicy.md`](TolerancePolicy.md) **§8** — it is the plan for the next item, written
   after steps 1-4 and specifically about how to do step 5 without breaking something else. Then
   §3, which it refines, and §5 for what each finished step actually did.
3. [`CodeReview_Fable.md`](CodeReview_Fable.md) **Section 13** (why `BucketLink2` is closed and
   what replaced it) and the "Still open" list at the end of Section 12.

## Do, in this order

**1. Phase 1 item 4, step 5a — measure the gap, decide nothing.** Rim-based closure produces
`maxGap`, `unmatchedRimLength` and per-rim counts; `closed` and the edge counts that
`NavigationReliability` derives from stay exactly as they are. Report `maxGap` in `Print()`, in the
harness per-part line and in `--json` under `navigation`. Recipe, including the cheap way to get
rims without touching all seven surface classes, in `TolerancePolicy.md` §8.

**The gate must come out bit-identical.** That is not a formality here: if it moves, rim matching
has leaked into a verdict, which is 5b's job and carries 5b's obligations.

**2. Phase 1 item 4, step 5b — let it decide.** Switch `NavigationReliability` to the rim verdict,
run the §4.2 direction-independence sweep **in the same commit**, and reconcile
`containsByParity`'s `Reliable` gating with what the sweep says. The Section 4.4 ambiguity-band
refinement of `containsByParity` belongs here too — it was deliberately left out of step 4 because
it changes which points get re-shot and needs the same sweep.

**3. Then pick from these three, in whichever order the evidence favours.** None blocks the others.
- **Green's-theorem capacity for wire-trimmed quadrics.** `integrateOverCurveTrim` is a fixed-128
  midpoint rule over a characteristic function, so it converges at O(1/N) and the gate's 1e-6
  relative capacity band is unreachable at any practical N. The whole capacity column of the Bagger
  gate is this and nothing else — measured, see Section 13. The integrand does not depend on `h`
  for a cylinder, so an antiderivative in phi exists in closed form and the area integral collapses
  to a contour integral around the trim wire, evaluated by the same Gauss-Legendre the planar
  `signedAreaContribution` already runs. Same reduction for cone, sphere and torus. This would also
  make `capacityIsExact()` mean something.
- **Gate ray-category soundness.** The sample generator assigns a ray to `distout` or `distin`
  using the tessellated reference, which on `BucketLink2` is not watertight and is half a
  centimetre out of place. The oracle should report its own containment for each ray origin, and
  the harness should re-label — or at minimum decline to score — a ray whose category the oracle
  contradicts. The rule is sound because it never consults the candidate. Expect the Bagger gate to
  improve when this lands; that is a *gate* change and must be recorded as one, not read as a
  kernel improvement.
- **K4 and K6**, both untouched. K4: `bsplineSampleRecursive`'s degenerate-chord recursion probes
  only the parametric midpoint. K6: cancellation in the naive quadratic formula, and absolute,
  scale-dependent tolerances in the cone-degeneracy and torus-quartic branches. Note that
  `BucketLink2` is **no longer evidence for either** — if you work on them, find a new reproducer
  first.

Run the oracle gate **before and after every item** so each effect is attributable, and update the
docs and commit as you go.

## Baseline to beat

- `ctest -R BVHSurfaceSolid` green, **53 cases**.
- Oracle gate: fixtures **6/9** converting parts (`oblique_cut_cyl` still does not convert),
  Bagger **4/12**.
- `contains` disagreements outside tolerance: **0** on fixtures, **2** on Bagger.
- BVH == `_Loop` bit-identical everywhere.

Gate totals and `contains` counts are separate numbers. Do not quote one without the other.

## Two things that will not be obvious

**The fast path is licensed by a measurement, not by an argument.** `Contains` takes a single
parity shot when the solid is `Reliable` and votes over five directions otherwise. That is
justified by a sweep over the Phase 0 corpus *plus* the fact that the closure check under-reports
`Reliable`. Making closure succeed more often moves solids onto that path, so **any change to the
closure criterion must re-run the direction-independence sweep** (`TolerancePolicy.md` §4.2). This
is the entire reason step 5 is split into 5a and 5b. If a newly-`Reliable` part shows direction
disagreements, fix the closure criterion or drop the gating — do not leave it.

**A bit-identical gate is sometimes the expected answer and sometimes a red flag.** Steps 1 and 3
had to be bit-identical. Steps 2 and 4 came out bit-identical because nothing on this corpus was
ever near the thresholds they changed — that is worth stating, but it is not evidence the change
was inert, and those steps are pinned by tests that were checked to fail against the old code.
Step 5a must be bit-identical or something is wrong. Know which case you are in before you report
a result.

## Traps, each of which cost time before

- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or `ctest` and the harness
  silently use stale installed libraries.
- pythonOCC needs the alibuild **python3.10**, not the system 3.12. `runOracleGate.py` sets this
  itself; a hand-run `occtOracle.py` does not — export `PYTHONPATH` to
  `$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/...` and
  `LD_LIBRARY_PATH` to `$SW/OCCT/latest/lib:$SW/Python/latest/lib`.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`; use a scratch folder.
- For a quick kernel probe, compile a standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase`
  (recipe in `TolerancePolicy.md` §6, and add `-lGeom` — `root-config --libs` alone does not pull
  it in and the link fails on `TGeoBBox`). A ROOT macro trips over unrelated O2 headers here.
- `--skip-convert` reuses `<workdir>/db`, which saves minutes — but only when the *converter* did
  not change. Step 3 touched the sidecar writer and needed a full re-convert.
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
