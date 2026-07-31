# NEXT — session-start instruction for the BVHSurfaceSolid work

This file is the current hand-over. Whoever finishes a session should **rewrite it** for the next
one. Last updated 2026-07-31, after Phase 1 items 1-3 (commits `63d1d08119`, `29e8322f79`,
`f809b38dd2`, `346d4675d0`, `3eef31dfbe`, `e05797a8b9`).

Branch `swenzel/bvhsurfacesolid`. The tree is clean; everything through Phase 1 item 3 is
committed.

## Read first, in this order

1. The last handoff note at the end of [`BVHSurfaceSolid.md`](BVHSurfaceSolid.md).
2. [`TolerancePolicy.md`](TolerancePolicy.md), in full — it is the plan for the next item.
3. [`CodeReview_Fable.md`](CodeReview_Fable.md) Sections 5, 6 and 12.

## Do, in this order

**0. Diagnose `Bagger/BucketLink2_0_1_1_7` — and timebox it.** The only Bagger part that is
*navigable* and still fails the gate: 24 missed `distout` crossings, 48 `distin`, 6.3% capacity
drift. No current theory covers it. Separate the suspects (K6, K4) by measurement; write a
throwaway probe and do not plan a fix before you have evidence. If it has not yielded in an hour
or two, record what you found and move on.

**1. Phase 1 item 4** — the five commits in `TolerancePolicy.md` §5: `parametricMetric` on all six
surface classes; the wire join checks in kernel and IO; sidecar v2; the on-boundary band;
rim-based closure with a `maxGap`.

**2. K4 and K6**, if step 0 did not already resolve them.

Run the oracle gate **before and after every item** so each effect is attributable, and update the
docs and commit as you go.

## Baseline to beat

- `ctest -R BVHSurfaceSolid` green, **48 cases**.
- Oracle gate: fixtures **6/9** converting parts (`oblique_cut_cyl` still does not convert),
  Bagger **4/12**.
- `contains` disagreements outside tolerance: **0** on fixtures, **2** on Bagger.
- BVH == `_Loop` bit-identical everywhere.

Gate totals and `contains` counts are separate numbers. Do not quote one without the other.

## The one thing that will not be obvious

`Contains` now takes a **single parity shot when the solid is `Reliable`** and votes over five
directions otherwise. That fast path is licensed by a measurement *plus* the fact that the closure
check under-reports `Reliable`. Making closure succeed more often moves solids onto that path, so
**any change to the closure criterion must re-run the direction-independence sweep**
(`TolerancePolicy.md` §4.2). If a newly-`Reliable` part shows direction disagreements, fix the
closure criterion or drop the gating — do not leave it.

## Traps, each of which cost time before

- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or `ctest` and the harness
  silently use stale installed libraries.
- pythonOCC needs the alibuild **python3.10**, not the system 3.12.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`; use a scratch folder.
- For a quick kernel probe, compile a standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase`
  (recipe in `TolerancePolicy.md` §6). A ROOT macro trips over unrelated O2 headers on this
  machine.

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
