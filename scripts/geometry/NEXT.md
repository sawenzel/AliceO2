# NEXT — session-start instruction for the BVHSurfaceSolid work

This file is the current hand-over. Whoever finishes a session should **rewrite it** for the next
one. Last updated 2026-08-01, after the four follow-up items of the previous hand-over (commits
`7046188aef`, `d60584b276`, `01f64b4708`).

Branch `swenzel/bvhsurfacesolid`. The tree is clean; everything below is committed.

**All four items of the previous hand-over are done.** Three landed as code; the fourth was a
diagnosis, and it retired its own premise. What follows is not a menu this time — the gate has
collapsed onto a single question and there is one thread to pull.

## The one result that reframes everything

**Every part that passes the oracle gate is navigable, and every part that fails is not.** On both
corpora, after this session: fixtures 8/9 with only `tube_window` failing, Bagger 5/12 with exactly
the seven non-navigable parts failing. The capacity column stopped being a quadrature artifact and
the distance columns stopped being a categorisation artifact, so what is left in every failing
column is downstream of the surface set being open.

So: **the closure check is the whole remaining subject.** And the closure check is now known to be
reporting something other than what it is read as.

## Read first, in this order

1. [`TolerancePolicy.md`](TolerancePolicy.md) **§12** — item 4's measurements. It corrects §9.1,
   which the previous hand-over told people to trust, so read §12 before §9.
2. [`CodeReview_Fable.md`](CodeReview_Fable.md) **§19** (the same conclusion, with the next steps)
   and **§16** (the trim-boundary band, which is the acceptance test Phase 2 will want).
3. The last handoff note at the end of [`BVHSurfaceSolid.md`](BVHSurfaceSolid.md).

## Do — in this order, and the order matters

**1. Expose the rims. (Small. Everything else waits on it.)** The kernel reports rim *counts* and
aggregate lengths and nothing per rim, so the offending rim cannot be named — every conclusion in
§12 had to be inferred from counts, lengths and a hand-rolled B-spline evaluator over
`GetSurfaceRecords()`. Add a per-rim accessor or a `--json` block: face index, wire index, length,
matched/boundary/non-manifold state, and the position of the worst chord. Then `BoomCylinderInner`'s
one unmatched rim of 3.837 cm can be pointed at instead of deduced.

**2. Read the pairing, with the rim named.** §12.2 measured that on `BoomCylinderInner` the two
faces of the junction carry trims tracing the *same* curve — Hausdorff 2.4e-7 cm, each on the
other's carrier to 4e-8 cm, both 3.837 cm long, which is exactly the unmatched length. And §12.3
measured that loosening the matching epsilon to 1e-4 cm changes nothing. A rim whose partner
provably exists, in the same place, cannot be unmatched for a tolerance reason. Two hypotheses, and
the dump from step 1 separates them in minutes:

- **One of the two rims is never emitted.** The suspicious shape is a wire that is a *single closed
  curve*: `BoomCylinderInner`'s hole and `BucketCylinderOuter`'s face 9 are both exactly that, and
  dropping `BoomCylinderInner`'s hole rim would put its rim count at 11, which is what it reports.
  This class has bitten before — `BSplineHoleInCylinderWall` regresses a bug where a closed
  B-spline flattened to two coincident points and every polyline query saw an empty curve.
- **Both are emitted and the matcher fails to pair them.** Matching probes chord *midpoints*; two
  faces sampling one curve at different phases put their midpoints up to a chord apart, which is
  0.0103 cm here — well above the 1e-8 epsilon and well below the 2.4e-7 the curves actually agree
  to. Note those two numbers are on opposite sides of the epsilon, which is worth thinking about.

**3. Rename or redefine `maxGap`.** It is the largest distance from a chord to the nearest chord of
a *different face* — "how isolated is the loneliest chord". It has been read as "how far apart are
the two faces at the seam" in `CloseShape`'s error text, the harness line, `--json`, and three
documents including the previous hand-over. It has no dependence on the matching epsilon at all,
which is the cleanest proof that it is not the seam separation. Either make it that quantity or
call it what it is.

**4. Only then re-scope Phase 2 (adjacency-based exact trims).** It is still the right destination
— it is what would remove the sliver of §11 rather than label it — but the evidence that motivated
doing it *for these six parts* has evaporated, and the reason they fail may be an ordinary bug two
levels below it. When it is scoped, `RayHit::onTrimBoundary` is its ready-made acceptance test: on
a model with exact adjacency, that flag should never fire.

**K4 and K6 remain untouched, and K6 has no reproducer at all now.** Item 1 took its only one away
— the `cyl_cross_cyl` point was a trim-boundary defect and the quadratic roots there are exact to
the last printed digit. K4 (`bsplineSampleRecursive`'s degenerate-chord recursion probes only the
parametric midpoint) never had one.

**One deliberate non-decision to inherit.** `capacityIsExact()` still returns false for any wire
trim, although the contour capacity now lands at 1e-12 on every closed part. In this codebase
"exact" means analytically exact, and Gauss-Legendre on `sin(u(t))` is at machine precision but not
that; for a B-spline trim the pcurve approximates the true seam however the integral is done.
Whether to introduce a weaker, honest predicate is a real question, left open rather than settled by
flipping a flag the measurements do not license.

Run the oracle gate **before and after every item** so each effect is attributable, and update the
docs and commit as you go.

## Baseline to beat

- `ctest -R BVHSurfaceSolid` green, **59 cases**.
- Oracle gate: fixtures **8/9** (`tube_window` fails), Bagger **5/12**.
- `contains` disagreements outside tolerance: **0** on fixtures, **2** on Bagger.
- Capacity: **≤1e-11 relative on every navigable part**; 1e-4 on the six open Bagger cylinders and
  3e-5 on `tube_window`, all of which are open-surface-set.
- BVH == `_Loop` bit-identical everywhere.
- Direction-independence sweep: **0** direction-split points in 143000 over the 13 `Reliable` parts.

Gate totals and `contains` counts are separate numbers. Do not quote one without the other.

## Five things that will not be obvious

**The §4.2 direction-independence sweep is blind to the re-shoot.** Its metric is
`ContainsAlongDirection`, which documents that it bypasses the policy `Contains` applies. So the
previous session's 1-in-143000 and this session's 0-in-143000 are two samples of the same rate, not
a before and an after. To measure the re-shoot, compare `Contains(p)` against
`ContainsAlongDirection(p, kFixed)` — on a `Reliable` solid the latter *is* the old `Contains`.

**The gate's timing column cannot resolve anything below a few per cent.** Its `contains` numbers
moved by −29% to +46% between two runs of the same code path on this machine. Never quote it for a
micro-optimisation; write a dedicated loop with a fixed point set instead.

**`maxGap` has a resolution and it is printed next to it, and it also does not mean what it says.**
The resolution caveat still holds — a rim is a polyline at `kArcSamples = 24` per turn, so two faces
sampling one shared *curved* edge at different phases differ by up to the chord sagitta
`r(1-cos(π/24)) ≈ 8.6e-3·r` cm, which `rimChordResolution` reports. On top of that, see item 3
above for what the number actually is.

**The fast path is licensed by a measurement, not by an argument.** `Contains` takes a single parity
shot when the solid is `Reliable`, votes over five directions otherwise, and now also votes when the
single shot rested on a trim-boundary tie-break. **Any change to the closure criterion must re-run
the §4.2 sweep**, because a criterion that succeeds more often moves solids onto that path. The
probe is a hundred lines against the built library (§6 recipe); rewrite it rather than hunt for it.

**Bit-identical is sometimes the answer and sometimes a red flag.** Item 1's gates *had* to be
bit-identical and were, which is what said the re-shoot had not leaked into anything the gate can
see; its real effect had to be measured separately, and was. Item 3's gate change moved columns and
had to be labelled as a gate change so no one reads it as a kernel improvement.

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
  `GetSurfaceRecords()` is public, so a probe can inspect a solid's exact trim data without
  touching the kernel — that is how §12 was measured.
- `--skip-convert` reuses `<workdir>/db`, which saves minutes — but only when the *converter* did
  not change. Copy a finished workdir and re-run with `--skip-convert` to re-measure a kernel
  change; that is the whole before/after loop.
- The gate exits non-zero whenever any part fails, which is the normal state. Read the counts, not
  the exit code.

## Commands

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH ctest -R BVHSurfaceSolid

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
```

To compare two gate runs, diff `<workdir>/gate.json` with the timing fields (`timing*`,
`*Seconds` — including `closeShapeSecondsMesh`/`closeShapeSecondsSurface`, which are easy to miss)
stripped: everything else, including the offender lists and the query checksums, is deterministic
and must match when a change is meant to be inert.

To see the rim numbers for one part without a full gate run:

```bash
$B/stage/bin/o2-bench-detectorsbase-solid-harness --db <workdir>/db --parts <part> \
    --points 1 --rays 1 | grep -E 'navigation|rim gap'
```

To check that no face was dropped — the first thing to run on any new model, and a two-line script:
compare `nFaces` in `<workdir>/oracle/answers_<model>_<part>.json` against `nSurfaces` in
`<workdir>/gate.json`.
