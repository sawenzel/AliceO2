# NEXT — session-start instruction for the BVHSurfaceSolid work

This file is the current hand-over. Whoever finishes a session should **rewrite it** for the next
one. Last updated 2026-08-01, after items 1-3 of the previous hand-over.

Branch `swenzel/bvhsurfacesolid`. The tree is clean; everything below is committed.

**Items 1-3 are done, and item 4 (re-scope Phase 2) is done as a scoping, below.** Item 2 turned
out not to be a pairing bug at all: the rim was never emitted, because the B-spline sampler
replaced a 179-pole curve with a straight line. That is finding **K4**, which the last hand-over
recorded as having no reproducer.

## The one result that matters

**Every oracle disagreement outside tolerance, on both corpora, on every column, is now zero — 35
→ 0.** Including on the seven parts that still fail. `contains` 2 → 0, `distin` 17 → 0, `distout`
16 → 0. Those parts fail on two things now and neither is a wrong navigation answer:

1. their closure check, on a residual open boundary of 0.3-0.9 cm per part;
2. their capacity, off by 1e-4 to 7e-4 relative.

And **the capacity column did not move by a single bit** through all of this, while those parts'
faces went from 0.75 cm apart to 5e-05 cm and their navigation became exact. So the baseline's
attribution of the capacity error to "the surface set is open" is void. It is the only numeric
column still failing anywhere, and it is now unexplained.

## Read first, in this order

1. [`TolerancePolicy.md`](TolerancePolicy.md) **§13** — this session's measurements, including
   §13.3 (why the rim vanished), §13.4 (why the fix looked like a regression), and §13.7 (the two
   premises that died).
2. [`CodeReview_Fable.md`](CodeReview_Fable.md) **§20** — the same, with the next steps.
3. The last handoff note at the end of [`BVHSurfaceSolid.md`](BVHSurfaceSolid.md).

Note that §12's *conclusion* (the defect is in rim extraction, not the geometry) survived and was
right; its *inference about which rim* did not. §9.1 is still corrected by §12.3.

## Do — the two threads, in either order

**1. Phase 2 — adjacency-based exact trims. It finally has its evidence, and it is now the only
thing that can close the remaining 0.3-0.9 cm.**

§13.8 swept `kBSplineFlatness` on `BoomCylinderInner` and the answer is not what it looks like.
The disagreement between the two faces *is* mostly chording — the isolation falls from 4.8e-05 to
2.6e-06 cm as the flattening tolerance goes 1e-5 → 1e-8 — but **the open length gets worse, from
0.325 cm to 0.984 cm**, because the match band is built from the sagitta and shrinks faster than
the disagreement does.

So tightening the flattening is not the lever, and the honest reading of §13.4's band is that it
under-estimates the polyline-to-polyline disagreement by about an order of magnitude: the sagitta
bounds each polyline against *its own* curve and says nothing about the other face's curve, which
is the term that dominates. It works at the shipped tolerance because the two happen to be within
an order of magnitude — a coincidence of scale, not a criterion. Do not paper over that with a
fudge factor on the band.

What removes the term is the two faces deriving their trims from one shared object. Phase 2
(`CodeReview_Fable.md` §4.2) does that by construction: both faces reference the *same* neighbour
surface, so their trims agree identically. Every Bagger seam is cylinder∩cylinder, which §4.2 lists
as exactly solvable. The floor the sweep does not go below — about 2.6e-06 cm — is a candidate for
the genuine difference between the two pcurves, and is the part no amount of sampling would ever
have removed.

`RayHit::onTrimBoundary` remains its ready-made acceptance test: on a model with exact adjacency
that flag should never fire. And there is now a second one — on a Phase 2 model the closure check
should pass at the *declared* tolerance without any sagitta widening at all.

**2. The capacity of the six cylinders — the last failing number, and now unattributed.**

+4.9e-4 on `BoomCylinderInner`, +6.9e-4 on `StickCylinderInner`, -3.8e-4 on `StickCylinderOuter`,
and -3.0e-5 on `tube_window`; ≤1e-11 on every part that passes. `integrateOverCurveTrimByParts`
runs Gauss-Legendre on the curves themselves and never touches the flattened polyline, which is
why nothing this session moved it — and also why the sampler bug cannot have caused it. Candidates
worth measuring before theorising, in cost order:

- The **seam bridges**. `integrateOverCurveTrimByParts` bridges each join with a straight segment
  only when `deltaV != 0`. Face 0 of `BoomCylinderInner` is a wire whose curves meet across a full
  turn in u; check whether a join with `deltaV == 0` but a non-zero **u** jump is being skipped.
- The **contour of a wire that wraps the whole periodic domain**. Face 0's junction B-spline runs
  u = 0 → 2pi. Green's theorem on a periodic chart needs the loop to close in the chart, not just
  in 3D.
- Compare against `integrateOverCurveTrim` (the midpoint-rule integrator kept precisely as the
  independent check) on these six parts at high N. If the two agree, the defect is in the trim
  region, not the integrator.

`GetSurfaceRecords()` is public, so a standalone probe can rebuild a wire and integrate it both
ways without touching the kernel — that is how §13.2-13.3 were measured (§6 recipe, and add
`-I$HOME/alisw/O2/Detectors/Base/src` to reach the kernel header directly).

## Baseline to beat

- `ctest -R BVHSurfaceSolid` green, **61 cases**.
- Oracle gate: fixtures **8/9** (`tube_window` fails), Bagger **6/12**.
- **Oracle disagreements outside tolerance: 0 on every column of both corpora.** This is the number
  to defend; a change that improves a gate total while breaking this is a regression.
- Capacity: **≤1e-11 relative on every navigable part**; 1e-4 to 7e-4 on the six open Bagger
  cylinders and 3e-5 on `tube_window`.
- Open boundary: `BoomCylinderInner` 0.325 of 62.8 cm; `tube_window` 0.276 of 58.1 cm; worst rim
  isolation on any failing part **6.9e-05 cm**.
- BVH == `_Loop` bit-identical everywhere.
- Direction-independence sweep: **0** direction-split points in 154000 over the 14 `Reliable` parts,
  and 0 points where `Contains` differs from a single fixed-direction shot.

Gate totals and disagreement counts are separate numbers. Do not quote one without the other.

## Six things that will not be obvious

**A bit-identical numeric column can accompany a real change, and did.** The sampler fix moved
`cyl_cross_cyl` and `cyl_inter_cyl` from closed to open while every checksum, every `contains`,
every distance and the capacity stayed bit-identical. Reading the gate totals alone would have
called that a regression and reverted the session's main result. Diff the *columns*, not the
verdicts.

**The match band is no longer the declared tolerance, and `SetModelTolerance` is now nearly
inert on curved seams.** A chord is matched within `rimEpsilon + own sagitta + partner sagitta`.
On a 24-sample arc the sagitta is ~8.6e-3·r cm, which dwarfs any declared tolerance. That is
deliberate (§13.4) and it is also the reason §12.3's epsilon sweep saw nothing move.

**The non-manifold test is deliberately on a different, stricter band.** Do not "unify" them.
Widening it to the sagitta turned `Base`, `Boom`, `Stick` and `BucketLink2` non-manifold, because a
corner is where a third face's rim legitimately comes within a chord's length. Neither corpus has
a coincident-face part, so the test has no measured true positive to calibrate against.

**`maxRimIsolation` is not a seam width, and never was.** It is how alone the loneliest chord is.
It was called `maxGap` and read as a separation in four places for two sessions. Renamed
everywhere, including the `--json` key; nothing else in the report changed.

**The gate's timing column cannot resolve anything below a few per cent**, and safety timings moved
by 3x between runs of identical code this session. Never quote it for a micro-optimisation; write a
dedicated loop with a fixed point set instead.

**Any change to the closure criterion must re-run the §4.2 sweep**, because a criterion that
succeeds more often moves solids onto `Contains`'s single-shot fast path — `BucketLink1` moved onto
it this session. The probe is a hundred lines against the built library (§6 recipe); rewrite it
rather than hunt for it.

## Traps, each of which cost time before

- **The gate needs the O2 environment in the same shell.** `runOracleGate.py` sets the OCC python
  itself but not `LD_LIBRARY_PATH` for the harness; run `eval "$(alienv printenv ...)"` and export
  `LD_LIBRARY_PATH` in the *same* command, or step 2 dies on `libarrow_acero.so`.
- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or `ctest` and the harness
  silently use stale installed libraries.
- The test binary is at `$B/stage/tests/`, not `$B/stage/bin/`; the harness is at `$B/stage/bin/`.
- pythonOCC needs the alibuild **python3.10**, not the system 3.12.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`; use a scratch folder.
- For a quick kernel probe, compile a standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase`
  (recipe in `TolerancePolicy.md` §6, and add `-lGeom`). Add
  `-I$HOME/alisw/O2/Detectors/Base/src` and the probe can build `Curve2D`/`CurveWire` objects
  directly and call the samplers — that is how §13.3 was measured.
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
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
ctest -R BVHSurfaceSolid

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
```

To compare two gate runs, diff `<workdir>/gate.json` with the timing fields (`timing*`,
`*Seconds` — including `closeShapeSecondsMesh`/`closeShapeSecondsSurface`, which are easy to miss)
stripped: everything else, including the offender lists, the query checksums and
`navigation.rimDetail`, is deterministic and must match when a change is meant to be inert.

To name a part's rims — which loop of which face is open, and where its worst chord is:

```bash
$B/stage/bin/o2-bench-detectorsbase-solid-harness --db <workdir>/db --parts <part> \
    --points 1 --rays 1 --rims | grep -E 'navigation|rim isolation|    rim'
```

To check that no face was dropped — the first thing to run on any new model: compare `nFaces` in
`<workdir>/oracle/answers_<model>_<part>.json` against `nSurfaces` in `<workdir>/gate.json`.
