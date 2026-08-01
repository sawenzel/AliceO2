# NEXT — session-start instruction for the BVHSurfaceSolid work

This file is the current hand-over. Whoever finishes a session should **rewrite it** for the next
one. Last updated 2026-08-01, after a review session plus wave 0 of the new plan.

Branch `swenzel/bvhsurfacesolid`. The tree is clean; everything below is committed (`ccb919a877`).

**The plan changed this session.** Read [`Workstreams.md`](Workstreams.md) first — the work is now
six streams meant to run in parallel, and v1's Phase 2 (adjacency trims) is superseded. Wave 0
(Stream C) is done. Waves 1-3 are open.

## The one result that matters

**The capacity column is solved, and it was quadrature after all.** The residual on the six Bagger
cylinders had survived five sessions and was recorded as "unexplained". It is Gauss-Legendre
spanning a B-spline's knots, and it hid because a *second* defect — the interval count taken from
the curve's endpoints, which is zero for every closed hole loop — made the obvious experiment
incapable of detecting it. Full record: [`Stream_C_Hygiene.md`](Stream_C_Hygiene.md).

Capacity relative to OpenCascade: `BoomCylinderInner` 4.88e-04 → **2.8e-07**; `StickCylinderInner`,
`BucketCylinderOuter` and the `tube_window` fixture now **pass outright**; the other three
1.1e-04…3.8e-04 → **1.3e-06**. Raising the quadrature order to 40 leaves those bit-identical, so
what is left is the trim data, not the integration.

**Every part that now fails, on either corpus, fails on navigability and nothing else.** The gate
has collapsed onto one question and it is the closure check's.

## Read first, in this order

1. [`Workstreams.md`](Workstreams.md) — the six streams, the contract every stream obeys, the
   file-ownership matrix, the launch order. §1 is not optional.
2. [`CodeReview_Fable_v2.md`](CodeReview_Fable_v2.md) — the second deep review: findings N1-N6,
   the feasibility re-scoring, and §5 on the CSG direction.
3. [`CSG_Pipeline.md`](CSG_Pipeline.md) — B-rep → CSG: prior art, the four tiers, the acceptance
   test, and §6's measurement of what ALICE3 actually contains.
4. [`Stream_C_Hygiene.md`](Stream_C_Hygiene.md) — what wave 0 did and the two traps it left.

`CodeReview_Fable.md` (v1) remains the register of record for the K/S/C findings. **Its §9-10 plan
is superseded** — do not execute it.

## Do — wave 1, three streams, near-zero overlap

Run them in parallel, in separate worktrees. Briefs are in `Workstreams.md`.

- **Stream A — the CSG / simplification pipeline.** New Python package. Start with the *census*
  (`CSG_Pipeline.md` §8 step 1): which parts are quadric-only, how many faces, how many concave
  edges, which match a `TGeoTube`/`TGeoBox`/`TGeoPcon` template. **Do not build the emitter before
  that table exists** — it decides how much of tiers 1-3 is worth building.
- **Stream E — scale and the gate at ALICE3 size.** First item is the translate/scale sweep, and it
  is the cheapest experiment on the board with the largest possible consequence: every measurement
  on this branch was taken a few centimetres from the origin, against absolute constants
  (`kBVHBoxTolerance` 1e-3 cm, `kWireJoinTolerance` 1e-5 cm). If the ladder is not
  column-identical at (0,0,400), a lot of recorded numbers need re-qualifying.
- **Stream F — sidecar v3 edge identity.** The closure check decides a topological question with a
  proximity query, and the previous hand-over reached the end of that road honestly. The converter
  has the `TopoDS_Edge` identity and throws it away. This is what makes "watertight" checkable
  rather than tuned, and it is what the six remaining failing parts are waiting for.

Wave 2 (Streams B and D) starts after A's census and E's sweep, because both can change what is
worth building.

## Baseline to beat

- `ctest -R BVHSurfaceSolid` green, **64 cases**.
- Oracle gate: fixtures **8/9** (`tube_window`), Bagger **6/12**.
- **Oracle disagreements outside tolerance: 0 on every column of both corpora.** Still the number
  to defend; a change that improves a gate total while breaking this is a regression.
- Capacity: ≤1e-11 on every part without a B-spline trim; ≤2.8e-07 on three of the six cylinders
  and 1.3e-06 on the other three.
- Open boundary and rim isolation unchanged from the previous hand-over — wave 0 did not touch them.
- BVH == `_Loop` bit-identical everywhere.

Gate totals and disagreement counts are separate numbers. Do not quote one without the other.

## Things that will not be obvious

**A refuted experiment is not a refuted hypothesis.** The `kContourMaxSpanU` sweep said "no" to the
right answer, and a whole section of the review was rewritten around that "no" before the localiser
showed the knob could never have reached the curve. When an experiment fails to move a number,
check that it was *capable* of moving it.

**`BRepGProp::VolumeProperties` on a single planar face returns 0.** Not documented in OCCT. Any
per-face comparison must restrict itself to curved faces and account for planes by difference —
`BasePin` is the self-check (20.943951 curved + 10.471976 caps = 31.415927).

**A bit-identical numeric column can accompany a real change, and has, twice.** Diff the columns,
not the verdicts. Strip `timing*`, `*Seconds` and `nsPerCall` from `gate.json` and everything else
must match when a change is meant to be inert.

**The match band is not the declared tolerance**, and `SetModelTolerance` is nearly inert on curved
seams — a chord is matched within `rimEpsilon + own sagitta + partner sagitta`. Deliberate; see the
previous hand-over and `TolerancePolicy.md` §13.4.

**The non-manifold test is on a different, stricter band. Do not "unify" them.**

**Any change to the closure criterion must re-run the §4.2 direction-independence sweep**, because
a criterion that succeeds more often moves solids onto `Contains`'s single-shot fast path.

**The gate's timing column cannot resolve anything below a few per cent** and moved 3x between runs
of identical code. Write a dedicated fixed-point loop instead.

## Traps, each of which cost time before

- **The gate needs the O2 environment in the same shell.** `runOracleGate.py` sets the OCC python
  itself but not `LD_LIBRARY_PATH` for the harness; `eval "$(alienv printenv ...)"` and export
  `LD_LIBRARY_PATH` in the *same* command, or step 2 dies on `libarrow_acero.so`.
- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or `ctest` and the harness
  silently use stale installed libraries.
- The test binary is at `$B/stage/tests/`, the harness at `$B/stage/bin/`.
- pythonOCC needs the alibuild **python3.10**; OCCT and pythonOCC live under
  `$ALIBUILD_WORK_DIR/ubuntu2404_aarch64/{OCCT,pythonOCC}/latest`.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`; use a scratch folder.
- `--skip-convert` reuses `<workdir>/db` — valid only when the *converter* did not change. Copy a
  finished workdir and re-run with it; that is the whole before/after loop for a kernel change.
- The gate exits non-zero whenever any part fails, which is the normal state. Read the counts.
- For a kernel probe, compile a standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase -lGeom`
  plus `-I$HOME/alisw/O2/Detectors/Base/src` (`TolerancePolicy.md` §6). Wave 0's per-face and
  Monte-Carlo probes were both written that way in a few minutes each.

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

To name a part's rims — which loop of which face is open, and where its worst chord is:

```bash
$B/stage/bin/o2-bench-detectorsbase-solid-harness --db <workdir>/db --parts <part> \
    --points 1 --rays 1 --rims | grep -E 'navigation|rim isolation|    rim'
```

To check that no face was dropped — the first thing to run on any new model: compare `nFaces` in
`<workdir>/oracle/answers_<model>_<part>.json` against `nSurfaces` in `<workdir>/gate.json`.
