# Stream X — the sub-patch BVH, and safety that may stop before any leaf

Date: 2026-08-09. Follows [`Stream_S_SafetyBVH.md`](Stream_S_SafetyBVH.md) (the safety traversal
this refines) and [`Stream_P_RepresentationBench.md`](Stream_P_RepresentationBench.md) §2.2 (the
per-candidate cost that says why fewer candidates is the lever). Companion to
[`Tutorial.md`](Tutorial.md).

**In one paragraph.** The BVH used to hold **one conservative box per surface**, and for a swept
quadric that box is enormous — a full cylinder's box is the box of its two full rim circles, a
sphere's is the whole ball — so every ray through the box paid an analytic trim-wire intersection
that mostly reported nothing. Each surface now contributes **several exact sub-patch cover boxes**
(`BoundedSurface::appendCoverBoxes`, chunks of ≤ π/4 per angle, boxed by closed-form sinusoid
extremisation — no sampling slack), and every traversal dedups surfaces through a per-query epoch
stamp so each is still tested **at most once**. Measured on the three loadable large ALICE3 exact
parts, the candidate sets shrink by **half to two thirds** — `ST1829909_002` (965 patches):
Safety 8.63 → **3.94** patches/call, `DistFromOutside` 4.2 → **1.7** (unpruned 10.4 → 3.4) —
and wall-clock follows where geometry dominates: transport −8 %, `DistFromInside` −11 % on that
part, Safety −24 % on `ST0923290_019`; `ST1829909_003` is flat within noise because the surviving
candidates are exactly the expensive near patches. Every answer is **bit-identical** to the loop
twins (0 disagreements: 1024-point kernels on all three parts, 112 unit cases, fixtures gate
10/10 surface parts oracle-clean, Bagger gate — see §5). On top of the exact path, Safety gained
two things from the latest VecGeom master's tessellated solid: an **anchor seed** (the upper
bound from 24 precomputed on-patch points opens the traversal, VecGeom's `fTestPoints`), and an
opt-in **approximate mode** (`SetSafetySlack`) whose best-first traversal may return the running
*lower bound* as soon as it is within the slack of the seed — a far point is answered by the
root box alone: **39 µs → 127 ns, zero patches evaluated**, and never below (1 − slack) of the
exact answer by construction, 0 violations in 18 000 checks.

Everything below is measured on `swenzel/bvhsurfacesolid` (uncommitted at measurement time),
single-threaded, on the same 10-core aarch64 box as Streams P/S, with other tenants; the
before/after pairs are a same-session A/B via `git stash`, not a comparison against Stream S's
older absolute numbers.

---

## 1. What changed, and why it is sound

### 1.1 Cover boxes per family

`BoundedSurface::appendCoverBoxes(std::vector<CoverBox>&)` — a new virtual whose boxes must
jointly contain **two** point sets:

* the **trimmed patch** (every `appendIntersections` hit and every on-surface point) — what makes
  skipping a surface the ray crosses no box of *complete*;
* the **realisation set of `distanceSqToPatch`** — what makes pruning a subtree on its box
  distance *sound for Safety*.

The second set is the larger one and differs per family, which is the whole design constraint:

| family | realisation set (verified in the kernel) | cover emitted |
| --- | --- | --- |
| planar / curved planar | the wire region itself | the single `conservativeBounds` box (default) |
| cylinder / cone | **inside the sweep window**: the query's own azimuth when in sweep (an on-axis query realises the same value at every azimuth), else one of the two seam edges | sweep window in π/4 phi chunks, each the exact box of its two rim arcs (a generator is a convex combination of its rim points) |
| sphere | the **whole sphere** — `(|p−C|−R)²` projects radially and never reads the trim | full (θ, φ) grid, 4 × 8 chunks |
| torus | the **whole torus** — the meridian distance never reads the trim | full (ring, tube) grid, 8 × 8 chunks; spindle torus (r > R) falls back to the single box |

The per-chunk boxes are **exact**: one helper, `sinusoidRange(a, b, t0, t1)` (endpoints widened
by the crest/trough where they fall inside the interval), composed per coordinate — two nested
applications for the sphere (θ after φ) and torus (tube after ring, using R + r cos v ≥ 0). No
sagitta fudge and no sampling, so the unit test can assert containment with a 1e-9 slack and the
lower-bound property with 1e-18.

**Why the sphere must cover the ball even for a small cap:** its `distanceSqToPatch` realises at
the radial projection, which for a query behind the cap is on the far side of the sphere. Boxing
only the cap would prune the subtree that realises the minimum and Safety would come out **too
large** — the one forbidden direction. `StreamX_CoverBoxesAreATightLowerBoundEnvelopePerFamily`
pins exactly this on a polar cap.

### 1.2 One surface, once: the dedup is correctness, not optimisation

With several boxes per surface a ray can enter the same surface's leaves more than once, and a
second `appendIntersections` call duplicates its hits — **flipping crossing parity** and
corrupting the graze clustering. Every traversal (`visitRayCandidates`, `nearestCrossing`,
`visitPointCandidates`, `nearestPatchDistanceSq`, the approximate kernel) therefore opens a
`SurfaceVisitMarker` — an epoch-stamped thread_local array, no clearing and no allocation on the
hot path — and hands each surface on exactly once. Traversals never nest, which is what makes one
stamp array per thread enough. `StreamX_CurvedFixturesStayIdenticalToTheLoop` pins the crossing
*lists* (same multiset, element by element) on the fixtures whose surfaces own many boxes.

The candidate counters and `CountBVHRayCandidates` now count **distinct surfaces handed on**, the
same quantity as before (one box per surface made the two readings coincide), so the numbers stay
comparable.

### 1.3 The anchor seed (VecGeom master, tessellated solid, stage 1)

`TessellatedImplementation` in current VecGeom master seeds its BVH safety query with the
distance to N precomputed surface points (`fTestPoints`) before traversing. The same trick here:
24 display-mesh vertices (each exactly on its patch) are kept at `CloseShape`, and
`nearestPatchDistanceSq` starts its running best from the nearest anchor instead of infinity, so
pruning works from the first node instead of after the first evaluated patch. The seed is
inflated by `(1 + 1e-12) + 1e-10 cm` — above every rounding of the anchor distance and of
`distanceSqToPatch`, and far below the `kBVHBoxTolerance = 1e-3` the winner-visit argument rests
on — so the winner is always still visited and the returned value **and index** are bit-identical
to the loop. The existing 2000-point × 9-fixture sweeps, the sabotage control and the tie-break
case all still pass unchanged (the sabotage control in fact got *stronger*: a single-patch sphere
now has 32 leaves to prune, so all 9 fixtures are sensitive to it, not 7 — the test was updated
to say so).

### 1.4 The approximate mode: stopping before any leaf

`O2BVHSurfaceSolid::SetSafetySlack(s)` (default 0 = exact; clamped to [0, 0.9]; only `Safety()`
consults it — `ComputeNormal` stays exact because a normal from the wrong patch is wrong in kind).
The kernel is a best-first traversal (the shape of `O2Tessellated::SafetyKernel`): the heap pops
nodes in increasing bound order, so the first popped bound that reaches `(1 − s)² ×` the current
upper bound (anchor seed or best evaluated patch) underestimates **everything not yet expanded**,
and the traversal returns it on the spot. What is returned is always a *lower bound* on the exact
answer — under-estimating is the direction a safety is allowed to err — and never less than
`(1 − s)` of it, because that is the stop condition itself. A far point stops at the root:
**zero patches evaluated**.

```
ST1829909_002 (965 patches), 2000 points per regime, contract checked against exact per point
regime               slack 0.00          slack 0.05          slack 0.20
near (inside box)    21939 ns  3.89 p/c  17974 ns  2.03 p/c  14429 ns  1.31 p/c
mid (3 box radii)    31962 ns  4.67 p/c  22800 ns  1.55 p/c   5814 ns  0.29 p/c
far (50 box radii)   38997 ns  5.40 p/c    127 ns  0.00 p/c     61 ns  0.00 p/c
violations of [ (1-s)·exact, exact ]: 0 of 18000
```

The far-regime exact column is the pathology the mode exists for: far away, *every* top-level
box bound looks similar and the exact traversal wanders before it can prove the minimum; the
approximate mode answers from the root in 127 ns — 307× — while a navigator loses at most 5 % of
step length.

---

## 2. What it costs

* **Build**: one BVH over ~10–60× more primitives for quadric-heavy solids (`ST1829909_002`:
  965 surfaces → 5.4 k cover boxes). `CloseShape` on that part: 9.97 s before, 9.63 s after —
  unchanged within noise; the closure validation dominates it, not the BVH build.
* **Memory**: the tree grows with the primitive count; the measured heap delta of the whole load
  + close on `ST1829909_002` moved from ~14.4 MB to ~14.8 MB.
* **Per-candidate cost rises** where candidates fall: the survivors are the truly-hit patches,
  which take the expensive trim-wire branch (`Stream_P` §2.2's `Curve2D::closestPoint`). That is
  why `ST1829909_003`'s wall-clock is flat although its candidate sets halved — the eliminated
  candidates were the cheap misses. The candidate counters, not ns/call, are the honest measure
  of what the *tree* does; ns/call is what the user feels, and both are reported.

## 3. Measured: the same three ALICE3 parts, before → after (same session, git-stash A/B)

`o2-bench-detectorsbase-xray --surfaces … --perf --raster 6 --beams 3 --perf-points 1024
--perf-rays 1024 --perf-passes 5`, median of 5 passes, warm cache.

| | `ST1829909_002` (965 p) | `ST1829909_003` (236 p) | `ST0923290_019` (44 p) |
| --- | --- | --- | --- |
| Contains ns | 12268 → **11312** (−7.8 %) | 12754 → 12925 (+1.3 %) | 2552 → 2555 (0) |
| Safety ns | 21712 → **20879** | 45491 → 45903 (+0.9 %) | 6439 → **4896** (−24 %) |
| Safety candidates/call | 8.63 → **3.94** | 4.92 → **3.82** | 3.47 → **1.96** |
| DistFromOutside ns | 8784 → 8578 | 7645 → 8102 | 4345 → 4213 |
| distout candidates/call (unpruned) | 4.2 → **1.7** (10.4 → 3.4) | 2.2 → **1.5** (3.7 → 2.6) | 2.9 → **1.9** (5.8 → 3.6) |
| DistFromInside ns | 9747 → **8667** (−11 %) | 20632 → 20566 | 3901 → 3819 |
| transport ns/ray | 19773 → **18202** (−7.9 %) | 25916 → 26517 (+2.3 %) | 15489 → **14872** (−4 %) |
| kernel disagreements vs _Loop | 0 → 0 | 0 → 0 | 0 → 0 |

(`ST1829909_01`, the largest sidecar, does not load before or after: its sidecar has a
pre-existing wire-join gap of 5.4e-6 cm over the 1e-6 tolerance, rejected by `LoadSurfaceSolid`
validation. Not this stream's doing; noted for whoever owns the converter.)

## 4. Tests

Four new cases (`StreamX_*`, `testBVHSurfaceSolid.cxx`), written before the implementation:

* `CoverBoxesAreATightLowerBoundEnvelopePerFamily` — per family, on skew frames, partial sweeps
  and wire trims: >1 box for the curved families, every trimmed-patch sample inside a box, and
  the nearest-box distance never above `distanceSqToPatch` at 400 random points ×3 scales.
* `RaysThroughEmptyBoxRegionsReachNoPatch` — corner ray past a sphere, axial ray through a tube,
  ray behind a quarter cylinder's back: **0 candidates** where the old single box forced a paid
  intersection (the tube's axial ray: 3 → 2, only the caps remain).
* `CurvedFixturesStayIdenticalToTheLoop` — distances, containment and the crossing multisets on
  sphere/torus/capsule/cone/wire-trim fixtures.
* `ApproximateSafetyIsBoundedCheapAndOptIn` — the [(1−s)·exact, exact] contract on all nine
  nearest-patch fixtures, exactness restored at slack 0, and the far-point zero-leaf claim.

One updated case: `StreamS_BreakingThePruningBoundIsCaught` — see §1.3.

## 5. Reproducing it

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-xray
ctest -R BVHSurfaceSolid                                  # 112 cases
$B/stage/tests/o2-test-detectorsbase-BVHSurfaceSolid --run_test='StreamX_*' --log_level=message
$B/stage/bin/o2-bench-detectorsbase-xray --self-test      # 34 checks

# the gates (surface totals must stay 0/0/0/0)
cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/x_fix --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/x_bag \
    --model scripts/geometry/STEP_examples/Bagger.step

# the per-part numbers (ALICE3 sidecars, e.g. from a previous /tmp/S_a3-style conversion)
$B/stage/bin/o2-bench-detectorsbase-xray --surfaces <dir>/surfaces_ST1829909_002_0_1_1_36.bin \
    --perf --raster 6 --beams 3 --perf-points 1024 --perf-rays 1024 --perf-passes 5
```

The approximate-safety probe used for §1.4 is a standalone `.cxx` against the stage lib (the
NEXT.md quick-probe pattern): load a sidecar, `CloseShape`, time `Safety()` at slacks
{0, 0.05, 0.2} over three point regimes and check the contract per point.
