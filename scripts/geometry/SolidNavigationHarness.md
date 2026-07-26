# Solid navigation harness, test-part database, and BVH distance navigation

Plan and (once implemented) usage documentation for the reusable validation / performance
comparison harness of `O2BVHSurfaceSolid` versus `O2Tessellated`, the CAD-derived test-part
database it runs on, and the BVH implementations of `DistFromOutside` / `DistFromInside`.

**Why this is a separate document from `BVHSurfaceSolid.md`.** That file is the durable milestone
checklist of the exact-surface project: a long list of surface/converter capabilities, each closed
once. This work is different in kind — it produces a *tool* that stays in daily use, so the document
has to become permanent usage documentation (how to build the part DB, how to run the benchmark, how
to read the output, how to profile it under `perf`), not a closed checkbox. Keeping it here also
keeps `BVHSurfaceSolid.md` readable.

It covers three of that file's open milestones, which are cross-linked from there but planned here:

- "Reusable test harness for validation and performance measurement as part of the repo"
- "Implement `DistFromOutside`" and "Implement `DistFromInside`"
- "Code optimization pass" (the concrete optimization wishes recorded in that milestone)

Status: **Steps 0-6 implemented and measured (2026-07-26).** Results in the section at the end.

## Why these three are one piece of work

They are not independent. The harness is a measuring instrument and needs a subject; the distance
functions are the subject; the optimization wishes are only verifiable with the instrument. Done in
one pass they reinforce each other:

```
part DB  ->  harness  ->  baseline measurement  ->  BVH DistFrom*  ->  re-measure
```

The rename of today's loop-based `DistFromOutside` to `DistFromOutside_Loop` yields, from a single
edit, both a correctness oracle and the performance baseline the new implementation is measured
against.

The risk of combining them is that a bug in the harness masks a bug in navigation. Two guards, both
already established practice in this project, are non-negotiable:

1. the harness is exercised first against **analytic** references (`TGeoBBox`, `TGeoTube`, ...) in
   the unit test, where the truth is known independently of any mesh;
2. **`BVH == _Loop`** is compared on every sample. This is independent of the tessellated reference
   entirely and separates BVH/traversal bugs from surface-kernel bugs.

## Ground rules for validation (do not re-litigate these)

- **The tessellated solid is a reference, not the truth.** It is a chorded mesh and it is
  systematically *inscribed*, so on curved parts `DistFromOutside` measured from outside is
  systematically *longer* on the exact solid, by roughly the chord sagitta. That is expected,
  quantifiable, and is a result to report — not a bug.
- **Classify, do not just count, mismatches.** For every disagreeing sample, measure `Safety` at the
  point and report whether it falls inside the mesh-precision band. This is the method that made
  "the reference is wrong, not the exact solid" a claim rather than a hope in the 2026-07-26 ALICE3
  sweep; keep it.
- **Never compare `Safety` for equality against another shape.** Safety is only contracted to be a
  lower bound of the true distance, so two correct implementations may differ arbitrarily. Check the
  contract (`safety <= true distance`, `safety >= 0`) and compare against analytic distances in unit
  tests.
- **Compare against the accelerated `O2Tessellated`, never ROOT `TGeoTessellated`** — the latter
  fills concavities on thin non-convex shells (~27% false-inside was measured on ALICE3).
- Fixed seeds everywhere. A performance or accuracy number that cannot be reproduced exactly is not
  a measurement.

## Deliverables

| # | Artifact | Kind | Status |
| --- | --- | --- | --- |
| 0 | `scripts/geometry/checkSurfaceSidecars.C` renamed to `.macro` | rename + doc fix | done |
| 1 | `scripts/geometry/makeTestPartDB.py` | new, Python | done |
| 2 | `Detectors/Base/include/DetectorsBase/O2SolidHarness.h`, `Detectors/Base/src/O2SolidHarness.cxx` | new, C++ | done |
| 3 | `Detectors/Base/test/runSolidHarness.cxx` -> `o2-bench-detectorsbase-solid-harness` | new, C++ | done |
| 4 | `LoadFacetSolid(...)` in `O2SurfaceSolidIO.h/.cxx` | extension | done |
| 5 | BVH `DistFromOutside` / `DistFromInside` in `O2BVHSurfaceSolid.cxx` (+ `_Loop` twins) | change | done |
| 6 | Measurement results, recorded in this file; milestone checkboxes in `BVHSurfaceSolid.md` | doc | done |

**2026-07-26 implementation note (Steps 0-4).** All four built and smoke-tested against a live
`makeTestPartDB.py` run (`Bagger.step` + `as1-oc-214.stp` + `oTOF System V3-R92cm.step`, 19/19
parts paired, all 19 sidecars load via `checkSurfaceSidecars.macro`). `o2-bench-detectorsbase-solid-harness`
ran validation + timing across all 19 parts with `--loop-crosscheck` without a crash; sample
result on `as1_oc_214/bolt`: Contains 99.84% agree (2/1250 mismatches, both mesh-band, 0
unexplained), Safety lower-bound contract holds for both shapes (100%), BVH==`Contains_Loop` on
1250/1250 points. DistFromOutside/DistFromInside show the expected 100%-mesh-band mismatch pattern
against the tessellated reference (at the time of this note `O2BVHSurfaceSolid::DistFromOutside`/
`DistFromInside` were still the all-surfaces loop; superseded by the Steps 5-6 note below, and the
numbers here are therefore historical -- the ones to quote are in the Results section). One unrelated
latent build break was hit and fixed the same way as Step 0: a stray untracked
`scripts/geometry/foo.C` (an old ad-hoc conversion output from a prior session) also tripped
`O2ReportNonTestedMacros`; renamed to `foo.macro`. Editing `Detectors/Base/CMakeLists.txt` again
triggered the documented incremental-`ninja`/rootcling GLIBCXX failure; fixed by a full
`aliBuild build O2` per the existing build notes below -- `ctest -R BVHSurfaceSolid` stayed green
throughout.

**2026-07-26 implementation note (Steps 5-6).** Both distance queries now traverse the BVH via
`bvh->intersect<false, /*robust*/ true>` with a leaf lambda, sharing one `Impl::nearestCrossing<bool>`
template (the entering/exiting test is the template argument, as planned) and tightening `ray.tmax`
from inside the lambda. `DistFromOutside_Loop` / `DistFromInside_Loop` are public and share
`Impl::nearestCrossingLoop<bool>`, so oracle and implementation cannot drift apart in anything but
traversal order. Scratch buffers are `static thread_local` inside the templated functions, which
gives each of the four entry points (`Contains`, `Contains_Loop`, entering, exiting) its own buffer
by construction; `Contains` no longer allocates either. New public switches
`SetRayTMaxPruning`/`GetRayTMaxPruning` and `ResetRayCandidateCounter`/`GetRayCandidateCount` exist
purely so the benchmark can price the pruning; the counter is `thread_local`, the flag is not and
is documented as a between-passes knob.

Five test cases were added to `testBVHSurfaceSolid.cxx` (`DistanceBVHMatchesLoopOnAllFixtures`,
`DistanceSweepsMatchRootPrimitives`, `DistanceHardCases`, `RayTMaxPruningIsOptimizationOnly`,
`RayTMaxPruningKeepsNearTies`), and the harness front-end gained `--pruning-ab`, `_Loop` baseline
timings, and the extended `--loop-crosscheck`. `ctest -R BVHSurfaceSolid` green (**35 cases**, 0.80 s).

Two things worth carrying forward, both recorded rather than papered over:

- **Mutation-testing the new tests.** Scaling the pruning bound by 0.5 is caught loudly by the
  sweeps; scaling it by 0.999 is caught by *nothing* in the suite. That is a property of the
  algorithm, not a hole in the tests: a node is culled only when the ray enters its *box* beyond
  `tmax`, and a box is entered no later than its patch is hit, so a slightly-too-tight bound can
  only lose a candidate whose box entry falls in the razor-thin window between the bad bound and
  the true best. The correctness argument is therefore structural (see the comment on
  `nearestCrossing`), and `RayTMaxPruningKeepsNearTies` pins the geometry rather than the bound.
- **`Contains` vs `Contains_Loop` disagrees on non-manifold parts** (301 / 142500 over the DB, 295
  of them in the two `oTOF` parts). Pre-existing, unrelated to Step 5, detail in the Results
  section below.

## Step 0 — clear a latent build break first

`BUILD_TEST_ROOT_MACROS` is `ON` in this build, and `cmake/O2ReportNonTestedMacros.cmake` globs
**every** `*.C` in the source tree and hard-fails on any that is neither registered via
`o2_add_test_root_macro` nor listed in `cmake/O2RootMacroExclusionList.cmake`.
`scripts/geometry/checkSurfaceSidecars.C` (added 2026-07-26) is neither. It has not bitten yet only
because CMake has not reconfigured since it was added — and adding the Step 3 target *forces* a
reconfigure.

Fix: rename it to `checkSurfaceSidecars.macro` and update the references in `BVHSurfaceSolid.md`.
`.macro` is the established house convention for exactly this (40 such files already exist, e.g.
`Detectors/gconfig/src/*.macro`, `run/SimExamples/**/*.macro`). Any future ROOT macro added under
`scripts/geometry/` follows the same rule.

Acceptance: CMake reconfigures without `Macro ... should be tested`.

## Step 1 — test-part database

`scripts/geometry/makeTestPartDB.py` builds a directory of CAD parts held in **both**
representations, which is the harness's precondition: the same CAD solid as an exact surface sidecar
and as a triangle mesh. The converter already emits both in one run (`surfaces_<VOL>_<LID>.bin` and
`facets_<VOL>_<LID>.bin`), so the script mostly orchestrates and indexes rather than computes.

Behaviour:

- For each selected model, run `O2_CADtoTGeo.py <model> --exact-surfaces auto --mesh
  --surface-report <out>/surface_report.json` into `<db>/<model-slug>/`.
- Pair the two artifact families by `<VOL>_<LID>`; a part enters the DB only when both exist.
- Write `<db>/manifest.json`: per part `{id, model, volume, lid, surfaces, facets, nTriangles,
  bbox}`, plus a per-model summary (`volumes total / exact`) and the converter command line used.
- CLI: `--models`, `--output`, `--skip-existing`, `--force`.

**MVP scope: the three fast models only** — `Bagger.step`, `as1-oc-214.stp`,
`oTOF System V3-R92cm.step` (~19 parts, minutes rather than the ~30 min ALICE3 needs). They are a
genuinely good MVP set, not just a cheap one:

| model | why it earns its place |
| --- | --- |
| `Bagger.step` | 12/13 exact and fully analytic across plane / cylinder / sphere / torus |
| `as1-oc-214.stp` | 5/5 exact, entirely *recognized* NURBS-as-cylinder — exercises the newest converter path |
| `oTOF System V3-R92cm.step` | a 1505-surface all-planar part: the scaling case |

`ALICE3_CAD_pure.step` is one `--models` flag away and is added in Step 6, *after* the harness is
trusted — so the expensive sweep becomes the harness's first real job rather than a prerequisite for
writing it.

Generated binaries stay out of git, per project convention; the generator is what gets committed.

Acceptance: `manifest.json` lists ~19 parts; every listed sidecar loads via `LoadSurfaceSolid`
(cross-check with `checkSurfaceSidecars.macro`).

## Step 2 — harness core

`o2::base::harness` in `DetectorsBase/O2SolidHarness.h` / `src/O2SolidHarness.cxx`. Typed on
`TGeoShape*` throughout, deliberately: the same code then serves surface-vs-tessellated comparison
*and* surface-vs-ROOT-primitive checks in the unit test.

**Sampling** — `SampleSet`, deterministic from a seed and a bounding box:

- bulk points, uniform in the inflated bbox;
- boundary-band points, concentrated within a configurable band of the surface (these are where
  disagreements live, and uniform sampling finds far too few of them);
- inside-biased points, accepted by the reference `Contains`;
- rays: isotropic directions from outside points, **plus** directions aimed at random points inside
  the bbox, so `DistFromOutside` hit rates are not degenerate. A benchmark of a function that
  almost always misses measures the BVH's reject path and nothing else.

**Validation** — `validateContains`, `validateDistFromOutside`, `validateDistFromInside`,
`validateSafety`, each returning agreement counts, mismatch counts split into
*within-mesh-band* / *unexplained*, worst deviation, and a bounded list of worst offenders with the
point/direction printed so a failure is directly reproducible. `validateSafety` checks the
lower-bound contract only (see ground rules).

**Timing** — warmup, repeat count, identical sample order for both shapes, and a checksum
accumulated from the results so the optimizer cannot elide the calls. Reports ns/call and the
candidate/reference ratio.

Also reported per part, because they are half the story: surface count vs triangle count,
`CloseShape`/BVH build time, and BVH ray-candidate counts.

## Step 3 — front-end tool

`Detectors/Base/test/runSolidHarness.cxx`, registered with
`o2_add_executable(solid-harness COMPONENT_NAME DetectorsBase IS_BENCHMARK ...)`, giving
`o2-bench-detectorsbase-solid-harness`.

Options: `--db <dir>` (or explicit `--surfaces`/`--facets`), `--parts <pattern>`, `--points`,
`--rays`, `--seed`, `--only contains,distout,distin,safety`, `--loop-crosscheck`, `--pruning-ab`,
`--json <out>`, `--warmup`, `--repeat`.

The two cross-cutting flags:

- `--loop-crosscheck` also runs `Contains_Loop` / `DistFromOutside_Loop` / `DistFromInside_Loop` and
  requires **exact** agreement with the BVH paths (not agreement within a tolerance: both minimise
  over the same hits from the same kernels, so any difference is a traversal bug). This is the
  correctness guard that does not involve the mesh at all, and it is what the headline result in the
  Results section rests on.
- `--pruning-ab` re-runs the distance kernels with `SetRayTMaxPruning(false)` and reports the leaf
  candidate counts and ns/call both ways, plus how many answers changed (which must be none).

The `--only` single-kernel mode is the `perf record` entry point, which is the stated purpose of the
harness as a foundation for profile-guided optimization:

```
perf record -g o2-bench-detectorsbase-solid-harness --db <db> --parts oTOF --only distout --rays 200000
```

## Step 4 — facet loader in C++

`bool LoadFacetSolid(const std::string& file, O2Tessellated& solid)` in `O2SurfaceSolidIO.h/.cxx`,
mirroring the `LoadFacets` helper that `O2_CADtoTGeo.py`'s generated macro already contains, so the
harness does not re-implement the format. `O2SurfaceSolidIO` is the right home: it already owns
"read the converter's binary artifacts".

Amended during Step 6: a facet that `O2Tessellated::AddFacet` rejects as degenerate (collapsed
vertices) is skipped and counted in a warning rather than failing the load. Real CAD tessellations
contain a few slivers, and the first ALICE3 run lost four whole parts — one of them 211554
triangles — to a single one each. Truncated and unreadable files remain fatal.

## Step 5 — `DistFromOutside` and `DistFromInside` on the BVH

Both are implemented in this pass. They are structural twins (the entering/exiting test is a sign
flip on `dot(normal, direction)`), and distance validation is lopsided with only one of them.

- Today's all-surfaces loops become **`DistFromOutside_Loop` / `DistFromInside_Loop`**, public and
  kept, mirroring the existing `Contains_Loop`: oracle and baseline in one.
- New implementations traverse with `bvh->intersect<false, /*robust*/ true>` and a leaf lambda —
  **the existing bvh2 traversal with lambdas, as `O2Tessellated` does; no hand-written traversal.**
- **Ray tightening (the pruning gap noted in `O2Tessellated`).** `Node::intersect_robust` re-reads
  `ray.tmax` on every node test (`bvh/v2/node.h:117`, verified), so mutating `ray.tmax` from inside
  the leaf lambda as the best hit shrinks genuinely prunes subsequent nodes. Because leaf boxes are
  inflated by `kBVHBoxTolerance` and traversal is `float`, the new `tmax` must be rounded **up**
  conservatively (the `truncate_roundup` pattern in `O2Tessellated`) so a nearer hit can never be
  pruned away.
- Cheap early reject of the ray against the BVH root box under `stepmax`, as `O2Tessellated` does.
- **No allocation inside navigation.** The `appendIntersections` scratch buffer becomes a
  `thread_local` reused across surfaces and calls (capacity persists; `clear()` per surface). One
  distinctly named buffer per entry point, so no function can be re-entered through another's
  buffer.
- Fast-math hygiene: the direction is unit by TGeo contract and is not renormalized; comparisons are
  done on squared quantities where a square root is otherwise needed; small hot helpers are marked
  for inlining.
- `stepmax` is respected both as the initial traversal `tmax` and in the return value.

The same `thread_local` scratch treatment is applied to `Contains`, which allocates today — it is
the identical one-line change on an already-BVH-accelerated function, and it is measured by the same
harness run.

Unit tests in `Detectors/Base/test/testBVHSurfaceSolid.cxx`:

- `BVH == _Loop` for both functions on the existing box / tube / hollow-tube / cone / sphere /
  torus / capsule fixtures;
- comparison against `TGeoBBox` / `TGeoTube` / `TGeoCone` / `TGeoTorus` for outside and inside
  points over several directions;
- hard cases: grazing/tangent rays, rays through shared edges and vertices, rays starting exactly on
  a surface, and `stepmax` honoured (including a hit just beyond `stepmax`).

## Step 6 — measure, then expand

1. Run the tool over the three-model DB. Record: agreement and band classification per function,
   ns/call for surface vs tessellated, the ratio, memory/primitive-count reduction, BVH candidate
   counts **with and without** `tmax` pruning (a one-flag A/B that directly prices the optimization),
   and `_Loop` vs BVH speedup.
2. Only then extend the DB to `ALICE3_CAD_pure.step` and re-run, so the large-part numbers come from
   an instrument already trusted on small parts.
3. Write the results into the "Results" section below, tick the three milestones in
   `BVHSurfaceSolid.md` with a pointer to this file, and add a dated handoff note there.

### Expectation to set before measuring

The exact solid is likely to be **slower per query** than the tessellated one on large parts: an
analytic patch test costs far more than a ray/triangle test, and every curved patch is currently a
single loose conservative AABB (the open "Subdivision-BVH acceleration for curved surface patches"
milestone). The wins expected to be demonstrable are exactness and a large reduction in primitive
count and memory. Quantifying the gap — and thereby sizing that subdivision-BVH milestone with
measured numbers instead of intuition — is a legitimate headline result of this work, not a
disappointment. Report whatever the harness says.

## Model / effort guidance per step

Steps 0–4 are mechanical against a fully specified design. Step 5 carries the numerical subtlety
(tolerances, float BVH pruning, tangency and boundary cases) and Step 6 the interpretation. A
reasonable split is a fast model for 0–4 and a stronger one for 5–6, but the guards in
"Why these three are one piece of work" matter more than the model choice: they fail loudly rather
than silently.

**On automating the switch: a document cannot do it.** Claude Code takes its model from the user
(`/model`), from `settings.json`, or — for a subagent — from the `model:` field in the agent
definition or the `model` parameter of the `Agent` tool. Nothing reads an instruction out of a
markdown file and switches models on it, so a line like "use Opus for Step 5" in this document is a
note to the human, not an instruction to the harness. Two ways to get the effect:

- **Manual checkpoint (recommended here).** The step boundaries below are natural stopping points
  anyway, because each ends in a build-and-test gate. Switch with `/model` at the boundary.
- **Dispatch steps as subagents.** An orchestrating session that spawns each step with the `Agent`
  tool can set `model` per spawn, which *is* automatic. The cost is that every subagent starts cold
  and re-derives context this project has a lot of; for a five-step sequence with shared state, the
  manual checkpoint is usually the better trade.

Suggested checkpoints: after Step 1 (DB exists and load-checks), after Step 4 (harness builds and
runs against an analytic fixture), after Step 5 (`ctest -R BVHSurfaceSolid` green).

## Build and environment notes

Full detail is in the `BVHSurfaceSolid.md` handoff notes; the essentials:

- Toolchain: `cd /data/swenzel && eval $(./alibuild/alienv printenv O2/latest,ninja/latest,pythonOCC/latest,CMake/latest)`,
  then build in `/data/swenzel/sw/BUILD/O2-latest/O2`. pythonOCC needs the `PYTHONPATH` and OCCT
  `LD_LIBRARY_PATH` fixes; `alienv` needs `ALIBUILD_WORK_DIR=/data/swenzel/sw` when run from another
  directory.
- Incremental: `ninja Detectors/Base/all`; focused test: `ctest -R BVHSurfaceSolid --output-on-failure`.
  The `testMatBudLUT` failure in the same label is unrelated to this work.
- **Steps 2–4 change `Detectors/Base/CMakeLists.txt`.** A previous CMakeLists edit made the
  incremental `ninja` fail in rootcling dictionary regeneration; the fix was a full
  `./alibuild/aliBuild build O2`. Budget for that.
- Interpreted ROOT macros resolve `libO2DetectorsBase` from the *installed* prefix, not the ninja
  stage — `export LD_LIBRARY_PATH=/data/swenzel/sw/BUILD/O2-latest/O2/stage/lib64:$LD_LIBRARY_PATH`
  before running one against freshly built code. The compiled harness does not have this problem,
  which is one more reason it is the primary front-end.

## Results

### 2026-07-26 — three-model DB

Reproduce with:

```bash
python3 scripts/geometry/makeTestPartDB.py --output <db>          # 19 parts, 3 models
o2-bench-detectorsbase-solid-harness --db <db> --points 3000 --rays 3000 --seed 1 \
    --loop-crosscheck --pruning-ab --json step6_report.json
```

DB: `Bagger.step` (12 parts) + `as1-oc-214.stp` (5) + `oTOF System V3-R92cm.step` (2) = 19 parts,
1851 analytic surface patches against 48703 reference triangles. Seed 1. Single-threaded, one
machine, `--warmup 1 --repeat 3`.

**Correctness, independent of the mesh (`BVH == _Loop`).** This is the guard that isolates the
Step 5 traversal from everything else, and it is the headline correctness result:

| Query | agreement | worst deviation |
| --- | --- | --- |
| `DistFromOutside` | 57000 / 57000 | 0 (bit-identical) |
| `DistFromInside` | 28500 / 28500 | 0 (bit-identical) |
| `Contains` | 142199 / 142500 | — |

The distance twins agree bit-for-bit on every ray. The 301 `Contains` disagreements (0.21%) are
**not** from this work: 295 of them are in the two `oTOF` parts, which `CloseShape` reports as
non-manifold (32 and 9 non-manifold edges, i.e. coincident/duplicated faces). Parity-based
containment is genuinely order-dependent on such input — `oddCrossingParity` clusters near-equal
hits, and the BVH hands it the same hits in a different order — so this is a property of the
containment algorithm on degenerate geometry, not of the traversal. Worth its own investigation;
it is recorded as an open item in `BVHSurfaceSolid.md` rather than fixed here.

**Accuracy against the tessellated reference.** As predicted by the ground rules, the mesh is the
thing that is wrong wherever the two differ:

| Query | agree | mismatch, mesh-band | mismatch, unexplained |
| --- | --- | --- | --- |
| `Contains` | 95.74% | 1482 | 4588 |
| `DistFromOutside` | 75.10% | 13941 | 254 |
| `DistFromInside` | 47.52% | 14739 | 218 |

The low `DistFromInside` agreement is the inscribed-mesh effect stated up front, not a regression:
inside rays exit the exact solid later than they exit its inscribed chording, on essentially every
curved patch. The `Contains` unexplained count is dominated by `oTOF/oTOF_v2` (4063 of 4588) —
again the non-manifold part. Excluding it, unexplained `Contains` mismatches are 525 / 135000.

**Speed.** Per-call medians over the 19 parts (ns/call):

| Query | exact surface solid | tessellated | ratio (median, range) | BVH vs `_Loop` |
| --- | --- | --- | --- | --- |
| `Contains` | 1919 | 610 | 2.44x (0.24 – 51.98) | — |
| `DistFromOutside` | 1391 | 620 | 1.34x (0.21 – 23.71) | **2.03x** |
| `DistFromInside` | 1368 | 746 | 1.63x (0.25 – 45.17) | **1.21x** |
| `Safety` | 10132 | 2261 | 4.74x (0.27 – 141.04) | — |

So the exact solid costs ~1.3–2.4x per navigation query at the median, with a very wide spread in
both directions — it is *faster* than the mesh on the parts whose patch count is far below the
triangle count, and much slower on `oTOF/oTOF_v2` (1505 patches). This is the expected outcome
stated before measuring, and it sizes the open "Subdivision-BVH acceleration for curved surface
patches" milestone with numbers rather than intuition.

`Safety` is the clear outlier and the cheapest remaining win: it is still a plain loop over every
patch with no BVH at all (`O2BVHSurfaceSolid::Safety`), which is why it is 4.7x the mesh's cost and
by far the most expensive kernel in absolute terms. A BVH-guided nearest-patch search (the
`SafetySqToNode` kernel already exists in `bvh2_extra_kernels.h` and `O2Tessellated` already uses a
priority queue for exactly this) is the obvious next optimization.

**Ray `tmax` pruning A/B** (`--pruning-ab`, `DistFromOutside`, 57000 rays):

- answers identical in 57000 / 57000 cases — the switch is a cost knob, never a semantic one;
- leaf candidates handed to the patch intersector: 202946 pruned vs 276537 unpruned, i.e. the
  tightening removes **26.6%** of the analytic patch tests (1.36x reduction);
- wall time: median speedup **1.01x** (range 0.90 – 1.30).

That gap between work removed and time saved is the honest result: on parts this small the
per-patch test is cheap enough that a quarter fewer of them does not show up above the traversal
and measurement noise. The optimization is worth keeping — it is free at runtime and its value
grows with patch count and with patch cost — but on this DB it is not where the time goes. The
same A/B on a part with expensive patches (a subdivision BVH over curved patches would make each
leaf test dearer, not cheaper) is the measurement that would price it properly.

**Memory / primitive count.** 1851 analytic patches replace 48703 triangles: a 26x reduction
overall, median 184x per part (range 2x – 280x). `CloseShape` cost is comparable between the two
representations (surface 0.0001 – 0.0993 s, mesh 0.0001 – 0.0761 s per part) with no systematic
winner.

### 2026-07-26 — ALICE3, the large-part re-run

Step 6's second half: the same instrument, now on real detector geometry, run only after the small
parts had validated it.

```bash
python3 scripts/geometry/makeTestPartDB.py --models ALICE3_CAD_pure.step --output <db>
o2-bench-detectorsbase-solid-harness --db <db> --points 2000 --rays 2000 --seed 1 \
    --loop-crosscheck --pruning-ab --json step6_alice3.json
```

`ALICE3_CAD_pure.step`: 55 volumes, 36 exact-eligible, **20 paired** into the DB, 19 of which load
(one `surfaces_*.bin` is rejected by `LoadSurfaceSolid` — a trim wire whose edges do not join, a
converter-side defect, not a harness one). 2656 analytic patches against 934096 triangles, i.e. an
order of magnitude more geometry per part than the three-model DB, including four parts of
~200k triangles each.

**Correctness.** Cleaner than on the small DB, and on the two axes that matter most:

| Query | `BVH == _Loop` | vs mesh: agree / band / unexplained |
| --- | --- | --- |
| `DistFromOutside` | 38000 / 38000, deviation 0 | 65.34% / 13172 / **0** |
| `DistFromInside` | 19000 / 19000, deviation 0 | 46.19% / 10224 / **0** |
| `Contains` | 95000 / 95000 | 99.60% / 375 / **1** |

Two things to read off this. The distance twins are again bit-identical, now over parts with up to
965 patches — the traversal and the pruning are exercised far harder here than by any unit-test
fixture. And **every single disagreement with the tessellated reference except one falls inside the
mesh band**: on geometry this size the "the mesh is the approximation, not the exact solid" claim
stops being an argument and becomes a measurement. `Contains` also agrees with `Contains_Loop`
everywhere here, which confirms the 301 disagreements on the three-model DB really are specific to
its two non-manifold `oTOF` parts.

**Speed.** Per-call medians over the 19 parts (ns/call):

| Query | exact surface solid | tessellated | ratio (median, range) | BVH vs `_Loop` (median, max) |
| --- | --- | --- | --- | --- |
| `Contains` | 5282 | 518 | 12.09x (2.15 – 72.15) | — |
| `DistFromOutside` | 4715 | 559 | 8.96x (2.46 – 33.01) | **2.19x**, up to 30.27x |
| `DistFromInside` | 4492 | 667 | 10.20x (1.29 – 42.86) | **1.91x**, up to 19.28x |
| `Safety` | 44749 | 2408 | 17.68x (0.70 – 441.70) | — |

The gap to the mesh widens exactly as predicted — 9–12x here against 1.3–2.4x on small parts —
because patch count per part grows while each curved patch is still one loose conservative AABB.
That is the quantified case for the subdivision-BVH milestone, and it is the number to beat.

Conversely the BVH is worth much more here than on small parts: 2.19x over the all-surfaces loop at
the median and **30x** on the largest, where the loop's cost is linear in 965 patches. `Safety`,
which has no BVH at all, is 17.7x the mesh and 44.7 µs/call — an order of magnitude more expensive
than any other kernel, and on one part 442x the mesh. Whatever else happens next, `Safety` gets a
BVH.

**Ray `tmax` pruning A/B** (38000 rays): answers identical 38000 / 38000; leaf candidates 164981
pruned vs 227892 unpruned, i.e. **27.6% of the patch tests removed** (1.38x); wall time median
**1.02x** (0.98 – 1.25). Essentially the same picture as on the small DB — the work removed is real
and reproducible, the time saved is not yet visible. The conclusion stands: keep it (free, correct,
grows with patch cost), but it is not where the time goes today.

**Primitive count.** 2656 analytic patches replace 934096 triangles — **352x fewer primitives**
overall, median 137x per part (range 19x – 1293x). `CloseShape` is where the exact representation
wins outright at this scale: 0.0024 – 2.70 s for the surface solid against 0.0003 – **40.54 s** for
the mesh.

### 2026-07-26 (later) — navigation-reliability labelling, and the canonical-trim null result

Two changes landed after the runs above; both are documented in
[`ExactTrimTopology.md`](ExactTrimTopology.md) (items 4 and 3).

**Every part is now labelled with whether it is navigable at all.** The harness prints a
`navigation:` line per part, marks a non-navigable one next to its accuracy columns, carries the
state in `--json` under `navigation` (`reliability`, `navigable`, `boundaryEdges`,
`nonManifoldEdges`, `reversedEdges`), and repeats the list of unnavigable parts at the end of the
run. This is what the "Reading the numbers" caveat below used to have to be taken on trust:

| DB | reliable | open-surface-set | non-manifold |
| --- | --- | --- | --- |
| three-model (19 parts) | 8 | 9 | 2 |
| ALICE3 (19 loaded) | 5 | 14 | 0 |

So **11 of 19 three-model parts and 14 of 19 ALICE3 parts are not navigable.** That is the same
population the "Reading the numbers" caveat below described as "not closed manifolds", now
measured rather than asserted, and it is the honest denominator for every accuracy row above.

An `unexplained` count on a part that is not `reliable` is not a measurement of the exact solid's
accuracy — the solid's own answer is undefined in the shadow of each gap — and the output now says
so where the number is printed.

**Canonical trim-curve recognition changes no answer, by design.** Re-running both DBs against
sidecars produced by the new converter, with identical seeds and sample counts:

| DB | `contains` unexplained | `distout` | `distin` |
| --- | --- | --- | --- |
| three-model | 4588 -> 4588 | 254 -> 254 | 218 -> 218 |
| ALICE3 | 1 -> 1 | 0 -> 0 | 0 -> 0 |

Bit-identical per part, with `BVH == _Loop` unchanged and no part losing its exact conversion.
Recognition claims to change the representation and not the geometry, so this null result is the
acceptance criterion, not a disappointment. What it *did* change is how much B-spline is left for
the kernel to flatten: stored B-spline trim edges drop 88 -> 50 on the three-model DB and
**15034 -> 4528 on ALICE3**, whose sidecars shrink from 6.09 MB to 3.65 MB.

One side effect worth knowing: the three ALICE3 parts `ST1829909_002/003/004` moved from
`non-manifold` to `open-surface-set` (ALICE3 goes from 5/11/3 to 5/14/0 across
reliable/open/non-manifold). Nothing about their navigation changed — `unexplained` stayed 0 for
all three — but recognising a trim curve as a line changes how `appendDirectedEdges` chords it, so
coincidences in the closure half-edge map shift. Both states are unnavigable, so this reorders a
diagnostic rather than fixing or breaking anything. The three-model DB's labels are unchanged.

**Do not read timings from that comparison.** Those runs shared the machine with a converter run
and the per-part `ns/call` swings by up to 2x in both directions as a result. The one clean
measurement taken on an idle machine, on `BoomCylinderOuter`, is `contains` 3090 -> 3000 ns/call:
a real but small gain, consistent with only 4 of that part's trim curves being recognised.

### Reading the numbers

Two caveats travel with every row of the three-model table. First, 11 of the 19 parts are reported by `CloseShape` as
not closed manifolds (missing faces where a B-spline surface did not convert — the standing
`bspline` ceiling recorded in `BVHSurfaceSolid.md`), so their `Contains` and distance answers are
answers about an incomplete solid. The `BVH == _Loop` cross-check is unaffected by this, which is
precisely why it is the guard that carries the correctness claim. Second, the ratios have a very
wide range; the medians are the summary, not the story, and per-part numbers are in the JSON
report.
