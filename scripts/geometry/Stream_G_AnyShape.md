# Stream G — gating any `TGeoShape`, and the first measurement of the tessellated fallback

Date: 2026-08-01. Companion to [`Tutorial.md`](Tutorial.md) §5.3 and §6 item 1, which name this as
*the* blocking prerequisite for the MVP. Written by the session that built it; the integrating
session folds it.

**In one paragraph.** `runOracleGate.py` and `o2-bench-detectorsbase-solid-harness` could score
exactly one thing: an `O2BVHSurfaceSolid` loaded from a `surfaces_<part>.bin` sidecar. They can now
score **any `TGeoShape`**, and they score **every representation a part has, side by side, against
the same oracle answers**, one column per representation in `gate.json`. Two consequences arrived
with it: a CSG-emitted part is now verifiable *before* the emitter exists (proved end to end on a
`TGeoBBox` and a `TGeoCompositeShape`, both **0 disagreements on all four columns**), and the
tessellated fallback has been measured against the oracle for the first time — §4, and it is worse
than the float32 story predicted, for a different reason.

---

## 1. The convention: `shape_<part>.root`

This is the contract the CSG emitter is written against. It is stated in
`DetectorsBase/O2SolidHarness.h` next to the loader, and repeated here.

| | |
| --- | --- |
| **file** | `shape_<VOL>_<LID>.root`, in the same converter output directory as `surfaces_<VOL>_<LID>.bin`, `facets_<VOL>_<LID>.bin` and `brep_<VOL>_<LID>.brep` |
| **content** | exactly one object inheriting from `TGeoShape`, stored under the key name **`shape`** |
| **units** | centimetres |
| **frame** | the part's own **local** frame — the leaf solid as the converter emits it, **no placement matrix applied**. The same frame as the sidecar, the mesh and the `.brep`. |
| **composites** | a `TGeoCompositeShape` is written whole; its `TGeoBoolNode` and component shapes stream with it. No `TGeoManager` is needed on either side (ROOT creates a default one as a side effect of constructing a bool node; that is harmless — `TGeoShape::~TGeoShape` de-registers itself, so ownership is not shared). |

Write one with `o2::base::harness::saveShapeToRootFile(path, shape)`; read one with
`o2::base::harness::loadShapeFromRootFile(path)`. Both live in `O2SolidHarness.h/.cxx` **so that
producer and consumer cannot drift**: the unit test, the fixture generator and the harness all go
through the same pair. A file whose `shape` key is missing is searched for the first
`TGeoShape`-derived key — that exists only so a hand-made file (`root -e '...'`) still loads;
emitters must write `shape`.

**How the harness finds it.** Two ways, both live:

1. `makeTestPartDB.py` indexes `shape_<suffix>.root` into `manifest.json` as `"shape"` when the
   file exists, exactly as it already does for `"brep"`. That entry wins when present.
2. Failing that, the harness **derives** the path from the sidecar's:
   `surfaces_<suffix>.bin` → `shape_<suffix>.root` in the same directory. This is the one that
   matters in practice, because it means dropping a shape file into a finished workdir and
   re-scoring with `--skip-convert` works without re-indexing — which is the loop an emitter
   author will actually run.

Ad-hoc mode gained `--shape <file>` alongside `--surfaces` / `--facets`.

### The frame, and how it is guaranteed rather than assumed

A `TGeoShape` answers in its own local frame; the oracle answers about the `.brep` in the frame the
converter wrote it in. Get that wrong and you get a full table of plausible-looking nonsense.

So the harness **measures** it. Per representation, `bboxDeviationFromOracle` is the max deviation
in cm between that representation's own bounding box and the oracle's `bboxMin`/`bboxMax` (which
`occtOracle.py` already wrote and nothing read until now). It is reported in `gate.json` and in the
scorecard.

**Positive control, because a check that cannot fail is not a check.** Writing the `box` fixture
deliberately in the wrong frame — a `TGeoBBox` centred on the origin instead of corner-at-origin,
i.e. shifted by (−1, −1.5, −2) cm — moves the number and the columns:

| | correct frame | wrong frame |
| --- | --- | --- |
| `bboxDeviationFromOracle` | 1.0e-07 cm | **2.0 cm** |
| disagreements, all four columns | 0 | **8904** |
| `capacityRelativeDeviation` | 1.5e-16 | **1.5e-16** |

That third row is the reason the bbox check exists at all: **capacity is completely blind to a
frame error**, so a capacity-only acceptance would have passed a shape sitting 2 cm from the part
it claims to be.

Read the number with its scale in mind, and read it per representation:

| representation | typical `bboxDeviationFromOracle` | why |
| --- | --- | --- |
| `mesh` | 2.8e-07 – 9.6e-04 cm | a chorded mesh is inscribed; the box is tight to the chords |
| `shape` (analytic ROOT shape) | 1.0e-07 cm | `TGeoBBox`/`TGeoCompositeShape` compute a tight box |
| `surface` | 1e-07 cm on planar parts, but **0.13 – 16.1 cm** on parts with curved faces | every curved patch is one loose conservative AABB (the open "subdivision-BVH" milestone), and the union of those is far larger than the solid |

So the check is sharp exactly where it needs to be — on the new `shape` path, where the frame
convention is the thing that can be got wrong — and coarse on the surface representation, where
it can only catch an order-of-part-size error. That is stated rather than papered over: the
`surface` column's 16.1 cm on `Bagger/Stick` is a property of the conservative AABBs, not a frame
problem, and anyone reading this table needs to know that before drawing a conclusion from it.

---

## 2. What changed, and what deliberately did not

**Abstracted.** The scoring path is now `scoreAgainstOracle(const TGeoShape* candidate, …)` in
`runSolidHarness.cxx` — one function, called once per representation. Everything specific to
`O2BVHSurfaceSolid` (closure, rims, `NavigationReliability`, the `_Loop` twins, the BVH candidate
counters, the `tmax` pruning A/B) hangs off a `surfaceSolid` pointer that is null for every other
representation and is reported only where it means something.

**Three representations, as parallel columns and not as a flag:**

| name | source | class |
| --- | --- | --- |
| `surface` | `surfaces_<part>.bin` | `O2BVHSurfaceSolid` |
| `mesh` | `facets_<part>.bin` | `O2Tessellated` |
| `shape` | `shape_<part>.root` | any `TGeoShape` |

`gate.json` gains a `representations` array per part; each entry carries `name`, `source`,
`shapeClass`, `primitives`/`primitiveKind`, `capacityMethod`, `capacityComparable`,
`bboxDeviationFromOracle`, `closureApplicable`, `disagreements`, and the full four-column `oracle`
block. `runOracleGate.py` prints a `=== REPRESENTATION SCORECARD ===` table and per-representation
totals; `o2-bench-…-solid-harness` prints a compact one-line-per-part version.

**Not changed.** `partJson["oracle"]` is still the exact-surface representation's block, produced
by the same code, and the gate's pass/fail verdict and exit code still score that representation
only. The `representations[0]` entry duplicates it — deliberately, so the scorecard is uniform and
a consumer never has to special-case the surface.

### Two decisions where the question is meaningless for a representation

**Closure, rims and `NavigationReliability` are `O2BVHSurfaceSolid` concepts.** A
`TGeoCompositeShape` has no rims and no closure check; asking whether it is "navigable" is a
category error, and reporting either `true` or `false` would be a lie of a different kind.

Decision: **the keys are absent**, and `closureApplicable: false` says why. There is no default.
The scorecard prints `-  (not applicable to this representation)`. The mesh reports its own,
**differently named** statement, `meshClosedBody` (`O2Tessellated::IsClosedBody()`) — deliberately
not called `navigable`, because it is a property of a triangle soup decided by half-edge counting
over chords and is not the same claim. The part-level `navigation` block in `gate.json` is
untouched and still describes the surface solid.

**`TGeoShape::Capacity()` is Monte-Carlo for composites.** Verified in ROOT 6.36:
`TGeoCompositeShape::Capacity()` (`TGeoCompositeShape.cxx:282`) throws random points into the
bounding box until 10000 land inside, and it is the **only** `Capacity()` in `geom/geom/src` that
touches `gRandom`. Its relative error is ~1e-2 — four orders of magnitude above the 1e-6 gate band.

Decision: **report it, mark it, never gate on it.** `capacityMethod` is one of
`exact-divergence` (the surface solid: divergence theorem in closed form), `mesh-divergence`
(`O2Tessellated`: exact signed-tetrahedra sum over its own triangles, deterministic, and therefore
a real measurement — of the chording deficit), `root-analytic`, or `root-montecarlo`.
`capacityComparable` is false only for the last, and the scorecard prints `n/a` there.

The measured justification, from the end-to-end check below: the `box_minus_cyl`
`TGeoCompositeShape` is **geometrically exact** — 0 disagreements on all four oracle columns — yet
ROOT reports its capacity as 55.9392 cm³ against OCCT's 55.95752, a relative deviation of
**3.3e-04**. A capacity gate would have failed a shape that is exactly right, by 330×. The
project's intended acceptance for CSG parts remains OCCT symmetric-difference volume; that is a
later step and is not this one.

The small residual capacity deviations on the three Bagger parts (1.39/1.31/1.39e-06) were
explicitly out of scope and are untouched.

---

## 3. Inertness of the existing path

The rule is *diff columns, not verdicts*. `gate.json` was compared field by field with `timing*`,
`*Seconds` and `nsPerCall` stripped (`/tmp/diffgate.py`, a flatten-and-compare over every leaf),
against the baseline captured immediately before the change on the same machine.

| corpus | leaf fields removed | leaf fields **changed** | added |
| --- | --- | --- | --- |
| fixtures (9 parts) | **0** | **0** | 3953, all under `representations` |
| Bagger (12 parts) | **0** | **0** | 5707, all under `representations` |

So: **no pre-existing field moved.** The only new top-level key is `representations`.

Gate totals and disagreement counts, before → after, quoted together as the contract requires:

| | before | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 75 cases, green | **78 cases, green** (3 appended, §5) |
| fixtures gate | 9/9 | **9/9** |
| Bagger gate | 9/12 (capacity 1.39/1.31/1.39e-06) | **9/12, the same three, the same values** |
| unexplained disagreements, surface, fixtures | contains 0, distout 0, distin 0, safety 0 | **0, 0, 0, 0** |
| unexplained disagreements, surface, Bagger | contains 0, distout 0, distin 0, safety 0 | **0, 0, 0, 0** |

One build-system change: `ROOT::RIO` was added explicitly to the `solid-harness` target in
`Detectors/Base/CMakeLists.txt`. It turned out to be a no-op in effect — CMake reported
`ninja: no work to do` and `ldd` shows `libRIO` was already linked transitively through
`O2::DetectorsBase` — but the dependency is now declared where it is used. **The documented
incremental-`ninja`/rootcling failure did not recur**; the reconfigure was clean.

---

## 4. The prize: how wrong the tessellated fallback actually is

Nobody had ever measured this. The old acceptance criterion compared *against* a mesh, which is
exactly why it could not certify anything; the mesh had therefore never been on the *candidate*
side of the oracle. Now it has, on all 21 parts of both corpora, with the same sample sets and the
same band (the model's declared tolerance, 1.0e-07 cm on both corpora).

Each cell is `disagreements / scored (worst finite deviation, cm)`. The parenthesised number is the
worst deviation among the recorded offenders **excluding total misses** — a ray where one side
finds a crossing and the other finds none has deviation 1e30 and is counted separately, in the
`missedSurface` column of the totals table below. A blank means every recorded offender for that
cell was a total miss.

| part | triangles | contains | distout | distin | safety | \|ΔV\|/V |
| --- | ---: | --- | --- | --- | --- | ---: |
| `box` | 12 | **0**/4999 | **0**/2000 | **0**/1000 | **0**/3000 | 1.5e-16 |
| `box_union_box` | 20 | **0**/5000 | **0**/2000 | **0**/1000 | **0**/3000 | 0 |
| `box_minus_cyl` | 520 | 10/5000 (2.3e-04) | 411/2000 (1.7e-03) | 159/1000 (3.6e-03) | 569/3000 (2.5e-04) | 6.0e-05 |
| `cyl_inter_cyl` | 604 | 36/4999 (2.5e-04) | 1200/2000 (4.5e-03) | 999/1000 (3.3e-03) | 699/3000 (5.4e-04) | 3.7e-04 |
| `cyl_cross_cyl` | 1612 | 16/5000 (4.7e-04) | 1056/2000 (2.7e-03) | 855/1000 (3.8e-03) | 772/3000 (5.0e-04) | 3.0e-04 |
| `tube_window` | 1248 | 18/5000 (2.9e-04) | 1072/2000 (1.9e-02) | 818/1000 (4.9e-03) | 856/3000 (4.7e-04) | 3.4e-04 |
| `cyl_plus_cone` | 1760 | 47/5000 (1.1e-03) | 1059/2000 (1.1e-02) | 895/1000 (4.9e-03) | 653/3000 (1.1e-03) | 7.1e-04 |
| `sphere_minus_cyl` | 7048 | 127/5000 (1.8e-03) | 1136/2000 (4.5e-02) | 1000/1000 (2.1e-02) | 1141/3000 (1.9e-03) | 1.9e-03 |
| `torus_union_cyl` | 23432 | 29/5000 (1.1e-03) | 1002/2000 (1.8e-02) | 778/1000 (5.4e-03) | 871/3000 (1.3e-03) | 5.5e-04 |
| `Bagger/BasePin` | 500 | 19/5000 (3.1e-04) | 1011/2000 (2.5e-02) | 922/1000 (6.4e-03) | 573/3000 (3.1e-04) | 4.1e-04 |
| `Bagger/Base` | 4364 | 3/4999 (0) | 136/2000 (3.3e-02) | 100/1000 (1.7e-02) | 232/3000 (7.5e-04) | 2.8e-05 |
| `Bagger/Boom` | 4520 | 2/5000 (1.3e-04) | 92/2000 (1.9e-03) | 67/1000 (1.8e-03) | 158/3000 (1.1e-03) | 1.1e-05 |
| `Bagger/Stick` | 4412 | 1/5000 (0) | 452/2000 (1.2e-03) | 232/1000 (1.7e-03) | 221/3000 (4.2e-04) | 2.0e-06 |
| `Bagger/BoomCylinderInner` | 1612 | 7/4999 (1.3e-04) | 651/2000 (6.7e-03) | 950/1000 (2.1e-03) | 1101/3000 (7.6e-04) | 2.8e-04 |
| `Bagger/BoomCylinderOuter` | 2244 | 5/5000 (2.5e-04) | 708/2000 (7.0e-03) | 928/1000 (6.1e-03) | 1491/3000 (4.6e-04) | 2.7e-04 |
| `Bagger/StickCylinderInner` | 1612 | 3/5000 (4.8e-05) | 657/2000 (2.1e-03) | 943/1000 (2.4e-03) | 1076/3000 (3.7e-04) | 2.8e-04 |
| `Bagger/StickCylinderOuter` | 2244 | 5/4999 (1.6e-05) | 679/2000 (7.0e-03) | 933/1000 (4.2e-03) | 1510/3000 (4.7e-04) | 2.7e-04 |
| `Bagger/BucketCylinderInner` | 1636 | 1/5000 (0) | 552/2000 (8.9e-03) | 936/1000 (4.8e-03) | 1148/3000 (2.8e-04) | 2.8e-04 |
| `Bagger/BucketCylinderOuter` | 2178 | 2/5000 (1.4e-04) | 844/2000 (1.3e-02) | 922/1000 (5.4e-03) | 1393/3000 (5.9e-04) | 3.5e-04 |
| `Bagger/BucketLink1` | 4316 | 4/5000 (1.3e-04) | 479/2000 (3.3e-03) | 500/1000 (2.4e-03) | 619/3000 (5.3e-04) | 1.4e-04 |
| `Bagger/BucketLink2` | 4628 | **366**/5000 (**2.2**) | 460/2000 | 540/1000 | 777/3000 (**0.70**) | **9.0e-02** |

**Totals over all 21 parts:**

| column | disagreements | scored | rate | of which "missed a surface entirely" |
| --- | ---: | ---: | ---: | ---: |
| `contains` | 701 | 104995 | 0.67 % | 0 |
| `distout` | 13657 | 42000 | **32.5 %** | 152 |
| `distin` | 13477 | 21000 | **64.2 %** | 71 |
| `safety` | 15860 | 63000 | **25.2 %** | 0 |

Compare against the surface solid on the same rows: **0, 0, 0, 0** on every column of every part.

### What the numbers say

1. **A mesh of an all-planar part is exact.** `box` and `box_union_box` score **0 on all four
   columns** — the only two parts in either corpus with no curved face. The tessellated fallback
   is not approximate *per se*; it is approximate *where there is curvature*. That is the single
   most useful fact here for the `auto`-mode fallback policy: the policy can be conditioned on
   whether the part has any non-planar face, not on a global "meshes are bad".

2. **A dead premise, reported as a result.** The brief states — correctly — that facets are stored
   float32, half-ulp 1.5e-05 cm at 400 cm, and that this is a floor on what the mesh can achieve.
   It is a floor, and it is **irrelevant**: the measured worst deviations are 1e-03 to 4.5e-02 cm,
   with the Bagger parts sitting at ~60 cm coordinates where the float32 half-ulp is ~3.6e-06 cm.
   The mesh is **10³–10⁴ times worse than its own storage precision**. Chording is the entire
   error budget; converting `facets_*.bin` to float64 would buy nothing measurable. Any effort
   aimed at mesh precision has to go into the deflection parameters, not the number format.

3. **`distin` is the worst column, by construction, and the ordering is stable.** The mesh is
   inscribed, so an inside ray exits the exact solid later than it exits the chording, on
   essentially every curved patch — 64 % disagreement against 32 % for `distout`. This reproduces
   the long-standing qualitative claim in `SolidNavigationHarness.md` with the *oracle* rather than
   with the mesh on the reference side, which is the first time it has been more than an argument.

4. **`contains` is nearly right and that is misleading.** 0.67 % overall; the mesh's error is
   sub-millimetre and a random point is rarely that close to a wall. A fallback policy set on
   `contains` alone would look ten times better than the thing it is choosing.

5. **One part is in a different class: `Bagger/BucketLink2`.** 366 `contains` disagreements with a
   worst deviation of **2.2 cm**, safety off by 0.70 cm, volume off by 9.0e-02 relative — a
   thousand times the next-worst part on every column. It is also **the only part in either corpus
   whose mesh is not a closed body** (`meshClosedBody=false`, 20/21). This is the same part whose
   non-watertight mesh is already documented in `O2SolidHarness.h` (finding H1: its left plate sits
   0.2 cm from the BREP's, which is why ray re-labelling by the oracle's origin classification had
   to be introduced). The oracle now measures that defect directly instead of it being inferred
   from offender listings. **A `meshClosedBody=false` part must not be allowed to ship as a
   tessellated fallback**; that is a concrete, measured input to Stream E step 3.

6. **152 + 71 total misses.** In 223 rays the mesh and the oracle disagree about whether there is
   any crossing at all. Those are the failure mode that matters for navigation (a lost track, not
   a slightly wrong step), and they are all on curved geometry.

---

## 5. The synthetic `TGeoShape` end-to-end check

Two ladder fixtures are, by construction, *exactly* a ROOT shape. `runOracleGate.py --fixture-shapes`
writes them in the sidecar convention and the ordinary gate then scores them:

```bash
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py \
    --workdir /tmp/gate --fixtures --fixture-shapes
```

| fixture | emitted as | contains | distout | distin | safety | bboxDev | capacity |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| `box` | `TGeoBBox(1, 1.5, 2)` at origin (1, 1.5, 2) | **0** | **0** | **0** | **0** | 1.0e-07 | 1.5e-16 |
| `box_minus_cyl` | `TGeoCompositeShape` = `TGeoBBox(2,2,2)` − `TGeoTube(0, 0.8, 2.5)` | **0** | **0** | **0** | **0** | 1.0e-07 | *n/a (MC)* |

**The any-shape path works end to end — convention, loader, frame, scoring — and a
`TGeoCompositeShape` built as ROOT booleans of primitives is oracle-clean on all four navigation
columns.** That is the thing the MVP was blocked on, demonstrated before the emitter exists.

`box` is deliberately the *offset* case: `TGeoBBox` is centred on `fOrigin` and this fixture is
corner-at-origin, so it is the cheapest way for a frame convention to be silently wrong, and it is
the one the negative control in §1 was run on.

Three unit cases were appended to `Detectors/Base/test/testBVHSurfaceSolid.cxx` in a delimited
block ending `// --- Stream G: gating any TGeoShape ---`:

- `ShapeSidecarRoundTripsAnyTGeoShape` — an offset `TGeoBBox` and a `TGeoCompositeShape` through
  `saveShapeToRootFile` → `loadShapeFromRootFile`, requiring bit-identical `Contains`, `Safety`,
  `DistFromOutside` and `DistFromInside` over 9 probe points × 5 directions;
- `ShapeSidecarRefusesWhatIsNotAShape` — a missing file, and a well-formed ROOT file whose `shape`
  key holds a `TNamed`. Both must be refused with a message, so an emitter that writes the wrong
  object does not look like an emitter that wrote nothing;
- `OracleValidatorsScoreAPlainRootShape` — the abstraction's own positive **and negative** control.
  A `TGeoTube` is scored against analytically-constructed oracle columns (`contains` and `safety`
  from closed form, independent of both shapes) and must show 0 disagreements; a second tube with
  its radius 0.05 cm too large must show more than 0, on every column. Without the negative half,
  a validator that is structurally unable to report a disagreement would pass.

---

## 6. Reproducing everything here

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
ctest -R BVHSurfaceSolid                       # 78 cases

cd $HOME/alisw/O2
# the ladder, including the two hand-built TGeoShape sidecars
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures --fixture-shapes
# Bagger
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
```

Both print `=== REPRESENTATION SCORECARD ===` after the gate summary, and the per-column
disagreement totals for the surface representation immediately after the gate total — because
those two numbers are never to be quoted apart.

To score a shape of your own against a part that is already in a workdir, with no reconversion:

```bash
# write <workdir>/db/<slug>/shape_<VOL>_<LID>.root, then
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir <workdir> --skip-convert
```

---

## 7. What this leaves open

- **`nsPerCall` for the `shape` representation is not reported.** The timing block still times the
  surface against the mesh only. Adding it is mechanical; nobody has needed it yet, and the CSG
  timing question (`CSG_Pipeline.md`'s composite-depth probe) wants a different experiment anyway.
- **The scorecard is a report, not a verdict.** Nothing in the exit code depends on the mesh or
  shape columns, deliberately: acceptance for a CSG part is OCCT symmetric-difference volume plus
  an oracle-clean gate in its emitted form, and only the second half exists today.
- **`--fixture-shapes` knows two fixtures.** It is a smoke test, not a recogniser. Extending it is
  Stream A's job and its output will arrive through `makeTestPartDB.py`'s `"shape"` index rather
  than through this flag.
- **`shape` sidecars are not produced by anything in the converter.** That is the emitter, and it
  is the next item.
