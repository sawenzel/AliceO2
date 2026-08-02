# Stream N — placed primitives: the self-union composite is gone

> **Name collision, since resolved.** This document and `Stream_O_ImplicitTrims.md` were written
> concurrently by two sessions and both were originally lettered **N** — a briefing mistake, not a
> disagreement. The trim census was renamed to **O** once both had landed. Cite streams by file
> name, never by the letter; the code comments in this change all name the file.

Date: 2026-08-02. Removes the *degenerate* `TGeoCompositeShape` — a recognised primitive unioned
with an identical copy of itself under the same matrix — that every placed CSG part used to be
written as. Genuine multi-leaf booleans are untouched and still emit a composite.

**In one paragraph.** A recognised primitive is now emitted **in its own canonical frame at the
origin** and the rigid transform travels beside it, as a `TGeoHMatrix` under the key `placement` in
`shape_<VOL>_<LID>.root` and as a 3×4 array in `csg_<part>.json`, `csg_report.json`,
`manifest.json` and `gate.json`. **An artefact that records no placement means the identity**, so
every file and every `db` directory written before this change still loads and still scores exactly
as it did. Consumers compose it where they used to rely on the shape already being in the part
frame: the harness and the X-ray benchmark transform points and rays into the shape's frame, and
`geom.C` places the volume with `partPlacement * shapePlacement`. On `Bagger.step` the one affected
part, `BasePin`, goes from `TGeoCompositeShape` with a Monte-Carlo `Capacity()` **2.045e-04** off
the oracle and `capacityComparable=false`, to `TGeoTube` with an analytic `Capacity()`
**6.785e-16** off and `capacityComparable=true`. Nothing else moved: `ctest` 92 → **97** cases
green, both gates keep **0/0/0/0** unexplained disagreements on `surface` and on `shape`, and a
field-by-field `gate.json` diff accounts for every one of its 115 moved fields.

---

## 1. Why this was worth doing, and what it cost before

`csg/primitives.py` used to state the constraint as permanent:

> **No ROOT shape class can carry a rigid transform.** … So a recognised tube that is not already
> on the z axis of the part's own frame can only be written as a composite.

The premise is still true of ROOT 6.36 — `TGeoBBox` has `fOrigin` and nothing else has anything;
the only `TGeoShape` holding a `TGeoMatrix` is `TGeoCompositeShape`, through its `TGeoBoolNode`,
which needs two operands. The **conclusion** was wrong, because it assumed the transform had to
live inside the artefact's single object. It does not have to live inside the shape at all.

What the degenerate composite cost, measured rather than asserted:

| cost | before | after |
| --- | --- | --- |
| `Capacity()` | Monte-Carlo (`TGeoCompositeShape` throws 10000 points into its bbox) — `Bagger/BasePin` **−2.045e-04** relative, and **different on every call** | closed form — **−6.785e-16** relative |
| the gate's capacity column | `capacityComparable=false`; the whole criterion unavailable for that part | `capacityComparable=true` |
| what `gate.json` says the part *is* | `TGeoCompositeShape` for a recognised plain tube | `TGeoTube` |
| navigation | a boolean tree: two `Contains` and two matrix applications per query | the primitive's own kernel, one matrix application by the caller |
| objects per part | 3 shapes + 2 matrices + 1 bool node | 1 shape + 1 matrix |

The user's instruction was exactly this: *"Maybe you could move the thing to the origin, save it,
and remember this transformation when finally placing the part in the assembly."*

---

## 2. The artefact format

One description, two spellings of the same number, because no single interpreter can read both
(the ROOT file is unreadable from an OCC-only Python, the JSON is unreadable from the C++ harness).
Both are produced by **one** function, `csg.primitives.placement_for_candidate()`, which needs no
ROOT and is therefore also what the deferred `emit.py --from-json` path uses.

```
placement := 3x4 row-major [R | t],  part = R * canonical + t
             R's COLUMNS are the leaf frame's basis vectors in the part frame
             null / absent  ==  identity
```

| artefact | where the placement lives | written by |
| --- | --- | --- |
| `shape_<VOL>_<LID>.root` | a `TGeoHMatrix` under the key **`placement`**, beside the `TGeoShape` under `shape` | `csg/emit.py`, `harness::saveShapeToRootFile` |
| `csg_<VOL>_<LID>.json` | `"placement"`: 3×4 or `null` | `csg/hook.py` |
| `csg_report.json` | `parts[].shapePlacement` | `csg/hook.py` |
| `manifest.json` | `parts[].shapePlacement` (mirrored from `csg_report.json`) | `makeTestPartDB.py` |
| `gate.json` | `representations[].placement`: 3×4 or `null`, for every representation | `runSolidHarness.cxx` |

Three properties were deliberate:

1. **Absence means identity, and the identity is never written.** `saveShapeToRootFile` refuses to
   store a matrix for which `IsIdentity()` holds, so "no key" stays the one and only spelling of
   the identity and a reader has one case fewer to get wrong. Unit test:
   `AbsentPlacementMeansIdentity`.
2. **An axis-aligned box keeps `fOrigin`.** `TGeoBBox` is the one ROOT primitive that carries its
   own translation, it already had an analytic capacity, and leaving it alone keeps every artefact
   for a box byte-identical to before. `placement_for_candidate()` returns `None` there.
3. **A genuine multi-leaf boolean has no placement.** The composite is already expressed in the
   part frame and its `TGeoBoolNode` carries the leaves' matrices itself. Unit-asserted in
   `csg/emit.py --self-test` ("a genuine two-leaf union is still an unplaced
   `TGeoCompositeShape`"), so "the composite is gone" can never quietly become true of the cases
   that legitimately need one.

**`build_occ()` is unchanged.** It still builds the solid in the *part* frame, so the OCCT
symmetric-difference acceptance keeps measuring the **placed** solid against the original CAD
solid. That is the one thing that could have silently changed meaning, and the reason it did not is
that the change is confined to `build_root()`. `crosscheck_contains()` — the sharp ROOT-vs-CAD check
— now composes the placement before asking the shape, which is what makes a wrong placement a
disagreement against the CAD body rather than a silent pass.

---

## 3. The composition order, and how it was proved

```
node matrix in the assembly  =  partPlacement * shapePlacement
```

A point travels shape → part → parent, so the part placement is on the left. In `geom.C`:

```cpp
TGeoHMatrix *tr_..._placed = new TGeoHMatrix(*tr_...);   // partPlacement
tr_..._placed->Multiply(shapePlace__0_1_1_2);            // * shapePlacement
asm__0_1_1_1->AddNode(vol__0_1_1_2, 1, tr_..._placed);
```

`TGeoHMatrix::Multiply(right)` computes `this = this * right`, so the part placement is the one
that gets copied and the shape placement is the right operand. The composition is emitted **per
`AddNode`**, so a prototype placed *n* times gets *n* correctly composed matrices rather than one
shared, mutated one. A CSG volume that is itself a *top* volume has no `AddNode`, so it is wrapped
in a one-node `TGeoVolumeAssembly` — the only place the matrix could otherwise have gone.

It is proved three ways, and every one of them carries a control that fails.

### 3.1 In the unit test, by navigating (`NodeMatrixIsPartPlacementTimesShapePlacement`)

A `TGeoManager` per candidate node matrix; the reference is built **without ever forming the
product** — the point is carried into the part frame by the part placement, then into the shape's
frame by the shape placement, and the tube membership is evaluated there. Over 24 389 lattice
points (of which 296 inside):

| node matrix | disagreements |
| --- | ---: |
| `partPlacement * shapePlacement` | **0** |
| `shapePlacement * partPlacement` (reversed) | > 0 |
| `partPlacement * shapePlacementᵀ` (transposed rotation) | > 0 |
| `partPlacement` alone (placement dropped) | > 0 |

### 3.2 In `O2_CADtoTGeo.py --self-test`, against OpenCascade

Same idea, but the reference is `BRepClass3d_SolidClassifier` on the CAD solid, reached by undoing
the *part* placement only — so the reference has never heard of the shape placement. 3000 probes,
163 of them inside the CAD body:

```
  [ok ] geom.C's node matrix partPlacement * shapePlacement puts the solid where the CAD body is
        -- 0 disagreement(s) over 3000 navigated points
  [ok ] the reversed product does move the count            -- 163 disagreement(s)
  [ok ] a transposed shape rotation does move the count     -- 229 disagreement(s)
  [ok ] dropping the shape placement does move the count    -- 163 disagreement(s)
```

### 3.3 End to end, on the real emitted macro

The strongest of the three, because it tests the macro the converter actually wrote rather than a
restatement of the rule. `geom.C` from a **pre**-change Bagger conversion and from a **post**-change
one were both built in ROOT, and 200 000 points of a *fixed* box were classified by
`TGeoManager::FindNode`:

| pair | points in a volume | disagreements |
| --- | ---: | ---: |
| pre-change `geom.C` vs post-change `geom.C` | 9477 / 9477 | **0** |
| post-change vs the same macro with the `Multiply` lines deleted | 9477 / 9450 | **165**, every one of them `BasePin` |

The control is what makes the zero mean something, and it also names the affected part: 165
disagreements, all `('BasePin', '')` or `('', 'BasePin')` or `('Boom', 'BasePin')`, i.e. exactly the
one Bagger part with a non-identity placement.

**A methodological trap, recorded because it cost a wrong conclusion first.** The first version of
that probe derived its sampling box from the top volume's own bounding box. A misplaced part
changes a `TGeoVolumeAssembly`'s bbox, so the control run sampled *different points* and every
volume in the model appeared to move — 7253 disagreements spread over parts that had not changed at
all. A control is only a control if the questions are identical. The box is now passed in.

---

## 4. The consumers, and the one design choice worth arguing about

### 4.1 The harness: transform the queries, do not wrap the shape

`runSolidHarness.cxx` transforms the sample points and rays into the shape's own frame
(`MasterToLocal` / `MasterToLocalVect`) and asks the shape directly. The alternative — wrapping the
shape in something that carries the matrix, or querying it through a placed `TGeoVolume` — was
rejected, and the reason is not taste:

* **the object the gate scores must be the object the converter emitted.** With a wrapper,
  `shapeClass` in `gate.json` is the wrapper's class, `Capacity()` is the wrapper's, and the
  wrapper's bounding box is the inflated corner hull again — i.e. every quantity this change was
  made to recover would be lost to the instrument that is supposed to measure it;
* a rigid transform preserves the inside/outside relation **and every distance along a ray**, so
  the oracle's answers need no correction at all — only the questions move;
* it costs nothing when there is no placement: `toLocal()` returns an empty vector and the
  unplaced case selects the caller's own vector by reference.

`bboxDeviationFromOracle` carries the shape's box into the part frame by transforming its eight
corners. For a rotated body that hull is strictly larger than the body, so the number is
conservative — **exactly as conservative as it already was for a composite**, whose
`TGeoBoolNode::ComputeBBox` does the same thing internally. Measured: unchanged at 1.000e-07 cm for
`BasePin`, and unchanged for all six rams (0.13 – 0.62 cm). It remains a *frame* check, not a
tightness measurement.

### 4.2 The X-ray benchmark: one of each, on purpose

Mode **(a)** (bare shape API) transforms each ray by hand. Mode **(b)** (the real `TGeoNavigator`)
puts the matrix on the **node** — `mWorld->AddNode(mPart, 1, placement)` — and leaves the rays in
the part frame, which is what a navigator is for and what production does. The two are therefore
independent statements about the same transform rather than one statement made twice, and
`modeAvsB` is a real cross-check of the placement. The raster window and the navigator world are
both computed from the *placed* box (`placedBox()`); rastering a rotated tube against the box of
the tube at the origin would be a window that misses the part.

---

## 5. The census the brief asked for: how common is a single placed primitive?

Counted from `csg_<part>.json` across every corpus the converter can process, with `--csg auto`.
"Recognised" here means *accepted by both acceptance tests* — the parts that actually become a
`shape_*.root`.

| corpus | leaf solids | accepted as CSG | **single primitive, unplaced** | **single primitive, placed** | **genuine multi-leaf** | rejected by dV_sym | declined |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| boolean fixtures (10 models) | 10 | 2 | 1 (`box`, `TGeoBBox`) | 0 | 1 | 2 | 6 |
| `Bagger.step` | 13 | 7 | 0 | **1** (`BasePin`, `TGeoTube`) | 6 | 0 | 6 |
| ALICE3 `CAD_noETA.stp` | 55 | 2 | 0 | **2** (`TGeoTubeSeg`) | 0 | 0 | 53 |
| `as1-oc-214.stp` | 5 | 0 | 0 | 0 | 0 | 0 | 5 |
| `oTOF System V3-R92cm.step` | 3 | 0 | 0 | 0 | 0 | 0 | 3 |
| **total** | **86** | **11** | **1** | **3** | **7** | **2** | **73** |

Multi-leaf detail — the distribution the brief asked for:

| leaf count | parts | operation | leaf types |
| ---: | ---: | --- | --- |
| 2 | 7 | `union` (all of them) | `TGeoTube` × 14 |
| ≥3 | 0 | — | — |

So: **4 of the 11 accepted parts (36%) are single primitives**, 3 of those 4 are placed and become
plain `TGeoTube`/`TGeoTubeSeg` under this change; **7 (64%) are genuine two-leaf unions** and keep
their composite. There is no part anywhere in these corpora with more than two leaves, and the only
operation the recogniser ever emits is `union` — `op` is declared as `primitive | union` and
`_root_composite` raises on anything else.

**Three things in that table contradict what the documents currently say, and they are the
findings, not the counts.**

1. **`Stream_A_CSG.md`'s "62560 placed six-plane boxes in oTOF" is not reachable through the
   converter.** `O2_CADtoTGeo.py` sees **3** leaf solids in `oTOF System V3-R92cm.step`
   (`oTOF_v2_v1`, `A_Side`, `Module_v1`), not 20 prototypes, and recognises none of them. Without
   `--mesh` it does not even complete: `triangulate_asbbox()` dies with
   `Standard_ConstructionError: Bnd_Box is void`. The census in `csg/census.py` walks the STEP
   assembly differently from the converter's own leaf-solid extraction, and only the converter's
   view can produce artefacts. **This is a pre-existing defect, unrelated to Stream N**, but it
   means the corpus with by far the most placed boxes contributes nothing today.
2. **ALICE3 yields 2 CSG parts, not the "36 eligible solids" figure.** Those two are
   `ST1A38494_002` and `_003`, both `tier1-tubeseg`, both placed. The 36/55 number in `Tutorial.md`
   §5.2b is about *surface-recognition eligibility*, which is a different question from whether the
   whole solid matches a CSG template; the two have been quoted near each other and are easy to
   confuse.
3. **`as1-oc-214` converts 0 of 5 as CSG.** `Tutorial.md` §5.2 says it goes "0/5 before Tier 0 and
   **5/5** after" — that is about the *exact-surface* path, not this one.

The honest summary of the census is therefore: **the change is unambiguously correct and unlocks
the capacity gate for every single-primitive part, but on today's reachable corpora it applies to
three parts.** Its value is in what it removes as a permanent constraint, and in the corpora
(oTOF-like assemblies of placed boxes and tubes) that the converter cannot yet reach.

---

## 6. Acceptance — every number run for this change

### 6.1 The four self-tests

| | result |
| --- | --- |
| `ctest -R BVHSurfaceSolid` | **Passed**, `Running 97 test cases... *** No errors detected` (was 92) |
| `python3 scripts/geometry/O2_CADtoTGeo.py --self-test` | **26 checks, 0 failures** (18 recognition + 8 placement; was 18) |
| `O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --self-test` | **17/17 verdict self-checks passed** |
| `python3 scripts/geometry/csg/emit.py --self-test` | **33/33** (11 acceptance + 22 recognise/emit) |

The five new `ctest` cases are `ShapeSidecarRoundTripsAPlacement`, `AbsentPlacementMeansIdentity`,
`PlacedPrimitiveAnswersExactlyLikeTheSelfUnionComposite`,
`PlacedPrimitiveRecoversTheAnalyticCapacity`, `NodeMatrixIsPartPlacementTimesShapePlacement`.

The tolerance the brief asked to be stated: the placed tube's `Capacity()` agrees with the closed
form to **3.02e-16** relative, and the placed tube *segment*'s agrees with OCCT's
`BRepGProp::VolumeProperties` on the CAD solid to **0.00e+00** — bit-identical, `19.6349540849` both
ways. The negative control (radius 1% too large) reports **2.68e-02**.

### 6.2 The gates — totals *and* disagreement counts, together

**Fixtures** (10 models, 9 scored parts), before and after, identical:

```
surface  contains=0  distout=0  distin=0  safety=0   (9/9 part(s) with zero disagreements)
mesh     contains=283  distout=6936  distin=5504  safety=5561   (2/9 part(s) with zero disagreements)
shape    contains=0  distout=0  distin=0  safety=0   (2/2 part(s) with zero disagreements)
```
9/9 parts PASS on the representation they ship in.

**`Bagger.step`** (13 leaf solids, 12 scored), before and after, identical:

```
surface  contains=0  distout=0  distin=0  safety=0   (12/12 part(s) with zero disagreements)
mesh     contains=418  distout=6721  distin=7973  safety=10299   (0/12 part(s) with zero disagreements)
shape    contains=0  distout=0  distin=0  safety=0   (7/7 part(s) with zero disagreements)
```
12/12 scored parts pass on the representation they ship in; 9/12 on the surface representation
(the historical number, unchanged); 1 leaf solid (`Bucket`) unscored, as before.

**Unexplained oracle disagreements are 0/0/0/0 on `surface` and 0/0/0/0 on `shape`, on both
corpora, before and after.**

### 6.3 The win, before and after, per part

`gate.json`, `representations[shape]`:

| part | before | after |
| --- | --- | --- |
| `Bagger/BasePin` | `TGeoCompositeShape`, `root-montecarlo`, comparable **false**, relDev **−2.0452e-04** | `TGeoTube`, `root-analytic`, comparable **true**, relDev **−6.7852e-16** |
| the six Bagger rams | `TGeoCompositeShape`, `root-montecarlo`, comparable false, relDev −2.4e-02 … +6.0e-03 | unchanged in kind; relDev redrawn (−2.1e-02 … +1.9e-02) because the Monte-Carlo is redrawn |
| `fixtures/box` | `TGeoBBox`, `root-analytic`, comparable true, relDev +1.4803e-16 | **byte-identical** |
| `fixtures/cyl_cross_cyl` | `TGeoCompositeShape`, comparable false | unchanged in kind |

`bboxDeviationFromOracle` is **unchanged for every part on every representation**, including
`BasePin`'s 1.0000e-07 cm — the frame check did not loosen.

### 6.4 `compareGateRuns.py`, Bagger pre vs post: 115 differences, all accounted for

| moved field | count | why |
| --- | ---: | --- |
| `representations.*.source`, `verdict.shipped*.source` | 43 | the two runs live in different workdirs; absolute paths |
| `representations.*.placement` (null) | 30 | the **new column**, `null` for every unplaced representation |
| `representations.2.oracle.capacityCandidate` / `capacityRelativeDeviation` | 14 | 1 × `BasePin` (Monte-Carlo → analytic, the win) + 6 × ram composites (Monte-Carlo redrawn) |
| `representations.2.placement.i.j` | 12 | `BasePin`'s placement matrix elements |
| `representations.2.shapeClass`, `verdict.shippedVerdict.shapeClass` | 2 | `BasePin`: `TGeoCompositeShape` → `TGeoTube` |
| `representations.2.capacityMethod`, `capacityComparable` | 2 | `BasePin` |

**Nothing else moved.** Not one oracle disagreement count, not one `bboxDeviationFromOracle`, not
one navigation field, not one verdict, on any part.

`BasePin`'s recorded placement, read straight out of the diff, is
`R = [[0,1,0],[1,0,0],[0,0,−1]]`, `t = (0, 5.916, 5)` — a proper rotation (det +1) taking local z to
`(0, 0, −1)`, which is precisely the axis `Stream_H_CSGEmitter.md` §2.2 states for that part.

### 6.5 X-ray transport

See §7 for the numbers; run with `--beams 96` on the fixtures corpus and additionally on
`Bagger.step`, because **the fixtures corpus contains no placed primitive** (its only single
primitive is the axis-aligned `box`, which keeps `fOrigin` and has no placement). A fixtures-only
X-ray run would have been a negative control with nothing to control — the same failure mode
`NEXT.md` warns about under *"a refuted experiment is not a refuted hypothesis"*.

---

## 7. X-ray transport results

`runXRayBench.py --beams 96` — 96 Fibonacci directions × 48² rays = 221 184 rays per part, both
stepping modes, per representation. Not an axis raster: a parallel beam is direction-poor and three
axis beams are three directions however many rays are fired (`NEXT.md`).

### 7.1 Fixtures (9 parts) — the no-regression run

```
    shape    (a) shape  steps=1083028  zeroStep=0 noAdv=0 unstick=0 capHit=0 unterm=0 oddList=0
                        nonAlt=0 dupXing=0 parity=0 parityNB=0 noTrans=0 outWorld=0 orgIn=0
    shape    (b) nav    steps=1083028  (every counter 0, as above)
             mode (a) vs mode (b): 0 of 442368 rays disagree
    shape    vs OCCT:   442368/442368 rays identical, LOST=0 extra=0 displaced=0 kind=0
                        worst dt=2.416e-13 cm   (2/2 part(s) fully clean)

    surface  (a) shape  steps=4916790  unterm=2 oddList=2 parity=2, every other counter 0
    surface  (b) nav    steps=4916790  (identical)
             mode (a) vs mode (b): 0 of 1990656 rays disagree
    surface  vs OCCT:   1990651/1990656 rays identical, LOST=0 extra=2 displaced=3 kind=0
                        worst dt=8.815e-06 cm   (8/9 part(s) fully clean)

    mesh     vs OCCT:   1193451/1990656 rays identical, LOST=552 extra=92 displaced=1545875
```

The two `surface` parity mismatches are the **pre-existing** ×1 residual `Tutorial.md` §5.7 already
records ("the parity audit … fires even at the shipped ×1 size"). They cannot have been caused by
this change: only the `shape` representation can carry a placement, `toShapeFrame(nullptr, …)` is
an exact copy, and `placedBox(box, nullptr, …)` reproduces the previous expression term for term.

**And this run exercises nothing new**, which is why it is not the whole story: the fixtures corpus
contains no placed primitive. Its only single primitive is the axis-aligned `box`, which keeps
`fOrigin` and has no placement; `cyl_cross_cyl` is a genuine union. The fixtures run is a
*no-regression* check, and quoting it as evidence that the placement works would be exactly the
mistake `NEXT.md` warns about under *"a refuted experiment is not a refuted hypothesis"*. The
Bagger run below is the one that transports through a placed shape.

### 7.2 `Bagger.step` (12 parts) — the run that exercises the placement

`BasePin` is the placed `TGeoTube`; the six rams are unplaced composites.

```
    shape    (a) shape  steps=2564290  zeroStep=0 noAdv=0 unstick=0 capHit=0 unterm=0 oddList=0
                        nonAlt=0 dupXing=0 parity=0 parityNB=0 noTrans=0 outWorld=0 orgIn=0
    shape    (b) nav    steps=2564290  (every counter 0, as above)
             mode (a) vs mode (b): 0 of 1548288 rays disagree
    shape    (a) vs OCCT: 1548288/1548288 rays identical, LOST=0 extra=0 displaced=0 kind=0
                          worst dt=1.852e-10 cm   (7/7 part(s) fully clean)
    shape    (b) vs OCCT: 1548288/1548288 rays identical, LOST=0 extra=0 displaced=0 kind=0
                          worst dt=1.856e-10 cm   (7/7 part(s) fully clean)

    surface  2654208/2654208 rays identical, LOST=0 extra=0 displaced=0, every counter 0,
             12/12 parts fully clean, both modes, 0 of 2654208 rays disagree between them
    mesh     1998259/2654208 identical, LOST=59925 extra=10429 displaced=1309104, 0/12 clean
```

**Mode (a) and mode (b) implement the placement differently and agree on every one of 1 548 288
rays.** Mode (a) transforms each ray into the shape's frame by hand; mode (b) leaves the rays in
the part frame and puts the matrix on the `TGeoNode`, so ROOT performs the transform. A wrong
placement in either would show up as a mode-(a)-vs-(b) disagreement *and* as LOST crossings against
OpenCascade. Both are zero.

### 7.3 `BasePin` alone, 96 beams × 24² — the placed part, isolated

Run against the post-change gate DB with `--reuse-db --parts BasePin`, so the three
representations of one part are directly comparable:

| repr | class | rays identical vs OCCT | LOST / extra / displaced | worst dt | robustness counters | chord V (cm³) | Capacity vs exact |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: |
| `surface` | `O2BVHSurfaceSolid` | 55296/55296 | 0 / 0 / 0 | 5.15e-13 cm | all zero | 31.463495 | −9.047e-16 |
| **`shape`** | **`TGeoTube`** | **55296/55296** | **0 / 0 / 0** | **4.72e-13 cm** | **all zero** | **31.463495** | **−6.785e-16** |
| `mesh` | `O2Tessellated` | 19781/55296 | 8 / 0 / 63625 | 2.56e-02 cm | all zero | 31.449424 | −4.144e-04 |

Mode (a) vs mode (b): **0 of 55296 rays disagree** on every representation. The exact solid and the
placed primitive agree on the chord volume to **7.3e-15** against OCCT's own chord integral over
the same rays — i.e. the two independent exact representations of this part now land on the same
number, which they could not be said to do while one of them was a composite with a sampled
capacity.

---

## 8. What is now wrong in the other documents

Corrected in this change:

* `csg/primitives.py`'s module docstring described the self-union as a permanent constraint. It now
  describes the placement, and states that `build_occ()` is deliberately unchanged.
* `DetectorsBase/O2SolidHarness.h`'s statement of the `shape_<part>.root` convention said the shape
  is "in the part's own LOCAL frame … with no placement matrix applied". It now documents the
  optional `placement` key and that its absence means the identity.
* `NEXT.md` and `Tutorial.md` both carried the "no `TGeoShape` carries a rigid transform, so a
  placed primitive is a composite" bullet as a live constraint. Updated.
* `testBVHSurfaceSolid.cxx`'s `makePlacedTube()` comment described the self-union as what the
  emitter does. It is now what the emitter *did*, and the helper survives as the reference the new
  equivalence case measures against.

Left alone deliberately: `Stream_H_CSGEmitter.md` §2.2 is a record of what that session measured
and is still true as history; a pointer to here has been added rather than rewriting it.

## 9. What this does not do

* It does not touch genuine booleans. A multi-leaf union is still a `TGeoCompositeShape` with
  Monte-Carlo `Capacity()` and `capacityComparable=false`, and its acceptance is still the OCCT
  symmetric difference. 7 of the 11 accepted parts across all corpora are in this class.
* It does not make more parts recognisable. The recogniser's scope is unchanged.
* It does not fix the oTOF converter failure (§5, finding 1), which is the corpus where this change
  would pay best.
