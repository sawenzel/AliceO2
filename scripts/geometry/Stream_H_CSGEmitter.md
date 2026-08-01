# Stream H — the CSG emitter: Bagger as a mixture of CSG, exact surfaces and mesh

Date: 2026-08-01. Delivers steps 2, 3 and 4 of the MVP in [`Tutorial.md`](Tutorial.md) §6, on top
of Stream G's any-`TGeoShape` gate (step 1) and Stream A's census. Written by the session that
built it; the integrating session folds it.

**In one paragraph.** Bagger now converts as a *mixture*: **7 leaf solids as native ROOT CSG
shapes, 5 as `O2BVHSurfaceSolid`, 1 as a mesh**, chosen by one hook in `O2_CADtoTGeo.py` and
recorded per part with the evidence that admitted it. The seven CSG parts are `BasePin` (a
`TGeoTube`) and all six "ram" parts (`TGeoTube ∪ TGeoTube`), and every one of them has an OCCT
symmetric-difference volume against the CAD solid of **exactly 0 cm³** and **0 disagreements with
the oracle on all four navigation columns**. The two acceptance tests never disagreed. Nothing
pre-existing moved: `ctest -R BVHSurfaceSolid` 78 → **80** cases green, fixtures **9/9**, Bagger
**9/12** with the same three capacity failures at the same values, unexplained oracle
disagreements on the surface representation **0 on every column of both corpora**, and a
field-by-field `gate.json` diff shows **zero non-path changes** on either corpus.

---

## 1. The result

### 1.1 The tiered table — every Bagger leaf solid, and what carries it

13 leaf solids. (The census counts 13 prototypes; the gate scores 12, because `Bucket` has no
exact-surface sidecar and therefore never enters the part database. Both numbers appear below,
which is the only way to keep them straight.)

| leaf solid | carried by | what it is | dV_sym (cm³) | band (cm³) | oracle: contains / distout / distin / safety |
| --- | --- | --- | ---: | ---: | --- |
| `BasePin` | **CSG** | `TGeoTube(0, 1, 5)` | **0** | 6.91e-06 | **0 / 0 / 0 / 0** |
| `BoomCylinderInner` | **CSG** | `TGeoTube(0, 0.6, 8.245) ∪ TGeoTube(0.7, 1.2, 0.75)` | **0** | 8.57e-06 | **0 / 0 / 0 / 0** |
| `BoomCylinderOuter` | **CSG** | `TGeoTube(0.6, 1, 7.991) ∪ TGeoTube(0.7, 1.5, 1.5)` | **0** | 2.11e-05 | **0 / 0 / 0 / 0** |
| `StickCylinderInner` | **CSG** | `TGeoTube(0, 0.6, 8.245) ∪ TGeoTube(0.7, 1.2, 0.75)` | **0** | 8.57e-06 | **0 / 0 / 0 / 0** |
| `StickCylinderOuter` | **CSG** | `TGeoTube(0.6, 1, 7.991) ∪ TGeoTube(0.7, 1.5, 1.5)` | **0** | 2.11e-05 | **0 / 0 / 0 / 0** |
| `BucketCylinderInner` | **CSG** | `TGeoTube(0, 0.35, 8.35) ∪ TGeoTube(0.5, 0.9, 0.75)` | **0** | 5.34e-06 | **0 / 0 / 0 / 0** |
| `BucketCylinderOuter` | **CSG** | `TGeoTube(0.35, 0.6, 8.1) ∪ TGeoTube(0.4, 1, 0.75)` | **0** | 1.15e-05 | **0 / 0 / 0 / 0** |
| `Base` | surface | 44 faces | — | — | 0 / 0 / 0 / 0 (surface) |
| `Boom` | surface | 31 faces | — | — | 0 / 0 / 0 / 0 (surface) |
| `Stick` | surface | 24 faces | — | — | 0 / 0 / 0 / 0 (surface) |
| `BucketLink1` | surface | 22 faces | — | — | 0 / 0 / 0 / 0 (surface) |
| `BucketLink2` | surface | 23 faces | — | — | 0 / 0 / 0 / 0 (surface) |
| `Bucket` | mesh | 97 faces, incl. 4 spherical + 2 toroidal | — | — | not in the part DB (no exact sidecar) |

**Tiers: CSG 7, exact surfaces 5, tessellated 1.** That is the honest shape of the answer — a
single coverage fraction cannot describe a system with three representations, and the project
asked for the tiered form explicitly.

Each CSG part also carries a third, cheaper measurement made at emission time: **ROOT-vs-CAD
`Contains`, 0 disagreements out of 4000** random points classified against the original OCCT solid
with `BRepClass3d_SolidClassifier`, on every one of the seven. That check exists because it is the
only one that tests the *ROOT* realisation against the *CAD body*; see §3.

### 1.2 Why the five declines are declines, per part

The recogniser's scope is Tier 1 (whole-part primitives) plus the one measured Tier-2 case. Every
decline names the structure it found:

| leaf solid | reason |
| --- | --- |
| `Base` | 7 axis clusters (44 faces: 16 cylinder, 28 plane) |
| `Boom` | 6 axis clusters (31 faces: 14 cylinder, 17 plane) |
| `Stick` | 6 axis clusters (24 faces: 12 cylinder, 12 plane) |
| `BucketLink2` | 3 axis clusters (23 faces: 11 cylinder, 12 plane) |
| `BucketLink1` | 2 clusters, but a planar face that is neither a cap nor a wedge of either (22 faces) |
| `Bucket` | a toroidal face — out of the recogniser's scope |

These are Tier-3 bodies. The census predicts exactly this: it found **one** Bagger prototype
matching a whole-part template and named the six rams as the two-cluster case, and those are
precisely the seven that converted. Nothing here contradicts the census; this is the census's
prediction, executed.

`BucketLink1` is the one worth a second look for a future session: it *is* two axis clusters, and
only one plane keeps it out. It is the cheapest next part in the corpus.

### 1.3 The boolean ladder, run through the same pipeline

The nine synthetic fixtures are the second corpus, and they are the harder test of *judgement*
because most of them are deliberately not primitives.

| fixture | outcome | detail |
| --- | --- | --- |
| `box` | **CSG** | `TGeoBBox(1, 1.5, 2)`, dV_sym **0**, emitted as a **bare `TGeoBBox`** — independently reproducing the fixture Stream G had to hand-build, including its analytic capacity (1.48e-16 relative) |
| `cyl_cross_cyl` | **CSG** | `TGeoTube(0,1,3) ∪ TGeoTube(0,1,3)`, dV_sym **0** |
| `cyl_inter_cyl` | rejected | proposed as a union; dV_sym 1.90 cm³ vs band 1.6e-06 |
| `tube_window` | rejected | proposed as a union; dV_sym 6.04 cm³ vs band 8.0e-06 |
| `box_minus_cyl` | declined | a planar face that is neither cap nor wedge (it is a difference) |
| `box_union_box` | declined | 10 planar faces: not a six-plane box |
| `cyl_plus_cone` | declined | mixed cone and cylinder on one axis |
| `sphere_minus_cyl` | declined | a sphere with additional faces (no theta/phi cuts in scope) |
| `torus_union_cyl` | declined | toroidal face |

Gate: **shape 3/3 parts with zero disagreements** on all four columns (the two above plus Stream
G's hand-built `box_minus_cyl`, which this recogniser declines and therefore left untouched).

---

## 2. The trap that had to be solved first, and what it cost

### 2.1 Can the converter write a `TGeoShape` from Python? Yes — measured, not assumed

The brief flagged this as possibly the first real obstacle. It was checked before anything was
designed around it, and the answer has two halves:

* **The alibuild Python 3.10 that pythonOCC is built against *is* the interpreter the O2
  environment provides.** With pythonOCC's site-packages prepended to `PYTHONPATH` and OCCT's
  `lib` prepended to `LD_LIBRARY_PATH`, one process imports **both** `OCC` and `ROOT`
  (`OCC 7.9.0`, `ROOT 6.36.04`), and `TFile::WriteTObject(shape, "shape")` round-trips a
  `TGeoTube` immediately. So the converter, run from an O2 shell, writes its own `.root` files
  and no intermediate step is needed.
* **Inside `runOracleGate.py` it does not hold.** That script builds the converter's environment
  by *replacing* `PYTHONPATH` and `LD_LIBRARY_PATH` with OCC-only values, which drops ROOT; the
  failure is not a missing module but `ImportError: libffi.so.6` from `cppyy`'s `ctypes` import,
  three levels down. Diagnosing that from the symptom would have cost the session the brief
  predicted.

So the hook does both, and says which it did. It **always** writes the description as
`csg_<VOL>_<LID>.json` — the small intermediate description the brief suggested — and writes
`shape_<VOL>_<LID>.root` when ROOT is importable. `csg/emit.py --from-json <dir>` completes a
deferred set afterwards, and a part whose `.root` file was not written is **never** referenced by
`geom.C`: the converter does not emit a reference to a file it did not write.

`csg/occ_env.py` now *prepends* to those two variables instead of replacing them. That one-line
change is what makes the single-interpreter path work at all, and it is a strict improvement:
prepending keeps OCCT's own libraries winning any conflict while leaving the rest of the
environment intact.

### 2.2 No `TGeoShape` in ROOT can carry a rigid transform

This is the real finding of the emission step, and it shapes the output format.

Bagger's leaf solids are **not** at the origin and their axes are **not** coordinate axes:
`BasePin` is a cylinder on `(0, 0, −1)` through `(0, 5.916, 0)`; `BoomCylinderInner`'s rod runs
along `(0, 0.707, −0.707)`. A recognised `TGeoTube` therefore has to be *placed*. In ROOT 6.36:

* `TGeoBBox` carries a translation, through `fOrigin`, and no rotation;
* every other primitive is fixed to the origin with its axis along z;
* the only `TGeoShape` subclass holding a `TGeoMatrix` at all is `TGeoCompositeShape`, through its
  `TGeoBoolNode` — which needs **two** operands.

So a placed primitive is written as **the union of the primitive with an identical copy of itself
under the same matrix**. That is set-theoretically the primitive, and it was measured to be
exactly that before anything was built on it: 20000 random probes against the same tube queried in
its own frame gave **0 disagreements on `Contains`, `Safety`, `DistFromOutside` and
`DistFromInside`**, bit-identical, with a negative control (`rmax` 1.00 vs 1.05) that does move
the count. The same claim is now a unit test, `CsgSelfUnionCarriesARigidTransformExactly`.

Two consequences, both reported rather than hidden:

1. **The shape's `Capacity()` becomes Monte-Carlo**, because it is a composite. Measured spread
   over repeated calls on the ram fixture: **1.8e-02 relative**, four orders of magnitude above
   the gate's 1e-06 band. This is exactly Stream G's rule — *report it, mark it, never gate on
   it* — and it is why the acceptance test for a CSG part is a symmetric-difference volume.
2. **An axis-aligned box escapes all of this.** A six-plane box is symmetric under flipping and
   permuting its own axes, so when those axes *are* the coordinate axes the frame is relabelled
   into the identity and a **bare `TGeoBBox` with its own origin** is emitted — one object instead
   of three plus a matrix, and an analytic `Capacity()`. That is not cosmetic for the corpus at
   large: the census counts **62560 placed six-plane boxes in oTOF alone**.

### 2.3 The first version of the frame relabelling was wrong, and the bounding box did not say so

An early `_match_box` computed the box centre from a point difference and got the sign of the
plane pairing wrong; it declined every box with "separations are not positive". Fixed by stating
each pair as *offsets along the outward normal*, where the face carrying that normal is by
definition the one at the larger offset, so the half-thickness is positive whichever face is found
first.

Worth recording because of what *did not* catch it: the OCCT-vs-ROOT bounding-box comparison. See
§3.

---

## 3. The two acceptance tests, and how they were kept independent

| | test | what it measures | what it cannot see |
| --- | --- | --- | --- |
| 1 | **OCCT symmetric-difference volume** — `BRepAlgoAPI_Cut` both ways, `BRepGProp` on each | the *OCCT* realisation of the description against the CAD solid | anything about the ROOT realisation |
| 2 | **the ordinary oracle gate** on `shape_<part>.root` | the *ROOT* realisation against OpenCascade's answers to the harness's own sample set | volume; it samples |

They are independent by construction: **one description, two builders** (`csg/primitives.py`).
A transposed rotation in the ROOT builder alone leaves test 1 perfect and breaks test 2, which is
the point of not sharing an object between them.

**Did they ever disagree? No.** On both corpora, every part that passed test 1 scored 0
disagreements on all four columns in test 2, and no part passed test 2 while failing test 1.

### 3.1 The band, and how sharp it is

The residual is a *volume*; the model tolerance is a *length* (1.0e-07 cm on both corpora). They
are related through the boundary — displacing a surface of area `A` by `t` moves `A·t` — so the
criterion applied is

    dV_sym  ≤  modelTolerance × area(original)

and the ratio `dV_sym / volume(original)` is reported alongside. Measured knife-edge, on a
2 × 3 × 4 cm box whose band is 5.2e-06 cm³:

| candidate | dV_sym | verdict |
| --- | ---: | --- |
| the same box | 0 | accepted |
| 0.5 × model tolerance too long | 0 (below OCCT's own boolean resolution) | accepted |
| 10 × model tolerance too long | 6.0e-06 (predicted 2·3·1e-06 = 6e-06) | **rejected** |
| 1e-04 cm too long | 6.0e-04 | **rejected** |
| a sphere of *equal volume* | 13.7 | **rejected** |

The last row is the reason this is a symmetric difference and not a volume comparison.

### 3.2 The guard fired on real geometry, twice, unprompted

The recogniser is deliberately greedy — it reads carriers and never inspects a trim curve —
because acceptance is exact. On the boolean ladder that asymmetry stopped being an argument:

| fixture | proposal | dV_sym | band | verdict |
| --- | --- | ---: | ---: | --- |
| `cyl_inter_cyl` | `TGeoTube ∪ TGeoTube` | **1.8997** | 1.6e-06 | **rejected** |
| `tube_window` | `TGeoTube ∪ TGeoTube` | **6.0354** | 8.0e-06 | **rejected** |

Both are *intersections* / single cells, not unions — exactly as the census says (`Stream_A_CSG.md`
§3, Tier 2) — and the recogniser proposed a union for both because their carrier structure is two
axis clusters. Neither was ever going to be caught by a heuristic. The volume caught both, by six
orders of magnitude.

That is also the answer to "recognition must never silently change geometry": on this corpus it
tried twice and was stopped twice.

### 3.3 A number that looks like a failure and is not

`bbox(ROOT vs OCCT)` is 0.13–0.62 cm on the six rams, and 0 on `BasePin`, the `box` fixture and
`cyl_cross_cyl`. It is **not** a frame error. `TGeoBoolNode::ComputeBBox` takes each operand's own
axis-aligned box, transforms its eight corners and unions the results, so a leaf whose axis is not
a coordinate axis contributes a box strictly larger than itself. It is the same effect Stream G
reports for the surface representation's conservative per-face AABBs — and indeed the gate's
`bboxDev` column shows the *same* values for the `surface` and `shape` representations of those
parts, to three digits.

Because that check is blunt for composites, the sharp one was added: **`crosscheck_contains`**,
4000 random points classified against the original CAD solid with `BRepClass3d_SolidClassifier`
and against the emitted `TGeoShape`. 0 disagreements on all seven Bagger parts and all CSG
fixtures. This is the check that would catch a transposed rotation immediately, and the bounding
box would not.

### 3.4 Traps the acceptance test is written around

`BRepGProp::VolumeProperties` returns exactly **0** for a single planar face, and OCCT's boolean
operations do not report an empty or degenerate result as an error. So `dV_sym = 0` is *also* what
a failed comparison produces. Three explicit guards therefore run before the verdict is believed:
both cuts must report `IsDone()`; the original's volume and area must be strictly positive; and
the candidate's own volume must be within a factor of 4 of the original's *before* the cuts are
attempted. `accept.self_test()` asserts all of it, with rejecting cases.

---

## 4. What was built

Everything is in `scripts/geometry/csg/`, plus **one** hook in `O2_CADtoTGeo.py`.

| file | what it is |
| --- | --- |
| `primitives.py` | the intermediate description (placed `TGeoBBox`/`TGeoTube`/`TGeoTubeSeg`/`TGeoCone`/`TGeoSphere`, or a union of two) and the two builders that realise it, in OCCT and in ROOT |
| `recognise.py` | carriers → axis clusters → a placed primitive, or a decline carrying the part's structure |
| `accept.py` | the symmetric-difference volume, its band, and 11 self-checks |
| `emit.py` | the driver: `--db`, `--brep`, `--from-json`, `--self-test`; writes `shape_*.root`, `csg_*.json` and the per-part report |
| `hook.py` | the converter's single integration point, the cascade report, and `geom.C`'s `LoadShape` |
| `occ_env.py` | now prepends rather than replaces (§2.1) |

`O2_CADtoTGeo.py` gains `--csg off|auto|required` and `--csg-report PATH`; **default `off`, so a
conversion that does not ask for CSG is unchanged.** With `auto` the per-part choice becomes
**CSG → exact surfaces → tessellated**, `geom.C` loads a native `TGeoShape` per CSG part, and
`csg_report.json` records the tier and the evidence for every leaf solid.

**The other representations are still written for a CSG part.** The gate scores all three side by
side (Stream G), and dropping the sidecar would quietly shrink the corpus that defends this
project's zero-disagreement invariant. The cascade decides what `geom.C` *builds*, which is the
decision that matters in production.

### 4.1 The emitted macro was actually built

Not just written: `geom.C` from a `--csg auto` conversion was loaded in ROOT and exported.

| | baseline `geom.C` | `--csg auto` |
| --- | --- | --- |
| builds and exports | yes | yes |
| illegal overlaps / extrusions | 0 | **0** |
| exported `geom.root` | 181 403 B | **127 420 B** |

(Both fail to compile under ACLiC — `geom.C`'s prelude does not `#include <TGeoTessellated.h>`.
That is pre-existing, identical with and without `--csg`, and untouched here.)

### 4.2 The two entry points agree

The converter hook recognises the in-memory leaf solid; `csg/emit.py --db` recognises the
`.brep` the converter wrote. Comparing the 12 shared descriptions field by field: **23 leaf fields
differ across 12 parts, worst |Δ| = 2.3e-13 cm** — the `.brep` file's decimal round trip, ~1e-15
relative on 60 cm coordinates. No verdict, recogniser or leaf count differs.

---

## 5. Self-checks: 28, each positive case paired with a negative one

`csg/emit.py --self-test` runs both suites and refuses to emit if any fails.

**Acceptance (11):** a shape against itself; an independently built description against an OCCT
box; the knife-edge table of §3.1; a sphere of equal volume; a tube built two ways; a solid
cylinder rejected as a tube; a rotated and translated tube against its placed description; and the
same tube *without* the rotation, rejected.

**Recognise/emit (17):** one per primitive in scope (box, solid cylinder, tube, tube segment,
cone, sphere), a placed tube, the rod-and-eye union; four negative controls that must decline or
be rejected (a blind bore, an L-plate, a torus, a cylinder with a milled flat); the emitted placed
tube's `Contains` against the closed form over 20000 points **plus a control proving that check
can fail**; a two-leaf union through the file convention; the OCCT/ROOT bounding-box agreement;
and the bare-`TGeoBBox` emission path.

Two defects were found *by* these checks rather than by inspection: the box frame sign error
(§2.3), and a negative control on a 0.05 cm shell that a coarse lattice missed entirely — it
reported 0 disagreements and would have advertised a control that could not fire. The probes are
now placed inside the shell and mapped out to the master frame.

Two C++ cases were appended to `Detectors/Base/test/testBVHSurfaceSolid.cxx` in a block ending
`// --- Stream H: the CSG emitter's ROOT-side claims ---`. They pin the two properties of *ROOT*
that the emitted file silently depends on — if a future ROOT changes either, every CSG part this
project writes becomes wrong geometry that still loads:

* `CsgSelfUnionCarriesARigidTransformExactly` — the §2.2 idiom, over 2000+ probes on all four
  queries, with a negative control;
* `CsgTwoLeafUnionRoundTripsAndMatchesTheClosedForm` — a barrel-plus-lug union against a
  hand-written membership function, before and after `saveShapeToRootFile`/`loadShapeFromRootFile`,
  plus the assertion that its `Capacity()` **scatters** between calls (a spread > 1e-4 proves the
  method is sampled, so no capacity criterion can ever apply to a shape written this way).

---

## 6. Before and after — totals and disagreement counts, quoted together

| | before | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 78 cases, green | **80 cases, green** (2 appended) |
| fixtures gate | 9/9 | **9/9** |
| Bagger gate | 9/12 (capacity 1.39/1.31/1.39e-06) | **9/12, the same three, the same values** |
| unexplained disagreements, **surface**, fixtures | contains 0, distout 0, distin 0, safety 0 | **0, 0, 0, 0** |
| unexplained disagreements, **surface**, Bagger | contains 0, distout 0, distin 0, safety 0 | **0, 0, 0, 0** |
| disagreements, **shape**, fixtures | 0/0/0/0 (2 hand-built parts) | **0/0/0/0 (3 parts)** |
| disagreements, **shape**, Bagger | — (no parts) | **0/0/0/0 (7 parts)** |

**Column diff**, `gate.json` compared leaf by leaf with `timing*`, `*Seconds` and `nsPerCall`
stripped, baseline versus a full reconversion with the modified converter (`--csg` defaulting to
off):

| corpus | leaf fields | removed | added | changed | of which non-path |
| --- | ---: | ---: | ---: | ---: | ---: |
| fixtures (9 parts) | 8315 → 8315 | 0 | 0 | 20 | **0** |
| Bagger (12 parts) | 14571 → 14571 | 0 | 0 | 24 | **0** |

Every changed field is a `representations[..].source` path, differing only because the two runs
used different scratch directories. **No numeric field moved on either corpus.**

---

## 7. The one thing the MVP text expects that did not happen, and why

`Tutorial.md` §6 says success would mean "all 12 passing the oracle gate, with the six rams exact
as booleans — which would also retire the last three capacity failures rather than tuning them
away."

The six rams **are** exact as booleans — symmetric difference 0 cm³, which is a far stronger
statement than any capacity band. But **the gate total is still 9/12**, and it will stay there
until someone decides to change it, because Stream G fixed the gate's pass/fail verdict on the
*surface* representation deliberately (`Stream_G_AnyShape.md` §2, "Not changed"). The three
failures are `capacityRelativeDeviation` of the `O2BVHSurfaceSolid` form of
`BoomCylinderOuter`, `StickCylinderOuter` and `BucketCylinderInner` — a property of that
representation's trim data, untouched by this stream and explicitly not a priority.

So the honest statement is: **those three parts are now carried by a representation that is exact,
and the number that reports them as failing describes a representation they no longer ship in.**
Closing that gap means making the gate's verdict representation-aware — deciding per part which
representation is the shipped one and gating *that*. It is a small change to `runOracleGate.py`
and a policy decision about what "the gate" means, which is why this stream reported it rather
than making it. Note also that a capacity criterion cannot simply be moved across: the CSG form of
those parts is a composite, whose `Capacity()` is Monte-Carlo (§2.2).

---

## 8. Reproducing everything here

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
ctest -R BVHSurfaceSolid                       # 80 cases

cd $HOME/alisw/O2/scripts/geometry
python3 csg/emit.py --self-test                # 28 checks

# the whole cascade, from the STEP file (needs ROOT and OCC in one interpreter):
SW=$HOME/alisw/sw/ubuntu2404_aarch64
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$LD_LIBRARY_PATH
python3 O2_CADtoTGeo.py STEP_examples/Bagger.step --output-folder /tmp/csgconv -o geom.C \
    --exact-surfaces auto --mesh --dump-brep --csg auto

# or, on a finished gate workdir, which is the loop Stream G designed for:
cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate \
    --model scripts/geometry/STEP_examples/Bagger.step
python3 scripts/geometry/csg/emit.py --db /tmp/gate/db --report /tmp/gate/csg_report.json
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --skip-convert
```

**One trap in that last loop, which cost a run here:** `manifest.json` stores **absolute** paths.
Copying a finished workdir and emitting into the copy scores the *original*'s files and silently
reports no `shape` representation at all. Either emit into the workdir itself, or rewrite the
paths in `manifest.json` after copying.

---

## 9. What this leaves open

- **`BucketLink1`** is two axis clusters and one awkward plane away from converting. Cheapest next
  part in the corpus.
- **The gate's verdict is representation-blind** (§7). Until that changes, the tiered scorecard is
  where the CSG result lives and the gate total describes something else.
- **Tier 0 is untouched and is still the highest-value item in the programme** (`Stream_A_CSG.md`
  §2.3): 1004 of ALICE3's B-spline surfaces are quadrics in disguise. This recogniser reads
  `BRepAdaptor_Surface` types directly, so it would see none of them today — on `as1-oc-214`,
  whose five prototypes are *all* quadric after Tier 0, it currently recognises nothing.
- **`primitive − N holes`** is not built. Two ladder fixtures (`box_minus_cyl`, `sphere_minus_cyl`)
  and `Bucket`'s bores want it, and it needs a subtraction node, which this stream had no measured
  case for on Bagger.
- **The tessellated tier is one part wide and that part is the worst one.** `Bucket` falls to the
  mesh because of 2 toroidal and 4 spherical faces. A `TGeoTorus`/`TGeoSphere`-aware Tier 1 would
  not help — `Bucket` is 97 faces — but it is worth noting that the single mesh part in the whole
  MVP is there for a reason the exact-surface path, not the CSG path, would have to fix.
- **No timing.** `CSG_Pipeline.md` §7's composite-depth probe is still unrun; these composites are
  two leaves deep, which is the shallow end of that question.
