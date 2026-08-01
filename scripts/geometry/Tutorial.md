# CAD → TGeo: what exists today, how to drive it, and what is missing

Orientation document, written 2026-08-01 after wave 1. If you have been away from this branch,
**read this one and nothing else first.** The other documents are deep records of single
investigations; this is the map.

Everything below is measured on this branch (`swenzel/bvhsurfacesolid`) unless marked *planned*.

---

## 1. The goal, and the three ways of reaching it

Give a CAD part the **best possible TGeo representation** — however that is achieved. "Best"
means, in priority order: exact geometry, watertight for navigation, cheap to query, and small.

There are three representations, and they are complementary rather than competing:

| # | representation | what it is | exact? | covers |
| --- | --- | --- | --- | --- |
| 1 | **CSG / primitives** | `TGeoTube`, `TGeoBBox`, … and `TGeoCompositeShape` booleans of them | yes, analytically | mechanical CAD: quadric-only parts |
| 2 | **`O2BVHSurfaceSolid`** | the part's real B-rep faces (trimmed analytic surfaces) carried into TGeo, with a BVH for ray queries | yes, per face | anything whose faces are plane/cylinder/cone/sphere/torus, however trimmed |
| 3 | **`O2Tessellated`** | a triangle mesh | no | everything, as the fallback |

The production behaviour is a **cascade**: try 1, else 2, else 3, and record per part which one
accepted it. All three now exist: `--csg auto --exact-surfaces auto --mesh` gives the full
cascade, and `csg_report.json` records the choice and its evidence per part.

Measured on `Bagger.step` — **CSG 7, exact surfaces 5, tessellated 1** of 13 leaf solids, with
every CSG part exact (`dV_sym = 0`) and oracle-clean on all four columns.

---

## 2. The pipeline

```
   STEP file
      │
      │  O2_CADtoTGeo.py            (pythonOCC / OpenCascade)
      ▼
   per leaf solid:
      ├── surfaces_<part>.bin       exact B-rep sidecar   → O2BVHSurfaceSolid
      ├── facets_<part>.bin         triangle mesh         → O2Tessellated
      ├── brep_<part>.brep          the OCCT solid itself → the oracle's input
      └── geom.C                    ROOT macro that builds the TGeo volumes
      │
      ▼
   ROOT / o2-sim
```

And, alongside it, the thing that makes any of this trustworthy:

```
   brep_<part>.brep ──► occtOracle.py ──► answers_<part>.json   (ground truth, OpenCascade)
                                                │
   surfaces_<part>.bin ──► O2BVHSurfaceSolid ───┤
                                                ▼
                                    o2-bench-...-solid-harness  ──► gate.json  (PASS / FAIL)
```

### Artefacts, in one line each

| file | written by | meaning |
| --- | --- | --- |
| `surfaces_<part>.bin` | converter, `--exact-surfaces` | the exact sidecar: surface records + trim wires + (v3) edge identity |
| `facets_<part>.bin` | converter, `--mesh` | `uint32 nTriangles` then 9 × float32 per triangle |
| `brep_<part>.brep` | converter, `--dump-brep` | the OCCT solid the sidecar was extracted from — the oracle's input, **not** used at run time |
| `shape_<part>.root` | converter, `--csg` | one object inheriting `TGeoShape` under key `shape`, in cm, in the part's local frame → the CSG representation |
| `csg_<part>.json`, `csg_report.json` | converter, `--csg` | the recognised description and its acceptance evidence, written whether or not the part was accepted |
| `geom.C` | converter | ROOT macro exporting `get_builder_hook_unchecked()`; loads the binaries relative to itself |
| `surface_report.json` | converter, `--surface-report` | per-face extraction result and why anything failed |
| `manifest.json` | `makeTestPartDB.py` | pairs each part's sidecar with its mesh — the harness's index |
| `samples_<part>.json` | harness, `--dump-samples` | the seeded point/ray set; **nothing outside the harness can regenerate it** |
| `answers_<part>.json` | `occtOracle.py` | OpenCascade's answers to exactly those samples |
| `gate.json` | `runOracleGate.py` | the full per-part scorecard |
| `crossings_<part>.json`, `xray.json` | `runXRayBench.py` | per-ray crossing lists from OpenCascade, and the transport scorecard |

---

## 3. How to drive it

Environment first — this has cost several sessions an hour each:

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH   # stage FIRST
```

`eval` and `export` must be in the **same shell command** as whatever you run next, or the
harness dies on `libarrow_acero.so`.

**Convert a model.**

```bash
python3 scripts/geometry/O2_CADtoTGeo.py my.step --output-folder out -o geom.C \
    --exact-surfaces auto --mesh --dump-brep --surface-report out/surface_report.json
```
`--exact-surfaces auto` emits a sidecar when **every** face of a leaf solid extracts exactly, and
silently keeps the mesh otherwise. `required` fails instead of falling back.

**Run the acceptance gate** — this is the real interface to the project:

```bash
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
```
It chains convert → dump samples → oracle → score. It exits non-zero whenever any part fails,
which is the normal state; **read the counts, not the exit code**.

Useful flags: `--skip-convert` (re-score without reconverting — valid only if the *converter* did
not change), `--transform translate:0,0,4000` / `scale:0.1` with `--load-samples` (the
position/scale sweep).

**Interrogate one part.**

```bash
$B/stage/bin/o2-bench-detectorsbase-solid-harness --db /tmp/gate/db --parts BoomCylinderInner \
    --points 1 --rays 1 --rims --edge-identity
```
`--rims` names which trim loop of which face is open and where; `--edge-identity` prints the
sidecar-v3 topology block; `--loop-crosscheck` re-runs every query through the non-BVH twin and
requires bit-identical answers.

**Transport geantinos through it** (§5.7) — the only tool here that *steps* rather than asking
single-shot questions, and the only one that runs on ALICE3:

```bash
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xray --fixtures --beams 96
```

---

## 4. The principles this project runs on

These are not style preferences; each was bought with a lost session.

1. **The oracle is the acceptance criterion, not a mesh.** A part passes only if it is a closed
   navigable manifold *and* agrees with OpenCascade outside the model's own declared tolerance.
   Comparing against a tessellation "within a chording band" cannot certify exactness and
   demonstrably hid whole classes of defect.
2. **Gate totals and disagreement counts are separate numbers. Never quote one without the
   other.** The invariant to defend is *zero unexplained disagreements on every column*; a change
   that improves a total while breaking that is a regression.
3. **Diff columns, not verdicts.** A bit-identical numeric column can accompany a real change,
   and has, twice. Strip `timing*`, `*Seconds`, `nsPerCall` from `gate.json` and compare the rest.
4. **Evidence before theory.** Build the instrument that *names* the object, then read it. The
   last five advances came that way; several confident hypotheses died on contact with a
   measurement, and the documents record the death rather than explaining it away.
5. **A refuted experiment is not a refuted hypothesis.** A sweep once said "no" to the correct
   answer while being structurally incapable of saying "yes". When an experiment fails to move a
   number, check that it *could* have. Corollary: any sweep reporting "no change" needs a
   positive control.
6. **A question asked twice becomes a harness flag.** That is where `--rims`, `--edge-identity`
   and `--load-samples` came from.
7. **Topology is decided by identity, never by proximity** (since wave 1 — §5.4).

---

## 5. Where each capability actually stands

### 5.1 `O2BVHSurfaceSolid` — working, and the furthest along

A `TGeoBBox` subclass carrying the part's real faces. Six surface kinds: planar polygon, curved
planar, cylindrical, spherical, conical, toroidal. Trim wires of lines, arcs and rational
B-splines. Closure, capacity (divergence theorem, closed form), `Contains` by ray parity, ray
distances, `Safety`, a BVH plus a `_Loop` twin that must agree bit-for-bit.

**Measured today:** ctest 75 green; ladder fixtures **9/9**; Bagger **9/12**; all 21 parts across
both corpora `reliable` and navigable; **zero unexplained oracle disagreements** on `contains`,
`distin`, `distout`, `safety`, both corpora.

**Missing for completeness:**
- **Free-form (B-spline) *surfaces* are not supported at all.** This is the ALICE3 coverage lever
  — after Tier-0 decoding (§5.2) it is 19 of 55 solids and 1373 of 9266 faces.
- Three Bagger parts fail **capacity** at 1.3–1.4e-06 relative against a 1e-06 band. Quadrature is
  converged; the residual is trim data.
- **A real correctness bug in the quartic solver** (torus intersections only) — §5.5.
- Never exercised at ALICE3 scale (9266 faces, metre coordinates).

### 5.2 B-rep → CSG — census **and** a working Bagger-scope emitter

Status: **recogniser, emitter and acceptance all exist and are wired into the converter**
(`--csg auto|required|off`, default `off`). Scope is deliberately narrow: whole-part primitive
templates plus the two-coaxial-cluster `tube ∪ tube`. General decomposition (Tier 3) is *not*
built and remains gated on a cell-count table that does not exist.

Each candidate must pass **two independent tests**: OCCT symmetric-difference volume
(`BRepAlgoAPI_Cut` both ways) inside the model tolerance, *and* the ordinary oracle gate on the
emitted shape. They have never disagreed; the volume test twice rejected a greedy union proposal
(`cyl_inter_cyl`, `tube_window`) that no heuristic would have caught — both are single cells, not
unions.

The census that scoped all of this is complete and it changed the plan substantially:

- **1004 of ALICE3's 2377 "B-spline" surfaces are quadrics in disguise** — cylinders, cones,
  spheres, planes written as NURBS by the exporter — recoverable at a max gap of 3.4e-8 mm.
  Decoding them ("Tier 0") takes quadric-only solids from **15/55 to 36/55**. This is the single
  highest-value item in the whole programme and it benefits representation 2 as much as CSG.
- Bagger is **13/13 quadric-only**. `as1-oc-214` is 0/5 before Tier 0 and **5/5** after.
  oTOF is 20 prototypes in 62628 placements, of which **19 are exact `TGeoBBox`es**.
- The property that matters for decomposition is **single-cell**, not convexity — a through hole
  has zero concave edges and is not convex.

**Missing for completeness:** Tier 0 (decoding NURBS-encoded quadrics — the highest-value item on
the board, and it helps representation 2 as much as CSG); Tier 3 for the parts that decline
(`Base`, `Boom`, `Stick`, `BucketLink1/2` — 2 to 7 axis clusters each); and coverage beyond
Bagger. Design in `CSG_Pipeline.md`, census in `Stream_A_CSG.md`, emitter in
`Stream_H_CSGEmitter.md`.

**Two things anyone touching this must know.** (i) The alibuild python3.10 that pythonOCC needs
*is* the O2 interpreter, so one process can import both `OCC` and `ROOT` — but only if
`PYTHONPATH`/`LD_LIBRARY_PATH` are **prepended** to. `runOracleGate.py` *replaces* them, and
PyROOT then dies on `libffi.so.6`; that is why the hook always writes `csg_<part>.json` and
completes deferred shapes via `emit.py --from-json`. (ii) **No `TGeoShape` in ROOT 6.36 can carry
a rigid transform** — `TGeoBBox` has `fOrigin`, nothing else has anything, and only
`TGeoCompositeShape` holds a matrix. A *placed* primitive is therefore emitted as its union with
an identical copy under the same matrix, which is why a recognised plain tube can appear as a
`TGeoCompositeShape` in `gate.json`.

**Important:** for Bagger-class geometry the emitter needs **no new O2 class**. `TGeoTube`,
`TGeoBBox`, `TGeoCone` and `TGeoCompositeShape` already exist in ROOT. A dedicated `O2CSGSolid`
(flat DNF) is only needed if native composite trees turn out too slow or too deep — that is a
later, gated decision, not a prerequisite.

### 5.3 The OCCT oracle and gate — working, and the most valuable asset here

`occtOracle.py` answers the harness's exact sample set from the `.brep` using
`BRepClass3d_SolidClassifier`, `IntCurvesFace_ShapeIntersector`, `BRepExtrema` and `BRepGProp`.
`runOracleGate.py` chains the four stages and scores per part, per column, with explicit
mismatch classes (`withinBand`, `missedSurface`, `unexplained`).

It scores **every representation a part has, side by side, against one oracle answer set** —
`surface`, `mesh` and `shape` columns in `gate.json`, each with its own disagreement counts, plus
a `bboxDeviationFromOracle` frame check per representation.

**Missing for completeness:**
- **The pass/fail verdict is still computed on the surface representation alone.** A part now
  shipping as exact CSG can therefore be reported as failing on a column belonging to a
  representation it no longer uses — see §6.
- No parallelism and no sample budget → an ALICE3-sized model is not yet gateable in bounded time.
- `makeTestPartDB.py` writes **absolute** paths into `manifest.json`, so copying a finished
  workdir and re-running with `--skip-convert` silently reads the *original* directory. It does
  not error; it quietly scores the wrong files. Rewrite the paths, or reconvert.

### 5.4 Sidecar v3 / edge identity — landed in wave 1

Closure used to decide a *topological* question (do these two faces share an edge?) with a
*proximity* query against a tuned band. It now decides on **identity**: the converter writes each
face's `TopoDS_Edge` ids, and closure requires every edge to appear exactly twice with opposite
sense — zero tolerance, zero sampling. The old proximity number survives as a *reported*
measurement, `GetMaxSharedEdgeDeviation()`.

Effect: parts failing on navigability went **7/21 → 0/21**, `tube_window` passes, Bagger 6/12 →
9/12. And the justification is not aesthetic: **ALICE3's edge tolerances are 4.7e-04 to 4.3e-03 cm
against Bagger's uniform 1.0e-07** — a geometric band would have had to be ~1e-3 cm, wide enough
to swallow real features. The band was never going to survive the next model.

Caveat that must always be quoted with the verdict: identity certifies that the source B-rep's
*topology* survived conversion. It does **not** certify each face's geometry — a mis-trimmed face
still reads closed.

### 5.5 Scale and position robustness — measured in wave 1, one real bug

- **Position-independent: yes.** The whole ladder at (0, 0, 400) cm is column-identical.
- **Scale-independent: no.** At ×0.1 the torus fixture starts disagreeing with the oracle.

Cause: `solveQuarticReal` guards its branches with `kTolerance` (**1e-9 cm, a length**) against
quantities of dimension **L²** and **L³**. Three guards share the defect. Consequences, both
confirmed on this branch:
- the resolvent guard fails ⇒ **the ray silently misses** a torus it does cross;
- the `termQ` guard misroutes an asymmetric quartic into the biquadratic branch ⇒ it returns
  **two confidently wrong roots** instead of four right ones (worse than a miss — a miss at least
  leaves a visible gap).

**How much should you care?** It is scale-dependent, not scale-exclusive: at the shipped fixture
size 1 ray in 5162 is already lost, which is why 2000-ray gates never saw it. Bagger has **2
toroidal faces out of 288**, so it is *not* an MVP blocker. ALICE3 has **350**, and small
toroidal features (fillets, bellows knuckles) are exactly where it bites. Fix it before ALICE3,
not before the MVP.

### 5.6 Tessellated fallback — works, and is what will not scale

`O2Tessellated` + `facets_*.bin`. Two measured limits: facets are stored **float32** (half-ulp
1.5e-05 cm at 400 cm), and meshing uses an **absolute** deflection for both linear and angular
terms, so converting a single 2 m sphere reached **22.9 GB resident** before being killed.
Tessellation, not the exact path, is the scaling problem.

Also still open: in `auto` mode an exact part whose reliability is not `Reliable` **ships anyway,
with only a warning**. That is the one live production risk on the branch.

### 5.7 X-ray / geantino transport — the newest instrument, and the only one that *transports*

`runXRayBench.py` + `o2-bench-detectorsbase-xray`. A structured beam is fired through a part and
each ray's **ordered crossing list** is produced by *stepping* — two independent ways — and
compared against OpenCascade **as a list, not as a rate**, so a missing crossing localises to a
face. Mode **(a)** is a bare shape-API loop (`Contains`, then alternating `DistFromOutside` /
`DistFromInside`); mode **(b)** is the real `TGeoNavigator`. Everything else in this project asks
single-shot questions from sampled points; only this composes steps, which is where navigation
actually fails.

```bash
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xray --fixtures
$B/stage/bin/o2-bench-detectorsbase-xray --self-test     # 19 checks, no DB, no oracle
```

Measured: **zero stalls of any kind over 1.15M steps** on fixtures and Bagger — the feared failure
mode (zero-length steps, ping-ponging) does not occur for the surface or CSG representations, and
mode (a) and (b) never disagree there. The mesh is where transport breaks: 6 rays enter
`cyl_inter_cyl`'s tessellation and **never leave**, and `O2Tessellated::Contains()` contradicts its
own `DistFrom*` on 534 + 15 intervals.

**Use `--beams N`, not the axis raster.** A parallel beam is *direction-poor* — three axis beams
are three directions however many rays are fired — and the known torus quartic defect (§5.5) is
invisible to them and visible to a Fibonacci fan. The benchmark's acceptance test is that it
detects that defect: at ×0.1, 96 directions × 24² rays give **LOST = 4 (in pairs — enter *and*
exit) and 2 parity mismatches**. The `parity` audit — `Contains()` at interval midpoints, the one
check independent of the stepping — is the more sensitive of the two and fires even at the
**shipped ×1 size** (2 in 262144 rays).

**Volume by chord integration** is a **1e-4 to 1e-5** instrument at practical densities and is
*non-monotone in N*, so it is quoted as an envelope at a stated density and never extrapolated. It
**cannot** resolve the 1.3e-06 capacity residuals and the tool says so before printing. What it is
for: gross errors, and **composites** — where it is the only independent volume a CSG part has.
Measured: the CSG composite and the exact surface solid agree to **2.6e-14** on `cyl_cross_cyl` and
≤1.6e-12 on all seven Bagger rams, while ROOT's Monte-Carlo `TGeoCompositeShape::Capacity()` is off
by up to 1.45e-02.

**It is the first instrument that runs on ALICE3** — 18 parts in ~32 s / 371 MB at N=48, because
*loading* dominates and stepping is nearly free. And it immediately found what nothing else can
see: **4 of 18 parts lose 418 crossings, invert the sense of 236, leave 70 transports
unterminated and contradict `Contains()` on 336 intervals — while all 18 report `reliable` and
`navigable` with zero non-manifold edges.** That is §5.4's stated caveat — identity certifies the
*topology* survived conversion, not each face's *geometry* — measured for the first time, on real
geometry, with a number. Three of the four are the `ST0923290` family. Reported, not diagnosed.

Assembly-level transport (several parts in one world, where leaking between volumes would show) is
**not** built; `Stream_J_XRay.md` §9 scopes what it needs, including a **leaking** counter that
none of the current metrics can express.

---

## 6. The MVP: Bagger as CSG + BVHSurfaceSolid + Tessellated — **done**

```bash
python3 scripts/geometry/O2_CADtoTGeo.py scripts/geometry/STEP_examples/Bagger.step \
    --csg auto --exact-surfaces auto --mesh --dump-brep --output-folder out -o geom.C
# -> tiers: CSG 7, exact surfaces 5, tessellated 1  (of 13 leaf solids)
```

| carried by | parts | evidence |
| --- | --- | --- |
| **CSG** | `BasePin`, and the six rams (`Boom`/`Stick`/`BucketCylinder` × `Inner`/`Outer`) | `dV_sym = 0` exactly; oracle 0/0/0/0 |
| **exact surfaces** | `Base`, `Boom`, `Stick`, `BucketLink1`, `BucketLink2` | oracle 0/0/0/0 |
| **tessellated** | `Bucket` (97 faces, 4 spherical + 2 toroidal) | no exact sidecar |

The six rams are `TGeoTube ∪ TGeoTube` — **unions**. Parts decline CSG for measured reasons: 2 to
7 axis clusters (Tier-3 bodies), one plane that is neither cap nor wedge, one toroidal face.

**The one prediction in this document that failed.** It said converting the rams as booleans would
retire the three capacity failures. They *are* exact now — `dV_sym = 0` is far stronger than any
capacity band — but the Bagger gate still reads 9/12, because the verdict is computed on the
**surface** representation. Those three parts are exact in the representation they now ship in,
while the number calling them failures describes one they no longer use. Closing that is a policy
decision in `runOracleGate.py`: make the verdict representation-aware. It is the natural next step
and it is not a geometry problem.

### The original scoping, kept for reference

This is the right next target, and Bagger is the right specimen: 13 prototypes, all quadric-only,
6 of them the tube-pair rams that CSG represents exactly and that are *precisely* the three parts
still failing the gate on capacity.

**What already exists:** the exact path (2), the tessellated path (3), the cascade between them,
the oracle, the gate, and the census that says what Bagger is.

**What is missing, in dependency order:**

1. **Gate any `TGeoShape`, not just `O2BVHSurfaceSolid`.** Without this, a CSG-emitted part cannot
   be scored, and "MVP" would mean "unverified". Everything else depends on it.
2. **A Bagger-sufficient recogniser**: whole-part primitive templates (tube / box / cone) plus the
   two-cluster `tube ∪ tube` case. Census-measured as enough for Bagger; deliberately *not* the
   general Tier-3 decomposition.
3. **An emitter to native ROOT shapes** — `TGeoTube`/`TGeoBBox`/`TGeoCompositeShape`. No new O2
   class (§5.2).
4. **Acceptance per part**: OCCT symmetric-difference volume ≤ model tolerance, *and* the part
   passes the ordinary oracle gate in its emitted form.
5. **Converter dispatch and policy**: try CSG → exact surfaces → tessellated, with the choice and
   its evidence recorded per part.
6. **A tiered scorecard in `gate.json`**: which representation accepted each part, and what it
   fell back to. Coverage then becomes a tiered number instead of a single fraction — the only
   honest way to report a system with three representations.

**What success looks like:** every Bagger part carried by the *strongest* representation that
accepts it, all 12 passing the oracle gate, with the six rams exact as booleans — which would also
retire the last three capacity failures rather than tuning them away.

**Explicitly not in the MVP:** free-form surfaces (Stream B), general CSG decomposition (Tier 3),
ALICE3 scale-up, the quartic fix, revolved/extruded profile solids. Each has its own brief.

---

## 7. Where to read further

| document | what it is |
| --- | --- |
| `BVHSurfaceSolid.md` | the solid's own design and format documentation |
| `SolidNavigationHarness.md` | the harness: every option, how to read the output, `perf` entry point |
| `TolerancePolicy.md` | every tolerance constant, what it governs, why it has that value |
| `CSG_Pipeline.md` | the B-rep → CSG design: prior art, the four tiers, acceptance |
| `Workstreams.md` | the six parallel streams, the contract, file ownership |
| `CodeReview_Fable.md`, `_v2.md` | the two deep reviews; the register of findings |
| `Stream_A_CSG.md`, `Stream_C_Hygiene.md`, `Stream_E_Scale.md`, `Stream_F_EdgeIdentity.md`, `Stream_G_AnyShape.md`, `Stream_H_CSGEmitter.md` | wave 0/1 and MVP investigation records |
| `NEXT.md` | the session-to-session hand-over; rewritten by whoever finishes |
