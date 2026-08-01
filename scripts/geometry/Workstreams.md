# Parallel workstreams — the CAD → TGeo programme, decomposed for independent execution

Date: 2026-08-01. Author: Claude (Fable 5). Companion to
[`CodeReview_Fable_v2.md`](CodeReview_Fable_v2.md) (what is wrong and why) and
[`CSG_Pipeline.md`](CSG_Pipeline.md) (the B-rep → CSG direction in detail).

Purpose: cut the remaining work into streams that can be executed **concurrently by separate
agents or sessions** without stepping on each other, and state for each one its goal, the files it
owns, its steps, its gate, and what it hands back. Section 1 is the contract every stream obeys;
Section 2 is the conflict map; Section 3 is the launch order; Sections A-F are the briefs.

**The strategic frame, in one paragraph.** There are three levers on CAD → TGeo coverage and
quality and they are complementary, not competing. (i) *Simplification* — recognise primitives,
booleans, extrusions and revolutions and emit exact simple solids; this is where mechanical CAD,
Bagger, oTOF and the beampipe live, and it removes the trim problem rather than solving it.
(ii) *Free-form surface support* — the measured coverage lever on ALICE3, where 36 of 55 solids and
2377 of 9266 faces are B-spline or NURBS, and where neither CSG nor better trims can help.
(iii) *Kernel and gate quality* — the residual defects and the ability to measure at scale. Going in
all three directions at once is the right call, and they decompose cleanly because they touch
mostly disjoint code.

---

## 1. The contract — every stream, no exceptions

**The invariant to defend.** Over both corpora (9 converted fixtures, 12 Bagger parts):
**zero disagreements with the OpenCascade oracle outside tolerance, on every column.** A change that
improves a gate *total* while breaking this is a regression. Gate totals and disagreement counts are
separate numbers; never quote one without the other.

**Diff columns, not verdicts.** A bit-identical numeric column can accompany a real change, and has:
the K4 sampler fix moved two fixtures from closed to open while every checksum stayed identical, and
reading the totals alone would have reverted the session's main result. Compare
`<workdir>/gate.json` with the timing fields stripped (`timing*`, `*Seconds`, including
`closeShapeSecondsMesh` / `closeShapeSecondsSurface`); everything else is deterministic and must
match when a change is meant to be inert.

**Evidence before theory.** This project's last four advances came from building the instrument that
names the object and then reading it — not from acting on the most plausible hypothesis. Two of the
last three hand-overs recorded premises that died on contact with a measurement. Measure first, and
report a dead premise as a result.

**Before and after, every time.** `ctest -R BVHSurfaceSolid` (61 cases at branch point) plus both
oracle gates (fixtures 8/9, Bagger 6/12 at branch point) before and after each committable step.
State which columns moved and why, or state that none did and why that is the expected outcome.

**Doc ownership is single-writer.** `NEXT.md`, `CodeReview_Fable.md`, `CodeReview_Fable_v2.md` and
`TolerancePolicy.md` are **not** to be edited by a stream. Each stream writes its own
`Stream_<X>_<topic>.md` and the integrating session folds them and rewrites `NEXT.md`. This is the
single most important rule for concurrent work here, because those four documents are the project's
memory and a merge conflict in them loses reasoning, not just text.

**Environment.** (Traps that have each cost time before.)

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH   # stage FIRST, or stale libs win
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
ctest -R BVHSurfaceSolid

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
```

- The gate needs the O2 environment **in the same shell**; `runOracleGate.py` sets the OCC python
  but not `LD_LIBRARY_PATH` for the harness.
- pythonOCC needs the alibuild **python3.10**, not the system 3.12. OCCT/pythonOCC live under
  `$ALIBUILD_WORK_DIR/ubuntu2404_aarch64/{OCCT,pythonOCC}/latest`.
- Test binary is at `$B/stage/tests/`, harness at `$B/stage/bin/`.
- `--skip-convert` reuses `<workdir>/db` — valid only when the *converter* did not change. Copy a
  finished workdir and re-run with it to re-measure a kernel change.
- The gate exits non-zero whenever any part fails, which is the normal state. Read the counts.
- **Never** write generated artifacts into `scripts/geometry/STEP_examples/`. Use a scratch folder.
- Quick kernel probe: standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase -lGeom`, plus
  `-I$HOME/alisw/O2/Detectors/Base/src` to reach the kernel header directly
  (`TolerancePolicy.md` §6).

---

## 2. Conflict map

The bottleneck is `Detectors/Base/src/BoundedSurface.h` (5076 lines), which four streams want.

| file | A CSG | B free-form | C hygiene | D profiles | E scale/gate | F edge-ID |
| --- | --- | --- | --- | --- | --- | --- |
| `BoundedSurface.h` | — | new class (append) | `contourIntegralAlongCurve` | new classes (append) | — | `measureRimClosure`, `ClosureReport` |
| `O2BVHSurfaceSolid.{h,cxx}` | — | `Add*Surface` | `containsByVote` | `Add*Surface` | — | report plumbing |
| `O2SurfaceSolidIO.cxx` | — | new record type | `edgeEndpoints` | new record type | — | edge-ID parse |
| `O2SolidHarness.*`, `runSolidHarness.cxx` | — | — | — | — | **owns** | — |
| `testBVHSurfaceSolid.cxx` | — | append | append | append | — | append |
| `O2_CADtoTGeo.py` | one hook | extractor | — | emitter | — | sidecar writer |
| `occtOracle.py`, `runOracleGate.py`, `make_boolean_fixtures.py` | reads | — | — | — | **owns** | — |
| `scripts/geometry/csg/**` (new) | **owns** | — | — | — | — | — |

Rules that make this work:

1. **New surface/solid classes are appended at the documented insertion point**, never interleaved.
   B and D therefore conflict only in whitespace.
2. **Tests are appended in a delimited block** ending with a marker comment naming the stream.
   Conflicts become trivial.
3. **C lands first and alone** — it is a few hours, it closes the last failing numeric column, and
   everyone else rebases onto it rather than around it.
4. Each agent works in its **own git worktree** and rebases before handing back.
5. *Optional enabler, not a prerequisite:* splitting `BoundedSurface.h` into
   `curves/`, `surfaces/`, `trim/`, `closure/` headers would remove the bottleneck entirely and was
   already recommended in v1 §5 ("~⅓ duplication across the five surface classes"). If it is done,
   it must be **one mechanical commit, alone, verified bit-identical** on both gates and ctest —
   never mixed with behaviour.

---

## 3. Launch order

| wave | streams | why together |
| --- | --- | --- |
| **0** | **C** (+ the `scripts/geometry/` cleanup) | Serial, half a day. Closes the last failing numeric column and clears a directory that five stale script variants make unreadable. Everyone rebases on it. |
| **1** | **A**, **E**, **F** — 3 agents | Near-zero overlap: A is a new Python package, E owns the harness and gate scripts, F owns the sidecar/closure path. All three produce *measurements* that waves 2-3 need. |
| **2** | **B**, **D** — 2 agents | Both add new kernel classes; developed in parallel, landed sequentially. Start only after wave 1's census (A step 1) and scale sweep (E step 1) are in, because both can change what is worth building. |
| **3** | Tier 3 / `O2CSGSolid` | Gated on A step 7's cell-count table. Do not pre-commit to it. |

---

## Stream A — the CSG / simplification pipeline

**Goal.** Recognise and emit exact simple solids: primitives, booleans of primitives, and
canonicalised swept surfaces. Removes the trim problem for the parts it covers rather than solving
it.

**Full brief:** [`CSG_Pipeline.md`](CSG_Pipeline.md), whose §8 is the step list. Summary: census →
acceptance harness (OCCT symmetric-difference volume) → Tier 0 swept canonicalisation → Tier 1
primitives → composite timing probe → Tier 2 booleans → Tier 3 spike.

**Owns.** New `scripts/geometry/csg/` package; one integration hook in `O2_CADtoTGeo.py`;
`Stream_A_CSG.md`.

**Gate.** Per tier: every eligible ladder fixture recognised; symmetric-difference volume ≤ model
tolerance; `runOracleGate.py` passes on the emitted shape; census table updated with the achieved
rate. Tier 2's headline: the six Bagger cylinder parts and `tube_window` convert as booleans →
**G3 becomes 12/12 and the tube-tube problem is finished.**

**Hands back.** The census table (which is what waves 2-3 are scoped on), the achieved recognition
rate per model, and the Tier-3 go/no-go table.

**Do not** touch the trim machinery, the closure check, or the kernel.

---

## Stream B — free-form surface support

**Goal.** The measured coverage lever on ALICE3: 36 of 55 solids and 2377 of 9266 faces are
B-spline or rational NURBS surfaces. Neither CSG nor better trims can reach them.

**Why now, against v1's ordering.** v1 deprioritised the B-spline-surface milestone as "largest
effort, do last", on the working assumption that trims were the blocker. This session's measurement
of `CAD_noETA.stp` (see `CSG_Pipeline.md` §6) says otherwise for this corpus: the 15 quadric-only
solids are essentially the 15 the converter already handles, so **everything else on the board adds
zero ALICE3 coverage.** The deprioritisation should be revisited on that evidence — but revisited,
not reversed by assumption: step 1 below is the check.

**Steps.**

1. **Scope it honestly first.** For the 36 free-form ALICE3 solids: how many faces each, are the
   free-form faces *trimmed* by analytic curves, are they rational, what degrees and pole counts,
   and how many are the "hybrid band" (mostly quadric with a few free-form faces — `#33027` is 996
   faces with 69 free-form, `#4411` is 245 with 40, `#123` is 50 with 4). The hybrid band may be
   reachable far more cheaply than the general case. *Gate:* the table exists.
2. Kernel: a `BSplineBoundedSurface` — ray intersection (Newton with interval/Bézier-clipping
   bracketing, or subdivision-BVH over Bézier patches), normal, `Safety` bound, conservative AABB,
   capacity contribution. This is the piece v1 called the largest effort and it is.
3. Sidecar record type + IO + converter extractor.
4. Honest exactness accounting: it is *not* exact and must say so — `capacityIsExact()`,
   `NavigationReliability`, and the tolerance it does achieve, stated per solid.

**Owns.** New surface class in `BoundedSurface.h` (appended), its IO record, its converter
extractor, `Stream_B_FreeForm.md`.

**Gate.** New ladder fixtures with a genuine free-form face; oracle-clean within a *declared*
tolerance; ALICE3 coverage re-measured and reported as a tiered number, not a fraction.

---

## Stream C — kernel hygiene (wave 0, serial)

**Goal.** Close the last failing numeric column and two latent defects. Small, sharp, first.

**Steps.**

1. **N1 — the capacity residual. Build the instrument; do not carry a hypothesis in.**
   `CodeReview_Fable_v2.md` §3/N1 records this being measured on 2026-08-01 and **four hypotheses
   dying**, including that document's own headline one: the quadrature is converged (a 64× finer
   `kContourMaxSpanU` moves nothing beyond the last digit), the chart cut segments are present and
   the contour closes to 7e-12, the `deltaV == 0` bridge skip is exactly correct, and the
   `Contains` Monte-Carlo witness is 3× too coarse to separate the 1.1e-2 cm³ in question. What
   survives is a perfect predictor — capacity fails **iff** a trim wire carries a B-spline with more
   than one interior knot span, 21 of 21 parts across both corpora — and no mechanism.

   The one instrument not yet built is **per-face surface area**: a different functional of the same
   region, free from OCCT per face (`BRepGProp.SurfaceProperties`), and computable on our side by
   `integrateOverCurveTrimByParts` with `F = r·u` — same contour, same wires, same quadrature,
   different integrand. Areas agreeing while volume disagrees indicts the integrand or the surface
   frame; areas disagreeing indicts the trim region and measures the mis-trim in cm². Add it as a
   `--face-areas` harness flag and an oracle column (the project's own rule: a question asked twice
   becomes a harness flag). *Gate:* the table exists and says which side is wrong. Only then fix.

   Independently of that, two small defects in `contourIntegralAlongCurve` stand on their own
   merits: `spanU` is taken from the endpoints only, so it reports **zero** for a closed curve and
   gives one interval for a whole loop; and the bridge `Curve2D` is constructed per seam per call.
   *Gate for those:* every column bit-identical.
2. **N2 — the loader's B-spline endpoints.** `O2SurfaceSolidIO.cxx:149-157` reads endpoints off the
   first/last poles, assuming clamping, while the kernel evaluates them (K1's fix). Nothing
   validates clamping. Use `curve.pointAt(0.)`/`pointAt(1.)`, and add clamping to `Curve2D::valid()`
   as a *reported* property, not a rejection, so the two layers cannot diverge silently again.
   *Gate:* a test with an unclamped input where loader and kernel previously disagreed.
3. **N4 — ambiguity-aware voting.** `containsByVote` calls `parityAlong` without the ambiguity
   pointer, so a point can be decided by a majority of tie-breaks. Collect `(answer, ambiguous)` per
   direction, decide among unambiguous shots when there are any, and report when there are none.
   *Gate:* measure the incidence on both corpora — expected ≈ 0, and the point is to have measured
   it rather than assumed it.
4. **Directory cleanup** (v1 §7, still open): five stale `O2_TGeoToCAD*` variants, four
   `O2_CADtoTGeo_with_*`, `*.py~`, `#...#`, `._*`, `a.out`, `tgeo_to_cad.o`, `__pycache__`. One
   commit, `git rm` plus `.gitignore`. The TGeo→CAD scripts matter again for Stream E and it is
   currently impossible to tell which of the five is live.

**Owns.** `contourIntegralAlongCurve`, `edgeEndpoints`, `containsByVote`, their tests, the cleanup,
`Stream_C_Hygiene.md`.

---

## Stream D — exact profile solids: revolved and extruded

**Goal.** Two small exact solids that cover a large share of mechanical CAD and, specifically, the
beampipe.

**Why they are not covered by anything else.** `TGeoPcon` and `TGeoXtru` carry **polygonal** profiles
only. Every fillet, bellows knuckle and rounded slot is an arc, and polygonising it is exactly the
exactness loss this project exists to avoid. `O2RevolvedSolid` (a closed `(r,z)` profile of lines
and arcs, swept through a φ range) and `O2ExtrudedSolid` (a closed `(x,y)` profile, extruded along
a slab) fix that with a small amount of new code, because **the 2D machinery already exists and is
validated**: `Curve2D` carries lines and arcs, `CurveWire` does winding and Green's-theorem areas,
and ray intersection reduces to a 2D problem per segment whose roots are the plane/cone/sphere/torus
solvers already in the file.

Both are watertight by construction (a profile loop that closes in 2D closes in 3D), both have
closed-form capacity (Pappus / Green), both are trivially AOT-codegen-able, and `O2RevolvedSolid`
makes most of the beampipe exact in ~10 parameters per section.

**Steps.** Detection in the converter (all carriers coaxial ⇒ revolution; all carriers ∥/⊥ to one
direction with two cap planes ⇒ extrusion) → profile construction by intersecting the solid with a
half-plane / a plane → kernel class → sidecar record → tests.

**Owns.** Two new classes in `BoundedSurface.h` (appended), their IO records, their converter
detectors/emitters, new ladder fixtures, `Stream_D_Profiles.md`.

**Gate.** New fixtures `pcon_with_fillets` and `extruded_bracket_with_holes`; exact capacity in
closed form; oracle-clean at machine precision; `tube_window` (the one remaining failing ladder
fixture) exact.

---

## Stream E — scale, robustness and the gate at ALICE3 size

**Goal.** Make it possible to measure a 55-solid, 9266-face model, and find out what breaks away
from the origin. This stream produces the evidence waves 2-3 are steered by.

**Steps, in priority order.**

1. **The translate/scale sweep — do this first, it can re-qualify recorded results.** Every
   measurement on this branch was taken on parts a few centimetres across near the origin. Run the
   whole G1 ladder translated by (0, 0, 400) and scaled ×10 and ×0.1. **Every column must be
   unchanged.** If it is not, the absolute constants are the cause — `kBVHBoxTolerance` 1e-3 cm,
   `kWireJoinTolerance` 1e-5 cm, `kRayTolerance` 1e-9, `sameIntersection`'s relative 1e-7·|t| — and
   every "0 disagreements" number on the branch becomes a statement about geometry near the origin
   only, to be re-qualified in v1 §20 and `NEXT.md`. This is v1's S12, still open, and it is the
   single cheapest experiment with the largest potential consequence.
2. **The ALICE3 corpus gate.** `CAD_noETA.stp`: 55 solids, 9266 faces. Parallelise `runOracleGate.py`
   per part; add a `--sample-budget` so a large model scores in bounded time; record per-part
   wall-clock and resident memory in `gate.json` so scaling laws accumulate for free. Establish
   `CloseShape` cost vs face count on the 2034-, 1052- and 996-face solids before assuming it scales.
3. **G4 — the `auto`-mode fallback policy.** An exact part whose `NavigationReliability` is not
   `Reliable` currently ships and only warns. This is the one live production risk on the branch and
   it is ~20 lines. *Gate:* a test that a deliberately-unreliable part converts to tessellated under
   `auto` and throws under `required`.
4. **The beampipe round-trip.** O2 already contains a beampipe in `Detectors/Passive/src/Pipe.cxx`
   and `PipeRun4.cxx`, hand-written from `TGeoTube`, `TGeoPcon`, `TGeoTorus` and
   `TGeoCompositeShape`. Round-trip it TGeo → STEP → shape and compare against the **original TGeo**.
   That gives three things nothing else does: an oracle that is *exact* rather than tolerant; a
   known-right answer to score Stream A's recognition against (we know the engineer drew a tube, a
   pcon, a torus); and a real specimen at metre coordinates with 0.08 cm features — a 10⁴ dynamic
   range against those absolute constants. The `O2_TGeoToCAD*.py` scripts exist but are stale (see
   Stream C step 4); check first that one can round-trip a single `TGeoCompositeShape` faithfully
   before budgeting the rest.
5. **Ladder extensions** for what the other streams need: `pcon_with_fillets`,
   `extruded_bracket_with_holes`, `tube_with_bolt_circle`, and one deliberately **non**-CSG,
   **non**-analytic part to pin the fallback path.
6. **A tiered coverage scorecard** in the gate output: per part, which representation accepted it
   (primitive / boolean / profile / exact surfaces / free-form / tessellated fallback), the
   acceptance evidence, and what it fell back to. Coverage becomes a measured tiered number instead
   of a single fraction — which is the honest way to report a programme with five representations.

**Owns.** `occtOracle.py`, `runOracleGate.py`, `make_boolean_fixtures.py`, `O2SolidHarness.*`,
`runSolidHarness.cxx`, the revived TGeo→CAD script, `Stream_E_Scale.md`.

---

## Stream F — sidecar v3: edge identity, and closure by topology

**Goal.** Stop deciding a topological question with a proximity query. This is
`CodeReview_Fable_v2.md` §3/N3.

**The problem.** `measureRimClosure` decides whether two faces share an edge by searching for a
nearby chord and comparing against `rimEpsilon + own sagitta + partner sagitta`. `NEXT.md` has
already reached the end of that road: "the sagitta bounds each polyline against *its own* curve and
says nothing about the other face's, which is the term that dominates… it works at the shipped
tolerance because the two happen to be within an order of magnitude — a coincidence of scale, not a
criterion." It is not fixable by choosing a better band, because the quantity thresholded is not the
quantity decided: tightening `kBSplineFlatness` improves the geometry and makes the *verdict* worse
(§13.8: isolation 4.8e-5 → 2.6e-6 while open length went 0.325 → 0.984 cm).

**The information exists and is discarded.** The converter walks `TopoDS_Edge`s; one
`TopExp.MapShapesAndAncestors` call gives both faces of every edge; the sidecar writes
`{curve type, params}` and no identity.

**Steps.** Converter emits a per-model edge table and a `(edgeId, sense)` per trim curve → sidecar
v3 → loader carries them → `validateClosure` decides on **identity** (every edge ID appears exactly
twice, opposite sense: zero tolerance, zero sampling, zero band) and *reports*
`maxSharedEdgeDeviation` in cm as a measured number rather than a verdict.

**Gate.** On every currently-`Reliable` part the verdict is unchanged and the new deviation number
is below the old match band. On the six Bagger cylinders the verdict becomes closed-by-topology with
a stated deviation — the first defensible "how far apart are the two faces" number the project has
had. `CloseShape` stops depending on `kBSplineFlatness`, so tightening the flattening becomes a pure
improvement again. Cost: 5 bytes per trim curve, nothing at run time.

**Owns.** `write_surfaces_bin` and the edge-table extraction in `O2_CADtoTGeo.py`, the v3 parse in
`O2SurfaceSolidIO.cxx`, `measureRimClosure` / `ClosureReport` / `SurfaceRim` in `BoundedSurface.h`,
`Stream_F_EdgeIdentity.md`.

**Note.** Worth doing whatever is decided about CSG: it is what makes any watertightness claim
checkable rather than tuned, and the same edge table is the input Stream A's Tier-3 decomposition
wants for its concavity test.

---

## 4. What is deliberately not a stream

- **Phase 2, adjacency-based exact trims** (v1 §4.2, §10). Superseded — see
  `CodeReview_Fable_v2.md` §5.7. Every case it was scoped for is a case Stream A covers with less
  new kernel code and a stronger acceptance test. Revisit only if Streams A and D leave an
  analytic-analytic residual that matters.
- **Phase 3, shared-edge sampling for free-form edges.** A coverage item, not a correctness one, and
  Stream B changes what it would apply to.
- **Further tuning of the rim match band.** Stream F makes it moot.
- **`Safety` BVH and the rest of G5.** After correctness, and after Stream A has changed what is
  being optimised.
- **Silent unreliable shipping.** Not deferred — it is Stream E step 3, and it is the only live
  production risk on the branch.
