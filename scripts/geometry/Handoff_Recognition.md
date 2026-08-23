# Handoff — the recognition programme: prisms, revolved profiles, then cells

**Written 2026-08-23 for a fresh session to execute, BEFORE the closure test**
([`Handoff_ClosureTest.md`](Handoff_ClosureTest.md) waits on this). Sandro's decision, from the
discussion: the round-trip studies showed that most parts declining CSG are not booleans at all
but single TGeo primitives our recogniser has no template for; building those templates raises
the closure test's "everything CSG" variant from a cascade toward its name.

Read first: [`Tutorial.md`](Tutorial.md) (the map), [`NEXT.md`](NEXT.md) (state + traps),
`Stream_AA_FlatCSG.md` (the measured ordering this follows), `CSG_Pipeline.md` (tiers and
acceptance design), `Stream_K_Tier0.md` (**the vacuous-cone lesson — read before writing any
scoring code**), and the round-trip studies `Stream_AD/AF/AG/AH/AI` (the demand tables and the
requirement list below).

## 1. Why this is unusually well-prepared work

Recognition work is normally starved of ground truth. Here every input carries its answer:
`O2_TGeoToCAD.py` (the branch's own TGeo→STEP writer, self-test 105) emits a report naming the
**source TGeo class and parameters of every part it writes**. Six modules round-trip cleanly
(PIPE, ITS, TPC, TRD, MAG, ABSO — ~7.4 M `Contains` samples, one disagreement), so hundreds of
(B-rep, known source shape) pairs exist and are regenerable on demand. Recognition development
here is genuine TDD against known-right answers.

## 2. The rungs, with measured demand

**Rung 1 — the revolved-profile recogniser → `TGeoPcon` (and cone/stack cases).**
Cluster coaxial caps and laterals, reconstruct the (r, z) profile, emit the native class.
Demand (decline counts): PIPE **46**, TPC **~19**, ITS **18**, ABSO **16**.
Requirements from the studies: mixed cone+cylinder laterals on one axis (`IBCYSSCone`);
z-steps (duplicate z planes — legal, and the writer now builds them validly); `rmin` stepping
through 0; half-turn and partial-phi wedges; Pgon-vs-Pcon disambiguation (Pgon laterals are
planes at the **apothem** radius — `conv_pgon` documents the convention).

**Rung 2 — the prism-family recogniser → `TGeoTrd1`/`Trd2`/`Arb8`/`Xtru`/`Pgon`.**
Planar-face graph → parameters. Demand: TPC **~44**, ITS **47** (Xtru-heavy: 23), TRD **149**
real decliners, MAG **14**. Requirements: slanted opposite faces (the current box rule's
"no opposite partner"); degenerate sections (a `Trd1` with `dx1 = 0`, TRD's `cutOnB045`);
twisted `Arb8` (ruled faces — decide and document whether in scope or declined with reason).

**Rung 3 — gated, only if 1+2 land with margin: the single-cell emitter** for true composites
(`Stream_AA` step 2): retires `tube_window` and `cyl_inter_cyl` exactly — and with them the two
known sliver rays of the fixtures gate. Composite *trees* (Tier 3 proper) stay out of scope.

All emission targets are **native ROOT classes — no new O2 class**; a recognised placed
primitive travels with its transform per `Stream_N_PlacedPrimitives.md`.

## 3. The acceptance, non-negotiable

Every candidate passes **three** independent tests, or it is declined with a recorded reason:
1. **OCCT symmetric difference** (`csg/accept.py`): `dV_sym` ≤ model tolerance, cut both ways.
2. **The oracle gate** on the emitted shape (`runOracleGate.py` scores any `TGeoShape`).
3. **NEW — the known-source check**, wherever the corpus provides it: the emitted shape against
   the writer-report's original TGeo shape (class + parameters where directly comparable;
   otherwise Capacity to 1e-9 and a seeded `Contains` cross-check against the original
   `TGeoShape` — the round-trip studies' instrument A, reversed).

**Scoring rule (the Stream_K lesson, it has already burned this project once): one measured
quantity decides for every candidate — achieved geometric gap relative to patch/solid scale —
never per-class criteria, never an angle, and every new recogniser ships with negative controls
in the self-test (a near-miss that MUST decline).**

## 4. Regression floors (before → after, every rung)

- `csg/emit.py --self-test` 33/33; `O2_CADtoTGeo.py --self-test` **54** (48 + the 6-check
  multibody suite); `O2_TGeoToCAD.py --self-test` **105**; `runOracleGate.py --self-test` 17/17;
  ctest `BVHSurfaceSolid` 113 + `BVHAssembly` 22.
- Bagger gate 13/13 exit 0 with CSG 7 unchanged **bit-identical minus timing**.
- Fixtures gate: the two sliver `distout` rays are the known state (Review Appendix A); rung 3,
  if reached, is what removes them — nothing else may move.
- The already-recognised sets stay recognised identically: ITS 128, TPC 82 (post-WSEG-fix
  baseline, `Stream_AI` addendum). Re-run the ITS and TPC reverse conversions and diff the
  decline histograms: only the *intended* classes may move.
- Update `website_data/decline_reasons.json` (via `csg/decline_catalogue.py`) and note the new
  coverage in the report — the website shows these reasons per part.

## 5. Corpora and regeneration (scratch is volatile — regenerate, do not hunt old dirs)

- Fixtures: `runOracleGate.py --workdir <w> --fixtures` (or `make_boolean_fixtures.py`).
- Bagger / ALICE3 / oTOF STEPs: in `scripts/geometry/STEP_examples/` and
  `scripts/geometry/ALICE_3_example/` (read-only — NEVER write there).
- Detector corpora with ground truth: sim env →
  `$B/stage/bin/o2-sim-serial -n 0 -g boxgen -m <MOD>` (with `export O2_ROOT=$B/stage`), then
  converter env → `O2_TGeoToCAD.py o2sim_geometry.root out.step --report report.json`, then
  `O2_CADtoTGeo.py out.step --exact-surfaces auto --csg auto` (no `--mesh`). PIPE, TPC, ITS
  first (they carry the demand); TRD if time (big).

## 6. Environment (verbatim; eval+command in ONE shell invocation)

```
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
# converter (needs ROOT *and* OCC — a bare OCC shell silently yields CSG 0):
SW=$HOME/alisw/sw/ubuntu2404_aarch64
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
# geometry regeneration: export O2_ROOT=$B/stage ; use o2-sim-serial (NEVER parallel o2-sim)
```
Standing traps (all in `NEXT.md`): two parallel `--csg auto` conversions race — run conversions
serially; detach >2 min runs to unique logs; `rm` is hook-blocked (use `mv`); one `ninja` at a
time (none should be needed); `manifest.json` stores absolute paths.

## 7. Files owned, and the staging

Owns: `scripts/geometry/csg/**` (recogniser/emitter/acceptance), the CSG hook surface of
`O2_CADtoTGeo.py`, new fixture generators if needed, `Stream_AJ_Recognition.md` (the report).
Does NOT touch: the kernel (`Detectors/Base/**`), `O2_TGeoToCAD.py` (writer is frozen ground
truth for this task), the website, the shared docs (`NEXT.md`/`Tutorial.md`/`INDEX.md`/
roadmaps — the integrating session folds).

Stages, commit per stage in the house style (title states what is now true; trailer exactly
`Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`):
1. Baselines: regenerate PIPE/TPC/ITS corpora, record the decline histograms as the before.
2. Rung 1 with self-tests + negative controls; re-run corpora; coverage table.
3. Rung 2 likewise.
4. Rung 3 only if 1+2 are green with time to spare, same discipline.
5. `Stream_AJ_Recognition.md`: per-corpus before → after tier tables, per-class recognition
   rates, every surviving decline with its reason, dead premises recorded.

**Definition of done**: the three-test acceptance holds for every new recognition; the demand
numbers of §2 measurably convert (state exactly how many of each); every floor in §4 holds; a
part the recogniser is UNSURE about still declines loudly — coverage never buys correctness.
