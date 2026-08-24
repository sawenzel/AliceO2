# NEXT — session-start instruction for the CAD → TGeo work

This file is the current hand-over. Whoever finishes a session should **rewrite it**.
Last rewritten 2026-08-24, at the end of the recognition-programme session (all three rungs
landed: revolved profiles, the prism family, the single-cell emitter).

Branch `swenzel/bvhsurfacesolid`. Everything below is committed unless marked otherwise.

## The next major step (queued 2026-08-23, for a FRESH session)

**The closure test**: real `o2-sim -m PIPE ITS TPC MAG` physics through the round-tripped
geometry in three representations, scored against the original C++ TGeo as its own oracle.
→ [`Handoff_ClosureTest.md`](Handoff_ClosureTest.md)

**Decided 2026-08-24 (Sandro): the flat-CSG programme runs in any case** — finish composite
recognition (torus cells, TGeoEltu, Tier-0, splitter decomposition; measured start 56/183 = 31 %
of composite-sourced parts recognised), then build `TGeoBVHCSG` as the optimisation layer for
complex composites (general DNF + cell BVH, in the spirit of Geant4's G4MultiUnion but not
union-only). → [`Handoff_FlatCSG.md`](Handoff_FlatCSG.md). Ordering vs the closure test is
Sandro's call; its R1/R2 (torus + Eltu) would lift PIPE to ~145/176 csg first.

The recognition programme that was queued before it is **DONE** (2026-08-24):
→ [`Stream_AJ_Recognition.md`](Stream_AJ_Recognition.md) is the record. Headline: every native
primitive class the writer emits is recognised at 100 % on PIPE/ITS/TPC/ABSO/TRD; 1655/1655 CSG
parts agree with their source `TGeoShape` (zero containment disagreements); the fixtures gate now
**exits 0** with `tube_window`, `cyl_inter_cyl` and `oblique_cut_cyl` shipped as csg — the two
sliver distout rays are out of every shipped verdict. The closure test's "everything CSG" variant
is now much closer to its name.

## Read this first

[`Tutorial.md`](Tutorial.md) is still the map; [`Review_2026-09.md`](Review_2026-09.md) is the
deep review with the verification appendix; [`INDEX.md`](INDEX.md) orders every document;
[`Plan_Presentation.md`](Plan_Presentation.md) is the plan of record for the WG talk;
[`Stream_AJ_Recognition.md`](Stream_AJ_Recognition.md) for everything about CSG recognition.

## Where the branch stands (recognition rows re-verified 2026-08-24; the rest 2026-08-22)

| | |
| --- | --- |
| `ctest -R BVHSurfaceSolid` / `BVHAssembly` | **113 + 22 cases**, green |
| `O2_CADtoTGeo.py --self-test` | **54 checks** = 18+8+10+12+6 — never quote the last line alone |
| `csg/emit.py --self-test` | **184** (11 acceptance + 173 recognise/emit, incl. three candidate-digest tables) |
| `checkKnownSource.py --self-test` (new) | **15/15** — the third acceptance test: emitted shape vs the source `TGeoShape` |
| `runOracleGate.py --self-test` / xray `--self-test` | 17/17 · clean |
| Bagger gate | **13/13, CSG 7 / surfaces 6 / tessellated 0, exit 0**, bit-identical through all three rungs |
| fixtures gate | **exit 0, 9/10 csg** (only `torus_union_cyl` declines); the two sliver rays remain only in the side-by-side surface columns, no shipped verdict carries one |
| detector corpora (regenerated, converted, known-source-checked) | PIPE **128**/25/23 of 176 · ITS **207**/52/0 of 259 · TPC **146**/26/0 of 172 · ABSO **26**/3/0 of 29 · TRD **1148**/22/0 of 1170 — every accepted CSG at dV_sym = 0, known-source 1655/1655 |
| **Geant integration demo** | works: IRIS + Bagger as sensitive external detectors, real hits, zero stuck tracks / nav errors (`Stream_Z_IntegrationDemo.md`) |
| material budget, exact vs tessellated | 0.039 % aggregate over 512 Fibonacci geantino rays; sole 8192-ray divergence is `BucketLink2`, the not-closed mesh |
| costs, exact vs tessellated | 1.94× transport, +42 s one-off build (IRIS CloseShape), +174 MB; exact is *smaller on disk* |
| website (`website/`) | mesh viewer, charts, JS exact raytracer bit-true vs the kernel, event-display player, self-check 30/30; `website_data/decline_reasons.json` regenerated 2026-08-24 (Bagger + ALICE3 + fixtures) |
| ray bridge (`tgeoRayService.py`) | the real kernel serving pixels: 10k rays in 4 ms |
| oTOF converts | 20 prototypes / 62 628 placements, 20/20 exact surfaces, 19/20 CSG `TGeoBBox` at dV_sym = 0 (`Stream_AC_OTOFTraversal.md`) |
| TGeo → STEP round trip | `O2_TGeoToCAD.py` self-test **105** (frozen through the recognition work); six-module study: ~7.4 M Contains samples, ONE disagreement (`Stream_AD/AF/AG/AH/AI`) |
| benchmark JSON (`website_data/`) | five hero parts complete; timing under load flagged `timingPreliminary` |

**Quote it correctly (unchanged):** `cyl_inter_cyl`'s mesh is closed as a triangle set; the
X-ray's lost rays are a **navigation** loss in `O2Tessellated`'s stepping, not holes.

## Open, in the order I would take them

1. **The closure test** (`Handoff_ClosureTest.md`) — the queued next step, fresh session.
2. **Converter: the JIT namespace bug.** `geom.C` emits `LoadSurfaceSolid` as a namespace-scoped
   forward declaration; the JIT wrapper makes it resolve wrongly and **the simulation continues
   silently without the module**. Workaround committed (`integration_demo/patch_exact_macro.py`);
   fix the converter's macro emission and make a failed module load loud.
3. **Converter: parallel `--csg auto` runs race** (two parallel conversions lost shapes):
   serialize/lock the deferred emit. And score **oTOF through the oracle gate** (converts, never
   scored).
4. **Recognition follow-ups — now folded into [`Handoff_FlatCSG.md`](Handoff_FlatCSG.md)** (details in `Stream_AJ_Recognition.md` §7): `TGeoEltu`
   (PIPE's 22-part elliptical population — needs the face reader to learn the elliptic
   carrier); run + score the **MAG corpus** (demand 14 prisms, never converted); the cell leaf
   budget (8) vs ITS's two 10-halfspace connector blocks; Xtru parameter comparison in
   `checkKnownSource.py`; a face-count cap before the cell matcher's boolean build if a corpus
   ever makes its ~12 min/1170 parts painful.
5. **Kernel-confirmed: `ST1829909_01` leaks parity through real inter-face slits** (8/15 299
   rays; `maxSharedEdgeDeviation` ≈ the declared model tolerance). Needs the per-face localiser,
   then widened crossing acceptance near identified shared edges or trim-snapping.
6. **Kernel: loose bounding boxes.** `Stick` claims 40.1 cm in z around an 11.0 cm solid (3.6×).
   Use the cover-box union in `ComputeBBox`.
7. **Writer (still frozen, now with a measured customer):** degenerate prism sections block the
   last TRD demand case (`B045cut`, a `TGeoTrd1` with `dx1 = 0` — declined as "[2, 4] distinct
   vertices"); plus the standing items: bare depth-32 chain constant, STEP-writer segfault
   bracket (38 676 fine / 74 601 crash), `TGeoPara`, `TGeoHalfSpace` (732 TPC placements),
   instruments B/C + ITS/MAG re-runs post-fix, loud-refusal path for genuine scales.
8. **Hand-written-geometry findings on the record** (unchanged, Sandro's call on upstreaming):
   PIPE bellows duplicates FIXED on fork branch (`4fe6b285e0` here); **`RB26s3Bellow` has ZERO
   daughters — missing steel**, restoring it changes the material budget; TPC prepreg z-typo
   FIXED (`3755c83277` here); ROOT `TGeoShapeAssembly` defects (see item 11).
9. **Track 3b, the real data:** MCStepLogger→`events.json` exporter, then the website's event
   tab replays the real IRIS/Bagger transports.
10. **Materials matching:** 26/55 IRIS volumes matched, ten by accidental prefix. Anchored
    part-number matching.
11. **`O2BVHAssembly` landed** (`Stream_AE_BVHAssembly.md`): flat oTOF Contains 4.3×,
    Safety(out) 28.7×, transport 6.5×; honest limit at ≤68 daughters; `FindNode` untouched
    (§8's three hook options — Sandro to pick). Two measured ROOT defects
    (`TGeoShapeAssembly::DistFromOutside` Big() outside the bbox of a voxelized assembly;
    `Safety` exceeding the true minimum) — upstream-report decision is Sandro's.
12. **The talk (Track 4):** assemble from Review + Stream_Z + Stream_AJ + the website; re-run
    `timingPreliminary` numbers on a quiet box; single-file website bundle.
13. **The face-normal gate column** and **the `auto`-mode unreliable-shipping policy** — both
    pre-corpus items, now five hand-overs old.
14. **Standing, unchanged:** free-form surfaces; models-are-not-legal overlaps and the broken
    `CheckOverlaps` on our shape; `Curve2D::closestPoint` as the kernel hot spot; mesh healing;
    composite trees (Tier 3 proper) remain deliberately not built.

## Traps in the environment (all still current)

- Env `O2/latest-swenzel-bvhsurfacesolid-o2`, build dir
  `B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2` (stack of 2026-08-22).
- **pythonOCC needs `--no-system SWIG`** when rebuilding (system 4.2.0 fails the required range).
- **Keep the converter env and the sim env separate.** The pythonOCC `PYTHONPATH` prepends make
  `o2-sim` segfault at startup; the converter needs O2/ROOT *in addition* to OCC — a bare OCC
  shell silently loses CSG (the deferred-emit WARN now says so per part).
- External-detector hits live in `o2sim.root` (`IRISHit`, `BAGRHit`) under `o2-sim-serial`.
- Testing a rebuilt detector library needs `export O2_ROOT=$B/stage`; use `o2-sim-serial` from
  `$B/stage/bin`.
- Reconfiguring CMake needs Clang on the prefix path
  (`export CMAKE_PREFIX_PATH=$HOME/alisw/sw/ubuntu2404_aarch64/Clang/v20.1.7-local1:$CMAKE_PREFIX_PATH`);
  new `.C` macros under scripts/ must be listed in `O2RootMacroExclusionList.cmake`.
- `rm` is blocked by a repo hook (`.claude/hooks/deny-deletions.py`); move files aside instead.
- Everything else in the 2026-08-09 list still applies: eval+command in one shell; stage lib
  first; prepend-never-replace; one ninja at a time; detach long runs to unique `--out` paths;
  never write into `STEP_examples/` or `ALICE_3_example/`; convert ALICE3 without `--mesh`;
  `manifest.json` stores absolute paths; run `--csg auto` conversions strictly serially.

## Commands

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ctest -R 'BVHSurfaceSolid|BVHAssembly'

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step

# converter env additions (on top of the O2 env above):
SW=$HOME/alisw/sw/ubuntu2404_aarch64
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
python3 scripts/geometry/csg/emit.py --self-test
# detector corpus with ground truth (sim env first, then converter env; strictly serial):
#   $B/stage/bin/o2-sim-serial -n 0 -g boxgen -m <MOD>   (with export O2_ROOT=$B/stage)
#   python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root <MOD>.step --report <MOD>_writer_report.json
#   python3 scripts/geometry/O2_CADtoTGeo.py <MOD>.step -o geom.C --exact-surfaces auto --csg auto
#   python3 scripts/geometry/checkKnownSource.py --original o2sim_geometry.root \
#       --writer-report <MOD>_writer_report.json --converted .

# the website (serve locally, then open the printed URL)
cd scripts/geometry/website && ./fetch_testdata.sh <gate-workdir> && python3 -m http.server 8231
# the ray bridge (inside the O2 env), for the website's RemoteEngine / engine-diff view
python3 scripts/geometry/tgeoRayService.py --port 8077
# the integration demo, end to end
scripts/geometry/integration_demo/   # see its README / convert_all.sh, run scripts
```
