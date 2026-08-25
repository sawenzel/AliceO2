# NEXT — session-start instruction for the CAD → TGeo work

This file is the current hand-over. Whoever finishes a session should **rewrite it**.
Last rewritten 2026-08-25, at the close of the flat-CSG programme's rung R5 (`o2::base::O2FlatCSG`
built, shipped, measured and recorded).

Branch `swenzel/bvhsurfacesolid`. Everything below is committed unless marked otherwise.

## Where the work stands, in one paragraph

The converter recognises **163 of 188 composite-sourced detector parts (87 %)** and every native
primitive class the writer emits. Parts whose cells would make a `TGeoCompositeShape` wider than
64 leaves now ship as **`o2::base::O2FlatCSG`** — a union of intersection-cells over signed
implicit halfspaces with a BVH over sub-boxes of the cells — and that representation is measured
to be 3–34× faster on `Contains`, 4.5–30× on `Safety` and 0.83–7.1× on transport against the same
part emitted as a plain composite. The two things still queued are the **closure test** (real
physics through the round trip) and **R6, the breadth census** over the remaining Run 3 modules.

## The two queued next steps, both for a FRESH session

1. **R6 — the breadth census.** Convert and score MAG, TOF, MFT, MCH, MID, FT0, FV0, ZDC, EMC,
   PHS, CPV, HMP — one command each now. It is what turns "87 % of composite-sourced parts" and
   the flat-CSG crossover from statements about five detectors into statements about the geometry.
   → [`Handoff_FlatCSG.md`](Handoff_FlatCSG.md) §2 (R6).
2. **The closure test** — real `o2-sim -m PIPE ITS TPC MAG` physics through the round-tripped
   geometry in three representations, scored against the original C++ TGeo as its own oracle.
   → [`Handoff_ClosureTest.md`](Handoff_ClosureTest.md). Unaffected by R5: it scores the
   representations that exist, and there is now one more of them.

**The flat-CSG programme's R1–R5 are DONE** (2026-08-24/25).
→ [`Stream_AK_FlatCSG.md`](Stream_AK_FlatCSG.md) is the record and is where every R5 number lives;
[`Design_FlatCSGSolid.md`](Design_FlatCSGSolid.md) is the spec of the class;
[`Handoff_FlatCSG.md`](Handoff_FlatCSG.md) §1.1 is the rung ledger.
Before it, the recognition programme is DONE (2026-08-24) →
[`Stream_AJ_Recognition.md`](Stream_AJ_Recognition.md).

## Read this first

[`Tutorial.md`](Tutorial.md) is still the map; [`Review_2026-09.md`](Review_2026-09.md) is the
deep review with the verification appendix; [`INDEX.md`](INDEX.md) orders every document;
[`Plan_Presentation.md`](Plan_Presentation.md) is the plan of record for the WG talk;
[`Stream_AJ_Recognition.md`](Stream_AJ_Recognition.md) for CSG recognition and
[`Stream_AK_FlatCSG.md`](Stream_AK_FlatCSG.md) for the flat halfspace solid.

## Where the branch stands (flat-CSG rows re-verified 2026-08-25; the rest 2026-08-22)

| | |
| --- | --- |
| `ctest -R 'FlatCSG\|BVHSurfaceSolid\|BVHAssembly'` | **31 + 113 + 22 cases**, green — `FlatCSG` is new in R5 and is part of the floor |
| `csg/emit.py --self-test` | **353** (16 acceptance + 337 recognise/emit, incl. six candidate-digest tables — all six must stay green) |
| `O2_CADtoTGeo.py --self-test` | **54 checks** = 18+8+10+12+6 — never quote the last line alone |
| `checkKnownSource.py --self-test` | **17/17** — the third acceptance test: emitted shape vs the source `TGeoShape` |
| `O2_TGeoToCAD.py --self-test` | **107** |
| `runOracleGate.py --self-test` / xray `--self-test` | 17/17 · clean |
| Bagger gate | **13/13, CSG 7 / surfaces 6 / tessellated 0, exit 0**, bit-identical through every rung |
| fixtures gate | **exit 0, 10/10 csg**; the two sliver rays remain only in the side-by-side surface columns |
| detector corpora (csg / surface / mesh) | PIPE **166**/9/1 of 176 · ITS **255**/9/0 of 264 · TPC **165**/7/0 of 172 · ABSO **29**/0/0 · TRD **1170**/0/0 (100 %) — every accepted CSG at `dV_sym = 0`, known-source **1785/1785** |
| composite-sourced recognition | **163/188 (87 %)** |
| flat-CSG parts shipped | **10**, 12–47 cells / 72–282 halfspaces; measured crossover **45 composite leaves** against a routing threshold of 64 |
| **Geant integration demo** | works: IRIS + Bagger as sensitive external detectors, real hits, zero stuck tracks / nav errors (`Stream_Z_IntegrationDemo.md`) |
| material budget, exact vs tessellated | 0.039 % aggregate over 512 Fibonacci geantino rays; sole 8192-ray divergence is `BucketLink2`, the not-closed mesh |
| costs, exact vs tessellated | 1.94× transport, +42 s one-off build (IRIS CloseShape), +174 MB; exact is *smaller on disk* |
| website (`website/`) | mesh viewer, charts, JS exact raytracer bit-true vs the kernel, event-display player, self-check 30/30; `website_data/decline_reasons.json` regenerated 2026-08-25 post-R5 |
| ray bridge (`tgeoRayService.py`) | the real kernel serving pixels: 10k rays in 4 ms |
| oTOF converts | 20 prototypes / 62 628 placements, 20/20 exact surfaces, 19/20 CSG `TGeoBBox` at dV_sym = 0 (`Stream_AC_OTOFTraversal.md`) |
| TGeo → STEP round trip | six-module study: ~7.4 M Contains samples, ONE disagreement (`Stream_AD/AF/AG/AH/AI`) |
| benchmark JSON (`website_data/`) | five hero parts complete; timing under load flagged `timingPreliminary` |

**Quote it correctly (unchanged):** `cyl_inter_cyl`'s mesh is closed as a triangle set; the
X-ray's lost rays are a **navigation** loss in `O2Tessellated`'s stepping, not holes.

**Quote the flat solid correctly:** its `Capacity` is OCCT's `GProp` on the CAD pieces, not an
independent measurement of the shipped solid; its `Safety` is a structural lower bound, not a
distance, and returns `0.` for 78–100 % of interior points.

## Open, in the order I would take them

1. **R6, the breadth census** (`Handoff_FlatCSG.md` §2) — cheap, one command per module, and it is
   what makes every claim above a claim about the geometry rather than about five detectors.
2. **The closure test** (`Handoff_ClosureTest.md`) — fresh session.
3. **The flat solid's three named follow-ups** (`Stream_AK_FlatCSG.md` §12), in order: the
   1-Lipschitz per-halfspace `Safety`; front-to-back ordering and early exit in
   `DistFromOutside` (the one kernel the flat solid can lose); and split at every carrier crossing
   **together with** raising `decompose.PART_MAX_CELLS` — priced, and neither half is worth doing
   alone: `IBCYSSFlangeA` with the budget raised decomposes into 157 cells in 27 s and then
   declines on the boundary gap instead. The streamer question is **closed, not open**: the
   `#pragma read` rule in `DetectorsBaseLinkDef.h` fires on every read path (measured,
   `Stream_AK_FlatCSG.md` §9.1).
4. **Converter: the JIT namespace bug.** `geom.C` emits `LoadSurfaceSolid` as a namespace-scoped
   forward declaration; the JIT wrapper makes it resolve wrongly and **the simulation continues
   silently without the module**. Workaround committed (`integration_demo/patch_exact_macro.py`);
   fix the converter's macro emission and make a failed module load loud.
5. **Converter: parallel `--csg auto` runs race** (two parallel conversions lost shapes):
   serialize/lock the deferred emit. And score **oTOF through the oracle gate** (converts, never
   scored).
6. **Kernel-confirmed: `ST1829909_01` leaks parity through real inter-face slits** (8/15 299
   rays; `maxSharedEdgeDeviation` ≈ the declared model tolerance). Needs the per-face localiser,
   then widened crossing acceptance near identified shared edges or trim-snapping.
7. **Kernel: loose bounding boxes.** `Stick` claims 40.1 cm in z around an 11.0 cm solid (3.6×).
   Use the cover-box union in `ComputeBBox`.
8. **Writer:** degenerate prism sections block the last TRD demand case (`B045cut`, a
   `TGeoTrd1` with `dx1 = 0`) and ITS's `IBGammaConvWireOuterSupport` — declined as "[2, 4]
   distinct vertices". Still open: STEP-writer segfault bracket (38 676 fine / 74 601 crash),
   `TGeoPara`, `TGeoHalfSpace` (732 TPC placements), instruments B/C + ITS/MAG re-runs, and a
   loud-refusal path for genuine scales.
9. **Small recognition items** from `Stream_AJ_Recognition.md` §7: run + score the **MAG corpus**
   (demand 14 prisms, never converted — folded into R6); the cell leaf budget (8) vs ITS's two
   10-halfspace connector blocks; Xtru parameter comparison in `checkKnownSource.py`; a face-count
   cap before the cell matcher's boolean build if a corpus ever makes its ~12 min/1170 parts
   painful.
10. **Hand-written-geometry findings on the record** (unchanged, Sandro's call on upstreaming):
    PIPE bellows duplicates FIXED on fork branch (`4fe6b285e0` here); **`RB26s3Bellow` has ZERO
    daughters — missing steel**, restoring it changes the material budget; TPC prepreg z-typo
    FIXED (`3755c83277` here); ROOT `TGeoShapeAssembly` defects (see item 13).
11. **Track 3b, the real data:** MCStepLogger→`events.json` exporter, then the website's event
    tab replays the real IRIS/Bagger transports.
12. **Materials matching:** 26/55 IRIS volumes matched, ten by accidental prefix. Anchored
    part-number matching.
13. **`O2BVHAssembly` landed** (`Stream_AE_BVHAssembly.md`): flat oTOF Contains 4.3×,
    Safety(out) 28.7×, transport 6.5×; honest limit at ≤68 daughters; `FindNode` untouched
    (§8's three hook options — Sandro to pick). Two measured ROOT defects
    (`TGeoShapeAssembly::DistFromOutside` Big() outside the bbox of a voxelized assembly;
    `Safety` exceeding the true minimum) — upstream-report decision is Sandro's.
14. **The talk (Track 4):** assemble from Review + Stream_Z + Stream_AJ + Stream_AK + the website;
    re-run `timingPreliminary` numbers on a quiet box; single-file website bundle.
15. **The face-normal gate column** and **the `auto`-mode unreliable-shipping policy** — both
    pre-corpus items, now six hand-overs old.
16. **Standing, unchanged:** free-form surfaces; models-are-not-legal overlaps and the broken
    `CheckOverlaps` on our shape; `Curve2D::closestPoint` as the kernel hot spot; mesh healing;
    composite trees (Tier 3 proper) remain deliberately not built.

## Traps in the environment (all still current)

- Env `O2/latest-swenzel-bvhsurfacesolid-o2`, build dir
  `B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2` (stack of 2026-08-22).
- **pythonOCC needs `--no-system SWIG`** when rebuilding (system 4.2.0 fails the required range).
- **Keep the converter env and the sim env separate.** The pythonOCC `PYTHONPATH` prepends make
  `o2-sim` segfault at startup; the converter needs O2/ROOT *in addition* to OCC — a bare OCC
  shell silently loses CSG (the deferred-emit WARN now says so per part). `runOracleGate.py`
  handles this itself: it *prepends* OCC to the inherited environment rather than replacing it.
- **An `O2FlatCSG` read off a ROOT file IS closed** — its boxes and BVH are transient by design,
  but the `#pragma read` rule in `DetectorsBaseLinkDef.h` calls `CloseShape()` on every object
  ROOT reads back, and it was measured firing on `TFile::Get` and on `TGeoManager::Import`, from
  C++ and from PyROOT (`Stream_AK_FlatCSG.md` §9.1). Three readers re-close explicitly anyway
  (`checkKnownSource.py`, `harness::loadShapeFromRootFile`, `geom.C`); that is belt-and-braces and
  is idempotent. **A shape you build by hand still owes itself a `CloseShape()`.**
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
ctest -R 'FlatCSG|BVHSurfaceSolid|BVHAssembly'

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

# the flat halfspace solid, scored against the SAME part as a plain composite (O2 env):
#   two subjects, one raster, one sample set; ratios > 1 mean the flat solid is faster
$B/stage/bin/o2-bench-detectorsbase-xray --perf --raster 32 \
    --shape  <dir>/shape_<part>.root \
    --flatcsg <dir>/flatcsg_<part>.bin
#   and the two split knobs, for sweeping them from outside the class:
#     --flat-split-depth N   --flat-min-box-fraction X

# the website (serve locally, then open the printed URL)
cd scripts/geometry/website && ./fetch_testdata.sh <gate-workdir> && python3 -m http.server 8231
# the ray bridge (inside the O2 env), for the website's RemoteEngine / engine-diff view
python3 scripts/geometry/tgeoRayService.py --port 8077
# the integration demo, end to end
scripts/geometry/integration_demo/   # see its README / convert_all.sh, run scripts
```
