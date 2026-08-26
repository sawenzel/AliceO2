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
part emitted as a plain composite. **R6, the breadth census, is now done**: all seventeen Run 3
modules round-trip at **22 736 native CSG of 22 901 leaf solids (99.3 %)**, every one agreeing with
its source. The one thing still queued is the **closure test** (real physics through the round
trip).

## Two findings from 2026-08-25 that change how numbers here should be read

**1. The flat solid's speed-up is against the round-tripped composite, not against what O2 ships.**
`Stream_AK_FlatCSG.md`'s 3.2–34x on `Contains` is measured against the same cells emitted as a
plain `TGeoCompositeShape`. Measured instead against the **source** `TGeoShape` each part was made
from, on the same 4000 points in the same kernel, the picture depends entirely on how deep that
source tree is: TRD `BREF1`'s original is 3 leaves at boolean depth **2** and beats the flat solid
on all four kernels (0.9x / 0.4x / 0.2x / 0.24x); ITS `IBCYSSFlangeC`'s is 36 leaves at depth **35**
and loses 6.9x on `Contains` and 5.7x on `Safety`; TRD `VolTOFrail` (22 leaves, depth 8) is in
between. This does not undercut the representation — its purpose is CAD-native geometry where no
original exists — but any quoted ratio must now say *against what*. The baseline is on the website
per part, and `exportSourceShapes.py` is how it gets there.

**2. `tgeoRayService.py` was corrupting every composite it traced, and is fixed.**
`TGeoBoolNode`'s per-thread scratch state is indexed by `TGeoManager::ThreadId()`, which returns 0
for every thread unless the manager is multi-threaded — and that needs `SetMaxThreads()`, which
needs a closed geometry, which this service does not have. On `BREF1` with 8 threads the composite
missed 13 033 of 35 420 hits, invented phantom hits and returned normals up to 90° wrong on ~19 %
of pixels; at one thread every one of those is **0**. Boolean shapes are now traced single-threaded
(`raysvc::isThreadSafe`) and `/load` reports it. **No `--perf`, gate, X-ray or `/bench` number is
affected** — all of those are single-threaded by construction — but any *picture* of a composite
from this bridge taken before 2026-08-25 was wrong. It is a harness bug, not a ROOT defect: ROOT
documents the `SetMaxThreads` requirement and this service never met it.

**3. An exact tessellation is not a fallback, and 82 % of ITS+TPC+TRD has one.**
A triangulation of a planar polygon IS that polygon, so a part whose every face is a plane with
straight trim edges has a mesh that is the *same solid*, not an approximation. Measured: 30 random
all-planar parts, 600 000 points, **zero** `Contains` disagreements against the exact surface solid.
The census is 143/264 ITS, 72/172 TPC, 1107/1170 TRD — **1322 of 1606** — at a median of 12 facets,
because 92 % of them are boxes. `csg/planar.py` is the predicate (the same `PlanarPolygon` vs
`CurvedPlanar` rule `LoadSurfaceSolid` applies) and `csg_report.json` now carries
`tessellationExact` per part plus a census line in the tier table. **It is an annotation and
routes nothing**, deliberately: exactness says the mesh costs no accuracy for such a part, not that
it is the right thing to emit — a box is a `TGeoBBox` and should stay one. Two measurements sit
beside it for whenever emission is revisited: preferring the mesh over the exact `surface` solid
would move **zero** parts (across all seven corpora every one of the 25 `surface` parts has curved
faces), while against the **flat solid** the exact mesh wins every kernel on all nine all-planar
flat parts (median 2.0x `Contains`, 7.7x `DistFromOutside`, 9.8x `DistFromInside`, 3.5x `Safety`) at
72–324 facets. Sandro's call, 2026-08-25: no decision now; what a part emits is ultimately a small
per-part benchmark's job, and `/bench` already times every subject of one part on one shared sample
set. `MeshHealing.md`'s caveat — a mesh can be *invalid*, not merely inaccurate — is the argument on
the other side that none of these numbers address.

## R6, the breadth census: DONE (2026-08-25)

All seventeen Run 3 modules round-tripped and scored — **22 901 leaf solids**, ABSO CPV EMC FT0
FV0 HMP ITS MAG MCH MFT MID PHS PIPE TOF TPC TRD ZDC. The document is
`roundTripReport.py --corpus <root>`; the numbers:

| | |
| --- | --- |
| native CSG tree | **22 736** (99.28 %) |
| `O2FlatCSG` | 20 |
| exact surfaces | 130 |
| tessellated | 15 |
| agrees with its source `TGeoShape` | **22 747 / 22 756** |
| tessellation is exact (all-planar) | 19 622 (86 %) |
| declined by the *writer*, never reaching the converter | 31 volumes |

**Every native primitive class round-trips at 100 %** — `TGeoBBox` 18 069, `TGeoTube` 754,
`TGeoTrd1` 133, `TGeoTrap` 126, `TGeoPcon` 123, `TGeoTubeSeg` 100, `TGeoXtru` 45, `TGeoPgon` 33,
`TGeoEltu` 26, `TGeoArb8` 9, `TGeoCtub` 8 — except `TGeoCone` (41/42) and `TGeoTorus` (302/308).
**Everything that declines is a `TGeoCompositeShape`**: 3125 of them, 2969 tree / 18 flat / 123
surface / 15 mesh. Coverage is a question about booleans and about nothing else.

The writer's 31 declines are the other half, and they are invisible if only converted parts are
counted: **`TGeoPara` is not mapped at all** (13 volumes), six composites are *unbounded* (a
half-space has no B-rep body), six `TGeoTorus` are hollow, and the rest are prism section-count
mismatches in `TGeoTrd1`/`TGeoArb8`/`TGeoPgon` plus two `TGeoCompositeShape(union)`.

**The nine "disagreements" are not defects.** All nine are MFT bodies of a multi-body label where
the body is ~6 ppm of the label's volume (2.0e-5 cm^3 inside 3.49 cm^3), and
`checkKnownSource.contains_crosscheck` samples the *source's* bounding box, so 20 000 uniform
points never land in the body and the instrument correctly refuses to pass an empty comparison.
Re-scored in the emitted body's own box, 200 000 points each: 3342–41 868 points inside the body,
**zero** outside the source. The fix is to sample the emitted body's box when `one_way` is set;
it is not made, because it changes an acceptance instrument.

**MFT is an outlier and it is not our doing.** 16 965 of its 19 185 solid-carrying volumes are
geometric duplicates: `capacitor`, `welding0` and `welding1` are **5144 separate `TGeoVolume`
objects with 5144 separate `TGeoBBox` objects each**, all one box, one medium. 25 853 volumes
where ~2000 would do. It costs the manager, the STEP file, ~90 min of conversion (each duplicate
gets a full 4000-point acceptance test) and presumably `o2-sim` load time and navigation. Worth a
JIRA against the MFT geometry; the report has a section for it per module.

## The queued next step, for a FRESH session
1. **The closure test** — real `o2-sim -m PIPE ITS TPC MAG` physics through the round-tripped
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
| website (`website/`) | mesh viewer, charts, JS exact raytracer bit-true vs the kernel, event-display player, self-check **99/99** over 29 parts; `website_data/decline_reasons.json` regenerated 2026-08-25 post-R5. **Six subjects** are now first class (2026-08-25): the ten R5 flat-CSG parts are in `testdata/`, `tgeoRayService.py` loads every artefact by magic, and the live `/bench` measures each on the *same* sample points. Two of the six are not conversion products — the **Original TGeo** baseline (`exportSourceShapes.py`) and the cells-tree comparison — and a `role` field keeps both out of every "ships" statement |
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

1. **The closure test** (`Handoff_ClosureTest.md`) — fresh session. The only queued item left.
2. **The `one_way` sampling box in `checkKnownSource.py`** — a one-way containment comparison
   samples the *source's* bounding box, so a body that is parts-per-million of a multi-body label
   is never hit and the part reports as a failure it cannot be scored for. Nine MFT parts, proven
   clean on 200 000 points each in the body's own box. Sample the emitted body's box when
   `one_way` is set.
3. **MFT's duplicate logical volumes** — 16 965 of 19 185. Written up as
   [`MFT_deduplication_task.md`](MFT_deduplication_task.md) for its own session; the complication
   is that `MFTSensor` is one of the duplicated families and hit scoring resolves the sensor from
   four levels of copy numbers above it.
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

# the TGeo -> STEP -> TGeo report over a corpus (one directory per module, each with
# o2sim_geometry.root, <MOD>_writer_report.json and a conv/ holding csg_report.json)
python3 scripts/geometry/roundTripReport.py --corpus <root> --out report.md
python3 scripts/geometry/roundTripReport.py --corpus <root> --out report.html --html  # print to PDF
python3 scripts/geometry/roundTripReport.py --corpus <root> --part BREF1              # one part
# and the round trip that fills such a corpus, per module:
#   o2-sim-serial -n 0 -g boxgen -m $MOD
#   O2_TGeoToCAD.py o2sim_geometry.root $MOD.step --report ${MOD}_writer_report.json
#   O2_CADtoTGeo.py $MOD.step -o geom.C --output-folder conv --exact-surfaces auto --csg auto --mesh
#   checkKnownSource.py --original o2sim_geometry.root --writer-report ${MOD}_writer_report.json \
#       --converted conv --json conv/knownsource.json

# the website (serve locally, then open the printed URL)
cd scripts/geometry/website && ./fetch_testdata.sh <gate-workdir> && python3 -m http.server 8231
# the ray bridge (inside the O2 env), for the website's RemoteEngine / engine-diff view.
# It now loads all four representations and dispatches on the file's own magic, so a bridge left
# running from before 2026-08-25 will refuse flatcsg_*.bin: restart it after pulling.
python3 scripts/geometry/tgeoRayService.py --port 8077
# the integration demo, end to end
scripts/geometry/integration_demo/   # see its README / convert_all.sh, run scripts
```
