# Website benchmark data — schema and provenance

This directory is consumed directly by the presentation website
(`scripts/geometry/website/`, Track 3, committed separately) via `index.json` — the file list its
`js/data.js` reads, since the static-hosting/Artifact model that site is built for has no directory
listing. It holds one JSON file per hero part, a cross-part `summary.json`, and that `index.json`.
Every number in every per-part file was produced by one of three existing instruments —
`runOracleGate.py`, `runXRayBench.py`, and `o2-bench-detectorsbase-xray --perf` — never invented,
never hand-edited. `makeWebsiteBench.py` (in the parent directory) reads those tools' own JSON
output and reshapes it into the schema below; it runs no conversion, no gate, no benchmark itself.

Read `SolidNavigationHarness.md` and `Stream_J_XRay.md`/`Stream_P_RepresentationBench.md` before
quoting any of these numbers in a talk: the ground rules there (mesh is a reference, not the
truth; an aggregate says *that*, never *where*; quote the sampling with every result) apply here
unchanged.

## The five hero parts

| file | part | corpus | why it is here |
| --- | --- | --- | --- |
| `BoomCylinderOuter.json` | `Bagger/BoomCylinderOuter_0_1_1_9` | Bagger | the CSG tube-union story: ships as `TGeoCompositeShape` |
| `Bucket.json` | `Bagger/Bucket_0_1_1_6` | Bagger | 97 faces, spheres + tori; the last part on this corpus to ship tessellated, made exact by the ellipse trim (`Stream_Q_EllipseTrim.md`) |
| `cyl_inter_cyl.json` | `cyl_inter_cyl/cyl_inter_cyl_0_1_1_1` | synthetic fixtures | the Steinmetz solid; the mesh-watertightness story |
| `torus_union_cyl.json` | `torus_union_cyl/torus_union_cyl_0_1_1_1` | synthetic fixtures | the torus/quartic story |
| `ST1829909_002.json` | ALICE3 leaf solid `ST1829909_002` (`CAD_noETA.stp`) | ALICE3 | scale: 965 patches, the sub-patch BVH part (`Stream_X_SubPatchBVH.md`) |

Plus `summary.json` (the cross-part headline table, see below) and `index.json` (`{"files": [...]}`, the list `js/data.js` reads to find the five files above -- regenerate it by re-running `makeWebsiteBench.py` rather than hand-editing it if a file is added or renamed).

## Schema

```
{
 "meta": {
   "part": "...",                    // the file's own hero-part name
   "partId": "...",                  // the id the underlying tools use (manifest/gate/xray key)
   "model": "...",                   // source CAD file (or "synthetic fixture" for the ladder)
   "story": "...",                   // one sentence: why this part is on the website
   "generated": "2026-08-22",
   "machine": "aarch64 10-core shared (other tenants running o2-sim concurrently)",
   "loadAverageDuringTiming": <float or null>,   // 1-min `uptime` load average taken immediately
                                                  // before the --perf run that produced this
                                                  // file's "functions"/"transport" numbers
   "timingPreliminary": <bool>,       // true whenever loadAverageDuringTiming > ~2 (see below)
   "samples": { ... }                 // exact points/rays/seed/raster/beams for every instrument,
                                       // see "Samples" below
 },
 "representations": [
   {
     "name": "surface" | "mesh" | "shape",
     "class": "o2::base::O2BVHSurfaceSolid" | "o2::base::O2Tessellated" | "TGeoCompositeShape" | ...,
     "primitives": N,                 // patches / triangles / composite leaves
     "primitiveKind": "patches" | "triangles" | "leaves",
     "memoryBytes": N or null,        // --perf's structural memory model (see below), null if
                                       // --perf did not report this representation
     "memoryFormula": "...",          // the formula --perf used to arrive at memoryBytes
     "sidecarBytes": N or null,       // size on disk of the sidecar this representation loads from
     "closeSeconds": N or null,       // CloseShape()/BVH-build wall time, from --perf
     "meshClosedBody": true|false|null,  // O2Tessellated's own watertightness statement (mesh only)
     "reliability": "reliable"|... or null,     // O2BVHSurfaceSolid navigability (surface only)
     "navigable": true|false or null,
     "functions": {
       "contains":         {"nsPerCall": N, "nsPerCallMin": N, "nsPerCallMax": N},
       "distFromOutside":  {...},
       "distFromInside":   {...},
       "safety":           {...}
     } or null,
     "transport": {"nsPerRay": N, "nsPerCrossing": N} or null,
     "accuracy": {
       "capacityRelativeDeviation": N or null,   // null when Capacity() is not comparable (below)
       "capacityComparableNote": "..." or null,  // why it is null, when it is
       "disagreements": {"contains": n, "distout": n, "distin": n, "safety": n},
       "unexplained":   {"contains": n, "distout": n, "distin": n, "safety": n},
       "oracleTolerance": N,          // the model's own declared BRep tolerance, cm
       "oracleValid": true|false      // BRepCheck_Analyzer on the reference .brep
     } or null,
     "xray": {
       "lost": n, "extra": n, "displaced": n,
       "unterminated": n, "parity": n, "parityNearBoundary": n, "oddCrossingLists": n,
       "raysTotal": n, "raysIdentical": n, "worstDeltaTCm": N
     } or null
   }, ...
 ]
}
```

`summary.json` carries, per part: `ships` (the representation the converter's cascade actually
emits), `speedRatioSurfaceVsFastest` (per function, the surface solid's ns/call divided by the
fastest of `shape`/`mesh` — this is the number Track 2's plan calls "the CSG story, 40-145x"), and
`meshLeak`/`surfaceLeak` (the mesh's and the surface solid's own X-ray robustness counters, so the
watertightness claim is visible without opening the per-part file).

### Two fields that are deliberately `null`, not a made-up number

- **`accuracy.capacityRelativeDeviation` is `null` for every `shape` representation built from a
  `TGeoCompositeShape`.** `TGeoShape::Capacity()` is Monte-Carlo sampled in ROOT for a composite
  (~1e-2 relative noise against the 1e-6 band the gate uses everywhere else), so a capacity
  deviation number there would fail a correct shape by two orders of magnitude for no reason. The
  gate's own volume criterion for a CSG part is `dV_sym`, the OCCT symmetric-difference volume
  computed at conversion time (reported in the converter's `csg_report.json`, not repeated here
  because it is a build-time, not a query-time, number).
- **`ST1829909_002.json` has no `shape` representation at all** — the part is declined for CSG
  recognition (a toroidal face is out of the recogniser's scope, `csg_report.json`) and it ships
  `surface`. Its `representations` array therefore has exactly two entries, `surface` and `mesh`.

## Samples — quoted with every result, never assumed

Per the project rule ("an aggregate says *that*, never *where*; a counter is a statement about a
ray lattice, not a solid"), every number here came from a specific, seeded sample or ray set. The
three instruments use three independent sample sets, on purpose (a single-shot accuracy query and
a transport step are different failure modes) — this is why `meta.samples` has three sub-blocks:

| sub-block | instrument | fields |
| --- | --- | --- |
| `gate*` | `runOracleGate.py` | `gatePoints`, `gateRays`, `gateSeed` — points/rays per category fed to the harness, then answered by `occtOracle.py` |
| `xray*` | `runXRayBench.py` | `xrayRaster` (N×N rays per beam), `xrayBeams` (Fibonacci-spiral beam directions; **never the axis raster alone** — a parallel beam is direction-poor, see `Stream_J_XRay.md`) |
| `perf*` | `o2-bench-detectorsbase-xray --perf` | `perfPoints`, `perfRays`, `perfPasses` (timed passes, 2 untimed warmup passes before them) |

Fixtures and Bagger both used `--beams 96` at the default raster `48` (48×48×96 = 221 184 rays per
part). **ST1829909_002 used `--raster 16` instead of the default 48** (16×16×96 = 24 576 rays) —
a deliberate cost decision, not a default: at raster 48 the OCCT crossing-list oracle
(`xrayOracle.py`) costs on the order of 15-35 ms per ray against a 965-face BRep (it classifies
every candidate interval's midpoint with `BRepClass3d_SolidClassifier`), which put the full
221 184-ray run at multiple hours on a shared box. Raster 16 is the same value
`Stream_L_ALICE3Defect.md` §2.5 already used for this model for the same reason, so this is a
reproduction of an established choice, not a new one. It is recorded per-part in `meta.samples` so
no reader compares `ST1829909_002`'s ray counts against the other four parts' without noticing the
raster differs.

## Timing: single-threaded, and the load average is reported not hidden

Per the environment rules for this session (another agent was running `o2-sim` jobs on the same
10-core box throughout), all conversion, accuracy (`--ref-answers`) and X-ray transport
(`--ref-crossings`) work was done first, and every `--perf` (timing) run was taken last,
single-threaded, with the 1-minute `uptime` load average recorded immediately before the run and
written into `meta.loadAverageDuringTiming`. `meta.timingPreliminary` is `true` whenever that load
average exceeds ~2 (i.e. more than about one core of contention beyond this session's own timing
process) — read the `functions`/`transport` numbers in such a file as directionally correct, not as
a clean single-tenant measurement, and prefer the corresponding numbers already published in
`Stream_P_RepresentationBench.md`/`Stream_X_SubPatchBVH.md` (measured on a quieter machine) if the
two disagree by more than the stated ratio uncertainty there.

## Exact commands (reproduce from scratch)

All of this needs the environment blocks from the top of this Track's task (O2 build env for the
gate/harness/xray binaries, **plus** the pythonOCC env, prepended never replacing, for every
Python conversion/oracle step). `<W>` below is any scratch workdir outside
`scripts/geometry/{STEP_examples,ALICE_3_example,IRIS}`.

**Fixtures corpus** (`cyl_inter_cyl`, `torus_union_cyl`, and the rest of the synthetic ladder):

```bash
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir <W>/gate_fixtures --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py  --workdir <W>/xray_fixtures \
    --reuse-db <W>/gate_fixtures/db --beams 96
$B/stage/bin/o2-bench-detectorsbase-xray --db <W>/gate_fixtures/db --perf \
    --json <W>/perf_fixtures/perf.json
```

**Bagger corpus** (`BoomCylinderOuter`, `Bucket`):

```bash
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir <W>/gate_bagger \
    --model scripts/geometry/STEP_examples/Bagger.step
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py  --workdir <W>/xray_bagger \
    --reuse-db <W>/gate_bagger/db --beams 96
$B/stage/bin/o2-bench-detectorsbase-xray --db <W>/gate_bagger/db --perf \
    --json <W>/perf_bagger/perf.json
```

**ALICE3 `ST1829909_002`** — converted alone via `--include-name`, so the run costs 10-14 s and
under 1.1 GB resident rather than repeating the whole-model conversion `NEXT.md` warns is unsafe at
the converter's default mesh precision. This single part meshes safely at the converter's default
`--mesh-prec` (11.8 s, 1.1 GB peak RSS, measured with `/usr/bin/time -v`) — the documented 22.9 GB
failure is a whole-55-solid-model number (`Stream_P_RepresentationBench.md` §6.1: the whole model
is also safe at `--mesh-prec 0.5` or coarser; this run used one part at the default precision,
which is smaller still):

```bash
# converter env (prepended to the O2 env)
python3 scripts/geometry/O2_CADtoTGeo.py scripts/geometry/ALICE_3_example/CAD_noETA.stp \
    --include-name "ST1829909_002" --exact-surfaces auto --csg auto --mesh --dump-brep \
    --surface-report <W>/alice3_st002/db/surface_report.json \
    --output-folder <W>/alice3_st002/db -o <W>/alice3_st002/db/geom.C

# O2 env: samples + OCCT single-shot oracle + scoring (mirrors runOracleGate.py's own steps,
# run by hand because the gate's own driver always converts the WHOLE model)
$B/stage/bin/o2-bench-detectorsbase-solid-harness \
    --surfaces <W>/alice3_st002/db/surfaces_ST1829909_002_0_1_1_36.bin \
    --facets   <W>/alice3_st002/db/facets_ST1829909_002_0_1_1_36.bin \
    --points 3000 --rays 3000 --seed 1 --dump-samples <W>/alice3_st002/oracle
# converter env
python3 scripts/geometry/occtOracle.py --brep <W>/alice3_st002/db/brep_ST1829909_002_0_1_1_36.brep \
    --samples <W>/alice3_st002/oracle/samples_adhoc.json \
    --out <W>/alice3_st002/oracle/answers_adhoc.json
# O2 env
$B/stage/bin/o2-bench-detectorsbase-solid-harness \
    --surfaces <W>/alice3_st002/db/surfaces_ST1829909_002_0_1_1_36.bin \
    --facets   <W>/alice3_st002/db/facets_ST1829909_002_0_1_1_36.bin \
    --points 3000 --rays 3000 --seed 1 --loop-crosscheck \
    --ref-answers <W>/alice3_st002/oracle \
    --json <W>/alice3_st002/oracle/ref_answers_run.json

# X-ray transport, at raster 16 (see "Samples" above for why)
$B/stage/bin/o2-bench-detectorsbase-xray \
    --surfaces <W>/alice3_st002/db/surfaces_ST1829909_002_0_1_1_36.bin \
    --facets   <W>/alice3_st002/db/facets_ST1829909_002_0_1_1_36.bin \
    --raster 16 --beams 96 --dump-rays <W>/alice3_st002/xray
# converter env
python3 scripts/geometry/xrayOracle.py --brep <W>/alice3_st002/db/brep_ST1829909_002_0_1_1_36.brep \
    --rays <W>/alice3_st002/xray/xrays_adhoc.json \
    --out  <W>/alice3_st002/xray/crossings_adhoc.json
# O2 env
$B/stage/bin/o2-bench-detectorsbase-xray \
    --surfaces <W>/alice3_st002/db/surfaces_ST1829909_002_0_1_1_36.bin \
    --facets   <W>/alice3_st002/db/facets_ST1829909_002_0_1_1_36.bin \
    --raster 16 --beams 96 --ref-crossings <W>/alice3_st002/xray \
    --json <W>/alice3_st002/xray/xray.json

# timing (single-threaded, last, load average recorded from `uptime` immediately before)
$B/stage/bin/o2-bench-detectorsbase-xray \
    --surfaces <W>/alice3_st002/db/surfaces_ST1829909_002_0_1_1_36.bin \
    --facets   <W>/alice3_st002/db/facets_ST1829909_002_0_1_1_36.bin \
    --perf --json <W>/alice3_st002/perf/perf.json
```

**Generate the website JSON** from the nine files above:

```bash
python3 scripts/geometry/makeWebsiteBench.py \
    --gate-fixtures <W>/gate_fixtures/gate.json --xray-fixtures <W>/xray_fixtures/xray.json \
    --perf-fixtures <W>/perf_fixtures/perf.json --load-avg-fixtures <uptime 1-min avg> \
    --gate-bagger   <W>/gate_bagger/gate.json   --xray-bagger   <W>/xray_bagger/xray.json \
    --perf-bagger   <W>/perf_bagger/perf.json   --load-avg-bagger   <uptime 1-min avg> \
    --gate-alice3   <W>/alice3_st002/oracle/ref_answers_run.json \
    --xray-alice3   <W>/alice3_st002/xray/xray.json \
    --perf-alice3   <W>/alice3_st002/perf/perf.json --load-avg-alice3 <uptime 1-min avg> \
    --gate-points-alice3 3000 --gate-rays-alice3 3000 --xray-raster-alice3 16 \
    --out-dir scripts/geometry/website_data
```

This writes the five `<part>.json` files, `summary.json`, and `index.json` (the file list
`scripts/geometry/website/js/data.js` reads to discover what is in this directory -- there is no
directory listing in the static-hosting/Artifact model that site is built for).

## Known residuals, reported rather than hidden

- **`cyl_inter_cyl` fails its own single-shot oracle gate**: 1 of 2000 `distout` rays disagrees
  with OpenCascade outside the model's declared tolerance (worst deviation 5.65e-06 cm — a grazing
  ray at machine-precision scale). `oblique_cut_cyl` (not a hero part here, but in the same gate
  run) shows the identical pattern. Neither is papered over: it is the number the gate printed on
  this run, with this seed. See `NEXT.md` item 8 for the standing, not-yet-localised class of
  small single-ray residuals this belongs to.
- **`cyl_inter_cyl`'s X-ray transport (96 beams) shows `surface`: extra=2, unterminated=2,
  parity=2 of 221 184 rays** — this reproduces `NEXT.md`'s recorded number exactly. Its `mesh`
  representation shows a much larger `lost` count at this raster/beam density; the specific
  `unterminated` count on the mesh is **known to move run-to-run at the last-ULP level**
  (`Stream_J_XRay.md` §"A robustness property that fell out by accident": the same ray lattice
  computed by two algebraically identical but differently-rounded expressions moved the mesh's LOST
  count from 14 to 40 on this exact fixture corpus, while `surface`/`shape` stayed bit-identical).
  Whatever this run's JSON says is what was measured on this run, with this seed; do not expect it
  to reproduce a specific historical count for the mesh column bit-for-bit.
- **`ST1829909_002`'s `contains` oracle check has 2 unexplained mismatches** (of 7472 scored
  points) and `safety` has 65 (of 5500) at ~1e-5 to 1e-6 cm deviation — both are recorded in the
  file's `accuracy` block rather than rounded away; they sit right at the edge of the declared
  model tolerance (4.3e-4 cm) and have not been individually localised to a face.
- **`Bucket`'s X-ray transport (96 beams) shows `surface`: extra=1, unterminated=1, parity=1 of
  221 184 rays**, even though its single-shot oracle gate is 0/0/0/0 clean (`gate.json`,
  `runOracleGate.py`) -- the two instruments measure different failure modes (a sampled query vs.
  a stepped transport) and this is the first time a transport-only residual on this part has been
  recorded rather than assumed clean by extension from the gate. Not yet localised to a face.
- **`ST1829909_002`'s X-ray transport (96 beams, raster 16) shows `surface`: extra=6, parity=1 of
  24 576 rays, 0 lost.** The `mesh` representation on the same rays shows `lost=13, extra=17,
  displaced=6499` -- consistent with the corpus-wide finding in `Stream_P_RepresentationBench.md`
  §6 that ALICE3's coarse-precision mesh "loses walls" that the exact solid does not.
