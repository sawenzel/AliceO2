# Plan — the WG presentation, and the two weeks in front of it

Written 2026-09-20, on return to the branch after six weeks of silence. The presentation is a
~30 minute collaboration/WG talk in two weeks. This file is the plan of record for getting there;
it was agreed in session and the decisions below are Sandro's, not inferred.

## Decisions taken

- **Venue**: WG/collaboration meeting, ~30 min. Technical audience; live demo plausible.
- **Demo models**: **IRIS + Bagger**. oTOF stays blocked (XCAF traversal, `Handoff_IntegrationTest.md`
  §5) and appears in the talk as a named, measured blocker — not silently absent, and not fixed
  under time pressure.
- **Visualization**: a website with an in-browser raytracer of the exact solids is the ambition,
  with offline renders as the guaranteed fallback. It is **one track, not the center of the plan**:
  the review and the working simulation carry the talk; the website makes them visible.
- **Review first**: a deep review of the codebase (mathematical, algorithmic) precedes everything.
  It has two halves with different prerequisites: the **reading half** needs no build and starts
  now; the **verification half** (gates, ctest, self-tests re-run) waits for the alidist build in
  progress and re-confirms every number the review or the talk quotes.
- **Overlaps**: the "models are not legal" finding (`Stream_T_AssemblyOracle.md`) is presented
  honestly as one slide. The source models are not repaired for the demo.

## Tracks

### Track 0 — the deep review (now → build-done + 1 day)

The centerpiece of week 1, and what "review" means here: read the code, not the documents about
the code.

- **Reading half, no build needed**: `BoundedSurface.h` (6017 lines — the curves, trim
  containment, the six surface kernels, the quartic, capacity by divergence theorem, closure),
  `O2BVHSurfaceSolid.{h,cxx}` (navigation, BVH, sub-patch covers, safety bounds),
  `O2SurfaceSolidIO.cxx` (sidecar v3), `O2Tessellated.cxx`, the converter's recognition and
  extraction math in `O2_CADtoTGeo.py`, `csg/` and the oracle scripts.
- **Deliverable**: `Review_2026-09.md` — what was achieved mathematically and algorithmically,
  assessed from the source; new findings if the reading surfaces any; the honest gaps
  (free-form surfaces, arrangement-cell trims, mesh validity); and a verified-state appendix
  filled in once the build lands.
- **Verification half, build needed**: `ctest -R BVHSurfaceSolid` (113 cases), fixtures gate
  10/10, Bagger 13/13 exit 0, the four self-tests. Numbers that move get investigated before
  anything else proceeds.

### Track 0b — document order (1 day, alongside Track 0)

- `INDEX.md`: reading order (`Tutorial.md` → `Review_2026-09.md` → per-topic records), each
  stream document marked **landed / superseded / open**.
- Execute the Stream-C directory cleanup at last: `git rm` the stale `O2_TGeoToCAD*` variants,
  `O2_CADtoTGeo_with_*`, `*.py~`, `#…#`, `._*`, `a.out`, `*.o`, `__pycache__`; extend
  `.gitignore`. One mechanical commit.

### Track 1 — the demo simulation, IRIS + Bagger (build-done → ~day 8)

Execute `Handoff_IntegrationTest.md` with Bagger substituted for oTOF.

- IRIS as an **`externalDetector`** — sensitive volumes producing real hits
  (`o2sim_Hits*.root`), materials from `IRIS/IRIS_MATERIALS.csv` + the NIST JSON.
- Bagger as the second module, its `Bucket` sensitive: a 30-second memorable demo, and Bagger is
  the one model that is exact end to end.
- Each converted **twice** — exact cascade (`--csg auto --exact-surfaces auto`) and pure
  tessellated — and compared on steps (MCStepLogger, per volume), secondaries, CPU, memory, and
  a geantino **material-budget fan** (Fibonacci directions, never axis rasters), which is the
  strongest equivalence test available.
- Anchor/placement per the handoff §4 (`barrel`, y = +30); geantinos first, then real particles;
  deterministic seeds; everything reproducible from a script in the repo.
- Deliverable: `Stream_Z_IntegrationDemo.md` with the numbers, failures stated as plainly as
  successes.

### Track 2 — single-part benchmark data (parallel with Track 1)

Re-run the representation benchmark and X-ray transport on the hero parts, emitting JSON the
website consumes directly:

- a ram (`TGeoTube ∪ TGeoTube`; the CSG story, 40–145×),
- `Bucket` (97 faces, spheres + tori; the last mesh part, made exact by one ellipse trim),
- `ST1829909_002` (ALICE3 scale; the sub-patch BVH numbers),
- `cyl_inter_cyl` (the part whose **mesh leaks** — 6 rays enter the tessellation and never leave),
- one torus-bearing part for the quartic story.

Per part: per-function ns (`Contains`, `DistFromOutside`, `DistFromInside`, `Safety`), memory,
oracle deviation, X-ray counters, for all representations present.

### Track 3 — website + JS raytracer (days ~5–12, timeboxed)

An artifact page in three layers, in order of certainty:

1. **Guaranteed**: three.js viewer of `facets_*.bin` meshes; charts from Track 2's JSON;
   side-by-side coarse/fine tessellation against the exact silhouette.
2. **The ambition**: a JS raytracer loading `surfaces_*.bin` — analytic ray–patch intersection
   (quadratics for plane/cylinder/sphere/cone, the quartic for torus), trim containment by 2D
   parity, progressive rendering in a worker. Same-camera exact-vs-mesh renders with a
   difference view where faceting and leaks light up.
3. **Known-hard 10%, decided in advance**: exact parity against *rational B-spline trim curves*
   in JS. Fallback is densely sampled trim wires, labeled as such — the surfaces stay exact,
   which is what the eye sees.

If layer 2 slips, layer 1 plus offline raytraces from the C++ kernel is the talk's material.
The website never becomes the critical path: the timebox is hard.

### Track 4 — the talk (days 12–14)

~30 min: the problem (CAD → simulation, what tessellation costs); the three-representation
cascade; the oracle methodology as the credibility spine; the review's headline numbers; live
website demo; the o2-sim demo with hits; the one honest overlaps slide; gaps and roadmap
(free-form, oTOF, arrangement cells, Embree).

## Schedule and gates

| days | work | gate to pass before moving on |
| --- | --- | --- |
| 1–2 | Track 0 reading half, Track 0b | `Review_2026-09.md` drafted; INDEX committed |
| build-done | Track 0 verification half | recorded numbers re-confirmed or divergences named |
| 3–8 | Track 1 (owns the build), Track 2 (Python-only, parallel) | o2-sim transports both modules without stuck tracks; benchmark JSON exists |
| 5–12 | Track 3, timeboxed | layer 1 done by day 9 or layer 2 is cut |
| 12–14 | Track 4 | dry run against the clock |

## Risks, named

- **The Geant integration has never been executed.** Biggest unknown; scheduled earliest of the
  build-dependent work for exactly that reason.
- **The build in progress** may land late or broken; the reading review and Track 0b are
  deliberately build-independent so the first days cannot stall.
- **`--mesh-prec` on IRIS**: never default 0.1 (a 2 m sphere reached 22.9 GB); pick, justify,
  record.
- **The JS raytracer** is genuinely new code; its failure mode is bounded by the layer structure
  above.
- Shared machine, one `ninja` at a time, long runs detached to `--out` — all standing rules from
  `NEXT.md` apply.
