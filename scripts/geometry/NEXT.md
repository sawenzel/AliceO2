# NEXT — session-start instruction for the CAD → TGeo work

This file is the current hand-over. Whoever finishes a session should **rewrite it**.
Last rewritten 2026-08-22, at the end of the return-to-project session (deep review,
stack re-verification, and the three parallel presentation tracks).

Branch `swenzel/bvhsurfacesolid`. Everything below is committed unless marked otherwise.

## Read this first

[`Tutorial.md`](Tutorial.md) is still the map; [`Review_2026-09.md`](Review_2026-09.md) is the
deep review of the mathematics and algorithms with the verification appendix;
[`INDEX.md`](INDEX.md) orders every document; [`Plan_Presentation.md`](Plan_Presentation.md) is
the plan of record for the WG talk (~1.5 weeks out from this rewrite).

## Where the branch stands (all re-verified 2026-08-22 on the rebuilt stack)

| | |
| --- | --- |
| `ctest -R BVHSurfaceSolid` | **113 cases**, green |
| `O2_CADtoTGeo.py --self-test` | **48 checks** = four suites 18+8+10+12 — never quote the last line alone |
| `runOracleGate.py --self-test` / `csg/emit.py --self-test` / xray `--self-test` | 17/17 · 33/33 · clean |
| Bagger gate | **13/13, CSG 7 / surfaces 6 / tessellated 0, exit 0**, all columns clean |
| fixtures gate | 10/10 scored; **2 `distout` sliver rays**, probed and explained (Review Appendix A) |
| **Geant integration demo** | **works**: IRIS + Bagger as sensitive external detectors, real hits, geantinos/e/π through both representations, zero stuck tracks / nav errors / aborts (`Stream_Z_IntegrationDemo.md`) |
| material budget, exact vs tessellated | **0.039 %** aggregate over 512 Fibonacci geantino rays, with both positive controls; the sole 8192-ray divergence is `BucketLink2`, the not-closed mesh |
| costs, exact vs tessellated | 1.94× transport, +42 s one-off build (IRIS CloseShape), +174 MB; exact is *smaller on disk* |
| the talk visual | exact IRIS hits at exactly r = 5.460 cm; tessellated smeared 5.427–5.460 (330 µm sagitta) |
| website (`website/`) | mesh viewer, charts, **JS exact raytracer validated bit-true against the kernel** (0 disagreements, max Δt = 0 in float32, 6 parts × 141k rays), event-display player, self-check 30/30 |
| ray bridge (`tgeoRayService.py`) | the real kernel serving pixels: 10k rays in 4 ms |
| benchmark JSON (`website_data/`) | five hero parts complete; ram shape 141× surface on `Contains`; ALICE3 part mesh 12–62× faster than surface; timing under load flagged `timingPreliminary` |

**Sharpened during Track 3, quote it correctly:** `cyl_inter_cyl`'s mesh is *closed as a
triangle set* (0 odd-parity rays over 60k). The X-ray's 6 lost rays are a **navigation** loss in
`O2Tessellated`'s stepping. "The mesh leaks" means navigation, not holes — which is the stronger
statement for the talk.

## Open, in the order I would take them

1. **Converter: the JIT namespace bug.** `geom.C` emits `LoadSurfaceSolid` as a namespace-scoped
   forward declaration; `CADGeometryUtils.cxx` wraps the macro in `o2_cadgeom_<TAG>_<N>`, the
   symbol resolves against the wrong namespace, the macro fails to compile, and **the simulation
   continues silently without the module**. This is why the exact path had never run in o2-sim.
   Workaround committed (`integration_demo/patch_exact_macro.py`); fix it in the converter's
   macro emission (emit an `#include`), and consider making a failed module load loud.
2. **Converter: `--csg auto` outside the gate env yields CSG 0, silently.** The deferred shape
   completion needs ROOT; without it every candidate declines, and the six Bagger cylinders
   decline with the empty reason `declined CSG: None`. Error loudly when ROOT is unreachable and
   fix the empty reason. (This produced the false "baseline moved" in `Stream_Z` §2, reconciled
   there.)
3. **The CSG path forward is now measured** (`Stream_AA_FlatCSG.md`, 2026-08-22): the cell-count
   table Tier 3 was gated on exists at last, and the feared blow-up did not happen — Bagger's
   decliners decompose to 7–28 cells and 15/16 tessellated-but-analytic ALICE3 solids to 6–14,
   by OCCT's own splitter at converter time, with every failure traced to a NURBS carrier
   (making Tier-0 canonicalisation the measured prerequisite). Order of work: Tier-0
   canonicalisation → the risk-free single-cell emitter (retires `cyl_inter_cyl` and
   `tube_window` exactly — and with them the two known sliver rays) → balanced
   `TGeoCompositeShape` trees for the 16 ALICE3 solids under the existing two-test acceptance.
   A flat `TGeoBVHCSG` is legitimate but deferred behind falsifiable criteria (Stream_AA §5):
   build it only if the trees measure too slow on the X-ray bench or at AOT-codegen time.
   Per-part decline reasons now ship in `csg_report.json` (`whyNotCSG`) and
   `surface_report.json` (`why_not_surface`), joined in `website_data/decline_reasons.json`
   (78 parts, 0 missing); the empty-reason and silent-no-ROOT defects of old item 2 are fixed.
4. **NEW, kernel-confirmed: `ST1829909_01` leaks parity through real inter-face slits.**
   The website's JS raytracer flagged odd-parity rays; a kernel probe confirms **8 of 15 299
   hitting rays** (BVH ≡ Loop, `onTrimBoundary` = 0 on all of them — not the sliver class), and
   the mechanism is measured: `maxSharedEdgeDeviation` = 4.70e-4 cm against the model's own
   declared tolerance 4.723e-4 cm, i.e. the faces genuinely sit up to ~5 µm apart along shared
   edges. Identity closure calls this closed *by design* (topology, not geometry — the documented
   caveat, now measured biting). The slit is wider than the trim boundary band, so no ambiguity
   flag fires and the vote never triggers. Fix directions to evaluate: widen crossing acceptance
   near identified shared edges by the sidecar's own measured deviation, or snap paired trims to
   one canonical curve per edge identity (the per-edge deviation data already exists). Needs the
   per-face localiser first. Probes: scratchpad `probe_parity.cxx` / `probe_dev.cxx`.
5. **Kernel: loose bounding boxes.** `Stick` claims 40.1 cm in z around an 11.0 cm solid (3.6×),
   verified against OCCT `Bnd_Box` and the mesh (<0.1 mm apart). Costs broad-phase efficiency
   everywhere and confuses CheckOverlaps. Likely the conservative per-family boxes (full rim
   circles / full torus) propagating into `ComputeBBox`; the sub-patch cover boxes are already
   tighter — use their union.
6. **Track 3b, the real data:** an exporter MCStepLogger trees → `events.json` (schema in
   `website/README.md`; tree layout in `integration_demo/data/README`), then the website's event
   tab shows the real IRIS/Bagger transports. Data is preserved (uncommitted) in
   `integration_demo/data/` and the scratch.
7. **Materials matching:** 26/55 IRIS volumes matched, 29 fell back to vacuum `Default`, and ten
   of the 26 match only by accidental string prefix (`ST0923290_01` ⊂ `ST0923290_010…_019`).
   Needs anchored/exact part-number matching.
8. **Gate: credit the sliver.** Teach the relabel class to explain a distout mismatch whose
   candidate crossing is `onTrimBoundary`-flagged within the trim band (Review Appendix A).
9. **The face-normal gate column** (fourth hand-over in a row) and **the `auto`-mode
   unreliable-shipping policy** (~20 lines) — both still open, both pre-corpus items.
10. **The talk (Track 4):** assemble from Review + Stream_Z + website; re-run the
   `timingPreliminary` numbers on a quiet box; build the single-file website bundle for
   publishing (all data inlined — the Artifact CSP allows nothing external, and the bridge is
   local-only).
11. **Standing, unchanged:** oTOF XCAF traversal; free-form surfaces; the models-are-not-legal
   overlap finding and the broken `CheckOverlaps` on our shape; `Curve2D::closestPoint` as the
   kernel hot spot; mesh healing.

## Traps in the environment (changed ones first)

- **The stack was rebuilt 2026-08-22 and the names changed**: env is
  `O2/latest-swenzel-bvhsurfacesolid-o2`, build dir
  `B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2`.
- **pythonOCC needs `--no-system SWIG`** when rebuilding: system SWIG 4.2.0 fails the required
  4.2.1…4.4.1 range and the recipe's system check does not catch it.
- **Keep the converter env and the sim env separate.** The pythonOCC `PYTHONPATH` prepends make
  `o2-sim` segfault at startup. And the converter needs O2/ROOT *in addition* for CSG emission
  (item 2 above): the gate composes both; a bare OCC shell silently loses CSG.
- **External-detector hits live in `o2sim.root`** (`IRISHit`, `BAGRHit` branches), not in
  per-detector `o2sim_Hits*.root` files, under `o2-sim-serial`.
- `rm` is blocked by a repo hook (`.claude/hooks/deny-deletions.py`); move files aside instead.
- Everything else in the 2026-08-09 list still applies: eval+command in one shell; stage lib
  first; prepend-never-replace; one ninja at a time; detach long runs to unique `--out` paths;
  never write into `STEP_examples/`; convert ALICE3 without `--mesh` (or `--include-name` one
  part); `manifest.json` stores absolute paths.

## Commands

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ctest -R BVHSurfaceSolid

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step

# the website (serve locally, then open the printed URL)
cd scripts/geometry/website && ./fetch_testdata.sh <gate-workdir> && python3 -m http.server 8231
# the ray bridge (inside the O2 env), for the website's RemoteEngine / engine-diff view
python3 scripts/geometry/tgeoRayService.py --port 8077
# the integration demo, end to end
scripts/geometry/integration_demo/   # see its README / convert_all.sh, run scripts
```
