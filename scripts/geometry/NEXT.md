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
| **oTOF converts** (fixed 2026-08-22, `Stream_AC_OTOFTraversal.md`) | 20 prototypes / **62 628 placements**, 20/20 exact surfaces, **19/20 CSG `TGeoBBox` at dV_sym = 0**, 44 s / 753 MB, 3 964 `AddNode`; converter self-test now 48+6 |
| **TGeo → STEP round trip closes** (`Stream_AD_TGeoToStep.md`, 2026-08-23) | `O2_TGeoToCAD.py`: 16/17 Run 3 shape classes, self-test 71/71; PIPE: 1.69M Contains samples, **168/169 parts zero disagreements** (1 point open), CSG-vs-source 620k/0, capacity median 1.5e-13; full geometry 21 150 solids (99.7 %) in 3 min — the STEP *writer* segfaults on the full model |
| benchmark JSON (`website_data/`) | five hero parts complete; ram shape 141× surface on `Contains`; ALICE3 part mesh 12–62× faster than surface; timing under load flagged `timingPreliminary` |

**Sharpened during Track 3, quote it correctly:** `cyl_inter_cyl`'s mesh is *closed as a
triangle set* (0 odd-parity rays over 60k). The X-ray's 6 lost rays are a **navigation** loss in
`O2Tessellated`'s stepping. "The mesh leaks" means navigation, not holes — which is the stronger
statement for the talk.

## Open, in the order I would take them

1. **Converter: the JIT namespace bug.** `geom.C` emits `LoadSurfaceSolid` as a namespace-scoped
   forward declaration; the JIT wrapper makes it resolve wrongly, the macro fails to compile, and
   **the simulation continues silently without the module**. Workaround committed
   (`integration_demo/patch_exact_macro.py`); fix the converter's macro emission and make a
   failed module load loud.
2. **Converter: `--csg auto` without ROOT yields CSG 0 loudly-but-recoverably; parallel runs
   race.** The per-part WARN exists now; still open: serialize/lock the deferred emit (two
   parallel conversions lost shapes), and score **oTOF through the oracle gate** (converts,
   never scored).
3. **The CSG path forward is measured** (`Stream_AA_FlatCSG.md`): Tier-0 canonicalisation → the
   single-cell emitter (ships `tube_window`/`cyl_inter_cyl` as CSG, retires the two sliver rays)
   → composite trees for the 16 ALICE3 solids; `TGeoBVHCSG` deferred behind §5's falsifiable
   criteria. Decline reasons ship in the reports and `website_data/decline_reasons.json`.
4. **Kernel-confirmed: `ST1829909_01` leaks parity through real inter-face slits** (8/15 299
   rays, BVH ≡ Loop, no trim flags; `maxSharedEdgeDeviation` 4.70e-4 cm ≈ the declared model
   tolerance — identity closure is gap-blind by design). Needs the per-face localiser, then
   widened crossing acceptance near identified shared edges or trim-snapping to canonical
   curves. Probes in scratch: `probe_parity.cxx` / `probe_dev.cxx`.
5. **Kernel: loose bounding boxes.** `Stick` claims 40.1 cm in z around an 11.0 cm solid (3.6×,
   verified vs OCCT). Use the cover-box union in `ComputeBBox`.
6. **TGeo→STEP follow-ups** (`Stream_AD_TGeoToStep.md`): bisect the full-geometry
   `STEPCAFControl_Writer` segfault / try per-detector `--top` writes; map `TGeoPara` (13
   volumes); reflected placements of subtree-carrying volumes (98 dropped); localise the
   1-in-1.69M Contains point. On the record: the **beam pipe places 81 volumes exactly on
   identical copies of themselves** (TGeo-side; `--dedup-world` handles it; independently
   reproduced by a second walk, and `CheckOverlaps` is *structurally blind* to it — verified:
   1 overlap reported in all of PIPE, none a plie — because the coincident discs sit in sibling
   `TGeoShapeAssembly` wiggles, whose daughters the checker never cross-compares; materially
   benign since the pairs are identical volumes, which is why a decade of physics never noticed); `BRepAlgoAPI_Fuse`
   can return `IsDone()` with a silently dropped operand (guarded by a volume invariant);
   **46/106 PIPE CSG declines are exactly "a TGeoPcon"** — PIPE's 58 Pcons are the measured
   corpus for Roadmap (h)'s revolved-profile detector.
7. **Track 3b, the real data:** the MCStepLogger→`events.json` exporter (schema in
   `website/README.md`; tree layout in `integration_demo/data/README`), then the website's event
   tab replays the real IRIS/Bagger transports.
8. **Materials matching:** 26/55 IRIS volumes matched, 29 on vacuum `Default`, ten of the 26 by
   accidental string prefix. Anchored part-number matching.
9. **Gate: credit the sliver** — relabel a distout mismatch whose candidate crossing is
   `onTrimBoundary`-flagged within the trim band (Review Appendix A).
10. **The face-normal gate column** (fourth hand-over in a row) and **the `auto`-mode
    unreliable-shipping policy** (~20 lines) — both pre-corpus items.
11. **The talk (Track 4):** assemble from Review + Stream_Z + the website; re-run the
    `timingPreliminary` numbers on a quiet box; build the single-file website bundle for
    publishing (the Artifact CSP allows nothing external; the bridge stays local-only).
12. **In flight:** `O2BVHAssembly` (Roadmap (j)) — an Opus agent is implementing it
    (baseline-first, TGeoShapeAssembly contract, oTOF acceptance); fold its
    `Stream_AE_BVHAssembly.md` on landing.
13. **Standing, unchanged:** free-form surfaces; the models-are-not-legal overlap finding and the
    broken `CheckOverlaps` on our shape; `Curve2D::closestPoint` as the kernel hot spot; mesh
    healing.

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
