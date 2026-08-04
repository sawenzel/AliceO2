# Stream X — the Geant integration test: converted CAD geometry under real transport

**2026-08-04/05.** Executes [`Handoff_IntegrationTest.md`](Handoff_IntegrationTest.md): ALICE3 IRIS
and oTOF converted to their own `geom.C` with materials, injected through
`o2::passive::ExternalModule`, transported with Geant4 (`geomRoot`, i.e. **TGeoNavigator** — the
representation actually under test), once purely tessellated and once exact
(CSG + `O2BVHSurfaceSolid`), and compared. Everything below was run, not relayed; commands are in
`scripts/geometry/IntegrationTest/`.

**Definition-of-done scorecard (Handoff §2):** ① both detectors convert with materials — **done,
after fixing the oTOF traversal blocker**; ② both transport under `o2-sim` with **zero stuck
tracks, zero navigation errors, zero aborts** — done; ③ tessellated vs exact compared on steps,
secondaries, CPU, memory and material budget — done; ④ reproducible from committed scripts —
done; ⑤ this document — done. Failures and caveats are in §8; they are real but none blocks the
verdict.

**The one-line verdict:** the exact representation survives a real Geant4 transport
indistinguishably from the mesh on robustness, agrees with it on the material budget to
`8.6e-05 %X₀` mean over 21 600 directions, and costs **2.0× in geometry-dominated stepping and
~0 in a physics-dominated run** — not the 40–145× of the per-call microbenchmark.

---

## 1. Stage 0 — the oTOF blocker was one traversal assumption, now fixed (`685851ce76`)

The converter mapped **one XCAF simple-shape label to one leaf solid**. oTOF stores its entire
geometry in three simple-shape labels whose TopoDS shapes are *compounds*: `oTOF v2 v1` (3
solids, 1505 faces), `Module v1` (17 solids), and `A Side` (**empty** — 0 solids, 0 faces). So
the converter saw 3 "solids" where the STEP file places 62 628, and `triangulate_asbbox()` died
on the empty label with `Bnd_Box is void` — the handoff's suspicion that the crash was a symptom
of the traversal defect was correct.

Fix: a multi-solid compound expands to one leaf per contained solid under a synthetic assembly
node; a single-solid compound unwraps; an empty label is skipped on every visit (including
re-visits through the visited-defs early return, which would otherwise emit a placement edge to a
volume that was never registered). Pinned by a fifth `--self-test` block (8 checks; suite now
**56 checks, 0 failures**).

Measured after the fix, matching the plain-reader census exactly:

| | |
| --- | --- |
| oTOF | **20 leaf prototypes, 62 628 placements, all distinct** (17 Module solids × 3672 + 3 shell solids × 68) |
| ALICE3 `CAD_noETA.stp` | 55 definitions / 103 placements — unchanged |
| Bagger oracle gate | exit 0, surface 13/13, shape 7/7 — unchanged |

## 2. Two integration defects the instruments could never see (`848d952605`)

Both were found by loading the converted geometry **the way `ExternalModule` actually loads it**
(`buildCADVolumeFromMacro`), which the tessellated-only IRIS prior art never exercised:

1. **The exact `geom.C` did not survive namespace wrapping.** The emitted macro re-opened
   `namespace o2 { namespace base {...} }` to declare `LoadSurfaceSolid` by prototype.
   `buildCADVolumeFromMacro` wraps the macro body in a unique namespace, so that block became
   `<wrapper>::o2::base` and *every* later `o2::base::` reference resolved to it — the JIT failed
   on every exact-surface macro. The prototype now comes from the public
   `DetectorsBase/O2SurfaceSolidIO.h`; no namespace is re-opened.
2. **`--exact-surfaces auto` shipped sidecars the loader rejects.** Tutorial §5.6's "one live
   production risk", materialized: ALICE3's `ST1829909_01` sidecar fails `LoadSurfaceSolid`
   (wire-join gap 5.4e-06 cm against the fixed 1e-06 tolerance — Stream_L mechanism 3, the known
   open item), and `geom.C` then throws at build time, killing the whole module. The converter now
   **load-validates every emitted sidecar through the real O2 loader** (exactly the check
   `geom.C` runs) and demotes failures to the tessellated fallback. On ALICE3 exactly the two
   known parts demote (`ST1829909_01`, `ST1829909_004`); `required` mode fails instead.

## 3. Stage 1+2 — conversions, `IntegrationTest/convert.sh`

Mesh precision **0.5** everywhere meshing happens (measured-affordable — 0.1 reached 22.9 GB on
one sphere and was killed; sweep is a roadmap item, not this task). One deliberate deviation from
the handoff's flag table, recorded: the exact configs also pass `--mesh`, because a part the
cascade cannot carry exactly would otherwise fall back to a **bounding-box placeholder**, and
ALICE3 has 35 such parts — the comparison would have been meaningless.

| config | tiers (CSG / surface / mesh) | wall | peak RSS |
| --- | --- | --- | --- |
| ALICE3 exact | 2 / **16** (18 emitted − 2 demoted) / 37 | 30 s | 1.14 GB |
| ALICE3 mesh | 0 / 0 / 55 | 20 s | 0.86 GB |
| oTOF exact | **19 / 1 / 0 — fully exact** | 54 s | 0.76 GB |
| oTOF mesh | 0 / 0 / 20 | 71 s | 0.42 GB |

oTOF's 19 CSG parts are plain `TGeoBBox`es (`dV_sym = 0` exactly); the 20th (`oTOF v2 v1_s3`,
1493 planar faces) is one `O2BVHSurfaceSolid`. **The corpus NEXT.md item 6 called "the one that
would pay best" is now fully exact.**

**Materials.** ALICE3 via `IRIS/IRIS_MATERIALS.csv` + the G4 NIST DB: the BOM parses as intended
(21 entries, header rows skipped); **26 of 55 leaves carry a BOM material** (Al alloys, stainless,
carbon fiber, CuBe — 24 NIST-resolved, plus 2 on `St. Steel EN 1.4306 (304L)` which does **not**
resolve and ships as a named dummy with a FIXME), **29 fall back to `Default`** (A=0, Z=0, ρ=0 —
transparent). oTOF: one material for the whole module per the brief — **POLYSTYRENE**
(ρ 1.06, X₀ 41.3 cm; low-Z/few-secondaries, not air, not a stopper, physically plausible for a
TOF), resolved to `G4_POLYSTYRENE`, on all 20 leaves. The `G4_Si` contrast variant was not run
(time went to the oTOF fix); it remains the suggested second material.

## 4. Stage 3 — placement, measured not assumed (`56c4de4a4b`)

`IntegrationTest/cadModuleProbe.cxx` loads a `geom.C` exactly as `ExternalModule` does and
reports the world-frame bbox under a candidate placement:

- **ALICE3 IRIS is authored beam-along-+Y** (836 cm local Y extent — IRIS *with beam pipe*), CAD
  origin = IP. `rotation_deg [90,0,0]` + the barrel `y=+30` shift → ALICE frame
  x [−12.9, 8.7], y [−7.7, 34.2], z [−398.1, 439.6] cm. The transverse asymmetry is the CAD's own.
- **oTOF is authored beam-along-+X with a corner-anchored frame** (axis through local
  (y, z) = (85.330, 0)) — the handoff's "do not take RotateX(90) on faith" warning was earned:
  neither model needed the worked example's rotation unmodified. `rotation_deg [0,90,0]` +
  translation (0, 30−85.330, +173.805) → a centred barrel: x,y [−96.6, 96.6], z [−174.0, 174.0].
- Both fit far inside `barrel` (R 790.5, z ±714.6): the handoff's "might not fit" worry does not
  bite. The two modules' material shells do not intersect (oTOF starts at r≈85; ALICE3 ends at 34).

`o2-sim-serial -n 0`: both modules log `Configured external module`, media remap to the O2
MaterialManager, and the world closes at **66 697 nodes / 102 volume UIDs — identical for both
representations** (only the shapes differ). `CheckOverlaps(0.1 cm)`, as a hint only: mesh world
**0 overlaps in 21.9 s** — although OCCT instruments measured real 2.3 mm CAD interpenetrations in
this very model (Stream_T), so that 0 says more about the checker than the model; exact world
**does not finish in 300 s** (the known starvation, Handoff §8.3, unchanged). Overlaps stay
deprioritised per Sandro; reported, moved on.

## 5. Stage 4 — transport (seed 424242, TGeant4 + geomRoot, BoxGun π± ×10/event, p 0.1–5 GeV, |η|<1)

Every run: **exit 0, zero stuck tracks, zero looping-killed tracks, zero G4Exceptions, zero
navigation errors** — grepped, both representations, all eight logs.

| run | init (s) | event loop (s) | peak RSS (MB) |
| --- | --- | --- | --- |
| geomcheck exact / mesh | 26.0 / 11.2 | — | 925 / 743 |
| pions ×5 exact / mesh | 24.7 / 10.6 | 0.38 / 0.36 | 908 / 743 |
| steplog ×50 exact / mesh | 24.6 / 10.5 | 0.59 / 0.51 | 909 / 743 |
| matbudget (21 600 geantinos) exact / mesh | 24.8 / 10.5 | **3.77 / 1.90** | 908 / 743 |

The ~14 s init difference is the exact module build itself (sidecar loads + BVH construction +
JIT; probe: 16.6 s vs 2.3 s for ALICE3), paid once. Memory: +166 MB for the exact world.

## 6. Stage 5 — the comparison

**Material budget (the headline).** `o2-sim-evalmat`, one event = a geantino scan over
**60 η × 360 φ = 21 600 directions** (η ∈ [−1.5, 1.5] — a direction fan, not an axis raster; the
direction-poor-beam trap is real and bitten this project twice). Identical grid, identical seed,
only the representation differing — `compare_matbudget.py`:

| map | mean (exact) | mean abs diff | max rel diff | bins >1 % | bins >10 % |
| --- | --- | --- | --- | --- | --- |
| X₀ (hradl, %) | 4.5613 | **8.6e-05** | 6.5 % | 1443 / 21600 | **0** |
| g/cm² (hgcm2) | 109.71 | 3.4e-03 | 7.9 % | 1620 / 21600 | **0** |

The two representations present the same material to the same rays at the few-1e-5 level in the
mean; 6.7 % of individual directions differ by 1–8 %, which is the mesh's chordal error at
precision 0.5 made visible under transport. The instrument can show a difference — it did — and
none exceeds 10 % anywhere. (Quirks, same in both files: `habso` is all-zero and `rzA` holds NaNs;
`MaterialBudgetMap` artifacts, not representation effects.)

**Steps and secondaries** (50 events, 500 π, MCStepLogger; per-volume table via
`summarize_steplog.sh`):

- total steps: **82 965 (exact) vs 82 894 (mesh) — −0.09 %**; tallied-in-module volumes 40 306 vs
  40 240. Per volume, the large volumes agree to a few %: `oTOF v2 v1_s1/s2` within 1 %,
  `ST1A38494_002/003/004` within 2 %.
- secondaries: 1972 vs 2148 (+8.9 % mesh) and tracks 2278 vs 2450. The bulk of the difference sits
  in `ST1A38494_01` (1081 vs 1280) — a volume that is **the same mesh in both runs**, so this is
  shower divergence (tiny boundary differences decorrelate the RNG streams), not a representation
  effect. The same caution applies to `ST2487458_01` (0 vs 47 steps; same mesh both runs).
- the one *genuinely* different-representation volume with a visible difference is
  `oTOF v2 v1_s3` (the `O2BVHSurfaceSolid`): steps 3253 vs 3102 (+4.9 % exact), secondaries 56 vs
  10 — small absolute numbers; with N=500 primaries this is indicative, not significant. A
  higher-statistics run is cheap (event loop is 0.6 s) if a number is ever quoted from it.

**Cost.** Where geometry dominates (geantino scan), exact/mesh = 3.77/1.90 = **2.0×**. Where
physics dominates (π transport), 0.59/0.51 = **1.16×**, i.e. within noise of free. The
representation benchmark's 40–145× per-call gap does not survive contact with a real transport,
where geometry queries are one cost among many and the navigator caches aggressively. **No
optimisation was performed while measuring.**

## 7. What this unlocks (and did not measure)

This harness is the instrument Roadmap (b)/(c) need — the step-length cost of a loose safety and
the sub-surface BVH payoff are now *measurable* (compare step counts/CPU under a knobbed
`Safety`). Neither knob exists yet; nothing was measured on them. Assembly-level transport under
`TGeoNavigator` now *runs* (the Roadmap's "will be exercised in the Geant test automatically"),
but a leaking *counter* still does not exist — zero robustness incidents is evidence, not proof,
of no leaks.

## 8. Failures, caveats, and what was not done — stated plainly

- **29 of 55 ALICE3 leaves transport with a `Default` (A=0, Z=0, ρ=0) material** — they are
  transparent. The material-budget equivalence is unaffected (both representations share the
  gap), but the absolute X₀ map understates the real detector. The 304L steel string resolves to
  nothing in the NIST matcher and ships as a named dummy on 2 leaves.
- **Mesh validity remains unhealed** (MeshHealing.md): the mesh world drops 137–175 degenerate
  facets at load, and `ST1A38494_01` reports not-fully-connected facets. No transport incident
  followed at these statistics, but "no incident at N=500" is weak evidence; the leaking counter
  (§7) is the real answer.
- **`CheckOverlaps` on the exact world does not terminate in 5 minutes** — unchanged from the
  handoff; the WIP in `00cc7d9eb3` was not touched. The OCCT-side census (Stream_T) remains the
  only trustworthy overlap instrument, and the full ALICE3 1699-pair census remains unrun.
- **The steplog totals ("did N steps") and the per-volume tally sums differ** (82 965 vs 40 306):
  the printed tally covers the volumes the logger lists per event; the totals line is
  authoritative for whole-world sums. The per-volume table is still valid volume-by-volume.
- **The oTOF `G4_Si` second variant and the mesh-precision sweep were not run**; the pion sample
  is 500 primaries (adequate for robustness and totals, thin for per-volume secondaries);
  MCStepLogger's ROOT-tree output (`mcStepAnalysis`) was not used — text tallies sufficed.
- **Environment findings that cost time here** (all reproduced, none guessed): the converter's
  pythonOCC/`PYTHONPATH` prepends **segfault o2-sim at startup** — keep converter and sim
  environments separate (`convert.sh` header vs `run_sim.sh` header); a stale `o2-sim-serial`
  against a freshly rebuilt `libO2DetectorsBase` segfaults even on `--help` — rebuild executables
  with the library; `ninja … | tail` hides a failed build behind tail's exit 0 — the first sim
  build died silently on a **full /tmp** (the machine's root filesystem hit 100 %; ~3.4 GB of
  regenerable gate scratch from Aug 1–2 was deleted to proceed: `/tmp/{xa3*,q_xa_*,tier0,streamE,
  S_xbag,J_bigfan*,q_xt_before,q_gate2_before,I_off_bag,gateF_*,csg_bag2,new_fix}`).

## 9. What this stream invalidates elsewhere (left for the next session to reconcile)

- `NEXT.md` item 6 ("the oTOF corpus is unreachable from the converter") — **fixed**; oTOF is now
  the best-covered corpus on the branch (fully exact).
- `NEXT.md`'s "two sidecars that do not load" (Stream_L mechanism 3) — the *tolerance defect*
  stands, but its blast radius changed: auto mode now demotes instead of shipping a
  build-time throw. Any text saying ALICE3 "emits 20" should read "emits 20, ships 18" until the
  tolerance item is fixed.
- Tutorial §5.6's "exact part ships with only a warning" production risk — **closed** by
  load-validation (for load failures; the reliability-flag half of that sentence is untouched).
- The Roadmap deferral "assembly-level transport … will be exercised in the Geant test
  automatically" — now true; the missing piece is only the leaking counter.
- `--self-test` baseline: 48 → **56 checks**.

## 10. Reproduce

```bash
# converter env (pythonOCC); see Handoff §11 — prepend, never replace
bash scripts/geometry/IntegrationTest/convert.sh <workdir>
bash scripts/geometry/IntegrationTest/make_sim_configs.sh <workdir>
# sim env: alienv + stage-first LD_LIBRARY_PATH + VMCWORKDIR=<srcdir> +
# ROOT_INCLUDE_PATH=<srcdir>/Detectors/Base/include — and WITHOUT the converter's prepends
for r in exact mesh; do for m in geomcheck pions steplog matbudget; do
  bash scripts/geometry/IntegrationTest/run_sim.sh <workdir> $m $r
done; done
python3 scripts/geometry/IntegrationTest/compare_matbudget.py \
  <workdir>/sim/matbudget_exact/o2sim_matbudget.root <workdir>/sim/matbudget_mesh/o2sim_matbudget.root
bash scripts/geometry/IntegrationTest/summarize_steplog.sh \
  <workdir>/sim/steplog_exact/run.log <workdir>/sim/steplog_mesh/run.log
```
