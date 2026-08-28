# Plan — the ALICE Upgrade Week talk, Sardinia, mid-September 2026

Written 2026-08-28, agreed in session. `Plans.md` next to this file is Sandro's statement of
*what* he wants to present; this file is the plan of record for *how* the supporting artefacts get
made, and who each slide's evidence is. It supersedes `../Plan_Presentation.md` for this talk —
that file described an earlier ~30 min WG slot and is kept for its track structure.

## Decisions taken

- **Slot**: ~20 min + 5 discussion. That is **15 slides**, ~75 s each; one claim and one artefact
  per slide. Everything else is backup.
- **All flagship artefacts are wanted** — renders, tables, benchmarks and the ALICE3 example.
  There is no single must-have, so the plan is ordered by risk instead: the risky work starts on
  day 1 and degrades gracefully.
- **Date**: mid-September, i.e. **~2 weeks** from 2026-08-28.
- **oTOF materials come from AliceO2 source**, not from a Mattermost round-trip. Checked in
  session: `Detectors/Upgrades/ALICE3/IOTOF/simulation/src/Detector.cxx::createMaterials()`
  defines exactly two — a four-component `AIR$` mixture and a plain `SILICON$` material
  (A=28.086, Z=14, rho=2.33, X0=9.36, lambda=999) — so this is an afternoon, not a blocker.
- **The hero render's subset is not pre-chosen.** W3 produces a census of candidate subsets with
  their representation histograms and cheap previews, and Sandro picks. See W3 step 1.
- **Deliverable is artefacts, not a deck.** Every figure, table and number as files, plus a
  slide-by-slide script. Sandro assembles the deck.

## Two findings that already de-risk the plan

1. **The corpus report survives.** `~/roundtrip_report/alice_roundtrip.json` (30 MB) carries
   `parts` and `modules` records for all 22 901 leaf solids of the R6 breadth census, alongside
   `alice_roundtrip.md/.html` and the 17 module `.step` files. The recognition scoreboard is
   therefore a **re-formatting** job, not a 17-module re-run (MFT alone was ~90 min). The
   per-module `conv/` directories are gone; anything not in that JSON does need a re-run, and
   W2 step 1 establishes exactly what is and is not in there before any slide depends on it.
2. **oTOF's materials are trivial** (above). The remaining oTOF risk is the XCAF traversal and the
   fact that it converts but has **never been scored through the oracle gate**
   (`Stream_AC_OTOFTraversal.md`), not the physics inputs.

## The 15 slides

| # | slide | the one claim | artefact | workstream |
| --- | --- | --- | --- | --- |
| 1 | Title | — | — | — |
| 2 | Why CAD -> simulation | quicker R&D, no hand-written C++ geometry | recap figure from the TB talk | — |
| 3 | What changed since the TB talk | one row per capability gained | delta table | W2 |
| 4 | The representation cascade | tessellated -> exact surfaces -> CSG / FlatCSG, chosen per part | concept figure | W1 |
| 5 | `O2BVHSurfaceSolid` | a solid is its trimmed surfaces plus a BVH; navigation is exact | simple-case render | W1 |
| 6 | `O2FlatCSG` | union of intersection-cells over signed halfspaces, BVH over sub-boxes | simple-case render | W1 |
| 7 | Truth by inversion, and recognition | TGeo -> CAD first, so every part has a known answer; recognition in five sentences | pipeline figure | W2 |
| 8 | Recognition scoreboard | 22 736 / 22 901 native CSG (99.3 %); 163/188 composite-sourced (87 %); every primitive class 100 % | three tables | W2 |
| 9 | **Hero render, before / after** | all-grey tessellated -> coloured and labelled by representation | two Blender renders | W3 |
| 10 | Full-system closure test | Run 3 PIPE/ITS/TPC/MAG round-tripped: median per-ray x/X0 difference 2.6e-12 | plot | W4 |
| 11 | What exactness costs | all-tessellated vs mostly-exact: transport, memory, disk, observables | chart | W5 |
| 12 | Composing a detector, and the website | many CAD sources via JSON, sensitive actions via JIT macros; the website benchmarks representations | code snippet + screenshots | W6 |
| 13 | **ALICE3 end to end** | beam pipe + IRIS + oTOF -> hits in IRIS and oTOF | hit counts, event display | W7 |
| 14 | Side story: MFT | 16 965 of 19 185 volumes are duplicates; a geometry doctor should catch this | figure | W8 |
| 15 | Gaps and roadmap | media / field / cuts assignment is the open pain point; free-form surfaces; arrangement cells | text | W9 |

**Backup slides**: the safety benchmark (W10); the "models are not legal worlds" overlaps finding
(`Stream_T_AssemblyOracle.md`); the tolerance policy; `TGeoTessellated` does not navigate and what
that cost (`NEXT.md`, closure test).

## Workstreams

Each one is a self-contained investigation with its own gate. A workstream that fails its gate
degrades its slide to an honest statement of where it stopped; it never silently disappears.

### W2 — recognition scoreboard (1 day, build-independent) — slides 3, 7, 8

1. Establish what `alice_roundtrip.json` actually contains: schema of `parts` and `modules`, and
   which of the quoted headline numbers are reconstructible from it alone.
2. Emit three slide-sized tables to `tables/`: the module summary; the primitive-class feature
   matrix (source class vs what the round trip made of it); the decline census split into
   **converter** declines and **writer** declines (`TGeoPara` unmapped 13, six unbounded
   composites, six hollow tori, prism section-count mismatches).
3. Emit the slide-3 delta table: capability by capability, TB talk vs now.
4. Every number lands in `numbers.json` with the file and record it came from.

**Gate**: no number on a slide that is not traceable to a record in `numbers.json`. Numbers that
cannot be reconstructed are listed as needing a re-run, and the re-run is scheduled or the claim
is dropped.

### W9 — the media / field / cuts pain point (1 day, build-independent) — slide 15

Read `remapCADMedia` in the converter, the media handling in the closure test
(`ClosureTest/check_media.py`, `remap_cuts.py`) and `Stream_AL_ClosureTest.md`'s media result
(746/746 volumes carried their source medium, 95/95 media agreed on cuts and processes), and write
down precisely what works today, what a CAD-only model has no source for, and what the user must
supply by hand. Deliverable `notes/media_gap.md` plus three bullets for the slide.

**Gate**: the slide names a mechanism, not a worry.

### W7 — ALICE3 end to end (4 days, highest risk, starts day 1) — slide 13

1. **Materials** (day 1, build-independent): lift oTOF's two media from `Detector.cxx` into
   `oTOF_MATERIALS.csv` in the shape of `IRIS/IRIS_MATERIALS.csv`.
2. Convert `STEP_examples/oTOF System V3-R92cm.step` and **score it through the oracle gate** —
   never done before.
3. Convert the ALICE3 beam pipe and IRIS (`ALICE_3_example/VTX_noETA.root` route; convert ALICE3
   **without** `--mesh`; never default `--mesh-prec 0.1`, a 2 m sphere reached 22.9 GB).
4. Compose all three as `externalDetectors` in one JSON (README.md section) with a JIT sensitive
   action per sensitive volume.
5. `o2-sim-serial`, fixed seed; count `IRISHit` and the oTOF hits in `o2sim.root`.

**Gates, in order**: the gate scores oTOF -> all three convert -> `o2-sim` transports with zero
stuck tracks and no navigation errors -> both detectors produce hits. Stopping at gate 3 means
slide 13 shows the composed geometry and states the blocker.

**Traps**: keep the converter env and the sim env separate; run `--csg auto` conversions strictly
serially (parallel runs race and lose shapes); the JIT namespace bug means a failed module load is
silent — apply `integration_demo/patch_exact_macro.py` and check the module is really there.

### W1 — concept renders for the two new solids (2 days) — slides 4, 5, 6

Two synthetic parts chosen to be legible at slide size: a cylinder-intersect-cylinder for the BVH
solid, drawn as trimmed patches with the BVH boxes over them; and one flat-CSG part drawn as
coloured cells with the sub-box BVH. Rendered through `tgeoRayService.py` (the **real** kernel) or
the website raytracer. Plus the cascade figure for slide 4.

**Trap**: restart the ray bridge — a bridge left running from before 2026-08-25 traces composites
multi-threaded and every picture it made of one is wrong.

**Gate**: legible at one third of slide width, and the picture matches the kernel (the website's
engine-diff view is bit-true against it).

### W3 — the hero render (4 days, the visual centrepiece) — slide 9

1. **Candidate census**: from `alice_roundtrip.json` and W7's ALICE3 output, list subsets of 20-60
   parts whose representation mix is visually interesting (surface-exact, CSG-exact, tessellated
   all present), with the histogram and a cheap preview per candidate. **Sandro picks.**
2. **Placement recovery** for the chosen subset via the STEP route of
   `O2DPGAgentic/doc/detector-cad-rendering.md` — the mesh dump's repeated leaf volumes are
   unplaced by design, and matching by volume then clustering by radius is the technique that
   works.
3. **BEFORE**: every part tessellated, uniform grey, camera chosen so faceting shows on
   silhouettes.
4. **AFTER**: same camera, colour by representation, surface-attached labels via
   `render_its.py`'s centroid-projected labeller.

**Gate**: identical camera in both images; every label readable without zooming; the colour legend
has at most four entries.

### W4 — closure-test plot (1 day) — slide 10

One figure from `Stream_AL_ClosureTest.md`'s existing data: per-ray relative x/X0 difference over
2000 Fibonacci rays on a log axis, the five rays above 1 % visible and named. Caption states the
honest part — hit-by-hit agreement of charged tracks is not achievable, and why.

### W5 — the cost benchmark (2 days, needs an idle machine) — slide 11

Re-run `integration_demo`'s exact-vs-tessellated pair and chart transport factor, one-off build
time, memory, file size, and the material-budget agreement. The existing numbers (1.94x transport,
+42 s build, +174 MB, 0.039 % over 512 rays) are flagged `timingPreliminary` in `NEXT.md` and must
be re-measured before they go on a slide.

**Gate**: nothing else running on the box; one `ninja` at a time rule applies.

### W6 — website and composition (1 day) — slide 12

Serve `website/`, restart the ray bridge, capture the per-part feature matrix and a live `/bench`,
and build a single-file bundle so the demo survives a bad network. Plus a minimal worked
`externalDetectors` JSON and JIT sensitive-action macro as a code snippet.

### W8 — MFT duplication (half a day) — slide 14

One figure from the numbers already in `NEXT.md` / `MFT_deduplication_result.md`, and the framing
sentence about an automatic geometry doctor.

### W10 — safety benchmark (backup, 2 days) — cut first

Accurate safety vs fast-return plus safety caching. The one item that may need kernel work rather
than measurement, which is why it is backup.

## Schedule

| days | work | gate before moving on |
| --- | --- | --- |
| 1-2 | W2, W9, W7 step 1 (all build-independent) | tables emitted and traceable; `oTOF_MATERIALS.csv` exists |
| 2-5 | W7 steps 2-5 (serial conversions), W1 alongside | oTOF scored; three modules composed |
| 4-8 | W3 (long pole; needs W7 if an ALICE3 subset is picked) | subset picked by day 5 |
| 6-8 | W5 on a quiet box, W4, W8 | benchmark numbers no longer `timingPreliminary` |
| 8-10 | W6, artefact ledger | single-file website bundle runs offline |
| 10-12 | `Script.md`, dry run against the clock | fits 20 min |

**Cut order** if time runs short: W10, then W6's live demo (screenshots remain), then W1's second
concept render, then W7 degrades to geometry-without-hits. **W2 and W3 are never cut** — they are
the talk.

## Deliverables

```
TalkUpgradeWeek2026/
  Plan.md          this file
  Plans.md         Sandro's statement of what to present
  Script.md        slide by slide: claim, artefact filename, speaker note
  numbers.json     every headline number with the file and record it came from
  figs/            PNG and PDF, one per slide artefact
  tables/          markdown and CSV
  notes/           per-workstream investigation records
```
