# Result — the capacitor/welding/X7R0402 half of `MFT_deduplication_task.md`

The "cheap 15 432-volume win" this task doc called out — `capacitor`, `welding0`,
`welding1` — is done, verified, and pushed, plus a follow-on addendum: the `X7R0402`
assembly containing them was duplicated the same way (5144 fresh `TGeoVolumeAssembly`s)
and missed by this doc's own table, since it only tabulated plain box-shaped volumes. The
sensor chain (`MFTSensor`/`MFT_C`/`MFT_G`) is **still untouched**, exactly as this doc's
item 4 scoped it: `GeometryTGeo::extractMatrixSensor()` and
`Detector::addAlignableVolumesChip()` both hard-code `MFT_C`'s volume *name* as
`MFT_C_<half>_<disk>_<ladder>` for navigation and CCDB alignment attachment — deduplicating
that family (or the ladder assembly `MFT_L`, same problem) needs that path construction
rewritten first, and is a separate task.

## Where

- **PR open**: [`AliceO2Group/AliceO2#15723`](https://github.com/AliceO2Group/AliceO2/pull/15723)
  against `dev`, squashed to one commit, head `1610003e` on fork branch
  **`swenzel/mft-electric-component-dedup`** (note the `swenzel/` prefix). The original
  working branch `mft-electric-component-dedup` (no prefix, 3 separate commits ending
  `65ec4bc9c9`) is superseded — **branch further work off the PR's head, not that one.**
- Full plan, verification scripts and raw measurement data:
  `/home/swenzel/alisw/agent-workspace/mft-tgeo-dedup/` (this machine) — see `RESULT.md`
  there for the full writeup, including which further targets are ready to implement,
  which need real design work first, and which need a scope decision.

## The numbers

`capacitor`/`welding0`/`welding1`/`X7R0402`: 5144 volumes/shapes each → 1 each. `MFT_G`/
`MFT_C`/`MFTSensor` unchanged at 280/280/280. Total `TGeoManager` volumes: 25853 → 10424
(phase 1) → 5281 (phase 2).

Verified after every commit: a 200-event seeded run (`SimCutParams.trackSeed=true`) gives
bit-for-bit identical `MFTHit`s (8850/8850, all fields) before/after. Containment-checked
for all three deduplicated leaf families (200 sampled points each), with a negative control
proving the check can actually fail — which caught a real unrelated segfault in the
verification script itself along the way.

Bonus, not asked for in the original doc but worth knowing: geometry-initialization wall time
(MFT-only `o2-sim-serial -n 0`) drops from a median 54.2 s to 41.6 s at phase 1 — about 23%
faster, confirmed with a run-order-reversed control batch. Not re-measured after phase 2's
further ~50% volume-count cut; plausibly faster still.

## Still open

- Fixing the sensor-chain half (`MFTSensor`/`MFT_C`/`MFT_G`) or the ladder assembly
  (`MFT_L`) needs the `MFT_C`/`MFT_L` navigation paths in `GeometryTGeo.cxx`/`Detector.cxx`
  rewritten to not depend on the volume's name embedding `<half>_<disk>_<ladder>` — see the
  Global Constraints section of the full plan
  (`/home/swenzel/alisw/agent-workspace/mft-tgeo-dedup/2026-08-25-mft-tgeovolume-dedup-plan.md`)
  for exactly where those call sites are.
- `connectord`/`boxconnectord`/`varnishlayer`/`flex`/`lineslayer`/`alulayer`/`kaptonlayer`
  are confirmed safe (same check as everything above), just not yet turned into a task —
  roughly another 2218 volumes off the current 5281, down to ~3063. `MetalStack` is also
  safe but shared with ITS via `AlpideChip.cxx` — needs an explicit scope decision, not
  just an MFT-local one; another ~279 on top, to ~2784.

Whoever picks any of this up next: read `RESULT.md` first — the verification recipe
(geometry summary + hit diff + containment check, all scripted) should carry over directly.
