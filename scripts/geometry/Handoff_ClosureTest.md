# Handoff — the closure test: real physics through the round-tripped geometry

**Written 2026-08-23 for a fresh session to execute.** Sandro's charge, verbatim:

> *"We have a running ALICE run3 simulation using C++ TGeo. It is producing hits from events.
> Say we take: `o2-sim -m PIPE ITS TPC MAG --seed FIXED + trackSeeding=ON`. What about the
> following: (a) we convert the exact geometry to STEP. (b) we do the backconversion in 3
> possible ways: everything tessellated, everything CSG, everything surfacesolid (c) we use the
> external modules/detectors mechanism and simulate the same setup as above but with the new
> geometry. The goal is to see if we get the same hits (ITS), if the step number is under
> control, if everything goes through etc.; We can then also gain some realistic CPU benchmarks."*

This is the experiment the whole programme has been building toward: the same physics through
the same geometry expressed four ways, with the original C++ TGeo as its own exact oracle.
Everything below is the analysis of how to do it honestly; nothing here is executed.

Prerequisite state (all on the branch, all measured): the TGeo→STEP writer round-trips PIPE,
ITS, TPC, MAG individually at ~6.2 M `Contains` samples with one disagreement total
(`Stream_AD/AF/AG/AH/AI`); writer defects (i)–(iii) and the WSEG collinear face are fixed
(self-test 105); the external modules/detectors mechanism is exercised end to end with hits
(`Stream_Z_IntegrationDemo.md`); MCStepLogger instrumentation and the comparison discipline
exist there too.

---

## 1. The central methodological point: what "the same hits" can mean

Bit-identical hits across representations is **not the right acceptance for charged particles in
material**, and chasing it would fail for a reason that is not geometry: Geant4 consumes the RNG
stream per step (multiple scattering, δ-ray and interaction sampling), and a representation that
takes even one extra boundary step — which tessellation does by construction — shifts every
subsequent draw of that track. Two geometrically identical worlds then produce diverging
trajectories after the first stochastic sampling. Per-track seeding contains the divergence to
the individual track (that is why the charge says trackSeeding=ON) but does not remove it within
a track.

So the acceptance is **layered**, strongest claim first:

- **L0, exactness where exactness is definable — geantinos.** No physics, no RNG: the ordered
  volume-crossing sequence and the Σx/X₀ integral per ray must agree between baseline and each
  representation within stated tolerance (the `Stream_Z` fan discipline: Fibonacci directions,
  never axis rasters; both positive controls — a bit-identical repeat, and a deliberately
  degraded mesh that must move the number).
- **L1, first-contact identity — charged primaries up to the first sensor crossing.** With
  per-track seeds and identical media/cuts/field, a primary's path to its FIRST ITS hit is
  deterministic transport in identical material; position/entry momentum of that first hit
  should match to navigation tolerance for tracks whose path crosses only
  representation-identical geometry. Diverging tracks must be *attributable* (a crossed region
  where representations differ, or an upstream stochastic process) — count and attribute, do
  not average away.
- **L2, distributional equality — full events.** Hit counts per chip/layer, dE/dx spectra,
  secondaries, step counts per volume: statistically compatible, with the comparison powered
  enough to detect the deliberately-degraded control.
- **Robustness gate throughout**: zero stuck tracks, navigation errors, aborts; step numbers
  "under control" = per-volume MCStepLogger tallies with named outliers, not one aggregate.

## 2. The three back-conversions, named honestly

"Everything CSG" is not currently achievable and the hand-over must not pretend it is: ITS
recognises 128/235, TPC 82/172, and the two missing recognisers (sloped prisms, Pcon revolved
profile — `Stream_AA`/round-trip studies) plus `TGeoHalfSpace` (5 TPC composites, 732
placements) bound today's ceiling. The honest variants:

| variant | cascade | what it measures |
| --- | --- | --- |
| T | tessellated only (`--mesh`, precision chosen and justified per §5.6) | the fallback everyone else uses |
| C | CSG → surface → mesh (`--csg auto --exact-surfaces auto --mesh`) | the shipped cascade, maximally native-ROOT |
| S | surface → mesh (`--exact-surfaces auto --mesh`) | the exact-surface path under real transport |

Optionally, building the two recognisers FIRST (they have six exactly-known corpora now) would
raise C's coverage substantially — a scoping decision for the fresh session, not a prerequisite.

## 3. Pain point (i) — tracking media, and in-field/out-field

The medium parameters (`isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin`) live on
each `TGeoMedium` and control per-medium stepping in field; the field itself is global
(`o2::field::MagneticField`), so "in-field vs out-field" is carried by the medium flags, not by
geometry. **Route: extract, do not re-derive.** Extend the TGeo→STEP step (a) to dump, per
volume, the medium + material record verbatim (name, Z/A/ρ/radlen/intlen, the 8 medium
parameters) into a sidecar JSON keyed by the emitted STEP part name (the writer already walks
every volume; this is bookkeeping, not geometry). The generated `geom.C` gains a mode that
builds media directly from this sidecar — exact reproduction, no matching heuristics.

## 4. Pain point (ii) — materials and the CSV route

Two routes, used for two different purposes:
- **Exact-dump route (the fair comparison)**: §3's sidecar. The physics comparison is only fair
  if the media are bit-equal; NIST matching noise would masquerade as geometry effects.
- **CSV/NIST route (the user-facing path test)**: additionally generate the BOM CSV from the
  same dump and run the existing `--materials-csv` + NIST machinery once, reporting match rate
  and worst density/composition deviations — exercising the route Sandro wants exercised,
  *without* letting its deviations contaminate L0–L2 (those run on the exact-dump build).
  Known bug to avoid or fix first: the accidental-prefix matching (NEXT list).

## 5. Pain point (iii) — Geant cuts and processes

**The mechanism exists and is verified present**: `MaterialManager::writeCutsAndProcessesToJSON`
/ `loadCutsAndProcessesFromJSON`, driven by `MaterialManagerParam.outputFile` / `.inputFile`
(`--configKeyValues`). Plan: baseline run dumps the cuts/processes JSON; the external-geometry
runs load it. **The open question to resolve early**: keying. Cuts are registered per (module,
medium); the external modules carry different module names, so either the loader tolerates a
name map, the dump is rewritten with the external names, or global+special cuts are applied at
medium granularity. Acceptance: after loading, re-dump from the external run and diff — the two
JSONs must agree modulo the intended renaming, or the comparison is not fair and must not run.

## 6. Additional pain points the charge did not list

- **(iv) Hit semantics.** The real ITS produces `o2::itsmft::Hit` with chip IDs through
  `GeometryTGeo` and ITS's own stepping logic; an external detector produces generic
  `o2::ext::Hit` with its own sensor indexing and the built-in entrance/exit action. Identical
  *physics* ≠ identical *digitised* hit records. Comparison layer: match on position/energy/
  track, or port the ITS stepping logic into a custom `sensitiveMacro` (the JIT hook exists).
  TPC sensitivity is a different animal (per-electron ionisation) — keep TPC passive in v1, ITS
  is the hit detector, exactly as the charge says.
- **(v) World identity.** Anchoring (`barrel`, y = +30), rotation, and the mother/air volumes:
  verify by bbox/placement-fidelity walk (the Stream_AH fourth instrument) BEFORE any transport;
  a module that does not land where the original detector sat invalidates everything downstream
  silently.
- **(vi) Baseline self-coincidences.** The branch carries the PIPE and TPC geometry fixes; the
  baseline must be run from THIS branch so both sides share them. `RB26s3Bellow` remains empty
  on both sides (recorded, material-changing, untouched) — fine for fairness, note it.
- **(vii) Writer scale.** Combined PIPE+ITS+TPC+MAG is ~42 k components — inside the known
  segfault bracket (38 676 fine / 74 601 crash). **Write per-module STEPs** (four files, four
  external modules); it avoids the bracket entirely and matches the mechanism.
- **(viii) Tessellation precision vs sensor scale.** Facets are float32 and `--mesh-prec` is an
  angular knob in disguise (`Stream_P` §5); ITS sensors are ~50 µm structures. Pick the
  precision against the L0 material-budget acceptance, per module, and state it; do not reuse a
  cosmetic value.
- **(ix) Determinism control first.** Before comparing representations, run the baseline twice
  with the same seed and require bit-identical hits — that verifies the seeding mechanism
  (`SimConfig::startSeed`; the exact per-track flag to be confirmed on the running version) and
  is the positive control everything else stands on.
- **(x) CPU benchmark hygiene.** Same cuts/media/field established first (else the ratio
  measures configuration, not geometry); single-threaded, load average recorded,
  `timingPreliminary` discipline; geometry build/close time reported separately from transport.

## 7. Staging for the fresh session (suggested order)

0. Determinism control (§6-ix); baseline run with cuts dump (§5) and media dump (§3).
1. Writer: media/material sidecar + `geom.C` exact-media mode (small, testable, self-test).
2. Per-module STEP writes (a) + three back-conversions (b) + placement-fidelity gate (§6-v).
3. Cuts round trip (§5) with the re-dump diff as the gate.
4. L0 geantino closure, with both controls.
5. L1/L2 with ITS hits via external detector; step tallies; robustness counters.
6. CPU benchmarks last, on the surviving fair configuration.
7. `Stream_*_ClosureTest.md` with the numbers; partial-with-named-blockers beats
   unfinished-everything (the standing rule).

**Definition of done for v1**: L0 passes on all three variants for all four modules; L1
attribution exists (even if imperfect); one honest L2 table; the pain points (i)–(iii) resolved
by mechanism, not by hand-editing; every number reproducible from committed scripts.

## 8. Why this is worth it beyond closure

If L0–L2 pass, the sentence for the talk and the paper is: *"the converted geometry is
physics-equivalent to the hand-written one, measured through the full simulation chain"* — which
no per-solid gate can say. The CPU table (T vs C vs S vs native) is then the first fair,
realistic cost comparison of the representations under production-like transport. And the
media/cuts extraction machinery built for it is exactly what any future external user of the
CAD pipeline needs anyway.
