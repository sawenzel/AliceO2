# Stream AL — the closure test: real physics through the round-tripped geometry

**2026-08-26/27.** Executes [`Handoff_ClosureTest.md`](Handoff_ClosureTest.md): `o2-sim -m PIPE
ITS TPC MAG` run twice, once on the hand-written C++ TGeo geometry and once on that same geometry
taken out to STEP and back, with a fixed seed and `SimCutParams.trackSeed=true`, and compared.
Everything below was run, not relayed; the scripts are in `scripts/geometry/ClosureTest/`.

**The one-line verdict:** the round-tripped geometry presents the same material to the same rays
as the hand-written one — **mean x/X₀ 27.904215 vs 27.904907 over 2000 rays to r = 800 cm, a
median per-ray relative difference of 2.6e-12, and the same 246.5 volumes crossed per ray** — and
hit-by-hit agreement of charged tracks is *not* achievable, for the reason the handoff predicted
rather than any defect.

---

## 1. What had to be built first

Three things the round trip did not carry, each of which would have silently invalidated the test.

**Media (handoff pain point i).** Every converted volume got `TGeoMaterial("Default", 0, 0, 0)` —
transparent. `O2_TGeoToCAD.py --media-json` now dumps each volume's medium and material verbatim,
keyed by the emitted STEP part name, and `O2_CADtoTGeo.py --media-json` rebuilds them.
`ClosureTest/check_media.py` is the acceptance test. **746 of 746 volumes carry their source
medium, with radiation- and interaction-length deviation exactly 0.**

Two traps inside that. `TGeoMaterial::SetRadLen` **recomputes** from the recipe instead of storing
what it is handed, so writing the source's own X₀ back makes it *wrong* (35.7585 cm became 36.6601
on `MAG_WATER`); carrying the recipe alone reproduces both derived lengths bit-for-bit. And a
mother's material lives in its `__body` leaf, not on its assembly label, so the table needs both or
every mother returns transparent.

**Cuts and processes (pain point iii).** `loadCutsAndProcessesFromJSON` resolves an entry by
`getMediumID(module, local index)` and skips a mismatch **silently**. The handoff's suggested
module-key rename therefore cannot work: the CAD side numbers its media in its own traversal order.
`ClosureTest/remap_cuts.py` maps by medium *name*, which both sides preserve exactly, and the driver
runs a CAD probe pass first to learn those indices. **95 media matched by name, 95 identical in cuts
and processes, 0 unmatched.**

**Controls (pain point ix).** Two baselines at one seed give 1243 ITS hits with every Δ exactly 0.
A different seed moves the mean |Δr| to 34.7 cm, so the instrument can see a difference.

## 2. Four defects that only real transport could find

Each was invisible to the per-solid gate, and NEXT.md's own caveat says why: *every "0
disagreements" number on this branch is a statement about solids in isolation.*

1. **Mixtures were flattened.** `remapCADMedia` called `MaterialManager::Material` for everything,
   so a `TGeoMixture` lost its element composition and kept only an effective A and Z. **55 of the
   90 media in these four modules are mixtures** — the TPC drift gas, every beam-pipe alloy,
   `CAVE_Air`. Fixed (`98bcc4469e`).
2. **One material per medium.** The same function registered a material per medium, so a
   field-free `_NF` twin sharing its material tripped `MaterialManager`'s uniqueness assert and
   aborted the run. Fixed (`98bcc4469e`).
3. **The JIT namespace defect** (NEXT.md item 4): `geom.C` declared its loaders inside
   `namespace o2`, which the JIT wrapper nests and shadows, so the module did not load and
   `o2-sim` exited 0 anyway. Fixed (`3ef0e5f048`).
4. **Anchoring and nesting** — §3.

## 3. Anchoring and nesting: what the round trip drops

**A module has no single mother.** It adds nodes straight into the hall, so "hook the top volume in
at the right place" has no single answer: PIPE attaches at **10** places across `cave`, `barrel` and
`caveRB24`, ITS at **33**, TPC and MAG at **1**. `ClosureTest/module_anchors.py` reads those
attachment points and matrices out of the source geometry. Everything under `barrel` is one
conversion placed back into the real `barrel` with the identity; PIPE's `RB26Pipe` and `RB24C` are
converted from their own roots and placed with their own matrices (both a plain 180° about y, each
verified by rebuilding the Euler triple through ROOT and comparing).

Anchoring a whole module at `cave` instead does not work, and fails *silently*: its detectors are
then daughters of `cave` while sitting geometrically inside the native `barrel` solid, and the
navigator, having descended into `barrel`, never looks at them. Measured — a ray from the IP crossed
`barrel`, then `cave`, and nothing else, and transport fell to **161 steps per event against the
baseline's 36 004**.

**STEP cannot express a solid that contains other solids.** XCAF makes a label either a simple shape
or an assembly, so the writer emits a mother as an assembly holding its own `__body` leaf *beside*
its children. Read back literally that is a full-size solid overlapping its own daughters. TGeo
carves a mother implicitly — a daughter takes precedence over its mother's solid — so the converter
has to put the nesting back, and now does, from a `bodyOfAssembly` table in the sidecar. Measured
without it: a ray crossed `ITSUWrapVol0__body` from r = 2.1 to r = 16.4 cm and met **no ITS layer at
all**, and the run produced zero hits.

**`--carve-mothers` is not a substitute** and was tried first. It subtracts daughter *solids*, and
an assembly daughter has none; `ITSUWrapVol0`'s only daughter is an assembly, so it came back
uncarved. Nesting is also exact and costs no boolean, where carving costs one per mother and lost
four of ITS's outright (`SpaceFrameVolumeLay3..6`, "OCCT result contains no solid").

## 4. The result

**Material budget, the closure measure.** `ClosureTest/matbudget_diff.py` integrates x/X₀ and
x/λ along a fixed Fibonacci ray set through two geometry files. No RNG, so it asks the geometry
question and nothing else. Directions are a Fibonacci sphere, never an axis raster — that trap has
reported a real defect as clean twice in this project.

| 2000 rays, r ≤ 800 cm | mean x/X₀ | median rel. diff | rays > 1 % | rays > 10 % | volumes crossed |
| --- | --- | --- | --- | --- | --- |
| baseline | 27.904215 | — | — | — | 246.5 |
| CSG cascade, `TGeoTessellated` | 27.904907 (+0.002 %) | 2.6e-12 | 5 / 2000 | 0 | 246.5 |
| **CSG cascade, `O2Tessellated`** | **27.903624** | **2.6e-12** | **3 / 2000** | **0** | **246.5** |
| tessellated only, `TGeoTessellated` | 91.780165 (+229 %) | 7.4e-01 | 1994 / 2000 | 1923 | 81.3 |
| **tessellated only, `O2Tessellated`** | **27.903700** | **1.9e-05** | **4 / 2000** | **0** | **249.8** |

The two `TGeoTessellated` rows are the defect of §5, kept because they are what a naive
tessellated fallback costs. The `O2Tessellated` rows are what the pipeline ships.

Restricted to the ITS region (500 rays, r ≤ 45 cm) the cascade agrees to a mean absolute difference
of **1.3e-10** on 0.0647, i.e. floating point. ITS sensor placements are identical detector-wide:
108 / 144 / 180 / 2688 / 3360 / 8232 / 9408 for layers 0–6 on both sides.

**Hits.** Both worlds transport cleanly — zero stuck tracks, zero G4Exceptions, zero navigation
errors, zero aborts — at 24 689 steps and 3437 secondaries per event against 36 004 and 3347. The
primaries are bit-identical. Where both sides record a crossing the positions agree to ~0.03 cm.
But **only 10 of 198 primaries keep every hit within 1 cm**, and that is the handoff §1 answer, not
a defect: Geant draws from the RNG once per step, the round-tripped world has a different boundary
structure, and a track that takes one extra boundary step decorrelates from there on. Per-track
seeding confines the divergence to the track; it cannot remove it.

Two things confound a naive hit comparison and must be stated separately:

- **The two sides do not define a hit the same way.** ITS's own `ProcessHits` versus the external
  detector's built-in entrance/exit action gives a **1.2–2.7× multiplicity difference per layer**,
  largest in the inner barrel. Handoff pain point (iv), as predicted.
- **A track index is not a shared name.** The two runs keep different secondaries, so matching on it
  compares unrelated particles — a median 23 cm residual between geometries whose layers agree to
  the millimetre. `compare_hits.py` gained primary-ordinal keying and nearest-neighbour matching so
  both traps are visible rather than silent.

## 5. The tessellated-only variant, and the defect it exposed

Asked for as a benchmark. In its first form it did not close at all -- 3.3x the radiation
length while crossing 81.3 volumes per ray instead of 246.5, and in the ITS region 3.6 against
169.2, i.e. the inner detector was simply invisible. The cause is neither chordal error nor, as
this document first claimed, any inability to express a cavity:

> **`TGeoTessellated` does not navigate.** It derives from `TGeoBBox` and overrides none of
> `Contains`, `DistFromInside`, `DistFromOutside` or `Safety`, so every volume built on one is
> navigated as its filled **bounding box**. ROOT says so in the class doc -- *"The class does not
> provide navigation functionality, it just wraps the data for the composing faces"* -- and real
> navigation exists only through `TGeoVGShape`, which needs a VecGeom build this install does not
> have.

The minimal reproducer is a **tetrahedron**: four facets, correct outward winding,
`IsClosedBody()` true, and `Contains(0.9, 0.9, 0.9)` true for a point nowhere near it. A convex,
cavity-free, closed body is already wrong, which retires the cavity hypothesis outright. (A
cube-with-a-cavity is a bad reproducer twice over: a cube *is* its bounding box, and the obvious
hand-winding of the inner shell is wrong.)

**The emitted mesh was never at fault.** `IP_PIPE`: 1008 facets over 504 unique vertices, all 3024
directed edges present exactly once with their reverse, Euler V-E+F = 0 as a genus-1 shell must be,
0 degenerate facets, and a divergence-theorem signed volume of 82.988 cm^3 against the analytic
83.0227 -- a 0.04 % chordal deficit. An inward-wound inner shell would have read 1931.17 and a
missing one 1007.09. The decisive measurement: over 200 000 random points,
`TGeoTessellated::Contains` and a `TGeoBBox` of the same half-lengths and origin disagreed **0
times**.

**This reached the simulation, not only this document's instrument.** O2 runs Geant4 with
`G4Params::navmode = kTGeo`, and `TG4RootSolid` forwards `Inside`/`DistanceToIn`/`DistanceToOut`
straight to the `TGeoShape` virtuals, so Geant4 saw the same boxes.

**Why no gate caught it.** `runOracleGate.py`, `tgeoRayService.py` and `O2SolidHarness` all score
`o2::base::O2Tessellated`, which does navigate. The emitted `geom.C` used ROOT's class. The
instrument and the product were testing different shapes -- a sharper form of NEXT.md's standing
caveat that a per-solid number says nothing about an assembled world.

**The fix.** `O2_CADtoTGeo.py --mesh-solid {o2,tgeo}`, default `o2`, emits
`o2::base::O2Tessellated` through the existing `LoadFacetSolid`; `tgeo` keeps ROOT's class for a
macro that must load outside O2 and warns that every such volume becomes its bounding box.
`TGeoGeometryUtils::TGeoShapeToTGeoTessellated` carries the same warning; it has no callers today.

**It also moved the cascade.** The one mesh-tier part in the whole test is `ARB8`, a twisted
`TGeoArb8` with four free-form faces: 26.73 cm^3 of M55J6K carbon navigated as an 875.49 cm^3 box,
a factor 32.7. Fixing it cut the cascade's mean absolute difference 3.2x and its worst ray 6.8x
(0.7007 -> 0.1027). That was this document's own open item 4.

**Cost.** 2000 rays take 1.47 s through the hand-written geometry, 1.66 s through the cascade and
1.76 s through the tessellated-only world -- **+20 % wall, +17 % per boundary crossing** over
hand-written CSG. `o2sim_geometry.root` grows 16.7 -> 38.6 MB, because `O2Tessellated` streams its
outward normals.

Two incidental findings. `LoadFacetSolid` reports 2 degenerate facets across the PIPE barrel meshes
(`RB24ValveMA1` index 1516, `RB24VMABCPirani` index 3827) that `TGeoTessellated::AddFacet` had been
dropping silently; both solids still close, so they are benign. And ROOT's
`TGeoTessellated::FlipFacets()` iterates `for (auto facet : fFacets)` **by value**, so it is a
no-op -- a bug `O2Tessellated` inherited verbatim when it was copied.

## 6. Reproduce

```bash
S=<studydir>                       # holds env_o2.sh and env_converter.sh
for m in MAG TPC PIPE ITS; do
  python3 scripts/geometry/ClosureTest/roundtrip_module.py $S $m csg,mesh
done
python3 scripts/geometry/ClosureTest/make_configs.py $S                       # cascade
python3 scripts/geometry/ClosureTest/make_configs.py $S --variant mesh \
        --out-prefix mesh_ --name CADMESH
bash    scripts/geometry/ClosureTest/run_closure.sh $S 20 424242
python3 scripts/geometry/ClosureTest/matbudget_diff.py \
        $S/run/base1/o2sim_geometry.root $S/run/cad/o2sim_geometry.root \
        --rays 2000 --rmax 800
```

## 7. Open after this

1. **Report `TGeoTessellated`'s missing navigation to ROOT** — as a feature request, not a bug:
   ROOT documents the limitation, it simply does not shout when a `TGeoTessellated` is put in a
   `TGeoVolume` without VecGeom, and silence there costs a factor 32 on a real part. ROOT's master
   already carries an implementation similar to `O2Tessellated`. Also worth reporting: the
   `FlipFacets()` by-value loop, which needs fixing in `O2Tessellated` too.
2. **A geantino/`o2-sim-evalmat` cross-check** of §4's numbers through the real transport, so the
   claim rests on two independent instruments rather than on `matbudget_diff.py` alone.
3. **Hit semantics** — a `sensitiveMacro` reproducing ITS's own `ProcessHits` would make the hit
   records comparable one-to-one; today only positions are.
4. **Tree closeness.** The round trip is material-equivalent but the tree is not identical: a
   mother with daughters returns as an assembly plus a `__body` volume where the source had one
   volume, and each module sits under one extra hall wrapper (338 347 nodes / 649 volume UIDs
   against 374 285 / 875). Collapsing an assembly that has exactly one `__body` child at the
   identity into a single volume should reproduce the source tree one-for-one. Until then the extra
   boundaries are extra steps, which is one of the reasons per-track hits cannot correspond.
5. **The remaining outlier rays** — 3 of 2000 above 1 % in the cascade, 4 in the tessellated world.
   A per-ray localiser would name the volume each crosses.
6. **CPU comparison under real transport** — §5 measures ray tracing only. The configuration is
   fair now, so a transport-level number can be taken.
