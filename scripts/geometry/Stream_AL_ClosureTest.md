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
| **CSG cascade** | **27.904907** (+0.002 %) | **2.6e-12** | 5 / 2000 | **0** | **246.5** |
| tessellated only | 91.780165 (+229 %) | 7.4e-01 | 1994 / 2000 | 1923 | 81.3 |

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

## 5. The tessellated-only variant, and why it fails

Asked for as a benchmark. It does not close, and the reason is **structural, not chordal**:

> **`TGeoTessellated` has no notion of an internal cavity.** A closed two-shell body reads as
> filled.

Measured on the converted beam pipe: `IP_PIPE` comes back with 1008 facets and 504 vertices whose
radii run 1.82 to 1.90 cm — both shells present — and `IsClosedBody()` is true, yet
`Contains(0,0,0)` and `Contains(1,0,0)` are **both true**, inside the vacuum bore. Reproduced from
scratch on a hand-built 24-facet cube-with-a-cubic-cavity: `Contains` is true inside the cavity.

Every hollow volume therefore becomes solid, which is why the tessellated world crosses **81.3
volumes per ray instead of 246.5** while carrying **3.3× the radiation length**. In the ITS region
it crosses **3.6** volumes per ray against 169.2 — the inner detector is simply invisible. Mesh
precision is not the lever; the default 0.1 was used and refining it cannot add a cavity the shape
class does not represent.

This is `MeshHealing.md`'s "a mesh can be *invalid*, not merely inaccurate" in its sharpest form,
and it bounds the tessellated fallback for any TGeo-derived geometry: it is usable for solids
without cavities and not otherwise.

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

1. **`TGeoTessellated` cavities** — decide whether to report upstream, and meanwhile refuse to emit
   a tessellated fallback for a solid with an inner shell rather than emitting one that is wrong.
2. **A geantino/`o2-sim-evalmat` cross-check** of §4's numbers through the real transport, so the
   claim rests on two independent instruments rather than on `matbudget_diff.py` alone.
3. **Hit semantics** — a `sensitiveMacro` reproducing ITS's own `ProcessHits` would make the hit
   records comparable one-to-one; today only positions are.
4. **The five outlier rays** above 1 % in the cascade column are unnamed. A per-ray localiser would
   say which volume they cross, and one of them is probably PIPE's single tessellated part.
5. **CPU comparison** — not measured. The configuration is fair now, so it can be.
