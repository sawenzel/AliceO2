# Stream Z — the Geant integration demo: IRIS and Bagger through `o2-sim`

This is Track 1 of [`Plan_Presentation.md`](Plan_Presentation.md), executing
[`Handoff_IntegrationTest.md`](Handoff_IntegrationTest.md) with Bagger substituted for the
blocked oTOF. Two CAD models are converted twice each — once through the exact cascade
(CSG → `O2BVHSurfaceSolid` → tessellated fallback) and once purely tessellated — injected into
`o2-sim` as **sensitive** external detectors, and compared under real Geant4 transport.

It works. Both models build inside `barrel`, Geant4 transports geantinos, electrons and pions
through them with no stuck track, no navigation error and no aborted event, and both detectors
produce real hits. The two representations agree on the material budget to **0.039 %** over 512
Fibonacci directions, and every one of those rays takes the identical number of boundary
crossings in both. Pushed to 8192 rays and 302 757 steps the agreement breaks in exactly one
place — `BucketLink2`, whose mesh ROOT itself calls not closed. The exact path costs **1.94× in
transport time**, **+42 s of one-off geometry build** and **+174 MB resident**.

Getting there required one workaround, because the exact path could not be loaded into `o2-sim`
at all before today.

Everything below was measured on this machine on 2026-08-22, branch `swenzel/bvhsurfacesolid`,
and is reproducible from `scripts/geometry/integration_demo/`.

---

## 1. The blocker: exact geometry could not reach `o2-sim`

`O2_CADtoTGeo.py::emit_cpp_prelude(exact_surfaces=True)` writes the loader prototype as a
namespace-scoped forward declaration:

```cpp
namespace o2 { namespace base {
bool LoadSurfaceSolid(const std::string& file, O2BVHSurfaceSolid& solid);
} }
```

`o2::base::loadCADGeometryHook` (`Detectors/Base/src/CADGeometryUtils.cxx`) JITs the macro
inside a unique namespace, hoisting only lines that begin with `#` to global scope. That turns
the declaration above into `o2_cadgeom_IRIS_0::o2::base`, and every subsequent
`o2::base::O2BVHSurfaceSolid` in the macro then resolves against that nested `o2` instead of the
global one:

```
input_line_120:45:61: error: no type named 'O2BVHSurfaceSolid' in namespace 'o2_cadgeom_IRIS_0::o2::base'
[ERROR] Failed to JIT external geometry macro .../iris_exact/geom.C
[ERROR] No geometry could be built for external detector IRIS
```

The simulation then **continues silently with the module missing** — a one-line `[ERROR]` in a
250 000-line log, and a geometry file that is 11 KB instead of 3.9 MB.

`integration_demo/patch_exact_macro.py` replaces the declaration with
`#include "DetectorsBase/O2SurfaceSolidIO.h"`. Being a `#` line it is hoisted to global scope,
where `o2::base` is the real one; the header is installed and the macro already adds its include
path. It is a **workaround in the demo pipeline, not a fix** — the converter is not mine to
change on this branch, and the right repair is one of two one-line changes upstream (emit the
`#include`, or make the wrapper hoist the declaration too).

Two smaller mechanism notes, both measured:

- The passive path logs `External geometry JSON ... must contain an 'externalModules' array` on
  every run that configures only `externalDetectors`. It is harmless, but it is an `[ERROR]`,
  and it sits three lines above the two `Configured external detector` lines that actually
  matter.
- `o2-sim-serial` does **not** migrate external-detector hits to their own file.
  `macro/migrateSimFiles.C` iterates DetIDs the GRP marks as read out and looks their branch
  names up in `SimTraits::DETECTORBRANCHNAMES`, which knows nothing about external detectors.
  The hits are there — branches `IRISHit` and `BAGRHit` in `o2sim.root` — but no
  `o2sim_HitsTST.root` appears. In parallel mode `O2HitMerger` writes them normally.

## 2. Conversion, four ways plus one control

`integration_demo/convert_all.sh` produces all of them. IRIS is
`ALICE_3_example/CAD_noETA.stp`, Bagger is `STEP_examples/Bagger.step`.

| output | flags | wall | peak RSS | tiers (CSG / exact / mesh) | on disk |
| --- | --- | --- | --- | --- | --- |
| `iris_exact` | `--csg auto --exact-surfaces auto --mesh --mesh-prec 0.5` | 24.6 s | 839 MB | 0 / 20 / 35 of 55 | 3.7 MB sidecars + 12 MB meshes |
| `iris_tess` | `--mesh --mesh-prec 0.5` | 15.8 s | 844 MB | — | 12 MB, 345 517 triangles |
| `bagger_exact` | `--csg auto --exact-surfaces auto` | 2.3 s | 226 MB | 0 / 13 / 0 of 13 | 220 KB sidecars |
| `bagger_tess` | `--mesh --mesh-prec 0.1` | 1.1 s | 251 MB | — | 1.6 MB, 43 822 triangles |
| `iris_coarse`, `bagger_coarse` | `--mesh --mesh-prec 2.0` | — | — | — | 85 206 / 2 960 triangles |

**The IRIS "exact" run is a cascade, not a pure exact conversion, and it needs `--mesh`.**
Only 20 of 55 leaf definitions emit an exact sidecar; without `--mesh` the remaining 35 fall
back to `triangulate_asbbox`, i.e. to their bounding boxes, which is not a geometry anyone can
transport through. Both IRIS arms therefore carry the *same* 35 meshes at the *same* precision,
and only 20 volumes actually differ between them. That is the honest scope of the comparison,
and it shows up cleanly in the per-volume tallies in §7.

**`--mesh-prec 0.5` for IRIS**, chosen and not swept. The default 0.1 is unaffordable (a 2 m
sphere reached 22.9 GB and was killed, `NEXT.md`); 0.5 is the finest value previously measured
to fit, and it did: 839–844 MB peak, under 25 s. Bagger is small enough for the default 0.1.
`--mesh-prec 2.0` is the deliberately-degraded control of §8.

**Two recorded baselines have moved.** `Handoff_IntegrationTest.md` §13 records Bagger as
*CSG 7 / surfaces 6 / tessellated 0*; today's cascade gives **CSG 0 / surfaces 13 /
tessellated 0**. Every part is declined by the CSG recogniser, and for the six cylinders the
declared reason is the empty string (`declined CSG: None`), which is at least a reporting
defect. The same file records `--self-test` at 48 checks; it now runs **12**. Neither is mine to
chase — both are Track 0's to reconcile — but both are quoted numbers that no longer hold.

*(Reconciled by the integrating session, same day: neither baseline moved. The 48-check
self-test is four suites (18+8+10+12) and 12 is the last suite's own line — re-verified at 48
this morning. The CSG 0 is an environment artifact: this run's converter shell had the pythonOCC
env but not the O2/ROOT env, so the CSG hook's deferred shape completion (`emit.py --from-json`)
could not run and every candidate declined; the gate pipeline, which composes both envs, measured
Bagger at CSG 7 with all shape columns clean the same day. Two genuine items survive this
reconciliation and are on `NEXT.md`'s list: the empty decline reason (`declined CSG: None`), and
that a direct `--csg auto` invocation silently yields zero CSG instead of erroring when ROOT is
unreachable.)*

## 3. Materials

**IRIS** uses `IRIS/IRIS_MATERIALS.csv` with `g4_nist_database/G4_NIST_DB.json`. The CSV parses:
21 BOM entries are read from a file whose first ten rows are report headers. Of the 55 logical
volumes, **26 receive a real medium and 29 fall back to `Default`** — a `TGeoMaterial` with
A = Z = ρ = 0, i.e. vacuum. Media emitted: Alu EN AW-5083 (H116) ×7, Alu EN AW-6082 (T6) ×2,
Alu EN AW-5083 (O-H111) ×1, Stainless Steel ×11, St. Steel EN 1.4306 (304L) ×2,
Cu Be C17410 (TH02) ×2, Carbon Fiber ×1.

Of those 26, **ten are matched by accident.** `build_volume_to_material_map` falls back to a
substring test on the part number, and the BOM's `ST0923290_01` is a substring of the volume
names `ST0923290_010` … `ST0923290_019` but not of `ST0923290_002` … `_009` or `_020` … `_028`.
All twenty-seven are sub-solids of the same gate-valve assembly; ten of them got Stainless Steel
and seventeen got vacuum, decided by nothing but string prefixing. Steel is plausible for all of
them, so the effect here is benign — but the mechanism is not, and it is silent.

The BOM also lists `ST2513437_01, SILICON SENSORS IRIS TRACKER-B0-B1-B2, Silicium - Silicon`,
which does **not** appear in `CAD_noETA.stp`. There is no silicon anywhere in the converted
model; it is the mechanical and beam-pipe assembly.

**Bagger** gets one material for the whole model via `integration_demo/Bagger_MATERIALS.csv`, a
thirteen-row BOM in the converter's own format naming `Stainless Steel`, which resolves to
`G4_STAINLESS-STEEL` (ρ = 8 g/cm³, three elements). Steel is what an excavator is made of, and
for a few events at 1 GeV it is neither transparent nor a stopper.

## 4. Placement, measured rather than assumed

`Cave.cxx` places `barrel` at cave (0, −30, 0), so a module that is to land on the ALICE origin
carries translation y = +30. That part is as the handoff says.

The rotation is not. The handoff suggests `RotateX(90)` on the assumption that the CAD is
Z-up. Both models are in fact authored with their long axis along CAD **y**: IRIS measures
dx = 21.8, dy = 837.7, dz = 43.1 cm, and its beam-pipe volume `SOLID` is centred on CAD x = z = 0
and spans y ∈ [−390, +390], so CAD y **is** the beam and CAD y = 0 **is** the interaction point.
`TGeoCombiTrans::RotateX(+90)` maps local +y → master +z and local +z → master −y (verified
numerically, not read off), so `rotation_deg [90,0,0]` is right — for the right reason.

Placements, and where the modules actually land, read back out of the assembled world
(`check_geometry.C` on `o2sim_geometry.root`):

| module | placement in `barrel` | world box, ALICE frame | leaves |
| --- | --- | --- | --- |
| IRIS | translation (0, 30, 0), rotation (90, 0, 0) | x [−12.93, 8.89] y [−8.92, 34.22] z [−398.10, 439.61] | 103 |
| BAGR | translation (120.93, 102.06, 20.33), rotation (90, 0, 0) | x [94.82, 105.19] y [−36.08, 36.08] z [−37.05, 37.05] | 13 |

IRIS straddles the origin with its beam axis on (0, 0) as the CAD says. Its y extent is
asymmetric because the IRIS support base hangs to one side in the CAD; that is the model, not
the placement. The Bagger's translation was chosen to centre its bounding box on ALICE
(100, 0, 0): clear of IRIS, well inside the barrel, and close enough that a box generator at the
origin reaches it.

**Both fit.** `barrel` is a `TGeoPgon` of radius 790.5 cm spanning z ∈ [−714.6, 714.6] minus two
end pockets, the nearer of which starts at z = −505. IRIS reaches z = 439.6 and r = 36.6; the
Bagger reaches r = 111. The "might not fit" worry does not bite.

## 5. Does the geometry assemble

Identical structure in both arms: **90 volumes, 138 nodes**, max depth 7.

| shape class | exact | tessellated |
| --- | --- | --- |
| `o2::base::O2BVHSurfaceSolid` | 33 | 0 |
| `TGeoTessellated` | 35 | 68 |
| `TGeoShapeAssembly` | 19 | 19 |
| `TGeoBBox` / `TGeoPcon` / `TGeoCompositeShape` | 1 / 1 / 1 | 1 / 1 / 1 |

## 6. Overlaps — reported, not chased

`TGeoManager::CheckOverlaps(0.1)` on the assembled world:

| representation | overlaps |
| --- | --- |
| exact | 4 |
| tessellated | 0 |

All four are inside one IRIS sub-assembly, `ST1829909_005`, between `ST1829909_01` and its
`_002` / `_003` / `_004` siblings; the largest is 2.20 cm. **None** between the two modules, and
none against any ALICE volume.

Read this as a hint. The checker is documented-unreliable on `O2BVHSurfaceSolid`
(`Handoff_IntegrationTest.md` §8: starved by default, and wrong when fed), and there is a second
candidate explanation that has nothing to do with the checker — a mesh at precision 0.5 shrinks
its part, so a genuine interpenetration of the CAD solids can simply vanish from the tessellated
arm. `Stream_T_AssemblyOracle.md` already found these source models are not legal. Per Sandro's
standing instruction these were reported and left alone.

## 7. Transport

`integration_demo/run_all.sh`, `o2-sim-serial`, seed 42, single-threaded, MCStepLogger attached.
Five events each of 50 geantinos, 20 electrons at 1 GeV, 20 π⁺ at 1 GeV, plus one event of a
512-ray Fibonacci geantino fan.

**No stuck track, no navigation error, no aborted event, in any run, in either representation.**

| run | init s | transport s | peak RSS MB | steps | secondaries | IRIS hits | BAGR hits |
| --- | --- | --- | --- | --- | --- | --- | --- |
| geantino exact | 45.21 | 0.438 | 773.2 | 9 410 | 0 | 0 | 0 |
| geantino tess | 3.60 | 0.419 | 599.4 | 9 410 | 0 | 0 | 0 |
| electron exact | 46.04 | 0.435 | 773.6 | 6 765 | 519 | 85 | 2 |
| electron tess | 3.64 | 0.419 | 599.0 | 6 431 | 459 | 88 | 0 |
| pion exact | 46.39 | 0.531 | 773.8 | 6 506 | 420 | 94 | 0 |
| pion tess | 3.62 | 0.417 | 599.4 | 6 044 | 341 | 89 | 0 |
| matfan 512 exact | 46.31 | 0.464 | 835.1 | 18 917 | 0 | 0 | 0 |
| matfan 512 tess | 3.58 | 0.417 | 783.9 | 18 917 | 0 | 0 | 0 |
| bigfan 8192 exact | 44.35 | **1.488** | 837.3 | 302 757 | 0 | 0 | 0 |
| bigfan 8192 tess | 3.42 | **0.767** | 788.1 | 302 761 | 0 | 0 | 0 |

The two 8192-ray rows were timed without MCStepLogger; their step counts come from a repeat of
the same two runs with it attached, which costs about 25 % in transport (1.812 s and 1.079 s) and
50 MB more resident, and is why the timing pair and the counting pair are separate runs.

**Geometry build.** The exact cascade costs **+42 s** once per job (45.6 s against 3.6 s,
averaged over the six runs of the table). Timestamped inside the log,
IRIS's builder hook takes **43 s** against **3 s** tessellated, and Bagger's takes about 1 s
either way — so essentially all of it is IRIS's 20 exact solids being constructed and closed,
about 2 s each. The tessellated arm builds 55 meshes and 13 more in 3.4 s.

**Transport.** On the 8192-ray fan, where geantinos do nothing but cross boundaries, the exact
arm takes **1.94×** the tessellated arm's transport time (1.488 s vs 0.767 s for ≈302 760 steps:
4.9 µs vs 2.5 µs per step). **55 %** of those steps are in `cave` and `barrel`, which are the
same volumes in both arms, leaving 137 762 steps in CAD volumes; attributing the whole difference
to those gives **+5.2 µs per CAD step**, and if a `cave`/`barrel` step costs about what a mesh
step costs, roughly **3× per CAD step**. That second number rests on an assumption and is
labelled as an estimate.

Either way this is the point the representation benchmark could not answer: `O2BVHSurfaceSolid`
measured **40–145×** a `TGeoCompositeShape` on isolated queries
(`Stream_P_RepresentationBench.md`), and under a real transport that becomes **1.94×** on the
job and about **3×** on the steps that touch it.

**Memory.** +174 MB resident for the exact geometry (773 vs 599 MB) in the small runs. On disk
the exact representation is *smaller*: Bagger is 220 KB of sidecars against 1.6 MB of triangles.

**Hits are real.** With `sensitiveMedia: ["Carbon Fiber"]` IRIS resolves to exactly one volume,
`ST2487721_01` (`CENTRAL PIPE-IRIS-2ND VACUUM`), which is one of the 20 that convert exactly —
so the hits come out of the exact kernel. With `sensitiveVolumes: ["Bucket"]` the Bagger's
substring match selects the whole bucket group (`Bucket`, `BucketLink1`, `BucketLink2`,
`BucketCylinderInner`, `BucketCylinderOuter`). Geantinos produce no hits, correctly: the
built-in action records charged tracks only.

The hit radii are the most direct picture of the difference anywhere in this document. In the
electron run:

```
exact       IRISHit n=85  r[ 5.460, 5.460]
tessellated IRISHit n=88  r[ 5.427, 5.460]
```

Every exact hit sits at r = 5.460 cm, the true cylinder radius. The tessellated hits are spread
over 5.427–5.460 cm — the sagitta of the faceted polygon, **330 µm** on a 5.46 cm vertex-detector
pipe at mesh precision 0.5.

### Per-volume steps

Geantinos, 5 × 50, exact vs tessellated:

| volume | representation | steps exact | steps tess | path exact cm | path tess cm |
| --- | --- | --- | --- | --- | --- |
| `barrel` | (same) | 4 750 | 4 750 | 195 600.808 | 195 600.686 |
| `ST1A38494_01` | mesh in **both** | 1 367 | 1 367 | 23.5170 | 23.5170 |
| `ST1A38494_004` | exact vs mesh | 494 | 494 | 9.0794 | 9.0603 |
| `ST1A38494_002` | exact vs mesh | 482 | 482 | 8.5693 | 8.5316 |
| `ST1A38494_003` | exact vs mesh | 460 | 460 | 8.2134 | 8.1749 |
| `ST2487455_01` | exact vs mesh | 460 | 460 | 13.6166 | 13.5794 |
| `ST2487721_01` | exact vs mesh | 422 | 422 | 12.5307 | 12.7848 |
| `SOLID` | mesh in **both** | 422 | 422 | 12.9032 | 12.9032 |
| `cave` | (same) | 371 | 371 | 249 863.338 | 249 863.338 |
| `ST2487457_01` | mesh in **both** | 176 | 176 | 1.4697 | 1.4697 |
| `Boom` (Bagger) | exact vs mesh | 6 | 6 | 3.9567 | 3.9577 |

This is the internal consistency check the comparison needed: **the step counts are identical
volume by volume**, and the path lengths differ in exactly the volumes whose representation
changed and are bit-identical in the ones that did not. The three IRIS volumes that are meshes
in both arms (`ST1A38494_01`, `SOLID`, `ST2487457_01`) agree to the last digit, as they must.

Electrons and pions **do** diverge — +5.2 % and +7.6 % more steps in the exact arm, +13 % and
+23 % more secondaries. This is **not** a geometry-equivalence signal and should not be read as
one: a charged track in a 5 kG field that is deflected by a 330 µm surface difference on its
first crossing is a different track from then on, and the two runs stop being the same
experiment. Their value here is only that transport completes cleanly for both.

## 8. The material budget — the headline, with its controls

`integration_demo/fibonacci_geantinos.macro` fires 512 geantinos along Fibonacci-sphere
directions with **no random numbers at all**, so ray *i* is the same direction in every run.
`matbudget.C` integrates step length over the radiation length of each step's own medium,
taken from the geometry that run actually used.

| | 512-ray fan | 50-geantino sample |
| --- | --- | --- |
| rays whose direction differs between the two runs | 0 | 0 |
| **rays whose boundary-crossing count differs** | **0** | **0** |
| Σ x/X₀ exact | 42.385502 | 22.004069 |
| Σ x/X₀ tessellated | 42.402216 | 22.012385 |
| aggregate difference | **−0.0394 %** | −0.0378 % |
| per-ray \|Δ\|/X₀ median | 0.022 % | 0.036 % |
| per-ray \|Δ\|/X₀ worst | 3.30 % | 0.47 % |
| worst absolute Δ | 0.0019 X₀ | 0.0015 X₀ |

Two representations of the same CAD present the same material to 0.04 % in the aggregate and
0.02 % on the median ray, and — the stronger statement — **not one of 512 rays sees a different
number of surfaces**. That is the leak test: a mesh with `meshClosedBody = false` lets a track
walk through a wall, and that shows up as a missing crossing pair, not as a small path-length
error. There were none at 512 rays. The sign is consistent (tessellation slightly thicker in both
samples), which is what a faceted approximation of these particular parts gives.

**At 8192 rays it does break, in one place.** The same fan with sixteen times the directions
gives 302 757 steps exact against 302 761 tessellated, and the per-volume tally localises the
whole difference:

| volume | exact | tess | Δ | path exact cm | path tess cm |
| --- | --- | --- | --- | --- | --- |
| `BucketLink2` | 6 | 4 | **+2** | 1.9657 | 1.4501 |
| `barrel` | 152 936 | 152 942 | −6 | 6 287 703.483 | 6 287 700.067 |
| every other volume | | | 0 | | |

`BucketLink2` is the one Bagger part whose mesh `TGeoTessellated::Check` calls **not fully
connected** (§9). One ray in 8192 crosses it in the exact solid and does not in the mesh. Two
explanations fit and this measurement does not separate them: a genuine leak through the open
mesh, or a grazing ray that misses an inscribed facet it would have clipped on the true surface.
The accompanying six-step difference in `barrel` is not accounted for by either, and is recorded
without an explanation. The relative size of the whole effect is 1.3 × 10⁻⁵ of the steps.

The lesson is the one the handoff warned about: a per-ray equivalence result is only as strong as
the number of rays, and the place it eventually fails is the part that was already known to have
an invalid mesh.

**A refuted difference is worth nothing unless the instrument could have found one, so both
controls were run.**

- **Determinism.** `matfan_exact` repeated with the same seed reproduces Σ x/X₀ = 42.385502 and
  18 917 steps **bit for bit**. The 0.0394 % is not run-to-run noise.
- **Sensitivity.** The same fan through the mesh-precision-2.0 control gives
  Σ x/X₀ = **47.912467**, **+13.0 %** against exact, with 19 039 steps instead of 18 917 and 13
  volumes touched instead of 12. The instrument moves by 330× the observed difference when the
  geometry is genuinely degraded, so the 0.04 % agreement at precision 0.5 is a result and not a
  blind spot.

The coarse control also says something on its own. Its excess is not spread evenly: the path
length inside `ST1A38494_01` goes from 46.66 cm to **70.45 cm**, +51 %, in a single volume. That
volume is one of the two whose mesh `TGeoTessellated::Check` calls **not fully connected** (§9).
A broken mesh does not merely lose material; here it invents it.

## 9. Mesh validity, measured inside `o2-sim`

Loading the geometry, ROOT itself reports:

| run | degenerate facets rejected | solids not fully connected |
| --- | --- | --- |
| exact (IRIS fallbacks only) | 121 | `ST1A38494_01` |
| tessellated | 175 | `ST1A38494_01`, `BucketLink2` |
| coarse (prec 2.0) | 175 | `BucketLink2` |

`ST1A38494_01` — the single busiest CAD volume in every geantino run, 1 367 of 9 410 steps — is
**not a closed body** at mesh precision 0.5, in *both* arms, because it is one of the 35 IRIS
parts that the exact cascade does not reach. And note the last row: at precision 2.0 that same
part happens to close, while `BucketLink2` does not. Mesh validity is not monotone in mesh
precision, which is exactly the trap [`MeshHealing.md`](MeshHealing.md) names — chordal accuracy
does not detect it and refining does not reliably cure it.

That no ray leaked at 512 directions is therefore a statement about *those* rays, not a clean
bill of health for the mesh — and at 8192 directions one ray did diverge, on `BucketLink2`,
which is one of the two flagged parts (§8).

## 10. `O2BVHSurfaceSolid` bounding boxes are loose, by up to 3.6×

Not something this task went looking for; it fell out of comparing the two Bagger conversions'
world boxes, which differ by 12 cm in z on a 72 cm model.

Taking OCCT's own `Bnd_Box` (written as the 12-triangle `facets_*.bin` of the no-mesh run) as the
oracle, and the precision-0.1 mesh as a second witness, the two agree to **< 0.1 mm** in local
coordinates on every part. The `O2BVHSurfaceSolid`'s own `TGeoBBox` does not:

| part | axis | true extent cm | `O2BVHSurfaceSolid` box cm | excess |
| --- | --- | --- | --- | --- |
| `Stick` | z | [25.856, 36.839] | [9.766, 49.855] | 40.1 cm for an 11.0 cm solid — **3.6×** |
| `BucketLink2` | z | [17.982, 21.737] | [10.138, 21.738] | 11.6 cm for 3.8 cm — 3.1× |
| `Boom` | y | [−60.538, −24.138] | [−64.009, −24.138] | +3.47 cm |
| `StickCylinderOuter` | z | [2.084, 15.769] | [1.465, 15.769] | +0.62 cm |
| `BoomCylinderOuter` | y | [−27.223, −12.925] | [−27.224, −12.304] | +0.62 cm |

The pattern — one axis, one end, on parts built from cylinders whose trimmed patch is a fraction
of the full surface — points at a box taken from untrimmed surface extents rather than from the
trim. It is conservative, so it is not a correctness bug: TGeo only requires box ⊇ solid. It
costs navigation (a looser box is a weaker rejection test) and it inflates `CheckOverlaps`, which
is a candidate explanation for §6's asymmetry. Recorded, not fixed.

## 11. What is not done

- **oTOF stays blocked**, as planned. Not attempted here.
- **`--mesh-prec` was chosen, not swept.** 0.5 for IRIS with 2.0 as a degraded control; the sweep
  remains a roadmap item.
- **The charged-particle comparison is not a controlled experiment** past the first surface
  crossing (§7). Making it one would need the field off and single scattering suppressed.
- **Parallel mode (`-j > 1`) was not exercised.** Everything here is `o2-sim-serial`, which is
  also why no `o2sim_Hits*.root` appears (§1).
- **Only one IRIS volume is sensitive**, because the model contains no silicon. A silicon
  sensitive layer would need `ST2513437_01`, which is in the BOM but not in `CAD_noETA.stp`.
- **CPU numbers are wall-clock on a shared 10-core box**, single-threaded and `nice`d. The 1.94×
  transport ratio was measured back to back; treat the third digit as noise.

## 12. Reproducing it

```bash
cd $HOME/alisw/O2/scripts/geometry
OUT=/some/scratch/dir

# 1. four conversions plus the JIT workaround (pythonOCC env, handled inside)
integration_demo/convert_all.sh $OUT

# 2. the full run matrix: 3 species x 2 representations + the 512-ray fan (O2 env, inside)
integration_demo/run_all.sh $OUT

# 3. reduce the MCStepLogger trees and compute the per-ray material budget
integration_demo/analyse_all.sh $OUT

# 4. the tables
python3 integration_demo/summarise_runs.py $OUT --per-volume

# geometry check and overlaps on any assembled world
root -l -b -q "integration_demo/check_geometry.C(\"$OUT/runs/matfan_exact/o2sim_geometry.root\",1)"
# external-detector hits
root -l -b -q "integration_demo/count_hits.C(\"$OUT/runs/electron_exact/o2sim.root\")"
```

`integration_demo/env_converter.sh` and `env_o2.sh` hold the two environments. They must never be
mixed: the pythonOCC `PYTHONPATH` prepends make `o2-sim` segfault at startup.

Raw run products for Track 3b are in `integration_demo/data/` (uncommitted; see the README there
for the tree and branch layout) and under
`/tmp/claude-1000/-home-swenzel-alisw-O2-scripts-geometry/b1a39244-c5c1-4627-a5e6-334a46b9188a/scratchpad/track1/runs/`.

## 13. What this invalidates

Left for the next session to reconcile, per the handoff's rule:

- `Handoff_IntegrationTest.md` §13: Bagger is now CSG 0 / surfaces 13 / tessellated 0, not
  7 / 6 / 0; `--self-test` runs 12 checks, not 48.
- `Handoff_IntegrationTest.md` §9 says `ExternalModule` is passive and produces no hits. True of
  `ExternalModule`; `ExternalDetector` exists and does produce hits, and this demo uses it.
- The claim that a converted CAD model can be injected into `o2-sim` on the exact path was never
  true before the workaround in §1.
