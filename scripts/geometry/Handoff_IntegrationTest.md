# Handoff — the Geant integration test for converted CAD geometry

**Written 2026-08-04 for a fresh session to execute autonomously.** You will have a long stretch of
unattended time. Everything you need to start is in this file; read it fully before touching
anything. Where it says *measured*, the number was taken on this machine and can be trusted; where
it says *unverified*, check it yourself before you build on it.

Branch `swenzel/bvhsurfacesolid`. Read [`Tutorial.md`](Tutorial.md) for what the conversion pipeline
is, [`NEXT.md`](NEXT.md) for the current state and the environment traps, and
[`Roadmap.md`](Roadmap.md) for what is deliberately *not* in scope.

---

## 1. What this is for

Every correctness result on this branch so far comes from instruments we built ourselves — the
oracle gate, the X-ray transport benchmark, the representation benchmark. They are good, but they
all ask geometric questions in isolation. **This task puts the converted geometry under a real Geant
transport and asks whether particles get through it.** That is a different question and it is the
one that decides whether any of this is usable.

It is also the harness that unlocks two roadmap items which cannot be measured any other way: the
step-length cost of a loose safety, and the real per-event cost of the exact path against a mesh.

Sandro's framing, close to verbatim:

> *"It needs to be based on realistic examples, say ALICE3 (inner vertex detector) + oTOF. The test
> harness is simply a conversion to separate geom.C files (including a material assignment: for
> ALICE3 present, for TOF use something which doesn't kill particles immediately but also not just
> air and which also produces few secondaries). The goal is verify that simulating few particles
> goes through via the externalModule mechanism. Note that you need to rotate and shift the detector
> so that it aligns with the 0 of ALICE. Compare pure tessellation mode with the case where
> CSG/BVHSurface solid is used. Study step numbers etc. Make it systematic and comprehensive.
> Probably you need to place it inside the centralBarrel virtual container (pitfall: might not
> fit?? -> choose something else)."*

## 2. Definition of done

1. Both detectors convert to their **own** `geom.C`, each with materials assigned.
2. Both are injected through `ExternalModule` and `o2-sim` transports particles through them
   **without a stuck track, a navigation error, or an abort**.
3. The same is done twice — once **pure tessellated**, once **exact (CSG + BVHSurfaceSolid)** — and
   the two are compared systematically on step counts, secondaries, CPU, memory and material budget.
4. Everything is reproducible from a script in the repo, not from your shell history.
5. A stream document records the numbers, with the failures stated as plainly as the successes.

**Partial completion is fine and expected.** If oTOF cannot be made to convert (see §5), do the
whole thing with ALICE3 alone and say so. A complete, honest single-detector result beats a
half-finished two-detector one.

---

## 3. The mechanism, concretely

`o2::passive::ExternalModule` loads a ROOT macro and hooks its top volume into an existing volume.
Two ways to drive it; **prefer the JSON**, it needs no recompilation.

**JSON** — `ExternalModule::createFromJSON`, wired to `o2-sim --extGeomFile <file>`. Schema, read
from `Detectors/Passive/src/ExternalModule.cxx` (verified against the source, not the docs):

```json
{
  "externalModules": [
    {
      "name":  "A3VTX",
      "title": "ALICE 3 IRIS vertex detector",
      "macro": "/abs/path/to/alice3/geom.C",
      "anchor": "barrel",
      "placement": { "translation": [0, 30, 0], "rotation_deg": [90, 0, 0] }
    }
  ]
}
```
`name`, `macro` and `anchor` are required; an entry missing any of them is **skipped with an error
and the simulation continues** — so check the log for `Configured external module` rather than
assuming it loaded. A module is only added when its `name` is in the active module list, so it must
also appear in the detector-list JSON (§7).

**C++** — the `build_geometry.C` pattern in
[`simulating_CAD_modules.md`](simulating_CAD_modules.md) §4. Use only if the JSON path proves
insufficient; note that document predates most of this branch and its "Known Limitations" section is
substantially out of date (it says the tool converts *only* to tessellated solids, which has not
been true for some time).

## 4. Coordinates — the shift is real and here is where it comes from

`Detectors/Passive/src/Cave.cxx` builds the container and places it **off-centre**:

```cpp
TGeoVolume* voBarrel = new TGeoVolume("barrel", shCaveTR, kMedAir);
cavevol->AddNode(voBarrel, 1, new TGeoTranslation(0., -30., 0.));
```

So `barrel` sits 30 cm **below** the ALICE origin, and everything already inside it is placed at
`y = +30` to compensate (see `Pipe.cxx`, `Absorber.cxx`). **A module anchored to `barrel` must
therefore be placed at `y = +30` to land on the ALICE origin.** That is the `SetDy(30)` in the
worked example, and it is the shift Sandro means.

The rotation in that example is `RotateX(90)`, i.e. the CAD is authored Z-up and ALICE wants Z along
the beam. **Do not take that on faith for these two models** — verify it from each model's own
bounding box (a vertex detector is long in the beam direction and narrow across it, so the right
axis is obvious from the extents) and state what you found.

**Does it fit?** `barrel` is a `TGeoTube`-based composite of radius **790.5 cm** spanning roughly
`z = -714.6 … +714.6`, minus two end pockets. ALICE3's IRIS and an oTOF at R = 92 cm are both far
inside that, so `barrel` should be the right anchor and the "might not fit" worry probably does not
bite. **Verify anyway** by comparing each converted model's bounding box against the barrel shape
before running, because a module that does not fit produces an overlap rather than an error. If one
genuinely does not fit, `cave` is the fallback anchor — and then the y = +30 shift does **not**
apply, since `cave` is the frame `barrel` is displaced within.

---

## 5. Stage 0 — the blocker, and it is oTOF

**`oTOF System V3-R92cm.step` does not currently convert, and this is the single biggest risk to the
task.** Both halves of the problem are measured:

- **The converter sees 3 leaf solids.** A plain `STEPControl_Reader` sees **62 628 placed solids,
  all distinct** (multiplicity histogram `{1: 62628}`, measured 2026-08-04 — no coincident
  duplication of the kind ALICE3 had). So the STEP file is fine and the defect is in the converter's
  **XCAF traversal**, which is a different code path from the plain reader.
- **Without `--mesh` it does not even complete** — `triangulate_asbbox()` dies with
  `Standard_ConstructionError: Bnd_Box is void`.

62 628 placements over ~20 prototypes is exactly the shape TGeo handles well — one volume, many
nodes — and exactly the case the placed-primitive work (`Stream_N_PlacedPrimitives.md`) made cheap.
So this is worth fixing properly rather than working around.

**Do this first, timeboxed.** If after a genuine attempt oTOF still will not convert, **stop, record
precisely where it fails, and proceed with ALICE3 alone.** Do not fabricate a substitute geometry to
have two detectors; a one-detector result with a named blocker is the more useful outcome.

A note on the `Bnd_Box is void` failure: it comes from the no-mesh path building a bounding box for a
solid it could not read, so it is probably a *symptom* of the traversal problem rather than an
independent bug. Fix the traversal first and re-check.

## 6. Stage 1 — convert, twice each

Four conversions in total. Keep them in separate output folders and **never** write into
`STEP_examples/` or `ALICE_3_example/`.

| | model | mode | flags |
| --- | --- | --- | --- |
| A | ALICE3 IRIS | exact | `--exact-surfaces auto --csg auto` |
| B | ALICE3 IRIS | tessellated | `--mesh --mesh-prec <p>` |
| C | oTOF | exact | `--exact-surfaces auto --csg auto` |
| D | oTOF | tessellated | `--mesh --mesh-prec <p>` |

**ALICE3 IRIS is `scripts/geometry/ALICE_3_example/CAD_noETA.stp`** — note that IRIS *is* the inner
vertex detector Sandro means: its root child is `ST2487728_01`, "IRIS 3 SECTORS ASSY. WITH BEAM
PIPE". Current state, measured: **103 leaf placements over 55 definitions**, of which **20 emit an
exact sidecar**; the rest fall back to tessellation, which is the standing and accepted bargain.

**On `--mesh-prec`:** it currently sets linear *and* angular deflection to the same value, and
`Stream_P_RepresentationBench.md` §5 measured that it behaves as an **angular** knob — a 20× change
in linear deflection at fixed angular moved nothing. Do **not** use the default 0.1 on ALICE3: one
2 m sphere reached **22.9 GB** resident and was killed. Coarser values are affordable — the full
model meshes at 679–841 MB at precision 1.0/0.5. Pick a value, justify it, and record it; the
precision *sweep* is a roadmap item, not this task.

**Prior art worth reading before you start:** `scripts/geometry/IRIS/` already contains a worked
run of much of this — `IRIS_MATERIALS.csv`, `detList.json`, several `TGeoOutput_p*` folders at
different mesh precisions, a `stepsG3/SimpleStepAnalysis` tool, and `baz_*` simulation outputs.
Reuse it rather than reinventing; but re-derive any number you intend to quote, because those files
predate everything on this branch.

## 7. Stage 2 — materials

**ALICE3: materials exist.** Use `scripts/geometry/IRIS/IRIS_MATERIALS.csv` with
`--materials-csv`, plus `--g4-nist-json scripts/geometry/g4_nist_database/G4_NIST_DB.json`. That CSV
is a Bill-of-Materials export with header rows above the real table — check the converter parses it
as intended rather than assuming, and report how many parts got a material and how many fell back.

**oTOF: choose one, and justify it.** Sandro's constraints are: not air, does not stop particles
immediately, few secondaries. That means **low Z** (few secondaries: less bremsstrahlung and pair
production, long radiation length) at **moderate density** (not air, but not a stopper).

Recommended default: **`G4_POLYSTYRENE`** (ρ ≈ 1.06 g/cm³, ⟨Z⟩ ≈ 5.6, X₀ ≈ 42 cm). It is physically
plausible for a TOF (plastic scintillator), it is unambiguously not air, and it is about as quiet as
anything of that density gets. `G4_PLEXIGLASS` (ρ 1.19) is an equally defensible alternative.

Consider **`G4_Si`** as a *second* variant if time allows — it is what an LGAD-based TOF would really
be, and at X₀ ≈ 9.4 cm it produces visibly more secondaries, which makes it a useful contrast rather
than a competitor. Do not make it the default.

Assign one material to the whole oTOF module; per-part materials there are not the point of this
test.

## 8. Stage 3 — place and check before simulating

Build the geometry **without transporting anything first**, and check:

1. Both modules appear in the log (`Configured external module ...`).
2. Bounding boxes land where intended — the module's world-frame extent should straddle the ALICE
   origin as the CAD says, not sit 30 cm off.
3. `TGeoManager::CheckOverlaps` on the assembled world. **Caveat, measured:** that checker is
   currently unreliable on `O2BVHSurfaceSolid` — starved by default because `GetPointsOnSegments`
   returns without filling its buffer, and wrong when fed (it reported a provably 0.09 cm gap as a
   0.41 cm overlap, and one of its Bagger entries was a 24-gon's chord error). A partial fix is
   committed as **unverified WIP** in `00cc7d9eb3`; treat its output as a hint, not a verdict, and
   say so in your write-up. `scripts/geometry/overlapCensus.py` and `assemblyOracle.py` are the
   trustworthy OCCT-side instruments if you need a real answer.
4. Overlaps *between* the two modules and with existing ALICE volumes. Sandro has explicitly
   **deprioritised overlaps in general** (TGeoNavigator tolerates them, the ALICE geometry has had
   them for years, and the target is the navigator plugin, not native Geant4) — so **report them and
   move on**. Do not spend the night fixing overlaps.

## 9. Stage 4 — simulate

Few particles, and make it deterministic: fix the seed, fix the generator, keep event counts small
enough that a full run is minutes not hours. Start with **geantinos** (they transport without
interacting, so any failure is purely geometric), then move to real particles.

Note that `ExternalModule` is a **passive** module — no sensitive volumes, so **no hits**. The
observables are steps, secondaries, timing and material budget. That is exactly what Sandro asked
for; do not go looking for hit trees.

The instrument for step counting is **MCStepLogger** (present in the stack at
`$ALIBUILD_ARCH_ROOT/MCStepLogger`). It produces per-volume step tallies, which is what makes the
tessellated-vs-exact comparison meaningful rather than a single aggregate. `IRIS/stepsG3/` shows a
previous analysis of this kind.

## 10. Stage 5 — the comparison, which is the actual deliverable

Same particles, same seed, same everything, **only the representation differing**. Per detector,
per representation, record at minimum:

- **steps per event**, total and **per volume** (an aggregate says *that*, never *where*)
- **secondaries produced**
- **CPU time per event**, and geometry **build/close** time separately from transport
- **resident memory**, and the TGeo volume/node counts
- **transport robustness**: stuck tracks, navigation errors, aborted events, tracks leaving through a
  wall they should not have
- **material budget along fixed directions**, via a geantino scan

**The material-budget scan is the strongest equivalence test available here** and I would treat it
as the headline. Two representations of the same CAD solid must present the same X₀ integral along
the same ray. It is the physics-level analogue of the X-ray benchmark that has driven this whole
branch, and unlike step counts it has a *right answer* the two representations must agree on.

Two cautions carried over from measurements on this branch:

- **A parallel or axis-aligned beam is direction-poor** — three axis beams are three directions
  however many rays you fire. Use a Fibonacci fan of directions for any scan. This trap has bitten
  this project twice, including once inside an instrument built to detect it.
- **A mesh can be invalid, not merely inaccurate.** `meshClosedBody = false` on even one part means
  tracks walk through walls, and chordal accuracy does not detect it — see
  [`MeshHealing.md`](MeshHealing.md). So if the tessellated run shows fewer steps, check whether it
  is faster or simply **leaking** before reporting it as a win.

Expect the exact path to be **slower per query** — the representation benchmark measured
`O2BVHSurfaceSolid` at 40–145× a `TGeoCompositeShape` doing identical work, with 61% of `Contains`
in `Curve2D::closestPoint`. Whether that survives contact with a real transport, where geometry is
one cost among many, is precisely what this test exists to find out. **Do not optimise anything
while measuring**; a speedup found mid-run is a finding to record, not a change to make.

---

## 11. Environment, and the traps that have cost this project time

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
SW=$HOME/alisw/sw/ubuntu2404_aarch64
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid \
      o2-bench-detectorsbase-solid-harness o2-bench-detectorsbase-xray
ctest -R BVHSurfaceSolid      # 108 cases, green
```

For the converter (needs pythonOCC, and does **not** self-re-exec):
```bash
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
O2_BUILD_DIR=$B $SW/Python/latest/bin/python3.10 scripts/geometry/O2_CADtoTGeo.py --self-test
```

- The `eval` and the `export` must be in the **same shell command** as whatever runs next, or
  anything linking O2 dies on `libarrow_acero.so`. This applies even to `--help`.
- **Prepend, never replace** `PYTHONPATH`/`LD_LIBRARY_PATH`. Replacing kills PyROOT on
  `libffi.so.6`. `scripts/geometry/csg/occ_env.py` encodes the correct form.
- **`which o2-sim` resolves to the *installed* O2**, not this branch's build. Verified: it returns
  `$SW/O2/dev-local2/bin/o2-sim`. Put `$B/stage/lib:$B/stage/lib64` first on `LD_LIBRARY_PATH` and
  use `$B/stage/bin/o2-sim`, or you will be testing code you did not change.
- Default Bash tool timeout is 120 s and will kill real runs — pass `timeout: 600000`, and **detach**
  anything longer, writing to a uniquely-named `--out` path rather than trusting stdout. The scratch
  directory is shared between concurrent sessions; an hour-long census was lost exactly that way.
- Only **one `ninja` at a time** — the build directory is shared.
- `manifest.json` stores **absolute** paths; a copied workdir silently scores the original.
- The oracle gate exits non-zero whenever any part fails. Normal — read the counts. Bagger is the
  exception and now exits 0.
- Machine hygiene: shared interactive box, 10 cores. Keep timing runs single-threaded and `nice`
  anything long, or you will corrupt your own measurements.

## 12. How to work

- **Test-driven, stage by stage, committing after each**, in this branch's style: the commit subject
  states what is now true, not what you did. Match `git log --oneline -10`.
- **Verify, do not relay.** Every number you report should be one you ran. This branch has a history
  of claims that did not survive an independent check, including several of mine.
- **A refuted experiment is not a refuted hypothesis.** If a comparison shows no difference, prove
  the instrument *could* have shown one before concluding anything.
- **Never `git stash`. Never `git add -A` / `git add .`** — the tree has many untracked files that
  are not yours; stage explicit paths and commit in one step.
- Write `scripts/geometry/Stream_X_IntegrationTest.md`. **Do not edit `NEXT.md`, `Tutorial.md` or
  `Roadmap.md`** — list what your work invalidates and leave the reconciliation to the next session.
- Leave the repo green: `ctest -R BVHSurfaceSolid` at 108 cases, the fixtures gate 10/10 with
  `surface` 0/0/0/0, the Bagger gate 13/13 exiting 0. If you break one of those, fix it or revert
  before you stop.

## 13. Current baseline, so you can tell what you changed

| | |
| --- | --- |
| `ctest -R BVHSurfaceSolid` | 108 cases, green |
| `O2_CADtoTGeo.py --self-test` | 48 checks, 0 failures |
| `runOracleGate.py --self-test` / `csg/emit.py --self-test` | 17/17 · 33/33 |
| `o2-bench-detectorsbase-xray --self-test` | 34 checks, 0 failures |
| ladder fixtures | 10/10 scored and passing, `surface` 0/0/0/0 |
| `Bagger.step` | 13/13 passing, CSG 7 / surfaces 6 / **tessellated 0**, gate exits 0 |
| ALICE3 `CAD_noETA.stp` | 103 leaf placements / 55 definitions, 20 exact sidecars |
| oTOF | **does not convert** — see §5 |
