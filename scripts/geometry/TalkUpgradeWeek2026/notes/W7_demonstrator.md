# W7 steps 4-5 — the demonstrator runs, and IRIS produces hits

2026-08-28. First end-to-end ALICE3 run: the March 2026 IRIS assembly (which carries its own beam
pipe) converted from CAD, injected into `o2-sim` as a sensitive external detector, transported, and
writing real hits. Artefacts in `../artefacts/demo/`.

## What ran

```
o2-sim-serial -n 3 -g boxgen --seed 42 \
    --detectorList EXTCAD:detectorlist.json \
    --extGeomFile externalGeometry.json \
    --configKeyValues 'BoxGun.number=20;BoxGun.pdg=211;BoxGun.prange[0]=1.0;BoxGun.prange[1]=2.0'
```

One `externalDetectors` entry: macro from the `--exact-surfaces auto --csg auto --mesh-prec 0.25
--exclude-name '^SOLID\b' --in-field` conversion, anchored to `barrel`, `detID` **TRK** (the
ALICE3 tracker slot, semantically right where the older demo borrowed `TST`), sensitive volume
`ST2513437_01` — "SILICON SENSORS IRIS TRACKER-B0-B1-B2", the part the cascade recognises exactly
as `TGeoTubeSeg` at dV_sym = 0.

**Placement was verified, not inherited.** `integration_demo/make_configs.py` uses translation
`[0, 30, 0]` and rotation `[90, 0, 0]` — `barrel` sits at cave `(0, -30, 0)`, and the CAD long axis
is CAD *y*, so `RotateX(+90)` puts it on the beam. That was measured on the *older* model, so the
new one's top bounding box was re-measured: `dx = 11.68, dy = 498.55, dz = 21.57`. Long axis still
CAD y, so the recipe carries over.

## The result

**178 hits over 3 events** — 59, 60 and 59, from 118/124/118 sensitive steps. The `IRISHit` branch
is in `o2sim.root`, as expected for `o2-sim-serial`. No stuck tracks, no navigation errors.

## `--in-field` is confirmed live, by a number the seed cannot produce

All nine IRIS media come out of the run with `ifield = 2` and **`fieldm = 15`**, none left at 0.
The seed compiled into the macro is `2,10`; 15 is what the loaded `o2::field::MagneticField`
reports for `Max()`. So the macro really did query the live field at geometry-construction time —
had the query silently failed, every medium would read 10.

That is the check worth keeping: **assert on `fieldm`, not on `ifield`**, because `ifield = 2` is
also the seed and proves nothing on its own.

## One defect found on the way, and how it was caught

The first implementation had the macro call `o2::base::Detector::initFieldTrackingParams`
directly — it is public and static, so this looked right. It works inside `o2-sim`, where
`libO2DetectorsBase` is loaded. It **segfaults a bare `root -l geom.C` session**: including
`DetectorsBase/Detector.h` drags in `FairDetector`, Cling reports "'FairDetector' is an incomplete
type" while parsing the dictionary payload, and then dies — breaking the `build_and_export()`
route the README documents and the closure test uses.

String-matching the emitted macro passed all of this happily. Only building the geometry for real
caught it. The fix is a small helper in the macro prelude over `Field/MagneticField.h` +
`TVirtualMC.h`, which Cling parses standalone; it mirrors `initFieldTrackingParams`'s body and is
commented as such. `iris_geom.root` now builds in a bare ROOT session, 36.6 MB, exit 0.

## Known, and not hidden

- **117 degenerate facets** are dropped by `O2Tessellated::AddFacet` at `--mesh-prec 0.25`
  ("Triangular facet at index N degenerated. Not adding."). Worth a look before the talk: it is a
  mesh-quality signal on the 29 mesh-tier parts, not on the exact ones.
- **18 illegal overlaps/extrusions** on `CheckOverlaps` — the known "the models are not legal
  worlds" finding (`Stream_T_AssemblyOracle.md`), and already a planned honest slide.
- `o2-sim` logs one harmless error that the JSON has no `externalModules` array. It does not: it
  has `externalDetectors`. Worth softening upstream.

## Next

Add oTOF as a second `externalDetectors` entry (`detID` **TF3**, the ALICE3 TOF slot, sensitive
volume `Module v1`) and confirm hits in both. Its own orientation still has to be measured the way
IRIS's was.

---

# Both detectors: IRIS + oTOF, hits in each

2026-08-28, same session. `../artefacts/demo/externalGeometry_two.json`, `sim_two.log`.

## Result

| | hits over 3 events | per event |
| --- | --- | --- |
| IRIS (`detID` TRK, `ST2513437_01`) | **177** | 59 / 58 / 60 |
| oTOF (`detID` TF3, `Module v1`) | **70** | 22 / 24 / 24 |

Ten external media, all carrying a field mode, none left at `ifield = 0`. Two macros JIT'd into
their own namespaces in one run, exactly as the loader is designed for.

## The placement was derived, then confirmed by the physics

oTOF's orientation was measured rather than guessed, the same way IRIS's was, and it is a
different axis: scanning every placement's world translation gives the barrel axis along CAD **x**
(0.02 to 334.90 cm), radius 92.19 in the CAD y–z plane, centred at y = 85.33. So the transform is
`rotation_deg [0, 90, 0]` — verified numerically to map local +x to master −z — with translation
`[0, -55.330, 167.460]`, which carries the measured barrel centre (167.46, 85.33, 0) onto the
barrel-frame point (0, +30, 0), i.e. the ALICE origin.

**The hits then confirm it independently.** oTOF hits come out at radius **91.98–92.50 cm** — and
the model is called `oTOF System V3-R92cm`. That radius was never an input to the transform; it
was derived from the placement scan alone, so the model's own name is an independent check that
landed. IRIS hits sit at radius **0.50–2.52 cm**, the vertex-detector barrels.

That pair of numbers is the slide: a vertex detector at 2 cm and a TOF barrel at 92 cm, both
authored in CAD, both transported, both producing hits, in one `o2-sim` run.

## Traps worth recording

- `build(false)` throws "gGeoManager is null" — the macro's `build_and_export()` creates the
  manager, `build()` does not. Create a `TGeoManager` first when calling `build()` directly.
- `build_and_export()` calls `CheckOverlaps()` unconditionally, ignoring its own `check` argument.
  On oTOF's 62 628 placements that is ~15 minutes; use `build(false)` for anything that only needs
  the volumes.
- A ROOT macro named `extent.C` collides with `std::extent` — Cling reports "reference to 'extent'
  is ambiguous" and the macro never runs. Name probe macros distinctively.
