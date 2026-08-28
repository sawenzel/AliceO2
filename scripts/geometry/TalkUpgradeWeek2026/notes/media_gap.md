# W9 — media, magnetic field and physics cuts on a CAD-authored part

Read 2026-08-28 from the source, for slide 15. The question is what a part that exists **only** as
CAD — no TGeo original to copy from — actually gets for its magnetic-field tracking parameters and
its physics cuts, and the answer is: nothing, silently.

## The three routes, and which one the closure test exercised

The converter emits media by one of two routes, and `remapCADMedia` then re-registers whatever it
finds with `MaterialManager`.

1. **`emit_media_sidecar_cpp`** — used when the geometry came out of TGeo in the first place (the
   round trip, the closure test). The media are known exactly, so they are rebuilt field for field:
   the material recipe (Z, A, density, or every element of a mixture with its weight) **and all
   eight Geant medium parameters**, `isvol ifield fieldm tmaxfd stemax deemax epsil stmin`, written
   out as explicit `SetParam` calls. `ifield` is what distinguishes a field-free `_NF` twin from
   its in-field original, and this route carries it.
2. **`emit_materials_cpp`** — the genuine-CAD route: a BOM material name matched against a Geant4
   NIST database, effective densities derived where the CSV gives a fill fraction. It emits
   `new TGeoMedium(name, id, mat)` and **no `SetParam` at all**. ROOT's three-argument
   `TGeoMedium` constructor zeroes all twenty parameter slots
   (`TGeoMedium.cxx:50`), so every one of the eight is **0**.
3. `remapCADMedia` (`Detectors/Base/src/CADGeometryUtils.cxx:106`) walks the volumes, registers each
   material once — a `TGeoMixture` through `Mixture()` so the element composition survives, a plain
   material through `Material()` — reads the eight parameters back off the `TGeoMedium` and hands
   them to `MaterialManager::Medium`.

**So the closure test's "95 of 95 media agree on cuts and processes" is a statement about route 1
and says nothing about route 2.** That distinction has to be on the slide, or the result reads as
a guarantee it is not.

## What is actually missing for a CAD-only part

- **Magnetic field tracking**: `ifield = 0`, `fieldm = 0`. An ordinary O2 detector calls
  `o2::base::Detector::initFieldTrackingParams()` in its `createMaterials()`, which asks the live
  `o2::field::MagneticField` for its integration method and maximum field and falls back to
  `integration = 2, maxfield = 10` (`Detectors/Base/src/Detector.cxx:144`). Neither
  `ExternalDetector` nor `ExternalModule` calls it — grep is empty in both.
- **Step control**: `tmaxfd stemax deemax epsil stmin` are all 0 rather than a detector's chosen
  values.
- **Special physics cuts and process flags**: `Detector::SetSpecialPhysicsCuts()` reads
  `$O2_ROOT/share/Detectors/<NAME>/simulation/data/simcuts.dat`
  (`Detectors/Base/src/Detector.cxx:125`). For a CAD module called `IRIS` or `oTOF` that file does
  not exist, and nothing in `ExternalDetector` / `ExternalModule` calls the method anyway, so the
  media fall through to the global defaults. This is the same failure mode already recorded in
  another context: a tracking medium with no `simcuts` line collapses to the 1 MeV default.

None of this is loud. There is no warning, no refusal — the simulation runs and the numbers look
plausible.

## The three bullets for slide 15

- A CAD model carries a **material** (via its BOM, matched to the Geant4 NIST database) but it
  cannot carry a **medium**: field integration, step limits and production cuts are simulation
  choices, not properties of a part, and no CAD format has a place to put them.
- Today they silently default to zero for a CAD-authored part, because the emitter writes a
  three-argument `TGeoMedium` and the external-module wrappers never call the two methods
  (`initFieldTrackingParams`, `SetSpecialPhysicsCuts`) that every hand-written O2 detector calls.
- The fix is small and worth naming as the next step: let the media CSV carry the eight parameters
  and a cuts block per medium, default them the way `initFieldTrackingParams` does rather than to
  zero, and make a missing entry **loud**.

## Follow-ups this opens (not for the talk, for the roadmap)

- Extend `IRIS_MATERIALS.csv` / the media CSV schema with the eight medium parameters and an
  optional per-medium cuts block; have `emit_materials_cpp` write them.
- Have `ExternalDetector` / `ExternalModule` call `initFieldTrackingParams()` and offer a
  `simcuts` path through the external-geometry JSON.
- A check in `remapCADMedia` that logs at warning level when a medium arrives with `ifield == 0`
  and a zero step-control block, so the default is visible in the log.

---

## Resolved 2026-08-28: `--in-field`, and it asks the live field

Sandro's call: a switch that says some modules are transported in the field, and it must be
**live**, not an asserted pair. Implemented in `O2_CADtoTGeo.py`.

With `--in-field`, the emitted macro carries

```cpp
int   cad_ifield = 2;
float cad_fieldm = 10;
o2::base::Detector::initFieldTrackingParams(cad_ifield, cad_fieldm);
```

and every medium then takes `SetParam(1, cad_ifield)` and `SetParam(2, cad_fieldm)`.
`initFieldTrackingParams` is a **public static** member of `o2::base::Detector`
(`Detectors/Base/include/DetectorsBase/Detector.h:236`) which queries the loaded
`o2::field::MagneticField` for its `Integral()` and `Max()`, and the macro is JIT'd inside
`o2-sim` during geometry construction — the same moment every hand-written detector calls it from
its own `createMaterials()`. So a CAD module is no longer treated differently from a coded one, and
nothing is baked in: `2,10` is only the seed those variables carry into the call, which is exactly
what the function itself falls back to when no field is loaded. `--in-field 1,5.5` overrides the
seed.

Two things this deliberately does **not** do:

- It does not touch the media-sidecar route, which already carries the source geometry's own eight
  parameters verbatim — and where `ifield = 0` is sometimes a real assertion, because ALICE builds
  field-free `_NF` twins of media. A blanket "substitute the live field wherever ifield is 0" in
  `remapCADMedia` would wrongly put those twins in the field; the converter is the right place to
  declare intent because only there is the route known.
- It does not invent step control. `tmaxfd`, `stemax`, `deemax`, `epsil` and `stmin` stay 0, the
  transport default. Those are per-detector physics choices and a CAD file carries no such
  statement.

Sixteen self-test checks cover it, including two negative controls: without the flag **no
`SetParam` is written at all** (so no existing converted module changes its transport), and
`DetectorsBase/Detector.h` is not pulled into the macro prelude.

Terminology, corrected on the record: `ifield = 2` is *live per-step evaluation* of a non-uniform
field, not a constant one — `ifield = 1` is the uniform case — and `fieldm` is only the maximum
field used for step estimation. The defect this flag fixes was never "the field is constant", it
was `ifield = 0`, meaning no field tracking at all.
