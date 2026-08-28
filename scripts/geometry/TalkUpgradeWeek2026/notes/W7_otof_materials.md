# W7 step 1 — oTOF materials, deduced from AliceO2 source

Done 2026-08-28. Product: `../oTOF_MATERIALS.csv` (in `scripts/geometry/`, beside
`IRIS/IRIS_MATERIALS.csv`).

## What the AliceO2 model actually says

`Detectors/Upgrades/ALICE3/IOTOF/simulation/src/Detector.cxx::createMaterials()` defines exactly
**two** media for the whole of IOTOF:

- `AIR$` — a four-component mixture (C, N, O, Ar at 0.000124 / 0.755267 / 0.231781 / 0.012827,
  rho = 1.20479e-3 g/cm3);
- `SILICON$` — a plain material, A = 28.086, Z = 14, rho = 2.33 g/cm3, X0 = 9.36 cm,
  lambda = 999 cm.

`Layer.cxx` then shows where each one goes: the **sensor** volume and the **chip** volume are both
built with `medSi`, and `medAir` is used only for the **layer envelope** that contains them
(`Layer.cxx:143-147`, and the same pattern at :281 and :457).

## The mapping to the CAD model, and why it is defensible

The CAD file `STEP_examples/oTOF System V3-R92cm.step` contains 208 distinct `PRODUCT` names in
six families — `Plate 1` (67), `Plate 2` (67), `Component1` (67), plus `Module v1`, `Component`,
`Plate` — and three pure assemblies (`oTOF System V3-R92cm`, `A Side`, `oTOF v2`).

**Every solid body gets silicon.** That is not a simplification imposed for the demo: in the
AliceO2 model every volume that has a *body* is silicon, and air appears only as an envelope, which
in the CAD world has no counterpart — the mother is the world volume. So the CSV assigns `Silicon`
to all 205 body-bearing products and does not name air at all.

## Verified, not assumed

Run against the converter's own code (converter env, `read_bom_csv` and `resolve_bom_material`):

- the CSV parses to **205 BOM entries**, one per body-bearing product, one material (`Silicon`);
- `Silicon` resolves to Geant4 NIST **`G4_Si`**: rho **2.33 g/cm3**, X0 **9.3661 cm**, one element
  Si at Z = 14, A = 28.0854 — token score 0.667, combined 0.500, comfortably over the 0.35
  acceptance threshold, and unambiguous;
- those values reproduce AliceO2's hand-written `SILICON$` (2.33, 9.36, Z = 14, A = 28.086) to
  three or four digits. The NIST route is therefore not a substitute for the source values, it
  **is** them.

Mass is left blank on purpose: with no mass the converter takes the NIST density rather than
deriving an effective one from mass over volume, which is what we want for a solid silicon plate.

## What this does not yet give

The eight Geant medium parameters and the physics cuts — see `media_gap.md`. AliceO2 sets
`ifield = 2, fieldm = 10` via `initFieldTrackingParams()`, `stemax = 0.0075 cm` and
`tmaxfd = 0.1 deg` for silicon; the CAD route emits a three-argument `TGeoMedium` and all eight
land at zero. For the demo that is transport with default step control, and slide 15 says so.

## Next (W7 step 2)

Convert the STEP with `--materials-csv oTOF_MATERIALS.csv --g4-nist-json
g4_nist_database/G4_NIST_DB.json`, then score it through `runOracleGate.py` — oTOF converts today
(20 prototypes, 62 628 placements, 20/20 exact surfaces, 19/20 CSG `TGeoBBox` at dV_sym = 0) but
has never been scored. Convert **without** `--mesh`, and strictly serially.

---

# W7 step 2 — oTOF scored, and converted with its materials

Run 2026-08-28 on an idle box. Artefacts under `../artefacts/otof/`; the table is
`../tables/T4_otof_gate.md`.

## The gate, run on oTOF for the first time

```
python3 scripts/geometry/runOracleGate.py --workdir <wd> \
    --model scripts/geometry/STEP_examples/oTOF\ System\ V3-R92cm.step
```

**20 of 20 parts pass on their shipped representation, no failures.** 19 ship as CSG — every one a
`TGeoBBox` found by `tier1-box` at **dV_sym = 0** — and one as an exact `O2BVHSurfaceSolid`. Summed
over parts, against the OCCT oracle:

| representation | contains | distFromOutside | distFromInside | safety | parts with zero disagreements |
| --- | --- | --- | --- | --- | --- |
| surface | 0 | 0 | 0 | 0 | **20/20** |
| shape (CSG) | 0 | 0 | 0 | 0 | **19/19** |
| mesh | 0 | 1950 | 5 | 86 | 0/20 |

This confirms `Stream_AC_OTOFTraversal.md`'s claim and upgrades it from *converted* to *scored*.

**Every one of the 20 has an exact tessellation** — 100 %, every face a planar polygon — so for
oTOF the mesh is not an approximation of the part, it is the part. That is a better line for the
talk than any speed ratio: on this detector the three representations are not a trade-off between
accuracy and cost, they are three exact descriptions of the same solid.

The single part that declines CSG, `oTOF_v2_v1_0_1_1_3_b3`, does so for a stated reason worth
quoting: 1493 planar faces, no six-plane box, no prism axis, 1993 trusted concave edges of 4330 —
and the cell decomposition stops at the budget of 64 after producing 62 cells. It is not a
free-form-surface failure, it is a budget.

## The conversion with materials

```
O2_CADtoTGeo.py "oTOF System V3-R92cm.step" -o geom.C --output-folder <dir> \
    --exact-surfaces auto --csg auto \
    --materials-csv oTOF_MATERIALS.csv --g4-nist-json g4_nist_database/G4_NIST_DB.json
```

206 BOM entries, **226 CAD logical volumes matched**, tiers CSG 19 / exact surfaces 1 /
tessellated 0. The emitted macro carries

```cpp
TGeoMaterial *mat_Silicon = new TGeoMaterial("Silicon", 28.08536146, 14, 2.33);
mat_Silicon->SetRadLen(9.366070292, 45.66030737);
TGeoMedium   *med_Silicon = new TGeoMedium("Silicon", 2, mat_Silicon);
```

and **all 20 leaf volumes sit on it, none on `Default`**.

### One mistake worth recording

The first CSV excluded `oTOF v2` along with the file root and `A Side`, on the assumption it was an
assembly. It is not: XCAF label `0:1:1:3` is a *compound* holding three solid bodies — exactly the
leaf-granularity distinction Stream AC was written about. The three bodies silently fell back to
`med_Default`, and the only visible symptom was `med_Default` appearing three times in `geom.C`.
Corrected: the exclusion list is now the two labels that carry **no solid body at all** (the file
root, and `A Side` at 0 bodies per Stream AC's table).

The lesson generalises to the check the converter should make: a leaf volume that ends on
`Default` is almost always an unmatched BOM row, and it is currently silent.

### And it confirms `media_gap.md` on a real product

`geom.C` contains **no `SetParam` call anywhere**. Every one of oTOF's media therefore enters
Geant with `ifield = 0`, `fieldm = 0` and a zero step-control block, where AliceO2's own IOTOF sets
`ifield = 2`, `fieldm = 10`, `stemax = 0.0075 cm`, `tmaxfd = 0.1 deg` for silicon. Slide 15's claim
is now demonstrated, not read off the source.

## Gate 1 of W7 is passed

Next: convert the ALICE3 beam pipe and IRIS, compose all three as `externalDetectors`, and run
`o2-sim-serial`.
