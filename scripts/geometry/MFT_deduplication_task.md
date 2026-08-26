# Task — MFT builds tens of thousands of duplicate logical volumes

MFT's geometry creates a separate `TGeoVolume` **and** a separate `TGeoShape` for every instance of
a component, where one volume placed many times would do. This is a task description only; the fix
belongs to another session.

Found 2026-08-25 while running the TGeo → STEP → TGeo report over the full Run 3 corpus
(`roundTripReport.py`), where MFT stood out against sixteen other modules.

## The measurement

Read off `o2sim_geometry.root` for `-m MFT`, comparing `TGeoVolume` objects, `TGeoShape` objects and
actual shape parameters:

| family | volumes | shape objects | distinct box parameters | media |
| --- | ---: | ---: | ---: | ---: |
| `capacitor` | 5144 | 5144 | **1** — `TGeoBBox(0.05, 0.025, 0.025)` | 1 |
| `welding0` | 5144 | 5144 | **1** — `TGeoBBox(0.0125, 0.025, 0.0125)` | 1 |
| `welding1` | 5144 | 5144 | **1** | 1 |
| `MFT_G` | 280 | 280 | **1** — `TGeoBBox(1.485, 0.735, 0.005)` | 1 |
| `MFT_C` | 280 | 280 | **1** — `TGeoBBox(0.75, 0.0025, 1.5)` | 1 |
| `MFTSensor` | 280 | 280 | **1** — `TGeoBBox(0.75, 0.0015, 1.5)` | 1 |

Across the module: **16 965 of 19 185 solid-carrying volumes are geometric duplicates** of another
volume in their own family. There are 2220 distinct (family, shape class, capacity) signatures. The
manager holds **25 853 volumes** where roughly 2000 would do.

Two shapes of it are worth telling apart. `capacitor` and the weldings are 5144 *separately named*
volumes (`capacitor_<i>_<j>_<k>`), so nothing in ROOT could have shared them. `MFTSensor`, `MFT_C`
and `MFT_G` are 280 volumes that all carry **one name** and are still 280 distinct objects — those
are the ones a reader would assume are already shared, and they are not.

## What it costs

- **The geometry manager**: 25 853 volumes and shapes, with the voxelisation that goes with them.
- **The STEP export**: `MFT.step` is **511 MB** with 18 871 B-rep bodies, against ~2200 distinct
  solids. Every other module is 0.3–31 MB. It is at the edge of what a CAD viewer will open.
- **The CAD round trip**: ~90 minutes for MFT alone, because each of the 5144 identical capacitors
  gets its own full acceptance test at 4000 `Contains` points.
- **Presumably `o2-sim` itself** — load time, memory, navigation. Not measured here; worth measuring
  before anyone decides how much this matters.

## The complication: sensitive volumes and hit scoring

This is the part that makes it more than a bookkeeping change, and it must be settled before any
deduplication is attempted.

MFT has exactly **one sensitive volume, `MFTSensor`** (volume id 173, medium `ALPIDE_SI$`), and it
is one of the duplicated families: 280 identical volume objects, each placed exactly once, every
placement with copy number 1. None of `capacitor`, `welding0`, `welding1`, `MFT_G` or `MFT_C` is
sensitive.

`Detector::ProcessHits` (`Detectors/ITSMFT/MFT/simulation/src/Detector.cxx`) identifies which
sensor was hit by walking **four levels of copy numbers** up the node path:

```cpp
fMC->CurrentVolOffID(++level, sensorID);
fMC->CurrentVolOffID(++level, ladderID);
fMC->CurrentVolOffID(++level, diskID);
fMC->CurrentVolOffID(++level, halfID);
Int_t sensorIndex = mGeometryTGeo->getSensorIndex(halfID, diskID, ladderID, sensorID);
```

So the sensor's identity does **not** live in the sensor volume — it lives in the copy numbers of
the placements above it, and `MFT_C` (280 duplicates) sits in that chain. Replacing 280 distinct
`MFTSensor` volumes with one volume placed 280 times is only safe if all four copy numbers still
come back the same for every sensor. Anything that changes the placement structure has to be
checked against `getSensorIndex`, and against `GeometryTGeo`'s own path/symbolic-name assumptions.

Ladder-level mothers are named `MFT_C_<half>_<disk>_<ladder>`, which suggests identity is partly
carried by *name* as well as by copy number; that needs establishing before, not after.

## What a fix has to demonstrate

1. The sensor index resolved for every hit is unchanged — same `halfID/diskID/ladderID/sensorID` for
   the same track, on a seeded `o2-sim` run compared hit for hit.
2. Alignment and `GeometryTGeo` still resolve: symbolic names, matrices, sensor count.
3. The geometry is the same solid: capacities and a containment cross-check per volume, not just a
   volume count that went down.
4. Passive volumes (`capacitor`, `welding*`) can be done independently of the sensor chain and are
   the cheap 15 432-volume win; the sensor chain is the risky part and can be left alone.

## Reproducing the measurement

```bash
o2-sim-serial -n 0 -g boxgen -m MFT          # writes o2sim_geometry.root
```

then group `gGeoManager->GetListOfVolumes()` by name with trailing `_<digits>` stripped, and per
family count distinct `TGeoShape` addresses against distinct shape parameters. The
`roundTripReport.py` "Repeated logical volumes" section does this per module from the writer report
and is the quickest way to see MFT against the other sixteen.

Corpus and report from the run that found it: `~/roundtrip_report/` (the STEP files are in
`~/roundtrip_report/step/`).
