# W7 step 3 — the ALICE3 beam pipe and IRIS, converted

Run 2026-08-28 on an idle box, immediately after the oTOF work. Artefacts under
`../artefacts/alice3/`; the 53 MB conversion output (62 facet files) stays in the session
scratchpad at `.../scratchpad/iris_conv`.

## What the models actually are

Established from the STEP `PRODUCT` histograms and the BOM, not assumed.

**`IRIS/ST2487728_01-03032026.stp` already contains the beam pipe, so no separate beam-pipe
conversion is needed.** The BOM's own title says so — "IRIS 3 SECTORS ASSY. WIITH BEAM PIPE" — and
three of its 21 line items are pipe:

| part number | name | material |
| --- | --- | --- |
| `ST1782525_01` | BEAM PIPE C SIDE-DOUBLE VACUUM-VERSION 072023 | St. Steel EN 1.4306 (304L) |
| `ST2487736_01` | EXTERNAL BP SECONDARY VACUUM-3 SECTORS-DECENTRE | Alu EN AW-5083 (H116) |
| `ST2487721_01` | CENTRAL PIPE-IRIS-2ND VACUUM-3 sectors-DECENTRE | Carbon Fiber |

The sensitive detector is `ST2513437_01`, "SILICON SENSORS IRIS TRACKER-B0-B1-B2" — three barrels,
and they come out of the cascade as exact `TGeoTubeSeg` at dV_sym = 0 (rmax 2.52 and 1.22 cm,
dz 25 cm). That matters for step 4: **the volumes we need hits from are the exactly-recognised
ones.**

The STEP holds 26 `PRODUCT` entities against the BOM's 21 rows, with **20 exact name matches**.
One BOM row (`ST1873216_01`, a DN100CF gate valve) has no counterpart in the STEP, and six STEP
products have no BOM row: the top assembly `ST2487728_01`, `ST2487730_01`, `ST1A38469_01`,
`ST1A38517_01`, and two `Gathered ST2487458_01 on CircPattern.N` labels (which do match, correctly,
by substring to the bellows row).

**`ALICE_3_example/VTX_noETA.root` is not an input — it is an old output, and it is the talk's
BEFORE picture.** It holds 73 volumes: **55 `TGeoTessellated` plus 18 assemblies, on a single
`Default` medium**. That is precisely the "every part is a tessellated blob with no material" state
slide 9 is meant to contrast against, and it is worth using as the literal BEFORE render rather
than re-manufacturing one. `ALICE_3_example/CAD_noETA.stp` (73 products) is its STEP source and an
older, more expanded export of the same assembly; `CAD.stp` (24 products) is older still. The March
2026 `IRIS/ST2487728_01-03032026.stp` is the current model and the only one with a matching BOM, so
it is the one converted here.

## The `--mesh` decision, which came first

A probe conversion **without** `--mesh` (cheap, four minutes) returned the cascade result:

```
tiers: CSG 6, exact surfaces 26, tessellated 30  (of 62 leaf solids)
```

**30 of 62 leaf solids land on the tessellated tier.** Without `--mesh` the converter does not
triangulate at all — `triangulate_asbbox()` emits a 12-triangle *bounding box* per part
(`O2_CADtoTGeo.py:4023-4025`) — so those 30 parts would have shipped as boxes. That is the same
failure the closure test caught the hard way, where one part navigated as its bounding box was
32.7x too large. NEXT.md's standing "convert ALICE3 without `--mesh`" trap is about the 22.9 GB
blow-up at the default precision, not about wanting boxes; on a model where half the parts fall to
the mesh tier, `--mesh` is mandatory and the precision is what has to be controlled.

## `--mesh-prec 0.25`, and why

Never the default 0.1. The choice is anchored on the four IRIS conversions a previous session left
on disk, measured on the dominant part `ST2487458_01` (the bellows, 2034 faces):

| `--mesh-prec` | that part's facet file | whole directory |
| --- | --- | --- |
| `TGeoOutput5_high` (≈ 0.1) | 294 MB | 420 MB |
| 0.2 | 54 MB | 84 MB |
| **0.25** | **30 MB** | **49 MB** |
| 0.6 | 4.1 MB | 11 MB |
| 1.2 | 1.3 MB | 4.4 MB |

Size goes as roughly `p^-2.35` over that range. **0.25** is chosen because it is the finest
precision a previous session already ran to completion on this model, it is 2.5x coarser than the
dangerous default, and for the 30 parts on the mesh tier the mesh *is* the shipped geometry, so
accuracy there is not a luxury — 0.6 and 1.2 trade it away for size we do not need.

The run confirms the choice exactly: the dominant part came out at **30 MB** and the whole
conversion at **49 MB of facets over 62 files**, matching the `TGeoOutput6_p0.25` precedent to the
megabyte.

## The conversion

```
O2_CADtoTGeo.py IRIS/ST2487728_01-03032026.stp -o geom.C --output-folder <dir> \
    --exact-surfaces auto --csg auto --mesh --mesh-prec 0.25 \
    --materials-csv IRIS/IRIS_MATERIALS.csv --g4-nist-json g4_nist_database/G4_NIST_DB.json
```

62 leaf solids: **CSG 6, exact surfaces 26, tessellated 30**. The 31 surface solids and 30
navigable `o2::base::O2Tessellated` volumes are emitted with `--mesh-solid o2`, the class that
actually navigates.

**Not one IRIS part has an exact tessellation.** The predicate was evaluated for 31 of the 62 (the other 31 carry no verdict in `csg_report.json`), and of those 31 the count of exact tessellations is **zero**. That is the sharpest contrast in the whole
talk and it should be on a slide: oTOF is 20/20 exact-tessellation and Run 3 is 86 %, but IRIS is
**0 %** — every part carries curved faces, so here the mesh really is an approximation and the
exact representations really are buying something. The three detectors make the point in three
different ways and none of them was arranged.

Why the 56 non-CSG parts declined, bucketed:

| reason | parts |
| --- | --- |
| a partial torus is not a whole torus | 15 |
| a sphere with extra faces is out of scope | 14 |
| free-form faces | 13 |
| too many axis clusters for the cell matcher | 12 |
| boundary gap over tolerance | 1 |
| a planar face is neither a cap nor a wedge | 1 |

**13 free-form parts, against exactly one in the whole of Run 3.** The bellows `ST2487458_01` is
924 free-form faces of 2034; the transition foil `ST2487195_01` is 118 of 182. The gap that the
Run 3 corpus makes look negligible is a real gap on ALICE3 CAD, and slide 15 should say that
plainly rather than quoting the Run 3 number alone.

## The media reality — better than the record, with one honest caveat

NEXT.md records 26 of 55 IRIS volumes matching, ten of them by accidental prefix. **That is no
longer the state**, because this is the March 2026 model with a BOM written for it:

- **61 of 62 leaf volumes carry a real medium; one falls to `med_Default`.**
- **Nine distinct media** are emitted — three aluminium alloys, two coppers, two steels, carbon
  fibre and silicon — with **zero unresolved (`FIXME`) materials**. The two steels come out as
  three-element `TGeoMixture`s.
- The one default volume is named **`SOLID`**: an unnamed body in the STEP for which no BOM row
  exists. It is not a matching failure to be fixed by tuning, it is missing source data, and it is
  not invented here.

The caveat that must travel with the 61/62: **28 of those 62 leaf volumes — 45 % — take their
material from a single BOM row**, `ST0923290_01`, the bought-in UHV gate valve. All 28 are bodies
of that one valve compound and all match its label exactly, so this is not the accidental-prefix
problem; it is the BOM treating a purchased assembly as one line item. Every one of those 28 bodies
is therefore stainless steel by assumption, including whatever inside the valve is not. For a
material-budget number that is a real, quotable uncertainty and it should be stated rather than
buried.

## The gate was not run, deliberately

`makeTestPartDB.py` builds the gate's part database by calling the converter with `--mesh` and
**no `--mesh-prec`** (`makeTestPartDB.py:131-141`), so it always meshes at the converter's default
**0.1** — the one value this model is known to be dangerous at. Extrapolating the measured size
series to p = 0.1 gives **275 MB for the dominant part alone**, which the existing
`TGeoOutput5_high` directory independently confirms at 294 MB; the whole model would be roughly
**480 MB of facets**, about ten times the run above, with the oracle's mesh column scoring
proportionally more triangles.

Two ways forward, neither taken here because both need a decision:

1. **Add a `--mesh-prec` passthrough to `makeTestPartDB.py`.** One argument, but it changes an
   acceptance instrument, and this branch has a standing habit of not doing that casually — the
   `one_way` sampling-box fix in `checkKnownSource.py` is still deliberately unmade for the same
   reason. With it, an IRIS gate at p = 0.25 should cost roughly what the numbers below say.
2. **Run it as it stands at p = 0.1** and accept ~480 MB and the longer mesh scoring.

Cost estimate, stated as an estimate: the oTOF gate proper took about six minutes on 20 parts and
1607 faces with small meshes. IRIS is 62 parts and 4177 faces, so the surface and CSG columns
should be well under an hour; the mesh column is the unknown and scales with triangle count, which
is why route 1 is worth the one-line change. A scoped-down alternative exists if the full run is
unattractive: the converter's `--include-name` can restrict the gate to the parts step 4 actually
needs — the three `ST2513437_01` sensor barrels and the three pipe parts — which is six parts of
low face count and would run in minutes.

## What blocks step 4

Nothing in the geometry. The conversion is done, the media are assigned, the sensitive volumes are
the exactly-recognised `TGeoTubeSeg` sensor barrels, and the beam pipe is already inside the same
model. Step 4 needs three decisions rather than more conversion:

1. **The gate question above** — score IRIS before or after the o2-sim run, and by which route.
2. **The `SOLID` volume on `Default`** — leave it (a zero-density placeholder material transports
   as nothing, which may be the honest thing for an unidentified body) or assign it explicitly.
3. **`media_gap.md` applies here in full.** As with oTOF, `geom.C` carries no `SetParam` call, so
   all nine media enter Geant with `ifield = 0` and a zero step-control block. On IRIS this matters
   more than on oTOF: this model is meant to sit in a magnetic field, and a silicon tracker's step
   control is not a detail. Slide 15's fix — default the eight parameters the way
   `initFieldTrackingParams()` does instead of to zero — should be applied before any material-budget
   or hit-multiplicity number from this model is quoted as physics.

---

# The gate on IRIS, run twice — and it is the first model that does not score clean

Run 2026-08-28, after `runOracleGate.py` gained a `--mesh-prec` passthrough. Artefacts:
`../artefacts/alice3/gate_prec0.25.json` and `gate_prec0.10.json`.

## The result

31 parts carry a scoreable representation (the other 30 ship as mesh). **The CSG tier is perfect,
5 of 5. The exact-surface tier passes 8 of 26.** Bagger scores 13/13 and oTOF 20/20 on the same
code, so this is not a regression in the solid — it is the first corpus hard enough to expose the
extraction.

The 18 failures split in two, and the split is the interesting part:

- **Three fail on capacity alone**, with zero disagreements on all four navigation kernels:
  `ST1A38495_01`, `ST2487455_01_b1`, `ST2487462_01`, off by 1.2e-6 to 5.8e-5 relative. A capacity
  deviation is a statement about *what geometry was extracted*, not about how it navigates.
- **Fifteen carry kernel disagreements**, always small against the sample counts. The largest
  single contributor is `ST1829909_01`'s four bodies at 22–51 `Safety` disagreements each, and
  that part is already `NEXT.md` open item 6: parity leaking through real inter-face slits at
  roughly the model's own declared tolerance. Also geometry fidelity, not algorithm.

What the gate tests is worth restating, because it is easy to read this as a verdict on
`O2BVHSurfaceSolid` alone. `--dump-brep` writes the **original OCCT solid read from the STEP
file**, scaled to cm — not a re-export of anything the converter built — so the oracle is the CAD
file itself and the gate scores the whole chain: extraction, recognition, and the shipped shape's
four kernels plus its capacity. On IRIS the evidence points at extraction on genuinely hard CAD.

## The control, and why it had to be run

The sample generator rejection-samples through the tessellated reference, so `--mesh-prec` moves
where the sample points land — and borderline `Contains`/`Safety` calls near a surface are exactly
what a coarser reference would perturb. Quoting an 8/26 obtained at 0.25 without checking that
would have been quoting our own sampling.

The whole model was therefore re-gated at the default 0.1:

| | parts | CSG | exact surfaces | verdicts changed | cost |
| --- | --- | --- | --- | --- | --- |
| `--mesh-prec 0.25` | 31 | 5/5 | **8/26** | — | 53 MB, ~15 min |
| `--mesh-prec 0.10` (default) | 31 | 5/5 | **8/26** | **0** | 472 MB, ~75 min |

**Identical verdicts, and the same three capacity-only failures.** Per-part disagreement counts
move by at most ±2 (e.g. `ST1829909_01_b1` 29 → 27 `Safety`, `ST2500409_01_b5` 2 → 4 `distout`).
The IRIS result is real, and it is not an artefact of the precision choice.

That also settles the gate route Sandro left open: **0.25 is the right setting for this model** —
same answer, a tenth of the disk, a fifth of the time. The passthrough defaults to the converter's
own 0.1, so every gate result on record is reproduced unchanged; Bagger re-run after the change is
still 13/13 surface and 7/7 shape, and the gate self-test is still 17/17.

## What this means for the talk

It is a better slide than a clean sweep would have been. Three detectors, three different answers
from one instrument: Bagger exact end to end, oTOF 20/20 with an exact tessellation for every
part, and IRIS — real, current, industrially-authored ALICE3 CAD — passing 8 of 26 on the exact
surface tier with 13 free-form parts. The honest claim is that the pipeline is exact where the CAD
is clean and measurably not where it is not, and that the instrument says which is which per part.
