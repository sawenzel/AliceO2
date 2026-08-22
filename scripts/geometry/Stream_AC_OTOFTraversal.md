# Stream AC — the oTOF XCAF traversal

`oTOF System V3-R92cm.step` now converts. The converter used to see **3 leaf solids** in a file
that holds **62 628 placed solids over 20 prototypes**, and without `--mesh` it died in
`triangulate_asbbox()` with `Standard_ConstructionError: Bnd_Box is void`. Both were one defect,
and the defect was **leaf granularity**, not placement expansion.

Measured 2026-08-22 on `swenzel/bvhsurfacesolid`. Everything below is a number taken on this
machine.

---

## 1. The mechanism

An XCAF *simple shape* is a **label**, not a **body**. `XCAFDoc_ShapeTool::IsSimpleShape(label)`
says only "this label is not an assembly of components"; the `TopoDS_Shape` hanging off it may
still be a `TopAbs_COMPOUND` holding any number of solid bodies. CAD exporters do this routinely,
and `oTOF System V3-R92cm.step` does it for every one of its leaves:

| leaf label | XCAF name | shape type | solid bodies | shells | faces |
| --- | --- | --- | --- | --- | --- |
| `0:1:1:3` | `oTOF v2 v1` | `TopAbs_COMPOUND` | **3** | 3 | 1505 |
| `0:1:1:4` | `Module v1` | `TopAbs_COMPOUND` | **17** | 17 | 102 |
| `0:1:1:7` | `A Side` | `TopAbs_COMPOUND` | **0** | 0 | 0 |

`expand_definition`'s leaf branch registered **one logical volume per label**, whatever the shape
held. So the whole `Module v1` compound — 17 separate bodies — became a single
`TGeoTessellated`/bounding box, and 16 of its 17 bodies vanished into the union's bbox. The empty
`A Side` compound is where `Bnd_Box is void` came from: a bounding box of nothing.

For comparison, every leaf of the other three corpora is a bare `TopAbs_SOLID` with exactly one
body — Bagger 13/13, `as1-oc-214.stp` 5/5, ALICE3 55/55. That is why the defect only ever showed
on oTOF.

## 2. Placement expansion was already correct — the reconciliation

Sandro's reading of the symptom was that this is the difference between logical and placed
volumes, and that there really are only three logical volumes. Half of that is right: XCAF really
does declare only three *leaf labels*. But those three labels are not three volumes, and the
placement side was never broken. Driving the production traversal (`expand_free_shapes`) over
oTOF, before any change, with meshing stubbed out so the empty leaf could not abort it:

```
logical_volumes (leaf DEFINITIONS): 3        0:1:1:3, 0:1:1:4, 0:1:1:7
assemblies:                       205
placement edges:                 3945
leaf OCCURRENCES:                3741        0:1:1:3 x 68, 0:1:1:4 x 3672, 0:1:1:7 x 1
Placement check: 3741 leaf placement(s), all at distinct world transforms.
```

The four counting bases that were floating around, reconciled:

| basis | count | what it counts | who reports it |
| --- | --- | --- | --- |
| XCAF free shapes | **1** | `0:1:1:1` `oTOF System V3-R92cm v2`, one root compound | `GetFreeShapes` |
| leaf *labels* | **3** | labels passing `IsSimpleShape` | the converter, before the fix |
| solid *bodies* in those labels | **20** | 3 + 17 + 0 | `csg/census.py`, as "prototypes" |
| label placements | **3 741** | 68 + 3672 + 1 | the converter's own placement expansion |
| body placements | **62 628** | 3x68 + 17x3672 + 0x1 | plain `STEPControl_Reader`, `Stream_W` |

`3x68 + 17x3672 = 204 + 62 424 = 62 628`, exactly the plain reader's figure and exactly the
census's. Nothing was lost in the placement expansion; the converter was placing three *bags* of
bodies 3 741 times instead of twenty *bodies* 62 628 times.

Why `csg/census.py` got it right is now also plain, and it is not a better walk: it never uses the
label walk for geometry at all. `load_step_solids` explodes each free shape topologically
(`TopExp_Explorer(shape, TopAbs_SOLID)`) and uses the label walk only to *attach names* by
`IsPartner`. Exploding to `TopAbs_SOLID` is precisely the step the converter's leaf branch was
missing.

## 3. The fix

`O2_CADtoTGeo.py`, `expand_definition`'s `IsSimpleShape` branch:

- A leaf label whose shape holds **more than one** solid body becomes an **assembly**, and each
  body becomes its own logical volume, placed once under it at the identity. The body's own
  `TopLoc_Location` is already baked into the `TopoDS_Solid` the explorer returns, so the identity
  placement is exact — no transform has to be recovered or composed.
- Body definitions are keyed `<label entry>#b<n>` (`0:1:1:4#b7`), which sanitises to
  `0_1_1_4_b7` in filenames and C++ identifiers. The display name stays the parent label's name
  unchanged, so BOM matching by exact part number behaves exactly as before.
- A leaf label with **no geometry at all** is dropped with a warning instead of being handed to
  `triangulate_asbbox`. That is the `Bnd_Box is void` failure, fixed at its source rather than
  worked around: it was a symptom of leaf granularity, as §5 of `Handoff_IntegrationTest.md`
  suspected.
- A leaf label with exactly **one** body keeps the bare label entry as its definition key and its
  original shape object, untouched. This is what makes the change free for every other corpus.

Stream_W's deduplication semantics are unchanged and untouched: the rule still keys on
(definition, world transform), it now simply sees the twenty real definitions. On oTOF it
suppresses nothing, because all 62 628 placements are at distinct world transforms — which is what
`Stream_W_DoublePlacement.md` measured of this file in the first place.

Two helpers were factored out for it: `solid_bodies_of(shape)` and `_register_leaf_shape(...)`,
the latter holding the volume/shape/triangles bookkeeping that both paths share.

### Geometry is preserved exactly

The decisive check is the world bounding box, computed two ways — from the CAD's own free shape,
and from the emitted graph by composing every placement and transforming every definition:

```
CAD world bbox (mm):   -1.509 -108.719 -962.019  3477.616 1815.319  962.019
graph world bbox (mm): -1.509 -108.719 -962.019  3477.616 1815.319  962.019
max |delta| (mm): 0.0
```

62 628 placed bodies, bit-identical extent.

## 4. Regression evidence

The traversal is shared by every corpus, so all three gates were re-run.

**Self-test.** The four existing suites are unmoved at **18 + 8 + 10 + 12 = 48 checks, 0
failures**. A fifth suite, `run_multibody_leaf_self_test` (**6 checks**), was added rather than
folded into an existing one, precisely so the 48 stays readable as 48. It covers a two-body leaf
instanced twice (2 volumes, 4 placements, nothing suppressed), the single-body control (one
volume, keyed on the bare label entry, no `#b` suffix), and an empty leaf label (dropped, siblings
still convert, no exception).

**Bagger and ALICE3, reconverted** (`--exact-surfaces auto --csg auto`, no `--mesh`), before and
after the fix, into separate folders:

| | Bagger before | Bagger after | ALICE3 before | ALICE3 after |
| --- | --- | --- | --- | --- |
| leaf definitions | 13 | 13 | 55 | 55 |
| leaf placements declared -> emitted | 13 -> 13 | 13 -> 13 | 206 -> 103 | 206 -> 103 |
| `surfaces_*.bin` | 13 | 13 | 20 | 20 |
| `facets_*.bin` | 13 | 13 | 55 | 55 |
| `shape_*.root` | 7 | 7 | 2 | 2 |
| `csg_*.json` | 13 (+report) | 13 (+report) | 55 (+report) | 55 (+report) |
| tiers | CSG 7 / surf 6 / tess 0 | same | CSG 2 / surf 18 / tess 35 | same |

`geom.C`, `csg_report.json` and **every** `.bin` sidecar are byte-identical modulo the output
path; the `.root` shape files differ only in ROOT's embedded timestamps. `diff -rq` reports
nothing else.

**Bagger oracle gate**, `runOracleGate.py --model STEP_examples/Bagger.step`:

```
surface  contains=0  distout=0  distin=0  safety=0   (13/13 parts with zero disagreements)
shape    contains=0  distout=0  distin=0  safety=0   (7/7 parts with zero disagreements)
mesh     contains=419 distout=7364 distin=8512 safety=10883   (0/13)
EXIT=0
```

Every column identical to the state recorded in `NEXT.md`.

## 5. oTOF, converted

`--exact-surfaces auto --csg auto`, **no** `--mesh`, one process:

| | |
| --- | --- |
| leaf definitions (prototypes) | **20** |
| leaf placements | **62 628**, multiplicity `{1: 62628}` — nothing suppressed |
| dropped | 1 empty leaf label (`A Side`), with a warning |
| exact surfaces | **20/20** extracted, 0 fall back to tessellation |
| CSG | **19/20** accepted, all `tier1-box` `TGeoBBox`, `dV_sym = 0 cm^3` on every one |
| the 20th | `oTOF v2 v1 #b3`: *declined CSG: 1493 planar faces: not a six-plane box*; ships exact-surface |
| files | 20 `surfaces_*.bin`, 20 `facets_*.bin`, 19 `shape_*.root`, 20 `csg_*.json` + `csg_report.json` |
| `geom.C` | 2.3 MB, **3 964 `AddNode` calls** |
| wall time / peak RSS | **44.4 s / 753 MB** |

With `--mesh --mesh-prec 0.1` on the whole model: **44.1 s / 775 MB**, output 3.4 MB. Meshing
twenty small prototypes costs nothing; it is placements, not prototypes, that oTOF has a lot of.

Two things worth reading off that table. The 19 exact `TGeoBBox` prototypes are exactly what
`Stream_A_CSG.md` predicted from the census (19/20 one-cell) and exactly the case
`Stream_N_PlacedPrimitives.md` made cheap — the placed-primitive path fires on the whole corpus.
And **3 964 `AddNode` for 62 628 placed bodies** is TGeo doing what it is good at: the assembly DAG
is shared, so a 62 628-body detector is a 2.3 MB macro rather than a 30 MB one. This is the
one-volume-many-nodes case, and it is now reachable.

The model is not centred on the ALICE origin — x runs 0 to 347.8 cm, z is +/-96.2 cm — so the
rotate-and-shift step of `Handoff_IntegrationTest.md` §6 still applies before it can be placed in
`centralBarrel`.

## 6. Website hand-off

Five representative prototypes were copied into `website/testdata/` with
`fetch_testdata.sh --append --all --group "oTOF (...)"`, which merges by part name and left the 17
existing parts untouched (manifest 17 -> 22). No file under `website/` was edited.

| part | what it is | ships |
| --- | --- | --- |
| `oTOF_v2_v1_0_1_1_3_b1` | `TGeoBBox(173.9, 0.01, 2.656)`, the long support strip | CSG shape + exact + mesh |
| `oTOF_v2_v1_0_1_1_3_b3` | 1493 planar faces, the scaling case | exact + mesh, CSG declined |
| `Module_v1_0_1_1_4_b1` | `TGeoBBox(1.6, 0.015, 1.25)`, placed 3 672 times | CSG shape + exact + mesh |
| `Module_v1_0_1_1_4_b6` | `TGeoBBox(1.6, 0.015, 0.1)`, the thin layer | CSG shape + exact + mesh |
| `Module_v1_0_1_1_4_b15` | `TGeoBBox(6.43, 0.005, 2.51)`, the carrier | CSG shape + exact + mesh |

`fetch_testdata.sh`'s `partname()` strips four trailing numeric segments, and `_b1` is not
numeric, so these keep their full stem as the part name. Ugly, correct, and unique — the website
script was deliberately not touched to prettify it.

## 7. What this does not fix, and one environment note

- **The 1493-face part is still the scaling case.** `oTOF v2 v1 #b3` has 1993 concave edges
  (`Stream_A_CSG.md` §2.5) and is not a CSG candidate; it ships as an exact surface solid and is
  the natural stress test for `CloseShape` cost against face count.
- **oTOF is not in the oracle gate.** The gate was run on Bagger only, as the regression control.
  Putting a 20-prototype oTOF through `runOracleGate.py` is the obvious next measurement, and the
  corpus is now reachable for it.
- **No transport was attempted.** `Handoff_IntegrationTest.md` stages 1-3 for oTOF are still open;
  what is unblocked is stage 0.
- **Do not run two `--csg auto` conversions at once.** Two converters started in parallel produced
  `[WARN] ...: accepted as CSG but NOT emitted -- ROOT unavailable` in one of them and zero
  `shape_*.root` files, while the serial re-run wrote all seven. The loud warning added for
  `NEXT.md` item 2 is doing its job, but the contention itself is real and silent-looking in a
  log tail.
