# W2 — the recognition scoreboard, and what is traceable

Run 2026-08-28. Source: `~/roundtrip_report/alice_roundtrip.json`, the R6 breadth census written by
`roundTripReport.py` on 2026-08-25. Regenerate everything with
`python3 make_scoreboard.py` (writes `tables/` and `numbers.json`).

## What the corpus JSON contains

`parts` is a list of 22 901 records, one per leaf solid, each carrying its module, part and source
volume name, the **source** `TGeoShape` class, boolean depth and leaf-class histogram, which tier it
**ships** as (`csg` / `flatcsg` / `surface` / `mesh`), the recogniser that claimed it, the
acceptance numbers (`dVsym`, `band`, `capacityRelativeDeviation`, `containsMismatches` out of
`containsPoints`), `whyNotCSG` when it declined, the exact-tessellation predicate and its reason,
a surface census, and the artefact sizes. `modules` is 17 records with the per-module tier counts,
the writer's per-shape-class converted/declined tallies with reasons, and the known-source score.

That is enough for every table on slides 3, 7 and 8 without re-running anything. The per-module
`conv/` directories are **gone** — only `step/` survives next to the report — so anything not in
this JSON does need a re-run.

## Every headline number, recomputed

| claim | recomputed | agrees with `NEXT.md` |
| --- | --- | --- |
| leaf solids in the census | 22 901 | yes |
| native CSG tree | 22 736 (99.28 %) | yes |
| `O2FlatCSG` | 20 | yes |
| exact surfaces | 130 | yes |
| tessellated | 15 | yes |
| CSG of any kind | 22 756 (99.37 %) | consistent — `NEXT.md`'s 99.3 % is the native figure alone |
| agrees with its source `TGeoShape` | 22 747 / 22 756 | yes |
| tessellation is exact (all-planar) | 19 622 (85.7 %) | yes |
| writer declines, never reaching the converter | 31 volumes | yes |
| writer declines by class | `TGeoPara` 13, `TGeoTorus` 6, `TGeoCompositeShape` 12 | yes |

## One number on the wanted-list is NOT reconstructible: 163/188 (87 %)

`NEXT.md` and `Plans.md` quote **163 of 188 composite-sourced detector parts (87 %)**. That
population is not in this corpus: it is the older five-module study (PIPE/ITS/TPC/ABSO/TRD),
counted over distinct parts, from before R6. The corpus-wide equivalent, recomputed here, is
**2 987 of 3 125 composite-sourced parts (95.6 %)** over 1 963 distinct source volumes.

**Decision for the talk: quote 2 987 / 3 125 = 95.6 %, not 87 %.** It is the stronger number, it is
over the whole of Run 3 rather than five modules, and it traces to a record. The 87 % figure is
retired from the slides rather than re-derived.

## Where the remaining 138 declines go

`tables/T3a_converter_declines.md` groups them; `T3a_converter_declines_raw.md` keeps every reason
verbatim. The shape of it: 57 parts exceed the cell matcher's axis-cluster budget, 33 are partial
tori, 28 disagree with their own loop twin, and the long tail is single digits. Exactly **one** part
in all of Run 3 declines for genuinely free-form surfaces — the gap that sounds largest in the
design documents is the smallest one in practice, and that is worth saying out loud on slide 15.

## Caveats that must travel with these tables

- The tier a part **ships** as is a routing decision, not a statement that no better representation
  exists. 19 622 parts have an exact tessellation and still ship as CSG, deliberately — a box is a
  `TGeoBBox`.
- The nine known-source "failures" are all MFT bodies of a multi-body label at ~6 ppm of the
  label's volume, where the instrument samples the source's bounding box and cannot score an empty
  comparison. Re-scored in the body's own box they are clean. Slide 8 says 22 747 / 22 756 **and**
  names this, or it says nothing about the nine.
- Any speed ratio quoted anywhere must say *against what*: the flat solid's 3-34x is against the
  same cells as a plain composite, not against what O2 ships.
