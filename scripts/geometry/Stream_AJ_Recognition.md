# Stream AJ — the recognition programme: revolved profiles, the prism family, single cells

Executed 2026-08-23/24 from [`Handoff_Recognition.md`](Handoff_Recognition.md), in three rungs,
each landed and reviewed before the next began. Commits, in order: `4d1780035c` (rung 1, the
revolved-profile recogniser → `TGeoPcon`), `e278709fd1` (`checkKnownSource.py`, the third
acceptance test), `a27568ec3c` (rung 1 fix round: the acceptance retry and the native-cone
canonicalisation), `0a013b95f5` (rung 2, the prism family → `TGeoTrd1`/`Trd2`/`Arb8`/`Xtru`/
`Pgon`), `389ed4befb` (known-source one-way scoring for multi-body labels), `c66a46b991` +
`b69715766f` (rung 3, the single-cell emitter). Detailed engineering reports for each rung were
written during the work; this document is the programme record: what converted, on what evidence,
and which premises died on the way.

## 1. The headline

**Every native primitive class the writer emits is now recognised at 100 % on all five detector
corpora.** PIPE 58/58 `TGeoPcon`, ITS 23/23 `TGeoXtru`, TPC 32/32 `TGeoTrd1` and 8/8 `TGeoPgon`,
TRD 47/47 `TGeoTrap`, and every `TGeoTube`/`TGeoBBox`/`TGeoTubeSeg`/`TGeoCone` there was. What
still declines is exactly the honest remainder: genuine boolean composites (Tier 3, deliberately
not built), PIPE's 22 `TGeoEltu` (elliptical cylinders arrive as free-form carriers — a possible
future rung the handoff never demanded), toroidal parts, and free-form CAD.

Three acceptance tests gate every conversion, none waived anywhere: the OCCT symmetric difference
(`dV_sym` within the model-tolerance band — **0 exactly on every part this programme added**),
the oracle gate, and the new known-source check — **1655 of 1655 CSG parts across the five
corpora agree with the original `TGeoShape` they were generated from, zero containment
disagreements** (20 000 seeded points per part).

And the fixtures gate, which carried two sliver `distout` rays since Review Appendix A, **now
exits 0**: `tube_window` (`TGeoTube - TGeoTube`), `cyl_inter_cyl` (`TGeoTube ^ TGeoTube`) and
`oblique_cut_cyl` (`TGeoTube ^ TGeoBBox`) all ship as csg with clean shape columns, so no sliver
ray is left in any shipped verdict.

## 2. Per-corpus tiers, before → after each rung

csg / surface / mesh per corpus. "Before" is the regenerated 2026-08-23 baseline
(`o2-sim-serial -n 0 -g boxgen -m <MOD>` → `O2_TGeoToCAD.py` → `O2_CADtoTGeo.py
--exact-surfaces auto --csg auto`), not the handoff's stale numbers (§6.1).

| corpus | leaf solids | before | after rung 1 | after rung 2 | after rung 3 (final) |
| --- | ---: | --- | --- | --- | --- |
| PIPE | 176 | 62 / 91 / 23 | 120 / 33 / 23 | 122 / 31 / 23 | **128 / 25 / 23** |
| ITS  | 259 | 152 / 107 / 0 | 171 / 88 / 0 | 197 / 62 / 0 | **207 / 52 / 0** |
| TPC  | 172 | 82 / 90 / 0 | 100 / 72 / 0 | 145 / 27 / 0 | **146 / 26 / 0** |
| ABSO | 29  | 9 / 20 / 0 | 24 / 5 / 0 | 26 / 3 / 0 | **26 / 3 / 0** |
| TRD  | 1170 | — (corpus new at rung 2) | — | 1148 / 22 / 0 | **1148 / 22 / 0** |

Across the four baseline corpora the CSG population went **305 → 507** (+202); with TRD the five
final corpora carry **1655** CSG parts. Every tier change was one way, `surface`/`mesh` → `csg`;
**zero parts left csg and zero previously accepted parts changed a single byte of candidate,
evidence or placement at any rung** — asserted by part-by-part report diffs at every stage and by
three SHA-256 digest tables in `csg/emit.py --self-test` (8 whole-part, 26 prism, 5 cell
candidates frozen).

Final recogniser breakdown over the 1655: tier1-box 1063, tier1-tube 182, rung2-trd1 110,
revolved-pcon 109, rung2-arb8 57, tier1-tubeseg 54, rung2-xtru 47, cell-intersection 17,
tier1-cone 12, revolved-cone 2, rung2-pgon 2.

## 3. Demand versus delivered (the handoff's §2 numbers)

| rung | demand (decline counts, handoff) | delivered |
| --- | --- | --- |
| 1 revolved → Pcon | PIPE 46, TPC ~19, ITS 18, ABSO 16 | PIPE **+58**, TPC **+18**, ITS **+19**, ABSO **+15** |
| 2 prism family | TPC ~44, ITS 47 (23 Xtru), TRD 149, MAG 14 | TPC **+45**, ITS **+26** (24 Xtru), TRD **141 in-corpus** (§6.4), MAG not run (§6.5) |
| 3 single cells (gated) | `tube_window` + `cyl_inter_cyl` | both, plus `oblique_cut_cyl`, `box_minus_cyl`, `sphere_minus_cyl`, and **17 corpus parts** |

Requirement bullets, each landed on real geometry: mixed cone+cylinder laterals (`IBCYSSCone` →
`TGeoPcon(nz=6, phi1=180, dphi=180)`); z-steps (duplicate z planes, `TPC_M__body` nz=30); `rmin`
through 0; half-turn and partial-phi wedges (`OBConeSideA` phi1=1.278 dphi=177.4); Pgon-vs-Pcon
disambiguation (a 48-edge hollow polygon — apothem planes within 0.2 % of a cylinder — comes out
`rung2-pgon`, never `revolved-pcon`, asserted both ways); slanted opposite faces (TPC's 44 "no
opposite partner" declines, including `TPC_IRB1`'s 0.5 % slant that an angular criterion would
have called a box); non-convex Xtru polygons (18 of ITS's 23); annular Pgon sections
(`TPC_Strip`, a 0.1 mm wall on an 85 cm radius, 18 edges).

The 17 rung-3 corpus parts deserve their own sentence: **all seventeen were `TGeoCompositeShape`
in the O2 geometry to begin with** — `ITSServicesWaterOB` is `TGeoTube - TGeoTube - TGeoTube`,
`ConeAStep0` is `TGeoBBox ^ TGeoBBox ^ TGeoTube` — so the single-cell emitter is not inventing
booleans, it is recovering the ones O2 itself wrote, through a STEP file that had thrown the
boolean tree away. The `Stream_AA` §5.2 falsifier was never hit: on 22 genuinely single-cell parts
measured, none was rejected by volume or gate.

## 4. The machinery that now exists

- `csg/recognise.py`: `_match_revolved` (one axis, any number of z sections → `TGeoPcon`, with
  2-section full-turn profiles canonicalised to `TGeoCone`/`TGeoTube`), `_match_prism` (all-planar
  face graph → section stack, templates most-specific-first `trd1, trd2, pgon, arb8, xtru`, corner
  order propagated through the side-face adjacency — non-convex sections cannot be read by angle
  sorting), `_match_single_cell` (no trusted concave edge → one intersection of the part's own
  oriented halfspaces, folded to at most `_CELL_MAX_LEAVES = 8` bounded native leaves). Every
  matcher runs strictly after everything that existed before it, and only on declines.
- Every matcher decides by **one measured quantity** (the Stream_K rule): a geometric gap in cm
  over the part's bounding-box diagonal against `REL_TOL = 1e-6` — `_profile_gap` (point to
  profile polygon), `_point_set_gap` (symmetric Hausdorff on corners *and edge midpoints* — the
  midpoints are what catch a hexahedron read in the wrong corner order), `_boundary_gap` (samples
  against the *other solid's boundary*, because two constructions of the same transcendental rim
  parametrise it differently). No per-class criteria, no angles. Each instrument is self-tested
  to report exactly a known displacement, so it can measurably say "no".
- `emit.process_solid` retries `recognise_revolved` then `recognise_single_cell` **only after the
  acceptance test rejects** a candidate — a part accepted first time never reaches a retry, so no
  accepted part can change. This is how an all-cone stack (proposed wrongly as one `TGeoCone`)
  and `tube_window`/`cyl_inter_cyl` (proposed wrongly as two-tube unions) become reachable at all.
- `csg/primitives.py`: six new leaf types with array-valued parameters (per-type length groups —
  `TGeoXtru` carries an independent polygon count and section count), an `outside` leaf flag
  (`op: "intersection"` folds to `TGeoIntersection`/`TGeoSubtraction` chains), and **both**
  builders realising every description from the same data — `BRepPrimAPI_MakeRevol` for the Pcon,
  sewn planar faces for the prisms, bounded halfspace leaves for cells (`TGeoHalfSpace` is
  deliberately not used: it has no OCCT counterpart and would break the one-description/
  two-builders invariant).
- `scripts/geometry/checkKnownSource.py` — the programme's third acceptance test, new: for every
  part the cascade ships as CSG, find the source `TGeoVolume` through the writer report
  (`emittedName`, `bodyComponent` for `__body` solids, `name#N` disambiguation by class+capacity,
  `__mirrored` Z-reflection), then compare class, `Capacity()` (to 1e-9 where both analytic), the
  Pcon profile where directly comparable, and a seeded 20k-point `Contains` cross-check through
  the recorded `shapePlacement`. One body of a multi-body CAD label is scored one-way (inside the
  body ⇒ inside the label), with an encloses-nothing guard so the weaker rule cannot pass empty.
  Its own self-test (15 checks) includes deliberately wrong shapes that must be caught.

Self-test growth over the programme: `csg/emit.py --self-test` 33 → **184** (11 acceptance +
173); `checkKnownSource.py --self-test` — → **15**; `O2_CADtoTGeo.py --self-test` 54 unchanged;
`O2_TGeoToCAD.py --self-test` 105 **untouched** (the writer stayed frozen ground truth);
`runOracleGate.py --self-test` 17/17; ctest `BVHSurfaceSolid` 113 + `BVHAssembly` 22 green
throughout. Bagger gate 13/13 exit 0 with its CSG 7 bit-identical at every stage.

Negative controls now in the suite, because coverage never buys correctness: near-miss profiles
and prisms displaced by ten model tolerances (refused by gap or volume, both legs exercised); a
V-notch ladder at 2e-3/1e-5/1e-6/1e-7 rad showing which of the three tests catches a two-cell
body at each depth — and that a notch one model tolerance deep is correctly *accepted*; Pgon
pseudo-revolutions; a cylinder with a coaxial hexagonal section; a twisted hexahedron; illegal
descriptions refused before either builder sees them.

## 5. Every surviving decline, and what it is

The remaining non-csg populations, from the final reports (each part carries its full recorded
reason; `website_data/decline_reasons.json` regenerated for Bagger + ALICE3 + the fixtures):

| corpus | surface+mesh left | composition |
| --- | ---: | --- |
| PIPE | 48 | 22 `TGeoEltu` (elliptical carrier = free-form to the face reader), 16 toroidal (bellows plies), 23-part mesh population unchanged, plus a handful of multi-axis composites and one sliver-band near-miss |
| ITS | 52 | space-frame/cage boolean composites (16–74 planes), multi-cluster service composites, the two `IBConnectorBlockTail*` blocks refused by the 8-leaf cell budget at 10 halfspaces |
| TPC | 26 | `IRCAP`/`ORCAP`/`OCGEM` composites (a cap with two holes is out of the cell emitter's scope), `TPC_hook`, multi-cluster bodies |
| ABSO | 3 | two mixed-plane composites, one 6-cap stack of several z-groups |
| TRD | 22 | 18 two-body BTOF labels (each body a genuine boolean), `BREF1`, `VolTOFrail`, `B077__body`, `barrel__body` |

Every one is a genuine multi-body or boolean composite or a carrier class with no template —
Tier-3 territory the handoff deliberately excluded, or the `TGeoEltu`/torus classes it never
demanded. None is a missed single primitive.

## 6. Dead premises, recorded

1. **The handoff's regression floor "ITS 128, TPC 82" was stale for ITS**: the regenerated
   baseline says 152 (the writer fixes since `Stream_AI` moved it). The floor was applied against
   the regenerated baseline — the diff is what matters, and it was zero at every rung.
2. **"Most CSG decliners are single primitives" — confirmed, and undercounted for PIPE**: the
   demand table said 46 Pcons; the corpus held 58. For ITS's prism demand the opposite: "47"
   counted boolean composites among the decliners; the true single-prism population was 26 (the
   remaining ~23 are genuine space-frame composites).
3. **"Retiring `tube_window` and `cyl_inter_cyl` retires the two sliver rays"** — almost: three
   fixtures carried sliver rays across the surface columns (`oblique_cut_cyl` was the third). All
   three now ship as csg, and the gate exits 0.
4. **TRD's "149 real decliners"**: 141 prism-family parts existed in the corpus; the remainder
   the frozen writer never emitted, chiefly the degenerate-section `TGeoTrd1` (`dx1 = 0`,
   `B045cut`) which `O2_TGeoToCAD._prism_from_rings` declines as "[2, 4] distinct vertices" — a
   writer work item, unreachable by any recogniser until the writer is unfrozen.
5. **MAG (demand 14) was never converted**: the handoff's own corpus list (§5: "PIPE, TPC, ITS
   first; TRD if time") never scheduled it. Its demand classes (`Trd1`/`Pcon`) are the same
   templates now at 100 % elsewhere; running the corpus is a one-command follow-up.
6. **Three negative controls flipped to conversions** as the rungs landed (the blind bore →
   polycone, the L-plate → Xtru, the milled flat → cell): a control written when nothing could
   convert a shape stops being a control once something can. Each flip was ruled on explicitly,
   each candidate is digest-frozen, and no check count shrank.
7. **`concaveEdgesTrusted == 0` is necessary but not sufficient** for a single cell: a V notch
   below the census's near-tangential filter (`NEAR_TANGENTIAL_SIN = 1e-3`) is a two-cell body
   the trust filter reports as one. The measured gap catches it at 1e-5 rad, the symmetric
   difference at 1e-6 rad; nothing wrong can ship, but `Stream_AA` §2's `concave (trusted)`
   column has that blind band in it and the volume is what closes it.
8. **Twisted `TGeoArb8` is out of scope by measurement**: the writer emits its sides as ruled
   `brepfill` B-spline patches, so `_face_records` declines the solid before any matcher sees it
   (reason names the case). Zero corpus parts affected; the description format already carries
   the vertices, so ruling it in later is a recogniser change only.

## 7. Follow-ups this programme leaves on the table

- **`TGeoEltu`** — PIPE's largest surviving single-class population (22 parts): an elliptical
  revolved rung, needing the face reader to learn the elliptic-cylinder carrier first.
- **MAG corpus** run + score (one command each; §6.5).
- The **cell leaf budget** (8) refuses ITS's two 10-halfspace connector blocks — raise it only
  with query-cost evidence in hand.
- `checkKnownSource.py` **class-specific parameter comparison** for the prism classes (an Xtru
  polygon deviation in the shape of `pcon_profile_deviation`); capacity + containment carry the
  verdict today.
- `crosscheck_bbox` on a rotated cell composite reports 1–3 cm — `TGeoIntersection::ComputeBBox`
  taking the axis-aligned hull of a rotated operand's box, a ROOT property, not a geometry
  defect (the sharp checks are 0/4000 and 0/20000 on every such part). A tighter box needs a
  tighter-than-box halfspace leaf.
- The prism/cell matchers run on every declining part; TRD (1170 leaves) converts in ~12 min. If
  a corpus ever makes that painful, a face-count cap of the census's kind
  (`carrier_face_cap = 40`) placed before the boolean build is the cheap fix.
- The **writer's degenerate prism sections** (already NEXT's open list) now measurably block the
  one demand case this programme could not reach (§6.4).
