# Handoff — the flat-CSG programme: finish composite recognition, then TGeoBVHCSG

**Written 2026-08-24, Sandro's decision at the close of the recognition programme
(`Stream_AJ_Recognition.md`).** The flat-CSG representation is pushed **in any case** — the
decision is no longer gated on per-detector demand, because ITS and TPC alone are too small a
sample to decide a representation on, and the substrate serves every future CAD-native geometry.
The stated goal: **100 % CSG recognition for every part that came from CSG in the first place**,
then `TGeoBVHCSG` as the optimisation layer for complex composites — in the sense that Geant4's
`G4MultiUnion` already attempts, but general: a flat DNF (union of intersection-cells over
halfspace primitives) with a BVH over cell AABBs, not union-of-solids only.

Read first: `Stream_AJ_Recognition.md` (what exists), `Stream_AA_FlatCSG.md` (the measured
plan this executes — its §5 gate 2, the single-cell emitter, is DONE), `Stream_K_Tier0.md`
(the scoring rule), `NEXT.md` (state + traps). The three-test acceptance discipline of
`Handoff_Recognition.md` §3 binds every rung below unchanged.

## 1. The measured starting point (2026-08-24, rung-3 corpora)

Composite-sourced parts recognised: PIPE 7/31, ITS 12/69, TPC 1/27, ABSO 0/3, TRD 36/58 —
**56/188, 30 %** (ITS regenerated 2026-08-24 after the depth-guard fix; its five deep cages
are now in the corpus). What the remaining 127 are, from the reports and the source-tree measurement:

- **~20 torus-carrying single cells** (PIPE's 16 bellows plies at 2–7 faces, `RB24ValveMA2`,
  RB24/RB26 plies) — blocked only by the missing torus carrier, not by decomposition;
- **deep connected unions**: ITS space-frame cages at depth 18–60 / up to 61 leaves, TPC
  `rodl`/`rods` (depth 13/11), TRD `TOFrail` (22 leaves), TRD's 18 two-body BTOF barrels;
- **multi-hole caps** (TPC `IRCAP`/`ORCAP`/`OCGEM` — a cap with >1 hole is outside the current
  cell emitter), multi-axis service composites, and the two ITS connector blocks refused by the
  8-leaf cell budget;
- PIPE's one genuinely twisted `ARB8` composite (ruled sides — stays out, recorded).

Separately, PIPE's 22 `TGeoEltu` are not composite-sourced but are the largest surviving
primitive population anywhere and belong in this programme's first stage.

## 1.1 Progress

**R1 + R2 are DONE** (2026-08-24, commits `9b4708b4c5`, `67ac71a315`, `d138313796`): PIPE
128 → **162 csg / 13 surface / 1 mesh** (+22 `tier1-eltu`, +1 `tier1-torus`, +11 bellows plies as
torus-carrying cells, all at `dV_sym = 0`, class and capacity matching their sources to 1.1e-14);
zero movement on ITS/TPC/ABSO/TRD; known-source 1689/1689; fixtures and Bagger gates exit 0
unchanged. The composite-sourced rate moved 56/188 → **67/188 (36 %)**. ALICE3 gains one
conversion (`ST2487455_01`, a torus-carrying cell) and now declines its fillet blends readably —
a crash on its self-intersecting tori forced `InvalidDescription` (a typed "the numbers are not a
legal solid") to be declined at every proposal site, with a bad matcher still raising. The
same-name resolution in `checkKnownSource.py` now uses the exact bounding box, never a
Monte-Carlo capacity. **R3 + R4 are DONE** (2026-08-24/25, commits `5e57015d4f`, `09ffc512b3`, `7d3a203504`): Tier-0
canonicalises 1082 ALICE3 faces (1014 at ≤1e-7 relative — gate 1 passed); the decomposition rung
converts **+86 detector parts** (union-of-cells, every one at dV_sym exactly 0, Contains
0/343 999) — corpora now PIPE 166/9/1, ITS 250/14/0, TPC 163/9/0, ABSO **29/0/0**, TRD 1167/3/0,
known-source **1775/1775**; ALICE3 3 → 6; the fixtures gate's shape column is **10/10** (all ten
fixtures ship as csg, `torus_union_cyl` included). Composite-sourced rate **153/188 (81 %)**. A
false accept was caught and closed (`accept.contains_disagreements`: OCCT `Cut` can return
`IsDone()` with zero solids both ways — `ST0923290_01#b19`). R5's demand is measured: 21 parts
over the cell/leaf budgets (task-r4-report §7), plus 24 declines where a dihedrally-convex piece
is not a carrier-arrangement cell — the split-at-every-carrier-crossing trade is R5's to price.
Next: **R5 (`TGeoBVHCSG`)** and/or **R6 (the census)** — R6 is cheap and sharpens R5's demand.

## 2. The rungs, in order

**R1 — torus carriers (small).** `TGeoTorus` leaf in `csg/primitives.py` with BOTH builders
(OCCT has the primitive; ROOT class is native; the kernel already solves quartics for the
surface representation); a tier-1 whole-torus template; the torus admitted as a carrier in
`_halfspace_carriers`/`_cell_leaf` (a torus halfspace is already bounded — like the sphere it
needs no extension). Retires PIPE's plies exactly. Demand: ~20 parts, PIPE-dominated.

**R2 — `TGeoEltu` (small).** Recognise the extruded-ellipse lateral (the carrier arrives
free-form today — read the basis curve), two caps → native `TGeoEltu(a, b, dz)`. Moves 22 PIPE
parts mesh → csg. An elliptic patch type in the kernel (mesh → surface fallback for whatever
declines) is a separate, optional work item — not this rung.

**R3 — Tier-0 canonicalisation** (`Stream_AA` §5 gate 1, still open). Pure converter work,
prerequisite for the splitter on real CAD. Falsifier as stated there: fewer than ~1000 ALICE3
faces decoding at ≤1e-7 relative gap, or any existing conversion moving, stops the rung.

**R4 — decomposition → union-of-cells emission.** `BRepAlgoAPI_Splitter` on trusted concave
edges (the ~150-line probe loop `Stream_AA` §2 already measured at 15/16 on ALICE3), then each
cell through the EXISTING `_match_single_cell`/fold machinery, emitted as a balanced
union-of-cells `TGeoCompositeShape` (`primitives` already has N-leaf union and the `outside`
flag). The per-part cell budget replaces the per-cell leaf budget; plain-composite emission is
the honest baseline the next rung must beat. Acceptance unchanged: symmetric difference, oracle
gate, known-source (composite sources score by containment; capacity is Monte-Carlo and stays
flagged not-comparable). Known risks, all recorded with evidence: OCCT splitter volume drift
(2/16 on ALICE3 — the acceptance rejects those, they fall one tier, nothing loosens); the
near-tangential blind band of the trust filter (the notch ladder — the volume closes it);
the writer's depth-32 chain constant is FIXED (`6241b9b118`, `MAX_BOOLEAN_DEPTH = 512`):
ITS's depth-35-to-60 cages convert to B-rep, ship as surface solids, and decline CSG with
measured cell structure (EndWheelCBasis* at only 4 trusted concave edges each — cheap first
splitter targets; IBCYSSFlangeA at 47).

**R5 — `TGeoBVHCSG`.** The new C++ shape class, only now: flat DNF cell blob (fixed-length
halfspace coefficient blocks — deliberately the AOT-codegen-friendly layout), BVH over cell
AABBs (the substrate validated three times on this branch: patch BVH, sub-patch BVH,
`O2BVHAssembly`), `Contains`/distances by per-cell interval clipping (no parity, no trims,
no CloseShape), `Safety`, analytic `Capacity` as the sum of disjoint cell volumes, ROOT IO,
gate integration (the gate scores any `TGeoShape` — nothing to change there), X-ray bench
scoring against the same parts emitted as plain composites, and a `TVirtualGeoConverter` seam
for Geant4 (`G4MultiUnion` where the structure is a pure union; the general case needs its own
mapping and is part of this rung's design work). Emission policy: cells above a measured
crossover (Stream_AA derives ≥10; verify) ship as `TGeoBVHCSG`, below it as plain composite.
This is the weeks-scale rung; it is the point of the programme, not a maybe.

**R6 — breadth: the full composite census.** Convert and score the remaining Run 3 modules
(MAG, TOF, TRD is done, MFT, MCH, MID, FT0, FV0, ZDC, EMC, PHS, CPV, HMP — whatever
`o2-sim-serial -m` accepts), one command each now, so the 100 %-of-composite-sourced claim and
the BVHCSG crossover are measured on the whole geometry rather than on two detectors. This rung
can run early and cheaply in parallel with R1/R2 (conversions strictly serial as always).

## 3. Definition of done

- The composite-sourced recognition rate (§1's 56/188 baseline, extended by R6's census) is
  **100 % minus an explicitly named remainder** — every surviving decline carries a recorded
  reason a reader can act on (twisted ruled sides, scaled shapes, whatever R6 surfaces).
- Every conversion passes the three tests; every floor of `Handoff_Recognition.md` §4 holds
  (self-tests grown never shrunk, digest tables green, Bagger bit-identical, byte-identity of
  previously accepted parts at every rung).
- `TGeoBVHCSG` exists, is gate-scored and X-ray-benched against plain-composite emission of the
  same parts, with the crossover measured and the emission policy recorded.
- `Stream_AJ` gets a successor stream document; `NEXT.md` is rewritten at each session close.

## 4. Environment (verbatim; eval+command in ONE shell invocation)

```
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
SW=$HOME/alisw/sw/ubuntu2404_aarch64
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
# sim env (geometry regeneration) stays separate: export O2_ROOT=$B/stage, o2-sim-serial, no OCC paths
```
Standing traps: all of `NEXT.md`'s (serial conversions, deletion hook, one ninja, detached
long runs). R5 touches `Detectors/Base/**` — the kernel freeze of the recognition programme
does NOT apply here; ctest floors do.

## 5. Ordering against the closure test

The closure test (`Handoff_ClosureTest.md`) remains queued and is unaffected: it scores the
representations that exist. Sandro's call whether it runs before R1 or after R1/R2 — R1+R2
would let the closure test's "everything CSG" variant carry PIPE at ~145/176 csg.
