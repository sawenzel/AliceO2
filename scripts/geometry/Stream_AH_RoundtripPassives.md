# Stream AH — the TGeo → STEP → TGeo round trip on TRD, MAG and ABSO

`Stream_AD_TGeoToStep.md` closed the round trip on the beam pipe and left the bench pointed at
the rest of the experiment. This document runs the same bench on three more modules — the TRD
detector, the L3 magnet (`MAG`) and the front absorber (`ABSO`) — because each was expected to
stress a different part of the two converters: TRD is `TGeoTrd1`-heavy and gives the polyhedral
face-sewing path its first hard workout, and MAG and ABSO are `TGeoPgon` and `TGeoPcon`/`TGeoCone`
stacks, which is exactly the corpus Roadmap item (h) asks about.

The method is Stream AD's, unchanged: the source `TGeoShape` is an exact oracle because the STEP
was generated from it, so every disagreement is a defect on one side or the other and there is no
tolerance argument to have about the model. What is new here is a **fourth instrument** — a
placement-fidelity walk (§2, instrument D) — because the three Stream AD instruments all score a part in its
own frame and are structurally blind to a part that arrives at the wrong place in the world.

Measured 2026-08-23 on `swenzel/bvhsurfacesolid`, aarch64, OCCT/pythonOCC 7.9.0, ROOT 6.36.10.
Every number below was run on this machine.

**Verdict up front.** Both converters carry the geometry correctly where nothing is mirrored —
**1 325 548 Contains samples** across the three modules, **zero disagreements**, and capacity
agreement at the double-precision floor. But four defects came out of it that PIPE could not have
shown, all in `O2_TGeoToCAD.py`, and three of them silently change the geometry:

1. the definition cache is keyed on the volume **name**, and TRD reuses names for
   volumes with different shapes — **43.5 % of TRD's 637 374 placements would be filled with
   the wrong solid** (§3.2);
2. a reflected placement of a volume **with daughters** is dropped, and the orphaned subtree
   survives as a free root shape at the identity — for ABSO that is the **whole front absorber**,
   27 of its 30 placements at the wrong world transform (§5.3);
3. baking a reflection through `BRepBuilderAPI_GTransform` turns every analytic face into a
   **B-spline** — 26 076 of TRD's 33 017 faces — which costs 8 505 CSG acceptances (§3.5);
4. and that same bake is **not volume-preserving on curved carriers**: every mirrored
   `TGeoTube` in TRD comes back **0.86 % too large**, while mirrored planar solids stay exact to
   1e-12 (§3.8).

None of these is a defect in the source geometry. Unlike PIPE, whose 81 coincident plie placements
were a real defect in the beam pipe, **all three modules are clean on the TGeo side**: zero
coincident (volume, world transform) duplicates in any of them.

---

## 1. The three source geometries

```bash
# the sim env only -- NEVER with the pythonOCC prepends, they segfault o2-sim
export VMCWORKDIR=$B/stage/share
$B/stage/bin/o2-sim -n 0 -g boxgen -m <MOD> --configKeyValues "align-geom.mDetectors=none"
```

The module names are the `PassiveBase` / detector registry names and were used verbatim:
`-m TRD`, `-m MAG`, `-m ABSO`. No variant spelling was needed. `-n 0` writes the geometry without
transporting anything: 4.04 s for TRD, 1.02 s for MAG, 1.56 s for ABSO. Everything below uses the
ideal `o2sim_geometry.root`, not the aligned one.

Each of the three drags in the cave and the barrel, so `cave`, `barrel` and `caveRB24` appear in
all three censuses; `-m TRD` additionally pulls in the space frame (`BSEGMO*`, `B077`, `BTOFS*`),
which is where its reflections live.

### 1.1 Census

**The PIPE census script is wrong on TRD, and that is the first finding.**
`scratchpad/tgeo2step/census.py` keys its reachability walk on `vol.GetName()`. A `TGeoManager`
does not require volume names to be unique and TRD's is not: **199 names cover 16 751 of TRD's
17 678 volumes**. The name-keyed walk therefore stops at the first volume of each name and reports
1 126 "reachable volumes" and 15 121 DAG edges, when the identity-keyed answer is 17 678 and
23 006. Every number in this section comes from an identity-keyed rewrite (`census2.py`).

| | TRD | MAG | ABSO | PIPE (Stream AD) |
| --- | ---: | ---: | ---: | ---: |
| logical volumes in the file | 17 678 | 186 | 34 | 208 |
| reachable from `cave`, by identity | **17 678** | **186** | **33** | 198 |
| distinct *names* among them | 1 126 | 19 | 33 | 198 |
| DAG placement edges (`AddNode`) | 23 006 | 691 | 32 | 333 |
| fully expanded placements | **637 374** | **1 026** | **33** | 1 475 |
| max depth | 13 | 7 | 6 | 6 |

| shape class | TRD | MAG | ABSO |
| --- | ---: | ---: | ---: |
| `TGeoTube` | 16 065 | — | 1 |
| `TGeoBBox` | 1 437 | 1 | 1 |
| `TGeoTrd1` | 48 | — | — |
| `TGeoCompositeShape` | 41 | 1 | 4 |
| `TGeoShapeAssembly` | 40 | 2 | 3 |
| `TGeoTrap` | 36 | — | 2 |
| `TGeoXtru` | 4 | 4 | — |
| `TGeoArb8` | 2 | — | — |
| `TGeoTubeSeg` | 2 | — | — |
| `TGeoPgon` | 2 | **177** | — |
| `TGeoPcon` | 1 | 1 | **13** |
| `TGeoCone` | — | — | 9 |
| **total** | **17 678** | **186** | **33** |

MAG is a stack of 8-edge `TGeoPgon`s — the L3 octagon — plus four `TGeoXtru` doors and plugs.
ABSO is thirteen `TGeoPcon`s and nine `TGeoCone`s, i.e. the r–z profile stack the revolved-profile
detector is meant for. TRD is where the polyhedra are: 48 `Trd1` + 36 `Trap` + 4 `Xtru` +
2 `Arb8` + 2 `Pgon` = 92 volumes that are neither boxes nor solids of revolution.

### 1.2 The coincidence walk

The PIPE 81-plie pattern, checked on all three with `scratchpad/tgeo2step/dupcheck.py`:

| | leaf-solid placements | mother-body placements | distinct (volume, world transform) | coincident duplicates |
| --- | ---: | ---: | ---: | ---: |
| TRD | 426 545 | 210 022 | 636 567 | **0** |
| MAG | 517 | 507 | 1 024 | **0** |
| ABSO | 22 | 8 | 30 | **0** |
| PIPE | 1 152 | 71 | 1 142 | 81 |

All three source geometries are clean. The beam pipe's 81 coincident plies remain the only such
defect found so far. TRD's reverse conversion *does* report coincident placements (§3.6) — but
those are made by the exporter, not present in the source, which is the point of running this walk
on the TGeo side first.

## 2. The instruments

The three from Stream AD §4.3, plus one new one.

**A — the STEP read back, against the source `TGeoShape`.** Uniform points in the source shape's
own bbox and frame, `BRepClass3d_SolidClassifier` on the shape OCCT reads out of our own file
against `TGeoShape::Contains`, points within 1e-6 cm of a surface skipped as unfair
(`TGeoShape::Safety`).

**B — the recognised CSG, against the source `TGeoShape`.** For every part accepted as CSG, the
`TGeoShape` the reverse converter emitted into `shape_*.root`, put back in the part frame with the
converter's own `shapePlacement`, against the original shape.

**C — capacity.** `TGeoShape::Capacity()` of the source against `BRepGProp` on the STEP read-back.

**D — placement fidelity, new (`placecmp.py`).** A and B score a part in its own frame, so a part
that is geometrically perfect but landed at the wrong world transform scores clean. D walks both
trees accumulating world matrices and compares the multiset of (leaf name, world matrix) pairs.
A reflecting TGeo placement is excluded from the comparison and counted separately, because the
exporter bakes the reflection into the solid and the STEP component location is then the identity
by design.

Sampling: 10 000 points per part for MAG (22 parts) and ABSO (30 parts); TRD has 9 608 parts, so a
stratified sample keeping every part of a rare source class and thinning the common ones, at
2 000 points per part under a 1 000 000-point budget. The scoring script is `rtscore.py`, an
adaptation of Stream AD's `rtcompare.py`; both live in this session's scratch and are not
committed.

**A caveat that matters for TRD.** Both the exporter and the scorer look a source volume up by
name. Where a name carries several distinct volumes they both get the same one — the prototype the
exporter converted — so A, B and C are *structurally blind* to the name-collision damage. That
damage is measured separately, in §3.2, by replaying the exporter's own build order.

---

## 3. TRD

`-m TRD`. 17 678 volumes, 637 374 placements, and the module that broke three assumptions.

### 3.1 The write

```
1126 logical volumes, 15521 components
1085 solids, 467 volumes with daughters, 40 pure assemblies, 15521 components, 1 volume declined
capacity check: max relative deviation 1.101e-02, median 0.000e+00
8531 reflecting placements baked as mirrored copies
```

140.35 MB, **35.1 s**, peak RSS **2.10 GB**. The STEP writer, which dies with SIGSEGV on the full
Run 3 model (Stream AD §5), completes here without complaint at 15 521 components — so the
writer's limit is somewhere above TRD and below the full geometry's 74 601.

One volume declined, and it is a genuine mapping gap: `B045cut`, a `TGeoSubtraction` whose right
operand `cutOnB045` is a `TGeoTrd1` with **`dx1 = 0`**. Its −z section is a line segment, so it has
two distinct vertices where the +z section has four, and `_prism_from_rings` refuses:

> `TGeoTrd1: sections have [2, 4] distinct vertices; a prism needs the same count in every section`

ABSO declines one volume for the same underlying reason (§5.1). See §6.1.

The 1.1e-2 capacity outlier is `BSEGMO14`, a `TGeoCompositeShape` — `TGeoCompositeShape::Capacity()`
is itself a Monte Carlo, established in Stream AD §4.3 as the noisy side of that comparison. The
median over all 1 085 solids is exactly 0.

### 3.2 The name-collision defect — 43.5 % of TRD's placements

`O2_TGeoToCAD.py` caches definitions as `self.definitions[str(vol.GetName())]`. TRD has 199
colliding names, and **11 of them carry more than one distinct shape**:

| name | volumes | distinct shapes | placements filled with the prototype instead | median \|ΔV\|/V |
| --- | ---: | ---: | ---: | ---: |
| `UTCP` | 7 752 | 7 | 134 260 | **29.1** |
| `UTCH` | 7 752 | 7 | 134 260 | **29.1** |
| `UTPL` | 456 | 6 | 6 644 | 1.28e-1 |
| `UTP3` | 48 | 3 | 504 | 2.00 |
| `UTC3` | 48 | 2 | 432 | 1.88 |
| `UTC4` | 48 | 2 | 432 | 1.88 |
| `UTG4` | 41 | 6 | 146 | 1.35e-1 |
| `UTG3` | 41 | 6 | 145 | 1.35e-1 |
| `USCR` | 34 | 6 | 123 | 8.68e-1 |
| `UTP1` | 48 | 2 | 36 | 4.00e-2 |
| `UTA2` | 3 | 2 | 36 | 5.66e-1 |

Replaying the exporter's depth-first build order and then walking the true tree
(`collide.py`): of TRD's 637 374 expanded placements, **277 018 — 43.5 % — sit on a volume whose
shape differs from the one the exporter converted for that name**, 142 181 of them leaves. The
worst offenders are the cooling pipes `UTCP`/`UTCH`, where the prototype is on average **29× the
volume** of the shape it stands in for.

This is silent. Nothing in the report, nothing in the log, no capacity check fires — the capacity
check compares the converted solid against *its own* volume's `Capacity()`, and that volume is the
prototype. It is also invisible to instruments A, B and C for the same reason.

The fix is one line of keying: `ROOT.addressof(vol)` instead of `vol.GetName()`, with the XCAF
label name disambiguated. MAG has one colliding name (`L3CD` ×168) but all 168 carry the same
shape, so its damage is 0; ABSO has no collisions at all. PIPE had none, which is why Stream AD
never saw this.

### 3.3 Reflections — 8 531 of them

TRD's space frame mirrors its halves, and it is where nearly all of the full Run 3 geometry's
9 309 reflecting placements come from: **8 531 in this one module**, of which 8 505 are baked into
private mirrored copies of the child solid and **26 are dropped** — the four `BTOFS0*` TOF
supermodule assemblies, which have daughters and so have no single solid to bake.

That is Stream AD's known "98 reflected placements of subtree-carrying volumes are dropped" item,
and TRD shows what it costs: TGeo has 636 567 leaf + mother-body placements, the STEP declares
**499 307**, so **137 260 placements — 21.6 % — are lost** to those 26 dropped edges.

### 3.4 What came back

`O2_CADtoTGeo.py --exact-surfaces auto --csg auto`, no `--mesh`: **3 min 21 s**, peak RSS
**2.14 GB**.

| | |
| --- | ---: |
| leaf solid definitions in the STEP | **9 608** |
| of which baked mirrored copies | 8 505 |
| leaf placements declared | 499 307 |
| leaf placements emitted | 499 293 |
| tiers | **CSG 952 / exact surface 8 656 / tessellated 0** |

Zero tessellated fallbacks: every solid in TRD is represented exactly.

**The faces, as the reverse converter classifies them:**

| plane | bspline | cylinder |
| ---: | ---: | ---: |
| 6 889 | **26 076** | 52 |

### 3.5 The mirror-bake carrier defect

**All 26 076 B-splines are faces of the 8 505 baked mirrored copies**, and the reverse converter's
canonical-form pre-pass recovers every one of them: 17 666 planes and 8 410 cylinders, "recovered
from stored bspline", at exact tolerance. So the exact-surface tier is unharmed —
9 608/9 608 solids represented exactly. The **CSG recogniser is not**, because it reads the stored
carrier type:

| | parts | CSG | surface | tessellated |
| --- | ---: | ---: | ---: | ---: |
| non-mirrored | 1 103 | **952 (86.3 %)** | 151 | 0 |
| baked mirrored | 8 505 | **0** | 8 505 | 0 |
| all | 9 608 | 952 (9.9 %) | 8 656 | 0 |

Every one of the 8 505 declines with the same reason, *free-form faces: N of N (surface kind
outside plane/cylinder/cone/sphere)* — including `UTCP__mirrored`, which is a **cylinder** and
whose un-mirrored twin is accepted as `TGeoTube` at dV_sym = 0.

The cause is `apply_tgeo_matrix`, which sends a reflecting matrix through
`BRepBuilderAPI_GTransform`. A `gp_GTrsf` is an arbitrary affine map, so OCCT degrades every
analytic carrier to a B-spline rather than assume it is preserved. A pure mirror is an isometry
and `gp_Trsf::SetMirror` expresses it exactly, which `BRepBuilderAPI_Transform` would then apply
carrier-preserving. This is the single largest CSG-acceptance lever measured anywhere in this
study: it is worth 8 505 parts on TRD alone.

The second, cheaper half of the fix is on the reverse side: the CSG recogniser should consult the
recognition pre-pass's recovered carrier, which already says "plane" and "cylinder" at
6.8e-14 cm, instead of the stored B-spline.

### 3.6 The coincidence report is the exporter's, not TRD's

The reverse converter warns:

```
WARNING: this CAD model DECLARES 14 leaf solid placement(s) that coincide exactly with another
placement of the same solid.
  Leaf placements: 499307 declared (multiplicity {1: 499287, 2: 2, 4: 4}) -> 499293 distinct.
  Suppressed 8 placement edge(s) ... dropped BTOFS00 -> BTOFS1 [coincident-occurrence] ...
```

On PIPE the same warning was arbitrated onto the TGeo side and turned out to be a real defect in
the beam pipe. Here the TGeo-side walk (§1.2) says **0 coincident duplicates**, and every
suppressed edge is a `BTOFS*` — the same four assemblies whose reflected placements were dropped
in §3.3. The orphaned subtrees land as free root shapes at the identity, on top of each other.
**The coincidence is manufactured by the exporter.** Running the cheap TGeo-side arbiter first is
what makes the two cases distinguishable, and it cost one script.

### 3.7 CSG acceptance, and why not more

952 of 1 103 non-mirrored parts accepted, **all 952 at symmetric difference exactly 0 cm³**, and
the converter's own ROOT-vs-CAD Contains check clean over **3 808 000 points**:

| recognised as | n |
| --- | ---: |
| `TGeoBBox` | 922 |
| `TGeoTube` | 28 |
| `TGeoTubeSeg` | 2 |

The 151 non-mirrored declines, by reason:

| n | reason |
| ---: | --- |
| 80 | *a box face has no opposite partner* |
| 44 | *N planar faces: not a six-plane box* |
| 25 | *the three plane pairs are not mutually perpendicular* |
| 1 | *a planar face is neither a cap nor a wedge of any axis cluster* |
| 1 | *N cap plane(s) perpendicular to the axis, expected 2* |

**The three top classes are one missing recogniser.** They are 149 of the 151, and they are TRD's
92 polyhedra plus the mother bodies over them: a `TGeoTrd1` or `TGeoTrap` has six planar faces of
which the sloped pairs are neither parallel nor perpendicular, so the box recogniser rejects them
three different ways depending on which test fires first. This is the answer to Roadmap (h) for
`Trd1`: **no, and the gap is a sloped-prism recogniser, not three separate gaps.** The last two
declines are the shared `barrel__body` composite and the shared `caveRB24` `TGeoPcon`.

### 3.8 The mirror bake is not volume-preserving on curved carriers

Instrument C separates cleanly by whether a part is a baked mirror:

| | n | median | max |
| --- | ---: | ---: | ---: |
| non-mirrored, non-composite | 227 | **1.8e-13** | 8.99e-13 |
| non-mirrored `TGeoCompositeShape` | 22 | 2.3e-03 | 1.35e-02 |
| mirrored `TGeoBBox` / `Trd1` / `Trap` | 75 | **7.7e-13** | 8.92e-12 |
| mirrored `TGeoTube` | 139 | **8.6154e-03** | 8.85e-03 |

The un-mirrored solids are exact — `UTCP`'s forward record reads `relDev 0.0`, `UTCH` 1.6e-16 —
and so are the *planar* mirrored copies. But every mirrored cylinder comes back **larger**, by the
same 8.6154e-03 relative whatever its size (`UTCP` 6.39 → 6.44505 cm³, `UTPL` 2.84 → 2.86447,
`UTCO` 127.235 → 128.331; the one outlier, `BTSHT2_AM`, is 8.8496e-03 and has a different
`rmin/rmax` ratio). `BRepBuilderAPI_GTransform` re-approximates the cylindrical carrier and the
approximation sits outside the true surface.

So the mirror bake is not merely a carrier annoyance that the reverse recogniser can undo: on
curved faces it is a **0.9 % material-budget error**, present in 8 505 of TRD's parts. The fix in
§6.4 item 2 removes it, because a `gp_Trsf` mirror moves the exact analytic carrier.

The 7.2e-01 outlier the raw table shows for `UTG3__mirrored` is not this effect: it is the name
collision of §3.2 leaking into the instrument. For five of TRD's 1 126 names — `UTG3`, `UTG4`,
`USCR`, `UTC1`, `UTP1` — the exporter's depth-first prototype is a *different volume* from the one
`geo.GetVolume(name)` returns, and for `UTG3` and `UTG4` the two differ in volume by 72 % and
81 %. Instrument C therefore catches the collision by accident, from the other side.

### 3.9 Scoring

Stratified sample, 500 of 9 608 parts × 2 000 points (1 000 000-point budget), keeping every part
of a rare source class and thinning the common ones.

| instrument | result |
| --- | --- |
| A, STEP read-back vs source `TGeoShape` | 249 parts, **497 774 points, 0 disagreements** |
| B, recognised CSG vs source `TGeoShape` | 134 parts, **267 774 points, 0 disagreements** |
| C | §3.8 |
| D, placement fidelity | below |

Instrument D (43 s, 2.12 GB):

```
free root shapes in the STEP: 32                     (should be 1)
TGeo non-reflecting placements: 352 693  ->  352 115 at the right world transform (99.84 %)
   578 missing (UTCO__body x288, UTCL x288, B051__body, B052), 63 extra (all BTOFS*/BTSHT*)
TGeo reflecting placements: 283 874  (44.6 % of TRD lies under a reflection)
```

The 32 free root shapes are the orphaned subtrees of §3.3, and the 63 extra placements are all
inside them. The 578 missing non-reflecting placements are two patterns of 288, `UTCO__body` and
`UTCL`, which is the §3.2 collision again — those names' subtrees are not the prototype's.

The reflecting half is not placement-comparable by construction, since the reflection is baked into
the solid and the component location becomes the identity; the loss on that side is the 137 260
placements counted in §3.3.

## 4. MAG

`-m MAG`. The L3 magnet: 186 volumes, all but nine of them 8-edge `TGeoPgon`s.

### 4.1 The round trip

| | forward | reverse |
| --- | --- | --- |
| wall | **0.51 s** | **1.81 s** |
| peak RSS | 744 MB | 715 MB |
| size | 2.01 MB STEP, 196 components | |
| coverage | 17 solids, 8 with daughters, 2 pure assemblies, **0 declined** | |
| reflections | 5, all baked | |
| capacity check | max 5.87e-3, median **9.3e-16** | |
| tiers | | **CSG 1 / exact surface 21 / tessellated 0** |

The 5.87e-3 is the shared `barrel` composite; the median over all 17 solids is at the
double-precision floor. Every `TGeoPgon` and `TGeoXtru` in the module converts.

### 4.2 Faces, and CSG

| plane | bspline | cylinder |
| ---: | ---: | ---: |
| 410 | **98** | 5 |

The 98 B-splines are, again, exactly the faces of the five baked mirrored copies (`L3CR`,
`L3DoorR`, `L3DoorL`, `L3PlugSPR`, `L3PlugSPL`); the recognition pre-pass recovers all 98 as
planes at max gap **6.8e-14 cm**, so 22/22 solids are still exact surfaces.

CSG: **1 of 22**, and the one is the `cave` box. The 21 declines:

| n | reason |
| ---: | --- |
| 14 | *N planar faces: not a six-plane box* |
| 5 | *free-form faces* (the five mirrored copies) |
| 1 | *a planar face is neither a cap nor a wedge of any axis cluster* (`barrel__body`) |
| 1 | *N cap plane(s) perpendicular to the axis, expected 2* (`caveRB24`, a `TGeoPcon`) |

**This is the direct answer to Roadmap (h) for `TGeoPgon`: no.** All ten of MAG's distinct Pgons
and all four Xtrus land in one class, *"18 / 42 / 50 planar faces: not a six-plane box"* — every
individual face is recognised as a plane, and the whole-part rule has nothing between "six planes
in three perpendicular pairs" and "give up". A prism recogniser (a closed planar section swept
along an axis) takes all 14, and it is the same recogniser TRD's 149 declines want.

### 4.3 Scoring

| instrument | result |
| --- | --- |
| A, STEP read-back vs source `TGeoShape` | 17 parts × 10 000 = **170 000 points, 0 disagreements** |
| B, recognised CSG vs source `TGeoShape` | 1 part × 10 000 = **10 000 points, 0 disagreements** |
| D, placement fidelity | 1 free root shape; **1 019/1 019 non-reflecting placements at the right world transform**, 5 reflecting baked as expected, 0 extra |

C, capacity of the STEP read-back against `TGeoShape::Capacity()`, all 22 parts:

| source class | n | median | max |
| --- | ---: | ---: | ---: |
| `TGeoPgon` | 11 | 2.9e-13 | 2.0e-09 |
| `TGeoXtru` | 8 | 4.1e-13 | 4.3e-13 |
| `TGeoPcon` | 1 | — | 3.0e-14 |
| `TGeoBBox` | 1 | — | 0.0 |
| **all non-composite** | **21** | **2.9e-13** | **2.0e-09** |
| `TGeoCompositeShape` | 1 | — | 5.9e-03 |

The 2.0e-09 is `L3CCO__body`, a thin-walled 4-plane Pgon coil casing — `BRepGProp`'s integration
floor over a shell whose volume is 4e4 cm³ inside a 5e8 cm³ envelope, not a geometry difference.
The five mirrored parts, scored on capacity only (their source shape is not mirrored), sit at
2.1e-15 to 4.1e-13 — but that is because all five are *planar* (`Pgon`, `Xtru`). On a curved
carrier the same bake loses 0.9 % of the volume; TRD §3.8 has the measurement.

MAG is the clean module of the three: no name collisions with divergent shapes, no dropped
placements, no declines on the forward side, nothing tessellated, and perfect placement fidelity.

## 5. ABSO

`-m ABSO`. The front absorber: 33 volumes, 13 `TGeoPcon` and 9 `TGeoCone`.

### 5.1 The round trip

| | forward | reverse |
| --- | --- | --- |
| wall | **0.19 s** | **2.37 s** |
| peak RSS | 732 MB | 714 MB |
| size | 0.61 MB STEP, 38 components | |
| coverage | 29 solids, 11 with daughters, 3 pure assemblies, **1 declined** | |
| reflections | 2, of which **1 baked and 1 dropped** | |
| capacity check | max 5.87e-3, median **1.5e-16** | |
| tiers | | **CSG 9 / exact surface 21 / tessellated 0** |

The declined volume is `AFassCentral`, a subtraction whose left operand `FassCone` is a
`TGeoPgon` with `rmin = 0` on its first two z planes and `rmin = 177` on its last two:

> `TGeoPgon: inner sections do not match the outer count`

The inner ring at `rmin = 0` collapses to a single point, so it has one distinct vertex where the
outer ring has eight. This is the same underlying gap as TRD's `cutOnB045` (§3.1) — a section that
degenerates to a point or a line — reached from the other direction. See §6.1.

### 5.2 Faces, and CSG

| plane | cylinder | cone | bspline |
| ---: | ---: | ---: | ---: |
| 113 | 52 | 64 | **6** |

The six B-splines are the one baked mirrored copy `AFassUMFlange__mirrored`; recovered as six
planes at 6.2e-14 cm.

CSG: **9 of 30**, all at symmetric difference exactly 0 and the converter's own ROOT-vs-CAD
Contains clean over 36 000 points — 7 × `TGeoCone`, 1 × `TGeoTube`, 1 × `TGeoBBox`.

The 21 declines:

| n | reason |
| ---: | --- |
| 12 | *N cap plane(s) perpendicular to the axis, expected 2* |
| 2 | *a planar face is neither a cap nor a wedge of any axis cluster* |
| 2 | *N coaxial cone faces, expected 1 or 2* |
| 2 | *mixed lateral surface kinds `['cone', 'cylinder']` on one axis* |
| 2 | *a box face has no opposite partner* (the two `TGeoTrap` flanges) |
| 1 | *free-form faces* (the mirrored flange) |

**This reproduces PIPE's headline exactly.** PIPE's largest decline class was 46 of 106 "more than
two cap planes perpendicular to the axis" — a `TGeoPcon`. Here it is 12 of 21, plus the two
*"4 coaxial cone faces, expected 1 or 2"* which are the same Pcon shape failing on the cone count
first (`AFaPbCone`, a 3-plane stack with `rmin > 0`), and the two *"mixed lateral surface kinds"*
which are one-band Pcons whose inner wall is a cylinder while the outer is a cone. **Sixteen of
ABSO's 21 declines are one missing revolved-profile detector**, on a module that is nothing but
revolved profiles. Roadmap (h) for `TGeoPcon`: no, still, and with a second exactly-known corpus
behind it — PIPE's 58 Pcons and ABSO's 13.

### 5.3 The dropped reflection — the whole absorber

ABSO records two reflecting placements, one baked and one not. The one that is not is
`barrel/AFA_1`, and `AFA` is the whole front absorber assembly:

```cpp
// Detectors/Passive/src/Absorber.cxx:270, :1039
TGeoRotation* rotxz = new TGeoRotation("rotxz", 90., 0., 90., 90., 180., 0.);   // diag(1, 1, -1)
barrel->AddNode(voFA, 1, new TGeoCombiTrans(0., 30., -90., rotxz));
```

`rotxz` is a pure z-mirror, determinant −1. `AFA` is a `TGeoVolumeAssembly`, so it has no solid of
its own to bake the reflection into, and `build()` records
*"reflected placement of a volume with daughters cannot be baked"* and **skips the component**.
The subtree has already been built, so it survives in the XCAF document unreferenced and is
written as a free root shape at the identity.

Instrument D measures the consequence:

```
free root shapes in the STEP: 3          (should be 1)
reflecting placements in TGeo (baked, no location to compare): 26
placements in TGeo but not in the STEP: 1
placements in the STEP but not in TGeo: 26
```

**27 of ABSO's 30 placements are at the wrong world transform** — only `cave__body`,
`barrel__body` and `caveRB24`, which sit above `AFA`, are right. Every part is geometrically
perfect — instrument A scores 29 parts × 10 000 = 290 000 points with zero disagreements — and the
whole absorber sits 90 cm and one mirror away from where it belongs. This is precisely the class
of defect that A, B and C cannot see, and it is why instrument D exists.

It is also the reason to take Stream AD's item 1b seriously: "98 reflected placements of volumes
with daughters are dropped in the full geometry" reads like a rounding error until one of the 98
turns out to be the entire front absorber.

### 5.4 Scoring

| instrument | result |
| --- | --- |
| A, STEP read-back vs source `TGeoShape` | 29 parts × 10 000 = **290 000 points, 0 disagreements** |
| B, recognised CSG vs source `TGeoShape` | 9 parts × 10 000 = **90 000 points, 0 disagreements** |
| D, placement fidelity | **27 of 30 placements at the wrong world transform**; 3 free root shapes (§5.3) |

C, capacity:

| source class | n | median | max |
| --- | ---: | ---: | ---: |
| `TGeoPcon` | 13 | 3.6e-14 | 4.8e-13 |
| `TGeoCone` | 9 | 6.3e-14 | 5.0e-12 |
| `TGeoTrap` | 3 | 1.2e-13 | 1.5e-13 |
| `TGeoTube` | 1 | — | 5.5e-15 |
| `TGeoBBox` | 1 | — | 0.0 |
| **all non-composite** | **27** | **6.3e-14** | **5.0e-12** |
| `TGeoCompositeShape` | 3 | 2.5e-03 | 5.9e-03 |

Thirteen `TGeoPcon`s at a median of 3.6e-14 relative is the strongest statement in this document
about the revolve construction: the r–z profile revolve reproduces `TGeoPcon::Capacity()`, which is
analytic, to fourteen digits, over stacks of up to 16 z planes.

## 6. What the three modules add up to

### 6.1 One gap behind two forward declines

TRD's `cutOnB045` (`TGeoTrd1` with `dx1 = 0`) and ABSO's `FassCone` (`TGeoPgon` with `rmin`
stepping through 0) are the same defect in `_prism_from_rings`: a section that **degenerates to a
point or a line** has fewer distinct vertices than its partner, and the sewer requires equal
counts. Both are legal TGeo shapes and both are used in production geometry. The fix is to allow a
degenerate ring, sewing triangles to the collapsed vertex or edge instead of quads — the same thing
`BRepPrimAPI_MakeCone` does for a cone that closes to an apex.

### 6.2 The forward converter is right, and its keying is not

Across the three modules the *shape mapping* is not the problem: 1 131 of 1 133 shaped volumes
convert (99.8 %), the capacity medians are 0, 9.3e-16 and 1.5e-16, and 1 130 000 Contains samples
found zero disagreements. Every defect this study found is in the **bookkeeping around** the
shapes — which volume a name means, where a reflected subtree goes, and what carrier a mirrored
face keeps.

### 6.3 Cross-module comparison

| | PIPE (AD) | **TRD** | **MAG** | **ABSO** |
| --- | ---: | ---: | ---: | ---: |
| reachable volumes (identity) | 198 | 17 678 | 186 | 33 |
| expanded placements | 1 475 | 637 374 | 1 026 | 33 |
| coincident duplicates in TGeo | 81 | 0 | 0 | 0 |
| **forward** | | | | |
| solids built / shaped volumes | 170/170 | 1 085/1 086 | 17/17 | 29/30 |
| declined | 0 | 1 | 0 | 1 |
| components | 1 450 | 15 521 | 196 | 38 |
| reflections baked / dropped | 0 / 0 | 8 505 / 26 | 5 / 0 | 1 / 1 |
| STEP size | 4.36 MB | 140.35 MB | 2.01 MB | 0.61 MB |
| wall / peak RSS | 3.5 s | 35.1 s / 2.10 GB | 0.51 s / 744 MB | 0.19 s / 732 MB |
| capacity, median (max non-composite) | 1.5e-13 (1.8e-08) | 0 (—) | 9.3e-16 (2.0e-09) | 1.5e-16 (5.0e-12) |
| **reverse** | | | | |
| leaf solids | 172 | 9 608 | 22 | 30 |
| wall / peak RSS | 40 s | 3 min 21 s / 2.14 GB | 1.8 s / 715 MB | 2.4 s / 714 MB |
| CSG | 62 (36 %) | 952 (9.9 %) — **86.3 % of non-mirrored** | 1 (4.5 %) | 9 (30 %) |
| exact surface | 87 | 8 656 | 21 | 21 |
| tessellated | 23 | **0** | **0** | **0** |
| non-quadric faces | 26 of 1 093 | 26 076 of 33 017 (all mirror-bake) | 98 of 513 (all mirror-bake) | 6 of 235 (all mirror-bake) |
| **scoring** | | | | |
| A, Contains points / disagreements | 1 689 997 / 1 | 497 774 / **0** | 170 000 / **0** | 290 000 / **0** |
| B, CSG-vs-source points / disagreements | 620 000 / 0 | 267 774 / **0** | 10 000 / **0** | 90 000 / **0** |
| D, non-reflecting placements right | not run | 352 115/352 693 | **1 019/1 019** | **3/4** (26 more under the dropped reflection) |
| D, free root shapes (should be 1) | not run | **32** | 1 | **3** |
| top decline reason | Pcon caps, 46/106 | mirror B-splines, 8 505/8 656 | Pgon prisms, 14/21 | Pcon caps, 12/21 |

### 6.4 The single ranked list of what to fix

1. **Key the definition cache on volume identity** (`O2_TGeoToCAD.py`). Silently wrong geometry on
   43.5 % of TRD's placements; zero cost to fix.
2. **Mirror as a `gp_Trsf`, not a `gp_GTrsf`.** Worth 8 505 CSG acceptances on TRD, it removes
   26 180 spurious B-spline faces across the three modules, and it removes a **0.86 % volume
   error** on every mirrored cylinder (§3.8) — the only finding here that changes a material
   budget.
3. **Mirror a whole subtree** instead of dropping the component. Costs ABSO its entire absorber
   placement and TRD 137 260 placements. If it stays unfixed, it must at minimum be *loud*: today
   it is one line in the report's `volumes` array.
4. **A sloped-prism / polygonal-prism recogniser** on the reverse side. 149 of TRD's 151
   non-mirrored declines and 14 of MAG's 21 — one recogniser, two modules.
5. **The revolved-profile detector** (Stream D, `O2RevolvedSolid`). 16 of ABSO's 21 declines, 46 of
   PIPE's 106. Now backed by two independent exactly-known corpora.
6. **Degenerate sections in `_prism_from_rings`.** Two forward declines, two modules, one cause.
7. **Let the CSG recogniser read the recovered carrier**, not the stored one. Cheap, and it makes
   fix 2 unnecessary for the reverse side even if the STEP keeps its B-splines.

### 6.5 Dead premises, recorded

- **"The PIPE census script generalises."** It keys on volume names, and ALICE volume names are
  not unique. On TRD it under-reports reachable volumes by 15×. Anything keyed on a TGeo volume
  name is suspect until the collision count is checked.
- **"The reverse converter's coincidence warning means the source geometry has duplicates."** True
  on PIPE, false on TRD, where the exporter manufactured them by dropping reflected components.
  The cheap TGeo-side arbiter is what tells the two apart, and it must be run first every time.
- **"Contains agreement plus capacity agreement means the round trip is faithful."** It does not.
  ABSO scores 290 000 Contains points and 27 capacities at 6e-14 with the entire absorber in the
  wrong place. A part-frame instrument cannot see a world-frame defect; instrument D exists
  because of this module.
- **"`TGeoPara` will show up in ABSO."** It does not. The full geometry's 13 `TGeoPara` volumes are
  all TPC's (`Detectors/TPC/simulation/src/Detector.cxx`); none of TRD, MAG or ABSO uses one, so
  the last unmapped shape class was not exercised here.
- **"A mirrored solid is the same solid."** MAG's five mirrors agree in volume to 4e-13, and that
  was nearly the conclusion — until TRD, where the mirrors are cylinders and come back 0.86 % too
  large. OCCT's general transform rewrites every plane and cylinder as a B-spline, which is the
  same failure mode as the `ThruSections` lofts Stream AD retired, arriving by a different route;
  on a plane the rewrite is exact, on a cylinder it is not. Checking the volume on a *planar*
  corpus would have licensed the wrong general statement.

## 7. Reproducing it

```bash
# 1. the geometry (sim env; NEVER with the pythonOCC prepends)
export VMCWORKDIR=$B/stage/share
$B/stage/bin/o2-sim -n 0 -g boxgen -m TRD --configKeyValues "align-geom.mDetectors=none"

# 2. census (identity-keyed) and the coincidence walk -- ROOT only, no lock needed
python3 census2.py  o2sim_geometry.root census2.json
python3 dupcheck.py o2sim_geometry.root
python3 collide.py  o2sim_geometry.root        # the name-collision damage

# 3/4. the converter env is the O2 env PLUS the OCC prepends, in one process
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root TRD.step --report TRD_report.json
python3 scripts/geometry/O2_CADtoTGeo.py TRD.step --output-folder rt -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report rt/surface_report.json --csg-report rt/csg_report.json

# 5. score it
python3 rtscore.py  o2sim_geometry.root TRD.step rt --points 2000 --max-total-points 1000000
python3 placecmp.py o2sim_geometry.root TRD.step
```

`--dedup-world` was **not** used on any of the three: all three source geometries have zero
coincident (volume, world transform) pairs, so it would have dropped nothing.

`census2.py`, `collide.py`, `rtscore.py` and `placecmp.py` live in this session's scratch and are
not committed; every number they produced is quoted above with the sample size it came from.

## 8. What is open

1. The seven items of §6.4, in that order.
2. The **exporter's STEP-writer limit** is now bracketed: 15 521 components write fine, 74 601
   segfault. Bisecting it is still Stream AD item 2, and TRD gives it a lower bound.
3. `TGeoPara` remains the only unmapped Run 3 shape class and lives in TPC, not in any module
   studied here.
4. **TRD's 578 missing non-reflecting placements** (`UTCO__body` ×288, `UTCL` ×288) are attributed
   to the name collision but not localised further; they should fall out once §6.4 item 1 lands.
