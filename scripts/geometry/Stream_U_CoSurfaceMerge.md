# Stream U — merging co-surface faces: the exporter really does fragment surfaces, and it buys one solid

Date: 2026-08-03. Branch `swenzel/bvhsurfacesolid`. **Measurement only — no production code was
changed, nothing was built, no sidecar was written.** Prerequisites and companions:
[`Stream_R_CoSurfaceTrims.md`](Stream_R_CoSurfaceTrims.md) §4/§5.3/§9.4, whose residual analysis
this acts on; [`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md), the edge census underneath
it; [`Stream_Q_EllipseTrim.md`](Stream_Q_EllipseTrim.md), the outcome this was tested against.

---

## 0. The answer first

**The hypothesis does not hold. After dissolving every patch seam, exactly 1 of the 16 ALICE3
solids has every remaining boundary edge iso in the recognised (φ, h) domain.** It is
`ST0923290_014`, and it would emit through the shipping `_recognized_quadric_wire_block` with no
change to that function at all. The other 15 do not, and the reason is a single measured fact:

> Of the **443** rejected boundary edges on the 16 solids, **108 (24%) are seams** — shared with a
> face lying on the same analytic surface — and **335 (76%) are true boundaries**, shared with a
> face on a *different* analytic surface. Dissolving a seam cannot touch a true boundary.

The premise of the hypothesis is nevertheless **confirmed, and more strongly than
`Stream_R` §4 put it.** The exporter does fragment surfaces, everywhere:

| corpus | analytic faces | distinct analytic surfaces | faces per surface |
| --- | ---: | ---: | ---: |
| ALICE3 `CAD_noETA.stp`, all 55 leaf solids | 7802 | 4051 | **1.93** |
| ALICE3, the 16 target solids | 711 | 370 | **1.92** |
| `Bagger.step` | 288 | 193 | 1.49 |
| `as1-oc-214.stp` | 53 | 39 | 1.36 |
| ladder fixture `cyl_inter_cyl` | 6 | 2 | **3.00** |

and on the 16 target solids **393 of the 427 faces that reach the wire block (92%)** have a
neighbouring face on their own surface. `Stream_R` §4's "380 of 678 (56%)" reproduces here exactly —
**395 of 711 (56%)** of all analytic faces — and the 92% is the same fact restricted to the faces
that matter.

So the diagnosis was right and the inference from it was wrong. **Fragmentation is nearly universal
and almost all of it is on edges that are already iso.** Seams are 1666 of the 4434 boundary edges
on ALICE3's wire-block faces (37.6%) but only 240 of the 763 rejected ones (31.5%): they are
slightly *under*-represented in the population that blocks emission, not concentrated in it.

**What this closes off.** Together with `Stream_R`'s refutation of the naive conjunction (1 of 15)
and its finding that the bounded extension needs an exact cell enumeration that does not exist,
this was the last cheap route to the 16 solids. There is no longer a one-line vocabulary gap to
find here. `ST0923290_014` is a real 1-solid win available for a bounded price (§7); everything
else on that list needs either the exact cell enumeration (`Stream_R` §9.3) or tessellation.

---

## 1. The claim under test, stated so it can fail

`Stream_R` §5.3 predicted, explicitly as a prediction and not a measurement:

> *"merging adjacent NURBS patches recognised as the same analytic surface into one face before
> trimming would remove the coincident neighbour entirely and, on this evidence, take R2 from 13/15
> to 15/15"*

and §9.4 asked for the census. The stronger reading — the one worth testing, because it costs
nothing if true — is that a merge pre-pass makes the *existing parametric* path work:

> **H** — if the exporter split one analytic surface into several faces, then an edge that is "not
> iso" is often an internal seam, not a boundary. Merge the co-surface faces first and the
> remaining true boundary is iso, so `_recognized_quadric_wire_block` already carries it: no new
> format, no kernel change, no approximation.

**H is falsified by counting.** A solid emits only if *every* face extracts, so H needs every
rejected edge of a solid to be either a seam or already iso. Measured: on 15 of the 16 solids it is
not, and the shortfall is not marginal — `ST0923290_022`…`_026` have **0** of their 18 blocking
edges dissolved by any merge.

---

## 2. The instrument, and how to re-run every number here

One new probe, [`probes/coSurfaceMerge.py`](probes/coSurfaceMerge.py), standalone, no build products
involved. It imports the converter's own `extract_graph` (so the leaf solids are the converter's
`def_shapes`) and `probes/trimEdgeCensus.py`'s `stored_model` / `recognize` / `recognized_model` /
`projector_for` / `wire_block_census`, whose face verdicts `Stream_O` §3 checked against the
shipping `recognize_and_extract_face` over 998 ALICE3 faces with **0 mismatches**.

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2 ; SW=$HOME/alisw/sw/ubuntu2404_aarch64
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
cd $HOME/alisw/O2/scripts/geometry

# the instrument's own control -- must be 5 checks, 0 failures, before anything else is believed
python3 probes/coSurfaceMerge.py --synthetic

# ALICE3, every leaf solid, ~4 min, no --mesh
python3 probes/coSurfaceMerge.py --model ALICE3:ALICE_3_example/CAD_noETA.stp \
    --coincidence-check --negative-control --json /tmp/u/alice3.json

# the controls
python3 probes/coSurfaceMerge.py --model Bagger:STEP_examples/Bagger.step \
    --model as1:STEP_examples/as1-oc-214.stp --fixtures \
    --coincidence-check --negative-control --json /tmp/u/controls.json

# the threshold sweep of §3.3 -- the answer is flat over four orders
for f in 0.5 2 10 100 ; do
  python3 probes/coSurfaceMerge.py --model ALICE3:ALICE_3_example/CAD_noETA.stp \
      --solids ST0923290_002,ST0923290_003,ST0923290_004,ST0923290_005,ST0923290_011,ST0923290_014,ST0923290_016,ST0923290_017,ST0923290_021,ST0923290_022,ST0923290_023,ST0923290_024,ST0923290_025,ST0923290_026,ST0923290_027,ST0923290_028 \
      --merge-factor $f ; done
```

### 2.1 How "the same analytic surface" is decided — geometrically, never by parameter tuples

Comparing two recognised parameter tuples needs a tolerance per parameter and a way to turn an
angle into a length. The probe does none of that. For a candidate pair (A, B) it

1. samples a 5 × 5 grid of points on face B's own **stored** surface over its (u, v) rectangle;
2. projects them onto B's **ideal** recognised surface in closed form (plane, cylinder, cone,
   sphere, torus each have one), so the projected points lie on B *exactly*, to double precision,
   and the patch-versus-ideal recognition residual is removed from the comparison;
3. measures their largest distance to A's implicit surface;
4. symmetrises, and calls the result the **surface separation** — a length in cm, the largest gap
   between the two surfaces over the region either patch occupies.

The threshold is the **model's own declared BRep tolerance** (the larger of the two faces'),
scaled by `--merge-factor`, default 2.0 — the same yardstick and the same default as
`Stream_O` §3's bucket-B threshold.

### 2.2 What "merge" means, and the distinction that matters

The probe keeps **two** groupings and they answer different questions:

- **surfaces** — every face lying on one analytic surface, adjacent or not. This is the
  fragmentation census of §3.
- **merge components** — the connected components of *"shares a boundary edge **and** lies on the
  same surface"*. Only these may be fused. Six windows cut in one tube are one surface but not one
  face, and fusing them would invent a face whose region is disconnected. On ALICE3 this is not a
  hypothetical: **620 surfaces carry more than one disconnected patch set**, `ST2487457_01` has
  four cylinders with **70 faces each**, and `ST0923290_016` faces 30/34/38/42/46/50 are six
  disjoint patches of one cylinder with **zero** shared edges between them.

Two merge **policies** are scored throughout, because merging geometry that already works is a risk
and not a gain:

| policy | what it fuses |
| --- | --- |
| **P-need** | only a component containing a face the shipped wire block declines today |
| **P-all** | every adjacent co-surface component |

---

## 3. Three self-checks, run before anything below was believed

### 3.1 The per-edge iso verdict is the shipped one

The probe builds its own (φ, other) frame per merge group — the reference direction is arbitrary,
and it samples each edge in the 3D curve's own direction rather than the wire's, so that the merged
boundary can be chained into loops through `FirstVertex`/`LastVertex`. Neither can move an iso
verdict (the test reads `max − min` of each run), and that is checked rather than argued: for every
*unmerged* face the probe re-runs `wire_block_census` with the shipped frame and the shipped
traversal and compares edge for edge.

```
   ALICE3   0 mismatches over 4434 edges
   as1      0 mismatches over  112 edges
```

### 3.2 The instrument can fail — a deliberately fragmented cylinder, and a negative control beside it

`--synthetic` builds four quadrant patches of **one** cylinder in-process, NURBS-converts them so
the recogniser rather than the stored type has to decide, and puts a genuinely different cylinder
beside them:

```
  [ok ] the split cylinder really does produce four separately-stored patches -- found 4, stored types ['bspline']
  [ok ] the control cylinder produces a recognisable lateral face -- found 1
  [ok ] the four patches of ONE cylinder are recognised as the same surface
        -- worst separation 3.997e-15 cm vs threshold 2.000e-07 cm
  [ok ] a genuinely DIFFERENT cylinder is not merged with them (negative control)
        -- nearest separation 1.955e+01 cm vs threshold 2.000e-07 cm
  [ok ] the two populations are separated by orders of magnitude, not by the threshold
        -- ratio 4.891e+15

5 checks, 0 failure(s)
```

### 3.3 The threshold is over-determined, and the answer does not depend on it

`--negative-control` reports, per solid, the nearest same-kind pair the probe did **not** merge.

| corpus | merged pairs, max separation | nearest NON-merged pair |
| --- | ---: | ---: |
| ALICE3, all 55 solids | **0.215** declared tolerances | **1.0e+04** declared tolerances |
| ALICE3, the 16 targets | 0.0186 | 1.48e+06 |
| `Bagger.step` | 2.64e-05 | **3.54** ← the tightest in the whole corpus |
| `as1-oc-214.stp` | 0.00115 | 2.5e+07 |
| ladder fixtures | 5.33e-08 | 8.44e+07 |

Bagger's is the one to look at: `Bucket` f#22 and f#23 are two planes **3.545e-08 cm** apart against
a declared tolerance of 1e-08 cm — genuinely distinct, 0.35 nm apart, and correctly kept apart. The
threshold at 2.0 therefore sits, **across the whole corpus**, in an interval bounded below by 0.215
(ALICE3's worst merged pair) and above by 3.54 (Bagger's nearest split pair): a factor of 16 of
empty room globally, and 5 orders on ALICE3's own targets. `Stream_O`'s standard was 1.024 versus
1.39, a factor of 1.36; this is wider.

Swept directly, the decisive number does not move at all:

| `--merge-factor` | 0.5 | 2.0 (default) | 10 | 100 |
| --- | ---: | ---: | ---: | ---: |
| adjacent co-surface components on the 16 targets | 451 | 451 | 451 | 451 |
| rejected edges: seam / true boundary | 108 / 335 | 108 / 335 | 108 / 335 | 108 / 335 |
| **solids that would emit** | **1 / 16** | **1 / 16** | **1 / 16** | **1 / 16** |

### 3.4 The seam verdict, cross-checked against an instrument that shares no code with it

`Stream_R` §4 calls a neighbour **coincident** when its implicit function is identically zero over
the face's own surface samples — a statement about function values on a (u, v) grid. The grouping
above is a statement about closed-form projections between two ideal surfaces. `--coincidence-check`
runs both on every boundary edge and compares:

```
   ALICE3   4116 / 4116 agree,  0 disagree
            max |f_neighbour| on the face's own surface:
              seam edges           max  1.53e-09 cm
              true-boundary edges  min  2.23e-03 cm      -- six orders apart, no overlap
   as1       112 /  112 agree,  0 disagree
              seam max 1.36e-11 cm, true-boundary min 0.281 cm
```

Not one edge is classified by a tolerance decision. The two populations are separated by six orders
of magnitude and the threshold could be moved anywhere between them.

---

## 4. Fragmentation — the premise, and it is confirmed

### 4.1 ALICE3, all 55 leaf solids

```
   7802 analytic faces on 4051 distinct surfaces          1.93 faces per surface
   faces per surface:  1:2444  2:970  3:292  4:207  5:17  6:68  7:4  8:15  11:1  12:3
                       14:1  16:1  18:4  20:7  23:2  24:2  36:1  46:2  47:6  70:4
   ... of which ADJACENT co-surface components:  6088
                       1:4943  2:731  3:268  4:142  5:1  6:1  7:2
   surfaces carrying more than one disconnected patch set:  620
```

The worst cases, named: **`ST2487457_01`, four cylinders of 70 faces each** (f#3/9/15/21/27/33 …,
f#410/634/635/…, f#412/708/…, f#413/424/…); **`ST2487458_01`, six planes of 47 faces each**
(f#465/473/487/499/511/523 …). Neither is on the target list — and the gap between the two
histograms is the point: 70 faces on one cylinder collapse to components of size ≤ 7 once adjacency
is required, because they are separate windows in one tube, not four quadrants of one patch.

### 4.2 The 16 target solids

```
   711 analytic faces on 370 distinct surfaces            1.92 faces per surface
   faces per surface:              1:179  2:66  3:115  4:4  5:2  6:3  11:1
   adjacent co-surface components: 1:316  2:14  3:118  4:2  5:1
   395 of 711 analytic faces (56%) sit in a component of size >= 2   <-- Stream_R §4's 380/678
   393 of 427 wire-block faces (92%) have a co-surface neighbour
```

The dominant motif is unmistakable and it is exactly what `Stream_R` §4 described: **115 surfaces
carrying exactly three faces, and 118 adjacent components of exactly three** — one 120° third each.
(118 > 115 because three of the larger surfaces break into three-patch components.)
`ST0923290_014` f#105…f#110+ is eleven coplanar faces on one plane; `ST0923290_016`
f#16/17/18/22/25 is five patches of one cone.

---

## 5. Seam versus true boundary — the measurement that decides

### 5.1 The split

| population | edges | seam | true boundary |
| --- | ---: | ---: | ---: |
| ALICE3, **rejected** boundary edges | **763** | **240 (31.5%)** | **523 (68.5%)** |
| ALICE3, accepted (control) | 3671 | 1426 (38.8%) | 2245 (61.2%) |
| the 16 target solids, **rejected** | **443** | **108 (24.4%)** | **335 (75.6%)** |
| the 16 target solids, accepted (control) | 1521 | 644 | 877 |

The control row is what refutes the hypothesis rather than merely failing to support it. Seams are
**37.6%** of all 4434 boundary edges on ALICE3's wire-block faces and only **31.5%** of the rejected
ones. A seam is no more likely to be non-iso than any other edge — most seams are constant-φ
generators, which the shipped test accepts. The fragmentation is real, it is everywhere, and it is
overwhelmingly on edges that already work.

### 5.2 What the 335 true boundaries actually are

| neighbour's surface | edges |
| --- | ---: |
| cylinder | 168 |
| plane | 154 |
| cone | 13 |

These are `Stream_O`'s bucket B — genuine quadric∩quadric intersection curves, `plane ∩ sphere` and
`cylinder ∩ cylinder` at the head of the list. They are not artefacts of the exporter's patch
layout and no amount of merging removes them.

### 5.3 Per solid, and what blocks each one

`nonIso` is the number of boundary edges still failing the shipped iso test after every seam in
every component that needed merging has been dissolved.

| solid | faces | rejected | seam | true | merged groups (need/all) | seams dissolved | nonIso after | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `ST0923290_002` | 24 | 18 | 4 | 14 | 1 / 3 | 3 | 14 | blocked |
| `ST0923290_003` | 24 | 16 | 4 | 12 | 1 / 3 | 3 | 12 | blocked |
| `ST0923290_004` | 24 | 18 | 4 | 14 | 1 / 3 | 3 | 14 | blocked |
| `ST0923290_005` | 24 | 16 | 4 | 12 | 1 / 3 | 3 | 12 | blocked |
| `ST0923290_011` | 193 | 200 | 52 | 148 | 20 / 51 | 56 | 148 | blocked |
| **`ST0923290_014`** | 138 | **12** | **12** | **0** | **2 / 35** | **6** | **0** | **EMITS** |
| `ST0923290_016` | 59 | 6 | 0 | 6 | 1 / 7 | 5 | 6 | blocked |
| `ST0923290_017` | 24 | 15 | 0 | 15 | 1 / 3 | 3 | 15 | blocked |
| `ST0923290_021` | 33 | 0 | 0 | 0 | 0 / 5 | 0 | 1 | blocked, see §5.5 |
| `ST0923290_022` … `_026` | 24 | 22 | 4 | 18 | 1 / 3 | 3 | 18 | blocked |
| `ST0923290_027`, `_028` | 24 | 16 | 4 | 12 | 1 / 3 | 3 | 12 | blocked |
| **total** | | **443** | **108** | **335** | | **103** | | **1 / 16** |

The blocking merged group is one per solid on all but `ST0923290_011`, and it is always the same
object — the spherical or conical end cap:

```
ST0923290_002/003/004/005/017/022..026/027/028   f#8,9,10   sphere, 3 patches
     21 edges -> 3 seams dissolved -> 18 boundary, 12 to 18 of them still non-iso
ST0923290_016                                    f#16,17,18,22,25  cone, 5 patches
     23 edges -> 5 seams dissolved -> 18 boundary,  6 still non-iso
ST0923290_011   twelve groups, e.g. f#72,73,74 (cylinder, 23 boundary, 14 non-iso),
     f#131,132,133,134 (cylinder, 23 boundary, 23 non-iso), f#107,109,111,142 (cone, 18/13)
```

The sphere caps are `Stream_O` §8.1's and `Stream_R` §5.2's faces: a spherical patch bounded by
circles in five or six *different* planes, which no single polar axis makes iso and which merging
does not touch, because none of those circles is a seam.

### 5.4 The one that works, in full

`ST0923290_014`'s twelve rejected edge records are **six distinct edges seen from both sides**, and
every one is a seam between two patches of the same cone:

```
   f#116 <-> f#117   max |f_neighbour| on the face  2.7e-12 cm
   f#116 <-> f#118                                  1.0e-13 cm
   f#117 <-> f#118                                  2.8e-12 cm
   f#119 <-> f#120 <-> f#121   likewise, 8.8e-14 .. 7.8e-12 cm
```

Two countersinks, each written as three 120° cone patches, each patch's two radial seams slanted in
(φ, h) because a cone's generator through a patch corner is not a constant-φ line in the frame the
recogniser picks for *that patch*. Merging each triple leaves **one** boundary loop of three rim
arcs, all iso in `h`, winding a full 2π, plus the degenerate apex — three edges, one outer wire,
sweep exactly 2π. `_recognized_quadric_wire_block` accepts that unchanged.

### 5.5 `ST0923290_021` is worse than `Stream_O` §7 recorded, and it is worth its own line

`Stream_O` §7 says `_021` is declined by `recognized quadric trim wire has fewer than 3 edges` and
that this "is not a geometry problem … an implementation limit worth its own line in the ledger,
not a representation gap". Measured directly on the face in question:

```
   ST0923290_021 f#20, recognised as a cone, residual 3.4e-13
     wire 0, outer, 2 edges:
       bspline   spanPhi 1.534 rad   spanOther 2.73e-11 cm   length 1.150 cm   <- iso rim
       bspline   spanPhi 1.534 rad   spanOther 3.436e-01 cm  length 1.283 cm   <- NOT iso
```

The two-edge wire also contains a **genuinely non-iso** boundary edge. Relaxing the three-edge guard
alone would not make `_021` emit; it would move it from one rejection to another. The census could
not see this because `wire_block_census` returns before testing any edge when a wire is too short.

---

## 6. Merged-domain sanity — and the cost of merging what already works

For each merged face the probe chains its surviving boundary edges into loops through shared
vertices, unwraps φ continuously around each loop, and asks whether the merged domain is a plain
rectangle in (φ, other), which is what the sidecar's `phiStart / phiSweep / hLo / hHi` expresses.

**ALICE3, P-all, every merged group of ≥ 2 faces (293 groups):**

```
   rectangle                     215
   wraps a full 2pi              245   of which a FULL BAND (two disjoint rim loops)  205
   with holes (inner wires)       11
   boundary not a clean loop set   0   -- every merged boundary chained into closed loops
```

**This is where merging a working face costs something.** Fusing the three 120° patches of a *full*
cylinder produces a band whose boundary is two disjoint rim circles — two "outer" wires — and
`_recognized_quadric_wire_block` rejects `n_outer != 1`. Measured on ALICE3's 20 solids that emit
today:

| policy | solids that would stop emitting |
| --- | ---: |
| **P-need** | **0** — nothing in an already-emitting solid is touched, by construction |
| P-all | **4**: `ST0923290_013`, `_018`, `_019`, `_020` |

All four are recoverable in principle without any new format, because a rectangular merged domain
needs **no trim wire at all** — it is the parameter-only quadric record `extract_cylindrical_face`
has always emitted for a natively-analytic cylinder, and `recognize_and_extract_face` simply never
takes that branch. But that is a second change to justify for zero coverage gain, and P-need makes
it unnecessary. **The measured recommendation is P-need: merge only where the shipped path
declines.**

Restricted to the 35 groups P-need would actually form on the 16 targets: 21 wrap 2π, 7 are full
bands, 11 have holes, 0 have a boundary that is not a clean loop set, and **10 of 35 come out with
every boundary edge iso** — of which the two on `ST0923290_014` are the only ones whose solid then
emits.

### 6.1 Edge identity would survive, and that is measured, not argued

`applyEdgeIdentityClosure` (`BoundedSurface.h:5591`) demands every model edge id appear exactly
twice with opposite sense. A dissolved seam must therefore have been used exactly twice inside its
group — once by each side — so that removing both occurrences leaves every surviving id still at
two. Over ALICE3's **833 dissolved seams under P-all**:

```
   dissolved seams NOT used exactly twice inside their group        0
   surviving boundary edges used more than once inside their group  0
```

Closure stays balanced by construction, with no tolerance and no sampling. (`as1-oc-214`: 28
dissolved, both counts 0.) The count that would break is exactly the one that is zero, which is why
it is reported rather than asserted.

---

## 7. What is actually on the table, and its price

`ST0923290_014` is a real, bounded, exact win: **ALICE3 20 → 21 sidecars.** What it would cost:

1. A per-solid pre-pass in `extract_surfaces_for_shape`, since `recognize_and_extract_face` is
   per-face today and merging needs solid-level context.
2. Merging entirely in the recognised (φ, h) domain — collect the component's boundary edges,
   discard the seams, chain the rest. It does **not** need `ShapeUpgrade_UnifySameDomain`, which
   would not fire here anyway: the patches' *stored* surfaces are distinct B-splines and only the
   recogniser knows they are one cylinder.
3. One surface record for the component instead of one per face, carrying the union of `edge_refs`
   minus the dissolved seams (§6.1).
4. Gated on P-need, so no shipping geometry moves.

Whether one ALICE3 solid is worth that is a judgement, not a measurement, and it is not mine to
make. What is a measurement: it is one solid, not sixteen, and `ST0923290_014` has never emitted, so
a new sidecar is a new chance to hit `Stream_L_ALICE3Defect.md` mechanism 3 — ALICE3 already has 2
emitted sidecars that do not load.

---

## 8. Being adversarial about this result

**"You measured the wrong merge."** The probe fuses only components that are *edge-connected* and
co-surface. Fusing all 70 faces of one `ST2487457_01` cylinder would be a larger merge — and an
inadmissible one, because the region would be disconnected and the resulting face would have six
outer loops. The looser merge is measured too, as P-all, and it is worse: it breaks 4 emitting
solids and rescues no additional one.

**"Your seam test is too strict."** It cannot be: the sweep in §3.3 runs it 50× looser and 4×
tighter with byte-identical results, and the independent coincidence instrument agrees on 4116 of
4116 edges with six orders of separation between the populations.

**"`ST0923290_016` shows 0 seam edges among its 6 rejected, yet 5 seams dissolved."** Correct and
consistent: its five-patch cone f#16/17/18/22/25 has seams, but those seams are already *iso* and
therefore not in the rejected population. The six rejected edges are its true boundary. This is the
whole refutation in one solid.

**The 92% is not a bigger version of `Stream_R`'s 56%.** They are different denominators of the same
fact and both reproduce here: 395 of 711 analytic faces (56%), 393 of 427 wire-block faces (92%).

**`Stream_R` §5.3's mechanistic prediction survives; its inference does not.** All seven of its
residual R2 failures — `ST0923290_011` f#98/f#134/f#109/f#112/f#111 and `ST0923290_016` f#16/f#18 —
are in merge components of size ≥ 2, so merging *does* remove their coincident neighbour, and the
coincident-neighbour count for a merged face is zero by construction (a co-surface adjacent
neighbour is necessarily inside the component). Nothing here scores R2 itself; whether that takes
R2 from 13/15 to 15/15 remains untested. What is now measured is that it does **not** take the
parametric path from 1/16 to anything else.

**A face-count reduction is a side effect, not a result.** P-need would replace **106 faces by 35
records** on the 16 targets; P-all replaces **855 wire-block faces by 293 records** across ALICE3.
The sidecar-size consequence of that was not measured, and it is worth nothing on its own — 15 of
those solids still emit no sidecar at all.

---

## 9. What this licenses and does not license

### It licenses

1. **Closing the "merge the patches first" route as a coverage lever.** 1 of 16, measured two ways,
   with the threshold swept over four orders and an independent cross-check at 4116/4116.
2. **Recording `ST0923290_014` as an available exact 1-solid gain** with the scope in §7, if anyone
   wants it.
3. **`ST0923290_021` needing more than the three-edge guard** (§5.5). Its two-edge wire also carries
   a non-iso edge; `Stream_O` §7's "implementation limit, not a representation gap" is half right.
4. **P-need over P-all**, if a merge is ever built for any reason: P-all breaks 4 of ALICE3's 20
   emitting solids (§6).
5. **"These solids ship tessellated"** as the honest, evidence-backed answer for the 15, unless and
   until the exact cell enumeration of `Stream_R` §9.3 exists.

### It does not license

1. **Any statement about the kernel, transport or navigation.** No sidecar was written, no gate was
   run, no ray was traced. Nothing here has been through `O2BVHSurfaceSolid`.
2. **Any statement about R2 / implicit trims.** Not scored here. `Stream_R` §5.5's caveat — a
   sampled cell set is inadmissible — is untouched by this work.
3. **A coverage number for any other model.** `Stream_O` §7's lesson stands: coverage is a property
   of the decomposition. Every figure here names ALICE3's 55-solid decomposition.
4. **Anything about merging *planar* faces.** ALICE3's planes reach `extract_planar_face`, not the
   wire block, and its 11-face coplanar group on `ST0923290_014` and Bagger's 7-face coplanar group
   on `Bucket` were censused but their merge was not simulated (§10).
5. **A claim that the merge is free.** §6 measures its cost on working geometry, and it is not zero.

---

## 10. What could not be measured

- **The planar merge path.** The merge simulation covers only faces reaching
  `_recognized_quadric_wire_block`. Bagger's and ALICE3's coplanar groups are counted in §4 but
  their merged domains are not analysed, because `extract_planar_face` expresses a trim differently.
- **Anything downstream of the converter** — see §9.
- **Whether merging changes the recognition residual.** The probe uses one member's recognition for
  the whole component; refitting the merged sample set could give a different (probably better)
  quadric, and that was not tried.
- **Torus merges.** `torus_union_cyl` f#0/f#2 group as one torus at 5.33e-08 declared tolerances,
  which exercises the torus branch of the identity test and nothing else: no torus face is on the
  16 target solids, and no torus face reaches the wire block anywhere in the corpus.
- **IRIS and `oTOF System V3-R92cm.step`.** IRIS would gain nothing in any case: `Stream_O` §7 shows
  its decomposition puts every recognised face inside one 1734-face solid.

---

## 11. Claims in other documents that this work bears on

I have edited none of these; they are for central reconciliation.

- **[`Stream_R_CoSurfaceTrims.md`](Stream_R_CoSurfaceTrims.md) §5.3's prediction** and **§9.4's
  "patch merging as a prerequisite worth measuring"** — measured. The census it asked for: 56% of
  analytic faces and 92% of wire-block faces on the 16 solids are in a co-surface component; the
  coincident-neighbour count does go to zero for a merged face; and it is worth **1 of 16** solids
  on the parametric path. §5.3's "on this evidence, take R2 from 13/15 to 15/15" is neither
  confirmed nor refuted here — R2 was not scored.
- **`Stream_R` §4's "380 of the 678 faces (56%)"** — reproduces exactly, as 395 of 711 (56%) on the
  16 solids including `ST0923290_021`.
- **[`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md) §7's `ST0923290_021` row** — "it is not
  a geometry problem: one recognised cone's trim wire has fewer than three edges … an implementation
  limit worth its own line in the ledger, not a representation gap". Half right: that wire's other
  edge is genuinely non-iso (§5.5), so `_021` needs both.
- **[`NEXT.md`](NEXT.md) item 3** — already amended by `Stream_R`. It should now also say that the
  patch-merge route was measured and is worth **1 of 16**, that `ST0923290_014` is the one, and that
  the gating work item for the rest is unchanged: the exact cell enumeration of `Stream_R` §9.3.
- **[`Tutorial.md`](Tutorial.md) §5.1's list of what is missing for completeness** — `Stream_R` §11
  already records that it omits the trim vocabulary as a coverage lever. Nothing here changes that;
  it removes one candidate from the list of cheap levers.
