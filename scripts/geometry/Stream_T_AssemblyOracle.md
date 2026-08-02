# Stream T — the assembly oracle, and whether this geometry is legal at all

Date: 2026-08-02. Companion to [`Stream_J_XRay.md`](Stream_J_XRay.md) (which scoped this work in
§9 and called the oracle side "the small half") and to [`Tutorial.md`](Tutorial.md) (the map).
Written by the session that built it; the integrating session folds it.

**Read §0 first, because it changes what the next item on the board should be.**

---

## 0. The verdict

**Neither `Bagger.step` nor ALICE3's `CAD_noETA.stp` is a legal TGeo/Geant4 geometry as it stands.**
Both contain pairs of placed solids that genuinely interpenetrate — not touch, *interpenetrate*,
with positive shared volume. Overlapping volumes are forbidden by both engines and produce silently
wrong transport: a track in an undefined volume, energy deposited in the wrong material, and a
navigator whose behaviour is not specified.

**Bagger — 3 of 78 pairs:**

| pair | shared volume | fraction of smaller part | shared slab thickness | sampled depth |
| --- | ---: | ---: | ---: | ---: |
| `Base` ∩ `BoomCylinderOuter` | 8.664e-02 cm³ | 1.80e-03 | 5.41e-02 cm | 2.25e-02 cm |
| `Stick` ∩ `Bucket` | 7.274e-02 cm³ | 1.25e-03 | 2.19e-01 cm | 4.92e-02 cm |
| `BucketLink2` ∩ `BucketLink1` | 3.563e-02 cm³ | 2.90e-03 | 8.20e-01 cm | 3.30e-03 cm |

**ALICE3, `ST0923290` sub-assembly alone — 7 of 378 pairs, and much worse:**

| pair | shared volume | fraction of smaller part | sampled depth |
| --- | ---: | ---: | ---: |
| `ST0923290_013` ∩ `ST0923290_{006,007,008,009}` (4 pairs) | 2.979e-01 cm³ each | **5.74e-02** | **2.29e-01 cm** |
| `ST0923290_019` ∩ `ST0923290_{022,023,024}` (3 pairs) | 3.393e-02 cm³ each | **3.18e-01** | 1.15e-01 cm |

**2.3 mm deep, and one part has 32 % of its volume inside a neighbour.** That is not a rounding
slip. And these are the *only 28 parts I could census in the time available* — the full 206-instance
run is §3.4.

**The other assembly in the corpus is clean, and that is a measurement:** `as1-oc-214.stp`, 18
placed solids, 153 pairs, **0 interpenetrating**, 32 pairs touching with zero shared volume. So the instrument
does not report overlaps everywhere, and Bagger's 3 and ALICE3's 7 can be believed. §3.3 shows the
same clean model turning dirty the moment a deliberate 0.2 cm interpenetration is injected — and
the instrument recovering the injected displacement to six digits.

### And the tool that should have caught this is broken on our own representation

`TGeoManager::CheckOverlaps` is what a Geant4 user reaches for, and on `O2BVHSurfaceSolid` it does
not work:

* at **ROOT's default settings** it prints `GetPointsOnSegments: You should require at least N
  points` for every surface solid in the world — the sampler is starved (default 1000 mesh points
  against 1694–6172 display vertices), and `TGeoShape::GetPointsOnSegments` **returns without
  filling its buffer**. Whatever verdict comes out is not computed from points on the shapes.
* raise `SetNmeshPoints` above the largest mesh and the errors stop and **the verdict changes** —
  and two of its five entries contradict an exact distance computation. It reports
  `BucketLink2 | BucketCylinderInner` as a **0.41 cm overlap** when `BRepExtrema_DistShapeShape`
  says those two solids are **0.09 cm apart**, and `BasePin | Base` as a **0.87–0.91 cm** overlap
  when their boolean intersection has **zero volume**.

Details and the full tables in §3.5. **This is a second, independent blocking item**: even once the
CAD is fixed, there is currently no working overlap check on the representation this branch ships.

---

## 1. What was missing, and what these two scripts are

Everything this branch measures is about **one leaf solid at a time**. `Stream_J_XRay.md` §9 named
the gap precisely: the X-ray benchmark's world holds one `TGeoVolume`, so its crossing list is a
two-state (in-part / out-of-part) sequence, and **leaking** — a track exiting A and the navigator
not reporting it entering B — cannot be expressed by any counter it has.

| file | what it is |
| --- | --- |
| [`assemblyOracle.py`](assemblyOracle.py) | ground truth for a *compound*: per ray, the ordered crossing list annotated with **which volume the track is in after each crossing** |
| [`overlapCensus.py`](overlapCensus.py) | whether the placed parts compose into a world TGeo and Geant4 will accept, pair by pair, **by name** |

```bash
cd $HOME/alisw/O2
python3 scripts/geometry/assemblyOracle.py --self-test      # 26 checks, no model, no build dir
python3 scripts/geometry/overlapCensus.py  --self-test      # 14 checks, no model, no build dir

python3 scripts/geometry/overlapCensus.py --step scripts/geometry/STEP_examples/Bagger.step
python3 scripts/geometry/overlapCensus.py --step .../CAD_noETA.stp --parts ST0923290_01,...
python3 scripts/geometry/assemblyOracle.py --step scripts/geometry/STEP_examples/Bagger.step \
        --parts Stick,Bucket --beams 64 --raster 24
```

Both re-exec themselves under the pythonOCC interpreter via `csg/occ_env.py`, so neither needs the
O2 environment and neither can be caught by the `PYTHONPATH`-replacement trap `NEXT.md` records.
Neither touches C++ or the build directory.

---

## 2. The oracle: occupancy, not a boolean

`xrayOracle.py` classifies the midpoint of every interval between candidate crossings and keeps a
**boolean** — inside the part, or not. `assemblyOracle.py` keeps a **set of occupants**, obtained by
running one `BRepClass3d_SolidClassifier` per part at the same midpoint. That single generalisation
is the whole design, and every assembly-level situation falls out of it:

| situation | occupancy sequence | what it means for a navigator |
| --- | --- | --- |
| touching parts | `{A} → {B}` at **one** distance | report the transition; never emit a vacuum step |
| a genuine gap | `{A} → {} → {B}` | two boundaries, with the vacuum run's length stated |
| a part nested in another | `{A} → {A,B} → {A}` | legal only if declared mother/daughter |
| interpenetration | `{A} → {A,B} → {B}` | **nothing correct exists**; the oracle refuses to pick |
| a ray starting inside | segment 0's occupancy is non-empty | reported as `s0` |

The emitted JSON, one entry per crossing, from a real Bagger ray:

```
t=655.341823  part=Stick        sense=+1  occupancy-after=['Stick']        group=0
t=669.581836  part=Stick        sense=-1  occupancy-after=[]               group=1
t=691.489256  part=BucketLink1  sense=+1  occupancy-after=['BucketLink1']  group=2
...
seg = [[0, 655.34, []], [655.34, 669.58, ['Stick']], [669.58, 691.49, []], ...]
```

**Distances are in the model's native units (mm here), not cm**; `scaleToCm` travels in the
document. That is deliberate — a `TopLoc_Location` must be rigid, so placing a part by its XCAF
location is free while *scaling* it would mean rebuilding every B-rep, which is what makes ALICE3's
206 placements of 55 prototypes affordable at all.

Each crossing carries `{t, part, sense, occupancy-after, group}`. The **occupancy-after** field is
the one that makes leaking detectable: a later C++ round can step a `TGeoNavigator` through the same
assembly and compare not just *where* the boundaries are but *which volume the track is in between
them*. Crossings that share a distance share a `group`, which is how a touching transition is
expressed as one event rather than two with an undefined gap between them.

### The one design decision that is load-bearing

**Candidate positions are merged ACROSS parts before the intervals are cut.** A's exit face and B's
entry face at a shared boundary land at the same ray parameter to within round-off; merging them
makes the answer one transition `{A} → {B}`. Not merging them would put a vacuum sliver of
floating-point width into the ground truth, and *every* navigator in the world would then look wrong
against it. The merge tolerance travels with the answer.

That cuts both ways, so the self-test pins both directions: a **1e-6 cm gap is resolved** as a real
vacuum run and counted as `thin`, and — the control that makes the first statement mean anything —
**at a merge tolerance of 1e-4 the same gap merges away and the pair becomes a touching
transition**. Without the second check, "resolved" would be a property of the assertion rather than
of the instrument.

A vacuum run shorter than `--thin-vacuum` is counted, named and flagged rather than smoothed away,
because a sliver that thin is itself a modelling defect and not something to reproduce faithfully.

### The synthetic assembly, and its results

`assemblyOracle.py --self-test` builds an eight-part assembly whose every answer is arithmetic — no
CAD file, no model, no build directory:

```
x:  0     2 2      4  5      7  7+1e-6   9
    |--A--|--B-----|  |---C--|  |----D---|
     touching face     1 cm gap  1e-6 cm gap

E = [12,18]^3 with F = [14,16]^3 nested wholly inside it
G = [20,24] and H = [23,27] interpenetrating over [23,24]
```

| case | what is checked | result |
| --- | --- | --- |
| **touching** | 4 crossings at 1, 3, 3, 5; occupancy `A, B, B, vacuum`; one contact transition; **no vacuum segment between A and B** | **pass** (5 checks) |
| **gap** | a vacuum run of exactly 1 cm; occupancy after exiting B is vacuum | **pass** (2) |
| **thin gap** | the 1e-6 cm vacuum is resolved, counted `thin`, and D is entered not skipped | **pass** (3) |
| **thin-gap CONTROL** | at merge tolerance 1e-4 the same gap **merges away** and C\|D becomes a touching transition | **pass** (1) |
| **nesting** | occupancy `{}, {E}, {E,F}, {E}, {}`; entering F does **not** exit E; crossings at 2, 4, 6, 8 | **pass** (4) |
| **interpenetration** | occupancy `{}, {G}, {G,H}, {H}, {}`; shared slab exactly `[23,24]`; reported ambiguous rather than resolved to one occupant | **pass** (3) |
| **ray starting inside** | `s0 = {A}`; first crossing is exit A at t = 1 | **pass** (2) |
| **grazing a shared face / edge** | no interior crossing invented; the list still alternates per part | **pass** (2) |
| **Fibonacci fan, 1152 rays** | every part's crossings alternate enter/exit; every crossing's `occ` equals its segment's occupancy | **pass** (2) |
| **negative + positive control** | a touching-only pair reports **no** multiple occupancy; nudging one box 0.1 cm into the other **does** fire the flag | **pass** (2) |

**26 checks, 0 failures.** They passed on the first run, which is precisely when to distrust them —
hence the two mutation controls (the coarse merge tolerance and the 0.1 cm nudge), which are the
only reason the other 24 mean anything.

**What OCCT does with a grazing ray, measured rather than assumed.** A ray running exactly *along*
the shared face `x = 2`, and one along the shared *edge* `x = 2, z = 2`, both give `amb = True` with
empty occupancy and **zero invented crossings**. That is the honest answer — OCCT classifies `ON`
and has no opinion — and a consumer must exclude those rays rather than score them, exactly as
`Stream_J_XRay.md` §1 excludes `amb` rays and the sample gate excludes `nNoVerdict` points.

### It finds the real overlaps independently of the boolean

This is the cross-validation that matters: the census (§3) uses `BRepAlgoAPI_Common`; the oracle uses
ray classification. They share no code path beyond loading the assembly.

**Bagger**, 64 Fibonacci directions × 24 × 24 = 36864 rays, restricted to each pair:

| pair | rays with **ambiguous occupancy** | touching transitions | crossings |
| --- | ---: | ---: | ---: |
| `Stick`, `Bucket` | **16** / 36864 | 39 | 9140 |
| `Base`, `BoomCylinderOuter` | **54** / 36864 | 172 | 14054 |
| `BucketLink2`, `BucketLink1` | **390** / 36864 | 417 | 23656 |

**ALICE3**, 48 directions × 16 × 16 = 12288 rays:

| pair | rays with **ambiguous occupancy** | crossings |
| --- | ---: | ---: |
| `ST0923290_013`, `ST0923290_006` | **10** / 12288 | 2478 |
| `ST0923290_019`, `ST0923290_022` | **52** / 12288 | 11840 |

**Every interpenetration the boolean found is independently visible as a ray segment with two
occupants.** Two instruments, two algorithms, the same pairs.

### A warning about ray density — the same trap as `Stream_J_XRay.md` §6.2

The *first* whole-Bagger run — 24 Fibonacci directions × 6 × 6 = 864 rays over the assembly's whole
bounding sphere — reported **0 rays with ambiguous occupancy**. The overlaps are there; the raster
could not see them. At Bagger's ~100 cm scale a 6 × 6 impact grid has ~30 cm cells and the thickest
of the three defects is 0.22 cm.

**A clean assembly-oracle run at low density is not evidence of a legal geometry.** For the legality
question the census is what should be trusted; the oracle's value is that it answers in the
transport's own terms and can be pointed at a named pair. This is the third time in this project's
history that a whole-model raster has reported a known real defect as clean.

---

## 3. The census: is this geometry legal?

### 3.1 Method, and what each stage costs

1. **AABB rejection.** Free, and what makes N² affordable. Boxes are inflated by `--pad` (default
   0.1 cm), which does **not** change which pairs are found to overlap — an overlapping pair's boxes
   overlap at pad 0 — but decides **which disjoint pairs get their separation measured**. It is a
   stated scope ("every pair within 1 mm"), not a tolerance.
2. **`BRepExtrema_DistShapeShape`** on the survivors. A *positive* distance settles a pair with no
   boolean at all: disjoint, by this much. That number is worth as much as the overlaps —
   **a 1e-9 cm gap and a 1e-9 cm overlap are the same modelling intent and completely different for
   a navigator**, and only a measured separation tells them apart.
3. **`BRepAlgoAPI_Common`** only where the distance is zero. Its **volume** is the discriminator:
   ≈ 0 → **coincident faces** (touching; the *normal* case for an assembly, and legal); > 0 →
   **real interpenetration** (illegal); ≈ the smaller part's volume → **containment** (legal only as
   mother/daughter, which a flat CAD→TGeo conversion does not produce).
4. **Penetration depth**, for interpenetrating pairs only: the max over sampled interior points of
   the shared region of `min(dist to A's surface, dist to B's surface)`. A **lower bound**, reported
   as one, with the shared region's bounding-box diagonal beside it as the upper bound.

**Measured cost:** Bagger 20 pairs in 2.7 s; `as1-oc-214` 44 pairs in 10.6 s; ALICE3 58 pairs in
110 s (≈ 1.9 s/pair — the B-reps are much larger).

### 3.2 Bagger — 13 placed solids, and it is not legal

```
13 placed solids -> 78 pairs; 20 survive the AABB rejection (25.64 %)
disjoint 6 | coincident faces 11 | INTERPENETRATING 3 | contained 0 | failed 0      (2.7 s)
```

The three interpenetrating pairs are in §0. The rest of what the census learned is worth as much:

**11 coincident-face pairs (touching, zero shared volume — legal):** `BasePin`\|`Base`,
`Base`\|`Boom`, `Boom`\|`Stick`, `Stick`\|`StickCylinderInner`, `Stick`\|`BucketCylinderOuter`,
`Stick`\|`BucketLink2`, `Bucket`\|`BucketLink1`, `BoomCylinderOuter`\|`BoomCylinderInner`,
`StickCylinderInner`\|`StickCylinderOuter`, `BucketCylinderInner`\|`BucketCylinderOuter` (+1).

**6 disjoint pairs, tightest gaps:** `BucketLink1`\|`BucketCylinderInner` **0.05 cm**,
`BucketLink2`\|`BucketCylinderInner` **0.09 cm**, `Boom`\|`BoomCylinderInner` 0.75 cm,
`Boom`\|`StickCylinderInner` 0.81 cm, `Boom`\|`BoomCylinderOuter` 0.81 cm,
`Stick`\|`BucketCylinderInner` 2.04 cm.

Note the shape of it. The excavator's joints are modelled either as *exactly touching* (11 pairs) or
with a *clean sub-millimetre clearance* (0.05, 0.09 cm — round numbers, i.e. designed), and all three
defects are at joints where the modeller evidently meant one of those two and got neither.
`Stick`∩`Bucket` and `BucketLink2`∩`BucketLink1` are both bucket-linkage joints;
`Base`∩`BoomCylinderOuter` is the boom ram's mount. **This is what CAD-authored interpenetration
looks like: it is at the joints, it is tenths of a millimetre, and it was never wrong for any
purpose the model previously had.**

### 3.3 `as1-oc-214.stp` — the negative control, on a real model, and the injection that turns it

```
18 placed solids -> 153 pairs; 44 survive the AABB rejection (28.76 %)
disjoint 12 | coincident faces 32 | INTERPENETRATING 0 | contained 0 | failed 0     (10.6 s)
```

32 pairs of bolts, nuts, brackets and a plate share boundaries exactly; 12 pairs are cleanly apart
with the tightest gap 0.865 cm; nothing interpenetrates.

*A refuted experiment is not a refuted hypothesis*, so that 0 is only worth having if the instrument
would have said otherwise. `--inject NAME:DX,DY,DZ` translates one part before the census:

| injection | result |
| --- | --- |
| `plate` by (0, 0, 0.2) cm | 0 → **2 interpenetrating**: `plate`∩`l-bracket` and `plate`∩`l-bracket#1`, 9.5288 cm³ each |
| `plate` by (0, 0.2, 0) cm | 0 → **6 interpenetrating**, the plate now eating into six bolts, 0.3973 cm³ each |
| `plate` by (0.2, 0, 0) cm | 0 → **6 interpenetrating**, 0.3973 cm³ each |

and the first is checkable on paper:

| quantity | measured | expected |
| --- | ---: | --- |
| shared slab's smallest extent | **0.200000 cm** | **exactly the injected displacement** |
| sampled penetration depth | **0.083333 cm** | 0.2 × 2.5/6 = **0.083333** — the deepest point a 6-sample grid can reach in a 0.2 cm slab |
| implied contact area (V / 0.2) | 47.6438 cm² | the plate/bracket mating face |

**The instrument recovers the injected penetration exactly, in two independent quantities.**

### 3.4 ALICE3 — and a correction to how this model has been counted

**The converter's "55 solids" is a count of leaf *definitions*. The world it builds contains 206
*placed* volumes.** `Stream_J_XRay.md` §7 says "55 solids / 9266 faces" and scores "18 parts"; both
are per-definition numbers and both are correct as such. But legality is a property of *instances*:
`ST2487458_01` alone is placed **12** times, `ST1A38476_01` and the four `ST1A38494_*` six times
each, and two placements of one prototype can perfectly well interpenetrate each other. The census
is therefore over **206 instances and 21115 pairs**, not 55 and 1485.

My loader reproduces the converter's own counts exactly — **55 unique definitions / 206 instances**
for ALICE3, **13 / 13** for Bagger — which is the cross-check that it walks the XCAF assembly the
same way `O2_CADtoTGeo.py` does.

```
206 placed solids -> 21115 pairs; 1699 survive the AABB rejection (8.05 %)
```

**The full 1699-pair run did not finish inside this session** (≈ 1.9 s/pair ⇒ ≈ 54 min, and it was
still running when this was written). What *is* complete is the `ST0923290` sub-assembly — 28
instances, all 378 of its internal pairs:

```
28 placed solids -> 378 pairs; 58 survive the AABB rejection (15.34 %)
disjoint 22 | coincident faces 29 | INTERPENETRATING 7 | contained 0 | failed 0    (110.0 s)
```

| pair | shared volume | fraction of smaller | sampled depth |
| --- | ---: | ---: | ---: |
| `ST0923290_013` ∩ `ST0923290_006` | 2.9785e-01 cm³ | 5.738e-02 | 2.293e-01 cm |
| `ST0923290_013` ∩ `ST0923290_007` | 2.9785e-01 | 5.738e-02 | 2.293e-01 |
| `ST0923290_013` ∩ `ST0923290_008` | 2.9785e-01 | 5.738e-02 | 2.293e-01 |
| `ST0923290_013` ∩ `ST0923290_009` | 2.9785e-01 | 5.738e-02 | 2.293e-01 |
| `ST0923290_019` ∩ `ST0923290_022` | 3.3929e-02 | **3.177e-01** | 1.146e-01 |
| `ST0923290_019` ∩ `ST0923290_023` | 3.3929e-02 | **3.177e-01** | 1.146e-01 |
| `ST0923290_019` ∩ `ST0923290_024` | 3.3929e-02 | **3.177e-01** | 1.146e-01 |

**It is two defects, each replicated by the sub-assembly's own symmetry** — one part (`_013`)
against a four-fold family, one (`_019`) against a three-fold family. The identical volumes to five
digits are the signature of that symmetry, and they are the reason this is two fixes rather than
seven.

And `_019` is 0.02 cm (0.2 mm, a designed clearance) clear of `_027` and `_028` while sitting
0.11 cm *inside* `_022`, `_023` and `_024`. **That part is misplaced relative to its neighbours, not
mis-modelled.**

**Note which parts these are.** `ST0923290_013`, `_018` and `_019` are exactly the parts
`Stream_J_XRay.md` §7 singled out as the only ALICE3 solids whose transport was not clean (68/144/192
LOST crossings, 68 unterminated transports on `_013`). Those single-solid defects were later fixed
(`Stream_L_ALICE3Defect.md`, 18/18 clean). **The same parts also interpenetrate their neighbours,
and that has not been fixed and is a different problem.** Whether that is coincidence or a common
cause in how this sub-assembly was authored is not established here.

**22 disjoint pairs, tightest gaps:** `_019`\|`_027` and `_019`\|`_028` at **0.02 cm**,
`_{006,007,008,009}`\|`_011` at 0.0309 cm, `_018`\|`_020` at 0.05 cm, `_{006,007,008}`\|`_010` at
0.06 cm. **29 coincident (touching) pairs**, including `_018` touching `_019`, `_022`, `_023`,
`_024`, `_027` and `_028` — see §3.5 for what actually touches in them.

### 3.5 What actually touches, in the pairs the census calls `coincident`

"Coincident" as the census defines it is *distance zero and empty boolean intersection*. That is
not the same as "these two share a face", and it is worth knowing which, because **a plane
tessellates exactly and a curved surface does not** — a touching pair on a curved face is the pair
that turns into an overlap the moment anything in the pipeline chords it.

`BRepExtrema_DistShapeShape` names the sub-shapes that realise the minimum, so this is a direct
question rather than an inference. (The first attempt — testing whether a *face centroid* lies on
the contact — returned "no face centroid on the boundary" for all 11 Bagger pairs, which is itself
the finding that these are **not** whole-face contacts.)

| corpus | pairs | contact realised on |
| --- | ---: | --- |
| Bagger | 11 | **10 × plane** + edges/vertices; **1 × cylinder** — `BasePin` \| `Base`, and nothing else |
| ALICE3 `ST0923290` | 29 | **14 × plane** only; **14 × plane + B-spline**; **1 × B-spline** only |

**Bagger's touching is almost entirely plane-on-plane and therefore robust.** ALICE3's is not:
**15 of its 29 contacts involve a free-form face**, which is both the tessellation-fragile case and
exactly the population that has no exact sidecar (`NEXT.md` open item 4).

And note where Bagger's single curved contact is: `BasePin` \| `Base` — **the same pair
`TGeoManager::CheckOverlaps` reports a 0.87–0.91 cm phantom overlap for** (§3.6). That is a
suggestive coincidence and it is reported as one; the mechanism is not established, and a chord
sagitta on a pin of that size does not obviously account for 9 mm.

### 3.6 `TGeoManager::CheckOverlaps` — a third instrument, and it is broken here

The converted Bagger world (`O2_CADtoTGeo.py --exact-surfaces auto`, **13/13 leaf solids emitted as
exact `O2BVHSurfaceSolid`, 0 tessellated fallbacks**) run through ROOT's own overlap checker — the
one whose verdict a Geant4 user would actually act on, and an algorithm sharing nothing with either
instrument above.

**At ROOT's default `SetNmeshPoints(1000)`:**

```
Info in <TGeoNodeMatrix::CheckOverlaps>: Number of illegal overlaps/extrusions : 4
  [ 0] OVERLAP  Base        | BoomCylinderOuter   overlap = 0.0443384 cm
  [ 1] OVERLAP  Stick       | Bucket              overlap = 0.0145601 cm
  [ 2] OVERLAP  BucketLink2 | BucketLink1         overlap = 0.01      cm
  [ 3] OVERLAP  BasePin     | Base                overlap = 0.00855504 cm
```

**The first three are exactly the three pairs the OCCT census found**, which is genuine
corroboration from a third algorithm. But the run also prints, for every surface solid in the world:

```
Error in <o2::base::O2BVHSurfaceSolid::GetPointsOnSegments>: You should require at least 5236 points
Error in ...: at least 1802 / 1694 / 5012 / 6172 points
```

`GetPointsOnSegments` is not overridden by `O2BVHSurfaceSolid`; the base
`TGeoShape::GetPointsOnSegments` refuses when `npoints < GetNmeshVertices()` and **returns without
filling the caller's array**. The shipped display meshes have 1694–6172 vertices against ROOT's
default budget of 1000. **The default-configuration verdict is therefore not computed from points on
these shapes**, and its agreement with the census cannot be taken as evidence of anything.

Why the vertex counts are that large is visible in `CloseShape`: the display mesh is assembled by
calling `appendDisplayMesh` once **per bounded surface**, each appending its own vertices, so it is
a *soup of per-face patches with no vertices shared between surfaces*. `GetNmeshVertices` counts
the duplicates. That explains the starvation; it does not explain the phantom overlaps, and I could
not go further without touching C++.

**Feed the sampler properly and the verdict changes:**

| `SetNmeshPoints` | overlaps | entries |
| ---: | ---: | --- |
| 1000 (default) | 4 | as above, with errors on every shape |
| 8000 | **5** | `BasePin`\|`Base` **0.869474**, `BucketLink2`\|`BucketCylinderInner` **0.41**, `Stick`\|`Bucket` 0.0706829, `Base`\|`BoomCylinderOuter` 0.0443384, `BucketLink2`\|`BucketLink1` 0.01 |
| 20000 | **5** | `BasePin`\|`Base` **0.909091**, `BucketLink2`\|`BucketCylinderInner` **0.41**, `Stick`\|`Bucket` 0.0731736, `Base`\|`BoomCylinderOuter` 0.0443384, `BucketLink2`\|`BucketLink1` 0.01 |

Two of those five are refuted by exact computation:

| ROOT says | OCCT says |
| --- | --- |
| `BucketLink2` \| `BucketCylinderInner` overlap **0.41 cm** | `BRepExtrema_DistShapeShape` = **0.09 cm apart** — disjoint |
| `BasePin` \| `Base` overlap **0.87–0.91 cm** | boolean intersection has **zero volume** — touching |

A pair 0.9 mm apart cannot overlap by 4.1 mm, and a 9 mm interpenetration of a pin would have an
enormous common volume rather than none. `BRepExtrema_DistShapeShape` is the exact-distance
workhorse and the census self-test verifies it recovers a **1e-7 cm** gap exactly.

**Conclusion: `TGeoManager::CheckOverlaps` is not a trustworthy instrument on `O2BVHSurfaceSolid`,
in either configuration.** Starved by default, and wrong when fed. It is also the *only* overlap
check most users will ever run. That makes it a blocking item of its own, independent of the CAD
defects — and a cheap one, since the fix is an override of `GetPointsOnSegments` (and possibly
`GetNmeshVertices`) on the surface solid. **I did not attempt it: this session must not modify
C++.**

*(One caveat I could not eliminate: `SetNsegments` at 20/60/180 changed nothing at all, which is
expected — these shapes carry a fixed sidecar mesh — but it means I have no lever that isolates
"ROOT's sampling density" from "the shape's mesh" and the 0.87 → 0.91 cm drift with `nmesh` is
unexplained.)*

---

## 4. The fixture ladder has no census, and that is the right answer

The brief asked for the ladder fixtures. They cannot be censused, and inventing a number would be
worse than saying so. `make_boolean_fixtures.py` writes **one solid per STEP file**, ten separate
files — measured, not assumed:

| fixture | placed solids | bounding box (cm) |
| --- | ---: | --- |
| `box` | 1 | x[0, 2] y[0, 3] z[0, 4] |
| `box_minus_cyl` | 1 | x[−2, 2] y[−2, 2] z[−2, 2] |
| `box_union_box` | 1 | x[0, 4] y[0, 3] z[0, 4] |
| `cyl_cross_cyl` | 1 | x[−3, 3] y[−1, 1] z[−3, 3] |
| `cyl_inter_cyl` | 1 | x[−1, 1] y[−1, 1] z[−1, 1] |
| `cyl_plus_cone` | 1 | x[−1, 1] y[−1, 1] z[0, 5] |
| `oblique_cut_cyl` | 1 | x[−1.2, 1.2] y[−1.2, 1.2] z[0, 4.578] |
| `sphere_minus_cyl` | 1 | x[−2, 2] y[−2, 2] z[−1.908, 1.908] |
| `torus_union_cyl` | 1 | x[−3.572, 3.572] y[−3.572, 3.572] z[−2, 2] |
| `tube_window` | 1 | x[−1.502, 1.5] y[−1.5, 1.5] z[−3, 3] |

Every one is modelled about the origin. **If all ten were placed in one world at their modelled
positions, 45 of 45 pairs would have overlapping bounding boxes** — not because the ladder has a
defect but because it was never an assembly. It is ten one-solid worlds, and the correct census
result for a one-solid world is "no pairs".

---

## 5. What a reviewer should challenge

1. **The ALICE3 number is incomplete and the headline leans on 28 of 206 parts.** `ST0923290`'s 7
   interpenetrating pairs are complete and reproducible for that sub-assembly, but the full
   1699-pair run had not finished. **If the remaining 1641 pairs are clean the verdict does not
   change** (7 > 0 already), but the *scale* of the problem is unknown, and the 12-fold-placed
   `ST2487458_01` in particular has not been checked against itself. Re-run:
   `overlapCensus.py --step .../CAD_noETA.stp --out alice3.json` and quote the counts.
2. **`--pad 0.1 cm` bounds the gap table, not the overlap table.** Every overlapping pair is found
   at pad 0; a *disjoint* pair further apart than the pad never appears. So "tightest gaps" tables
   are scoped to 1 mm and nothing in them is a claim about pairs beyond that.
3. **The penetration depth is a sampled lower bound.** A 6³ grid over the shared region plus its
   vertices, capped at 120 points. For a thin slab it is short by up to one grid half-cell — proven
   in the self-test, where the true 0.25 cm reads 0.208333. The bracketing upper bound (the shared
   region's bbox diagonal) is far too loose to be useful. **A pair whose shared region is a thin
   sheet with one deep spike would be under-reported.** The shared volume and the shared bbox
   extents are the trustworthy numbers.
4. **`coincident` is defined by a volume threshold of 1e-12 cm³**, and OCCT's `BRepAlgoAPI_Common`
   of two solids meeting on a face returns an **empty shape**, not a zero-thickness sheet — so
   "coincident" is really "distance zero and empty intersection". A genuine interpenetration
   thinner than the parts' own B-rep tolerance would classify as coincident. Nothing here measures
   how close to that line the 11 Bagger and 29 ALICE3 coincident pairs sit.
5. **Coincident contacts are legal but fragile, and §3.5 only half-answers it.** The chording
   sagitta reaches **2.9e-02 cm** (`Stream_J_XRay.md` §4) — *larger* than two of Bagger's three
   real penetration depths (3.3e-03 and 2.25e-02 cm). So a pair that touches exactly in CAD can
   interpenetrate after tessellation **by more than the defects this census reports**, and the
   census on the STEP source is a lower bound on the trouble in any tessellated world. §3.5 says
   Bagger's touching is 10/11 planar (robust) and ALICE3's is 15/29 free-form (not) — but I never
   *measured* a post-tessellation census to confirm the mechanism. Converting with `--mesh` and
   re-running the census would settle it, and would cost one run.
6. **The oracle has been run only on pairs, never on a whole assembly at useful density.** The
   whole-Bagger run at 864 rays saw nothing. The occupancy-annotated crossing list is validated on
   the synthetic assembly and on named real pairs; it has never been asked to produce a full
   assembly's ground truth, and its cost at a density that would find a 0.2 cm defect across a 100 cm
   assembly is not measured.
7. **`load_assembly` uses `TopLoc_Location` and never rebuilds geometry**, so N placements of a
   prototype share one B-rep. That is what makes ALICE3 affordable at all, but it means I have
   *assumed* the XCAF locations I compose are the same ones `O2_CADtoTGeo.py` writes into `geom.C`.
   The instance/definition counts match exactly (55/206 and 13/13) and ROOT's checker reproduces
   `Base`\|`BoomCylinderOuter` at 0.0443384 cm on the converted world, which is strong evidence —
   but it is evidence, not a proof, and a per-placement matrix diff would settle it.
8. **The ROOT `CheckOverlaps` section is a finding I did not set out to make and did not fully
   diagnose.** In particular I never established *how* the default-configuration run produced three
   correct pairs from an unfilled buffer. That is worth one hour from someone who can read
   `TGeoChecker.cxx`.

---

## 6. What this invalidates elsewhere

I do not edit `NEXT.md` or `Tutorial.md`. These are the claims in them this stream moves:

* **`NEXT.md` "Where the branch stands"** — the table's "12/12 shipped, 0/0/0/0 unexplained
  disagreements" and ALICE3's "18/18 parts clean" are all **per leaf solid** and remain true. They
  should now be read next to the statement that **neither corpus composes into a legal world**.
  Every correctness number on this branch is about parts; none was ever about the assembly.
* **`NEXT.md` open item 5, "Assembly-level transport"** — the oracle half is done
  (`assemblyOracle.py`, 26 analytic checks), and the estimate in `Stream_J_XRay.md` §9 that "the
  oracle change is the small half" holds. What that item did not anticipate is that **the geometry
  it would be run on is illegal**, so the C++ half now has a prerequisite: there is no correct
  navigator behaviour to compare against inside an overlap region.
* **`Stream_J_XRay.md` §7's "55 solids"** — correct as a count of *definitions*; the world has
  **206 placed instances**. Any assembly-level statement must use the second number.
* **`Stream_J_XRay.md` §9's list of what an assembly needs** — item (iv), the **leaking** counter, is
  now expressible: it is a crossing whose `occupancy-after` is non-empty in the oracle and empty in
  the navigator, or vice versa. Item (ii), "the volume-id crossing list is the real work", is
  unchanged and is still C++.
* **`Tutorial.md` §5.4's caveat on edge identity** (topology survives conversion; face *geometry* is
  not certified) gains a sibling: **nothing in the pipeline certifies that two certified parts do
  not occupy the same space.**
* **New, and not previously recorded anywhere**: `TGeoManager::CheckOverlaps` does not work on
  `O2BVHSurfaceSolid` (§3.6). This belongs in the traps list.

---

## 7. Files

| path | what |
| --- | --- |
| `scripts/geometry/assemblyOracle.py` | the occupancy-annotated crossing-list oracle; `--self-test` = 26 analytic checks |
| `scripts/geometry/overlapCensus.py` | the pairwise legality census; `--self-test` = 14 analytic checks; `--inject` = the positive control on any model |

Neither writes into `STEP_examples/`; both take `--out` for JSON. All artefacts of this session were
written to a scratch directory and none are committed.
