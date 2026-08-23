# Stream AI — the three writer fixes

`Stream_AD` built the TGeo → STEP bench and `Stream_AF`/`AG`/`AH` ran it over six modules. Those
studies were read-only and left three defects in `O2_TGeoToCAD.py` on the record, ranked. This
stream fixes them, in that order, one commit each, and re-runs the bench to price every claim.

The three, as they were stated:

1. **the definition cache is keyed on the volume name**, and TRD reuses one name for volumes with
   different shapes — 43.5 % of its placements get another volume's solid, silently;
2. **mirror-baking through `BRepBuilderAPI_GTransform`** rewrites every analytic face as a B-spline
   and moves the volume by 0.5–1.1 % on curved carriers;
3. **a reflected placement of a volume with daughters is dropped**, with the orphaned subtree
   surviving as a free root shape at the identity.

Measured 2026-08-23 on `swenzel/bvhsurfacesolid`, aarch64, OCCT/pythonOCC 7.9.0, ROOT 6.36.10.
Every number below was run on this machine, on geometry regenerated today — the branch has since
received the RB26/2 bellow fix and the TPC prepreg fix, so the baselines are re-measured rather
than quoted from the earlier streams.

**The one-line verdict.** All three are fixed and measured: TRD's 277 018 wrongly-filled placements
go to 0, the mirrored `TPC_CDCE` goes from 0.887 % to 6.1e-15 against its analytic `Capacity()`,
every mirrored face in TPC, TRD and ABSO is analytic again (26 182 → 0 non-quadric), and the
placement-fidelity walk goes from 3/30 to 29/30 on ABSO, 20 873/37 315 to 36 583/37 315 on TPC and
499 307/636 567 to 636 566/636 567 on TRD. Every placement still missing is a volume the *forward
mapping* declines for a reason of its own — Stream AH's item 6, degenerate prism sections, and
`TGeoHalfSpace` — and nothing is missing for a placement reason any more. Self-test 71 → 101 checks,
0 failures, with at least one negative control per fix.

---

## 1. What was measured, and against what

The instruments are the four of `Stream_AD` §4.3 and `Stream_AH` §2. Two of them had to be
rewritten, because the scripts the earlier streams used lived in their session scratch and were
never committed:

- **the placement-fidelity walk (instrument D)**, which is the one that matters here. It walks the
  source `TGeoManager` accumulating a world matrix for every leaf-solid and mother-body placement,
  walks the written STEP's XCAF tree accumulating the same, and compares the two multisets. Unlike
  Stream AH's version it can score the *reflected* placements too, because after fix (iii) a
  mirrored prototype is placed with a proper location `L` and stands for `Z·V`, so its true world
  transform is `L·Z` and the walk undoes that.
- **the mirror instrument**, which reads every `__mirrored` leaf definition back out of the STEP
  and compares its volume and its analytic surface types against its un-mirrored twin in the same
  file.
- **instrument A in the world frame**, which is new. Stream AD's A scores a part in its own frame
  and cannot see where it lands; D compares world transforms and never looks at the solid. This one
  draws points in the world, asks the STEP's own shape through its XCAF location and the source
  `TGeoShape` through `TGeoManager`'s world matrix whether each point is inside, and counts the
  disagreements. A mirrored prototype is scored exactly like any other part, which is the point of
  having it.

Two lessons about the *instrument* came out of writing it, both worth keeping:

- **Matching world transforms on rounded keys does not work.** TGeo composes doubles down the tree
  and OCCT composes `gp_Trsf`s, and the two drift by ~1e-5 mm over a deep chain — 1e-8 of the ALICE
  cave. A key rounded to 1e-5 mm reported 304 TPC placements as both missing and extra. The walk
  now quantises only the rotation and matches translations with a 1e-3 mm tolerance.
- **A sorted merge is the wrong way to compare two such multisets.** Two placements a few
  nanometres apart can sort in different orders on the two sides, and the merge then calls both
  unmatched — 156 more phantom TPC mismatches. It is a spatial hash now.

Neither was a geometry finding; both would have been reported as one.

## 2. Fix (i) — the definition cache is keyed on volume identity

`O2_TGeoToCAD.py` cached definitions as `self.definitions[str(vol.GetName())]`. A `TGeoManager`
does not require volume names to be unique and ALICE's are not.

### What it cost, replayed

`collide.py` replays the exporter's depth-first build order and prices it by expanded placements:

```
volumes(identity)=17678  expanded placements=637373
names=1126  colliding=199  with DIVERGENT shapes=11
   UTCP 7, UTCH 7, USCR 6, UTG3 6, UTG4 6, UTPL 6, UTP3 3, UTP1 2, UTC3 2, UTC4 2, UTA2 2
OLD name-keyed cache: 15914 volumes get another volume's solid
   -> 277018 of 637373 expanded placements (43.5 %), 142181 of them leaves
```

which reproduces `Stream_AH` §3.2 exactly, including the per-name shape counts.

### The fix

The cache is keyed on volume identity (`ROOT.addressof(vol)`), and two volumes share one definition
only when they agree **by value**:

- a **leaf** definition is keyed on `(name, shape signature)`, where the signature is the shape's
  class and its own parameters — every mapped class has one written out by hand, and a class the
  signature does not know by value (`TGeoCompositeShape` among them) is keyed on the shape's
  address, which shares with nothing else. That is the conservative direction on purpose.
- an **assembly** definition is keyed on `(name, own shape signature, the ordered list of
  (child definition, placement) pairs)`, so two same-named volumes are one definition only if
  their own shape *and* their whole placed content agree. The child keys are interned as integers,
  or the top volume's key would be a nested copy of the whole tree.

This keeps the dedup that mattered: TRD's 17 678 volumes still cost **1 152 solids**, not 17 644.
Distinct *shape objects* in TRD number 17 644 — the shapes are not shared in the source at all —
but they take only 570 distinct values, and only 1 164 distinct (name, value) pairs.

Where one TGeo name ends up covering several definitions the emitted STEP names are disambiguated
in build order as `name`, `name#2`, `name#3`, …, and the mapping goes into the report as
`nameDisambiguation`. On TRD that is **12 names: the 11 divergent ones listed above, plus `UTC1`**,
which is two structurally distinct assemblies that happen to share both a name and their own
`TGeoBBox` — a case `collide.py` cannot see, because it compares shapes and this one differs in its
content. On PIPE it is 4 (the half-plie volumes that share two names — `NEXT.md` item 6's "24 false
positives of any name-keyed walk"); on TPC 2 (`TPC_ENDCAP`, `TPC_SECT` — `TPC_Strip`'s twins are
geometrically identical and correctly stay one definition); on ABSO 0.

### The self-check that would have caught it

Every volume that *reuses* a definition is now priced against its own `TGeoShape::Capacity()`, and
the worst deviation goes into the report as `sharedDefinitionMaxRelDev`. The defect used to be
invisible because the capacity check compared the converted solid against *its own* volume's
capacity, and that volume was the prototype.

| module | volumes | `sharedDefinitionMaxRelDev` | worst volume |
| --- | ---: | ---: | --- |
| TRD | 17 678 | **2.7e-16** | `UTCH` |
| PIPE | 202 | 4.1e-14 | `RB24B1PT` |
| TPC | 168 | 1.8e-12 | `TPC_OGRO` |
| ABSO | 33 | 0.0 | — |

`TGeoCompositeShape` is excluded from that check, because its `Capacity()` is a Monte Carlo and
disagrees with itself.

### What it changed in the output

| | PIPE before | PIPE after | TPC before | TPC after | TRD before | TRD after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| solids | 170 | **174** | 130 | 130 | 1 085 | **1 123** |
| components | 1 450 | **1 474** | 38 677 | 38 677 | 15 521 | **15 545** |
| coincident placements dropped (`--dedup-world`) | 24 | **0** | 0 | 0 | — | — |
| wall | 3.5 s | 3.7 s | 6.1 s | 6.2 s | 33.2 s | 35.3 s |

The 24 PIPE drops were the fix's first visible consequence and are worth stating plainly: they were
**not** coincidences. Four distinct half-plie `TGeoVolume`s share two names, and the coincidence key
was the name, so the dedup threw away 24 real placements of real material. The key is now the
*definition*, which is what `O2_CADtoTGeo.py` actually checks; PIPE's coincidence count is 0 and its
placement fidelity goes from 1 142/1 166 to **1 166/1 166**.

The PIPE census reconciles exactly, and it is worth writing out because three different numbers get
called "leaf placements". On the post-bellow-fix geometry:

| | |
| --- | ---: |
| leaf-solid placements in TGeo | **1 095** |
| of those, distinct by (**name**, world transform) | **1 071** — the 24 half-plie collisions |
| of those, distinct by (**definition**, world transform) | **1 095** — after fix (i), none coincide |
| mother-body placements | 71 |
| total placements the writer emits | **1 166** |
| leaf occurrences `O2_CADtoTGeo.py` then counts | **1 170** |

The last +4 is `BPS_supportBarCarbon`, a `TGeoCompositeShape` whose OCCT result is a compound of
three disjoint bodies that the reverse converter's multibody-leaf rule correctly splits (Stream AD
§4.2). Before the fix the writer dropped 24 of the 1 095 as name-keyed coincidences and the reverse
side saw 1 146.

## 3. Fix (ii) — bake a reflection with `gp_Trsf`, not `gp_GTrsf`

`Stream_AG` §3 proved the direction in nine lines. The finding stands exactly as stated: a
`gp_GTrsf` is an arbitrary affine map, so OCCT degrades every analytic carrier to a B-spline rather
than assume it survives, and on a curved carrier the approximation sits outside the true surface.

What this stream adds is that **`gp_Trsf::SetValues` accepts an improper orthogonal matrix
directly** — OCCT models it as a uniform scale of −1 — so there is no need to decompose the
reflection at all. `apply_tgeo_matrix` builds the exact `gp_Trsf` for any isometry and calls
`BRepBuilderAPI_Transform`; `gp_GTrsf` is kept only for a genuine non-uniform scale, which is what
it is for.

Measured on the self-test's hollow cylinder:

```
tube {'plane': 2, 'cylinder': 2} -> gp_Trsf mirror {'plane': 2, 'cylinder': 2}, rel 0.000e+00
the retired gp_GTrsf route: {'bspline': 4}, rel 8.850e-03
```

and on the geometry, comparing every `__mirrored` leaf against its un-mirrored twin **in the same
STEP file**:

| | TPC before | TPC after | TRD before | TRD after | ABSO before | ABSO after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| distinct mirrored definitions | 3 | 41 | 42 | 100 | 1 | 25 |
| non-quadric faces on them | 60 of 60 | **0 of 377** | 207 of 207 | **0 of 470** | 6 of 6 | **0 of 196** |
| non-quadric faces in the whole file | 60 of 1 091 | **0 of 1 408** | 207 of 7 148 | **0 of 7 250** | 6 of 235 | **0 of 229** |
| worst mirrored-vs-twin volume | 1.063e-02 | **1.134e-04** | 8.850e-03 | **3.805e-13** | 1.369e-13 | 0.0 |

The face counts are per *distinct definition name*, which is the unit this instrument works in.
Stream AH's 26 076 counts the same TRD faces once per baked **copy**, of which there were 8 505 over
those 42 names; both numbers go to zero.

The analytic oracle `Stream_AG` asked for, `TPC_CDCE`, a `TGeoPcon` whose `Capacity()` is exact:

| | volume | vs the un-mirrored twin | faces |
| --- | ---: | ---: | --- |
| before | 1085.736636 | **8.870e-03** | 10 bspline |
| after | 1076.191048 | **6.127e-15** | 6 cylinder, 4 plane |

**The one number that is not at the floor is honest and is not the mirror.** TPC's worst after the
fix is `TPC_rods__mirrored` at 1.134e-04 against its twin. Both are built from the same
`TGeoCompositeShape` and at build time they are the *same double*, 84.42487965394021 cm³, recorded
for both definitions in `TPC_report.json`. The 1.1e-04 appears only after the STEP write and
read-back, where `BRepGProp` integrates the two 23-face composites to 84.429735 and 84.420162 — it
straddles the build-time value and is the integrator, not the geometry. Both carry the identical
face table (17 cylinder, 6 plane).

This path also gained the volume invariant it never had. An isometry cannot change a volume, so
`apply_isometry` prices every bake against it and declines with the number if it fails, the same
discipline `_boolean` has carried since Stream AD. The band is 1e-6 relative.

### The orthogonality band, and what happens inside it

Requiring a placement matrix to be orthogonal at 1e-9 turned two TRD placements into baked copies
and lost the subtree under one of them. They are not scaled: `BM49/B051_1` is a `TGeoCombiTrans`
whose cosines were written out to nine digits, and its `MᵀM` departs from the identity by
**6.658e-08**:

```
+0.681268213  0  +0.732033940
 0            1   0
-0.732033894  0  +0.681268164
```

That is 0.3 µm over the whole cave (`BM49` is reached both plainly and under a reflection, so the
report lists the same placement twice). The band is therefore dimensionless and set at **1e-6**, which
separates "a rotation with hand-written constants" from a genuine non-uniform scale (O(0.1)) with
four orders of margin on each side. Three things happen around it, and none of them is silent:

- **Inside the band the matrix is snapped to the nearest rotation** before anything uses it —
  polar decomposition by Newton's iteration `R ← (R + R⁻ᵀ)/2`, which converges to the orthogonal
  factor and preserves the sign of the determinant, so a sloppy reflection stays a reflection. On
  `BM49/B051_1` that takes `|MᵀM − I|` from **6.658e-08** to 8.9e-16 and moves the rotation entries
  by **2.437e-08** — the self-test uses the nine-digit copy printed above and gets 6.721e-08 and
  2.460e-08 from it. Snapping is what keeps the rest exact: `gp_Trsf`, `apply_isometry`'s volume invariant and
  every composed world transform are then built from a matrix that really is an isometry, rather
  than one that is nearly one.
- **The snap is recorded, not absorbed.** The report carries `orthonormalisedPlacements`,
  `maxOrthogonalityDeviation`, `maxRotationCorrection`, the placement it came from, and the first
  50 individually. Corrections below 1e-12 are not recorded: composing a rotation from angles in
  double precision moves the entries by ~1e-16 and every placement in the geometry would otherwise
  be listed.
- **Outside the band the matrix is refused as a placement, loudly.** It is baked through the general
  transform, counted as `scaledPlacementsBaked`, listed in `scaledPlacements` with its deviation,
  and each one prints a `[WARN]` naming the placement, the number and — if the child has daughters —
  the fact that they cannot follow it. There is no isometry prototype for a non-uniform scale
  (`S·C·S⁻¹` is not rigid), so that limitation is real and is stated rather than hidden.

Had the band stayed at 1e-9, TRD would have lost 8 placements to a fix meant to recover 137 260.
Across PIPE, ABSO, TPC and TRD exactly **two** placements are snapped and **none** is outside the
band, so the loud path has no occurrences in these four modules and is untested against a real
scaled matrix.

## 4. Fix (iii) — reflected subtrees, as mirrored prototypes

### The choice, and why

Write `Z = diag(1, 1, −1)` and `V^ = Z·V` for a volume's *mirrored prototype*. A reflecting
placement `M` of `V` is then

```
M·V = (M·Z)·(Z·V) = (M·Z)·V^
```

and `M·Z` is proper, because `det(M·Z) = (−1)(−1) = +1`. So a reflection never needs a general
transform and never needs a solid to bake into: it becomes a *rigid placement of the child's
prototype*, which is exactly what a STEP component location can express. Applying the same identity
to `Z·M` pushes the reflection one level further down, so for a daughter with local matrix `C` under
a parent of parity `m_p`:

```
A   = Z^m_p · C
m_c = 1 if det A < 0 else 0
L   = A · Z^m_c        (proper, by construction)
```

That is the whole rule. Two consequences make it the right one of the two options on the table:

- **Every volume has at most two prototypes**, itself and `Z·itself`, shared by every reflected use
  of it. Per-leaf baking would have kept one private copy per reflecting *placement* — which is
  what the old code did, and it is why TRD's STEP carried 8 505 mirrored solids. TRD now carries
  **100 mirrored leaf definitions** and the file goes from **140.35 MB to 32.76 MB**.
- **It works for an assembly**, which per-leaf baking structurally cannot: an assembly has no solid,
  and that is the whole defect. Pushing the parity down means the reflection is only ever applied
  to a solid at the leaves, where it is an exact isometry (fix (ii)).

A non-uniform scale has no such prototype — `S·C·S⁻¹` is not rigid — so a scaling placement is
still baked, is counted separately in the report as `scaledPlacementsBaked`, and would still lose a
subtree if one ever appeared under it. In the four modules measured here, none does.

### What it recovered

Instrument D, the placement-fidelity walk, over every leaf-solid and mother-body placement:

| | free root shapes (should be 1) | | placements at TGeo's own world transform | |
| --- | ---: | ---: | ---: | ---: |
| | before | after | before | after |
| **PIPE** | 1 | 1 | 1 142 / 1 166 (97.94 %) | **1 166 / 1 166 (100 %)** |
| **ABSO** | **3** | **1** | 3 / 30 (10.00 %) | **29 / 30 (96.67 %)** |
| **TPC** | **38** | **1** | 20 873 / 37 315 (55.94 %) | **36 583 / 37 315 (98.04 %)** |
| **TRD** | **32** | **1** | 352 117 / 636 567 (55.31 %) | **636 566 / 636 567 (100.00 %)** |

and **zero placements in the STEP that are not in TGeo**, in all four — the orphaned free root
shapes that used to sit on top of each other at the identity are gone, and with them the
manufactured coincidences they caused downstream (§6).

The "extra" column is the one that says how bad it was. Before the fixes the TRD STEP carried
**147 172 placements that TGeo does not have** — 135 592 `UTCP` at the wrong transform from fix
(i)'s collision, the rest `BTSHT*`/`BTOFS*` orphans from fix (iii) — against 284 450 of TGeo's own
that it did not carry. Both are zero and one now.

Every placement still missing is a volume the forward *mapping* declines, not a placement the
writer drops:

| module | still missing | why |
| --- | ---: | --- |
| PIPE | 0 | — |
| ABSO | 1 | `AFassCentral`, a `TGeoPgon` with `rmin` stepping through 0 (Stream AH §6.1) |
| TPC | 732 | the five `TGeoHalfSpace` composites: `TPC_RR_CU` ×660, `TPC_OMH5` ×34, `TPC_IHSS` ×18, `TPC_OHS5` ×18, `TPC_MMHC` ×2 |
| TRD | 1 | `B045cut`, a `TGeoTrd1` with `dx1 = 0` (Stream AH §6.1) |

ABSO is the clearest single statement: `barrel/AFA_1` is a pure z-mirror of the whole front
absorber assembly, it used to be dropped with 27 of ABSO's 30 placements landing at the identity,
and it now lands at its true world position with the one declined `TGeoPgon` the only thing absent.

### And the solid is right, not only the transform

The world-frame instrument settles the other half of the question — a part can be at the right world
matrix and still be the wrong point set, which is exactly the risk a mirrored prototype carries.
Points drawn in the world, the STEP's shape through its XCAF location against the source
`TGeoShape` through TGeo's world matrix:

| | parts scored | points | disagreements |
| --- | ---: | ---: | ---: |
| ABSO, before | 3 | 30 000 | 0 |
| **ABSO, after** | **29** | **290 000** | **0** |

Before the fix only three of ABSO's parts could be matched to a TGeo world transform at all — the
three above `AFA` — so there was nothing to score. After it, every part of the reflected absorber
scores, and none of the 290 000 points disagrees.

## 5. The self-test — 101 checks, 0 failures

```
--- primitives vs TGeoShape::Capacity(), band 1e-9: 28 checks, 0 failures
--- composites vs an independent MC of TGeo (N=200k, 4 sigma): 7 checks, 0 failures
--- Contains agreement, TGeo vs OCCT classifier (3000 pts each): 3 checks, 0 failures
--- placement transforms: 5 checks, 0 failures
--- negative controls (a wrong parameter must be rejected): 4 checks, 0 failures
--- analytic surface types of every face: 16 checks, 0 failures
--- XCAF assembly document, written and read back: 8 checks, 0 failures
--- definition cache keyed on volume identity: 9 checks, 0 failures
--- mirror baking: exact isometry, analytic carriers: 12 checks, 0 failures
--- reflected subtrees: mirrored prototypes, placed: 9 checks, 0 failures

101 checks, 0 failures
```

Never quote the last line alone; the ten suites are 28 + 7 + 3 + 5 + 4 + 16 + 8 + 9 + 12 + 9. The
original 71 are unchanged and green.

**Suite 8, the definition cache (9 checks).** Four price the value signature itself: two equal tubes
share a signature, a wrong `rmin` does not (negative control), a wrong `Pcon` inner radius does not
(negative control), and two structurally identical `TGeoCompositeShape`s are *not* shared by value
(the conservative direction). Five price the cache on an in-memory geometry with two volumes called
`dup` whose volumes differ by 3.00 relative and two called `same` that agree: two definitions come
out, named `dup` and `dup#2`, the mapping is in the report, `same` is still emitted once and placed
twice, and every definition matches its own volume's `Capacity()` at 2.2e-16. The last check is the
control on the control — the two `dup` shapes really do differ by 3.00, so a name-keyed cache would
have been caught.

**Suite 9, the mirror bake (12 checks).** A mirrored tube keeps its volume (achieved 0.000e+00) and
its analytic faces; the retired `gp_GTrsf` route on the same tube is wrong by 8.850e-03 with all
four faces B-spline (negative control); a mirrored `Pcon` matches its analytic `Capacity()` and
keeps its faces; a real reflecting `TGeoRotation` taken from `Absorber.cxx` is exact; the mirrored
solid is not inside out; and a transform that is *not* an isometry is rejected by the volume
invariant (negative control on the invariant). Four more price the orthogonality band: TRD's
`BM49/B051_1` matrix, quoted verbatim to nine digits, is snapped from 6.721e-08 to 8.88e-16 by a
correction of 2.460e-08; an exact rotation is left alone to 1e-16; the snap keeps a reflection a reflection; and a
genuine `diag(1, 1, 2)` is still refused as a rigid placement (negative control).

**Suite 10, reflected subtrees (9 checks).** An in-memory geometry places one assembly twice, once
plainly and once through a z-mirror, and the assembly contains a tube, a rotated box and a
*third daughter that is itself placed reflected*. `TGeoManager::cd(path)` +
`GetCurrentMatrix()` is the oracle. The checks: exactly one free shape (no orphan), all 7 leaf
occurrences emitted, four named paths land where TGeo puts them (worst |Δ| 2.7e-13 mm), the
un-conjugated convention — placing the prototype at `M` instead of `M·Z` — would be wrong by
1.73 mm (negative control), a reflection under a reflection comes back as the plain volume, and
3 mirrored solid definitions carry 4 mirrored components rather than one bake per placement.

## 6. The regression gates

### 6.1 Forward, all four modules

| | PIPE (`--dedup-world`) | | ABSO | | TPC (`--dedup-world`) | | TRD | |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| | before | after | before | after | before | after | before | after |
| volumes visited, by identity | 198 | **202** | 33 | 33 | 165 | **168** | 1 126 | **17 678** |
| solids | 170 | 174 | 29 | 29 | 130 | 171 | 1 085 | 1 152 |
| mothers / pure assemblies | 78 / 28 | 78 / 28 | 11 / 3 | 11 / 3 | 69 / 30 | 87 / 41 | 467 / 40 | 492 / 40 |
| components | 1 450 | 1 474 | 38 | 39 | 38 677 | 38 714 | 15 521 | 15 607 |
| declined | 0 | 0 | 1 | 1 | 5 | 6 | 1 | 1 |
| reflecting placements | 0 | 0 | 2 | 2 | 74 | 74 | 8 531 | 8 540 |
| mirrored definitions / components | — | — | 1 / 1 | **19 / 29** | 37 / 37 | **41 / 16 540** | 8 505 / 8 505 | **62 / 16 377** |
| coincident placements dropped | 24 | **0** | — | — | 0 | 0 | — | — |
| names disambiguated | 0 | 4 | 0 | 0 | 0 | 2 | 0 | 12 |
| matrices snapped to a rotation | — | 0 | — | 0 | — | 0 | — | **2** |
| capacity check, median | 3.1e-16 | 3.3e-16 | 1.5e-16 | 1.5e-16 | 2.0e-15 | 9.1e-16 | 0.0 | 0.0 |
| STEP size | 4.36 MB | 4.41 MB | 0.61 MB | 0.68 MB | 37.06 MB | **34.21 MB** | 140.35 MB | **32.76 MB** |
| wall / peak RSS | 4.7 s | 4.7 s | 1.2 s | 1.2 s | 7.2 s / 1.16 GB | 8.6 s / 1.16 GB | 34.5 s / 2.10 GB | **24.6 s / 1.15 GB** |

The declined count rises by one on TPC because a mirrored prototype of a declining volume is its own
definition and declines too. The solid counts rise where mirrored prototypes are now first-class
definitions rather than anonymous baked copies. `TPC_Strip`'s two same-named twins are geometrically
identical and correctly stay one definition, which is why TPC disambiguates 2 names and not 3.

TRD is **faster and four times smaller** after the fixes, which is the mirror arithmetic: 8 505
private B-spline copies became 100 shared analytic ones.

### 6.2 Reverse, all four modules

`O2_CADtoTGeo.py --exact-surfaces auto --csg auto`, no `--mesh`. Both columns were run today; the
TRD "before" reproduces `Stream_AH` §3.4 exactly (9 608 definitions, 499 307 declared placements,
14 coincidences, CSG 952 / surface 8 656 / mesh 0):

| | PIPE | | ABSO | | TPC | | TRD | |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| | before | after | before | after | before | after | before | after |
| leaf solid definitions | 172 | **176** | 30 | 29 | 168 | **172** | 9 608 | **1 170** |
| leaf placements declared | 1 146 | **1 170** | 30 | 29 | 36 619 | **36 619** | 499 307 | **636 584** |
| coincidence warning | none | none | none | none | **1 535** | **none** | **14** | **none** |
| CSG / exact surface / mesh | 62 / 87 / 23 | 62 / **91** / 23 | 9 / 21 / 0 | 9 / 20 / 0 | 73 / 94 / 1 | **82** / 88 / 2 | 952 / 8 656 / 0 | **1 006** / 164 / 0 |
| declines *free-form faces* | 23 | 23 | 1 | **0** | 37 | **0** | 8 505 | **0** |
| wall | 33.3 s | 32.9 s | 2.3 s | 2.3 s | 2 min 31 s | 2 min 19 s | 3 min 10 s | **1 min 58 s** |

**PIPE** reproduces Stream AD §4 with the two corrections the fixes make to it: 176 leaf solids and
1 170 leaf placements, all at distinct world transforms, against 172 / 1 146. The tier counts are
unchanged for CSG and mesh; the four recovered half-plie definitions land in the exact-surface tier.

**TPC's 1 535 coincident placements were the exporter's own** — the orphaned `TPC_mrod` roots piling
on top of each other at the identity, exactly as Stream AG §4 diagnosed. So were **TRD's 14**, all of
them `BTOFS*`. Both are gone, and both models now come back with every declared placement at a
distinct world transform.

**The free-form decline class is gone everywhere.** On TPC it was 37 parts, all of them the mirrored
copies and only them; on TRD 8 505; on ABSO 1. Those parts now decline, or not, on their own
geometry: TPC's CSG count rises 73 → 82 and TRD's 952 → 1 006. TRD's 8 656 exact-surface parts
collapse to 164 because 8 505 of them were duplicate mirrored copies of a hundred shapes.

TPC's mesh count goes 1 → 2: `TPC_WSEG`, whose two degenerate two-edge planar wires Stream AG §4
named, now has a mirrored prototype that tessellates for the same reason.

## 7. What is still open

1. **Degenerate prism sections** (`Stream_AH` §6.1) are now the *only* reason a placement does not
   reach the STEP in PIPE, ABSO and TRD: `AFassCentral` and `B045cut`, one each.
2. **`TGeoHalfSpace`** is the only such reason in TPC, and it is 732 placements over five
   composites. Still the cheapest coverage win on the board.
3. **A non-uniform scale on a volume with daughters** would still lose its subtree. There is no
   isometry prototype for it and the report counts it, but nothing in PIPE, ABSO, TPC or TRD
   exercises it. Whether the full Run 3 geometry does is unmeasured.
4. **`TGeoPara`**, the last unmapped Run 3 shape class, is untouched.
5. **The full-geometry STEP write.** The bracket is unchanged in kind but the numbers moved: TPC
   writes 38 714 components and TRD 15 607, and the full model's 74 601 still segfaults. The fixes
   make the full model *smaller* (the mirror bake was 9 211 private copies), so it is worth
   re-running, and was not.
6. **The reverse recogniser's two missing detectors** — the sloped prism and the revolved profile —
   are where the decline histograms still point, unchanged by anything here.

### What was not run, and should be said

- Stream AD's **instruments B and C in full** — the recognised CSG scored against the source shape,
  and the whole capacity table per shape class — were not re-run. The forward converter's own
  per-definition capacity check is in every report and its medians are in §6.1; nothing in these
  three fixes touches the shape *mapping*, which is what B and C measure.
- **ITS and MAG**, the two remaining modules of the six-module corpus, were not re-run. MAG's five
  reflections are all planar and all of leaf volumes, so fix (ii) applies to it and fix (iii) does
  not; ITS was not exercised at all.
- **The full Run 3 geometry** was not run, forward or backward.
- The **loud path for a genuinely scaled matrix** has no occurrences in the four modules measured
  here, so it is exercised only by the self-test's `diag(1, 1, 2)` control.

## 8. Dead premises, recorded

- **"Keying on identity will explode the definition count."** It looked certain: TRD's 17 678
  volumes hold 17 644 *distinct shape objects*, so identity keying naively means 16× the solids.
  They take 570 distinct **values** and 1 164 distinct (name, value) pairs, so keying on value with
  identity as the fallback costs 1 152 solids against the old 1 085 — 6 %, not 1 600 %.
- **"A reflection has to be baked, because a STEP location must be proper."** The premise is true
  and the conclusion does not follow: `M = (M·Z)·Z` moves the reflection out of the location and
  into the definition. The old code accepted the premise and drew the conclusion, and that is the
  whole of defect (iii).
- **"`gp_Trsf::SetValues` will refuse a determinant of −1."** It accepts it, reports
  `IsNegative()`, and `BRepBuilderAPI_Transform` then produces a bit-exact mirrored solid with its
  carriers intact and its orientation right way out. The decomposition into a canonical mirror
  composed with a proper rotation was written and then not needed for the bake — only for the
  *placement*, where it is the point.
- **"Tightening the orthogonality test is free."** It cost TRD 8 placements, because a hand-written
  ALICE rotation matrix is only orthogonal to 6.7e-08. A band has to be priced against the input it
  will actually see, not against the defect it is aimed at.
- **"304 TPC placements at the wrong world transform."** They were at the right one; the instrument
  rounded world transforms to 1e-5 mm and TGeo and OCCT compose doubles differently. Then 156 more,
  from comparing two multisets by sorted merge. Both were instrument defects reported as geometry
  defects for as long as it took to look.
- **"The mirrored composite that is still 1e-4 off is the mirror."** Both definitions hold the same
  double at build time. The 1.1e-04 is `BRepGProp` integrating two 23-face composites after a STEP
  round trip, and it straddles the build-time value.

## 9. Reproducing it

```bash
# 1. the geometry (sim env; NEVER with the pythonOCC prepends; O2_ROOT for the stage libs)
export O2_ROOT=$B/stage VMCWORKDIR=$B/stage/share
$B/stage/bin/o2-sim-serial -n 0 -g boxgen -m PIPE --configKeyValues "align-geom.mDetectors=none"

# 2. the converter env is the O2 env PLUS the OCC prepends, in one process
python3 scripts/geometry/O2_TGeoToCAD.py --self-test                     # 101 checks
python3 scripts/geometry/O2_TGeoToCAD.py o2sim_geometry.root PIPE.step \
        --dedup-world --report PIPE_report.json

# 3. back again
python3 scripts/geometry/O2_CADtoTGeo.py PIPE.step --output-folder rt -o geom.C \
        --exact-surfaces auto --csg auto \
        --surface-report rt/surface_report.json --csg-report rt/csg_report.json

# 4. the instruments
python3 placecmp.py    o2sim_geometry.root PIPE.step    # D, placement fidelity
python3 worldscore.py  o2sim_geometry.root PIPE.step --points 10000   # A, world frame
python3 mirrorcheck.py PIPE.step                        # mirrored vs un-mirrored twins
python3 collide.py     o2sim_geometry.root              # the name-collision damage
```

`placecmp.py`, `worldscore.py`, `mirrorcheck.py` and `collide.py` live in this session's scratch and
are not committed; every number they produced is quoted above with the sample size it came from.
