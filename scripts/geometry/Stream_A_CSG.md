# Stream A — the recognition census

Date: 2026-08-01. Stream A of the parallel programme in
[`Workstreams.md`](Workstreams.md). Full brief: [`CSG_Pipeline.md`](CSG_Pipeline.md), whose §8
step 1 this document delivers and whose §7 says, in bold, *do not build the emitter before this
table exists*. Nothing was emitted; nothing outside `scripts/geometry/csg/` and this file was
touched.

Instrument: `scripts/geometry/csg/census.py` (+ `csg/occ_env.py`), pythonOCC 7.9.0 / OCCT 7.9.3.

**Read §2 first.** Four premises in `CSG_Pipeline.md` §6 and `Workstreams.md` Stream B did not
survive the measurement, and two of them change what is worth building.

---

## 0. The instrument, and why it can be believed

`census.py --self-test` builds ten solids whose every census column is known in closed form and
asserts the census against them: **43 checks, all passing.** It runs automatically before every
census unless `--no-self-test` is given. The checks that earned their place:

| check | what it caught / guards |
| --- | --- |
| per-face `BRepGProp.VolumeProperties` sums to **0**, not the solid volume | the documented trap this project already paid for once; asserted, so it can never be re-learned |
| a box has 12 **convex** edges | the *sign* of `(n1 x n2) . t`. Had it been inverted, every concave count in this document would be its complement and nothing else would have shown it |
| an L-shape has exactly 1 concave edge; a blind hole 1; a groove 2; a through hole **0** | the concavity test against four hand-computable answers, including the one that killed a premise (§2.1) |
| a NURBS-converted cylinder **is** recognised as a cylinder (gap 1.8e-15) | positive control for the Tier-0 claim |
| a genuine free-form saddle is **not** recognised as any quadric | **negative control.** Without it the "1004 B-splines are secretly quadrics" headline would be worthless |
| geometric halfspace side agrees with the `ORIENTATION` flag on every curved face | the two are computed independently; a disagreement would invalidate every exterior-halfspace count |
| `hash(TShape())` induces the same classes as `IsPartner` | the prototype key, on which every "55 solids" number depends |
| face-type histogram sums to the face count (every solid, always) | assertion inside `census_solid` |
| placement counts sum to the bodies found (every model) | assertion inside `census_model` |

Two defects in the instrument were found and fixed *by* these checks rather than by inspection:
a prototype flag that compared a prototype id against a list index (it made ALICE3 look like it
had 730 B-spline faces instead of 2377), and a `templateMatched` roll-up that counted
size-skipped solids as matches (5/55 became 10/55). Both are noted in the code.

**The strongest external validation.** `CSG_Pipeline.md` §6 measured `CAD_noETA.stp` by parsing
the STEP text directly. This census measures it through OCCT's topology, sharing no code and no
assumption with that. The two agree **exactly, on every cell**: 55 solids, 9266 faces, plane
3321, cylinder 2724, cone 409, torus 350, B-spline 2377, linear extrusion 73, revolution 12, and
15 quadric-only solids. Where this document contradicts §6 below, it is therefore not a
disagreement about what is in the file.

---

## 1. The census

Machine-readable per-model JSON is cached (`--cache`); the tables below are `--markdown` output.
Reproduce with:

```bash
cd $HOME/alisw/O2/scripts/geometry
python3 csg/census.py --self-test
python3 csg/census.py --cache /tmp/csgcache \
    --model STEP_examples/Bagger.step \
    --model STEP_examples/as1-oc-214.stp \
    --model "STEP_examples/oTOF System V3-R92cm.step" \
    --model ALICE_3_example/CAD_noETA.stp \
    --model IRIS/ST2487728_01-03032026.stp
python3 csg/census.py --report --markdown --cache /tmp/csgcache
python3 csg/census.py --ladder            # the boolean ladder fixtures, rebuilt in-process
```

Any interpreter works: the script re-execs itself under the alibuild Python 3.10.

### 1.1 Prototypes versus placements — a distinction the earlier numbers did not make

An assembly places the same body many times. oTOF is **62628 placed solids drawn from 20
distinct bodies**; ALICE3 is 206 placements of 55; IRIS 126 of 62. Every census column except
the bounding box is a property of the prototype, so the census measures each prototype once and
records its placement count. **Both roll-ups are reported and they are very different numbers** —
recognition work scales with prototypes, conversion output with placements.

### 1.2 Wall-clock

| model | size | prototypes | placements | proto faces | load | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Bagger.step` | 0.5 MB | 13 | 13 | 288 | 0.1 s | **0.3 s** |
| `as1-oc-214.stp` | 0.4 MB | 5 | 18 | 53 | 0.1 s | **0.1 s** |
| `oTOF System V3-R92cm.step` | 5.6 MB | 20 | 62628 | 1607 | 3.2 s | **5.9 s** |
| `ALICE_3_example/CAD_noETA.stp` | 32 MB | 55 | 206 | 9266 | 7.8 s | **21.5 s** |
| `IRIS/ST2487728_01-03032026.stp` | 35 MB | 62 | 126 | 9460 | 7.1 s | **18.1 s** |

Zero errors, zero skipped solids, zero non-manifold edges on all five models. The whole corpus
censuses in under a minute, so this is cheap to re-run after any converter change.

### 1.3 Summary roll-up, per prototype

| model | protos | faces | quadric-only | **after Tier 0** | 1 CSG cell | primitive template | any template | halfspaces / faces |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Bagger | 13 | 288 | **13/13** | 13/13 | 1 | 1 (`TGeoTube`) | 1 | 187 / 288 |
| as1 | 5 | 53 | **0/5** | **5/5** | 3 | 0 | 0 | — |
| oTOF | 20 | 1607 | **20/20** | 20/20 | 19 | **19 (`TGeoBBox`)** | 19 | 114 / 114 |
| ALICE3 | 55 | 9266 | **15/55** | **36/55** | 9 | 2 (`TGeoTubeSeg`) | 5 | 530 / 812 |
| IRIS | 62 | 9460 | **26/62** | **47/62** | 12 | 2 (`TGeoTubeSeg`) | 9 | 890 / 1281 |

"any template" adds `revolution/TGeoPcon-like` and `extrusion/TGeoXtru-like` (Stream D's
targets). Solids above 400 faces are not template-tested — 5 prototypes in ALICE3, 5 in IRIS, 1
in oTOF, 0 in Bagger and as1 — and are counted as neither matched nor unmatched; they are far
above any sane Tier-1/3 budget in any case.

### 1.4 Edge and concavity structure, per prototype

| model | edges | convex | concave | mixed | tangential | seam+degenerate | **concave, trusted** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Bagger | 705 | — | 139 | 0 | 130 | — | **139** |
| as1 | 126 | — | 3 | 0 | 28 | — | **3** |
| oTOF | 4558 | — | 1993 | 0 | 0 | — | **1993** |
| ALICE3 | 22835 | — | 4230 | 929 | 7645 | — | **3189** |
| IRIS | 22929 | — | 4086 | 711 | 8714 | — | **3159** |

"Trusted" excludes concave/mixed verdicts taken on edges where `|n1 x n2| < 1e-3` (0.057°), i.e.
where the sign of the dihedral test is decided by a cross product at its own noise floor.
**On ALICE3 that is 1970 of 5159 concave verdicts (38%)**, on IRIS 1638 of 4797 (34%), and on
Bagger and oTOF it is zero. This distinction matters: `#ST2487458_01`, ALICE3's largest solid,
reports 1611 concave edges of which only **32** survive the trust filter — it is a
blend-saturated body, not a body with 1611 genuine reflex edges, and a decomposition told to
split on all 1611 would produce nonsense.

### 1.5 Concave-edge distribution (the Tier-3 input), per prototype

| concave edges | Bagger | as1 | oTOF | ALICE3 | IRIS |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 (one CSG cell) | 1 | 3 | 19 | 9 | 12 |
| 1-2 | 5 | 2 | 0 | 0 | 0 |
| 3-10 | 4 | 0 | 0 | 6 | 11 |
| 11-50 | 2 | 0 | 0 | 25 | 24 |
| 51-200 | 1 | 0 | 0 | 10 | 11 |
| > 200 | 0 | 0 | 1 | 5 | 4 |

### 1.6 Tier-0 canonicalisation, per prototype

| model | non-quadric faces | canonicalise | to what | max gap |
| --- | ---: | ---: | --- | ---: |
| Bagger | 0 | 0 | — | — |
| as1 | 28 B-spline | **28 (100%)** | all cylinder | 7.3e-12 mm |
| oTOF | 0 | 0 | — | — |
| ALICE3 | 2377 B-spline + 85 swept | **1004 (41%)** | 786 cyl, 174 cone, 36 sph, 8 plane | 3.4e-8 mm |
| IRIS | 2337 B-spline + 73 swept | **1004 (43%)** | identical breakdown (shared parts) | 3.4e-8 mm |

**Of the 85 ALICE3 swept faces, zero canonicalise.** All 73 linear extrusions and all 12
revolutions are swept from a **B-spline basis curve**, not a line or a circle (measured column
`basisCurves`). See §2.2.

Full per-solid tables (one row per prototype, 20 columns) are in the cached JSON and are
regenerated by `--report --markdown`; they are not inlined here because ALICE3 and IRIS are 55
and 62 rows wide.

---

## 2. Premises that did not survive

### 2.1 "A solid with zero concave edges is convex and is a pure intersection of halfspaces" — half wrong, and the wrong half matters

Measured, in `--self-test`: a plate with a **through hole** has **zero** concave edges. At each
rim the material fills a 90° quadrant, so both rims are convex. The solid is obviously not a
convex region.

What zero concave edges does buy is the property the decomposition actually needs: the solid is
**one CSG cell** — the intersection of its own faces' oriented halfspaces, where "halfspace"
includes the *exterior* of a quadric, which is not a convex set. A blind hole, by contrast, has
exactly one concave edge (the bottom rim), and that is precisely the witness Shapiro's split
loop consumes.

Consequences, and they are not cosmetic:

- The census column is named `singleCell`, not `convex`, and this document never says "convex".
- **Tier 3's reach is much wider than a convexity framing suggests.** A bracket riddled with
  through-holes and bores is a *one-cell* solid, not a decomposition problem. Measured
  one-cell counts: oTOF 19/20, ALICE3 9/55, IRIS 12/62, as1 3/5, Bagger 1/13.
- The census therefore also reports `exteriorHalfspaces` per solid — the count of faces whose
  material lies outside its own carrier. Bagger: 43 of 288 faces. ALICE3: 1481 of 9266.

### 2.2 "ALICE3 has 85 swept faces that are almost all trivially canonical" (§4 Tier 0, §3.3) — false: none of them are

`CSG_Pipeline.md` §4 scopes Tier 0 as "convert `SURFACE_OF_{LINEAR_EXTRUSION,REVOLUTION}` of a
line/circle into the equivalent plane/cylinder/cone/torus. ALICE3: 85 faces across 8 solids, of
which 4 solids become quadric-only. Half a day."

Measured: all 73 linear-extrusion faces and all 12 revolution faces have a **B-spline basis
curve**. Extruding or revolving a B-spline yields a genuinely free-form surface. The structural
test finds nothing to canonicalise, and OCCT's own `ShapeAnalysis_CanonicalRecognition` also
declines every one of the 85. **The rescue rate on the faces Tier 0 was scoped for is 0/85.**

### 2.3 …but Tier 0 is worth far *more* than §4 claims, by a mechanism nobody wrote down

**1004 of ALICE3's 2377 B-spline *surfaces* are quadrics in disguise** — 786 cylinders, 174
cones, 36 spheres, 8 planes — recognised by OCCT at a maximum deviation of **3.4e-8 mm
(3.4e-9 cm)**, i.e. machine precision, six orders of magnitude inside any tolerance this project
uses. On `as1-oc-214.stp` it is **every one — 28 of 28 prototype faces, 70 of 70 placed, at
7.3e-12 mm**: that model has *no* analytic curved surface at all and is nonetheless entirely
quadric once decoded. This is an exporter
artefact — some CAD systems write every surface as NURBS — and it is invisible to a recogniser
that switches on the STEP entity type.

Effect on the coverage number that gates the whole programme:

| ALICE3, per prototype | §6 as written | measured |
| --- | ---: | ---: |
| quadric-only solids | 15 / 55 | 15 / 55 ✓ |
| + rescued by Tier 0 | 15 + 4 = **19** | 15 + **21** = **36** |
| solids genuinely bounded by free-form | **36** | **19** |
| faces genuinely free-form | **2377 (25.7%)** | **1373 (14.8%)** |

Tier 0 is still "half a day", still pure converter work, still no new kernel — and it takes the
analytically-representable share of ALICE3 from **15/55 to 36/55**. On IRIS, 26/62 to 47/62.

### 2.4 "36 of 55 solids and 2377 of 9266 faces are B-spline or NURBS. Neither CSG nor better trims can reach them." (`Workstreams.md`, Stream B) — the number is 42% smaller

Free-form remains ALICE3's largest single coverage item, and Stream B remains justified. But the
honest figure after Tier 0 is **19 of 55 solids and 1373 of 9266 faces**, not 36 and 2377. 17
solids move from "no representation can reach them" to "reachable by a converter change costing
half a day". Stream B's step-1 scoping table should be taken **after** Tier 0, or it will scope
`BSplineBoundedSurface` against surfaces that a canonicaliser removes for free — including the
three "hybrid band" specimens §6.3 names, whose free-form counts this census confirms exactly
(`ST2487457_01` 996 faces / 69 free-form, `ST2487455_003` 245 / 40, `ST2487736_01` 50 / 4).

### 2.5 "oTOF: 100% planar … a 1505-face planar tile, CSG is hopeless and the BVH solid is the natural answer" (§6, §9.2) — true of one body out of twenty

`oTOF System V3-R92cm.step` is 62628 placed solids from 20 prototypes, 100% planar (confirmed:
1607 planar faces, no other type). **19 of the 20 prototypes are six-plane boxes and match
`TGeoBBox` outright** — 62560 of 62628 placed bodies, 99.89%. The 1493-face support frame is the
twentieth, it has 1993 concave edges, and for *it* the §9.2 judgement is exactly right.

So oTOF is not a hard case for the constructive route with one awkward part; it is 62560
`TGeoBBox`es and one awkward part. Tier 1 alone converts it, with no boolean algebra, no trims,
no closure check and no new C++ — and it is the single largest body count in the corpus.

---

## 3. The verdict: how much of Tiers 0-3 is worth building

### Tier 0 — swept canonicalisation → **build it, but re-scope it. Highest value per hour in the stream.**

Not for the reason it was scoped (that mechanism recovers nothing, §2.2) but for the one
measured in §2.3. The implementation is *simpler* than planned: `ShapeAnalysis_CanonicalRecognition`
already does it, is already installed, and reports its own deviation, so the converter can accept
a canonicalisation on evidence (`GetGap` below a declared bound) rather than on a type switch.
Payoff: ALICE3 15→36 of 55, IRIS 26→47 of 62, as1 0→5 of 5. Keep the swept-surface handling as
well — it is cheap and other models will have line/circle sweeps — but do not budget coverage
against it.

### Tier 1 — whole-part primitive recognition → **build it. Unambiguously.**

Measured whole-part primitive matches, per prototype: oTOF **19/20**, Bagger 1/13, ALICE3 2/55,
IRIS 2/62. Per *placed body* that is **62560 + 1 + 12 + 6 = 62579 bodies** across the corpus, the
overwhelming majority in oTOF. Add Stream D's profile templates (`revolution/TGeoPcon-like`,
`extrusion/TGeoXtru-like`) and the prototype match rate rises to 5/55 on ALICE3 and 9/62 on IRIS.
No new C++, exact output, ROOT- and Geant4-native. The acceptance harness (§8 step 2) must land
first, exactly as planned — recognising a box is easy, *proving* it equals the CAD body is the
part that must not be skipped.

### Tier 2 — recognised booleans → **build it; the headline claim holds, with one correction**

**The commissioned question: do the six Bagger cylinder parts convert as booleans?**
Measured, and the answer is **yes** — but they are **unions, not differences.**

| part | faces | halfspaces | concave | axis cluster A | axis cluster B | plane dirs |
| --- | ---: | ---: | ---: | --- | --- | ---: |
| `BoomCylinderInner` | 6 | 6 | 1 | cyl r=7, r=12 | cyl r=6 | 2 |
| `StickCylinderInner` | 6 | 6 | 1 | cyl r=7, r=12 | cyl r=6 | 2 |
| `BucketCylinderInner` | 6 | 6 | 1 | cyl r=5, r=9 | cyl r=3.5 | 2 |
| `BoomCylinderOuter` | 8 | 7 | 2 | cyl r=7, r=15, r=15 | cyl r=6, r=10 | 2 |
| `StickCylinderOuter` | 8 | 7 | 2 | cyl r=7, r=15, r=15 | cyl r=6, r=10 | 2 |
| `BucketCylinderOuter` | 10 | 8 | 6 | cyl r=4, r=10, r=10 | cyl r=3.5, r=6 | 2 |

Every one is exactly **two coaxial clusters plus caps**: a barrel `TGeoTube(rmin, rmax)` and a
lug `TGeoTube` (or solid cylinder) on a second, non-parallel axis. In `TGeoCompositeShape` terms
that is `TGeoUnion(TGeoTube, TGeoTube)` with one transformation and **two leaves, not a chain** —
so §5's objection to deep composites ("a DNF with 20 cells is 20 nested boolean nodes") does not
apply here at all. The bores are `rmin` on the leaves, so no difference node is needed.

Corrections to `CSG_Pipeline.md` §4/§8 step 6, which names the pattern `tube − tube`:

1. It is `tube ∪ tube`, and the difference form is not present in these parts.
2. The union being exact is *not* established by the census. It is highly likely — the barrel
   bore is a single unsplit face in every Outer part, so the lug does not intrude into it — but
   the only thing that settles it is §8 step 2's symmetric-difference volume. Do not skip it.

**`tube_window`, and the ladder, measured** (`census.py --ladder`, which rebuilds the fixture
geometry in-process; `make_boolean_fixtures.py` is Stream E's file and was not touched):

| fixture | faces | halfspaces | concave | one cell? | Tier-1 template |
| --- | ---: | ---: | ---: | :---: | --- |
| `cyl_inter_cyl` (Steinmetz) | 6 | **2** | 0 | **Y** | none |
| `tube_window` | 4 | **4** | 0 | **Y** | none |
| `cyl_plus_cone` | 4 | 4 | 0 | Y | **`TGeoCone`** |
| `cyl_cross_cyl` (union) | 8 | 6 | 7 | . | none |

This is the sharpest result in the census. **`tube_window` needs no boolean tree at all**: it is
a *single cell* of four oriented halfspaces — `{r_z ≤ 15} ∩ {z ≥ −30} ∩ {z ≤ 30} ∩ {r_x ≥ 8}` —
because a through hole contributes no concave edge (§2.1). And `cyl_inter_cyl`, the archetype of
the quadric-quadric seam that six sessions of trim work could not represent exactly, is a single
cell of **two** halfspaces.

That is `CSG_Pipeline.md` §1's argument reduced to two numbers, and it says something about the
build order: the one remaining failing ladder fixture is reachable by the *simplest* form of the
flat cell representation, not by the boolean-tree emission Tier 2 was scoped around.

Beyond Bagger, the Tier-2 vocabulary is thin on the detector models. Of ALICE3's 15 quadric-only
prototypes, 4 are already one cell (Tier-1 territory — 3 of those 4 do match a template) and the
other 11 carry 4 to 904 concave edges over up to 167 distinct halfspaces; none is a
two-primitive boolean. IRIS is the same shape of result (10 of 26 one-cell). **Tier 2's measured
value is Bagger, `as1` (5/5 quadric after Tier 0; 3 of 5 already one cell; the nuts and bolts
carry 1-2 concave edges) and mechanical CAD generally** — which is what §6 already says, and the
census confirms rather than extends it.

### Tier 3 — general convex decomposition → **do the time-boxed spike, gate it hard, expect it to be narrow**

The evidence, per prototype, on the ≤40-face gate §6 proposes:

- **Bagger**: 7 of 13 are ≤10 faces with ≤8 halfspaces and 0-6 concave edges. Cheap. But Tier 1
  and Tier 2 already cover all 7 (1 tube + the 6 rams), so Tier 3's *marginal* Bagger yield on
  the ≤40-face gate is the four 22-44-face links and arms, at 6-18 concave edges each.
- **ALICE3 / IRIS**: of the 36 (ALICE3) and 47 (IRIS) prototypes that are quadric-only after
  Tier 0, the concave-edge distribution is dominated by the 11-200 band (35 of 55 ALICE3
  prototypes, 35 of 62 IRIS). A split per concave edge is a cell count in the tens to hundreds
  per body. That is above any budget worth defending, and §6's own conclusion ("above any sane
  Tier-3 decomposition budget") is confirmed and in fact understated.
- **The trust filter changes the picture in one specific way and only one.** On the 5 ALICE3
  prototypes with >200 concave edges the raw counts are 505, 904, 207, 461 and 1611 and the
  trusted counts are 484, 904, 207, 443 and **32** — one of the five (`ST2487458_01`) collapses
  by a factor of 50. If a Tier-3 spike is done, it must split on
  *trusted* concave edges only, and it should measure whether the tangential-band threshold
  moves the cell count, because on this corpus it moves it by 38%.
- **oTOF**: 19 of 20 already one cell (Tier 1 handles them); the 20th has 1993 concave edges.

**Recommendation.** Tier 3 as a research spike producing a *table*, exactly as §8 step 7 scopes
it — with the gate tightened from "≤40 faces" to "≤40 faces **and** ≤12 trusted concave edges",
which on this corpus is roughly 7 Bagger prototypes, 5 as1 prototypes, 19 oTOF prototypes and a
handful of ALICE3 / IRIS ones. Do **not** pre-commit to `O2CSGSolid`.

**But note the one-cell subset, which is not a spike and is not research.** A solid with zero
trusted concave edges *is* a flat cell with no decomposition step at all: no
`BRepAlgoAPI_Splitter`, none of the robustness risk §3.2 flags, no arrangement, no search. That
subset is oTOF 19/20, IRIS 12/62, ALICE3 9/55, as1 3/5, Bagger 1/13 — **and `tube_window` and
`cyl_inter_cyl`**. If any part of Tier 3 is built early, it should be this: the single-cell
emitter, gated on `concaveEdgesTrusted == 0` and accepted by symmetric-difference volume. It is
a strict subset of the general algorithm and carries none of its risk.

The measured order of value is **Tier 0 > Tier 1 > single-cell Tier 3 ≈ Tier 2 >> general
Tier 3**, and Tier 0 is the one that was mis-scoped, in both directions at once.

### One thing the census says that no tier covers

`distinctCarriers` versus faces: Bagger **187 halfspaces for 288 faces**, ALICE3 530 for 812
(on the solids small enough to test), IRIS 890 for 1281. CAD splits one carrier into two, four or
more faces routinely. Any cost model — Tier 3 cell counts, but also the BVH solid's per-solid
surface count — that is stated in faces overstates by 25-35%.

---

## 4. What this stream did *not* do

Per the brief, work stopped at the census. Not built: the acceptance harness (§8 step 2), any
tier, any emitter, any converter hook. Nothing outside `scripts/geometry/csg/` and this file was
modified; `NEXT.md`, `CodeReview_Fable*.md`, `TolerancePolicy.md`, `Workstreams.md` and
`CSG_Pipeline.md` are untouched and the corrections in §2 are recorded here for the integrating
session to fold.

The next step is unchanged from §8: **the acceptance harness, before any emitter.** Every claim
in §3 about a template "matching" is a claim about the carrier multiset, and the census does not
and cannot prove that the matched primitive equals the CAD body. `BRepAlgoAPI_Cut` both ways plus
`BRepGProp` does, and it is the same test all four tiers need.

---

## Appendix — the per-solid tables

One row per prototype; `n` is the placement count. Generated by
`csg/census.py --report --markdown`. Columns: `halfsp` = distinct carrier halfspaces (§3, last
paragraph); `concave` = raw concave + mixed dihedrals; `trusted` = the same excluding
near-tangential verdicts (§1.4); `Tier-0` = non-quadric faces that canonicalise to a quadric.
A `-` in `halfsp` means the solid is not quadric-only or is above the 400-face analysis cap.

### `Bagger.step`

Unit `mm` (x0.1 to cm); 13 prototype solids in 13 placements; load 0.1 s, census total 0.3 s. One row per prototype.

| # | n | part | faces | halfsp | plane | cyl | cone | sph | tor | free-form | swept | quadric-only | edges | concave | trusted | 1 cell? | Tier-0 | template | volume | bbox diag |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: | ---: | :---: | ---: | --- | ---: | ---: |
| 0 | 1 | `/Assembly/BasePin001` | 3 | 3 | 2 | 1 | 0 | 0 | 0 | 0 | 0 | Y | 3 | 0 | 0 | Y | 0 | TGeoTube | 3.142e+04 | 103.9 |
| 1 | 1 | `/Assembly/Base001` | 44 | 26 | 28 | 16 | 0 | 0 | 0 | 0 | 0 | Y | 122 | 16 | 16 | . | 0 | none | 2.413e+05 | 182.6 |
| 2 | 1 | `/Assembly/Boom001` | 31 | 20 | 17 | 14 | 0 | 0 | 0 | 0 | 0 | Y | 87 | 6 | 6 | . | 0 | none | 1.205e+06 | 485.4 |
| 3 | 1 | `/Assembly/Stick001` | 24 | 19 | 12 | 12 | 0 | 0 | 0 | 0 | 0 | Y | 66 | 8 | 8 | . | 0 | none | 3.339e+05 | 399.9 |
| 4 | 1 | `/Assembly/Bucket001` | 97 | 53 | 69 | 22 | 0 | 4 | 2 | 0 | 0 | Y | 247 | 68 | 68 | . | 0 | none | 5.831e+04 | 186.6 |
| 5 | 1 | `/Assembly/BucketLink003` | 23 | 14 | 12 | 11 | 0 | 0 | 0 | 0 | 0 | Y | 57 | 18 | 18 | . | 0 | none | 1.708e+04 | 99.37 |
| 6 | 1 | `/Assembly/BucketLink004` | 22 | 12 | 12 | 10 | 0 | 0 | 0 | 0 | 0 | Y | 54 | 10 | 10 | . | 0 | none | 1.23e+04 | 74.78 |
| 7 | 1 | `/Assembly/BoomCylinderOuter001` | 8 | 7 | 3 | 5 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 2 | 2 | . | 0 | none | 4.82e+04 | 200.2 |
| 8 | 1 | `/Assembly/BoomCylinderInner001` | 6 | 6 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | Y | 9 | 1 | 1 | . | 0 | none | 2.299e+04 | 195.9 |
| 9 | 1 | `/Assembly/StickCylinderInner001` | 6 | 6 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | Y | 9 | 1 | 1 | . | 0 | none | 2.299e+04 | 193.4 |
| 10 | 1 | `/Assembly/StickCylinderOuter001` | 8 | 7 | 3 | 5 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 2 | 2 | . | 0 | none | 4.82e+04 | 196.6 |
| 11 | 1 | `/Assembly/BucketCylinderInner001` | 6 | 6 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | Y | 9 | 1 | 1 | . | 0 | none | 9046 | 190.2 |
| 12 | 1 | `/Assembly/BucketCylinderOuter001` | 10 | 8 | 5 | 5 | 0 | 0 | 0 | 0 | 0 | Y | 18 | 6 | 6 | . | 0 | none | 1.601e+04 | 189 |

Prototype roll-up: {"basisCurves": {}, "byType": {"cylinder": 110, "plane": 172, "sphere": 4, "torus": 2}, "canonicalisableByType": {}, "canonicalisableTo": {}, "carriersVsFaces": [187, 288], "concaveHistogram": {"0": 1, "1-2": 5, "11-50": 2, "3-10": 4, "51-200": 1}, "concaveHistogramTrusted": {"0": 1, "1-2": 5, "11-50": 2, "3-10": 4, "51-200": 1}, "concaveNearTangential": 0, "concaveTotal": 139, "concaveTotalTrusted": 139, "edgeErrors": 0, "edgesTotal": 705, "errors": 0, "exteriorHalfspaces": 43, "faces": 288, "maxCanonicalGap": null, "mixedEdges": 0, "mixedNearTangential": 0, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 1, "quadricOnly": 13, "quadricOnlyAfterTier0": 13, "singleCell": 1, "singleCellAndQuadric": 1, "singleCellAndQuadricAfterTier0": 1, "singleCellTrusted": 1, "skipped": 0, "solids": 13, "tangentialEdges": 130, "templateMatched": 1, "templateNotAttempted": 0, "templates": {"TGeoTube": 1, "none": 12}, "tier0Rescues": 0}

Placement roll-up: {"basisCurves": {}, "byType": {"cylinder": 110, "plane": 172, "sphere": 4, "torus": 2}, "canonicalisableByType": {}, "canonicalisableTo": {}, "carriersVsFaces": [187, 288], "concaveHistogram": {"0": 1, "1-2": 5, "11-50": 2, "3-10": 4, "51-200": 1}, "concaveHistogramTrusted": {"0": 1, "1-2": 5, "11-50": 2, "3-10": 4, "51-200": 1}, "concaveNearTangential": 0, "concaveTotal": 139, "concaveTotalTrusted": 139, "edgeErrors": 0, "edgesTotal": 705, "errors": 0, "exteriorHalfspaces": 43, "faces": 288, "maxCanonicalGap": null, "mixedEdges": 0, "mixedNearTangential": 0, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 1, "quadricOnly": 13, "quadricOnlyAfterTier0": 13, "singleCell": 1, "singleCellAndQuadric": 1, "singleCellAndQuadricAfterTier0": 1, "singleCellTrusted": 1, "skipped": 0, "solids": 13, "tangentialEdges": 130, "templateMatched": 1, "templateNotAttempted": 0, "templates": {"TGeoTube": 1, "none": 12}, "tier0Rescues": 0}

### `CAD_noETA.stp`

Unit `mm` (x0.1 to cm); 55 prototype solids in 206 placements; load 7.8 s, census total 21.5 s. One row per prototype.

| # | n | part | faces | halfsp | plane | cyl | cone | sph | tor | free-form | swept | quadric-only | edges | concave | trusted | 1 cell? | Tier-0 | template | volume | bbox diag |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: | ---: | :---: | ---: | --- | ---: | ---: |
| 0 | 2 | `87728_01/ST2487730_01/ST2487736_01/SOLID` | 50 | - | 8 | 28 | 10 | 0 | 0 | 4 | 0 | . | 104 | 28 | 28 | . | 0 | none | 2.276e+06 | 7803 |
| 1 | 2 | `/ST1A38517_01/ST2487455_004/ST2487455_01` | 66 | 44 | 33 | 31 | 2 | 0 | 0 | 0 | 0 | Y | 198 | 25 | 25 | . | 0 | none | 1.203e+06 | 4872 |
| 2 | 2 | `ST1A38517_01/ST2487455_004/ST2487455_002` | 6 | 4 | 2 | 1 | 2 | 0 | 1 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | revolution/TGeoPcon-like | 548 | 138.1 |
| 3 | 2 | `ST1A38517_01/ST2487455_004/ST2487455_003` | 245 | - | 41 | 24 | 20 | 0 | 120 | 40 | 0 | . | 732 | 40 | 0 | . | 0 | none | 180.3 | 46.06 |
| 4 | 2 | `d/ST2487728_01/ST1A38517_01/ST1A38495_01` | 65 | 43 | 20 | 36 | 9 | 0 | 0 | 0 | 0 | Y | 177 | 9 | 9 | . | 0 | none | 3.316e+04 | 159.6 |
| 5 | 2 | `d/ST2487728_01/ST1A38517_01/ST1A38526_01` | 53 | 35 | 20 | 24 | 9 | 0 | 0 | 0 | 0 | Y | 141 | 9 | 9 | . | 0 | none | 2.69e+04 | 158.7 |
| 6 | 6 | `1/ST1A38517_01/ST1A38469_01/ST1A38476_01` | 148 | - | 51 | 46 | 0 | 0 | 11 | 19 | 21 | . | 381 | 85 | 57 | . | 0 | none | 6379 | 162.5 |
| 7 | 6 | `/ST1A38469_01/ST1A38494_005/ST1A38494_01` | 22 | - | 8 | 0 | 0 | 0 | 0 | 0 | 14 | . | 60 | 0 | 0 | Y | 0 | none | 2.287e+04 | 704.4 |
| 8 | 6 | `ST1A38469_01/ST1A38494_005/ST1A38494_002` | 6 | 6 | 4 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoTubeSeg | 7687 | 502.1 |
| 9 | 6 | `ST1A38469_01/ST1A38494_005/ST1A38494_003` | 6 | 6 | 4 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoTubeSeg | 3589 | 500.5 |
| 10 | 6 | `ST1A38469_01/ST1A38494_005/ST1A38494_004` | 10 | 10 | 8 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 24 | 0 | 0 | Y | 0 | extrusion/TGeoXtru-like | 1724 | 500.1 |
| 11 | 6 | `1/ST1A38517_01/ST1A38469_01/ST1A38486_01` | 148 | - | 51 | 46 | 0 | 0 | 11 | 19 | 21 | . | 381 | 85 | 57 | . | 0 | none | 6379 | 162.5 |
| 12 | 6 | `1/ST1A38517_01/ST1A38469_01/ST2487457_01` | 996 | - | 570 | 357 | 0 | 0 | 0 | 69 | 0 | . | 2544 | 70 | 0 | . | 0 | not-attempted(size) | 2509 | 700.7 |
| 13 | 2 | `/ST1A38517_01/ST2487459_005/ST2487459_01` | 202 | 167 | 23 | 91 | 6 | 0 | 82 | 0 | 0 | Y | 474 | 6 | 6 | . | 0 | none | 5024 | 159.5 |
| 14 | 2 | `ST1A38517_01/ST2487459_005/ST2487459_002` | 10 | - | 2 | 2 | 2 | 0 | 2 | 0 | 2 | . | 20 | 0 | 0 | Y | 0 | none | 275.1 | 23.75 |
| 15 | 2 | `ST1A38517_01/ST2487459_005/ST2487459_003` | 10 | - | 2 | 2 | 2 | 0 | 2 | 0 | 2 | . | 20 | 0 | 0 | Y | 0 | none | 275.1 | 23.8 |
| 16 | 2 | `ST1A38517_01/ST2487459_005/ST2487459_004` | 10 | - | 2 | 2 | 2 | 0 | 2 | 0 | 2 | . | 20 | 0 | 0 | Y | 0 | none | 275.1 | 23.74 |
| 17 | 2 | `d/ST2487728_01/ST1A38517_01/ST2487461_01` | 128 | - | 17 | 38 | 0 | 0 | 4 | 46 | 23 | . | 298 | 99 | 0 | . | 0 | none | 3963 | 159.1 |
| 18 | 2 | `d/ST2487728_01/ST1A38517_01/ST2487721_01` | 32 | 14 | 4 | 20 | 8 | 0 | 0 | 0 | 0 | Y | 76 | 4 | 4 | . | 0 | none | 1.648e+05 | 968.4 |
| 19 | 2 | `d/ST2487728_01/ST1A38517_01/ST2487462_01` | 80 | 46 | 30 | 23 | 13 | 0 | 14 | 0 | 0 | Y | 219 | 4 | 4 | . | 0 | none | 1.465e+06 | 3585 |
| 20 | 2 | `/ST1A38517_01/ST1829909_005/ST1829909_01` | 1052 | - | 383 | 511 | 142 | 0 | 16 | 0 | 0 | Y | 2748 | 505 | 484 | . | 0 | not-attempted(size) | 1.606e+06 | 285 |
| 21 | 2 | `ST1A38517_01/ST1829909_005/ST1829909_002` | 965 | - | 535 | 308 | 102 | 0 | 20 | 0 | 0 | Y | 2400 | 904 | 904 | . | 0 | not-attempted(size) | 7.997e+05 | 342.6 |
| 22 | 2 | `ST1A38517_01/ST1829909_005/ST1829909_003` | 236 | 128 | 112 | 80 | 24 | 0 | 20 | 0 | 0 | Y | 558 | 207 | 207 | . | 0 | none | 1.192e+05 | 192.5 |
| 23 | 2 | `ST1A38517_01/ST1829909_005/ST1829909_004` | 720 | - | 363 | 268 | 56 | 0 | 33 | 0 | 0 | Y | 1824 | 461 | 443 | . | 0 | not-attempted(size) | 1.953e+05 | 292 |
| 24 | 2 | `d/ST2487728_01/ST1A38517_01/ST1782525_01` | 50 | 27 | 18 | 20 | 0 | 0 | 12 | 0 | 0 | Y | 144 | 28 | 28 | . | 0 | revolution/TGeoPcon-like | 1.202e+05 | 167.4 |
| 25 | 4 | `/ST1A38517_01/ST0923290_029/ST0923290_01` | 352 | - | 78 | 0 | 0 | 0 | 0 | 274 | 0 | . | 886 | 172 | 166 | . | 246 | none | 5.174e+05 | 358.5 |
| 26 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_002` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 27 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_003` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 28 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_004` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 29 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_005` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 30 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_006` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 18 | 15 | . | 24 | none | 5191 | 123.1 |
| 31 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_007` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 20 | 15 | . | 24 | none | 5191 | 123.1 |
| 32 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_008` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 19 | 15 | . | 24 | none | 5191 | 123.1 |
| 33 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_009` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 19 | 15 | . | 24 | none | 5191 | 123.1 |
| 34 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_010` | 86 | - | 57 | 0 | 0 | 0 | 0 | 29 | 0 | . | 225 | 63 | 60 | . | 23 | none | 2.986e+04 | 130.6 |
| 35 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_011` | 193 | - | 32 | 0 | 0 | 0 | 0 | 161 | 0 | . | 541 | 77 | 72 | . | 161 | none | 2.363e+05 | 188.1 |
| 36 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_012` | 10 | - | 6 | 0 | 0 | 0 | 0 | 4 | 0 | . | 24 | 0 | 0 | Y | 4 | none | 4118 | 129.5 |
| 37 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_013` | 20 | - | 11 | 0 | 0 | 0 | 0 | 9 | 0 | . | 51 | 5 | 5 | . | 9 | none | 1.189e+05 | 190.4 |
| 38 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_014` | 138 | - | 37 | 0 | 0 | 0 | 0 | 101 | 0 | . | 351 | 68 | 68 | . | 101 | none | 7.761e+04 | 138 |
| 39 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_015` | 246 | - | 61 | 0 | 0 | 0 | 0 | 185 | 0 | . | 620 | 52 | 43 | . | 117 | none | 1.136e+05 | 121.6 |
| 40 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_016` | 59 | - | 18 | 0 | 0 | 0 | 0 | 41 | 0 | . | 150 | 20 | 20 | . | 41 | none | 2.932e+04 | 102.9 |
| 41 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_017` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 42 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_018` | 48 | - | 21 | 0 | 0 | 0 | 0 | 27 | 0 | . | 138 | 22 | 22 | . | 27 | none | 3915 | 45.26 |
| 43 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_019` | 44 | - | 14 | 0 | 0 | 0 | 0 | 30 | 0 | . | 102 | 24 | 24 | . | 30 | none | 1.71e+04 | 51.13 |
| 44 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_020` | 37 | - | 20 | 0 | 0 | 0 | 0 | 17 | 0 | . | 93 | 19 | 19 | . | 17 | none | 1.562e+04 | 51.15 |
| 45 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_021` | 33 | - | 17 | 0 | 0 | 0 | 0 | 16 | 0 | . | 79 | 17 | 17 | . | 16 | none | 7198 | 52.92 |
| 46 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_022` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 47 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_023` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 48 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_024` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 49 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_025` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 50 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_026` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 51 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_027` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 243 | 32.57 |
| 52 | 4 | `ST1A38517_01/ST0923290_029/ST0923290_028` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 243 | 32.57 |
| 53 | 4 | `d/ST2487728_01/ST1A38517_01/ST2487195_01` | 182 | - | 20 | 24 | 0 | 0 | 0 | 138 | 0 | . | 388 | 88 | 40 | . | 8 | none | 1798 | 103.7 |
| 54 | 12 | `d/ST2487728_01/ST1A38517_01/ST2487458_01` | 2034 | - | 374 | 736 | 0 | 0 | 0 | 924 | 0 | . | 4440 | 1611 | 32 | . | 0 | not-attempted(size) | 1.135e+04 | 110.8 |

Prototype roll-up: {"basisCurves": {"extrusion(bspline)": 73, "revolution(bspline)": 12}, "byType": {"bspline": 2377, "cone": 409, "cylinder": 2724, "extrusion": 73, "plane": 3321, "revolution": 12, "torus": 350}, "canonicalisableByType": {"bspline": 1004}, "canonicalisableTo": {"bspline->cone": 174, "bspline->cylinder": 786, "bspline->plane": 8, "bspline->sphere": 36}, "carriersVsFaces": [530, 812], "concaveHistogram": {"0": 9, "11-50": 25, "3-10": 6, "51-200": 10, ">200": 5}, "concaveHistogramTrusted": {"0": 12, "11-50": 27, "3-10": 6, "51-200": 6, ">200": 4}, "concaveNearTangential": 1085, "concaveTotal": 5159, "concaveTotalTrusted": 3189, "edgeErrors": 0, "edgesTotal": 22835, "errors": 0, "exteriorHalfspaces": 1481, "faces": 9266, "maxCanonicalGap": 3.389199711702714e-08, "mixedEdges": 929, "mixedNearTangential": 885, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 2, "quadricOnly": 15, "quadricOnlyAfterTier0": 36, "singleCell": 9, "singleCellAndQuadric": 4, "singleCellAndQuadricAfterTier0": 5, "singleCellTrusted": 12, "skipped": 0, "solids": 55, "tangentialEdges": 7645, "templateMatched": 5, "templateNotAttempted": 5, "templates": {"TGeoTubeSeg": 2, "extrusion/TGeoXtru-like": 1, "none": 45, "not-attempted(size)": 5, "revolution/TGeoPcon-like": 2}, "tier0Rescues": 21}

Placement roll-up: {"basisCurves": {"extrusion(bspline)": 346, "revolution(bspline)": 48}, "byType": {"bspline": 16934, "cone": 818, "cylinder": 14676, "extrusion": 346, "plane": 14438, "revolution": 48, "torus": 788}, "canonicalisableByType": {"bspline": 4016}, "canonicalisableTo": {"bspline->cone": 696, "bspline->cylinder": 3144, "bspline->plane": 32, "bspline->sphere": 144}, "carriersVsFaces": [1148, 1712], "concaveHistogram": {"0": 36, "11-50": 92, "3-10": 14, "51-200": 44, ">200": 20}, "concaveHistogramTrusted": {"0": 46, "11-50": 110, "3-10": 14, "51-200": 28, ">200": 8}, "concaveNearTangential": 11290, "concaveTotal": 29346, "concaveTotalTrusted": 8938, "edgeErrors": 0, "edgesTotal": 113358, "errors": 0, "exteriorHalfspaces": 7738, "faces": 48048, "maxCanonicalGap": 3.389199711702714e-08, "mixedEdges": 9550, "mixedNearTangential": 9118, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 12, "quadricOnly": 42, "quadricOnlyAfterTier0": 126, "singleCell": 36, "singleCellAndQuadric": 20, "singleCellAndQuadricAfterTier0": 24, "singleCellTrusted": 46, "skipped": 0, "solids": 206, "tangentialEdges": 43464, "templateMatched": 22, "templateNotAttempted": 24, "templates": {"TGeoTubeSeg": 12, "extrusion/TGeoXtru-like": 6, "none": 160, "not-attempted(size)": 24, "revolution/TGeoPcon-like": 4}, "tier0Rescues": 84}

### `ST2487728_01-03032026.stp`

Unit `mm` (x0.1 to cm); 62 prototype solids in 126 placements; load 7.1 s, census total 18.1 s. One row per prototype.

| # | n | part | faces | halfsp | plane | cyl | cone | sph | tor | free-form | swept | quadric-only | edges | concave | trusted | 1 cell? | Tier-0 | template | volume | bbox diag |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: | ---: | :---: | ---: | --- | ---: | ---: |
| 0 | 1 | `7728_01/Product1.1/Part641.1/=>[0:1:1:5]` | 50 | - | 8 | 28 | 10 | 0 | 0 | 4 | 0 | . | 104 | 28 | 28 | . | 0 | none | 2.718e+06 | 9402 |
| 1 | 1 | `2487728_01/ST1A38517_01.1/ST1A92654_01.1` | 95 | 47 | 29 | 62 | 4 | 0 | 0 | 0 | 0 | Y | 234 | 38 | 38 | . | 0 | none | 2.734e+06 | 5673 |
| 2 | 1 | `2487728_01/ST1A38517_01.1/ST1A92654_01.1` | 8 | 3 | 0 | 2 | 4 | 0 | 2 | 0 | 0 | Y | 16 | 0 | 0 | Y | 0 | revolution/TGeoPcon-like | 1096 | 174.6 |
| 3 | 1 | `2487728_01/ST1A38517_01.1/ST1709494_01.1` | 65 | 43 | 20 | 36 | 9 | 0 | 0 | 0 | 0 | Y | 177 | 9 | 9 | . | 0 | none | 3.289e+04 | 159.6 |
| 4 | 1 | `2487728_01/ST1A38517_01.1/ST1A38526_01.1` | 53 | 35 | 20 | 24 | 9 | 0 | 0 | 0 | 0 | Y | 141 | 9 | 9 | . | 0 | none | 2.663e+04 | 158.7 |
| 5 | 3 | `38517_01.1/ST1A38469_01.1/ST1703820_01.1` | 148 | - | 51 | 46 | 0 | 0 | 14 | 19 | 18 | . | 381 | 81 | 58 | . | 0 | none | 6636 | 162.5 |
| 6 | 3 | `38517_01.1/ST1A38469_01.1/ST1719125_01.1` | 22 | - | 8 | 0 | 0 | 0 | 0 | 0 | 14 | . | 60 | 0 | 0 | Y | 0 | none | 2.287e+04 | 704.4 |
| 7 | 3 | `38517_01.1/ST1A38469_01.1/ST1719140_01.1` | 148 | - | 51 | 46 | 0 | 0 | 14 | 19 | 18 | . | 381 | 81 | 58 | . | 0 | none | 6636 | 162.5 |
| 8 | 3 | `38517_01.1/ST1A38469_01.1/ST2486458_01.1` | 996 | - | 570 | 357 | 0 | 0 | 0 | 69 | 0 | . | 2544 | 70 | 0 | . | 0 | not-attempted(size) | 2509 | 700.7 |
| 9 | 3 | `38517_01.1/ST1A38469_01.1/ST1563021_01.4` | 100 | 99 | 36 | 28 | 16 | 0 | 20 | 0 | 0 | Y | 208 | 8 | 8 | . | 0 | none | 1.304e+05 | 4485 |
| 10 | 3 | `38517_01.1/ST1A38469_01.1/ST1563021_01.4` | 78 | 78 | 22 | 20 | 18 | 0 | 18 | 0 | 0 | Y | 156 | 0 | 0 | Y | 0 | none | 8.906e+04 | 5104 |
| 11 | 3 | `38517_01.1/ST1A38469_01.1/ST1563021_01.4` | 78 | 77 | 22 | 20 | 18 | 0 | 18 | 0 | 0 | Y | 156 | 0 | 0 | Y | 0 | none | 5.712e+04 | 5100 |
| 12 | 3 | `38517_01.1/ST1A38469_01.1/ST1564028_01.4` | 34 | 19 | 8 | 18 | 8 | 0 | 0 | 0 | 0 | Y | 76 | 8 | 8 | . | 0 | none | 2.964e+04 | 574.2 |
| 13 | 3 | `38517_01.1/ST1A38469_01.1/ST1564028_01.4` | 6 | 6 | 4 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | extrusion/TGeoXtru-like | 374 | 15.09 |
| 14 | 3 | `38517_01.1/ST1A38469_01.1/ST1564028_01.4` | 22 | 15 | 10 | 12 | 0 | 0 | 0 | 0 | 0 | Y | 48 | 8 | 8 | . | 0 | none | 1248 | 43.56 |
| 15 | 3 | `38517_01.1/ST1A38469_01.1/ST1564028_01.4` | 46 | 24 | 2 | 28 | 0 | 0 | 16 | 0 | 0 | Y | 92 | 6 | 6 | . | 0 | none | 7960 | 4505 |
| 16 | 3 | `38517_01.1/ST1A38469_01.1/ST1564028_01.4` | 46 | 24 | 2 | 28 | 0 | 0 | 16 | 0 | 0 | Y | 92 | 6 | 6 | . | 0 | none | 7956 | 4505 |
| 17 | 3 | `38517_01.1/ST1A38469_01.1/ST2513437_01.1` | 10 | 10 | 8 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 24 | 0 | 0 | Y | 0 | extrusion/TGeoXtru-like | 1233 | 500.1 |
| 18 | 3 | `38517_01.1/ST1A38469_01.1/ST2513437_01.1` | 6 | 6 | 4 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoTubeSeg | 5135 | 502 |
| 19 | 3 | `38517_01.1/ST1A38469_01.1/ST2513437_01.1` | 6 | 6 | 4 | 2 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoTubeSeg | 2405 | 500.5 |
| 20 | 1 | `2487728_01/ST1A38517_01.1/ST1A92680_01.1` | 202 | 167 | 23 | 91 | 6 | 0 | 82 | 0 | 0 | Y | 474 | 6 | 6 | . | 0 | none | 5024 | 159.5 |
| 21 | 1 | `2487728_01/ST1A38517_01.1/ST1A92680_01.1` | 10 | 6 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | Y | 18 | 0 | 0 | Y | 0 | revolution/TGeoPcon-like | 275.1 | 23.75 |
| 22 | 1 | `2487728_01/ST1A38517_01.1/ST1A92680_01.1` | 10 | 6 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | Y | 18 | 0 | 0 | Y | 0 | revolution/TGeoPcon-like | 275.1 | 23.8 |
| 23 | 1 | `2487728_01/ST1A38517_01.1/ST1A92680_01.1` | 10 | 6 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | Y | 18 | 0 | 0 | Y | 0 | revolution/TGeoPcon-like | 275.1 | 23.74 |
| 24 | 1 | `2487728_01/ST1A38517_01.1/ST1A92678_01.1` | 128 | - | 17 | 38 | 0 | 0 | 4 | 46 | 23 | . | 298 | 84 | 0 | . | 0 | none | 3963 | 159.1 |
| 25 | 1 | `/ST2487728_01/ST1A38517_01.1/Part27.1` | 32 | 14 | 4 | 20 | 8 | 0 | 0 | 0 | 0 | Y | 76 | 4 | 4 | . | 0 | none | 1.648e+05 | 968.4 |
| 26 | 1 | `2487728_01/ST1A38517_01.1/ST1A92729_01.1` | 80 | 46 | 30 | 23 | 13 | 0 | 14 | 0 | 0 | Y | 219 | 4 | 4 | . | 0 | none | 1.771e+06 | 4385 |
| 27 | 1 | `2487728_01/ST1A38517_01.1/ST1829909_01.1` | 1052 | - | 383 | 511 | 142 | 0 | 16 | 0 | 0 | Y | 2748 | 493 | 472 | . | 0 | not-attempted(size) | 1.606e+06 | 285 |
| 28 | 1 | `2487728_01/ST1A38517_01.1/ST1829909_01.1` | 965 | - | 535 | 308 | 102 | 0 | 20 | 0 | 0 | Y | 2400 | 888 | 888 | . | 0 | not-attempted(size) | 7.997e+05 | 342.6 |
| 29 | 1 | `2487728_01/ST1A38517_01.1/ST1829909_01.1` | 236 | 128 | 112 | 80 | 24 | 0 | 20 | 0 | 0 | Y | 558 | 191 | 191 | . | 0 | none | 1.192e+05 | 192.5 |
| 30 | 1 | `2487728_01/ST1A38517_01.1/ST1829909_01.1` | 720 | - | 363 | 268 | 56 | 0 | 33 | 0 | 0 | Y | 1824 | 454 | 436 | . | 0 | not-attempted(size) | 1.953e+05 | 292 |
| 31 | 1 | `2487728_01/ST1A38517_01.1/ST1342189_01.1` | 48 | 25 | 16 | 20 | 0 | 0 | 12 | 0 | 0 | Y | 96 | 28 | 28 | . | 0 | revolution/TGeoPcon-like | 1.923e+05 | 177.3 |
| 32 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 352 | - | 78 | 0 | 0 | 0 | 0 | 274 | 0 | . | 886 | 172 | 166 | . | 246 | none | 5.174e+05 | 358.5 |
| 33 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 34 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 35 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 36 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 818.3 | 34.41 |
| 37 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 21 | 16 | . | 24 | none | 5191 | 123.1 |
| 38 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 19 | 15 | . | 24 | none | 5191 | 123.1 |
| 39 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 21 | 15 | . | 24 | none | 5191 | 123.1 |
| 40 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 45 | - | 16 | 0 | 0 | 0 | 0 | 29 | 0 | . | 112 | 19 | 15 | . | 24 | none | 5191 | 123.1 |
| 41 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 86 | - | 57 | 0 | 0 | 0 | 0 | 29 | 0 | . | 225 | 63 | 60 | . | 23 | none | 2.986e+04 | 130.7 |
| 42 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 193 | - | 32 | 0 | 0 | 0 | 0 | 161 | 0 | . | 541 | 75 | 73 | . | 161 | none | 2.363e+05 | 188.1 |
| 43 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 10 | - | 6 | 0 | 0 | 0 | 0 | 4 | 0 | . | 24 | 0 | 0 | Y | 4 | none | 4118 | 129.5 |
| 44 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 20 | - | 11 | 0 | 0 | 0 | 0 | 9 | 0 | . | 51 | 5 | 5 | . | 9 | none | 1.189e+05 | 190.4 |
| 45 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 138 | - | 37 | 0 | 0 | 0 | 0 | 101 | 0 | . | 351 | 68 | 68 | . | 101 | none | 7.761e+04 | 138 |
| 46 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 246 | - | 61 | 0 | 0 | 0 | 0 | 185 | 0 | . | 620 | 52 | 43 | . | 117 | none | 1.136e+05 | 121.6 |
| 47 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 59 | - | 18 | 0 | 0 | 0 | 0 | 41 | 0 | . | 150 | 20 | 20 | . | 41 | none | 2.932e+04 | 102.9 |
| 48 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 49 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 48 | - | 21 | 0 | 0 | 0 | 0 | 27 | 0 | . | 138 | 22 | 22 | . | 27 | none | 3915 | 45.26 |
| 50 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 44 | - | 14 | 0 | 0 | 0 | 0 | 30 | 0 | . | 102 | 24 | 24 | . | 30 | none | 1.71e+04 | 51.13 |
| 51 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 37 | - | 20 | 0 | 0 | 0 | 0 | 17 | 0 | . | 93 | 19 | 19 | . | 17 | none | 1.562e+04 | 51.15 |
| 52 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 33 | - | 17 | 0 | 0 | 0 | 0 | 16 | 0 | . | 79 | 17 | 17 | . | 16 | none | 7198 | 52.92 |
| 53 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 54 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 55 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 56 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 57 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 106.8 | 17.24 |
| 58 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 243 | 32.57 |
| 59 | 2 | `2487728_01/ST1A38517_01.1/ST0923290_01.1` | 24 | - | 15 | 0 | 0 | 0 | 0 | 9 | 0 | . | 60 | 23 | 23 | . | 9 | none | 243 | 32.57 |
| 60 | 2 | `/ST2487728_01/ST1A38517_01.1/Part649.1` | 182 | - | 20 | 24 | 0 | 0 | 0 | 138 | 0 | . | 388 | 80 | 40 | . | 8 | none | 1798 | 103.7 |
| 61 | 6 | `2487728_01/ST1A38517_01.1/ST2487011_01.1` | 2034 | - | 374 | 736 | 0 | 0 | 0 | 924 | 0 | . | 4440 | 1320 | 0 | . | 0 | not-attempted(size) | 1.135e+04 | 110.8 |

Prototype roll-up: {"basisCurves": {"extrusion(bspline)": 73}, "byType": {"bspline": 2337, "cone": 453, "cylinder": 2888, "extrusion": 73, "plane": 3378, "torus": 331}, "canonicalisableByType": {"bspline": 1004}, "canonicalisableTo": {"bspline->cone": 174, "bspline->cylinder": 786, "bspline->plane": 8, "bspline->sphere": 36}, "carriersVsFaces": [890, 1281], "concaveHistogram": {"0": 12, "11-50": 24, "3-10": 11, "51-200": 11, ">200": 4}, "concaveHistogramTrusted": {"0": 15, "11-50": 26, "3-10": 11, "51-200": 7, ">200": 3}, "concaveNearTangential": 950, "concaveTotal": 4797, "concaveTotalTrusted": 3159, "edgeErrors": 0, "edgesTotal": 22929, "errors": 0, "exteriorHalfspaces": 1497, "faces": 9460, "maxCanonicalGap": 3.389199711702714e-08, "mixedEdges": 711, "mixedNearTangential": 688, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 2, "quadricOnly": 26, "quadricOnlyAfterTier0": 47, "singleCell": 12, "singleCellAndQuadric": 10, "singleCellAndQuadricAfterTier0": 11, "singleCellTrusted": 15, "skipped": 0, "solids": 62, "tangentialEdges": 8714, "templateMatched": 9, "templateNotAttempted": 5, "templates": {"TGeoTubeSeg": 2, "extrusion/TGeoXtru-like": 2, "none": 48, "not-attempted(size)": 5, "revolution/TGeoPcon-like": 5}, "tier0Rescues": 21}

Placement roll-up: {"basisCurves": {"extrusion(bspline)": 173}, "byType": {"bspline": 8427, "cone": 573, "cylinder": 7814, "extrusion": 173, "plane": 7488, "torus": 563}, "canonicalisableByType": {"bspline": 2008}, "canonicalisableTo": {"bspline->cone": 348, "bspline->cylinder": 1572, "bspline->plane": 16, "bspline->sphere": 72}, "carriersVsFaces": [1618, 2145], "concaveHistogram": {"0": 27, "11-50": 45, "3-10": 22, "51-200": 23, ">200": 9}, "concaveHistogramTrusted": {"0": 37, "11-50": 49, "3-10": 22, "51-200": 15, ">200": 3}, "concaveNearTangential": 5060, "concaveTotal": 12906, "concaveTotalTrusted": 4357, "edgeErrors": 0, "edgesTotal": 58453, "errors": 0, "exteriorHalfspaces": 4069, "faces": 25038, "maxCanonicalGap": 3.389199711702714e-08, "mixedEdges": 3545, "mixedNearTangential": 3489, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 6, "quadricOnly": 48, "quadricOnlyAfterTier0": 90, "singleCell": 27, "singleCellAndQuadric": 22, "singleCellAndQuadricAfterTier0": 24, "singleCellTrusted": 37, "skipped": 0, "solids": 126, "tangentialEdges": 26285, "templateMatched": 17, "templateNotAttempted": 12, "templates": {"TGeoTubeSeg": 6, "extrusion/TGeoXtru-like": 6, "none": 97, "not-attempted(size)": 12, "revolution/TGeoPcon-like": 5}, "tier0Rescues": 42}

### `as1-oc-214.stp`

Unit `mm` (x0.1 to cm); 5 prototype solids in 18 placements; load 0.1 s, census total 0.1 s. One row per prototype.

| # | n | part | faces | halfsp | plane | cyl | cone | sph | tor | free-form | swept | quadric-only | edges | concave | trusted | 1 cell? | Tier-0 | template | volume | bbox diag |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: | ---: | :---: | ---: | --- | ---: | ---: |
| 0 | 8 | `/as1/rod-assembly_1/nut_1` | 8 | - | 6 | 0 | 0 | 0 | 0 | 2 | 0 | . | 18 | 0 | 0 | Y | 2 | none | 664.4 | 28.44 |
| 1 | 1 | `/as1/rod-assembly_1/rod_1` | 4 | - | 2 | 0 | 0 | 0 | 0 | 2 | 0 | . | 6 | 0 | 0 | Y | 2 | none | 1.571e+04 | 201.2 |
| 2 | 6 | `et-assembly_1/nut-bolt-assembly_1/bolt_1` | 7 | - | 3 | 0 | 0 | 0 | 0 | 4 | 0 | . | 12 | 2 | 2 | . | 4 | none | 3201 | 49.94 |
| 3 | 2 | `/as1/l-bracket-assembly_1/l-bracket_1` | 16 | - | 8 | 0 | 0 | 0 | 0 | 8 | 0 | . | 42 | 1 | 1 | . | 8 | none | 9.686e+04 | 127.9 |
| 4 | 1 | `/as1/plate_1` | 18 | - | 6 | 0 | 0 | 0 | 0 | 12 | 0 | . | 48 | 0 | 0 | Y | 12 | none | 5.306e+05 | 235.2 |

Prototype roll-up: {"basisCurves": {}, "byType": {"bspline": 28, "plane": 25}, "canonicalisableByType": {"bspline": 28}, "canonicalisableTo": {"bspline->cylinder": 28}, "carriersVsFaces": [0, 0], "concaveHistogram": {"0": 3, "1-2": 2}, "concaveHistogramTrusted": {"0": 3, "1-2": 2}, "concaveNearTangential": 0, "concaveTotal": 3, "concaveTotalTrusted": 3, "edgeErrors": 0, "edgesTotal": 126, "errors": 0, "exteriorHalfspaces": 0, "faces": 53, "maxCanonicalGap": 7.27674612609663e-12, "mixedEdges": 0, "mixedNearTangential": 0, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 0, "quadricOnly": 0, "quadricOnlyAfterTier0": 5, "singleCell": 3, "singleCellAndQuadric": 0, "singleCellAndQuadricAfterTier0": 3, "singleCellTrusted": 3, "skipped": 0, "solids": 5, "tangentialEdges": 28, "templateMatched": 0, "templateNotAttempted": 0, "templates": {"none": 5}, "tier0Rescues": 5}

Placement roll-up: {"basisCurves": {}, "byType": {"bspline": 70, "plane": 90}, "canonicalisableByType": {"bspline": 70}, "canonicalisableTo": {"bspline->cylinder": 70}, "carriersVsFaces": [0, 0], "concaveHistogram": {"0": 10, "1-2": 8}, "concaveHistogramTrusted": {"0": 10, "1-2": 8}, "concaveNearTangential": 0, "concaveTotal": 14, "concaveTotalTrusted": 14, "edgeErrors": 0, "edgesTotal": 354, "errors": 0, "exteriorHalfspaces": 0, "faces": 160, "maxCanonicalGap": 7.27674612609663e-12, "mixedEdges": 0, "mixedNearTangential": 0, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 0, "quadricOnly": 0, "quadricOnlyAfterTier0": 18, "singleCell": 10, "singleCellAndQuadric": 0, "singleCellAndQuadricAfterTier0": 10, "singleCellTrusted": 10, "skipped": 0, "solids": 18, "tangentialEdges": 70, "templateMatched": 0, "templateNotAttempted": 0, "templates": {"none": 18}, "tier0Rescues": 18}

### `oTOF System V3-R92cm.step`

Unit `mm` (x0.1 to cm); 20 prototype solids in 62628 placements; load 3.2 s, census total 5.9 s. One row per prototype.

| # | n | part | faces | halfsp | plane | cyl | cone | sph | tor | free-form | swept | quadric-only | edges | concave | trusted | 1 cell? | Tier-0 | template | volume | bbox diag |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: | ---: | :---: | ---: | --- | ---: | ---: |
| 0 | 68 | `em V3-R92cm v2/Component1:1/oTOF v2 v1:1` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 3.695e+04 | 3478 |
| 1 | 68 | `em V3-R92cm v2/Component1:1/oTOF v2 v1:1` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 3.77e+04 | 3478 |
| 2 | 68 | `em V3-R92cm v2/Component1:1/oTOF v2 v1:1` | 1493 | - | 1493 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 4330 | 1993 | 1993 | . | 0 | not-attempted(size) | 4.556e+05 | 3481 |
| 3 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 4 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 5 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 6 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 7 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 8 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 9 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 10 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 11 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 12 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 13 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 240 | 40.64 |
| 14 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 15 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 16 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 17 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 645.6 | 138.1 |
| 18 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |
| 19 | 3672 | `em V3-R92cm v2/Component1:1/Module v1:61` | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | Y | 12 | 0 | 0 | Y | 0 | TGeoBBox | 19.2 | 32.07 |

Prototype roll-up: {"basisCurves": {}, "byType": {"plane": 1607}, "canonicalisableByType": {}, "canonicalisableTo": {}, "carriersVsFaces": [114, 114], "concaveHistogram": {"0": 19, ">200": 1}, "concaveHistogramTrusted": {"0": 19, ">200": 1}, "concaveNearTangential": 0, "concaveTotal": 1993, "concaveTotalTrusted": 1993, "edgeErrors": 0, "edgesTotal": 4558, "errors": 0, "exteriorHalfspaces": 0, "faces": 1607, "maxCanonicalGap": null, "mixedEdges": 0, "mixedNearTangential": 0, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 19, "quadricOnly": 20, "quadricOnlyAfterTier0": 20, "singleCell": 19, "singleCellAndQuadric": 19, "singleCellAndQuadricAfterTier0": 19, "singleCellTrusted": 19, "skipped": 0, "solids": 20, "tangentialEdges": 0, "templateMatched": 19, "templateNotAttempted": 1, "templates": {"TGeoBBox": 19, "not-attempted(size)": 1}, "tier0Rescues": 0}

Placement roll-up: {"basisCurves": {}, "byType": {"plane": 476884}, "canonicalisableByType": {}, "canonicalisableTo": {}, "carriersVsFaces": [375360, 375360], "concaveHistogram": {"0": 62560, ">200": 68}, "concaveHistogramTrusted": {"0": 62560, ">200": 68}, "concaveNearTangential": 0, "concaveTotal": 135524, "concaveTotalTrusted": 135524, "edgeErrors": 0, "edgesTotal": 1045160, "errors": 0, "exteriorHalfspaces": 0, "faces": 476884, "maxCanonicalGap": null, "mixedEdges": 0, "mixedNearTangential": 0, "nonManifoldEdges": 0, "orientationDisagreements": 0, "primitiveMatched": 62560, "quadricOnly": 62628, "quadricOnlyAfterTier0": 62628, "singleCell": 62560, "singleCellAndQuadric": 62560, "singleCellAndQuadricAfterTier0": 62560, "singleCellTrusted": 62560, "skipped": 0, "solids": 62628, "tangentialEdges": 0, "templateMatched": 62560, "templateNotAttempted": 68, "templates": {"TGeoBBox": 62560, "not-attempted(size)": 68}, "tier0Rescues": 0}

