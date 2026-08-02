# Stream R — trimming a face by its neighbours' half-spaces: design, and the measurement that decides it

Date: 2026-08-02. Branch `swenzel/bvhsurfacesolid`. **Design and measurement only — no C++ was
changed, nothing was built, no sidecar was written.** Companion and prerequisite:
[`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md), whose census this acts on;
[`Stream_Q_EllipseTrim.md`](Stream_Q_EllipseTrim.md), the other outcome of that census;
[`Stream_F_EdgeIdentity.md`](Stream_F_EdgeIdentity.md) §8 for what closure does and does not certify.

---

## 0. The answer first

**The naive conjunction is not sufficient. On the 15 ALICE3 solids `Stream_O` says implicit trims
would rescue, it carries the faces that block emission on exactly one of them.** Restricted to the
101 faces that fail today — the ones carrying `Stream_O`'s 691 bucket-B edges — the rule
"inside iff on the correct side of every neighbouring surface" is right on 50 of 101 faces, and
**1 of 15 solids has all of its blocking faces right**. The disagreements are not tolerance
artefacts: 89 of 89 sampled false positives were independently confirmed off the face by
`BRepExtrema_DistShapeShape`, and the worst under-acceptance is **3.09 cm** deep on
`ST0923290_011` f#98 and **0.66 cm** on the spherical faces of `ST0923290_002/003/004/005`.

**A specific, bounded extension does carry it: 13 of 15.** Replace the single sense per surface by a
small set of admissible *sign vectors* — a disjunction of conjunctions over the same half-spaces,
which is exactly "the trimmed region is a union of cells of the arrangement". Median cell count is
**1** (where it degenerates to the conjunction), p90 **3**, max **19**. It fixes 44 of the 51 faces
the conjunction gets wrong and costs one `uint32` mask lookup per query.

**The 7 faces it still gets wrong all fail for one identified reason, and it is not the co-surface
idea.** Every one of them has a boundary edge shared with another face lying on the *same* analytic
surface — a NURBS patch seam, which the exporter creates when it writes one cylinder as two or four
patches. Such an edge has **no** co-surface representation at all: the neighbour's implicit function
is identically zero on the whole face, so its sign carries no information. Where that seam is iso it
is carried by the quadric record's own (φ, h) box; where it is not, nothing carries it, and the face
over-accepts by up to **0.87 cm**.

**The control result is the one that should stop a wholesale replacement.** On geometry that is
known-good and ships today, the conjunction is *worse* than what exists:

| corpus | ships today | solids whose every face the conjunction gets right | … the cell set gets right |
| --- | --- | ---: | ---: |
| ladder fixtures | 10/10 | **10 / 10** | 10 / 10 |
| `as1-oc-214.stp` | 5/5 | 4 / 5 (`l-bracket`) | **5 / 5** |
| `Bagger.step` | 13/13 | **1 / 13** | 2 / 13 |
| ALICE3, the 6 solids that emit | 6 eligible | 2 / 6 | 3 / 6 |
| ALICE3, the 15 `Stream_O` solids | 0/15 | 0 / 15 (all faces) · **1 / 15** (blocking faces) | 12 / 15 (all faces) · **13 / 15** (blocking faces) |

Bagger's worst face is `Boom` f#0, a planar face with 15 trimming surfaces whose region spans 8
cells: the conjunction is wrong by **15.1 cm**. `as1-oc-214`'s single failure is `l-bracket` — the
literal L-shaped planar face, which is the textbook case and is present in the corpus.

**So an implicit trim must be a per-*face* option that the converter verifies for that face, never a
replacement for the parametric trim.** Written that way it is a real gain: it is the only exact
representation available for the 691 edges, it makes shared-edge geometry identical by construction
rather than by measurement, and — contrary to the concern that motivated this brief — it is
**cheaper** per ray than what ships, not 10× slower (§7).

---

## 1. The question, stated as a set equality

`Stream_O` measured that 691 of ALICE3's 763 rejected boundary edges lie, to within the source
B-rep's own declared tolerance, on **both** their own face's analytic surface and a neighbouring
face's analytic surface. That licenses describing the boundary without any parametric curve.

It does **not** license the containment rule that suggests itself. Write `f_k` for the signed
implicit function of the k-th distinct neighbouring surface. The proposed rule is

> **R1** — a point `p` on the face's surface is inside the patch iff `sense_k · f_k(p) ≥ 0` for
> every `k`.

That is a claim of **set equality** between the trimmed face and an intersection of half-spaces, and
it is right exactly when the trimmed region is a single cell of the arrangement of those surfaces
*and* that cell is the face. Both halves can fail, independently, and this document measures both.

Three rules are scored throughout:

| | rule | what it can express |
| --- | --- | --- |
| **R0** | the conjunction alone | an intersection of half-spaces, unbounded where no neighbour bounds it |
| **R1** | R0 **∧ the record's own (u, v) box** | R0 confined to `[φStart, φStart+φSweep] × [hMin, hMax]`, which the sidecar already stores per quadric — and which is *itself* a conjunction of at most four more half-spaces (§6.3) |
| **R2** | a set of admissible sign vectors | any union of arrangement cells — a DNF over the same half-spaces |

R1 is the honest baseline, because nothing has to be added to the format to get the box: it is
already there. R0 is measured separately (§5.4) only to show that the box is load-bearing.

---

## 2. The instrument, and how to re-run every number here

Two probes in [`probes/`](probes/), standalone, no build products involved.

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2 ; SW=$HOME/alisw/sw/ubuntu2404_aarch64
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
cd $HOME/alisw/O2/scripts/geometry

# ALICE3, every leaf solid: 7773 faces, 8.0M surface samples, ~5 min, ~1.2 GB
python3 probes/implicitTrimValidate.py --model ALICE3:ALICE_3_example/CAD_noETA.stp \
    --grid 32 --far-grid 16 --far-check 300 --verbose --json /tmp/r/alice3.json

# the controls: Bagger, as1-oc-214, the ten ladder fixtures
python3 probes/implicitTrimValidate.py --model Bagger:STEP_examples/Bagger.step \
    --model as1:STEP_examples/as1-oc-214.stp --fixtures \
    --grid 32 --far-grid 16 --far-check 300 --json /tmp/r/controls.json

# every table below
python3 probes/implicitTrimValidateReport.py /tmp/r/*.json --per-solid --split-population

# the two negative controls (§3.2) -- both must move the counts
python3 probes/implicitTrimValidate.py --fixtures --grid 24 --far-grid 0 --flip-sense 0
python3 probes/implicitTrimValidate.py --fixtures --grid 24 --far-grid 0 --perturb-radius 0.001

# the independent ground-truth check (§3.3) and the density sweep (§5.5)
python3 probes/implicitTrimValidate.py --model ALICE3:ALICE_3_example/CAD_noETA.stp \
    --solids ST0923290_011,ST0923290_016 --grid 64 --far-grid 0 --verify 40 \
    --solid-crosscheck 200 --json /tmp/r/verify.json
```

### 2.1 The method, in one paragraph

For every face carrying an analytic surface — **stored** in the STEP file, or **recognised** by the
shipping `_recognize_analytic_surface`, reached through `probes/trimEdgeCensus.py`'s own helpers so
the recogniser is byte-for-byte the shipped one — walk the wires, take each boundary edge's adjacent
face through `TopExp::MapShapesAndAncestors(EDGE, FACE)` (identity, never proximity), take that
face's analytic surface, and **deduplicate geometrically**: two neighbours that lie on the same
plane contribute one half-space. Then sample the face's *surface* on an `N × N` grid over its own
`(u, v)` rectangle, evaluated on the real `Geom_Surface`, and classify every sample two ways.

### 2.2 Why `BRepTopAdaptor_FClass2d` is the honest ground truth

`FClass2d` answers exactly the question asked — *is this `(u, v)` inside this face's trim wires* —
with no projection, no proximity band, and no dependence on the solid being closed or on any other
face. `BRepClass3d_SolidClassifier` answers a **different** question (is this point inside the
solid); it is undefined ON the boundary, which is where every one of these samples sits by
construction; and it would conflate a mis-trimmed face with a neighbouring face happening to cover
the same region. It is therefore used only as a cross-check (§3.3), not as the criterion.

Sampling over the face's **own** `(u, v)` rectangle is also deliberate: the rectangle contains the
trimmed region *and* the parts of the same surface just outside it, which is precisely where a
conjunction over-accepts, and it needs no extrapolation of a NURBS surface beyond its knot range.
It is also, exactly, the region the kernel's `[φStart, φSweep] × [hMin, hMax]` box admits — so the
near-field measurement *is* the score of R1, with the box included and nothing assumed.

### 2.3 The measurement that is independent of any sense convention

For each sample compute the **sign vector** `σ(p) = (sign f_1, …, sign f_K)` — the cell of the
arrangement the point lies in. Then

- **`cellsIn`** = the number of distinct sign vectors over the truth-IN samples. `cellsIn == 1` is
  *equivalent* to "some conjunction of one sense per surface can express this region". `cellsIn == n`
  means the smallest exact implicit description over these surfaces is a DNF with `n` terms.
- **`leak`** = truth-OUT samples whose sign vector is one of the IN cells. These are points **no**
  DNF over the same surfaces can reject: the arrangement cell genuinely straddles the boundary.
  This is the number that decides whether the idea works at all.

Neither quantity depends on how a sense is chosen, so neither can be tuned. R1's own senses are read
off the truth-IN sample furthest from the boundary — one bit per surface, exactly what the converter
would store — and R1 is then scored as false positives (rule IN, truth OUT) and false negatives
(rule OUT, truth IN), each with the 3D distance from the sample to the face's boundary wire and with
the depth on the wrong side of the deciding surface, in cm.

---

## 3. Three self-checks, run before anything below was believed

### 3.1 The recogniser and the neighbour walk are the shipped ones

The probe imports `O2_CADtoTGeo`'s `extract_graph` (so the leaf solids are the converter's own
`def_shapes`) and `probes/trimEdgeCensus.py`'s `stored_model` / `recognize` / `recognized_model`,
whose face verdicts `Stream_O` §3 checked against the shipping `recognize_and_extract_face` over 998
ALICE3 faces with **0 mismatches**. The census was re-run here and reproduces `Stream_O` cell for
cell: 4533 edges, 763 non-iso, **B 691 / C 56 / A 16**, rollup 15 of 16 solids covered.

### 3.2 The instrument can fail — two negative controls, both move the count

A validation that cannot fail has not passed. On the ten ladder fixtures, where the unperturbed
answer is **0 false positives, 0 false negatives, 0 leak on all 56 faces**:

| control | faces disagreeing | FP | FN | leak |
| --- | ---: | ---: | ---: | ---: |
| none | **0 / 56** | 0 | 0 | 0 |
| `--flip-sense 0` (invert one stored sense per face) | **56 / 56** | 3342 | 28482 | 0 |
| `--perturb-radius 0.05` (move every trimming surface 500 µm) | 39 / 56 | 720 | 440 | 1520 |
| `--perturb-radius 0.001` (**10 µm**) | 4 / 56 | 16 | 0 | 16 |

A 10 µm displacement of a trimming surface moves the count off zero. The instrument resolves far
below anything this document argues about.

### 3.3 The ground truth itself, checked against a second OCCT algorithm

Every false positive is a claim that a point is on the face's surface and *not* on the face.
`BRepExtrema_DistShapeShape` decides that in 3D against the trimmed face as a shape and shares no
code with `FClass2d`. Over `ST0923290_011/014/016/018` at grid 32 and over the 15 target solids at
grid 64:

```
   grid 32:  179 of 179 sampled false positives confirmed off the face
   grid 64:   89 of  89 sampled false positives confirmed off the face
   minimum confirmed distance   1.2e-05 cm      median per face   5.2e-03 … 8.0e-01 cm
```

Not one false positive is an `FClass2d` artefact. The one face where the confirmed distance is
1.2e-05 cm (`ST0923290_011` f#112, a single sample) is exactly the tolerance artefact the brief asks
to be separated from a broken design, and the instrument separates it: it sits 1.3e-05 cm from the
boundary, while `ST0923290_011` f#98's 284 leaking samples sit up to **0.46 cm** from it.

`BRepClass3d_SolidClassifier`, run on 200 samples per face over the same four solids, agrees with the
`FClass2d` verdict on **88731** of 89790 samples and disagrees on 1059 (1.2%). That residual is
expected and is why it is not the criterion: it answers a different question, it is undefined on the
boundary, and a point outside face A but on face B reads `ON` for the solid while `OUT` for A.

---

## 4. What the design population actually is

`Stream_O` classified only the *rejected* edges. A face is carried by an implicit trim only if
**every** one of its boundary edges is, so the populations here are wider and must be named:

| population | definition | faces (all corpora) |
| --- | --- | ---: |
| **P0** | every face with an analytic surface | 8170 |
| **P1** | P0, and no boundary edge has a free-form neighbour | 6475 |
| **P2** | P1, and no seam edge and no neighbour lying on the face's own surface | 4066 |

**All 678 faces of the 15 `Stream_O` solids are in P1** — no face of theirs has a free-form
neighbour, which independently confirms `Stream_O` §7's population. Only 298 of them are in P2.

The gap between P1 and P2 is the single most consequential structural fact in this document:

> **380 of the 678 faces (56%) have a neighbouring face lying on their own analytic surface.**

That happens because the ALICE3 exporter writes one cylinder, cone or sphere as two or four NURBS
patches. Adjacent patches are recognised as the *same* analytic surface, so the shared edge's
"co-surface trim" is a function that is identically zero over the entire face — measured max |f| of
**0.0 cm** — and its sign is pure floating-point noise. Such an edge is not trimmed by a
neighbour at all. Where it is iso (a constant-φ generator, a constant-h rim) the record's own
(φ, h) box carries it and R1 is exact; where it is not, **nothing** carries it, and that is the sole
mechanism behind every one of the 7 remaining R2 failures (§5.3).

Two smaller mechanisms, both counted and both harmless here: 77 faces across the corpora have a
periodic **seam** edge (a face that is its own neighbour), and 153 have a **degenerate** edge (a cone
apex, a sphere pole) which is a point and needs no trim. Zero faces in the whole corpus have an edge
shared by more than two faces.

---

## 5. The measurement

### 5.1 The number that decides — the faces that block emission today

A solid emits only if every face extracts. The faces that block the 15 `Stream_O` solids are the
101 carrying at least one of the 691 bucket-B edges; the rest of each solid already extracts through
the parametric path. Cross-referencing `probes/trimEdgeCensus.py`'s per-face `nRejected` with this
probe, face by face:

| | grid 32 | grid 64 |
| --- | ---: | ---: |
| blocking faces R1 gets exactly right | 53 / 101 | **50 / 101** |
| blocking faces R2 gets exactly right | 97 / 101 | **94 / 101** |
| **solids whose blocking faces are ALL right, R1** | **1 / 15** | **1 / 15** |
| **solids whose blocking faces are ALL right, R2** | **13 / 15** | **13 / 15** |

The solid-level verdict is identical at both densities. The one solid R1 carries is
**`ST0923290_014`** (6 blocking faces, all single-cell). The two R2 misses are **`ST0923290_011`**
and **`ST0923290_016`**.

### 5.2 Every blocking face R1 gets wrong, named (grid 64)

51 faces, and they fall into exactly three families.

| family | faces | K | cells | worst depth | leak |
| --- | ---: | ---: | ---: | ---: | ---: |
| **spherical patches** — `_002/003/004/005/017/022…026/027/028` f#8, f#9, f#10 | 36 | 5–6 | 5–8 | **0.663 cm** | **0** |
| **`ST0923290_011`** f#72/73/98/106/107/109/111/112/115/134/141/173 | 12 | 1–6 | 1–3 | **3.094 cm** | 5 faces |
| **`ST0923290_016`** f#16/17/18 | 3 | 4–6 | 3–6 | **0.288 cm** | 2 faces |

Over all **678** faces of the 15 solids (not only the blocking ones) R1 is wrong on **89**, and every
one is accounted for: **85** fail through non-convexity — the region spans more than one arrangement
cell — and **4** through a coincident neighbour. There is no third mechanism and no residue.

The spherical family is the cleanest statement in the whole measurement: **36 faces, 0 false
positives, up to 2821 false negatives each, and leak exactly 0.** The conjunction *under*-accepts —
the spherical patch is a region of the sphere bounded by 5 or 6 surfaces whose union of cells is 5
to 8, not 1. Nothing about the surfaces is wrong; the region is simply not convex in the half-space
sense. R2 carries all 36 exactly. These are `Stream_O` §8.1's faces — the same
`ST0923290_002` f#8 whose 8 boundary edges no single polar axis can make iso.

### 5.3 Every face R2 also gets wrong, named, and the one mechanism behind them

| solid | face | kind | K | cells | leak | max distance from the boundary | coincident neighbour? |
| --- | --- | --- | ---: | ---: | ---: | ---: | :-- |
| `ST0923290_011` | f#98 | cylinder | 6 | 3 | 284 / 478 | **4.58e-01 cm** | **yes** (also 1 hole) |
| `ST0923290_011` | f#134 | cylinder | 2 | 2 | 482 / 996 | 1.02e-01 cm | **yes** |
| `ST0923290_011` | f#109 | cone | 1 | 2 | 185 / 185 | 1.62e-03 cm | **yes** |
| `ST0923290_011` | f#112 | cylinder | 2 | 1 | 6 / 1353 | 1.27e-05 cm | **yes** |
| `ST0923290_011` | f#111 | cone | 2 | 1 | 1 / 1332 | 1.56e-06 cm | **yes** |
| `ST0923290_016` | f#16 | cone | 4 | 3 | 65 / 1781 | **8.74e-01 cm** | **yes** |
| `ST0923290_016` | f#18 | cone | 4 | 3 | 105 / 1688 | **8.34e-01 cm** | **yes** |

**Seven of seven carry a neighbour lying on their own surface.** Two of them (f#111, f#112) leak by
1.6e-06 and 1.3e-05 cm and are tolerance artefacts by any reading. The other five are real, and the
cause is §4's: a non-iso patch seam that no co-surface trim can express and the (φ, h) box does not
cover.

**Prediction, not a measurement:** merging adjacent NURBS patches recognised as the same analytic
surface into one face before trimming would remove the coincident neighbour entirely and, on this
evidence, take R2 from 13/15 to 15/15. Nothing here tests it. It is a converter-side operation of
its own and it would also shrink the recognised-face count substantially.

### 5.4 The (u, v) box is load-bearing — R0 alone over-accepts by centimetres

The far-field sampling puts a grid over the face's whole analytic surface — the full 2π of a
cylinder's or cone's azimuth, **both nappes** of a cone, twice the face's axial span — asks R0, and
then checks each accepted point against OCCT. On the ten ladder fixtures, where R1 is perfect:

```
   18 of 56 faces accept points that are not on the face
   720 spurious points, up to 1.994 cm from the boundary
   e.g. cyl_inter_cyl f#0 : K=1, 80 of 256 far-field points accepted, 56 of them off the face,
                            max 1.893 cm
        box_union_box f#1 : K=3, 96 of 256 accepted, 32 off the face, max 1.131 cm
```

`box_union_box`'s planar faces are bounded by four neighbours of which **one is coplanar with the
face** (the union seam), leaving three half-planes to bound a rectangle. `cyl_inter_cyl` f#0 and
f#2 have **K = 1**. Neither is a defect of the corpus; both are what happens when a boundary edge
has no co-surface representation, and both are cured by the box. **R0 is not a candidate design.**

### 5.5 The instrument can prove a face wrong; it cannot prove a face right

Sweeping the grid over 16 / 24 / 32 / 48 / 64 on `ST0923290_002`, `_011` and `_016` (276 faces):

```
   cellsIn identical across all five densities   266 / 276  (96.4%)
   leak verdict identical                        270 / 276  (97.8%)
   R1 verdict identical                          270 / 276  (97.8%)
```

Ten faces change their cell count with density (`ST0923290_011` f#12: 17 / 17 / 19 / 23 / 23) and six
change their **leak verdict** — `ST0923290_011` f#109 reads leak-free at grid ≤ 32 and leaking at
grid 48 and 64. Over the 15 solids, R1-exact faces fall 604 → 589 and R2-exact 665 → 649 between
grid 32 and grid 64.

This has a hard consequence and it is the most important caveat in this document:

> **A cell set derived by sampling is an approximation with an unbounded error.** A lobe the
> sampling missed becomes a false negative in production — a hole in the solid, in the shadow of
> which parity containment is undefined. R2 as a *sampled* construction does not meet this
> project's standing bargain (*exact, or tessellated, never quietly approximate*).

What R2 needs to be admissible is an **exact** cell enumeration: intersect the face's surface with
each of the K trimming surfaces (OCCT does this analytically for quadric pairs), build the
arrangement on the surface, and read the cells the face occupies off the topology the B-rep already
carries. That is real work and it is *not* what "trim by half-spaces" promised. It is, however,
bounded and converter-side — see §9.

The asymmetry is why every number above is quoted at two densities: a face the probe calls **wrong**
is wrong, because the disagreement is a located point that OCCT independently confirms is off the
face. A face the probe calls **right** is only right at the density measured.

### 5.6 The controls, per solid

| model | solid | faces | R1 exact | R2 exact | worst depth |
| --- | --- | ---: | ---: | ---: | ---: |
| fixtures | all ten | 56 | **56** | 56 | 0 |
| `as1` | `bolt`, `nut`, `plate`, `rod` | 37 | 37 | 37 | 0 |
| `as1` | **`l-bracket`** | 16 | 14 | 16 | **4.91 cm** |
| Bagger | `BasePin` | 3 | 3 | 3 | 0 |
| Bagger | `Boom` | 31 | 25 | 29 | **15.14 cm** |
| Bagger | `Bucket` | 97 | 67 | 85 | 6.74 cm |
| Bagger | `Base` | 44 | 30 | 32 | 3.28 cm |
| Bagger | `Stick` | 24 | 18 | 22 | 6.62 cm |
| Bagger | `BucketLink1` | 22 | 19 | 19 | 0.65 cm |
| Bagger | `BucketLink2` | 23 | 15 | 15 | 0.24 cm |
| Bagger | the six rams | 44 | 38 | 39 | 0.85 cm |

`l-bracket` is the textbook failure and it is in the corpus: an **L-shaped planar face** whose region
is the union of two cells, which no intersection of half-spaces can be. `Boom` f#0 and f#20 are
planar faces with K = 15 whose region spans 8 cells; the conjunction misses 171 interior samples by
up to 15.1 cm.

R2 is never worse than R1 per face, and cannot be: `cellsIn == 1` with no false positive *is* zero
leak. But it barely helps on the six rams (38 → 39) or on `BucketLink1`/`BucketLink2` (no change at
all), because their cylindrical faces are split by a patch seam whose implicit function is
identically zero — §4's mechanism again, and the reason the cell set has to be built from topology
rather than from samples.

---

## 6. The design

Everything below is a proposal. None of it is built.

### 6.1 What is stored — sidecar v4, one appended per-surface block

The format is versioned and has been extended twice by strict append (v1→v2 header tail, v2→v3
header tail plus a per-record tail), so the mechanism exists. The natural slot is a new per-face
block **after** the v3 edge-ref block; putting it inside the wire block would fight
`wireToCurves`'s endpoint-and-join contract (`O2SurfaceSolidIO.cxx:259`), and an implicit trim has
no endpoints.

```
  uint32   nTrimSurfaces  K
    per trimming surface:
      uint32   surfaceType     1 plane · 2 cylinder · 3 cone · 4 sphere · 5 torus
      uint32   nParams ; float64 params[nParams]      the canonical face-record layout, verbatim
      int32    sense           +1 : f >= 0 is inside · -1 : f <= 0 is inside
      uint32   nEdgeIds ; uint32 edgeId[]             which MODEL edges this surface carries
  float64    trimBand          the acceptance band, in cm
  uint32     nCells            0 = pure conjunction (R1) ; otherwise a DNF (R2)
    per cell: uint32 mask      bit k set = "sense_k * f_k >= 0" required.  K <= 32.
  float64    verifiedDensity           samples per patch diagonal the converter checked at
  uint32     verifiedDisagreements     and what it found.  Must be 0 to emit.
```

`kFlagImplicitTrim`, one of the 31 unused bits of the existing per-surface `flags` word
(`O2SurfaceSolidIO.cxx:444`), selects it; `nWires` is then 0. `kSidecarVersionMax` (`O2SurfaceSolidIO.cxx:48`) goes 3 → 4.

Three fields are not decoration.

- **`edgeId[]` is what makes closure work** (§6.5). It is the same `uint32` identity the v3 block
  already uses, so the two blocks agree by construction.
- **`trimBand` is a length in cm, and that is the design's single largest simplification.** Every one
  of the five implicit functions below is a first-order signed distance, so `|f|` *is* a distance.
  Today `CurveWire::boundaryBand` (`BoundedSurface.h:2164`) has to convert `kTolerance` into
  parametric units through `ParametricMetric::maxScale` — the exact class of defect
  `Stream_E_Scale.md` and `Stream_M_Quartic.md` document. Here there is nothing to convert. The band
  must be the **model's own declared tolerance**, up to 8.1e-05 cm on ALICE3 per `Stream_O` §8.3,
  not `kTolerance`'s 1e-09 cm.
- **`verifiedDensity` / `verifiedDisagreements` are the acceptance evidence**, and §5.5 is why they
  are mandatory: the converter must run this probe's own comparison against `FClass2d` on the face
  it is about to emit, and emit only on zero disagreements at a recorded density. That converts a
  global gamble into a per-face, self-checked decision — the same discipline as the CSG emitter's
  `dV_sym` test and the recogniser's measured-gap test.

### 6.2 What the kernel evaluates

```cpp
bool implicitTrimContains(const Vec3& p, bool* boundary) const {
  uint32_t mask = 0;
  double   worst = std::numeric_limits<double>::max();
  for (int k = 0; k < mK; ++k) {
    const double f = mSense[k] * implicitValue(mTrim[k], p);
    if (f >= 0.) mask |= 1u << k;
    worst = std::min(worst, std::abs(f));
    if (mNCells == 0 && f < -mBand) { if (boundary) *boundary = false; return false; }
  }
  if (boundary) *boundary = (worst <= mBand);
  return mNCells == 0 ? true : cellSetContains(mask);
}
```

`implicitValue` is five closed forms the kernel **already writes down**, today only along a ray:

| kind | value | cost |
| --- | --- | --- |
| plane | `dot(p − o, n)` | 3 sub, 3 mul, 2 add — this is `PlanarBoundedSurface::planeDistance`, `BoundedSurface.h:3006`, already present and already correct, merely never asked about a *neighbour* |
| sphere | `‖p − c‖ − R` | + 1 `sqrt` |
| cylinder | `‖perp(p − o)‖ − R` | 1 dot, 3 fma, 1 `sqrt` — the same `x² + y² − r²` as `quadraticC`, `BoundedSurface.h:3603` |
| cone | `r·cos a − h·sin a` | 1 dot, 1 `sqrt`, 2 mul |
| torus | `hypot(ρ − R, z) − r` | 2 `sqrt` — the implicit form the quartic is built from, `BoundedSurface.h:4746` |

**No transcendental. No division. No branch per curve. No cached polyline.**

`cellSetContains` is a linear scan of a sorted `uint32` array; measured `nCells` is 1 at the median
and 19 at the maximum on the target solids, so it is a handful of compares. `K ≤ 32` is a real
constraint, and **14 faces of the whole 8170-face corpus exceed it** — `ST1829909_01` f#101 (K = 291,
392 boundary edges), `ST2487457_01` f#984/f#985 (K = 211), f#12/f#18 (K = 142), `ST1829909_01`
f#778/f#779 (K = 141), `ST2487459_01` f#3/f#168 (K = 53) and six more. All are big planar, conical or
cylindrical faces on solids that are out of scope for other reasons; none is on the 15 target solids,
whose maximum is K = 33 and whose *blocking* faces are all K ≤ 6. A face over the limit declines the
implicit trim and keeps its wire.

### 6.3 The (u, v) box is not a special case — it is four more half-spaces

The measurement in §5.4 says the box is required. It needs no new machinery:

- `heightMin` / `heightMax` on a cylinder or cone ≡ two planes perpendicular to the axis;
- `thetaMin` / `thetaMax` on a sphere ≡ two planes perpendicular to the polar axis;
- `phiStart` / `phiSweep` with sweep ≤ π ≡ two planes containing the axis; with sweep > π it is the
  *union* of the two complements, which is one extra cell in R2's vocabulary and needs nothing new.

So the box costs at most four plane evaluations — 20 flops — and it disappears into the same loop.
This also removes the last reason to form `(φ, h)` from the hit point at all.

### 6.4 Ray-surface acceptance — structurally unchanged, and one transcendental cheaper

In `appendIntersections`, `pointInTrim(hitPhi, hitHeight, &onTrimBoundary)` becomes
`implicitTrimContains(hitPoint, &onTrimBoundary)`. Nothing above the surface class changes:
`crossingSense`, `forEachCrossingCluster`, `oddCrossingParity`, `nearestCrossing`, the BVH and every
TGeo entry point see the same `RayHit{distance, normal, onTrimBoundary}`.

What *disappears* is the parameter formation that exists only to feed the current test — 1 `atan2`
plus 1 `hypot` per root on a cylinder or cone, 1 `acos` plus 1 `atan2` plus 1 `hypot` on a sphere,
2 `atan2` plus two `unwrapAngleInto` on a torus. On aarch64 one `atan2` is worth several `sqrt`s, so
for a cylinder with the measured median K = 2 the whole implicit test plausibly costs less than the
coordinates the current test needs before it starts.

### 6.5 Closure and edge identity — a *stronger* check, and one real casualty

**Unchanged and still the verdict:** `applyEdgeIdentityClosure` (`BoundedSurface.h:5591`). Every
`edgeId` must appear exactly twice with opposite sense, zero tolerance, zero sampling. The v4 block
carries `edgeId[]` per trimming surface, so its precondition — every surface names its edges — is
met.

**New, and it is what the brief predicted:** the two faces sharing an edge would be trimmed by
*literally the same* implicit surface. Face A's own surface record and face B's trimming-surface
record for that edge are the same recognition of the same `TopoDS_Face`, so they are the same
doubles. That makes a new, exact, tolerance-free load-time check available:

> **Co-surface reciprocity** — for every `edgeId` shared by faces A and B and carried implicitly on
> both sides, A's trimming-surface record for that edge must equal B's own surface record, and B's
> must equal A's, bit for bit.

A converter that violates it has made two faces disagree about a shared edge — the exact defect
`Stream_F_EdgeIdentity.md` exists to remove — and it would be caught at load, not measured
afterwards. Correspondingly `measureSharedEdgeDeviation` (`BoundedSurface.h:5226`) returns
**structurally zero** for such an edge: the two faces describe it by the same pair of equations
`{f_A = 0} ∩ {f_B = 0}`, so there is nothing left to sample. That also retires, for these faces,
`Stream_Q_EllipseTrim.md` §6's polyline-flattening artefact, which is a property of *sampling* two
realisations of one curve at different phases.

**The casualty:** `sampleTrimCurve` (`BoundedSurface.h:2822`), `appendDirectedEdges` (2907) and
`appendRims` (2915) all demand an ordered boundary polyline per edge, and an implicit trim has no
natural ordering. Three consumers, in decreasing importance: the rim measurement (still reported),
the directed half-edge chord count (already superseded as a verdict — the header says so just below `validateClosure`),
and `measureSharedEdgeDeviation` (structurally zero, above). The cheapest honest route is to keep
writing the boundary polyline into v4 as a **reported side-channel that no containment test may
consult** — the converter already has it, it is what `_quadric_trim_wire` produces, and it costs
bytes rather than correctness. If it is ever re-consulted for containment the exactness claim
reverts, so it should be stored under a name that says so.

### 6.6 Coexistence — per face, and it must be per face

A face carries either a wire trim or an implicit trim, chosen by the flag. Both kinds coexist inside
one solid, inside one BVH, with no interaction: the trim test is entirely local to the surface
object. §0's control table is why this is not optional — on Bagger the implicit trim is *worse* than
what ships on 12 of 13 parts, and a wholesale replacement would destroy a gated corpus.

The converter's rule should be: attempt the parametric wire first (it is what is proven); attempt the
implicit trim only where the wire declines; and emit the implicit trim only when the face's own
verification is clean at a recorded density. `surface_report.json` records the choice and the
evidence per face.

---

## 7. The cost model

### 7.1 Trimming surfaces per face

| population | p50 | p90 | p99 | max |
| --- | ---: | ---: | ---: | ---: |
| the 15 `Stream_O` solids, recognised faces | **2** | 4 | 6 | **6** |
| the 15 `Stream_O` solids, all faces | **2** | 4 | 8 | 33 |
| all corpora, P1 | 4 | 6 | 13 | 291 |
| boundary edges with a neighbour, the 15 solids | 4 | 8 | 18.7 | 59 |

K is bounded by the number of distinct neighbour *surfaces*, which is always ≤ the number of boundary
edges and on this corpus is typically **half** of it — one cap plane bounds four edges of a NURBS
patch quartet. Cells: p50 **1**, p90 3, p99 8, max 19 on the target solids.

### 7.2 Against what ships

`CurveWire::classify` (`BoundedSurface.h:2174`) is **two full passes** over the wire's N curves with
no bounding-box prefilter, no spatial index and **no early exit in the crossing pass**. Per curve:
a line is ~15–20 flops plus a division; an arc is 1 `atan2` plus up to 8 break parameters, a
`std::sort`, a `sqrt` and 2 trig per monotone sub-arc; a B-spline walks its cached polyline of M
chords, twice. `curveTrimContains` short-circuits only when the outer wire says `Outside`.

| | today, median target face | implicit, median target face |
| --- | --- | --- |
| trim test | 2 passes × ~4 curves = ~8 curve visits, arcs with `atan2`/`sin`/`cos` | **2** signed evaluations, ~2 `sqrt`, no transcendental |
| parameter formation per root | 1 `atan2` + 1 `hypot` (cyl/cone), 1 `acos` + 1 `atan2` (sphere), 2 `atan2` (torus) | **none** |
| early exit | crossing pass: none | first violated half-space, typically the first one for a ray outside the patch |
| p99 face | ~38 curve visits | 8 evaluations |
| worst face in the corpus | `ST1829909_01` f#101: 392 curves × 2 = 784 visits | K = 291 with early exit — one of the 14 faces over the 32-surface limit, declines |

**The implicit trim is cheaper, not 10× slower, and the reason is structural rather than
constant-factor:** it replaces an O(N) two-pass scan with transcendentals by an O(K) single pass with
early exit and no transcendentals, and it deletes the coordinate transform that today precedes the
test. The place a reviewer should look for a regression is not throughput but the **band**: a
co-surface trim's band is the model's declared tolerance (up to 8.1e-05 cm), four orders looser than
`kTolerance` and one order looser than `kBSplineFlatness`, so `onTrimBoundary` will fire on many more
hits and `containsByVote`'s ambiguity accounting will see more of them. That is a real cost and it is
unmeasured.

### 7.3 Size

A trimming surface is one type word, ≤ 15 doubles, a sense and its edge ids — about 140 bytes. It
replaces, for the edges it carries, `Curve2D` records of 4 doubles (line) to several kilobytes (the
179-pole rational B-spline of `Stream_F`'s reproducer). ALICE3's sidecar is 3.6 MB after canonical
recognition; implicit trims would shrink it further, not grow it.

### 7.4 Tangency — the conditioning cost, quantified

99 of the 678 faces on the target solids (14.6%) have a trimming surface **tangent** to the face
along a shared edge — a milled flat tangent to a fillet is the common case, and it is measured as
`min |sin(angle between the two surfaces' normals)|` over samples that genuinely lie on the trimming
surface. At a tangency `|f|` grows as `ε²/(2R)` in the arc distance ε rather than as ε, so a band `b`
resolves the boundary only to `ε = sqrt(2Rb)`: with the declared tolerance `b = 8.1e-05 cm` and a
1 cm radius that is **130 µm** of positional uncertainty in where the boundary sits, against a trim
whose own accuracy is 0.8 µm.

**Tangency is not the mechanism behind any observed failure.** Over all 678 faces of the 15 solids at
grid 64, R1 is wrong on 89; **85** of them fail through non-convexity (`cellsIn > 1`) and the
remaining **4** through a coincident neighbour (`ST0923290_011` f#10, f#15, f#111, f#112 — every one
of them). Eight of the 99 tangent faces are among the 89, but none of the 89 is unaccounted for by
those two mechanisms. So the tangency figure above is a quantified *risk*, not an observed defect,
and it is a reason the band must be per-face and recorded rather than a global constant.

---

## 8. Being adversarial about this result

**The conjunction's failures are false *negatives*, mostly.** 191042 of the 205027 disagreements on
P2 are the rule saying OUT where the truth says IN. A false negative is a hole in a face — a ray
passes through material — and in a parity containment scheme it is exactly as bad as a false
positive, so nothing is saved by the asymmetry. It does mean the naive rule fails *conservatively* in
the sense that it never claims material that is not there.

**The 15/16 rollup in `Stream_O` was never a coverage claim and this does not make it one.** It was a
converter-side edge count. This document says what happens to those edges when they become a
containment rule, and the answer for the naive version is 1 of 15.

**`ST0923290_021` is untouched.** It is `Stream_O`'s 16th solid, it has zero rejected edges, and it
is declined by `recognized quadric trim wire has fewer than 3 edges`. Measured here it is
perfect — 33 faces, 0 false positives, 0 false negatives, every region a single cell — which is a
reminder that it never needed a trim representation at all.

**The cell counts are lower bounds.** A cell present in the true region but missed by the grid
lowers `cellsIn`; §5.5 shows this happening on 10 of 276 faces between grid 16 and 64. Every
`cellsIn` in this document is therefore "at least", and every "R2 carries this face" is "at the
density stated".

**The far-field ground truth is weaker than the near-field one.** It projects onto the face's
surface with `ShapeAnalysis_Surface` and tests the residual before calling `FClass2d`, so a
projection landing in a local minimum would mis-report. It is used only for §5.4's qualitative claim
that the box is required, where the effect is 1.97 cm and the mechanism (`K = 1`; three half-planes
for a rectangle) is independently visible in the record.

**One thing was not tried and might change the verdict for R1.** Nothing here attempts to *choose*
the trimming set. The probe takes every distinct neighbour surface; a converter is free to add
half-spaces that are not neighbours at all — a separating plane through two boundary vertices, say —
and turn a multi-cell region into a single cell. That is a synthesis problem (find the smallest set
of half-spaces whose intersection is the face) and it is genuinely open. It would, however, break
the property that makes this attractive: a synthesised plane is not a neighbour's surface, so
co-surface reciprocity (§6.5) no longer holds for the edge it bounds.

---

## 9. What this licenses and does not license

### It licenses

1. **Retiring R1 — the naive conjunction — as the design.** It carries 1 of the 15 target solids and
   destroys 12 of Bagger's 13. It should not be built.
2. **Building the per-face, converter-verified implicit trim with a cell set (R2)** for the faces the
   parametric wire declines, *provided* the cell set is enumerated exactly rather than sampled.
   Scope, measured: 101 blocking faces on 15 solids; K ≤ 6 on all of them; cells ≤ 8; seven
   surface-pair kinds.
3. **The exact cell enumeration as its own work item.** OCCT intersects quadric pairs analytically;
   the face's own wire topology already says which arcs bound which lobe. This is the piece that
   turns R2 from an approximation into a representation, and it is the gate on everything else.
4. **Patch merging as a prerequisite worth measuring.** 56% of the target faces have a neighbour on
   their own surface, and all 7 of R2's residual failures do. Merging same-surface NURBS patches
   would remove the mechanism outright.
5. **Co-surface reciprocity as a load-time closure check** (§6.5), independent of whether the trim
   itself is implicit — it is a converter-consistency assertion that costs nothing.

### It does not license

1. **Any statement about the kernel, transport or navigation.** Nothing here has been through
   `O2BVHSurfaceSolid`, the oracle gate or the X-ray benchmark. No sidecar was written. The cost
   model in §7 is an operation count, not a measurement.
2. **"Implicit trims give exact geometry."** They give geometry as exact as the source B-rep, exactly
   as `Stream_O` §9 says — up to 8.1e-05 cm — *and* only as exact as the cell enumeration, which is
   currently a sampling.
3. **Replacing the parametric trim.** §5.6. On three gated corpora the implicit trim is worse on
   most parts.
4. **A coverage number.** `Stream_O` §7's lesson applies unchanged: coverage is a property of how the
   assembly decomposes into leaf solids. Every figure here names ALICE3's 55-solid decomposition.
5. **Anything about torus trims.** 208 torus-trimmed-by-torus and 539 torus-trimmed-by-plane
   surface pairs appear in the corpus and the implicit function is written down, but no torus face is
   in the 15 target solids and none was scrutinised individually.

### The order I would take it in

1. **The exact cell enumeration** (§5.5, §9.3). Everything else is an approximation until it exists.
   Its acceptance test already exists: this probe, at rising density, must find zero disagreements.
2. **Patch merging** (§9.4), measured first as a census — how many recognised faces merge, and does
   the coincident-neighbour count go to zero.
3. **The v4 record and the kernel evaluation** (§6.1, §6.2), gated per face on
   `verifiedDisagreements == 0`.
4. **Co-surface reciprocity in the loader** (§6.5). Small, independent, and it is the one part of
   this design that strengthens an existing guarantee rather than adding a new representation.

---

## 10. What could not be measured

- **Anything downstream of the converter.** No sidecar, no gate, no X-ray. The 15 solids still
  produce nothing that has ever been loaded.
- **Whether the exact cell enumeration is tractable** for the seven surface-pair kinds. §9.3 asserts
  OCCT can intersect them; nothing here ran that.
- **The `onTrimBoundary` cost** of a band four orders looser than `kTolerance` (§7.2). It will change
  how often `containsByVote` reports an ambiguous shot and nothing here predicts by how much.
- **`oTOF System V3-R92cm.step`**, the fourth corpus, and IRIS. IRIS would in any case gain nothing:
  `Stream_O` §7 shows its decomposition puts every recognised face inside one 1734-face solid.
- **Faces with K > 32.** Fourteen in the corpus; all declined by the proposed limit, none examined
  individually.
- **Whether a synthesised (non-neighbour) half-space could rescue R1** (§8, last paragraph).

---

## 11. Claims in other documents that this work bears on

I have edited none of these; they are for central reconciliation.

- **[`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md) §9's licence (2)** — "investigating an
  implicit / co-surface trim representation as the replacement, with the measured scope in hand" —
  stands, but its §10 open question ("whether an implicit trim can be *evaluated* cheaply in the
  kernel … nor about how such a trim orders its boundary, nor about how `CloseShape` would count its
  edges") is now answered: cheaply, yes (§7); its boundary has no order and three consumers want one
  (§6.5); closure is unchanged because it already decides on identity, and the geometric measurement
  becomes structurally zero.
- **`Stream_O` §7's rollup, "15 of 16 fully covered"**, is a converter-side edge count and remains
  correct as such. It must not be read as "15 solids would emit": with the naive conjunction the
  number is **1**, and with a sampled cell set **13**.
- **[`Stream_F_EdgeIdentity.md`](Stream_F_EdgeIdentity.md) §8's caveat** — identity certifies that
  the topology survived, not each face's geometry — gets a partial answer for co-surface-trimmed
  faces: reciprocity (§6.5) is a geometric check that identity currently lacks, alongside
  `Stream_L_ALICE3Defect.md`'s outward-normal criterion.
- **[`ExactTrimTopology.md`](ExactTrimTopology.md)'s note** that "the exact intersection curve of two
  quadrics is a degree-4 algebraic curve and is representable exactly in *neither* face's parametric
  domain with the current `Curve2D` vocabulary" is the premise this whole stream acts on, and it is
  correct. What this document adds is that not representing the curve does not by itself give a
  containment test.
- **`Tutorial.md` §5.1's list of what is missing for completeness** does not mention the trim
  vocabulary as a coverage lever at all; `Stream_O` §5.2b and this document both say it is the
  lever. That is a `Tutorial.md` gap, not an error.
