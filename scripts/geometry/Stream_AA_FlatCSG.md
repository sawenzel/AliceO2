# Stream AA — a flat TGeoBVHCSG: the cell-count table, and the verdict

Date: 2026-08-22. Branch `swenzel/bvhsurfacesolid`. This answers the commissioned question
(`Roadmap.md` (g)):

> *"a subagent should review if we've already done everything we could to make a CSG conversion
> (via a possible new TGeoBVHCSG ... without recursion ... flat) possible? I'd like to understand
> if this could have its own place next to BVHSurfaceSolid or if it is strictly not worth it."*

Instruments: `probes/cellCountProbe.py` (new; 7 self-checks, all passing — it runs the
Shapiro-style split loop with `BRepAlgoAPI_Splitter` and counts the pieces) and the per-part
decline catalogue (`csg/decline_catalogue.py` → `website_data/decline_reasons.json`), both built
this session. `CSG_Pipeline.md` §8 step 7 gated the flat solid on "the cell-count and
success-rate table over the eligible ALICE3 and Bagger solids"; that table did not exist. It now
does (§2), and it was measured, not estimated.

---

## 0. The answer first

**A flat, non-recursive, BVH-indexed CSG solid has a legitimate place next to
`O2BVHSurfaceSolid` — but not yet, and not as the next step.** The measurements below support a
sharper statement than either "build it" or "strictly not worth it":

1. **The decomposition it needs is tractable — a refuted fear.** Every CSG-declining Bagger part
   and 15 of the 16 ALICE3 solids that today ship tessellated despite being analytic decompose
   completely, in under a second each, into **3–37 cells** (median ~11), with volume conserved
   to ≤1.3e-6 relative and zero splitter failures on the completed parts. The standing
   projection ("a cell count in the tens to hundreds per body", `Stream_A_CSG.md` §3 Tier 3) is
   refuted on this corpus for everything under ~250 faces, and confirmed only for the four
   ST1829909 giants (§2.3).
2. **What it buys is coverage and exactness, not speed of anything that ships today.** Be
   precise about the motivation: `Stream_P` measured `TGeoCompositeShape` **40–145× faster**
   than the surface solid on the shallow composites the corpus actually has, so speed of
   existing CSG parts needs no new class. The prize is the **16 ALICE3 solids** (29% of the
   model) that are analytic after recognition and ship as meshes because their trims are
   transcendental in the (φ, h) chart — plus `cyl_inter_cyl` and `tube_window`, the last
   fixture-gate residual, which are **single cells of 2 and 4 halfspaces** (§3.1).
3. **The same cells serve the other open research route.** A union of arrangement cells *is* a
   flat recursion-free CSG: `Stream_R`'s 13/15 result and this decomposition are two shadows of
   one object, the arrangement of the part's own carriers, and the splitter just enumerated it
   in 3D for the very solids Stream_R could not finish in 2D (§3.4). One enumeration, built once
   at converter time, can be spent either as implicit trims for the surface solid or as a flat
   solid.
4. **A strictly cheaper first step exists and should come first**: Tier-0 canonicalisation
   (measured prerequisite — every observed splitter failure is a NURBS-encoded carrier the
   probe could not extend), then emission of the decomposition as plain balanced
   `TGeoCompositeShape` trees, gated by the existing symmetric-difference + oracle acceptance.
   That ships the coverage with **no new C++ class**. Build `TGeoBVHCSG` only if the measured
   cost of those trees on the X-ray benchmark fails the criteria in §5 — the numbers say it
   plausibly will for the torus-family parts, but that is exactly the measurement that must
   decide, not this document.

---

## 1. The decline catalogue this assessment stands on (deliverable B)

Every number in §2 is about parts that *decline* one representation or the other, so the first
deliverable was making the converter say, per part, why. That is now shipped, not just written
here:

- `csg_report.json` rows carry a top-level machine-readable **`whyNotCSG`** (None for a part
  that ships as CSG). The free-form and toroidal declines carry counts ("toroidal faces: 2 of
  97"), the structural declines carry the cluster structure ("7 axis clusters: beyond the
  recogniser's scope [44 faces: 16 cylinder, 28 plane]"), and the acceptance rejections carry
  the measured symmetric difference.
- `surface_report.json` volumes carry **`why_not_surface`**, distilled from the per-face
  classification and extraction reasons ("3 face(s): bspline face extraction not implemented
  yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not
  axis-aligned in (phi, h/theta)").
- The two defects `NEXT.md` item 2 records are fixed: a `--csg auto` run whose interpreter has
  no PyROOT now warns **loudly per part** and records `"csg deferred: ROOT unavailable"` instead
  of the empty `declined CSG: None`.
- `csg/decline_catalogue.py` joins the two reports into
  `website_data/decline_reasons.json` — **78 parts** (10 fixtures + 13 Bagger + 55 ALICE3),
  schema `{"parts":[{name, model, shipsAs, whyNotCSG, whyNotSurface}]}`, **0 rows with a missing
  reason**. The gates did not move: Bagger 13/13 and fixtures with the same two known distout
  sliver rays; a field-by-field `gate.json` diff against a pre-change run shows only scratch-dir
  paths and the two intended reason strings; `O2_CADtoTGeo.py --self-test` 48,
  `csg/emit.py --self-test` 33/33, `probes/cellCountProbe.py --self-test` 7/7.

The catalogue's shape, in one paragraph: of the 78 parts, 11 ship as CSG (with `whyNotCSG` null
and the dV_sym evidence attached), 32 as surface solids, 35 as meshes — all 35 in ALICE3. Of
those 35 mesh parts, 19 carry genuinely free-form surfaces (Stream B's remit, no constructive
representation can reach them) and **16 are the analytic-after-recognition solids blocked only
by their trims**. Those 16 are the coverage question this stream was asked about, and they are
all in the ST0923290 (torus + NURBS-encoded-quadric) family.

---

## 2. The cell-count table (deliverable A1)

Method: split at the sharpest *trusted* concave/mixed edge (the `Stream_A` §1.4 near-tangential
filter, sin ≥ 1e-3), extending the carrier of one adjacent face — preferring the plane — to a
full surface; recurse; a piece with no trusted concave edge is one cell. Budgets (128 cells, 256
splits, 300 s) turn a blow-up into a reported lower bound. `halfsp` is distinct carriers, not
faces (CAD splits one carrier into 1.5–3 faces, `Stream_A` §3); `Σcell-h` is the summed per-cell
distinct-halfspace count — the total AND-length of the flat DNF, i.e. its size. Where a part
still carries stored-NURBS faces the halfspace columns fall back to face counts and therefore
**overstate by the 25–35%** the census measured.

### 2.1 Bagger — every CSG decliner (all ship as exact surface solids today)

| part | faces | halfsp | clusters | concave (trusted) | **cells** | splits | failures | Σcell-h | max cell | vol conserved |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| `Base` | 44 | 26 | 7 | 16 | **10** | 8 | 0 | 79 | 10 | Y |
| `Boom` | 31 | 20 | 6 | 6 | **7** | 6 | 0 | 59 | 13 | Y |
| `Stick` | 24 | 19 | 6 | 8 | **7** | 6 | 0 | 56 | 10 | Y |
| `BucketLink1` | 22 | 12 | 2 | 10 | **7** | 4 | 0 | 36 | 8 | Y |
| `BucketLink2` | 23 | 14 | 3 | 18 | **11** | 6 | 0 | 45 | 7 | Y |
| `Bucket` | 97 | 53 | 13 | 68 | **28** | 14 | 0 | 185 | 10 | Y |

All six complete in ≤0.44 s, zero failures, volume conserved to ≤1e-6 relative, zero residual
untrusted concave edges inside any cell. `Bucket` — 97 faces, 68 concave edges, tori and
spheres — was the emblem of "Tier 3 territory", and it is 28 cells.

### 2.2 ALICE3 — the quadric-only CSG decliners (all ship as surface solids today)

| part | faces | halfsp | concave (trusted) | **cells** | splits | failures | Σcell-h | max cell | vol conserved (rel) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `ST2487455_01` | 66 | 44 | 25 | **22** | 17 | 0 | 122 | 12 | 1.4e-05 |
| `ST2487455_002` | 6 | 4 | 0 | **1** | 0 | 0 | 4 | 4 | 0 |
| `ST1A38495_01` | 65 | 43 | 9 | **13** | 8 | 0 | 108 | 17 | 1.4e-07 |
| `ST1A38526_01` | 53 | 35 | 9 | **10** | 9 | 0 | 88 | 14 | 1.4e-08 |
| `ST1A38494_004` | 10 | 10 | 0 | **1** | 0 | 0 | 10 | 10 | 0 |
| `ST2487459_01` | 202 | 167 | 6 | **10** | 3 | 0 | 197 | **164** | 1.1e-04 |
| `ST2487721_01` | 32 | 14 | 4 | **3** | 2 | 0 | 20 | 8 | 7.3e-08 |
| `ST2487462_01` | 80 | 46 | 4 | **3** | 2 | 0 | 51 | 41 | 6.4e-08 |
| `ST1782525_01` | 50 | 27 | 28 | **15** | 14 | 0 | 90 | 6 | 2.5e-13 |

### 2.3 ALICE3 — the four ST1829909 giants: the projection is confirmed *here* and only here

| part | faces | halfsp | concave (trusted) | cells | failures | residual untrusted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ST1829909_01` | 1052 | 564 | 484 | **≥129** (budget) | 30 | 36 |
| `ST1829909_002` | 965 | 484 | 904 | **≥129** (budget) | 8 | 21 |
| `ST1829909_003` | 236 | 128 | 207 | **≥130** (budget) | 19 | 55 |
| `ST1829909_004` | 720 | 325 | 443 | **≥129** (budget) | 0 | 21 |

These blow the 128-cell budget within seconds and would keep going. They are exactly the bodies
`CSG_Pipeline.md` §6 said the BVH surface solid is uniquely good at — **and all four already
ship as exact surface solids**, so nothing is lost by excluding them. The ≤40-face-and-≤12-
trusted-concave gate `Stream_A` proposed is, on this evidence, too tight (Bucket at 97/68 and
ST1782525_01 at 50/28 decompose in under half a second); a workable gate is a *cell budget*, not
a face budget — stop at ~64 cells and fall back.

### 2.4 ALICE3 — the ST0923290 family: **the 16 tessellated-but-analytic solids, measured**

This is the population a flat CSG would actually add coverage on (§1). All 16, plus the rest of
their family, went through the probe. Halfspace columns are face-count fallbacks here (stored
NURBS), so read them as upper bounds.

| part (16 = in the target set) | faces | concave (trusted) | **cells** | failures | Σcell-h (≤) | vol conserved (rel) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 `ST0923290_002/003/004/005` | 24 | 23 | **11** each | 0 | 82 | 3.4e-08 |
| 16 `ST0923290_011` | 193 | 72 | **≥10** | **4** | 195 | 1.8e-07 |
| 16 `ST0923290_014` | 138 | 68 | **6** | 0 | 374 | 6.2e-08 |
| 16 `ST0923290_016` | 59 | 20 | **8** | 0 | 85 | 3.0e-07 |
| 16 `ST0923290_017` | 24 | 23 | **13** | 0 | 95 | 1.3e-06 |
| 16 `ST0923290_021` | 33 | 17 | **7** | 0 | 61 | 6.3e-08 |
| 16 `ST0923290_022/023/024` | 24 | 23 | **11** each | 0 | 81 | 1.2e-06 |
| 16 `ST0923290_025/026` | 24 | 23 | **14** each | 0 | 100 | 1.3e-06 |
| 16 `ST0923290_027/028` | 24 | 23 | **11** each | 0 | 82 | 8.5e-08 |
| — `ST0923290_01` | 352 | 166 | ≥24 | 4 | 381 | 3.8e-07 |
| — `ST0923290_006/007/008/009` | 45 | 15 | ≥2 | 1 each | 8 | 5.9e-07 |
| — `ST0923290_010` | 86 | 60 | **37** | 0 | 336 | 8.1e-09 |
| — `ST0923290_012/013/015/018/019/020` | 10–246 | 0–43 | 1–14 | 0–1 | 26–149 | ≤9.1e-08 |

**Headline: 15 of the 16 decompose completely — 6 to 14 cells each, zero splitter failures,
volume conserved to ≤1.3e-6 relative.** The 16th (`_011`) reaches 10 cells with 4 unresolved
pieces. **Every splitter failure observed anywhere in the family has one measured cause: the
witness edge's both faces are NURBS-encoded carriers the probe cannot extend.** Those are
precisely the 1004 faces `Stream_A` §2.3 measured as quadrics in disguise at 3.4e-8 mm — so
Tier-0 canonicalisation is not merely nice-to-have for a flat CSG; it is the measured
prerequisite for enumerating its cells at all (and for writing any cell's coefficient blocks).

---

## 3. What a flat CSG would buy, quantified (deliverable A2)

### 3.1 Coverage — the honest count

- **Bagger: zero new coverage.** All 13 parts already ship exactly (7 CSG + 6 surface). A flat
  CSG changes their *representation quality* (§3.2, §3.3), not their coverage.
- **ALICE3: up to 16 solids** (55 → potentially 51 of 55 analytically exact; the remaining 19
  are genuinely free-form and belong to Stream B; the 4 giants and a handful of swept-B-spline
  parts stay with the surface solid or the mesh). 15 of the 16 are measured tractable today
  modulo Tier-0; `_011` needs Tier-0 first even to finish enumerating.
- **Fixtures: the last residual.** `cyl_inter_cyl` (1 cell, 2 halfspaces) and `tube_window`
  (1 cell, 4 halfspaces) — the two parts carrying the fixture gate's only remaining
  disagreements (2 distout sliver rays) and the ladder's most expensive surface solids
  (22.2 and 17.2 µs/ray, `Stream_P` §4). A **single-cell** emitter — the no-decomposition,
  no-risk subset `Stream_A` §3 already recommended — covers both, plus 9+12 more one-cell
  prototypes on ALICE3/IRIS.
- **oTOF: none.** 19 of 20 prototypes are Tier-1 boxes; the blocker there is XCAF traversal,
  not representation.

### 3.2 Exactness and capacity

A cell decomposition's cells are **disjoint**, which buys two things a `TGeoCompositeShape`
cannot offer: `Capacity()` is the *sum of per-cell volumes*, computable exactly by OCCT at
conversion time and asserted at load (against the composite's Monte-Carlo capacity with its
measured 1.8e-2 spread, `Stream_H` §2.2); and distance queries reduce to per-cell interval
clipping followed by a union of abutting intervals — no parity, no trims, no `CloseShape`, no
join tolerance. The entire N3/closure problem class does not exist in this representation
(`CSG_Pipeline.md` §5), which on the 16 target solids is the whole point: their *only* defect
today is the trim.

Two honesty notes. The splitter is OCCT's tolerant boolean core, not an exact arrangement
enumeration: two parts show volume drift above their model-tolerance band (`ST2487455_01`
1.4e-5, `ST2487459_01` 1.1e-4 relative), so per-part acceptance must remain the existing
symmetric-difference + oracle pair, which rejects exactly such cases. And splitting only on
*trusted* concave edges leaves blend-seam reflexes inside cells (0 on Bagger and on 14 of the
16; 1–8 on `_010`/`_011`/`_01`; 21–55 on the giants) — a cell containing one is not exactly the
intersection of its halfspaces, the error being bounded by the blend's deviation from tangency.
Also that must be caught by the symmetric difference, not assumed small.

### 3.3 Query cost — anchored on Stream_P, with derived parts labelled as such

Measured anchors (`Stream_P`): a bare `TGeoTube::Contains` is 1.8–2.5 ns, i.e. ~1 ns per
halfspace sign test; `TGeoCompositeShape` costs 5.7 ns per leaf on `Contains` and 24–36 ns per
leaf on `DistFromOutside`, exhaustively, with a 7× entry cost; the surface solid costs 135–337
ns per *candidate patch* with line/arc trims and 2.3–5.9 µs with curved-intersection trims, and
`Safety` is unaccelerated (0.8 ms on a 965-patch part).

Derived for a flat DNF (estimate, not a measurement): `Contains` without pruning is Σcell-h sign
tests — **~60–200 ns** on the 16 target solids (Σcell-h 61–100, ≤374) — and a BVH over cell
AABBs cuts the OR to the 2–4 candidate cells a point can be near, i.e. **tens of ns**; that BVH
is the substrate this branch has already built and validated twice (patch BVH, sub-patch BVH).
Distances add a root solve per lateral halfspace (quadratic; quartic for the torus — solvers
already exist and are validated in `BoundedSurface.h`). For comparison, the same parts today: as
meshes, 0.4–1.2 µs/ray with LOST rays and validity risk; as surface solids (their Bagger-class
analogues), 0.4–2.1 µs on `Contains` and 2.6–10.8 µs on `Safety`.

The plain-composite alternative, from the same anchors: emitting a 11-cell decomposition as a
balanced union of intersection-cells is ~80–100 `TGeoBoolNode` leaves, i.e. **~0.5 µs
`Contains`, ~2–3.5 µs `DistFromOutside`** by the measured linear law — roughly at parity with
the surface solids of that size and clearly better than the mesh only on correctness, not on
speed. So: **as pure speed, the flat solid beats the composite emission by an estimated order of
magnitude at ≥10 cells; but neither is needed for speed on anything that ships today** — the
corpus's real composites are 2-leaf and already 40–145× faster than the surface solid.

### 3.4 The arrangement-cell identity (Stream_R §9.3), made explicit

`Stream_R` reached 13/15 by replacing one sign per neighbouring surface with a small *set* of
admissible sign vectors — "the trimmed region is a union of cells of the arrangement" — and then
gated everything on an exact cell enumeration it judged research-grade, because sampling cannot
prove a cell set complete. **These are the same cells.** A face's 2D arrangement cells are the
boundary traces of the 3D cells of the carrier arrangement; a union of 3D arrangement cells over
oriented quadric halfspaces *is* a flat recursion-free CSG — DNF, two levels, no tree. What this
stream adds is the measurement that the enumeration Stream_R lacked is performable by
`BRepAlgoAPI_Splitter` at converter time on precisely Stream_R's target population (15/16
complete, volume-checked), with exact acceptance available afterwards. The two open research
routes — "exact trims via arrangement cells" for the surface solid and "flat CSG cells" — are
therefore one work item with two spend paths, and whichever representation consumes the cells,
the enumeration and its Tier-0 prerequisite are shared. This is the strongest structural
argument that a flat solid is *not* a parallel universe next to `O2BVHSurfaceSolid` but the same
substrate with a second leaf kind, exactly as `CSG_Pipeline.md` §9.3 framed it.

---

## 4. What it costs (deliverable A3)

- **The decomposition algorithm — no longer the hard part on this corpus.** §2 is the evidence:
  the split loop is ~150 lines of probe code, runs in seconds, and fails only where Tier-0 has
  not been done. What remains genuinely hard: the halfspace orientation walk per cell (the
  census's `halfspace_side` already computes and cross-checks it per face), degenerate/tangent
  carriers (OCCT's classic weakness — not yet hit on this corpus, but the corpus is small), and
  the near-tangential residuals of §3.2. The mitigation is structural and already the project's
  design: greedy proposal, exact acceptance, fall back one tier on any failure.
- **Tier-0 canonicalisation first** — pure converter work, measured at "half a day" scope,
  needed by both spend paths and by nothing less than every observed splitter failure.
- **A new C++ shape class + IO + gate integration** — a real stream: cell blob format, loader,
  `Contains`/distances/`Safety`/`Capacity`, the cell-AABB BVH, `TVirtualGeoConverter` for G4.
  The gate itself needs *nothing*: it scores any `TGeoShape` (Stream G), and the acceptance
  harness (`csg/accept.py`) and the X-ray bench are representation-blind. Realistic effort is
  the same order as the sub-patch BVH stream, i.e. weeks not days — and it is the natural AOT
  codegen target recorded as the preferred end state (fixed-length coefficient blocks,
  straight-line code, no virtual dispatch).
- **The composite-tree alternative costs almost nothing new** — the emitter learns to write a
  union-of-cells `TGeoCompositeShape` (balanced, per `Stream_P` §3 reading 4), `primitives.py`
  gains intersection/subtraction nodes, and both acceptance tests run unchanged. Its running
  cost is §3.3's µs-range estimate and Monte-Carlo capacity.

---

## 5. Verdict (deliverable A4)

**Not strictly worth it today; a place next to `O2BVHSurfaceSolid` is real but must be earned
through two cheaper, measurable gates.** The order:

1. **Tier-0 canonicalisation** in the converter (already recommended by three documents;
   now also the measured unblocking condition for cell enumeration). Falsifier: fewer than
   ~1000 ALICE3 faces decode at ≤1e-7 relative gap, or existing conversions move — then stop.
2. **The single-cell emitter** (`concaveEdgesTrusted == 0`, no splitter, no new class if emitted
   as a composite of halfspace-primitives — or deferred until step 3 decides the target),
   retiring `cyl_inter_cyl` and `tube_window` exactly. Falsifier: the symmetric difference or
   the oracle gate rejects a single-cell part — that would break the §2.1 premise itself.
3. **Decomposition → balanced `TGeoCompositeShape` emission for the 16**, cell budget ~64,
   accepted per part by symmetric difference + oracle gate, scored on the X-ray bench alongside
   the meshes they replace. This is the coverage prize, with zero new kernel risk.
4. **Build `TGeoBVHCSG` only on measured need**, i.e. if step 3 ships parts whose composite form
   is measurably too slow where it matters (transport ns/ray against the tessellated baseline on
   the X-ray bench; or `Safety`, which is the composite's worst measured law) — the §3.3
   estimates say the torus-family parts plausibly land there, at ~80–100 boolean leaves —
   **or** when the AOT-codegen phase begins, whose natural input the flat cell list is either
   way. If step 3's trees prove fast enough, record that as the decision and stop; that outcome
   is a success, not a failure of this analysis.

What would change this verdict in the *other* direction (build the flat solid sooner): a corpus
where multi-cell quadric bodies dominate placements (an oTOF-like mechanical assembly once XCAF
traversal is fixed, or the beam pipe), because there the per-part decision repeats thousands of
times and the composite's 7× node entry cost and Monte-Carlo capacity are paid at scale.

---

## 6. What was run, and what was not

Reproduce (converter env; the gate composes O2 + OCC):

```bash
cd $HOME/alisw/O2/scripts/geometry
$OCCPY probes/cellCountProbe.py --self-test                       # 7/7
$OCCPY probes/cellCountProbe.py --db <bagger-gate-workdir>/db \
    --include Base_ Boom_ Stick_ BucketLink Bucket_ --markdown    # sec 2.1
$OCCPY probes/cellCountProbe.py --model ALICE_3_example/CAD_noETA.stp \
    --include ST2487455 ST1A38495_01 ST1A38526_01 ST2487459_01 ST2487721_01 \
              ST2487462_01 ST1829909 ST1782525_01 --markdown      # sec 2.2/2.3
$OCCPY probes/cellCountProbe.py --model ALICE_3_example/CAD_noETA.stp \
    --include ST0923290 --timeout 120 --markdown                  # sec 2.4
python3 csg/decline_catalogue.py --gate-db <fixtures-workdir>/db \
    --run Bagger=<bagger-workdir>/db/Bagger --run ALICE3=<alice3-conv-dir> \
    --out website_data/decline_reasons.json                       # sec 1
```

Not measured, and why:

- **Any flat-cell query cost.** No flat solid exists; §3.3's flat-DNF numbers are derived from
  Stream_P's measured anchors and are labelled estimates throughout.
- **The composite emission of a decomposition** (§5 step 3). The emitter does not exist yet;
  its predicted cost comes from the measured K-ladder law, which was built from unions of tubes,
  not intersections with subtractions — treat the µs figures as scaling estimates.
- **IRIS.** Time-boxed out; its census profile tracks ALICE3's (shared parts), and its
  decomposition would add placements, not new mechanisms.
- **Exactness of any produced cell.** The probe counts cells and checks volume conservation; it
  never wrote a cell as sign conditions nor accepted one against the CAD body. That is
  deliberately left to the emitter stage, where the existing two-test acceptance applies.
- **Splitter behaviour on tangential carriers at scale.** The known OCCT cliff; the corpus
  here did not exercise it beyond the near-tangential filter, and a larger mechanical corpus
  might.

---

## Appendix - the per-part decline table (part x {ships-as, whyNotCSG, whyNotSurface})

Generated by `csg/decline_catalogue.py`; the machine-readable original, with untruncated
reasons, is `website_data/decline_reasons.json`. A `-` means the part ships in that
representation.

| model | part | faces | ships as | why not CSG | why not SurfaceSolid |
| --- | --- | ---: | --- | --- | --- |
| Bagger | `BasePin` | 3 | **csg** | — | — |
| Bagger | `Base` | 44 | **surface** | 7 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [44 faces: 16 cylinder, 28 plane; 7 axis cluster(s)] | — |
| Bagger | `Boom` | 31 | **surface** | 6 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [31 faces: 14 cylinder, 17 plane; 6 axis cluster(s)] | — |
| Bagger | `Stick` | 24 | **surface** | 6 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [24 faces: 12 cylinder, 12 plane; 6 axis cluster(s)] | — |
| Bagger | `Bucket` | 97 | **surface** | toroidal faces: 2 of 97 (out of the recogniser's scope) | — |
| Bagger | `BucketLink2` | 23 | **surface** | 3 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [23 faces: 11 cylinder, 12 plane; 3 axis cluster(s)] | — |
| Bagger | `BucketLink1` | 22 | **surface** | a planar face is neither a cap nor a wedge of any axis cluster [22 faces: 10 cylinder, 12 plane; 2 axis cluster(s)] | — |
| Bagger | `BoomCylinderOuter` | 8 | **csg** | — | — |
| Bagger | `BoomCylinderInner` | 6 | **csg** | — | — |
| Bagger | `StickCylinderInner` | 6 | **csg** | — | — |
| Bagger | `StickCylinderOuter` | 8 | **csg** | — | — |
| Bagger | `BucketCylinderInner` | 6 | **csg** | — | — |
| Bagger | `BucketCylinderOuter` | 10 | **csg** | — | — |
| ALICE3 | `SOLID` | 50 | **mesh** | free-form faces: 4 of 50 (surface kind outside plane/cylinder/cone/sphere) | 2 face(s): bspline face extraction not implemented yet |
| ALICE3 | `ST2487455_01` | 66 | **surface** | 8 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [66 faces: 2 cone, 31 cylinder, 33 plane; 8 axis cluster(s)] | — |
| ALICE3 | `ST2487455_002` | 6 | **surface** | toroidal faces: 1 of 6 (out of the recogniser's scope) | — |
| ALICE3 | `ST2487455_003` | 245 | **mesh** | free-form faces: 40 of 245 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 120 of 245 (out of the recogniser's scope) | 40 face(s): bspline face extraction not implemented yet |
| ALICE3 | `ST1A38495_01` | 65 | **surface** | 19 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [65 faces: 9 cone, 36 cylinder, 20 plane; 19 axis cluster(s... | — |
| ALICE3 | `ST1A38526_01` | 53 | **surface** | 13 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [53 faces: 9 cone, 24 cylinder, 20 plane; 13 axis cluster(s... | — |
| ALICE3 | `ST1A38476_01` | 148 | **mesh** | free-form faces: 40 of 148 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 11 of 148 (out of the recogniser's scope) | 19 face(s): bspline face extraction not implemented yet; 18 face(s): extrusion face extraction not implemented yet; 3 face(s): revolution face extract... |
| ALICE3 | `ST1A38494_01` | 22 | **mesh** | free-form faces: 14 of 22 (surface kind outside plane/cylinder/cone/sphere) | 14 face(s): extrusion face extraction not implemented yet |
| ALICE3 | `ST1A38494_002` | 6 | **csg** | — | — |
| ALICE3 | `ST1A38494_003` | 6 | **csg** | — | — |
| ALICE3 | `ST1A38494_004` | 10 | **surface** | a planar face is neither a cap nor a wedge of any axis cluster [10 faces: 2 cylinder, 8 plane; 1 axis cluster(s)] | — |
| ALICE3 | `ST1A38486_01` | 148 | **mesh** | free-form faces: 40 of 148 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 11 of 148 (out of the recogniser's scope) | 19 face(s): bspline face extraction not implemented yet; 18 face(s): extrusion face extraction not implemented yet; 3 face(s): revolution face extract... |
| ALICE3 | `ST2487457_01` | 996 | **mesh** | free-form faces: 69 of 996 (surface kind outside plane/cylinder/cone/sphere) | 69 face(s): bspline face extraction not implemented yet |
| ALICE3 | `ST2487459_01` | 202 | **surface** | toroidal faces: 82 of 202 (out of the recogniser's scope) | — |
| ALICE3 | `ST2487459_002` | 10 | **mesh** | free-form faces: 2 of 10 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 2 of 10 (out of the recogniser's scope) | 2 face(s): revolution face extraction not implemented yet |
| ALICE3 | `ST2487459_003` | 10 | **mesh** | free-form faces: 2 of 10 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 2 of 10 (out of the recogniser's scope) | 2 face(s): revolution face extraction not implemented yet |
| ALICE3 | `ST2487459_004` | 10 | **mesh** | free-form faces: 2 of 10 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 2 of 10 (out of the recogniser's scope) | 2 face(s): revolution face extraction not implemented yet |
| ALICE3 | `ST2487461_01` | 128 | **mesh** | free-form faces: 69 of 128 (surface kind outside plane/cylinder/cone/sphere); toroidal faces: 4 of 128 (out of the recogniser's scope) | 46 face(s): bspline face extraction not implemented yet; 23 face(s): extrusion face extraction not implemented yet |
| ALICE3 | `ST2487721_01` | 32 | **surface** | 7 axis clusters: beyond the recogniser's scope (Tier 3 territory, deliberately not built) [32 faces: 8 cone, 20 cylinder, 4 plane; 7 axis cluster(s)] | — |
| ALICE3 | `ST2487462_01` | 80 | **surface** | toroidal faces: 14 of 80 (out of the recogniser's scope) | — |
| ALICE3 | `ST1829909_01` | 1052 | **surface** | toroidal faces: 16 of 1052 (out of the recogniser's scope) | — |
| ALICE3 | `ST1829909_002` | 965 | **surface** | toroidal faces: 20 of 965 (out of the recogniser's scope) | — |
| ALICE3 | `ST1829909_003` | 236 | **surface** | toroidal faces: 20 of 236 (out of the recogniser's scope) | — |
| ALICE3 | `ST1829909_004` | 720 | **surface** | toroidal faces: 33 of 720 (out of the recogniser's scope) | — |
| ALICE3 | `ST1782525_01` | 50 | **surface** | toroidal faces: 12 of 50 (out of the recogniser's scope) | — |
| ALICE3 | `ST0923290_01` | 352 | **mesh** | free-form faces: 274 of 352 (surface kind outside plane/cylinder/cone/sphere) | 68 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-al... |
| ALICE3 | `ST0923290_002` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_003` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_004` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_005` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_006` | 45 | **mesh** | free-form faces: 29 of 45 (surface kind outside plane/cylinder/cone/sphere) | 6 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-ali... |
| ALICE3 | `ST0923290_007` | 45 | **mesh** | free-form faces: 29 of 45 (surface kind outside plane/cylinder/cone/sphere) | 6 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-ali... |
| ALICE3 | `ST0923290_008` | 45 | **mesh** | free-form faces: 29 of 45 (surface kind outside plane/cylinder/cone/sphere) | 6 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-ali... |
| ALICE3 | `ST0923290_009` | 45 | **mesh** | free-form faces: 29 of 45 (surface kind outside plane/cylinder/cone/sphere) | 6 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-ali... |
| ALICE3 | `ST0923290_010` | 86 | **mesh** | free-form faces: 29 of 86 (surface kind outside plane/cylinder/cone/sphere) | 6 face(s): bspline face extraction not implemented yet |
| ALICE3 | `ST0923290_011` | 193 | **mesh** | free-form faces: 161 of 193 (surface kind outside plane/cylinder/cone/sphere) | 28 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-al... |
| ALICE3 | `ST0923290_012` | 10 | **surface** | free-form faces: 4 of 10 (surface kind outside plane/cylinder/cone/sphere) | — |
| ALICE3 | `ST0923290_013` | 20 | **surface** | free-form faces: 9 of 20 (surface kind outside plane/cylinder/cone/sphere) | — |
| ALICE3 | `ST0923290_014` | 138 | **mesh** | free-form faces: 101 of 138 (surface kind outside plane/cylinder/cone/sphere) | 6 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-ali... |
| ALICE3 | `ST0923290_015` | 246 | **mesh** | free-form faces: 185 of 246 (surface kind outside plane/cylinder/cone/sphere) | 68 face(s): bspline face extraction not implemented yet; 12 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as... |
| ALICE3 | `ST0923290_016` | 59 | **mesh** | free-form faces: 41 of 59 (surface kind outside plane/cylinder/cone/sphere) | 5 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric boundary edge is not axis-ali... |
| ALICE3 | `ST0923290_017` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_018` | 48 | **surface** | free-form faces: 27 of 48 (surface kind outside plane/cylinder/cone/sphere) | — |
| ALICE3 | `ST0923290_019` | 44 | **surface** | free-form faces: 30 of 44 (surface kind outside plane/cylinder/cone/sphere) | — |
| ALICE3 | `ST0923290_020` | 37 | **surface** | free-form faces: 17 of 37 (surface kind outside plane/cylinder/cone/sphere) | — |
| ALICE3 | `ST0923290_021` | 33 | **mesh** | free-form faces: 16 of 33 (surface kind outside plane/cylinder/cone/sphere) | 1 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as cone but recognized quadric trim wire has fewer than 3 ed... |
| ALICE3 | `ST0923290_022` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_023` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_024` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_025` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_026` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_027` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST0923290_028` | 24 | **mesh** | free-form faces: 9 of 24 (surface kind outside plane/cylinder/cone/sphere) | 3 face(s): bspline face extraction not implemented yet; recognition attempted: recognized as sphere but recognized quadric boundary edge is not axis-a... |
| ALICE3 | `ST2487195_01` | 182 | **mesh** | free-form faces: 138 of 182 (surface kind outside plane/cylinder/cone/sphere) | 138 face(s): bspline face extraction not implemented yet |
| ALICE3 | `ST2487458_01` | 2034 | **mesh** | free-form faces: 924 of 2034 (surface kind outside plane/cylinder/cone/sphere) | 924 face(s): bspline face extraction not implemented yet |
| box | `box` | 6 | **csg** | — | — |
| box_minus_cyl | `box_minus_cyl` | 7 | **surface** | a planar face is neither a cap nor a wedge of any axis cluster [7 faces: 1 cylinder, 6 plane; 1 axis cluster(s)] | — |
| box_union_box | `box_union_box` | 10 | **surface** | 10 planar faces: not a six-plane box [10 faces: 10 plane; 0 axis cluster(s)] | — |
| cyl_cross_cyl | `cyl_cross_cyl` | 8 | **csg** | — | — |
| cyl_inter_cyl | `cyl_inter_cyl` | 6 | **surface** | symmetric difference 1.8997 cm^3 exceeds the band 1.6e-06 cm^3 (= 1e-07 cm x area 16 cm^2); extra 1.8997, missing 0 | — |
| cyl_plus_cone | `cyl_plus_cone` | 4 | **surface** | mixed lateral surface kinds ['cone', 'cylinder'] on one axis [4 faces: 1 cone, 1 cylinder, 2 plane; 1 axis cluster(s)] | — |
| oblique_cut_cyl | `oblique_cut_cyl` | 3 | **surface** | a planar face is neither a cap nor a wedge of any axis cluster [3 faces: 1 cylinder, 2 plane; 1 axis cluster(s)] | — |
| sphere_minus_cyl | `sphere_minus_cyl` | 2 | **surface** | a sphere with additional faces is out of scope (no theta/phi cuts) [2 faces: 1 cylinder, 1 sphere; 1 axis cluster(s)] | — |
| torus_union_cyl | `torus_union_cyl` | 6 | **surface** | toroidal faces: 2 of 6 (out of the recogniser's scope) | — |
| tube_window | `tube_window` | 4 | **surface** | symmetric difference 6.03536 cm^3 exceeds the band 8.04457e-06 cm^3 (= 1e-07 cm x area 80.4457 cm^2); extra 6.03536, missing 0 | — |
