# Design — `O2FlatCSG`, the flat-DNF halfspace solid with a sub-cell BVH

**Written 2026-08-24, approved by Sandro at the close of the R4 review.** This is the design
document for rung R5 of the flat-CSG programme (`Handoff_FlatCSG.md` §2). It fixes what gets
built before any of it is built; the implementation plan derives from it.

Read first: `Handoff_FlatCSG.md` (the programme and R5's measured demand), `Stream_AA_FlatCSG.md`
(the plan this executes, its §5 step 4), `Stream_AJ_Recognition.md` (what the recogniser and the
cell emitter already do), `NEXT.md` (state and environment traps).

---

## 1. The problem this solves

R4 left 21 parts that decompose into cells perfectly well and then decline, because the union of
those cells would have to ship as a `TGeoCompositeShape` — a binary boolean tree — and the tree
is refused at `recognise._PART_MAX_LEAVES = 64` halfspaces or `decompose.PART_MAX_CELLS = 64`
cells. The measured table is `Handoff_FlatCSG.md`'s R5 block: TRD `BREF1` at 47 cells and 66+
leaves, ITS `IBCYSSFlangeA` at ≥59 cells, TPC's two-hole GEM cap `TPC_OCGEM` at 12 cells, four
`ST1829909_01` bodies at ≥57 cells each, and fourteen more.

The budget is not arbitrary. A 66-leaf boolean tree is 66 nested `TGeoBoolNode`s, each of which
answers a query by recursing into both children, so a `Contains` costs 66 virtual calls whatever
the point is and a `DistFromOutside` costs considerably more. The budget is an honest statement
that ROOT's representation stops being usable there.

There is a second problem, quieter and arguably worse, and it is a *fidelity* problem rather than
a cost one. A cell is an intersection of **halfspaces**, and ROOT has no usable halfspace, so
`recognise._cell_leaf` bounds each one into a padded native primitive sized to the part's
bounding box grown by `_CELL_MARGIN` (`_CellBox.window`). The shipped solid is therefore the
intersection of *padded boxes, tubes and cones*, which equals the cell on the part's neighbourhood
and is guaranteed to equal it nowhere else. It is correct, it is tested, and it carries a margin
constant that has to be right.

`O2FlatCSG` removes both. It stores the halfspaces themselves, so there is no padding and no
margin; and its cost grows with the number of halfspaces that are *locally undecided* rather than
with the number the part has.

## 2. Scope

**In:** the C++ shape class, its sidecar IO, the converter path that emits it, the raised budget
so the 21 parts ship, the three-test acceptance on those parts, the `_Loop` twins and their
ctest cases, and the crossover measurement against plain-composite emission.

**Out, and recorded as such:**

- **Geant4.** Sandro, at the design review: "Geant4 is not important to me at this moment.
  `G4MultiUnion` is merely an inspiration." No `TVirtualGeoConverter` seam in this rung. §11
  keeps the one-paragraph correspondence for whoever picks it up.
- **The 24 boundary-gap declines** — the parts where a dihedrally-convex splitter piece is not a
  cell of the carrier arrangement (2–23 % off, median 2 %). Fixing them means splitting at every
  carrier crossing rather than at trusted concave edges, which is a change to `csg/decompose.py`
  and an independent risk. They stay declined. This class is what makes that change *affordable*
  later, since its whole point is that the multiplied cell count stops being a problem.
- **R6, the census** over the remaining Run 3 modules. Sandro chose to go straight into R5.

## 3. The representation

A part is a **union of cells**; a cell is an **intersection of halfspaces**; a halfspace is a
signed implicit surface. Depth is two, always. That is the whole data model.

### 3.1 A halfspace

Every carrier `recognise._halfspace_carriers` produces is one of five kinds — plane, cylinder,
cone, sphere, torus — plus the elliptic cylinder that R2 taught the reader. Four of the five and
the ellipse are **quadrics**, so they collapse into a single fixed-length block:

```
Q(x) = xᵀ A x + 2 bᵀ x + c        A symmetric 3×3
```

stored as ten doubles `(a00, a01, a02, a11, a12, a22, b0, b1, b2, c)` with a sign `s ∈ {+1, −1}`.
The material side of the halfspace is `s · Q(x) ≤ 0`.

This is the "mathematically generic" instruction taken literally: there is **one** surface
evaluator, **one** ray intersector, and no per-kind branch on the hot path. It is also the
AOT-codegen-friendly layout `Stream_AA` asked for, because a cell is then a contiguous run of
identical blocks and specialising a cell is unrolling a loop over them.

The carrier-to-quadric map, for the record:

| carrier | `A` | `2b` | `c` |
| --- | --- | --- | --- |
| plane, outward unit `n` through `p` | `0` | `n` | `−n·p` |
| sphere, centre `p`, radius `r` | `I` | `−2p` | `\|p\|² − r²` |
| cylinder, axis `(p, d)`, radius `r` | `I − ddᵀ` | `−2Ap` | `pᵀAp − r²` |
| cone, axis `(p, d)`, ref radius `r`, semi-angle `α`, `k = tan α` | `I − (1+k²) ddᵀ` | `−2Ap − 2rk d` | `pᵀAp + 2rk (p·d) − r²` |
| elliptic cylinder, axes `x̂, ŷ`, semi-axes `a, b` | `x̂x̂ᵀ/a² + ŷŷᵀ/b²` | `−2Ap` | `pᵀAp − 1` |

The **torus** is the one carrier that is not a quadric. It gets a second block type carrying its
canonical parameters `(p, d, R, r)` — eight doubles — because both things the class needs
from it are cleaner in canonical form than in expanded quartic coefficients: the ray intersection
is a quartic the kernel already knows how to solve (it solves them for the surface representation),
and the range bound of §4.2 wants the exact signed distance, which the canonical form gives and
the quartic does not. No reference direction is stored: a carrier torus is a full
surface of revolution, and whatever bounds it in phi is a plane carrier of the same cell.

Both block types are fixed length. A part's halfspaces live in one flat array; a cell is a
`(first, count)` pair into it.

### 3.2 Why this is more faithful, not just faster

Stated once more because it is the part that is easy to miss: the composite emission has to
choose a finite extent for every halfspace and does so from `_CellBox`, whose `margin` is
`_CELL_MARGIN` times the part diagonal. `O2FlatCSG` has no such constant. The solid it ships is
the intersection of the halfspaces the faces actually carry, everywhere, not on a neighbourhood.

### 3.3 The sign convention, and how it is checked

Two things decide the material side: the carrier's own orientation (a plane's normal is already
flipped for `TopAbs_REVERSED` in `_halfspace_carriers`) and the `side` field that
`census.halfspace_side` decides. They compose into the single `s` above.

Getting that composition wrong is the most likely way to ship a wrong solid, and it would be a
*silent* wrongness — an inverted halfspace still produces a solid. So the emitter does not get to
be trusted about it. The existing, corpus-validated `_cell_leaf` is the oracle: for every cell of
every part on the flat path, a test asserts that the flat halfspace conjunction and
`_cell_leaf`'s padded-primitive conjunction classify an identical sample of points inside the
part's bounding box the same way. `_cell_leaf` has shipped 1775 known-source-clean parts; it is
the right thing to be measured against.

## 4. The sub-cell BVH

Sandro's addition at the design review: bake the sub-cell split in from the start rather than
retrofit it. The BVH primitives are **sub-boxes of cells**, not cells.

### 4.1 Why not cells

A cell's AABB is a bad proxy for the cell on exactly the bodies this class is for. A long
diagonal rod, a curved cell following a torus, an L-shaped cell — each has an AABB mostly full of
nothing, so a cell-level BVH lets rays and points into cells they cannot possibly be in, and then
pays the full halfspace list to find out. The parts in the demand table are 10 to 59 cells of
precisely this character.

### 4.2 The build

For each cell, start from its AABB and recursively split (median of the longest axis). At each
box, classify every halfspace still active in the parent using a **rigorous range bound** of the
implicit function over the box:

- **quadric.** With box centre `m` and half-extents `h`, write `u = x − m`, so
  `Q(x) = Q(m) + 2 g·u + uᵀAu` with `g = Am + b`. Then
  `|Q(x) − Q(m)| ≤ L` with `L = 2 Σᵢ |gᵢ| hᵢ + Σᵢⱼ |Aᵢⱼ| hᵢ hⱼ`,
  giving `Q ∈ [Q(m) − L, Q(m) + L]`. Cheap, rigorous, and it tightens **quadratically** as the
  boxes shrink, which is exactly the behaviour recursive subdivision wants.
- **torus.** Its exact signed distance `f(x) = √((ρ − R)² + z²) − r` (with `ρ` the distance from
  the axis and `z` the coordinate along it) is 1-Lipschitz, so `f ∈ [f(m) − ‖h‖, f(m) + ‖h‖]`
  is rigorous and tight.

Then, per halfspace and box:

- if `min_box s·Q > 0` the box lies **wholly outside** the halfspace, hence wholly outside the
  cell: **drop the box**;
- if `max_box s·Q ≤ 0` the halfspace holds **everywhere** in the box: **drop it from that box's
  active list**;
- otherwise it stays active and the box is undecided.

A box whose active list empties is **wholly inside** the cell and is marked `solid`. Recursion
stops on a depth cap, a minimum box size relative to the part diagonal, or an active-list length
that has stopped shrinking (§9 measures where those go). The surviving boxes across all cells go
into one `bvh::v2::Bvh` — the same substrate validated three times on this branch: the patch BVH,
the sub-patch BVH (`Stream_X`), and `O2BVHAssembly` (`Stream_AE`).

The pruning is conservative in the only direction that is safe: an over-wide range bound loses
pruning, never correctness.

### 4.3 What it buys, in three places at once

1. **Tight boxes** on the long, diagonal and curved cells that motivate the class.
2. **Short active lists at the leaves.** This is the real win. A 70-halfspace part evaluates the
   two to four halfspaces that are undecided in the box a query actually lands in, not seventy.
3. **A rigorous, wholly generic `Safety`** — see §5.4. No closed-form point-to-quadric distance
   is needed anywhere, which matters because for a general quadric there isn't one.

### 4.4 The correctness crux

An active list is a correct description of the cell **only inside its box**. Every query must
therefore clip to the box before using its list. For a point query this is automatic. For a ray
it is not: the traversal walks boxes in `t` order and, within each box, clips the ray to that
box's slab interval before running §5.2's interval clipping over that box's active list. Any
implementation that gathers active lists across boxes and then clips once is wrong. This is the
single most important invariant in the class and the `_Loop` twins of §6 exist mainly to catch a
violation of it.

## 5. The queries

### 5.1 `Contains(p)`

BVH point query to the boxes containing `p`; for each, `p` is in the part if every halfspace in
that box's active list satisfies `s·Q(p) ≤ 0` (a `solid`-marked box answers true immediately).
Cells are disjoint by construction, so the first box that says yes is the answer.

### 5.2 Distances, by generic interval clipping

Within a box, along a ray clipped to `[t₀, t₁]`: collect the roots of every active halfspace.
A quadric gives `α t² + 2β t + γ` with `α = dᵀAd`, `β = dᵀ(Ao + b)`, `γ = Q(o)`, so at most two
roots — and the degenerate `α ≈ 0` case is a plane, or a ray parallel to a cylinder axis, and is
handled as the linear equation it is. The torus gives at most four. Sort the roots, and classify
the midpoint of each resulting sub-interval by the §5.1 test.

The result is the cell's occupancy as a set of intervals, with **no convexity assumption
anywhere** — which is required, not merely elegant, because a complemented carrier (`side ==
"exterior"`) makes cells routinely non-convex, and every part with a hole in it has one.

`DistFromOutside` is then the first entry point at `t > 0` across the boxes the ray meets, and
the traversal can stop at the first box that produces one because boxes are visited in `t` order.
`DistFromInside` is the far end of the occupancy interval containing `t = 0`, which needs the
union across cells — a point can leave one cell into an adjacent one, and the part's boundary is
where the union ends, not where a cell does. Both obey ROOT's usual `iact`/`step` contract and
`TGeoShape::Tolerance()` conventions.

### 5.3 `Capacity`

The sum of the per-cell volumes, which the converter writes into the sidecar from OCCT's `GProp`
on the piece each cell came from. The cells are disjoint — that is what the decomposition
guarantees and what `decompose`'s volume guard checks — so there is no inclusion–exclusion.

This is the honest reading of the handoff's "analytic `Capacity`": the number is exact to OCCT's
accuracy on the source solid rather than sampled, which is what makes it better than
`TGeoCompositeShape`'s Monte-Carlo. A general quadric cell has no closed-form volume and the
class does not pretend otherwise.

### 5.4 `Safety`

Entirely from the box structure, with no surface-distance formula:

- **outside:** the distance from `p` to the nearest non-empty box is a lower bound on the
  distance to the solid, because every point of the solid is in some box.
- **inside:** if `p` is in a `solid`-marked box, the distance from `p` to that box's faces is a
  lower bound on the distance to the boundary. If `p` is in an undecided box, fall back to the
  distance to that box's faces, which is still a bound but a poor one — §9 measures whether
  refining undecided boxes around a point pays, and the fallback is always the legal `0`.

Box-to-point distance is exact and trivial, so this is rigorous by construction. §11 records that
1-Lipschitz forms exist for the plane, sphere, cylinder and torus and would tighten it; that is a
measured improvement for later, not the shipped rule.

### 5.5 The rest of the `TGeoShape` contract

`ComputeBBox` from the union of the retained boxes — which is tighter than the union of the cell
AABBs, and worth noting because `NEXT.md` item 6 records loose bounding boxes as a live defect on
the surface class. `ComputeNormal` from the gradient of whichever active halfspace the point is
closest to satisfying with equality; `∇(s·Q) = 2s(Ax + b)`, normalised. `GetPointsOnSegments`,
`InspectShape`, and the drawing stubs follow `O2Tessellated`'s precedent.

## 6. Falsifiability: the `_Loop` twins

Every accelerated query gets a twin that walks **all cells and all halfspaces**, with no BVH and
no active-list pruning: `Contains_Loop`, `DistFromOutside_Loop`, `DistFromInside_Loop`,
`Safety_Loop`. The twin is the definition of the answer. The tests require bit identity between
each query and its twin over the fixture corpus, and `Safety`'s twin is required to be a
*sound* bound rather than an identical number, since the box structure is where its value comes
from.

This is the discipline `O2BVHAssembly` already carries and the reason to keep it here is concrete:
§4.4's clip-inside-the-box invariant and §4.2's two pruning rules are the three places this class
can be silently wrong, and all three are exactly what a twin catches.

## 7. Persistence

A **binary sidecar**, following `O2SurfaceSolidIO`'s pattern: the generated `geom.C` constructs
the shape by name and a `LoadFlatCSG(file, solid)` fills it from `flatcsg_*.bin`, after which the
caller calls the class's build entry point. Reasons: `geom.C` stays small; Cling's compile limits
do not get tested by a 59-cell part; and it reuses a loader shape the surface path has already
proven. `NEXT.md` item 2 records a live JIT defect on the surface path's macro emission — this
design deliberately does not add a second, larger consumer of that mechanism.

The sidecar carries the halfspace blocks, the cell `(first, count)` table, the per-cell volumes,
and the part bounding box. It does **not** carry the BVH: that is built on load, because it is
fast to build, because a stored BVH is a second thing that can disagree with the data it
describes, and because §9 will change the split parameters.

The class additionally gets a normal ROOT streamer, so a shape that reached a `TGeoManager`
survives being written to a `.root` geometry file without the sidecar. The blob is small — a
60-cell part at ~10 halfspaces per cell is under 50 kB — so there is no reason for it to be
anything cleverer.

## 8. The converter side

`csg/decompose.py` is **unchanged**: the same split loop, the same connectivity-first rule, the
same volume guard, the same `PART_MAX_CELLS = 64` (every part in the demand table is under it).

What changes:

- a **flat emitter** that reads `_halfspace_carriers`' records directly and maps them to quadric
  and torus blocks per §3.1. It does *not* go through `_cell_leaf`, which exists only to bound a
  halfspace for tree emission — but §3.3's bridge test measures the two against each other.
- the **budget**. `_PART_MAX_LEAVES = 64` refuses a boolean tree that wide and there is no tree
  any more, so the flat path gets its own bound, sized to what the sidecar and the BVH build
  carry rather than to what a `TGeoCompositeShape` can be evaluated through. The tree path keeps
  its budget unchanged.
- the **emission policy**: a part ships flat above the measured crossover and as a plain
  composite below it. `Stream_AA` predicts the crossover at ≥10 cells; §9 measures it. Until it
  is measured the policy is not written down, and no part that converts today changes
  representation.
- a **single cell is admissible.** `primitives.union_of_cells` requires two or more cells, and
  that stays true of the DNF *description*; but the C++ class accepts one, because the ITS
  connector blocks that the 8-leaf per-cell budget refuses are one cell each and belong here.

## 9. What R5 must measure to be done

1. **The 21 parts convert** at `dV_sym = 0`, pass the oracle gate, and pass
   `checkKnownSource.py` against their original `TGeoShape`.
2. **The crossover.** The same parts emitted both ways — flat, and as plain composites wherever
   the tree budget permits — and X-ray-benched against each other, so the crossover is a
   measurement and not `Stream_AA`'s prediction. Report it in cells and in halfspaces, since it
   is not obvious in advance which is the better predictor.
3. **The split knobs.** Depth cap, minimum box size, and stopping rule against leaf active-list
   length, box count, memory, and query cost. Chosen on data and recorded, the way
   `Stream_AE` recorded `O2BVHAssembly`'s honest limit at ≤68 daughters.
4. **`Safety` quality**, since §5.4's bound is structural: measured against the twin and against
   the true distance on the fixtures, because a sound but weak safety is a transport-speed
   defect that no correctness test sees.

## 10. Acceptance and floors

The three-test discipline of `Handoff_Recognition.md` §3 binds unchanged: symmetric difference,
oracle gate, known-source. The false-accept guard `accept.contains_disagreements` runs on this
path too — the R4 report is explicit that `BRepAlgoAPI_Cut` can report `IsDone()` with zero
solids in both directions (`ST0923290_01#b19`), and the flat path builds its OCCT candidate the
same way, so it inherits the defect and must inherit the guard.

Floors at handover, none of which may shrink: corpora PIPE 166/9/1 · ITS 250/14/0 · TPC 163/9/0 ·
ABSO 29/0/0 · TRD 1167/3/0; known-source 1775/1775; fixtures gate 10/10 csg exit 0; Bagger 13/13
+ 7/7 exit 0 with CSG 7 bit-identical; `emit.py --self-test` 301 with its six digest tables;
`checkKnownSource.py` 17; `O2_TGeoToCAD.py` 107; `O2_CADtoTGeo.py` 54; ctest `BVHSurfaceSolid`
113 and `BVHAssembly` 22. R5 adds a `testFlatCSG` target; it does not touch the others.

R5 touches `Detectors/Base/**`, which the recognition programme's kernel freeze forbade. That
freeze does not apply here (`Handoff_FlatCSG.md` §4); the ctest floors do.

## 11. Risks and recorded remainders

1. **The sign convention (§3.3)** is the likeliest silent wrongness. Mitigated by the bridge test
   against `_cell_leaf`, which is the strongest oracle available for it.
2. **The clip-inside-the-box invariant (§4.4)** is the likeliest silent *acceleration* bug, and
   an acceleration bug that only bites on some rays is the worst kind. Mitigated by the `_Loop`
   twins, which is what they are for.
3. **The range bound (§4.2) is conservative, not exact.** On a cell whose halfspaces are nearly
   tangent over a large region the bound stays undecided for many levels, so the box count grows
   and the benefit shrinks. That is a cost, never a correctness failure, and §9's knob
   measurement is where it becomes visible.
4. **`Capacity` depends on OCCT's `GProp`** on the source piece rather than on the shipped solid.
   If the two disagree, the acceptance tests catch it as a symmetric-difference failure first —
   but the number is inherited, and the report should say so wherever it is quoted.
5. **Geant4.** For whoever picks it up: a cell whose halfspaces are all `side == "interior"` and
   whose part is a pure union maps onto `G4MultiUnion` directly. The general case — complemented
   halfspaces, hence non-convex cells — has no `G4` counterpart and needs either a
   `G4BooleanSolid` tree per cell (which reintroduces exactly the cost this class removes) or a
   `G4VSolid` subclass mirroring this one. The second is the honest option and is its own rung.
6. **The 24 boundary-gap declines** stay declined (§2). This class is the thing that makes the
   split-at-every-carrier-crossing fix affordable, because the cell-count multiplication it
   causes is no longer what refuses a part.
