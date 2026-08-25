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

**The plane row's scaling is an obligation on the emitter, not a convention of convenience.** A
plane must be stored with `2b = n` for a *unit* outward normal `n` — that is, `b = n/2` — and not
with any other multiple. `s·Q ≤ 0` describes the same halfspace for any positive rescaling of
`(A, b, c)`, so nothing in the geometry, in `Contains`, or in the acceptance tests would notice a
rescaled plane. What it costs is **bit identity between the accelerated queries and their `_Loop`
twins** (§6), which is this class's whole self-check discipline.

The mechanism, since it is not guessable from the row: a cell's bounding box routinely has a face
lying exactly on one of that cell's own axis-aligned plane halfspaces, and there the box face and
the halfspace surface are the same plane, crossed at one parameter. The BVH traversal reaches it
as the slab bound `(v − o_k) / d_k`; the interval clipping reaches it as the quadratic's root
`−½γ/β`. Under `2b = n` those two divide a numerator and a denominator each scaled by exactly one
half, and IEEE division returns the identical double for both, so the endpoint is one value. Any
**power-of-two** rescaling is equally safe, for the same reason. A **non-power-of-two** one is not:
an unnormalised carrier normal `(3,0,0)` gives `b_x = 1.5`, `fl(1.5·d_x)` rounds, the root moves an
ulp off the slab bound, and the interval endpoint the accelerated query returns is the slab bound
where the twin returns the root.

One ulp on an exit distance is navigationally irrelevant, so this is not a correctness rule — it is
what keeps the twins usable as an exact oracle. The `testFlatCSG` case
`the_accelerated_distances_track_their_twins_when_a_plane_is_rescaled` measures what a ×3 rescale
actually costs, so the size of the loss is on the record rather than assumed.

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
stops on a minimum box size relative to the part diagonal, an active-list length that has stopped
shrinking (§9 measures where those go), or one of two separate level budgets running out —
**depth** and **cubify**, split this way because a level cap denominated purely in tree depth
charges full price for a split that only equalises aspect ratio. `SplitBox` always halves the
current *longest* axis, so on a cell whose extents are not already close to a cube — a 20×2×2 arm
is the motivating case — the first several levels are spent making the box roughly cubic before
any of them can start detaching a leaf from a face, and a depth-only cap can be exhausted entirely
by that alone.

So a split is charged to **cubify** rather than **depth** whenever the box, before that split, is
still far from cubic (`longest > 2 · shortest` over its three extents). `depth` is untouched until
the box is within a factor of two on every axis, at which point every further split is charged to
`depth` exactly as if there were only one budget — which is also what makes the split
**aspect-ratio-neutral**: a cell that starts near-cubic never pays into `cubify` and its box tree
is identical to what a depth-only rule would have produced. This follows from an invariance in the
splitting rule itself: for extents `a ≤ b ≤ c` with `c ≤ 2a` (i.e. already within the near-cubic
threshold), halving the longest extent `c` gives `c⁄2 ≤ a ≤ b`, so the new ratio is `2b⁄(c⁄2) ≤ 2b⁄a
≤ 2` — ratio `≤ 2` is preserved by every longest-axis halving, so a box that starts there can never
become far-from-cubic. A cube sits at ratio exactly 2 every second level and, since the test is a
strict `>`, still counts as near-cubic, so the invariance holds for any cubify budget, not a
specifically tuned one.

`cubify` is a hard ceiling, not a target: a **per-path** budget (`kMaxCubifySplits` in the C++,
currently 10) on the aspect-ratio-equalising splits one root-to-leaf branch may spend before the
box is kept however far from cubic it still is. It only ever decreases going down the recursion, so
it bounds one branch of a cell's tree, not the cell as a whole — the worst case for one cell is
`2^kMaxCubifySplits` leaves along its widest branch. How far it reaches depends on the cell's
*shape*, not only its worst ratio: a rod `(r, 1, 1)` keeps the same axis longest every split, so `N`
splits buy the full `2^N` reduction; a plate `(r, r, 1)` alternates between its two long axes, so
only every other split reduces either one, covering `2^(N⁄2)` — four orders of magnitude short of
the rod case at `N = 10`. The ceiling exists purely so a pathological cell cannot recurse without
bound; an ordinary cell resolves in a handful of splits long before it matters.

One more case the arithmetic above does not fall out of automatically: a cell bbox with a genuinely
zero (or merely tiny) extent on one axis is not something `CloseShape`'s bbox validation rejects —
it checks only for unset, inverted or non-finite bounds — but it pins `shortest` at (near) zero on
an axis that is never the longest and so never gets split, which would otherwise make
`longest > 2 · shortest` permanently true and burn the entire cubify ceiling on a cell a depth-only
rule would have resolved cheaply. The implementation floors `shortest` at `minSize` for this test
only: once the other two axes have shrunk to the size floor there is nothing left worth cubifying
towards, and the box falls through to the ordinary `longest ≤ minSize` stop instead.

The surviving boxes across all cells go into one `bvh::v2::Bvh` — the same substrate validated
three times on this branch: the patch BVH, the sub-patch BVH (`Stream_X`), and `O2BVHAssembly`
(`Stream_AE`).

The pruning is conservative in the only direction that is safe: an over-wide range bound loses
pruning, never correctness.

**The cell AABB the converter supplies must contain the cell, and that is a correctness
obligation.** An intersection of halfspaces does not bound itself, so `SetCellBBox` is the only
thing that says where a cell ends, and the split starts from it: no box is ever built outside it.
A cell that reaches past its own declared box therefore has material the accelerated `Contains`
cannot find, because it looks only inside boxes, and that the accelerated distances cut short --
while the `_Loop` twins, which know nothing of boxes, still see it. The disagreement is silent and
one-sided, which is the shape of defect §4.4 and §6 exist to prevent, so the converter owes an
**outer** bound: erring wide costs only boxes that prune to nothing, erring narrow loses material.

**Which of the two answers is wrong depends on where the box came from, and for this converter it
is the twin** (added 2026-08-25, task 9). `recognise._flat_cell_box` hands in the bounding box of
the CAD *piece* the cell was read off, so a point outside that box is outside the piece and
outside the part: a cell still holding out there is **larger than the part**, the accelerated
`Contains` is right and `Contains_Loop` is inventing material. Widening the box would silence the
twins by shipping the phantom material, so it is not a repair. The sound response is to refuse the
part, which `recognise._flat_box_holds_cell` does before a sidecar is written — it refused ITS
`BPSupportLowerCollar` and TPC `TPC_ORH` on the first corpus run. The obligation itself is
discharged positively rather than by that sampler alone: `emit.crosscheck_contains` compares
`Contains` against `Contains_Loop` on the shipped shape, and a disagreement drops the part a tier.

`CloseShape` already refuses a bbox that is unset, inverted or non-finite, with an `Error` and an
all-or-nothing build. Containment it cannot decide cheaply in general, so it carries a
**debug-build check** instead: each face is sampled on a 5x5 grid offset outward by `1e-6` of the
part diagonal, and every sample must fail `CellContains`. That catches the ordinary way this goes
wrong -- a bbox measured on the wrong piece, or in the wrong frame -- without putting a sampling
loop on the release path.

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

`DistFromOutside` is then the first entry point at `t > 0` across the boxes the ray meets. An
earlier draft of this section said the traversal may stop at the first box that produces an entry,
because boxes are visited in `t` order. **That is unsound as written, and the implementation does
not do it.** The entry only counts if its interval's exit clears `TGeoShape::Tolerance()` -- a point
already on the boundary is inside, not entering -- so the first box to produce an entry can produce
a sub-tolerance sliver the answer must reject, and the real entry then lies in a later box. The same
tolerance rule is why a cell's pieces have to be rejoined across the boxes that tile it *before* the
rule is applied: a box boundary crossed within the tolerance of the ray origin would otherwise cut
the real interval into a stub the rule discards.

The shipped traversal is therefore **order-independent**: it visits every box the ray meets, each
box contributes its own intervals over its own clipped window, the pieces are rejoined per cell, and
the rule is applied to the rejoined intervals. That is exactly what `DistFromOutside_Loop` computes,
which is the point. A sound early exit does exist -- a box whose slab entry already exceeds the
standing best cannot improve it -- and §9 measures whether it is worth the code.
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
  lower bound on the distance to the boundary. If `p` is in an undecided box, that box's active
  list may hold boundary anywhere within it, so its face distance is NOT a bound at all — the
  sound answer there is `0`. §9 measures how often that fallback fires and whether refining
  undecided boxes around a point would pay for itself.

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
same volume guard, the same `PART_MAX_CELLS = 64`.

> **Corrected 2026-08-25, at the close of the rung (task 9).** This paragraph originally claimed
> "every part in the demand table is under it". That is false, and it is what keeps ITS
> `IBCYSSFlangeA` declined after the rung landed. `decompose.split_into_cells` tests
> `len(cells) + len(pending) + len(unresolved) > max_cells` at the top of each iteration
> (`csg/decompose.py:305`), so the bound is on the whole working set and not on the terminal cell
> count: 59 terminal cells with a non-empty queue trips it. `PART_MAX_CELLS` — not either of the
> flat budgets, which the largest shipped part uses 47/256 and 282/1024 of — is now what decides
> this rung's coverage. Raising it stays out of scope here because this section freezes
> `decompose.py`; it is a real question for §9's measurement rung.

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
  that stays true of the DNF *description*; the C++ class and `primitives.flat_cells` both accept
  one.

  > **Corrected 2026-08-25 (task 9).** The named justification — "the ITS connector blocks that
  > the 8-leaf per-cell budget refuses are one cell each and belong here" — did not survive the
  > corpus run. Those blocks are no longer one cell: `IBConnectorBlockBodyASide` decomposes into
  > 17 pieces and is refused because cell 3's halfspaces do not close it up (see the ownership
  > rule below and `recognise._flat_box_holds_cell`), and `IBConnectorBlockBodyCSide` is refused
  > on the boundary gap. The rung therefore does **not** deliver the single-cell case this bullet
  > names; one cell stays legal in the description and in the class, and nothing in the five
  > corpora exercises it.

- **ownership is not reopened by running last.** The flat path runs only after every matcher
  above it declines, and a *one-piece* decomposition keeps the two whole-part guards
  `_match_single_cell` applies: an all-planar body belongs to the prism family's templates, and a
  one-carrier body to the tier-1 ones (`_cell_leaves(..., whole_part=len(pieces) == 1)`). Running
  last is a licence to take what nobody else can state, not a licence to overrule a matcher that
  looked at a part and judged it not one of its own. Concretely: a three-section prism 1e-5 cm out
  of similarity, which the prism family refuses because it is not a `TGeoXtru`, would otherwise
  ship here as an exact ten-plane cell — a body the recogniser had already judged unlike its
  template, accepted on a technicality. Multi-cell parts are unaffected; every part this rung
  ships is 12 to 47 cells.

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
