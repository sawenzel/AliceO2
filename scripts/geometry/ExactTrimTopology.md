# Exact trim topology: shared edges, analytic point-in-trim, loud failure

Plan for the four work items opened by the 2026-07-26 finding that
`O2BVHSurfaceSolid::Contains` can be wrong by centimetres on CAD parts whose faces meet along
B-spline seams. The finding, its evidence and the milestone entry live in
[`BVHSurfaceSolid.md`](BVHSurfaceSolid.md); the harness that found it and how to run it are in
[`SolidNavigationHarness.md`](SolidNavigationHarness.md). **Read the finding there first** — this
document assumes it.

Status: **superseded in its diagnosis (2026-07-26).** Items 4 and 3 are done. Items 1 and 2 are
**not recommended**: the defect that motivated this whole document turned out to be a kernel bug —
a closed B-spline trim curve flattened to nothing — not the lost shared-edge topology this plan
assumes. Fixing it cut the reproducer's `contains` error by two thirds and every Bagger cylinder's
`distin` error to zero, with no converter change at all. Read
[What the 2026-07-26 pass established](#what-the-2026-07-26-pass-established) and
[Status after the fix](#status-after-the-fix) before acting on anything below them.

## The one-paragraph version

> **This paragraph is wrong and is kept only because the rest of the document argues from it.** It
> is a plausible causal story that was never tested before two of the four items were built. The
> measured cause is a kernel bug — a closed B-spline trim curve flattened to two coincident points
> and vanished — and the shared-edge disagreement it blames was measured at `1.3e-5` model units,
> three orders of magnitude too small to produce the observed error. See
> [What the 2026-07-26 pass established](#what-the-2026-07-26-pass-established).

In a BREP, watertightness is a *topological* guarantee: two adjacent faces reference the same
`TopoDS_Edge`, so there is exactly one 3D curve between them. Each face *also* stores a **pcurve**,
that edge's image in that face's own `(u, v)` domain — and the two pcurves of one edge are
independently fitted representations carried with their own tolerance, not derived from each other.
`O2_CADtoTGeo.py`'s `_quadric_trim_wire` extracts the per-face pcurve
(`BRep_Tool.CurveOnSurface(edge, face)`, `O2_CADtoTGeo.py:1012`) and never records which edge it
came from. **The shared-edge identity is discarded at extraction time**, so the two patches disagree
about where their common boundary is, and nothing downstream can detect or repair it. The kernel
then compounds this by flattening each B-spline pcurve to a polyline *per face* for the point-in-trim
winding test. The result is a sliver gap; parity containment casts one fixed skew ray, so every
point whose ray threads that gap loses a crossing and flips to "inside". The error is not the size
of the gap — it is the gap's shadow.

## What the 2026-07-26 pass established

Items 4 and 3 are implemented and measured. The measurements moved two of this plan's load-bearing
assumptions, so read this before working on items 1 or 2.

**Item 4 (fail loudly) — done.** `O2BVHSurfaceSolid` now carries a
`NavigationReliability` state (`Undetermined` / `Reliable` / `ReversedFaces` / `OpenSurfaceSet` /
`NonManifold`, ordered by severity, worst defect wins) with `IsNavigable()`,
`GetNavigationReliabilityName()` and the raw `GetBoundaryEdgeCount()` /
`GetNonManifoldEdgeCount()` / `GetReversedEdgeCount()`. `CloseShape` raises the three closure
defects from `Warning` to `Error` and each message now states the consequence rather than the
count alone. The harness prints the state per part, flags a non-navigable part next to its accuracy
columns, carries it in `--json` under `navigation`, and repeats the list at the end of a run.
Test: `NavigationReliabilityIsQueryable`. Measured: of the three-model DB's 19 parts **8 are
reliable, 9 open-surface-set, 2 non-manifold**; of ALICE3's 19 loaded parts **5 reliable, 14
open-surface-set**. So 11 of 19 and 14 of 19 parts respectively were answering navigation queries
that nothing in the output previously flagged as undefined.

**Item 3 (canonical trim curves) — done, and it is semantics-preserving.**
`analyze_surface_geometry.py --trim-curves` measures the coverage; `_recognize_canonical_curve` in
`O2_CADtoTGeo.py` acts on it from both the quadric path (`_quadric_trim_wire`) and the planar path
(`extract_planar_face`). A straight line is recognised from *control-polygon collinearity*, which
is a proof rather than a fit — the curve lies in the convex hull of its poles — and a circle from a
dense sample at `1e-9` relative residual. Over the three-model DB stored B-spline trim edges drop
**88 -> 50** (28 became arcs, 10 lines) and **every `unexplained` count is bit-identical**, which
is the result an exact recognition must produce: it changes the representation, never the geometry.

**Two things the measurement did not support.**

1. *The reproducer is not fixed by item 3.* `BoomCylinderOuter` still reports 699 boundary edges
   and the same `contains` `unexplained` count. What item 3 recognised there were its four
   **periodic cylinder seams** (`phi = 0` and `phi = 2*pi`, written by the kernel as 25-pole
   B-splines and exactly straight in `(phi, h)`). The seams that actually cause the defect are the
   two cylinder-cylinder intersection curves (179 and 157 poles), which are genuinely free-form.
2. *The polyline flattening is not the dominant error at these seams.* `bsplineSampleInto` is
   adaptive to a chord flatness of `1e-5` in the parametric domain, and the two pcurves of one
   shared edge were measured to disagree by at most `1.3e-5` model units on `Bagger.step` (mean
   `7.3e-8`; 6 of 705 shared edges above `1e-6`) and `2.9e-5` on `as1-oc-214.stp`. Both are far
   below the observed error. Item 2 is still right for *cost* and for removing the last sampling
   step, but it will not on its own close the gap.

**What actually went wrong: a kernel bug, not lost topology.** *(Corrected 2026-07-26; an earlier
version of this section blamed near-tangency amplification. That was inferred, not measured, and it
is wrong — the two cylinders meet at 59-60 degrees, a healthy transversal crossing.)*

Enumerating the parity ray's crossings at a wrong point (via `DistFromOutside_Loop` /
`DistFromInside_Loop`, which are "nearest entering" and "nearest exiting" independently of
containment) shows the failure is not a *missing* crossing but a **doubled** one: two consecutive
`ENTER`s about 0.04 cm apart. Locating both in their faces' domains:

- crossing 1 lies exactly on the boom tube's outer wall (`R = 1`);
- crossing 2 lies exactly on the fat tube's outer wall (`R = 1.5`) at `r = 0.9745` from the boom
  axis — i.e. **0.026 cm inside the boom footprint**, squarely inside that face's hole, which must
  therefore exclude it.

Reconstructing the hole curve directly from the sidecar (de Boor, 179 poles, rational) and running
an independent winding test confirms the curve is closed to `5.3e-12`, has exactly the intended
extent, and *does* contain the point. **The converter data was correct and the kernel disagreed
with it** — three orders of magnitude beyond any shared-edge disagreement (measured max `1.3e-5`).

The bug: `bsplineSampleRecursive` ended its recursion when the chord `p0 -> p1` fell below the
flatness scale. A **closed** curve has `p0 == p1` exactly, so a full circle flattened to two
coincident points and every polyline-based query saw an empty curve. It survived because
`signedAreaContribution` integrates by Gauss-Legendre rather than from the polyline, so the wire
still validated and still reported the right area — every check that could have caught it used the
analytic path. Fixed, with `BSplineHoleInCylinderWall` as the regression test.

**Consequence for this whole plan.** Its premise — "the converter discards shared-edge identity,
therefore adjacent patches disagree" — did not describe the reproducer. Items 1 and 2 address a
problem that was not the one causing the observed error. See
[Status after the fix](#status-after-the-fix) for what remains.

For the record, if item 1 is ever revisited: the exact intersection curve of two quadrics is a
degree-4 algebraic curve and is representable exactly in *neither* face's parametric domain with
the current `Curve2D` vocabulary, so item 1's headline representation ("shared 3D curve, mapped per
face analytically") is exact only where that image is already canonical — which is what item 3
covers. The achievable version is this plan's own ground rule, consistency rather than exactness:
derive both faces' trims from one shared sample set of the single `TopoDS_Edge` 3D curve. Note that
the closure check re-chords a cylinder's parametric edges by phi-span at `kArcSamples` per turn, so
shared samples must span less than `2*pi/24` in phi or one side subdivides again.

## Status after the fix

Three-model DB, `unexplained` counts before -> after the closed-B-spline fix:

| | `contains` | `distout` | `distin` |
| --- | --- | --- | --- |
| DB total | 4588 -> 4430 | 254 -> 251 | 218 -> **114** |
| excluding the non-manifold `oTOF` part | 525 -> **367** | | |
| `BoomCylinderOuter` (the reproducer) | 51 -> **16** | 0 | 30 -> **0** |
| `BucketCylinderInner` | 84 -> **13** | 0 | 9 -> **0** |
| `StickCylinderOuter` | 53 -> **17** | 0 | 28 -> **0** |
| `BucketCylinderOuter` | 29 -> **5** | 0 | 0 |

Every Bagger cylinder part's `distin` count goes to zero, and two parts whose `Contains` disagreed
with `Contains_Loop` now agree on all 7500 points. ALICE3 is unchanged in every column — no
regression, and its trims arrive through the recognised-NURBS path as multiple *open* edges rather
than one closed one, so it never hit this.

The closure diagnostic also became honest: `BoomCylinderInner`, `BucketCylinderInner` and
`StickCylinderInner` used to report 0 boundary edges and "navigable", but only because a whole
face's wire had vanished. The DB-wide navigable count drops 8/19 -> 5/19 as a result. That is the
diagnostic working.

**What is left, and it is no longer this document's four items.** `contains` still has 367
unexplained points outside the non-manifold part, and no part is a closed manifold yet. Before
resuming any converter work, find out what the *remaining* disagreements are: the one clean lesson
of this pass is that a plausible causal story written into a plan is not evidence, and two of the
four items were executed before anyone checked whether the premise held. Diagnose first.

## Ground rules

- **No sampling where an analytic answer exists.** This is the project's premise and the explicit
  instruction behind this plan. Flattening a curve to a polyline to answer a containment question is
  the thing being removed, not a thing to be tuned.
- **Consistency is the requirement, not exactness-on-surface.** A shared trim curve does not have to
  lie exactly on either surface — STEP tolerances mean it generally does not. Parity only needs both
  sides to agree on where the boundary is. Do not spend effort projecting curves onto surfaces.
- **The regression oracle already exists.** Every change here must be checked with
  `o2-bench-detectorsbase-solid-harness --loop-crosscheck`, and the number that must go to zero is
  the `unexplained` column of `contains` (and `distin`). `BVH == _Loop` must stay bit-identical
  throughout: it is unaffected by any of this work, so if it breaks, the change did something
  unintended.
- **`Bagger/BoomCylinderOuter_0_1_1_9` is the reproducer.** 8 faces (5 cylinders + 3 annular planes,
  four carrying general B-spline trim wires), 699 boundary edges, 61/10000 sampled points wrong by
  up to 1.71 cm. Small, fast, and it fails today. Use it as the unit of progress.

## Item 4 first — fail loudly — **DONE 2026-07-26**

**Do this before the other three, even though it fixes nothing.** It is the smallest change, and it
is the reason the defect went unnoticed: `CloseShape` *already knew* (it printed "699 boundary
edges") and the solid then went on to answer navigation queries that were wrong by centimetres. A
warning that a caller can miss is not a diagnostic.

- Give `O2BVHSurfaceSolid` a queryable state — e.g. `NavigationReliability` / `IsNavigable()` —
  set by `CloseShape` from the closure report, distinguishing "closed" from "has boundary edges"
  from "non-manifold".
- Make the `CloseShape` message state the consequence, not just the count: an open surface set means
  parity containment is undefined in the shadow of every gap.
- Decide (this is a judgement call for the human) whether an open solid should refuse to be used in
  a geometry at all, or only report. The conservative default is to keep answering but make the
  state impossible to ignore.
- The harness should print this state per part and carry it in `--json`, so every future measurement
  is labelled with whether its subject was navigable.

Acceptance: running the harness over the three-model DB, every part whose surface set is open is
visibly flagged as such in the output, and no future reader can attribute its `unexplained` column
to mesh chording. **Met.** The judgement call was resolved the conservative way the paragraph above
recommends: an unnavigable solid keeps answering queries, but the state is queryable, the closure
messages are `Error` rather than `Warning`, and the harness both annotates each part and repeats
the list of unnavigable parts at the end of the run.

## Item 1 — preserve shared edges (converter)

The actual fix. Everything else is mitigation.

- In `O2_CADtoTGeo.py`, build the edge->faces map once per solid
  (`TopExp.MapShapesAndAncestors` with `TopAbs_EDGE` / `TopAbs_FACE`), so each trim edge is
  identified by its `TopoDS_Edge`, not by the face that happens to be looking at it.
- Extract each shared edge's trim geometry **once** and give both adjacent patches the same curve.
  Two candidate representations, and the choice matters:
  - *Shared 3D curve, mapped per face.* Take `BRep_Tool.Curve(edge)` and map it into each face's
    `(u, v)` analytically using that face's known analytic parametrisation (cylinder, cone, sphere,
    torus — all invertible in closed form). Both faces then derive their trim from one curve, so
    they agree by construction. This is the principled version and it composes with item 3.
  - *Shared pcurve, one side wins.* Cheaper: pick one face's pcurve as authoritative and map it to
    the other. Removes the disagreement but bakes in one face's fitting error.
- The sidecar format needs to carry the sharing, or at least not destroy it. Decide whether
  `surfaces_*.bin` gains an edge-identity field (so the kernel can know two wires are the same
  curve) or whether agreement at write time is sufficient. **Agreement at write time is probably
  sufficient** for parity and is far less invasive — start there and only add identity to the format
  if a later item needs it.
- Watch for seam edges of periodic surfaces (a full cylinder's `u = 0 / u = 2*pi` seam), which are
  shared by a face *with itself* and need the periodic unwrapping the kernel already does.

### Concrete route (revised 2026-07-26, after items 3 and 4)

Take the *second* bullet's spirit but drive it from the 3D edge, because the first bullet's
"map the shared 3D curve into `(u, v)` analytically" is exact only for canonical curves and those
are already handled by item 3. Steps, in the existing converter structure:

1. Build the edge -> faces map once per solid (`TopExp.MapShapesAndAncestors`, `TopAbs_EDGE` /
   `TopAbs_FACE`) and key it on the edge *ignoring orientation* — the two faces see the same edge
   with opposite orientations.
2. For each shared edge, sample `BRep_Tool.Curve(edge)` once at N canonical parameters, uniform
   over the edge's own `[first, last]`. This sample set is the shared object; cache it by edge.
3. Per adjacent face, turn each 3D sample into that face's `(u, v)` with
   `ShapeAnalysis_Surface(BRep_Tool.Surface(face)).ValueOfUV(point, tol)`, then push it through the
   `map_uv` each quadric extractor already passes to `_quadric_trim_wire`. This reuses the existing
   affine machinery instead of writing five per-surface inverse maps.
4. Emit the trim edge as a polyline of line segments through those points. Both faces then agree
   *exactly* at every shared vertex, and the residual sliver between vertices is a chord error that
   shrinks with N rather than an uncontrolled fitting difference.

Three traps, all found while diagnosing this:

- **Periodic unwrapping.** `ValueOfUV` returns `u` in an arbitrary `2*pi` window; it must be
  unwrapped into the window the face's own UV bounds use, or a seam-crossing wire folds over.
- **Sampling density is not free to choose.** `CloseShape` re-chords a cylinder's parametric edges
  by phi-span at `kArcSamples` (24) per turn, so the shared samples must be dense enough that each
  segment spans less than `2*pi/24` in phi — otherwise one side subdivides again and the vertices
  stop matching, which is the very thing this item exists to fix.
- **Self-shared seam edges.** A full cylinder's `u = 0 / u = 2*pi` seam is shared by a face *with
  itself*; it must not be forced through the two-faces path.

Acceptance: `BoomCylinderOuter` converts with **0 boundary edges** reported by `CloseShape`, and its
`contains` `unexplained` count in the harness is 0. Then re-run the full three-model and ALICE3 DBs
and confirm the DB-wide `unexplained` totals drop (three-model: 4588 for `contains` today, though
most of that is the separate non-manifold `oTOF` issue — see below).

## Item 2 — analytic point-in-trim (kernel)

Removes the last sampling step from containment, and the cost that dominates these parts.

- `BoundedSurface`'s B-spline wire containment currently flattens the curve to a cached polyline and
  runs a winding test. Replace with an exact 2D crossing count: split each B-spline span into Bézier
  segments and root-find the intersections of the test ray with each segment (Bézier clipping, or
  de Casteljau subdivision with a convex-hull reject), which converges to machine precision.
- The convex-hull property gives a cheap early reject per segment, so this should be *faster* than
  the polyline for typical wires, not just more accurate — the polyline is O(number of samples)
  unconditionally.
- Keep the flattened polyline for the display mesh and for `CloseShape`'s closure diagnostics; it is
  fit for those purposes. Only containment and distance must stop using it.
- Tangency and endpoint hits need the same treatment the 3D kernels already use: a crossing exactly
  at a segment endpoint must be counted once, not twice, and a tangential touch must not flip
  parity. `sameIntersection` / the cluster logic in `oddCrossingParity` is the existing precedent.

Acceptance: `ctest -R BVHSurfaceSolid` green with new kernel-level cases for a B-spline wire that is
exactly a circle (compare against the analytic circle), a wire with a cusp, and a ray through a
segment endpoint. Harness `contains` timing on `BoomCylinderOuter` should improve from today's
3102 ns/call; if it does not, the polyline was not the bottleneck and that assumption needs
re-testing before continuing.

## Item 3 — canonical recognition of trim curves (converter) — **DONE 2026-07-26**

Cheapest of the four, but the measured coverage is far below what "highest value per unit effort"
suggested, and it does **not** replace item 1 — on the reproducer it fixes nothing at all.

- CAD kernels routinely write an exact circle or an exact straight line as a B-spline. That is not
  an approximation — it is the same curve in a heavier representation. Recognising it and storing it
  as a circle/line is exact recognition, not fitting.
- This is the explicitly-deferred "cheaper half" of the canonical-form recognition milestone
  (`BVHSurfaceSolid.md`); the surface-side half is already implemented and its residual-threshold
  discipline is the model to copy: relative, tight, and a curve that is *almost* a circle stays a
  B-spline. Silent geometry changes are worse than a slow trim.
- Payoff is double: both sides of a canonical seam agree analytically (so item 1's disagreement never
  arises for that seam), and the trim test drops to the existing cheap circle/line path (so item 2's
  Bézier machinery is never entered).
- `analyze_surface_geometry.py` already does the analogous job for *surfaces* and reports the
  coverage forecast; extend it to trim curves first, so the payoff is measured before the converter
  is changed.

Acceptance: a coverage number ("N of M B-spline trim curves are exactly circles/lines") measured
over the three-model and ALICE3 DBs, then the converter change, then a re-run showing the
`unexplained` and timing improvements attributable to it. **Met, with a null result on
`unexplained` that is the correct outcome** — see below.

### Measured coverage (`analyze_surface_geometry.py --trim-curves`)

Counted only over faces whose *surface* the exact converter can represent (a trim curve on a face
that falls back to the mesh anyway would inflate the number). "Real geometry" is the
machine-precision classification of curves the file stores as B-splines.

| model | B-spline pcurves | exactly line | exactly circle | free-form |
| --- | --- | --- | --- | --- |
| `Bagger.step` | 48 | 6 | 0 | 42 |
| `as1-oc-214.stp` | 210 | 70 | 70 | 70 |

The same tool also measures the item-1 defect directly, by pushing each shared edge's *two*
pcurves through their own faces and comparing the resulting 3D polylines:

| model | shared edges | max disagreement | mean | above `1e-6` |
| --- | --- | --- | --- | --- |
| `Bagger.step` | 705 | 1.26e-5 | 7.3e-8 | 6 |
| `as1-oc-214.stp` | 354 | 2.86e-5 | 5.1e-6 | 70 |

Both models run in seconds. `oTOF System V3-R92cm.step` and `ALICE3_CAD_pure.step` were **not**
completed with this mode: the shared-edge comparison evaluates two 48-point polylines per shared
edge in a Python loop, and on models with ~9000 and ~40000 edges that has not finished in an hour.
If the number is wanted for those two, either sample fewer points or vectorise
`_max_point_to_polyline` across edges first. It was not needed here, because the converter's own
output answers the coverage question directly and over the whole model (table below).

### Result of the converter change

Counted on the converter's own output, which is the number that matters — how many trim edges the
kernel still has to flatten:

| DB | B-spline trim edges | became lines | became arcs | sidecar bytes |
| --- | --- | --- | --- | --- |
| three-model | 88 -> **50** (-43%) | +10 | +28 | 720,936 -> 708,968 (-2%) |
| ALICE3 | 15034 -> **4528** (-70%) | +10337 | +169 | 6,088,072 -> **3,645,848** (-40%) |

ALICE3 is where this pays: it is full of NURBS surfaces recognised as cylinders, whose trim edges
are overwhelmingly straight in the recognised `(phi, h)` domain, and 10506 of its 15034 B-spline
trim edges turn out to be exactly a line or a circle.

And the correctness result — measured with identical seeds and sample counts on both DBs, against
sidecars from the old converter:

| DB | `contains` unexplained | `distout` | `distin` |
| --- | --- | --- | --- |
| three-model | 4588 -> 4588 | 254 -> 254 | 218 -> 218 |
| ALICE3 | 1 -> 1 | 0 -> 0 | 0 -> 0 |

Bit-identical per part, `BVH == _Loop` unchanged, and no part lost its exact conversion (19/19 and
20/20 paired both ways). That null result is the acceptance criterion for this item rather than a
disappointment: recognition claims to change the representation and not the geometry, so any
movement in those columns would have been a bug in the recogniser.

One side effect: three ALICE3 parts (`ST1829909_002/003/004`) moved from `non-manifold` to
`open-surface-set`, because recognising a trim curve as a line changes how `appendDirectedEdges`
chords it and so shifts coincidences in the closure half-edge map. Their `unexplained` counts
stayed 0 and both states are unnavigable, so this reorders a diagnostic rather than changing
behaviour.

Recognition thresholds follow the surface-side discipline. Straightness is decided on the
*control polygon*, not on samples: collinear poles put the whole curve inside their convex hull, so
it is a proof that the curve is the segment, and a curve whose samples are straight but whose poles
are not is rejected. A curve is also rejected if it doubles back along its own chord. The circle
test is a single least-squares solve at `1e-9` relative residual; an ellipse and a line bowed by
`1e-6` are both correctly refused.

## Do not conflate with the other open item

`Contains` also disagrees with `Contains_Loop` on 301/142500 points over the three-model DB, 295 of
them in the two `oTOF` parts that `CloseShape` reports as **non-manifold** (coincident/duplicated
faces, not gaps). That is a *different* defect — order-dependent clustering in `oddCrossingParity` on
duplicated hits — and none of the four items above will fix it. It has its own entry in
`BVHSurfaceSolid.md`. Keep the two apart when reading harness output: gaps show up as `unexplained`
against the mesh with `BVH == _Loop` intact; duplicated faces show up as `BVH != _Loop`.

## Suggested order and checkpoints

1. ~~Item 4 (fail loudly)~~ — **done 2026-07-26**.
2. ~~Item 3's *measurement*~~ — **done 2026-07-26**; it sized items 1-3 and ruled item 3 out as a
   fix for the reproducer.
3. ~~Item 3's converter change~~ — **done 2026-07-26**, semantics-preserving as required.
4. **Item 1 (shared edges) — do not start.** The premise it rests on was not what caused the
   observed error. If the remaining 367 `contains` disagreements are ever traced back to genuine
   shared-edge disagreement, revisit it then, with that evidence in hand.
5. **Item 2 (Bézier point-in-trim) — not the fix either.** The flattening is adaptive to `1e-5`;
   what was broken was that a closed curve was not flattened *at all*. It remains a legitimate
   *performance* item and would remove the last sampling step, but rank it against the `Safety`
   BVH milestone rather than treating it as a correctness fix.
6. **Diagnose the residual first.** 367 `contains` points outside the non-manifold part, and no
   part closing yet. Enumerate the crossing list at a failing point, as was done here — it took one
   throwaway probe to move from a wrong three-item theory to a one-line fix.

Checkpoints are the natural build-and-test gates: after each, `ctest -R BVHSurfaceSolid` green and a
harness run over the three-model DB with `--loop-crosscheck`, recording the `unexplained` totals in
`SolidNavigationHarness.md`'s Results section so the trend is visible.

## Environment

Unchanged from `SolidNavigationHarness.md`, and the two traps that cost time this session:

- The `o2-bench-detectorsbase-solid-harness` on `$PATH` comes from the *installed* prefix and is
  stale after an incremental build. Use `<BUILD>/O2-latest/O2/stage/bin/` with
  `LD_LIBRARY_PATH=<BUILD>/O2-latest/O2/stage/lib64:$LD_LIBRARY_PATH`. The same applies to
  interpreted ROOT macros, which resolve `libO2DetectorsBase` from the installed prefix.
- `makeTestPartDB.py` needs pythonOCC on `PYTHONPATH` and OCCT on `LD_LIBRARY_PATH` on top of
  `alienv`:
  `export PYTHONPATH=$ALIBUILD_ARCH_PREFIX/pythonOCC/latest/lib/python3.10/site-packages:$PYTHONPATH`
  and `export LD_LIBRARY_PATH=$ALIBUILD_ARCH_PREFIX/OCCT/latest/lib:$LD_LIBRARY_PATH`.
  Converting `ALICE3_CAD_pure.step` takes roughly 10 minutes; the harness over its 20 parts takes
  roughly 30.
