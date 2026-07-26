# Exact trim topology: shared edges, analytic point-in-trim, loud failure

Plan for the four work items opened by the 2026-07-26 finding that
`O2BVHSurfaceSolid::Contains` can be wrong by centimetres on CAD parts whose faces meet along
B-spline seams. The finding, its evidence and the milestone entry live in
[`BVHSurfaceSolid.md`](BVHSurfaceSolid.md); the harness that found it and how to run it are in
[`SolidNavigationHarness.md`](SolidNavigationHarness.md). **Read the finding there first** — this
document assumes it.

Status: **planned, not started (2026-07-26).**

## The one-paragraph version

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

## Item 4 first — fail loudly

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
to mesh chording.

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

## Item 3 — canonical recognition of trim curves (converter)

Cheapest of the four and highest value per unit effort, but partial coverage, so it does **not**
replace item 1.

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
`unexplained` and timing improvements attributable to it.

## Do not conflate with the other open item

`Contains` also disagrees with `Contains_Loop` on 301/142500 points over the three-model DB, 295 of
them in the two `oTOF` parts that `CloseShape` reports as **non-manifold** (coincident/duplicated
faces, not gaps). That is a *different* defect — order-dependent clustering in `oddCrossingParity` on
duplicated hits — and none of the four items above will fix it. It has its own entry in
`BVHSurfaceSolid.md`. Keep the two apart when reading harness output: gaps show up as `unexplained`
against the mesh with `BVH == _Loop` intact; duplicated faces show up as `BVH != _Loop`.

## Suggested order and checkpoints

1. Item 4 (fail loudly) — small, unblocks honest measurement of everything after it.
2. Item 3's *measurement* only (how many trim curves are canonical?) — cheap, and it sizes items 1-3.
3. Item 1 (shared edges) — the real fix; `BoomCylinderOuter` closure goes to 0 boundary edges.
4. Item 3's converter change — folds in the cheap analytic path.
5. Item 2 (Bézier point-in-trim) — the remaining free-form seams and the last sampling step.

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
