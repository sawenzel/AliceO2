# Roadmap — ideas and deferred work for the CAD → TGeo path

Distinct from [`NEXT.md`](NEXT.md), which is the hand-over for *what is in flight*. This file holds
things we have decided to do **later**, so that deferring them is a recorded decision rather than a
thing that quietly fell off the list. Nothing here is scheduled. Nothing here is started.

Each item records who raised it, what is already measured that bears on it, and what would have to
be true for it to pay. **Do not treat the annotations as commitments** — several of them are
inferences from tonight's benchmark and are labelled as such.

---

## Raised by Sandro, 2026-08-03

These are his, recorded close to verbatim, with measurement notes added underneath.

### (a) A `TGeoOCCTSolid` — OCCT itself as a shape, as fallback or as inspiration

> *"Since we have OCCT as oracle, would it make sense to make a TGeoOCCTSolid as another fallback or
> inspiration for performance?"*

**Why it is attractive.** It would make "everything is representable" trivially true. The oracle
already answers exactly the questions a `TGeoShape` must — `BRepClass3d_SolidClassifier` for
`Contains`, `IntCurvesFace_ShapeIntersector` for the distance kernels, `BRepExtrema` for `Safety` —
and it is the thing every acceptance test on this branch is scored against. A part that no tier can
represent would stop being a coverage gap and become a performance decision instead.

**What is already known.** Nothing about its speed: the oracle has only ever been run offline, in a
separate process, on sampled points. It is almost certainly far slower than anything we ship —
tonight's benchmark put `O2BVHSurfaceSolid` at 2.6 µs for `Contains` and `TGeoCompositeShape` at
18 ns, and OCCT's classifier is not built for per-step queries. Two things to check before any
design: **thread safety** (parts of OCCT are not re-entrant, and a simulation will be threaded) and
the memory cost of holding B-reps resident for a whole detector.

**The "inspiration" half may be the more valuable half** — if a `TGeoOCCTSolid` existed even as a
debug-only shape, it would let the X-ray benchmark and the gate compare *in-process* rather than
across a subprocess boundary, which is currently the slowest part of the gate.

### (b) Safety: stop early at the BVH, and cache

> *"we could use the BVH as far as possible and also return early (never talking to any final exact
> face). This was recently done in VecGeom's master in TessellatedSolid and it gains heaps. Of
> course, the safeties will be less exact which has the downside of allowing less free travel of
> particles steps. So we have to instrument the simulation and find a sweet spot. But the idea
> **needs to be implemented** in BVHSurfaceSolid. Secondly, we could employ safety caching as in
> https://github.com/AliceO2Group/AliceO2/pull/15128 (with the necessary thread-local treatment
> added)."*

**This lands exactly on what tonight's measurement left behind, and it is safe by construction.**
`Stream_S_SafetyBVH.md` made `Safety` use the BVH and got 835 µs → 22.4 µs on ALICE3's 965-patch
solid, *bit-for-bit identical* to the old loop. But it still descends to `distanceSqToPatch` on the
surviving candidates, and **that is where the remaining time is**: ~8.6 candidates on that solid,
each costing 2.4–8.6× an average patch, ≈ 9 µs of real geometry. Returning the traversal's own
bound instead of descending would remove essentially all of it.

And the correctness argument is already written down and already measured. `Safety` must never
exceed the true distance; the BVH's point-to-AABB distance is a **genuine lower bound** on the
distance to any patch in the subtree (`Stream_S` §"why the bound is a lower bound"), verified on
7 surface families × 4000 points, 28 000 checks, 0 violations. So an early return is automatically a
*valid* safety — it is only a *loose* one. That converts the whole question from correctness to the
step-length trade-off Sandro describes, which is the right place for it to live.

**What is missing is the instrument, not the idea.** Nothing here can currently measure the cost of
a loose safety, because that cost is paid in extra navigation steps during transport, not in the
shape. That is what the integration/Geant harness is for. Suggested knob: a depth or
tightness threshold, with exact-at-the-patch as one end of the range, so the sweet spot is a
measured parameter rather than a rewrite.

Safety caching (PR 15128) is orthogonal and composes with it; note the thread-local requirement.

### (c) BVH over *sub*-meshes and *sub*-surfaces — flagged as a priority

> *"We have a super-fast BVH and we should employ it as much as possible --> So I think we need to
> investigate BVH on top of *sub-meshes* or *sub-surfaces* as a priority. This likely improves
> intersection/ray-tracing as well as the safety calculations."*

**The benchmark independently points at the same place.** `Stream_P_RepresentationBench.md` §2.2
found that the surface solid's cost splits by **trim-wire complexity, not patch count**: 135–337 ns
for a line/arc trim against **2253–5908 ns** for a B-spline/ellipse trim, with nothing in between,
and 61% of `Contains` inside `Curve2D::closestPoint`. The BVH is already doing its job — 1.7
candidate patches per query out of 8. The cost is *per surviving candidate*.

Subdividing a patch attacks that from both sides: tighter boxes mean fewer candidates survive, and a
sub-patch carries a shorter trim wire, so each survivor is cheaper. It is the most direct route to
the 40–145× gap against `TGeoCompositeShape`, and it does not change any representation — the same
exact surfaces, indexed more finely. Note it would also give `Safety`'s early return (item b) a much
tighter bound to return.

### (d) Intel Embree as an optional BVH engine

> *"For tessellated solid, I recently showed that we can use Intel Embree as the BVH engine, gaining
> another factor of 2 due to SIMD and auto-CPU architecture tuning. I plan to integrate Intel Embree
> in the software stack and we should optionally be able to use it."*

Blocked on the stack integration, which Sandro is doing. Design implication for us: the BVH should
sit behind an interface thin enough that the engine is swappable, and any Embree path must be
**optional** — this branch is developed on aarch64, where the SIMD story differs. Worth keeping in
mind while doing (c), since sub-surface indexing is a change to the same layer.

### (e) Two mesh-precision knobs, and then vary them by physics relevance

> *"We can easily make it separate for linear and angular and have 2 separate knobs. Then we can keep
> in mind a physics perspective: Precision should be greater where it matters for us: Close to the
> interaction point or within some small pseudo-rapidity. A tessellated far away in the experimental
> pit does not need to be precise."*

**The first half is half-done and the measurement behind it is the interesting part.** `--mesh-prec`
currently sets `lin_defl` **and** `ang_defl` to the same value, and
`Stream_P_RepresentationBench.md` §5 decoupled them experimentally (via a `remeshFromBrep.py` probe)
and found that **`--mesh-prec` is an angular knob in disguise**: a 20× change in `lin_defl` at fixed
`ang_defl` moved nothing measurable on Bagger. So splitting the flag is not cosmetic — one of the
two was never doing anything, and anyone who has ever "tuned the mesh precision" on this project was
tuning the angle.

That matters for the second half, because the two knobs scale differently and the physics wants the
one that is currently inert:

- `lin_defl` is an **absolute** sagitta bound in cm. It means the same thing everywhere.
- `ang_defl` is an **angle**, so it refines *relative* to curvature — on a large-radius surface far
  from the IP it permits a large absolute error, which is exactly the regime this idea is about.

So a physics-driven scheme most likely wants **`lin_defl` set per volume** as a function of position
(distance from the IP, or |η|), with `ang_defl` held merely loose-enough to avoid degenerate
faceting. Setting it the other way round would make the far-field *worse* than intended.

**It would also fix the meshing memory problem directly.** ALICE3 cannot be meshed at the default
0.1 — one **2 m sphere** reached 22.9 GB resident and was killed. That sphere is precisely a
"far away, does not need precision" volume, so a spatially-varying scheme is not only a physics
optimisation but the thing that makes meshing the full model tractable at all.

**Two cautions, both measured, for whoever builds it:**

1. **Mesh validity is non-monotone in precision.** Refining `BucketLink2` from 0.1 to 0.05 took it
   from 6600 to **10697** LOST crossings and 673 → 1843 unterminated. "Coarser far away" therefore
   cannot be assumed safe either — a per-volume precision scheme has to be **validated per volume**,
   not signed off globally from a single deviation number. See [`MeshHealing.md`](MeshHealing.md).
2. **Far from the IP, chordal deviation is the wrong acceptance criterion.** Those volumes are there
   for their **material budget**, not their position, so what must be bounded is the volume error —
   which the gate already computes as `capacityRelativeDeviation`, and which a coarse mesh degrades
   in a way a sagitta bound does not describe. Suggest: near field gated on deviation, far field
   gated on capacity, with the crossover stated explicitly rather than implied by a formula.

---

## Deferred, with the reason

- **Assembly-level transport under `TGeoNavigator`** — the oracle exists (`assemblyOracle.py`,
  per-ray occupancy sets, 28 analytic self-checks green); the navigator side and a `leaking` counter
  do not. Deferred because it *"will be exercised in the Geant test automatically"*: the integration
  harness loads the whole translated geometry as `geom.C` does and operates on that.
- **`Curve2D::closestPoint`** — 61% of `Contains`, and on the nearest-patch path that `Safety` and
  `ComputeNormal` now use, so it is the shared hot spot of both. A parallel investigation, later.
  Overlaps heavily with (c), which may be the better handle on it.
- **Mesh healing** — [`MeshHealing.md`](MeshHealing.md). A mesh can be *invalid* rather than merely
  inaccurate, and chordal accuracy does not detect it. Later.
- **Overlaps in general** — demoted deliberately. `TGeoNavigator` tolerates them, the ALICE
  simulation geometry has had them for years, and the target is the navigator plugin rather than
  native Geant4; Bagger's three 0.1 mm joint overlaps are a toy model's cosmetics. **Not** to be
  confused with the coincident *double placement*, which is a different defect and is being fixed.
- **Algorithmic overlap repair at the STEP level** — *"if there is a way to algorithmically fix
  overlaps directly on the STEP level we should of course do it. Just maybe not now."* Mechanically
  straightforward with OCCT (`BRepAlgoAPI_Cut` the intersection out of one part of each pair); the
  real question is *which* part yields, which is a modelling decision and probably needs a rule per
  assembly rather than a global one.
- **Exact arrangement-cell enumeration for trims** — the only remaining exact route for ALICE3's
  unrepresented solids (`Stream_R_CoSurfaceTrims.md` §9.3). Reaches 13/15 but cannot be derived by
  sampling. Research-grade; absent it those solids ship tessellated, which is the standing bargain.
- **Free-form surfaces** — 19 of 55 ALICE3 solids, 1373 faces. Largest effort on the board.
- **The oTOF corpus** — unreachable from the converter (3 leaf solids seen, not 20 prototypes in
  62628 placements; `triangulate_asbbox()` dies on `Bnd_Box is void`). It is ~62560 *placed boxes*,
  the case placed primitives made cheap, so it would pay well.

## Next up, and the reason it is next

**Integration test: exercise the whole pipeline, assign materials, and simulate with the Geant
setup.** This is the next harness, and the point of it is that it *instruments things nothing else
can* — the step-length cost of a loose safety (item b), the real per-event cost of the surface solid
against a mesh, and assembly-level navigation, all under a workload that is the actual use case
rather than a benchmark of it.

## Raised by Sandro, 2026-08-22

### (f) The in-browser event display, layer 2: a live o2-sim loop

> *"What about **actually** simulating a part in the visualization in the browser with the o2-sim
> Geant4 backend and display some transported traces (either pixels on a screen behind the object)
> ... pointclouds within the object or things like this --> A sort of event-display. This will push
> **part testing** to a whole new level **and** we can sell it nicely. To have the impression of
> real-time, we might need to use MCStepLogger ... or hook into the O2HitMerger FairMQ channel."*

Layer 1 (batch replay from MCStepLogger, embedded in the website) is in
`Plan_Presentation.md` Track 3b. What is deferred here is the **live loop**:

- a local bridge (WebSocket/SSE) with a "fire N events" button spawning runs on demand;
- o2-sim **service mode** to keep geometry-initialized workers warm, so a shot costs the batch
  and not the Geant4 init — whether MCStepLogger flushes usably per batch in that mode is the
  gate to measure first;
- the **O2HitMerger FairMQ channel** as an alternative live tap: despite its name it sees **all
  MCTracks** (Sandro's correction, 2026-08-22), so it can feed track-level display (vertices and
  kinematics; line/helix rendering) without the step logger — less granular than steps, but
  already streaming in parallel mode.
