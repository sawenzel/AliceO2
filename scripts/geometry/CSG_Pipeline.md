# B-rep → CSG: the simplification pipeline — algorithms, prior art, and the staged plan

Date: 2026-08-01. Author: Claude (Fable 5). Companion to
[`CodeReview_Fable_v2.md`](CodeReview_Fable_v2.md) §5, which argues *why*; this document is *how*.

Status: **not previously planned.** Nothing in `BVHSurfaceSolid.md`, `ExactTrimTopology.md`,
`TolerancePolicy.md` or the v1 review proposes emitting boolean solids; the closest prior thought is
v1 §4.2's "semi-algebraic membership: a candidate hit is inside the face iff a small set of sign
conditions holds", which is the *local* form of the same idea, applied per face instead of per
solid. This document is the first plan for it.

Every capability claim about OCCT below was **verified on this machine** (OCCT 7.9.3, pythonOCC
7.9.0) by inspecting the installed headers and bindings; every number about ALICE3 was measured
from `ALICE_3_example/CAD_noETA.stp` this session. Nothing was run against the converter or the
kernel.

---

## 1. The problem, stated so the algorithm is obvious

A B-rep says: *here are the faces, here is where each one stops.* The "where it stops" is the pain —
it needs trim curves, and the trim curve of two quadrics is transcendental in either face's chart,
so the two faces disagree, forever, by an amount nobody can bound from inside one chart. Six
sessions of this branch are that sentence.

A CSG says: *here are the carriers, here are the signs.* A point is in the solid iff its sign
vector against the carriers matches. **There is no "where it stops", so there is nothing to
disagree about.** The transcendental curve is still there in the geometry — it is just never
*represented*, only implied by the conjunction of two sign tests, each of which is a polynomial
evaluation exact to machine precision.

That is the whole of the argument, and it is why the conversion is worth doing even though the
recognition side is harder than the trimming side.

---

## 2. Prior art — this is a solved problem with a name

**It is called B-rep → CSG conversion, and it has been production practice for fifteen years in a
neighbouring field.** Monte-Carlo neutron transport codes (MCNP, TRIPOLI, Serpent) accept *only*
CSG, so the fusion community had to solve exactly our problem — take vendor STEP, emit exact CSG,
at ITER scale.

**Theory.**

- **V. Shapiro & D. Vossler**, *Construction and optimization of CSG representations*, CAD 23(1),
  1991; and *Separation for boundary to CSG conversion*, ACM TOG 12(1), 1993. The foundational
  result: a solid bounded by natural quadrics admits an exact CSG representation over the
  halfspaces of its own faces, **provided the halfspace set is "separating"**; where it is not, they
  give a construction for the auxiliary halfspaces that make it so. They also treat *optimisation* —
  removing redundant halfspaces from each term — which matters because the naive form carries every
  face of the solid in every term.
- **Buchele & Crawford**, *Three-dimensional halfspace constructive solid geometry tree
  construction from implicit boundary representations*, CAD 36(11), 2004. The practical 3D
  construction.
- The relevant *negative* result to know: exact CSG conversion is worst-case exponential in the
  number of halfspaces, and undecidable in general for free-form carriers. Both are handled by
  bounding the input (Section 6) rather than by cleverness.

**Practice — and this is the part worth copying.**

- **McCAD** (KIT, ~2010-) — CAD → MCNP. The first to make the decomposition approach work on
  reactor-scale assemblies.
- **GEOUNED** (UNED / F4E / CERN, open source, 2023-) — STEP → MCNP / OpenMC / Serpent / PHITS.
  **Built on OCCT via FreeCAD**, i.e. on the identical kernel this project already installs and
  already drives from Python. It handles planes, cylinders, cones, spheres and tori — our exact
  five families minus none. It is the closest existing thing to what we want, it is readable, and
  its published papers describe the decomposition heuristics that took years to find. Reading it
  before writing ours is worth several days.

**What they converged on** is not pattern-matching and not search. It is one algorithm:

> **Decompose by extending the carriers.** Find a face whose carrier witnesses a non-convexity;
> cut the solid with that carrier extended to infinity; recurse on the pieces. Terminate when every
> piece is convex. A convex piece is *by definition* the intersection of the oriented halfspaces of
> its own faces. The solid is the union of the pieces:
>
>     S  =  ⋃ᵢ ⋂ⱼ Hᵢⱼ
>
> A two-level disjunctive normal form over oriented quadric halfspaces. No trims, no pcurves, no
> wires, no shared edges, no closure check, no tolerance policy.

The termination argument is the useful part: **a solid is convex iff it has no concave edge**, and
every split at a concave edge's carrier strictly reduces the concave-edge count of at least one
piece. So the loop is driven by an explicit, cheap, local test rather than by a global search.

---

## 3. The core primitives, each of which is small

### 3.1 Concavity test (the engine)

For an edge shared by faces `F₁`, `F₂` with outward normals `n₁`, `n₂` at a point on the edge, and
edge tangent `t` oriented along `F₁`:

    convex   if  (n₁ × n₂) · t  ≥ 0
    concave  otherwise

evaluated at a few parameters along the edge (a curved edge can change character). Tangential edges
(`n₁ × n₂ ≈ 0`) are neither and are skipped — they are also where the whole method is delicate.

A solid with zero concave edges is convex; its CSG term is `⋂ Hⱼ` over its faces' carriers, with
each halfspace oriented by the face's `same_sense`/`REVERSED` flag. **This is ~50 lines with
`BRepGProp_Face::Normal` and `BRepAdaptor_Curve`.**

### 3.2 Splitting

`BRepAlgoAPI_Splitter` — **verified present in the pythonOCC bindings on this machine.** Feed the
solid as argument and the extended carrier as tool; take the resulting pieces. `BOPAlgo_Splitter`
is the lower-level entry point if the API one proves fragile.

Splitting is the risky operation: OCCT's boolean core is the least robust part of the kernel, and
tangential or near-coincident carriers are its classic failure mode. Mitigation is structural —
every failure falls back one tier, never propagates (Section 5).

### 3.3 Carrier extraction and canonicalisation

The converter already recognises plane/cylinder/cone/sphere/torus at machine precision and this is
credited as sound (v1 §2). Two additions are worth making regardless of CSG:

- **`ShapeAnalysis_CanonicalRecognition` — verified present in OCCT 7.9.3 and exposed in pythonOCC
  on this machine** (`IsPlane`, `IsCylinder`, `IsCone`, `IsSphere` with a tolerance). This is
  OCCT's own recogniser, it accepts a *set* of faces as well as one, and it is the natural
  cross-check for — and possibly replacement of — the hand-rolled recogniser that is currently
  duplicated between `O2_CADtoTGeo.py` and `analyze_surface_geometry.py` with a subtle policy
  difference (v1 §7). Use it as a second opinion first; adopt only if it agrees.
- **Swept-surface canonicalisation** (see §4, Tier 0): `SURFACE_OF_LINEAR_EXTRUSION` of a line is a
  plane, of a circle is a cylinder; `SURFACE_OF_REVOLUTION` of a line is a cone (or cylinder, or
  plane), of a circle is a torus (or sphere). ALICE3 has 85 such faces that are currently
  unsupported and are almost all trivially canonical.

### 3.4 Redundant-halfspace removal

A convex cell built from its own faces carries one halfspace per face, and many are implied by the
others. Removing them matters for both speed and readability. Cheap sound test: a halfspace `H` is
redundant in cell `C` iff no face of `C` lies on `H` — which is already known from the
construction — plus a sampling check. Shapiro's optimisation is the rigorous version; start with
the cheap one and measure.

### 3.5 Acceptance — the reason this is feasible now

    volume( candidate − original ) + volume( original − candidate )  ==  0   (to model tolerance)

computed with `BRepAlgoAPI_Cut` both ways and `BRepGProp`. **This is a proof of equality, not a
sample-based argument** — nothing on the trimmed-surface path has ever had an acceptance test this
strong. Plus, as a second and independent gate, the shape goes through `runOracleGate.py` unchanged:
the gate scores a `TGeoShape` and does not care how it is implemented.

So the pipeline can be **greedy, heuristic and incomplete on the proposal side, because acceptance
is exact.** That asymmetry is the design.

---

## 4. The tiers — four, not three, and Tier 0 is free

Ordered by (value ÷ risk). Each is independently shippable and independently verified.

**Tier 0 — swept-surface canonicalisation.** Convert `SURFACE_OF_{LINEAR_EXTRUSION,REVOLUTION}` of
a line/circle into the equivalent plane/cylinder/cone/torus. Pure converter work, no new kernel, no
CSG. ALICE3: 85 faces across 8 solids, of which 4 solids become quadric-only. Half a day.

**Tier 1 — whole-part primitive recognition.** The part's entire face set matches one ROOT
primitive: `TGeoTube(Seg)`, `TGeoCone(Seg)`, `TGeoSphere`, `TGeoTorus`, `TGeoBBox`, `TGeoPcon`,
`TGeoPgon`, `TGeoArb8`. Detection is a closed-form match against the recognised carrier multiset;
emission is one TGeo constructor. **No new C++ at all**, and the output is exact, fast,
ROOT-native, Geant4-native and visualisable.

**Tier 2 — recognised boolean of primitives.** `tube − tube` (the Bagger junction), `tube ∪ tube`,
`primitive − N holes`, `primitive ∩ primitive`. Emission: `TGeoCompositeShape` over `TGeoBoolNode`
— **but measure first** (§7), because ROOT's boolean navigation may be the wrong target and the
alternative (§5) is already needed for Tier 3. Note `TGeoBoolNode` accepts *any* `TGeoShape`, so an
`O2BVHSurfaceSolid` can be a leaf: "one awkward body minus four bolt holes" is expressible with the
awkward body kept in the existing representation. That hybrid is underrated and may be the highest
value form on real assemblies.

**Tier 3 — general convex decomposition** (§2). The general answer for quadric-bounded solids, and
where the research risk lives. Do not start it until Tiers 0-2 are measured, because they may leave
little for it.

---

## 5. The emission target: `O2CSGSolid`

For Tier 3 — and possibly Tier 2, pending §7 — do **not** emit a deep `TGeoCompositeShape`. A DNF
with 20 cells is 20 nested boolean nodes, each re-transforming the point and re-dispatching
virtually.

Emit a flat solid instead: a `std::vector<Cell>`, each cell a contiguous run of oriented quadric
coefficient blocks. Then

- `Contains` = `OR` over cells of `AND` over sign tests — a few dozen FLOPs, no traversal, no
  allocation, no BVH.
- `DistFromOutside` / `DistFromInside` = per cell, intersect the ray's parameter interval with each
  halfspace's inside-interval (Cyrus-Beck generalised to quadrics: a plane gives a half-line, a
  quadric a root pair), then min/max over cells. **Every root solver this needs already exists and
  is already validated in `BoundedSurface.h`** — this is reuse, not new mathematics.
- `Safety` = min over cells of min over halfspaces of the exact point-to-quadric distance, already
  implemented.
- `Capacity` = closed form for Tier 0-2; divergence theorem over cell boundaries for Tier 3, or
  simply the OCCT value recorded at conversion time and asserted at load.

Two further properties worth stating because they are the strategic payoff:

- **Watertight is not a property that can fail.** There is no surface set, so there is no closure
  check, no rim matching, no join tolerance, no `kBSplineFlatness`, no sagitta band. The entire
  §3.3/N3 problem class does not exist in this representation.
- **It is the natural AOT-codegen target** (the direction recorded in v1 §10 Phase 4 and the user's
  stated preference). A cell list of fixed length with fixed coefficients generates as straight-line
  code with no virtual dispatch and no heterogeneous loop — far more so than a trimmed-surface
  solid.

---

## 6. Where this does *not* work — measured, not guessed

`CAD_noETA.stp`, measured this session by direct STEP parsing (55 solids, 9266 faces):

| carrier | faces | share |
| --- | --- | --- |
| plane | 3321 | 35.8% |
| cylinder | 2724 | 29.4% |
| cone | 409 | 4.4% |
| torus | 350 | 3.8% |
| **analytic total** | **6804** | **73.4%** |
| B-spline (with knots) | 1259 | 13.6% |
| rational NURBS (complex entity) | 1118 | 12.1% |
| **free-form total** | **2377** | **25.7%** |
| surface of linear extrusion | 73 | 0.8% |
| surface of revolution | 12 | 0.1% |

Per solid:

- **15 solids are bounded by natural quadrics only** → CSG-eligible in principle.
- **4 more are quadrics + swept only** → CSG-eligible after Tier 0.
- **36 contain free-form surfaces** → not CSG-convertible at the top level, at any tier.

Size distribution: 8 solids ≤10 faces, 17 ≤40, 14 ≤100, 11 ≤400, **5 >400** (2034, 1052, 996, 965,
720 faces).

**The uncomfortable conclusion, stated plainly: on ALICE3 the CSG route adds essentially no
coverage.** The 15 quadric-only solids are, to within the recognition rescues, the same 15 the
exact-surface converter already handles (`BVHSurfaceSolid.md` records "coverage 15/55 already
analytic"). What CSG buys *there* is quality — exactness at the seams, no closure question, less
memory, more speed, portable output — not new parts.

Three consequences that should shape the whole programme:

1. **ALICE3 coverage is gated by free-form *surfaces*, not by trims and not by CSG.** 36 of 55
   solids and 2377 faces. The B-spline-surface milestone that v1 deprioritised as "largest effort,
   do last" is, on this corpus, *the* coverage lever. That deprioritisation was made when the
   working assumption was that trims were the blocker; the measurement says otherwise.
2. **Three of the five biggest solids are quadric-only and huge** (1052, 965, 720 faces). Those are
   above any sane Tier-3 decomposition budget and should never be attempted as one — they are
   exactly what the BVH solid is good at. **The two representations are complementary in precisely
   the right way: CSG scales badly in face count and beautifully in boolean complexity; the BVH
   solid is the reverse.** Gate Tier 3 at ≤40 faces (17 ALICE3 solids qualify, plus 8 tiny ones)
   and by a cell-count cap.
3. **The hybrid band is real and is the ALICE3-specific opportunity**: `#33027` has 996 faces of
   which 69 are free-form (7%); `#4411` has 245 with 40; `#123` has 50 with 4. Those cannot be pure
   CSG and cannot be pure quadric B-rep — they want a quadric body with a small free-form patch
   set, i.e. exactly the Tier-2 hybrid leaf.

Where CSG *is* dominant, by contrast: `Bagger.step` (12/13 parts, every seam quadric∩quadric),
`oTOF` (100% planar), the beampipe (§ companion doc), and mechanical/engineering CAD generally.
Those are also the models where the current representation is failing today. So the tiering is
right; the coverage expectation on ALICE3 specifically must be honest.

Further limits, for the record: variable-radius blends and lofts are free-form and out (constant
fillets are torus/cylinder halfspaces and are fine). Tangential carriers are the decomposition's
robustness cliff. Acceptance is against OCCT's tolerant model, so "equal" means equal to model
tolerance — the same bar as today, neither better nor worse.

---

## 7. The measurement that must come before the code

**Recognition rate before emission.** One Python session, no emission, no C++: walk Bagger, as1,
oTOF and all 55 ALICE3 solids; report per solid whether it is quadric-only, its face count, its
concave-edge count, and — for the quadric-only ones — whether a Tier-1 template matches. That
single table decides how much of Tiers 1-3 is worth building, and it is a few hours. **Do not build
the emitter before this table exists.**

**`TGeoCompositeShape` timing.** Build `cyl_inter_cyl` by hand as a composite and run it through
`runSolidHarness` against `O2BVHSurfaceSolid` and `O2Tessellated`. One afternoon; decides Tier 2's
emission target empirically. Use a dedicated fixed-point loop, not the gate's timing column, which
"cannot resolve anything below a few per cent" and moved 3× between identical runs (`NEXT.md`).

---

## 8. Plan — Stream A

Each step ends in a committable state with its own gate. Steps 1-2 are the decision point.

1. **Recognition census** (§7). Deliverable: `scripts/geometry/csg/census.py` + a table in this
   document. *Gate:* the table exists and the numbers are reproducible.
2. **Acceptance harness** — `csg/accept.py`: symmetric-difference volume between a candidate CSG
   (built in OCCT) and the original solid. *Gate:* self-test on hand-built pairs, including a
   deliberately-wrong candidate that must be rejected.
3. **Tier 0** — swept-surface canonicalisation in the converter. *Gate:* the 85 ALICE3 swept faces
   are recognised; existing conversions bit-identical; `--surface-report` numbers move only where
   expected.
4. **Tier 1** — primitive recognition + TGeo emission. *Gate:* every G1 ladder fixture that is a
   primitive is recognised and accepted; symmetric-difference volume ≤ model tolerance; the oracle
   gate passes on the emitted shape; census table updated with the achieved rate.
5. **`TGeoCompositeShape` timing probe** (§7) → decide Tier 2's target.
6. **Tier 2** — `tube − tube` and `primitive − N holes`. *Gate:* the six Bagger cylinder parts and
   `tube_window` convert as booleans and pass both acceptance tests. **This is where G3 becomes
   12/12 and the tube-tube problem is finished.**
7. **Tier 3 spike** — time-boxed, ≤40-face solids only, `BRepAlgoAPI_Splitter` + concavity test,
   emitting the DNF as data (no C++ yet). *Gate:* the cell-count and success-rate table over the
   eligible ALICE3 and Bagger solids. Decide `O2CSGSolid` on that table, not before.

Explicitly **not** in this stream: `O2CSGSolid` itself (a C++ stream, gated on step 7),
`O2RevolvedSolid`/`O2ExtrudedSolid` (separate stream — they are exact *solids*, not CSG), and any
change to the existing trim machinery.

---

## 9. Relation to `BVHSurfaceSolid` — substrate, not rival

Three clarifications, because "CSG versus BVH" is the wrong axis and reasoning on it leads to the
wrong build order.

### 9.1 `TGeoBox`/`TGeoTube` recognition is *not* a consequence of CSG — it is a prerequisite, and it stands alone

Tier 1 needs **no boolean algebra at all**. "These six planes are pairwise parallel, mutually
orthogonal, and their signed offsets are ±(dx, dy, dz)" is a template match on the carrier multiset
the converter already recognises at machine precision. Same for `TGeoTube` (two coaxial cylinders +
two planes ⊥ the axis), `TGeoTubeSeg` (+ two planes through the axis), `TGeoCone`, `TGeoTorus`,
`TGeoPcon` (coaxial cones/cylinders/planes, ordered by z), `TGeoArb8`.

So Tier 1 ships **before** any CSG machinery exists, carries none of Tier 3's robustness risk, and
its output is the most portable thing in the whole programme: exact, ROOT-native, Geant4-native,
visualisable, overlap-checkable, ~40 bytes. It is listed inside this document only because it shares
the census and the acceptance harness with the tiers above it — not because it depends on them.

The only reason it is step 4 and not step 1 is that the acceptance harness (step 2) must exist
first. Recognising a tube is easy; *proving* the recognised tube equals the CAD part is the part
that must not be skipped, and it is the same `BRepAlgoAPI_Cut` symmetric-difference test every other
tier uses.

### 9.2 Which is faster is a real question with a face-count answer, and it must be measured

The intuition that `BVHSurfaceSolid` may outperform CSG is not wrong — it is *conditionally* right,
and the condition is face count.

| | flat CSG DNF | `BVHSurfaceSolid` |
| --- | --- | --- |
| `Contains` | `OR` over cells of `AND` over sign tests. A point outside must be rejected by **every** cell (the `AND` short-circuits, the `OR` does not), so cost is ~O(cells), no spatial pruning | ray parity: BVH descent O(log N) + a ray/surface root solve per candidate + a point-in-trim per hit, and a wire trim walks a polyline |
| distances | per-cell interval clipping over **all** cells — no pruning at all in the flat form | ray-BVH traversal prunes hard; tmax tightening prunes further |
| build | none | sidecar parse + BVH build + `CloseShape` |
| scaling | **badly in face count** (cells grow with the arrangement), beautifully in boolean complexity | **badly in boolean complexity** (every seam is a trim curve someone must fit), beautifully in face count — O(log N) |

So: for a 6-face box or a 3-primitive junction, CSG wins by a wide margin — no traversal, no trims,
no allocation, no `CloseShape`. For oTOF's 1505-face planar tile, CSG is hopeless and the BVH solid
is the natural answer. The crossover is somewhere in between and **nobody should guess where**; it
is a measurement, and the harness already times any `TGeoShape`.

The important structural note: the flat DNF's lack of pruning is fixable, and the fix is the thing
already built — **a BVH over cell AABBs**. Once that is done, the two representations have converged
architecturally: one traversal, one set of quadric root solvers, two kinds of leaf. Which is the
next point.

### 9.3 The end state is one navigator with several leaf kinds, and a per-part selection policy

The right mental model is not "CSG solid vs surface solid". It is:

- **A substrate** — the BVH, the exact quadric root solvers, the `Curve2D` 2D machinery, the
  harness, the oracle gate, the reliability reporting. All of it is already built, validated, and
  shared.
- **Several leaf kinds** on top of it: a trimmed analytic patch (today), a halfspace cell (Tier 3),
  an exact 2D profile swept or revolved (Stream D), a free-form patch (Stream B). Plus whole-part
  shortcuts that bypass the substrate entirely when they apply: a ROOT primitive (Tier 1), a shallow
  boolean of primitives (Tier 2).
- **A selection policy** that picks per part.

That policy is what "the best possible TGeo representation for a CAD part, however the mechanism"
actually means operationally, and it should be explicit and measured rather than implicit in the
order of `if` statements:

1. Try each applicable representation, cheapest and most portable first: primitive → boolean of
   primitives → revolved/extruded profile → exact surfaces (BVH) → free-form (BVH) → tessellated.
2. **Accept on evidence, not on recognition**: symmetric-difference volume for the constructive
   ones, the oracle gate for all of them, `NavigationReliability` for the surface ones.
3. Among those that pass, choose by *measured* per-query cost on a fixed sample set — not by tier
   rank. A recognised `TGeoPcon` and a 200-face BVH solid may both be exact; the pcon is better, and
   a case will eventually exist where it is not.
4. **Record the decision and the runners-up** in the manifest, so coverage is reported as a tiered
   scorecard with reasons rather than as a single fraction (Stream E step 6).

Nothing in this document argues for retiring `BVHSurfaceSolid`. It argues for taking its hardest
inputs away from it — the quadric-quadric seams it cannot represent without approximation — and
letting it do what it is uniquely good at: many faces, arbitrary trims, and, once Stream B lands,
the free-form 26% of ALICE3 that no constructive representation can reach.
