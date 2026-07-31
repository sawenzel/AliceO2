# BVHSurfaceSolid — comprehensive review, conceptual analysis, and reset plan

Date: 2026-07-31. Reviewer: Claude (Fable 5), on branch `swenzel/bvhsurfacesolid`.
Method: full read of `BVHSurfaceSolid.md`, `ExactTrimTopology.md` and the code; four parallel
deep-review passes (kernel `BoundedSurface.h`; solid/IO/tests; converter `O2_CADtoTGeo.py`;
OCCT-oracle feasibility on this machine); plus targeted manual verification of the parity and
trim-classification code paths. No code was changed.

This document answers the questions posed for this session:

1. Why does the code fail on "easy" examples (Bagger tube-tube combinations) although all the
   surfaces involved are plain cylinders and planes?
2. Can such cases be handled *exactly* at all — and if not exactly, what is the correct target?
3. What do the concepts (patches, trims, curves, Bézier/B-spline, BREP solids) require, with
   OpenCascade as the reference?
4. Should we build an OCCT-based reference solid first?
5. What are the new success criteria, and what does a proper test-driven plan look like?

---

## 1. Verdict (executive summary)

**The kernel and the per-surface mathematics are in good shape; the architecture of *trimming*
is what fails, and it fails for a fundamental reason, not an incidental one.** The project's
"easy" failures (tube-tube intersections in `Bagger.step`) are in fact the mathematically hard
part of the whole problem: the intersection curve of two cylinders is *transcendental in each
cylinder's parametric chart*, so the per-face 2D trim curves ("pcurves") that both the CAD file
and our sidecar carry are inherently independent approximations. Two adjacent faces therefore
never agree exactly about their common boundary, the solid can never pass the watertight closure
check, and one fixed-direction ray-parity `Contains` turns every seam sliver into an error
*shadow* that extends arbitrarily far from the boundary. Everything observed on Bagger — the 367
residual `unexplained` points, "no part closes", the 1.7 cm containment errors before the
closed-B-spline fix — is this one design decision playing out.

The good news is threefold:

- For **analytic∩analytic edges** (cylinder-cylinder, cylinder-plane, cone-plane, …, i.e. all of
  Bagger and as1 and much of ALICE3) an *exact* trim representation exists — just not as a 2D
  curve in the parametric domain. Represent the trim boundary as "the neighbor surface itself"
  and every membership question becomes a closed-form root/sign computation (Section 4.4). Both
  faces then agree *by construction*. This is the single most important architectural change to
  make.
- For genuinely free-form edges, **consistency, not exactness**, is the achievable and sufficient
  goal (parity only needs both faces to agree on the boundary), via shared-edge sampling — the
  route already sketched in `ExactTrimTopology.md` item 1, which was correctly diagnosed as "not
  the cause of the observed bug" but remains the right *data model*.
- Independently of the data, single-ray parity must become **self-checking** (ambiguity band on
  trim hits + re-shoot in a different direction). This is cheap, catches every seam problem at
  query time instead of silently misclassifying, and is what OCCT, Geant4 and VecGeom all do in
  some form.

A trusted oracle is missing: today's reference (`O2Tessellated` within a "mesh band") is too weak
to certify exactness — worse, the review found two soundness holes in the harness's band
classification that *systematically hide* the worst bug classes (a candidate that misses an
entire wall is always "explained by mesh chording"), so the residual "367 unexplained" is a
floor, not a ceiling (Section 6, S9). OCCT 7.9.3 + pythonOCC are installed and working **on this
machine**; a Python-side OCCT oracle feeding the existing harness is ~2 days of work and should
come first, because every subsequent claim is judged by it.

Two further findings deserve headline status: **ROOT persistence is actively harmful** (a
streamed-and-read solid comes back as an empty shape that reports itself `Reliable` — S1), and
the documented "order-dependent clustering" explanation for the non-manifold `Contains`
disagreements is **refuted by the code** (the clusterer is deterministic; the hit multisets must
differ between BVH and loop traversal — S6). Both change what the next debugging session should
do.

The recommended reset is therefore *not* "go back and rewrite" but: freeze coverage work, build
the oracle and a synthetic Boolean fixture ladder (tube-tube first), fix the identified kernel
bugs, then change the trim representation for analytic-analytic edges. Sections 9-10 define the
gates and the order.

---

## 2. What is solid and should be kept

Credit where due — the following parts were reviewed and found sound, some of them very carefully
engineered:

- **Analytic surface kernels**: ray/quadric intersections with stable handling of tangent double
  roots (parity-safe suppression), exact divergence-theorem capacities (incl. the torus closed
  form), oriented normals, the frame conventions. Verified consistent across the five families.
- **The Curve2D/CurveWire winding machinery for lines and arcs**: v-monotone arc splitting,
  half-open canonical-endpoint seam convention, Green's-theorem areas. This is genuinely robust
  2D computational geometry.
- **BVH integration**: conservative AABBs with documented tolerance and outward float rounding,
  `max_leaf_size = 1` with a measured justification, tmax-tightening with a documented pruning
  argument, `_Loop` twins bit-identical over 85k rays.
- **The measurement culture**: `SolidNavigationHarness`, the part DB, fixed seeds, the honest
  handoff notes, and `NavigationReliability`/fail-loudly. The 2026-07-26 root-cause session
  ("diagnose before planning") is the right template for everything below.
- **Canonical-form recognition** (surfaces + trim curves) in the converter: model *selection* at
  machine precision with a bit-identical-geometry acceptance criterion. The as1 lesson ("the
  stored type is a statement about the exporter, not the geometry") is important and correctly
  institutionalized.

The test suite (36 cases) is real and green, but it validates almost exclusively *synthetic solids
built through the C++ API against ROOT primitives* — it never sees the failure class that
dominates real CAD (independently approximated adjacent trims). That is why everything passed
while Bagger was wrong by centimetres.

---

## 3. BREP anatomy — the concepts, precisely

Terms used throughout (OpenCascade vocabulary, since it is the reference implementation):

- **Surface**: an untrimmed analytic or NURBS carrier, a map S(u,v) → R³ (infinite plane, full
  cylinder, …).
- **Face**: a surface plus a finite region of its (u,v) domain.
- **Wire**: a closed loop of edges bounding that region (one outer, holes inner).
- **Edge**: a *topological* entity shared by exactly two faces in a closed manifold shell. It
  carries: a 3D curve C(t) → R³, a tolerance, and **per adjacent face** a *pcurve*
  P_F(t) → (u,v)_F, the edge's image in that face's parametric domain.
- **Shell/Solid**: an oriented, closed collection of faces; watertightness is a *topological*
  statement (every edge has exactly two incident faces), not a geometric one.
- **The crucial fine print**: C(t) and the pcurves are *independently fitted* representations.
  OCCT guarantees only that S_F(P_F(t)) stays within the edge tolerance of C(t)
  ("SameParameter"). The three curves of one edge (3D + two pcurves) never coincide exactly
  except in special (canonical) cases. A BREP is exact in its surfaces and *tolerant* in its
  boundary curves — by design, in every CAD kernel.

Why must it be so? Because closed-form parametric representations of intersection curves do not
exist in general:

**The tube-tube computation.** Cylinder C1 (axis z, radius r), cylinder C2 (axis x, radius R),
r < R. A point of C1 is (r cos φ, r sin φ, h). Membership on C2 requires
(r sin φ)² + h² = R², so the intersection curve in C1's own chart is

    h(φ) = ± sqrt(R² − r² sin²φ)

This is not a line, not an arc, not a polynomial or rational B-spline in (φ, h) — no exact
`Curve2D` exists for it. (The 3D curve is a degree-4 algebraic space curve; its image in an
*angular* chart is transcendental.) The same holds for cylinder-cone, cone-cone, anything-torus,
and — the Bagger `Bucket` case — a plane cutting a cylinder viewed from the *cylinder's* chart
(from the plane's chart it is an exact ellipse). So:

> **The trim curves of the "easy" Boolean combinations are exactly the objects that have no exact
> representation in the parametric-wire vocabulary.** The current failure is not a bug in the
> B-spline code; it is the data model reaching its mathematical ceiling.

This directly answers the session's fear ("if this cannot be done exactly … probably nothing can
be done"): it *can* be done exactly — but only by changing what a trim *is* (Section 4), not by
improving the curves.

---

## 4. What "exact" can mean — three regimes and the way out

### 4.1 Regime E1 — exact carrier surfaces (achieved)

Ray ∩ surface, normals, capacities are exact for plane/cylinder/cone/sphere/torus. This is the
real win over tessellation and it is done and validated. Nothing below touches it.

### 4.2 Regime E2 — analytic∩analytic edges: exactness is possible, but not with pcurves

Where the trim boundary of face F (carrier S1) is "the intersection with neighbor surface S2",
membership can be decided without any 2D curve:

- The 2D winding test in F's chart casts a scanline (say h = h0, φ increasing). On the carrier
  this scanline is a *3D circle* (or generator line, meridian, …). Its crossings with the trim
  boundary are exactly the solutions of {scanline} ∩ S2 — for quadric S2 a **quadratic** (torus:
  quartic) with closed-form roots. Count crossings with the right branch/interval bookkeeping and
  the point-in-face answer is exact to machine rounding.
- Equivalently, and even simpler for ray hits: a candidate hit p on S1 is inside the face iff a
  small set of *sign conditions* f_S2(p) ≶ 0 … holds (semi-algebraic membership). For locally
  convex trim regions this is just "AND over neighbors"; in general the winding formulation above
  is the robust one.
- **Adjacency is the representation.** Store per trim edge: (neighbor surface, parameter
  interval, branch selector) instead of a fitted 2D curve. The two faces sharing the edge
  reference the *same* neighbor surface, so their boundaries agree **identically** — watertight
  by construction, which no join tolerance can ever achieve. The closure check degenerates into a
  structural check on edge identity (as in a real BREP kernel).

STEP delivers the required adjacency (each `TopoDS_Edge` knows both faces; the converter already
has `TopExp.MapShapesAndAncestors` on its shelf), and the surface pair for each edge is known
after canonical recognition. Applicability on the current examples:

| model | blocking edges | covered by E2? |
| --- | --- | --- |
| Bagger tube-tube seams | cyl∩cyl (transcendental pcurve) | **yes — exact** |
| Bagger `Bucket` ellipses | plane∩cyl | **yes — exact** (from the plane chart it is also an exact arc/ellipse) |
| as1 (after recognition) | cyl∩plane | **yes — exact** |
| ALICE3 non-iso recognized trims (373 faces) | mostly quadric∩quadric/plane | **yes — exact** |
| genuinely free-form edges (B-spline surface involved) | — | no → regime E3 |

This is the single highest-leverage change available: it makes the "easy" examples *actually
exact*, closes the watertightness gap for them, removes the polyline flattening from navigation,
and shrinks the sidecars (an edge record becomes a surface reference + interval).

Cost assessment (honest): a new trim-edge kind in the kernel (`Curve2D` gains an implicit-curve
variant, or better: the trim domain object gains "edge = carrier ∩ neighbor" entries), scanline
crossing code per carrier/neighbor pair (closed-form; quadric pairs are quadratics, torus pairs
quartics — all solvers already exist in the file), converter emission of neighbor-surface records
+ sidecar format bump. Estimate 1-2 weeks of focused work including tests — comparable to what
the B-spline trim milestone took, with a far larger correctness payoff.

### 4.3 Regime E3 — free-form edges: consistency, not exactness

Where a B-spline *surface* is involved, no exact trim exists in any vocabulary; the correct goal
(as `ExactTrimTopology.md` already concluded) is that both faces derive their trims from **one
shared sample set of the single 3D edge curve**, with endpoints pinned to the shared vertex
coordinates. Then the residual disagreement is a *known, coherent* chord error instead of an
uncontrolled fitting difference, and the closure check can identify paired half-edges by ID with
a declared epsilon. `ExactTrimTopology.md` item 1's four-step recipe (edge map → shared samples →
per-face closed-form projection → polyline emission) stands; its "do not start" verdict was about
*that bug*, not about the data model — with E2 in place, E3 is the fallback for the shrinking
remainder.

### 4.4 Robustness — parity must be self-checking regardless of data

Verified in code: `pointInTrim` is a hard binary (`BoundedSurface.h:2402-2406` and siblings),
`oddCrossingParity` (`O2BVHSurfaceSolid.cxx:101-127`) clusters near-equal t only, and `Contains`
uses one fixed direction (`kContainsTestDirection`) with no retry. Consequences: any hit landing
within ε of a trim boundary is classified by floating-point coin flip; a lost crossing
misclassifies the *entire shadow* of the seam; nothing at query time notices.

The standard remedy (OCCT `BRepClass3d_SolidClassifier` re-shoots on degenerate hits; Geant4 and
VecGeom use surface-tolerance bands):

1. Each patch reports, per hit, whether the (u,v) hit point lies within a declared band of the
   trim boundary (the 2D distance machinery already exists — `closestPoint`/`distanceSq`).
2. If any hit in a query's list is boundary-ambiguous, or an unmatched near-cluster appears,
   **re-shoot along a different direction** (2-3 attempts, then decide by `Safety`).
3. Result: seam defects cost a bounded retry near edges instead of silent wrong answers far away.

This should be implemented *before* the data-model work, because (i) it very likely removes most
of the residual 367 `unexplained` points immediately, (ii) it converts future converter/kernel
regressions from silent to loud, and (iii) it is what makes the solid safe to use in production
even while coverage grows.

### 4.5 What OCCT itself does (and what to copy vs. not)

- OCCT never tries to make boundary curves agree exactly; it makes tolerance a first-class,
  per-edge/per-vertex attribute and requires all classification to be tolerance-aware. **Copy the
  concept**: the sidecar must carry the model tolerance (today it carries none — the C++ closure
  check literally cannot know what ε to glue with), and per-edge identity.
- OCCT's classifier (`BRepClass3d`) is a ray classifier with per-face 2D classification, ON
  states, and ray restart. **Copy the algorithm shape** (Section 4.4), not the implementation —
  it is ms-per-call, allocation-heavy and effectively single-threaded, which is precisely why
  this project exists.
- OCCT's `ShapeAnalysis`/`BRepCheck` machinery is the model for our closure diagnostics: verdicts
  based on topology + tolerances, not on exact float equality of resampled polylines.

---

## 5. Findings — kernel (`Detectors/Base/src/BoundedSurface.h`, 3894 lines)

Confirmed/strongly-suspected defects, ranked. (Line numbers from the current branch state.)

| # | Where | Defect | Consequence |
| --- | --- | --- | --- |
| K1 | `Curve2D::valid()` 1000-1038 vs `startPoint/endPoint` 1058-1077 | Knot-vector **clamping assumed, never validated**; a periodic/unclamped B-spline (exactly what OCC writes for tube-tube curves before `SetNotPeriodic`) returns off-curve poles as endpoints | wire falsely `Open` → whole face silently dropped (see K7); or off-curve canonical endpoint corrupts winding |
| K2 | `buildCurveTrim` 1681-1684 | Full-turn trim rejection uses the **pole hull** span, not the curve span; a closed intersection curve with overshooting poles reads as "> 2π" | legitimate through-hole host faces rejected → face dropped |
| K3 | `kWireJoinTolerance` (line 63) applied to **mixed-unit** (φ[rad], h[cm]) distances | closure acceptance depends on radius: small holes falsely open, big cylinders accept real 3D gaps; the known `ST1829909_01` loader rejection is this | faces dropped / bad wires accepted; needs a per-domain metric |
| K4 | `bsplineSampleRecursive` 964-986 | The closed-curve fix is not fully robust: degenerate-chord recursion probes only the parametric midpoint | a slit-like/self-touching span can still collapse; same failure class as the historical bug |
| K5 | `CurveWire::classify` 1517-1519 | Boundary band 1e-9 tested against a polyline that is only ~1e-5 accurate; also winding and distance use two slightly different polylines (canonical endpoints vs raw cache) | `Boundary` state unreachable for B-spline trims; in/out flips within ±1e-5 of the true curve |
| K6 | quadratic roots 2488/2815/3150; cone branch switch 3139; torus quartic 660-732 | Naive `(-b±√)/2a` (cancellation); cone degeneracy switch and torus Ferrari tolerances are **absolute** and scale-dependent; 2 Newton steps can't rescue near-double roots | wrong/lost roots on near-axis cone rays and toroidal fillets → parity flips (live suspect for part of the 367) |
| K7 | `O2BVHSurfaceSolid.cxx` `Add*Surface` (e.g. 349-352, 471-474) | A face that fails to build is logged and **silently omitted** from the parity solid | one dropped face = wrong `Contains` in its entire shadow; cheapest hypothesis to check against the 367 |
| K8 | `oddCrossingParity` 108-125 + `sameIntersection` 222-226 | **Transitive** clustering with a relative window: N hits each 1e-7·t apart chain into one cluster spanning N·1e-7·t | thin features at large t merge into "mixed" clusters → even parity → misclassification growing with distance |
| K9 | closure half-edge check (sampling conventions 570-574, 1698-1706, 3842) | Chord *counts* are synchronized but sample *phases* only match for quarter-turn-related frames; B-spline edges resampled per-face | "no real CAD part closes" is **structural**, not a tuning issue — the check cannot pass on imported CAD as designed |
| K10 | `CurveWire::initialize` 1437-1473 | no self-touch/self-intersection validation (the polygon path has it, 377-384) | meaningless signed area on dirty input → wrong orientation "normalization" |
| K11 | sphere `directionInTrim` 2766-2770 | pole over-accept with wire trims | minor, measure-zero |
| K12 | polygon wires join with `kToleranceSq` (417) but curve wires with `kWireJoinToleranceSq` (1453) | same extractor feeds both | inconsistent acceptance |

**Approximation inventory** (everything that enters navigation, beyond the CAD data itself):
B-spline→polyline flattening (the effective trim boundary, ~1e-5 param units); canonical seam
endpoint substitution (≤1e-5 shifts near every seam); pole-hull unwrap windows; tangent-pair
suppression window scaling with |t|; `kRayTolerance` 1e-9 lower cut; 1e-9 on-surface bands that
pretend more precision than the data has. Capacity quadrature and display meshes are correctly
quarantined. **Architectural debts**: no shared-edge entity (the root cause), no exact
point-in-trim (Bézier clipping/root isolation would retire the polyline from navigation), no
coherent tolerance policy (four absolute constants across three incompatible unit systems), one
3900-line header with ~⅓ duplication across the five surface classes (a `TrimDomain2D` value type
plus a surface-of-revolution base would remove most of it).

---

## 6. Findings — solid / IO / tests / harness

| # | Where | Defect | Consequence |
| --- | --- | --- | --- |
| S1 | `Streamer` cxx:1184-1192, `//! fImpl` | **ROOT persistence is broken and actively harmful**: a deserialized solid comes back empty, `CloseShape(false)` then *zeroes* the streamed bbox, and an empty `ClosureReport` defaults to `closed=true` → the husk reports `NavigationReliability::Reliable` | any `TGeoManager::Export/Import` silently converts every solid into an "authoritatively reliable" empty point. Fix first: stream the sidecar-record blob and replay it, or fail loudly / return `Undetermined` |
| S2 | `Contains` cxx:987-997 | bbox pre-check runs *before* the documented `Contains_Loop` fallback; pre-CloseShape `fDX=0` | `Contains` is effectively disabled before `CloseShape` — the documented fallback path is unreachable |
| S3 | boundary policy: `kRayTolerance` filter + `isWantedCrossing` | a point exactly on a face is "inside" for `Contains`, but its t≈0 exit is invisible to `DistFromInside` → returns `Big` | navigator tunneling from a boundary; ROOT primitives snap \|t\|<tol to 0. Unify the convention |
| S4 | `nearestCrossing` cxx:212-247 | distance queries classify each hit independently — the mixed-cluster (tangential graze) cancellation that `Contains` has is missing | phantom entering crossing at a reflex/concave edge graze: `DistFromOutside` returns t, parity says outside → zero-step ping-pong. No concave fixture exists |
| S5 | cxx:86-94 vs cxx:114 | tangency tolerance differs between parity (strict `<0`) and distance classification (`\|dot\|<=kTolerance`) | Contains and distances disagree exactly on tangency-stressed configurations |
| S6 | `oddCrossingParity` (whole) | **The documented "order-dependent clustering" diagnosis is wrong**: on a given hit multiset the clusterer is deterministic (sort + consecutive-value chaining). The observed BVH-vs-loop `Contains` disagreements on non-manifold parts must come from *differing hit multisets* (float traversal / candidate sets) | the `CloseShape` error text and enum doc encode a wrong root cause; instrument (dump both sorted hit lists on divergence) before "fixing" clustering |
| S7 | `ComputeNormal` cxx:1144-1150 | nearest patch chosen by `distanceSqToPatch`, which is only a *lower bound* for wire-trimmed quadrics | wrong-face normals possible near trimmed-away regions — matters for optics/reflection |
| S8 | `validateClosure` BoundedSurface.h:3838-3890 | exact quantized (1e-7) vertex matching, no neighbor-bin handling, of chords produced through per-face frames consistent only to ~1e-6; B-spline edges resampled per independent pcurve | **the closure criterion structurally cannot say "closed" on real CAD** (independently confirms K9); counts are chord-inflated (~`kArcSamples` per unmatched rim). Needs topology-level matching + a quantitative gap metric |
| S9 | harness `O2SolidHarness.cxx:238-248` | two soundness holes: (i) candidate-returns-`Big` is probed at the *reference's own hit* → refSafety≈0 → auto-"within band"; (ii) tunneling to a parallel far wall also lands on the reference mesh → "within band"; `meshBand=1e-2` is a fixed guess | **a candidate that misses a whole wall is always "explained by mesh chording"** — the 367 unexplained points are a floor, not a ceiling. Add a `missedSurface` category and classify by \|dc−dr\|, not probe-point proximity |
| S10 | `O2SurfaceSolidIO.cxx` | `kJoinTolerance=1e-5` mixed-unit (same as K3, at IO layer); `readDoubles` resizes to an unvalidated `uint32` (corrupt file → 32 GB alloc attempt, exception escapes); B-spline clamping assumed, never validated (K1's IO twin) | loader fragility; cap counts against file size, use `angularTolerance(r)` for the φ coordinate |
| S11 | `DistFrom*` cxx:1039-1040 | O(N) `Safety` computed whenever `safe!=nullptr` regardless of `iact`; `validateClosure` re-runs the 128² capacity quadrature on every `CloseShape` even with `check=false`; per-query stack allocations in `visitPointCandidates`/`nearestCrossing` | avoidable hot-path cost |
| S12 | scale dependence | `kBVHBoxTolerance=1e-3` cm and relative `sameIntersection` merging (1e-7·\|t\|) untested away from the origin; no fixture at 1e3-1e4 cm offsets or with ≲1e-3 cm faces | unknown behavior at detector-scale coordinates — must be tested, ALICE geometries live at metres |

**TGeoShape support matrix (honest):** `Contains`/`DistFromOutside`/`DistFromInside` implemented
with BVH + `_Loop` oracles (bit-identical, the branch's strongest asset); `Safety` O(N) loop, no
BVH (measured worst kernel), underestimate-safe; `ComputeNormal` implemented with the S7 caveat;
`Capacity` exact for planar/untrimmed quadrics, 128² quadrature for wire trims (exactness flag
not exposed through `Capacity()`); `ComputeBBox` + visualization implemented;
`DistancetoPrimitive` stub; `GetPointsOnSegments` **inherits TGeoBBox** (wrong geometry for
ROOT's overlap checker sampling — worth an override or a loud comment); streaming **broken**
(S1). `iact` semantics ignored in distance queries.

**Test-coverage gaps, ranked:** (1) persistence round-trip (would catch S1 immediately);
(2) non-navigable-solid semantics pinned in a test, incl. the known BVH-vs-loop divergence as a
regression anchor; (3) adversarial distance fixtures — concave/reflex edges, rim grazes from
outside, on-boundary points vs ROOT primitives, near-tangential `dot` band; (4) randomized
metamorphic sweeps (parity invariant across several test directions;
`Contains(origin+d·(DistFromOutside+ε))`; chord consistency) with fixed seeds; (5) scale/offset
robustness (S12); (6) Safety/ComputeNormal on wire-trimmed quadrics; (7) IO fuzzing (counts,
NaNs, unclamped knots); (8) concurrency (thread_local buffers, global pruning flag). The
historical closed-B-spline bug class is now covered; **no fixture anywhere has
independently-approximated adjacent trims** — the dominant real-CAD failure class enters the test
suite only via Phase 0's Boolean ladder. No Python-side tests exist at all.

## 7. Findings — converter (`scripts/geometry/O2_CADtoTGeo.py`, exact-surface path ~1450 lines)

| # | Where | Defect | Consequence |
| --- | --- | --- | --- |
| C1 | `recognize_and_extract_face` L1679 | `inner_wall = (Orientation()==REVERSED)` for **recognized** faces — but the flag is relative to the *stored* NURBS normal, whose sign vs. the recognized frame is exporter-arbitrary; the reassuring comment is false | **containment-inverting** when wrong; the recognizer already computes the needed normal field — fold its sign in. Explains the known `ST0923290_013` winding warnings |
| C2 | quadric extractors L1254/1299/1337/1381 | `inner_wall` ignores **indirect (mirrored) gp_Ax3**: OCC's natural normal points inward there, so `reversed == inner_wall` inverts (planes handle this; quadrics don't) | containment-inverting on mirrored parts |
| C3 | `_quadric_trim_fills_uv_box` L1223-1229 | compares **signed** area to +box area → REVERSED faces likely never take the scalar fast path | correctness OK but doubled trim cost; or (if winding assumption differs) the test can't distinguish winding — needs a unit test either way |
| C4 | `_quadric_trim_wire` L1179-1188, `extract_planar_face` L1083-1090 | line-line junctions are chained tight, but **arc/B-spline edges keep their own endpoints** → vertex-tolerance gaps at line↔curve joints | parity flips through the gaps; cheap fix: snap all edge endpoints to shared vertex coordinates |
| C5 | `_recognized_quadric_wire_block` L1614-1626 | wire starting with a degenerate (apex/pole) edge anchors φ=0 arbitrarily; >π wedges across the apex pick the wrong unwrap branch | spurious polygon legs / inverted regions |
| C6 | `face_supported` L589-619 vs extractors | report gates on **3D** curve types, extractors on **pcurve** types — disagreements in both directions | eligibility numbers are neither an upper nor lower bound; keep extraction as source of truth and say so |
| C7 | sphere recognition L1750-1752 | polar axis fixed to global ẑ | recognized sphere caps cut by tilted planes rejected unnecessarily |
| C8 | `_quadric_phi_range` L1114 | silent clamp of >2π UV ranges | wrong trims on pathological seam-crossing faces |

Structural: the exact-surface subsystem (44% of a 3316-line script) is coupled to module globals
and has **zero automated tests**; it should become a package (`sidecar.py`, `extract.py`,
`recognize.py`, `report.py`) with in-process `BRepPrimAPI` fixtures (FORWARD and REVERSED and
mirrored variants of every family — that alone settles C1-C3). The recognizer is duplicated
between `O2_CADtoTGeo.py` and `analyze_surface_geometry.py` with a subtle policy difference
(first-at-machine-precision vs smallest-residual). Directory hygiene (five stale script variants,
editor backups, `a.out`) will bite the first outside contributor.

**The tube-tube walkthrough (what the converter emits today)** — face A and face B each go:
area test fails → `_quadric_trim_wire` → own pcurve (`BRep_Tool.CurveOnSurface`) → exact affine
pole map → sidecar B-spline in own (φ,h) chart; recognition finds no canonical form (correctly).
The only *shared* object, the edge's 3D curve, is dropped; the two emitted trims disagree by up
to 2× edge tolerance along their whole length, with independent parametrizations and re-derived
endpoints. The C++ closure check then honestly reports hundreds of boundary edges. No tolerance
tuning fixes this; only Section 4.2 (exact, for this case) or 4.3 (consistent) can.

---

## 8. The OCCT reference solid — build it, Python-first

Feasibility was verified **on this machine** (aarch64, alibuild): OCCT 7.9.3 with all needed
toolkits and headers (`BRepClass3d_SolidClassifier`, `IntCurvesFace_ShapeIntersector`,
`BRepExtrema_DistShapeShape`, `BRepGProp`, `BRepCheck_Analyzer`, `ShapeAnalysis`, STEP I/O);
pythonOCC 7.9.0 imports and answers containment/ray/volume queries correctly using
`$SW/Python/latest/bin/python3.10` with `PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages`
and `LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib` (system python3 is 3.12 — will
not import). The local `alidist/o2.sh` already lists OCCT as a dependency, and the untracked
`Utilities/TGeo_Utils/CMakeLists.txt` contains a working `find_package(OpenCASCADE)` pattern.

Recommended: **variant (b), the Python-side dataset oracle, first** (~2 days):

1. `O2_CADtoTGeo.py`/`makeTestPartDB.py`: also write `part_<VOL>_<LID>.brep` per exact part
   (`breptools.Write`, *after* the mm→cm scaling so units match the sidecar) and index it in the
   manifest.
2. Harness: `--dump-samples <json>` (the internal mt19937 samples are otherwise irreproducible).
3. New `occtOracle.py`: BREP + samples in → answers out (`BRepClass3d_SolidClassifier` for
   Contains with ON treated as "no verdict"; `IntCurvesFace_ShapeIntersector` min-positive-W for
   distances; `BRepExtrema_DistShapeShape` for the Safety upper bound; `brepgprop` volume;
   `BRepCheck_Analyzer` as per-part validity gate).
4. Harness: `--ref-answers <json>` — the existing `Offender`/validation machinery keeps working,
   only the reference source changes. Tolerance semantics: classify against a per-part band
   derived from the model's max face/edge tolerance (which the oracle can also report), **not**
   from mesh chord length.

The C++ `O2OCCTSolid` (TGeoBBox-derived, in-process, optional `find_package(OpenCASCADE)`) is a
later convenience (~2-3 days) once a live in-process reference is worth having; keep it out of
timing tables (BRepClass3d is ~ms/call, not thread-safe).

Note its value honestly: the oracle certifies *agreement with OCCT's tolerant classification*,
i.e. correctness up to model tolerance — exactly the right bar for real CAD, and a far stronger
bar than "inside the mesh band".

---

## 9. New success criteria

The current de-facto criterion ("`Contains` matches accelerated `O2Tessellated` outside the mesh
band, on parts that happen to convert") is retired. New gates, in order; each is a hard gate for
the work behind it:

- **G0 — kernel integrity.** `ctest -R BVHSurfaceSolid` green, plus new adversarial cases:
  periodic/unclamped knot input (K1), full-turn trims with overshooting poles (K2), small/large
  radius join metric (K3), degenerate-chord closed curves (K4), near-axis cone rays and torus
  fillet grazes (K6), seam-aligned rays, on-boundary points. Property held everywhere:
  `Contains(BVH) == Contains_Loop`.
- **G1 — the synthetic Boolean ladder (the new TDD backbone).** pythonOCC-generated STEP
  fixtures, in increasing order: box∪box (shared face), box−cyl (hole), cyl⊥cyl union ("plus"),
  cyl∩cyl, cyl−cyl (tube window — the Bagger reproducer minimized), obliquely cut cylinder
  (ellipse cap), cyl+cone, sphere−cyl, torus∪cyl. For EVERY fixture: `--exact-surfaces required`
  succeeds; `CloseShape`: **0 boundary edges, 0 non-manifold, Reliable**; vs OCCT oracle:
  **0 disagreements** for all samples farther than δ from the boundary, δ = the *declared* model
  tolerance (not mesh band); `Safety` ≤ true distance; capacity within stated tolerance;
  BVH == Loop bit-identical.
- **G2 — as1-oc-214**: 5/5 exact, all Reliable, oracle-clean per G1.
- **G3 — Bagger**: **13/13 exact, all Reliable, oracle-clean.** This is the headline target the
  user set for this phase; with Section 4.2 it is achievable *exactly* (every Bagger edge is
  analytic∩analytic).
- **G4 — oTOF + ALICE3**: coverage reported honestly; **an unreliable exact part must never ship
  silently** — `auto` mode falls back to tessellated when `NavigationReliability != Reliable`
  (policy change: today it ships and only warns).
- **G5 — performance** (after correctness): `Safety` gets the priority-queue BVH; per-query
  medians within ~2x of `O2Tessellated` on the DB; numbers recorded in
  `SolidNavigationHarness.md` as now.

A part "converts" only if it passes its gate; eligibility numbers without the Reliable+oracle
qualifier are no longer quoted.

---

## 10. Recommended plan (ordered)

**Phase 0 — oracle + fixtures (≈3-4 days).** Build the OCCT oracle (Section 8) and the G1
fixture generator (`scripts/geometry/make_boolean_fixtures.py`); fix the harness classification
holes (S9: `missedSurface` category, classify by |dc−dr|, band derived from actual deflection).
Run all of it on the current code and record the honest baseline per fixture. *No solid/converter
changes in this phase.* Also instrument, not fix, the non-manifold divergence: dump both sorted
hit lists on any BVH-vs-loop disagreement (S6) — evidence before theory, per the project's own
2026-07-26 lesson.

**Phase 1 — kernel robustness (≈1-2 weeks).**
1. Ambiguity band + re-shoot parity (Section 4.4) — expected to clear a large share of residual
   `unexplained` points and make all later work loud instead of silent.
2. Fix the safety-critical solid bugs: persistence (S1 — stream the sidecar records or fail
   loudly), pre-CloseShape `Contains` (S2), boundary-policy unification (S3/S5), cluster-aware
   distance queries sharing one helper with parity (S4).
3. Fix K1-K8 (each with its G0 regression test); make face-drop (K7) fatal in `required` mode
   and downgrade-to-fallback in `auto`.
4. Unify the tolerance policy: per-domain metrics (angular ↔ length via radius) applied in both
   kernel and IO (K3/S10), one documented ε family; sidecar format v2 carries the model
   tolerance (and, for Phase 3, edge IDs). Replace the closure criterion with topology-level
   matching + a quantitative gap metric (K9/S8) so "closed" becomes achievable and meaningful.

**Phase 2 — exact analytic-analytic trims (≈1-2 weeks).** Implement neighbor-surface trim edges
(Section 4.2) kernel-side and converter-side, driven strictly by the G1 ladder: first box−cyl,
then cyl⊥cyl, then the full ladder, then G2/G3. This is the phase that makes Bagger exact and
watertight, and it is where the TDD discipline the user asked for lives.

**Phase 3 — shared-edge consistency for the remainder (≈1 week).** `ExactTrimTopology.md` item
1's recipe for free-form edges (edge map, shared samples, per-face closed-form projection,
pinned vertices), now justified by evidence and needed only where E2 does not apply.

**Phase 4 — resume coverage and performance.** Converter fixes C1-C8; ellipse planar trims
(subsumed by E2 for plane∩cyl but keep the 3-line fallback); `Safety` BVH; subdivision-BVH;
then, and only then, re-open the B-spline *surface* milestone with the honest
exactness-vs-tessellation trade-off already recorded in `BVHSurfaceSolid.md`. Performance
direction decided 2026-07-31: since the geometry is static once converted, prefer **AOT
geometry-specific code generation** (converter-emitted specialized C++, VecGeom
specialized-navigator spirit) over runtime JIT, with cling available if in-process generation is
ever needed, and clad-style automatic differentiation if Newton-based surface intersection is
ever implemented — details in the `BVHSurfaceSolid.md` optimization milestone.

**Explicitly de-prioritized**: further coverage sweeps and ALICE3 milestones until G1-G3 are
green — measuring coverage of an unreliable representation is motion, not progress.

---

## 11. Phase 0 executed (2026-07-31) — the gate exists, and it separates the ladder exactly

Phase 0 of Section 10 is implemented and run. Everything below is measured, not projected.

**What was built.**
- `scripts/geometry/occtOracle.py` — the OpenCascade oracle. Answers `Contains`
  (`BRepClass3d_SolidClassifier`, with an explicit ON = *no verdict* state), `DistFrom{Outside,
  Inside}` (`IntCurvesFace_ShapeIntersector`; deliberately one computation for both, since entry
  versus exit is a property of the origin, not of the intersector), the exact boundary distance
  (`BRepExtrema`, the upper bound `Safety` must not exceed), capacity (`BRepGProp`), validity
  (`BRepCheck_Analyzer`) and the model's own max sub-shape tolerance — which is what makes the
  comparison band principled instead of guessed. It carries a `--self-test` mode that checks
  every kernel against closed-form geometry; **it found two defects in itself before it judged
  anything** (distances were being measured against the solid rather than its shell, so every
  interior reference distance was 0). Runtime is ~0.3 s for 1250 points + 750 rays.
- `scripts/geometry/make_boolean_fixtures.py` — the G1 ladder: 10 synthetic CAD fixtures in mm,
  all `BRepCheck`-valid, volumes cross-checked against closed form to 1e-11…1e-16 (the two
  without closed form verified by independent quadrature).
- `scripts/geometry/runOracleGate.py` — one command for the whole round trip (convert → dump
  samples → oracle → score) with a pass/fail verdict per part.
- `--dump-brep` in `O2_CADtoTGeo.py` (verified purely additive; units confirmed by volume
  agreement to 1e-16…7e-14 relative, where a missing mm→cm scale would show as 1000×) and a
  `brep` field in `makeTestPartDB.py`'s manifest.
- Harness: `--dump-samples` / `--ref-answers`, the S9 classification fix (a distinct
  `missedSurface` category that can never be explained away, plus incidence-scaled band
  `meshBand/|cos θ|`), and `O2BVHSurfaceSolid::DescribeContainsCrossings` — the S6 instrumentation
  that prints both hit lists when BVH and loop parity disagree, so that divergence gets diagnosed
  rather than theorized about.

**The measured ladder baseline** (2000 points / 2000 rays, seed 1; band = the model's own 1e-7
tolerance). A part passes only if it is navigable *and* clean on all four kernels *and* its
capacity matches:

| fixture | junction curve | navigable | contains | distout | distin | safety | capacity relDev | gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `box` | — | yes | 100% | 100% | 100% | 100% | 1.5e-16 | **PASS** |
| `box_union_box` | shared face | yes | 100% | 100% | 100% | 100% | 0 | **PASS** |
| `box_minus_cyl` | plane ∩ cyl (circle) | yes | 100% | 100% | 100% | 100% | 2.4e-14 | **PASS** |
| `cyl_plus_cone` | coaxial (circle) | yes | 100% | 100% | 100% | 100% | 7.7e-14 | **PASS** |
| `sphere_minus_cyl` | coaxial (circle) | yes | 100% | 100% | 100% | 100% | -8.3e-14 | **PASS** |
| `torus_union_cyl` | coaxial (circle) | yes | 100% | 100% | 100% | 100% | 9.3e-14 | **PASS** |
| `cyl_cross_cyl` | **cyl ∩ cyl** | yes | 100% | 100% | 100% | 100% | **2.3e-4** | FAIL |
| `cyl_inter_cyl` | **cyl ∩ cyl** | yes | 100% | 99.90% (**2 missed**) | 100% | 100% | **3.5e-4** | FAIL |
| `tube_window` | **cyl ∩ cyl** | **no (1418 boundary edges)** | 99.96% (2) | 99.55% (**8 missed**) | 100% | 100% | **-4.4e-3** | FAIL |
| `oblique_cut_cyl` | plane ∩ cyl (**ellipse**) | — does not convert: `planar boundary edge is a ellipse curve` — | | | | | | FAIL |

**The predicted failure set is exactly the observed failure set.** Every fixture whose junction
curves are circles (including the sphere and torus cases — they are coaxial, so their
intersections *are* circles) is exact to machine precision on all four kernels. Every fixture
with a genuine cylinder-cylinder intersection curve fails, and the severity tracks how much of
the boundary that curve carries: `cyl_cross_cyl` (two small seams) still navigates correctly and
only leaks into the numerically-integrated capacity; `cyl_inter_cyl` (whose entire boundary is
intersection curve) starts missing ray crossings; `tube_window` is not a closed manifold at all.
This is Section 4's argument, reproduced on 3-to-8-face fixtures that run in seconds — the
minimal reproducers Phase 2 needs, replacing an 8-face part inside a 13-volume CAD model.

`oblique_cut_cyl` reproduces the Bagger `Bucket` ellipse blocker in a 3-face fixture.

**Real CAD, same gate** (`Bagger.step`, 12 converted parts): **4/12 pass**. The failures are the
same three families: not-navigable parts (`BucketLink1`, `StickCylinderInner/Outer`, …), missed
crossings (`BucketLink2`: 24 missed in `distout`), and capacity drift (0.26%…6.2%). Notably the
one part checked in detail against both references — `BasePin` — shows **4 "within band"
mismatches against the mesh but 100% agreement with the oracle** on all four kernels, capacity
matching to 9e-16 relative: there the mesh was wrong, not the solid. That is the reason the
oracle exists.

**Gate status: G0 green** (`ctest -R BVHSurfaceSolid`, 36 cases, after every change in this
phase). **G1 6/10, G3 4/12** — the honest starting line for Phase 1 and Phase 2.

## 12. Phase 1 executed (2026-07-31) — items 1-3 of four

Phase 1's first three items are implemented and measured; item 4 (tolerance policy and closure
criterion) is not started. Everything below is measured, not projected. The gate was run before
and after every item, so each effect is attributable.

Commits, in order: `63d1d08119` (S1), `29e8322f79` (re-shoot), `f809b38dd2` (S2-S5),
`346d4675d0` (K1/K2/K8). `ctest -R BVHSurfaceSolid` green throughout, **48 cases** (from 37).

**Headline: `contains` is now clean.** Disagreements with the OpenCascade oracle outside
tolerance, per model:

| | fixtures | Bagger |
| --- | --- | --- |
| Phase 0 baseline | 2 | 56 |
| after re-shoot | **0** | **2** |

Gate *totals* are unchanged at **G1 6/10, G3 4/12** — every part that still fails does so on
navigability, on the distance queries, or on capacity, none of which these items touch. That is
the honest reading: Phase 1 removed the containment error class, and what is left is Phase 2's.

### What each item did

**1. Re-shoot parity (Section 4.4), `29e8322f79`.** Diagnosed before fixing. The new hook
`ContainsAlongDirection()` evaluates parity along an explicit direction, which turns the solid
into its own oracle — on a closed 2-manifold parity is a topological invariant, so *any*
direction dependence is a defect, with no reference shape involved. Measured over the whole
Phase 0 corpus (~11k points per part):

- every part the closure check calls `Reliable` has **zero** disagreements between shooting
  directions; every part that disagrees (0.12%-3.7% of its points) is one the closure check
  already rejects. The correlation is exact on this corpus.
- of the 55 points where the fixed direction disagrees with the oracle, **not one is wrong in
  every direction**. This corrects the recorded picture: a gap's shadow belongs to the
  *direction*, not to the point, so it is escapable.
- majority vote is right on 50/55 with three spread directions, 53/55 with five, 55/55 with
  thirteen. A hand-picked triple scored only 42/55 — two of its directions turned out correlated
  and agreed on the wrong answer, which is why the implementation uses a golden-angle spiral.

So `Reliable` solids keep the single shot and unreliable ones vote over five directions, stopping
at three agreements. Cost: **1.0-1.5x** on Reliable parts (tens of ns absolute), **2.9-4.7x** on
the already-unreliable ones. The closure check's own error is in the safe direction for this
policy — it under-reports `Reliable` rather than over-reporting it (S8/K9).

**2. S1-S5, `63d1d08119` + `f809b38dd2`.**

- **S1** (done first, being actively harmful): the solid streamed nothing but its `TGeoBBox` base,
  so a read-back solid was empty and — because an empty `ClosureReport` defaults to closed and
  consistently oriented — reported itself `Reliable`. It now persists the `Add*Surface` call
  sequence as `BVHSurfaceRecord`s and replays them on read, recomputing the closure diagnostics
  rather than trusting them; a record that fails to rebuild discards the whole solid. `CloseShape`
  refuses an empty surface set unconditionally and leaves the shape `Undetermined`.
- **S2**: the bbox pre-check ran ahead of the "no BVH yet, use the loop" fallback, and before
  `CloseShape` the box is all zeros — `Contains` was disabled on any unclosed solid. Order fixed.
- **S3**: a point on a face was inside for `Contains` but infinitely far from the wall for
  `DistFromInside`. Distance queries now admit hits from just behind the origin and clamp to zero
  (ROOT's own convention); parity keeps its positive lower bound so it never re-counts the surface
  the point sits on.
- **S5**: parity used a bare sign test and the distance queries a tolerance band. Both now go
  through one `crossingSense()` with `Entering`/`Exiting`/`Tangential`.
- **S4**: a ray that only *touches* the boundary has not crossed it. `Contains` knew this; the
  distance queries classified each hit alone and reported the touch as a crossing, handing the
  navigator a step to a point where it then finds itself outside. Both now share one
  `forEachCrossingCluster()` walk. Because a hit's meaning depends on its neighbours, the BVH
  query collects hits and classifies afterwards; tmax tightening still shrinks to the nearest
  candidate with a margin wide enough to keep its cluster partners, and if the candidate turns out
  to be a graze the query is redone unpruned. Cost **1.00x-1.10x**.

The **concave fixture the review found missing** now exists: an L-shaped prism with a genuine
reflex edge, the only configuration where a ray can touch the boundary from inside and stay
inside. It carries the full sweep battery.

**3. K1, K2, K8, `346d4675d0`.** K1 (endpoints of an unclamped B-spline are now evaluated instead
of read off the poles), K2 (the full-turn rejection re-measures on the curve before refusing, not
on the conservative pole hull), K8 (cluster membership tested against the cluster's first member,
so transitive chaining can no longer merge a thin far-away feature away). All three are latent on
this corpus — gate columns bit-identical — because the converter's canonical recognition converts
the offending inputs away before the kernel sees them. Each has a regression test; K2's was
verified to fail against the pre-fix code.

### Corrections to this document

- **K7 does not hold.** A face that fails to build is *not* silently omitted on any production
  path: `Add*Surface` returns false, `LoadSurfaceSolid` turns that into a whole-file rejection,
  and the generated macro turns that into an exception. The true, weaker statement is that the
  return value is the only signal. Pinned by a test rather than "fixed".
- **The shadow of a gap is direction-bound.** Earlier notes describe a gap as flipping `Contains`
  over its whole shadow, which is right, but read as though the point were lost. It is not: 55/55
  of the affected points are recoverable by aiming elsewhere.

### Still open

*(Updated later the same day, after Phase 1 item 4 steps 1-4 and the `BucketLink2` diagnosis of
Section 13. The list below is the state as of then.)*

- **Phase 1 item 4 is four fifths done.** Steps 1-4 are committed (`1bc5c4fbc9`, `9f45887ef7`,
  `cacd64e4a5`, `f612f895a9`): every surface reports a first fundamental form; the wire join checks
  in kernel *and* loader judge a gap as a length in cm against one constant, closing **K3, K12 and
  S10**; sidecar version 2 carries the model's declared tolerance; and the on-boundary band is
  sized from the representation's real accuracy with winding and distance sharing one polyline,
  closing **K5**. All four gate bit-identical (fixtures 6/9, Bagger 4/12) — for steps 1 and 3 that
  is the required outcome, and for steps 2 and 4 it is the expected one, since the measured worst
  join residual on this corpus is 4.06e-11 cm and a 1e-5-wide boundary shell is not something
  bbox-spread samples land in. All four are pinned by tests instead, and step 2's were verified to
  fail against the pre-fix rule.
- **Step 5 — rim-based closure and `maxGap` (K9/S8) — is not started**, and it is the one with
  design content. [`TolerancePolicy.md`](TolerancePolicy.md) §8 splits it in two so that building
  the gap measurement does not, by itself, invalidate the measurement licensing `Contains`'s
  single-shot fast path. Read that before starting.
- **K4** (degenerate-chord recursion probes only the parametric midpoint) and **K6** (cancellation
  in the naive quadratic formula; absolute, scale-dependent tolerances in the cone-degeneracy and
  torus-quartic branches) are untouched. `BucketLink2` is no longer evidence for either.
- **`BucketLink2` is resolved and is no longer a lead** — see Section 13. The two items that
  replace it are gate ray-category soundness and Green's-theorem capacity for wire-trimmed
  quadrics.
- The distance and capacity columns are otherwise untouched by Phase 1 by construction.

## 13. `BucketLink2` diagnosed (2026-07-31) — the kernel was never wrong here

The "best remaining lead" of Section 12 is resolved, and it resolves *away* from the kernel. Three
throwaway probes; no fix was needed in navigation. Two new findings replace it, both precisely
scoped. Numbers below are from the Phase-1 baseline tree (nothing was changed to obtain them).

`BucketLink2` failed the gate on three counts: 24 *missed* `distout` crossings, 48 `distin`
disagreements, 6.3% capacity drift. It has **12 planes and 11 cylinders, no cone and no torus**,
which already rules out two of K6's three branches by construction.

**H1 — the 24 + 48 distance disagreements are an artifact of the gate, not of the solid.**
Probe 1 asked OpenCascade to classify each of the 20 recorded offender ray origins and to list
every crossing with its face index. Every single `distout` origin is `IN` and every single `distin`
origin is `OUT` — that is, the *opposite* of the category the sample generator assigned. Probe 2
then asked the exact solid and the tessellated solid the same question:

- the exact solid's classification matches OpenCascade on **20 of 20**;
- the tessellated solid is opposite on **20 of 20**;
- and where the category is right, the exact solid's distance matches the oracle to every printed
  digit — e.g. `distin#0` exact `DistFromOutside` = `0.143949` against oracle `0.143949`,
  `distout#0` exact `DistFromInside` = `0.870144` against oracle `0.870144`, and so for all 20.

The cause is that `sampleSets()` categorises rays with the **mesh**, and this mesh is wrong. A
`Contains` scan along x at three of the offenders' (y, z) puts the mesh's left plate at
`x ∈ [-2.10, -1.90]` where the BREP's is `[-1.92, -1.42]`; `O2Tessellated::Check` reports several
hundred facets with only two neighbours, i.e. the mesh **is not watertight**, so its own parity
classification is undefined in the shadow of its holes. The gate then asks the candidate for
`DistFromOutside` from a point that is inside, compares it against an oracle number that is an
*exit* distance, and books the difference as a missed surface.

This is S9's failure mode one layer earlier: S9 says a wrong candidate answer can be excused by
the mesh band; this says a *right* candidate answer can be condemned by the mesh classification.
The fix belongs in the gate: the oracle should report its own containment for each ray origin, and
the harness should re-label — or at minimum decline to score — a ray whose category the oracle
contradicts. That rule is sound because it never consults the candidate.

**H2 — the capacity column is entirely a quadrature artifact.** `Contains` is the query the oracle
already agrees with 5000/5000 on this part, so probe 3 used it as the witness: a 4M-point Monte
Carlo of the exact solid gives **17.061 ± 0.052 cm³** against OpenCascade's **17.079** (0.35σ),
while `Capacity()` returns **16.004**. The solid is not missing material; the integrator is wrong.

`integrateOverCurveTrim` (`BoundedSurface.h:1690`) is a midpoint rule over the trim's parametric
bounding box with a hard in/out test per cell, at a fixed `samplesPerAxis = 128`. Integrating a
characteristic function that way converges at **O(1/N)**, not O(1/N²), because every boundary cell
is booked whole or not at all. Making `samplesPerAxis` temporarily tunable shows exactly that,
oscillating as the staircase re-phases:

| samples/axis | 128 | 256 | 512 | 1024 | 2048 | OCCT |
| --- | --- | --- | --- | --- | --- | --- |
| `Capacity()` | 16.004 | 17.710 | 16.927 | 17.244 | 17.032 | 17.079 |

And it accounts for the whole column, across the model:

| part | 128 | 1024 | OCCT |
| --- | --- | --- | --- |
| `BasePin` | 31.415927 | 31.415927 | 31.415927 |
| `Base` | 241.280940 | 241.280940 | 241.280940 |
| `Boom` | 1204.575249 | 1204.575249 | 1204.575249 |
| `Stick` | 333.942179 | 333.942179 | 333.942179 |
| `BoomCylinderInner` | 23.029627 | 22.927630 | 22.989411 |
| `BoomCylinderOuter` | 48.150825 | 48.243639 | 48.197441 |
| `BucketCylinderInner` | 9.050783 | 9.042326 | 9.045578 |
| `BucketCylinderOuter` | 16.338176 | 15.896115 | 16.011732 |
| `BucketLink1` | 13.065863 | 12.292618 | 12.304424 |
| `BucketLink2` | 16.003995 | 17.243583 | 17.079260 |
| `StickCylinderInner` | 23.061991 | 22.926726 | 22.989411 |
| `StickCylinderOuter` | 48.071733 | 48.243921 | 48.197441 |

The four parts that agree with OpenCascade to all printed digits are exactly the four with **no**
wire-trimmed quadric — their contribution is the closed form and is identical at both resolutions.
Every part that deviates has wire-trimmed quadrics and moves toward the oracle under refinement.
**No geometry defect is involved anywhere in the capacity column.**

The real fix is not a bigger N — at O(1/N) the gate's 1e-6 relative band is unreachable — but
Green's theorem, which is already the mechanism the planar wires use. For a cylinder the integrand
`(X·n)·r/3 = (r/3)(C·u cos φ + C·v sin φ + r)` does not depend on `h`, so an antiderivative in φ
exists in closed form and the area integral collapses to a contour integral around the trim wire,
evaluated by the same Gauss-Legendre the planar `signedAreaContribution` already runs. That makes
the wire-trimmed quadric capacity exact for line and arc trims and spectrally accurate for
B-splines — and it would let `capacityIsExact()` finally mean something. The same reduction exists
for the cone, sphere and torus.

Until then, note the honest statement: `Capacity()` on a wire-trimmed quadric is a ~1e-2 relative
number that **does not say so** — `capacityIsExact()` exists but is not exposed through
`Capacity()`, which Section 6's support matrix already flagged.

**Consequences for the record.** K6 and K4 remain untouched, but `BucketLink2` is no longer
evidence for either of them, and the "distinct defect that no current theory covers" is retired.
The two items that replace it are above: gate ray-category soundness, and Green's-theorem capacity
for wire-trimmed quadrics.

## 14. Phase 1 item 4 completed (2026-08-01) — the gap is a length now, and it is not what we thought

Step 5 of item 4, in two commits (5a measure, 5b decide), as split by
[`TolerancePolicy.md`](TolerancePolicy.md) §8. Both gates bit-identical on both commits (fixtures
6/9, Bagger 4/12, `contains` 0 and 2); `ctest -R BVHSurfaceSolid` green at **57 cases**, from 53.
Every measurement is in `TolerancePolicy.md` §9 and §10; the headlines and the corrections are
here. **K9 and S8 are closed, and Phase 1 is finished.**

**The closure check is now a curve comparison and a number in centimetres.** Each face emits its
boundary as rims — one ordered polyline per trim loop — and each rim chord is matched, at its
midpoint, against the chords of every other face. `ClosureReport` gains `maxGap`,
`unmatchedRimLength`, `totalRimLength`, the matching epsilon (the model's own declared tolerance
from sidecar v2, else a documented constant) and per-rim counts;
`NavigationReliability`/`IsClosed()` read them, `Print()`, the harness line and `--json` report
them, and `CloseShape`'s errors quote the gap in cm and the open fraction of the boundary instead
of a chord count. The per-chord counters remain as a diagnostic and decide nothing.

**Correction to the recorded expectation, by four orders of magnitude.** `TolerancePolicy.md` §3.3
predicted, from the 2026-07-26 shared-edge pcurve measurement, that the Bagger cyl-cyl parts would
stay open at ~1.3e-5 cm. Measured: they are open by **0.25 to 0.75 cm**, over 4-15% of their rim
length. That is a face missing or trimmed to the wrong curve, not two pcurves disagreeing in the
fifth decimal, and it changes what Phase 2 is for on these parts — they are not nearly closed and
then spoiled by tolerance. `tube_window`, the motivating case for "the counts are chord-inflated",
went from **1418 boundary edges** to seven rims of which three are open, 9.94 cm of 53 cm.

**Correction to §1.3 of that document, in the other direction.** It predicted that vertex matching
fails on real CAD because the two faces of a shared edge sample it at different phases. The
structural defect is real — a unit test builds a box whose last face is sampled at twice the chord
count and the old check calls that closed box open — but **no part of either corpus exercises it**:
every part the old check called closed is rim-matched and vice versa, part for part. So 5b's
verdict switch is a no-op here, and the gate says so. Take that as "this corpus does not test the
thing", not as "the thing was not a problem".

**The §4.2 sweep, and a new lead.** 21 parts, 11000 points each, 13 golden-spiral directions,
before and after the switch: the same 13 `Reliable` parts and **one** disagreement in 143000
points. The parts that do disagree between directions are exactly the ones the closure check
rejects, at 0.55%-7.0%. So the licence for `Contains`'s single-shot fast path survives the change
intact. The one offender — `cyl_cross_cyl`, 1 of 13 directions, `Safety` 0.0992 cm, 12 directions
including the fixed one agreeing — is **not** a gap shadow (a gap costs a point most of its
directions and puts it behind the gap). It is one unlucky ray near the curve where two cylinders
cross, and it is now the cheapest known reproducer for **K6**.

**One thing deliberately not shipped, with the measurement that says so.** §2.4's ambiguity band
was deferred into 5b. The free form of it — re-shoot whenever a parity cluster held more than one
coincident hit, i.e. whenever the cancellation rule rather than the geometry decided — was built
and measured: **0.2-1.3% of `Contains` on seven fixtures, +15.8% on `box_minus_cyl`, and not one
of the 143000 sweep points moved**, including the offender it was aimed at. It fires where the
cancellation rule is already right and stays silent where the residual defect is. Dropped, with the
numbers recorded in `containsByParity` so nobody rebuilds it. The instrument §2.4 actually asks for
— a band around the trim curve in the parametric domain, sized through the surface metric — is
unbuilt, and §10.1's headroom of one point in 143000 sets its priority behind diagnosing that
point.

**Three implementation findings worth carrying.** The plan's "cheap way in" — chain runs of
consecutively emitted edges — does not work, because the parametric-rectangle quadrics interleave
their two rims; chaining by matching endpoints does, and still needs no change in any of the seven
surface classes. A full-turn patch whose `fullSweep()` does not fire emits its seam twice in
opposite directions, which cancels in the half-edge check and so had never been noticed, but chains
into a two-point rim straddling the patch and produced a spurious boundary rim on five of nine
fixtures until reversed duplicates were cancelled first. And the sampling noise floor must be
estimated from the *turn angle*: the obvious estimator, a vertex's deviation from the chord joining
its neighbours, calls a box corner 2.4 cm of noise.

## 15. Environment notes (this machine)

- Build tree `/home/swenzel/alisw/sw/BUILD/O2-latest/O2` builds this checkout (source symlink);
  the branch's test target `o2-test-detectorsbase-BVHSurfaceSolid` appears after a CMake re-run.
  (Not built in this session per the user's choice.)
- pythonOCC/OCCT working invocation: see Section 8 (alibuild `python3.10`, not system 3.12).
- `alienv` is at `/home/swenzel/alisw/alibuild/alienv`; `ALIBUILD_WORK_DIR=/home/swenzel/alisw/sw`.

## 16. Phase 1 follow-up item 1 (2026-08-01) — the direction-dependent point is K5, not K6

NEXT.md's item 1 asked to diagnose the single disagreement the section 4.2 sweep left, on the
theory that it was the cheapest known reproducer for **K6** (cancellation in the naive quadratic
formula). It is not, and K6 loses its only reproducer. Full measurements in
[`TolerancePolicy.md`](TolerancePolicy.md) section 11; the headlines are here.

**The roots are fine; the trim is not.** On `cyl_cross_cyl` the offending ray gains a third
crossing on the x-cylinder at a point that lies 3.7e-6 cm *inside* the z-cylinder — interior to
the fused solid, where that patch should have been trimmed away. Both cylinders are hit exactly
where they should be, to the last printed digit. The two hits are 5.4e-5 apart in relative terms,
so no clustering rule could ever have merged them, which is also why section 10.2's coincident-hit
trigger was silent on the one point it was aimed at.

**The overhang is systematic and one-sided, and its floor names its cause.** Walking the analytic
seam of the two cylinders and bisecting where each patch actually ends: **all 1440** sampled
positions overhang the true seam, **none** undercuts, worst 1.95e-5 cm, floor **1.00e-5 cm** — and
1e-5 is `kBSplineFlatness`, i.e. `CurveWire::boundaryBand` for a B-spline trim. `curveTrimContains`
resolves the `Boundary` state as *inside the trim*. That tie-break is one-sided by construction, so
every wire-trimmed patch keeps a sliver of the band's width past its own trim curve, and on a
Boolean seam that sliver sits in the solid's interior.

So this is **K5's other half**. Step 4 widened the band correctly — a 1e-9 band around a 1e-5
polyline was measuring noise — but never said what a hit *inside* the band means to a ray, and
"accept it silently" converts an undecidable point into a wrong answer that leaves no trace.

**The fix is section 2.4's ambiguity band, and it was far cheaper than section 2.4 feared.** No new
geometry and no per-class trim mathematics: `curveTrimContains` already computed the flag and every
caller was passing `nullptr` to it. `RayHit` gains `onTrimBoundary`, the five `appendIntersections`
that can produce a wire-trimmed hit set it, `parityAlong` reports it, and `containsByParity`
re-shoots through the existing vote when a `Reliable` solid's single shot rests on one.

**Measured on `cyl_cross_cyl` against the analytic truth**: uniform sampling, 5e6 points — the
re-shoot fires on 48 shots (9.6e-6), changes 24 answers, **corrects 24, breaks 0**, and `Contains`
goes from 24 wrong to **0**. Aimed at the seam, 2e6 points — 453008 wrong before, **0 after, 0
broken**. Cost: one shot in 104000, at most five extra parity shots each, ~5e-5 amortised.

**Both gates bit-identical** (fixtures 6/9, Bagger 4/12, `contains` 0 and 2) and `ctest` green at
**58 cases**. The gate cannot see this: its `contains` column already agreed 5000/5000 on every
`Reliable` part and the defect is at the 1e-5 rate. Two traps recorded for whoever reads those
numbers next:

- **The section 4.2 sweep is blind to this change by construction.** Its metric is
  `ContainsAlongDirection`, which documents that it bypasses the re-shoot policy. The previous
  session's 1-in-143000 and this session's 0-in-143000 are two samples of the same rate, not a
  before and an after.
- **The gate's timing column cannot resolve anything at this scale.** Its `contains` numbers moved
  by −29% to +46% between two runs of the same code path on this machine.

**Scope, stated plainly.** Parity is now self-checking about its own tie-breaks. The distance
queries are not: `DistFromOutside`, `DistFromInside` and `ComputeNormal` still consume a flagged
hit silently, and `nearestCrossing` has no equivalent policy — which is S4 seen from the distance
side. And the sliver itself remains; labelling it is not removing it. Removing it needs the seam
known better than either face's chart knows it, which is Phase 2 — for which `onTrimBoundary`
becomes a ready-made acceptance test: on a model with exact adjacency it should never fire.
