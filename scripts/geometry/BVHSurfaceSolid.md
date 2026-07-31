The O2_CADtoTGeo tool translates a CAD geometry to a TGeo geometry, where all parts
are represented as a TGeoTessellated solid. This constitutes an approximation in most cases.

The goal of this project is to make the translation exact. The idea is to keep the surface
modeling of surface patches exact whenever this is possible (planar, second-order, etc.) and
to map the CAD parts to a new BVHSurfaceSolid. This is most likely like a BREP solid, but
the goal is to make a new, efficient implementation.

The design goals are:

(a) offer classes and modelling of bounded surfaces
(b) make the final BVHSurfaceSolid on top of bounded surfaces using the existing BVH acceleration structures, already used in O2Tessellated solid
(c) the final BVHSurfaceSolid needs to implement TGeoVolume/TGeoShape navigation functions (just like O2TessellatedSolid)
(d) the O2_CADtoTGeo translation should be modified to map to the new solid whenever possible (might not be possible in some complicated BezierCurve cases, etc.); in this case, keep the TessellatedSolid fallback
(e) unit tests should be constructed testing all elements of BVHSurfaceSolid

Next to the source code of O2TessellatedSolid, inspect the source code of VecGeom at
https://gitlab.cern.ch/VecGeom/VecGeom, which has a similar concept of SurfaceBounded
geometries, although in a slightly different setting.

Geant4 previously had a BREPSolid but abandoned it due to complicated maintenance. Now is
the time to make a good and fast implementation again.

# Implementation plan

This section is meant to be a durable handoff checklist. Each agent session should pick a
single unchecked item or a small cluster of adjacent unchecked items, keep the change focused,
run the listed validation, and update the checkbox plus a short note when done.

## Current anchors

- Existing accelerated tessellated shape: `Detectors/Base/include/DetectorsBase/O2Tessellated.h` and `Detectors/Base/src/O2Tessellated.cxx`.
- Existing BVH glue: `Detectors/Base/src/bvh2_third_party.h` and `Detectors/Base/src/bvh2_extra_kernels.h`.
- Build and ROOT dictionary integration point: `Detectors/Base/CMakeLists.txt`.
- CAD conversion entry point: `scripts/geometry/O2_CADtoTGeo.py`.
- Existing STEP examples and generated outputs: `scripts/geometry/STEP_examples/`.

## Proposed target shape

Use a new `TGeoBBox`-derived shape in `DetectorsBase`, tentatively named `O2BVHSurfaceSolid`.
The public shape should follow the `O2Tessellated` pattern: keep heavy BVH implementation out of
the public header, build all acceleration structures in `CloseShape`, and override the same TGeo
navigation methods:

- `Contains`
- `DistFromOutside`
- `DistFromInside`
- `Safety`
- `ComputeNormal`
- `Capacity`
- `ComputeBBox`
- display mesh methods (`GetBuffer3D`, `MakeBuffer3D`, `SetPoints`, `SetSegsAndPols`) using a visualization triangulation only

The key difference from `O2Tessellated` is that BVH leaf primitives are bounded analytic surface
patches instead of triangles. A bounded patch should be modeled as an infinite analytic surface plus
one or more trim wires, following the OpenCascade face/wire/edge idea: the surface provides the
analytic kernel, wires define the finite domain, and edges/curves define each wire. Triangles remain
only for visualization and fallback paths.

## Adopted design principles (2026-07-31)

Decided after the comprehensive review in [`CodeReview_Fable.md`](CodeReview_Fable.md) (read its
Sections 3-4 for the full derivation). These principles govern all subsequent milestones:

1. **Trims are adjacency, not curves (regime E2).** For an edge where both incident faces have
   analytic carriers, the trim boundary is represented as a reference to the *neighbor surface*
   plus one discrete side/branch bit and a parameter interval — never as a fitted 2D curve. The
   face-adjacency graph (the shared `TopoDS_Edge` identity that STEP already carries) is the
   authoritative object; membership tests evaluate closed-form sign/root computations against the
   neighbor's implicit equation. Both faces then place their common boundary at the identical
   point set, so watertightness holds *by construction* and the closure check becomes a
   structural graph check (every edge has exactly two faces). This is the same move that made the
   tessellated solid watertight: agreement through shared topology instead of tolerance-matched
   independent geometry. Transcendental pcurves (tube-tube intersections) never need to be
   stored or flattened. Fitted/sampled trim curves remain only as the fallback where a free-form
   surface is involved (regime E3: one shared sample set per edge, pinned vertices).
2. **Parity must be self-checking.** Point classification reports trim-boundary-ambiguous hits
   and re-shoots along a different direction instead of silently guessing; degenerate
   configurations cost a bounded retry, never a silent wrong answer.
3. **Simplify to the exactly-equivalent simplest form, at every level.** Canonical recognition
   (surfaces, curves) is one instance of a general principle: STEP delivers numeric expression
   soup (poles/weights/knots, transforms, swept/revolved constructions), and wherever the stored
   expression is *exactly* (machine-precision residual) a simpler analytic object, the simpler
   form is adopted — e.g. an extrusion/revolution whose basis B-spline is an exact line/circle
   *is* a plane/cylinder/cone/torus. Near-zero is not zero: a rewrite that would move the
   geometry beyond machine precision is refused (the existing recognition contract).

## Data model milestones

- [x] Decide final public names before creating files.
	- Proposed C++ names: `O2BVHSurfaceSolid`, `BoundedSurface`, `SurfaceWire`, `SurfaceEdge`, `PlanarBoundedSurface`, `CylindricalBoundedSurface`, `ConicalBoundedSurface`, `SphericalBoundedSurface`.
	- Deliverable: update this section with accepted names.
	- Validation: none beyond agreement, but keep names consistent with ROOT dictionary style.
	- Done 2026-07-22: first implementation uses public `O2BVHSurfaceSolid`; planar/wire implementation types are private in `Detectors/Base/src/O2BVHSurfaceSolid.cxx`.

- [x] Add a minimal compiling skeleton for the new shape.
	- Files: `Detectors/Base/include/DetectorsBase/O2BVHSurfaceSolid.h`, `Detectors/Base/src/O2BVHSurfaceSolid.cxx`, `Detectors/Base/CMakeLists.txt`.
	- Implement constructor, destructor, `ComputeBBox`, empty `CloseShape`, ROOT `ClassDefOverride`, and dictionary registration.
	- Keep copy assignment and copy construction deleted unless a real deep-copy policy is implemented.
	- Validation: build the `DetectorsBase` target and load the dictionary in ROOT.
	- Done 2026-07-22: skeleton plus first planar implementation added and object-compiled with the configured O2 build flags; full CMake target build is blocked by the unrelated non-tested macro report noted below.

- [x] Introduce private bounded-surface, wire, and edge interfaces used by the solid implementation.
	- Keep it private to `Detectors/Base/src` unless another package has a clear need for it.
	- `BoundedSurface` owns one analytic support surface plus trim wires in that surface's parametric domain.
	- `SurfaceWire` represents one closed oriented boundary loop: one outer wire plus zero or more inner wires for holes.
	- `SurfaceEdge` represents one bounded curve segment in the surface parameter domain and, when needed, its 3D embedding.
	- Required surface methods: conservative AABB, ray intersection candidates filtered by wires, oriented normal at an intersection, signed or unsigned distance for safety, visualization triangulation, and exact/approximate capacity contribution flag.
	- Deliverable: interface plus one trivial test/dummy surface used only by unit tests.
	- Validation: build and a focused unit test that constructs the solid with one dummy surface.
	- Partial 2026-07-22: private concrete `SurfaceWire` and `PlanarBoundedSurface` types exist; a generic bounded-surface interface and dummy test surface remain open.
	- Done 2026-07-23: extracted a private header `Detectors/Base/src/BoundedSurface.h` (namespace `o2::base::surface`) with `SurfaceEdge`, `SurfaceWire` (role + status), the abstract `BoundedSurface` interface (`conservativeBounds`, `intersectRay` returning the oriented hit normal, `distanceSqToPatch`, `normalAt`, `containsPointOnSurface`, `capacityContribution` + `capacityIsExact`, `appendDisplayMesh`, `appendDirectedEdges`), `PlanarBoundedSurface`, and a trivial `DummyBoundedSurface`. `O2BVHSurfaceSolid` now stores `std::vector<std::unique_ptr<BoundedSurface>>` and navigates polymorphically. Test `DummyBoundedSurfaceInterface` exercises the interface.

- [x] Define numerical conventions in code and tests.
	- Choose tolerances for ray `t`, boundary checks, duplicate intersection clustering, and point-on-surface classification.
	- Record these constants in one implementation-local namespace, not scattered literals.
	- Validation: unit tests for near-boundary points and near-tangent rays.
	- Partial 2026-07-22: implementation-local constants were introduced; dedicated near-boundary and near-tangent tests remain open.
	- Done 2026-07-23: added a `NumericalConventions` test exercising near-boundary point classification (inside/outside/boundary across `kTolerance`), near-tangent ray behaviour on a `PlanarBoundedSurface` (grazing miss vs. steep hit, and a hit at the `minDistance` boundary), and `sameIntersection` clustering across `kIntersectionTolerance`.

- [x] Implement the initial wire/edge data model before concrete analytic surfaces.
	- Start with oriented closed wires made from line-segment edges in a 2D parameter domain.
	- Operations needed by all surfaces: closure validation, orientation/area, point-in-wire classification, edge distance/projection, parametric AABB, duplicate vertex cleanup within tolerance, and visualization sampling.
	- Keep wire classification independent of planar surfaces so cylinders, spheres, cones, and future surfaces reuse the same trimming model.
	- Acceptance test: outer square wire, square-with-hole wires, reversed wire, open wire, and point-on-edge cases.
	- Partial 2026-07-22: line-segment closed wires support cleanup, area, point classification, and edge distance; the dedicated wire fixture tests are still pending.
	- Done 2026-07-23: added the remaining reusable, surface-independent primitives to `SurfaceEdge`/`SurfaceWire` (`SurfaceEdge::closestPoint` returning the projected point + clamped parameter, `SurfaceWire::parametricBounds` for the 2D AABB, `SurfaceWire::sampledBoundary` returning the ordered closed boundary polyline as a hook for future curved edges and mesh-independent visualization). Added a `WireDataModel` test covering the outer square (area/orientation/AABB/boundary sampling), a square-with-hole `PlanarBoundedSurface` (point in hole outside, in material inside, on hole edge boundary), reversed wire, open edge list, and point-on-edge / edge-projection cases.

## Surface representation milestones

- [x] Implement bounded planar surfaces first.
	- Model: infinite plane frame plus one outer trim wire and optional inner trim wires in local 2D coordinates.
	- First supported boundary curves: straight line segments. More exact curve types can follow.
	- Required kernels: AABB, ray-plane intersection, wire-domain containment, closest point/distance to wire edges, normal, visualization triangulation.
	- Acceptance test: build an exact box from six planar surfaces and compare navigation against `TGeoBBox`.
	- Done 2026-07-22: initial planar kernels, display triangulation, loop-based navigation, and exact box comparison test added. BVH acceleration remains a later milestone.

- [x] Implement wire orientation and closure validation for planar surfaces.
	- Detect open wires, self-evident degeneracies, zero area, and orientation inconsistent with solid outward normals.
	- Do not silently repair complex invalid input unless the repair is simple and logged.
	- Acceptance test: intentionally reversed face, missing face, and duplicated edge fixtures.
	- Done 2026-07-23: `SurfaceWire::initialize` now takes a `WireRole` and returns a `WireStatus`, rejecting non-finite, too-few-vertex/open, zero-area and self-touching (pinched) wires, and normalizing orientation (outer CCW, inner CW) as a logged simple repair (`O2BVHSurfaceSolid` warns when a wire was re-oriented). Solid-level `validateClosure` runs in `CloseShape` using a directed half-edge (manifold) check plus signed-volume sanity, exposed via `IsClosed()` / `IsOrientationConsistent()` and emitted as warnings. Tests `WireValidationAndOrientation` and `SolidClosureDetectsMissingAndReversedFaces` cover reversed/open/zero-area/pinched wires and missing/reversed face fixtures.

- [x] Add curve classes for trimmed surface boundaries.
	- Start with line and circle arc boundaries, because these cover many mechanical STEP faces.
	- Required operations: endpoint, tangent, bounding box contribution, point projection, robust 2D winding/crossing support.
	- Acceptance test: planar disk or annulus represented by circular trim curves.
	- Done 2026-07-23: added a surface-independent `Curve2D` value type in `BoundedSurface.h` carrying either a straight line segment or a circular arc (centre/radius/signed angular sweep, full circle = +/-2pi) with all required exact operations: `startPoint`/`endPoint`, `tangentAt`, exact `extendBounds` (arc includes the cardinal axis-extreme points inside the sweep), `closestPoint`/`distanceSq`, exact `signedAreaContribution` (Green's theorem) and robust `rightwardCrossings` (arcs are split into v-monotonic sub-arcs so tangent extrema and shared endpoints are not double-counted). Added a `CurveWire` closed-loop type reusing the wire role/status/orientation conventions with exact `signedArea`, `parametricBounds`, winding-based `classify`, and `sampledBoundary` (arc chord sampling). Added a `TrimmedCurveBoundaries` unit test covering line/arc endpoint/tangent/bbox/projection, a full-circle disk (exact pi*r^2 area, inside/outside/boundary), a re-oriented clockwise outer circle, an annulus (outer CCW + inner CW hole, net pi*(R^2-r^2) area, material vs. hole classification), a mixed line+arc half-disk, and an open-loop rejection. Curve areas and classifications verified numerically.

- [x] Implement cylindrical bounded surfaces.
	- Model: axis frame, radius, parametric `u` angle and `v` height/domain, trim wires in `(u, v)`.
	- Required kernels: analytic ray-cylinder intersection, trim-domain check, normal, distance estimate, AABB, visualization mesh.
	- Acceptance test: exact closed cylinder using two planar caps plus one cylindrical surface, compared against `TGeoTube`.
	- Done 2026-07-23: `CylindricalBoundedSurface` with an orthonormal frame (centre + axis + reference `U` axis), radius, and a parametric-rectangle trim (`phi` sweep x height range) matching `TGeoTube`/`TGeoTubeSeg`; general wires on the periodic `(phi, h)` domain are deferred. Kernels: quadratic ray intersection reporting both roots (tangent double-roots suppressed for parity), exact separable patch distance (2D point/segment in the `(rho, h)` half-plane in-sweep, 3D seam-segment distance off-sweep), signed outer/inner-wall normal, exact closed-form divergence-theorem capacity, rim-circle AABB, grid display mesh, and oriented sampled boundary edges. Tests: `CylindricalSurfaceKernels`, `ClosedCylinderMatchesTGeoTube` (with disk caps), `HollowCylinderMatchesTGeoTube` (outer + inner wall + annulus caps).

- [x] Implement spherical bounded surfaces.
	- Model: center, radius, parametric domain and trim wires.
	- Required kernels: analytic ray-sphere intersection, normal, distance estimate, AABB, visualization mesh.
	- Acceptance test: full sphere or sphere section compared against `TGeoSphere` where applicable.
	- Done 2026-07-23: `SphericalBoundedSurface` with centre, polar axis frame, radius and a `TGeoSphere`-style parametric trim (polar `theta` range x `phi` sweep); the full sphere is self-closing (no boundary edges). Kernels: quadratic ray intersection with trim filtering, radial signed normal, exact divergence-theorem capacity for arbitrary sections, centre+-radius AABB, lat/long display mesh, oriented boundary edges for sections (pole rims and full-sweep seams skipped). Patch distance is exact for the full sphere and for points whose direction lies in the trim; otherwise the full-sphere distance is returned as a documented conservative lower bound (keeps `Safety` safe). Tests: `SphereMatchesTGeoSphere` (closed-solid comparison) and `SphericalSectionKernels` (hemisphere trim/capacity/intersections).

- [x] Implement conical bounded surfaces.
	- Model: apex or reference point, axis, opening angle/radii, parametric domain and trim wires.
	- Required kernels: analytic ray-cone intersection, normal, distance estimate, AABB, visualization mesh.
	- Acceptance test: cone or truncated cone compared against an equivalent ROOT shape if available.
	- Done 2026-07-23: `ConicalBoundedSurface` with a linear radius law `r(h)` (radii given at the two height bounds, zero allowed at one end for apex cones; slope 0 degenerates to a cylinder) and the same frame/trim/inner-wall conventions as the cylinder. Kernels: quadratic ray intersection with mirror-nappe rejection and linear fallback for degenerate quadratics, slanted signed normal, exact separable patch distance against the straight generator segment, exact divergence-theorem capacity (verified against the truncated-cone volume formula), rim AABB, display mesh and boundary edges with apex-rim degeneration handled (an apex cone closes against a single cap). Tests: `TruncatedConeMatchesTGeoCone` and `ApexConeClosesWithSingleCap`. Known measure-zero limitation for later hard-case work: a ray exactly along the axis through the apex is a tangent double root and reports no lateral crossing.

- [ ] Decide how to handle torus, Bezier, BSpline, and other high-order CAD surfaces.
	- Conservative default: leave unsupported surfaces on the existing tessellated fallback path.
	- Optional later path: add bounded torus analytically; keep BSpline/Bezier tessellated unless a robust exact kernel is designed.
	- Deliverable: converter support matrix documented in this file and in the script help.
	- Decision 2026-07-24 (informed by the ALICE3_CAD_pure.step real-CAD test): the umbrella above is
	  split into three concrete, independently-shippable milestones below, ordered by value/cost as measured
	  on that assembly (55 unique unplaced parts; only 3/55 convert today). The measured blocker ranking:
	  B-spline dominates (44/55 volumes touch it, 34/55 are blocked *only* by it), then torus and
	  quadric-curved trims (15 volumes each). Crucially the B-spline story separates into **trim curves**
	  (cheap, surface stays exact) vs **surfaces** (expensive, exactness-vs-tessellation trade-off) — do the
	  first before the second.

- [x] **B-spline trim curves on analytic surfaces** (highest value/cost ratio — do first). SCOPED PLAN below.
	- Done 2026-07-24: implemented all five scoped steps. Kernel (`BoundedSurface.h`): a third `Curve2D`
	  kind `BSpline` carrying degree/poles/weights/clamped-flat-knots, with de Boor (`bsplineEval` via
	  Piegl DersBasisFuns, rational-aware) evaluation, clamped endpoints, control-hull `extendBounds`,
	  exact Gauss-Legendre-per-knot-span `signedAreaContribution` (non-rational exact, rational numeric),
	  robust `rightwardCrossings` via an adaptively-flattened on-curve polyline reusing the existing
	  half-open canonical-endpoint seam convention, polyline-based `closestPoint`/`distanceSq`, and
	  `reverseInPlace` (poles/weights reversed, knots complemented). A **lazily-cached** flattened
	  polyline (`bsplineCache`) makes the dense per-point winding/distance queries (numeric capacity
	  integration, point-in-wire grids) pay the de Boor sampling only once — without it `CloseShape`'s
	  capacity integration re-sampled every bspline per grid point and hung. `CurvedPlanarBoundedSurface`
	  reports `capacityIsExact()==false` when any wire is a bspline. Public POD `PlanarBoundaryCurve`
	  gained a `makeBSpline` case (reused by the planar and quadric wire overloads). Sidecar: `curveType=2`
	  (`O2SurfaceSolidIO.cxx` `parseBSplineEdge`); relaxed the reader `kJoinTolerance` to `1e-5` AND added
	  a matching kernel `kWireJoinTolerance = 1e-5` for `CurveWire` closure (the reader relax alone was
	  insufficient — the kernel's `1e-9` `CurveWire` closure was the real gate for mixed line/bspline
	  loops, incl. the pre-existing D-cap bug). Python (`O2_CADtoTGeo.py`): `extract_planar_face` emits
	  bspline/Bezier boundary edges (3D poles projected into the plane frame); the three quadric
	  extractors call a new `_quadric_trim_wire` that keeps line edges as lines and converts every curved
	  pcurve (circle/ellipse/Bezier/bspline) to a bspline whose poles are transformed by the affine
	  `(u,v)->(phi,height/theta)` map (exact; the old arc-on-quadric ellipse gate is thereby closed too);
	  `face_supported` / `_SUPPORTED_*_CURVES` reconciled. Tests: `BSplineTrimCurveKernels`,
	  `BSplineSidecarRoundTrip`, `BSplineWindowInCylinderWall` (26 cases, `ctest -R BVHSurfaceSolid`
	  green). **End-to-end**: `ST1A38495_01` and `ST1A38526_01` (previously blocked purely by bspline
	  trims) now convert under `--exact-surfaces auto`; each loads via `LoadSurfaceSolid` and matches the
	  accelerated `O2Tessellated` `Contains` with **0 mismatches over 20k random points** (and BVH==loop),
	  capacity within 0.17-0.24% of the divergence-theorem mesh volume. **ALICE3 coverage 3/55 -> 7/55.**
	  Known: `IsClosed()` may warn (two exact faces sharing one 3D bspline edge carry different pcurve
	  parametrizations, so equal-count chord sampling gives slightly different 3D points — the documented
	  scope caveat; navigation stays exact). Remaining ceiling is bspline *surfaces* (34/55 blocked only by
	  them) and torus.
	- Motivation: on real CAD, whole all-analytic parts are blocked purely by *trim* B-splines, not
	  B-spline surfaces. Example studied: `ST1A38495_01` (lid `0:1:1:10`, 65 faces) and `ST1A38526_01`
	  (lid `0:1:1:11`, 53 faces) — every surface is plane/cylinder/cone, yet each is blocked by ~20 faces
	  whose trim wires contain a B-spline: bspline planar boundary loops (`{line:6, bspline:3}`), and
	  cylinders/cones bounded by a bspline pcurve (inner hole `{bspline:2, circle:2}`, outer `{bspline:2,
	  line:2}`). **0** of the blocking quadric faces there are pure line/arc — so the earlier "quadric
	  arc-pcurve" fix does not rescue them; a bspline *trim-edge* type does.
	- Key insight #1: the ray/surface hit stays **fully analytic** (plane/quadric). Only the 2D trimming test
	  and boundary sampling need the B-spline. Point-in-wire reduces to a scalar root-find `C_v(t) - v0 = 0`
	  in parameter space, **not** a 3D ray/surface intersection.
	- Key insight #2 (de-risks the quadric path AND subsumes the old arc-on-quadric gap): the quadric
	  `(U,V)→(phi,h)`/`(phi,theta)` remap is **affine** for all three families (cylinder `h=V`; cone
	  `h=V·cosα`; sphere `theta=π/2−V`; `phi=±U+offset`). A B-spline is **closed under affine maps** — just
	  transform its control-point poles (same degree/knots/weights). So a bspline pcurve maps to a bspline in
	  `(phi,h)` *exactly* (no ellipse problem — that problem was specific to the *arc* `Curve2D` type, which an
	  anisotropic map turns into an ellipse). Bonus: convert *arc* pcurves on quadrics to bsplines
	  (`Geom2dConvert`) and transform poles too, so this one edge type also closes the 2297-face
	  "quadric trim edge pcurve not a straight line" fallback — verify the affine claim by round-tripping a
	  sampled 3D point.

	SCOPED IMPLEMENTATION PLAN (bottom-up; each step builds + tests before the next):

	1. **Kernel — 2D B-spline `Curve2D` variant** in `Detectors/Base/src/BoundedSurface.h`.
	   - Add a third `Curve2D` kind (line=0, arc=1, bspline=2) carrying degree `p`, poles `(u,v)[]`, weights
	     `w[]` (all 1 ⇒ non-rational), and a clamped flat knot vector (length `nPoles+p+1`). Bézier is the
	     single-span special case.
	   - Implement the ops `CurveWire` already calls, matching existing conventions:
	     - `startPoint`/`endPoint`: de Boor evaluation at the parameter ends.
	     - `extendBounds`: conservative AABB = control-point convex-hull box (exact bound; matches the
	       conservative-box philosophy of the BVH).
	     - `signedAreaContribution`: `∮(u dv − v du)/2` by Gauss–Legendre per Bézier span
	       (`ceil((2p)/2)+1` points is exact for non-rational degree `p`; rational ⇒ numeric, flagged).
	     - **`rightwardCrossings(point)`** (the core): decompose to Bézier spans (knot insertion); per span
	       cull by control-point convex hull (all `v` one side of `v0` ⇒ skip); solve `C_v(t)=v0` on survivors
	       (cubic ⇒ Cardano/3×3 companion; general ⇒ recursive de Casteljau / Bézier clipping); count roots
	       with `C_u(t*)>u0` using the sign of `dC_v/dt` and the **existing half-open canonical-endpoint
	       convention** in `CurveWire::classify` (reuse the shared-loop-endpoint trick to avoid seam
	       double-counts).
	     - `closestPoint`/`distanceSq`: conservative **lower bound** (distance to the control-point convex
	       hull), consistent with the wire-trimmed-quadric `Safety` policy; exact boundary distance deferred.
	     - `sampledBoundary`/chord sampling: adaptive de Casteljau to a flatness tolerance (visualization +
	       closure half-edge check).
	   - Unit test `BSplineTrimCurveKernels` (isolated): a planar bspline outer loop — `classify`
	     inside/outside/boundary, area vs a fine-polygon reference, crossings incl. a horizontal-tangent case.

	2. **Public API** — extend the `PlanarBoundaryCurve` POD (reused by both planar and quadric `(u,v)`
	   overloads) with a bspline case (`int kind`, `int degree`, `std::vector<std::array<double,2>> poles`,
	   `std::vector<double> weights, knots`). No new `Add*Surface` entry points needed — the existing
	   `AddCurvedPlanarSurface` and the quadric wire overloads already take `PlanarBoundaryCurve` wires.

	3. **Sidecar + reader** (`O2SurfaceSolidIO.cxx`, format section above). Add `curveType = 2` (bspline2d),
	   record = `degree, nPoles` then poles `2*nPoles`, weights `nPoles`, flat knots `nPoles+degree+1` — packed
	   in the self-describing `curveParams[]` (bump not required; `nCurveParams` stays generic). `wireToCurves`
	   parses it into the kernel `Curve2D`. **Also relax `kJoinTolerance`** from `1e-9` here: the `{line:…,
	   bspline:…}` mixed loops in these parts need line/bspline endpoints to join within the extractor's
	   ~`1e-6`, and this is the same pre-existing D-cap bug — fold the one-line relaxation in as a prerequisite.
	   Test `BSplineSidecarRoundTrip`: write a bspline wire block, read via `LoadSurfaceSolid`, compare
	   `Contains`/area to the in-memory kernel; confirm line/arc round-trips stay byte-identical.

	4. **Python extractor** (`O2_CADtoTGeo.py`). Emit bspline pcurve records from `Geom2d` in **both** paths:
	   - `extract_planar_face`: bspline planar boundary edges (fixes `plane with unsupported boundary curves`).
	   - `_quadric_line_wire`: bspline (and arc, via `Geom2dConvert` → bspline) pcurves, applying the **affine
	     `(U,V)→(phi,h)/(phi,theta)` pole transform** (same map used for the scalar bounds, incl. the
	     left-handed-frame `phi` mirror). Get the pcurve via `BRep_Tool.CurveOnSurface`
	     (`Geom2d_BSplineCurve`: degree, poles, weights, knots+mults → flat vector; Bezier via
	     `Geom2dConvert`). Remove the "arc pcurve / bspline falls back" gates for the trim path.
	   - Reconcile `_SUPPORTED_PLANAR_CURVES` / `face_supported` so the eligibility report matches the extractor
	     (bspline trims now supported; bspline *surfaces* still not).

	5. **End-to-end validation** (reuse the ALICE3 harness). Convert `ST1A38495_01` and `ST1A38526_01` under
	   `--exact-surfaces auto`; load each via `LoadSurfaceSolid`; compare `Contains` against the accelerated
	   `O2Tessellated` (NOT ROOT `TGeoTessellated`) over a grid + fixed-seed random set, and capacity against
	   the divergence-theorem mesh volume; expect boundary-band-only mismatches. Add a synthetic
	   `BSplineHoleInCylinder` fixture (cylinder with a bspline-bounded window) with an analytic/tessellated
	   reference. Re-run the ALICE3 sweep and record the new coverage number.

	SCOPE BOUNDARIES (explicitly out of Point 1 — keep the change focused):
	- B-spline *surfaces* are untouched (that is Point 3, and carries the exactness-vs-tessellation caveat).
	- **Capacity** of a bspline-trimmed face is numerically integrated → `capacityIsExact()==false` (same
	  policy already used for wire-trimmed quadrics); **`Safety`/distance** returns a conservative lower bound.
	- **Closure diagnostic caveat:** two exact faces sharing one 3D bspline edge each carry their *own* pcurve
	  (different parametrizations of the same 3D curve), so equal-count chord sampling can give slightly
	  different 3D points and the `IsClosed`/half-edge check may *warn*. Navigation (`Contains`) is exact
	  regardless — it uses the analytic hit + the exact root-find trim test, never the sampling. Accept the
	  warning initially; consistent shared-3D-edge sampling is a later refinement if the warnings are noisy.
	- Validation: `ST1A38495_01` converts under `--exact-surfaces auto` and matches `O2Tessellated` `Contains`
	  to the mesh-precision band; `ninja` + `ctest -R BVHSurfaceSolid` green; line/arc output byte-identical.

- [x] **Bounded torus** (analytic; second-largest single analytic blocker).
	- Motivation: 15 ALICE3 volumes touch a torus; near-miss `ST2487455_002` (lid `0:1:1:8`) is a single
	  torus face away from exact (cyl + 2 cone + 2 plane + 1 torus).
	- Approach: a `TorusBoundedSurface` with a `(phi_ring, phi_tube)` parametric domain and the same optional
	  `CurveWire` trim + inner-wall conventions as the quadrics. Ray/torus is a quartic — use a stable quartic
	  solver (closed-form Ferrari or a numerically-robust companion-matrix root finder); mind grazing/tangent
	  double roots for `Contains` parity, as already handled for cylinder/cone.
	- Done 2026-07-24: implemented `TorusBoundedSurface` in `BoundedSurface.h` (major/minor radius, frame,
	  scalar `phiRing x phiTube` rectangle trim, optional `CurveWire` trim + inner-wall). Kernel: a shared
	  `solveQuarticReal` (Ferrari via resolvent cubic + Newton polishing; `solveDepressedCubic` helper) —
	  the ray/torus leading coefficient is `|dir|^4 > 0` so it is always a genuine quartic; real roots are
	  clustered by `sameIntersection` and **even-sized clusters (tangencies) are dropped** so crossing parity
	  stays consistent (the torus analogue of the cylinder/cone double-root suppression). Signed radial
	  (away-from-tube-spine) normal; **exact** divergence-theorem capacity over the rectangle (closed form
	  below; verified `= 2*pi^2*R*r^2` for the full torus), numeric for a wire trim (`capacityIsExact()==false`).
	  `distanceSqToPatch` is the exact `(rho, z)` meridian distance to the tube circle for the full torus and a
	  conservative lower bound for a trimmed patch (keeps `Safety` safe — sphere policy). Full-torus AABB
	  (conservative for partial sweeps); grid display mesh; directed edges for the closure half-edge check
	  (self-closing when both sweeps are full). Both ring and tube angles are periodic, so `pointInTrim`
	  unwraps **both** into their wire windows (unlike the cylinder/cone/sphere, whose second parameter is not
	  periodic); a trim wrapping more than a full turn in either angle is rejected. Public
	  `AddToroidalSurface` (+ wire-trim overload) reusing the `PlanarBoundaryCurve` POD; sidecar `surfaceType=5`
	  (15 params) in `O2SurfaceSolidIO.cxx`; Python `extract_toroidal_face` (OCC `U=phiRing`, `V=phiTube`; the
	  affine `(u,v)->(phiRing,phiTube)` map is identity/mirror so it reuses `_quadric_trim_wire` and
	  `_quadric_trim_fills_uv_box`), `"torus"` added to `SURFACE_TYPE_ENUM`/`_SUPPORTED_SURFACE_TYPES`/
	  `_FACE_EXTRACTORS`. Tests (`testBVHSurfaceSolid.cxx`, now 30 cases): `ToroidalSurfaceKernels`
	  (4-crossing ray, tangent graze, normals, meridian distance, exact capacity, quarter-tube section),
	  `FullTorusMatchesTGeoTorus` (closed-solid grid + distances + capacity vs `TGeoTorus`),
	  `WireTrimmedTorusMatchesSection` (rectangle-vs-wire equivalence, numeric capacity), `TorusSidecarRoundTrip`
	  (`surfaceType=5` reader vs `TGeoTorus`). `ninja` + `ctest -R BVHSurfaceSolid` green. End-to-end: a
	  synthetic full-torus STEP converts under `--exact-surfaces required` (1/1 exact) and its Python-written
	  sidecar loads via `LoadSurfaceSolid` with **0 `Contains` mismatches over 9261 points vs `TGeoTorus`** and
	  exact capacity. Capacity closed form for the `(phiRing u, phiTube v)` rectangle (centre `C`,
	  `cu=C.U, cv=C.V, cw=C.W`): `contribution = mNormalSign/3 * [ du*Iv + du*cw*Iw + (cu*Su+cv*Cu)*Iuv ]`
	  where `Su=sin u1-sin u0`, `Cu=cos u0-cos u1`, `Iv = r*((R^2+r^2)*Sv + R*r*dv + R*r*C2v)`,
	  `Iw = r*(R*(cos v0-cos v1) + r*(cos2v0-cos2v1)/4)`, `Iuv = r*(R*Sv + r*C2v)`, `Sv=sin v1-sin v0`,
	  `C2v=dv/2+(sin2v1-sin2v0)/4`. **ALICE3 coverage re-measured 2026-07-26: 7/55 -> 15/55** (see the
	  2026-07-26 handoff note). All 8 newly-converting volumes contain torus faces, including the predicted
	  near-miss `ST2487455_002` (exactly the forecast 1 cyl + 2 cone + 2 plane + 1 torus).

- [x] **Canonical-form recognition: recover the exact analytic model behind a stored NURBS** (highest value on the board — measured +14/55 on ALICE3 and +5/5 on as1; do before the general NURBS milestone).
	- **Premise, established by measurement 2026-07-26 (see the handoff note): the surface type stored in a
	  STEP file is a statement about the exporter, not about the geometry.** CAD kernels routinely write an
	  exact cylinder, cone, sphere or circle as a *rational* B-spline/Bezier patch — that representation is
	  exact, not an approximation, because a rational quadratic/cubic reproduces conics exactly. The current
	  extractor dispatches on `BRepAdaptor_Surface::GetType()` and therefore throws away geometry it fully
	  supports. This is the single biggest source of unnecessary tessellated fallback we have found.
	- Measured prize (relative residual < 1e-9, i.e. machine precision):
	  - `as1-oc-214.stp`: **0/5 -> 5/5**; all 70 bspline faces are exact cylinders, zero free-form.
	  - `ALICE3_CAD_pure.step`: **15/55 -> 29/55**; 2889 of 8664 non-analytic faces are exact quadrics
	    (1419 cylinder, 1398 cone, 72 sphere). Lower bound — torus recognition is not yet in the classifier.
	- Design principle (as the request put it: *simplify to the smallest/best-fitting model*): this is **model
	  selection**, not fitting. Try candidate models in increasing parameter count — plane (3) < sphere (4) <
	  cylinder (5) < cone (6) < torus (7) < free-form — and accept the **first** whose residual is at machine
	  precision. Never accept a looser fit silently: an "almost cylinder" that is really free-form must stay
	  free-form, or the converter would silently change the geometry it was built to represent exactly.
	- Recommended method — **differential, not algebraic**: classify from the surface normal field, which is
	  cheap, needs no initial guess, and each test is a small linear solve on a sampled grid (validated in the
	  session scratchpad prototype `recognize.py`):
	  - plane: all unit normals parallel.
	  - sphere: normal lines concurrent — solve `P_i = C + r*N_i` for `(C, r)` in least squares.
	  - cylinder: normals coplanar; the axis is the smallest right singular vector of the normal matrix;
	    then project out the axis and fit a circle, requiring constant radius.
	  - cone: `N_i . (P_i - A) = 0` is **linear in the apex** `A`; solve, then require a constant half-angle
	    about the mean ruling direction.
	  - torus: needs a second step (fit the spine circle in the meridian half-plane); add it after the
	    quadrics, since it only tightens the reported bound.
	- Do the **same for curves** (the request explicitly says curves too, and it is the cheaper half): a bspline
	  trim edge that is an exact line/circle/ellipse should be recovered as one. Signature weights
	  (`1, 1/2, 1` for a 120-degree rational quadratic arc; `1, 1/3, 1/3, 1` for a rational cubic semicircle)
	  are a useful fast path, but a numeric fit is more general and should be the authority. Note this is
	  *demotion for exactness/compactness*, not enablement — the bspline `Curve2D` kernel already handles these
	  correctly, so the win here is smaller (exact capacity, tighter `Safety`, fewer closure warnings) and it
	  should not block the surface work.
	- Do **not** use OCCT's `BRepLib_CanonicalRecognition` (7.7+): verified 2026-07-26 that it is not exposed in
	  this pythonOCC build (`ImportError` from `OCC.Core.BRepLib`). Write our own numeric recognizer in
	  `O2_CADtoTGeo.py`; it is ~100 lines and we then control the exactness contract.
	- Placement: a pre-pass in the Python extractor that rewrites a recognized face's *surface* record to the
	  analytic type, keeping its existing trim wires (the pcurves are in `(u,v)` of the *stored* surface, so the
	  reparametrization to `(phi, h)`/`(phi, theta)` must be derived from the recognized frame — this is the
	  real work, and the affine-pole-transform machinery in `_quadric_trim_wire` is the right tool). **No C++
	  kernel, sidecar-format or reader change should be needed** — the output is an ordinary quadric record.
	  Verify that before coding rather than assuming it.
	- Report it: add the recognized-vs-stored breakdown to `--surface-report` so the eligibility report stops
	  under-reporting, and add a flag (suggested `--recognize-surfaces exact|off`, default `exact`) plus a
	  per-face residual in the JSON so a reviewer can audit what was reclassified and how well it fit.
	- Validation: `as1-oc-214.stp` reaches 5/5 and `ALICE3_CAD_pure.step` >= 29/55 under `--exact-surfaces auto`;
	  each newly-converting volume loads via `LoadSurfaceSolid` and matches the accelerated `O2Tessellated`
	  `Contains` (**not** ROOT `TGeoTessellated`) to the mesh-precision band; a unit fixture asserts that a
	  deliberately free-form patch is **rejected** (guarding the exactness contract, which matters more than the
	  positive cases); `ninja` + `ctest -R BVHSurfaceSolid` green.
	- Done 2026-07-26 (surface recognition pre-pass; curve recognition deferred — see gaps below).
	  Implemented entirely in `O2_CADtoTGeo.py` (no C++/kernel/sidecar/reader change — verified, not
	  assumed, exactly as scoped): `_recognize_analytic_surface` is the plane/sphere/cylinder/cone
	  differential-normal-field recognizer (ported from the validated `analyze_surface_geometry.py`
	  prototype, same `1e-9` relative-residual machine-precision gate, no torus test — matching the
	  documented lower bound). `recognize_and_extract_face` is the pre-pass: for a face whose *stored*
	  type has no direct extractor (bspline/bezier/revolution/extrusion), it recognizes the frame and
	  re-derives the trim **from the 3D boundary edges directly** (not the stored pcurves) via
	  `_recognized_quadric_wire_block` — sample each edge's own 3D curve, project through the recognized
	  frame's closed-form inverse (`point -> (phi, h)` for cylinder/cone, `(phi, theta)` for sphere), and
	  accept only edges that come out *exactly* axis-aligned (constant-`h`/`theta` = a rim/cap, or
	  constant-`phi` = a generator/meridian) in that domain. Every accepted edge is then a straight line
	  segment in `(phi, h)` by construction, so **the affine-pole-transform machinery in
	  `_quadric_trim_wire` turned out not to be reusable and was not needed**: that machinery reparametrizes
	  a *curved* pcurve under an affine map, but the map from a recognized bspline surface's own `(u, v)` to
	  `(phi, h)` is only ever *separable*, not affine (verified numerically — e.g. a NURBS circular arc's
	  parameter maps to `phi` via a Möbius-type nonlinear function, not a linear one) — the milestone's
	  suggestion to reuse it was checked and found not to hold; sampling in the recognized domain directly
	  sidesteps the question entirely and is simpler. A recognized plane reuses `extract_planar_face`
	  unchanged via a new `frame_override` parameter, since that function already builds its wires from 3D
	  edges (not pcurves), so no separate machinery was needed there at all.
	  Continuous per-sample phi-unwrapping (not a coarser vertex-to-vertex or grid unwrap) was required
	  for correctness: a first version unwrapped only at wire vertices and produced a **17*pi** unwrap
	  error on one as1 face, because a rim edge's own sweep can approach pi and a non-uniform rational
	  parametrization can place two *adjacent* vertices' raw delta-phi right at the +/-pi branch cut;
	  densely sampling *along* each edge keeps every consecutive step small and fixed it (also replaced a
	  separate coarse-grid "conservative window" computation with the wire's own dense-sample extent, more
	  robust and simpler). A cone/sphere apex or pole (a degenerate edge: single point, no 3D curve, common
	  on real CAD - countersinks, chamfers) is handled by carrying the running unwrapped phi through
	  unchanged at that point (phi is numerically indeterminate there - atan2 of a near-zero radial vector -
	  but the point trivially satisfies "iso" in both directions), rather than hard-rejecting the face.
	  Report: `classify_face`/`face_supported`/`build_surface_report` gained `recognize_surfaces` (default
	  on) and record `recognized_type`/`recognized_residual` per face plus a
	  `recognized_surface_counts`/`recognized_stored_type_counts` summary breakdown; new CLI flag
	  `--recognize-surfaces exact|off` (default `exact`) threads through both the report and
	  `--exact-surfaces auto|required`.
	  **Measured** (env: `ALIBUILD_WORK_DIR=/data/swenzel/sw`, usual pythonOCC `PYTHONPATH`/OCCT
	  `LD_LIBRARY_PATH` fixes; `python3 O2_CADtoTGeo.py <file> --exact-surfaces auto --mesh
	  --surface-report <path>`):
	  - `as1-oc-214.stp`: **0/5 -> 5/5**, exactly as forecast. All 28 bspline faces recognized as exact
	    cylinders. `Contains` vs `O2Tessellated`: 0/20000 mismatches on 4/5 volumes, 2/20000 (boundary-band,
	    max safety 9e-5) on the 5th; BVH == `Contains_Loop` everywhere.
	  - `ALICE3_CAD_pure.step`: **15/55 -> 20/55 extracted, 19/55 usable** — short of the 29/55 forecast.
	    1180 of the 8664 non-analytic faces recognized (786 cylinder, 358 cone, 36 sphere; the forecast's
	    2889 was measured on deduplicated `TopoDS` solids at a coarser residual gate and, per its own
	    documented caveat, was *"a forecast to be confirmed by an actual run"* — this is that run, and the
	    honest number is lower). The gap is almost entirely **non-iso trims**: a recognized face whose
	    boundary is not exactly a rim/generator in the recognized `(phi, h)` domain is rejected by design
	    (`373` -> `278+58+36+1 = 373` faces fall on this reason after the apex fix below) — out of scope for
	    this pass per the milestone's own stated boundary ("a genuinely slanted/curved cut... stays on the
	    tessellated fallback"). Fixed one sub-case discovered during measurement: cone-apex degenerate
	    edges were initially hard-rejected (99 faces); after the fix they no longer hard-fail, but in this
	    dataset **all 99 turned out to have a second, independent non-iso edge** elsewhere on the same face,
	    so net ALICE3 coverage was unchanged by that fix (20/55 before and after) — kept anyway since it is
	    correct and unit-relevant for CAD that does not also have this second issue (verified no regression
	    on as1 or the other three sweep files).
	  - `Bagger.step`/`oTOF System V3-R92cm.step`: unchanged (12/13, 2/3) — confirmed no regression; neither
	    file's fallback faces are recognizable quadrics-in-disguise (Bagger's gap is ellipse *trim curves* on
	    an already-recognized-correctly plane, ellipse remains a separate open milestone below).
	  - Sidecar load-check (`checkSurfaceSidecars.macro`) on the final ALICE3 output: 19/20 load (the 20th,
	    `ST1829909_01`, is the **pre-existing, unrelated** `kJoinTolerance` gap documented in the 2026-07-26
	    handoff note above — not touched by this session). 1/19 loaded volumes
	    (`ST0923290_013`, 9 recognized inner-wall cylinder holes) has `IsOrientationConsistent() == false`
	    (8 reversed boundary edges) — root-caused, not chased further, see gap below. All others closed or
	    warn only in the already-documented shared-3D-edge chord-sampling sense.
	  - `Contains` cross-validated (accelerated `O2Tessellated`, not ROOT `TGeoTessellated`) on all 5 newly
	    recognized ALICE3 volumes: 0/10000 mismatches on 4/5; the flagged `ST0923290_013` gave 3/200000
	    mismatches at a higher sample count, each within `5.6e-4` cm of the exact boundary via `Safety` —
	    the mesh-precision band, not a correctness break. `ninja` build (`DetectorsBase` target name changed
	    upstream; used `ctest -R BVHSurfaceSolid` directly) + focused ctest green (30 cases, unaffected as
	    expected — this session made no C++ changes).
	  **New gap discovered (root-caused, deliberately not fixed this session):**
	  `IsOrientationConsistent()` can be false for a recognized face even though `Contains` stays exact.
	  Root cause: a *stored* quadric's wire-traversal order (from OCC's `BRepTools_WireExplorer`) has a
	  reliable external handedness hint — `ax3.Direct()` — that the existing `_quadric_trim_wire` mirrors
	  `phi = -u` against when needed. A *recognized* face has no such hint: the underlying bspline surface's
	  own `(u, v)` handedness is unrelated to the recognized `(phi, h)` frame, and OCC's wire order is
	  relative to the *stored* surface, not the recognized one. Verified this is **not** an axis-sign
	  ambiguity in the recognizer (the radial-outward direction `rel - (rel . axis) * axis` is provably
	  axis-sign-invariant; measured `dot(radial, N)` is a *consistent* `-1.0` across all 9 affected faces,
	  matching their uniform `inner_wall = True`, i.e. individually correct) — it is a **per-face wire
	  winding-sense** question, independent of any per-face frame choice, that would need cross-face
	  consistency propagation across the whole solid to fix in general. Deferred: `Contains`/`Safety`/
	  `DistFrom*` do not depend on cross-face winding (each surface's own ray-parity test is independent of
	  its neighbors), only the closure *diagnostic* is affected — same acceptance already established for
	  the pre-existing "`IsClosed()` may warn on shared bspline edges" caveat.
	  **Remaining scope for a follow-up session:**
	  (1) non-iso trims on a recognized quadric (373 ALICE3 faces) — would need either a numeric
	  re-fit of the boundary curve in `(phi, h)` (e.g. a Bezier-clipping-style projection) or accepting a
	  bspline-in-`(phi,h)` via a *non-affine* numeric fit, both a materially bigger step than this pass;
	  (2) curve recognition (bspline trim edges that are exactly a line/circle/ellipse) — not started, the
	  milestone's "cheaper half"; (3) torus recognition — explicitly deferred by the milestone itself;
	  (4) the cross-face wire-winding consistency gap above. `g4Config.C` and generated
	  `scripts/geometry/STEP_examples/`/`facets_*`/`surfaces_*` artifacts remain intentionally out of the
	  commit (all measurement runs in this session used a scratch output folder, not the repo).

- [ ] **Ellipse (and remaining conic) trim curves on planar faces** (small, Python-only — cheapest open item).
	- Motivation: found by the 2026-07-26 `Bagger.step` sweep (see the handoff note). That model is 12/13
	  exact; the single fallback volume `Bucket` (97 faces: 69 plane, 22 cylinder, 4 sphere, 2 torus — the
	  only volume in the file exercising sphere *and* torus) is blocked **solely** by 8 **ellipse** planar
	  boundary curves. Closing this converts a whole volume and takes that file to 13/13. Elliptical planar
	  boundaries are common wherever a cylinder is cut obliquely, so this should recur on other real CAD.
	- Key insight: the work is almost certainly already done, in the wrong branch. On **quadrics** ellipse
	  pcurves are supported today — `_SUPPORTED_QUADRIC_CURVES` includes `ellipse`, and `_quadric_trim_wire`
	  converts every curved pcurve (circle/ellipse/Bezier/bspline) to a bspline and pushes the poles through
	  the affine `(u,v) -> (phi, h/theta)` map. Only the **planar** path lags: `extract_planar_face` rejects
	  anything outside `(GeomAbs_Line, GeomAbs_Circle, GeomAbs_BSplineCurve, GeomAbs_BezierCurve)`
	  (`O2_CADtoTGeo.py:898`), and `_SUPPORTED_PLANAR_CURVES = {line, circle, bspline, bezier}` matches it.
	- Approach: route ellipse edges into the **existing** bspline branch. `_planar_bspline_edge_params`
	  already calls `geomconvert.CurveToBSplineCurve(Geom_TrimmedCurve(...), Convert_TgtThetaOver2)`, and
	  that conversion is **exact** for conics (a conic is a rational quadratic NURBS; `Convert_TgtThetaOver2`
	  is precisely the conic parameterisation), then projects the 3D poles into the plane frame — an affine
	  map, under which a bspline is closed. So the expected change is only: add `GeomAbs_Ellipse` to the
	  accepted-type tuple, let it fall through to the `else` bspline branch, and add `"ellipse"` to
	  `_SUPPORTED_PLANAR_CURVES` so the eligibility report and the extractor stay reconciled.
	- **No kernel, sidecar-format, or reader change is expected**: the record emitted is an ordinary
	  `curveType=2` bspline. Verify that expectation before writing code rather than assuming it — if it
	  holds, this is a ~3-line Python diff.
	- Consider folding in `GeomAbs_Hyperbola` / `GeomAbs_Parabola` the same way (both are conics that
	  `CurveToBSplineCurve` handles), but only if a real fixture exercises them — an unbounded-branch conic
	  trimmed to a finite edge is fine, an untrimmed one is not. Do not add them speculatively.
	- Caveats inherited from the bspline-trim milestone (state them, do not re-litigate): capacity of an
	  ellipse-trimmed planar face is numerically integrated (`capacityIsExact() == false`), and `Safety` is a
	  conservative lower bound. Navigation stays exact — the surface is still an analytic plane and the trim
	  test is still the exact winding/root-find.
	- Validation: `Bucket` converts under `--exact-surfaces auto` and the whole of `Bagger.step` reaches
	  13/13; the sidecar loads via `LoadSurfaceSolid`; `Contains` compared against the accelerated
	  `O2Tessellated` (**not** ROOT `TGeoTessellated`) over a fixed-seed random set, expecting
	  boundary-band-only mismatches. Add a synthetic unit fixture (an obliquely cut cylinder, whose cap is an
	  exact ellipse, or a plate with an elliptical hole) with an analytic reference, and check line/arc/bspline
	  output stays byte-identical. `ninja` + `ctest -R BVHSurfaceSolid` green.

- [ ] **B-spline / NURBS bounded surfaces** (bulk coverage; largest effort — do last).
	- Motivation: the dominant blocker, and after the torus milestone it is very nearly the *only* one.
	  Re-measured on ALICE3 2026-07-26 (coverage now 15/55): **33 of the 40 remaining fallback volumes are
	  blocked purely by bspline**, and 36/55 touch it (down from 44/55 before the bspline-trim + torus work,
	  because those volumes moved into the exact set). The full remaining blocker list is only three entries
	  wide — bspline (2377 faces), extrusion (73 faces, 4 volumes), revolution (12 faces, 5 volumes) — and a
	  basis-curve probe on 2026-07-26 showed the last two are **bspline in disguise**: all 14 extrusion faces
	  of `ST1A38494_01` and both revolution faces of each of `ST2487459_002/003/004` sweep/revolve a *bspline*
	  basis curve, not a line or circle. So there is no cheap analytic win left to harvest, and **15/55 is the
	  hard analytic ceiling** — everything above it requires this milestone.
	- **CORRECTION 2026-07-26: the "15/55 is the hard analytic ceiling" claim above is WRONG, and so is the
	  reasoning behind it.** It equated *stored surface type* with *geometry*. A measurement that classifies
	  faces by their actual differential geometry instead (see the canonical-recognition milestone below)
	  shows **2889 of ALICE3's 8664 non-analytic faces (33%) are exact quadrics at machine precision**
	  (1419 cylinder, 1398 cone, 72 sphere; relative residual < 1e-9), and that recovering them takes the
	  assembly to **29/55, not 15/55**. The "bspline basis curve" probe was right about the *storage* and
	  wrong about the *geometry* — a rational bspline basis curve is routinely an exact circle. This
	  milestone is still needed, but for **26/55 genuinely free-form volumes**, not 40.
	- Cheaper reduced-dimension special case worth considering first (covers the 4-5 extrusion/revolution
	  volumes without a general NURBS intersector): a surface swept along a straight line reduces ray/surface
	  intersection to a 2D ray/curve problem in the cross-section plane, and a surface of revolution reduces to
	  a 2D problem in the meridian half-plane. Both reuse the existing 2D bspline `Curve2D` kernel. Low total
	  value (~4 volumes) but a much smaller step than a general trimmed-NURBS intersector, and a useful
	  proving ground for the trim/BVH plumbing.
	- Approach: a `BSplineBoundedSurface` with a genuine **iterative** ray/surface intersection (unlike the
	  closed-form quadrics). Grounded, well-established methods: Newton from a seed (fast, not globally
	  robust alone); **Bézier clipping** — Nishita, Sederberg & Kakimoto, *"Ray tracing trimmed rational
	  surface patches"*, SIGGRAPH 1990 (handles intersection **and** trimming together); or subdivision into a
	  BVH of Bézier sub-patches + Newton refine (fits our architecture best — see the subdivision-BVH
	  performance milestone, which this reuses). The 2D trim test is shared with the analytic surfaces.
	- Reference to validate against: OpenCASCADE's own `IntCurvesFace` / `Geom(2d)API` already ray-cast and
	  trim these — cross-check the intersector against it.
	- Exactness caveat (be honest before committing): an iterative intersector is exact-*to-a-tolerance*. For a
	  bspline *surface* the marginal `Contains` fidelity over a fine adaptive tessellation in the existing BVH
	  is much smaller than it was for quadrics (where exactness was cheap and free). Weigh against simply
	  keeping bspline surfaces on the tessellated fallback. The trim-curve milestone above does **not** carry
	  this caveat — there the surface stays exact and cheap, so the win is unambiguous.

## BVHSurfaceSolid navigation milestones

- [x] Implement `CloseShape` and BVH construction over bounded-surface AABBs.
	- Reuse the BVH2 pattern from `O2Tessellated::BuildBVH`.
	- Store primitive ids and centers exactly as `O2Tessellated` does, but primitive boxes come from surfaces.
	- Expand AABBs by a documented tolerance to avoid missed boundary hits after float conversion.
	- Validation: unit test proves multiple surfaces are traversed and root bounding box is correct.
	- Done 2026-07-23: `CloseShape` now builds a float BVH2 over per-surface `conservativeBounds` AABBs, expanded by the documented `kBVHBoxTolerance` (= 1e-3, must dominate all navigation length tolerances) and then rounded outward during the double->float conversion via `nextafterf`. `config.max_leaf_size = 1` because analytic patch intersections are far more expensive than node box tests and because the bvh2 `traverse_top_down` enters a leaf start node without a box test (a default single-leaf tree would defeat all pruning for small solids). Test hooks `HasBVH()`, `GetBVHRootBounds()` and `CountBVHRayCandidates()` were added; the `BVHConstructionAndTraversal` test checks the root box is conservative but tight (within 2x tolerance), that a through-ray traverses >= 2 leaf surfaces and a miss-ray traverses 0.

- [x] Implement `Contains` using ray parity over surface intersections.
	- Use a skew test direction, as in `O2Tessellated`, but cluster equal or near-equal `t` intersections so an edge or vertex hit does not double-count.
	- Boundary policy: points within tolerance of any surface should be considered inside unless ROOT convention requires otherwise.
	- Acceptance test: box, cylinder, hollow-cylinder, cylinder with added spherical endcaps, and points exactly on faces, edges, and vertices.
	- Done 2026-07-23: `Contains` now uses the BVH: a point-in-box traversal for the boundary policy (leaf boxes are expanded by `kBVHBoxTolerance` >= `kTolerance`, so a point within tolerance of a patch cannot be pruned) followed by a BVH ray traversal collecting all patch intersections along the skew test direction, evaluated by the existing mixed-cluster parity logic (factored into `oddCrossingParity`). A trivial full-loop `Contains_Loop` (named after `O2Tessellated::Contains_Loop`) is kept for debugging/cross-validation, and `Contains` falls back to it before `CloseShape`. Tests: pre-existing box/cylinder/hollow-cylinder/sphere/cone ROOT comparisons now exercise the BVH path; new `BVHConstructionAndTraversal` (two-disjoint-boxes fixture, BVH vs. loop vs. analytic on a grid) and `ContainsBoundaryPointsAndCapsule` (points exactly on faces/edges/vertices, capsule = cylinder + two spherical endcaps with closure check, grid cross-validation and exact capacity).

- [x] Implement `DistFromOutside`.
	- Traverse the BVH with a ray and ask candidate surfaces for intersections.
	- Cull BVH nodes further away than the current best intersection
	- Return the smallest positive entering intersection, respecting `stepmax` and normal orientation.
	- Acceptance test: compare against `TGeoBBox`/`TGeoTube` for outside points and several directions, including grazing misses.
	- Done 2026-07-26, planned and measured in [`SolidNavigationHarness.md`](SolidNavigationHarness.md) (Step 5). Shares `Impl::nearestCrossing<wantEntering>` with `DistFromInside`; `bvh->intersect<false, /*robust*/ true>` with a leaf lambda, `ray.tmax` tightened from inside the lambda, cheap bounding-box reject under `stepmax`, `thread_local` scratch. Public `DistFromOutside_Loop` kept as oracle and baseline: **bit-identical on 57000/57000 rays** over the 19-part DB, and 2.03x slower than the BVH path.

- [x] Implement `DistFromInside`.
	- Return the smallest positive exiting intersection, respecting normal orientation and avoiding immediate re-hit at `t = 0`.
	- Cull BVH nodes further away than the current best intersection
	- Acceptance test: compare against ROOT primitives for central, near-boundary, and oblique rays.
	- Done 2026-07-26 together with `DistFromOutside` (same template, opposite `dot(normal, dir)` sign). `DistFromInside_Loop` bit-identical on 28500/28500 rays; BVH 1.21x faster than the loop.


- [x] Reusable test harness for validation and performance measurement as part of the repo
    - Starting from a surface and tessellated representation obtained from the same CAD solid, create a (compiled or ROOT macro) harness to create test points and directions which can be used to validate Distance, Contains, Safety values and to obtain CPU performance comparisons between the 2 versions
	- This harness can be the foundation of further profile guided optimizations when run within perf
	- We can create a test database from all convertable parts extractable from CAD examples in folder STEP_EXAMPLES... usable on which this harness can run iteratively and repeatedly.
	- Done 2026-07-26. `DetectorsBase/O2SolidHarness.{h,cxx}` (sampling / validation / timing, typed on `TGeoShape*`), the `o2-bench-detectorsbase-solid-harness` front-end (`Detectors/Base/test/runSolidHarness.cxx`), and `scripts/geometry/makeTestPartDB.py` which builds the part DB from `STEP_examples/`. Full design, usage, `perf` entry point and results: [`SolidNavigationHarness.md`](SolidNavigationHarness.md).

- [ ] Code optimization pass
    - Not a real milestone, rather something that be revisited many times, but concrete ideas that are **important** to me:
	    - [x] use the existing BVH traversal functions with lambdas from bvh2 (as in O2Tessellated) instead of writing your own
		- [x] use BVH pruning via tightening the ray properties (not yet done correctly in O2Tessellated) — done in `O2BVHSurfaceSolid` 2026-07-26 and **measured**: removes 26.6% of the analytic patch tests, but only ~1% wall time on the current DB (`--pruning-ab`; numbers and interpretation in `SolidNavigationHarness.md`). `O2Tessellated` still does not do it.
		- [x] do not use std::vector allocations inside navigations functions; Reuse existing members (thread-local) or stay with space on the stack — `Contains`, `Contains_Loop` and both distance queries now use `static thread_local` scratch, one buffer per entry point.
		- [x] try to inline as much as possible
		- [x] think in fast-math terms: Avoid divisions, square roots, etc.
		- [ ] **AOT geometry-specific code generation** (decided 2026-07-31, do after the correctness
		  phases): the geometry is static once converted, so the converter (or a post-pass) can emit
		  specialized C++ per part/patch — constants (radius, frame, neighbor-trim surfaces) folded,
		  no virtual dispatch, branchless sign tests, batchable over rays — compiled once into a
		  shared library shipped with the geometry. This is ahead-of-time specialization in the
		  spirit of VecGeom's specialized navigators, not runtime JIT; if in-process generation is
		  ever wanted, ROOT's cling is already available. Prerequisite: adjacency-trims (regime E2)
		  make the inner loops small closed-form kernels, which is the shape that specializes well.
		  Measure against the harness before and after; template specialization over patch
		  archetypes is the cheaper first step and may capture most of the win.
		- [ ] If iterative surface intersection (Newton for B-spline surfaces) is ever implemented,
		  consider automatic differentiation (e.g. clad, which is clang/cling-based and already in
		  the ROOT ecosystem) for exact derivatives of composed/specialized kernels rather than
		  hand-coded Jacobians.
	- Kept open deliberately: the measurement identified the next two targets, both bigger than anything above. See "measured optimization targets" below.

- [ ] **Measured optimization targets (2026-07-26).** Both come from the first full harness run; see `SolidNavigationHarness.md` for the numbers behind them.
	- `Safety` is the most expensive kernel by a wide margin (10132 ns/call median, 4.74x the tessellated solid) because it is still a plain loop over every patch with no BVH. The priority-queue traversal in the `Safety` milestone below is therefore *the* optimization to do next, not a nicety; `bvh::v2::extra::SafetySqToNode` already exists for it.
	- Per-query cost of the exact solid is ~1.3–2.4x the mesh at the median, dominated by parts with many patches, because every curved patch is one loose conservative AABB. This is the quantified case for the "Subdivision-BVH acceleration for curved surface patches" milestone.

- [ ] **The converter drops BREP edge sharing, and that makes `Contains` genuinely wrong.** Found 2026-07-26 by drilling into the harness's "unexplained" column on `Bagger/BoomCylinderOuter_0_1_1_9`. This is the single most important defect currently known in the exact-surface path; it revises a caveat recorded on 2026-07-24 and it is the reason the two items after it exist.

	**Symptom.** 61 of 10000 sampled points get `Contains = inside` while lying **up to 1.71 cm from the nearest patch** (`Safety` says so). `Contains_Loop` agrees, so it is not the BVH. On a 40^3 grid, 61 of 4453 "inside" answers are >0.5 cm from any patch. `CloseShape` reports 699 boundary edges. The part is not a tube: 8 faces = 5 cylinders + 3 annular planes, four of them carrying general B-spline trim wires.

	**Root cause.** In a BREP the watertightness guarantee is *topological*: two adjacent faces reference the **same** `TopoDS_Edge`, hence the same 3D curve. Each face additionally stores a **pcurve** — that edge's image in *that face's* own (u,v) domain — and the pcurves are independently fitted representations carried with their own tolerance, not derived from one another. `_quadric_trim_wire` extracts `BRep_Tool.CurveOnSurface(edge, face)` (`O2_CADtoTGeo.py:1012`), i.e. the per-face pcurve, and never records which edge it came from. **The shared-edge identity — the thing that made the model watertight in the first place — is discarded at extraction time.** Two independent errors then stack on top of that: the two pcurves of one edge already differ within model tolerance, and the kernel additionally flattens each B-spline pcurve to a polyline *per face* for the winding test. The result is a sliver gap along the seam.

	**Why a sliver gap produces a 1.7 cm error.** Parity containment casts one *fixed* skew ray. Any point whose ray threads the gap loses one crossing and flips to "inside" — so the wrong region is not the size of the gap, it is the gap's **shadow**. Predicted and confirmed experimentally: the wrong region is a narrow tube aligned with the parity test direction, contiguous ~0.5 cm along it and gone within 0.25 cm perpendicular.

	**This supersedes** the 2026-07-24 note's "IsClosed() may warn on the shared-3D-bspline-edge sampling mismatch (documented caveat; navigation exact)". Navigation is **not** exact when this happens, and the framing there (a B-spline sampling artifact) understates it: the defect is the lost topology, and B-splines only make it visible because theirs are the pcurves we flatten.

	**Work items, in the order they should be attacked.** Planned in detail, with acceptance criteria, ordering and environment notes, in [`ExactTrimTopology.md`](ExactTrimTopology.md) — start there.
	1. **Preserve shared edges (converter).** Key trim curves by their `TopoDS_Edge` identity (`TopExp::MapShapesAndAncestors` gives edge -> faces) and hand both adjacent patches the *same* trim curve. Consistency is what parity needs — the curve does not have to lie exactly on either surface, both sides merely have to agree on where the boundary is. This is the fix that addresses the actual cause; everything else is mitigation.
	2. **Analytic point-in-trim, no polyline (kernel).** Replace the flattened-polyline winding test on B-spline wires with an exact 2D ray/curve crossing count: convert each span to Bézier form and root-find (Bézier clipping), which is exact to machine precision rather than sampled. Removes the last sampling step from containment, and removes the cost that dominates these parts.
	3. **Canonical recognition of trim *curves* (converter, cheapest).** ✅ **Done 2026-07-26.** `_recognize_canonical_curve` in `O2_CADtoTGeo.py`, driven from both `_quadric_trim_wire` and `extract_planar_face`; straightness proved from control-polygon collinearity (convex hull), circles at `1e-9` relative. Stored B-spline trim edges **88 → 50** (three-model) and **15034 → 4528** on ALICE3, whose sidecars shrink 6.09 MB → 3.65 MB; **every `unexplained` count bit-identical** on both DBs, which is what an exact recognition must do. **It does not fix the reproducer**: the four curves it recognised there are `BoomCylinderOuter`'s periodic cylinder seams, not the two cylinder-cylinder intersection curves that cause the defect.
	4. **Fail loudly (kernel).** ✅ **Done 2026-07-26.** `O2BVHSurfaceSolid::NavigationReliability` + `IsNavigable()` + boundary/non-manifold/reversed edge counts; the three `CloseShape` closure defects are now `Error`s whose text states the consequence; the harness prints the state per part, carries it in `--json`, and lists every unnavigable part at the end of a run. Test `NavigationReliabilityIsQueryable`. Three-model DB: 11 reliable, 6 open-surface-set, 2 non-manifold.

	**ROOT CAUSE FOUND 2026-07-26: this was a kernel bug, not lost topology — items 1 and 2 are not the fix and should not be started.** Enumerating the parity ray's crossings at a wrong point shows a **doubled** crossing, not a missing one: two consecutive `ENTER`s ~0.04 cm apart, the second lying `0.026 cm` inside the other face's hole, which must therefore exclude it. Reconstructing that hole curve from the sidecar independently (de Boor, 179 poles) shows it is closed to `5.3e-12`, has the intended extent, and *does* contain the point — **the converter data was right and the kernel disagreed with it**. `bsplineSampleRecursive` ended its recursion whenever the chord `p0 -> p1` fell below the flatness scale; a **closed** curve has `p0 == p1` exactly, so a full circle flattened to two coincident points and every polyline-based query (winding, closest point, boundary band, display mesh, closure check) saw an empty curve. It hid because `signedAreaContribution` integrates by Gauss-Legendre, so the wire still validated and still reported the correct area. Fixed; regression test `BSplineHoleInCylinderWall`. Result on the three-model DB: `contains` unexplained 4588 → 4430 (525 → **367** excluding the unchanged non-manifold `oTOF` part), `distin` 218 → **114**, the reproducer 51 → **16**, and every Bagger cylinder's `distin` → **0**; two parts whose `Contains` disagreed with `Contains_Loop` now agree everywhere; ALICE3 unchanged (no regression). *An earlier version of this note blamed near-tangency amplification — that was inferred rather than measured and is wrong: the two cylinders meet at 59-60°.*

- [ ] **Patch *cost*, not patch *count*, is what makes small parts slow.** Same part, measured 2026-07-26: 8 analytic patches against 2244 triangles, yet `Contains` is 1.25x slower than the mesh (3102 vs 2490 ns) and `DistFromInside` 1.58x. The BVH cannot help — with 8 patches it is only 1.18-1.20x over `_Loop`. The cost is per-patch: each query solves the quadric and then runs a winding test against a flattened B-spline polyline. Items 2 and 3 above are the lever, not the acceleration structure.
	- Counterpoint worth keeping: on this same part `Safety` is **1.9x faster** than the tessellated solid (4557 vs 8705 ns), because 8 patch-distance evaluations beat a priority-queue walk over 2244 triangles. The "Safety is the worst kernel" finding from the DB-wide medians is driven by parts with many patches; on low-patch parts the exact representation already wins.

- [ ] **`Contains` disagrees with `Contains_Loop` on non-manifold parts.** 301/142500 points over the 19-part DB, 295 of them in the two `oTOF` parts that `CloseShape` reports as non-manifold (32 and 9 non-manifold edges = coincident/duplicated faces). Parity containment is order-dependent on such input — `oddCrossingParity` clusters near-equal hits and the BVH supplies them in a different order than the loop — so this is the containment algorithm meeting degenerate geometry, not a traversal bug (the distance twins are bit-identical on the same DB). Decide whether to make the clustering order-independent or to reject non-manifold input at `CloseShape`.


- [ ] Implement `Safety` using priority-queue BVH traversal.
	- Start from the `O2Tessellated::SafetyKernel` pattern.
	- Node ordering uses distance to node AABB; leaves use surface distance or a conservative lower bound.
	- For unsupported exact distance on a surface, return a conservative estimate and record this in the surface support matrix.
	- Acceptance test: safety is never larger than the true distance for representative fixtures.

- [ ] Implement `ComputeNormal`.
	- Preferred source: nearest surface by safety query with face id return, then evaluate the exact oriented surface normal.
	- Direction convention should match ROOT expectations and `O2Tessellated` behavior.
	- Acceptance test: normals on box faces, cylinder side, caps, and edges.

- [ ] Implement `Capacity`.
	- Preferred exact path: integrate oriented surface patches using divergence theorem where each supported surface can provide its contribution.
	- Interim acceptable path: return visualization-mesh volume with a clear flag/comment that capacity is approximate until exact per-surface contributions are implemented.
	- Acceptance test: box and cylinder capacities within defined tolerance.

- [ ] Implement visualization mesh generation independent of navigation.
	- Each surface provides a display triangulation with configurable or fixed display quality.
	- Navigation must never depend on display triangles.
	- Acceptance test: `GetBuffer3D` returns nonzero vertices/segments/polygons and ROOT can draw a simple solid.

- [ ] **Subdivision-BVH acceleration for curved surface patches** (performance; revisit after correctness).
	- Motivation: today each surface is one BVH leaf with one *conservative* AABB. For curved patches this box
	  is loose — a trimmed-arc cylinder/sphere gets the full-cylinder / full-sphere envelope (observed in the
	  ALICE3 test: `ComputeBBox` over-covers the ~117° channel arc). Loose leaves inflate ray-candidate counts
	  and every candidate is an expensive analytic patch test.
	- Idea (wanted even for the *existing* analytic planes/quadrics, not only for bspline): subdivide each
	  curved patch in its parametric domain into sub-patches, bound each sub-patch tightly, and store the
	  sub-patches as BVH leaves referencing the same analytic surface plus a parametric sub-window. Navigation
	  stays **exact** — only the acceleration granularity changes.
	- Build it generically over `BoundedSurface` (a leaf = surface pointer + parametric sub-window) so every
	  curved surface opts in with one shared mechanism, and so the eventual `BSplineBoundedSurface` reuses the
	  exact same subdivision structure its intersector needs.
	- Acceptance test: ray-candidate counts drop materially on a trimmed-arc quadric fixture vs. the current
	  single-box leaf, with `Contains`/distance results bit-for-bit unchanged.

## Surface sidecar format

The exact-surface sidecar (`surfaces_<VOLNAME>_<LID>.bin`), written by
`write_surfaces_bin` in `scripts/geometry/O2_CADtoTGeo.py` and read by
`o2::base::LoadSurfaceSolid` (`Detectors/Base/include/DetectorsBase/O2SurfaceSolidIO.h`).
All integers are little-endian `uint32`, all geometry values are little-endian `float64`,
lengths in cm, angles in radians.

**Version 2** appends one `float64` to the fixed header and changes nothing else, so a version-1
file is a version-2 file that does not state its model's tolerance. The converter writes version 2;
the reader accepts both, and warns when it falls back for a version-1 file.

```
header:
  char[4]  magic        = "O2SS"
  uint32   version      = 2
  uint32   nSurfaces
  uint32   reserved     = 0
  float64  modelTolerance          # version 2 only; cm; 0 = not stated
per surface (nSurfaces times):
  uint32   surfaceType     1=plane 2=cylinder 3=cone 4=sphere 5=torus
  uint32   flags           bit 0: innerWall (quadrics/torus: outward normal towards the axis/center/tube spine)
  uint32   nParams
  float64  params[nParams]  fixed per-type layout, see below
  uint32   nWires
  per wire (nWires times):
    uint32   wireRole      0=outer 1=inner (hole)
    uint32   nEdges
    per edge (nEdges times):
      uint32   curveType     0=line segment 1=circular arc 2=bspline (2D, clamped)
      uint32   nCurveParams
      float64  curveParams[nCurveParams]
        line:    u0 v0 u1 v1
        arc:     cu cv radius phiStart phiSweep (signed sweep, full circle = +/-2*pi)
        bspline: degree nPoles poles[2*nPoles] weights[nPoles] knots[nPoles+degree+1]
                 (clamped flat knot vector; weights all 1 => non-rational; Bezier is the
                  single-span special case)
```

Per-type `params` layouts (matching the `O2BVHSurfaceSolid::Add*Surface` signatures):

- plane (9): origin xyz, axisU xyz, axisV xyz. The trim is carried in the wire block as
  general line/arc/bspline loops (polygon, disk/annulus, rounded rectangle, slot,
  spline-bounded plate, ...); exactly one outer wire, holes as inner wires. A pure line-segment
  loop loads through `AddPlanarSurface` (polygon); any arc or bspline edge loads through
  `AddCurvedPlanarSurface` (`CurvedPlanarBoundedSurface`, orthonormal frame). A disk is one
  full-circle outer wire; an annulus adds an inner circle.
- cylinder (14): centerPoint xyz, axis xyz, referenceAxisU xyz, radius, heightMin,
  heightMax, phiStart, phiSweep. An empty wire block means the trim is the scalar parametric
  rectangle. A non-empty wire block carries a general line/arc/bspline trim in the periodic
  `(u, v) = (phi[rad], h[cm])` domain (one outer wire, holes as inner wires); it is
  authoritative and the scalar params then act as the frame + a conservative window.
- cone (15): centerPoint xyz, axis xyz, referenceAxisU xyz, radiusAtMin, radiusAtMax,
  heightMin, heightMax, phiStart, phiSweep. Wire block as for the cylinder (domain
  `(phi[rad], h[cm])`); the scalar radii pin the linear radius law `r(h)`.
- sphere (14): center xyz, polarAxis xyz, referenceAxisU xyz, radius, thetaMin, thetaMax,
  phiStart, phiSweep. Wire block as for the cylinder, in the `(u, v) = (phi[rad], theta[rad])`
  domain.
- torus (15): centerPoint xyz, axis xyz, referenceAxisU xyz, majorRadius, minorRadius,
  phiStart, phiSweep, tubeStart, tubeSweep. `majorRadius` is the axis-to-tube-centre radius and
  `minorRadius` the tube radius. The two periodic angles are the ring angle
  `phiRing` (around the axis, from `referenceAxisU`) and the tube angle `phiTube` (around the
  tube, measured from the outer equator towards the +axis pole). An empty wire block means the
  trim is the scalar parametric rectangle `phiRing x phiTube`; a non-empty wire block carries a
  general line/arc/bspline trim in the `(u, v) = (phiRing[rad], phiTube[rad])` domain (both
  periodic). A full torus (both sweeps `2pi`) is self-closing.

The `nParams`/`nCurveParams` counts make every record self-describing, so future minor
extensions can add trailing parameters or new curve types without breaking version-1
readers that choose to skip them; incompatible changes bump `version`. A quadric trim wire
must be a non-wrapping loop within one phi window (`<= 2pi`; for the torus, also within one turn
in the tube angle); the writer emits it only when the face is not the plain parametric
rectangle. Line edges stay lines; every curved pcurve
(circle/ellipse/Bezier/bspline) is converted to a bspline and its poles are pushed through the
*affine* `(u, v) -> (phi, height/theta)` map (a bspline is closed under affine maps, so this is
exact - an anisotropic map merely turns a circle into an exactly-represented ellipse-bspline).
Consecutive edge endpoints may differ by the extractor's ~1e-6 projection precision; the reader
joins within `kJoinTolerance = 1e-5` and the kernel's `CurveWire` accepts the loop within the
matching `kWireJoinTolerance = 1e-5` (both looser than the `1e-9` boundary tolerance).

## CAD conversion milestones

- [x] Add exact-surface extraction probes to `O2_CADtoTGeo.py` without changing default output.
	- Use OpenCascade face/surface adaptors to classify faces by analytic surface type.
	- For each simple shape, collect: face type, orientation, surface parameters, boundary wires/edges, and whether all faces are supported.
	- Deliverable: a `--surface-report` or debug JSON path that reports exact-conversion eligibility per logical volume.
	- Validation: run on one STEP example and confirm existing `geom.C`/facet output is unchanged.
	- Done 2026-07-23: `--surface-report <path>` writes a JSON report built from a new `def_shapes` store (the possibly clipped leaf shape). `classify_face` (`BRepAdaptor_Surface` type + per-type parameters in cm + `TopAbs_REVERSED` orientation + UV bounds), wire/edge probing (`BRepAdaptor_Curve` types, degenerated-edge handling) and `face_supported` encode the current C++ support matrix: planes need line/circle boundaries, quadrics need a single parametric-rectangle trim (detected by sampling each edge's pcurve for iso-parametric u/v). Validated on `oTOF_V3_R92cm.step`: output folders byte-identical with/without the flag; report finds 1607 planar faces / 9116 line edges, 2/3 volumes eligible (the third is a faceless leaf).

- [x] Design and document the generated surface sidecar format.
	- It should be versioned and independent of Python object pickling.
	- Candidate: compact binary sidecar similar to `facets_*.bin`, with magic/version, surface count, type enum, orientation, parameters, wire count, edge count, edge curve records, and optional display mesh.
	- Deliverable: format description in this file and loader skeleton in generated C++ macro.
	- Validation: round-trip one planar box sidecar into an `O2BVHSurfaceSolid`.
	- Done 2026-07-23: format documented in the "Surface sidecar format" section above; Python writer `write_surfaces_bin` added to `O2_CADtoTGeo.py`. Deviation from the original deliverable: the loader is a shared library function `o2::base::LoadSurfaceSolid` (`Detectors/Base/include/DetectorsBase/O2SurfaceSolidIO.h` + `src/O2SurfaceSolidIO.cxx`) rather than code pasted into the generated macro, because the macro needs libO2DetectorsBase anyway to instantiate the solid and the unit test can then exercise the exact same parser; `emit_cpp_prelude` will only gain a thin call in the "generated C++ macro support" milestone. Validated by the `SurfaceSidecarRoundTrip` unit test (C++-written box sidecar vs `TGeoBBox`, cylinder+caps sidecar vs `TGeoTube`, bad-magic/truncation rejection) plus a cross-language parity run: Python-written box and tube sidecars loaded via `LoadSurfaceSolid` in a ROOT macro with 0 `Contains` mismatches and exact capacities.

- [x] Add generated C++ macro support for `O2BVHSurfaceSolid`.
	- Extend `emit_cpp_prelude` with a `LoadSurfaceSolid(...)` helper and include the new header.
	- Add `emit_surface_solid_cpp(...)` next to `emit_tessellated_cpp(...)`.
	- Keep `emit_tessellated_cpp` as the fallback path.
	- Validation: generated macro compiles and builds a ROOT geometry with exact surface solids for supported volumes.
	- Done 2026-07-23: `emit_cpp_prelude(exact_surfaces=False)` optionally appends an exact-surface block (`R__ADD_INCLUDE_PATH($O2_ROOT/include)`, `R__LOAD_LIBRARY(libO2DetectorsBase)`, textual include of the dictionary header `O2BVHSurfaceSolid.h`, a manual `LoadSurfaceSolid` prototype since `O2SurfaceSolidIO.h` is not a dictionary header, and a `LoadSurfaces(file, solid, check)` helper that loads + `CloseShape(check)` + throws on closure/orientation failure when checking). New `emit_surface_solid_cpp(lid, name, sidecar_path, medium)` mirrors `emit_tessellated_cpp`. `emit_root_macro` gained an optional `surface_files: Dict[lid, sidecar_path]` argument: listed volumes become `O2BVHSurfaceSolid`, everything else keeps the tessellated fallback; the default (`None`) leaves output byte-identical (re-verified on `oTOF_V3_R92cm.step`). No CLI wiring yet — that is the mode-flags milestone. Validated end-to-end on a generated centered-box STEP: pipeline run with `surface_files` produced a `geom.C` that loads in ROOT interpreted mode, `build_and_export` closes the geometry with 0 overlaps and exports, top shape is `o2::base::O2BVHSurfaceSolid` with exact capacity 48 and 0 `Contains` mismatches vs `TGeoBBox` on a 9^3 grid.

- [x] Add converter mode flags.
	- Proposed flags: `--exact-surfaces off|auto|required` and `--surface-report <path>`.
	- `off`: current behavior.
	- `auto`: exact where every face of a logical volume is supported; fallback to tessellated otherwise.
	- `required`: fail with a useful report when any selected logical volume cannot be represented exactly.
	- Validation: command-line help and one smoke run per mode.
	- Done 2026-07-23: `--exact-surfaces off|auto|required` (default `off`) added to `main()` and threaded into `emit_root_macro(exact_surfaces=...)`. In `auto`/`required` the converter runs `extract_surfaces_for_shape` over every leaf `def_shapes` entry after `extract_graph`; each fully-extracted solid gets a `surfaces_<VOLNAME>_<LID>.bin` sidecar and is emitted as `O2BVHSurfaceSolid`, otherwise it keeps the tessellated fallback (`auto`) or aborts with a per-volume/per-face reason report (`required`). `off` leaves output byte-identical to the pre-flag default (facets identical, no sidecars). The explicit `surface_files` test hook is now only honoured when `exact_surfaces == "off"`. Validated: `--help` shows the flag; box STEP smoke run per mode (off/auto/required), cylinder STEP in `auto` falls back (circular caps -> "not a line segment", wall -> "not implemented yet") and in `required` fails with the report.

- [x] Implement planar face extraction from OpenCascade.
	- Convert OpenCascade face wires and edges into the planar surface representation.
	- Reject unsupported boundary curves until curve classes are implemented.
	- Validation: generated exact box or box-like STEP uses `O2BVHSurfaceSolid`; unsupported faces still fallback in `auto` mode.
	- Done 2026-07-23: `extract_planar_face` builds a sidecar `plane` record (origin/axisU/axisV in cm, line-segment wires) from a `GeomAbs_Plane` face. It walks each wire with `BRepTools_WireExplorer` in connected order and projects the ordered edge-start vertices into the plane's local `(u, v)` frame. Outward normal is enforced by choosing the sign of `axisV` so `axisU x axisV` equals the face's outward normal (surface normal for a FORWARD face, its opposite for REVERSED; also robust to a left-handed OCC `ax3`). Non-line boundary edges (circles/arcs), degenerated edges, and non-plane faces are rejected (version-1 planar sidecar wires are line-only) so those solids fall back in `auto`. Dispatch table `_FACE_EXTRACTORS` currently holds only `plane`; the quadric extractors slot in at the next milestone. Validated: the generated box sidecar loaded via `LoadSurfaceSolid` closes with `IsClosed()`/`IsOrientationConsistent()` true, `Capacity` = 0.192 exactly, and 0 `Contains` mismatches vs `TGeoBBox` over 2e5 random points; cylinder faces still fall back.

- [x] Implement cylindrical, spherical, and conical face extraction.
	- Reuse the support matrix from the C++ surface milestones.
	- Preserve face orientation from OpenCascade so normals point outward.
	- Validation: generated cylinder/sphere/cone fixtures navigate like their ROOT primitive equivalents.
	- Done 2026-07-23: added `extract_cylindrical_face`, `extract_conical_face`, `extract_spherical_face` (dispatched from `_FACE_EXTRACTORS`), each emitting the quadric sidecar record consumed by `AddCylindricalSurface`/`AddConicalSurface`/`AddSphericalSurface`. Also generalized `extract_planar_face` to route circle-bounded planar faces to a new `_extract_planar_disk` (sidecar `disk`/annulus record -> `AddPlanarDiskSurface`), which is what actually closes cylinder/cone end caps — the milestone is only "useful" once the caps convert too. OCC->C++ mapping: U is always azimuth, V is height (cylinder), ruling-line distance (cone; `h=v·cos α`, `r=RefRadius+v·sin α`, endpoints reordered so `heightMin<heightMax`), or latitude (sphere; `θ=π/2−v`). `_quadric_phi_range` maps the OCC U-range to `(phiStart, phiSweep)` and mirrors it (`phiStart=−umax`) when the `gp_Ax3` is left-handed (`ax3.Direct()` false), so partial sweeps stay correct. `inner_wall = (face.Orientation()==REVERSED)` — the OCC canonical quadric normal points outward radially, so a REVERSED face is a hole wall. Quadric trims must be parametric rectangles (`_quadric_trim_is_rectangular` reuses `_edge_pcurve_is_iso`); non-iso trims fall back. Python validates ranges (positive height/sweep/theta, non-negative radii) so a record is only emitted when the C++ `initialize` will accept it. Validated: extracted sidecars for a solid cylinder (r0.3, h0.5), truncated cone (r0.4->0.2, h0.6), full sphere (r0.3) and a hollow tube (R0.4/r0.2, inner cylinder correctly flagged `inner_wall`) each `CloseShape` to `IsClosed()`/`IsOrientationConsistent()` true with **0** `Contains` mismatches vs analytic references over 4e5 random points; full converter CLI on `cyl.step --exact-surfaces required` reports "1/1 leaf solids represented exactly" and emits `O2BVHSurfaceSolid`; a torus still falls back in `auto` and aborts with a report in `required`.

	Session handoff note (2026-07-23, quadric + disk extraction): the exact-surface CAD path now covers all four analytic surface families the C++ kernel supports, so simple revolved/prismatic CAD solids (box, cylinder, cone, tube, sphere, and their combinations) convert end-to-end to exact `O2BVHSurfaceSolid` with exact navigation. Design decision worth remembering: **disk/annulus caps are extracted from *planar* faces** (circle-bounded), not from the quadric extractors — `extract_planar_face` now branches on boundary-curve kind (`{line}`->polygon, `{circle}`->disk, mixed/other->fallback), so a face that mixes straight and curved boundary edges (e.g. a D-shaped or rounded-rectangle plate) still falls back. Known gaps for the next session: (1) quadric faces are limited to *parametric-rectangle* trims — a cylinder whose end is cut by a slanted plane has a non-iso boundary and falls back; general curved-boundary trims on quadrics need the Curve2D/CurveWire sidecar wiring that already exists on the C++ side but is not emitted by the writer yet; (2) planar faces bounded by a mix of lines and arcs (rounded polygons) are still unsupported for the same reason; (3) `inner_wall` is derived purely from `face.Orientation()`, which assumes the OCC canonical normal is outward-radial — solid for `BRepPrimAPI`/STEP solids but worth a geometric double-check if a vendor STEP ever trips `IsOrientationConsistent()`. Next milestone in this file: "Keep clipping and name filtering compatible with exact output."

- [x] Support general line+arc planar bounded surfaces (retire the disk special-case).
	- Planar CAD faces whose boundary mixes straight segments and arcs (rounded rectangle, slot, D-shape) must convert exactly, not just polygons and circular disks. OpenCascade produces these routinely.
	- The kernel already models this (`Curve2D`/`CurveWire`/`CurvedPlanarBoundedSurface`); wire it through the public API, the sidecar reader and the Python extractor.
	- Validation: rounded-rectangle round-trip closes and navigates exactly; the quadric fixtures (caps now via the general path) still pass.
	- Done 2026-07-23: the gap was purely plumbing — the kernel and the sidecar *format* already carried general line/arc wires. Added a public POD `O2BVHSurfaceSolid::PlanarBoundaryCurve` (line/arc) + `AddCurvedPlanarSurface(origin, axisU, axisV, outerWire, innerWires)` translating to internal `Curve2D` and `CurvedPlanarBoundedSurface`; **removed** `AddPlanarDiskSurface` and the sidecar `kPlanarDisk` type (surfaceType 5). The reader's `kPlane` case now builds line/arc wires for every plane: a pure line loop keeps the `AddPlanarSurface` polygon path (byte-identical output, general-metric frame), any arc routes to `AddCurvedPlanarSurface`. Python `extract_planar_face` is now one general extractor emitting a mix of `line`/`arc` edges per wire; arcs recover the signed sweep by sampling the 3D edge in wire-traversal order (backwards when the edge is `REVERSED`) and unwrapping the projected polar angle — robust to full circles (single periodic edge -> +/-2pi), >pi arcs, and either winding. A disk is now just one full-circle arc wire; an annulus adds an inner circle. Unit tests: cylinder/cone/capsule cap fixtures reworked onto `AddCurvedPlanarSurface`, the `SurfaceSidecarRoundTrip` tube caps rewritten as arc-wire plane records, and a new `CurvedPlanarStadiumPrism` test (stadium caps = two lines + two semicircle arcs, extruded with flat walls + half-cylinders) checks `IsClosed`/`IsOrientationConsistent`, exact capacity, and `Contains` vs an analytic stadium-prism predicate on a grid. Python end-to-end: a rounded-rectangle plate STEP converts with `--exact-surfaces required` (no fallback) and matches an analytic reference with 0 `Contains` mismatches; cylinder/cone/sphere/tube unchanged.

	Session handoff note (2026-07-23, general planar wires): planar exact surfaces are now fully general (any line+arc boundary), which was the level a real detector needs. The `disk` special-case is gone from every layer (kernel API, sidecar type, extractor). Remaining known gaps: (1) **general curved trims on quadrics** — a cylinder/cone/sphere cut by a slanted or curved boundary still falls back, because `CylindricalBoundedSurface`/`ConicalBoundedSurface`/`SphericalBoundedSurface` trim by a parametric rectangle (scalar height/phi/theta bounds), not a `CurveWire` in their `(phi, h)`/`(phi, theta)` domain; closing this needs kernel work (a curve-wire trim on quadrics) plus the sidecar quadric wire block (reserved but unused) and matching extractor — the natural next surface milestone; (2) `inner_wall` on quadrics is still derived from `face.Orientation()` only (see prior note). Next milestone in this file: "Keep clipping and name filtering compatible with exact output."

- [x] Support general curved (line/arc `CurveWire`) trims on quadrics.
	- Give `CylindricalBoundedSurface`/`ConicalBoundedSurface`/`SphericalBoundedSurface` an optional general trim in their periodic `(u, v)` domain (`(phi, h)` for cylinder/cone, `(phi, theta)` for sphere) instead of only a scalar parametric rectangle, so quadric faces with holes/windows or non-rectangular line-bounded parametric regions convert exactly. Wire it through the public API, the (already-reserved) sidecar quadric wire block, and the Python extractor.
	- Scope: exactly representable trims are those whose parametric boundary edges are lines or arcs in `(u, v)`. A genuinely slanted/curved 3D planar cut maps to a sinusoid/small-circle in `(u, v)` and stays on the tessellated fallback — this is a fundamental `Curve2D` limitation, not a bug.
	- Done 2026-07-24: kernel — added an optional `CurveWire mTrimOuter` + inner holes to each quadric, a wire-taking `initialize(...)` overload, a `pointInTrim(u, v)` classifier (phi unwrapped into the wire window), and branched `containsPointOnSurface`/`appendIntersections`/`capacityContribution`/`appendDisplayMesh`/`appendDirectedEdges` on `mHasWireTrim`. Shared free helpers in `BoundedSurface.h`: `curveTrimContains`, `integrateOverCurveTrim` (numeric midpoint quadrature of `(1/3) X·n |X_u x X_v|` — wire-trimmed capacity is therefore *inexact*, flagged via `capacityIsExact()==false`; the untrimmed rectangle keeps its exact closed form), `buildCurveTrim` (validates + rejects a trim spanning `> 2pi` in phi, tightens the conservative window to the wire bounds), and `sampleCurveWireByU`/`appendCurveTrimMesh`/`appendCurveTrimEdges`. Key subtlety: a straight `(phi, h)` line maps to a *curved* 3D edge (a `h=const` rim is a circle), so `sampleCurveWireByU` chords each parametric edge by its phi-span at `kArcSamples`/turn, making a rim match adjacent chord-sampled caps in the closure half-edge check (a phi-const seam stays one segment). `distanceSqToPatch` returns the rectangle-window value as a documented conservative lower bound when trimmed (keeps `Safety` safe). Sphere's `(phi, theta)` map is orientation-reversed vs the outward normal, so its directed-edge sign is `-mNormalSign` (cylinder/cone: `+mNormalSign`). Public API: wire-taking overloads of `AddCylindricalSurface`/`AddConicalSurface`/`AddSphericalSurface` (reusing the `PlanarBoundaryCurve` POD for the `(u, v)` domain). Reader: `collectQuadricTrim` feeds a non-empty quadric wire block through the new overloads; an empty block keeps the scalar path (byte-identical). Python: the three quadric extractors keep the scalar path only when the trim *fills* its UV box (`_quadric_trim_fills_uv_box` — a shoelace-area check, strictly more correct than the old every-edge-iso test, which mis-classified L-shaped/notched all-iso faces as full rectangles and would have silently filled the notch); otherwise `_quadric_line_wire` emits a line-only wire block (arc pcurves fall back — they would become ellipses under the non-uniform remap). Tests (`testBVHSurfaceSolid.cxx`, now 23 cases): `WireTrimmedCylinderMatchesTube`/`WireTrimmedConeMatchesCone` (rectangle-as-wire equivalence vs `TGeoTube`/`TGeoCone`, incl. closure + capacity), `WireTrimmedQuadricKernels` (cylinder window/hole exclusion + ray filtering, an arc-in-`(phi,h)` trim, a sphere-section-as-wire vs the scalar section, and a `>2pi` rejection), `WireTrimmedSidecarRoundTrip` (line wire block through the reader vs `TGeoTube`). Python end-to-end: a plain cylinder keeps the scalar path (no wire block); a box-notched cylinder wall emits a 12-line wire block that round-trips through `LoadSurfaceSolid`.

	Session handoff note (2026-07-24, quadric curve-wire trims): quadric exact surfaces now accept general line/arc `CurveWire` trims (holes/windows, notched/polygonal parametric regions). Navigation (Contains/DistFrom*/normals) is exact; **capacity of a wire-trimmed quadric is numerically integrated (inexact, flagged)** and its `distanceSqToPatch`/`Safety` is a conservative lower bound — tighten later if needed (an exact Green's-theorem line integral is derivable per surface but was deprioritized). Discovered but **out of scope** (pre-existing, in the planar path, not touched here): a planar face whose boundary mixes a line and a partial arc (e.g. a D-shaped cap from a box-cut cylinder) can fail the sidecar reader's strict `1e-9` edge-join tolerance (`wireToCurves` in `O2SurfaceSolidIO.cxx`) because `extract_planar_face`'s arc endpoint (`center + r*cos(a0)`) and the neighbouring line's projected vertex differ by ~`1e-6`; the kernel's own `CurveWire` closure tolerance is looser, so relaxing `kJoinTolerance` in the reader is the likely fix. This blocks *whole* box-cut solids that contain such caps from converting, though the quadric wall face itself round-trips fine. Remaining quadric gaps: (1) arc pcurve trims on quadrics (need an ellipse/parametric-conic curve type or a similarity-only fast path); (2) exact wire-trim capacity; (3) `inner_wall` still from `face.Orientation()` only (see prior notes). Next milestone in this file: "Keep clipping and name filtering compatible with exact output."

- [ ] Keep clipping and name filtering compatible with exact output.
	- Current clipping can create new faces through boolean intersection; after clipping, rerun eligibility detection on the clipped shape.
	- If clipping creates unsupported curves/surfaces, fall back to tessellated output in `auto` mode.
	- Validation: existing clipped STEP example commands still produce output, with exact surfaces only where eligible.

## Test and validation milestones

- [x] Add a focused C++ unit test target for `O2BVHSurfaceSolid`.
	- File candidate: `Detectors/Base/test/testBVHSurfaceSolid.cxx`.
	- Register in `Detectors/Base/CMakeLists.txt` with label `detectorsbase`.
	- Test categories: primitive surface kernels, solid closure/orientation, TGeo navigation, visualization mesh, streaming if practical.
	- Done 2026-07-22: `BVHSurfaceSolid` Boost test registered under `DetectorsBase`; current coverage builds a six-planar-face box and compares `Contains`, distances, safety, normal, capacity, and mesh counts with `TGeoBBox`.

- [ ] Add randomized comparison tests against ROOT primitives for supported exact shapes.
	- For box, tube, sphere, and cone fixtures, sample points and directions with a fixed seed.
	- Compare `Contains`, `DistFromInside`, `DistFromOutside`, `Safety`, `ComputeNormal`, and `Capacity` within declared tolerances.

- [ ] Add regression tests for hard numerical cases.
	- Tangent rays, points on a boundary, rays through shared edges/vertices, very small faces, large coordinate offsets, repeated or nearly repeated intersections.

- [ ] Add converter smoke tests if pythonOCC is available in the test environment.
	- Keep tests optional or environment-gated if pythonOCC is not guaranteed in CI.
	- At minimum, test the pure Python eligibility/report functions without requiring full STEP import.

- [ ] Add performance benchmarks or timing checks outside regular unit tests.
	- Compare exact-surface BVH navigation against `O2Tessellated` for representative CAD parts.
	- Track build time, memory, `Contains`, distance, and safety throughput.

## Documentation and rollout milestones

- [ ] Document the exact-surface support matrix.
	- Include C++ support, converter support, distance exactness, capacity exactness, and fallback behavior.

- [ ] Document generated-file behavior in `O2_CADtoTGeo.py` help text.
	- Explain when `facets_*.bin` is emitted, when surface sidecars are emitted, and how `auto|required` modes behave.

- [ ] Add one small exact-surface example under `scripts/geometry/STEP_examples/` or document how to generate one.
	- Keep generated binary artifacts out of the repository unless existing conventions say otherwise.

- [ ] Add a final migration note once exact surfaces are stable.
	- Explain how users can switch existing `O2_CADtoTGeo.py` workflows from tessellated-only to exact-surface auto mode.

## Open questions and ambiguities

- What should the final public class name be: `BVHSurfaceSolid`, `O2BVHSurfaceSolid`, or something else matching ROOT naming plans?
- Should bounded-surface classes be public API in `include/DetectorsBase`, or private implementation details behind the final shape?
- Is exact `Capacity()` required for the first usable version, or is a documented visualization-mesh approximation acceptable initially?
- Which exact surface types are mandatory for the first milestone: planar only, planar plus cylinders, or the full planar/cylinder/cone/sphere set?
- What is the expected persistence format for generated exact CAD data: binary sidecar, JSON sidecar, direct generated C++ construction, or ROOT streaming only?
	- Decision 2026-07-23: ROOT persistence of the solid is deliberately deprioritized for now. The current implementation keeps all surface data in a transient private `Impl` (`fImpl //!`), so a streamed-then-read solid comes back empty; this is acceptable while the shape is under active development. Add a streamable persistent representation (a POD mirror of the surfaces streamed like `O2Tessellated`'s `fVertices`/`fFacets`, rebuilt via `CloseShape(false)` on read) once the surface/navigation functionality is complete.
- Should `Contains` classify boundary points as inside, outside, or follow ROOT primitive behavior exactly case by case?
- Are OpenCascade circle/ellipse/BSpline trim curves required for the first converter milestone, or can unsupported trims force tessellated fallback?
- Should the internal wire model mirror OpenCascade orientation exactly, or normalize to an O2 convention at import time and keep the original orientation only for diagnostics?
- Should this implementation aim for eventual upstreaming to ROOT, and if yes, should naming and dependencies avoid O2-specific choices from the start?

## Session handoff notes

- 2026-07-22: Plan created from the initial project note, current `O2Tessellated` implementation, `DetectorsBase` CMake integration, and the current `O2_CADtoTGeo.py` tessellated emission path. No implementation has started yet.
- 2026-07-22: Added first `O2BVHSurfaceSolid` implementation with private line-trimmed planar surfaces, display triangulation, loop-based TGeo navigation methods, ROOT dictionary wiring, and `Detectors/Base/test/testBVHSurfaceSolid.cxx`. Validated by object-compiling the new implementation and test with `/data/swenzel/sw/BUILD/O2-latest/O2/compile_commands.json`, then manually linking/running the focused Boost test against a temporary ROOT dictionary: `*** No errors detected`. Full `O2lib-DetectorsBase` build currently stops during CMake regeneration because `O2ReportNonTestedMacros` reports existing `scripts/geometry/STEP_examples/oTOF_TGeo_clip_xminus/geom.C` as untested.
- 2026-07-23: Confirmed a clean full build via `./alibuild/aliBuild build O2 --debug`. Reviewed the initial two milestones and committed the code, build wiring, and this plan on branch `swenzel/bvhsurfacesolid` (`g4Config.C` and the generated `scripts/geometry/STEP_examples/` artifacts intentionally left out of the commit). Decision recorded: ROOT persistence of the solid is deprioritized until the surface/navigation functionality is complete (see Open questions).
- 2026-07-23: Completed the private bounded-surface/wire/edge interface milestone and the planar wire orientation & closure validation milestone. Added `Detectors/Base/src/BoundedSurface.h` (abstract `BoundedSurface`, `SurfaceEdge`, role/status-aware `SurfaceWire`, `PlanarBoundedSurface`, `DummyBoundedSurface`, and `validateClosure`), refactored `O2BVHSurfaceSolid` to navigate polymorphically over `std::unique_ptr<BoundedSurface>` with entering/exiting decided by the per-hit oriented normal, added `IsClosed()` / `IsOrientationConsistent()`, and extended the unit test with dummy-surface, wire-validation, and missing/reversed-face fixtures. Both source and test object-compile with the configured O2 flags.
- 2026-07-23: Completed the remaining two data-model milestones (numerical conventions tests and the initial wire/edge data model). Added reusable surface-independent primitives (`SurfaceEdge::closestPoint`, `SurfaceWire::parametricBounds`, `SurfaceWire::sampledBoundary`) in `BoundedSurface.h` and new `NumericalConventions` and `WireDataModel` test cases in `testBVHSurfaceSolid.cxx` (near-boundary/near-tangent, square-with-hole, reversed/open wires, point-on-edge, edge projection). Full build and focused test confirmed passing by the user. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-23: Completed the "curve classes for trimmed surface boundaries" milestone. Added a surface-independent `Curve2D` (line + circular arc) value type and a `CurveWire` closed-loop type to `BoundedSurface.h` with exact endpoint/tangent/bounding-box/projection/signed-area operations and robust arc-aware winding (v-monotonic sub-arc splitting). Added a `TrimmedCurveBoundaries` unit test (line/arc primitives, disk, annulus, mixed line+arc half-disk, open-loop rejection). Header object-compiles with C++17 and the curve areas/classifications were verified numerically; the existing planar path is unchanged. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-23: Completed the cylindrical/spherical/conical bounded-surface milestones plus the supporting infrastructure they required.
	- Interface change: the single-hit virtual `BoundedSurface::intersectRay` was replaced by `appendIntersections(...)` returning every ray/patch intersection (a `RayHit` = distance + oriented normal) because quadric patches can be crossed twice per ray; a nearest-hit non-virtual `intersectRay` wrapper remains for callers that only need the closest hit, and `DistFromInside`/`DistFromOutside` prune with a shrinking `maxDistance`.
	- `Contains` parity was upgraded: near-equal intersection clusters that mix entering and exiting normals are tangential edge/corner grazes and now contribute even parity (found via a hollow-cylinder ray passing exactly through the inner rim corner).
	- Added `CurvedPlanarBoundedSurface` (plane with orthonormal frame, trimmed by `CurveWire` line/arc loops) so cylinder/cone caps are exact disks/annuli; exposed as `AddPlanarDiskSurface`. A planned shared `PlaneFrame` refactor was dropped: the curved class requires orthonormal axes, so the polygonal class's metric machinery would be dead weight.
	- Robustness fix in `CurveWire::classify`/`Curve2D::rightwardCrossings`: crossing counts now use canonical shared loop endpoints (each curve's start plus the next curve's start), because floating-point seam drift (e.g. `sin(2*pi) != 0` for a one-arc full circle) broke the half-open scanline convention for points aligned with the seam.
	- All curved-surface boundary sampling shares the `kArcSamples` (= 24, divisible by 4) full-circle chord count with `CurveWire::sampledBoundary`, so shared rims cancel exactly in the closure half-edge check; the cylinder/cone/sphere tests validate `IsClosed()`/`IsOrientationConsistent()` on mixed planar/quadric solids.
	- Testing policy per review feedback: `Safety` is never compared against ROOT shapes (TGeo safeties are legitimate underestimates and e.g. `TGeoCone::Safety` returns 0 on the axis of an `rmin = 0` cone); safety checks use analytic expected distances instead. `Contains`/distance functions are still compared against `TGeoTube`/`TGeoSphere`/`TGeoCone` on grids and ray sets.
	- Validation: full `DetectorsBase` ninja build in `/data/swenzel/sw/BUILD/O2-latest/O2` and the focused `ctest -R BVHSurfaceSolid` run pass (15 test cases). The `testMatBudLUT` failure in the same label is unrelated (missing `libO2TPCSimulation.so` in the partial build). `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-23: Completed the first two navigation milestones: `CloseShape`/BVH construction and BVH-accelerated `Contains` (details in the milestone notes above). Key implementation points: BVH lives in the transient `Impl` (`std::unique_ptr<Bvh>` built from per-surface conservative AABBs expanded by the new documented `kBVHBoxTolerance` and rounded outward in the float conversion); generic `visitRayCandidates` (bvh2 ray traversal) and `visitPointCandidates` (explicit stack, point-in-box) visitors on `Impl` will be reused by the upcoming `DistFromOutside`/`DistFromInside`/`Safety` milestones; `max_leaf_size = 1` (the bvh2 traversal does not box-test a leaf start node, and analytic patches are expensive - found via a diagnostic ray-candidate count of 6 instead of 0 on a miss ray when the whole 6-face box collapsed into one root leaf). `DistFromOutside`/`DistFromInside`/`Safety`/`ComputeNormal` still loop over all surfaces - they are the next milestones. Environment note for running builds/tests: use `eval "$(alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"` (comma-separated package list), then `ninja`/`ctest` in `/data/swenzel/sw/BUILD/O2-latest/O2`. Focused `ctest -R BVHSurfaceSolid` passes (17 test cases).
- 2026-07-23: Completed the first two CAD conversion milestones (surface-report probes and the sidecar format), enabling data-driven testing with real CAD input.
	- `O2_CADtoTGeo.py`: new exact-surface classification section (`classify_face`, `face_supported`, `build_surface_report`), a `def_shapes` leaf-shape store, the `--surface-report <path>` flag, and the sidecar writer `write_surfaces_bin` (not yet called by the pipeline; the extraction milestones will hook it up). Default output verified byte-identical.
	- Sidecar format v1 documented in this file ("Surface sidecar format"); C++ reader `o2::base::LoadSurfaceSolid` in `DetectorsBase` (new files `O2SurfaceSolidIO.h/.cxx`, registered in CMake, no dictionary entry needed). `testBVHSurfaceSolid.cxx` gained a `boxFaceFrame` helper refactor and the `SurfaceSidecarRoundTrip` case (18 test cases total, all pass).
	- Environment notes: pythonOCC needs a PYTHONPATH fix after `alienv` (`export PYTHONPATH=/data/swenzel/sw/slc9_x86-64/pythonOCC/v7.9.3-local1/lib/python3.10/site-packages:$PYTHONPATH`). An incremental ninja build after the CMakeLists change failed in rootcling dictionary regeneration (stale captured LD_LIBRARY_PATH); a full `./alibuild/aliBuild build O2` fixed it. For interpreted ROOT macros against the solid: use `R__ADD_INCLUDE_PATH($O2_ROOT/include)` + `R__LOAD_LIBRARY(libO2DetectorsBase)` inside the macro (exporting `ROOT_INCLUDE_PATH` breaks ROOT's C++ modules), and declare the `LoadSurfaceSolid` prototype manually since `O2SurfaceSolidIO.h` is not a dictionary header.
	- Support-matrix caveat surfaced by the probes: the C++ quadric API only accepts parametric-rectangle trims, so `auto`-mode conversion of quadric faces with general trim wires will fall back to tessellation until general `(u, v)` wires on quadrics are implemented; the report's `trim_kind` field tracks how often this occurs in real CAD data.
	- `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-23: Completed the "generated C++ macro support" milestone (details in the milestone note above). Python-only change to `O2_CADtoTGeo.py`: `emit_cpp_prelude` grew an `exact_surfaces` flag, `emit_surface_solid_cpp` sits next to `emit_tessellated_cpp`, and `emit_root_macro(surface_files=...)` selects per-lid between exact and tessellated emission (tessellated remains the fallback and the default). The generated exact block reuses the interpreted-mode recipe recorded in the previous note verbatim, so the macro works both interpreted and compiled. Facet files are still written even for exact volumes (harmless; the mode-flags/extraction milestones can revisit). Next milestone: `--exact-surfaces off|auto|required` CLI wiring, which will populate `surface_files` from `face_supported` eligibility plus the (still missing) face-extraction code that turns classified faces into `write_surfaces_bin` records.
- 2026-07-23: Completed the "converter mode flags" and "planar face extraction" milestones (details in the milestone notes above). Python-only change to `O2_CADtoTGeo.py`:
	- New extraction section: vector helpers, `extract_planar_face` (plane -> sidecar `plane` record via `BRepTools_WireExplorer` + outward-normal-preserving `axisV` sign choice), the `_FACE_EXTRACTORS` dispatch (only `plane` wired up), and `extract_surfaces_for_shape` (all-faces-or-nothing per leaf solid).
	- New `--exact-surfaces off|auto|required` CLI flag threaded into `emit_root_macro`; `auto`/`required` run the extractor over `def_shapes` after `extract_graph`, write per-volume sidecars, and drive the exact-vs-tessellated emission. `required` aborts with a per-volume/per-face reason report. `off` output is byte-identical to the pre-flag default. The explicit `surface_files` test hook is now honoured only when `exact_surfaces == "off"`.
	- Design note: extraction success (not the `--surface-report`/`face_supported` eligibility) is the source of truth for `auto`/`required`. Reconciled 2026-07-23 once the quadric extractors and general planar line+arc wires landed: `_SUPPORTED_PLANAR_CURVES = {line, circle}` now matches what the extractor actually emits (mixed line+arc planar wires), so the eligibility report and the emitter agree. The remaining superset is quadric faces with non-parametric-rectangle trims, which the report may still call eligible but the extractor rejects.
	- Validation: `--help`, one smoke run per mode on a generated centered box STEP, cylinder STEP `auto`-fallback / `required`-fail; the extracted box sidecar loads via `LoadSurfaceSolid`, closes with consistent orientation, `Capacity` = 0.192 exactly, 0 `Contains` mismatches vs `TGeoBBox` over 2e5 points. No C++ changes, so no rebuild needed. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
	- Next milestone: cylindrical/spherical/conical face extraction — add extractors to `_FACE_EXTRACTORS` reusing the C++ `Add{Cylindrical,Spherical,Conical}Surface` parametric-rectangle API, preserving OCC face orientation for outward normals.
- 2026-07-24: Real-CAD assessment on `scripts/geometry/STEP_examples/ALICE3_CAD_pure.step` (30.7 MiB, 103 leaf
	placements, **55 unique unplaced parts**) — no code change, planning only. Ran `--print-tree`,
	`--surface-report`, and a full `--exact-surfaces auto --mesh` sweep, and cross-checked `Contains` of the
	exact solid against the accelerated `O2Tessellated` (not ROOT `TGeoTessellated`, whose `Contains` fills
	concavities on thin shells and is unreliable). Findings that shaped the new milestones:
	- **Coverage 3/55 (~5.5%)** today (family `ST1A38494_002/003/004`, thin curved channels). Those 3 are
	  exact: `Contains` matches `O2Tessellated` to 1.4e-4 over 681k points (residual = mesh precision band)
	  and capacity matches the divergence-theorem mesh volume to <0.05%. So the kernel is correct; the ceiling
	  is the CAD *extractor's* surface/trim-type coverage.
	- **B-spline is the dominant blocker**: 44/55 volumes touch it, 34/55 are blocked *only* by it; all
	  analytic capabilities combined (quadric-curved trims + torus + revolution + extrusion) top out at
	  ~11/55. The bspline story splits into **trim curves** (cheap; surface stays exact) and **surfaces**
	  (expensive; exactness-vs-tessellation trade-off).
	- Concrete study of `ST1A38495_01`/`ST1A38526_01`: all-analytic surfaces (plane/cyl/cone), blocked purely
	  by bspline *trim* curves — 0 of the blocking quadric faces are pure line/arc, so the arc-pcurve fix does
	  not help them; a 2D bspline trim-edge type does.
	- Added three concrete surface milestones (bspline trim curves → torus → bspline surfaces) under "Surface
	  representation milestones" and a "Subdivision-BVH acceleration for curved surface patches" performance
	  milestone under the navigation section. **Recommended next implementation target: bspline trim curves**
	  (surgical, high value, surface stays exact). Also re-confirmed the still-open pre-existing planar
	  line+arc `kJoinTolerance` reader bug (did not trigger on this dataset because blocking faces are rejected
	  upstream, but remains a real one-line-ish fix for D-shaped caps).
- 2026-07-24: Completed "Support general curved (line/arc `CurveWire`) trims on quadrics" (details + design in the milestone note above). Quadric surfaces now take an optional line/arc trim wire in their `(phi, h)`/`(phi, theta)` domain (kernel `BoundedSurface.h`, public `Add*Surface` overloads, sidecar reader `collectQuadricTrim`, Python `_quadric_line_wire` + the stricter `_quadric_trim_fills_uv_box` scalar-path gate). Navigation stays exact; wire-trim capacity is numeric/inexact (flagged) and `Safety` a conservative lower bound. Full `ninja` build + `ctest -R BVHSurfaceSolid` pass (23 cases); pythonOCC smoke confirms scalar-path regression and a notched-cylinder line wire block round-tripping through `LoadSurfaceSolid`. Discovered a **pre-existing** planar line+arc join-tolerance issue in the sidecar reader (`kJoinTolerance` `1e-9` too strict vs `extract_planar_face`'s ~`1e-6` arc-endpoint precision) that blocks whole box-cut solids containing D-shaped caps — left for a focused planar fix. Env note: pythonOCC also needs `LD_LIBRARY_PATH+=/data/swenzel/sw/slc9_x86-64/OCCT/v7.9.3-local1/lib` (for `libTKFeat.so.7.9`, pulled in by the new `Geom2dAdaptor` import) on top of the usual `PYTHONPATH` fix. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-24: Completed the **B-spline trim curves on analytic surfaces** milestone (full details in the milestone note above). Kernel adds a `Curve2D::BSpline` kind (de Boor eval, GL-per-span exact area, polyline-based winding/distance with a lazily-cached flattened polyline that is essential for `CloseShape`/capacity performance), public `PlanarBoundaryCurve::makeBSpline`, sidecar `curveType=2`, and the affine-pole-transform quadric extractor `_quadric_trim_wire` (which also subsumes the old arc-on-quadric ellipse gap). Also folded in the long-standing planar line+arc join fix: reader `kJoinTolerance` and a new kernel `kWireJoinTolerance` are both `1e-5` (the kernel `CurveWire` `1e-9` closure was the actual gate, not just the reader). `ctest -R BVHSurfaceSolid` green (26 cases). End-to-end: `ST1A38495_01`/`ST1A38526_01` (previously bspline-trim-blocked) convert under `--exact-surfaces auto` and match `O2Tessellated` `Contains` with **0 mismatches / 20k points** (BVH==loop too), capacity within 0.17-0.24% of the mesh volume; **ALICE3 coverage 3/55 -> 7/55**. `IsClosed()` may warn on the shared-3D-bspline-edge sampling mismatch (documented caveat; navigation exact). **[Superseded 2026-07-26: this caveat is wrong. On `Bagger/BoomCylinderOuter_0_1_1_9` the same mismatch leaves real sliver gaps that flip `Contains` to "inside" up to 1.71 cm from any patch. See the "B-spline seam gaps" open item in the milestone list above.]** Remaining ceiling: bspline *surfaces* (largest effort, next-to-last milestone) and torus. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-24: Completed the **Bounded torus** milestone (full details + the capacity closed form in the milestone note above). New `TorusBoundedSurface` (kernel `BoundedSurface.h`) with a quartic ray intersection via a shared `solveQuarticReal`/`solveDepressedCubic` (Ferrari + Newton polish; even-sized root clusters dropped as tangencies for parity), exact rectangle capacity + numeric wire-trim capacity, meridian-distance `Safety`, and **both** ring/tube angles unwrapped in `pointInTrim` (the torus is the first surface with two periodic parameters). Public `AddToroidalSurface` (+ wire-trim overload), sidecar `surfaceType=5` (15 params), Python `extract_toroidal_face` (reuses the quadric `_quadric_trim_wire`/`_quadric_trim_fills_uv_box`; the `(u,v)->(phiRing,phiTube)` map is identity/mirror). `ninja` + `ctest -R BVHSurfaceSolid` green (**30 cases**: `ToroidalSurfaceKernels`, `FullTorusMatchesTGeoTorus`, `WireTrimmedTorusMatchesSection`, `TorusSidecarRoundTrip`). End-to-end: a synthetic full-torus STEP converts under `--exact-surfaces required` (1/1 exact) and its Python sidecar loads via `LoadSurfaceSolid` with **0 `Contains` mismatches over 9261 points vs `TGeoTorus`** and exact capacity `2*pi^2*R*r^2`. **Not** re-measured this session: the ALICE3 coverage number after adding torus (the 15 torus-touching volumes, incl. the near-miss `ST2487455_002`, should now improve on 7/55 — a follow-up `--exact-surfaces auto` sweep should confirm and update it). Remaining ceiling: bspline *surfaces* (the last, largest milestone). `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-26: **ALICE3 coverage re-measured after the torus milestone — measurement only, no source change** (this file is the only edit). Command: `python3 O2_CADtoTGeo.py STEP_examples/ALICE3_CAD_pure.step --output-folder <tmp> --exact-surfaces auto --mesh --surface-report <tmp>/surface_report.json`; env needs `ALIBUILD_WORK_DIR=/data/swenzel/sw` before `alienv` on top of the usual pythonOCC `PYTHONPATH` / OCCT `LD_LIBRARY_PATH` fixes.
	- **Coverage 7/55 -> 15/55 (27%).** The eligibility report and the extractor now agree exactly (15 eligible == 15 extracted), closing the long-standing "report is a superset of the extractor" gap.
	- The 8 new volumes are exactly the torus-containing ones: `ST1782525_01`, `ST1829909_01/_002/_003/_004`, `ST2487455_002`, `ST2487459_01`, `ST2487462_01`. The predicted near-miss `ST2487455_002` converts with precisely the forecast composition (1 cylinder + 2 cones + 2 planes + 1 torus). The 7 pre-existing volumes (`ST1A38494_002/003/004`, `ST1A38495_01`, `ST1A38526_01`, `ST2487455_01`, `ST2487721_01`) contain no torus at all, so the torus milestone alone accounts for the whole 7 -> 15 jump.
	- **Validation against the accelerated `O2Tessellated`** (not ROOT `TGeoTessellated`), fixed seed, points in the inflated mesh bbox: **14 of the 15 load and navigate correctly**. `Contains` mismatches run 0-17 per 20k points and *every* mismatching point sits within 1.7e-3 cm of the exact boundary (measured with `Safety`), i.e. purely the chorded-mesh precision band — there the reference is wrong, not the exact solid. The three largest loadable parts (236 / 720 / 965 surfaces) gave **0 mismatches over 5k points each**. `Contains` (BVH) == `Contains_Loop` on every point of every volume, and `IsOrientationConsistent()` is true everywhere. Capacity matches the divergence-theorem mesh volume to 0.001-0.27%, except `ST2487462_01` at 2.2% (bspline-trimmed faces -> numerically integrated capacity, already flagged `capacityIsExact()==false`).
	- **One real gap found: `ST1829909_01` (1052 surfaces) extracts but does NOT load.** `LoadSurfaceSolid` rejects it with "surface 1006: wire edge 1 end does not join the next edge start". Measured the actual drift over all 15 sidecars (17212 edge joins): only this volume exceeds `kJoinTolerance`/`kWireJoinTolerance = 1e-5`, and only at **6 joins, worst 2.96e-5** (~3x over), all on cylinder wire trims. So the honest number is *15/55 extracted, 14/55 usable*. Worth noting before anyone just loosens the constant again: the tolerance is applied to a **mixed-unit** `(phi[rad], h[cm])` metric, so the same absolute value means a different physical precision depending on the cylinder radius; the principled fix is a radius-aware (arc-length) join metric on quadrics. Deliberately **not** changed in this session — a looser tolerance also accepts genuinely broken wires, so it deserves its own focused decision plus a regression test.
	- **Remaining blockers re-ranked.** Face counts over the whole assembly: plane 3321, cylinder 2724, bspline 2377, cone 409, torus 350, extrusion 73, revolution 12. Of the 40 fallback volumes, 33 are blocked *purely* by bspline, 3 purely by revolution, 1 purely by extrusion, 3 by a mix. A basis-curve probe (downcast to `Geom_SurfaceOfRevolution` / `Geom_SurfaceOfLinearExtrusion`, then `BasisCurve()`) showed every one of those extrusion/revolution faces sweeps or revolves a **bspline** basis curve — they are not analytic-in-disguise, so there is no cheap +4 to harvest. **15/55 is the hard analytic ceiling**; bspline surfaces are the entire remaining milestone.
	- `ctest -R BVHSurfaceSolid` re-confirmed green (30 cases) against the installed build. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-26: **Coverage sweep over three non-ALICE CAD models — measurement only, no source change** (this file is the only edit; it adds the "Ellipse trim curves on planar faces" milestone the sweep motivated). Same command and environment as the ALICE3 sweep above, plus a load check of every emitted sidecar through `LoadSurfaceSolid` + `CloseShape` (no navigation comparisons this session).

	| model | volumes exact | sidecars loading | blocker |
	| --- | --- | --- | --- |
	| `Bagger.step` (503 KiB, mm, 13 vol / 288 faces) | **12/13 (92%)** | 12/12 | 1 volume, **ellipse** planar trim curves (8 curves on 2 faces) |
	| `oTOF System V3-R92cm.step` (5.3 MiB, mm, 3 vol / 1607 faces) | **2/3 (67%)** | 2/2 | 1 node with *zero* faces (empty assembly label, not a surface-type gap) |
	| `as1-oc-214.stp` (424 KiB, mm, 5 vol / 53 faces) | **0/5** | n/a | **bspline surfaces** (28 faces), exclusively |

	- These sit far above the ALICE3 assembly's 15/55 (27%), so **27% is not representative of mechanical CAD in general** — it reflects how ALICE3 was authored. On models that use analytic surfaces the converter is already near-complete.
	- **Bagger** is the useful one: `plane 172, cylinder 110, sphere 4, torus 2`; curves `line 902, circle 458, bspline 42, ellipse 8`. Fully analytic. Its one fallback (`Bucket`) is blocked *only* by the 8 ellipse curves — see the new milestone. Eligibility report and extractor agree exactly on all three models.
	- **oTOF** is 100% planar (`plane 1607`, `line 9116`, nothing else). Both real solids convert; `oTOF v2 v1` is 1505 surfaces. "2/3" understates it — every actual solid converts.
	- **as1** (the classic OCC demo assembly) initially read as the pure counter-example — the converter reports `plane 25, bspline 28` and 0/5, with no cylinder/cone/torus at all. **That reading was wrong, and chasing it produced the most important finding of the session** (see the follow-up note below): every one of those bspline faces is an *exact* cylinder. The lesson is recorded as its own milestone: **the stored STEP surface type is not the geometry.**
	- **FOLLOW-UP, same day — the stored surface type is not the geometry.** Prompted by disbelief at as1's 0/5 ("this is literally rectangles and rods"), the bspline faces were dumped and measured. Result: **as1 contains no free-form geometry at all.** Every bspline face is `degU=1, degV=3, poles=(2,4)` with weights `1, 1/3, 1/3, 1` over four poles on a rectangle — the textbook *exact* rational cubic Bezier for a **180-degree arc**. Fitting a cylinder to 401 samples of each gives R = 5.000000000 mm with max radial deviation **1.1e-11 mm (relative 2e-12)**, i.e. floating-point noise, not approximation error. Each hole is simply two exact half-cylinder patches. A generic geometry classifier (normal-field based: plane / sphere / cylinder / cone, relative residual < 1e-9) then measured both assemblies:
		- **as1: all 70 bspline faces are exact cylinders; 5/5 solids would convert. 0/5 -> 5/5.**
		- **ALICE3: 2889 of 8664 non-analytic faces (33%) are exact quadrics** — 1419 cylinder, 1398 cone, 72 sphere — and per solid, **15/55 already analytic + 14/55 rescued = 29/55**, with 26/55 genuinely free-form. The classifier's "already analytic" count reproduces the converter's measured 15/55 exactly, which cross-validates it.
		- So the previous session's **"15/55 is the hard analytic ceiling" conclusion is wrong** (corrected in place in the bspline-surface milestone). Its "bspline basis curve" probe was right about the storage and wrong about the geometry: a *rational* bspline basis curve is routinely an exact circle. Recognition is **exact recovery, not tolerance fitting**.
		- Caveat on the numbers: the classifier tests plane/sphere/cylinder/cone only — **no torus** — so 5775 "free-form" ALICE3 faces is an *upper* bound and 29/55 a *lower* bound. Also, per-face recognition was measured on `TopoDS` solids deduplicated by `TShape`; the exact converter works on XCAF logical volumes, so treat 29/55 as a forecast to be confirmed by an actual `--exact-surfaces auto` run.
		- Dead end worth not repeating: **OCCT's own `BRepLib_CanonicalRecognition` (7.7+) is NOT exposed in this pythonOCC build** (`ImportError` from `OCC.Core.BRepLib` in pythonOCC v7.9.3-local1). The recognizer has to be our own numeric one.
	- **Closure diagnostics (observed, not chased).** All 14 loaded solids are `IsOrientationConsistent() == true`; `IsClosed()` warns on 6. Three Bagger volumes (`Boom/Stick/BucketCylinderOuter`) carry bspline trims → the documented shared-3D-bspline-edge chord-sampling caveat. **Two cases are not explained by that caveat and are new**: (a) Bagger `BucketLink1` warns with 96 boundary edges but has *no* bspline (pure line+circle), while `BucketLink2` *does* have bsplines and does *not* warn; (b) both oTOF volumes warn with **0 boundary edges but 32 / 9 non-manifold edges** — a flavour not seen on ALICE3 (edges shared by >2 faces, typical of coincident faces or T-junctions in an all-planar tiled model). Neither blocks loading and neither touches navigation, but both should be understood before `IsClosed()` is trusted as a health signal on such models.
- 2026-07-26: **Implemented the canonical-form recognition pre-pass** (full details, measured numbers and
  the new gap in the milestone note above — this entry is the short version). First: reproduced all four
  2026-07-26 baselines fresh (`as1` 0/5, `Bagger` 12/13, `oTOF` 2/3, `ALICE3` 15/55) with the committed
  `--exact-surfaces auto --mesh --surface-report` sweep, confirming nothing had drifted; every sidecar
  load-checked via `checkSurfaceSidecars.macro` matched the prior session's counts exactly. Then implemented
  the recognizer in `O2_CADtoTGeo.py` only (`_recognize_analytic_surface`, `recognize_and_extract_face`,
  `_recognized_quadric_wire_block`; `extract_planar_face` gained a `frame_override` param; no C++ change,
  verified rather than assumed) and wired it into `extract_surfaces_for_shape`/`classify_face`/
  `build_surface_report` behind a new `--recognize-surfaces exact|off` flag (default on).
  **Key design finding, worth not re-deriving**: the milestone suggested reusing `_quadric_trim_wire`'s
  affine-pole-transform machinery for the reparametrization; checked and found it does **not** apply — the
  map from a recognized bspline surface's stored `(u, v)` to the recognized `(phi, h)` is only *separable*
  (each recognized coordinate a monotonic function of one stored parameter), not affine (a NURBS circular
  arc's own parameter maps to `phi` via a Möbius-type nonlinear function). The wire is instead rebuilt by
  sampling each 3D boundary edge directly and testing "is it exactly iso in the recognized domain", not by
  transforming the stored pcurve at all — simpler, and correct without needing the affine assumption.
  **Measured**: `as1-oc-214.stp` **0/5 -> 5/5**, exactly the forecast (`Contains` vs `O2Tessellated`: 0 or
  near-0 mismatches over 20k points/volume). `ALICE3_CAD_pure.step` **15/55 -> 20/55 extracted, 19/55
  usable** — short of the 29/55 forecast (that number was explicitly flagged in its own note as *"a
  forecast to be confirmed by an actual run"*; this was that run). The gap is almost entirely faces whose
  boundary is genuinely non-iso in the recognized frame (a real slanted/curved cut, out of scope for this
  pass by the milestone's own stated boundary), plus one pre-existing, unrelated `kJoinTolerance` rejection
  (`ST1829909_01`, already documented, not touched). `Bagger`/`oTOF` unchanged (no regression). All newly
  recognized ALICE3 volumes spot-checked against the accelerated `O2Tessellated`: 0 mismatches on 4/5,
  boundary-band-only (max `5.6e-4` cm) on the 5th at 200k points. `ctest -R BVHSurfaceSolid` green (30
  cases, unaffected — Python-only session). **New, deliberately-not-fixed gap found**: one ALICE3 volume's
  recognized faces load with `IsOrientationConsistent() == false` (cross-face wire-winding disagreement —
  root-caused to recognized faces lacking the `ax3.Direct()`-style external handedness hint that stored
  quadrics have; `Contains` is unaffected since it does not depend on cross-face winding). Remaining scope
  for next session: non-iso trims on a recognized quadric (the majority of the 29 vs 20 gap), curve
  recognition (bspline trim edges that are exactly line/circle/ellipse — the milestone's own "cheaper
  half", not started), torus recognition (explicitly deferred by the milestone), and the winding-consistency
  gap. `g4Config.C` and generated `scripts/geometry/STEP_examples/`/`facets_*`/`surfaces_*` artifacts
  remain intentionally out of the commit (all runs this session used a scratch output folder).- 2026-07-26: **Completed the navigation-harness / test-part-DB / BVH-distance block** — three milestones
  above ticked in one pass (`DistFromOutside`, `DistFromInside`, "Reusable test harness"), plus the five
  concrete wishes of the "Code optimization pass" milestone. Planned, executed and measured against
  [`SolidNavigationHarness.md`](SolidNavigationHarness.md), which is now the permanent usage documentation
  for the harness (how to build the part DB, how to run and read the benchmark, how to profile it under
  `perf`) and carries the full result tables; only the headlines are repeated here.
  New code: `DetectorsBase/O2SolidHarness.{h,cxx}` (deterministic sampling, validation with mesh-band
  classification, timing), `o2-bench-detectorsbase-solid-harness` (`Detectors/Base/test/runSolidHarness.cxx`,
  flags incl. `--loop-crosscheck` and `--pruning-ab`), `LoadFacetSolid` in `O2SurfaceSolidIO`,
  `scripts/geometry/makeTestPartDB.py`, and in `O2BVHSurfaceSolid` the BVH distance queries
  (`Impl::nearestCrossing<wantEntering>` — bvh2 lambda traversal, `ray.tmax` tightened from inside the
  leaf callback, bbox reject under `stepmax`, `thread_local` scratch) with public `DistFromOutside_Loop` /
  `DistFromInside_Loop` twins and the `SetRayTMaxPruning` / `GetRayCandidateCount` benchmark hooks.
  Five new test cases; `ctest -R BVHSurfaceSolid` green (**35 cases**, 0.80 s).
  Headline measurements over the 19-part three-model DB (seed 1, 3000 points / 3000 rays per part):
  **BVH == `_Loop` bit-identical on 57000/57000 + 28500/28500 distance rays**; ray-`tmax` tightening
  removes 26.6% of the analytic patch tests for ~1% wall time on parts this small (answers identical
  57000/57000); the exact solid costs ~1.3–2.4x the tessellated one per query at the median with a wide
  spread; 1851 analytic patches replace 48703 triangles (26x fewer primitives, median 184x per part).
  Following up on the harness's "unexplained" column then turned up the most important defect currently
  known in the exact-surface path, recorded as its own open item: **the converter discards BREP shared-edge
  identity**, trimming each face by its own pcurve (`_quadric_trim_wire` uses `BRep_Tool.CurveOnSurface`),
  so the topology that made the CAD model watertight is lost at extraction and adjacent patches disagree
  about where their common boundary is. On `Bagger/BoomCylinderOuter_0_1_1_9` the resulting sliver gaps
  flip `Contains` to "inside" **1.71 cm from the nearest patch** (confirmed: `Contains_Loop` agrees, so not
  the BVH; and the wrong region is a narrow tube aligned with the parity test direction, i.e. the gap's
  shadow, exactly as predicted). This supersedes the 2026-07-24 "navigation exact" caveat. Four ordered
  work items are recorded there and planned in detail in the new
  [`ExactTrimTopology.md`](ExactTrimTopology.md), headed by preserving shared edges in the converter.
  Two further findings were recorded as open items rather than fixed: `Safety` is the most expensive kernel
  (10132 ns/call, 4.74x the mesh) purely because it still has no BVH — the `Safety` milestone below is
  now the highest-value optimization — and `Contains` disagrees with `Contains_Loop` on 301/142500 points,
  295 of them in the two non-manifold `oTOF` parts (order-dependent parity clustering on duplicated faces;
  the distance twins are unaffected). Also recorded: mutation-testing showed the new sweeps catch a 2x
  over-prune loudly but not a 0.1% one, which is a property of the culling rule (a node is culled on its
  *box* entry, always <= its patch hit), so the pruning-bound guarantee rests on the argument documented
  at `nearestCrossing`, not on a test. Step 0 of the plan also cleared a latent `BUILD_TEST_ROOT_MACROS`
  hard-fail by renaming `checkSurfaceSidecars.C` -> `.macro`, and the first ALICE3 run showed
  `LoadFacetSolid` was too strict — one degenerate sliver triangle failed a whole 211k-triangle mesh —
  so degenerate facets are now skipped and counted rather than fatal. `g4Config.C` and generated
  `scripts/geometry/STEP_examples/`/`facets_*`/`surfaces_*` artifacts remain intentionally out of the
  commit (all runs this session used a scratch output folder).
- 2026-07-26: **Exact trim topology, items 4 and 3** (plan and full evidence in
  [`ExactTrimTopology.md`](ExactTrimTopology.md); milestone checkboxes in the "converter drops BREP edge
  sharing" open item above). Two of the four work items opened by that finding are done, and the
  measurements taken along the way revised the finding itself.
	- **Item 4, fail loudly (kernel + harness).** `O2BVHSurfaceSolid::NavigationReliability`
	  (`Undetermined`/`Reliable`/`ReversedFaces`/`OpenSurfaceSet`/`NonManifold`, ordered by severity, worst
	  defect wins) with `IsNavigable()`, `GetNavigationReliabilityName()` and the raw boundary /
	  non-manifold / reversed edge counts; `CloseShape`'s three closure defects are now `Error`s whose text
	  states the *consequence* rather than the count; `Print()` shows the state. The harness prints it per
	  part, marks a non-navigable part next to its accuracy columns, carries it in `--json` under
	  `navigation`, and repeats the list at the end of a run. Test `NavigationReliabilityIsQueryable`
	  (36 cases green). The conservative default was kept: an unnavigable solid still answers queries.
	  Measured: **11 of 19 three-model parts and 14 of 19 ALICE3 parts are not navigable** — that is the
	  honest denominator behind every accuracy row in `SolidNavigationHarness.md`.
	- **Item 3, canonical trim curves (converter).** `analyze_surface_geometry.py` gained `--trim-curves`
	  (classify trim curves by their real geometry; also measures how far apart the two faces of a shared
	  edge place their common boundary). `_recognize_canonical_curve` in `O2_CADtoTGeo.py` acts on it from
	  both `_quadric_trim_wire` and `extract_planar_face`. Straightness is decided on the **control
	  polygon** — collinear poles put the curve inside their convex hull, so it is a proof, not a fit — and
	  circles at `1e-9` relative; an ellipse and a line bowed by `1e-6` are correctly refused. Result:
	  B-spline trim edges **88 → 50** (three-model) and **15034 → 4528** (ALICE3, sidecars 6.09 MB →
	  3.65 MB), with **every `unexplained` count bit-identical** on both DBs and no part losing its exact
	  conversion. The null result *is* the acceptance criterion: recognition must change representation,
	  never geometry.
	- **What the measurements changed.** (a) Item 3 does **not** fix the reproducer — the curves it
	  recognises on `BoomCylinderOuter` are its four periodic cylinder seams, not the two
	  cylinder-cylinder intersection curves that cause the defect. (b) The polyline flattening is **not**
	  the dominant error at these seams (`bsplineSampleInto` is adaptive to `1e-5`, and the two pcurves of
	  a shared edge were measured to differ by at most `1.3e-5` model units on Bagger), so item 2 is a cost
	  item plus a small accuracy item, not the fix. (c) Enumerating the parity ray's crossings at a wrong
	  point shows the failure is a **doubled** crossing, not a missing one — two consecutive `ENTER`s ~0.04
	  cm apart — because the two cylinders meet **near-tangentially** there and amplify a `1e-5` parametric
	  disagreement into a `4e-2` cm one. The size of the error is set by the tangency, so tightening
	  tolerances cannot fix it. (d) The exact intersection curve of two quadrics is degree 4 and is
	  representable exactly in *neither* face's `(phi, h)` domain with the current `Curve2D` vocabulary, so
	  item 1 must take the "consistency, not exactness" route (one shared sample set per `TopoDS_Edge`,
	  both faces built from it); a concrete four-step recipe and three traps are written up in
	  `ExactTrimTopology.md`.
	- **Left open deliberately:** items 1 (shared edges — now the only remaining fix for this defect class)
	  and 2 (Bézier point-in-trim). Item 1 is a converter refactor with real risk of a half-finished
	  intermediate state, so it was scoped and handed over rather than started. Also not completed: the
	  `--trim-curves` sweep over `oTOF`/`ALICE3` (the shared-edge polyline comparison is a Python loop and
	  does not finish in an hour on those; vectorise `_max_point_to_polyline` first if that number is
	  wanted). `g4Config.C` and generated `scripts/geometry/STEP_examples/`/`facets_*`/`surfaces_*`
	  artifacts remain intentionally out of the commit.
- 2026-07-26: **Root cause of the "B-spline seam gaps" defect found and fixed — it was a one-line kernel
  bug, and `ExactTrimTopology.md`'s premise was a red herring.** `bsplineSampleRecursive` stopped
  subdividing whenever the chord `p0 -> p1` fell below the flatness scale. A **closed** B-spline has
  `p0 == p1` exactly, so a full circle — which is what a CAD kernel writes for a tube-tube intersection
  curve — flattened to two coincident points and disappeared from every polyline-based query: winding
  classification, closest point, boundary band, display mesh and the closure half-edge check. It hid for
  three sessions because `signedAreaContribution` integrates the curve by Gauss-Legendre rather than from
  the polyline, so the wire still validated and still reported its correct enclosed area; every check that
  could have caught it used the analytic path. The existing tests covered a *line* hole and a B-spline
  outer wire built from four *open* quarter arcs, so neither could see it.
	- **Measured effect** (three-model DB, `unexplained`): `contains` 4588 → 4430, and 525 → **367**
	  excluding the unchanged non-manifold `oTOF` part; `distin` 218 → **114**; the reproducer
	  `BoomCylinderOuter` 51 → **16** with its `distin` 30 → **0**. Every Bagger cylinder part's `distin`
	  goes to zero. Two parts whose `Contains` disagreed with `Contains_Loop` (`BucketCylinderInner`,
	  `BucketCylinderOuter`) now agree on all 7500 points. ALICE3 is unchanged in every column — its trims
	  arrive through the recognised-NURBS path as multiple *open* edges, so it never hit this.
	- **The closure diagnostic became honest.** `BoomCylinderInner`, `BucketCylinderInner` and
	  `StickCylinderInner` used to report 0 boundary edges and "navigable" — only because a whole face's
	  wire had vanished. DB-wide navigable drops 8/19 → 5/19. That is item 4 working, not a regression.
	- **Process lesson worth keeping.** Two of `ExactTrimTopology.md`'s four items were built before anyone
	  tested its central causal claim, which measurement then contradicted twice: the shared-edge pcurve
	  disagreement is `1.3e-5` model units (three orders too small), and the "near-tangency amplification"
	  offered as the explanation was inferred, not measured, and is wrong (the surfaces meet at 59-60°).
	  One throwaway probe enumerating the crossing list moved this from a three-item converter refactor to
	  a one-line fix. Diagnose before planning.
	- **Still open:** 367 `contains` disagreements outside the non-manifold part, and no part closing yet.
	  Diagnose those the same way before resuming any converter work. Items 1 and 2 of
	  `ExactTrimTopology.md` are marked "do not start" pending such evidence.
- 2026-07-31: **Comprehensive review + Phase 0 (test infrastructure) — commit `09f20bc32f`.**
  A full review of the branch is in [`CodeReview_Fable.md`](CodeReview_Fable.md); its Sections 3-4
  derive why per-face 2D trim curves can never make a tube-tube part watertight, Section 9 defines
  the new success gates, Section 10 the phased plan, and **Section 11 records the measured Phase 0
  baseline**. The durable design decisions it produced are at the top of this file under
  "Adopted design principles (2026-07-31)". Phase 0 built the missing reference and fixtures; it
  deliberately changed **no** navigation or converter semantics.

	**Everything is committed and `ctest -R BVHSurfaceSolid` is green (36 cases).** A new session can
	start Phase 1 from a clean tree.

	**How to run the gate** (all paths relative to the O2 checkout):
	```bash
	# build (incremental; the branch's targets appear after a cmake re-run)
	export ALIBUILD_WORK_DIR=$HOME/alisw/sw
	B=$HOME/alisw/sw/BUILD/O2-latest/O2
	cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
	ninja o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-solid-harness
	# NOTE: ctest and the harness resolve libO2DetectorsBase from the *installed* prefix unless the
	# build's stage dirs come first -- otherwise you get an undefined-symbol failure or, worse, a
	# silently stale measurement.
	LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH ctest -R BVHSurfaceSolid

	# the whole gate in one command (converts, dumps samples, runs the oracle, scores)
	O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
	O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
	    --model scripts/geometry/STEP_examples/Bagger.step
	# validate the oracle itself first if anything looks surprising -- it is the trust anchor
	$HOME/alisw/sw/ubuntu2404_aarch64/Python/latest/bin/python3.10 scripts/geometry/occtOracle.py --self-test
	```
	pythonOCC needs the alibuild Python 3.10 (the system `python3` is 3.12 and cannot import `OCC`);
	`runOracleGate.py` sets that environment itself. Never write generated artifacts into
	`scripts/geometry/STEP_examples/` -- use a scratch folder, as every session here has.

	**Baseline to beat (measured, seed 1, 2000 points / 2000 rays):** fixture ladder **6/10**,
	`Bagger.step` **4/12**. The three failing fixtures are exactly the three with a genuine
	cylinder-cylinder intersection curve (`cyl_cross_cyl`, `cyl_inter_cyl`, `tube_window`), and
	`oblique_cut_cyl` does not convert at all (the ellipse planar-trim blocker, i.e. Bagger's
	`Bucket` reduced to three faces). Full table in `CodeReview_Fable.md` Section 11.

	**Phase 1 is next**, in this order (details and file:line references in `CodeReview_Fable.md`
	Sections 5-6; each item wants its own regression test, and the gate above is the acceptance
	measurement -- run it before and after every item so the effect is attributable):
	1. **Ambiguity band + re-shoot parity** (Section 4.4). `pointInTrim` is a hard binary and
	   `Contains` casts one fixed direction with no retry, so a hit landing within epsilon of a trim
	   boundary is classified by coin flip and the error is the whole shadow of the seam. Expected to
	   move the `contains` columns the most.
	2. **The safety-critical solid bugs**: S1 (a streamed-and-read solid comes back *empty* and
	   reports itself `Reliable` -- fix first, it is actively harmful), S2 (`Contains` is disabled
	   before `CloseShape`), S3/S5 (boundary policy differs between `Contains` and the distance
	   queries), S4 (distance queries lack the mixed-cluster cancellation `Contains` has).
	3. **Kernel bugs K1-K8**, notably K1 (B-spline clamping assumed, never validated), K2 (full-turn
	   trims rejected on the pole hull), K3 (`kWireJoinTolerance` applied to mixed rad/cm units),
	   K7 (a face that fails to build is silently dropped from the parity solid).
	4. **Tolerance policy and the closure criterion** (K3/S10, K9/S8): per-domain metrics
	   (angular <-> length via radius) in kernel *and* IO, and a closure check matched at the
	   topology level with a quantitative gap metric, so "closed" becomes achievable and meaningful.

	**Diagnose before planning.** `DescribeContainsCrossings` now prints both parity hit lists when
	BVH and loop disagree (`--loop-crosscheck`); the recorded "order-dependent clustering"
	explanation for those disagreements is **not supported by the code** (the clusterer sorts first,
	so it is a function of the sorted distances alone -- see S6). Measure it before fixing it. Phase 2
	(adjacency-based exact trims, the actual fix for the tube-tube class) starts only after Phase 1.
- 2026-07-31: **Phase 1, items 1-3 — commits `63d1d08119`, `29e8322f79`, `f809b38dd2`,
  `346d4675d0`.** Full account, with every measurement, in
  [`CodeReview_Fable.md`](CodeReview_Fable.md) **Section 12**; only the headlines and the handover
  are here. Item 4 of Phase 1 is *not* started. `ctest -R BVHSurfaceSolid` green, **48 cases**
  (from 37). The tree is clean and a new session can start from item 4.

	**Headline: the containment error class is gone.** `contains` disagreements with the
	OpenCascade oracle outside tolerance: fixtures **2 -> 0**, Bagger **56 -> 2**. Gate *totals*
	are unchanged (**G1 6/10, G3 4/12**) because every part that still fails does so on
	navigability, on the distance queries or on capacity — none of which these items touch. Both
	statements are true at once and neither should be quoted without the other.

	**S1 was done first**, ahead of the recommended order, because it was the one defect that was
	actively harmful: `TGeoManager::Export/Import` silently replaced every one of these solids with
	an empty shape that reported itself `Reliable`. The solid now persists its `Add*Surface` call
	sequence and replays it on read.

	**The re-shoot was measured before it was built**, and the measurement corrected the recorded
	picture. A new hook `ContainsAlongDirection()` makes the solid its own oracle: on a closed
	2-manifold parity cannot depend on where the ray is aimed, so any direction dependence *is* a
	defect, with no reference shape involved. Over the whole Phase 0 corpus, every `Reliable` part
	has **zero** direction disagreements in ~11k points and every part that disagrees is one the
	closure check already rejects; and of the 55 points where the fixed direction disagrees with the
	oracle, **not one is wrong in every direction**. A gap's shadow belongs to the *direction*, not
	to the point. So `Reliable` solids keep the single shot (1.0-1.5x) and unreliable ones vote over
	five golden-spiral directions (2.9-4.7x). Use `ContainsAlongDirection` for the next diagnosis
	too — it is the cheapest oracle in the codebase.

	**Two corrections to `CodeReview_Fable.md`**, both recorded there: **K7 does not hold** (a face
	that fails to build is not silently dropped on any production path — the loader rejects the
	whole file and the generated macro throws), and the shadow of a gap is escapable rather than
	lost. Neither was "fixed"; both were pinned by tests, in the spirit of the review's own S6
	correction.

	**New in the test suite:** the **concave fixture the review found missing** — an L-shaped prism
	with a genuine reflex edge, the only configuration where a ray can touch the boundary from
	inside and stay inside — carrying the full sweep battery; a persistence round trip over every
	surface family and both trim flavours; and regression tests for the boundary policy, edge
	grazes, unclamped B-spline endpoints and the full-turn trim rejection.

	**Where to pick up, in order.**
	1. **Phase 1 item 4** (the only Phase 1 work left), planned in detail in
	   [`TolerancePolicy.md`](TolerancePolicy.md): per-domain metrics (angular <-> length via
	   radius) in kernel *and* IO (K3/S10), sidecar v2 carrying the model tolerance, and a closure
	   criterion matched at the topology level with a quantitative gap metric (K9/S8) so that
	   "closed" becomes achievable and meaningful. **K5 belongs with it** — a 1e-9 boundary band
	   tested against a ~1e-5 polyline needs exactly the same metric, and doing it separately means
	   building the metric twice.
	2. **Diagnose `BucketLink2` before anything else in Phase 2.** It is the one Bagger part that is
	   **navigable** and still has 24 *missed* crossings in `distout` (48 in `distin`, 6.3% capacity
	   drift). Every other failure is explained by an open surface set; this one is not, and nothing
	   in Phase 1 moved it. Probe it the way item 1 was probed — evidence, then theory.
	3. **K4 and K6** remain open and untouched (degenerate-chord recursion; cancellation in the
	   naive quadratic formula plus absolute, scale-dependent tolerances in the cone and torus
	   branches). K6 is still a live suspect for distance-query residue.
	4. Only then Phase 2 (adjacency-based exact trims), which is what the remaining distance and
	   capacity columns are actually made of.

	**Traps that cost time here, worth carrying forward.** Put `$B/stage/lib:$B/stage/lib64` first
	on `LD_LIBRARY_PATH` or `ctest` and the harness silently resolve stale installed libraries.
	pythonOCC needs the alibuild python3.10, not the system 3.12 (`runOracleGate.py` sets this
	itself). Never write generated artifacts into `scripts/geometry/STEP_examples/` — every run this
	session used a scratch folder. And a ROOT macro is the wrong vehicle for a quick kernel probe on
	this machine (cling trips over unrelated headers); compiling a 40-line `.cxx` against
	`$B/stage/lib -lO2DetectorsBase` takes a minute and just works.
- 2026-07-31: **`BucketLink2` diagnosed, and Phase 1 item 4 steps 1-4 of five** — commits
  `5716d0070d`, `1bc5c4fbc9`, `9f45887ef7`, `cacd64e4a5`, `f612f895a9`. Full account in
  [`CodeReview_Fable.md`](CodeReview_Fable.md) **Section 13** (the diagnosis) and
  [`TolerancePolicy.md`](TolerancePolicy.md) §5 (the four steps, with their measurements).
  `ctest -R BVHSurfaceSolid` green, **53 cases** (from 48). Tree clean.

	**The best remaining lead turned out not to be a kernel defect at all.** `BucketLink2` was the
	only navigable Bagger part still failing the gate — 24 missed `distout` crossings, 48 `distin`,
	6.3% capacity drift — and no theory covered it. Three throwaway probes settled it. Its 24 + 48
	distance disagreements are an artifact of the *gate*: the sample generator assigns a ray to the
	"outside" or "inside" category using the tessellated reference, and that mesh is not watertight
	(its own `Check` reports hundreds of two-neighbour facets) and puts the part's left plate half a
	centimetre off. OpenCascade classifies all 20 recorded offender origins opposite to their
	assigned category; the exact solid agrees with OpenCascade on **20 of 20**, and where the
	category is right its distances match the oracle to every printed digit. The gate was asking
	`DistFromOutside` of an interior point and comparing it against an exit distance. And the 6.3%
	capacity drift is **quadrature, not geometry**: a 4M-point Monte Carlo of the exact solid gives
	17.061 ± 0.052 cm³ against OpenCascade's 17.079 while `Capacity()` returns 16.004.
	`integrateOverCurveTrim` is a fixed-128 midpoint rule over a characteristic function, so it
	converges at O(1/N); across the whole model the only four parts whose capacity matches the
	oracle exactly are the only four with no wire-trimmed quadric. So `BucketLink2` is not evidence
	for K4 or K6 and is no longer a lead; two precisely-scoped items replace it (gate ray-category
	soundness, and Green's-theorem capacity for wire-trimmed quadrics).

	**Phase 1 item 4, steps 1-4.** Every surface now reports its first fundamental form
	(`parametricMetric`), so a parametric separation can be turned into the length it actually
	spans. The wire join checks in the kernel *and* in the loader go through it and compare against
	one constant in cm, closing **K3, K12 and S10**: three places used to decide the same question
	with three different rules — 1e-5, 1e-9 and a bare per-coordinate 1e-5 — none of them a
	distance, although one extractor feeds all three. The constant moved 1e-5 → 1e-6 cm, measured
	rather than assumed (690 joins across both corpora, worst residual **4.06e-11 cm**). Sidecar
	**version 2** carries the model's own declared tolerance (`Bagger.step`: 1e-8 cm) so the kernel
	stops guessing what epsilon two faces of an imported solid should agree to. And the on-boundary
	band is sized from the representation's real accuracy instead of optimism, closing **K5**: a
	1e-9 band around a 1e-5 polyline made `Boundary` unreachable for every B-spline trim, and
	winding and distance measured against two different polylines, which is now one.

	**All four steps gate bit-identical** (fixtures 6/9, Bagger 4/12, `contains` 0 and 2). For steps
	1 and 3 that is the required outcome. For steps 2 and 4 it is the *expected* one and not
	evidence that nothing changed — no face on this corpus was ever near either join threshold, and
	a 1e-5-wide boundary shell is not something bbox-spread samples land in. Each is pinned by tests
	instead, and step 2's two tests were verified to fail against the pre-fix rule (five assertions,
	each named in the test's own comments).

	**Where to pick up.** Step 5 — rim-based closure and `maxGap` (K9/S8) — is the only Phase 1 work
	left, and `TolerancePolicy.md` **§8** now splits it in two: build the rim gap measurement and
	report it while changing no verdict (the gate must stay bit-identical, which is a real check
	that it has not leaked into one), and only then let it decide navigability, in the same commit
	as the mandatory direction-independence sweep. That split exists because a closure criterion
	that succeeds more often moves solids onto `Contains`'s single-shot fast path, whose licence is
	a measurement taken while the closure check under-reported `Reliable`.
- 2026-08-01: **Phase 1 item 4 step 5, and with it Phase 1, is finished** — commits `213b8aca8c`
  (5a, measure) and `d32196029e` (5b, decide). Full account in
  [`TolerancePolicy.md`](TolerancePolicy.md) **§9 and §10** (every table) and
  [`CodeReview_Fable.md`](CodeReview_Fable.md) **Section 14**. `ctest -R BVHSurfaceSolid` green,
  **57 cases** (from 53). Both gates bit-identical on both commits — fixtures 6/9, Bagger 4/12,
  `contains` 0 and 2. Tree clean. **K9 and S8 are closed.**

	**"How far apart are the faces, in cm?" now has an answer, and it is not the recorded one.**
	Closure is a curve comparison: each face emits rims (one ordered polyline per trim loop) and
	each rim chord is matched at its midpoint against the chords of every other face, at the
	model's own declared tolerance from sidecar v2. `TolerancePolicy.md` §3.3 predicted the Bagger
	cyl-cyl parts would stay open at ~1.3e-5 cm, the shared-edge pcurve disagreement. **They are
	open by 0.25 to 0.75 cm**, over 4-15% of their rim length — a face missing or trimmed to the
	wrong curve, not two pcurves disagreeing in the fifth decimal. Phase 2's premise for those
	parts changes accordingly: they are not nearly closed and then spoiled by tolerance.
	`tube_window` went from **1418 boundary edges** to seven rims of which three are open, 9.94 cm
	of 53 cm.

	**The other prediction was wrong in the opposite direction.** §1.3 said vertex matching fails
	on real CAD because the two faces of a shared edge sample it at different phases. The
	structural defect is real — a unit test builds a box whose last face carries twice the chord
	count and the old check calls that closed box open — but no part of either corpus exercises
	it. Every part the old check called closed is rim-matched and vice versa, part for part, which
	is why 5b's verdict switch is a gate no-op. "This corpus does not test the thing" is the right
	reading, not "the thing was not a problem".

	**The mandatory sweep passed and left one lead.** 21 parts, 11k points, 13 golden-spiral
	directions, before and after the switch: the same 13 `Reliable` parts, **one** disagreement in
	143000 points, and the parts that do disagree between directions are exactly the ones the
	closure check rejects (0.55%-7.0%). So `Contains`'s single-shot fast path keeps its licence.
	The offender — `cyl_cross_cyl` at `(0.65334264649649720, 0.88394684996026007,
	0.97463122696308724)`, 1 of 13 directions, `Safety` 0.0992 cm — is not a gap shadow and is now
	the cheapest known reproducer for **K6**.

	**One piece deliberately not shipped, measured rather than assumed.** §2.4's ambiguity band:
	the free form — re-shoot whenever a parity cluster held more than one coincident hit — costs
	0.2-1.3% of `Contains` on seven fixtures and **15.8%** on `box_minus_cyl` and moved **not one**
	of the 143000 sweep points, including the offender it was aimed at. Dropped, with the numbers
	left in `containsByParity`. The parametric-domain band §2.4 actually asks for is unbuilt and
	is deferred behind diagnosing that point.

	**Three things the plan did not anticipate**, all pinned by tests. The parametric-rectangle
	quadrics interleave their two rims rather than emitting each loop consecutively, so §8's
	"chain consecutive runs" does not work — chaining by matching endpoints does, and still needs
	no change in any of the seven surface classes. A full-turn patch whose `fullSweep()` does not
	fire emits its seam twice in opposite directions; that cancels in the half-edge check, so it
	had never been noticed, but it chains into a two-point rim straddling the patch and produced a
	spurious boundary rim on five of nine fixtures. And the sampling noise floor has to be
	estimated from the *turn angle* — a vertex's deviation from the chord joining its neighbours
	calls a box corner 2.4 cm of noise. `maxGap` is reported next to that floor and should never
	be quoted without it.
