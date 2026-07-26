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
	  `C2v=dv/2+(sin2v1-sin2v0)/4`. Not yet re-measured: the ALICE3 coverage number (torus now available for
	  the 15 torus-touching volumes; a follow-up sweep should confirm `ST2487455_002` and update 7/55).

- [ ] **B-spline / NURBS bounded surfaces** (bulk coverage; largest effort — do last).
	- Motivation: the dominant blocker. On ALICE3, 34/55 volumes are blocked *only* by bspline and 44/55
	  touch it; every analytic capability combined tops out at ~11/55. Majority exact coverage of a real
	  assembly fundamentally requires a genuine bspline surface kernel (surface + trim).
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

- [ ] Implement `DistFromOutside`.
	- Traverse the BVH with a ray and ask candidate surfaces for intersections.
	- Return the smallest positive entering intersection, respecting `stepmax` and normal orientation.
	- Acceptance test: compare against `TGeoBBox`/`TGeoTube` for outside points and several directions, including grazing misses.

- [ ] Implement `DistFromInside`.
	- Return the smallest positive exiting intersection, respecting normal orientation and avoiding immediate re-hit at `t = 0`.
	- Acceptance test: compare against ROOT primitives for central, near-boundary, and oblique rays.

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

Version 1 of the exact-surface sidecar (`surfaces_<VOLNAME>_<LID>.bin`), written by
`write_surfaces_bin` in `scripts/geometry/O2_CADtoTGeo.py` and read by
`o2::base::LoadSurfaceSolid` (`Detectors/Base/include/DetectorsBase/O2SurfaceSolidIO.h`).
All integers are little-endian `uint32`, all geometry values are little-endian `float64`,
lengths in cm, angles in radians.

```
header:
  char[4]  magic        = "O2SS"
  uint32   version      = 1
  uint32   nSurfaces
  uint32   reserved     = 0
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
- 2026-07-24: Completed the **B-spline trim curves on analytic surfaces** milestone (full details in the milestone note above). Kernel adds a `Curve2D::BSpline` kind (de Boor eval, GL-per-span exact area, polyline-based winding/distance with a lazily-cached flattened polyline that is essential for `CloseShape`/capacity performance), public `PlanarBoundaryCurve::makeBSpline`, sidecar `curveType=2`, and the affine-pole-transform quadric extractor `_quadric_trim_wire` (which also subsumes the old arc-on-quadric ellipse gap). Also folded in the long-standing planar line+arc join fix: reader `kJoinTolerance` and a new kernel `kWireJoinTolerance` are both `1e-5` (the kernel `CurveWire` `1e-9` closure was the actual gate, not just the reader). `ctest -R BVHSurfaceSolid` green (26 cases). End-to-end: `ST1A38495_01`/`ST1A38526_01` (previously bspline-trim-blocked) convert under `--exact-surfaces auto` and match `O2Tessellated` `Contains` with **0 mismatches / 20k points** (BVH==loop too), capacity within 0.17-0.24% of the mesh volume; **ALICE3 coverage 3/55 -> 7/55**. `IsClosed()` may warn on the shared-3D-bspline-edge sampling mismatch (documented caveat; navigation exact). Remaining ceiling: bspline *surfaces* (largest effort, next-to-last milestone) and torus. `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.
- 2026-07-24: Completed the **Bounded torus** milestone (full details + the capacity closed form in the milestone note above). New `TorusBoundedSurface` (kernel `BoundedSurface.h`) with a quartic ray intersection via a shared `solveQuarticReal`/`solveDepressedCubic` (Ferrari + Newton polish; even-sized root clusters dropped as tangencies for parity), exact rectangle capacity + numeric wire-trim capacity, meridian-distance `Safety`, and **both** ring/tube angles unwrapped in `pointInTrim` (the torus is the first surface with two periodic parameters). Public `AddToroidalSurface` (+ wire-trim overload), sidecar `surfaceType=5` (15 params), Python `extract_toroidal_face` (reuses the quadric `_quadric_trim_wire`/`_quadric_trim_fills_uv_box`; the `(u,v)->(phiRing,phiTube)` map is identity/mirror). `ninja` + `ctest -R BVHSurfaceSolid` green (**30 cases**: `ToroidalSurfaceKernels`, `FullTorusMatchesTGeoTorus`, `WireTrimmedTorusMatchesSection`, `TorusSidecarRoundTrip`). End-to-end: a synthetic full-torus STEP converts under `--exact-surfaces required` (1/1 exact) and its Python sidecar loads via `LoadSurfaceSolid` with **0 `Contains` mismatches over 9261 points vs `TGeoTorus`** and exact capacity `2*pi^2*R*r^2`. **Not** re-measured this session: the ALICE3 coverage number after adding torus (the 15 torus-touching volumes, incl. the near-miss `ST2487455_002`, should now improve on 7/55 — a follow-up `--exact-surfaces auto` sweep should confirm and update it). Remaining ceiling: bspline *surfaces* (the last, largest milestone). `g4Config.C` and generated `scripts/geometry/STEP_examples/` artifacts remain intentionally out of the commit.