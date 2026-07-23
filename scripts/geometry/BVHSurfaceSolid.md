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

- [ ] Define numerical conventions in code and tests.
	- Choose tolerances for ray `t`, boundary checks, duplicate intersection clustering, and point-on-surface classification.
	- Record these constants in one implementation-local namespace, not scattered literals.
	- Validation: unit tests for near-boundary points and near-tangent rays.
	- Partial 2026-07-22: implementation-local constants were introduced; dedicated near-boundary and near-tangent tests remain open.

- [ ] Implement the initial wire/edge data model before concrete analytic surfaces.
	- Start with oriented closed wires made from line-segment edges in a 2D parameter domain.
	- Operations needed by all surfaces: closure validation, orientation/area, point-in-wire classification, edge distance/projection, parametric AABB, duplicate vertex cleanup within tolerance, and visualization sampling.
	- Keep wire classification independent of planar surfaces so cylinders, spheres, cones, and future surfaces reuse the same trimming model.
	- Acceptance test: outer square wire, square-with-hole wires, reversed wire, open wire, and point-on-edge cases.
	- Partial 2026-07-22: line-segment closed wires support cleanup, area, point classification, and edge distance; the dedicated wire fixture tests are still pending.

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

- [ ] Add curve classes for trimmed surface boundaries.
	- Start with line and circle arc boundaries, because these cover many mechanical STEP faces.
	- Required operations: endpoint, tangent, bounding box contribution, point projection, robust 2D winding/crossing support.
	- Acceptance test: planar disk or annulus represented by circular trim curves.

- [ ] Implement cylindrical bounded surfaces.
	- Model: axis frame, radius, parametric `u` angle and `v` height/domain, trim wires in `(u, v)`.
	- Required kernels: analytic ray-cylinder intersection, trim-domain check, normal, distance estimate, AABB, visualization mesh.
	- Acceptance test: exact closed cylinder using two planar caps plus one cylindrical surface, compared against `TGeoTube`.

- [ ] Implement spherical bounded surfaces.
	- Model: center, radius, parametric domain and trim wires.
	- Required kernels: analytic ray-sphere intersection, normal, distance estimate, AABB, visualization mesh.
	- Acceptance test: full sphere or sphere section compared against `TGeoSphere` where applicable.

- [ ] Implement conical bounded surfaces.
	- Model: apex or reference point, axis, opening angle/radii, parametric domain and trim wires.
	- Required kernels: analytic ray-cone intersection, normal, distance estimate, AABB, visualization mesh.
	- Acceptance test: cone or truncated cone compared against an equivalent ROOT shape if available.

- [ ] Decide how to handle torus, Bezier, BSpline, and other high-order CAD surfaces.
	- Conservative default: leave unsupported surfaces on the existing tessellated fallback path.
	- Optional later path: add bounded torus analytically; keep BSpline/Bezier tessellated unless a robust exact kernel is designed.
	- Deliverable: converter support matrix documented in this file and in the script help.

## BVHSurfaceSolid navigation milestones

- [ ] Implement `CloseShape` and BVH construction over bounded-surface AABBs.
	- Reuse the BVH2 pattern from `O2Tessellated::BuildBVH`.
	- Store primitive ids and centers exactly as `O2Tessellated` does, but primitive boxes come from surfaces.
	- Expand AABBs by a documented tolerance to avoid missed boundary hits after float conversion.
	- Validation: unit test proves multiple surfaces are traversed and root bounding box is correct.

- [ ] Implement `Contains` using ray parity over surface intersections.
	- Use a skew test direction, as in `O2Tessellated`, but cluster equal or near-equal `t` intersections so an edge or vertex hit does not double-count.
	- Boundary policy: points within tolerance of any surface should be considered inside unless ROOT convention requires otherwise.
	- Acceptance test: box, cylinder, and points exactly on faces, edges, and vertices.

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

## CAD conversion milestones

- [ ] Add exact-surface extraction probes to `O2_CADtoTGeo.py` without changing default output.
	- Use OpenCascade face/surface adaptors to classify faces by analytic surface type.
	- For each simple shape, collect: face type, orientation, surface parameters, boundary wires/edges, and whether all faces are supported.
	- Deliverable: a `--surface-report` or debug JSON path that reports exact-conversion eligibility per logical volume.
	- Validation: run on one STEP example and confirm existing `geom.C`/facet output is unchanged.

- [ ] Design and document the generated surface sidecar format.
	- It should be versioned and independent of Python object pickling.
	- Candidate: compact binary sidecar similar to `facets_*.bin`, with magic/version, surface count, type enum, orientation, parameters, wire count, edge count, edge curve records, and optional display mesh.
	- Deliverable: format description in this file and loader skeleton in generated C++ macro.
	- Validation: round-trip one planar box sidecar into an `O2BVHSurfaceSolid`.

- [ ] Add generated C++ macro support for `O2BVHSurfaceSolid`.
	- Extend `emit_cpp_prelude` with a `LoadSurfaceSolid(...)` helper and include the new header.
	- Add `emit_surface_solid_cpp(...)` next to `emit_tessellated_cpp(...)`.
	- Keep `emit_tessellated_cpp` as the fallback path.
	- Validation: generated macro compiles and builds a ROOT geometry with exact surface solids for supported volumes.

- [ ] Add converter mode flags.
	- Proposed flags: `--exact-surfaces off|auto|required` and `--surface-report <path>`.
	- `off`: current behavior.
	- `auto`: exact where every face of a logical volume is supported; fallback to tessellated otherwise.
	- `required`: fail with a useful report when any selected logical volume cannot be represented exactly.
	- Validation: command-line help and one smoke run per mode.

- [ ] Implement planar face extraction from OpenCascade.
	- Convert OpenCascade face wires and edges into the planar surface representation.
	- Reject unsupported boundary curves until curve classes are implemented.
	- Validation: generated exact box or box-like STEP uses `O2BVHSurfaceSolid`; unsupported faces still fallback in `auto` mode.

- [ ] Implement cylindrical, spherical, and conical face extraction.
	- Reuse the support matrix from the C++ surface milestones.
	- Preserve face orientation from OpenCascade so normals point outward.
	- Validation: generated cylinder/sphere/cone fixtures navigate like their ROOT primitive equivalents.

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