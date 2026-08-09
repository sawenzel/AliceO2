# Stream Y — the wire-join band honours the sidecar's declared model tolerance

Status: **done** (2026-08-09). One loader + kernel change, one pinned unit test
(`StreamY_LoaderHonoursTheDeclaredModelTolerance`), verified on all 20 ALICE3 sidecars.

The largest ALICE3 exact-surface sidecar, `surfaces_ST1829909_01_0_1_1_35.bin` (1052 surfaces,
1.46 MB), was rejected by `LoadSurfaceSolid` with

> surface 1006: wire edge 1 end does not join the next edge start (gap 5.41e-06 cm, tolerance 1e-06 cm)

while the other 19 sidecars of the same conversion loaded. This is the *same* face and the same
class of rejection as TolerancePolicy.md §1.1 (finding S10) — measured correctly this time (the
gap really is 5.41e-6 cm of 3D length, not a mixed rad/cm artifact), and still wrong to reject.

## 1. The evidence

Read straight from the file (python struct walk of the documented format, no C++ involved):

- **The sidecar is version 3 and declares its model tolerance: 4.7230e-4 cm.** That is the source
  CAD assembly's own statement of how closely its edges meet — written by the converter
  (`O2_CADtoTGeo.py`, `write_surfaces_bin`) from the OCC shape, and read by the loader into
  `O2BVHSurfaceSolid::SetModelTolerance` *before* any surface is added.
- **Surface 1006** is a cylinder patch: r = 0.29 cm, phi ∈ [−1.53e-5, 3.14161] rad (half a turn,
  *not* a full-turn seam face), v ∈ [0.109, 0.8] cm, one outer wire of 5 edges
  (bspline, bspline, line, line, line), with v3 edge identities present.
- **The offending join** is edge 1 (a degree-3, 25-pole B-spline in (u, v)) ending at
  (−1.53327e-5, 0.55094420538) against edge 2 (a line) starting at (−3.35886e-6, 0.55094835365):
  Δu = 1.1974e-5 rad, Δv = 4.1483e-6 cm. Through the cylinder's first fundamental form,
  gap = sqrt((r·Δu)² + Δv²) = sqrt((3.472e-6)² + (4.148e-6)²) = **5.41e-6 cm** — exactly the
  number in the error. So the loader's units and metric are right; the rejection is real
  measurement against the wrong band.
- The band it was judged against, `kWireJoinTolerance = 1e-6 cm`, is documented in both
  `BoundedSurface.h` and `O2SurfaceSolidIO.cxx` as the CAD **extractor's** precision — "a
  fallback, not a measurement of the model". This model measures itself at 4.7e-4 cm, and
  5.41e-6 cm is 87× *inside* that statement. The five ST1829909 / large-assembly sidecars all
  declare 4.2–4.7e-4 cm; the small parts declare 1e-8 to 9e-5 cm. Only this file happened to
  carry a join drift above the bare floor.

The kernel already applies exactly this reasoning elsewhere: `validateClosure` matches rims at
`modelTolerance > 0 ? modelTolerance : kRimMatchTolerance` — "a model that declares its own
tolerance gets that instead, because it knows better." The wire-join gate was the one place still
holding a declared-4.7e-4 model to the 1e-6 fallback.

## 2. The decision

**Accept a wire join up to `max(kWireJoinTolerance, declared model tolerance)`** — the declared
value when the model states one looser than the floor, the floor otherwise (no declaration can
make the extractor's endpoint sampling more precise than it is; a declaration of 0 means "not
stated" and changes nothing). This is *not* a widened constant: a v1 sidecar, a v2 file declaring
0, and every hand-built solid keep the 1e-6 band bit-for-bit, and a declared tolerance *below*
the gap still rejects (the model itself calls such a seam open).

The band must hold in **two** places or in neither: the loader's own join check, and the kernel's
`CurveWire::initialize`, which re-checks the same joins moments later inside `Add*Surface`. A
loader-only fix would pass the wire and then watch the kernel return `WireStatus::Open`.

What the residual gap costs afterwards: nothing new. `CurveWire::initialize` has always fixed one
canonical vertex per seam and given it to both curves (K5), so an accepted wire is exactly closed
downstream regardless of which band admitted it; and the rim-closure verdict already used the
declared tolerance before this change.

## 3. What changed

- `Detectors/Base/src/BoundedSurface.h`
  - new `constexpr double wireJoinToleranceFor(double modelTolerance)` beside the constant, with
    the reasoning above;
  - `CurveWire::initialize`, `SurfaceWire::initializeFromEdges`, `buildCurveTrim`, the four
    trimmed-quadric `initialize` overloads and `CurvedPlanarBoundedSurface::initialize` take a
    `joinTolerance = kWireJoinTolerance` parameter and judge joins against it. Defaulted, so
    every existing caller is unchanged.
- `Detectors/Base/src/O2BVHSurfaceSolid.cxx` — `AddCurvedPlanarSurface` and the trimmed
  `AddCylindricalSurface` / `AddConicalSurface` / `AddSphericalSurface` / `AddToroidalSurface`
  pass `wireJoinToleranceFor(fModelTolerance)`. The loader sets the model tolerance from the
  header before the surface loop, so the value is present when it matters.
- `Detectors/Base/src/O2SurfaceSolidIO.cxx` — `wireToCurves` / `collectQuadricTrim` take the same
  band, computed once per file, and the diagnostic now names it and its origin:
  `(gap 5e-06 cm, tolerance 2e-06 cm, declared by the model)` vs `(..., the extractor-precision
  fallback)`.
- `Detectors/Base/test/testBVHSurfaceSolid.cxx` — `StreamY_LoaderHonoursTheDeclaredModelTolerance`
  pins all three behaviours on a synthetic wire-trimmed cylinder with a 5e-6 cm mid-wire join gap:
  loads (and closes, and navigates like the gap-free tube) under a declared 1e-4; rejected as a
  v1 file (floor stands); rejected under a declared 2e-6 (model itself calls it open).

## 4. A caveat the synthetic test measured

The first version of the test put the gap at the trim rectangle's phi-wrap corner of a
**full-turn** cylinder — and the solid loaded but reported 12 cm of open rim. That is not the
join band's doing: a full-turn patch emits its u = 0 / u = 2π seam twice and the rim chaining
cancels the reversed duplicate pair only when the two runs coincide within `kTolerance`. A
within-model-tolerance offset on a seam-adjacent join breaks that exact coincidence, and both
seam runs then count as partnerless rims. The test now keeps the gap mid-wire, where the real
face's gap is (surface 1006 spans half a turn; no full-turn face in this corpus carries an
above-floor join drift — all 20 sidecars close, see §5). If a future model does combine both, the
symptom will be "loads but boundaryRims > 0 on a full-turn face", and the fix belongs in the
seam-pair cancellation band, not here.

## 5. Verification (all on this branch's build, 2026-08-09)

- `ctest -R BVHSurfaceSolid`: **113 of 113 cases, 323406 of 323406 assertions pass** (112
  pre-existing + the new one).
- `surfaces_ST1829909_01_0_1_1_35.bin` now loads: 1052 surfaces, model tolerance 4.723e-4 cm.
  `CloseShape` → **reliability `reliable`**, IsClosed = IsOrientationConsistent = IsNavigable =
  true, 1218/1218 rims matched.
- Standalone probe, 500 random rays/points over 1.5× the AABB: `DistFromOutside` ==
  `DistFromOutside_Loop` and `DistFromInside` == `DistFromInside_Loop` **exactly** on every ray
  (0 mismatches), `Safety` == `Safety_Loop` **exactly** on every point (0 mismatches).
- The same probe over **all 20** ALICE3 sidecars: every file loads, every solid `reliable`,
  0 Dist/Safety mismatches everywhere — the 19 previously-loading files are unaffected (their
  joins were already under the floor; the band only ever widens acceptance, never narrows it).
- One `Contains` vs `Contains_Loop` disagreement was observed on ST1829909_01, **at a point
  outside the solid's AABB** (z = −13.34 vs half-z 8.92, Safety 7.5 cm): `Contains()` answers
  false there by its documented box fast path (correctly — the point cannot be inside);
  `Contains_Loop`, which has no box gate, mis-answers on that one parity ray. The crossing dump
  shows **identical** BVH and Loop hit lists (so the BVH-vs-Loop machinery agrees); the list
  itself has odd parity because the ray double-enters two faces 2.0e-4 cm apart — a face-face
  overlap at the scale of the model's own 4.7e-4 cm tolerance, interior to both trims, so no
  `onTrimBoundary` re-shoot fires. That is a property of the CAD model, not of this change, and
  outside the box it is unreachable through the public `Contains()`.
