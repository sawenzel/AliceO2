# Stream F — sidecar v3: edge identity, and closure by topology

Date: 2026-08-01. Branch `swenzel/bvhsurfacesolid`. Brief:
[`Workstreams.md`](Workstreams.md) Stream F; finding
[`CodeReview_Fable_v2.md`](CodeReview_Fable_v2.md) §3/N3.

**Headline.** `measureRimClosure` decided a topological question — do these two faces share an
edge — with a proximity query. The converter already had the answer and threw it away. Sidecar v3
carries it (5 bytes per trim curve), and `validateClosure` now decides by counting edge
identities: **zero tolerance, zero sampling, zero band.** The geometric disagreement between the
two faces becomes a *reported* number, `maxSharedEdgeDeviation`, which is the first defensible
answer this project has had to "how far apart are these two faces really".

Gate: fixtures **8/9 → 9/9**, Bagger **6/12 → 9/12**, oracle disagreements outside tolerance
**0 → 0** on every column of both corpora, `ctest` **64 → 73** cases. Every part on both corpora is
now navigable; the three remaining Bagger failures are the capacity residual (1.3e-06 relative
against a 1e-07 declared tolerance) and nothing else.

---

## 1. Step 1 — the falsifier, run before any implementation

The brief's premise was that the converter's own edge walk sees every edge exactly twice, once
each way. That is a claim about the corpus, not a theorem, so it was measured first. Read-only
pythonOCC probe over the per-part BREPs the gate itself dumps (`<workdir>/db/**/brep_*.brep`, in
cm), plus the ALICE3 STEP read directly.

**Self-check first.** The `box` fixture must read 6 faces, 12 edges, 12 edges with two
occurrences, 0 seam, 0 degenerate, 12 opposite-sense. It does. Nothing below would mean anything
if it did not.

Columns: `=2 / =1 / >2` are edges by number of face occurrences; `seam` is
`BRepTools::IsReallyClosed`; `deg` is `BRep_Tool::Degenerated`; `wOpp / wSame / wSeam / wBad` are
the same census taken through **`BRepTools_WireExplorer`** — the exact walk `_quadric_trim_wire`,
`extract_planar_face` and `_recognized_quadric_wire_block` use — with `wSeam` counting the pairs
whose two occurrences are on *one* face; `walk!` is the number of edges the converter's walk
visits a different number of times than the face explorer does.

| part | faces | edges | =2 | =1 | >2 | seam | deg | wOpp | wSame | wSeam | wBad | walk! | edge tol (cm) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `box` | 6 | 12 | 12 | 0 | 0 | 0 | 0 | 12 | 0 | 0 | 0 | 0 | 1.0e-07 |
| `box_minus_cyl` | 7 | 15 | 15 | 0 | 0 | 1 | 0 | 15 | 0 | 1 | 0 | 0 | 1.0e-07 |
| `box_union_box` | 10 | 20 | 20 | 0 | 0 | 0 | 0 | 20 | 0 | 0 | 0 | 0 | 1.0e-07 |
| `cyl_cross_cyl` | 8 | 15 | 15 | 0 | 0 | 4 | 0 | 15 | 0 | 4 | 0 | 0 | 1.0e-07 |
| `cyl_inter_cyl` | 6 | 9 | 9 | 0 | 0 | 0 | 0 | 9 | 0 | 0 | 0 | 0 | 1.0e-07 |
| `cyl_plus_cone` | 4 | 5 | 5 | 0 | 0 | 2 | 0 | 5 | 0 | 2 | 0 | 0 | 1.0e-07 |
| `sphere_minus_cyl` | 2 | 4 | 4 | 0 | 0 | 2 | 0 | 4 | 0 | 2 | 0 | 0 | 1.0e-07 |
| `torus_union_cyl` | 6 | 9 | 9 | 0 | 0 | 4 | 0 | 9 | 0 | 4 | 0 | 0 | 1.0e-07 |
| `tube_window` | 4 | 8 | 8 | 0 | 0 | 3 | 0 | 8 | 0 | 3 | 0 | 0 | 1.0e-07 |
| `BasePin` | 3 | 3 | 3 | 0 | 0 | 1 | 0 | 3 | 0 | 1 | 0 | 0 | 1.0e-07 |
| `Base` | 44 | 122 | 122 | 0 | 0 | 6 | 0 | 122 | 0 | 6 | 0 | 0 | 1.0e-07 |
| `BoomCylinderInner` | 6 | 9 | 9 | 0 | 0 | 3 | 0 | 9 | 0 | 3 | 0 | 0 | 1.0e-07 |
| `BoomCylinderOuter` | 8 | 12 | 12 | 0 | 0 | 4 | 0 | 12 | 0 | 4 | 0 | 0 | 1.0e-07 |
| `Boom` | 31 | 87 | 87 | 0 | 0 | 7 | 0 | 87 | 0 | 7 | 0 | 0 | 1.0e-07 |
| `BucketCylinderInner` | 6 | 9 | 9 | 0 | 0 | 3 | 0 | 9 | 0 | 3 | 0 | 0 | 1.0e-07 |
| `BucketCylinderOuter` | 10 | 18 | 18 | 0 | 0 | 4 | 0 | 18 | 0 | 4 | 0 | 0 | 1.0e-07 |
| `BucketLink1` | 22 | 54 | 54 | 0 | 0 | 10 | 0 | 54 | 0 | 10 | 0 | 0 | 1.0e-07 |
| `BucketLink2` | 23 | 57 | 57 | 0 | 0 | 11 | 0 | 57 | 0 | 11 | 0 | 0 | 1.0e-07 |
| `StickCylinderInner` | 6 | 9 | 9 | 0 | 0 | 3 | 0 | 9 | 0 | 3 | 0 | 0 | 1.0e-07 |
| `StickCylinderOuter` | 8 | 12 | 12 | 0 | 0 | 4 | 0 | 12 | 0 | 4 | 0 | 0 | 1.0e-07 |
| `Stick` | 24 | 66 | 66 | 0 | 0 | 7 | 0 | 66 | 0 | 7 | 0 | 0 | 1.0e-07 |

**The falsifier passes, and more cleanly than the brief assumed.** 546 edges over 21 parts:

- every edge has **exactly two** face occurrences — not one, not three, on any part;
- every pair carries **opposite** orientation, on the face explorer *and* on the converter's own
  wire walk (`wOpp` equals the edge count for every part; `wSame` and `wBad` are zero everywhere);
- **seams need no special case.** A seam edge is the two visits of one face, still twice and still
  opposite, so "exactly twice, opposite sense" covers it uniformly. 63 of the 546 are seams;
- there is not a single degenerate edge on either corpus;
- the per-edge tolerance is uniformly 1.0e-07 cm — the model tolerance, with no spread at all.

The six failing Bagger cylinders specifically: 9, 12, 9, 18, 9, 12 edges, every one of them
non-degenerate and shared exactly twice with opposite sense. Identity-based closure can decide
them.

### 1.1 The converse trap, checked

An edge shared by exactly two faces in OCCT's map could still correspond to two *different* trim
curves in our sidecar after canonical recognition and splitting. It does not, and here is how that
was verified rather than assumed: the probe runs the **real extractors** (`_FACE_EXTRACTORS` and
`recognize_and_extract_face`, imported from `O2_CADtoTGeo.py`) over every face and compares, per
face, the number of edges `_face_wire_edges` yields against the number of sidecar curves the
record carries.

| | wire faces | param-rect faces | unextractable | curves == edges | curves != edges |
| --- | --- | --- | --- | --- | --- |
| 9 fixtures | 47 | 6 | 0 | 47 | 0 |
| 12 Bagger parts | 181 | 30 | 0 | 181 | 0 |

The correspondence is one-to-one on every face that carries a wire block. Canonical recognition
changes a curve's *type* (a B-spline pcurve recognised as a line or an arc) and never its count,
and no extractor splits an edge.

**The trap it did surface is a different one, and it changed the design.** 36 of the 264 faces
carry no wire block at all: a quadric whose trim is exactly its parametric UV rectangle stores the
trim in its scalar parameters, so there is no per-trim-curve slot to hang an edge id on. The
format therefore carries the edge list **per face**, alongside the wire block rather than inside
it, with an `anchored` bit saying whether entry *i* is also trim curve *i*. The verdict needs only
the identities (it is a count); the deviation measurement needs the anchoring. All six failing
Bagger cylinders and `tube_window` have **zero** param-rect faces, so the parts that mattered are
fully measurable.

### 1.2 ALICE3, spot-checked — two facts that do not hold on Bagger

`CAD_noETA.stp` read flat: **206 `TopoDS_Solid`s, 6 to 2034 faces** (the 55 in
`CSG_Pipeline.md` §6 is a count of leaf logical volumes through XCAF, a different quantity —
noted, not contradicted). Four solids across the size distribution:

| solid | faces | edges | =2 | =1 | >2 | seam | deg | wOpp | wSame | wBad | edge tol (cm) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| #66 | 24 | 60 | 60 | 0 | 0 | 0 | 0 | 60 | 0 | 0 | 4.7e-04 |
| #45 | 45 | 112 | 109 | 3 | 0 | 0 | **3** | 109 | 0 | 0 | 5.2e-04 |
| #109 | 148 | 381 | 381 | 0 | 0 | 0 | 0 | 381 | 0 | 0 | 7.6e-04 |
| #138 | 965 | 2400 | 2388 | 12 | 0 | 0 | **12** | 2388 | 0 | 0 | 4.3e-03 |

The falsifier passes there too, and:

1. **Degenerate edges are real.** Neither Bagger nor the fixtures contain one; ALICE3 does, and
   every `=1` above is exactly a degenerate edge. The `kEdgeDegenerate` flag and its exclusion from
   the incidence count are therefore load-bearing rather than defensive — but they are *untested on
   a converted part*, because no part on either gated corpus has one. First real one should be
   checked by hand.
2. **ALICE3 edge tolerances are 4700x to 43000x Bagger's.** 4.7e-04 to 4.3e-03 cm against a
   uniform 1.0e-07 cm. Any *geometric* closure band on that model would have to be around 1e-3 cm
   — four orders of magnitude wider than the one being argued over on Bagger, and comfortably wide
   enough to swallow real features. This is the strongest argument for the identity criterion that
   the measurement produced, and it was not anticipated: the band was never going to survive the
   next model.

---

## 2. What changed, per layer

### Converter — `scripts/geometry/O2_CADtoTGeo.py`

- `build_edge_table(shape)`: one `TopExp::MapShapes(shape, TopAbs_EDGE, ...)` per solid, giving a
  0-based `uint32` per `TopoDS_Edge`.
- `face_boundary_edge_refs(face, edge_id, anchored)`: the face's `(edgeId, flags)` in
  `_face_wire_edges` order — the same order the three trim extractors emit their curves in. Flags:
  `EDGE_FLAG_REVERSED` from `edge.Orientation()`, `EDGE_FLAG_DEGENERATE` from
  `BRep_Tool.Degenerated`, `EDGE_FLAG_ANCHORED` when the record carries a wire block.
- `extract_surfaces_for_shape` attaches `edge_refs` to every accepted record and returns the edge
  table size. It matters that a single unsupported face rejects the whole solid: an emitted
  sidecar always describes *all* of the shape's faces, so the reader's incidence count is a
  complete statement rather than a statement about the subset that converted.
- `write_surfaces_bin` writes version 3.

### Sidecar format — version 3

Documented in [`BVHSurfaceSolid.md`](BVHSurfaceSolid.md). Two additions, both appended:

```
header:  ... version = 3 ... float64 modelTolerance
         uint32  nModelEdges              # v3 only; size of the solid's edge table
per surface: ... existing record, unchanged, including the wire block ...
         uint32  nBoundaryEdges           # v3 only
         per boundary edge:
           uint32  edgeId                 # index into the model's edge table
           uint8   flags                  # bit0 reversed, bit1 degenerate, bit2 anchored
```

5 bytes per trim curve, as scoped. Nothing at run time.

### Loader — `Detectors/Base/src/O2SurfaceSolidIO.cxx`

- Accepts versions 1-3; parses the header field and the per-surface block.
- `reorderEdgeRefsToKernelOrder`: the sidecar lists wires in the order the CAD walk found them,
  each tagged outer or inner; the `Add*Surface` API takes the outer wire as one argument and the
  inner wires as another. Those orders coincide only when the outer wire happens to come first.
  Permuted, not hoped.
- Rejects an `edgeId` outside the stated edge table.
- **Hardening that was not in the brief but had to be done.** Every count in this format is a
  `uint32` read from the file. A v2 body behind a v3 header yields counts of billions, and
  `resize`-ing to one is not a parse error, it is an out-of-memory kill with no diagnostic at all
  — which is exactly what happened the first time the new reader met the old test fixture. Every
  allocation is now sized against the bytes the file actually holds.

### Kernel — `Detectors/Base/src/BoundedSurface.h`, `O2BVHSurfaceSolid.{h,cxx}`

- `BoundedSurface::BoundaryEdgeRef` + `setBoundaryEdges`/`boundaryEdges`, and a
  `sampleTrimCurve(index)` virtual implemented on all six concrete surfaces (returning false for a
  parametric-rectangle trim, which is the honest answer).
- `SurfaceWire::sourceEdge` / `CurveWire::sourceCurve`: **the trap that would have been silent.**
  Both wire types may reverse a loop in `initialize()` to fix its winding, and `SurfaceWire` may
  drop coincident vertices, so storage index *i* is in general not input index *i*. Pairing curve
  *i* of one face with curve *i* of another after one of them was reversed compares two different
  edges — and still produces a number, and a plausible-looking one. The mapping is recorded and
  is `-1` ("unmeasurable") whenever de-duplication broke the correspondence. Falsified rather than
  asserted: with the naive `return inputIndex`, `TrimCurveIdentitySurvivesWireReorientation`
  reports a deviation of **7.21 cm** instead of 0.
- `applyEdgeIdentityClosure`: applies **only when every** surface states its edges — a face that
  named none would otherwise look like a face with no missing neighbours, which is the exact
  failure being replaced. Rules: twice with opposite sense ⇒ shared; once ⇒ boundary; three or
  more ⇒ non-manifold; twice with the same sense ⇒ reversed; degenerate ⇒ excluded.
- `measureSharedEdgeDeviation`: both faces' realisations of one edge sampled to 33 points and
  compared by symmetric Hausdorff distance, point-to-*polyline* both ways, so the two
  parametrisations need not run the same way and no correspondence between samples is assumed.
  Evaluated on `Curve2D::pointAt`, i.e. on the curves and **not** on their flattened polylines,
  which is why the number is independent of `kBSplineFlatness` (§5).
- `ClosureReport` keeps every geometric measurement it had — `maxRimIsolation`,
  `unmatchedRimLength`, `rimChordResolution`, the per-rim records. They no longer decide anything;
  the four rim *counters* are re-derived from identity so the per-rim reports stay consistent with
  the verdict.
- `GetNavigationReliability()` reads the identity counters directly rather than their per-rim
  projection: a face with no trim loop of its own (a full sphere) contributes no rim record, so a
  defect on its edges would have nowhere to be reported.
- New public accessors: `HasEdgeIdentity`, `GetSourceEdgeCount` and its four siblings,
  `GetMaxSharedEdgeDeviation`, `GetMeasuredSharedEdgeCount`, `GetUnmeasuredSharedEdgeCount`.
  `BVHSurfaceRecord` carries the identity so a streamed solid does not come back deciding closure
  by a different rule than the one that was written.

---

## 3. Before and after

`ctest -R BVHSurfaceSolid`: **64 → 73** cases, green. Nine appended in a delimited block ending
`// --- Stream F: sidecar v3 edge identity ---`.

| | before | after |
| --- | --- | --- |
| fixtures gate | 8/9 (`tube_window` open) | **9/9** |
| Bagger gate | 6/12 | **9/12** |
| **oracle disagreements outside tolerance, both corpora, every column** | **0** | **0** |
| parts failing on navigability | 7 of 21 | **0 of 21** |
| parts failing on capacity | 3 of 21 | 3 of 21 (unchanged) |

The three residual failures are `BoomCylinderOuter`, `BucketCylinderInner` and
`StickCylinderOuter` at 1.31e-06 / 1.39e-06 relative capacity against the model's 1e-07 — the
residual wave 0 left, untouched here. **The gate has moved off navigability entirely.**
`NEXT.md`'s "every part that still fails, fails on navigability and nothing else" is now inverted:
no part fails on navigability, and the only column left is capacity.

### The inertness check, which is the v2 compatibility answer

Both baseline `db` folders (version-2 sidecars) re-scored with the new kernel, `--skip-convert`:
**`gate.json` bit-identical on both corpora** with `timing*`, `*Seconds`, `nsPerCall` and
`closeShapeSeconds*` stripped. A v2 sidecar states no identities, `applyEdgeIdentityClosure`
returns without touching anything, and the geometric verdict stands exactly as before — including
the seven parts it calls open. That is the intended and demonstrated behaviour, and it means an
existing `db` folder and `--skip-convert` keep working and keep their old answers.

No face was dropped: `nFaces` in every `oracle/answers_*.json` equals `nSurfaces` in `gate.json`,
21 of 21.

---

## 4. `maxSharedEdgeDeviation` — the new number

Per part, in cm, with `maxRimIsolation` beside it to show they are not the same quantity.
`meas/unm` is measured / unmeasurable shared edges (unmeasurable = one side is a
parametric-rectangle face).

| part | edges | shared | bnd | n-mf | rev | deg | **maxSharedEdgeDeviation** | maxRimIsolation | meas/unm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `box` | 12 | 12 | 0 | 0 | 0 | 0 | **0** | 0 | 12/0 |
| `box_minus_cyl` | 15 | 15 | 0 | 0 | 0 | 0 | **3.31e-13** | 2.07e-14 | 15/0 |
| `box_union_box` | 20 | 20 | 0 | 0 | 0 | 0 | **0** | 0 | 20/0 |
| `cyl_cross_cyl` | 15 | 15 | 0 | 0 | 0 | 0 | **3.81e-09** | 6.66e-06 | 15/0 |
| `cyl_inter_cyl` | 9 | 9 | 0 | 0 | 0 | 0 | **3.81e-09** | 6.66e-06 | 9/0 |
| `cyl_plus_cone` | 5 | 5 | 0 | 0 | 0 | 0 | **4.14e-13** | 2.78e-13 | 2/3 |
| `sphere_minus_cyl` | 4 | 4 | 0 | 0 | 0 | 0 | **9.02e-13** | 8.32e-13 | 1/3 |
| `torus_union_cyl` | 9 | 9 | 0 | 0 | 0 | 0 | **—** | 3.31e-13 | 0/9 |
| `tube_window` | 8 | 8 | 0 | 0 | 0 | 0 | **4.71e-09** | 2.26e-04 | 8/0 |
| `BasePin` | 3 | 3 | 0 | 0 | 0 | 0 | **—** | 9.49e-16 | 0/3 |
| `Base` | 122 | 122 | 0 | 0 | 0 | 0 | **1.12e-12** | 1.24e-12 | 82/40 |
| `BoomCylinderInner` | 9 | 9 | 0 | 0 | 0 | 0 | **2.46e-08** | 4.80e-05 | 9/0 |
| `BoomCylinderOuter` | 12 | 12 | 0 | 0 | 0 | 0 | **1.80e-07** | 6.92e-05 | 12/0 |
| `Boom` | 87 | 87 | 0 | 0 | 0 | 0 | **9.73e-12** | 9.72e-12 | 59/28 |
| `BucketCylinderInner` | 9 | 9 | 0 | 0 | 0 | 0 | **1.80e-08** | 3.29e-05 | 9/0 |
| `BucketCylinderOuter` | 18 | 18 | 0 | 0 | 0 | 0 | **2.12e-07** | 4.04e-05 | 18/0 |
| `BucketLink1` | 54 | 54 | 0 | 0 | 0 | 0 | **1.77e-11** | 1.63e-03 | 38/16 |
| `BucketLink2` | 57 | 57 | 0 | 0 | 0 | 0 | **4.06e-11** | 4.06e-11 | 46/11 |
| `StickCylinderInner` | 9 | 9 | 0 | 0 | 0 | 0 | **2.46e-08** | 4.80e-05 | 9/0 |
| `StickCylinderOuter` | 12 | 12 | 0 | 0 | 0 | 0 | **1.80e-07** | 6.92e-05 | 12/0 |
| `Stick` | 66 | 66 | 0 | 0 | 0 | 0 | **1.23e-11** | 1.33e-11 | 46/20 |

Read it:

- **The C++ edge counts reproduce the pythonOCC census of §1 exactly**, part for part, over an
  entirely separate code path — converter, format, loader, kernel. That agreement is the strongest
  self-check available and it was the first thing checked.
- **On every part that was already `Reliable`, the deviation is at or below 4.1e-11 cm**, except
  `cyl_cross_cyl` / `cyl_inter_cyl` at 3.8e-09. The gate asked for "below the old match band";
  those two had a band of `1e-8 + 2 x 6.7e-6` ≈ 1.3e-05 cm, so the deviation is 3500x under it,
  and every other part is under its band by many more orders. No previously-`Reliable` verdict
  moved.
- **`maxRimIsolation` is not this number and never was.** On `BucketLink1` they differ by 8 orders
  of magnitude (1.6e-03 against 1.8e-11) — the isolation is how alone the loneliest chord is, and
  on a solid whose faces meet to 18 picometres it is still a millimetre.
- **The six cylinders meet to 1.8e-08 … 2.1e-07 cm.** That is the number the project did not have.
  Note that the two outer cylinders' 1.8e-07 and 2.1e-07 sit *above* the model's own declared
  tolerance of 1.0e-07 cm — by a factor of two. Their faces are further apart than the CAD system
  says its geometry is defined to. That is a small, real and previously invisible defect, and it is
  now a number somebody can act on rather than a verdict somebody has to argue with.
- Three parts (`torus_union_cyl`, `BasePin`) measure nothing: every one of their shared edges has a
  parametric-rectangle face on at least one side. Closing that gap needs the converter to emit a
  wire block for such faces, which changes their geometry and is deliberately not done here.

---

## 5. `CloseShape` no longer depends on `kBSplineFlatness`

The demonstration the brief asked for. Same sidecars, kernel rebuilt with the flattening tolerance
at its shipped 1e-5 and at 1e-8, and both the v2 (geometric verdict) and v3 (identity verdict)
databases scored:

| part | `kBSplineFlatness` | old criterion: open length | old verdict | **maxSharedEdgeDeviation** | new verdict |
| --- | --- | --- | --- | --- | --- |
| `BoomCylinderInner` | 1e-5 | 0.325 cm of 62.8 | open | **2.458e-08** | closed |
| `BoomCylinderInner` | 1e-8 | **0.984 cm** of 62.8 | open | **2.458e-08** | closed |
| `StickCylinderOuter` | 1e-5 | 0.907 cm of 95.8 | open | **1.804e-07** | closed |
| `StickCylinderOuter` | 1e-8 | **3.089 cm** of 95.8 | open | **1.804e-07** | closed |
| `tube_window` | 1e-5 | 0.276 cm of 58.1 | open | **4.712e-09** | closed |
| `tube_window` | 1e-8 | 0.082 cm of 58.1 | open | **4.712e-09** | closed |

The inversion `TolerancePolicy.md` §13.8 recorded is reproduced exactly — 0.325 → 0.984 cm on
`BoomCylinderInner`, and a new instance of it, 0.907 → 3.089 cm on `StickCylinderOuter`, three
times *worse* geometry-for-geometry. The identity verdict does not move, and
`maxSharedEdgeDeviation` is **bit-identical** across the three decades, because it is evaluated on
the curves rather than on their polylines. Tightening the flattening is a pure improvement again.

## 6. The §4.2 direction-independence sweep, re-run

Mandatory here, because seven parts newly become `Reliable` and therefore newly move onto
`Contains`'s single-shot parity fast path — a closure improvement can degrade `Contains`
robustness without touching `Contains`.

Standalone probe (`TolerancePolicy.md` §6), 21 parts x 11000 bbox-spread points x 13 golden-spiral
directions = 231000 points, `Contains` against `ContainsAlongDirection`, same seed before and
after.

| | before (v2, geometric verdict) | after (v3, identity verdict) |
| --- | --- | --- |
| `Reliable` parts | 14 | **21** |
| points | 231000 | 231000 |
| direction-split points | **4** | **4** |
| points where the fixed shot differs from `Contains` | **0** | **0** |

**Identical.** The same four points, on the same two parts — `cyl_inter_cyl` (3) and
`torus_union_cyl` (1) — both of which were `Reliable` before and after. **The seven parts that
newly acquired the fast path contribute zero direction-split points and zero fixed-shot
differences.** The fast path's licence is untouched.

(This probe's seed and spread are its own, so its absolute count is not comparable with §10.1's
"one offender in 143000"; the before/after identity is what it is for.)

---

## 7. What a version-2 sidecar does under the new loader

Exactly what it did before, demonstrated rather than asserted:

- it loads (the reader accepts 1-3);
- it states no identities, so `HasEdgeIdentity()` is false and `applyEdgeIdentityClosure` returns
  immediately;
- the geometric rim verdict decides, unchanged — including calling the seven parts open;
- both baseline gates re-scored against v2 `db` folders are **bit-identical** on every non-timing
  column.

A version-1 sidecar is unaffected in the same way and still gets its documented tolerance
fallback and warning. A **version-3** file may also legitimately state no identities per face
(`nBoundaryEdges = 0`), and is then a v2 file in all but the header; both paths are pinned by
`SidecarV3EdgeIdentityRoundTrip`. Existing `db` folders and `--skip-convert` do not break; they
simply keep the old, weaker verdict until reconverted.

---

## 8. What this criterion does and does not certify — read this before quoting it

It certifies that **the topology of the source B-rep survived conversion**: the faces the sidecar
carries are bounded by edges that pair up exactly, so no face is missing, duplicated or reversed.
That is strictly more than the old criterion established, because the old one thresholded a
quantity that did not answer the question.

It does **not** certify that each face's *geometry* is right. If an extractor mis-trims a face,
identity still says closed. That is why `maxSharedEdgeDeviation` is reported next to the verdict
and why it must be quoted with it. The honest summary of a part is now two numbers, not one:
*closed by topology*, and *the two faces meet to X cm*.

Three consequences worth deciding elsewhere, none of them taken here:

1. **`NavigationReliability` could gate on the deviation** against the model tolerance — the two
  outer cylinders would then be flagged, at 2x their model's declared tolerance. That is a fallback
  *policy* question and Stream E owns it (its step 3, the `auto`-mode fallback).
2. **The harness does not print the new numbers.** `O2SolidHarness` is Stream E's file and was not
  touched, so `gate.json` carries no `maxSharedEdgeDeviation` column. The numbers in §4 come from a
  standalone probe. **Requested harness flag:** `--edge-identity`, emitting
  `HasEdgeIdentity`/`GetSourceEdgeCount` and siblings/`GetMaxSharedEdgeDeviation` into the
  `navigation` block of `--json`. `O2BVHSurfaceSolid::Print()` does carry them.
3. **The degenerate-edge path is untested on a converted part.** ALICE3 has them (§1.2); neither
  gated corpus does. The first converted part with one should be checked by hand.

---

## 9. Files touched

- `scripts/geometry/O2_CADtoTGeo.py` — `build_edge_table`, `face_boundary_edge_refs`,
  `write_surfaces_bin` (v3), `extract_surfaces_for_shape`.
- `scripts/geometry/BVHSurfaceSolid.md` — the version-3 format.
- `Detectors/Base/src/O2SurfaceSolidIO.cxx` — v3 parse, wire-order permutation, allocation bounds.
- `Detectors/Base/src/BoundedSurface.h` — `BoundaryEdgeRef`, `sampleTrimCurve`,
  `sourceEdge`/`sourceCurve`, `applyEdgeIdentityClosure`, `measureSharedEdgeDeviation`,
  `ClosureReport` fields.
- `Detectors/Base/include/DetectorsBase/O2BVHSurfaceSolid.h`,
  `Detectors/Base/src/O2BVHSurfaceSolid.cxx` — the public accessors, `SetSurfaceBoundaryEdges`,
  record persistence, `Print()`, `CloseShape`'s message, `GetNavigationReliability`.
- `Detectors/Base/test/testBVHSurfaceSolid.cxx` — nine cases, appended in the delimited Stream F
  block. Two pre-existing edits were unavoidable: `appendSidecarHeader` gained the v3 header field,
  and `SidecarModelToleranceRoundTrip` asserted that version 3 was an unknown version, which is now
  version 4.

No stream-owned file of any other stream was touched. `NEXT.md`, `CodeReview_Fable*.md`,
`TolerancePolicy.md`, `Workstreams.md`, `CSG_Pipeline.md`, `occtOracle.py`, `runOracleGate.py`,
`make_boolean_fixtures.py`, `O2SolidHarness.*` and `runSolidHarness.cxx` are untouched, and no
`scripts/geometry/csg/**` was created.
