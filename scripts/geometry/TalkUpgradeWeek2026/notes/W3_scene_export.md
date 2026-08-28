# W3 step 1 — the scene exporter

2026-08-28. `export_scene.py` turns a built CAD-converted TGeo geometry into a Blender-ready
bundle in which every placed part carries its identity and its representation. It makes no look
decisions: no colours, no materials, no cameras. Step 2 (Blender) is a separate session.

## What replaced what

`O2DPGAgentic/doc/detector-cad-rendering.md` had to recover placements out of the STEP file by
matching solids on volume and clustering by radius, because the mesh dump's repeated leaf volumes
are unplaced by design. **None of that is needed any more.** We build the converted geometry and
read the placements straight off the node tree. The join is:

    facets_<part>.bin   local triangles (cm)
  x built TGeo geometry global placement per node
  x csg_report.json     representation, recogniser, description

## Format: OBJ, one object per placement

Chosen over glTF and PLY because it is the only one of the three I could **verify without opening
Blender** — a text format I can re-parse and check face-index ranges, object names and bounding
boxes against the manifest, which is exactly what caught the defect below.

The manifest also carries, per placement, the absolute path of the source `.bin` and the corrected
3x4 world matrix. That is the fast route for IRIS: the existing `render_its.py` already reads this
`.bin` format directly, so a Blender script can build meshes itself and skip a 657 MB OBJ import.
Both routes describe the same scene.

## The identity problem, and how it is solved

A part's lid lives only in the C++ variable name in `geom.C` (`vol__0_1_1_2`); the TGeo volume
carries the display name alone, and display names repeat across the bodies of one multi-body
label — IRIS has **28** volumes called `ST0923290_01`. So the exporter parses `geom.C` for the
volume creation order, asserts that the built geometry's non-assembly volumes appear in exactly
that order and with those names, and then **renames each volume to its part key**. After that a
placement identifies its part unambiguously.

## The defect this found: a CSG part's facets were misplaced

The first bundle passed every count and still had every CSG part in the wrong place.

`geom.C` composes the shape placement into the node matrix (the `_placed` matrices:
`tr->Multiply(shapePlace__<lid>)`). So a CSG node's matrix maps the **canonical shape** frame to
the world. The facet file is in the part's **CAD-local** frame. Applying the node matrix to the
facets therefore applies `shapePlacement` twice. Surface and mesh parts have no shape placement,
which is why only the scene extent gave it away: IRIS came out `dx = 24.65` against the geometry's
own `dx = 11.68`.

Fixed by `world = M * shapePlacement^-1 * cadLocal`. Negative control on Bagger, per part:

| part | representation | excursion uncorrected | corrected |
| --- | --- | --- | --- |
| `BasePin_0_1_1_2` | csg | **5.916 cm** | 0.0000 |
| `Base_0_1_1_3` | surface | 0.000 | 0.0000 |
| the other 11 | both | < 0.07 (inside) | unchanged |

5.916 cm is exactly `shapePlacement`'s own translation for that part. Only `BasePin` has a
non-identity shape placement in Bagger, which is why the rest are unmoved.

**The lesson for the gates:** part counts cannot see a misplacement. What sees it is comparing the
triangles against the box of the shape TGeo actually navigates, and that check is now in the
exporter.

## Verification, both models

| gate | Bagger | IRIS |
| --- | --- | --- |
| representation counts vs `csg_report` | 7 csg / 6 surface / 0 mesh | 6 csg / 26 surface / 29 mesh |
| every part exactly once | 13 / 13 | 61 / 61 |
| facets inside the shape TGeo navigates | worst 0.0000 cm | worst **0.0007 cm** (`ST2500394_01_b1`) |
| scene inside the built geometry's own box | worst -0.0127 cm (inside) | worst **0.0000 cm** |
| triangles | 43 822 local = 43 822 world over 13 placements | 1 398 571 local, **6 001 720** world over 125 placements |
| scene half-lengths | 5.17, 33.80, 28.35 | 10.81, 497.96, 21.18 |
| top volume box (half) | 5.19, 37.05, 36.08 | 11.68, 498.55, 21.57 |

The scene is *inside* the top box rather than equal to it, which is the correct relation: an
assembly's box is the union of its daughters' shape boxes and those can be loose. IRIS's x minimum
(-12.928) touches the box exactly.

Re-parsing the written OBJ independently: Bagger 13 objects / 131 466 vertices / 43 822 faces, max
face index equal to the vertex count; IRIS 125 objects / 18 005 160 vertices / 6 001 720 faces, no
duplicate object names, bounding boxes agreeing with the manifest to 1e-2 cm.

## Label anchor

The centroid of the part's **largest-area triangle**, with that triangle's normal flipped to point
away from the part centroid. The biggest facet is a broad flat face, which is where a leader line
reads best. Verified on the surface of its own part for 13/13 Bagger and 61/61 IRIS parts, normals
unit length.

## `labelShort`

Sandro's call: dV_sym is detail the slide does not need. So each part carries a short form —
`<volume> · <class>` — beside the full `recogniser` and `description`, which stay for a caption:

- csg → the shape class, e.g. `ST0923290_01 · TGeoBBox`, `BasePin · TGeoTube`
- surface → `Base · exact surfaces`
- mesh → `ST0923290_01 · tessellated`

## What step 2 needs to know

- **Label ambiguity is real.** `ST0923290_01` is 28 separate bodies and they do not share a
  representation — the same volume name appears as csg, surface *and* mesh. `labelShort` alone
  cannot disambiguate them; the unique key is `part` (name + lid). Label a representative body, or
  label the group, but do not draw 28 identical labels.
- **Placement multiplicity** in IRIS: 16 parts placed once, 29 twice, 15 three times, 1 six times.
  Colour rotates per part (Sandro's call), so all placements of one part share its colour.
- The IRIS OBJ is 657 MB. Prefer the manifest's `facetFile` + `matrix3x4` route in Blender.
- ROOT warns `TGeoVoxelFinder::SortAll: Wrong bounding box for volume ST2487730_01` when IRIS is
  closed. That is one of the 14 assemblies, not a part, and it is the loose-bounding-box item
  already open in `NEXT.md`. It does not affect the export, which uses shape boxes per node.
- Units are cm and the frame is the built geometry's top volume — the same frame the earlier TB
  slide-20 pictures were taken in.

## Bundles

`scratchpad/scene/bagger` (4.4 MB) and `scratchpad/scene/iris` (627 MB). Bulk, not committed.
