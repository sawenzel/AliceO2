# The tessellated fallback is not automatically sound — a note for reference

Recorded 2026-08-02 at the user's request, before the representation benchmark runs. This file is a
**standing note, not a work item and not a design**. It exists because the project has just settled
on a policy — *exact, or tessellated* — whose second half is less safe than it sounds, and the
evidence for that was scattered across three stream documents.

## The policy this note qualifies

Everything that cannot be represented exactly ships as `O2Tessellated`. That is the right call and
it is not in question here. But it makes the mesh a **load-bearing** representation rather than a
convenience, and a mesh is the one representation on this branch that can be *invalid* rather than
merely *inaccurate*. An inaccurate representation is off by a distance you can quote. An unclosed
mesh loses a wall, and a particle walks through it.

## What is actually measured, and on which parts

Two different Bagger parts get conflated when this is recalled from memory. They are not the same
part and only one of them ships as a mesh:

| part | what it ships as | state of its mesh |
| --- | --- | --- |
| `Bagger/Bucket` | **mesh** — 97 faces, 4 spherical + 2 toroidal, no exact sidecar exists | not reported defective |
| `Bagger/BucketLink2` | **exact surfaces** (oracle 0/0/0/0) | **`meshClosedBody = false`** — the only part in either corpus whose mesh is not a closed body |

So the broken mesh is on a part whose *shipped* representation is exact. The defect costs us nothing
today. It tells us what the fallback would cost if that part had needed it.

What `BucketLink2`'s mesh costs, measured (`Stream_J_XRay.md` §2–3, `Stream_G_AnyShape.md` §4 item 5):

- **1430 LOST boundary crossings** — every one of Bagger's — on **2605 of 6912 rays**. A LOST
  crossing is a wall a transported particle walks straight through.
- **366 `contains` disagreements**, worst deviation **2.2 cm**.
- **534 parity events**: the same `O2Tessellated` object's `Contains()` disagrees with its own
  `DistFrom*` about which side of the surface a point is on. Zero on the surface solid, zero on the
  CSG shapes. `parityNB` is 0, so none of these is excused by near-tangency.

And it is not confined to one awkward Bagger part. On the ladder fixtures the meshes of
`cyl_cross_cyl` (16 LOST) and `cyl_inter_cyl` (24 LOST) also lose walls, and both produce
*unterminated* transports. Worse, they are **unstable**: two runs whose ray origins differed by at
most **1.8e-15 cm** (≈ 8 double ulp) moved the fixture meshes from **14 to 40** LOST crossings, while
`surface` and `shape` were bit-identical across the same two runs. Bagger's mesh did not move.

The pattern in all three cases is the same: **curved-surface intersections**. That is where the
tessellator's two patches are meshed independently and do not agree on their shared edge.

## Why this matters more than it looks

1. It is a **correctness** failure, not an accuracy one, so the "quantify the tessellation error in
   cm" number that the representation benchmark will produce does **not** cover it. A mesh can have a
   perfectly respectable chordal deviation and still be open.
2. `meshClosedBody` is already computed and already says `false` on the one bad part — so the
   detector exists. What does not exist is any **policy** attached to it. Nothing currently refuses to
   ship an unclosed mesh, and nothing repairs one.
3. It interacts with the `auto` cascade. If a part falls through CSG and surfaces to mesh, that is
   precisely the part most likely to have the curved intersections that break the mesh — the fallback
   is worst exactly where it is most needed.

## The healing family, for whoever picks this up

The user's instinct — foresee a healing stage — is right. Rough ordering by cost, cheapest first,
and none of these has been tried:

- **Sew before meshing, in OCCT, which we already have.** `BRepBuilderAPI_Sewing` and the
  `ShapeFix_*` family operate on the B-rep *before* triangulation, so they fix the cause (faces that
  do not share an edge) rather than the symptom (triangles that do not share a vertex). This is the
  first thing to try and it costs no new dependency.
- **Mesh both sides of a shared edge from the same discretisation.** OCCT's incremental mesher can
  be driven so that adjacent faces reuse the edge polygon. If the 8-ulp instability is what it looks
  like, this is the direct fix for it.
- **Repair the triangle soup.** CGAL's `Polygon_mesh_processing` (`stitch_borders`,
  `orient_polygon_soup`, `polygon_soup_to_polygon_mesh`, `remove_self_intersections`) or MeshLab /
  libigl equivalents. Powerful and well-tested, but a new external dependency in O2 is a real cost
  and should not be paid before the OCCT-side options are shown to be insufficient.
- **Refuse rather than repair.** The cheapest honest option: make `meshClosedBody = false` a hard
  failure of the converter, so an unclosed mesh is never shipped and the part is reported as
  unrepresentable instead of silently leaky. This is not a substitute for healing but it is a
  substitute for *not knowing*, and it can land in an afternoon.

Whichever is chosen, the acceptance test already exists and should be the one used: the X-ray
benchmark's LOST / unterminated / parity counters on the healed mesh, plus a re-run at an 8-ulp ray
perturbation to show the instability is gone and not merely re-rolled.

## What this note does not claim

- Not that `Bagger/Bucket` — the part that actually ships as a mesh — is defective. It is not
  reported as such. Nobody has specifically checked whether its curved faces produce the same
  shared-edge problem, and that check is worth doing before any healing work starts, because it
  decides whether this is a live defect in shipped output or a latent one in the fallback.

  **Addendum, same day, from the implicit-trim census (`Stream_O_ImplicitTrims.md`): `Bucket` may
  not stay a mesh part at all.** The census reports that `Bucket` has **zero** faces without an
  analytic surface, and that it falls back to mesh for an unrelated reason — `extract_planar_face`
  declines an *ellipse* boundary, on 4 edges whose deviation from the two implicit surfaces is
  ≤ 3e-11 cm. If that holds, `Bucket` becomes exact and Bagger ships **no** mesh part at all. The
  "4 spherical + 2 toroidal faces" explanation in `Tutorial.md` §6 and in the table above is then
  wrong about the *cause*, though still right that no sidecar exists today. Treat the table's
  `Bucket` row as provisional until the ellipse trim is either fixed or shown not to fix it. This
  does not weaken the note: it removes the one part that would have exercised the fallback in
  shipped output, which makes the fallback *less* tested, not more trustworthy.
- Not that CGAL is the answer. It is the best-known name in the family, which is why the user named
  it; it is also the most expensive to adopt.
- No estimate of effort. Nothing here has been prototyped.
