"""Split a part into single cells, so a union of them can be emitted.

This is `probes/cellCountProbe.py`'s loop, promoted to production rather than rewritten: the
probe measured it at 15 of 16 complete on ALICE3 and on every CSG-declining Bagger part, in
under a second each, with volume conserved to 1.3e-6 relative (`Stream_AA_FlatCSG.md` §2). The
probe now imports these functions, so the same code carries the probe's own seven closed-form
self-checks -- a box is one cell, an L-plate is two, a grooved block is three, a blind bore is
two, a through-hole plate is *one* -- and there is exactly one implementation to be wrong.

    while a piece has a trusted concave (or mixed) edge:
        extend the carrier of one of the edge's two faces to a full surface;
        split the piece with BRepAlgoAPI_Splitter;
        recurse on the pieces.
    a piece with no trusted concave edge is one CSG cell.

Connectivity comes first, and that is not cosmetic
--------------------------------------------------
R3 §6.1 measured the lesson the hard way: **zero trusted concave edges does not mean one cell.**
ALICE3's `ST0923290_01#b12` is four pins on four parallel axes; it has no concave edge because
its pieces share no edge at all, and read as one cell its halfspaces intersect to nothing. So
the loop starts from the shape's connected solids, which costs no boolean at all, and only then
splits on concave edges. Every split afterwards re-enumerates solids for the same reason.

What this module does NOT do is decide anything
-----------------------------------------------
It returns a report. `recognise._match_union_of_cells` turns a `stop` into a readable decline and
`csg/accept.py` still decides whether the emitted union is the part. Two failure modes are known
and are reported rather than hidden: OCCT's splitter is a tolerant boolean core and not an exact
arrangement, so the summed piece volumes can drift from the original (measured at 1.4e-5 and
1.1e-4 relative on two ALICE3 parts, `Stream_AA` §3.2) -- `volumeConserved` says so and the
caller declines; and splitting only on *trusted* concave edges leaves blend-seam reflexes inside
cells, which nothing here can see and only the part-level symmetric difference catches.

Carriers a CAD exporter wrote as B-splines
------------------------------------------
Every splitter failure the probe ever observed had one cause: the witness edge's faces were
NURBS-encoded quadrics it could not extend into a cutting tool (`Stream_AA` §2.4). That is what
Tier 0 was built for, and `carrier_tool_face` now asks `csg/tier0.py` whenever the adaptor has no
analytic surface of its own, so the tool is built from the canonical parameters.
"""

import math
import time

# The per-part cell budget. `Stream_AA_FlatCSG.md` §2.3 measured that a face budget is the wrong
# gate -- `Bucket` at 97 faces and 68 concave edges decomposes in under half a second -- and
# proposed a cell budget of about 64 instead. A part over it is declined naming the bound; it is
# never shipped as a tree that wide.
PART_MAX_CELLS = 64

# The wall clock and the split count, so a blow-up is a decline and never a hang. Both are the
# probe's, scaled down: a production conversion runs over a thousand parts and cannot spend the
# probe's 300 s on one of them.
MAX_SPLITS = 256
TIMEOUT_S = 60.0

# The splitter's volume guard: the pieces must still sum to the part. Relative, because the
# quantity is a volume ratio; 1e-6 is the probe's own bound, met by every part it completed on
# Bagger and by 15 of the 16 ALICE3 targets, and exceeded by exactly the two parts `Stream_AA`
# §3.2 names. A breach declines -- it is never accepted and never loosened.
VOLUME_REL_TOL = 1.0e-6

# How far a cutting tool has to reach, in bounding-box diagonals of the piece being split.
TOOL_EXTENT_DIAGONALS = 4.0


def _occ():
    from csg.occ_env import ensure_occ
    ensure_occ()


# ------------------------------------------------------------------------------------------
# connectivity
# ------------------------------------------------------------------------------------------

def solid_components(shape):
    """The shape's own `TopoDS_Solid` bodies, or the shape itself when it carries none.

    A CAD leaf that is a compound of disjoint bodies arrives here as exactly that, and separating
    them is free -- no boolean, no tolerance, no witness edge. `recognise` used to read such a
    shape as one cell and get an empty intersection out of it (R3 §6.1).
    """
    from OCC.Core.TopAbs import TopAbs_SOLID
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
    out = []
    walk = TopExp_Explorer(shape, TopAbs_SOLID)
    while walk.More():
        out.append(topods.Solid(walk.Current()))
        walk.Next()
    return out or [shape]


# ------------------------------------------------------------------------------------------
# finding the split witness: the sharpest trusted concave/mixed edge, with its two faces
# ------------------------------------------------------------------------------------------

def first_trusted_concave_edge(solid):
    """(edge, face1, face2) of a trusted concave or mixed dihedral, or None.

    "Trusted" is `Stream_A_CSG.md` §1.4's filter: a concave/mixed verdict whose |n1 x n2| stays
    below `NEAR_TANGENTIAL_SIN` is a blend seam at the noise floor and no decomposition should
    split on it. Preference order: the sharpest (largest max-sin) witness, so the split happens
    where the reflex is best conditioned.
    """
    from csg.census import NEAR_TANGENTIAL_SIN, FaceEdgeOrientations, edge_dihedral
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE
    from OCC.Core.TopExp import topexp
    from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCC.Core.TopoDS import topods

    amap = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, amap)
    orients = FaceEdgeOrientations(solid)
    best = None
    for i in range(1, amap.Size() + 1):
        edge = topods.Edge(amap.FindKey(i))
        if BRep_Tool.Degenerated(edge):
            continue
        faces = list(amap.FindFromIndex(i))
        distinct = []
        for f in faces:
            if not any(f.IsSame(g) for g in distinct):
                distinct.append(f)
        if len(distinct) != 2:
            continue
        f1, f2 = topods.Face(distinct[0]), topods.Face(distinct[1])
        verdict, max_sin = edge_dihedral(edge, f1, f2, orients)
        if verdict in ("concave", "mixed") and max_sin >= NEAR_TANGENTIAL_SIN:
            if best is None or max_sin > best[3]:
                best = (edge, f1, f2, max_sin)
    return None if best is None else best[:3]


def count_trusted_concave(solid):
    from csg import census
    counts = census.edge_census(solid)
    return (counts["concave"] + counts["mixed"]
            - counts["concaveNearTangential"] - counts["mixedNearTangential"])


# ------------------------------------------------------------------------------------------
# extending a face's carrier into a splitting tool
# ------------------------------------------------------------------------------------------

def carrier_tool_face(face, extent, scale=None):
    """A face covering the whole carrier of `face`, big enough to cut anything within `extent`.

    For a plane and the lateral quadrics the parametric window is widened to `extent` (the caller
    passes a few bounding-box diagonals); a sphere and a torus are already closed. A face whose
    stored surface is not analytic goes through `csg/tier0.py` and the tool is built from the
    canonical parameters -- which is the whole reason R3 came first, since every splitter failure
    the probe ever recorded was a NURBS-encoded quadric it could not extend. Returns None when
    the carrier is genuinely not one of the five families.
    """
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.GeomAbs import (GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane,
                                  GeomAbs_Sphere, GeomAbs_Torus)
    ad = BRepAdaptor_Surface(face, True)
    kind = ad.GetType()
    try:
        if kind == GeomAbs_Plane:
            return _tool_from_plane(ad.Plane(), extent)
        if kind == GeomAbs_Cylinder:
            return _tool_from_cylinder(ad.Cylinder(), extent)
        if kind == GeomAbs_Cone:
            return _tool_from_cone(ad.Cone(), extent)
        if kind == GeomAbs_Sphere:
            return _tool_from_sphere(ad.Sphere())
        if kind == GeomAbs_Torus:
            return _tool_from_torus(ad.Torus())
    except Exception:                                            # noqa: BLE001
        return None
    if scale is None:
        return None
    return _canonical_tool_face(face, ad, extent, scale)


def _canonical_tool_face(face, adaptor, extent, scale):
    """The tool a Tier-0 canonicalised face extends to, or None if the face is not canonical."""
    from csg import tier0
    carrier, _gap = tier0.canonicalise(face, adaptor, scale)
    if carrier is None:
        return None
    from OCC.Core.gp import (gp_Ax3, gp_Cone, gp_Cylinder, gp_Dir, gp_Pln, gp_Pnt, gp_Sphere,
                             gp_Torus)
    try:
        if carrier["kind"] == "plane":
            return _tool_from_plane(
                gp_Pln(gp_Pnt(*carrier["p"]), gp_Dir(*carrier["n"])), extent)
        frame = gp_Ax3(gp_Pnt(*carrier["p"]), gp_Dir(*carrier["d"]), gp_Dir(*carrier["x"]))
        if carrier["kind"] == "cylinder":
            return _tool_from_cylinder(gp_Cylinder(frame, carrier["r"]), extent)
        if carrier["kind"] == "cone":
            return _tool_from_cone(gp_Cone(frame, carrier["a"], carrier["r"]), extent)
        if carrier["kind"] == "sphere":
            return _tool_from_sphere(gp_Sphere(
                gp_Ax3(gp_Pnt(*carrier["p"]), gp_Dir(0.0, 0.0, 1.0)), carrier["r"]))
        if carrier["kind"] == "torus":
            return _tool_from_torus(gp_Torus(frame, carrier["r"], carrier["rt"]))
    except Exception:                                            # noqa: BLE001
        return None
    return None


def _face_of(surface, umin, umax, vmin, vmax):
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    return BRepBuilderAPI_MakeFace(surface, umin, umax, vmin, vmax, 1.0e-7).Face()


def _tool_from_plane(pln, extent):
    from OCC.Core.Geom import Geom_Plane
    return _face_of(Geom_Plane(pln), -extent, extent, -extent, extent)


def _tool_from_cylinder(cyl, extent):
    from OCC.Core.Geom import Geom_CylindricalSurface
    return _face_of(Geom_CylindricalSurface(cyl), 0.0, 2.0 * math.pi, -extent, extent)


def _tool_from_cone(cone, extent):
    from OCC.Core.Geom import Geom_ConicalSurface
    return _face_of(Geom_ConicalSurface(cone), 0.0, 2.0 * math.pi, -extent, extent)


def _tool_from_sphere(sphere):
    from OCC.Core.Geom import Geom_SphericalSurface
    return _face_of(Geom_SphericalSurface(sphere), 0.0, 2.0 * math.pi,
                    -0.5 * math.pi, 0.5 * math.pi)


def _tool_from_torus(torus):
    from OCC.Core.Geom import Geom_ToroidalSurface
    return _face_of(Geom_ToroidalSurface(torus), 0.0, 2.0 * math.pi, 0.0, 2.0 * math.pi)


def split_solid(piece, tool):
    """Split `piece` by `tool`; returns the list of solids, or None on failure."""
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Splitter
    from OCC.Core.TopAbs import TopAbs_SOLID
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopTools import TopTools_ListOfShape
    from OCC.Core.TopoDS import topods

    splitter = BRepAlgoAPI_Splitter()
    args = TopTools_ListOfShape()
    args.Append(piece)
    tools = TopTools_ListOfShape()
    tools.Append(tool)
    splitter.SetArguments(args)
    splitter.SetTools(tools)
    try:
        splitter.Build()
    except Exception:                                            # noqa: BLE001
        return None
    if not splitter.IsDone():
        return None
    out = []
    walk = TopExp_Explorer(splitter.Shape(), TopAbs_SOLID)
    while walk.More():
        out.append(topods.Solid(walk.Current()))
        walk.Next()
    return out or None


# ------------------------------------------------------------------------------------------
# the loop
# ------------------------------------------------------------------------------------------

def bbox_diagonal(shape):
    from csg import census
    box = census.bounding_box(shape)
    if box is None:
        return 1.0
    xmin, ymin, zmin, xmax, ymax, zmax = box
    return math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)


def split_into_cells(solid, max_cells=PART_MAX_CELLS, max_splits=MAX_SPLITS,
                     timeout_s=TIMEOUT_S, scale=None, verbose=False):
    """Split at trusted concave edges until every piece is one cell.

    Returns a report. `stop` is None when the decomposition finished within every budget and
    every piece reached a cell; otherwise it names the bound that was hit, in words a decline can
    be built from. `pieces` are the terminal single-cell solids, `unresolved` the pieces whose
    witness edge could not be cut by either of its faces' carriers.

    `scale` is the part's `max(bounding-box diagonal, 1 cm)`, handed to Tier 0 when a witness
    face has to be canonicalised before it can be extended; None turns that off.
    """
    _occ()
    from csg.census import volume_of
    start = time.time()
    diagonal = bbox_diagonal(solid)
    extent = TOOL_EXTENT_DIAGONALS * max(diagonal, 1.0)
    original_volume = volume_of(solid)
    pending = list(solid_components(solid))
    n_components = len(pending)
    cells, unresolved = [], []
    n_splits = n_split_failures = 0
    stop = None

    while pending:
        if len(cells) + len(pending) + len(unresolved) > max_cells:
            stop = f"the cell budget of {max_cells} was exceeded"
            break
        if n_splits >= max_splits:
            stop = f"the split budget of {max_splits} was exceeded"
            break
        if time.time() - start > timeout_s:
            stop = f"the {timeout_s:.0f} s decomposition timeout was exceeded"
            break
        piece = pending.pop()
        witness = first_trusted_concave_edge(piece)
        if witness is None:
            cells.append(piece)
            continue
        _edge, face1, face2 = witness
        parts = _split_at(piece, face1, face2, extent, scale)
        if parts is None:
            n_split_failures += 1
            unresolved.append(piece)
            continue
        n_splits += 1
        pending.extend(parts)
        if verbose:
            print(f"    split {n_splits}: {len(parts)} piece(s), {len(cells)} cell(s) so far, "
                  f"{len(pending)} pending")

    piece_volume = sum(volume_of(c) for c in cells) + sum(volume_of(u) for u in unresolved)
    conserved = None
    if stop is None:
        conserved = abs(piece_volume - original_volume) <= \
            VOLUME_REL_TOL * max(abs(original_volume), 1.0)
    return {
        "pieces": cells,
        "unresolved": unresolved,
        "pending": list(pending),
        "components": n_components,
        "splits": n_splits,
        "splitFailures": n_split_failures,
        "stop": stop,
        "volumeOriginal": original_volume,
        "volumePieces": piece_volume,
        "volumeConserved": conserved,
        "volumeDrift": (abs(piece_volume - original_volume) / max(abs(original_volume), 1e-30)
                        if stop is None else None),
        "seconds": round(time.time() - start, 2),
        "diagonal": diagonal,
    }


def _split_at(piece, face1, face2, extent, scale):
    """Cut `piece` on one of the witness edge's two carriers; the pieces, or None.

    The planar carrier is preferred as the knife because OCCT's boolean core is most robust on
    planes; the other face is the fallback, and a cut that returns a single solid counts as a
    failure of that tool rather than as a result.
    """
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.GeomAbs import GeomAbs_Plane
    faces = sorted((face1, face2),
                   key=lambda f: BRepAdaptor_Surface(f, True).GetType() != GeomAbs_Plane)
    for face in faces:
        tool = carrier_tool_face(face, extent, scale)
        if tool is None:
            continue
        parts = split_solid(piece, tool)
        if parts is not None and len(parts) > 1:
            return parts
    return None
