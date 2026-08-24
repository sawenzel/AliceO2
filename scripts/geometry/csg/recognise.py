"""Bagger-sufficient CSG recognition: whole-part primitives, and the two-cluster tube union.

Scope, deliberately narrow (`Stream_A_CSG.md` §3, `Tutorial.md` §6 item 2):

  * **Tier 1** — the part's whole face set is one ROOT primitive: `TGeoBBox`, `TGeoTube`,
    `TGeoTubeSeg`, `TGeoCone`, `TGeoSphere`;
  * **Tier 2, the one measured case** — two coaxial *clusters* of cylinders, i.e. a barrel and a
    lug on a second axis, emitted as `TGeoTube u TGeoTube`. The census measured all six Bagger
    "ram" parts as exactly this, and measured them as a **union**, not the difference the older
    plan text names: the bores are `rmin` on the leaves, so no subtraction node is needed.

Two whole-part families were added after that census, each strictly after the tier-1 matcher for
the same face set declines, so nothing recognised before them changed: the **revolved profile**
(one axis, any number of z sections, a `TGeoPcon`) and the **prism family** (an all-planar face
graph read as a stack of sections: `TGeoTrd1`, `TGeoTrd2`, `TGeoArb8`, `TGeoXtru`, `TGeoPgon`).

This is *not* the general Tier-3 decomposition, which the census says Bagger does not need and
which the project has gated on a cell-count table that does not exist yet.

Proposal is cheap, acceptance is exact
--------------------------------------
Nothing here proves anything. Every function below *proposes* a placed primitive from the carrier
structure, and `accept.symmetric_difference` then decides. That asymmetry is the design
(`CSG_Pipeline.md` §3.5), and it is why this module never inspects a trim curve: a proposal built
from an untrimmed carrier that happens to be wrong is rejected by a volume, not by a heuristic.

What it *must* not do is accept a part it has not actually recognised. Every threshold here is
relative to the part's own bounding-box diagonal and tight (1e-6), and every unhandled structure
returns a reason rather than a guess. A declined part costs coverage; a wrongly accepted one costs
correctness.

Where the extents come from
---------------------------
The radius of a cylinder is in its carrier; its *length* is not — that lives in the trim. This
module reads it from `BRepTools.UVBounds`, which is the exact parametric extent of the trimmed
face, `v` being signed axial distance for a cylinder and `v*cos(semiangle)` for a cone. Using the
face's own extent rather than the whole carrier matters for the lug case: the rod of a Bagger ram
is cut where it meets the eye's outer cylinder, and a rod emitted any longer than that would
intrude into the eye's bore and fill a hole that must stay open.
"""

import math

from csg import primitives as prim
from csg.primitives import _add, _cross, _dot, _norm, _scale, _sub, _unit

# Relative tolerance on directions, radii and offsets, scaled by the part's bounding-box
# diagonal. Same value and the same reasoning as the census: CAD arrives with ~1e-7 relative
# agreement between faces the engineer meant to be identical.
REL_TOL = 1.0e-6
ANG_TOL = 1.0e-6


class Declined(Exception):
    """Raised internally with the reason; recognise() turns it into a report entry."""


# ------------------------------------------------------------------------------------------
# face analysis
# ------------------------------------------------------------------------------------------

def _face_records(solid):
    """[{kind, ...carrier..., uv bounds}] for every face, or a reason why the solid is out."""
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.BRepTools import breptools
    from OCC.Core.GeomAbs import (GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane,
                                  GeomAbs_Sphere, GeomAbs_Torus)
    from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods

    records = []
    n_torus = 0
    n_freeform = 0
    n_faces = 0
    exp = TopExp_Explorer(solid, TopAbs_FACE)
    while exp.More():
        face = topods.Face(exp.Current())
        exp.Next()
        n_faces += 1
        ad = BRepAdaptor_Surface(face, True)
        umin, umax, vmin, vmax = breptools.UVBounds(face)
        t = ad.GetType()
        rec = {"uv": (umin, umax, vmin, vmax), "face": face,
               "reversed": face.Orientation() == TopAbs_REVERSED}
        if t == GeomAbs_Plane:
            pl = ad.Plane()
            n = _xyz(pl.Axis().Direction())
            if rec["reversed"]:
                n = _scale(n, -1.0)
            rec.update(kind="plane", n=n, p=_xyz(pl.Axis().Location()))
        elif t == GeomAbs_Cylinder:
            cy = ad.Cylinder()
            rec.update(kind="cylinder", d=_xyz(cy.Axis().Direction()),
                       p=_xyz(cy.Axis().Location()), x=_xyz(cy.Position().XDirection()),
                       r=cy.Radius())
        elif t == GeomAbs_Cone:
            co = ad.Cone()
            rec.update(kind="cone", d=_xyz(co.Axis().Direction()),
                       p=_xyz(co.Axis().Location()), x=_xyz(co.Position().XDirection()),
                       r=co.RefRadius(), a=co.SemiAngle())
        elif t == GeomAbs_Sphere:
            sp = ad.Sphere()
            rec.update(kind="sphere", p=_xyz(sp.Location()), r=sp.Radius())
        elif t == GeomAbs_Torus:
            # Counted rather than reported at first sight, so the decline says how far out of
            # scope the part is ("2 of 97 faces") instead of naming one face.
            n_torus += 1
            continue
        else:
            n_freeform += 1
            continue
        records.append(rec)
    if n_torus or n_freeform:
        found = []
        if n_freeform:
            found.append(f"free-form faces: {n_freeform} of {n_faces} "
                         "(surface kind outside plane/cylinder/cone/sphere; a twisted "
                         "TGeoArb8 side is one of these and is out of scope)")
        if n_torus:
            found.append(f"toroidal faces: {n_torus} of {n_faces} "
                         "(out of the recogniser's scope)")
        return None, "; ".join(found)
    if not records:
        return None, "no faces"
    return records, None


def _xyz(v):
    return (v.X(), v.Y(), v.Z())


def _bbox_diagonal(shape):
    from OCC.Core.Bnd import Bnd_Box
    from OCC.Core.BRepBndLib import brepbndlib
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)


# ------------------------------------------------------------------------------------------
# direction / axis predicates
# ------------------------------------------------------------------------------------------

def _parallel(a, b):
    return _norm(_cross(a, b)) <= ANG_TOL and _dot(a, b) > 0.0


def _collinear(a, b):
    return _norm(_cross(a, b)) <= ANG_TOL


def _perpendicular(a, b):
    return abs(_dot(a, b)) <= ANG_TOL


_COORDINATE_AXES = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def _snap_to_coordinate_axis(vec):
    """(index, sign) if `vec` is a coordinate axis to within ANG_TOL, else None.

    Worth the few lines because of what it buys downstream: a frame whose rotation is exactly the
    identity lets `primitives.build_root` emit a bare `TGeoBBox` (which carries its own origin)
    instead of a `TGeoCompositeShape`. That is one shape instead of three objects and a matrix,
    an analytic `Capacity()` instead of a Monte-Carlo one, and it is the common case in real
    mechanical CAD -- the census counts 62560 placed six-plane boxes in oTOF alone.
    """
    for index, axis in enumerate(_COORDINATE_AXES):
        dot = _dot(vec, axis)
        if abs(abs(dot) - 1.0) <= ANG_TOL and _norm(_cross(vec, axis)) <= ANG_TOL:
            return index, (1.0 if dot > 0.0 else -1.0)
    return None


def _on_axis(point, loc, direction, tol):
    delta = _sub(point, loc)
    return _norm(_sub(delta, _scale(direction, _dot(delta, direction)))) <= tol


# ------------------------------------------------------------------------------------------
# clustering
# ------------------------------------------------------------------------------------------

def _axial_extent(rec, axis_dir, axis_loc):
    """The face's [tmin, tmax] along the cluster axis, and its radii at those two ends."""
    umin, umax, vmin, vmax = rec["uv"]
    base = _dot(_sub(rec["p"], axis_loc), axis_dir)
    sign = 1.0 if _dot(rec["d"], axis_dir) > 0.0 else -1.0
    if rec["kind"] == "cylinder":
        t0, t1 = base + sign * vmin, base + sign * vmax
        r0 = r1 = rec["r"]
    else:                                                   # cone
        ca, sa = math.cos(rec["a"]), math.sin(rec["a"])
        t0, t1 = base + sign * vmin * ca, base + sign * vmax * ca
        r0, r1 = rec["r"] + vmin * sa, rec["r"] + vmax * sa
    if t0 > t1:
        t0, t1, r0, r1 = t1, t0, r1, r0
    return t0, t1, r0, r1


def _cluster_axial(records, tol):
    """Group cylinder/cone faces by the axis *line* they sit on."""
    clusters = []
    for rec in records:
        if rec["kind"] not in ("cylinder", "cone"):
            continue
        for cl in clusters:
            if _collinear(rec["d"], cl["dir"]) and _on_axis(rec["p"], cl["loc"], cl["dir"], tol):
                cl["members"].append(rec)
                break
        else:
            clusters.append({"dir": rec["d"], "loc": rec["p"], "x": rec["x"], "members": [rec]})
    for cl in clusters:
        spans = [_axial_extent(m, cl["dir"], cl["loc"]) for m in cl["members"]]
        cl["tmin"] = min(s[0] for s in spans)
        cl["tmax"] = max(s[1] for s in spans)
        cl["spans"] = spans
        cl["kinds"] = sorted({m["kind"] for m in cl["members"]})
    return clusters


def _distinct_radii(values, tol):
    out = []
    for v in sorted(values):
        if not out or abs(v - out[-1]) > tol:
            out.append(v)
    return out


def _split_planes(records, clusters, tol):
    """Assign every planar face to a cluster as a cap, or as a wedge face through its axis."""
    caps = {i: [] for i in range(len(clusters))}
    wedges = {i: [] for i in range(len(clusters))}
    for rec in records:
        if rec["kind"] != "plane":
            continue
        placed = False
        for i, cl in enumerate(clusters):
            if _collinear(rec["n"], cl["dir"]):
                caps[i].append(rec)
                placed = True
                break
            if _perpendicular(rec["n"], cl["dir"]) and abs(
                    _dot(_sub(rec["p"], cl["loc"]), rec["n"])) <= tol:
                wedges[i].append(rec)
                placed = True
                break
        if not placed:
            raise Declined("a planar face is neither a cap nor a wedge of any axis cluster")
    return caps, wedges


# ------------------------------------------------------------------------------------------
# Tier 1: whole-part primitives
# ------------------------------------------------------------------------------------------

def _match_box(records, tol):
    planes = [r for r in records if r["kind"] == "plane"]
    if len(planes) != len(records):
        return None
    if len(planes) != 6:
        raise Declined(f"{len(planes)} planar faces: not a six-plane box")
    used = [False] * 6
    # Per axis: the outward normal `n` of one face of the pair, the mid-plane offset along `n`,
    # and the half-thickness. Stating it in offsets along `n` rather than in point differences
    # removes the ordering question: the face whose *outward* normal is `n` is by definition the
    # one at the larger offset, so the half-thickness is positive whichever face is found first.
    axes = []
    for i in range(6):
        if used[i]:
            continue
        for j in range(i + 1, 6):
            if used[j] or _dot(planes[i]["n"], planes[j]["n"]) > 0.0:
                continue
            if _collinear(planes[i]["n"], planes[j]["n"]):
                used[i] = used[j] = True
                n = _unit(planes[i]["n"])
                di = _dot(planes[i]["p"], n)
                dj = _dot(planes[j]["p"], n)
                axes.append((n, (di + dj) / 2.0, (di - dj) / 2.0))
                break
        else:
            raise Declined("a box face has no opposite partner")
    if len(axes) != 3:
        raise Declined("the six planes do not form three opposite pairs")
    for a in range(3):
        for b in range(a + 1, 3):
            if not _perpendicular(axes[a][0], axes[b][0]):
                raise Declined("the three plane pairs are not mutually perpendicular")
    if any(half <= 0.0 for _n, _mid, half in axes):
        raise Declined("the plane pair separations are not positive (inverted orientations?)")
    # A box is symmetric under flipping any of its own axes and under permuting them, so when the
    # three axes *are* the coordinate axes the frame can be relabelled into the identity without
    # changing the solid -- and then the emitter writes a bare TGeoBBox rather than a composite.
    snapped = [_snap_to_coordinate_axis(n) for n, _mid, _half in axes]
    if all(s is not None for s in snapped) and len({s[0] for s in snapped}) == 3:
        ordered = [None, None, None]
        for (index, sign), (_n, mid, half) in zip(snapped, axes):
            ordered[index] = (_COORDINATE_AXES[index], sign * mid, half)
        axes = ordered
    elif _dot(_cross(axes[0][0], axes[1][0]), axes[2][0]) < 0.0:
        axes = [axes[0], axes[2], axes[1]]                       # keep the frame right-handed
    x, y, z = axes[0][0], axes[1][0], _cross(axes[0][0], axes[1][0])
    halves = [axes[k][2] for k in range(3)]
    origin = (0.0, 0.0, 0.0)
    for k, axis in enumerate((x, y, z)):
        origin = _add(origin, _scale(axis, axes[k][1]))
    frame = {"origin": [float(c) for c in origin], "x": list(x), "y": list(y), "z": list(z)}
    return prim.candidate("primitive", [prim.leaf(
        "TGeoBBox", {"dx": halves[0], "dy": halves[1], "dz": halves[2]}, frame)], "tier1-box")


def _match_axial_primitive(records, clusters, caps, wedges, tol):
    """One axis cluster + two caps: a tube, a tube segment or a cone."""
    cl = clusters[0]
    cap = caps[0]
    wedge = wedges[0]
    if len(cap) != 2:
        raise Declined(f"{len(cap)} cap plane(s) perpendicular to the axis, expected 2")
    if len(wedge) not in (0, 2):
        raise Declined(f"{len(wedge)} wedge plane(s) through the axis, expected 0 or 2")
    t_caps = sorted(_dot(_sub(c["p"], cl["loc"]), cl["dir"]) for c in cap)
    if t_caps[1] - t_caps[0] <= 0.0:
        raise Declined("the two caps are coincident")
    # The caps bound the solid; the lateral faces must not stick out of them.
    if cl["tmin"] < t_caps[0] - tol or cl["tmax"] > t_caps[1] + tol:
        raise Declined("a lateral face extends beyond the cap planes")
    dz = (t_caps[1] - t_caps[0]) / 2.0
    centre = _add(cl["loc"], _scale(cl["dir"], (t_caps[0] + t_caps[1]) / 2.0))

    if cl["kinds"] == ["cylinder"]:
        radii = _distinct_radii([m["r"] for m in cl["members"]], tol)
        if len(radii) > 2:
            raise Declined(f"{len(radii)} distinct coaxial radii, expected 1 or 2")
        rmin = radii[0] if len(radii) == 2 else 0.0
        rmax = radii[-1]
        outer = [m for m in cl["members"] if abs(m["r"] - rmax) <= tol]
        frame = prim.frame_from_axis(centre, cl["dir"], outer[0]["x"])
        if wedge:
            phi1, phi2 = _phi_range(outer, frame)
            return prim.candidate("primitive", [prim.leaf(
                "TGeoTubeSeg", {"rmin": rmin, "rmax": rmax, "dz": dz, "phi1": phi1,
                                "phi2": phi2}, frame)], "tier1-tubeseg")
        return prim.candidate("primitive", [prim.leaf(
            "TGeoTube", {"rmin": rmin, "rmax": rmax, "dz": dz}, frame)], "tier1-tube")

    if cl["kinds"] == ["cone"]:
        if wedge:
            raise Declined("a phi-cut cone is out of scope (TGeoConeSeg not emitted)")
        if len(cl["members"]) > 2:
            raise Declined(f"{len(cl['members'])} coaxial cone faces, expected 1 or 2")
        radii_at = []
        for member in cl["members"]:
            radii_at.append(_cone_radii_at(member, cl, t_caps[0], t_caps[1]))
        radii_at.sort(key=lambda rr: rr[0] + rr[1])
        if len(radii_at) == 2:
            (rmin1, rmin2), (rmax1, rmax2) = radii_at
        else:
            (rmax1, rmax2), = radii_at
            rmin1 = rmin2 = 0.0
        frame = prim.frame_from_axis(centre, cl["dir"], cl["members"][0]["x"])
        return prim.candidate("primitive", [prim.leaf(
            "TGeoCone", {"dz": dz, "rmin1": rmin1, "rmax1": rmax1, "rmin2": rmin2,
                         "rmax2": rmax2}, frame)], "tier1-cone")

    raise Declined(f"mixed lateral surface kinds {cl['kinds']} on one axis")


def _cone_radii_at(member, cl, t0, t1):
    sa, ca = math.sin(member["a"]), math.cos(member["a"])
    base = _dot(_sub(member["p"], cl["loc"]), cl["dir"])
    sign = 1.0 if _dot(member["d"], cl["dir"]) > 0.0 else -1.0
    out = []
    for t in (t0, t1):
        v = sign * (t - base) / ca if ca != 0.0 else 0.0
        out.append(abs(member["r"] + v * sa))
    return out[0], out[1]


def _phi_range(outer_faces, frame):
    """Absolute phi bounds, in degrees, of a wedge, measured in the emitted frame's x/y."""
    lo, hi = None, None
    for face in outer_faces:
        umin, umax, _v0, _v1 = face["uv"]
        # the face's own reference direction may differ from the frame's x
        offset = math.atan2(_dot(face["x"], frame["y"]), _dot(face["x"], frame["x"]))
        for u in (umin + offset, umax + offset):
            lo = u if lo is None else min(lo, u)
            hi = u if hi is None else max(hi, u)
    span = math.degrees(hi - lo)
    if span >= 360.0 - 1.0e-6:
        raise Declined("the wedge spans a full turn")
    return math.degrees(lo), math.degrees(hi)


def _match_sphere(records, tol):
    spheres = [r for r in records if r["kind"] == "sphere"]
    if not spheres:
        return None
    if len(spheres) != len(records):
        raise Declined("a sphere with additional faces is out of scope (no theta/phi cuts)")
    radii = _distinct_radii([s["r"] for s in spheres], tol)
    centre = spheres[0]["p"]
    for s in spheres[1:]:
        if _norm(_sub(s["p"], centre)) > tol:
            raise Declined("spherical faces are not concentric")
    if len(radii) != 1:
        raise Declined(f"{len(radii)} distinct concentric sphere radii, expected 1")
    return prim.candidate("primitive", [prim.leaf(
        "TGeoSphere", {"rmin": 0.0, "rmax": radii[0]},
        prim.identity_frame(centre))], "tier1-sphere")


# ------------------------------------------------------------------------------------------
# The revolved profile: one axis, any number of z sections -> TGeoPcon
# ------------------------------------------------------------------------------------------
#
# `_match_axial_primitive` above stops at two caps and one lateral kind, which is a tube, a tube
# segment or a cone and nothing else. Most of what the detector corpora decline is the next thing
# along: one axis, several z sections, cylinders and cones mixed on it -- a `TGeoPcon`. This
# rebuilds that (r, z) profile.
#
# How the profile is read off, and why it is read off this way
# ------------------------------------------------------------
# Every lateral face is a straight segment in the (r, z) half-plane (a cylinder is vertical, a
# cone is slanted) and its endpoints are exact, from `_axial_extent`. The z levels of the profile
# are therefore the union of every lateral endpoint and every cap plane. Between two adjacent
# levels the cross-section is an annulus, and *which* annulus is read at the interval's midpoint,
# where the answer cannot be confused by a face that ends exactly on a level. One radius there
# means rmin = 0; two mean rmin and rmax; three or more is not a polycone and is declined.
#
# Nothing about orientation is used. The first thing measured on a revolved solid built by
# `O2_TGeoToCAD.conv_pcon` is that its inner faces are *not* reversed and their axis direction is
# flipped instead, so a rule keyed on `TopAbs_REVERSED` would have read the bore as the barrel on
# half the corpus.
#
# The one quantity that decides
# -----------------------------
# `Stream_K_Tier0.md` §3: a recogniser must score every candidate by one measured geometric gap
# relative to the part's own scale, never by a per-class criterion and never by an angle. So the
# proposal above is *rebuilt back into* the (r, z) half-plane and `_profile_gap` measures the
# largest distance from a point that is actually on the solid's boundary -- every vertex, and
# samples along every lateral -- to the proposed profile. That number, over the bounding-box
# diagonal, is the acceptance criterion of this matcher. It is also why the matcher may be
# generous about structure: a reconstruction that drifted anywhere shows up as a distance.

_PROFILE_SAMPLES = (0.0, 0.25, 0.5, 0.75, 1.0)


def _merge_levels(values, tol):
    """Sorted distinct z levels from `(value, exact)` pairs, merging anything within `tol`.

    Within a merged group the value marked exact wins. A cap plane states its own z directly,
    whereas a lateral face's endpoint comes through a parametric bound and carries a few ulp of
    arithmetic with it; taking the cap's number keeps the emitted `z` array reading as the
    engineer wrote it instead of as -8.9e-16.
    """
    out = []
    for value, exact in sorted(values):
        if out and value - out[-1][0] <= tol:
            if exact and not out[-1][1]:
                out[-1] = (value, True)
            continue
        out.append((value, exact))
    return [value for value, _exact in out]


def _span_radius_at(span, t):
    """The lateral's radius at axial coordinate `t`; the segment is straight in (r, z)."""
    t0, t1, r0, r1 = span[0], span[1], span[2], span[3]
    if t1 - t0 <= 0.0:
        return r0
    return r0 + (r1 - r0) * (t - t0) / (t1 - t0)


def _point_segment_distance(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    length2 = dx * dx + dy * dy
    if length2 <= 0.0:
        return math.hypot(p[0] - a[0], p[1] - a[1])
    s = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2
    s = min(1.0, max(0.0, s))
    return math.hypot(p[0] - (a[0] + s * dx), p[1] - (a[1] + s * dy))


def _profile_gap(profile, samples):
    """The largest distance, in cm, from a boundary sample `(r, z)` to the profile's outline.

    `profile` is the closed (r, z) ring the candidate will be revolved from. This is the whole
    acceptance criterion of `_match_revolved` and it is a distance, so it is comparable across
    parts once divided by the part's scale -- which is the property `Stream_K_Tier0.md` §3 says a
    recogniser's score must have and which an angle or a per-class residual does not.
    """
    worst = 0.0
    n = len(profile)
    for sample in samples:
        best = float("inf")
        for i in range(n):
            best = min(best, _point_segment_distance(sample, profile[i], profile[(i + 1) % n]))
            if best <= 0.0:
                break
        worst = max(worst, best)
    return worst


def _solid_vertices(solid):
    """Every vertex of the solid, in the part frame.

    Taken from the topology rather than from the carriers the profile was built out of, which is
    what makes `_profile_gap` a measurement and not a restatement.
    """
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.TopAbs import TopAbs_VERTEX
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
    out = []
    seen = set()
    exp = TopExp_Explorer(solid, TopAbs_VERTEX)
    while exp.More():
        pnt = BRep_Tool.Pnt(topods.Vertex(exp.Current()))
        exp.Next()
        key = (round(pnt.X(), 9), round(pnt.Y(), 9), round(pnt.Z(), 9))
        if key in seen:
            continue
        seen.add(key)
        out.append((pnt.X(), pnt.Y(), pnt.Z()))
    return out


def _canonical_revolved_leaf(lf, origin, axis, tol):
    """Say a two-section full-turn profile in the native class it already is.

    A profile that came out with exactly two sections and a full turn is a `TGeoCone`, or a
    `TGeoTube` when neither radius changes. Emitting it as a two-section `TGeoPcon` would be the
    same solid, but ROOT's own class carries what it is, so `checkKnownSource.py` compares like
    with like and nothing downstream has to special-case a polycone that is a cone. A wedge or a
    z-step needs more than two sections' worth of description and stays a `TGeoPcon`.

    Returns `(leaf, recogniser tag)`.
    """
    p = lf["params"]
    z, rmin, rmax = p["z"], p["rmin"], p["rmax"]
    if len(z) != 2 or abs(p["dphi"] - 360.0) > 1.0e-9:
        return lf, "revolved-pcon"
    # TGeoTube and TGeoCone are centred on their own frame, so the frame's origin moves to the
    # middle of the section pair; the axis and the reference x are unchanged.
    frame = dict(lf["frame"])
    frame["origin"] = [float(c) for c in _add(origin, _scale(axis, 0.5 * (z[0] + z[1])))]
    dz = 0.5 * (z[1] - z[0])
    if abs(rmin[0] - rmin[1]) <= tol and abs(rmax[0] - rmax[1]) <= tol:
        return prim.leaf("TGeoTube", {"rmin": 0.5 * (rmin[0] + rmin[1]),
                                      "rmax": 0.5 * (rmax[0] + rmax[1]), "dz": dz},
                         frame), "revolved-tube"
    return prim.leaf("TGeoCone", {"dz": dz, "rmin1": rmin[0], "rmax1": rmax[0],
                                  "rmin2": rmin[1], "rmax2": rmax[1]},
                     frame), "revolved-cone"


def _match_revolved(solid, records, clusters, caps, wedges, tol, diag):
    """One axis cluster, any number of z sections: a `TGeoPcon`."""
    cl = clusters[0]
    cap = caps[0]
    wedge = wedges[0]

    # A coordinate axis is taken pointing the positive way, so a part whose axis is the global z
    # through the global origin comes out with an identity frame and no placement at all -- and
    # then the emitted TGeoPcon's own parameters are directly comparable with the source shape's.
    axis = cl["dir"]
    snapped = _snap_to_coordinate_axis(axis)
    if snapped is not None and snapped[1] < 0.0:
        axis = _scale(axis, -1.0)
    # The axial origin is the foot of the perpendicular from the part frame's origin, for the
    # same reason: it is the one choice that is a property of the axis rather than of whichever
    # face happened to be seen first.
    origin = _sub(cl["loc"], _scale(axis, _dot(cl["loc"], axis)))

    spans = []
    for member in cl["members"]:
        t0, t1, r0, r1 = _axial_extent(member, axis, origin)
        if t1 - t0 <= tol:
            raise Declined("a lateral face has no axial extent (a cone at its own apex?)")
        if min(r0, r1) < -tol:
            raise Declined("a lateral face reaches a negative radius")
        spans.append((t0, t1, max(r0, 0.0), max(r1, 0.0), member))

    t_caps = [_dot(_sub(c["p"], origin), axis) for c in cap]
    levels = _merge_levels([(s[0], False) for s in spans] + [(s[1], False) for s in spans]
                           + [(t, True) for t in t_caps], tol)
    if len(levels) < 2:
        raise Declined("the axial faces span fewer than two distinct z levels")

    sections = []                       # (z, rmin, rmax), in profile order
    for k in range(len(levels) - 1):
        lo, hi = levels[k], levels[k + 1]
        mid = 0.5 * (lo + hi)
        here = [s for s in spans if s[0] - tol <= mid <= s[1] + tol]
        radii = _distinct_radii([_span_radius_at(s, mid) for s in here], tol)
        if not radii:
            raise Declined(f"no lateral face covers the z range [{lo:.6g}, {hi:.6g}]: "
                           "the solid is not one connected polycone")
        if len(radii) > 2:
            raise Declined(f"{len(radii)} distinct coaxial radii between z = {lo:.6g} and "
                           f"{hi:.6g}, expected 1 or 2")
        ends = []
        for t in (lo, hi):
            at = [_span_radius_at(s, t) for s in here]
            ends.append((min(at) if len(radii) == 2 else 0.0, max(at)))
        if k == 0:
            sections.append((lo, ends[0][0], ends[0][1]))
        elif (abs(sections[-1][1] - ends[0][0]) > tol
              or abs(sections[-1][2] - ends[0][1]) > tol):
            # A z-step: TGeo states it as two sections sharing one z, which is legal and is what
            # the writer's STEP already contains as a cap annulus at that plane.
            sections.append((lo, ends[0][0], ends[0][1]))
        sections.append((hi, ends[1][0], ends[1][1]))

    # Every radial jump in the profile is a face of the solid, so it must be there. This is what
    # separates a polycone from an open shell that merely looks like one.
    for k, (z, rmin, rmax) in enumerate(sections):
        needs_cap = (rmax - rmin > tol) if k in (0, len(sections) - 1) else \
                    (k + 1 < len(sections) and abs(sections[k + 1][0] - z) <= tol)
        if needs_cap and not any(abs(tc - z) <= tol for tc in t_caps):
            raise Declined(f"the profile steps or ends at z = {z:.6g} with no cap plane there")

    if wedge:
        normals = []
        for w in wedge:
            if not any(_collinear(w["n"], n) for n in normals):
                normals.append(w["n"])
        if len(normals) > 2:
            raise Declined(f"{len(normals)} distinct half-planes through the axis: "
                           "not a single phi wedge")

    # phi is read only off laterals whose own axis runs *with* the frame's, because a face whose
    # carrier axis is flipped parametrises phi the other way round and would report the wedge
    # mirrored.
    oriented = [m for m in cl["members"] if _parallel(m["d"], axis)]
    # On a coordinate axis the frame is pinned to the coordinate frame, which makes it the
    # identity, drops the placement, and leaves phi1 stated absolutely -- so the emitted
    # TGeoPcon's own numbers are the source shape's numbers and `checkKnownSource.py` can compare
    # them directly. Off a coordinate axis there is no canonical x and a lateral supplies it.
    ref_x = None if _snap_to_coordinate_axis(axis) is not None else (
        oriented[0]["x"] if oriented else None)
    frame = prim.frame_from_axis(origin, axis, ref_x)
    if wedge:
        if not oriented:
            raise Declined("no lateral face runs with the axis, so the phi wedge cannot be read")
        lo_phi, hi_phi = _phi_range(oriented, frame)
        phi1, dphi = lo_phi, hi_phi - lo_phi
    else:
        phi1, dphi = 0.0, 360.0

    z = [s[0] for s in sections]
    rmin = [s[1] for s in sections]
    rmax = [s[2] for s in sections]
    try:
        lf = prim.leaf("TGeoPcon", {"phi1": phi1, "dphi": dphi, "z": z, "rmin": rmin,
                                    "rmax": rmax}, frame)
    except ValueError as bad:
        raise Declined(f"the reconstructed profile is not a legal TGeoPcon: {bad}") from bad

    # The one measured quantity. It is taken on the reconstructed *profile*, before the leaf is
    # canonicalised below, because the gap is a property of the reconstruction and not of the
    # ROOT class it ends up wearing.
    # The same ring `build_occ` will revolve, deduped by the same rule, so the gap is measured
    # against the candidate that the acceptance test will actually be handed.
    profile = prim.pcon_profile_rz(lf["params"])
    if len(profile) < 3:
        raise Declined("the reconstructed (r, z) profile has fewer than three corners")
    samples = []
    for point in _solid_vertices(solid):
        rel = _sub(point, origin)
        zc = _dot(rel, axis)
        samples.append((math.sqrt(max(_dot(rel, rel) - zc * zc, 0.0)), zc))
    for span in spans:
        for f in _PROFILE_SAMPLES:
            t = span[0] + f * (span[1] - span[0])
            samples.append((_span_radius_at(span, t), t))
    gap = _profile_gap(profile, samples)
    scale = max(diag, 1.0)
    if gap > REL_TOL * scale:
        raise Declined(f"the boundary is {gap:.3g} cm off the reconstructed profile "
                       f"({gap / scale:.3g} of the part's {diag:.6g} cm diagonal, "
                       f"over {REL_TOL:.0e})")

    lf, tag = _canonical_revolved_leaf(lf, origin, axis, tol)
    return prim.candidate("primitive", [lf], tag,
                          notes={"nz": len(z), "nCaps": len(cap), "nWedges": len(wedge),
                                 "nLaterals": len(cl["members"]),
                                 "profileGapCm": gap, "profileGapRelative": gap / scale})


# ------------------------------------------------------------------------------------------
# The prism family: an all-planar face graph -> Trd1 / Trd2 / Arb8 / Pgon / Xtru
# ------------------------------------------------------------------------------------------
#
# `_match_box` above stops at six planes in three perpendicular opposite pairs, which is a box and
# nothing else. What the detector corpora decline next is everything else the writer builds out of
# a *stack of sections* (`O2_TGeoToCAD._prism_from_rings`): a `TGeoTrd1` whose x faces slant, a
# `TGeoArb8`, a `TGeoXtru` swept from a general polygon, a `TGeoPgon` whose laterals are planes at
# the apothem radius. All five are one construction and this rebuilds it.
#
# How the sections are read off
# -----------------------------
# A prism has an axis: the direction of an opposite plane pair that no third plane shares. Every
# vertex of the solid then sits at one of a small number of *levels* along it, every face lies
# either wholly in one level (a cap) or across two adjacent ones (a side), and each side joins one
# edge of the section below to one edge of the section above. Walking that adjacency from the
# bottom cap's own wire propagates the corner *order* upwards, which is the part a purely
# geometric reading cannot supply: a non-convex ITS rib section is not recoverable by sorting its
# corners by angle, and half the corpus's Xtru polygons are non-convex.
#
# Nothing here decides which ROOT class the solid is. The templates are tried most-specific
# first -- Trd1, Trd2, Pgon, Arb8, Xtru -- and each *proposes* parameters that the same one
# measured quantity then scores.
#
# The one quantity that decides
# -----------------------------
# `Stream_K_Tier0.md` §3: one measured geometric gap, relative to the part's own scale, never a
# per-class criterion and never an angle. Each proposal is rebuilt into the rings `build_occ` will
# sew (`primitives.prism_samples`) and `_point_set_gap` measures the symmetric Hausdorff distance,
# in cm, between those corners and edge midpoints and the *solid's own* corners and edge midpoints
# read from its topology. Corners alone would not be enough -- a hexahedron whose corner order is
# wrong has the same eight corners and a different solid -- which is why the edge midpoints are in
# both sets. That distance over the bounding-box diagonal is this matcher's whole acceptance
# criterion, and the OCCT symmetric difference is still the judge afterwards.
#
# Twisted `TGeoArb8`: out of scope, and declined before this matcher
# -----------------------------------------------------------------
# A `TGeoArb8` whose four lateral corners are not coplanar is written by `_quad_face` as a ruled
# `brepfill` patch, which is a B-spline surface and not a plane. Measured on a genuinely twisted
# hexahedron: 4 of its 6 faces come back free-form, so `_face_records` declines the solid before
# any matcher sees it, and the reason says so. Ruling it in would mean recognising a bilinear
# patch inside a general B-spline surface, which is a recognition problem of a different kind from
# the one this module solves; a declined part costs coverage and a wrongly accepted one costs
# correctness. Neither corpus contains one: TPC's five `TGeoArb8` volumes and ABSO's two
# `TGeoTrap`s are all-planar and are recognised here.

# Most specific first. `Trd1` and `Trd2` are the most constrained; a `Pgon` states a polygon ring
# in three numbers and two radii per section where an `Arb8` would need sixteen coordinates, and
# it is the only one of the five that can carry an annular section at all; `Arb8` is the general
# hexahedron; `Xtru` is the general right prism. The template loop is the OUTER one and the axis
# loop the inner one, because a `TGeoTrd1` is also a perfectly good `TGeoArb8` about its own y
# axis -- measured, on the writer's own `TGeoTrd1` output -- and the class the solid *is* must win
# over the axis that happened to be enumerated first.
_PRISM_TEMPLATES = ("trd1", "trd2", "pgon", "arb8", "xtru")


def _face_wires(face):
    """Ordered corner points of each wire of a planar face, the outer wire first."""
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.BRepTools import BRepTools_WireExplorer, breptools
    from OCC.Core.TopAbs import TopAbs_WIRE
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
    outer = breptools.OuterWire(face)
    found = []
    exp = TopExp_Explorer(face, TopAbs_WIRE)
    while exp.More():
        wire = topods.Wire(exp.Current())
        exp.Next()
        pts = []
        walk = BRepTools_WireExplorer(wire, face)
        while walk.More():
            pnt = BRep_Tool.Pnt(walk.CurrentVertex())
            pts.append((pnt.X(), pnt.Y(), pnt.Z()))
            walk.Next()
        if pts:
            found.append((not wire.IsSame(outer), pts))
    found.sort(key=lambda item: item[0])
    return [pts for _is_hole, pts in found]


def _solid_samples(solid):
    """Every vertex and every edge midpoint of the solid, in the part frame.

    Read from the topology, never from the carriers a proposal was built out of, which is what
    makes `_point_set_gap` a measurement rather than a restatement.
    """
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.TopAbs import TopAbs_EDGE
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
    out = list(_solid_vertices(solid))
    seen = set()
    exp = TopExp_Explorer(solid, TopAbs_EDGE)
    while exp.More():
        edge = topods.Edge(exp.Current())
        exp.Next()
        curve, first, last = BRep_Tool.Curve(edge)
        if curve is None:
            continue
        pnt = curve.Value(0.5 * (first + last))
        key = (round(pnt.X(), 9), round(pnt.Y(), 9), round(pnt.Z(), 9))
        if key in seen:
            continue
        seen.add(key)
        out.append((pnt.X(), pnt.Y(), pnt.Z()))
    return out


def _point_set_gap(a, b):
    """The symmetric Hausdorff distance, in cm, between two point sets."""
    def one_way(u, v):
        worst = 0.0
        for p in u:
            best = float("inf")
            for q in v:
                d2 = ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2)
                if d2 < best:
                    best = d2
                    if best == 0.0:
                        break
            worst = max(worst, best)
        return math.sqrt(worst)
    return max(one_way(a, b), one_way(b, a))


class _Corners:
    """Nearest-corner lookup, so a corner reported twice by two faces is one corner."""

    def __init__(self, points, tol):
        self.points = list(points)
        self.tol = tol

    def find(self, p):
        best, best_d = None, self.tol
        for i, q in enumerate(self.points):
            d = _norm(_sub(p, q))
            if d <= best_d:
                best, best_d = i, d
        return best


def _prism_axis_candidates(records):
    """Directions the part could be a stack of sections along.

    A cap pair is an opposite plane pair that *no third plane shares*: on a hollow `TGeoPgon` the
    laterals do come in opposite pairs too, but each such direction carries four collinear planes
    (two outer, two inner), so the axis is not confused with them.
    """
    out = []
    for rec in records:
        direction = _unit(rec["n"])
        same = [j for j, other in enumerate(records) if _collinear(other["n"], direction)]
        if len(same) != 2:
            continue
        a, b = same
        if _dot(records[a]["n"], records[b]["n"]) > 0.0:
            continue
        if _dot(records[a]["p"], direction) < _dot(records[b]["p"], direction):
            direction = _scale(direction, -1.0)
        snapped = _snap_to_coordinate_axis(direction)
        if snapped is not None:
            direction = _COORDINATE_AXES[snapped[0]]
        if any(_collinear(direction, seen) for seen in out):
            continue
        out.append(direction)
    # z, then y, then x, then anything else: every ROOT prism class is stated about its own z, so
    # a part the writer wrote in its own frame is read back in that frame and its emitted
    # parameters are the source's own numbers.
    def order(direction):
        snapped = _snap_to_coordinate_axis(direction)
        return (1, 0) if snapped is None else (0, -snapped[0])
    out.sort(key=order)
    return out


def _canonical_ring_start(ring):
    """Rotate a ring to start at its lexicographically smallest corner.

    A ring's *order* is a property of the solid; where the order starts is a property of whichever
    wire OCCT happened to walk first. Pinning it makes the emitted `TGeoArb8` corner list and
    `TGeoXtru` polygon a function of the geometry alone, so the same CAD emits the same numbers
    whoever read it.
    """
    start = min(range(len(ring)), key=lambda i: (round(ring[i][0], 12), round(ring[i][1], 12),
                                                 round(ring[i][2], 12)))
    return ring[start:] + ring[:start]


def _ring_signed_area(ring, ex, ey):
    total = 0.0
    for i, a in enumerate(ring):
        b = ring[(i + 1) % len(ring)]
        total += _dot(a, ex) * _dot(b, ey) - _dot(b, ex) * _dot(a, ey)
    return 0.5 * total


def _prism_sections(solid, records, axis, tol):
    """`(levels, rings)` along `axis`: the ordered corner rings of every section.

    `rings[k]` is the list of wires at level `k` -- one for a solid section, two for the annular
    section of a hollow `TGeoPgon` -- each in corner order and counterclockwise about `axis`, with
    corner `i` of section `k` joined to corner `i` of section `k + 1`.
    """
    corners = _Corners(_solid_vertices(solid), tol)
    levels = _merge_levels([(_dot(v, axis), False) for v in corners.points], tol)
    if len(levels) < 2:
        raise Declined("every corner sits on one plane: the part has no extent along the axis")

    def level_of(point):
        t = _dot(point, axis)
        best = min(range(len(levels)), key=lambda k: abs(levels[k] - t))
        return best if abs(levels[best] - t) <= tol else None

    caps, sides = [], []
    for rec in records:
        wires = _face_wires(rec["face"])
        if not wires:
            raise Declined("a planar face has no wire")
        seen = {level_of(p) for wire in wires for p in wire}
        if None in seen:
            raise Declined("a face corner sits on no section plane of the axis")
        if len(seen) == 1:
            caps.append((wires, seen.pop()))
        elif len(seen) == 2 and max(seen) - min(seen) == 1 and len(wires) == 1:
            sides.append((wires[0], min(seen)))
        else:
            raise Declined(f"a planar face spans sections {sorted(seen)} "
                           "and is neither a cap nor a single prism side")
    if len(caps) != 2:
        raise Declined(f"{len(caps)} face(s) lie wholly in one section plane, expected 2 caps")
    caps.sort(key=lambda cap: cap[1])
    if caps[0][1] != 0 or caps[1][1] != len(levels) - 1:
        raise Declined("the two caps are not the outermost sections")
    if len(caps[0][0]) != len(caps[1][0]):
        raise Declined(f"the caps carry {len(caps[0][0])} and {len(caps[1][0])} wires")
    if len(caps[0][0]) > 2:
        raise Declined(f"a cap has {len(caps[0][0])} wires: more than one hole is out of scope")

    # Which corner of the section above each side reaches, keyed on the ordered corner pair it
    # spans below. Both traversal directions are registered, since a wire's own direction is a
    # property of the face's orientation and not of the ring order being propagated.
    step = {}
    for wire, k in sides:
        n = len(wire)
        if n not in (3, 4):
            raise Declined(f"a prism side has {n} corners, expected 3 or 4")
        at = [level_of(p) for p in wire]
        for a in range(n):
            b = (a + 1) % n
            if at[a] != k or at[b] != k:
                continue
            if n == 4:
                up_b, up_a = wire[(a + 2) % n], wire[(a + 3) % n]
            else:
                up_a = up_b = wire[(a + 2) % n]
            step[(k, corners.find(wire[a]), corners.find(wire[b]))] = (up_a, up_b)

    frame = prim.frame_from_axis((0.0, 0.0, 0.0), axis)
    ex, ey = tuple(frame["x"]), tuple(frame["y"])
    rings = [[_canonical_ring_start(
        list(wire) if _ring_signed_area(wire, ex, ey) > 0.0 else list(reversed(wire)))
        for wire in caps[0][0]]]
    for k in range(len(levels) - 1):
        above = []
        for ring in rings[k]:
            n = len(ring)
            up = [None] * n
            for i in range(n):
                j = (i + 1) % n
                key = (k, corners.find(ring[i]), corners.find(ring[j]))
                if key in step:
                    up_i, up_j = step[key]
                else:
                    key = (k, corners.find(ring[j]), corners.find(ring[i]))
                    if key not in step:
                        raise Declined(f"no prism side joins two corners of section {k}: "
                                       "the sections do not stack")
                    up_j, up_i = step[key]
                for pos, value in ((i, up_i), (j, up_j)):
                    if up[pos] is not None and _norm(_sub(up[pos], value)) > tol:
                        raise Declined("two prism sides disagree about a corner of the next "
                                       "section")
                    up[pos] = value
            above.append(up)
        rings.append(above)
    used = {corners.find(p) for level in rings for ring in level for p in ring}
    distinct = {corners.find(p) for p in corners.points}
    if used != distinct:
        raise Declined(f"{len(distinct - used)} corner(s) of the solid lie on no section ring")
    return levels, rings


def _ring_xy(ring, frame):
    origin, ex, ey = tuple(frame["origin"]), tuple(frame["x"]), tuple(frame["y"])
    return [(_dot(_sub(p, origin), ex), _dot(_sub(p, origin), ey)) for p in ring]


def _in_plane_x_candidates(axis, ring):
    """Reference x directions to try, coordinate axes first so an aligned part stays aligned."""
    out = []
    snapped = _snap_to_coordinate_axis(axis)
    if snapped is not None:
        for k in range(3):
            if k != snapped[0]:
                out.append(_COORDINATE_AXES[k])
    for i, a in enumerate(ring):
        edge = _sub(ring[(i + 1) % len(ring)], a)
        flat = _sub(edge, _scale(axis, _dot(edge, axis)))
        if _norm(flat) > 1.0e-12:
            out.append(_unit(flat))
    return out


def _prism_leaf_gap(leaf, samples):
    """The one measured quantity: how far the proposal's boundary is from the solid's, in cm."""
    return _point_set_gap(prim.prism_samples(leaf), samples)


def _try_trd(levels, rings, axis, tol):
    """A `TGeoTrd1` or `TGeoTrd2`: two rectangular sections sharing a centre line."""
    if len(levels) != 2 or any(len(level) != 1 for level in rings):
        raise Declined("a TGeoTrd needs exactly two single-wire sections")
    lower, upper = rings[0][0], rings[1][0]
    if len(lower) != 4 or len(upper) != 4:
        raise Declined(f"a TGeoTrd needs four corners per section, got "
                       f"{len(lower)} and {len(upper)}")
    centre = _scale(_add(_centroid(lower), _centroid(upper)), 0.5)
    dz = 0.5 * (levels[1] - levels[0])
    out = []
    for ref_x in _in_plane_x_candidates(axis, lower):
        frame = prim.frame_from_axis(centre, axis, ref_x)
        low, high = _ring_xy(lower, frame), _ring_xy(upper, frame)
        dx1, dy1 = max(abs(p[0]) for p in low), max(abs(p[1]) for p in low)
        dx2, dy2 = max(abs(p[0]) for p in high), max(abs(p[1]) for p in high)
        if abs(dy1 - dy2) <= tol:
            out.append(("rung2-trd1", "TGeoTrd1",
                        {"dx1": dx1, "dx2": dx2, "dy": 0.5 * (dy1 + dy2), "dz": dz}, frame))
        out.append(("rung2-trd2", "TGeoTrd2",
                    {"dx1": dx1, "dx2": dx2, "dy1": dy1, "dy2": dy2, "dz": dz}, frame))
    return out


def _try_arb8(levels, rings, axis, tol):
    """A `TGeoArb8`: two four-corner sections, eight corners, stated as they are."""
    if len(levels) != 2 or any(len(level) != 1 for level in rings):
        raise Declined("a TGeoArb8 needs exactly two single-wire sections")
    lower, upper = rings[0][0], rings[1][0]
    if len(lower) != 4 or len(upper) != 4:
        raise Declined(f"a TGeoArb8 needs four corners per section, got "
                       f"{len(lower)} and {len(upper)}")
    origin = _scale(axis, 0.5 * (levels[0] + levels[1]))
    frame = prim.frame_from_axis(origin, axis, _prism_ref_x(axis))
    vertices = []
    for ring in (lower, upper):
        for corner in _ring_xy(ring, frame):
            vertices.extend([corner[0], corner[1]])
    return [("rung2-arb8", "TGeoArb8",
             {"dz": 0.5 * (levels[1] - levels[0]), "vertices": vertices}, frame)]


def _try_pgon(levels, rings, axis, tol):
    """A `TGeoPgon`: every section a regular polygon ring at the same set of angles.

    ROOT's rmin/rmax are the **apothem** radii -- the laterals are planes tangent to the inscribed
    circle -- so the corners sit at `r / cos(dseg / 2)` and that is what is inverted here.
    """
    frame = prim.frame_from_axis((0.0, 0.0, 0.0), axis, _prism_ref_x(axis))
    radial = []
    for level in rings:
        here = []
        for ring in level:
            for x, y in _ring_xy(ring, frame):
                here.append((math.hypot(x, y), math.atan2(y, x)))
        radial.append(here)
    biggest = max((r for here in radial for r, _a in here), default=0.0)
    if biggest <= tol:
        raise Declined("the sections have no radial extent")
    angle_tol = max(tol / biggest, ANG_TOL)
    angles = _merge_angles([a for here in radial for r, a in here if r > tol], angle_tol)
    if len(angles) < 2:
        raise Declined(f"{len(angles)} distinct corner angle(s): not a polygon ring")
    gaps = [(angles[(i + 1) % len(angles)] - angles[i]) % (2.0 * math.pi)
            for i in range(len(angles))]
    widest = max(range(len(gaps)), key=lambda i: gaps[i])
    if max(gaps) - min(gaps) <= angle_tol:
        nedges, dphi = len(angles), 360.0
        phi1 = angles[0]
        dseg = 2.0 * math.pi / nedges
    else:
        ordered = angles[widest + 1:] + angles[:widest + 1]
        steps = [(ordered[i + 1] - ordered[i]) % (2.0 * math.pi) for i in range(len(ordered) - 1)]
        if max(steps) - min(steps) > angle_tol:
            raise Declined("the corner angles are not equally spaced: not a polygon ring")
        dseg = sum(steps) / len(steps)
        nedges = len(steps)
        phi1 = ordered[0]
        dphi = math.degrees(dseg * nedges)
    half = math.cos(dseg / 2.0)
    rmin, rmax = [], []
    for here in radial:
        radii = _distinct_radii([r for r, _a in here if r > tol], tol)
        if len(radii) > 2:
            raise Declined(f"{len(radii)} distinct corner radii in one section, expected 1 or 2")
        if not radii:
            raise Declined("a section has no corner off the axis")
        rmax.append(radii[-1] * half)
        rmin.append(radii[0] * half if len(radii) == 2 else 0.0)
    return [("rung2-pgon", "TGeoPgon",
             {"phi1": math.degrees(phi1), "dphi": dphi, "nedges": nedges,
              "z": list(levels), "rmin": rmin, "rmax": rmax}, frame)]


def _try_xtru(levels, rings, axis, tol):
    """A `TGeoXtru`: one polygon, per section an offset and an isotropic scale."""
    if any(len(level) != 1 for level in rings):
        raise Declined("a TGeoXtru section is one closed polygon, and this part's is not")
    frame = prim.frame_from_axis((0.0, 0.0, 0.0), axis, _prism_ref_x(axis))
    sections = [_ring_xy(level[0], frame) for level in rings]
    nv = len(sections[0])
    if any(len(s) != nv for s in sections):
        raise Declined("the sections do not all carry the same number of corners")
    base = sections[0]
    centre0 = _centroid2(base)
    spread = sum((p[0] - centre0[0]) ** 2 + (p[1] - centre0[1]) ** 2 for p in base)
    if spread <= 0.0:
        raise Declined("the polygon has no extent")
    xoff, yoff, scale = [], [], []
    for section in sections:
        centre = _centroid2(section)
        num = sum((section[i][0] - centre[0]) * (base[i][0] - centre0[0])
                  + (section[i][1] - centre[1]) * (base[i][1] - centre0[1])
                  for i in range(nv))
        s = num / spread
        if s <= 0.0:
            raise Declined("a section scales to zero or turns the polygon inside out")
        scale.append(s)
        xoff.append(centre[0] - s * centre0[0])
        yoff.append(centre[1] - s * centre0[1])
    return [("rung2-xtru", "TGeoXtru",
             {"x": [p[0] for p in base], "y": [p[1] for p in base], "z": list(levels),
              "xoff": xoff, "yoff": yoff, "scale": scale}, frame)]


_PRISM_TRIES = {"trd1": _try_trd, "trd2": _try_trd, "arb8": _try_arb8,
                "pgon": _try_pgon, "xtru": _try_xtru}


def _prism_ref_x(axis):
    snapped = _snap_to_coordinate_axis(axis)
    if snapped is not None:
        return _COORDINATE_AXES[(snapped[0] + 1) % 3]
    return None


def _centroid(points):
    total = (0.0, 0.0, 0.0)
    for p in points:
        total = _add(total, p)
    return _scale(total, 1.0 / len(points))


def _centroid2(points):
    return (sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points))


def _merge_angles(values, tol):
    """Distinct angles in [0, 2pi), merged within `tol`, the wrap included."""
    out = []
    for a in sorted(v % (2.0 * math.pi) for v in values):
        if out and min(a - out[-1], (out[0] + 2.0 * math.pi) - a) <= tol:
            continue
        out.append(a)
    while len(out) > 1 and (out[0] + 2.0 * math.pi) - out[-1] <= tol:
        out.pop()
    return out


def _match_prism(solid, records, tol, diag):
    """One axis, a stack of planar sections: the `Trd1`/`Trd2`/`Pgon`/`Arb8`/`Xtru` family."""
    axes = _prism_axis_candidates(records)
    if not axes:
        raise Declined("no opposite plane pair that no third plane shares: no prism axis")
    samples = _solid_samples(solid)
    scale = max(diag, 1.0)
    best_gap, best_tag = None, None
    reasons = []
    sections = {}
    proposals = {}
    for template in _PRISM_TEMPLATES:
        for index, axis in enumerate(axes):
            if index not in sections:
                try:
                    sections[index] = _prism_sections(solid, records, axis, tol)
                except Declined as why:
                    sections[index] = None
                    reasons.append(str(why))
            if sections[index] is None:
                continue
            levels, rings = sections[index]
            if (template, index) not in proposals:
                try:
                    made = _PRISM_TRIES[template](levels, rings, axis, tol)
                except Declined as why:
                    made = []
                    reasons.append(f"as a {template}: {why}")
                proposals[(template, index)] = [item for item in made
                                                if item[0].endswith(template)]
            for tag, kind, params, frame in proposals[(template, index)]:
                try:
                    leaf = prim.leaf(kind, params, frame)
                except ValueError as bad:
                    reasons.append(f"as a {template}: not a legal {kind}: {bad}")
                    continue
                gap = _prism_leaf_gap(leaf, samples)
                if best_gap is None or gap < best_gap:
                    best_gap, best_tag = gap, tag
                if gap <= REL_TOL * scale:
                    return prim.candidate(
                        "primitive", [leaf], tag,
                        notes={"nSections": len(levels), "nWires": len(rings[0]),
                               "nCorners": sum(len(r) for r in rings[0]),
                               "prismGapCm": gap, "prismGapRelative": gap / scale})
    if best_gap is not None:
        raise Declined(f"the boundary is {best_gap:.3g} cm off the closest prism template "
                       f"({best_tag}, {best_gap / scale:.3g} of the part's {diag:.6g} cm "
                       f"diagonal, over {REL_TOL:.0e})")
    raise Declined("; ".join(dict.fromkeys(reasons)) or "no prism template applies")


# ------------------------------------------------------------------------------------------
# Tier 2: the two-cluster union
# ------------------------------------------------------------------------------------------

def _match_two_cluster_union(records, clusters, caps, wedges, tol):
    if any(wedges[i] for i in range(len(clusters))):
        raise Declined("a wedge plane in a two-cluster part is out of scope")
    axes = [cl["dir"] for cl in clusters]
    if _collinear(axes[0], axes[1]):
        raise Declined("the two clusters are parallel: not the lug case")
    leaves = []
    for i, cl in enumerate(clusters):
        if cl["kinds"] != ["cylinder"]:
            raise Declined(f"cluster {i} has lateral kinds {cl['kinds']}, expected cylinders only")
        radii = _distinct_radii([m["r"] for m in cl["members"]], tol)
        if len(radii) > 2:
            raise Declined(f"cluster {i} has {len(radii)} distinct radii, expected 1 or 2")
        rmin = radii[0] if len(radii) == 2 else 0.0
        rmax = radii[-1]
        t0, t1 = cl["tmin"], cl["tmax"]
        for cap in caps[i]:
            t = _dot(_sub(cap["p"], cl["loc"]), cl["dir"])
            t0, t1 = min(t0, t), max(t1, t)
        if t1 - t0 <= 0.0:
            raise Declined(f"cluster {i} has no axial extent")
        centre = _add(cl["loc"], _scale(cl["dir"], (t0 + t1) / 2.0))
        outer = [m for m in cl["members"] if abs(m["r"] - rmax) <= tol]
        frame = prim.frame_from_axis(centre, cl["dir"], outer[0]["x"])
        leaves.append(prim.leaf("TGeoTube", {"rmin": rmin, "rmax": rmax,
                                             "dz": (t1 - t0) / 2.0}, frame))
    return prim.candidate("union", leaves, "tier2-tube-union",
                          notes={"nCaps": [len(caps[i]) for i in range(len(clusters))]})


# ------------------------------------------------------------------------------------------
# entry point
# ------------------------------------------------------------------------------------------

def recognise(solid):
    """Propose a CSG description for one leaf solid in cm. Returns (candidate|None, reason)."""
    records, reason = _face_records(solid)
    if records is None:
        return None, reason
    diag = _bbox_diagonal(solid)
    tol = REL_TOL * max(diag, 1.0)
    try:
        try:
            cand = _match_box(records, tol)
        except Declined as box_declined:
            # All planar and not a box: the prism family. It runs strictly *after* the box
            # matcher and only on what that declines, so no part recognised today changes tier
            # or candidate.
            try:
                return _match_prism(solid, records, tol, diag), None
            except Declined as prism_declined:
                raise Declined(f"{box_declined}; as a prism: {prism_declined}") from None
        if cand is not None:
            return cand, None
        cand = _match_sphere(records, tol)
        if cand is not None:
            return cand, None
        clusters = _cluster_axial(records, tol)
        if not clusters:
            raise Declined("no cylindrical or conical face to key on")
        # The cluster count is checked before the planes are assigned, so that a part well
        # outside the scope is reported by its structure ("7 axis clusters") rather than by
        # whichever plane happened to fail first.
        if len(clusters) > 2:
            raise Declined(f"{len(clusters)} axis clusters: beyond the recogniser's scope "
                           "(Tier 3 territory, deliberately not built)")
        caps, wedges = _split_planes(records, clusters, tol)
        if len(clusters) == 1:
            # The revolved matcher runs strictly *after* the whole-part primitive one and only on
            # what that declines, so no part that is recognised today changes tier or candidate.
            try:
                return _match_axial_primitive(records, clusters, caps, wedges, tol), None
            except Declined as primitive_declined:
                try:
                    return _match_revolved(solid, records, clusters, caps, wedges,
                                           tol, diag), None
                except Declined as revolved_declined:
                    raise Declined(f"{primitive_declined}; as a revolved profile: "
                                   f"{revolved_declined}") from None
        return _match_two_cluster_union(records, clusters, caps, wedges, tol), None
    except Declined as declined:
        return None, f"{declined} [{_structure(records, tol)}]"


def recognise_revolved(solid):
    """Propose a revolved profile for a solid, skipping the whole-part matchers entirely.

    `recognise()` reaches `_match_revolved` only where `_match_axial_primitive` *declines*. A
    whole-part proposal that is instead **rejected by the acceptance test** never gets there, and
    an all-cone stack is exactly that case: two caps and two cone faces look like one `TGeoCone`
    to the tier-1 matcher, which proposes a cone that is not the solid and is refused by a volume.
    `emit.process_solid` calls this after such a rejection, so the retry costs nothing on a part
    that was accepted the first time and cannot change one.

    Returns `(candidate|None, reason)`, and never raises on a mere mismatch.
    """
    records, reason = _face_records(solid)
    if records is None:
        return None, reason
    diag = _bbox_diagonal(solid)
    tol = REL_TOL * max(diag, 1.0)
    try:
        clusters = _cluster_axial(records, tol)
        if len(clusters) != 1:
            raise Declined(f"{len(clusters)} axis cluster(s): not a single revolved profile")
        caps, wedges = _split_planes(records, clusters, tol)
        return _match_revolved(solid, records, clusters, caps, wedges, tol, diag), None
    except Declined as declined:
        return None, f"{declined} [{_structure(records, tol)}]"


def _structure(records, tol):
    """A one-line structural summary, appended to every decline so the reason is readable."""
    kinds = {}
    for rec in records:
        kinds[rec["kind"]] = kinds.get(rec["kind"], 0) + 1
    try:
        n_clusters = len(_cluster_axial(records, tol))
    except Exception:                                            # noqa: BLE001
        n_clusters = -1
    breakdown = ", ".join(f"{n} {k}" for k, n in sorted(kinds.items()))
    return f"{len(records)} faces: {breakdown}; {n_clusters} axis cluster(s)"
