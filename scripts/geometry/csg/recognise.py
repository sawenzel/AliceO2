"""Bagger-sufficient CSG recognition: whole-part primitives, and the two-cluster tube union.

Scope, deliberately narrow (`Stream_A_CSG.md` §3, `Tutorial.md` §6 item 2):

  * **Tier 1** — the part's whole face set is one ROOT primitive: `TGeoBBox`, `TGeoTube`,
    `TGeoTubeSeg`, `TGeoCone`, `TGeoSphere`;
  * **Tier 2, the one measured case** — two coaxial *clusters* of cylinders, i.e. a barrel and a
    lug on a second axis, emitted as `TGeoTube u TGeoTube`. The census measured all six Bagger
    "ram" parts as exactly this, and measured them as a **union**, not the difference the older
    plan text names: the bores are `rmin` on the leaves, so no subtraction node is needed.

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
        rec = {"uv": (umin, umax, vmin, vmax),
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
                         "(surface kind outside plane/cylinder/cone/sphere)")
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
        cand = _match_box(records, tol)
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
            return _match_axial_primitive(records, clusters, caps, wedges, tol), None
        return _match_two_cluster_union(records, clusters, caps, wedges, tol), None
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
