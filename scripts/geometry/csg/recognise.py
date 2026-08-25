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

from csg import primitives as prim, tier0
from csg.primitives import _add, _cross, _dot, _norm, _scale, _sub, _unit

# Relative tolerance on directions, radii and offsets, scaled by the part's bounding-box
# diagonal. Same value and the same reasoning as the census: CAD arrives with ~1e-7 relative
# agreement between faces the engineer meant to be identical.
REL_TOL = 1.0e-6
ANG_TOL = 1.0e-6


class Declined(Exception):
    """Raised internally with the reason; recognise() turns it into a report entry."""


def _leaf(kind, params, frame, outside=False):
    """`primitives.leaf`, with an illegal *solid* turned into a decline.

    Every proposal in this module goes through here. A description the validators refuse is a
    statement about the part, not a bug in the matcher: ALICE3's fillet blends carry tori whose
    minor radius exceeds their major one, which is a self-intersecting torus and is not a
    `TGeoTorus` at all. Before this wrapper existed that `ValueError` escaped `recognise()` --
    which catches only `Declined` -- and killed the whole conversion on the first blend it met.

    A missing parameter or an unknown leaf type is a bug *here*, stays a plain `ValueError`, and
    is deliberately left to escape: wrapping at this altitude keeps a bad proposal quiet and a
    bad matcher loud.
    """
    try:
        return prim.leaf(kind, params, frame, outside)
    except prim.InvalidDescription as illegal:
        raise Declined(str(illegal)) from None


def _candidate(op, leaves, recogniser, notes=None):
    """`primitives.candidate`, with an illegal description turned into a decline."""
    try:
        return prim.candidate(op, leaves, recogniser, notes)
    except prim.InvalidDescription as illegal:
        raise Declined(str(illegal)) from None


# ------------------------------------------------------------------------------------------
# face analysis
# ------------------------------------------------------------------------------------------

class _LazyScale:
    """`max(bounding-box diagonal, 1 cm)` of a solid, measured on first use and then kept."""

    def __init__(self, solid):
        self._solid = solid
        self._value = None

    @property
    def value(self):
        if self._value is None:
            self._value = max(_bbox_diagonal(self._solid), 1.0)
        return self._value


def _face_records(solid):
    """[{kind, ...carrier..., uv bounds}] for every face, or a reason why the solid is out.

    A face whose *stored* surface is a B-spline may still be exactly a quadric -- that is what a
    CAD exporter does with a cylinder -- so a face no adaptor branch below claims goes to
    `csg/tier0.py` before it is written off as free-form. What comes back is a carrier in this
    same vocabulary plus the measured gap it was accepted on, and it is flagged `canonicalised`
    so that every report downstream can say when a conversion rested on Tier 0.
    """
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.BRepTools import breptools
    from OCC.Core.GeomAbs import (GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane,
                                  GeomAbs_Sphere, GeomAbs_Torus)
    from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods

    records = []
    n_freeform = 0
    n_faces = 0
    n_canonical = 0
    best_declined = None
    # Measured only if a face actually needs canonicalising, so a part whose faces are all
    # natively analytic -- which is every part of every detector corpus -- costs what it did.
    scale = _LazyScale(solid)
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
            to = ad.Torus()
            rec.update(kind="torus", d=_xyz(to.Axis().Direction()),
                       p=_xyz(to.Position().Location()), x=_xyz(to.Position().XDirection()),
                       r=to.MajorRadius(), rt=to.MinorRadius())
        else:
            ellipse = _extruded_ellipse(ad)
            if ellipse is not None:
                rec.update(kind="eltu", **ellipse)
            else:
                canonical, gap = tier0.canonicalise(face, ad, scale.value)
                if canonical is None:
                    n_freeform += 1
                    if gap is not None and (best_declined is None or gap < best_declined):
                        best_declined = gap
                    continue
                rec.update(**canonical)
                if rec["kind"] == "plane" and rec["reversed"]:
                    # The service returns the underlying surface's own normal, exactly as
                    # `ad.Plane().Axis()` does above, so the face's flag is applied here in the
                    # one place that applies it.
                    rec["n"] = _scale(rec["n"], -1.0)
                n_canonical += 1
        records.append(rec)
    if n_freeform:
        how_far = ("" if best_declined is None else
                   f"; the nearest canonical surface any of them proposes is "
                   f"{best_declined:.3g} cm away, {best_declined / scale.value:.3g} of the part, "
                   f"against {tier0.REL_TOL:.0e}")
        rescued = f"; {n_canonical} canonicalised" if n_canonical else ""
        return None, (f"free-form faces: {n_freeform} of {n_faces} "
                      "(surface kind outside plane/cylinder/cone/sphere/torus and not a quadric "
                      "in disguise; a twisted TGeoArb8 side is one of these and is out of "
                      f"scope){rescued}{how_far}")
    if not records:
        return None, "no faces"
    return records, None


def _extruded_ellipse(ad):
    """`{d, p, x, y, a, b}` if this surface is a linear extrusion of an exact ellipse, else None.

    `O2_TGeoToCAD.conv_eltu` builds a `TGeoEltu` as a prism over a `gp_Elips`, so its lateral
    arrives as `GeomAbs_SurfaceOfExtrusion` whose basis curve is `GeomAbs_Ellipse` -- measured, not
    assumed. The writer also swaps the major axis onto y when `a < b`, so what comes back is
    always major >= minor with the major direction saying which of the two the source called `a`;
    reconstructing that is `_eltu_frame`'s job, not this one's. Any other extrusion -- a swept
    B-spline racetrack, for instance -- returns None and stays free-form.
    """
    from OCC.Core.GeomAbs import GeomAbs_Ellipse, GeomAbs_SurfaceOfExtrusion
    if ad.GetType() != GeomAbs_SurfaceOfExtrusion:
        return None
    try:
        basis = ad.BasisCurve()
        if basis.GetType() != GeomAbs_Ellipse:
            return None
        el = basis.Ellipse()
    except Exception:                                            # noqa: BLE001
        return None
    return {"d": _unit(_xyz(ad.Direction())), "p": _xyz(el.Location()),
            "x": _xyz(el.Position().XDirection()), "y": _xyz(el.Position().YDirection()),
            "a": el.MajorRadius(), "b": el.MinorRadius()}


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
    return _candidate("primitive", [_leaf(
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
            return _candidate("primitive", [_leaf(
                "TGeoTubeSeg", {"rmin": rmin, "rmax": rmax, "dz": dz, "phi1": phi1,
                                "phi2": phi2}, frame)], "tier1-tubeseg")
        return _candidate("primitive", [_leaf(
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
        return _candidate("primitive", [_leaf(
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
    return _candidate("primitive", [_leaf(
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
        return _leaf("TGeoTube", {"rmin": 0.5 * (rmin[0] + rmin[1]),
                                      "rmax": 0.5 * (rmax[0] + rmax[1]), "dz": dz},
                         frame), "revolved-tube"
    return _leaf("TGeoCone", {"dz": dz, "rmin1": rmin[0], "rmax1": rmax[0],
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
        lf = _leaf("TGeoPcon", {"phi1": phi1, "dphi": dphi, "z": z, "rmin": rmin,
                                    "rmax": rmax}, frame)
    except Declined as bad:
        raise Declined(f"the reconstructed profile is not a legal TGeoPcon: {bad}") from None

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
    return _candidate("primitive", [lf], tag,
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
        # A degenerate edge -- a cone's apex, a sphere's pole -- has no 3D curve, and PyROOT's
        # binding then returns a 2-tuple rather than the usual (curve, first, last).
        span = BRep_Tool.Curve(edge)
        if span is None or len(span) < 3 or span[0] is None:
            continue
        curve, first, last = span[0], span[1], span[2]
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
                    leaf = _leaf(kind, params, frame)
                except Declined as bad:
                    reasons.append(f"as a {template}: not a legal {kind}: {bad}")
                    continue
                gap = _prism_leaf_gap(leaf, samples)
                if best_gap is None or gap < best_gap:
                    best_gap, best_tag = gap, tag
                if gap <= REL_TOL * scale:
                    return _candidate(
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
        leaves.append(_leaf("TGeoTube", {"rmin": rmin, "rmax": rmax,
                                             "dz": (t1 - t0) / 2.0}, frame))
    return _candidate("union", leaves, "tier2-tube-union",
                      notes={"nCaps": [len(caps[i]) for i in range(len(clusters))]})


# ------------------------------------------------------------------------------------------
# The single cell: one intersection of the part's own halfspaces -> TGeoCompositeShape
# ------------------------------------------------------------------------------------------
#
# `Stream_AA_FlatCSG.md` §5 step 2. A part with **no trusted concave edge** is one cell of the
# arrangement of its own faces' carriers: the intersection of one oriented halfspace per distinct
# carrier, and nothing else. That is not the same as convex -- a tube's bore and a drilled hole
# are halfspaces whose material lies *outside* the carrier, which is why `tube_window` (4
# halfspaces, one of them a hole wall) is one cell while looking nothing like a convex body.
#
# What is emitted, and why it is exact
# ------------------------------------
# One bounded native leaf per carrier, folded into a `TGeoCompositeShape` of intersections and
# subtractions. OCCT cannot build an unbounded halfspace and ROOT's `TGeoHalfSpace` has no OCCT
# counterpart, so a halfspace is realised on BOTH sides as a bounded primitive sized to cover the
# part's bounding box inflated by a generous margin. That is exact, not an approximation, and the
# argument is short:
#
#   * write `B` for the inflated box, `H_i` for the halfspaces and `L_i` for the leaves;
#   * every leaf satisfies `H_i n B  ==  L_i n B` -- an interior leaf is `H_i` clipped to
#     something containing `B`, an exterior leaf's primitive covers all of `H_i`'s complement
#     within `B`;
#   * the fold is `L_0 n L_1 n ... `, and `L_0` is bounded by construction, so the whole fold
#     lies inside `B` and equals `(n H_i) n B`;
#   * the cell is inside the part's own bounding box, hence inside `B`, so the fold is the cell.
#
# The artificial faces therefore never touch the part's neighbourhood, which is the only place
# the symmetric difference looks. Where the single-cell hypothesis is *false* the fold is strictly
# larger than the part and `dV_sym` says so by exactly that much -- which is the falsifier
# `Stream_AA` §5.2 names, and it is left to the volume rather than guessed at here.
#
# The one quantity that decides
# -----------------------------
# `Stream_K_Tier0.md` §3: `_boundary_gap`, the symmetric Hausdorff distance in cm from each
# solid's boundary samples to the *other solid's boundary*, over the bounding-box diagonal. It is
# the same kind of measurement `_profile_gap` and `_point_set_gap` make for the earlier rungs, and
# it is a length over a length. No per-class criterion, no angle. It is load-bearing rather than
# decorative: a V notch shallower than the census's trust filter can see is refused by this
# number and by nothing else (see the notch ladder in `emit.self_test`).

# How far past the part's own bounding box a halfspace leaf is built. A quarter of the diagonal
# puts every artificial face clear of the material while keeping the leaf near the part's size:
# `TGeoIntersection::ComputeBBox` is the running overlap of the operands' boxes, so a leaf sized
# to a sphere of one and a half diagonals -- the first version here -- handed ROOT a composite
# whose box was twice the solid on `oblique_cut_cyl`, and ROOT rejects rays against that box.
_CELL_MARGIN = 0.25

# The budget on what a part may ship as, in boolean leaves -- summed over all its cells, since
# that is the size of the tree ROOT will walk. It replaces the per-cell budget of 8 that rung 3
# carried, which was the right bound while a part was one cell and the wrong one the moment a
# part could be several.
#
# 64 is one doubling past `Stream_P_RepresentationBench.md` §3's measured ladder, taken along
# that ladder's own dead-linear law -- 5.7 ns per leaf on `Contains`, 24-36 ns on
# `DistFromOutside`, measured constant from K=2 to K=32. Extrapolated once: ~365 ns and ~1.5-2.3
# us at 64 leaves, against the 0.4-2.1 us and 2.6-10.8 us those parts measure TODAY as exact
# surface solids (`Stream_AA_FlatCSG.md` §3.3). So the bound is where the composite stops being
# the cheaper of the two representations, and the extrapolation is labelled as one.
#
# It is a dial, and it was measured before it was set: over the five detector corpora, of the 99
# declining parts whose decomposition passes its own boundary gap, a budget of 8 admits 33, 16
# admits 44, 32 admits 78, 64 admits 89 and 128 admits 98. The parts above 64 -- `BREF1` at 282
# leaves, `B077__body` at 122, `VolTOFrail` at 92 -- are precisely the demand `Stream_AA` §5 step
# 4 says a flat `TGeoBVHCSG` has to be justified on, so declining them here is what leaves that
# decision measurable instead of pre-empting it.
_PART_MAX_LEAVES = 64

# At most this many boundary samples per side feed the gap. The samples are strided rather than
# truncated so a part with many edges is still sampled all over.
_CELL_GAP_SAMPLES = 200


def _stride(items, most):
    if len(items) <= most:
        return items
    step = len(items) / float(most)
    return [items[int(i * step)] for i in range(most)]


def _point_to_shape_distance(point, shape):
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
    from OCC.Core.gp import gp_Pnt
    probe = BRepBuilderAPI_MakeVertex(gp_Pnt(*point)).Vertex()
    dist = BRepExtrema_DistShapeShape(probe, shape)
    dist.Perform()
    if not dist.IsDone():
        return float("inf")
    return dist.Value()


def _boundary_gap(a, b, most=_CELL_GAP_SAMPLES):
    """Symmetric Hausdorff distance, in cm, from each solid's boundary samples to the OTHER's
    boundary.

    Deliberately not `_point_set_gap`, which compares one sample *set* against the other. That is
    right for a prism, whose edges are straight and whose midpoints are therefore canonical, and
    wrong here: `tube_window`'s window rim is the transcendental cylinder-cylinder curve, and two
    constructions of the *same* rim parametrise it differently, so their midpoints differ by
    1.6e-03 cm on a curve where the two solids coincide exactly. Measuring a point against the
    other shape rather than against its samples removes that artefact without weakening anything
    -- a sample that is genuinely off the other boundary still reports its true distance.
    """
    samples = (_solid_samples(a), _solid_samples(b))
    if not samples[0] or not samples[1]:
        # One of the two has no boundary at all. In practice that is always the proposal, and it
        # means the halfspaces it was folded from have no common interior -- four pins on four
        # parallel axes read as one cell intersect to nothing. Measured on ALICE3's
        # `ST0923290_01#b12`, which has zero trusted concave edges because its pieces share no
        # edge, not because it is convex. Said here rather than left to the distance below, which
        # would report every sample as unmeasurable and blame OCCT for it.
        raise Declined("the proposal is empty: these carriers have no common interior, so the "
                       "part is not one cell")
    worst = 0.0
    for points, other in ((samples[0], b), (samples[1], a)):
        for point in _stride(points, most):
            distance = _point_to_shape_distance(point, other)
            if not math.isfinite(distance):
                # `BRepExtrema_DistShapeShape` gave up. That is a measurement that did not
                # happen, not a measurement of zero, so it declines and says which.
                raise Declined("OCCT could not measure a boundary sample against the "
                               "proposal, so the gap is unknown")
            worst = max(worst, distance)
    return worst


def _same_carrier(a, b, tol):
    """Do two faces sit on the same oriented carrier surface?"""
    if a["kind"] != b["kind"]:
        return False
    if a["kind"] == "plane":
        return (_collinear(a["n"], b["n"]) and _dot(a["n"], b["n"]) > 0.0
                and abs(_dot(_sub(a["p"], b["p"]), a["n"])) <= tol)
    if a["kind"] == "sphere":
        return _norm(_sub(a["p"], b["p"])) <= tol and abs(a["r"] - b["r"]) <= tol
    if a["kind"] == "torus":
        # Pinned by its centre, its axis, and both radii; two tori of the same R on one axis but
        # different tube radii are the barrel and the bore of a ply and must stay distinct.
        return (_collinear(a["d"], b["d"]) and _norm(_sub(a["p"], b["p"])) <= tol
                and abs(a["r"] - b["r"]) <= tol and abs(a["rt"] - b["rt"]) <= tol)
    if not (_collinear(a["d"], b["d"]) and _on_axis(b["p"], a["p"], a["d"], tol)):
        return False
    if a["kind"] == "cylinder":
        return abs(a["r"] - b["r"]) <= tol
    # A cone is pinned by its apex and its half-angle; the reference radius is chart-dependent.
    return (abs(abs(a["a"]) - abs(b["a"])) <= ANG_TOL
            and _norm(_sub(_cone_apex(a), _cone_apex(b))) <= tol)


def _cone_apex(carrier):
    slope = math.tan(carrier["a"])
    if abs(slope) < 1.0e-30:
        return carrier["p"]
    return _add(carrier["p"], _scale(carrier["d"], -carrier["r"] / slope))


def _halfspace_carriers(solid, tol):
    """The distinct oriented halfspaces of a solid's faces, with the material side of each.

    `census.halfspace_side` is the authority on which side the material is on, and it is reused
    rather than reimplemented: it reads the face's own normal field at the patch centre and
    compares it against the carrier's outward radial direction, and the census cross-checks that
    verdict against `TopAbs_REVERSED` in its own self-test. A rule keyed on the orientation flag
    alone has already been measured to be wrong on revolved solids (see `_match_revolved`).

    A face whose stored surface is not analytic goes through `csg/tier0.py` first, so a cylinder a
    CAD exporter wrote as a B-spline is a carrier here like any other; the side is then decided by
    the same census rule against the *canonical* carrier rather than one read off the adaptor.
    """
    from csg import census
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods

    scale = _LazyScale(solid)
    carriers = []
    exp = TopExp_Explorer(solid, TopAbs_FACE)
    while exp.More():
        face = topods.Face(exp.Current())
        exp.Next()
        ad = BRepAdaptor_Surface(face, True)
        kind = census.SURFACE_TYPE_NAME.get(ad.GetType(), "other")
        canonical = None
        if kind not in ("plane", "cylinder", "cone", "sphere", "torus"):
            canonical, gap = tier0.canonicalise(face, ad, scale.value)
            if canonical is None:
                how_far = ("" if gap is None else
                           f" (the nearest canonical surface it proposes is {gap:.3g} cm away, "
                           f"{gap / scale.value:.3g} of the part)")
                raise Declined(f"a {kind} face is outside the single-cell emitter's "
                               f"carriers{how_far}")
            kind = canonical["kind"]
        rec = {"kind": kind, "side": None}
        if canonical is not None:
            rec.update({k: v for k, v in canonical.items() if k != "uv"})
            if kind == "plane" and face.Orientation() == TopAbs_REVERSED:
                rec["n"] = _scale(rec["n"], -1.0)
        elif kind == "plane":
            axis = ad.Plane().Axis()
            normal = _xyz(axis.Direction())
            if face.Orientation() == TopAbs_REVERSED:
                normal = _scale(normal, -1.0)
            rec.update(n=_unit(normal), p=_xyz(axis.Location()))
        elif kind == "cylinder":
            cy = ad.Cylinder()
            rec.update(d=_unit(_xyz(cy.Axis().Direction())), p=_xyz(cy.Axis().Location()),
                       r=cy.Radius(), x=_xyz(cy.Position().XDirection()))
        elif kind == "cone":
            co = ad.Cone()
            rec.update(d=_unit(_xyz(co.Axis().Direction())), p=_xyz(co.Axis().Location()),
                       r=co.RefRadius(), a=co.SemiAngle(),
                       x=_xyz(co.Position().XDirection()))
        elif kind == "sphere":
            sp = ad.Sphere()
            rec.update(p=_xyz(sp.Location()), r=sp.Radius())
        else:
            to = ad.Torus()
            rec.update(d=_unit(_xyz(to.Axis().Direction())), p=_xyz(to.Position().Location()),
                       x=_xyz(to.Position().XDirection()), r=to.MajorRadius(),
                       rt=to.MinorRadius())
        rec["side"] = (tier0.carrier_side(face, ad, rec) if canonical is not None
                       else census.halfspace_side(face, ad, kind))
        if rec["side"] is None:
            raise Declined(f"a {kind} face's material side could not be decided")
        for existing in carriers:
            if _same_carrier(existing, rec, tol):
                if existing["side"] != rec["side"]:
                    raise Declined("one carrier bounds material on both sides: not one cell")
                break
        else:
            carriers.append(rec)
    if not carriers:
        raise Declined("no faces to read halfspaces from")
    return carriers


def _bbox_of(shape):
    from OCC.Core.Bnd import Bnd_Box
    from OCC.Core.BRepBndLib import brepbndlib
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    box.SetGap(0.0)
    return box.Get()


class _CellBox:
    """The part's bounding box, and the extent every halfspace leaf has to cover.

    Write `B` for this box grown by `margin`. Every leaf is built so that `H_i n B == L_i n B`,
    which is what makes the fold of the leaves equal to the intersection of the halfspaces on the
    part's neighbourhood, and hence equal to the cell.
    """

    def __init__(self, solid, diag):
        xmin, ymin, zmin, xmax, ymax, zmax = _bbox_of(solid)
        self.centre = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax))
        self.corners = [(x, y, z) for x in (xmin, xmax)
                        for y in (ymin, ymax) for z in (zmin, zmax)]
        self.margin = _CELL_MARGIN * max(diag, 1.0)

    def window(self, origin, direction):
        """`[lo, hi]`, measured from `origin` along `direction`, that a leaf must span."""
        reach = [_dot(_sub(corner, origin), direction) for corner in self.corners]
        return min(reach) - self.margin, max(reach) + self.margin


def _cell_leaf(carrier, box):
    """One bounded native leaf covering this halfspace over the part's neighbourhood."""
    outside = carrier["side"] == "exterior"
    if carrier["kind"] == "plane":
        # A box whose +z face lies exactly on the carrier plane and whose body fills the material
        # side. The frame's z points into the material, i.e. against the outward normal.
        normal = carrier["n"]
        into = _scale(normal, -1.0)
        foot = _sub(box.centre, _scale(normal, _dot(_sub(box.centre, carrier["p"]), normal)))
        oriented = prim.frame_from_axis(foot, into)
        depth = max(box.window(foot, into)[1], box.margin)
        half = [max(max(abs(_dot(_sub(c, foot), tuple(oriented[axis]))) for c in box.corners)
                    + box.margin, box.margin) for axis in ("x", "y")]
        frame = dict(oriented)
        frame["origin"] = [float(v) for v in _add(foot, _scale(into, 0.5 * depth))]
        return _leaf("TGeoBBox", {"dx": half[0], "dy": half[1], "dz": 0.5 * depth},
                         frame, outside)
    if carrier["kind"] == "sphere":
        # Already bounded: the halfspace is the ball itself, at its true radius.
        return _leaf("TGeoSphere", {"rmin": 0.0, "rmax": carrier["r"]},
                         prim.identity_frame(carrier["p"]), outside)
    if carrier["kind"] == "torus":
        # Bounded for the same reason as the sphere, and needing no extension: the halfspace
        # "within `rt` of the circle of radius R" *is* the solid torus. A ply's bore is the same
        # carrier with the `outside` flag.
        return _leaf("TGeoTorus", {"r": carrier["r"], "rmin": 0.0, "rmax": carrier["rt"],
                                       "phi1": 0.0, "dphi": 360.0},
                         prim.frame_from_axis(carrier["p"], carrier["d"], carrier["x"]), outside)
    lo, hi = box.window(carrier["p"], carrier["d"])
    if carrier["kind"] == "cylinder":
        frame = prim.frame_from_axis(
            _add(carrier["p"], _scale(carrier["d"], 0.5 * (lo + hi))), carrier["d"],
            carrier["x"])
        return _leaf("TGeoTube", {"rmin": 0.0, "rmax": carrier["r"],
                                      "dz": 0.5 * (hi - lo)}, frame, outside)
    # A cone's halfspace is r <= rref + u tan(a), which is empty beyond the apex, so clipping the
    # window there loses nothing and keeps the second nappe out of the leaf.
    slope = math.tan(carrier["a"])
    if abs(slope) < 1.0e-30:
        raise Declined("a conical carrier with a zero half-angle")
    apex = -carrier["r"] / slope
    lo, hi = (max(lo, apex), hi) if slope > 0.0 else (lo, min(hi, apex))
    if hi - lo <= 0.0:
        raise Declined("a conical carrier whose halfspace does not reach the part")
    frame = prim.frame_from_axis(_add(carrier["p"], _scale(carrier["d"], 0.5 * (lo + hi))),
                                 carrier["d"], carrier["x"])
    return _leaf("TGeoCone", {"dz": 0.5 * (hi - lo), "rmin1": 0.0, "rmin2": 0.0,
                                  "rmax1": max(carrier["r"] + lo * slope, 0.0),
                                  "rmax2": max(carrier["r"] + hi * slope, 0.0)},
                     frame, outside)


def _fold_cell_leaves(carriers, box, tol):
    """One leaf per carrier, except where several carriers already ARE a native primitive.

    Two exact groupings are worth taking, and they are the only two taken here:

      * an interior cylinder or cone capped by two planes perpendicular to its axis is a
        `TGeoTube` / `TGeoCone` with a real `dz`, not three halfspaces;
      * six interior planes in three mutually perpendicular opposite pairs are a `TGeoBBox`, and
        that test is `_match_box` itself rather than a second copy of it.

    Both are equalities, not approximations: `{r <= R} n {z <= z1} n {z >= z0}` *is* the tube.
    The point is the shipped representation -- `tube_window` folds from four halfspaces to
    `TGeoTube - TGeoTube` and `box_minus_cyl` from seven to `TGeoBBox - TGeoTube`, which is one
    boolean node instead of six on parts whose whole reason for existing is query cost
    (`Stream_P` measures a composite's node entry at about 7x a primitive's).
    """
    planes = [c for c in carriers if c["kind"] == "plane"]
    axials = [c for c in carriers if c["kind"] in ("cylinder", "cone")]
    consumed = set()
    leaves = []

    for carrier in axials:
        if carrier["side"] != "interior":
            continue
        open_lo, open_hi = box.window(carrier["p"], carrier["d"])
        ends, span = {}, {}
        for sign, key, opened in ((1.0, "hi", open_hi), (-1.0, "lo", open_lo)):
            here = [pl for pl in planes if id(pl) not in consumed
                    and _collinear(pl["n"], carrier["d"])
                    and sign * _dot(pl["n"], carrier["d"]) > 0.0]
            ends[key] = here[0] if len(here) == 1 else None
            # An end with no cap of its own is left open at the same extent an unfolded
            # halfspace leaf would use, so folding one cap in is still exact.
            span[key] = (_dot(_sub(ends[key]["p"], carrier["p"]), carrier["d"])
                         if ends[key] is not None else opened)
        if ends["hi"] is None and ends["lo"] is None:
            continue
        if span["hi"] - span["lo"] <= tol:
            continue
        folded = _capped_axial_leaf(carrier, span["lo"], span["hi"])
        if folded is None:
            continue
        for key in ("hi", "lo"):
            if ends[key] is not None:
                consumed.add(id(ends[key]))
        consumed.add(id(carrier))
        leaves.append(folded)

    loose = [pl for pl in planes if id(pl) not in consumed]
    if len(loose) == 6:
        try:
            as_box = _match_box(loose, tol)
        except Declined:
            as_box = None
        if as_box is not None:
            leaves.append(as_box["leaves"][0])
            consumed.update(id(pl) for pl in loose)

    for carrier in carriers:
        if id(carrier) in consumed:
            continue
        leaves.append(_cell_leaf(carrier, box))
    if not leaves:
        raise Declined("no halfspace leaf could be built")
    return leaves


def _capped_axial_leaf(carrier, lo, hi):
    """The bounded primitive an interior cylinder or cone plus its two caps already is."""
    frame = prim.frame_from_axis(_add(carrier["p"], _scale(carrier["d"], 0.5 * (lo + hi))),
                                 carrier["d"], carrier["x"])
    dz = 0.5 * (hi - lo)
    if carrier["kind"] == "cylinder":
        return _leaf("TGeoTube", {"rmin": 0.0, "rmax": carrier["r"], "dz": dz}, frame)
    slope = math.tan(carrier["a"])
    rmax1 = carrier["r"] + lo * slope
    rmax2 = carrier["r"] + hi * slope
    if min(rmax1, rmax2) < 0.0:
        return None                     # the apex is between the caps: not one frustum
    return _leaf("TGeoCone", {"dz": dz, "rmin1": 0.0, "rmin2": 0.0,
                                  "rmax1": rmax1, "rmax2": rmax2}, frame)


def _leaf_bbox_volume(lf):
    """A ranking key only: the leaf's own box, used to fold the tightest operand first."""
    p = lf["params"]
    if lf["type"] == "TGeoBBox":
        return 8.0 * p["dx"] * p["dy"] * p["dz"]
    if lf["type"] == "TGeoTube":
        return 8.0 * p["rmax"] ** 2 * p["dz"]
    if lf["type"] == "TGeoCone":
        return 8.0 * max(p["rmax1"], p["rmax2"]) ** 2 * p["dz"]
    return 8.0 * p["rmax"] ** 3


def _cell_leaves(solid, tol, diag, whole_part=True):
    """The ordered halfspace leaves of one cell, or a `Declined` saying why it is not one.

    Shared by the single-cell matcher, where the cell is the whole part, and by the multi-cell
    one, where it is one terminal piece of the decomposition. `whole_part` carries the two guards
    that are statements about a *part* and not about a cell: an all-planar body belongs to the
    prism family's templates, and a body with one carrier is not a composite. Neither is true of
    a piece -- an L-plate splits into two boxes, and a sphere cut off a body is one carrier and a
    perfectly good cell.
    """
    carriers = _halfspace_carriers(solid, tol)
    if whole_part and len(carriers) < 2:
        raise Declined(f"{len(carriers)} distinct carrier(s): not a composite")
    if whole_part and all(c["kind"] == "plane" for c in carriers):
        # An all-planar body is the prism family's, and rung 2 has real templates for it. Read as
        # halfspaces it would ship as a pile of boolean nodes where a TGeoXtru or a TGeoArb8 says
        # the same solid in one -- so the boundary between the two matchers is drawn here rather
        # than left to whichever runs first.
        raise Declined(f"{len(carriers)} planar carriers and nothing else: an all-planar solid "
                       "belongs to the prism family, not to the cell emitter")
    box = _CellBox(solid, diag)

    inside_leaves, outside_leaves = [], []
    for lf in _fold_cell_leaves(carriers, box, tol):
        (outside_leaves if lf.get("outside") else inside_leaves).append(lf)
    if not inside_leaves:
        raise Declined("every carrier's material lies outside it: the cell is unbounded")
    # Tightest first, so `TGeoIntersection::ComputeBBox`'s running overlap starts small and the
    # emitted composite reports a bounding box of the part's own size.
    inside_leaves.sort(key=_leaf_bbox_volume)
    return inside_leaves + outside_leaves, carriers, outside_leaves


def _match_single_cell(solid, records, tol, diag):
    """One intersection cell of the part's own halfspaces: a `TGeoCompositeShape`."""
    from csg import census
    counts = census.edge_census(solid)
    trusted = (counts["concave"] + counts["mixed"]
               - counts["concaveNearTangential"] - counts["mixedNearTangential"])
    if trusted:
        raise Declined(f"{trusted} trusted concave edge(s) of {counts['edges']}: the part is "
                       "more than one cell")
    if counts["nonManifold"] or counts["error"]:
        raise Declined(f"{counts['nonManifold']} non-manifold and {counts['error']} undecidable "
                       "edge(s): the cell test cannot be trusted here")

    leaves, carriers, outside_leaves = _cell_leaves(solid, tol, diag)
    # The fold can collapse the whole cell into one native primitive -- six coplanar-face
    # carriers that are a box, one cylinder and its two caps that are a tube. That is not an
    # intersection of anything and is described as the primitive it is; the tag says which of the
    # two happened, because the two answer differently to the query cost this stream exists for.
    if len(leaves) > _PART_MAX_LEAVES:
        raise Declined(f"the cell is {len(leaves)} halfspaces wide, over the part budget of "
                       f"{_PART_MAX_LEAVES}: it would ship as a boolean tree that deep")
    op = "primitive" if len(leaves) == 1 else "intersection"
    cand = _candidate(op, leaves,
                          "cell-primitive" if op == "primitive" else "cell-intersection",
                          notes={"nCarriers": len(carriers),
                                 "nOutside": len(outside_leaves),
                                 "concaveEdgesTrusted": trusted,
                                 "nLeaves": len(leaves),
                                 "marginDiagonals": _CELL_MARGIN})

    # The one measured quantity. Built here rather than left to the acceptance test because a
    # proposal that does not even build is a decline, not a rejection.
    try:
        realised = prim.build_occ(cand)
    except Exception as exc:                                     # noqa: BLE001
        raise Declined(f"the cell did not build in OCCT: {exc}") from None
    gap = _boundary_gap(solid, realised)
    scale = max(diag, 1.0)
    if gap > REL_TOL * scale:
        raise Declined(f"the cell's boundary is {gap:.3g} cm from the part's "
                       f"({gap / scale:.3g} of the part's {diag:.6g} cm diagonal, "
                       f"over {REL_TOL:.0e})")
    cand["notes"]["cellGapCm"] = gap
    cand["notes"]["cellGapRelative"] = gap / scale
    return cand


# ------------------------------------------------------------------------------------------
# The union of cells: the flat two-level DNF
# ------------------------------------------------------------------------------------------

def _cell_decomposition(solid, scale, max_cells, cache=None):
    """`csg.decompose.split_into_cells` plus the four guards both cell matchers share.

    Computed once per part and memoised on `cache`, because the flat path runs on exactly what
    the tree path declined and re-splitting a body that took a dozen booleans to split is the
    most expensive thing either matcher does. The guards are unchanged and their wording is
    unchanged: a decline a corpus report already carries must keep reading the same.
    """
    from csg import decompose as decomp
    key = ("split", max_cells)
    if cache is not None and key in cache:
        report = cache[key]
    else:
        report = decomp.split_into_cells(solid, max_cells=max_cells, scale=scale)
        if cache is not None:
            cache[key] = report
    if report["stop"]:
        raise Declined(f"the decomposition stopped: {report['stop']} after {report['splits']} "
                       f"split(s) into {len(report['pieces'])} cell(s)")
    if report["unresolved"]:
        raise Declined(f"{len(report['unresolved'])} piece(s) of {len(report['pieces']) + len(report['unresolved'])} "
                       "could not be cut at their own witness edge, so the decomposition is "
                       "incomplete")
    if not report["volumeConserved"]:
        raise Declined(f"the split moved {report['volumeDrift']:.3g} of the part's volume, over "
                       f"{decomp.VOLUME_REL_TOL:.0e}: OCCT's splitter did not conserve it and "
                       "the decomposition is not the part")
    return report


def _match_union_of_cells(solid, records, tol, diag, max_cells=None, max_leaves=None,
                          cache=None):
    """The part decomposed into cells and emitted as their union.

    `Stream_AA_FlatCSG.md` §5 step 3, and the last rung before a flat solid would be built for
    speed rather than for coverage. The decomposition is `csg/decompose.py` -- the probe's own
    measured loop -- and every terminal piece goes through the *existing* single-cell machinery,
    so a cell here is the same object the single-cell emitter has been shipping since rung 3 and
    nothing about how a cell is read is new.

    Four things can go wrong and all four decline rather than ship:

      * the decomposition hits a budget, or a piece cannot be cut at its own witness edge;
      * the split does not conserve volume -- OCCT's splitter is a tolerant boolean core and not
        an exact arrangement, and `Stream_AA` §3.2 measured two ALICE3 parts drifting above their
        band. Those are refused here and fall one tier, which is the whole point of measuring it;
      * a piece is not a cell after all, or the part would ship as a tree wider than the budget;
      * the realised union classifies a point differently from the part.

    That last one is not redundant with the boundary gap or with the symmetric difference, and
    rung 4 measured why. Two defects live in this path that neither of those instruments can
    see: the splitter can hand back a piece that is not a subset of the part, and a cell can
    claim material its own piece does not while its boundary still hugs the piece's -- and on
    the part where both happen, `BRepAlgoAPI_Cut` returns zero solids in both directions, so the
    symmetric difference reports 0.0 against a band of 6.7e-04. `Stream_AA` §3.2's honesty notes
    say never to assume this residue small; the corroboration is what stops it being assumed at
    all.
    """
    from csg import decompose as decomp
    scale = max(diag, 1.0)
    max_cells = decomp.PART_MAX_CELLS if max_cells is None else max_cells
    max_leaves = _PART_MAX_LEAVES if max_leaves is None else max_leaves
    report = _cell_decomposition(solid, scale, max_cells, cache)
    pieces = report["pieces"]
    if len(pieces) < 2:
        raise Declined(f"the decomposition is {len(pieces)} piece(s): not a union of cells")

    cells, total_leaves, n_carriers, n_outside = [], 0, 0, 0
    for index, piece in enumerate(pieces):
        piece_diag = decomp.bbox_diagonal(piece)
        try:
            leaves, carriers, outside = _cell_leaves(piece, tol, piece_diag, whole_part=False)
        except Declined as declined:
            raise Declined(f"cell {index + 1} of {len(pieces)}: {declined}") from None
        total_leaves += len(leaves)
        n_carriers += len(carriers)
        n_outside += len(outside)
        if total_leaves > max_leaves:
            raise Declined(f"{len(pieces)} cells of {total_leaves}+ halfspaces in total, over "
                           f"the part budget of {max_leaves}: it would ship as a boolean tree "
                           "that wide")
        cells.append(_cell(len(leaves), leaves, index, len(pieces)))

    cand = _union_of_cells(cells, "cells-union",
                           notes={"nCells": len(cells),
                                  "nComponents": report["components"],
                                  "nSplits": report["splits"],
                                  "nLeaves": total_leaves,
                                  "nCarriers": n_carriers,
                                  "nOutside": n_outside,
                                  "cellLeaves": [len(c["leaves"]) for c in cells],
                                  "volumeDriftRelative": report["volumeDrift"],
                                  "marginDiagonals": _CELL_MARGIN})
    gap = _measured_gap(solid, cand, diag, "the union of cells")

    # The containment corroboration, and it is not belt-and-braces. Both defects this rung
    # measured are invisible to a boundary distance AND to the symmetric difference:
    # `BRepAlgoAPI_Splitter` can hand back a piece that is not a subset of the part (OCCT's
    # tolerant boolean core, `Stream_AA_FlatCSG.md` §3.2's first honesty note), and a cell can
    # claim material its own piece does not while its boundary still hugs the piece's. On
    # ALICE3's `ST0923290_01#b19` both happen, `BRepAlgoAPI_Cut` returns zero solids in both
    # directions, and `dV_sym` reports 0.0 against a band of 6.7e-04. A classification sees it.
    disagreements, scored, worst = accept_module().contains_disagreements(
        solid, prim.build_occ(cand), accept_module().model_tolerance_cm(solid))
    if disagreements:
        raise Declined(f"the union of cells disagrees with the part about {disagreements} of "
                       f"{scored} classified point(s), the farthest {worst:.3g} cm from the "
                       "part's boundary: the decomposition is not the part, whatever the "
                       "symmetric difference says")
    cand["notes"]["cellGapCm"] = gap
    cand["notes"]["cellGapRelative"] = gap / scale
    cand["notes"]["containsScored"] = scored
    return cand


def accept_module():
    """`csg.accept`, imported lazily to keep this module importable without pythonOCC."""
    from csg import accept
    return accept


def _cell(n_leaves, leaves, index, total):
    """`primitives.cell`, with an illegal cell turned into a decline naming which cell it was."""
    try:
        return prim.cell("primitive" if n_leaves == 1 else "intersection", leaves)
    except prim.InvalidDescription as illegal:
        raise Declined(f"cell {index + 1} of {total}: {illegal}") from None
    except ValueError as illegal:
        raise Declined(f"cell {index + 1} of {total}: {illegal}") from None


def _union_of_cells(cells, recogniser, notes=None):
    """`primitives.union_of_cells`, with an illegal description turned into a decline."""
    try:
        return prim.union_of_cells(cells, recogniser, notes)
    except prim.InvalidDescription as illegal:
        raise Declined(str(illegal)) from None
    except ValueError as illegal:
        raise Declined(str(illegal)) from None


# ------------------------------------------------------------------------------------------
# The flat DNF: the same cells, shipped as halfspaces in `o2::base::O2FlatCSG`
# ------------------------------------------------------------------------------------------

# The flat path's budgets. Not the tree's: `_PART_MAX_LEAVES` bounds how wide a
# `TGeoCompositeShape` may be, and `O2FlatCSG` is not one. These bound the sidecar and the box
# build instead, and are set an order of magnitude above the measured demand -- the worst part in
# `Handoff_FlatCSG.md`'s R5 table is `IBCYSSFlangeA` at about 59 cells.
_PART_MAX_FLAT_CELLS = 256
_PART_MAX_FLAT_HALFSPACES = 1024

# How far a cell's declared bounding box is grown past the piece's own OCCT box.
#
# The box is a CORRECTNESS obligation (`Design_FlatCSGSolid.md` section 4.2): `O2FlatCSG` builds
# its sub-cell boxes strictly inside it, so a cell reaching past its box makes the shape disagree
# with its own `_Loop` twins out there -- the accelerated queries look only inside boxes and find
# nothing, the twins have no boxes and find material.
#
# **The box handed in here is the CAD piece's own bounding box**, so which of the two is right is
# not open: a point outside it is outside the piece and outside the part, the accelerated answer
# is the correct one, and a cell that reaches out there is LARGER than the part because its
# halfspaces do not close it up. The repair is never to widen the box -- that ships the phantom
# material -- it is to refuse the part, which is what `_flat_box_holds_cell` below does.
#
# The margin is not there to buy room for such a cell. It is there because the piece's box is a
# measurement: the acceptance below establishes that the realised cell's boundary is within
# `REL_TOL * scale = 1e-6 * scale` of the piece's everywhere and that it classifies no sampled
# point differently, and 1e-3 of the part diagonal sits three orders of magnitude above that so
# the declared box cannot clip the cell at its own surface. It also clears `CloseShape`'s
# debug-build face sampler, which probes 1e-6 of the diagonal outside each face.
_FLAT_BOX_MARGIN = 1.0e-3


def _flat_cell_box(piece, margin):
    """The outer bound of one cell: the piece's own OCCT box, grown by `margin` on every axis."""
    xmin, ymin, zmin, xmax, ymax, zmax = _bbox_of(piece)
    lo = [xmin - margin, ymin - margin, zmin - margin]
    hi = [xmax + margin, ymax + margin, zmax + margin]
    if not all(math.isfinite(v) for v in lo + hi):
        raise Declined("the cell's bounding box is not finite, so nothing can say where the "
                       "cell ends")
    return lo, hi


# The outward probe that checks the box actually holds the cell. A 3x3 grid on each of the six
# faces, pushed out by each of these multiples of the box's own diagonal -- the small one is
# `CloseShape`'s debug sampler moved onto the release path, the large ones follow a cell that
# leaves through a face and only widens (a cone's mirror nappe, a cylinder's far end).
_FLAT_BOX_PROBE_GRID = 3
_FLAT_BOX_PROBE_OFFSETS = (1.0e-6, 0.25, 1.0, 4.0)


def _flat_box_holds_cell(blocks, lo, hi):
    """`Declined` when a sampled point OUTSIDE the declared box is still inside the cell.

    The converter owes `SetCellBBox` an outer bound and `CloseShape` cannot decide containment
    (`Design_FlatCSGSolid.md` section 4.2): it builds no box outside the one it is given, so a
    cell reaching past its box makes `O2FlatCSG` and its own `_Loop` twins answer differently out
    there, which is what design section 6 exists to forbid.

    **A hit means the cell is bigger than the part, and the fix is never a bigger box.** `lo`/`hi`
    is the CAD piece's own bounding box grown by `_FLAT_BOX_MARGIN`, so a probe point outside it
    is outside the piece and outside the part. The accelerated `Contains` -- which says "no
    material" there -- is right; `Contains_Loop` is inventing material, because the cell's
    halfspaces do not close the cell up. Widening the box would silence the twins by *shipping*
    the phantom material. Refusing the part is the only sound response and is what this does.

    This cannot *prove* containment, and does not claim to. What it does is take the one way this
    actually goes wrong -- an unclosed cell that runs off through a face and keeps going -- and
    look for it where it would be: straight out of every face, near and far.
    """
    from csg import flat
    span = [hi[i] - lo[i] for i in range(3)]
    reach = math.sqrt(sum(v * v for v in span))
    if not (reach > 0.0):
        raise Declined("the cell's bounding box has no extent, so no box can hold the cell")
    steps = [[lo[i] + span[i] * (k + 0.5) / _FLAT_BOX_PROBE_GRID
              for k in range(_FLAT_BOX_PROBE_GRID)] for i in range(3)]
    for axis in range(3):
        u, v = (axis + 1) % 3, (axis + 2) % 3
        for face, base in ((0, lo[axis]), (1, hi[axis])):
            direction = -1.0 if face == 0 else 1.0
            for offset in _FLAT_BOX_PROBE_OFFSETS:
                for su in steps[u]:
                    for sv in steps[v]:
                        point = [0.0, 0.0, 0.0]
                        point[axis] = base + direction * offset * reach
                        point[u], point[v] = su, sv
                        if flat.flat_contains(blocks, tuple(point)):
                            raise Declined(
                                f"the cell's halfspaces still hold {offset * reach:.3g} cm past "
                                f"the CAD piece's own bounding box on axis {axis}: they do not "
                                "close the cell up, so the cell is LARGER than the part. Widening "
                                "the declared box would ship that phantom material and would make "
                                "O2FlatCSG disagree with its own _Loop twins; the part is refused "
                                "instead")


def _match_flat_cells(solid, records, tol, diag, max_cells=None, max_halfspaces=None, cache=None):
    """The same decomposition, emitted as signed implicit halfspaces for `o2::base::O2FlatCSG`.

    Rung R5 of `Handoff_FlatCSG.md`, and the reason the rung exists: a part whose cells are
    perfectly good cells and whose *tree* would be 66 `TGeoBoolNode`s wide declines on the union
    path and ships one tier down today. Nothing about the decomposition or about how a cell is
    read changes here -- `decompose.split_into_cells` and `_cell_leaves` are the same calls the
    union path makes -- only what the cells are turned into: `csg/flat.py`'s halfspace blocks
    instead of `_cell_leaf`'s padded native primitives.

    Three things are this matcher's own and none of them is optional:

      * **the cell bounding box.** An intersection of halfspaces does not bound itself, so the
        converter is the only thing that can say where a cell ends, and `SetCellBBox` believing a
        box that does not contain its cell is a silent one-sided loss of material. See
        `_FLAT_BOX_MARGIN`.
      * **`flat.check_cell_box`, per cell, with the box that is written.** An exterior cone
        carrier's stored quadric is the DOUBLE cone, and beyond the apex its mirror nappe carves
        out material that is really there. An interior cone is repaired exactly by
        `blocks_from_carriers`'s apex plane; an exterior one cannot be, because the complement of
        one nappe is a union, so the part is refused. The refusal needs the box, which is why it
        lives here and not in the emitter.
      * **the containment corroboration**, exactly as on the union path. `Stream_AA` measured
        `BRepAlgoAPI_Cut` reporting `IsDone()` with zero solids in BOTH directions on
        `ST0923290_01#b19`, so the symmetric difference read 0.0 against a band of 6.7e-04. The
        flat path builds its OCCT candidate the same way and inherits the defect (design
        section 10).
    """
    from csg import decompose as decomp, flat
    scale = max(diag, 1.0)
    max_cells = _PART_MAX_FLAT_CELLS if max_cells is None else max_cells
    max_halfspaces = _PART_MAX_FLAT_HALFSPACES if max_halfspaces is None else max_halfspaces
    report = _cell_decomposition(solid, scale, decomp.PART_MAX_CELLS, cache)
    pieces = report["pieces"]
    if len(pieces) > max_cells:
        raise Declined(f"{len(pieces)} cells, over the flat part budget of {max_cells} cells: "
                       "the sidecar and the sub-cell box build are sized to what a part is, not "
                       "to what OCCT can split")
    margin = _FLAT_BOX_MARGIN * scale

    cells, occ_cells, total_blocks, n_carriers, n_outside = [], [], 0, 0, 0
    # A one-piece decomposition IS the whole part, so it keeps the whole-part guards
    # `_match_single_cell` applies -- an all-planar body belongs to the prism family's templates
    # and a one-carrier body to the tier-1 ones, and which matcher owns a part is not a question
    # the flat path is entitled to reopen just by running last. Without this, a near-miss prism
    # that the prism family correctly refuses would ship here instead of declining, and
    # `csg/emit.py --self-test`'s negative control for exactly that says so.
    whole_part = len(pieces) == 1
    for index, piece in enumerate(pieces):
        piece_diag = decomp.bbox_diagonal(piece)
        try:
            leaves, carriers, outside = _cell_leaves(piece, tol, piece_diag,
                                                     whole_part=whole_part)
            lo, hi = _flat_cell_box(piece, margin)
            # The obligation, discharged with the SAME box that is written to the sidecar and
            # handed to `SetCellBBox`. `flat.check_cell_box` raises `Declined`; it is caught here
            # only to say which cell, exactly as every other per-cell decline is.
            flat.check_cell_box(carriers, lo, hi)
            blocks = flat.blocks_from_carriers(carriers)
            _flat_box_holds_cell(blocks, lo, hi)
        except Declined as declined:
            raise Declined(f"cell {index + 1} of {len(pieces)}: {declined}") from None
        total_blocks += len(blocks)
        n_carriers += len(carriers)
        n_outside += len(outside)
        if total_blocks > max_halfspaces:
            raise Declined(f"{len(pieces)} cells of {total_blocks}+ halfspaces in total, over "
                           f"the flat part budget of {max_halfspaces} halfspaces")
        cells.append({"blocks": blocks, "volume": decomp_volume(piece),
                      "lo": lo, "hi": hi})
        occ_cells.append(_cell(len(leaves), leaves, index, len(pieces)))

    cand = _flat_cells(cells, "flat-cells",
                       notes={"nCells": len(cells),
                              "nComponents": report["components"],
                              "nSplits": report["splits"],
                              "nHalfspaces": total_blocks,
                              "nCarriers": n_carriers,
                              "nOutside": n_outside,
                              "cellHalfspaces": [len(c["blocks"]) for c in cells],
                              "cellBoxMarginCm": margin,
                              "volumeDriftRelative": report["volumeDrift"],
                              "occCells": occ_cells})
    gap = _measured_gap(solid, cand, diag, "the flat cells")
    disagreements, scored, worst = accept_module().contains_disagreements(
        solid, prim.build_occ(cand), accept_module().model_tolerance_cm(solid))
    if scored <= 0:
        # `accept.contains_disagreements` answers `(0, 0, 0.0)` when it could not sample at all --
        # a null bounding box, an unclassifiable solid -- and that reads exactly like a clean
        # result. Design section 10 makes the corroboration mandatory on this path, so a
        # corroboration that scored nothing is a decline rather than a pass.
        raise Declined("the containment corroboration scored no point at all, so it corroborated "
                       "nothing: the flat cells are not admitted on an empty measurement")
    if disagreements:
        raise Declined(f"the flat cells disagree with the part about {disagreements} of "
                       f"{scored} classified point(s), the farthest {worst:.3g} cm from the "
                       "part's boundary: the decomposition is not the part, whatever the "
                       "symmetric difference says")
    cand["notes"]["cellGapCm"] = gap
    cand["notes"]["cellGapRelative"] = gap / scale
    cand["notes"]["containsScored"] = scored
    return cand


def decomp_volume(piece):
    """The cell's own volume, from OCCT's `GProp` on the piece it came from.

    `O2FlatCSG::Capacity()` is the sum of these, which is exact because the decomposition's
    volume guard has already established that the pieces are disjoint and sum to the part. The
    number is OCCT's, and design section 11 item 4 says to say so wherever it is quoted.
    """
    from csg.census import volume_of
    return abs(volume_of(piece))


def _flat_cells(cells, recogniser, notes=None):
    """`primitives.flat_cells`, with an illegal description turned into a decline."""
    try:
        return prim.flat_cells(cells, recogniser, notes)
    except prim.InvalidDescription as illegal:
        raise Declined(str(illegal)) from None
    except ValueError as illegal:
        raise Declined(str(illegal)) from None


# ------------------------------------------------------------------------------------------
# Tier 1: the elliptic cylinder
# ------------------------------------------------------------------------------------------


def _eltu_frame(centre, axis, major_dir, minor_dir, major, minor):
    """The frame and the `(a, b)` pair to state an elliptic cylinder in.

    An ellipse is the same ellipse under four labellings -- either semi-axis may be `a`, and
    either sense of it may be `x` -- so the choice is free and is made to match the source. The
    writer puts the major axis on y when the source had `a < b`, so blindly pinning `x` to the
    major axis would report a `TGeoEltu(3, 1.5)` rotated a quarter turn where the source wrote
    `TGeoEltu(1.5, 3)`. Preferring the labelling whose frame is the identity gives the source's
    own two numbers back whenever the part sits square in its own frame, which is the case
    `checkKnownSource.py` compares parameters on; off-axis, the major axis takes `x` and the
    placement carries the rest.
    """
    options = [(major_dir, major, minor), (_scale(major_dir, -1.0), major, minor),
               (minor_dir, minor, major), (_scale(minor_dir, -1.0), minor, major)]
    fallback = None
    for ref_x, a, b in options:
        frame = prim.frame_from_axis(centre, axis, ref_x)
        if prim.frame_is_identity_rotation(frame):
            return frame, a, b
        if fallback is None:
            fallback = (frame, a, b)
    return fallback


def _match_eltu(solid, records, tol, diag):
    """One extruded-ellipse lateral between two perpendicular caps: a `TGeoEltu`."""
    laterals = [r for r in records if r["kind"] == "eltu"]
    planes = [r for r in records if r["kind"] == "plane"]
    other = [r for r in records if r["kind"] not in ("eltu", "plane")]
    if other:
        kinds = sorted({r["kind"] for r in other})
        raise Declined(f"an elliptic lateral together with {kinds} is not a whole TGeoEltu")
    axis = laterals[0]["d"]
    snapped = _snap_to_coordinate_axis(axis)
    if snapped is not None and snapped[1] < 0.0:
        axis = _scale(axis, -1.0)
    for lateral in laterals[1:]:
        if not (_collinear(lateral["d"], axis)
                and abs(lateral["a"] - laterals[0]["a"]) <= tol
                and abs(lateral["b"] - laterals[0]["b"]) <= tol
                and _on_axis(lateral["p"], laterals[0]["p"], axis, tol)):
            raise Declined(f"{len(laterals)} elliptic laterals that are not one carrier")
    if len(planes) != 2:
        raise Declined(f"{len(planes)} planar face(s) on an elliptic lateral, expected 2 caps")
    for plane in planes:
        if not _collinear(plane["n"], axis):
            raise Declined("a planar face of an elliptic cylinder is not perpendicular to it")
    caps = sorted(_dot(_sub(plane["p"], laterals[0]["p"]), axis) for plane in planes)
    if caps[1] - caps[0] <= tol:
        raise Declined("the two caps of an elliptic cylinder are coincident")
    centre = _add(laterals[0]["p"], _scale(axis, 0.5 * (caps[0] + caps[1])))
    frame, a, b = _eltu_frame(centre, axis, laterals[0]["x"], laterals[0]["y"],
                              laterals[0]["a"], laterals[0]["b"])
    try:
        lf = _leaf("TGeoEltu", {"a": a, "b": b, "dz": 0.5 * (caps[1] - caps[0])}, frame)
    except Declined as bad:
        raise Declined(f"the elliptic cylinder is not a legal TGeoEltu: {bad}") from None
    cand = _candidate("primitive", [lf], "tier1-eltu",
                      notes={"semiAxisRatio": min(a, b) / max(a, b)})
    gap = _measured_gap(solid, cand, diag, "the elliptic cylinder")
    cand["notes"]["eltuGapCm"] = gap
    cand["notes"]["eltuGapRelative"] = gap / max(diag, 1.0)
    return cand


# ------------------------------------------------------------------------------------------
# Tier 1: the whole torus
# ------------------------------------------------------------------------------------------
#
# A torus carrier used to end the cascade before it began: `_face_records` counted toroidal faces
# and declined the part. It is now a carrier like any other, which unblocks two things at once --
# this whole-part template, and the torus as a halfspace in `_match_single_cell`, which is what
# converts a bellows ply (a few tori, a cylinder or two and some planes).
#
# Nothing that converts today has a toroidal face, so admitting the kind cannot move an accepted
# part. What it can move is a part that declined: its reason changes, or it converts. The corpus
# diffs say which, per part.


def _match_torus(solid, records, tol, diag):
    """All-toroidal laterals on one axis, optionally phi-cut: a `TGeoTorus`."""
    tori = [r for r in records if r["kind"] == "torus"]
    planes = [r for r in records if r["kind"] == "plane"]
    other = [r for r in records if r["kind"] not in ("torus", "plane")]
    if not tori:
        raise Declined("no toroidal face to key on")
    if other:
        kinds = sorted({r["kind"] for r in other})
        raise Declined(f"a torus together with {kinds} is not a whole torus")

    axis = _unit(tori[0]["d"])
    snapped = _snap_to_coordinate_axis(axis)
    if snapped is not None and snapped[1] < 0.0:
        axis = _scale(axis, -1.0)
    centre = tori[0]["p"]
    for t in tori[1:]:
        if not _collinear(t["d"], axis):
            raise Declined("the toroidal faces do not share one axis")
        if _norm(_sub(t["p"], centre)) > tol:
            raise Declined("the toroidal faces are not concentric")
        if abs(t["r"] - tori[0]["r"]) > tol:
            raise Declined(f"{len(tori)} toroidal faces with different major radii")

    minors = _distinct_radii([t["rt"] for t in tori], tol)
    if len(minors) > 2:
        raise Declined(f"{len(minors)} distinct tube radii on one torus, expected 1 or 2")
    rmin = minors[0] if len(minors) == 2 else 0.0
    rmax = minors[-1]

    for plane in planes:
        if not (_perpendicular(plane["n"], axis)
                and abs(_dot(_sub(plane["p"], centre), plane["n"])) <= tol):
            raise Declined("a planar face of a torus is not a wedge through its axis")
    if planes:
        normals = []
        for plane in planes:
            if not any(_collinear(plane["n"], n) for n in normals):
                normals.append(plane["n"])
        if len(normals) > 2:
            raise Declined(f"{len(normals)} distinct half-planes through the torus axis")

    # phi is read only off tori whose own axis runs *with* the frame's, for the reason
    # `_match_revolved` gives: a flipped carrier axis parametrises phi the other way round.
    oriented = [t for t in tori if _parallel(t["d"], axis) and abs(t["rt"] - rmax) <= tol]
    ref_x = None if _snap_to_coordinate_axis(axis) is not None else (
        oriented[0]["x"] if oriented else None)
    frame = prim.frame_from_axis(centre, axis, ref_x)
    if planes:
        if not oriented:
            raise Declined("no toroidal face runs with the axis, so the phi wedge cannot be read")
        lo_phi, hi_phi = _phi_range(oriented, frame)
        phi1, dphi = lo_phi, hi_phi - lo_phi
    else:
        phi1, dphi = 0.0, 360.0

    try:
        lf = _leaf("TGeoTorus", {"r": tori[0]["r"], "rmin": rmin, "rmax": rmax,
                                     "phi1": phi1, "dphi": dphi}, frame)
    except Declined as bad:
        raise Declined(f"the torus is not a legal TGeoTorus: {bad}") from None

    cand = _candidate("primitive", [lf], "tier1-torus",
                      notes={"nTori": len(tori), "nWedges": len(planes)})
    gap = _measured_gap(solid, cand, diag, "the torus")
    cand["notes"]["torusGapCm"] = gap
    cand["notes"]["torusGapRelative"] = gap / max(diag, 1.0)
    return cand


def _measured_gap(solid, cand, diag, what):
    """The one measured quantity for a whole-part proposal, in cm over the part's diagonal.

    Same instrument as the single cell: the symmetric Hausdorff distance from either solid's
    boundary samples to the other solid's boundary. Measured against the *realised* proposal
    rather than against the carrier equations, so a wrong phi wedge, a swapped rmin or a
    transposed pair of semi-axes shows up as a distance rather than as nothing at all.
    """
    try:
        realised = prim.build_occ(cand)
    except Exception as exc:                                     # noqa: BLE001
        raise Declined(f"{what} did not build in OCCT: {exc}") from None
    gap = _boundary_gap(solid, realised)
    scale = max(diag, 1.0)
    if gap > REL_TOL * scale:
        raise Declined(f"{what}'s boundary is {gap:.3g} cm from the part's "
                       f"({gap / scale:.3g} of the part's {diag:.6g} cm diagonal, "
                       f"over {REL_TOL:.0e})")
    return gap


# ------------------------------------------------------------------------------------------
# entry point
# ------------------------------------------------------------------------------------------

def _with_tier0_notes(cand, records):
    """Record on the candidate what Tier-0 canonicalisation the part rested on.

    Attached only where at least one carrier was canonicalised, so a part built entirely from
    natively-analytic faces keeps the candidate it had before this rung, byte for byte -- which is
    what the digest tables in `csg/emit.py --self-test` assert.

    The number is the worst gap over the part's canonicalised *faces*, which is an upper bound on
    the worst over the carriers a given template actually used; a reader who sees it knows the
    conversion is no better than that, which is the direction a bound has to err in.
    """
    canonical = [r for r in records if r.get("canonicalised")]
    if cand is None or not canonical:
        return cand
    cand["notes"]["tier0Faces"] = len(canonical)
    cand["notes"]["tier0WorstGapCm"] = max(r["tier0GapCm"] for r in canonical)
    cand["notes"]["tier0WorstGapRelative"] = max(r["tier0GapRelative"] for r in canonical)
    return cand


def recognise(solid):
    """Propose a CSG description for one leaf solid in cm. Returns (candidate|None, reason)."""
    records, reason = _face_records(solid)
    if records is None:
        return None, reason
    diag = _bbox_diagonal(solid)
    tol = REL_TOL * max(diag, 1.0)
    cand, reason = _cascade(solid, records, tol, diag)
    return _with_tier0_notes(cand, records), reason


def _cascade(solid, records, tol, diag):
    """The matcher ladder itself, in order of increasing generality."""
    try:
        if any(r["kind"] == "eltu" for r in records):
            # Same reasoning as the torus below: an extruded ellipse was a free-form decline
            # before this rung, so no matcher underneath ever saw one.
            return _match_eltu(solid, records, tol, diag), None
        if any(r["kind"] == "torus" for r in records):
            # A part with a toroidal face reached none of the matchers below before this rung,
            # so routing it straight to the torus template and, on a decline, to the cell
            # emitter cannot change any decision they ever made -- and it keeps them from
            # proposing a tube for a body whose torus they cannot see.
            return _match_torus(solid, records, tol, diag), None
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
        # Everything above has declined, which is exactly the condition the single-cell emitter
        # runs under: no part recognised by any earlier matcher can reach it.
        try:
            return _match_single_cell(solid, records, tol, diag), None
        except Declined as cell_declined:
            # And the decomposition runs only on what the single cell declines, so it cannot
            # change a decision any matcher above it ever made. It is the last rung: a part that
            # reaches it has been refused by every whole-part template and by the one-cell read.
            cache = {}
            try:
                return _match_union_of_cells(solid, records, tol, diag, cache=cache), None
            except Declined as union_declined:
                # THE FLAT PATH IS TRIED ONLY HERE, after the tree path has declined, so no part
                # that converts today changes representation: everything above has already
                # refused this part and the alternative to a flat solid is a surface solid or a
                # mesh. The decomposition itself is shared through `cache`, so the flat attempt
                # costs the halfspace emission and the acceptance, not a second split.
                try:
                    return _match_flat_cells(solid, records, tol, diag, cache=cache), None
                except Declined as flat_declined:
                    return None, (f"{declined}; as a single cell: {cell_declined}; as a union of "
                                  f"cells: {union_declined}; as flat cells: {flat_declined} "
                                  f"[{_structure(records, tol)}]")


def recognise_single_cell(solid):
    """Propose one intersection cell for a solid, skipping every earlier matcher.

    The counterpart of `recognise_revolved`: `recognise()` reaches `_match_single_cell` only where
    everything above it *declines*, and a proposal that is instead **rejected by the acceptance
    test** never gets there. `cyl_inter_cyl` and `tube_window` are exactly that case -- two axis
    clusters read as a two-tube union, which is a superset of the part and is refused by a volume.
    `emit.process_solid` calls this after such a rejection.

    Returns `(candidate|None, reason)`, and never raises on a mere mismatch.
    """
    records, reason = _face_records(solid)
    if records is None:
        return None, reason
    diag = _bbox_diagonal(solid)
    tol = REL_TOL * max(diag, 1.0)
    try:
        return _with_tier0_notes(_match_single_cell(solid, records, tol, diag), records), None
    except Declined as declined:
        return None, f"{declined} [{_structure(records, tol)}]"


def recognise_union_of_cells(solid, max_cells=None, max_leaves=None):
    """Propose a union of cells for a solid, skipping every matcher above it.

    The counterpart of `recognise_single_cell`: `recognise()` reaches `_match_union_of_cells`
    only where everything above it *declines*, and a proposal that is instead **rejected by the
    acceptance test** never gets there. `emit.process_solid` calls this after such a rejection,
    so a part whose whole-part proposal was a superset of it -- the case the retry chain exists
    for -- gets the decomposition tried on it too.

    Returns `(candidate|None, reason)`, and never raises on a mere mismatch.
    """
    records, reason = _face_records(solid)
    if records is None:
        return None, reason
    diag = _bbox_diagonal(solid)
    tol = REL_TOL * max(diag, 1.0)
    try:
        return _with_tier0_notes(
            _match_union_of_cells(solid, records, tol, diag, max_cells, max_leaves),
            records), None
    except Declined as declined:
        return None, f"{declined} [{_structure(records, tol)}]"


def recognise_flat_cells(solid, max_cells=None, max_halfspaces=None):
    """Propose a flat halfspace DNF for a solid, skipping every matcher above it.

    The counterpart of `recognise_union_of_cells`, and it is reached the same two ways: from
    `recognise()`'s cascade, strictly after the union path declines, and from
    `emit.process_solid`'s retry chain after a whole-part proposal is *rejected* by the
    acceptance test rather than declined by a matcher.

    Returns `(candidate|None, reason)`, and never raises on a mere mismatch.
    """
    records, reason = _face_records(solid)
    if records is None:
        return None, reason
    diag = _bbox_diagonal(solid)
    tol = REL_TOL * max(diag, 1.0)
    try:
        return _with_tier0_notes(
            _match_flat_cells(solid, records, tol, diag, max_cells, max_halfspaces),
            records), None
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
        return _with_tier0_notes(
            _match_revolved(solid, records, clusters, caps, wedges, tol, diag), records), None
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
    canonical = [r for r in records if r.get("canonicalised")]
    tier0_note = ""
    if canonical:
        worst = max(r["tier0GapRelative"] for r in canonical)
        tier0_note = (f"; {len(canonical)} canonicalised at a worst gap of {worst:.3g} "
                      "of the part")
    return f"{len(records)} faces: {breakdown}; {n_clusters} axis cluster(s){tier0_note}"
