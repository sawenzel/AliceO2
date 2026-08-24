"""The intermediate CSG description, and the two builders that realise it.

Stream H. A recognised part is described here as data — a small JSON-serialisable tree of
*placed primitives* — and exactly one description is then realised twice:

  * `build_occ()`   -> a `TopoDS_Shape`, which is what the OCCT symmetric-difference acceptance
                       test measures against the original CAD solid;
  * `build_root()`  -> a `TGeoShape`, which is what `shape_<VOL>_<LID>.root` carries and what the
                       oracle gate scores.

Keeping the two builders behind one description is the whole point: the two acceptance tests of
`CSG_Pipeline.md` §3.5 are only independent if they are fed the same *description* rather than the
same *object*. A transposed rotation in one builder and not the other is then a disagreement
between the tests, which is a finding, rather than a silent pass.

Frames
------
Every leaf carries an explicit orthonormal frame `(origin, x, y, z)` in the part's local frame, in
cm. The primitive's own parameters are stated in that frame with `z` as its axis, which is both
ROOT's convention for `TGeoTube`/`TGeoCone`/`TGeoSphere` and OCCT's for `gp_Ax2`.

The ROOT trap this module absorbs, and how
------------------------------------------
**No ROOT shape class can carry a rigid transform.** `TGeoBBox` has `fOrigin` (a translation and
nothing more); every other primitive is fixed to the origin with its axis along z; and the only
`TGeoShape` in ROOT 6.36 that holds a `TGeoMatrix` is `TGeoCompositeShape`, through its
`TGeoBoolNode`.

That is still true of ROOT. What is no longer true is the conclusion the first version of this
module drew from it — that a placed primitive must therefore be written as a `TGeoCompositeShape`
(the primitive unioned with an identical copy of itself under the same matrix). **The shape is now
emitted in its own canonical frame at the origin and the rigid transform travels beside it**, as a
`placement`: a 3x4 row-major `[R | t]` with `part = R * canonical + t`, i.e. `R`'s columns are the
leaf frame's basis vectors. `build_root()` returns `(shape, placement)`; `placement is None` means
identity, which is what every artefact written before this change means too, so nothing older has
to be rewritten. See `Stream_N_PlacedPrimitives.md`.

Consumers compose it where they used to rely on the shape already being in the part frame:
`shape_<VOL>_<LID>.root` carries a `TGeoHMatrix` under the key `placement`; the harness and the
X-ray benchmark transform their points and rays into the shape's local frame; `geom.C` places the
volume with `partPlacement * shapePlacement`, in that order.

`build_occ()` is unchanged and still builds the solid **in the part frame** — the OCCT
symmetric-difference acceptance therefore keeps measuring the *placed* solid against the original
CAD solid, which is the property that would otherwise have silently changed meaning.

What this buys, measured: a placed primitive is a `TGeoTube`/`TGeoTubeSeg`/`TGeoCone`/... again, so
its `Capacity()` is analytic instead of `TGeoCompositeShape`'s Monte-Carlo sampling, and the gate's
capacity column becomes a real measurement (`capacityComparable=true`) for those parts instead of
being unavailable. Genuine multi-leaf booleans still emit a `TGeoCompositeShape` and still report
`capacityComparable=false`; for those, acceptance remains the symmetric-difference volume
(Stream G §2).
"""

import math

# Two frames are the same frame below this; a pure-translation or identity fast path is only
# taken when the rotation is *exactly* representable, so this is a strict test, not a tolerance
# on the geometry itself.
_IDENTITY_EPS = 1.0e-12

# Below this relative difference a cone's two radii are the same radius, and OCCT wants a
# cylinder rather than a cone. See `_occ_frustum`.
_CONE_DEGENERATE_EPS = 1.0e-12


def identity_frame(origin=(0.0, 0.0, 0.0)):
    return {"origin": [float(c) for c in origin],
            "x": [1.0, 0.0, 0.0], "y": [0.0, 1.0, 0.0], "z": [0.0, 0.0, 1.0]}


def frame_from_axis(origin, axis_z, ref_x=None):
    """An orthonormal right-handed frame with `z` along `axis_z`, `x` along `ref_x` if given."""
    z = _unit(axis_z)
    if ref_x is not None:
        x = _sub(ref_x, _scale(z, _dot(ref_x, z)))
        if _norm(x) < 1.0e-9:
            x = None
        else:
            x = _unit(x)
    else:
        x = None
    if x is None:
        # any vector not parallel to z
        seed = (1.0, 0.0, 0.0) if abs(z[0]) < 0.9 else (0.0, 1.0, 0.0)
        x = _unit(_sub(seed, _scale(z, _dot(seed, z))))
    y = _cross(z, x)
    return {"origin": [float(c) for c in origin], "x": list(x), "y": list(y), "z": list(z)}


def frame_is_identity_rotation(frame):
    return (abs(frame["x"][0] - 1.0) < _IDENTITY_EPS and abs(frame["x"][1]) < _IDENTITY_EPS
            and abs(frame["x"][2]) < _IDENTITY_EPS and abs(frame["y"][1] - 1.0) < _IDENTITY_EPS
            and abs(frame["y"][0]) < _IDENTITY_EPS and abs(frame["y"][2]) < _IDENTITY_EPS
            and abs(frame["z"][2] - 1.0) < _IDENTITY_EPS and abs(frame["z"][0]) < _IDENTITY_EPS
            and abs(frame["z"][1]) < _IDENTITY_EPS)


def frame_is_identity(frame):
    return frame_is_identity_rotation(frame) and all(abs(c) < _IDENTITY_EPS
                                                     for c in frame["origin"])


# ------------------------------------------------------------------------------------------
# tiny vector helpers (this package must not depend on numpy: the converter does not)
# ------------------------------------------------------------------------------------------

def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    n = _norm(a)
    if n == 0.0:
        raise ValueError("cannot normalise a zero vector")
    return (a[0] / n, a[1] / n, a[2] / n)


# ------------------------------------------------------------------------------------------
# the description
# ------------------------------------------------------------------------------------------

LEAF_TYPES = ("TGeoBBox", "TGeoTube", "TGeoTubeSeg", "TGeoCone", "TGeoSphere", "TGeoPcon",
              "TGeoTrd1", "TGeoTrd2", "TGeoArb8", "TGeoXtru", "TGeoPgon", "TGeoTorus", "TGeoEltu")

_REQUIRED_PARAMS = {
    "TGeoBBox": ("dx", "dy", "dz"),
    "TGeoTube": ("rmin", "rmax", "dz"),
    "TGeoTubeSeg": ("rmin", "rmax", "dz", "phi1", "phi2"),
    "TGeoCone": ("dz", "rmin1", "rmax1", "rmin2", "rmax2"),
    "TGeoSphere": ("rmin", "rmax"),
    "TGeoPcon": ("phi1", "dphi"),
    "TGeoTorus": ("r", "rmin", "rmax", "phi1", "dphi"),
    "TGeoEltu": ("a", "b", "dz"),
    "TGeoTrd1": ("dx1", "dx2", "dy", "dz"),
    "TGeoTrd2": ("dx1", "dx2", "dy1", "dy2", "dz"),
    "TGeoArb8": ("dz",),
    "TGeoXtru": (),
    "TGeoPgon": ("phi1", "dphi", "nedges"),
}

# Per-type array-valued parameters, stored as plain lists so the description stays
# JSON-serialisable. By default all arrays of one leaf type share a length -- that is the only
# structural rule the generic half enforces; anything else a type needs it states in
# `_LEAF_VALIDATORS`.
_REQUIRED_ARRAY_PARAMS = {
    "TGeoPcon": ("z", "rmin", "rmax"),
    "TGeoPgon": ("z", "rmin", "rmax"),
    "TGeoArb8": ("vertices",),
    "TGeoXtru": ("x", "y", "z", "xoff", "yoff", "scale"),
}

_MIN_ARRAY_LENGTH = {
    "TGeoPcon": 2,
    "TGeoPgon": 2,
    "TGeoArb8": 16,
}

# The one type whose arrays do *not* all share a length: a `TGeoXtru` is one polygon of `nvert`
# corners swept through `nz` sections, and the two counts are independent. Stated as groups, each
# with its own minimum; a type absent from here keeps the single-group rule above unchanged.
_ARRAY_LENGTH_GROUPS = {
    "TGeoXtru": ((("x", "y"), 3), (("z", "xoff", "yoff", "scale"), 2)),
}


class InvalidDescription(ValueError):
    """The numbers do not describe a legal solid of that class.

    Distinct from a plain `ValueError` on purpose. This one is a statement about the *part* -- a
    fillet blend whose torus has a minor radius larger than its major one is a self-intersecting
    torus and is not a `TGeoTorus`, and a recogniser that runs into that must decline readably.
    A missing parameter, an unknown leaf type or a leaf count that does not fit the op is a bug in
    the caller instead, stays a plain `ValueError`, and is left to escape.
    """


def _validate_eltu(p):
    for key in ("a", "b", "dz"):
        if p[key] <= 0.0:
            raise InvalidDescription(f"TGeoEltu: {key} = {p[key]} is not positive")


def _validate_torus(p):
    if p["r"] <= 0.0:
        raise InvalidDescription(f"TGeoTorus: the major radius {p['r']} is not positive")
    if p["rmax"] <= 0.0:
        raise InvalidDescription(f"TGeoTorus: rmax {p['rmax']} is not positive")
    if p["rmin"] < 0.0:
        raise InvalidDescription(f"TGeoTorus: rmin {p['rmin']} is negative")
    if p["rmin"] >= p["rmax"]:
        raise InvalidDescription(f"TGeoTorus: rmin {p['rmin']} is not below rmax {p['rmax']}")
    if p["rmax"] > p["r"]:
        # A tube radius over the major radius is a self-intersecting torus: ROOT accepts the
        # numbers and the two representations then disagree about the overlap, so it is refused
        # here rather than measured later. CAD fillet blends produce these routinely -- ALICE3
        # carries them on small-radius features -- so the message says what the part is, not just
        # which inequality failed.
        raise InvalidDescription(
            f"TGeoTorus: rmax {p['rmax']} exceeds the major radius {p['r']}, so this is a "
            "self-intersecting torus (a fillet blend) that TGeoTorus cannot state")
    if not 0.0 < p["dphi"] <= 360.0 + 1.0e-9:
        raise InvalidDescription(f"TGeoTorus: dphi {p['dphi']} is not in (0, 360]")


def _validate_pcon(p):
    if not 0.0 < p["dphi"] <= 360.0 + 1.0e-9:
        raise InvalidDescription(f"TGeoPcon: dphi {p['dphi']} is not in (0, 360]")
    z, rmin, rmax = p["z"], p["rmin"], p["rmax"]
    for i in range(len(z)):
        if rmin[i] < 0.0:
            raise InvalidDescription(f"TGeoPcon: rmin[{i}] = {rmin[i]} is negative")
        if rmin[i] > rmax[i]:
            raise InvalidDescription(
                f"TGeoPcon: rmin[{i}] = {rmin[i]} exceeds rmax[{i}] = {rmax[i]}")
    for i in range(1, len(z)):
        if z[i] < z[i - 1]:
            raise InvalidDescription(f"TGeoPcon: z is not non-decreasing at section {i} "
                             f"({z[i]} < {z[i - 1]})")
    if z[-1] <= z[0]:
        raise InvalidDescription("TGeoPcon: the profile has no axial extent")
    for i in range(2, len(z)):
        if z[i] == z[i - 1] == z[i - 2]:
            raise InvalidDescription(f"TGeoPcon: three sections share z = {z[i]}")


def _validate_pgon(p):
    _validate_pcon(p)
    if p["nedges"] < 1 or abs(p["nedges"] - round(p["nedges"])) > 1.0e-9:
        raise InvalidDescription(f"TGeoPgon: nedges {p['nedges']} is not a positive whole number")


def _validate_trd1(p):
    if p["dy"] <= 0.0 or p["dz"] <= 0.0:
        raise InvalidDescription(f"TGeoTrd1: dy {p['dy']} and dz {p['dz']} must both be positive")
    if min(p["dx1"], p["dx2"]) < 0.0 or max(p["dx1"], p["dx2"]) <= 0.0:
        raise InvalidDescription(f"TGeoTrd1: dx1 {p['dx1']}, dx2 {p['dx2']} do not bound a solid")


def _validate_trd2(p):
    if p["dz"] <= 0.0:
        raise InvalidDescription(f"TGeoTrd2: dz {p['dz']} must be positive")
    for a, b in (("dx1", "dx2"), ("dy1", "dy2")):
        if min(p[a], p[b]) < 0.0 or max(p[a], p[b]) <= 0.0:
            raise InvalidDescription(f"TGeoTrd2: {a} {p[a]}, {b} {p[b]} do not bound a solid")


def _validate_arb8(p):
    if p["dz"] <= 0.0:
        raise InvalidDescription(f"TGeoArb8: dz {p['dz']} must be positive")
    if len(p["vertices"]) != 16:
        raise InvalidDescription(f"TGeoArb8: needs 16 vertex coordinates, got {len(p['vertices'])}")
    for half, name in ((p["vertices"][:8], "-dz"), (p["vertices"][8:], "+dz")):
        corners = [(half[2 * i], half[2 * i + 1]) for i in range(4)]
        if len({(round(c[0], 12), round(c[1], 12)) for c in corners}) < 3:
            raise InvalidDescription(
                f"TGeoArb8: the {name} face has fewer than three distinct corners")


def _validate_xtru(p):
    z, scale = p["z"], p["scale"]
    for i in range(1, len(z)):
        if z[i] <= z[i - 1]:
            raise InvalidDescription(f"TGeoXtru: z is not strictly increasing at section {i} "
                             f"({z[i]} <= {z[i - 1]})")
    for i, s in enumerate(scale):
        if s <= 0.0:
            raise InvalidDescription(f"TGeoXtru: scale[{i}] = {s} is not positive")
    corners = {(round(a, 12), round(b, 12)) for a, b in zip(p["x"], p["y"])}
    if len(corners) != len(p["x"]):
        raise InvalidDescription("TGeoXtru: the polygon repeats a corner")


_LEAF_VALIDATORS = {
    "TGeoEltu": _validate_eltu,
    "TGeoTorus": _validate_torus,
    "TGeoPcon": _validate_pcon,
    "TGeoPgon": _validate_pgon,
    "TGeoTrd1": _validate_trd1,
    "TGeoTrd2": _validate_trd2,
    "TGeoArb8": _validate_arb8,
    "TGeoXtru": _validate_xtru,
}


def leaf(kind, params, frame, outside=False):
    """One placed primitive. `outside` marks a halfspace whose material is *outside* it.

    An intersection of halfspaces (`op: "intersection"`) may include a halfspace stated as the
    complement of a bounded primitive -- the bore of a tube, the wall of a drilled hole. ROOT
    writes that as a `TGeoSubtraction` node and OCCT as a `BRepAlgoAPI_Cut`, from this one flag.
    """
    if kind not in LEAF_TYPES:
        raise ValueError(f"unknown leaf type {kind!r}")
    arrays = _REQUIRED_ARRAY_PARAMS.get(kind, ())
    missing = [k for k in _REQUIRED_PARAMS[kind] + arrays if k not in params]
    if missing:
        raise ValueError(f"{kind}: missing parameter(s) {missing}")
    out = {k: float(params[k]) for k in _REQUIRED_PARAMS[kind]}
    for k in arrays:
        out[k] = [float(v) for v in params[k]]
    if arrays:
        groups = _ARRAY_LENGTH_GROUPS.get(kind,
                                          ((arrays, _MIN_ARRAY_LENGTH.get(kind, 1)),))
        for names, want in groups:
            lengths = {len(out[k]) for k in names}
            if len(lengths) != 1:
                raise InvalidDescription(
                    f"{kind}: array parameters {list(names)} have unequal lengths "
                    + ", ".join(f"{k}={len(out[k])}" for k in names))
            n = lengths.pop()
            if n < want:
                raise InvalidDescription(
                    f"{kind}: needs at least {want} of {list(names)}, got {n}")
    validator = _LEAF_VALIDATORS.get(kind)
    if validator is not None:
        validator(out)
    described = {"type": kind, "params": out, "frame": frame}
    if outside:
        # Only written when true, so every leaf recorded before halfspaces existed keeps its
        # bytes and the frozen digests of the self-test stay meaningful.
        described["outside"] = True
    return described


def placement_from_frame(frame):
    """The frame as a 3x4 row-major `[R | t]`, with `part = R * canonical + t`.

    `R`'s **columns** are the frame's basis vectors expressed in the part frame, which is the same
    convention `TGeoRotation::SetMatrix` takes (local -> master) and the same one
    `_root_matrix()` builds. Stated once here so that the JSON description, the `TGeoHMatrix` in
    `shape_<part>.root` and the C++ consumers cannot drift apart.
    """
    x, y, z, o = frame["x"], frame["y"], frame["z"], frame["origin"]
    return [[x[0], y[0], z[0], o[0]],
            [x[1], y[1], z[1], o[1]],
            [x[2], y[2], z[2], o[2]]]


def placement_to_local(placement, point):
    """`R^T (p - t)`: a point in the part frame expressed in the shape's own frame."""
    if placement is None:
        return tuple(float(c) for c in point)
    d = (point[0] - placement[0][3], point[1] - placement[1][3], point[2] - placement[2][3])
    return tuple(sum(placement[r][c] * d[r] for r in range(3)) for c in range(3))


def placement_direction_to_local(placement, direction):
    """`R^T d`: a direction in the part frame expressed in the shape's own frame."""
    if placement is None:
        return tuple(float(c) for c in direction)
    return tuple(sum(placement[r][c] * direction[r] for r in range(3)) for c in range(3))


def placement_for_candidate(cand):
    """The rigid transform `build_root()` will hand back beside the shape, or None for identity.

    Computable without ROOT, which matters: the hook writes `csg_<part>.json` in an interpreter
    that may not have PyROOT, and `emit.py --from-json` completes the `.root` file later. Both
    must agree about the placement, so exactly one function decides it and `build_root()` calls
    this one rather than repeating the rule.
    """
    if cand["op"] != "primitive":
        # A genuine multi-leaf boolean is still a TGeoCompositeShape, whose TGeoBoolNode carries
        # the leaves' matrices itself; the composite is already in the part frame.
        return None
    lf = cand["leaves"][0]
    frame = lf["frame"]
    if frame_is_identity(frame):
        return None
    if lf_is_box(lf) and frame_is_identity_rotation(frame):
        # TGeoBBox carries a pure translation itself, through fOrigin. Leaving it there keeps
        # every artefact written for an axis-aligned box byte-identical to before this change.
        return None
    return placement_from_frame(frame)


def candidate(op, leaves, recogniser, notes=None):
    """A described solid: `primitive`, `union`, or `intersection` (of halfspaces).

    `intersection` is the single-cell form. Its leaves are folded left to right; a leaf marked
    `outside` contributes the complement of its primitive, so the fold is a chain of
    intersections and subtractions and not a separate op. The first leaf sets the region the
    others cut down, so it cannot be a complement -- an intersection of complements alone is
    unbounded and is not a cell of anything.
    """
    if op not in ("primitive", "union", "intersection"):
        raise ValueError(f"unknown op {op!r}")
    if op == "primitive" and len(leaves) != 1:
        raise ValueError("op 'primitive' takes exactly one leaf")
    if op == "union" and len(leaves) < 2:
        raise ValueError("op 'union' takes at least two leaves")
    if op == "intersection":
        if len(leaves) < 2:
            raise ValueError("op 'intersection' takes at least two leaves")
        if leaves[0].get("outside"):
            raise ValueError("op 'intersection': the first leaf cannot be a complement")
    if op != "intersection" and any(lf.get("outside") for lf in leaves):
        raise ValueError(f"op {op!r} has no meaning for a complemented leaf")
    return {"op": op, "leaves": leaves, "recogniser": recogniser, "notes": notes or {}}


def describe(cand):
    """One line, for reports."""
    parts = []
    for lf in cand["leaves"]:
        p = lf["params"]
        if lf["type"] in ("TGeoTube", "TGeoTubeSeg"):
            parts.append(f"{lf['type']}(rmin={p['rmin']:.4g}, rmax={p['rmax']:.4g}, "
                         f"dz={p['dz']:.4g})")
        elif lf["type"] == "TGeoBBox":
            parts.append(f"TGeoBBox({p['dx']:.4g}, {p['dy']:.4g}, {p['dz']:.4g})")
        elif lf["type"] == "TGeoCone":
            parts.append(f"TGeoCone(dz={p['dz']:.4g}, {p['rmin1']:.4g}/{p['rmax1']:.4g} -> "
                         f"{p['rmin2']:.4g}/{p['rmax2']:.4g})")
        elif lf["type"] == "TGeoTrd1":
            parts.append(f"TGeoTrd1(dx {p['dx1']:.4g} -> {p['dx2']:.4g}, dy={p['dy']:.4g}, "
                         f"dz={p['dz']:.4g})")
        elif lf["type"] == "TGeoTrd2":
            parts.append(f"TGeoTrd2(dx {p['dx1']:.4g} -> {p['dx2']:.4g}, "
                         f"dy {p['dy1']:.4g} -> {p['dy2']:.4g}, dz={p['dz']:.4g})")
        elif lf["type"] == "TGeoArb8":
            v = p["vertices"]
            parts.append(f"TGeoArb8(dz={p['dz']:.4g}, x {min(v[0::2]):.4g}..{max(v[0::2]):.4g}, "
                         f"y {min(v[1::2]):.4g}..{max(v[1::2]):.4g})")
        elif lf["type"] == "TGeoXtru":
            parts.append(f"TGeoXtru(nvert={len(p['x'])}, nz={len(p['z'])}, "
                         f"z {p['z'][0]:.4g}..{p['z'][-1]:.4g}, "
                         f"scale {min(p['scale']):.4g}..{max(p['scale']):.4g})")
        elif lf["type"] == "TGeoPgon":
            parts.append(f"TGeoPgon(nedges={int(round(p['nedges']))}, nz={len(p['z'])}, "
                         f"phi1={p['phi1']:.4g}, dphi={p['dphi']:.4g}, "
                         f"z {p['z'][0]:.4g}..{p['z'][-1]:.4g}, "
                         f"rmin {min(p['rmin']):.4g}..{max(p['rmin']):.4g}, "
                         f"rmax {min(p['rmax']):.4g}..{max(p['rmax']):.4g})")
        elif lf["type"] == "TGeoEltu":
            parts.append(f"TGeoEltu(a={p['a']:.4g}, b={p['b']:.4g}, dz={p['dz']:.4g})")
        elif lf["type"] == "TGeoTorus":
            parts.append(f"TGeoTorus(r={p['r']:.4g}, rmin={p['rmin']:.4g}, "
                         f"rmax={p['rmax']:.4g}, phi1={p['phi1']:.4g}, dphi={p['dphi']:.4g})")
        elif lf["type"] == "TGeoPcon":
            parts.append(f"TGeoPcon(nz={len(p['z'])}, phi1={p['phi1']:.4g}, "
                         f"dphi={p['dphi']:.4g}, z {p['z'][0]:.4g}..{p['z'][-1]:.4g}, "
                         f"rmin {min(p['rmin']):.4g}..{max(p['rmin']):.4g}, "
                         f"rmax {min(p['rmax']):.4g}..{max(p['rmax']):.4g})")
        else:
            parts.append(f"TGeoSphere(rmin={p['rmin']:.4g}, rmax={p['rmax']:.4g})")
    if cand["op"] == "union":
        return " u ".join(parts)
    if cand["op"] == "intersection":
        out = [parts[0]]
        for lf, text in zip(cand["leaves"][1:], parts[1:]):
            out.append((" - " if lf.get("outside") else " ^ ") + text)
        return "".join(out)
    return parts[0]


# ------------------------------------------------------------------------------------------
# builder 1: OCCT (the acceptance test's candidate side)
# ------------------------------------------------------------------------------------------

def build_occ(cand):
    """Realise the description as a `TopoDS_Shape` in OCCT. Requires pythonOCC."""
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    leaves = cand["leaves"]
    out = _occ_leaf(leaves[0])
    for lf in leaves[1:]:
        nxt = _occ_leaf(lf)
        if cand["op"] == "union":
            maker, what = BRepAlgoAPI_Fuse, "BRepAlgoAPI_Fuse"
        elif lf.get("outside"):
            maker, what = BRepAlgoAPI_Cut, "BRepAlgoAPI_Cut"
        else:
            maker, what = BRepAlgoAPI_Common, "BRepAlgoAPI_Common"
        op = maker(out, nxt)
        op.Build()
        if not op.IsDone():
            raise RuntimeError(f"{what} failed while building the candidate")
        out = op.Shape()
    return out


def _occ_ax2(frame, along_z=0.0):
    from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt
    o = _add(tuple(frame["origin"]), _scale(tuple(frame["z"]), along_z))
    return gp_Ax2(gp_Pnt(*o), gp_Dir(*frame["z"]), gp_Dir(*frame["x"]))


def _occ_cut(outer, inner):
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
    op = BRepAlgoAPI_Cut(outer, inner)
    op.Build()
    if not op.IsDone():
        raise RuntimeError("BRepAlgoAPI_Cut failed while building the candidate")
    return op.Shape()


def _dedupe_ring(pts, tol=1.0e-12):
    """Drop consecutive duplicates in a closed (r, z) ring, the wrap included.

    Same rule as `O2_TGeoToCAD._dedupe_ring`, and it is what makes a `TGeoPcon` whose profile
    pinches to the axis (rmin = rmax, or a duplicated z plane with no radial jump) build at all:
    the duplicated corner would otherwise enter the wire as a zero-length edge.
    """
    out = []
    for pt in pts:
        if out and abs(pt[0] - out[-1][0]) < tol and abs(pt[1] - out[-1][1]) < tol:
            continue
        out.append(pt)
    while len(out) > 1 and abs(out[0][0] - out[-1][0]) < tol and abs(out[0][1] - out[-1][1]) < tol:
        out.pop()
    return out


def pcon_profile_rz(params, tol=1.0e-12):
    """The closed (r, z) profile of a `TGeoPcon`, outer chain then inner chain reversed.

    Exactly the ring `O2_TGeoToCAD.conv_pcon` revolves, so the candidate this package builds and
    the CAD the writer produced from the same numbers are the same construction, not two
    constructions that happen to agree.
    """
    z, rmin, rmax = params["z"], params["rmin"], params["rmax"]
    nz = len(z)
    outer = [(rmax[i], z[i]) for i in range(nz)]
    if all(r <= tol for r in rmin):
        inner = [(0.0, z[nz - 1]), (0.0, z[0])]
    else:
        inner = [(rmin[i], z[i]) for i in range(nz - 1, -1, -1)]
    return _dedupe_ring(outer + inner, tol)


def _occ_pcon(lf):
    """Revolve the (r, z) profile face: true cone/cylinder/plane faces, nothing tessellated.

    The acceptance test cuts this against the original B-rep, so an approximation could not
    pass -- `BRepPrimAPI_MakeRevol` of a polygonal profile is the exact construction, and it is
    the one the writer uses in the other direction.
    """
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeRevol
    from OCC.Core.gp import gp_Ax1, gp_Dir, gp_Pnt
    p, frame = lf["params"], lf["frame"]
    pts = pcon_profile_rz(p)
    if len(pts) < 3:
        raise ValueError("TGeoPcon: degenerate (r, z) profile "
                         f"({len(pts)} distinct corner(s))")
    # OCCT sweeps from the profile's own half-plane, so the profile is laid out at phi1 and the
    # revolution covers dphi -- the same convention `_occ_leaf` uses for a TGeoTubeSeg.
    phi1 = math.radians(p["phi1"])
    xr = _add(_scale(tuple(frame["x"]), math.cos(phi1)),
              _scale(tuple(frame["y"]), math.sin(phi1)))
    origin, zax = tuple(frame["origin"]), tuple(frame["z"])
    poly = BRepBuilderAPI_MakePolygon()
    for (r, zz) in pts:
        poly.Add(gp_Pnt(*_add(origin, _add(_scale(xr, r), _scale(zax, zz)))))
    poly.Close()
    if not poly.IsDone():
        raise RuntimeError("TGeoPcon: could not build the (r, z) profile wire")
    face = BRepBuilderAPI_MakeFace(poly.Wire())
    if not face.IsDone():
        raise RuntimeError("TGeoPcon: the (r, z) profile is not a valid planar face")
    rev = BRepPrimAPI_MakeRevol(face.Face(), gp_Ax1(gp_Pnt(*origin), gp_Dir(*zax)),
                                math.radians(p["dphi"]))
    rev.Build()
    if not rev.IsDone():
        raise RuntimeError("TGeoPcon: revolution of the (r, z) profile failed")
    return rev.Shape()


# ------------------------------------------------------------------------------------------
# the prism family: Trd1 / Trd2 / Arb8 / Xtru / Pgon
# ------------------------------------------------------------------------------------------
#
# All five are the same construction -- a stack of closed sections, corresponding corner by
# corner, with the corners of section k joined to the corners of section k+1 -- and they differ
# only in how the sections are *stated*. So exactly one function turns a description into rings
# (`prism_rings`), and everything else reads those: `build_occ` sews them into a solid, the
# recogniser's score measures against them, and `build_root` hands ROOT the parameters they came
# from. `O2_TGeoToCAD._prism_from_rings` is the same construction in the other direction, which is
# why the acceptance test can be exact rather than approximate.

_PRISM_TYPES = ("TGeoTrd1", "TGeoTrd2", "TGeoArb8", "TGeoXtru", "TGeoPgon")


def _dedupe_ring3(pts, tol=1.0e-9):
    """Drop consecutive duplicate corners in a closed 3-D ring, the wrap included.

    The same rule and the same tolerance as `O2_TGeoToCAD._dedupe_ring3`, so a section that
    pinches -- a `TGeoPgon` wedge closing on its own axis -- enters the wire as the writer's CAD
    has it, without a zero-length edge.
    """
    out = []
    for q in pts:
        if out and max(abs(q[i] - out[-1][i]) for i in range(3)) < tol:
            continue
        out.append(tuple(float(c) for c in q))
    while len(out) > 1 and max(abs(out[0][i] - out[-1][i]) for i in range(3)) < tol:
        out.pop()
    return out


def _pgon_section_ring(r_apothem, z, phi1_deg, dphi_deg, nedges, full):
    """One `TGeoPgon` section polygon. ROOT's rmin/rmax are inscribed-circle radii.

    Byte-for-byte the convention of `O2_TGeoToCAD._pgon_ring`, which is the only statement of it:
    the laterals are planes at the **apothem** radius, so the circumscribed radius the corners sit
    on is `r / cos(dseg / 2)`.
    """
    dseg = math.radians(dphi_deg) / nedges
    radius = r_apothem / math.cos(dseg / 2.0)
    n = nedges if full else nedges + 1
    return [(radius * math.cos(math.radians(phi1_deg) + k * dseg),
             radius * math.sin(math.radians(phi1_deg) + k * dseg), z) for k in range(n)]


def pgon_rings(params):
    """`(outer_stack, inner_stack|None)` for a `TGeoPgon`, as `conv_pgon` builds them."""
    z, rmin, rmax = params["z"], params["rmin"], params["rmax"]
    phi1, dphi, nedges = params["phi1"], params["dphi"], int(round(params["nedges"]))
    full = abs(dphi - 360.0) < 1.0e-9
    hollow = any(r > 0.0 for r in rmin)
    if hollow and full:
        # An annular section is two disjoint rings, which no single wire can express: the outer
        # and the inner prism are separate stacks and the caps are annular.
        return ([_pgon_section_ring(rmax[i], z[i], phi1, dphi, nedges, True)
                 for i in range(len(z))],
                [_pgon_section_ring(max(rmin[i], 0.0), z[i], phi1, dphi, nedges, True)
                 for i in range(len(z))])
    rings = []
    for i in range(len(z)):
        outer = _pgon_section_ring(rmax[i], z[i], phi1, dphi, nedges, full)
        if hollow:
            inner = _pgon_section_ring(max(rmin[i], 0.0), z[i], phi1, dphi, nedges, full)
            rings.append(outer + list(reversed(inner)))
        elif full:
            rings.append(outer)
        else:
            rings.append(outer + [(0.0, 0.0, z[i])])
    return rings, None


def prism_rings(lf):
    """`(outer_stack, inner_stack|None)`: the leaf's sections, in the leaf's own frame.

    Every ring is a closed polygon in corner order, and corner `i` of section `k` is joined to
    corner `i` of section `k + 1`. Ring lengths agree across the stack by construction.
    """
    kind, p = lf["type"], lf["params"]
    if kind == "TGeoTrd1":
        dx1, dx2, dy, dz = p["dx1"], p["dx2"], p["dy"], p["dz"]
        return ([[(-dx1, -dy, -dz), (dx1, -dy, -dz), (dx1, dy, -dz), (-dx1, dy, -dz)],
                 [(-dx2, -dy, dz), (dx2, -dy, dz), (dx2, dy, dz), (-dx2, dy, dz)]], None)
    if kind == "TGeoTrd2":
        dx1, dx2, dy1, dy2, dz = p["dx1"], p["dx2"], p["dy1"], p["dy2"], p["dz"]
        return ([[(-dx1, -dy1, -dz), (dx1, -dy1, -dz), (dx1, dy1, -dz), (-dx1, dy1, -dz)],
                 [(-dx2, -dy2, dz), (dx2, -dy2, dz), (dx2, dy2, dz), (-dx2, dy2, dz)]], None)
    if kind == "TGeoArb8":
        v, dz = p["vertices"], p["dz"]
        return ([[(v[2 * i], v[2 * i + 1], -dz) for i in range(4)],
                 [(v[8 + 2 * i], v[8 + 2 * i + 1], dz) for i in range(4)]], None)
    if kind == "TGeoXtru":
        x, y, z = p["x"], p["y"], p["z"]
        xoff, yoff, sc = p["xoff"], p["yoff"], p["scale"]
        return ([[(xoff[k] + sc[k] * x[i], yoff[k] + sc[k] * y[i], z[k])
                  for i in range(len(x))] for k in range(len(z))], None)
    if kind == "TGeoPgon":
        return pgon_rings(p)
    raise ValueError(f"{kind} is not a prism-family leaf")


def _to_part(frame, q):
    return _add(tuple(frame["origin"]),
                _add(_scale(tuple(frame["x"]), q[0]),
                     _add(_scale(tuple(frame["y"]), q[1]), _scale(tuple(frame["z"]), q[2]))))


def prism_samples(lf):
    """Every corner and every edge midpoint of a prism-family leaf, in the **part** frame.

    This is what the recogniser's one measured quantity is taken against. Corners alone would
    not do: a hexahedron read out in the wrong corner order has the same eight corners and a
    different solid, and only the edges say so.
    """
    outer, inner = prism_rings(lf)
    frame = lf["frame"]
    out = []
    for stack in (outer, inner):
        if stack is None:
            continue
        rings = [_dedupe_ring3(r) for r in stack]
        for k, ring in enumerate(rings):
            n = len(ring)
            for i, q in enumerate(ring):
                out.append(_to_part(frame, q))
                nxt = ring[(i + 1) % n]
                out.append(_to_part(frame, _scale(_add(q, nxt), 0.5)))
                if k + 1 < len(rings) and len(rings[k + 1]) == n:
                    up = rings[k + 1][i]
                    out.append(_to_part(frame, _scale(_add(q, up), 0.5)))
    return out


def _occ_quad_face(b0, b1, t1, t0, tol=1.0e-7):
    """One lateral patch: planar when its corners are coplanar, ruled when they are not.

    The same rule as `O2_TGeoToCAD._quad_face`, including the Newell area test that drops a patch
    bounding nothing (a `TGeoPgon` section repeated at one z collapses its closure edges onto a
    line). Returning the ruled face rather than refusing it is what lets the description express a
    twisted `TGeoArb8` exactly, even though the recogniser never proposes one.
    """
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
    from OCC.Core.BRepFill import brepfill
    from OCC.Core.gp import gp_Pnt
    pts = _dedupe_ring3([b0, b1, t1, t0])
    if len(pts) < 3:
        return None
    nrm = [0.0, 0.0, 0.0]
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        nrm[0] += (a[1] - b[1]) * (a[2] + b[2])
        nrm[1] += (a[2] - b[2]) * (a[0] + b[0])
        nrm[2] += (a[0] - b[0]) * (a[1] + b[1])
    span = max(_norm(_sub(q, pts[0])) for q in pts[1:])
    if _norm(nrm) <= tol * span * span:
        return None
    if len(pts) == 3:
        return BRepBuilderAPI_MakeFace(_occ_polygon_wire(pts)).Face()
    n = _cross(_sub(b1, b0), _sub(t0, b0))
    nn = _norm(n)
    scale = max(_norm(_sub(b1, b0)), _norm(_sub(t0, b0)), 1.0e-30)
    off = abs(_dot(n, _sub(t1, b0))) / nn if nn > 0.0 else 0.0
    if nn > 1.0e-24 and off <= tol * scale:
        mf = BRepBuilderAPI_MakeFace(_occ_polygon_wire(pts))
        if mf.IsDone():
            return mf.Face()
    e1 = BRepBuilderAPI_MakeEdge(gp_Pnt(*b0), gp_Pnt(*b1)).Edge()
    e2 = BRepBuilderAPI_MakeEdge(gp_Pnt(*t0), gp_Pnt(*t1)).Edge()
    return brepfill.Face(e1, e2)


def _occ_polygon_wire(pts):
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
    from OCC.Core.gp import gp_Pnt
    poly = BRepBuilderAPI_MakePolygon()
    for q in pts:
        poly.Add(gp_Pnt(float(q[0]), float(q[1]), float(q[2])))
    poly.Close()
    if not poly.IsDone():
        raise RuntimeError("prism: could not build a section wire")
    return poly.Wire()


def _occ_prism(lf):
    """Sew a prism-family leaf out of explicit faces -- no tessellation, no approximation."""
    from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeSolid,
                                         BRepBuilderAPI_Sewing)
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.TopoDS import topods
    kind = lf["type"]
    outer, inner = prism_rings(lf)
    frame = lf["frame"]
    stacks = []
    for stack in (outer, inner):
        if stack is None:
            continue
        rings = [_dedupe_ring3([_to_part(frame, q) for q in ring]) for ring in stack]
        nv = len(rings[0])
        if nv < 3 or any(len(r) != nv for r in rings):
            raise ValueError(f"{kind}: sections carry "
                             f"{sorted({len(r) for r in rings})} distinct corner counts")
        stacks.append(rings)
    faces = []
    for rings in stacks:
        nv = len(rings[0])
        for k in range(len(rings) - 1):
            lo, hi = rings[k], rings[k + 1]
            for i in range(nv):
                j = (i + 1) % nv
                face = _occ_quad_face(lo[i], lo[j], hi[j], hi[i])
                if face is not None:
                    faces.append(face)
    for idx in (0, -1):
        mf = BRepBuilderAPI_MakeFace(_occ_polygon_wire(stacks[0][idx]))
        if len(stacks) == 2:
            mf.Add(topods.Wire(_occ_polygon_wire(stacks[1][idx]).Reversed()))
        if not mf.IsDone():
            raise ValueError(f"{kind}: could not build a cap face")
        faces.append(mf.Face())
    extent = max(abs(c) for rings in stacks for r in rings for q in r for c in q) or 1.0
    sew = BRepBuilderAPI_Sewing(1.0e-7 * extent)
    for face in faces:
        sew.Add(face)
    sew.Perform()
    shell = sew.SewedShape()
    if shell is None or shell.IsNull():
        raise ValueError(f"{kind}: sewing the sections produced nothing")
    ms = BRepBuilderAPI_MakeSolid(topods.Shell(shell))
    ms.Build()
    solid = ms.Solid()
    props = GProp_GProps()
    brepgprop.VolumeProperties(solid, props)
    if props.Mass() < 0.0:
        solid = topods.Solid(solid.Reversed())
    return solid


def _occ_eltu(lf):
    """An elliptic cylinder, built exactly as `O2_TGeoToCAD.conv_eltu` builds it.

    An ellipse wire, a planar face on it, and a prism along the axis -- so the acceptance test
    measures an extruded exact ellipse against an extruded exact ellipse, with no tessellation
    anywhere. OCCT's `gp_Elips` wants its major radius first, which is why the frame's x is not
    assumed to be the major axis here even though the recogniser prefers to put it there.
    """
    from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                         BRepBuilderAPI_MakeWire)
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
    from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Elips, gp_Pnt, gp_Vec
    p, frame = lf["params"], lf["frame"]
    base = _sub(tuple(frame["origin"]), _scale(tuple(frame["z"]), p["dz"]))
    if p["a"] >= p["b"]:
        major_dir, major, minor = frame["x"], p["a"], p["b"]
    else:
        major_dir, major, minor = frame["y"], p["b"], p["a"]
    axis = gp_Ax2(gp_Pnt(*base), gp_Dir(*frame["z"]), gp_Dir(*major_dir))
    edge = BRepBuilderAPI_MakeEdge(gp_Elips(axis, major, minor)).Edge()
    face = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(edge).Wire())
    if not face.IsDone():
        raise RuntimeError("TGeoEltu: the ellipse wire is not a valid planar face")
    prism = BRepPrimAPI_MakePrism(face.Face(),
                                  gp_Vec(*_scale(tuple(frame["z"]), 2.0 * p["dz"])))
    prism.Build()
    if not prism.IsDone():
        raise RuntimeError("TGeoEltu: the prism failed")
    return prism.Shape()


def _occ_torus(lf):
    """The torus OCCT already has, built exactly as `O2_TGeoToCAD.conv_torus` builds it.

    Analytic on both sides, so the symmetric difference measures a torus against a torus. The
    hollow case is a Cut of the inner torus, which shares the outer's two wedge planes exactly;
    that is OCCT's fragile coincident-face case, so the inner one is swept a hair further in phi
    when there is a wedge to share. It changes nothing about the solid, because the material
    between the two wedge planes is removed either way.
    """
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
    p, frame = lf["params"], lf["frame"]
    phi1, dphi = math.radians(p["phi1"]), math.radians(p["dphi"])
    xr = _add(_scale(tuple(frame["x"]), math.cos(phi1)),
              _scale(tuple(frame["y"]), math.sin(phi1)))
    rotated = {"origin": frame["origin"], "x": list(xr), "y": frame["y"], "z": frame["z"]}

    def make(minor, sweep):
        maker = BRepPrimAPI_MakeTorus(_occ_ax2(rotated), p["r"], minor, sweep)
        maker.Build()
        if not maker.IsDone():
            raise RuntimeError("BRepPrimAPI_MakeTorus failed while building the candidate")
        return maker.Shape()

    outer = make(p["rmax"], dphi)
    if p["rmin"] > 0.0:
        full = dphi >= 2.0 * math.pi - 1.0e-12
        inner = make(p["rmin"], dphi if full else min(dphi + 1.0e-4, 2.0 * math.pi))
        outer = _occ_cut(outer, inner)
    return outer


def _occ_leaf(lf):
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,
                                      BRepPrimAPI_MakeSphere)
    from OCC.Core.gp import gp_Pnt
    kind, p, frame = lf["type"], lf["params"], lf["frame"]
    if kind == "TGeoTorus":
        return _occ_torus(lf)
    if kind == "TGeoEltu":
        return _occ_eltu(lf)
    if kind == "TGeoPcon":
        return _occ_pcon(lf)
    if kind in _PRISM_TYPES:
        return _occ_prism(lf)
    if kind == "TGeoBBox":
        corner = tuple(frame["origin"])
        for axis, half in (("x", p["dx"]), ("y", p["dy"]), ("z", p["dz"])):
            corner = _sub(corner, _scale(tuple(frame[axis]), half))
        ax2 = _occ_ax2({"origin": list(corner), "x": frame["x"], "y": frame["y"],
                        "z": frame["z"]})
        return BRepPrimAPI_MakeBox(ax2, 2 * p["dx"], 2 * p["dy"], 2 * p["dz"]).Shape()
    if kind in ("TGeoTube", "TGeoTubeSeg"):
        ax2 = _occ_ax2(frame, -p["dz"])
        if kind == "TGeoTubeSeg":
            # OCCT sweeps from the frame's own x direction, so rotate the reference direction to
            # phi1 and sweep by (phi2 - phi1); ROOT states the same wedge as two absolute angles.
            phi1 = math.radians(p["phi1"])
            xr = _add(_scale(tuple(frame["x"]), math.cos(phi1)),
                      _scale(tuple(frame["y"]), math.sin(phi1)))
            rotated = {"origin": frame["origin"], "x": list(xr), "y": frame["y"],
                       "z": frame["z"]}
            ax2 = _occ_ax2(rotated, -p["dz"])
            sweep = math.radians(p["phi2"] - p["phi1"])
            outer = BRepPrimAPI_MakeCylinder(ax2, p["rmax"], 2 * p["dz"], sweep).Shape()
            if p["rmin"] > 0.0:
                inner = BRepPrimAPI_MakeCylinder(_occ_ax2(rotated, -p["dz"] - _pad(p["dz"])),
                                                 p["rmin"], 2 * p["dz"] + 4 * _pad(p["dz"]),
                                                 sweep).Shape()
                outer = _occ_cut(outer, inner)
            return outer
        outer = BRepPrimAPI_MakeCylinder(ax2, p["rmax"], 2 * p["dz"]).Shape()
        if p["rmin"] > 0.0:
            # The inner cylinder is deliberately longer than the outer one: a Cut against two
            # exactly coincident planar caps is OCCT's classic fragile case, and the extension
            # changes nothing about the resulting solid.
            pad = _pad(p["dz"])
            inner = BRepPrimAPI_MakeCylinder(_occ_ax2(frame, -p["dz"] - pad), p["rmin"],
                                             2 * p["dz"] + 2 * pad).Shape()
            outer = _occ_cut(outer, inner)
        return outer
    if kind == "TGeoCone":
        outer = _occ_frustum(_occ_ax2(frame, -p["dz"]), p["rmax1"], p["rmax2"], 2 * p["dz"])
        if p["rmin1"] > 0.0 or p["rmin2"] > 0.0:
            pad = _pad(p["dz"])
            slope = (p["rmin2"] - p["rmin1"]) / (2 * p["dz"])
            inner = _occ_frustum(_occ_ax2(frame, -p["dz"] - pad),
                                 max(p["rmin1"] - slope * pad, 0.0),
                                 max(p["rmin2"] + slope * pad, 0.0),
                                 2 * p["dz"] + 2 * pad)
            outer = _occ_cut(outer, inner)
        return outer
    if kind == "TGeoSphere":
        o = tuple(frame["origin"])
        outer = BRepPrimAPI_MakeSphere(gp_Pnt(*o), p["rmax"]).Shape()
        if p["rmin"] > 0.0:
            inner = BRepPrimAPI_MakeSphere(gp_Pnt(*o), p["rmin"]).Shape()
            outer = _occ_cut(outer, inner)
        return outer
    raise ValueError(f"unhandled leaf type {kind!r}")


def _occ_frustum(ax2, r1, r2, height):
    """A cone frustum, or a cylinder when its two radii are the same.

    `BRepPrimAPI_MakeCone` raises `Standard_DomainError("cone with two identic radii")` rather
    than degenerating gracefully, and a `TGeoCone` with `rmax1 == rmax2` (a cylindrical barrel
    with a conical bore) or `rmin1 == rmin2` (a conical barrel with a cylindrical bore) is a
    perfectly ordinary shape -- both occur in ABSO. The non-degenerate call is unchanged, so
    nothing that builds today builds differently.
    """
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder
    if abs(r1 - r2) <= _CONE_DEGENERATE_EPS * max(abs(r1), abs(r2), 1.0):
        return BRepPrimAPI_MakeCylinder(ax2, 0.5 * (r1 + r2), height).Shape()
    return BRepPrimAPI_MakeCone(ax2, r1, r2, height).Shape()


def _pad(dz):
    return max(1.0e-3 * dz, 1.0e-6)


# ------------------------------------------------------------------------------------------
# builder 2: ROOT (what shape_<part>.root carries)
# ------------------------------------------------------------------------------------------

def build_root(cand, name="shape"):
    """Realise the description as `(TGeoShape, placement)`. Requires PyROOT.

    The shape is in its **own canonical frame**: a single recognised primitive comes back as the
    bare `TGeoTube`/`TGeoTubeSeg`/`TGeoCone`/`TGeoSphere`/`TGeoBBox` at the origin with its axis
    along z, and `placement` (3x4 `[R | t]`, or None for identity) says where that sits in the
    part frame. A genuine multi-leaf union is still a `TGeoCompositeShape` already expressed in
    the part frame, and its placement is None.

    Composing the two reproduces exactly what `build_occ()` builds -- that is the invariant the
    self-tests measure, and it is the only reason the two acceptance tests stay independent.
    """
    import ROOT
    placement = placement_for_candidate(cand)
    if cand["op"] == "primitive":
        lf = cand["leaves"][0]
        frame = lf["frame"]
        if placement is None and lf_is_box(lf) and not frame_is_identity(frame):
            # Axis-aligned box: TGeoBBox's own fOrigin is the placement.
            from array import array
            p = lf["params"]
            return ROOT.TGeoBBox(name, p["dx"], p["dy"], p["dz"],
                                 array("d", [float(c) for c in frame["origin"]])), None
        shape = _root_leaf(lf, name)
        return shape, placement
    shapes = [(_root_leaf(lf, f"{name}_l{i}"), lf["frame"])
              for i, lf in enumerate(cand["leaves"])]
    outside = [bool(lf.get("outside")) for lf in cand["leaves"]]
    return _root_composite(name, shapes, cand["op"], outside), placement


def root_placement_matrix(placement, name="placement"):
    """The placement as a `TGeoHMatrix` -- what `shape_<part>.root` carries under key `placement`.

    Returns None for an identity placement, which is what an artefact that records nothing means.
    """
    if placement is None:
        return None
    import ROOT
    # Built through TGeoRotation::SetMatrix / TGeoCombiTrans rather than by poking TGeoHMatrix's
    # arrays directly, because those two set the kGeoRotation / kGeoTranslation bits that decide
    # whether ROOT treats the matrix as anything other than the identity.
    combi = _root_matrix({"x": [placement[0][0], placement[1][0], placement[2][0]],
                          "y": [placement[0][1], placement[1][1], placement[2][1]],
                          "z": [placement[0][2], placement[1][2], placement[2][2]],
                          "origin": [placement[0][3], placement[1][3], placement[2][3]]}, name)
    matrix = ROOT.TGeoHMatrix(combi)
    matrix.SetName(name)
    ROOT.SetOwnership(matrix, False)
    return matrix


def placement_from_root_matrix(matrix):
    """The inverse of `root_placement_matrix()`, for reading an artefact back."""
    if matrix is None:
        return None
    rot = matrix.GetRotationMatrix()
    tr = matrix.GetTranslation()
    return [[rot[0], rot[1], rot[2], tr[0]],
            [rot[3], rot[4], rot[5], tr[1]],
            [rot[6], rot[7], rot[8], tr[2]]]


def lf_is_box(lf):
    return lf["type"] == "TGeoBBox"


def _root_matrix(frame, name):
    import ROOT
    from array import array
    rot = ROOT.TGeoRotation(name + "_r")
    # TGeoRotation::SetMatrix takes the local->master matrix row-major, i.e. the columns are the
    # local frame's basis vectors expressed in the part frame.
    m = array("d", [frame["x"][0], frame["y"][0], frame["z"][0],
                    frame["x"][1], frame["y"][1], frame["z"][1],
                    frame["x"][2], frame["y"][2], frame["z"][2]])
    rot.SetMatrix(m)
    combi = ROOT.TGeoCombiTrans(frame["origin"][0], frame["origin"][1], frame["origin"][2], rot)
    ROOT.SetOwnership(rot, False)
    ROOT.SetOwnership(combi, False)
    return combi


def _root_node_class(op, outside):
    import ROOT
    if op == "union":
        return ROOT.TGeoUnion
    if op == "intersection":
        return ROOT.TGeoSubtraction if outside else ROOT.TGeoIntersection
    raise ValueError(f"unhandled composite op {op!r}")


def _root_composite(name, shapes_and_frames, op, outside=None):
    """Left-fold the leaves into nested boolean nodes.

    The accumulated composite is already expressed in the part frame, so it enters the next node
    with a null matrix; only the fresh leaf carries one. `SetOwnership(.., False)` everywhere is
    not decoration: `TGeoBoolNode` deletes both its operands and both its matrices, so anything
    PyROOT still believes it owns would be freed twice.

    `op` is `union` or `intersection`; under `intersection` a leaf flagged in `outside` enters
    as a `TGeoSubtraction` instead, which is how a halfspace whose material lies outside its own
    primitive is written. The fold order is the leaves' order and nothing here reorders it: the
    composite's bounding box is `TGeoIntersection::ComputeBBox`'s running overlap, so a caller
    that puts its tightest leaf first gets a tight box, and one that does not still gets the
    right solid.
    """
    import ROOT
    flags = list(outside or [False] * len(shapes_and_frames))
    (s0, f0), (s1, f1) = shapes_and_frames[0], shapes_and_frames[1]
    ROOT.SetOwnership(s0, False)
    ROOT.SetOwnership(s1, False)
    node = _root_node_class(op, flags[1])(s0, s1, _root_matrix(f0, f"{name}_m0"),
                                          _root_matrix(f1, f"{name}_m1"))
    ROOT.SetOwnership(node, False)
    comp = ROOT.TGeoCompositeShape(f"{name}_c1", node)
    ROOT.SetOwnership(comp, False)
    for i, (shape, frame) in enumerate(shapes_and_frames[2:], start=2):
        ROOT.SetOwnership(shape, False)
        node = _root_node_class(op, flags[i])(comp, shape, ROOT.nullptr,
                                              _root_matrix(frame, f"{name}_m{i}"))
        ROOT.SetOwnership(node, False)
        comp = ROOT.TGeoCompositeShape(f"{name}_c{i}", node)
        ROOT.SetOwnership(comp, False)
    comp.SetName(name)
    return comp


def _root_leaf(lf, name):
    import ROOT
    kind, p = lf["type"], lf["params"]
    if kind == "TGeoBBox":
        return ROOT.TGeoBBox(name, p["dx"], p["dy"], p["dz"])
    if kind == "TGeoTube":
        return ROOT.TGeoTube(name, p["rmin"], p["rmax"], p["dz"])
    if kind == "TGeoTubeSeg":
        return ROOT.TGeoTubeSeg(name, p["rmin"], p["rmax"], p["dz"], p["phi1"], p["phi2"])
    if kind == "TGeoCone":
        return ROOT.TGeoCone(name, p["dz"], p["rmin1"], p["rmax1"], p["rmin2"], p["rmax2"])
    if kind == "TGeoSphere":
        return ROOT.TGeoSphere(name, p["rmin"], p["rmax"])
    if kind == "TGeoTorus":
        return ROOT.TGeoTorus(name, p["r"], p["rmin"], p["rmax"], p["phi1"], p["dphi"])
    if kind == "TGeoEltu":
        return ROOT.TGeoEltu(name, p["a"], p["b"], p["dz"])
    if kind == "TGeoPcon":
        shape = ROOT.TGeoPcon(name, p["phi1"], p["dphi"], len(p["z"]))
        for i, (zz, r0, r1) in enumerate(zip(p["z"], p["rmin"], p["rmax"])):
            shape.DefineSection(i, zz, r0, r1)
        return shape
    if kind == "TGeoPgon":
        shape = ROOT.TGeoPgon(name, p["phi1"], p["dphi"], int(round(p["nedges"])), len(p["z"]))
        for i, (zz, r0, r1) in enumerate(zip(p["z"], p["rmin"], p["rmax"])):
            shape.DefineSection(i, zz, r0, r1)
        return shape
    if kind == "TGeoTrd1":
        return ROOT.TGeoTrd1(name, p["dx1"], p["dx2"], p["dy"], p["dz"])
    if kind == "TGeoTrd2":
        return ROOT.TGeoTrd2(name, p["dx1"], p["dx2"], p["dy1"], p["dy2"], p["dz"])
    if kind == "TGeoArb8":
        from array import array
        return ROOT.TGeoArb8(name, p["dz"], array("d", [float(v) for v in p["vertices"]]))
    if kind == "TGeoXtru":
        from array import array
        shape = ROOT.TGeoXtru(len(p["z"]))
        shape.SetName(name)
        shape.DefinePolygon(len(p["x"]), array("d", [float(v) for v in p["x"]]),
                            array("d", [float(v) for v in p["y"]]))
        for k in range(len(p["z"])):
            shape.DefineSection(k, p["z"][k], p["xoff"][k], p["yoff"][k], p["scale"][k])
        return shape
    raise ValueError(f"unhandled leaf type {kind!r}")
