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

LEAF_TYPES = ("TGeoBBox", "TGeoTube", "TGeoTubeSeg", "TGeoCone", "TGeoSphere")

_REQUIRED_PARAMS = {
    "TGeoBBox": ("dx", "dy", "dz"),
    "TGeoTube": ("rmin", "rmax", "dz"),
    "TGeoTubeSeg": ("rmin", "rmax", "dz", "phi1", "phi2"),
    "TGeoCone": ("dz", "rmin1", "rmax1", "rmin2", "rmax2"),
    "TGeoSphere": ("rmin", "rmax"),
}


def leaf(kind, params, frame):
    if kind not in LEAF_TYPES:
        raise ValueError(f"unknown leaf type {kind!r}")
    missing = [k for k in _REQUIRED_PARAMS[kind] if k not in params]
    if missing:
        raise ValueError(f"{kind}: missing parameter(s) {missing}")
    return {"type": kind, "params": {k: float(params[k]) for k in _REQUIRED_PARAMS[kind]},
            "frame": frame}


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
    if op not in ("primitive", "union"):
        raise ValueError(f"unknown op {op!r}")
    if op == "primitive" and len(leaves) != 1:
        raise ValueError("op 'primitive' takes exactly one leaf")
    if op == "union" and len(leaves) < 2:
        raise ValueError("op 'union' takes at least two leaves")
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
        else:
            parts.append(f"TGeoSphere(rmin={p['rmin']:.4g}, rmax={p['rmax']:.4g})")
    return (" u ".join(parts)) if cand["op"] == "union" else parts[0]


# ------------------------------------------------------------------------------------------
# builder 1: OCCT (the acceptance test's candidate side)
# ------------------------------------------------------------------------------------------

def build_occ(cand):
    """Realise the description as a `TopoDS_Shape` in OCCT. Requires pythonOCC."""
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
    shapes = [_occ_leaf(lf) for lf in cand["leaves"]]
    out = shapes[0]
    for nxt in shapes[1:]:
        op = BRepAlgoAPI_Fuse(out, nxt)
        op.Build()
        if not op.IsDone():
            raise RuntimeError("BRepAlgoAPI_Fuse failed while building the candidate")
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


def _occ_leaf(lf):
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                                      BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere)
    from OCC.Core.gp import gp_Pnt
    kind, p, frame = lf["type"], lf["params"], lf["frame"]
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
        ax2 = _occ_ax2(frame, -p["dz"])
        outer = BRepPrimAPI_MakeCone(ax2, p["rmax1"], p["rmax2"], 2 * p["dz"]).Shape()
        if p["rmin1"] > 0.0 or p["rmin2"] > 0.0:
            pad = _pad(p["dz"])
            slope = (p["rmin2"] - p["rmin1"]) / (2 * p["dz"])
            inner = BRepPrimAPI_MakeCone(_occ_ax2(frame, -p["dz"] - pad),
                                         max(p["rmin1"] - slope * pad, 0.0),
                                         max(p["rmin2"] + slope * pad, 0.0),
                                         2 * p["dz"] + 2 * pad).Shape()
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
    return _root_composite(name, shapes, cand["op"]), placement


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


def _root_composite(name, shapes_and_frames, op):
    """Left-fold the leaves into nested `TGeoUnion`s.

    The accumulated composite is already expressed in the part frame, so it enters the next node
    with a null matrix; only the fresh leaf carries one. `SetOwnership(.., False)` everywhere is
    not decoration: `TGeoBoolNode` deletes both its operands and both its matrices, so anything
    PyROOT still believes it owns would be freed twice.
    """
    import ROOT
    if op != "union":
        raise ValueError(f"unhandled composite op {op!r}")
    (s0, f0), (s1, f1) = shapes_and_frames[0], shapes_and_frames[1]
    ROOT.SetOwnership(s0, False)
    ROOT.SetOwnership(s1, False)
    node = ROOT.TGeoUnion(s0, s1, _root_matrix(f0, f"{name}_m0"), _root_matrix(f1, f"{name}_m1"))
    ROOT.SetOwnership(node, False)
    comp = ROOT.TGeoCompositeShape(f"{name}_c1", node)
    ROOT.SetOwnership(comp, False)
    for i, (shape, frame) in enumerate(shapes_and_frames[2:], start=2):
        ROOT.SetOwnership(shape, False)
        node = ROOT.TGeoUnion(comp, shape, ROOT.nullptr, _root_matrix(frame, f"{name}_m{i}"))
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
    raise ValueError(f"unhandled leaf type {kind!r}")
