#!/usr/bin/env python3.10
"""Stream N probe: why is a recognised quadric's boundary edge rejected, and could an *implicit*
(co-surface) trim carry it exactly?

The question
------------
`O2_CADtoTGeo.py::_recognized_quadric_wire_block` accepts a boundary edge only if it is iso in the
recognised (phi, other) domain. ALICE3 has 36 surface-eligible leaf solids and emits 20 sidecars;
the 16 missing ones are declined by that one line. `Stream_K_Tier0.md` §2 proposed *fitting* the
rejected edges. This probe measures whether a fit is needed at all, by classifying **every rejected
edge** into exactly one bucket:

  A  iso after all -- iso in an equivalent/alternative parametrisation of the *same* recognised
     surface, or iso within the exporter's own declared edge tolerance and rejected only by the
     shipped test's fixed constants. Split into
       A1  iso in the shipped parametrisation at the level of the edge's own declared tolerance
           (a tolerance/constant issue in the shipped test, not a geometry fact), and
       A2  iso only in an alternative parametrisation. For a sphere the polar axis is free and the
           converter hard-codes (0,0,1); any *planar* edge on a sphere is a constant-theta circle
           about that plane's normal. Cylinders and cones have no such freedom (their axis is
           determined), so A2 is structurally sphere-only and the probe reports it as measured.
  B  exactly the intersection of two recognised analytic surfaces -- the edge is shared with a
     neighbouring face whose surface is also analytic (stored *or* recognised). Reported with the
     measured max distance of dense edge samples from BOTH implicit surfaces, absolutely and
     normalised by the edge length, the face patch diagonal and the edge's own BRep tolerance.
  C  the neighbour face is free-form: no exact co-surface trim exists.
  D  anything else -- named, never lumped (seam edge shared by a face with itself, an edge with
     more than one neighbour, an edge with no 3D curve, ...).

What it does NOT do
-------------------
Nothing here changes the converter, the sidecar or the kernel. It imports `O2_CADtoTGeo` and calls
its *shipping* functions: `extract_graph` (so the leaf solids are the converter's own `def_shapes`
-- a `STEPControl_Reader` load of the same file heals differently and gives a different count),
`_recognize_analytic_surface`, the `_FACE_EXTRACTORS`, and `recognize_and_extract_face`.

Self-checks it runs before reporting anything
---------------------------------------------
1. **Face-verdict agreement.** The probe re-implements the per-edge iso test so it can keep going
   past the first bad edge (the shipped function returns at the first one). Its predicted face
   verdict is compared against the shipped `recognize_and_extract_face` for every face; a
   disagreement is printed and counted. If that count is not 0 the census is not measuring the
   shipped criterion and nothing else in the output means anything.
2. **Bucket-B positive control.** Edges that the shipped test *accepts* (iso rims and generators)
   are measured with exactly the same instrument. A rim circle is the intersection of the
   recognised quadric with its neighbouring cap plane, so the accepted edges must land in B at a
   deviation of the same order. An instrument that only ever sees small numbers on the population
   it is arguing about has not been controlled.

Usage
-----
    export ALIBUILD_WORK_DIR=$HOME/alisw/sw
    B=$HOME/alisw/sw/BUILD/O2-latest/O2
    cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
    export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
    cd $HOME/alisw/O2/scripts/geometry
    python3 probes/trimEdgeCensus.py --model ALICE_3_example/CAD_noETA.stp --label ALICE3 \
        --json /tmp/n/alice3.json

`--label fixtures` is special: it builds the ladder fixtures with `make_boolean_fixtures.py` into a
scratch directory and censuses every one of them.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))

from csg.occ_env import ensure_occ  # noqa: E402

ensure_occ()

import numpy as np  # noqa: E402

from OCC.Core.BRep import BRep_Tool  # noqa: E402
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface  # noqa: E402
from OCC.Core.BRepTools import breptools  # noqa: E402
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED  # noqa: E402
from OCC.Core.TopExp import topexp  # noqa: E402
from OCC.Core.TopTools import (TopTools_IndexedMapOfShape,  # noqa: E402
                               TopTools_IndexedDataMapOfShapeListOfShape,
                               TopTools_ListIteratorOfListOfShape)
from OCC.Core.TopoDS import topods  # noqa: E402
from OCC.Core.Bnd import Bnd_Box  # noqa: E402
from OCC.Core.BRepBndLib import brepbndlib  # noqa: E402
from OCC.Extend.TopologyUtils import TopologyExplorer  # noqa: E402

import O2_CADtoTGeo as C  # noqa: E402


# ---------------------------------------------------------------------------------------------
# analytic surfaces: one model dict per face, one distance function for all of them
# ---------------------------------------------------------------------------------------------

def _vec(p):
    return np.array([p.X(), p.Y(), p.Z()], dtype=float)


def _unit(v):
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


# ---------------------------------------------------------------------------------------------
# the PRE-FIX recogniser, verbatim from `git show 237be7f81a^`, for one purpose only
# ---------------------------------------------------------------------------------------------
# `Stream_K_Tier0.md` §2's edge census (1891 edges, 1053 free-form) was measured before §5's
# acceptance-criterion fix, when 184 extra faces on `ST2487458_01` were still accepted as cones.
# Re-deriving §2's numbers therefore needs the recogniser §2 ran against. This is a *copy* used by
# `--legacy-recognizer` for the reconciliation and by nothing else; the production function is
# never modified and the default path calls it, not this.

def _recognize_analytic_surface_prefix(adaptor, uv_bounds):
    umin, umax, vmin, vmax = uv_bounds
    P, N = C._sample_surface_for_recognition(adaptor, umin, umax, vmin, vmax)
    if P is None:
        return None
    scale = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0)))
    if scale < 1e-12:
        return None
    dev = float(np.abs(np.abs(N @ N[0]) - 1.0).max())
    if dev < 1e-12:
        normal = N[0] / np.linalg.norm(N[0])
        out = {"kind": "plane", "residual": dev, "normal": normal, "point": P[0], "P": P, "N": N}
        out["gap"] = C._analytic_surface_gap("plane", out, P)
        out["gap_relative"] = out["gap"] / scale
        return out
    best = ("freeform", float("inf"), {})
    A = np.zeros((3 * len(P), 4))
    b = np.zeros(3 * len(P))
    for i in range(len(P)):
        A[3 * i:3 * i + 3, 0:3] = np.eye(3)
        A[3 * i:3 * i + 3, 3] = N[i]
        b[3 * i:3 * i + 3] = P[i]
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    centre, radius = sol[:3], abs(sol[3])
    res = float(np.abs(np.linalg.norm(P - centre, axis=1) - radius).max() / scale)
    if res < best[1]:
        best = ("sphere", res, {"centre": centre, "radius": radius})
    _u, _s, Vt = np.linalg.svd(N, full_matrices=False)
    axis = Vt[-1]
    if np.abs(N @ axis).max() < 1e-9:
        e1 = Vt[0]
        e2 = np.cross(axis, e1)
        x, y = P @ e1, P @ e2
        M = np.column_stack([x, y, np.ones_like(x)])
        D, E, F = np.linalg.lstsq(M, -(x ** 2 + y ** 2), rcond=None)[0]
        cx, cy = -D / 2, -E / 2
        r2 = cx * cx + cy * cy - F
        if r2 > 0:
            R = math.sqrt(r2)
            res = float(np.abs(np.hypot(x - cx, y - cy) - R).max() / scale)
            if res < best[1]:
                origin = cx * e1 + cy * e2
                best = ("cylinder", res, {"axis": axis, "refu": e1, "origin": origin, "radius": R})
    apex, *_ = np.linalg.lstsq(N, np.einsum('ij,ij->i', N, P), rcond=None)
    d = P - apex
    dn = np.linalg.norm(d, axis=1)
    ok = dn > 1e-12
    if ok.sum() > 10:
        u = d[ok] / dn[ok, None]
        res = float(np.abs(np.einsum('ij,ij->i', u, N[ok])).max())
        mean_dir = u.mean(axis=0)
        _u2, _s2, Vt2 = np.linalg.svd(u - mean_dir, full_matrices=False)
        ax2 = np.cross(Vt2[0], Vt2[1])
        n2 = np.linalg.norm(ax2)
        if n2 > 1e-12:
            ax2 = ax2 / n2
            if np.dot(mean_dir, ax2) < 0.0:
                ax2 = -ax2
            res = max(res, float(np.abs(u @ ax2).std()))
            if res < best[1]:
                ref = u[0] - np.dot(u[0], ax2) * ax2
                refn = np.linalg.norm(ref)
                if refn > 1e-9:
                    half_angle = float(np.arccos(np.clip(np.abs(u @ ax2), -1.0, 1.0)).mean())
                    best = ("cone", res, {"axis": ax2, "apex": apex, "refu": ref / refn,
                                          "half_angle": half_angle})
    kind, res, extra = best
    if res >= C._RECOGNIZE_TOL_EXACT:
        return None
    out = {"kind": kind, "residual": res, "P": P, "N": N}
    out.update(extra)
    out["gap"] = C._analytic_surface_gap(kind, out, P)
    out["gap_relative"] = out["gap"] / scale
    return out


def recognize(adaptor, uv_bounds):
    """The recogniser the census runs against -- the shipping one unless `--legacy-recognizer`."""
    if _LEGACY_RECOGNIZER:
        return _recognize_analytic_surface_prefix(adaptor, uv_bounds)
    return C._recognize_analytic_surface(adaptor, uv_bounds)


_LEGACY_RECOGNIZER = False


def stored_model(face):
    """The face's analytic surface, taken from what the STEP file actually stores. None for a
    bspline/bezier/revolution/extrusion/offset face."""
    adaptor = BRepAdaptor_Surface(face)
    kind = C._SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")
    try:
        if kind == "plane":
            axis = adaptor.Plane().Axis()
            return {"kind": "plane", "source": "stored",
                    "point": _vec(axis.Location()), "normal": _unit(_vec(axis.Direction()))}
        if kind == "cylinder":
            cyl = adaptor.Cylinder()
            axis = cyl.Axis()
            return {"kind": "cylinder", "source": "stored",
                    "origin": _vec(axis.Location()), "axis": _unit(_vec(axis.Direction())),
                    "radius": float(cyl.Radius())}
        if kind == "cone":
            cone = adaptor.Cone()
            semi = float(cone.SemiAngle())
            axis = _unit(_vec(cone.Axis().Direction()))
            if semi < 0.0:
                axis, semi = -axis, -semi
            return {"kind": "cone", "source": "stored",
                    "apex": _vec(cone.Apex()), "axis": axis, "half_angle": semi}
        if kind == "sphere":
            sph = adaptor.Sphere()
            return {"kind": "sphere", "source": "stored",
                    "centre": _vec(sph.Location()), "radius": float(sph.Radius())}
        if kind == "torus":
            tor = adaptor.Torus()
            axis = tor.Axis()
            return {"kind": "torus", "source": "stored",
                    "centre": _vec(axis.Location()), "axis": _unit(_vec(axis.Direction())),
                    "major": float(tor.MajorRadius()), "minor": float(tor.MinorRadius())}
    except Exception:
        return None
    return None


def recognized_model(rec):
    """The same dict shape, from `_recognize_analytic_surface`'s output."""
    kind = rec["kind"]
    if kind == "plane":
        return {"kind": "plane", "source": "recognized",
                "point": np.asarray(rec["point"], float), "normal": _unit(rec["normal"])}
    if kind == "cylinder":
        return {"kind": "cylinder", "source": "recognized",
                "origin": np.asarray(rec["origin"], float), "axis": _unit(rec["axis"]),
                "radius": float(rec["radius"])}
    if kind == "cone":
        return {"kind": "cone", "source": "recognized",
                "apex": np.asarray(rec["apex"], float), "axis": _unit(rec["axis"]),
                "half_angle": float(rec["half_angle"])}
    if kind == "sphere":
        return {"kind": "sphere", "source": "recognized",
                "centre": np.asarray(rec["centre"], float), "radius": float(rec["radius"])}
    return None


def surface_distance(model, P):
    """|distance| from every row of P (native CAD units) to the model's implicit surface."""
    P = np.atleast_2d(np.asarray(P, dtype=float))
    kind = model["kind"]
    if kind == "plane":
        return np.abs((P - model["point"]) @ model["normal"])
    if kind == "sphere":
        return np.abs(np.linalg.norm(P - model["centre"], axis=1) - model["radius"])
    if kind == "cylinder":
        axis = model["axis"]
        rel = P - model["origin"]
        perp = rel - np.outer(rel @ axis, axis)
        return np.abs(np.linalg.norm(perp, axis=1) - model["radius"])
    if kind == "cone":
        axis = model["axis"]
        rel = P - model["apex"]
        h = rel @ axis
        r = np.linalg.norm(rel - np.outer(h, axis), axis=1)
        a = model["half_angle"]
        return np.abs(r * math.cos(a) - h * math.sin(a))
    if kind == "torus":
        axis = model["axis"]
        rel = P - model["centre"]
        z = rel @ axis
        rho = np.linalg.norm(rel - np.outer(z, axis), axis=1)
        return np.abs(np.sqrt((rho - model["major"]) ** 2 + z * z) - model["minor"])
    return np.full(len(P), np.inf)


def face_diagonal(face):
    box = Bnd_Box()
    try:
        brepbndlib.Add(face, box)
    except Exception:
        return float("nan")
    if box.IsVoid():
        return float("nan")
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return float(math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2))


# ---------------------------------------------------------------------------------------------
# the shipped iso test, re-implemented so it does not stop at the first bad edge
# ---------------------------------------------------------------------------------------------

_N_SAMPLES = 9  # exactly `_recognized_quadric_wire_block`'s


def wire_block_census(face, project):
    """Mirror of `_recognized_quadric_wire_block`, returning a verdict for **every** edge.

    Returns (status, edges, window) where `status` is None on "the shipped function would have
    accepted this face's wires" and otherwise the shipped rejection reason for the *first*
    offending edge (so the two verdicts can be compared), and `edges` is one dict per boundary
    edge in `_face_wire_edges` order.
    """
    wires_edges = list(C._face_wire_edges(face))
    if not wires_edges:
        return "recognized quadric face has no wires", [], None

    per_wire = []
    all_other = []
    for _wire, is_outer, edges in wires_edges:
        if len(edges) < 3:
            return "recognized quadric trim wire has fewer than 3 edges", [], None
        edge_samples, degenerate_other, edge_objs = [], [], []
        for edge, start_vertex in edges:
            edge_objs.append(edge)
            if BRep_Tool.Degenerated(edge):
                _phi, other = project(BRep_Tool.Pnt(start_vertex))
                edge_samples.append(None)
                degenerate_other.append(other)
                all_other.append(other)
                continue
            degenerate_other.append(None)
            try:
                curve3d, first, last = BRep_Tool.Curve(edge)
            except Exception:
                curve3d = None
            if curve3d is None:
                return "recognized quadric boundary edge has no 3D curve", [], None
            reversed_edge = edge.Orientation() == TopAbs_REVERSED
            samples = []
            for k in range(_N_SAMPLES):
                tau = k / (_N_SAMPLES - 1.0)
                t = (1.0 - tau) if reversed_edge else tau
                phi, other = project(curve3d.Value(first + t * (last - first)))
                samples.append((phi, other))
                all_other.append(other)
            edge_samples.append(samples)
        per_wire.append((is_outer, edge_samples, degenerate_other, edge_objs))

    tol_other = 1e-6 * max(1.0, max(all_other) - min(all_other))
    tol_phi = 1e-7

    out_edges = []
    first_reject = None
    outer_window = None
    n_outer = 0
    for is_outer, edge_samples, degenerate_other, edge_objs in per_wire:
        prev_phi = None
        unwrapped_edges = []
        for samples in edge_samples:
            if samples is None:
                if prev_phi is None:
                    prev_phi = 0.0
                unwrapped_edges.append([prev_phi] * _N_SAMPLES)
                continue
            u_edge = []
            for phi_raw, _other in samples:
                if prev_phi is None:
                    phi_u = phi_raw
                else:
                    d = phi_raw - prev_phi
                    d -= 2.0 * math.pi * math.floor((d + math.pi) / (2.0 * math.pi))
                    phi_u = prev_phi + d
                u_edge.append(phi_u)
                prev_phi = phi_u
            unwrapped_edges.append(u_edge)

        all_phi_u, all_other_w = [], []
        for i, samples in enumerate(edge_samples):
            phis_u = unwrapped_edges[i]
            others = ([degenerate_other[i]] * _N_SAMPLES if samples is None
                      else [o for _p, o in samples])
            all_phi_u.extend(phis_u)
            all_other_w.extend(others)
            span_other = max(others) - min(others)
            span_phi = max(phis_u) - min(phis_u)
            iso_other = span_other <= tol_other
            iso_phi = span_phi <= tol_phi
            accepted = bool(iso_other or iso_phi)
            if not accepted and first_reject is None:
                first_reject = ("recognized quadric boundary edge is not axis-aligned in "
                                "(phi, h/theta)")
            out_edges.append({
                "edge": edge_objs[i],
                "is_outer": bool(is_outer),
                "degenerate": samples is None,
                "accepted": accepted,
                "iso_other": bool(iso_other),
                "iso_phi": bool(iso_phi),
                "span_other": float(span_other),
                "span_phi": float(span_phi),
                "tol_other": float(tol_other),
                "tol_phi": float(tol_phi),
                "uv": None if samples is None else list(zip(phis_u, others)),
            })
        if is_outer:
            n_outer += 1
            outer_window = (min(all_phi_u), max(all_phi_u) - min(all_phi_u),
                            min(all_other_w), max(all_other_w))

    if first_reject is not None:
        return first_reject, out_edges, outer_window
    if n_outer != 1:
        return (f"recognized quadric face has {n_outer} outer trim wires (expected exactly 1)",
                out_edges, outer_window)
    return None, out_edges, outer_window


def projector_for(rec, scale_to_cm):
    """The exact `project` closure `recognize_and_extract_face` builds for each recognised kind."""
    kind = rec["kind"]
    s = scale_to_cm
    if kind == "cylinder":
        axis = _unit(rec["axis"])
        refu = _unit(rec["refu"] - np.dot(rec["refu"], axis) * axis)
        e2 = np.cross(axis, refu)
        origin = np.asarray(rec["origin"], float)

        def project(pnt):
            rel = np.array([pnt.X(), pnt.Y(), pnt.Z()]) - origin
            return math.atan2(np.dot(rel, e2), np.dot(rel, refu)), float(np.dot(rel, axis)) * s
        return project, "h"
    if kind == "cone":
        axis = _unit(rec["axis"])
        refu = _unit(rec["refu"] - np.dot(rec["refu"], axis) * axis)
        e2 = np.cross(axis, refu)
        apex = np.asarray(rec["apex"], float)

        def project(pnt):
            rel = np.array([pnt.X(), pnt.Y(), pnt.Z()]) - apex
            return math.atan2(np.dot(rel, e2), np.dot(rel, refu)), float(np.dot(rel, axis)) * s
        return project, "h"
    if kind == "sphere":
        centre = np.asarray(rec["centre"], float)
        polar = np.array([0.0, 0.0, 1.0])
        refu = C._arbitrary_orthonormal_frame(polar)
        e2 = np.cross(polar, refu)
        radius = float(rec["radius"])

        def project(pnt):
            rel = (np.array([pnt.X(), pnt.Y(), pnt.Z()]) - centre) / radius
            theta = math.acos(max(-1.0, min(1.0, float(np.dot(rel, polar)))))
            return math.atan2(float(np.dot(rel, e2)), float(np.dot(rel, refu))), theta
        return project, "theta"
    return None, None


# ---------------------------------------------------------------------------------------------
# per-edge classification
# ---------------------------------------------------------------------------------------------

_DENSE = 129  # dense samples per edge for the deviation measurement

# How many declared BRep tolerances an edge may sit off its two implicit surfaces and still count
# as bucket B. This is NOT a free constant: it is calibrated by the probe's own positive control,
# the edges the shipped test *accepts* (rims and generators, which are exactly the intersection of
# the recognised quadric with a neighbouring cap plane). Those reach 1.39 declared tolerances on
# ALICE3, so a threshold below that would reject known-exact intersections. `--b-factor` moves it
# and the summary prints the whole sweep, so the choice is inspectable rather than asserted.
_B_FACTOR = 2.0


def dense_edge_samples(edge):
    try:
        curve3d, first, last = BRep_Tool.Curve(edge)
    except Exception:
        return None
    if curve3d is None:
        return None
    ts = np.linspace(first, last, _DENSE)
    return np.array([[p.X(), p.Y(), p.Z()] for p in (curve3d.Value(float(t)) for t in ts)])


def polyline_length(P):
    return float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())


def planarity(P):
    """Max distance of the samples from their own best-fit plane, and that plane's normal."""
    centre = P.mean(axis=0)
    _u, _s, vt = np.linalg.svd(P - centre, full_matrices=False)
    normal = vt[-1]
    return float(np.abs((P - centre) @ normal).max()), normal


def classify_edge(edge, entry, face_rec_model, face, neighbours, scale_to_cm,
                  other_name, radius_hint, sphere_alt_axis_ok=False):
    """Bucket one boundary edge. Returns a JSON-able dict (lengths in cm)."""
    s = scale_to_cm
    tol_edge = float(BRep_Tool.Tolerance(edge))          # native units
    P = dense_edge_samples(edge)
    info = {
        "accepted": entry["accepted"],
        "degenerate": entry["degenerate"],
        "spanOther": entry["span_other"],
        "spanPhi": entry["span_phi"],
        "tolOtherUsed": entry["tol_other"],
        "tolPhiUsed": entry["tol_phi"],
        "tolEdgeCm": tol_edge * s,
        "selfKind": face_rec_model["kind"],
    }
    if P is None:
        info["bucket"] = "D"
        info["bucketDetail"] = "edge has no 3D curve"
        return info

    edge_len = polyline_length(P) * s
    diag_self = face_diagonal(face) * s
    info["edgeLengthCm"] = edge_len
    info["patchDiagonalCm"] = diag_self

    # ---- physical size of the parametric spans, so a "near-iso" claim is dimensionally honest.
    # `other` is already cm for a cylinder/cone and radians for a sphere; `phi` is always radians.
    if other_name == "theta":
        span_other_cm = entry["span_other"] * radius_hint * s
    else:
        span_other_cm = entry["span_other"]
    span_phi_cm = entry["span_phi"] * radius_hint * s
    info["spanOtherCm"] = span_other_cm
    info["spanPhiCm"] = span_phi_cm

    # ---- A1: iso in the shipped parametrisation at the level of the edge's own tolerance
    tol_cm = tol_edge * s
    a1 = (span_other_cm <= tol_cm) or (span_phi_cm <= tol_cm)

    # ---- A2: iso only in an alternative but equally valid parametrisation of the same surface.
    # Sphere: the polar axis is free (the converter hard-codes z). A planar edge on a sphere is a
    # circle of constant theta about that plane's normal -- so planarity is the *per-edge* test.
    # It is not enough on its own: `_recognized_quadric_wire_block` needs ONE frame for the whole
    # face, so A2 is only claimed when `sphere_alternative_axis` finds a single axis that serves
    # every edge of the face. `edgeIsPlanar` is reported either way, because "each edge is planar
    # but no common axis exists" is the thing that has to be visible for the claim to be checkable.
    plane_res, _n = planarity(P)
    info["planarityCm"] = plane_res * s
    info["edgeIsPlanar"] = bool(plane_res * s <= tol_cm)
    a2 = (face_rec_model["kind"] == "sphere") and info["edgeIsPlanar"] and sphere_alt_axis_ok

    # ---- neighbours
    others = [g for g in neighbours if g is not None]
    info["nNeighbourFaces"] = len(others)
    dev_self = float(surface_distance(face_rec_model, P).max()) * s
    info["devSelfCm"] = dev_self

    if not entry["accepted"] and a1:
        info["bucket"] = "A"
        info["bucketDetail"] = "A1 iso within the edge's own declared tolerance"
    elif not entry["accepted"] and a2:
        info["bucket"] = "A"
        info["bucketDetail"] = "A2 iso for an alternative polar axis (sphere)"
    else:
        info["bucket"] = None  # decided below

    if len(others) != 1:
        detail = ("seam edge: the face is its own neighbour" if len(others) == 0
                  else f"edge shared by {len(others) + 1} faces")
        if info["bucket"] is None:
            info["bucket"] = "D"
            info["bucketDetail"] = detail
        info["neighbourKind"] = None
        return info

    gface, gmodel = others[0]
    info["neighbourKind"] = None if gmodel is None else gmodel["kind"]
    info["neighbourSource"] = None if gmodel is None else gmodel["source"]
    if gmodel is None:
        if info["bucket"] is None:
            info["bucket"] = "C"
            info["bucketDetail"] = "adjacent face is free-form"
        return info

    dev_nbr = float(surface_distance(gmodel, P).max()) * s
    diag_nbr = face_diagonal(gface) * s
    info["devNeighbourCm"] = dev_nbr
    info["neighbourDiagonalCm"] = diag_nbr
    dev = max(dev_self, dev_nbr)
    # The threshold: the *declared* tolerance of the three BRep entities involved. OCCT's own
    # contract is that an edge's 3D curve, its pcurves and its vertices agree to within the edge
    # tolerance, so "the edge lies on both surfaces to within what the exporter itself claims"
    # is exactly `dev <= tol`. The face tolerances are folded in because a face may declare a
    # looser one than its edges; the sweep in the summary shows how little that choice matters.
    tol_face_self = float(BRep_Tool.Tolerance(face)) * s
    tol_face_nbr = float(BRep_Tool.Tolerance(gface)) * s
    tol_ref = max(tol_cm, tol_face_self, tol_face_nbr)
    info["tolFaceSelfCm"] = tol_face_self
    info["tolFaceNeighbourCm"] = tol_face_nbr
    info["tolRefCm"] = tol_ref
    info["devCm"] = dev
    info["devOverEdgeLength"] = dev / edge_len if edge_len > 0 else float("inf")
    info["devOverPatchDiagonal"] = dev / diag_self if diag_self > 0 else float("inf")
    info["devOverEdgeTolerance"] = dev / tol_cm if tol_cm > 0 else float("inf")
    info["devOverTolRef"] = dev / tol_ref if tol_ref > 0 else float("inf")
    info["pair"] = "|".join(sorted([face_rec_model["kind"], gmodel["kind"]]))
    if info["bucket"] is None:
        within = dev <= _B_FACTOR * tol_ref
        info["bucket"] = "B" if within else "D"
        info["bucketDetail"] = (f"edge lies on both implicit surfaces within {_B_FACTOR:g}x the "
                                "declared BRep tolerance" if within else
                                "both surfaces analytic but the edge is not on both")
    return info


# ---------------------------------------------------------------------------------------------
# per-solid drive
# ---------------------------------------------------------------------------------------------

def _ancestor_faces(anc, edge):
    """The faces of the solid that reference this edge, as TopoDS_Face."""
    try:
        lst = anc.FindFromKey(edge)
    except Exception:
        return []
    out = []
    it = TopTools_ListIteratorOfListOfShape(lst)
    while it.More():
        out.append(topods.Face(it.Value()))
        it.Next()
    return out


def census_solid(lid, name, shape, scale_to_cm, verbose=False):
    faces = list(TopologyExplorer(shape).faces())
    face_map = TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape, TopAbs_FACE, face_map)
    anc = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, anc)

    model_cache = {}

    def model_of(face):
        idx = face_map.FindIndex(face)
        if idx in model_cache:
            return model_cache[idx]
        m = stored_model(face)
        if m is None:
            adaptor = BRepAdaptor_Surface(face)
            try:
                uv = breptools.UVBounds(face)
                rec = recognize(adaptor, uv)
            except Exception:
                rec = None
            m = recognized_model(rec) if rec is not None else None
        model_cache[idx] = m
        return m

    solid = {
        "lid": lid, "name": name, "nFaces": len(faces),
        "emitted": False, "faceReasons": [], "faces": [], "edges": [], "planarEdges": [],
        "verdictMismatches": 0,
    }

    n_unsupported = 0
    n_trim_declined = 0
    for i, face in enumerate(faces):
        adaptor = BRepAdaptor_Surface(face)
        surf_type = C._SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")
        extractor = C._FACE_EXTRACTORS.get(surf_type)
        if extractor is None:
            record, reason = None, f"{surf_type} face extraction not implemented yet"
        else:
            try:
                record, reason = extractor(face, scale_to_cm)
            except Exception as exc:
                record, reason = None, f"{surf_type} extractor raised {exc!r}"
        recognised_reason = None
        if record is None:
            try:
                rec_record, rec_reason = C.recognize_and_extract_face(face, scale_to_cm)
            except Exception as exc:
                rec_record, rec_reason = None, f"recognition raised {exc!r}"
            if rec_record is not None:
                record, reason = rec_record, None
            elif rec_reason is not None:
                recognised_reason = rec_reason
                reason = f"{reason}; recognition attempted: {rec_reason}"
        if record is None:
            solid["faceReasons"].append(reason)
            if recognised_reason is not None:
                n_trim_declined += 1
            else:
                n_unsupported += 1

        # ---- second census, different code path, same question. `extract_planar_face` declines a
        # planar face whose boundary carries a curve type outside {line, circle, bspline, bezier}.
        # On ALICE3 that never fires, but it is exactly what stops `oblique_cut_cyl` (the ladder
        # fixture that has never emitted) and Bagger's `Bucket`: an oblique cut through a cylinder
        # bounds its cap plane with an ELLIPSE -- the brief's own motivating example, and an exact
        # co-surface intersection. Censused separately so the two mechanisms are never added up.
        if record is None and reason and reason.startswith("planar boundary edge is a"):
            pmodel = stored_model(face)
            if pmodel is not None:
                for _wire, _outer, wedges in C._face_wire_edges(face):
                    for edge, _v in wedges:
                        if BRep_Tool.Degenerated(edge):
                            continue
                        try:
                            ctype = C._CURVE_TYPE_NAMES.get(BRepAdaptor_Curve(edge).GetType(),
                                                            "unknown")
                        except Exception:
                            ctype = "unknown"
                        if ctype in ("line", "circle", "bspline", "bezier"):
                            continue
                        entry = {"edge": edge, "accepted": False, "degenerate": False,
                                 "span_other": float("inf"), "span_phi": float("inf"),
                                 "tol_other": float("nan"), "tol_phi": float("nan"),
                                 "iso_other": False, "iso_phi": False, "is_outer": True,
                                 "uv": None}
                        nbrs = [(g, model_of(g)) for g in _ancestor_faces(anc, edge)
                                if not g.IsSame(face)]
                        info = classify_edge(edge, entry, pmodel, face, nbrs, scale_to_cm,
                                             "planar", 1.0)
                        info["face"] = i
                        info["solid"] = name
                        info["curveType"] = ctype
                        solid["planarEdges"].append(info)

        # ---- the census proper: only faces that reach the wire block
        if extractor is not None and record is not None and record.get("recognized") is None:
            continue   # a natively-analytic face that extracted on its own
        try:
            uv = breptools.UVBounds(face)
            rec = recognize(adaptor, uv)
        except Exception:
            rec = None
        if rec is None or rec["kind"] == "plane":
            continue
        model_cache[face_map.FindIndex(face)] = recognized_model(rec)
        project, other_name = projector_for(rec, scale_to_cm)
        if project is None:
            continue
        status, entries, _win = wire_block_census(face, project)

        # self-check 1: my verdict must be the shipped one
        shipped_ok = record is not None and record.get("recognized") is not None
        mine_ok = status is None
        if shipped_ok != mine_ok:
            # the shipped path can still decline an accepted wire block on the phi-sweep guard
            if not (mine_ok and not shipped_ok and recognised_reason and
                    "wraps more than a full turn" in recognised_reason):
                solid["verdictMismatches"] += 1
                if verbose:
                    print(f"    [mismatch] {name} f#{i}: shipped={shipped_ok} mine={mine_ok} "
                          f"status={status!r} reason={recognised_reason!r}")

        fmodel = recognized_model(rec)
        radius_hint = (float(rec["radius"]) if rec["kind"] in ("cylinder", "sphere")
                       else _cone_radius_hint(rec))
        n_rej = sum(1 for e in entries if not e["accepted"])
        alt_axis = None
        if rec["kind"] == "sphere" and n_rej:
            ok, axis, n_cand = sphere_alternative_axis(face, rec, scale_to_cm)
            alt_axis = {"ok": bool(ok), "axis": axis, "candidates": n_cand}
        solid["faces"].append({
            "index": i, "recognizedKind": rec["kind"], "storedType": surf_type,
            "nEdges": len(entries), "nRejected": n_rej,
            "emitted": bool(shipped_ok), "wireBlockStatus": status,
            "sphereAlternativeAxis": alt_axis,
        })
        for entry in entries:
            edge = entry["edge"]
            nbrs = []
            for g in _ancestor_faces(anc, edge):
                if g.IsSame(face):
                    continue
                nbrs.append((g, model_of(g)))
            info = classify_edge(edge, entry, fmodel, face, nbrs, scale_to_cm,
                                 other_name, radius_hint,
                                 sphere_alt_axis_ok=bool(alt_axis and alt_axis["ok"]))
            info["face"] = i
            info["solid"] = name
            solid["edges"].append(info)

    solid["nUnsupportedFaces"] = n_unsupported
    solid["nTrimDeclinedFaces"] = n_trim_declined
    # `extract_surfaces_for_shape` treats a shape with no faces as unsupported ("shape has no
    # faces"); IRIS carries exactly one such leaf (`COMPOUND`) and forgetting this makes the
    # emitted count 12 where the converter says 11.
    solid["emitted"] = bool(faces) and not solid["faceReasons"]
    return solid


def sphere_alternative_axis(face, rec, scale_to_cm):
    """Face-level test for bucket A2: is there ONE polar axis that makes *every* boundary edge of
    this recognised sphere iso?

    Per-edge A2 is necessary but not sufficient. `_recognized_quadric_wire_block` needs a single
    (phi, theta) frame for the whole face, so a patch bounded by circles in non-parallel planes is
    not rescued by any axis choice however planar each individual edge is. This is the honest
    version of the question and it is what the solid-level rollup uses.

    A circle of constant theta about `n` is the sphere cut by a plane perpendicular to `n`; a
    meridian (constant phi) is the sphere cut by a plane *containing* the centre and `n`. Both are
    tested directly in 3D against the edge's own declared tolerance, so no angular constant is
    invented.

    Returns (ok, axis_or_None, n_candidates_tried).
    """
    centre = np.asarray(rec["centre"], float)
    edges = []
    for _wire, _is_outer, wire_edges in C._face_wire_edges(face):
        for edge, _v in wire_edges:
            if BRep_Tool.Degenerated(edge):
                continue
            P = dense_edge_samples(edge)
            if P is None:
                return False, None, 0
            tol = float(BRep_Tool.Tolerance(edge))
            res, normal = planarity(P)
            edges.append((P, tol, res, normal))
    if not edges:
        return False, None, 0

    candidates = [np.array([0.0, 0.0, 1.0])]          # the shipped choice, as a control
    for _P, tol, res, normal in edges:
        if res <= tol:
            candidates.append(normal)
    # a lune (meridian edges only): the axis is the intersection of two edge planes
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            cross = np.cross(edges[i][3], edges[j][3])
            n = float(np.linalg.norm(cross))
            if n > 1e-6:
                candidates.append(cross / n)

    for cand in candidates:
        axis = _unit(cand)
        ok = True
        for P, tol, _res, _normal in edges:
            rel = P - centre
            h = rel @ axis
            iso_theta = (h.max() - h.min()) <= 2.0 * tol
            m = np.cross(axis, rel[0])
            mn = float(np.linalg.norm(m))
            iso_phi = False
            if mn > 1e-12:
                m = m / mn
                iso_phi = float(np.abs(rel @ m).max()) <= 2.0 * tol
            if not (iso_theta or iso_phi):
                ok = False
                break
        if ok:
            return True, axis.tolist(), len(candidates)
    return False, None, len(candidates)


def _cone_radius_hint(rec):
    """A representative radius for a recognised cone: the mean distance of its recognition
    samples from the axis. Used only to turn a phi span in radians into a length."""
    P = np.asarray(rec["P"], float)
    axis = _unit(rec["axis"])
    rel = P - np.asarray(rec["apex"], float)
    perp = rel - np.outer(rel @ axis, axis)
    return float(np.linalg.norm(perp, axis=1).mean())


# ---------------------------------------------------------------------------------------------

def run_model(step_path, label, verbose=False):
    unit = C.detect_step_length_unit(step_path)
    scale = C.step_unit_scale_to_cm(unit)
    C.extract_graph(step_path, meshparam=None, scale_to_cm=scale)
    solids = []
    for lid, shape in C.def_shapes.items():
        name = C.def_names.get(lid, lid)
        if verbose:
            print(f"  [{label}] {name}")
        solids.append(census_solid(lid, name, shape, scale, verbose=verbose))
    return {"model": step_path, "label": label, "unit": unit, "scaleToCm": scale,
            "solids": solids}


def _q(values, qs=(0.5, 0.9, 0.99, 1.0)):
    if not values:
        return {}
    a = np.sort(np.asarray(values, float))
    return {f"p{int(q * 100)}": float(a[min(len(a) - 1, int(math.ceil(q * len(a)) - 1))])
            for q in qs}


def summarise(result):
    solids = result["solids"]
    buckets = Counter()
    pairs = Counter()
    dev_rel_accepted, dev_rel_rejected = [], []
    dev_rel_analytic_nbr = []          # every rejected edge with an analytic neighbour
    n_edges_on_recognized_faces = 0
    n_accepted = 0
    for s in solids:
        for e in s["edges"]:
            n_edges_on_recognized_faces += 1
            if e["accepted"]:
                n_accepted += 1
                if "devOverTolRef" in e:
                    dev_rel_accepted.append(e["devOverTolRef"])
                continue
            buckets[e["bucket"]] += 1
            if "devOverTolRef" in e:
                dev_rel_analytic_nbr.append(e["devOverTolRef"])
            if e["bucket"] == "B":
                pairs[e["pair"]] += 1
                dev_rel_rejected.append(e["devOverTolRef"])
    sweep = {f"le{int(k)}x": sum(1 for v in dev_rel_analytic_nbr if v <= k)
             for k in (1, 2, 5, 10, 100, 1000)}

    # ---- the number that decides anything: the solid-level rollup.
    # Population: leaf solids that emit no sidecar and whose every declined face was declined for a
    # *trim* reason -- the recognised-quadric iso test, or the planar extractor's curve vocabulary.
    # A solid with even one genuinely unsupported surface is out of scope here and is listed
    # separately, because no trim work can rescue it.
    rollup = {"fullyCovered": [], "stillFails": [], "notTrimBlocked": []}
    for s in solids:
        if s["emitted"] or not s["nFaces"]:
            continue
        if s["nUnsupportedFaces"] or not s["nTrimDeclinedFaces"]:
            rollup["notTrimBlocked"].append(s["name"])
            continue
        rej = [e for e in s["edges"] if not e["accepted"]] + s["planarEdges"]
        bad = [e for e in rej if e["bucket"] not in ("A", "B")]
        other_face_reasons = sorted({
            f["wireBlockStatus"] for f in s["faces"]
            if f["wireBlockStatus"] and "not axis-aligned" not in f["wireBlockStatus"]})
        entry = {"solid": s["name"], "nRejected": len(rej),
                 "nPlanarTrimEdges": len(s["planarEdges"]),
                 "buckets": dict(Counter(e["bucket"] for e in rej)),
                 "nUncovered": len(bad),
                 "otherWireBlockFailures": other_face_reasons,
                 "sphereAxisFaces": sum(1 for f in s["faces"]
                                        if f.get("sphereAlternativeAxis")),
                 "sphereAxisFacesRescued": sum(
                     1 for f in s["faces"]
                     if f.get("sphereAlternativeAxis") and f["sphereAlternativeAxis"]["ok"])}
        if bad or other_face_reasons:
            rollup["stillFails"].append(entry)
        else:
            rollup["fullyCovered"].append(entry)
    return {
        "rollup": rollup,
        "nEligibleNotEmitted": len(rollup["fullyCovered"]) + len(rollup["stillFails"]),
        "bFactor": _B_FACTOR,
        "thresholdSweepRejectedWithAnalyticNeighbour": sweep,
        "nRejectedWithAnalyticNeighbour": len(dev_rel_analytic_nbr),
        "nSolids": len(solids),
        "nEmitted": sum(1 for s in solids if s["emitted"]),
        "nPlanarTrimEdges": sum(len(s["planarEdges"]) for s in solids),
        "planarTrimBuckets": dict(Counter(e["bucket"] for s in solids for e in s["planarEdges"])),
        "nEdgesOnRecognizedFaces": n_edges_on_recognized_faces,
        "nAcceptedEdges": n_accepted,
        "nRejectedEdges": n_edges_on_recognized_faces - n_accepted,
        "buckets": dict(buckets),
        "pairs": dict(pairs),
        "verdictMismatches": sum(s["verdictMismatches"] for s in solids),
        "devOverTolRefRejectedB": _q(dev_rel_rejected),
        "devOverTolRefAcceptedControl": _q(dev_rel_accepted),
        "devOverTolRefRejectedAnalyticNeighbour": _q(dev_rel_analytic_nbr),
    }


def main():
    global _B_FACTOR, _LEGACY_RECOGNIZER
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", action="append", default=[],
                    help="STEP file to census; repeatable. Use label:path to name it.")
    ap.add_argument("--fixtures", action="store_true",
                    help="also build and census the ladder fixtures")
    ap.add_argument("--json", help="write the full per-edge census here")
    ap.add_argument("--b-factor", type=float, default=_B_FACTOR,
                    help="how many declared BRep tolerances an edge may sit off its two implicit "
                         "surfaces and still count as bucket B (default %(default)s, calibrated "
                         "by the accepted-edge control -- see the module docstring)")
    ap.add_argument("--legacy-recognizer", action="store_true",
                    help="run the census against the PRE-FIX recogniser (git 237be7f81a^), which "
                         "is what Stream_K_Tier0.md §2's 1891/1053 figures were measured with")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    _B_FACTOR = args.b_factor
    _LEGACY_RECOGNIZER = args.legacy_recognizer

    jobs = []
    if args.fixtures:
        outdir = Path(tempfile.mkdtemp(prefix="trimcensus-fixtures-"))
        subprocess.run([sys.executable, str(_HERE.parent / "make_boolean_fixtures.py"),
                        "--outdir", str(outdir)], check=True)
        for f in sorted(outdir.glob("*.step")) + sorted(outdir.glob("*.stp")):
            jobs.append((f"fixture:{f.stem}", str(f)))
    for m in args.model:
        if ":" in m and not os.path.exists(m):
            label, path = m.split(":", 1)
        else:
            label, path = Path(m).stem, m
        jobs.append((label, path))

    results = []
    for label, path in jobs:
        print(f"== {label}: {path}")
        res = run_model(path, label, verbose=args.verbose)
        res["summary"] = summarise(res)
        results.append(res)
        s = res["summary"]
        print(f"   solids {s['nSolids']}, emitted {s['nEmitted']}, "
              f"edges on recognised faces {s['nEdgesOnRecognizedFaces']} "
              f"(accepted {s['nAcceptedEdges']}, rejected {s['nRejectedEdges']})")
        print(f"   buckets {s['buckets']}   verdict mismatches {s['verdictMismatches']}")
        print(f"   planar-trim edges (2nd mechanism) {s['nPlanarTrimEdges']} "
              f"{s['planarTrimBuckets']}")
        print(f"   pairs   {s['pairs']}")
        print(f"   dev/tol rejected-B          {s['devOverTolRefRejectedB']}")
        print(f"   dev/tol rejected w/ nbr     {s['devOverTolRefRejectedAnalyticNeighbour']}")
        print(f"   dev/tol accepted control    {s['devOverTolRefAcceptedControl']}")
        print(f"   threshold sweep             {s['thresholdSweepRejectedWithAnalyticNeighbour']}"
              f" of {s['nRejectedWithAnalyticNeighbour']}")
        ru = s["rollup"]
        print(f"   ROLLUP: eligible-but-not-emitted {s['nEligibleNotEmitted']}; "
              f"fully covered by A+B {len(ru['fullyCovered'])}; "
              f"still fails {len(ru['stillFails'])}")
        for e in ru["stillFails"]:
            print(f"     still fails: {e['solid']:16s} uncovered {e['nUncovered']:4d} of "
                  f"{e['nRejected']:4d}  {e['buckets']}  {e['otherWireBlockFailures']}")
        for e in ru["fullyCovered"]:
            print(f"     covered    : {e['solid']:16s} {e['nRejected']:4d} rejected "
                  f"{e['buckets']}  sphere-axis faces rescued "
                  f"{e['sphereAxisFacesRescued']}/{e['sphereAxisFaces']}")

    if args.json:
        out = Path(args.json)
        out.parent.mkdir(parents=True, exist_ok=True)
        def clean(o):
            if isinstance(o, dict):
                return {k: clean(v) for k, v in o.items() if k not in ("edge",)}
            if isinstance(o, list):
                return [clean(v) for v in o]
            if isinstance(o, (np.floating, np.integer)):
                return o.item()
            return o
        out.write_text(json.dumps(clean(results), indent=1))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
