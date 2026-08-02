#!/usr/bin/env python3.10
"""Stream R probe: is a face's trimmed region the conjunction of its neighbours' half-spaces?

The question
------------
`Stream_O_ImplicitTrims.md` measured that 691 of ALICE3's 763 rejected boundary edges are *exactly*
the intersection of the face's own analytic surface with a neighbouring face's analytic surface.
That licenses a trim representation with no parametric curve at all: a point on the face's surface
is inside the trimmed patch **iff it is on the correct side of every neighbouring surface**.

That rule is not automatically the same set as the true trimmed face. It is right exactly when the
trimmed region is a single cell of the arrangement of the neighbouring surfaces, and that cell is
the face. This probe measures, per face, on the real corpora, how often that holds -- and when it
does not, by how much and where.

The instrument
--------------
For every face with an analytic surface (stored *or* recognised by the shipping recogniser):

1. **Trimming set.** Walk the face's wires; for each boundary edge take the adjacent face through
   `TopExp::MapShapesAndAncestors(EDGE, FACE)` -- identity, never proximity -- and its analytic
   surface (stored, else recognised). Deduplicate *geometrically*: two neighbour faces lying on the
   same plane contribute one half-space. K = the number of distinct trimming surfaces.

2. **Sample the face's SURFACE, not its boundary.** A grid over the face's own (u, v) rectangle,
   evaluated on the real `Geom_Surface`. The rectangle is the bounding box of the trim wires in
   (u, v), so it contains the trimmed region *and* the parts of the surface just outside it -- which
   is exactly where a conjunction over-accepts.

3. **Ground truth: `BRepTopAdaptor_FClass2d`.** It answers precisely the question asked -- "is this
   (u, v) inside this face's trim wires" -- with no projection, no proximity band and no reliance on
   the solid being closed. `BRepClass3d_SolidClassifier` answers a *different* question (is this
   point inside the solid), is undefined ON the boundary, which is where every one of these samples
   sits, and would conflate a mis-trimmed face with a neighbouring face's coverage of the same
   region. It is used here only as an independent cross-check on a subsample (`--solid-crosscheck`).

4. **The sign vector.** For each sample compute sigma(P) = (sign f_1(P), ..., sign f_K(P)), the cell
   of the arrangement the point sits in. Then:

     cellsIn  = the number of *distinct* sign vectors over the truth-IN samples.
                cellsIn == 1  <=>  the trimmed region is a single cell  <=>  a conjunction of
                one sense per surface can express it. cellsIn == n means the smallest exact
                implicit description is a disjunction of n conjunctions (a DNF with n terms).
     leak     = truth-OUT samples whose sign vector is one of the IN cells. These are points the
                implicit description accepts and the CAD face does not, and **no** amount of DNF
                structure over the same K surfaces removes them: the arrangement cell genuinely
                straddles the boundary. This is the number that decides whether the idea works.

   `cellsIn` and `leak` are independent of any sense convention, so neither can be tuned.

5. **The naive conjunction, scored directly.** sense_k is read off the deepest interior sample
   (the truth-IN sample furthest from the boundary polyline) -- one bit per surface, exactly what
   the converter would store. The rule is `all(sense_k * f_k(P) >= 0)`. Reported as false
   positives (rule IN, truth OUT) and false negatives (rule OUT, truth IN), each with the 3D
   distance from the sample to the face's own boundary wire, because a disagreement 1e-12 cm from
   an edge is a tolerance artefact and one 3 cm inside is a broken design.

6. **The far field.** The near-field grid cannot see a cell that re-enters the surface outside the
   face's own (u, v) rectangle -- a cylinder's phi wrap, a cone's mirror nappe. A second sampling
   runs over the face's *analytic* surface across its full natural extent (full 2*pi in phi, both
   nappes of a cone, h over twice the face's span). Only the points the rule *accepts* are then
   checked against OCCT (`ShapeAnalysis_Surface` projection + residual + `FClass2d`), so the cost
   is paid only where it matters.

Positive control
----------------
`--flip-sense` inverts one stored sense on every face and `--perturb-radius R` moves every trimming
surface by R cm. Both must move the disagreement counts. A validation that cannot fail has not
passed.

Usage
-----
    python3 probes/implicitTrimValidate.py --model ALICE_3_example/CAD_noETA.stp \
        --solids ST0923290_002,ST0923290_011 --json /tmp/r/alice3.json
    python3 probes/implicitTrimValidate.py --fixtures --json /tmp/r/fixtures.json
    python3 probes/implicitTrimValidateReport.py /tmp/r/*.json
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))

from csg.occ_env import ensure_occ  # noqa: E402

ensure_occ()

import numpy as np  # noqa: E402

from OCC.Core.BRep import BRep_Tool  # noqa: E402
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCC.Core.BRepTools import breptools  # noqa: E402
from OCC.Core.BRepTopAdaptor import BRepTopAdaptor_FClass2d  # noqa: E402
from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier  # noqa: E402
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface  # noqa: E402
from OCC.Core.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_IN,  # noqa: E402
                             TopAbs_OUT, TopAbs_ON)
from OCC.Core.TopExp import topexp  # noqa: E402
from OCC.Core.TopTools import (TopTools_IndexedMapOfShape,  # noqa: E402
                               TopTools_IndexedDataMapOfShapeListOfShape,
                               TopTools_ListIteratorOfListOfShape)
from OCC.Core.TopoDS import topods  # noqa: E402
from OCC.Core.gp import gp_Pnt, gp_Pnt2d  # noqa: E402
from OCC.Extend.TopologyUtils import TopologyExplorer  # noqa: E402

import O2_CADtoTGeo as C  # noqa: E402

# the census probe already walks faces, neighbours and recognised surfaces on every model, and its
# face verdicts were checked against the shipped `recognize_and_extract_face` with 0 mismatches.
import trimEdgeCensus as TC  # noqa: E402


# =============================================================================================
# implicit functions -- one signed, distance-like scalar per surface kind
# =============================================================================================

def implicit(model, P):
    """Signed, approximately-signed-distance value of `model` at every row of P (native units).

    Sign convention: negative inside the bounded side where one exists (sphere, cylinder, torus),
    positive on the side the normal points to for a plane. A cone uses r*cos(a) - h*sin(a), which is
    negative strictly inside the *upper* nappe and positive everywhere on the mirror nappe -- that
    asymmetry is deliberate and is one of the things being measured (`coneMirrorNappe`).
    """
    P = np.atleast_2d(np.asarray(P, dtype=float))
    kind = model["kind"]
    if kind == "plane":
        return (P - model["point"]) @ model["normal"]
    if kind == "sphere":
        return np.linalg.norm(P - model["centre"], axis=1) - model["radius"]
    if kind == "cylinder":
        axis = model["axis"]
        rel = P - model["origin"]
        perp = rel - np.outer(rel @ axis, axis)
        return np.linalg.norm(perp, axis=1) - model["radius"]
    if kind == "cone":
        axis = model["axis"]
        rel = P - model["apex"]
        h = rel @ axis
        r = np.linalg.norm(rel - np.outer(h, axis), axis=1)
        a = model["half_angle"]
        return r * math.cos(a) - h * math.sin(a)
    if kind == "torus":
        axis = model["axis"]
        rel = P - model["centre"]
        z = rel @ axis
        rho = np.linalg.norm(rel - np.outer(z, axis), axis=1)
        return np.sqrt((rho - model["major"]) ** 2 + z * z) - model["minor"]
    return np.full(len(P), np.nan)


_EVALUABLE = ("plane", "cylinder", "cone", "sphere", "torus")


def gradient(model, P):
    """Unit gradient (the surface's own normal field) at every row of P."""
    P = np.atleast_2d(np.asarray(P, dtype=float))
    kind = model["kind"]
    if kind == "plane":
        return np.tile(TC._unit(model["normal"]), (len(P), 1))
    if kind == "sphere":
        d = P - model["centre"]
        return d / np.maximum(np.linalg.norm(d, axis=1)[:, None], 1e-300)
    if kind == "cylinder":
        axis = model["axis"]
        rel = P - model["origin"]
        perp = rel - np.outer(rel @ axis, axis)
        return perp / np.maximum(np.linalg.norm(perp, axis=1)[:, None], 1e-300)
    if kind == "cone":
        axis = TC._unit(model["axis"])
        rel = P - model["apex"]
        h = rel @ axis
        perp = rel - np.outer(h, axis)
        e = perp / np.maximum(np.linalg.norm(perp, axis=1)[:, None], 1e-300)
        a = model["half_angle"]
        return math.cos(a) * e - math.sin(a) * axis
    if kind == "torus":
        axis = TC._unit(model["axis"])
        rel = P - model["centre"]
        z = rel @ axis
        perp = rel - np.outer(z, axis)
        rho = np.linalg.norm(perp, axis=1)
        e = perp / np.maximum(rho[:, None], 1e-300)
        d = (rho - model["major"])[:, None] * e + z[:, None] * axis
        return d / np.maximum(np.linalg.norm(d, axis=1)[:, None], 1e-300)
    return np.full((len(P), 3), np.nan)


def perturb(model, dr):
    """Move a surface by `dr` in native units -- the negative control."""
    m = dict(model)
    if m["kind"] == "plane":
        m["point"] = np.asarray(m["point"], float) + dr * np.asarray(m["normal"], float)
    elif m["kind"] == "sphere":
        m["radius"] = float(m["radius"]) + dr
    elif m["kind"] == "cylinder":
        m["radius"] = float(m["radius"]) + dr
    elif m["kind"] == "cone":
        m["apex"] = np.asarray(m["apex"], float) + dr * np.asarray(m["axis"], float)
    elif m["kind"] == "torus":
        m["minor"] = float(m["minor"]) + dr
    return m


# =============================================================================================
# analytic parametrisation of the face's OWN surface -- used only for the far field
# =============================================================================================

def _frame(axis):
    axis = TC._unit(axis)
    e1 = TC._unit(C._arbitrary_orthonormal_frame(axis))
    return axis, e1, np.cross(axis, e1)


def far_field_points(model, pts_on_face, n_phi, n_h):
    """A grid over the face's analytic surface across its full natural extent.

    Deliberately generous: the whole 2*pi of a cylinder/cone/sphere azimuth, both nappes of a cone,
    and twice the face's own axial span. Returns (points, label) or (None, reason).
    """
    kind = model["kind"]
    if kind == "cylinder":
        axis, e1, e2 = _frame(model["axis"])
        rel = pts_on_face - model["origin"]
        h = rel @ axis
        h0, h1 = float(h.min()), float(h.max())
        pad = 0.5 * max(h1 - h0, 1e-9)
        hs = np.linspace(h0 - pad, h1 + pad, n_h)
        ph = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
        PH, H = np.meshgrid(ph, hs, indexing="ij")
        R = float(model["radius"])
        Q = (model["origin"] + R * (np.cos(PH)[..., None] * e1 + np.sin(PH)[..., None] * e2)
             + H[..., None] * axis)
        return Q.reshape(-1, 3), "cylinder(phi 0..2pi, h +/-50%)"
    if kind == "cone":
        axis, e1, e2 = _frame(model["axis"])
        rel = pts_on_face - model["apex"]
        h = rel @ axis
        h0, h1 = float(h.min()), float(h.max())
        pad = 0.5 * max(h1 - h0, 1e-9)
        # both nappes: mirror the face's own h window through the apex
        lo, hi = h0 - pad, h1 + pad
        hs = np.concatenate([np.linspace(lo, hi, n_h), np.linspace(-hi, -lo, max(4, n_h // 2))])
        ph = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
        PH, H = np.meshgrid(ph, hs, indexing="ij")
        t = math.tan(float(model["half_angle"]))
        Q = (model["apex"] + H[..., None] * axis
             + (np.abs(H) * t)[..., None] * (np.cos(PH)[..., None] * e1
                                             + np.sin(PH)[..., None] * e2))
        return Q.reshape(-1, 3), "cone(phi 0..2pi, both nappes)"
    if kind == "sphere":
        axis, e1, e2 = _frame(np.array([0.0, 0.0, 1.0]))
        ph = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
        th = np.linspace(1e-4, math.pi - 1e-4, n_h)
        PH, TH = np.meshgrid(ph, th, indexing="ij")
        R = float(model["radius"])
        Q = (model["centre"]
             + R * (np.sin(TH)[..., None] * (np.cos(PH)[..., None] * e1 + np.sin(PH)[..., None] * e2)
                    + np.cos(TH)[..., None] * axis))
        return Q.reshape(-1, 3), "sphere(full)"
    if kind == "plane":
        n = TC._unit(model["normal"])
        _a, e1, e2 = _frame(n)
        rel = pts_on_face - model["point"]
        a, b = rel @ e1, rel @ e2
        pa = 0.6 * max(float(a.max() - a.min()), 1e-9)
        pb = 0.6 * max(float(b.max() - b.min()), 1e-9)
        A = np.linspace(float(a.min()) - pa, float(a.max()) + pa, n_phi)
        B = np.linspace(float(b.min()) - pb, float(b.max()) + pb, n_h)
        AA, BB = np.meshgrid(A, B, indexing="ij")
        Q = model["point"] + AA[..., None] * e1 + BB[..., None] * e2
        return Q.reshape(-1, 3), "plane(+/-60%)"
    return None, f"no far-field parametrisation for {kind}"


# =============================================================================================
# geometry helpers
# =============================================================================================

def boundary_polyline(face, per_edge=65, max_segments=6000):
    """Dense 3D samples of every non-degenerate boundary edge, as segment endpoints."""
    A, B = [], []
    for _wire, _outer, wedges in C._face_wire_edges(face):
        for edge, _v in wedges:
            if BRep_Tool.Degenerated(edge):
                continue
            try:
                curve3d, first, last = BRep_Tool.Curve(edge)
            except Exception:
                curve3d = None
            if curve3d is None:
                continue
            ts = np.linspace(first, last, per_edge)
            P = np.array([[p.X(), p.Y(), p.Z()]
                          for p in (curve3d.Value(float(t)) for t in ts)])
            A.append(P[:-1])
            B.append(P[1:])
    if not A:
        return None, None
    A = np.concatenate(A)
    B = np.concatenate(B)
    if len(A) > max_segments:
        idx = np.linspace(0, len(A) - 1, max_segments).astype(int)
        A, B = A[idx], B[idx]
    return A, B


def dist_to_segments(P, A, B, chunk=4096):
    """Min distance from every row of P to the polyline segments A->B."""
    if A is None or len(A) == 0:
        return np.full(len(P), np.inf)
    d = B - A
    L2 = np.einsum('ij,ij->i', d, d)
    L2 = np.where(L2 > 0, L2, 1.0)
    out = np.empty(len(P))
    for s in range(0, len(P), chunk):
        Q = P[s:s + chunk]
        w = Q[:, None, :] - A[None, :, :]
        t = np.clip(np.einsum('ijk,jk->ij', w, d) / L2[None, :], 0.0, 1.0)
        proj = A[None, :, :] + t[..., None] * d[None, :, :]
        out[s:s + chunk] = np.linalg.norm(Q[:, None, :] - proj, axis=2).min(axis=1)
    return out


def _ancestor_faces(anc, edge):
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


def sign_vector_keys(F, eps):
    """Rows of F -> a hashable sign tuple, with |f| <= eps folded to the nearer non-zero sign."""
    S = np.where(F >= 0.0, 1, -1).astype(np.int8)
    return [tuple(int(x) for x in row) for row in S]


# =============================================================================================
# the per-face measurement
# =============================================================================================

def measure_face(face, fidx, anc, model_of, face_model, scale_to_cm, args, solid_shape=None):
    s = scale_to_cm
    surf = BRep_Tool.Surface(face)
    tol_face = float(BRep_Tool.Tolerance(face))

    rec = {
        "face": fidx,
        "selfKind": face_model["kind"],
        "selfSource": face_model["source"],
        "tolFaceCm": tol_face * s,
    }

    # ---- wires, edges, neighbours -----------------------------------------------------------
    wires = list(C._face_wire_edges(face))
    rec["nWires"] = len(wires)
    rec["nInnerWires"] = sum(0 if outer else 1 for _w, outer, _e in wires)
    nbr_models = []          # (model, wire_index, edge_index)
    n_edges = n_deg = n_seam = n_multi = n_freeform = 0
    nbr_kinds = Counter()
    for wi, (_wire, _outer, wedges) in enumerate(wires):
        for ei, (edge, _v) in enumerate(wedges):
            n_edges += 1
            if BRep_Tool.Degenerated(edge):
                n_deg += 1
                continue
            others = [g for g in _ancestor_faces(anc, edge) if not g.IsSame(face)]
            if len(others) == 0:
                n_seam += 1
                continue
            if len(others) > 1:
                n_multi += 1
            gm = model_of(others[0])
            if gm is None or gm["kind"] not in _EVALUABLE:
                n_freeform += 1
                nbr_kinds[None if gm is None else gm["kind"]] += 1
                continue
            nbr_kinds[gm["kind"]] += 1
            nbr_models.append((gm, wi, ei))
    rec.update({"nEdges": n_edges, "nDegenerateEdges": n_deg, "nSeamEdges": n_seam,
                "nMultiNeighbourEdges": n_multi, "nNonEvaluableNeighbours": n_freeform,
                "neighbourKinds": {str(k): v for k, v in nbr_kinds.items()}})

    if not nbr_models:
        rec["status"] = "no evaluable neighbour surface"
        return rec

    # ---- sample the face's SURFACE over its own (u, v) rectangle ------------------------------
    try:
        umin, umax, vmin, vmax = breptools.UVBounds(face)
    except Exception:
        rec["status"] = "UVBounds failed"
        return rec
    if not all(math.isfinite(x) for x in (umin, umax, vmin, vmax)) or umax <= umin or vmax <= vmin:
        rec["status"] = "degenerate (u, v) rectangle"
        return rec

    n = args.grid
    us = umin + (umax - umin) * (np.arange(n) + 0.5) / n
    vs = vmin + (vmax - vmin) * (np.arange(n) + 0.5) / n
    UU, VV = np.meshgrid(us, vs, indexing="ij")
    uvs = np.column_stack([UU.ravel(), VV.ravel()])

    P = np.empty((len(uvs), 3))
    for i, (u, v) in enumerate(uvs):
        p = surf.Value(float(u), float(v))
        P[i] = (p.X(), p.Y(), p.Z())

    clf = BRepTopAdaptor_FClass2d(face, max(tol_face, 1e-12))
    state = np.empty(len(uvs), dtype=np.int8)   # 1 IN, 0 ON, -1 OUT
    for i, (u, v) in enumerate(uvs):
        st = clf.Perform(gp_Pnt2d(float(u), float(v)))
        state[i] = 1 if st == TopAbs_IN else (0 if st == TopAbs_ON else -1)
    n_in = int((state == 1).sum())
    n_out = int((state == -1).sum())
    n_on = int((state == 0).sum())
    rec.update({"nSamples": len(uvs), "nIn": n_in, "nOut": n_out, "nOn": n_on})
    if n_in == 0:
        rec["status"] = "no interior sample found (patch too thin for this grid)"
        return rec

    # ---- boundary polyline, needed both for placement and for the tangency measurement --------
    A, B = boundary_polyline(face, per_edge=args.edge_samples)

    # ---- deduplicate trimming surfaces GEOMETRICALLY on the face's own samples ---------------
    raw = [m for m, _wi, _ei in nbr_models]
    Fraw = np.array([implicit(m, P) for m in raw])          # (Kraw, M)
    diag = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0)))
    dedup_tol = max(1e-9 * max(diag, 1.0), 1e-12)
    keep, groups = [], []
    for j in range(len(raw)):
        for gi, k in enumerate(keep):
            if float(np.abs(Fraw[j] - Fraw[k]).max()) <= dedup_tol:
                groups[gi].append(j)
                break
        else:
            keep.append(j)
            groups.append([j])
    allmodels = [raw[k] for k in keep]
    Fall = Fraw[keep]
    rec["nDistinctNeighbourSurfaces"] = len(allmodels)
    rec["nBoundaryEdgesWithNeighbour"] = len(raw)
    # how many distinct boundary edges each distinct surface bounds -- "a surface can bound the
    # face twice" is exactly `> 1` here
    rec["edgesPerNeighbourSurface"] = sorted((len(g) for g in groups), reverse=True)

    # ---- COINCIDENT neighbours: a neighbour lying on the face's OWN surface -------------------
    # A cylinder split by a seam, two coplanar faces joined by a union, a torus meeting a
    # tangent cylinder: the neighbour's implicit function is identically zero on the whole face,
    # so its sign carries no information at all and the edge it is supposed to trim has no
    # implicit representation. This is not a tolerance question -- it is a rank deficiency.
    coinc_tol = max(20.0 * tol_face, 1e-9 * max(diag, 1.0))
    maxabs = np.abs(Fall).max(axis=1)
    coincident = maxabs <= coinc_tol
    rec["nCoincidentTrimSurfaces"] = int(coincident.sum())
    rec["nEdgesOnCoincidentSurfaces"] = int(sum(len(groups[i]) for i in np.where(coincident)[0]))
    rec["coincidentMaxAbsFCm"] = float(maxabs[coincident].max() * s) if coincident.any() else 0.0

    # An edge is *unrepresented* by any co-surface trim when its neighbour has no evaluable
    # surface, when it is a seam (the face is its own neighbour), or when the neighbour lies on
    # the face's own surface. Across such an edge the implicit description simply does not
    # bound the face, so every over-acceptance there is explained before it is measured. A face
    # is `implicitComplete` only when that count is zero -- that is the population on which the
    # idea is even applicable, and every headline number is quoted for it separately.
    rec["nEdgesUnrepresented"] = int(n_freeform + n_seam + rec["nEdgesOnCoincidentSurfaces"])
    rec["implicitComplete"] = bool(rec["nEdgesUnrepresented"] == 0)

    eff = np.where(~coincident)[0]
    models = [allmodels[i] for i in eff]
    F = Fall[eff]                                           # (K, M)
    K = len(models)
    rec["nTrimSurfaces"] = K
    rec["trimSurfaceKinds"] = dict(Counter(m["kind"] for m in models))
    rec["trimSurfaceSources"] = dict(Counter(m["source"] for m in models))
    rec["edgesPerTrimSurface"] = sorted((len(groups[i]) for i in eff), reverse=True)
    rec["nSurfacesBoundingTwice"] = int(sum(1 for i in eff if len(groups[i]) > 1))
    if K == 0:
        rec["status"] = "every neighbour lies on the face's own surface"
        return rec

    # ---- transversality: the angle between the two surfaces along the shared edge -------------
    # The sign test's conditioning is |sin(angle)|; at a tangency it degenerates and the test is
    # asked precisely where it is worst.
    # Only samples that genuinely lie ON the trimming surface are admissible: a neighbour plane
    # parallel to the face never touches its boundary, and taking "the closest 5%" there would
    # report a tangency that does not exist. A surface with no such sample is reported as
    # unmeasured rather than as tangent.
    if A is not None:
        E = np.concatenate([A, B[-1:]])
        nface = gradient(face_model, E)
        on_tol = max(100.0 * tol_face, 1e-7 * max(diag, 1.0))
        sins, unmeasured = [], 0
        for k in range(K):
            fk = np.abs(implicit(models[k], E))
            on = fk <= on_tol
            if on.sum() < 3:
                unmeasured += 1
                sins.append(None)
                continue
            nk = gradient(models[k], E[on])
            c = np.abs(np.einsum('ij,ij->i', nface[on], nk)).clip(0.0, 1.0)
            sins.append(float(np.sqrt(1.0 - c ** 2).min()))
        got = [x for x in sins if x is not None]
        rec["minTransversalitySin"] = min(got) if got else None
        rec["nTransversalityUnmeasured"] = unmeasured
        rec["transversalitySin"] = [None if x is None else round(x, 6) for x in sins]

    if args.perturb_radius:
        models = [perturb(m, args.perturb_radius / s) for m in models]
        F = np.array([implicit(m, P) for m in models])

    # ---- the sign vector / arrangement-cell measurement --------------------------------------
    keys = sign_vector_keys(F.T, 0.0)
    inmask = state == 1
    outmask = state == -1
    cells_in = Counter(k for k, m in zip(keys, inmask) if m)
    cells_out = Counter(k for k, m in zip(keys, outmask) if m)
    rec["cellsIn"] = len(cells_in)
    rec["cellsInHistogram"] = sorted(cells_in.values(), reverse=True)
    # A cell holding a handful of samples can be a sliver produced by a trimming surface that
    # merely grazes the interior; a cell holding 1% of the interior is a real lobe. Both are
    # reported, because the first decides how many DNF terms an exact description needs and the
    # second decides whether the shape is genuinely non-convex.
    rec["cellsInMajor"] = int(sum(1 for v in cells_in.values() if v >= max(3, 0.01 * n_in)))
    leak_keys = set(cells_in) & set(cells_out)
    n_leak = int(sum(cells_out[k] for k in leak_keys))
    rec["nLeakSamples"] = n_leak
    rec["nLeakCells"] = len(leak_keys)

    # ---- boundary distance, so a disagreement can be placed ----------------------------------
    dbnd = dist_to_segments(P, A, B) * s

    if n_leak:
        leakmask = np.array([(k in leak_keys) for k in keys]) & outmask
        rec["leakMaxBoundaryDistCm"] = float(dbnd[leakmask].max())
        rec["leakMedBoundaryDistCm"] = float(np.median(dbnd[leakmask]))
    else:
        rec["leakMaxBoundaryDistCm"] = 0.0
        rec["leakMedBoundaryDistCm"] = 0.0

    # ---- the naive conjunction, with the sense read off the deepest interior sample -----------
    interior = np.where(inmask)[0]
    deepest = interior[int(np.argmax(dbnd[interior]))]
    rec["deepestInteriorDistCm"] = float(dbnd[deepest])
    sense = np.where(F[:, deepest] >= 0.0, 1.0, -1.0)
    if args.flip_sense is not None and K > 0:
        sense[args.flip_sense % K] *= -1.0
    rec["sense"] = [int(x) for x in sense]

    G = sense[:, None] * F                                   # (K, M), >= 0 means "inside side"
    rule_in = np.all(G >= 0.0, axis=0)
    fp = rule_in & outmask
    fn = (~rule_in) & inmask
    rec["nFalsePositive"] = int(fp.sum())
    rec["nFalseNegative"] = int(fn.sum())
    for tag, mask in (("fp", fp), ("fn", fn)):
        if mask.any():
            rec[tag + "MaxBoundaryDistCm"] = float(dbnd[mask].max())
            rec[tag + "MedBoundaryDistCm"] = float(np.median(dbnd[mask]))
            # depth: how far on the wrong side of the deciding surface, in cm
            if tag == "fp":
                depth = np.min(G[:, mask], axis=0) * s
            else:
                depth = -np.min(G[:, mask], axis=0) * s
            rec[tag + "MaxDepthCm"] = float(np.abs(depth).max())
        else:
            rec[tag + "MaxBoundaryDistCm"] = 0.0
            rec[tag + "MedBoundaryDistCm"] = 0.0
            rec[tag + "MaxDepthCm"] = 0.0

    # ---- independent verification of the ground truth on the disagreeing points --------------
    # `BRepTopAdaptor_FClass2d` decides in the face's (u, v) domain. `BRepExtrema_DistShapeShape`
    # decides in 3D against the trimmed face as a shape and shares no code with it. A point the
    # conjunction accepts and FClass2d calls OUT must be a measurable distance from the face; if
    # it is not, the disagreement is FClass2d's and not the rule's.
    if args.verify:
        from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        idx = np.where(fp)[0]
        if len(idx) > args.verify:
            idx = idx[np.linspace(0, len(idx) - 1, args.verify).astype(int)]
        dists = []
        for i in idx:
            v = BRepBuilderAPI_MakeVertex(gp_Pnt(*[float(x) for x in P[i]])).Vertex()
            ext = BRepExtrema_DistShapeShape(v, face)
            if ext.IsDone() and ext.NbSolution() > 0:
                dists.append(float(ext.Value()) * s)
        if dists:
            rec["fpVerifiedMinDistToFaceCm"] = float(min(dists))
            rec["fpVerifiedMedDistToFaceCm"] = float(np.median(dists))
            rec["fpVerifiedN"] = len(dists)
            # agreement: a genuine OUT point must be further from the face than the face's own
            # declared tolerance
            rec["fpVerifiedAgree"] = int(sum(1 for d in dists if d > 10.0 * tol_face * s))

    # ---- per-surface sense consistency: does any single surface cut the interior? -------------
    amb = 0
    amb_worst = 0.0
    for k in range(K):
        gk = sense[k] * F[k][inmask]
        if gk.min() < 0.0:
            amb += 1
            amb_worst = max(amb_worst, float(-gk.min()) * s)
    rec["nAmbiguousSenseSurfaces"] = amb
    rec["ambiguousWorstDepthCm"] = amb_worst

    # ---- tangency conditioning: |grad f| is 1 for all these kinds, so the conditioning number
    # is the angle between the face's surface and the trimming surface along the shared edge.
    # Measured as the smallest |f_k| gradient mismatch is not needed; instead measure how slowly
    # f_k grows away from the boundary, which is what an ill-conditioned sign test looks like.
    near = dbnd <= max(10.0 * tol_face * s, 1e-9)
    rec["nSamplesWithinTolOfBoundary"] = int(near.sum())

    # ---- far field ---------------------------------------------------------------------------
    if args.far_grid > 0:
        Q, flabel = far_field_points(face_model, P, args.far_grid, args.far_grid)
        rec["farFieldParam"] = flabel
        if Q is not None:
            FQ = np.array([implicit(m, Q) for m in models])
            GQ = sense[:, None] * FQ
            acc = np.all(GQ >= 0.0, axis=0)
            rec["nFarField"] = int(len(Q))
            rec["nFarFieldAccepted"] = int(acc.sum())
            # of the accepted ones, how many are genuinely NOT on the face?
            idx = np.where(acc)[0]
            if len(idx) > args.far_check:
                idx = idx[np.linspace(0, len(idx) - 1, args.far_check).astype(int)]
            sas = ShapeAnalysis_Surface(surf)
            bad = []
            for i in idx:
                p = gp_Pnt(float(Q[i][0]), float(Q[i][1]), float(Q[i][2]))
                uv = sas.ValueOfUV(p, max(tol_face, 1e-7))
                q = surf.Value(uv.X(), uv.Y())
                resid = math.dist((q.X(), q.Y(), q.Z()), (p.X(), p.Y(), p.Z()))
                if resid > max(1e-6 * max(diag, 1.0), 10.0 * tol_face):
                    onface = False
                else:
                    st = clf.Perform(gp_Pnt2d(uv.X(), uv.Y()))
                    onface = (st != TopAbs_OUT)
                if not onface:
                    bad.append(i)
            rec["nFarFieldChecked"] = int(len(idx))
            rec["nFarFieldSpurious"] = len(bad)
            if bad:
                db = dist_to_segments(Q[np.array(bad)], A, B) * s
                rec["farFieldSpuriousMaxBoundaryDistCm"] = float(db.max())
                rec["farFieldSpuriousMedBoundaryDistCm"] = float(np.median(db))
            else:
                rec["farFieldSpuriousMaxBoundaryDistCm"] = 0.0
                rec["farFieldSpuriousMedBoundaryDistCm"] = 0.0

    # ---- optional independent cross-check of the ground truth ---------------------------------
    if args.solid_crosscheck and solid_shape is not None:
        cls = BRepClass3d_SolidClassifier(solid_shape)
        step = max(1, len(P) // args.solid_crosscheck)
        agree = disagree = 0
        for i in range(0, len(P), step):
            p = gp_Pnt(*[float(x) for x in P[i]])
            cls.Perform(p, max(tol_face, 1e-9))
            st = cls.State()
            # a point on the trimmed face must read ON for the solid; a point on the surface but
            # outside the face must not.
            onface = state[i] == 1
            ison = (st == TopAbs_ON)
            if onface == ison:
                agree += 1
            else:
                disagree += 1
        rec["solidCrosscheckAgree"] = agree
        rec["solidCrosscheckDisagree"] = disagree

    rec["status"] = "ok"
    return rec


# =============================================================================================
# drive
# =============================================================================================

def measure_solid(lid, name, shape, scale_to_cm, args, census_faces=None):
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
        m = TC.stored_model(face)
        if m is None:
            adaptor = BRepAdaptor_Surface(face)
            try:
                uv = breptools.UVBounds(face)
                rec = TC.recognize(adaptor, uv)
            except Exception:
                rec = None
            m = TC.recognized_model(rec) if rec is not None else None
        model_cache[idx] = m
        return m

    out = {"lid": lid, "name": name, "nFaces": len(faces), "faces": []}
    for i, face in enumerate(faces):
        if args.faces is not None and i not in args.faces:
            continue
        fm = model_of(face)
        if fm is None or fm["kind"] not in _EVALUABLE:
            continue
        # which population is this face in?
        adaptor = BRepAdaptor_Surface(face)
        stored_type = C._SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")
        r = measure_face(face, i, anc, model_of, fm, scale_to_cm, args,
                         solid_shape=shape if args.solid_crosscheck else None)
        r["storedType"] = stored_type
        r["population"] = ("candidate" if fm["source"] == "recognized" else "control")
        if census_faces is not None:
            r["censusRejectedEdges"] = census_faces.get(i)
        out["faces"].append(r)
    return out


def run_model(step_path, label, args):
    unit = C.detect_step_length_unit(step_path)
    scale = C.step_unit_scale_to_cm(unit)
    C.extract_graph(step_path, meshparam=None, scale_to_cm=scale)
    solids = []
    for lid, shape in C.def_shapes.items():
        name = C.def_names.get(lid, lid)
        if args.solids and name not in args.solids:
            continue
        t0 = time.time()
        s = measure_solid(lid, name, shape, scale, args)
        s["seconds"] = time.time() - t0
        solids.append(s)
        if args.verbose:
            nf = len(s["faces"])
            bad = sum(1 for f in s["faces"]
                      if f.get("status") == "ok" and (f.get("nFalsePositive", 0)
                                                      or f.get("nFalseNegative", 0)))
            print(f"  [{label}] {name}: {nf} faces measured, {bad} with disagreements "
                  f"({s['seconds']:.1f} s)", flush=True)
    return {"model": step_path, "label": label, "unit": unit, "scaleToCm": scale,
            "solids": solids}


def build_fixtures(tmpdir):
    script = _HERE.parent / "make_boolean_fixtures.py"
    subprocess.run([sys.executable, str(script), "--outdir", str(tmpdir)],
                   check=True, stdout=subprocess.DEVNULL)
    return sorted(Path(tmpdir).glob("*.step")) + sorted(Path(tmpdir).glob("*.stp"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", action="append", default=[],
                    help="STEP file, optionally LABEL:path")
    ap.add_argument("--fixtures", action="store_true")
    ap.add_argument("--solids", default=None, help="comma-separated leaf-solid names")
    ap.add_argument("--faces", default=None, help="comma-separated face indices (single solid)")
    ap.add_argument("--grid", type=int, default=40, help="near-field grid, N x N samples per face")
    ap.add_argument("--far-grid", type=int, default=24, help="far-field grid (0 disables)")
    ap.add_argument("--far-check", type=int, default=120,
                    help="max far-field accepted points verified against OCCT per face")
    ap.add_argument("--edge-samples", type=int, default=65)
    ap.add_argument("--flip-sense", type=int, default=None,
                    help="NEGATIVE CONTROL: invert stored sense #k on every face")
    ap.add_argument("--perturb-radius", type=float, default=0.0,
                    help="NEGATIVE CONTROL: move every trimming surface by this many cm")
    ap.add_argument("--verify", type=int, default=0,
                    help="re-check up to this many false-positive samples per face against "
                         "BRepExtrema_DistShapeShape (an independent OCCT ground truth)")
    ap.add_argument("--solid-crosscheck", type=int, default=0,
                    help="cross-check the FClass2d ground truth against BRepClass3d_SolidClassifier "
                         "on this many samples per face")
    ap.add_argument("--json", default=None)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    args.solids = set(x for x in (args.solids or "").split(",") if x) or None
    args.faces = set(int(x) for x in (args.faces or "").split(",") if x) or None

    results = []
    tmp = None
    models = []
    for spec in args.model:
        if ":" in spec and not os.path.exists(spec):
            lab, path = spec.split(":", 1)
        else:
            lab, path = Path(spec).stem, spec
        models.append((lab, path))
    if args.fixtures:
        tmp = tempfile.mkdtemp(prefix="implicitTrimFixtures")
        for p in build_fixtures(tmp):
            models.append((f"fixture:{p.stem}", str(p)))

    for lab, path in models:
        if args.verbose:
            print(f"[{lab}] {path}", flush=True)
        results.append(run_model(path, lab, args))

    payload = {
        "args": {k: (sorted(v) if isinstance(v, set) else v) for k, v in vars(args).items()},
        "models": results,
    }
    if args.json:
        Path(args.json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json).write_text(json.dumps(payload, indent=1, default=str))
        print(f"wrote {args.json}")
    summarise(payload)


def summarise(payload):
    tot = Counter()
    worst = []
    for m in payload["models"]:
        for s in m["solids"]:
            for f in s["faces"]:
                tot["faces"] += 1
                st = f.get("status")
                tot["status:" + str(st)] += 1
                if st != "ok":
                    continue
                tot["measured"] += 1
                tot["samples"] += f.get("nSamples", 0)
                tot["in"] += f.get("nIn", 0)
                tot["out"] += f.get("nOut", 0)
                tot["fp"] += f.get("nFalsePositive", 0)
                tot["fn"] += f.get("nFalseNegative", 0)
                tot["leak"] += f.get("nLeakSamples", 0)
                tot["trimSurfaces"] += f.get("nTrimSurfaces", 0)
                if f.get("cellsIn", 1) > 1:
                    tot["facesMultiCell"] += 1
                if f.get("nLeakSamples", 0):
                    tot["facesLeaking"] += 1
                if f.get("nOut", 0) == 0:
                    tot["facesWithNoOutsideSample"] += 1
                if f.get("nInnerWires", 0):
                    tot["facesWithHoles"] += 1
                if f.get("nCoincidentTrimSurfaces", 0):
                    tot["facesWithCoincidentNeighbour"] += 1
                if f.get("nNonEvaluableNeighbours", 0):
                    tot["facesWithNonEvaluableNeighbour"] += 1
                if f.get("nSurfacesBoundingTwice", 0):
                    tot["facesWithASurfaceBoundingTwice"] += 1
                if f.get("nAmbiguousSenseSurfaces", 0):
                    tot["facesWithAmbiguousSense"] += 1
                if (f.get("minTransversalitySin") or 1.0) < 1e-3:
                    tot["facesWithTangentNeighbour"] += 1
                if f.get("nFalsePositive", 0) or f.get("nFalseNegative", 0):
                    tot["facesDisagreeing"] += 1
                    worst.append((max(f.get("fpMaxDepthCm", 0.0), f.get("fnMaxDepthCm", 0.0)),
                                  s["name"], f["face"], f.get("nFalsePositive", 0),
                                  f.get("nFalseNegative", 0), f.get("cellsIn"),
                                  f.get("nTrimSurfaces")))
                if f.get("nFarFieldSpurious", 0):
                    tot["facesFarFieldSpurious"] += 1
                    tot["farFieldSpurious"] += f.get("nFarFieldSpurious", 0)
    print("\n=== implicit-trim validation ===")
    for k in sorted(tot):
        print(f"  {k:34s} {tot[k]}")
    if worst:
        print("\n  worst faces by disagreement depth (cm):")
        for d, name, fi, fp, fn, cells, K in sorted(worst, reverse=True)[:25]:
            print(f"    {name:22s} f#{fi:<4d} depth {d:.3e}  FP {fp:<5d} FN {fn:<5d} "
                  f"cellsIn {cells} K {K}")


if __name__ == "__main__":
    main()
