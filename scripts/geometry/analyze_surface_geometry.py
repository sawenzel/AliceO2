#!/usr/bin/env python3
"""
Classify the *actual differential geometry* of every face in a STEP file, independently of the
surface type the file happens to store.

Why this exists
---------------
The stored STEP surface type is a statement about the exporter, not about the geometry. CAD
kernels routinely write an exact cylinder, cone or sphere as a *rational* B-spline/Bezier patch:
that representation is exact, not an approximation, because a rational quadratic/cubic reproduces
conics exactly. A converter that dispatches on `BRepAdaptor_Surface::GetType()` therefore throws
away analytic geometry it fully supports.

Measured on 2026-07-26 with this tool:
  - as1-oc-214.stp        : all 70 "bspline" faces are EXACT cylinders (R=5.000000000 mm, max
                            radial deviation 1.1e-11 mm). 0/5 solids convertible -> 5/5.
  - ALICE3_CAD_pure.step  : 2889 of 8664 non-analytic faces (33%) are exact quadrics
                            (1419 cylinder, 1398 cone, 72 sphere). 15/55 solids -> 29/55.

The classifier's "already analytic as stored" per-solid count reproduces the exact converter's
independently measured coverage, which cross-validates it.

Method
------
Model selection, not fitting: candidate models are tried in increasing parameter count and the
FIRST one that fits at machine precision is accepted.

    plane (3) < sphere (4) < cylinder (5) < cone (6) < free-form

Each test is a small linear solve on the sampled surface *normal field*, so no initial guess and
no iteration is needed:

    plane    : all unit normals parallel
    sphere   : normal lines concurrent  -> solve P_i = C + r*N_i for (C, r)
    cylinder : normals coplanar; axis = smallest right singular vector of the normal matrix,
               then project out the axis and fit a circle, requiring constant radius
    cone     : N_i . (P_i - A) = 0 is LINEAR in the apex A; solve, then require a constant
               half-angle about the mean ruling direction

Every candidate is scored by ONE quantity -- the largest distance from a sampled point to the
candidate surface, relative to the patch size -- so "exact" means machine precision. This is
deliberate: an "almost cylinder" that is really free-form must stay free-form, otherwise a
converter acting on this classification would silently change the geometry it exists to represent
exactly. It is also load-bearing rather than tidy: until 2026-08-02 the plane and the cone were
scored by an *angle* instead, and on ALICE3 that accepted 184 faces as cones that miss their own
surface by up to 79 cm (Stream_K_Tier0.md sec 3).

Known limitation: there is no torus test yet, so torus patches stored as NURBS are reported as
free-form. Every "free-form" count is therefore an UPPER bound and every analytic count a LOWER
bound.

Note: OCCT's own BRepLib_CanonicalRecognition (7.7+) is not exposed in the pythonOCC v7.9.3 build
used here (ImportError from OCC.Core.BRepLib), which is why this is a standalone numeric recognizer.
`ShapeAnalysis_CanonicalRecognition` *is* exposed, however, and `csg/census.py` uses it; the two
recognisers agree face for face on ALICE3's cylinders, cones and spheres once this one's acceptance
is a measured gap (Stream_K_Tier0.md sec 3.1), which is the strongest available cross-check on both.

Trim curves
-----------
`--trim-curves` asks the same question one dimension down, about the *boundary* curves, and
additionally measures the shared-edge defect written up in ExactTrimTopology.md: how far apart the
two faces of one `TopoDS_Edge` place their common boundary, given that each carries its own
independently fitted pcurve.

Measured on 2026-07-26 (faces the exact converter can represent; "real geometry" is the
machine-precision classification of curves the file stores as B-splines):
  - Bagger.step     : 48 B-spline pcurves -> 6 exactly lines, 0 circles, 42 free-form.
                      705 shared edges, max pcurve disagreement 1.3e-5 model units (mean 7.3e-8).
  - as1-oc-214.stp  : 210 B-spline pcurves -> 70 lines, 70 circles, 70 free-form.
                      354 shared edges, max pcurve disagreement 2.9e-5 (mean 5.1e-6).

Usage
-----
    python3 analyze_surface_geometry.py model.step [more.step ...]
    python3 analyze_surface_geometry.py --per-solid model.step      # coverage forecast
    python3 analyze_surface_geometry.py --trim-curves model.step    # trim-curve + shared-edge report
    python3 analyze_surface_geometry.py --json out.json model.step

Author:
  - Added 2026-07-26 alongside the BVHSurfaceSolid exact-surface work (see BVHSurfaceSolid.md,
    milestone "Canonical-form recognition").
"""

import argparse
import json
import math
import sys
from collections import Counter

import numpy as np

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Core.GeomAbs import (
    GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus,
    GeomAbs_BSplineSurface, GeomAbs_BezierSurface,
    GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion,
    GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse, GeomAbs_Hyperbola, GeomAbs_Parabola,
    GeomAbs_BSplineCurve, GeomAbs_BezierCurve, GeomAbs_OtherCurve,
)
from OCC.Core.gp import gp_Pnt, gp_Vec

# Surface types the exact-surface converter already handles natively.
ANALYTIC_TYPES = (GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus)

STORED_TYPE_NAMES = {
    GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Cone: "cone",
    GeomAbs_Sphere: "sphere", GeomAbs_Torus: "torus",
    GeomAbs_BSplineSurface: "bspline", GeomAbs_BezierSurface: "bezier",
    GeomAbs_SurfaceOfRevolution: "revolution", GeomAbs_SurfaceOfExtrusion: "extrusion",
}

CURVE_TYPE_NAMES = {
    GeomAbs_Line: "line", GeomAbs_Circle: "circle", GeomAbs_Ellipse: "ellipse",
    GeomAbs_Hyperbola: "hyperbola", GeomAbs_Parabola: "parabola",
    GeomAbs_BSplineCurve: "bspline", GeomAbs_BezierCurve: "bezier",
    GeomAbs_OtherCurve: "other",
}

# Curve types the converter's trim path already stores analytically (line/arc Curve2D kinds);
# everything else becomes a flattened B-spline.
CANONICAL_CURVE_TYPES = (GeomAbs_Line, GeomAbs_Circle)

# Relative residual below which a model is accepted as EXACT (machine precision).
TOL_EXACT = 1e-9
# A deliberately looser bound, reported alongside so the gap between "exact" and "close" is visible.
TOL_LOOSE = 1e-6

# Retired: the plane test used to be "normals within this of parallel", an ANGLE compared against
# the same bound as the distance-valued residuals of the other models. It is now the measured
# distance to the plane, like every other candidate (see `classify_surface`). Kept only because
# older JSON reports quote it.
TOL_PLANE_NORMALS = 1e-12
# Max |N . axis| for the normal field to count as coplanar (cylinder candidate).
TOL_CYLINDER_COPLANAR = 1e-9


def sample_surface(adaptor, n=9):
    """Sample an (n x n) parametric grid, returning (points, unit normals) or (None, None)."""
    u0, u1 = adaptor.FirstUParameter(), adaptor.LastUParameter()
    v0, v1 = adaptor.FirstVParameter(), adaptor.LastVParameter()
    if not all(map(math.isfinite, (u0, u1, v0, v1))):
        return None, None
    points, normals = [], []
    for u in np.linspace(u0, u1, n):
        for v in np.linspace(v0, v1, n):
            p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
            try:
                adaptor.D1(float(u), float(v), p, du, dv)
            except Exception:
                return None, None
            nrm = np.cross([du.X(), du.Y(), du.Z()], [dv.X(), dv.Y(), dv.Z()])
            length = np.linalg.norm(nrm)
            if length < 1e-14:      # parametric degeneracy (pole/seam), skip this sample
                continue
            points.append([p.X(), p.Y(), p.Z()])
            normals.append(nrm / length)
    if len(points) < 12:
        return None, None
    return np.array(points), np.array(normals)


def surface_gap(kind, model, P):
    """The largest distance from a sampled surface point to the candidate surface, in CAD units.

    ONE quantity for every candidate kind. See `classify_surface` for why that matters and
    `scripts/geometry/Stream_K_Tier0.md` for what scoring them differently cost.
    """
    if kind == "plane":
        n = np.asarray(model["normal"], dtype=float)
        n = n / np.linalg.norm(n)
        return float(np.abs((P - np.asarray(model["point"], dtype=float)) @ n).max())
    if kind == "sphere":
        return float(np.abs(np.linalg.norm(P - model["centre"], axis=1) - model["radius"]).max())
    if kind == "cylinder":
        a = np.asarray(model["axis"], dtype=float)
        a = a / np.linalg.norm(a)
        radial = P - model["origin"]
        radial = radial - np.outer(radial @ a, a)
        return float(np.abs(np.linalg.norm(radial, axis=1) - model["radius"]).max())
    if kind == "cone":
        a = np.asarray(model["axis"], dtype=float)
        a = a / np.linalg.norm(a)
        rel = P - model["apex"]
        h = rel @ a
        r = np.linalg.norm(rel - np.outer(h, a), axis=1)
        half = model["half_angle"]
        return float(np.abs(r * math.cos(half) - h * math.sin(half)).max())
    return float("inf")


def classify_surface(P, N):
    """Return (model_name, relative_gap) for the smallest model that fits the samples.

    Every candidate is scored by exactly one quantity: `surface_gap` divided by the sample
    bounding-box diagonal. The linear solves *propose* a model -- that is what they are good at,
    and no initial guess is needed -- but none of them decides, because their own residuals are
    not the same quantity: the sphere's and the cylinder's are distances, the plane's and the
    cone's were angles. Measured on ALICE3, scoring the cone by an angle accepted 184 faces whose
    recognized cone misses the real surface by up to 79 cm, at an internal residual of 6.7e-10
    against a 1e-9 bound: on a patch whose normal field is nearly rank-2 (a swept free-form
    profile), the least-squares apex runs off to 1e11 patch diagonals, the half-angle collapses to
    zero, and both cone tests then pass vacuously. `Stream_K_Tier0.md` sec 3 has the anatomy;
    `O2_CADtoTGeo.py --self-test` has the closed-form reproducer and the regression guard.
    """
    scale = np.linalg.norm(P.max(axis=0) - P.min(axis=0))
    if scale < 1e-12:
        return "degenerate", 0.0

    def score(kind, model):
        gap = surface_gap(kind, model, P)
        return gap / scale if math.isfinite(gap) else float(np.inf)

    # --- plane (3 parameters): all unit normals parallel
    plane = {"normal": N[0] / np.linalg.norm(N[0]), "point": P[0]}
    plane_res = score("plane", plane)
    if plane_res < TOL_EXACT:
        return "plane", float(plane_res)

    best = ("freeform", float(np.inf))

    # --- sphere (4): normal lines concurrent, P_i = C + r*N_i
    A = np.zeros((3 * len(P), 4))
    b = np.zeros(3 * len(P))
    for i in range(len(P)):
        A[3 * i:3 * i + 3, 0:3] = np.eye(3)
        A[3 * i:3 * i + 3, 3] = N[i]
        b[3 * i:3 * i + 3] = P[i]
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    res = score("sphere", {"centre": sol[:3], "radius": abs(sol[3])})
    if res < best[1]:
        best = ("sphere", res)

    # --- cylinder (5): normals coplanar; axis = null direction of the normal field
    _, _, Vt = np.linalg.svd(N, full_matrices=False)
    axis = Vt[-1]
    if np.abs(N @ axis).max() < TOL_CYLINDER_COPLANAR:
        e1 = Vt[0]
        e2 = np.cross(axis, e1)
        x, y = P @ e1, P @ e2
        M = np.column_stack([x, y, np.ones_like(x)])
        D, E, F = np.linalg.lstsq(M, -(x ** 2 + y ** 2), rcond=None)[0]
        cx, cy = -D / 2, -E / 2
        r2 = cx * cx + cy * cy - F
        if r2 > 0:
            res = score("cylinder", {"axis": axis, "origin": cx * e1 + cy * e2,
                                     "radius": math.sqrt(r2)})
            if res < best[1]:
                best = ("cylinder", res)

    # --- cone (6): N_i . (P_i - A) = 0 is linear in the apex A
    apex, *_ = np.linalg.lstsq(N, np.einsum('ij,ij->i', N, P), rcond=None)
    d = P - apex
    dn = np.linalg.norm(d, axis=1)
    ok = dn > 1e-12
    if ok.sum() > 10:
        u = d[ok] / dn[ok, None]
        mean_dir = u.mean(axis=0)
        _, _, Vt2 = np.linalg.svd(u - mean_dir, full_matrices=False)
        ax2 = np.cross(Vt2[0], Vt2[1])
        n2 = np.linalg.norm(ax2)
        if n2 > 1e-12:                       # a ruling axis exists
            ax2 = ax2 / n2
            if np.dot(mean_dir, ax2) < 0.0:
                ax2 = -ax2
            half_angle = float(np.arccos(np.clip(np.abs(u @ ax2), -1.0, 1.0)).mean())
            res = score("cone", {"axis": ax2, "apex": apex, "half_angle": half_angle})
            if res < best[1]:
                best = ("cone", res)

    return best


# ---------------------------------------------------------------------------------------------
# Trim curves
#
# The surface half of this tool asks "is a stored NURBS *surface* really a quadric". This half
# asks the same question one dimension down, about the *trim* curves, and it exists because of the
# defect written up in ExactTrimTopology.md: the converter stores every curved trim edge as a
# B-spline and the kernel then flattens it to a polyline to answer point-in-trim. When a trim edge
# is really a circle or a straight line, storing it as such is exact recognition (same curve,
# lighter representation) and makes both adjacent faces agree analytically -- so the number below
# sizes items 2 and 3 of that plan before either is written.
#
# Two classifications per edge, and they answer different questions:
#   - the *pcurve* in the face's own (u, v): what this face could store analytically today;
#   - the *3D curve*: what a shared-edge representation (item 1) could hand to both faces.
# ---------------------------------------------------------------------------------------------

def sample_curve(evaluate, first, last, n=25):
    """Sample `evaluate(t)` (returning a coordinate tuple) uniformly over [first, last]."""
    if not (math.isfinite(first) and math.isfinite(last)) or last <= first:
        return None
    try:
        return np.array([evaluate(float(t)) for t in np.linspace(first, last, n)])
    except Exception:
        return None


def classify_curve(S):
    """Return (model_name, relative_residual) for a sampled curve in 2D or 3D.

    Model selection in increasing parameter count, as for surfaces: line < circle < freeform.
    The residual is relative to the sample bounding-box diagonal, so "exact" means the stored
    representation *is* that curve rather than merely close to it -- a curve that is almost a
    circle must stay free-form for exactly the reason the surface half states: silently changing
    geometry is worse than a slow trim.
    """
    scale = float(np.linalg.norm(S.max(axis=0) - S.min(axis=0)))
    if scale < 1e-12:
        return "degenerate", 0.0

    # --- line: every sample on the chord through the endpoints
    direction = S[-1] - S[0]
    length = float(np.linalg.norm(direction))
    if length > 1e-12:
        unit = direction / length
        offsets = S - S[0]
        perpendicular = offsets - np.outer(offsets @ unit, unit)
        residual = float(np.linalg.norm(perpendicular, axis=1).max() / scale)
        if residual < TOL_EXACT:
            return "line", residual
        best = ("freeform", residual if residual < np.inf else np.inf)
    else:
        best = ("freeform", float(np.inf))

    # --- circle: planar (3D only) and at constant distance from a centre in that plane.
    # Solving |P - C|^2 = R^2 as the linear system 2 P.C + (R^2 - |C|^2) = |P|^2 keeps this a
    # single least-squares solve with no initial guess, in the plane of the samples.
    centred = S - S.mean(axis=0)
    if S.shape[1] == 3:
        _, singular, Vt = np.linalg.svd(centred, full_matrices=False)
        if singular[2] > TOL_EXACT * scale:
            return best                              # not planar: cannot be a circle
        basis = Vt[:2]
        flat = centred @ basis.T
    else:
        basis, flat = None, centred
    M = np.column_stack([2.0 * flat, np.ones(len(flat))])
    solution, *_ = np.linalg.lstsq(M, np.einsum('ij,ij->i', flat, flat), rcond=None)
    centre2d = solution[:-1]
    r2 = solution[-1] + float(centre2d @ centre2d)
    if r2 > 0:
        radius = math.sqrt(r2)
        residual = float(np.abs(np.linalg.norm(flat - centre2d, axis=1) - radius).max() / scale)
        if residual < best[1]:
            best = ("circle", residual)
    return best


def _pcurve_samples(edge, face):
    curve2d, first, last = BRep_Tool.CurveOnSurface(edge, face)
    if curve2d is None:
        return None, None
    stored = None
    try:
        stored = Geom2dAdaptor_Curve(curve2d).GetType()
    except Exception:
        pass
    def evaluate(t):
        p = curve2d.Value(t)
        return (p.X(), p.Y())
    return stored, sample_curve(evaluate, first, last)


def _pcurve_polyline_3d(edge, face, n=48):
    """The edge as *this face* sees it: its pcurve pushed through this face's own surface.

    This is the geometry the converter actually writes today, so the distance between the two
    polylines of one shared edge is precisely the seam gap ExactTrimTopology.md is about.
    """
    curve2d, first, last = BRep_Tool.CurveOnSurface(edge, face)
    if curve2d is None or not (math.isfinite(first) and math.isfinite(last)) or last <= first:
        return None
    surface = BRep_Tool.Surface(face)
    if surface is None:
        return None
    try:
        points = []
        for t in np.linspace(first, last, n):
            uv = curve2d.Value(float(t))
            p = surface.Value(uv.X(), uv.Y())
            points.append((p.X(), p.Y(), p.Z()))
        return np.array(points)
    except Exception:
        return None


def _max_point_to_polyline(points, polyline):
    """One-sided max distance from `points` to the segments of `polyline`."""
    starts, ends = polyline[:-1], polyline[1:]
    segments = ends - starts
    lengths2 = np.einsum('ij,ij->i', segments, segments)
    lengths2 = np.where(lengths2 > 1e-30, lengths2, 1.0)
    worst = 0.0
    for point in points:
        offsets = point - starts
        t = np.clip(np.einsum('ij,ij->i', offsets, segments) / lengths2, 0.0, 1.0)
        closest = starts + t[:, None] * segments
        worst = max(worst, float(np.linalg.norm(point - closest, axis=1).min()))
    return worst


def analyze_trim_curves(shape, convertible_surfaces_only=True):
    """Classify every trim edge of every face, and measure the shared-edge pcurve disagreement.

    With `convertible_surfaces_only` (the default) only faces whose *surface* the exact converter
    can represent are counted -- a trim curve on a face that falls back to the mesh anyway is not
    a trim the kernel will ever evaluate, so counting it would inflate the payoff.
    """
    pcurve_stored, pcurve_real = Counter(), Counter()
    curve3d_stored, curve3d_real = Counter(), Counter()
    residuals = []

    faces = [f for f in faces_of(shape)]
    keep = []
    for face in faces:
        if not convertible_surfaces_only:
            keep.append(face)
            continue
        adaptor = BRepAdaptor_Surface(face)
        if adaptor.GetType() in ANALYTIC_TYPES:
            keep.append(face)
            continue
        P, N = sample_surface(adaptor)               # recognized-as-quadric faces convert too
        if P is not None:
            kind, res = classify_surface(P, N)
            if res < TOL_EXACT and kind in ("plane", "cylinder", "cone", "sphere"):
                keep.append(face)

    for face in keep:
        explorer = TopExp_Explorer(face, TopAbs_EDGE)
        while explorer.More():
            edge = topods.Edge(explorer.Current())
            explorer.Next()
            if BRep_Tool.Degenerated(edge):
                continue

            stored2d, samples2d = _pcurve_samples(edge, face)
            name2d = CURVE_TYPE_NAMES.get(stored2d, "none" if stored2d is None else f"enum{stored2d}")
            pcurve_stored[name2d] += 1
            if stored2d in CANONICAL_CURVE_TYPES:
                pcurve_real[name2d] += 1            # already stored analytically
            elif samples2d is None:
                pcurve_real["unsampleable"] += 1
            else:
                kind, residual = classify_curve(samples2d)
                pcurve_real[kind if residual < TOL_EXACT else "freeform"] += 1
                residuals.append(residual)

            try:
                stored3d = BRepAdaptor_Curve(edge).GetType()
            except Exception:
                stored3d = None
            name3d = CURVE_TYPE_NAMES.get(stored3d, "none" if stored3d is None else f"enum{stored3d}")
            curve3d_stored[name3d] += 1
            if stored3d in CANONICAL_CURVE_TYPES:
                curve3d_real[name3d] += 1
            else:
                curve3d, first, last = BRep_Tool.Curve(edge)
                samples3d = None
                if curve3d is not None:
                    def evaluate(t, c=curve3d):
                        p = c.Value(t)
                        return (p.X(), p.Y(), p.Z())
                    samples3d = sample_curve(evaluate, first, last)
                if samples3d is None:
                    curve3d_real["unsampleable"] += 1
                else:
                    kind, residual = classify_curve(samples3d)
                    curve3d_real[kind if residual < TOL_EXACT else "freeform"] += 1

    # Shared-edge disagreement: for every edge carried by two of the kept faces, how far apart are
    # the two faces' own images of it? This is the sliver gap, measured in model length units.
    edge_to_faces = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_to_faces)
    kept_hashes = {hash(f) for f in keep}
    shared, disagreements = 0, []
    for index in range(1, edge_to_faces.Size() + 1):
        edge = topods.Edge(edge_to_faces.FindKey(index))
        if BRep_Tool.Degenerated(edge):
            continue
        adjacent = [topods.Face(f) for f in edge_to_faces.FindFromIndex(index)]
        adjacent = [f for f in adjacent if hash(f) in kept_hashes]
        if len(adjacent) != 2:
            continue
        shared += 1
        first_polyline = _pcurve_polyline_3d(edge, adjacent[0])
        second_polyline = _pcurve_polyline_3d(edge, adjacent[1])
        if first_polyline is None or second_polyline is None:
            continue
        gap = max(_max_point_to_polyline(first_polyline, second_polyline),
                  _max_point_to_polyline(second_polyline, first_polyline))
        disagreements.append(gap)

    disagreements_array = np.array(disagreements) if disagreements else np.zeros(0)
    return {
        "n_faces_considered": len(keep),
        "pcurve_stored_types": dict(pcurve_stored),
        "pcurve_real_geometry": dict(pcurve_real),
        "curve3d_stored_types": dict(curve3d_stored),
        "curve3d_real_geometry": dict(curve3d_real),
        "n_shared_edges": shared,
        "shared_edge_gap_max": float(disagreements_array.max()) if len(disagreements_array) else 0.0,
        "shared_edge_gap_mean": float(disagreements_array.mean()) if len(disagreements_array) else 0.0,
        "shared_edge_gap_over_1e-6": int((disagreements_array > 1e-6).sum()),
        "shared_edge_gap_over_1e-4": int((disagreements_array > 1e-4).sum()),
    }


def load_shape(path):
    reader = STEPControl_Reader()
    reader.ReadFile(path)
    reader.TransferRoots()
    return reader.OneShape()


def faces_of(shape):
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        f = topods.Face(exp.Current())
        exp.Next()
        yield f


def analyze_faces(shape):
    """Per-face pass: stored type counts plus the re-classification of non-analytic faces."""
    stored, exact, loose = Counter(), Counter(), Counter()
    residuals = []
    n_faces = 0
    for face in faces_of(shape):
        n_faces += 1
        adaptor = BRepAdaptor_Surface(face)
        stype = adaptor.GetType()
        stored[STORED_TYPE_NAMES.get(stype, f"enum{stype}")] += 1
        if stype in ANALYTIC_TYPES:
            continue
        P, N = sample_surface(adaptor)
        if P is None:
            exact["unsampleable"] += 1
            loose["unsampleable"] += 1
            continue
        kind, res = classify_surface(P, N)
        exact[kind if res < TOL_EXACT else "freeform"] += 1
        loose[kind if res < TOL_LOOSE else "freeform"] += 1
        residuals.append(res)
    return {
        "n_faces": n_faces,
        "stored_types": dict(stored),
        "recognized_exact": dict(exact),
        "recognized_loose": dict(loose),
        "n_nonanalytic": int(sum(v for k, v in stored.items()
                                 if k not in ("plane", "cylinder", "cone", "sphere", "torus"))),
    }


def analyze_solids(shape):
    """Per-solid pass: how many solids become fully analytic if recognition is applied?"""
    solids, seen = [], set()
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    while exp.More():
        s = exp.Current()
        exp.Next()
        key = hash(s.TShape())          # deduplicate repeated placements of one definition
        if key in seen:
            continue
        seen.add(key)
        solids.append(s)

    already, rescued, blocked = 0, 0, 0
    recovered = Counter()
    for solid in solids:
        n_nonanalytic, n_recognized = 0, 0
        kinds = Counter()
        for face in faces_of(solid):
            adaptor = BRepAdaptor_Surface(face)
            if adaptor.GetType() in ANALYTIC_TYPES:
                continue
            n_nonanalytic += 1
            P, N = sample_surface(adaptor)
            if P is None:
                continue
            kind, res = classify_surface(P, N)
            if res < TOL_EXACT and kind in ("plane", "cylinder", "cone", "sphere"):
                n_recognized += 1
                kinds[kind] += 1
        if n_nonanalytic == 0:
            already += 1
        elif n_recognized == n_nonanalytic:
            rescued += 1
            recovered.update(kinds)
        else:
            blocked += 1
    return {
        "n_solids": len(solids),
        "already_analytic": already,
        "rescued_by_recognition": rescued,
        "still_freeform": blocked,
        "recovered_surface_kinds": dict(recovered),
    }


def main():
    ap = argparse.ArgumentParser(
        description="Classify STEP faces by their real geometry, not their stored surface type.")
    ap.add_argument("step", nargs="+", help="Input STEP file(s)")
    ap.add_argument("--per-solid", action="store_true",
                    help="Also report, per solid, how many would become fully analytic with "
                         "recognition applied (the exact-conversion coverage forecast)")
    ap.add_argument("--trim-curves", action="store_true",
                    help="Also classify the *trim curves* by their real geometry (how many "
                         "B-spline trim edges are exactly circles or lines) and measure how far "
                         "apart the two faces of a shared edge place their common boundary")
    ap.add_argument("--all-faces", action="store_true",
                    help="With --trim-curves, count trim edges of every face rather than only of "
                         "faces whose surface the exact converter can represent")
    ap.add_argument("--json", default=None, metavar="PATH", help="Write the results as JSON")
    args = ap.parse_args()

    out = {}
    for path in args.step:
        shape = load_shape(path)
        result = analyze_faces(shape)
        print(f"file: {path}")
        print(f"  faces total:                {result['n_faces']}")
        print(f"  stored surface types:       {result['stored_types']}")
        print(f"  non-analytic faces:         {result['n_nonanalytic']}")
        print(f"  re-classified by geometry:")
        print(f"    exact (rel < {TOL_EXACT:g}):    {result['recognized_exact']}")
        print(f"    loose (rel < {TOL_LOOSE:g}):    {result['recognized_loose']}")
        if args.per_solid:
            solid = analyze_solids(shape)
            result.update(solid)
            total = solid["n_solids"]
            print(f"  unique solids:              {total}")
            print(f"    already fully analytic:   {solid['already_analytic']}/{total}")
            print(f"    rescued by recognition:   {solid['rescued_by_recognition']}/{total}"
                  f"  -> forecast coverage "
                  f"{solid['already_analytic'] + solid['rescued_by_recognition']}/{total}")
            print(f"    still genuinely freeform: {solid['still_freeform']}/{total}")
            print(f"    surface kinds recovered:  {solid['recovered_surface_kinds']}")
        if args.trim_curves:
            trim = analyze_trim_curves(shape, convertible_surfaces_only=not args.all_faces)
            result["trim_curves"] = trim
            scope = "all faces" if args.all_faces else "convertible faces only"
            print(f"  trim curves ({scope}, {trim['n_faces_considered']} faces):")
            print(f"    pcurve stored types:      {trim['pcurve_stored_types']}")
            print(f"    pcurve real geometry:     {trim['pcurve_real_geometry']}")
            print(f"    3D curve stored types:    {trim['curve3d_stored_types']}")
            print(f"    3D curve real geometry:   {trim['curve3d_real_geometry']}")
            print(f"  shared edges: {trim['n_shared_edges']}"
                  f"  max pcurve disagreement {trim['shared_edge_gap_max']:.3e}"
                  f"  mean {trim['shared_edge_gap_mean']:.3e}  (model units)")
            print(f"    edges disagreeing by >1e-6: {trim['shared_edge_gap_over_1e-6']}"
                  f"  >1e-4: {trim['shared_edge_gap_over_1e-4']}")
        print()
        out[path] = result

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"Wrote {args.json}")


if __name__ == "__main__":
    sys.exit(main())
