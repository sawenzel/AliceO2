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

Residuals are relative to the patch size, so "exact" means machine precision. This is deliberate:
an "almost cylinder" that is really free-form must stay free-form, otherwise a converter acting on
this classification would silently change the geometry it exists to represent exactly.

Known limitation: there is no torus test yet, so torus patches stored as NURBS are reported as
free-form. Every "free-form" count is therefore an UPPER bound and every analytic count a LOWER
bound.

Note: OCCT's own BRepLib_CanonicalRecognition (7.7+) is not exposed in the pythonOCC v7.9.3 build
used here (ImportError from OCC.Core.BRepLib), which is why this is a standalone numeric recognizer.

Usage
-----
    python3 analyze_surface_geometry.py model.step [more.step ...]
    python3 analyze_surface_geometry.py --per-solid model.step      # coverage forecast
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
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_SOLID
from OCC.Core.TopoDS import topods
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import (
    GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus,
    GeomAbs_BSplineSurface, GeomAbs_BezierSurface,
    GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion,
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

# Relative residual below which a model is accepted as EXACT (machine precision).
TOL_EXACT = 1e-9
# A deliberately looser bound, reported alongside so the gap between "exact" and "close" is visible.
TOL_LOOSE = 1e-6

# Normals within this of parallel count as a plane.
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


def classify_surface(P, N):
    """Return (model_name, relative_residual) for the smallest model that fits the samples."""
    scale = np.linalg.norm(P.max(axis=0) - P.min(axis=0))
    if scale < 1e-12:
        return "degenerate", 0.0

    # --- plane (3 parameters): all unit normals parallel
    dev = np.abs(np.abs(N @ N[0]) - 1.0).max()
    if dev < TOL_PLANE_NORMALS:
        return "plane", float(dev)

    best = ("freeform", float(np.inf))

    # --- sphere (4): normal lines concurrent, P_i = C + r*N_i
    A = np.zeros((3 * len(P), 4))
    b = np.zeros(3 * len(P))
    for i in range(len(P)):
        A[3 * i:3 * i + 3, 0:3] = np.eye(3)
        A[3 * i:3 * i + 3, 3] = N[i]
        b[3 * i:3 * i + 3] = P[i]
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    centre, radius = sol[:3], sol[3]
    res = float(np.abs(np.linalg.norm(P - centre, axis=1) - abs(radius)).max() / scale)
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
            R = math.sqrt(r2)
            res = float(np.abs(np.hypot(x - cx, y - cy) - R).max() / scale)
            if res < best[1]:
                best = ("cylinder", res)

    # --- cone (6): N_i . (P_i - A) = 0 is linear in the apex A
    apex, *_ = np.linalg.lstsq(N, np.einsum('ij,ij->i', N, P), rcond=None)
    d = P - apex
    dn = np.linalg.norm(d, axis=1)
    ok = dn > 1e-12
    if ok.sum() > 10:
        u = d[ok] / dn[ok, None]
        res = float(np.abs(np.einsum('ij,ij->i', u, N[ok])).max())
        _, _, Vt2 = np.linalg.svd(u - u.mean(axis=0), full_matrices=False)
        ax2 = np.cross(Vt2[0], Vt2[1])
        n2 = np.linalg.norm(ax2)
        if n2 > 1e-12:                       # constant half-angle about the ruling axis
            res = max(res, float(np.abs(u @ (ax2 / n2)).std()))
        if res < best[1]:
            best = ("cone", res)

    return best


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
        print()
        out[path] = result

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"Wrote {args.json}")


if __name__ == "__main__":
    sys.exit(main())
