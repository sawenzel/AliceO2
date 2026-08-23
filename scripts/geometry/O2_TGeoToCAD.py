#!/usr/bin/env python3
"""
O2_TGeoToCAD.py -- TGeo -> STEP (AP214) with XCAF assembly structure.

The inverse of `O2_CADtoTGeo.py`: it reads a ROOT geometry file (an `o2-sim`
`o2sim_geometry.root`, ideal or aligned), walks the `TGeoVolume` DAG, builds an
OCCT solid for every volume whose shape it can map, and writes one STEP file with
the assembly tree preserved.

Why: the ALICE Run 3 geometry is written in C++ and only exists as TGeo.  A STEP
export makes it CAD-editable, and -- the reason this lives on the CAD -> TGeo
branch -- it makes the *original TGeo an exact oracle* for the reverse
conversion: the right answer is known by construction, so a round trip
TGeo -> STEP -> `O2_CADtoTGeo.py` can be scored without a tolerance argument.

The mapping
-----------
    TGeoVolume with no daughters      -> one XCAF simple shape (a definition)
    TGeoVolume with daughters         -> one XCAF assembly label, one component
                                         per TGeoNode referring to the daughter's
                                         definition, carrying the node's TGeoMatrix
    TGeoVolume with daughters AND its
      own (non-assembly) shape        -> the above, plus one extra component
                                         `<name>__body` holding the mother's own
                                         solid at the identity
    TGeoVolumeAssembly                -> a pure XCAF assembly, no solid

A logical volume is converted **once** and referenced from every node that places
it, so a shared volume costs one solid, matching TGeo's own DAG.

Mother solids are exported **uncarved**: TGeo's containment semantics (a daughter
implicitly subtracts itself from its mother) are deliberately *not* applied, so
every part in the STEP is the shape the TGeo author wrote and per-part capacity
is directly comparable against `TGeoShape::Capacity()`.  `--carve-mothers` turns
the subtraction on for a CAD-facing export; it changes the geometry and is not
what the round-trip acceptance scores.

Units: TGeo is cm, STEP is written in mm, so every length and translation is
scaled by 10.

Usage
-----
    O2_TGeoToCAD.py INPUT.root OUTPUT.step [options]
    O2_TGeoToCAD.py --self-test

    --report FILE        per-volume JSON report (default: <output>.report.json)
    --top VOLNAME        start from this volume instead of the TGeoManager top
    --max-depth N        stop descending below depth N
    --include-name PAT   only convert volumes whose name matches this glob
                         (their ancestors are still emitted as assemblies)
    --no-mother-bodies   omit the `__body` component of volumes with daughters
    --skip-top-body      omit only the top volume's own solid (the `cave` box)
    --carve-mothers      subtract placed daughters from each mother solid
    --dedup-world        expand the tree per occurrence and drop any placement of a
                         volume that coincides exactly with another placement of the
                         same volume (see "coincident placements" below)
    --no-verify          skip the per-definition BRepGProp capacity check
    --no-step            build and report, but do not write the STEP
    --quiet

Coincident placements
---------------------
TGeo lets the same volume be placed twice at the same world transform, and the
ALICE Run 3 beam pipe does exactly that 81 times (adjacent bellows convolutions
share a "plie" disc).  `O2_CADtoTGeo.py` refuses such a model outright, because a
geometry with two volumes in the same place has no defined transport through it.
The default DAG export reproduces the duplicates faithfully -- that is the honest
thing for a STEP that claims to be the geometry.  `--dedup-world` instead expands
the assembly tree per occurrence (leaf definitions stay shared, so there is still
one prototype per TGeoVolume) and drops every repeated (volume, world transform),
counting them in the report.

Report
------
`tgeo2step_report.json`-style: one record per logical volume with
{name, shapeClass, converted, reason, capacity_cm3, occVolume_cm3, relDev, ...}
plus a summary keyed by shape class.  A volume that could not be mapped carries a
machine-readable `reason`, the same decline-reason discipline the reverse
converter uses.
"""

import argparse
import fnmatch
import json
import math
import os
import sys
import time

# --------------------------------------------------------------------------
# OCCT
# --------------------------------------------------------------------------

from OCC.Core.gp import (
    gp_Pnt, gp_Dir, gp_Vec, gp_XYZ, gp_Ax1, gp_Ax2, gp_Trsf, gp_GTrsf, gp_Mat,
    gp_Elips, gp_Pln,
)
from OCC.Core.GC import GC_MakeArcOfCircle
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Shape
from OCC.Core.TopAbs import TopAbs_SOLID, TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeCone,
    BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus, BRepPrimAPI_MakeRevol,
    BRepPrimAPI_MakePrism, BRepPrimAPI_MakeHalfSpace,
)
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire, BRepBuilderAPI_Transform, BRepBuilderAPI_GTransform,
    BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
)
from OCC.Core.BRepFill import brepfill
from OCC.Core.TopoDS import topods
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Common
from OCC.Core.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GProp import GProp_GProps
from OCC.Core.TDocStd import TDocStd_Document
from OCC.Core.TDataStd import TDataStd_Name
from OCC.Core.TDF import TDF_LabelSequence, TDF_Label
from OCC.Core.XCAFDoc import XCAFDoc_DocumentTool
from OCC.Core.STEPCAFControl import STEPCAFControl_Writer
from OCC.Core.Interface import Interface_Static
from OCC.Core.IFSelect import IFSelect_RetDone

SCALE_TO_MM = 10.0          # TGeo cm -> STEP mm
BOOLEAN_VOLUME_TOL = 1e-4   # relative slack on the boolean volume invariant
EPS = 1e-12


class ShapeDeclined(Exception):
    """A TGeo shape this mapper does not (or could not) convert. The message is
    the machine-readable decline reason that lands in the report."""


# --------------------------------------------------------------------------
# small OCCT helpers
# --------------------------------------------------------------------------

def _moved(shape, trsf):
    return BRepBuilderAPI_Transform(shape, trsf, True).Shape()


def _rotz(deg):
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), math.radians(deg))
    return t


def _translate(dx, dy, dz):
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(float(dx), float(dy), float(dz)))
    return t


def _ax2(z0, phi1_deg):
    ph = math.radians(phi1_deg)
    return gp_Ax2(gp_Pnt(0.0, 0.0, float(z0)), gp_Dir(0, 0, 1),
                  gp_Dir(math.cos(ph), math.sin(ph), 0.0))


def solid_volume_mm3(shape):
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    return abs(props.Mass())


def _has_solid(shape):
    if shape is None or shape.IsNull():
        return False
    return TopExp_Explorer(shape, TopAbs_SOLID).More()


def _check(shape, what):
    if shape is None or shape.IsNull():
        raise ShapeDeclined(f"{what}: OCCT returned a null shape")
    if not _has_solid(shape):
        raise ShapeDeclined(f"{what}: OCCT result contains no solid")
    return shape


def _dedupe_ring(pts, tol=1e-9):
    """Drop consecutive duplicates in a closed point ring, wrap included."""
    out = []
    for p in pts:
        if out and abs(p[0] - out[-1][0]) < tol and abs(p[1] - out[-1][1]) < tol:
            continue
        out.append(p)
    while len(out) > 1 and abs(out[0][0] - out[-1][0]) < tol and abs(out[0][1] - out[-1][1]) < tol:
        out.pop()
    return out


def _revolve_profile(pts_rz, phi1_deg, dphi_deg, what):
    """Revolve a closed (r, z) profile in the x>=0 half of the XZ plane about +Z.

    This is the exact route for every solid of revolution: one operation, no
    booleans, and rmin > 0 comes out as a real inner face rather than a cut.
    """
    pts = _dedupe_ring([(float(r), float(z)) for (r, z) in pts_rz])
    if len(pts) < 3:
        raise ShapeDeclined(f"{what}: degenerate r-z profile ({len(pts)} distinct points)")
    if min(p[0] for p in pts) < -1e-9:
        raise ShapeDeclined(f"{what}: negative radius in profile")
    poly = BRepBuilderAPI_MakePolygon()
    for (r, z) in pts:
        poly.Add(gp_Pnt(r, 0.0, z))
    poly.Close()
    if not poly.IsDone():
        raise ShapeDeclined(f"{what}: could not build the r-z profile wire")
    mf = BRepBuilderAPI_MakeFace(poly.Wire())
    if not mf.IsDone():
        raise ShapeDeclined(f"{what}: r-z profile is not a valid planar face")
    rev = BRepPrimAPI_MakeRevol(mf.Face(), gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                math.radians(dphi_deg))
    rev.Build()
    if not rev.IsDone():
        raise ShapeDeclined(f"{what}: revolution of the r-z profile failed")
    sh = rev.Shape()
    if abs(phi1_deg) > 1e-12:
        sh = _moved(sh, _rotz(phi1_deg))
    return _check(sh, what)


def _revolve_edges(elements, phi1_deg, dphi_deg, what):
    """Revolve a closed (r, z) profile made of line and arc elements about +Z.

    Elements are ("line", p1, p2) or ("arc", p1, pmid, p2), each point an (r, z)
    pair in the x >= 0 half of the XZ plane.  This is the sphere route: OCCT's
    BRepPrimAPI_MakeSphere cuts theta with *planes* (a spherical zone), while TGeo
    cuts it with *cones* through the centre (a spherical cone), so the primitive
    cannot be used for a theta-sectioned sphere at all.
    """
    def _p(rz):
        return gp_Pnt(float(rz[0]), 0.0, float(rz[1]))

    mw = BRepBuilderAPI_MakeWire()
    nedges = 0
    for e in elements:
        if e[0] == "line":
            p1, p2 = e[1], e[2]
            if math.hypot(p1[0] - p2[0], p1[1] - p2[1]) < 1e-9:
                continue
            mw.Add(BRepBuilderAPI_MakeEdge(_p(p1), _p(p2)).Edge())
        else:
            p1, pm, p2 = e[1], e[2], e[3]
            arc = GC_MakeArcOfCircle(_p(p1), _p(pm), _p(p2))
            if not arc.IsDone():
                raise ShapeDeclined(f"{what}: could not build a profile arc")
            mw.Add(BRepBuilderAPI_MakeEdge(arc.Value()).Edge())
        nedges += 1
    if nedges < 2 or not mw.IsDone():
        raise ShapeDeclined(f"{what}: could not close the r-z profile wire")
    mf = BRepBuilderAPI_MakeFace(mw.Wire())
    if not mf.IsDone():
        raise ShapeDeclined(f"{what}: r-z profile is not a valid planar face")
    rev = BRepPrimAPI_MakeRevol(mf.Face(), gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                math.radians(dphi_deg))
    rev.Build()
    if not rev.IsDone():
        raise ShapeDeclined(f"{what}: revolution of the r-z profile failed")
    sh = rev.Shape()
    if abs(phi1_deg) > 1e-12:
        sh = _moved(sh, _rotz(phi1_deg))
    return _check(sh, what)


def _signed_volume(shape):
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    return props.Mass()


def _polygon_wire(pts, what):
    poly = BRepBuilderAPI_MakePolygon()
    for (x, y, z) in pts:
        poly.Add(gp_Pnt(float(x), float(y), float(z)))
    poly.Close()
    if not poly.IsDone():
        raise ShapeDeclined(f"{what}: could not build a polygon wire")
    return poly.Wire()


def _dedupe_ring3(pts, tol=1e-9):
    out = []
    for p in pts:
        if out and max(abs(p[i] - out[-1][i]) for i in range(3)) < tol:
            continue
        out.append(p)
    while len(out) > 1 and max(abs(out[0][i] - out[-1][i]) for i in range(3)) < tol:
        out.pop()
    return out


def _quad_face(b0, b1, t1, t0, what, tol=1e-7):
    """One lateral patch of a prism: a planar face when the four corners are
    coplanar, a ruled face when they are not (a twisted TGeoArb8 side).

    Emitting a plane as a plane matters: the reverse converter's recogniser keys
    on the analytic surface type, and a lofted B-spline that happens to be flat is
    reported as a free-form face.
    """
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    pts = _dedupe_ring3([b0, b1, t1, t0])
    if len(pts) < 3:
        return None
    if len(pts) == 3:
        return BRepBuilderAPI_MakeFace(_polygon_wire(pts, what)).Face()
    n = cross(sub(b1, b0), sub(t0, b0))
    nn = math.sqrt(sum(c * c for c in n))
    scale = max(math.sqrt(sum(c * c for c in sub(b1, b0))),
                math.sqrt(sum(c * c for c in sub(t0, b0))), 1e-30)
    d = sub(t1, b0)
    off = abs(sum(n[i] * d[i] for i in range(3))) / nn if nn > 0 else 0.0
    if nn > 1e-24 and off <= tol * scale:
        mf = BRepBuilderAPI_MakeFace(_polygon_wire(pts, what))
        if mf.IsDone():
            return mf.Face()
    e1 = BRepBuilderAPI_MakeEdge(gp_Pnt(*b0), gp_Pnt(*b1)).Edge()
    e2 = BRepBuilderAPI_MakeEdge(gp_Pnt(*t0), gp_Pnt(*t1)).Edge()
    return brepfill.Face(e1, e2)


def _prism_from_rings(outer, inner=None, what="prism"):
    """Build a solid from a stack of closed sections by sewing explicit faces.

    `outer` (and the optional `inner`, which makes the caps annular) is a list of
    sections, each a list of (x, y, z) with the same vertex count and order.
    """
    outer = [_dedupe_ring3(r) for r in outer]
    if len(outer) < 2:
        raise ShapeDeclined(f"{what}: fewer than two sections")
    nv = len(outer[0])
    if nv < 3 or any(len(r) != nv for r in outer):
        raise ShapeDeclined(
            f"{what}: sections have {sorted(set(len(r) for r in outer))} distinct "
            "vertices; a prism needs the same count in every section")
    rings = [outer]
    if inner is not None:
        inner = [_dedupe_ring3(r) for r in inner]
        if any(len(r) != nv for r in inner):
            raise ShapeDeclined(f"{what}: inner sections do not match the outer count")
        rings.append(inner)

    faces = []
    for ring in rings:
        for k in range(len(ring) - 1):
            lo, hi = ring[k], ring[k + 1]
            for i in range(nv):
                j = (i + 1) % nv
                f = _quad_face(lo[i], lo[j], hi[j], hi[i], what)
                if f is not None:
                    faces.append(f)
    # caps, annular when there is an inner stack
    for idx in (0, -1):
        mf = BRepBuilderAPI_MakeFace(_polygon_wire(outer[idx], what))
        if inner is not None:
            mf.Add(topods.Wire(_polygon_wire(inner[idx], what).Reversed()))
        if not mf.IsDone():
            raise ShapeDeclined(f"{what}: could not build a cap face")
        faces.append(mf.Face())

    ext = max(abs(c) for r in outer for p in r for c in p) or 1.0
    sew = BRepBuilderAPI_Sewing(1e-7 * ext)
    for f in faces:
        sew.Add(f)
    sew.Perform()
    shell = sew.SewedShape()
    if shell is None or shell.IsNull():
        raise ShapeDeclined(f"{what}: sewing produced nothing")
    try:
        ms = BRepBuilderAPI_MakeSolid(topods.Shell(shell))
        ms.Build()
        solid = ms.Solid()
    except Exception as e:
        raise ShapeDeclined(f"{what}: faces did not sew into a closed shell ({e})")
    if _signed_volume(solid) < 0:
        solid = topods.Solid(solid.Reversed())
    return _check(solid, what)


def _unify(shape):
    """Merge co-planar / co-cylindrical neighbouring faces.

    A solid built by a chain of fuses carries a seam face at every merge, and
    OCCT's BOP can choke on those without saying so -- see `_boolean` below.
    """
    u = ShapeUpgrade_UnifySameDomain(shape)
    u.Build()
    return u.Shape()


def _run_boolean(op, a, b):
    algo = op(a, b)
    algo.Build()
    if not algo.IsDone():
        return None
    return algo.Shape()


def _boolean(op, a, b, what, lower=True):
    """A boolean with a volume invariant, because OCCT can fail silently.

    On one MFT part (`MF21`) `BRepAlgoAPI_Fuse` returned `IsDone() == True` and a
    valid solid that was **just the second operand** -- the 4.38 cm^3 first operand
    had vanished, and the enclosing volume came out 25 % wrong with nothing
    reporting an error. Every boolean is therefore priced against what its result
    must satisfy:

        fuse    max(vA, vB) <= v <= vA + vB
        cut     vA - vB     <= v <= vA
        common  0           <= v <= min(vA, vB)

    `lower=False` drops the lower bound, for the one case where the tool has no
    meaningful volume: a `BRepPrimAPI_MakeHalfSpace` solid is unbounded, so vB is
    not a number the bound can be built from.

    A violation is first treated as the seam-face problem it turned out to be
    (retry with the operands unified), and only then declined -- never shipped.

    The band is deliberately loose (`BOOLEAN_VOLUME_TOL`, 1e-4 relative): this is a
    net for a boolean that silently lost a whole operand, not an accuracy test.
    `BRepGProp` integrates a union of near-disjoint solids to about 1e-6 relative,
    and a tighter band fires on that instead.
    """
    try:
        va, vb = solid_volume_mm3(a), solid_volume_mm3(b)
    except Exception:
        va = vb = None

    def bounds(v):
        if va is None:
            return True, ""
        tol = BOOLEAN_VOLUME_TOL * max(va, vb, 1.0)
        if op is BRepAlgoAPI_Fuse:
            lo, hi = max(va, vb) - tol, va + vb + tol
        elif op is BRepAlgoAPI_Cut:
            lo, hi = va - vb - tol, va + tol
        else:
            lo, hi = -tol, min(va, vb) + tol
        if not lower:
            lo = -tol
        return lo <= v <= hi, f"{v:.6g} outside [{lo:.6g}, {hi:.6g}] mm^3"

    sh = _run_boolean(op, a, b)
    if sh is None:
        raise ShapeDeclined(f"{what}: OCCT boolean did not complete")
    _check(sh, what)
    ok, msg = bounds(solid_volume_mm3(sh))
    if ok:
        return sh
    retry = _run_boolean(op, _unify(a), _unify(b))
    if retry is not None and _has_solid(retry):
        ok2, msg2 = bounds(solid_volume_mm3(retry))
        if ok2:
            return retry
        msg = f"{msg}; after unifying the operands {msg2}"
    raise ShapeDeclined(
        f"{what}: OCCT's boolean returned a volume the operands cannot give "
        f"({msg}); the operation failed silently")


# --------------------------------------------------------------------------
# TGeoMatrix -> OCCT transform
# --------------------------------------------------------------------------

def tgeo_matrix_components(m):
    """(3x3 row-major matrix including any TGeoScale, translation in mm)."""
    if m is None:
        return [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0., 0., 0.]
    r = m.GetRotationMatrix()
    s = m.GetScale()
    t = m.GetTranslation()
    rot = [[float(r[3 * i + j]) for j in range(3)] for i in range(3)]
    sc = [float(s[j]) for j in range(3)]
    mat = [[rot[i][j] * sc[j] for j in range(3)] for i in range(3)]
    tr = [float(t[i]) * SCALE_TO_MM for i in range(3)]
    return mat, tr


def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def tgeo_matrix_to_trsf(m):
    """A rigid gp_Trsf, or None if the matrix is a reflection / non-uniform scale.

    STEP component placements are rigid by construction, so a non-rigid TGeoMatrix
    cannot be a placement -- it has to be baked into a copy of the solid instead.
    """
    mat, tr = tgeo_matrix_components(m)
    if abs(_det3(mat) - 1.0) > 1e-9:
        return None
    t = gp_Trsf()
    try:
        t.SetValues(mat[0][0], mat[0][1], mat[0][2], tr[0],
                    mat[1][0], mat[1][1], mat[1][2], tr[1],
                    mat[2][0], mat[2][1], mat[2][2], tr[2])
    except Exception:
        return None
    return t


def tgeo_matrix_to_gtrsf(m):
    mat, tr = tgeo_matrix_components(m)
    g = gp_GTrsf()
    g.SetVectorialPart(gp_Mat(mat[0][0], mat[0][1], mat[0][2],
                              mat[1][0], mat[1][1], mat[1][2],
                              mat[2][0], mat[2][1], mat[2][2]))
    g.SetTranslationPart(gp_XYZ(tr[0], tr[1], tr[2]))
    return g


def apply_tgeo_matrix(shape, m, what):
    """Move `shape` by a TGeoMatrix, using a general transform when it reflects."""
    t = tgeo_matrix_to_trsf(m)
    if t is not None:
        return _moved(shape, t)
    g = tgeo_matrix_to_gtrsf(m)
    algo = BRepBuilderAPI_GTransform(shape, g, True)
    if not algo.IsDone():
        raise ShapeDeclined(f"{what}: could not apply a reflecting/scaling matrix")
    return _check(algo.Shape(), what)


# --------------------------------------------------------------------------
# shape converters -- all output mm, centred as TGeo centres them
# --------------------------------------------------------------------------

def _phi_span(phi1, phi2):
    d = float(phi2) - float(phi1)
    while d <= 0:
        d += 360.0
    return float(phi1), min(d, 360.0)


def conv_box(sh, s):
    dx, dy, dz = sh.GetDX() * s, sh.GetDY() * s, sh.GetDZ() * s
    ox, oy, oz = (sh.GetOrigin()[0] * s, sh.GetOrigin()[1] * s, sh.GetOrigin()[2] * s)
    if min(dx, dy, dz) <= 0:
        raise ShapeDeclined("TGeoBBox: a half-length is zero or negative")
    box = BRepPrimAPI_MakeBox(2 * dx, 2 * dy, 2 * dz).Shape()
    return _moved(box, _translate(ox - dx, oy - dy, oz - dz))


def _tube_like(rmin, rmax, dz, phi1, dphi, what):
    if rmax <= 0 or dz <= 0:
        raise ShapeDeclined(f"{what}: rmax or dz is zero")
    if rmin >= rmax:
        raise ShapeDeclined(f"{what}: rmin >= rmax")
    if rmin <= EPS:
        cyl = BRepPrimAPI_MakeCylinder(_ax2(-dz, phi1), rmax, 2 * dz, math.radians(dphi))
        cyl.Build()
        if not cyl.IsDone():
            raise ShapeDeclined(f"{what}: BRepPrimAPI_MakeCylinder failed")
        return _check(cyl.Shape(), what)
    return _revolve_profile([(rmin, -dz), (rmax, -dz), (rmax, dz), (rmin, dz)],
                            phi1, dphi, what)


def conv_tube(sh, s):
    return _tube_like(sh.GetRmin() * s, sh.GetRmax() * s, sh.GetDz() * s,
                      0.0, 360.0, "TGeoTube")


def conv_tubeseg(sh, s):
    phi1, dphi = _phi_span(sh.GetPhi1(), sh.GetPhi2())
    return _tube_like(sh.GetRmin() * s, sh.GetRmax() * s, sh.GetDz() * s,
                      phi1, dphi, "TGeoTubeSeg")


def _cone_like(rmin1, rmax1, rmin2, rmax2, dz, phi1, dphi, what):
    if dz <= 0:
        raise ShapeDeclined(f"{what}: dz is zero")
    if max(rmax1, rmax2) <= 0:
        raise ShapeDeclined(f"{what}: both outer radii are zero")
    if rmin1 <= EPS and rmin2 <= EPS:
        cone = BRepPrimAPI_MakeCone(_ax2(-dz, phi1), rmax1, rmax2, 2 * dz, math.radians(dphi))
        cone.Build()
        if not cone.IsDone():
            raise ShapeDeclined(f"{what}: BRepPrimAPI_MakeCone failed")
        return _check(cone.Shape(), what)
    return _revolve_profile([(rmin1, -dz), (rmax1, -dz), (rmax2, dz), (rmin2, dz)],
                            phi1, dphi, what)


def conv_cone(sh, s):
    return _cone_like(sh.GetRmin1() * s, sh.GetRmax1() * s,
                      sh.GetRmin2() * s, sh.GetRmax2() * s, sh.GetDz() * s,
                      0.0, 360.0, "TGeoCone")


def conv_coneseg(sh, s):
    phi1, dphi = _phi_span(sh.GetPhi1(), sh.GetPhi2())
    return _cone_like(sh.GetRmin1() * s, sh.GetRmax1() * s,
                      sh.GetRmin2() * s, sh.GetRmax2() * s, sh.GetDz() * s,
                      phi1, dphi, "TGeoConeSeg")


def conv_pcon(sh, s):
    nz = int(sh.GetNz())
    if nz < 2:
        raise ShapeDeclined("TGeoPcon: fewer than two z planes")
    z = [sh.GetZ(i) * s for i in range(nz)]
    rmin = [sh.GetRmin(i) * s for i in range(nz)]
    rmax = [sh.GetRmax(i) * s for i in range(nz)]
    phi1, dphi = float(sh.GetPhi1()), float(sh.GetDphi())
    outer = [(rmax[i], z[i]) for i in range(nz)]
    if all(r <= EPS for r in rmin):
        inner = [(0.0, z[nz - 1]), (0.0, z[0])]
    else:
        inner = [(rmin[i], z[i]) for i in range(nz - 1, -1, -1)]
    return _revolve_profile(outer + inner, phi1, dphi, "TGeoPcon")


def _pgon_ring(r_apothem, z, phi1_deg, dphi_deg, nedges, full):
    """The polygon at one z plane. TGeo's rmin/rmax are inscribed-circle radii."""
    dseg = math.radians(dphi_deg) / nedges
    R = r_apothem / math.cos(dseg / 2.0)
    n = nedges if full else nedges + 1
    out = []
    for k in range(n):
        a = math.radians(phi1_deg) + k * dseg
        out.append((R * math.cos(a), R * math.sin(a), z))
    return out


def conv_pgon(sh, s):
    nz = int(sh.GetNz())
    nedges = int(sh.GetNedges())
    if nz < 2 or nedges < 1:
        raise ShapeDeclined("TGeoPgon: fewer than two z planes or no edges")
    phi1, dphi = float(sh.GetPhi1()), float(sh.GetDphi())
    full = abs(dphi - 360.0) < 1e-9
    z = [sh.GetZ(i) * s for i in range(nz)]
    rmin = [sh.GetRmin(i) * s for i in range(nz)]
    rmax = [sh.GetRmax(i) * s for i in range(nz)]
    hollow = any(r > EPS for r in rmin)
    rings = []
    for i in range(nz):
        outer = _pgon_ring(rmax[i], z[i], phi1, dphi, nedges, full)
        if hollow:
            inner = _pgon_ring(max(rmin[i], EPS), z[i], phi1, dphi, nedges, full)
            ring = outer + list(reversed(inner))
        elif full:
            ring = outer
        else:
            ring = outer + [(0.0, 0.0, z[i])]
        rings.append(ring)
    if hollow and full:
        # A full hollow polyhedra is an annulus in cross-section: two disjoint
        # rings, which no single lofted wire can express. Loft the outer and the
        # inner prism separately and subtract.
        outer_rings = [_pgon_ring(rmax[i], z[i], phi1, dphi, nedges, True) for i in range(nz)]
        inner_rings = [_pgon_ring(max(rmin[i], EPS), z[i], phi1, dphi, nedges, True) for i in range(nz)]
        return _prism_from_rings(outer_rings, inner_rings, what="TGeoPgon")
    return _prism_from_rings(rings, what="TGeoPgon")


def conv_sphere(sh, s):
    rmin, rmax = sh.GetRmin() * s, sh.GetRmax() * s
    th1, th2 = float(sh.GetTheta1()), float(sh.GetTheta2())
    phi1, dphi = _phi_span(sh.GetPhi1(), sh.GetPhi2())
    if rmax <= 0:
        raise ShapeDeclined("TGeoSphere: rmax is zero")
    if th2 <= th1:
        raise ShapeDeclined("TGeoSphere: theta2 <= theta1")

    def rz(r, theta_deg):
        a = math.radians(theta_deg)
        return (r * math.sin(a), r * math.cos(a))

    thm = 0.5 * (th1 + th2)
    p1, pm, p2 = rz(rmax, th1), rz(rmax, thm), rz(rmax, th2)
    elems = [("arc", p1, pm, p2)]
    if rmin > EPS:
        q1, qm, q2 = rz(rmin, th1), rz(rmin, thm), rz(rmin, th2)
        elems += [("line", p2, q2), ("arc", q2, qm, q1), ("line", q1, p1)]
    elif p1[0] < 1e-9 and p2[0] < 1e-9:
        elems += [("line", p2, p1)]                 # both poles: close on the axis
    else:
        elems += [("line", p2, (0.0, 0.0)), ("line", (0.0, 0.0), p1)]
    return _revolve_edges(elems, phi1, dphi, "TGeoSphere")


def conv_torus(sh, s):
    R = sh.GetR() * s
    rmin, rmax = sh.GetRmin() * s, sh.GetRmax() * s
    phi1, dphi = float(sh.GetPhi1()), float(sh.GetDphi())
    if rmax <= 0 or R <= 0:
        raise ShapeDeclined("TGeoTorus: R or Rmax is zero")

    def mk(r):
        m = BRepPrimAPI_MakeTorus(_ax2(0.0, phi1), R, r, math.radians(dphi))
        m.Build()
        if not m.IsDone():
            raise ShapeDeclined("TGeoTorus: BRepPrimAPI_MakeTorus failed")
        return _check(m.Shape(), "TGeoTorus")

    outer = mk(rmax)
    if rmin > EPS:
        return _boolean(BRepAlgoAPI_Cut, outer, mk(rmin), "TGeoTorus(hollow)")
    return outer


def conv_eltu(sh, s):
    a, b, dz = sh.GetA() * s, sh.GetB() * s, sh.GetDz() * s
    if a <= 0 or b <= 0 or dz <= 0:
        raise ShapeDeclined("TGeoEltu: a semi-axis or dz is zero")
    if a >= b:
        ax = gp_Ax2(gp_Pnt(0, 0, -dz), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
        maj, mnr = a, b
    else:
        ax = gp_Ax2(gp_Pnt(0, 0, -dz), gp_Dir(0, 0, 1), gp_Dir(0, 1, 0))
        maj, mnr = b, a
    edge = BRepBuilderAPI_MakeEdge(gp_Elips(ax, maj, mnr)).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    face = BRepBuilderAPI_MakeFace(wire).Face()
    pr = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, 2 * dz))
    pr.Build()
    if not pr.IsDone():
        raise ShapeDeclined("TGeoEltu: prism failed")
    return _check(pr.Shape(), "TGeoEltu")


def conv_trd1(sh, s):
    dx1, dx2 = sh.GetDx1() * s, sh.GetDx2() * s
    dy, dz = sh.GetDy() * s, sh.GetDz() * s
    return _prism_from_rings([
        [(-dx1, -dy, -dz), (dx1, -dy, -dz), (dx1, dy, -dz), (-dx1, dy, -dz)],
        [(-dx2, -dy, dz), (dx2, -dy, dz), (dx2, dy, dz), (-dx2, dy, dz)],
    ], what="TGeoTrd1")


def conv_trd2(sh, s):
    dx1, dx2 = sh.GetDx1() * s, sh.GetDx2() * s
    dy1, dy2 = sh.GetDy1() * s, sh.GetDy2() * s
    dz = sh.GetDz() * s
    return _prism_from_rings([
        [(-dx1, -dy1, -dz), (dx1, -dy1, -dz), (dx1, dy1, -dz), (-dx1, dy1, -dz)],
        [(-dx2, -dy2, dz), (dx2, -dy2, dz), (dx2, dy2, dz), (-dx2, dy2, dz)],
    ], what="TGeoTrd2")


def conv_arb8(sh, s):
    """TGeoArb8 and its subclasses (Trap): eight vertices, ruled lateral faces."""
    v = sh.GetVertices()
    dz = sh.GetDz() * s
    bot = [(v[2 * i] * s, v[2 * i + 1] * s, -dz) for i in range(4)]
    top = [(v[8 + 2 * i] * s, v[8 + 2 * i + 1] * s, dz) for i in range(4)]
    return _prism_from_rings([bot, top], what=sh.ClassName())


def conv_xtru(sh, s):
    nv, nz = int(sh.GetNvert()), int(sh.GetNz())
    if nv < 3 or nz < 2:
        raise ShapeDeclined("TGeoXtru: fewer than 3 vertices or 2 sections")
    x = [sh.GetX(i) for i in range(nv)]
    y = [sh.GetY(i) for i in range(nv)]
    rings = []
    for k in range(nz):
        z = sh.GetZ(k) * s
        x0, y0, sc = sh.GetXOffset(k) * s, sh.GetYOffset(k) * s, sh.GetScale(k)
        rings.append([(x0 + sc * x[i] * s, y0 + sc * y[i] * s, z) for i in range(nv)])
    return _prism_from_rings(rings, what="TGeoXtru")


def conv_ctub(sh, s):
    """A cut tube.

    TGeo's cut planes *replace* the +-dz end faces rather than being intersected
    with them: the solid is the tube's (r, phi) region between the two planes, and
    its z extent follows the planes.  So the tube has to be built long enough to
    reach past both planes before they cut it -- clipping at +-dz first would
    remove real material (measured: 5.283 cm^3 instead of 6.283 cm^3 on the
    self-test's slanted case).
    """
    rmin, rmax = sh.GetRmin() * s, sh.GetRmax() * s
    dz = sh.GetDz() * s
    phi1, dphi = _phi_span(sh.GetPhi1(), sh.GetPhi2())
    planes = []
    ext = 0.0
    for nvec, z0 in ((sh.GetNlow(), -dz), (sh.GetNhigh(), dz)):
        n = (float(nvec[0]), float(nvec[1]), float(nvec[2]))
        norm = math.sqrt(sum(c * c for c in n))
        if norm < EPS:
            raise ShapeDeclined("TGeoCtub: a cut normal is null")
        n = tuple(c / norm for c in n)
        if abs(n[2]) < 1e-9:
            raise ShapeDeclined("TGeoCtub: a cut plane is parallel to the axis")
        planes.append((n, z0))
        ext = max(ext, rmax * math.hypot(n[0], n[1]) / abs(n[2]))
    ext = ext * 1.5 + 1e-3 * max(rmax, dz)
    base = _tube_like(rmin, rmax, dz + ext, phi1, dphi, "TGeoCtub(base tube)")
    big = 4.0 * max(rmax, dz + ext) + 10.0
    for (n, z0) in planes:
        pl = BRepBuilderAPI_MakeFace(
            gp_Pln(gp_Pnt(0.0, 0.0, z0), gp_Dir(*n)), -big, big, -big, big)
        if not pl.IsDone():
            raise ShapeDeclined("TGeoCtub: could not build a cut plane")
        ref = gp_Pnt(n[0] * big, n[1] * big, z0 + n[2] * big)   # outside the solid
        hs = BRepPrimAPI_MakeHalfSpace(pl.Face(), ref)
        hs.Build()
        if not hs.IsDone():
            raise ShapeDeclined("TGeoCtub: half-space construction failed")
        base = _boolean(BRepAlgoAPI_Cut, base, hs.Solid(), "TGeoCtub", lower=False)
    return base


_BOOL_OPS = {
    "TGeoUnion": (BRepAlgoAPI_Fuse, "union"),
    "TGeoSubtraction": (BRepAlgoAPI_Cut, "subtraction"),
    "TGeoIntersection": (BRepAlgoAPI_Common, "intersection"),
}


def conv_composite(sh, s, depth=0):
    if depth > 32:
        raise ShapeDeclined("TGeoCompositeShape: boolean tree deeper than 32")
    bn = sh.GetBoolNode()
    if bn is None:
        raise ShapeDeclined("TGeoCompositeShape: no boolean node")
    op = _BOOL_OPS.get(bn.ClassName())
    if op is None:
        raise ShapeDeclined(f"TGeoCompositeShape: unknown boolean node {bn.ClassName()}")
    algo, opname = op
    left = shape_to_occ(bn.GetLeftShape(), s, depth + 1)
    right = shape_to_occ(bn.GetRightShape(), s, depth + 1)
    left = apply_tgeo_matrix(left, bn.GetLeftMatrix(), "composite left operand")
    right = apply_tgeo_matrix(right, bn.GetRightMatrix(), "composite right operand")
    return _boolean(algo, left, right, f"TGeoCompositeShape({opname})")


def conv_scaled(sh, s, depth=0):
    inner = shape_to_occ(sh.GetShape(), s, depth + 1)
    sc = sh.GetScale().GetScale()
    g = gp_GTrsf()
    g.SetVectorialPart(gp_Mat(sc[0], 0, 0, 0, sc[1], 0, 0, 0, sc[2]))
    algo = BRepBuilderAPI_GTransform(inner, g, True)
    if not algo.IsDone():
        raise ShapeDeclined("TGeoScaledShape: could not apply the scale")
    return _check(algo.Shape(), "TGeoScaledShape")


_DISPATCH = {
    "TGeoBBox": conv_box,
    "TGeoTube": conv_tube,
    "TGeoTubeSeg": conv_tubeseg,
    "TGeoCtub": conv_ctub,
    "TGeoCone": conv_cone,
    "TGeoConeSeg": conv_coneseg,
    "TGeoPcon": conv_pcon,
    "TGeoPgon": conv_pgon,
    "TGeoSphere": conv_sphere,
    "TGeoTorus": conv_torus,
    "TGeoEltu": conv_eltu,
    "TGeoTrd1": conv_trd1,
    "TGeoTrd2": conv_trd2,
    "TGeoArb8": conv_arb8,
    "TGeoTrap": conv_arb8,
    "TGeoXtru": conv_xtru,
}

# Shapes we know about and deliberately do not map, with the reason.
_KNOWN_DECLINES = {
    "TGeoHalfSpace": "unbounded solid: a half-space has no B-rep body of its own",
    "TGeoGtra": "twisted trapezoid: the lateral twist is not a ruled loft of the "
                "eight Arb8 vertices",
    "TGeoParaboloid": "quadric of revolution not mapped (no OCCT primitive; would "
                      "need a revolved parabola profile)",
    "TGeoHype": "hyperboloid of revolution not mapped",
    "TGeoPara": "parallelepiped not mapped",
    "TGeoTessellated": "already a mesh: STEP would carry facets, not a B-rep solid",
    "TGeoShapeAssembly": "assembly shape: emitted as a pure XCAF assembly, no solid",
}


def shape_to_occ(sh, s=SCALE_TO_MM, depth=0):
    """TGeoShape -> TopoDS_Shape in mm. Raises ShapeDeclined with a reason."""
    if sh is None:
        raise ShapeDeclined("volume has no shape")
    cls = sh.ClassName()
    if cls == "TGeoCompositeShape":
        return conv_composite(sh, s, depth)
    if cls == "TGeoScaledShape":
        return conv_scaled(sh, s, depth)
    fn = _DISPATCH.get(cls)
    if fn is None:
        raise ShapeDeclined(_KNOWN_DECLINES.get(cls, f"shape class {cls} is not mapped"))
    return fn(sh, s)


# --------------------------------------------------------------------------
# the XCAF builder
# --------------------------------------------------------------------------

class TGeoToStep:
    def __init__(self, opts):
        self.opts = opts
        self.doc = TDocStd_Document("O2_TGeoToCAD")
        self.shape_tool = XCAFDoc_DocumentTool.ShapeTool(self.doc.Main())
        self.definitions = {}      # volume name -> (label, occ solid or None)
        self.records = {}          # volume name -> report record
        self.ncomponents = 0
        self.nbaked = 0
        self.placed_world = set()      # (definition name, world key), --dedup-world only
        self.ndropped = 0
        self.dropped_examples = []
        self.reflected_nodes = []
        self.t0 = time.time()

    # ------------------------------------------------------------------

    def log(self, *a):
        if not self.opts.quiet:
            print(*a, file=sys.stderr, flush=True)

    def _record(self, vol, **kw):
        rec = self.records.setdefault(vol.GetName(), {
            "name": str(vol.GetName()),
            "shapeClass": vol.GetShape().ClassName() if vol.GetShape() else None,
            "ndaughters": int(vol.GetNdaughters()),
            "isAssembly": bool(vol.IsAssembly()),
            "converted": False,
            "reason": None,
            "capacity_cm3": None,
            "occVolume_cm3": None,
            "relDev": None,
        })
        rec.update(kw)
        return rec

    # ------------------------------------------------------------------

    def _solid_for(self, vol):
        """The OCCT solid of a volume's own shape, or (None, reason)."""
        sh = vol.GetShape()
        if sh is None or vol.IsAssembly() or sh.ClassName() == "TGeoShapeAssembly":
            return None, "pure assembly: no solid of its own, by design"
        try:
            occ = shape_to_occ(sh, SCALE_TO_MM)
        except ShapeDeclined as e:
            return None, str(e)
        except Exception as e:                                   # OCCT can throw
            return None, f"{sh.ClassName()}: OCCT raised {type(e).__name__}: {e}"
        return occ, None

    def _verify(self, vol, occ, rec):
        sh = vol.GetShape()
        try:
            cap = float(sh.Capacity())
        except Exception:
            cap = None
        rec["capacity_cm3"] = cap
        if not self.opts.verify:
            return
        try:
            v_cm3 = solid_volume_mm3(occ) / 1000.0
        except Exception as e:
            rec["occVolume_cm3"] = None
            rec["verifyError"] = str(e)
            return
        rec["occVolume_cm3"] = v_cm3
        if cap and abs(cap) > 0:
            rec["relDev"] = abs(v_cm3 - cap) / abs(cap)

    # ------------------------------------------------------------------

    def build(self, vol, depth=0):
        """Return the XCAF label for `vol`, building it (once) if needed."""
        name = str(vol.GetName())
        if name in self.definitions:
            return self.definitions[name][0]

        nd = int(vol.GetNdaughters())
        descend = nd > 0 and (self.opts.max_depth is None or depth < self.opts.max_depth)

        wanted = (self.opts.include_name is None
                  or fnmatch.fnmatch(name, self.opts.include_name))

        occ, reason = (None, None)
        if wanted:
            occ, reason = self._solid_for(vol)
        else:
            reason = "excluded by --include-name"
        rec = self._record(vol, converted=occ is not None, reason=reason)
        if occ is not None:
            self._verify(vol, occ, rec)

        if not descend:
            if occ is None:
                self.definitions[name] = (None, None)
                return None
            lab = self.shape_tool.AddShape(occ, False)
            TDataStd_Name.Set(lab, name)
            self.definitions[name] = (lab, occ)
            return lab

        # a volume with daughters becomes an assembly
        asm = self.shape_tool.NewShape()
        TDataStd_Name.Set(asm, name)
        self.definitions[name] = (asm, occ)
        ncomp0 = self.ncomponents

        emit_body = (occ is not None and self.opts.mother_bodies
                     and not (depth == 0 and self.opts.skip_top_body))
        placed_children = []
        for i in range(nd):
            node = vol.GetNode(i)
            child = node.GetVolume()
            clab = self.build(child, depth + 1)
            if clab is None:
                continue
            m = node.GetMatrix()
            t = tgeo_matrix_to_trsf(m)
            if t is None:
                # A reflecting or scaling placement cannot be a STEP component
                # location, so bake it into a private copy of the child solid.
                self.reflected_nodes.append(f"{name}/{node.GetName()}")
                csolid = self.definitions[child.GetName()][1]
                if csolid is None:
                    self._record(child, reason="reflected placement of a volume with "
                                               "daughters cannot be baked")
                    continue
                try:
                    baked = apply_tgeo_matrix(csolid, m, "reflected placement")
                except ShapeDeclined as e:
                    self._record(child, reason=str(e))
                    continue
                blab = self.shape_tool.AddShape(baked, False)
                TDataStd_Name.Set(blab, f"{child.GetName()}__mirrored")
                comp = self.shape_tool.AddComponent(asm, blab, TopLoc_Location(gp_Trsf()))
                self.nbaked += 1
            else:
                comp = self.shape_tool.AddComponent(asm, clab, TopLoc_Location(t))
                placed_children.append((self.definitions[child.GetName()][1], t))
            TDataStd_Name.Set(comp, str(node.GetName()))
            self.ncomponents += 1

        if emit_body:
            body = occ
            if self.opts.carve_mothers:
                body = self._carve(occ, placed_children, name) or occ
            blab = self.shape_tool.AddShape(body, False)
            TDataStd_Name.Set(blab, f"{name}__body")
            comp = self.shape_tool.AddComponent(asm, blab, TopLoc_Location(gp_Trsf()))
            TDataStd_Name.Set(comp, f"{name}__body")
            self.ncomponents += 1
            rec["bodyComponent"] = f"{name}__body"
        elif occ is not None:
            rec["bodyComponent"] = None
            rec["reason"] = "mother solid omitted (--no-mother-bodies/--skip-top-body)"
            rec["converted"] = False
        if self.ncomponents == ncomp0:
            # every child declined and there is no body: an empty XCAF label reads
            # back as a leaf holding an empty compound, so drop it instead.
            self.shape_tool.RemoveShape(asm)
            self.definitions[name] = (None, occ)
            return None
        return asm

    # ------------------------------------------------------------------

    @staticmethod
    def _world_key(mat, tr):
        return (tuple(round(mat[i][j], 9) for i in range(3) for j in range(3))
                + tuple(round(v, 6) for v in tr))

    @staticmethod
    def _compose(pmat, ptr, m):
        """Compose a parent world transform with a TGeoMatrix, in mm."""
        cmat, ctr = tgeo_matrix_components(m)
        mat = [[sum(pmat[i][k] * cmat[k][j] for k in range(3)) for j in range(3)]
               for i in range(3)]
        tr = [sum(pmat[i][k] * ctr[k] for k in range(3)) + ptr[i] for i in range(3)]
        return mat, tr

    def _leaf_definition(self, vol, depth):
        """The shared XCAF definition of a volume's own solid (built once)."""
        name = str(vol.GetName())
        hit = self.definitions.get(name)
        if hit is not None:
            return hit[0]
        occ, reason = self._solid_for(vol)
        rec = self._record(vol, converted=occ is not None, reason=reason)
        if occ is None:
            self.definitions[name] = (None, None)
            return None
        self._verify(vol, occ, rec)
        lab = self.shape_tool.AddShape(occ, False)
        TDataStd_Name.Set(lab, name)
        self.definitions[name] = (lab, occ)
        return lab

    def _body_definition(self, vol):
        name = str(vol.GetName()) + "__body"
        hit = self.definitions.get(name)
        if hit is not None:
            return hit[0]
        occ, reason = self._solid_for(vol)
        rec = self._record(vol, converted=occ is not None, reason=reason)
        if occ is None:
            self.definitions[name] = (None, None)
            return None
        self._verify(vol, occ, rec)
        rec["bodyComponent"] = name
        lab = self.shape_tool.AddShape(occ, False)
        TDataStd_Name.Set(lab, name)
        self.definitions[name] = (lab, occ)
        return lab

    def build_world(self, vol, depth=0, wmat=None, wtr=None, path=""):
        """Per-occurrence walk that drops coincident (volume, world transform) pairs.

        Assembly labels are one per *occurrence* so that a dropped duplicate deep in
        a shared subtree does not delete its needed siblings; leaf solids stay one
        definition per TGeoVolume, which is what the round trip is about.
        """
        if wmat is None:
            wmat = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
            wtr = [0., 0., 0.]
        name = str(vol.GetName())
        nd = int(vol.GetNdaughters())
        descend = nd > 0 and (self.opts.max_depth is None or depth < self.opts.max_depth)
        if not (self.opts.include_name is None
                or fnmatch.fnmatch(name, self.opts.include_name)):
            return None

        if not descend:
            key = (name, self._world_key(wmat, wtr))
            if key in self.placed_world:
                self.ndropped += 1
                if len(self.dropped_examples) < 50:
                    self.dropped_examples.append(path)
                return None
            lab = self._leaf_definition(vol, depth)
            if lab is None:
                return None
            self.placed_world.add(key)
            return lab

        # Gather the surviving children first: an occurrence whose every child was
        # dropped as a coincident duplicate must not leave an empty assembly label
        # behind, because XCAF hands an empty label to the reader as a leaf with an
        # empty compound and the reverse converter then warns about it, once per
        # occurrence.
        comps = []
        for i in range(nd):
            node = vol.GetNode(i)
            child = node.GetVolume()
            m = node.GetMatrix()
            cmat, ctr = self._compose(wmat, wtr, m)
            clab = self.build_world(child, depth + 1, cmat, ctr,
                                    f"{path}/{node.GetName()}")
            if clab is None:
                continue
            t = tgeo_matrix_to_trsf(m)
            if t is None:
                # A reflecting or scaling placement cannot be a STEP component
                # location, so bake it into a private copy of the child solid.
                self.reflected_nodes.append(f"{name}/{node.GetName()}")
                csolid = self.definitions.get(child.GetName(), (None, None))[1]
                if csolid is None:
                    self._record(child, reason="reflected placement of a volume with "
                                               "daughters cannot be baked")
                    continue
                try:
                    baked = apply_tgeo_matrix(csolid, m, "reflected placement")
                except ShapeDeclined as e:
                    self._record(child, reason=str(e))
                    continue
                blab = self.shape_tool.AddShape(baked, False)
                TDataStd_Name.Set(blab, f"{child.GetName()}__mirrored")
                self.nbaked += 1
                comps.append((blab, gp_Trsf(), str(node.GetName())))
                continue
            comps.append((clab, t, str(node.GetName())))

        if (self.opts.mother_bodies
                and not (depth == 0 and self.opts.skip_top_body)):
            key = (name + "__body", self._world_key(wmat, wtr))
            if key in self.placed_world:
                self.ndropped += 1
            else:
                blab = self._body_definition(vol)
                if blab is not None:
                    self.placed_world.add(key)
                    comps.append((blab, gp_Trsf(), f"{name}__body"))

        if not comps:
            return None
        asm = self.shape_tool.NewShape()
        TDataStd_Name.Set(asm, name)
        for (clab, t, cname) in comps:
            comp = self.shape_tool.AddComponent(asm, clab, TopLoc_Location(t))
            TDataStd_Name.Set(comp, cname)
            self.ncomponents += 1
        return asm

    # ------------------------------------------------------------------

    def _carve(self, mother, placed, name):
        cutters = [_moved(sh, t) for (sh, t) in placed if sh is not None]
        if not cutters:
            return mother
        try:
            tool = cutters[0]
            for c in cutters[1:]:
                tool = _boolean(BRepAlgoAPI_Fuse, tool, c, "carve fuse")
            return _boolean(BRepAlgoAPI_Cut, mother, tool, "carve cut")
        except ShapeDeclined as e:
            self.log(f"  [WARN] {name}: carving failed ({e}); keeping the uncarved mother")
            return mother

    # ------------------------------------------------------------------

    def write(self, path):
        self.shape_tool.UpdateAssemblies()
        Interface_Static.SetCVal("write.step.schema", "AP214IS")
        Interface_Static.SetCVal("write.step.unit", "MM")
        Interface_Static.SetCVal("write.step.product.name", "O2_TGeoToCAD")
        w = STEPCAFControl_Writer()
        w.SetNameMode(True)
        w.SetColorMode(False)
        w.SetLayerMode(False)
        if not w.Transfer(self.doc):
            raise RuntimeError("STEPCAFControl_Writer.Transfer failed")
        if w.Write(path) != IFSelect_RetDone:
            raise RuntimeError(f"STEP write failed for {path}")

    def report(self, source, out_step):
        by_class = {}
        npure = 0
        for r in self.records.values():
            pure = r["isAssembly"] or r["shapeClass"] == "TGeoShapeAssembly"
            c = by_class.setdefault(r["shapeClass"], {"converted": 0, "declined": 0,
                                                      "pureAssembly": 0, "reasons": {}})
            if r["converted"]:
                c["converted"] += 1
            elif pure:
                c["pureAssembly"] += 1
                npure += 1
            else:
                c["declined"] += 1
                key = (r["reason"] or "unknown").split(":")[0]
                c["reasons"][key] = c["reasons"].get(key, 0) + 1
        recs = sorted(self.records.values(), key=lambda r: r["name"])
        devs = [r["relDev"] for r in recs if r.get("relDev") is not None]
        return {
            "source": os.path.abspath(source),
            "output": os.path.abspath(out_step),
            "scaleToMm": SCALE_TO_MM,
            "generator": "O2_TGeoToCAD.py",
            "wallSeconds": round(time.time() - self.t0, 2),
            "definitions": sum(1 for r in recs if r["converted"]),
            "pureAssemblies": npure,
            "declined": sum(1 for r in recs if not r["converted"]
                            and not (r["isAssembly"] or r["shapeClass"] == "TGeoShapeAssembly")),
            "assemblies": sum(1 for r in recs if r["ndaughters"] > 0),
            "components": self.ncomponents,
            "bakedMirrorCopies": self.nbaked,
            "coincidentPlacementsDropped": self.ndropped,
            "coincidentPlacementExamples": self.dropped_examples,
            "reflectedPlacements": self.reflected_nodes[:50],
            "nReflectedPlacements": len(self.reflected_nodes),
            "maxRelDev": max(devs) if devs else None,
            "medianRelDev": sorted(devs)[len(devs) // 2] if devs else None,
            "byShapeClass": by_class,
            "volumes": recs,
        }


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def _cap_check(label, tgeo_shape, band, results, expect_fail=False, occ=None):
    """One capacity-parity check: BRepGProp on our solid vs TGeoShape::Capacity()."""
    try:
        if occ is None:
            occ = shape_to_occ(tgeo_shape, SCALE_TO_MM)
        v = solid_volume_mm3(occ) / 1000.0
        cap = float(tgeo_shape.Capacity())
        rel = abs(v - cap) / abs(cap) if cap else float("inf")
        ok = rel <= band
    except Exception as e:
        v, cap, rel, ok = None, None, None, False
        label = f"{label} [{type(e).__name__}: {e}]"
    passed = (ok != expect_fail)
    results.append((label, passed, cap, v, rel))
    return passed


def _print_suite(title, results):
    fails = [r for r in results if not r[1]]
    print(f"\n--- {title}: {len(results)} checks, {len(fails)} failures")
    for (label, ok, cap, v, rel) in results:
        mark = "ok  " if ok else "FAIL"
        if rel is None:
            print(f"  [{mark}] {label}")
        else:
            print(f"  [{mark}] {label:44s} TGeo {cap:14.6f}  OCC {v:14.6f}  rel {rel:.3e}")
    return len(fails)


def self_test():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    import array

    total = 0
    failures = 0

    import random
    import array as _arr
    rngc = random.Random(4242)

    def mc_volume(shape, n=200000):
        dx, dy, dz = shape.GetDX(), shape.GetDY(), shape.GetDZ()
        o = shape.GetOrigin()
        ox, oy, oz = o[0], o[1], o[2]
        vbox = 8.0 * dx * dy * dz
        pt = _arr.array("d", [0.0, 0.0, 0.0])
        hits = 0
        for _ in range(n):
            pt[0] = ox + rngc.uniform(-dx, dx)
            pt[1] = oy + rngc.uniform(-dy, dy)
            pt[2] = oz + rngc.uniform(-dz, dz)
            if shape.Contains(pt):
                hits += 1
        pf = hits / float(n)
        return pf * vbox, vbox * math.sqrt(max(pf * (1.0 - pf), 1e-15) / n)


    # ---- suite 1: primitives, analytic Capacity() ----
    band = 1e-9
    r1 = []
    _cap_check("TGeoBBox(1,2,3)", ROOT.TGeoBBox("b", 1, 2, 3), band, r1)
    _cap_check("TGeoTube(0,2,5)", ROOT.TGeoTube("t0", 0, 2, 5), band, r1)
    _cap_check("TGeoTube(1,2,5) rmin>0", ROOT.TGeoTube("t1", 1, 2, 5), band, r1)
    _cap_check("TGeoTubeSeg(1,2,5,30,150)", ROOT.TGeoTubeSeg("ts", 1, 2, 5, 30, 150), band, r1)
    _cap_check("TGeoTubeSeg(0,2,5,200,340)", ROOT.TGeoTubeSeg("ts2", 0, 2, 5, 200, 340), band, r1)
    _cap_check("TGeoCone(3,0,2,0,4)", ROOT.TGeoCone("c0", 3, 0, 2, 0, 4), band, r1)
    _cap_check("TGeoCone(3,1,2,0.5,4) rmin>0", ROOT.TGeoCone("c1", 3, 1, 2, 0.5, 4), band, r1)
    _cap_check("TGeoConeSeg(2,.5,1,.7,1.5,30,150)",
               ROOT.TGeoConeSeg("cs", 2, .5, 1, .7, 1.5, 30, 150), band, r1)
    _cap_check("TGeoEltu(2,3,4)", ROOT.TGeoEltu("e", 2, 3, 4), band, r1)
    _cap_check("TGeoTorus(10,0,2)", ROOT.TGeoTorus("to0", 10, 0, 2, 0, 360), band, r1)
    _cap_check("TGeoTorus(10,1,2) hollow", ROOT.TGeoTorus("to1", 10, 1, 2, 0, 360), band, r1)
    _cap_check("TGeoTorus(10,1,2,45,120) wedge",
               ROOT.TGeoTorus("to2", 10, 1, 2, 45, 120), band, r1)
    _cap_check("TGeoTrd1(1,2,3,4)", ROOT.TGeoTrd1("d1", 1, 2, 3, 4), band, r1)
    _cap_check("TGeoTrd2(1,2,3,4,5)", ROOT.TGeoTrd2("d2", 1, 2, 3, 4, 5), band, r1)
    _cap_check("TGeoSphere(0,2) full", ROOT.TGeoSphere("s0", 0, 2, 0, 180, 0, 360), band, r1)
    _cap_check("TGeoSphere(1,2) shell", ROOT.TGeoSphere("s1", 1, 2, 0, 180, 0, 360), band, r1)
    _cap_check("TGeoSphere(0,2,30,120) theta",
               ROOT.TGeoSphere("s2", 0, 2, 30, 120, 0, 360), band, r1)
    _cap_check("TGeoSphere(1,2,30,120,20,200)",
               ROOT.TGeoSphere("s3", 1, 2, 30, 120, 20, 200), band, r1)
    _cap_check("TGeoCtub(0,1,1) straight",
               ROOT.TGeoCtub("ct", 0, 1, 1, 0, 360, 0, 0, -1, 0, 0, 1), band, r1)
    arb8v = array.array("d", [-1, -1, -1, 1, 1, 1, 1, -1, -2, -2, -2, 2, 2, 2, 2, -2])
    _cap_check("TGeoArb8 (pyramid frustum)", ROOT.TGeoArb8("a8", 1.0, arb8v), band, r1)
    _cap_check("TGeoTrap(2,0,0,1,1,1,0,1,1,1,0)",
               ROOT.TGeoTrap("tp", 2, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0), band, r1)
    x = ROOT.TGeoXtru(2)
    x.DefinePolygon(4, array.array("d", [0, 0, 2, 2]), array.array("d", [0, 1, 1, 0]))
    x.DefineSection(0, -1, 0, 0, 1)
    x.DefineSection(1, 1, 0.5, 0, 2)
    _cap_check("TGeoXtru (4-gon, scaled+offset)", x, band, r1)
    pg = ROOT.TGeoPgon("pg", 0, 360, 6, 2)
    pg.DefineSection(0, -1, 0, 1)
    pg.DefineSection(1, 1, 0, 1)
    _cap_check("TGeoPgon(0,360,6) solid", pg, band, r1)
    pg2 = ROOT.TGeoPgon("pg2", 10, 90, 3, 2)
    pg2.DefineSection(0, -1, 0.5, 1)
    pg2.DefineSection(1, 1, 0.5, 1)
    _cap_check("TGeoPgon(10,90,3) hollow wedge", pg2, band, r1)
    pg3 = ROOT.TGeoPgon("pg3", 0, 360, 8, 3)
    pg3.DefineSection(0, -2, 0.5, 1)
    pg3.DefineSection(1, 0, 0.5, 2)
    pg3.DefineSection(2, 2, 0.8, 2)
    _cap_check("TGeoPgon(0,360,8) hollow stack", pg3, band, r1)
    pc = ROOT.TGeoPcon("pc", 0, 360, 3)
    pc.DefineSection(0, -1, 0, 1)
    pc.DefineSection(1, 0, 0, 2)
    pc.DefineSection(2, 1, 0, 2)
    _cap_check("TGeoPcon(0,360) rmin=0", pc, band, r1)
    pc2 = ROOT.TGeoPcon("pc2", 0, 360, 3)
    pc2.DefineSection(0, -1, 0.5, 1)
    pc2.DefineSection(1, 0, 0.5, 2)
    pc2.DefineSection(2, 1, 0.8, 2)
    _cap_check("TGeoPcon(0,360) rmin>0", pc2, band, r1)
    pc3 = ROOT.TGeoPcon("pc3", 20, 150, 4)
    pc3.DefineSection(0, -3, 0.5, 1)
    pc3.DefineSection(1, -1, 0.5, 2)
    pc3.DefineSection(2, -1, 1.2, 2)      # a zero-thickness radius jump
    pc3.DefineSection(3, 2, 1.2, 1.8)
    _cap_check("TGeoPcon(20,150) wedge + z-jump", pc3, band, r1)
    total += len(r1)
    failures += _print_suite("primitives vs TGeoShape::Capacity(), band 1e-9", r1)

    # ---- suite 2: composites, against an independent Monte-Carlo of TGeo itself ----
    # TGeoCompositeShape::Capacity() is itself a Monte-Carlo estimate, so it cannot
    # be a 1e-9 oracle. Instead we run our own MC with a *stated* N, quote its
    # sigma, and require the OCCT volume inside 4 sigma -- and we check that the
    # test could have failed by pricing the negative control the same way.
    r2 = []
    _keep = [ROOT.TGeoBBox("ca", 2, 2, 2), ROOT.TGeoTube("cb", 0, 1, 3)]
    tr = ROOT.TGeoTranslation("shift", 3, 0, 0)
    tr.RegisterYourself()
    rot = ROOT.TGeoRotation("rot90", 0, 90, 0)
    rot.RegisterYourself()
    composites = [
        # TGeoCtub's cut planes replace the end faces, so its z extent is not
        # +-dz and its bounding box is larger: MC over that box is the honest test.
        ("cut tube, slanted (TGeoCtub)",
         ROOT.TGeoCtub("ct2", 0, 1, 1, 0, 360, 0, -0.6, -0.8, 0, 0.6, 0.8)),
        ("box - tube (subtraction)", ROOT.TGeoCompositeShape("sub", "ca - cb")),
        ("box * tube (intersection)", ROOT.TGeoCompositeShape("inter", "ca * cb")),
        ("box + shifted tube (union)", ROOT.TGeoCompositeShape("uni", "ca + cb:shift")),
        ("(box - tube) + shifted tube (nested)",
         ROOT.TGeoCompositeShape("nest", "(ca - cb) + cb:shift")),
        ("box - rotated tube (rotated operand)",
         ROOT.TGeoCompositeShape("rotsub", "ca - cb:rot90")),
    ]
    for (label, cs) in composites:
        try:
            occ = shape_to_occ(cs, SCALE_TO_MM)
            v = solid_volume_mm3(occ) / 1000.0
            vmc, sig = mc_volume(cs)
            ok = abs(v - vmc) <= 4.0 * sig
            print(f"    {label:40s} OCC {v:10.6f}  MC {vmc:10.6f} +- {sig:.4f}"
                  f"  ({abs(v - vmc) / sig:.2f} sigma)")
        except Exception as e:
            ok = False
            label = f"{label} [{type(e).__name__}: {e}]"
        r2.append((label, ok, None, None, None))
    # the control on the control: a 1% wrong volume must be outside 4 sigma
    _, cs0 = composites[0]
    vmc0, sig0 = mc_volume(cs0)
    r2.append((f"a +2% wrong volume would be rejected ({0.02 * vmc0 / sig0:.1f} sigma)",
               abs(1.02 * vmc0 - vmc0) > 4.0 * sig0, None, None, None))
    total += len(r2)
    failures += _print_suite("composites vs an independent MC of TGeo (N=200k, 4 sigma)", r2)

    # ---- suite 2b: point-by-point Contains agreement, TGeo vs OCCT ----
    from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier
    from OCC.Core.TopAbs import TopAbs_IN
    r2b = []
    for (label, cs) in composites[:2] + [("pcon rmin>0", None)]:
        if cs is None:
            cs = ROOT.TGeoPcon("pcx", 0, 360, 3)
            cs.DefineSection(0, -1, 0.5, 1)
            cs.DefineSection(1, 0, 0.5, 2)
            cs.DefineSection(2, 1, 0.8, 2)
        occ = shape_to_occ(cs, SCALE_TO_MM)
        clf = BRepClass3d_SolidClassifier(occ)
        dx, dy, dz = cs.GetDX(), cs.GetDY(), cs.GetDZ()
        o = cs.GetOrigin()
        pt = _arr.array("d", [0.0, 0.0, 0.0])
        bad = skipped = 0
        ntot = 3000
        for _ in range(ntot):
            pt[0] = o[0] + rngc.uniform(-dx, dx)
            pt[1] = o[1] + rngc.uniform(-dy, dy)
            pt[2] = o[2] + rngc.uniform(-dz, dz)
            tin = bool(cs.Contains(pt))
            if cs.Safety(pt, tin) < 1e-6:      # on the surface: not a fair question
                skipped += 1
                continue
            clf.Perform(gp_Pnt(pt[0] * SCALE_TO_MM, pt[1] * SCALE_TO_MM,
                               pt[2] * SCALE_TO_MM), 1e-7)
            if (clf.State() == TopAbs_IN) != tin:
                bad += 1
        print(f"    {label:40s} {ntot - skipped} scored, {bad} disagreements"
              f" ({skipped} within 1e-6 cm of a surface)")
        r2b.append((f"Contains agrees, TGeo vs OCCT: {label}", bad == 0, None, None, None))
    total += len(r2b)
    failures += _print_suite("Contains agreement, TGeo vs OCCT classifier (3000 pts each)", r2b)

    # ---- suite 3: placement transforms ----
    r3 = []
    rng = random.Random(20260822)
    for k, m in enumerate([
        ROOT.TGeoTranslation("t", 1.5, -2.5, 3.5),
        ROOT.TGeoRotation("r", 30, 40, 50),
        ROOT.TGeoCombiTrans("ct", 1, 2, 3, ROOT.TGeoRotation("r2", 11, 22, 33)),
    ]):
        t = tgeo_matrix_to_trsf(m)
        worst = 0.0
        for _ in range(200):
            loc = [rng.uniform(-5, 5) for _ in range(3)]
            mas = array.array("d", [0, 0, 0])
            m.LocalToMaster(array.array("d", loc), mas)
            p = gp_Pnt(loc[0] * SCALE_TO_MM, loc[1] * SCALE_TO_MM, loc[2] * SCALE_TO_MM)
            p.Transform(t)
            worst = max(worst,
                        abs(p.X() - mas[0] * SCALE_TO_MM),
                        abs(p.Y() - mas[1] * SCALE_TO_MM),
                        abs(p.Z() - mas[2] * SCALE_TO_MM))
        ok = worst < 1e-9
        r3.append((f"{m.ClassName()} LocalToMaster vs gp_Trsf (200 pts, mm)", ok,
                   None, None, None))
        print(f"    worst |delta| = {worst:.3e} mm")
    # a reflection must be refused as a rigid placement and offered as a GTrsf
    refl = ROOT.TGeoRotation("refl")
    refl.ReflectZ(True)
    r3.append(("reflecting TGeoRotation refused as a gp_Trsf",
               tgeo_matrix_to_trsf(refl) is None, None, None, None))
    box = shape_to_occ(ROOT.TGeoBBox("rb", 1, 2, 3), SCALE_TO_MM)
    mirrored = apply_tgeo_matrix(box, refl, "reflection test")
    r3.append(("reflected box keeps its volume (baked GTrsf)",
               abs(solid_volume_mm3(mirrored) - solid_volume_mm3(box)) < 1e-6,
               None, None, None))
    total += len(r3)
    failures += _print_suite("placement transforms", r3)

    # ---- suite 4: negative controls ----
    r4 = []
    # each of these compares a deliberately WRONG TGeo shape against our solid for
    # the RIGHT one; the band must reject it, or the band proves nothing.
    good_tube = shape_to_occ(ROOT.TGeoTube("ngt", 1, 2, 5), SCALE_TO_MM)
    _cap_check("wrong rmin: Capacity(0,2,5) vs solid(1,2,5) must FAIL",
               ROOT.TGeoTube("ngt2", 0, 2, 5), 1e-9, r4, expect_fail=True, occ=good_tube)
    good_pcon = shape_to_occ(pc2, SCALE_TO_MM)
    pc2b = ROOT.TGeoPcon("pc2b", 0, 360, 3)
    pc2b.DefineSection(0, -1, 0.5, 1)
    pc2b.DefineSection(1, 0, 0.5, 2)
    pc2b.DefineSection(2, 1, 0.9, 2)          # rmin 0.8 -> 0.9
    _cap_check("wrong pcon rmin (0.8 -> 0.9) must FAIL",
               pc2b, 1e-9, r4, expect_fail=True, occ=good_pcon)
    good_pgon = shape_to_occ(pg, SCALE_TO_MM)
    pgb = ROOT.TGeoPgon("pgb", 0, 360, 7, 2)   # 6 -> 7 edges
    pgb.DefineSection(0, -1, 0, 1)
    pgb.DefineSection(1, 1, 0, 1)
    _cap_check("wrong pgon nedges (6 -> 7) must FAIL",
               pgb, 1e-9, r4, expect_fail=True, occ=good_pgon)
    # and a control on the control: the band accepts the right answer
    _cap_check("same pgon accepted (control on the control)", pg, 1e-9, r4, occ=good_pgon)
    total += len(r4)
    failures += _print_suite("negative controls (a wrong parameter must be rejected)", r4)

    # ---- suite 5: analytic surface types ----
    # The point of these STEPs is to be test material for the reverse converter,
    # whose recogniser keys on the analytic surface type. A flat B-spline is not a
    # plane to it, so "the volume is right" is not enough -- the carriers must be
    # the analytic surfaces TGeo meant.
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.GeomAbs import (
        GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus,
        GeomAbs_BezierSurface, GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
        GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface)
    _SNAME = {GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Cone: "cone",
              GeomAbs_Sphere: "sphere", GeomAbs_Torus: "torus",
              GeomAbs_BezierSurface: "bezier", GeomAbs_BSplineSurface: "bspline",
              GeomAbs_SurfaceOfRevolution: "revolution",
              GeomAbs_SurfaceOfExtrusion: "extrusion",
              GeomAbs_OffsetSurface: "offset", GeomAbs_OtherSurface: "other"}

    def face_types(shape):
        import collections as _c
        c = _c.Counter()
        ex = TopExp_Explorer(shape, TopAbs_FACE)
        while ex.More():
            c[_SNAME.get(BRepAdaptor_Surface(topods.Face(ex.Current())).GetType(), "?")] += 1
            ex.Next()
        return dict(c)

    pgf = ROOT.TGeoPgon("pgf", 0, 360, 6, 2)
    pgf.DefineSection(0, -1, 0.5, 1)
    pgf.DefineSection(1, 1, 0.5, 1)
    vtw = array.array("d", [-1, -1, -1, 1, 1, 1, 1, -1,
                            -1.5, -0.5, -0.5, 1.5, 1.5, 0.5, 0.5, -1.5])
    r6 = []
    for nm, tsh, want in [
        ("TGeoBBox", ROOT.TGeoBBox("fb", 1, 2, 3), {"plane": 6}),
        ("TGeoTube", ROOT.TGeoTube("ft", 1, 2, 5), {"plane": 2, "cylinder": 2}),
        ("TGeoTubeSeg", ROOT.TGeoTubeSeg("fts", 1, 2, 5, 30, 150),
         {"plane": 4, "cylinder": 2}),
        ("TGeoCone", ROOT.TGeoCone("fc", 3, 1, 2, 0.5, 4), {"plane": 2, "cone": 2}),
        ("TGeoPcon", pc2, {"plane": 2, "cylinder": 2, "cone": 2}),
        ("TGeoSphere", ROOT.TGeoSphere("fs", 1, 2, 30, 120, 20, 200),
         {"sphere": 2, "cone": 2, "plane": 2}),
        ("TGeoTorus", ROOT.TGeoTorus("fto", 10, 1, 2, 45, 120), {"torus": 2, "plane": 2}),
        ("TGeoEltu", ROOT.TGeoEltu("fe", 2, 3, 4), {"extrusion": 1, "plane": 2}),
        ("TGeoTrd1", ROOT.TGeoTrd1("fd1", 1, 2, 3, 4), {"plane": 6}),
        ("TGeoTrd2", ROOT.TGeoTrd2("fd2", 1, 2, 3, 4, 5), {"plane": 6}),
        ("TGeoXtru", x, {"plane": 6}),
        ("TGeoPgon hollow", pgf, {"plane": 14}),
        ("TGeoTrap", ROOT.TGeoTrap("ftp", 2, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0), {"plane": 6}),
        ("TGeoArb8 planar", ROOT.TGeoArb8("fa8", 1.0, arb8v), {"plane": 6}),
        ("TGeoArb8 twisted", ROOT.TGeoArb8("fa8t", 1.0, vtw), {"plane": 2, "bspline": 4}),
        ("TGeoCtub", ROOT.TGeoCtub("fct", 0, 1, 1, 0, 360, 0, -0.6, -0.8, 0, 0.6, 0.8),
         {"cylinder": 1, "plane": 2}),
    ]:
        try:
            got = face_types(shape_to_occ(tsh, SCALE_TO_MM))
            ok = got == want
            if not ok:
                nm = f"{nm} (got {got}, want {want})"
        except Exception as e:
            ok = False
            nm = f"{nm} [{type(e).__name__}: {e}]"
        r6.append((f"{nm}", ok, None, None, None))
    total += len(r6)
    failures += _print_suite("analytic surface types of every face", r6)

    # ---- suite 6: the XCAF document, written and read back ----
    r5 = []
    import tempfile
    mgr = ROOT.TGeoManager("stmgr", "self-test")
    vac = mgr.MakeBox("world", ROOT.nullptr, 50, 50, 50)
    mgr.SetTopVolume(vac)
    inner = mgr.MakeTube("innertube", ROOT.nullptr, 1, 2, 5)
    vac.AddNode(inner, 1, ROOT.TGeoTranslation(3, 0, 0))
    vac.AddNode(inner, 2, ROOT.TGeoTranslation(-3, 0, 0))
    grp = mgr.MakeVolumeAssembly("grp")
    leaf = mgr.MakeBox("leafbox", ROOT.nullptr, 1, 1, 1)
    grp.AddNode(leaf, 1, ROOT.TGeoTranslation(0, 4, 0))
    vac.AddNode(grp, 1, ROOT.TGeoTranslation(0, 0, 7))
    mgr.CloseGeometry()

    class _O:
        pass
    o = _O()
    o.quiet = True
    o.verify = True
    o.mother_bodies = True
    o.skip_top_body = False
    o.carve_mothers = False
    o.max_depth = None
    o.include_name = None
    o.dedup_world = False
    conv = TGeoToStep(o)
    conv.build(mgr.GetTopVolume())
    tmp = os.path.join(tempfile.mkdtemp(), "selftest.step")
    conv.write(tmp)
    rep = conv.report("in-memory", tmp)
    r5.append(("STEP file written", os.path.getsize(tmp) > 0, None, None, None))
    r5.append(("one definition per logical volume (3 shaped, 2 assemblies)",
               rep["definitions"] == 3 and rep["assemblies"] == 2, None, None, None))
    r5.append(("shared tube emitted once, placed twice",
               conv.ncomponents == 5, None, None, None))

    from OCC.Core.STEPCAFControl import STEPCAFControl_Reader
    d2 = TDocStd_Document("rb")
    rd = STEPCAFControl_Reader()
    rd.SetNameMode(True)
    r5.append(("STEP reads back", rd.ReadFile(tmp) == IFSelect_RetDone, None, None, None))
    rd.Transfer(d2)
    st2 = XCAFDoc_DocumentTool.ShapeTool(d2.Main())
    roots = TDF_LabelSequence()
    st2.GetFreeShapes(roots)
    r5.append(("exactly one free shape (the top assembly)",
               roots.Length() == 1, None, None, None))
    names = []
    leaves = []

    def walk(lb):
        ch = TDF_LabelSequence()
        st2.GetComponents(lb, ch)
        if ch.Length() == 0:
            leaves.append(lb)
            names.append(lb.GetLabelName())
            return
        for i in range(ch.Length()):
            c = ch.Value(i + 1)
            if st2.IsReference(c):
                ref = TDF_Label()
                st2.GetReferredShape(c, ref)
                walk(ref)
            else:
                walk(c)

    walk(roots.Value(1))
    r5.append(("names survive the write/read (innertube present)",
               "innertube" in names, None, None, None))
    r5.append(("the mother body is a named leaf (world__body)",
               "world__body" in names, None, None, None))
    r5.append((f"leaf occurrences == placements ({len(leaves)} == 4)",
               len(leaves) == 4, None, None, None))
    total += len(r5)
    failures += _print_suite("XCAF assembly document, written and read back", r5)

    print(f"\n{total} checks, {failures} failures")
    sys.stdout.flush()
    sys.stderr.flush()
    # PyROOT double-frees the loose TGeoShapes at interpreter teardown (they are
    # owned by gGeoManager as well), which aborts the process after every check
    # has already run. Leave before teardown so the exit status is the verdict.
    os._exit(1 if failures else 0)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def load_manager(path):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    geo = ROOT.TGeoManager.Import(path)
    if geo is None:
        raise RuntimeError(f"could not import a TGeoManager from {path}")
    return geo


def main(argv=None):
    ap = argparse.ArgumentParser(description="TGeo -> STEP (AP214) converter")
    ap.add_argument("input", nargs="?", help="ROOT geometry file")
    ap.add_argument("output", nargs="?", help="output .step file")
    ap.add_argument("--report", default=None)
    ap.add_argument("--top", default=None)
    ap.add_argument("--max-depth", type=int, default=None)
    ap.add_argument("--include-name", default=None)
    ap.add_argument("--no-mother-bodies", dest="mother_bodies", action="store_false")
    ap.add_argument("--skip-top-body", action="store_true")
    ap.add_argument("--carve-mothers", action="store_true")
    ap.add_argument("--dedup-world", action="store_true")
    ap.add_argument("--no-step", dest="write_step", action="store_false",
                    help="build every solid and write the report, but skip the STEP "
                         "write (which is where OCCT gives out on very large models)")
    ap.add_argument("--no-verify", dest="verify", action="store_false")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    opts = ap.parse_args(argv)

    if opts.self_test:
        return self_test()
    if not opts.input or not opts.output:
        ap.error("input and output are required (or use --self-test)")

    geo = load_manager(opts.input)
    if opts.top:
        vol = geo.GetVolume(opts.top)
        if vol is None:
            raise SystemExit(f"no volume named {opts.top}")
    else:
        vol = geo.GetTopVolume()

    conv = TGeoToStep(opts)
    conv.log(f"walking {vol.GetName()} ...")
    if opts.dedup_world:
        conv.build_world(vol)
    else:
        conv.build(vol)
    conv.log(f"  {len(conv.definitions)} logical volumes, "
             f"{conv.ncomponents} components")
    if opts.write_step:
        conv.log(f"  writing {opts.output}")
        conv.write(opts.output)

    rep = conv.report(opts.input, opts.output)
    rpath = opts.report or (os.path.splitext(opts.output)[0] + "_report.json")
    with open(rpath, "w") as f:
        json.dump(rep, f, indent=1)

    print(f"{rep['definitions']} solids, {rep['assemblies']} volumes with daughters, "
          f"{rep['pureAssemblies']} pure assemblies, {rep['components']} components, "
          f"{rep['declined']} volumes declined")
    if rep["maxRelDev"] is not None:
        print(f"capacity check: max relative deviation {rep['maxRelDev']:.3e}, "
              f"median {rep['medianRelDev']:.3e}")
    if rep["coincidentPlacementsDropped"]:
        print(f"{rep['coincidentPlacementsDropped']} coincident placement(s) dropped "
              f"(--dedup-world); e.g. {rep['coincidentPlacementExamples'][:2]}")
    if rep["nReflectedPlacements"]:
        print(f"{rep['nReflectedPlacements']} reflecting placements baked as mirrored copies")
    for cls, c in sorted(rep["byShapeClass"].items(), key=lambda kv: -(kv[1]["declined"])):
        if c["declined"]:
            print(f"  declined {cls}: {c['declined']} ({c['reasons']})")
    size = (f"{os.path.getsize(opts.output) / 1e6:.2f} MB"
            if opts.write_step else "no STEP written (--no-step)")
    print(f"report: {rpath}   ({rep['wallSeconds']} s, {size})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
