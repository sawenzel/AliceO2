#!/usr/bin/env python3
"""
A Python script, doing a deep STEP/XCAF -> ROOT TGeo conversion.
For now, all CAD solids are simply meshed. The ROOT geometry is build as a C++ ROOT macro
and facet data is stored in binary form to keep disc space minimal.

NEW (03/2026):
  - Optional material/medium emission from a BOM (bill of materials) CSV file.
    The CSV is expected to contain lines like:
      CAD, Mechanical/Part, <PartNumber>, <Rev>, <Name>, <Mass>, <Material>, ...
  - If both a part mass and a CAD volume are available, an effective density is computed
    and used in the emitted TGeoMaterial. Otherwise a reasonable default density is used
    for a few common materials, or 1.0 g/cm^3 as fallback.

Generates (into --output-folder):
  - geom.C (small ROOT macro)
  - facets_<VOLNAME>_<LID>.bin for each leaf logical volume (float32 triangles)
  - surfaces_<VOLNAME>_<LID>.bin for each exactly convertible leaf logical volume
    (with --exact-surfaces auto|required)
  - brep_<VOLNAME>_<LID>.brep, the OCCT BREP of the same leaf solid, scaled to cm
    (with --exact-surfaces auto|required --dump-brep); used as the reference oracle's
    input, see scripts/geometry/CodeReview_Fable.md Section 8

Facet file format (little-endian):
  uint32 nTriangles
  then nTriangles * 9 * float32:
    ax ay az bx by bz cx cy cz

VOLNAME is a filename-safe version of the XCAF label name when available (e.g. "nut"),
and LID is the XCAF label entry (e.g. "0:1:1:7" -> "0_1_1_7") to keep filenames unique.

Naming:
  - C++ variable names stay based on XCAF label entry (e.g. 0:1:1:7) for uniqueness.
  - ROOT object names (TGeoVolume / TGeoTessellated / TGeoVolumeAssembly) use the label's
    human name when available (e.g. "nut", "rod-assembly"), falling back to the entry.

Units:
  - By default, the script tries to detect the STEP LENGTH unit by scanning the STEP file
    header/contents (common patterns like .MILLI. / .CENTI. / .METRE. / INCH / FOOT).
  - You can override with --step-unit {auto,mm,cm,m,in,ft}. TGeo expects cm.

Author:
  - Sandro Wenzel, CERN (02/2026)
  - Material/BOM integration patch (03/2026)
"""

import warnings
warnings.filterwarnings("ignore", message=".*all to deprecated function.*", category=DeprecationWarning)

import argparse
import csv
import json
import math
import re
import struct
import sys
from dataclasses import dataclass
from pathlib import Path as _Path
from typing import Dict, List, Optional, Pattern, Tuple

import numpy as np

from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepTools import breptools, BRepTools_WireExplorer
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Core.Geom import Geom_TrimmedCurve
from OCC.Core.Geom2d import Geom2d_TrimmedCurve
from OCC.Core.GeomConvert import geomconvert
from OCC.Core.Geom2dConvert import geom2dconvert
from OCC.Core.Convert import Convert_TgtThetaOver2
from OCC.Core.GeomAbs import (
    GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus,
    GeomAbs_BezierSurface, GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
    GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface,
    GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse, GeomAbs_Hyperbola, GeomAbs_Parabola,
    GeomAbs_BezierCurve, GeomAbs_BSplineCurve, GeomAbs_OffsetCurve, GeomAbs_OtherCurve,
)
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopAbs import TopAbs_REVERSED, TopAbs_WIRE, TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Core.TopoDS import topods
from OCC.Extend.TopologyUtils import TopologyExplorer

from OCC.Core.STEPCAFControl import STEPCAFControl_Reader
from OCC.Core.TDocStd import TDocStd_Document
from OCC.Core.XCAFDoc import XCAFDoc_DocumentTool
from OCC.Core.IFSelect import IFSelect_RetDone

from OCC.Core.TDF import TDF_Label, TDF_LabelSequence, TDF_Tool
from OCC.Core.TCollection import TCollection_AsciiString
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Trsf

# volume properties for density calcs (may not be present in older pythonOCC builds)
try:
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop_VolumeProperties
    _HAS_VOLPROPS = True
except Exception:
    _HAS_VOLPROPS = False


# -------------------------------
# STEP/XCAF loading
# -------------------------------

def load_step_with_xcaf(path: str):
    doc = TDocStd_Document("pythonocc-doc")
    reader = STEPCAFControl_Reader()
    reader.SetColorMode(True)
    reader.SetNameMode(True)
    reader.SetLayerMode(True)

    status = reader.ReadFile(path)
    if status != IFSelect_RetDone:
        raise RuntimeError(f"STEP read failed for: {path}")

    reader.Transfer(doc)
    shape_tool = XCAFDoc_DocumentTool.ShapeTool(doc.Main())
    return doc, shape_tool


def label_id(label: TDF_Label) -> str:
    s = TCollection_AsciiString()
    TDF_Tool.Entry(label, s)
    return s.ToCString()


def label_name(label: TDF_Label) -> str:
    # Uses the XCAF/STEP name when present; can be empty.
    try:
        n = label.GetLabelName()
        if n:
            return str(n)
    except Exception:
        pass
    return ""


# -------------------------------
# Units
# -------------------------------

def step_unit_scale_to_cm(step_unit: str) -> float:
    step_unit = (step_unit or "auto").lower()
    if step_unit == "mm":
        return 0.1
    if step_unit == "cm":
        return 1.0
    if step_unit == "m":
        return 100.0
    if step_unit == "in":
        return 2.54
    if step_unit == "ft":
        return 30.48
    raise ValueError(f"Unknown --step-unit {step_unit} (use auto, mm, cm, m, in, ft)")


def detect_step_length_unit(step_path: str) -> str:
    """
    Heuristic unit detection by scanning STEP file text for common unit tokens.
    This avoids relying on OCCT APIs that can vary across pythonOCC builds.

    Returns one of: mm, cm, m, in, ft. Defaults to mm if uncertain.
    """
    p = _Path(step_path)
    # STEP can be huge: read only the first few MB; units are near the header.
    max_bytes = 4 * 1024 * 1024
    data = p.open("rb").read(max_bytes).decode("latin-1", errors="ignore").upper()

    if ".MILLI." in data:
        return "mm"
    if ".CENTI." in data:
        return "cm"
    if ".METRE." in data or ".METER." in data:
        return "m"
    if "INCH" in data:
        return "in"
    if "FOOT" in data or "FEET" in data:
        return "ft"

    # Conservative default for mechanical CAD STEP is mm
    return "mm"


@dataclass(frozen=True)
class ClipBox:
    xmin: float
    ymin: float
    zmin: float
    xmax: float
    ymax: float
    zmax: float

    @classmethod
    def from_values(cls, values: List[float]) -> "ClipBox":
        if len(values) != 6:
            raise ValueError("--clip-box expects 6 values: xmin ymin zmin xmax ymax zmax")
        xmin, ymin, zmin, xmax, ymax, zmax = (float(v) for v in values)
        if not (xmin < xmax and ymin < ymax and zmin < zmax):
            raise ValueError("--clip-box requires xmin<xmax, ymin<ymax, and zmin<zmax")
        return cls(xmin, ymin, zmin, xmax, ymax, zmax)

    def as_tuple(self) -> Tuple[float, float, float, float, float, float]:
        return (self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax)


@dataclass(frozen=True)
class NameFilter:
    include: Tuple[Pattern[str], ...]
    exclude: Tuple[Pattern[str], ...]

    @classmethod
    def from_patterns(cls, include: List[str], exclude: List[str], case_sensitive: bool = False) -> "NameFilter":
        flags = 0 if case_sensitive else re.IGNORECASE
        return cls(
            tuple(re.compile(pattern, flags) for pattern in include),
            tuple(re.compile(pattern, flags) for pattern in exclude),
        )

    @property
    def active(self) -> bool:
        return bool(self.include or self.exclude)

    @property
    def has_include(self) -> bool:
        return bool(self.include)

    def _text(self, lid: str, name: str) -> str:
        return f"{name} {lid}".strip()

    def matches_include(self, lid: str, name: str) -> bool:
        text = self._text(lid, name)
        return any(pattern.search(text) for pattern in self.include)

    def matches_exclude(self, lid: str, name: str) -> bool:
        text = self._text(lid, name)
        return any(pattern.search(text) for pattern in self.exclude)


# -------------------------------
# Triangulation helpers
# -------------------------------

def _scale_triangles(triangles, s: float):
    if s == 1.0:
        return triangles
    out = []
    for (a, b, c) in triangles:
        out.append((
            (a[0] * s, a[1] * s, a[2] * s),
            (b[0] * s, b[1] * s, b[2] * s),
            (c[0] * s, c[1] * s, c[2] * s),
        ))
    return out


def triangulate_asbbox(shape, scale_to_cm: float = 1.0):
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()

    p000 = (xmin, ymin, zmin)
    p001 = (xmin, ymin, zmax)
    p010 = (xmin, ymax, zmin)
    p011 = (xmin, ymax, zmax)
    p100 = (xmax, ymin, zmin)
    p101 = (xmax, ymin, zmax)
    p110 = (xmax, ymax, zmin)
    p111 = (xmax, ymax, zmax)

    triangles = [
        (p000, p100, p110), (p000, p110, p010),
        (p001, p111, p101), (p001, p011, p111),
        (p000, p101, p100), (p000, p001, p101),
        (p010, p110, p111), (p010, p111, p011),
        (p000, p010, p011), (p000, p011, p001),
        (p100, p101, p111), (p100, p111, p110),
    ]
    return _scale_triangles(triangles, scale_to_cm)


def triangulate_CAD_solid(my_solid, meshparam, scale_to_cm: float = 1.0):
    lin_defl = float(meshparam.get("lin_defl", 0.1))
    ang_defl = float(meshparam.get("ang_defl", 0.1))

    parallel = True
    try:
        BRepMesh_IncrementalMesh(my_solid, lin_defl, False, ang_defl, bool(parallel))
    except TypeError:
        BRepMesh_IncrementalMesh(my_solid, lin_defl, False, ang_defl)

    triangles = []
    for face in TopologyExplorer(my_solid).faces():
        loc = TopLoc_Location()
        triangulation = BRep_Tool.Triangulation(face, loc)
        if triangulation is None:
            continue

        trsf = loc.Transformation()
        reverse = (face.Orientation() == TopAbs_REVERSED)

        for i in range(1, triangulation.NbTriangles() + 1):
            tri = triangulation.Triangle(i)
            n1, n2, n3 = tri.Get()

            p1 = triangulation.Node(n1).Transformed(trsf)
            p2 = triangulation.Node(n2).Transformed(trsf)
            p3 = triangulation.Node(n3).Transformed(trsf)

            if reverse:
                p2, p3 = p3, p2

            triangles.append((
                (p1.X(), p1.Y(), p1.Z()),
                (p2.X(), p2.Y(), p2.Z()),
                (p3.X(), p3.Y(), p3.Z()),
            ))

    return _scale_triangles(triangles, scale_to_cm)


# -------------------------------
# Volume helpers (for density)
# -------------------------------

def volume_cm3_of_shape(shape, scale_to_cm: float) -> float:
    """Compute CAD solid volume in cm^3 (using STEP->cm scale)."""
    if _HAS_VOLPROPS:
        try:
            props = GProp_GProps()
            brepgprop_VolumeProperties(shape, props)
            # volume returned in STEP length units^3
            v = float(props.Mass())
            return v * (scale_to_cm ** 3)
        except Exception:
            pass

    # Fallback: bounding-box volume (rough but always defined)
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    dx, dy, dz = (xmax - xmin) * scale_to_cm, (ymax - ymin) * scale_to_cm, (zmax - zmin) * scale_to_cm
    return max(dx, 0.0) * max(dy, 0.0) * max(dz, 0.0)


# -------------------------------
# Naming helpers
# -------------------------------

def import_csg_hook():
    """The converter's one CSG integration point: `scripts/geometry/csg/hook.py`.

    Imported lazily and by path, so that the converter behaves exactly as before when `--csg` is
    off, and so that a missing or broken `csg` package can never affect a conversion that did not
    ask for it. See csg/hook.py for what the hook does and Stream_H_CSGEmitter.md for why.
    """
    here = str(_Path(__file__).resolve().parent)
    if here not in sys.path:
        sys.path.insert(0, here)
    from csg import hook
    return hook


def sanitize_cpp_name(s: str) -> str:
    safe = re.sub(r"[^0-9a-zA-Z]", "_", s)
    if not safe:
        safe = "x"
    if not (safe[0].isalpha() or safe[0] == "_"):
        safe = "_" + safe
    return safe


def sanitize_filename(s: str) -> str:
    safe = re.sub(r"[^0-9a-zA-Z]", "_", s)
    return safe or "x"


# -------------------------------
# Binary facet IO
# -------------------------------

def write_facets_bin(path: _Path, triangles):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(struct.pack("<I", len(triangles)))
        for (a, b, c) in triangles:
            f.write(struct.pack(
                "<9f",
                float(a[0]), float(a[1]), float(a[2]),
                float(b[0]), float(b[1]), float(b[2]),
                float(c[0]), float(c[1]), float(c[2]),
            ))


# -------------------------------
# Exact-surface classification probes (--surface-report)
# -------------------------------
#
# These probes classify each face of a leaf CAD solid by its analytic surface type and
# report whether the whole solid could be converted exactly to an O2BVHSurfaceSolid
# (see scripts/geometry/BVHSurfaceSolid.md). They never modify the emitted geometry.

_SURFACE_TYPE_NAMES = {
    GeomAbs_Plane: "plane",
    GeomAbs_Cylinder: "cylinder",
    GeomAbs_Cone: "cone",
    GeomAbs_Sphere: "sphere",
    GeomAbs_Torus: "torus",
    GeomAbs_BezierSurface: "bezier",
    GeomAbs_BSplineSurface: "bspline",
    GeomAbs_SurfaceOfRevolution: "revolution",
    GeomAbs_SurfaceOfExtrusion: "extrusion",
    GeomAbs_OffsetSurface: "offset",
    GeomAbs_OtherSurface: "other",
}

_CURVE_TYPE_NAMES = {
    GeomAbs_Line: "line",
    GeomAbs_Circle: "circle",
    GeomAbs_Ellipse: "ellipse",
    GeomAbs_Hyperbola: "hyperbola",
    GeomAbs_Parabola: "parabola",
    GeomAbs_BezierCurve: "bezier",
    GeomAbs_BSplineCurve: "bspline",
    GeomAbs_OffsetCurve: "offset",
    GeomAbs_OtherCurve: "other",
}

# Current C++ support matrix (keep in sync with O2BVHSurfaceSolid):
#  - planar faces: general line/arc/B-spline boundary wires
#  - cylinder/cone/sphere/torus: parametric-rectangle trims, or a general line/arc/B-spline trim in
#    the (phi, height/theta) or (phiRing, phiTube) domain (holes allowed); the trim must not wrap
#    more than a full turn in phi (for the torus, also not in the tube angle)
_SUPPORTED_SURFACE_TYPES = {"plane", "cylinder", "cone", "sphere", "torus"}
_SUPPORTED_PLANAR_CURVES = {"line", "circle", "bspline", "bezier"}
# Boundary curves whose 2D pcurve the quadric extractor can turn into an exact (phi, v) trim edge:
# a straight line stays a line; circle/ellipse/bezier/B-spline pcurves are converted to a B-spline
# and transformed by the affine (u, v) -> (phi, height/theta) map (a B-spline is closed under it).
_SUPPORTED_QUADRIC_CURVES = {"line", "circle", "ellipse", "bspline", "bezier"}


def _xyz(v, scale: float = 1.0) -> List[float]:
    return [v.X() * scale, v.Y() * scale, v.Z() * scale]


def _surface_params(adaptor: BRepAdaptor_Surface, surf_type: str, scale_to_cm: float) -> dict:
    """Extracts the analytic parameters (lengths in cm, angles in rad) for simple types."""
    s = scale_to_cm
    try:
        if surf_type == "plane":
            pln = adaptor.Plane()
            ax3 = pln.Position()
            return {
                "origin_cm": _xyz(ax3.Location(), s),
                "normal": _xyz(pln.Axis().Direction()),
                "axis_u": _xyz(ax3.XDirection()),
                "axis_v": _xyz(ax3.YDirection()),
            }
        if surf_type == "cylinder":
            cyl = adaptor.Cylinder()
            ax3 = cyl.Position()
            return {
                "origin_cm": _xyz(ax3.Location(), s),
                "axis": _xyz(cyl.Axis().Direction()),
                "ref_axis_u": _xyz(ax3.XDirection()),
                "radius_cm": cyl.Radius() * s,
            }
        if surf_type == "cone":
            cone = adaptor.Cone()
            ax3 = cone.Position()
            return {
                "origin_cm": _xyz(ax3.Location(), s),
                "axis": _xyz(cone.Axis().Direction()),
                "ref_axis_u": _xyz(ax3.XDirection()),
                "ref_radius_cm": cone.RefRadius() * s,
                "half_angle_rad": cone.SemiAngle(),
                "apex_cm": _xyz(cone.Apex(), s),
            }
        if surf_type == "sphere":
            sph = adaptor.Sphere()
            ax3 = sph.Position()
            return {
                "center_cm": _xyz(ax3.Location(), s),
                "polar_axis": _xyz(ax3.Direction()),
                "ref_axis_u": _xyz(ax3.XDirection()),
                "radius_cm": sph.Radius() * s,
            }
        if surf_type == "torus":
            tor = adaptor.Torus()
            ax3 = tor.Position()
            return {
                "center_cm": _xyz(ax3.Location(), s),
                "axis": _xyz(ax3.Direction()),
                "ref_axis_u": _xyz(ax3.XDirection()),
                "major_radius_cm": tor.MajorRadius() * s,
                "minor_radius_cm": tor.MinorRadius() * s,
            }
    except Exception as exc:
        return {"error": f"parameter extraction failed: {exc}"}
    return {}


def _edge_pcurve_is_iso(edge, face, uv_bounds) -> bool:
    """
    True when the edge's 2D pcurve on the face is an iso-parametric segment
    (u = const or v = const within tolerance). Used to detect parametric-rectangle
    trims on quadric faces, the only trim kind the C++ side supports today.
    """
    try:
        curve2d, first, last = BRep_Tool.CurveOnSurface(edge, face)
    except Exception:
        return False
    if curve2d is None:
        return False
    us, vs = [], []
    for i in range(5):
        t = first + (last - first) * i / 4.0
        p = curve2d.Value(t)
        us.append(p.X())
        vs.append(p.Y())
    umin, umax, vmin, vmax = uv_bounds
    tol_u = 1e-6 * max(1.0, abs(umax - umin))
    tol_v = 1e-6 * max(1.0, abs(vmax - vmin))
    return (max(us) - min(us) <= tol_u) or (max(vs) - min(vs) <= tol_v)


def classify_face(face, scale_to_cm: float, recognize_surfaces: bool = True) -> dict:
    """Classifies a single TopoDS face: analytic type, parameters, wires and edges.

    `recognize_surfaces`: for a face whose stored type has no direct extractor, also try the
    canonical-form recognizer (differential normal-field, machine-precision only) and, on success,
    record `recognized_type`/`recognized_residual` - the stored surface type is a statement about
    the exporter, not the geometry (see "Canonical-form recognition" below and BVHSurfaceSolid.md).
    This is a *surface*-only claim: it does not verify the boundary wires are extractable, so it is
    an optimistic upper bound on eligibility, same caveat as the existing `trim_kind` field.
    """
    adaptor = BRepAdaptor_Surface(face)
    surf_type = _SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")

    try:
        uv_bounds = list(breptools.UVBounds(face))
    except Exception:
        uv_bounds = [float("nan")] * 4

    record = {
        "type": surf_type,
        "orientation_reversed": face.Orientation() == TopAbs_REVERSED,
        "uv_bounds": uv_bounds,
        "params": _surface_params(adaptor, surf_type, scale_to_cm),
        "wires": [],
    }

    if recognize_surfaces and surf_type not in _SUPPORTED_SURFACE_TYPES and not any(math.isnan(x) for x in uv_bounds):
        rec = _recognize_analytic_surface(adaptor, uv_bounds)
        if rec is not None:
            record["recognized_type"] = rec["kind"]
            record["recognized_residual"] = rec["residual"]

    try:
        outer_wire = breptools.OuterWire(face)
    except Exception:
        outer_wire = None

    wx = TopExp_Explorer(face, TopAbs_WIRE)
    while wx.More():
        wire = topods.Wire(wx.Current())
        curve_types: Dict[str, int] = {}
        n_edges = 0
        n_degenerated = 0
        all_pcurves_iso = True

        ex = TopExp_Explorer(wire, TopAbs_EDGE)
        while ex.More():
            edge = topods.Edge(ex.Current())
            n_edges += 1
            if BRep_Tool.Degenerated(edge):
                # degenerate edges (sphere poles, cone apex) carry no 3D curve;
                # their pcurves are iso lines by construction
                n_degenerated += 1
            else:
                try:
                    ctype = _CURVE_TYPE_NAMES.get(BRepAdaptor_Curve(edge).GetType(), "unknown")
                except Exception:
                    ctype = "unknown"
                curve_types[ctype] = curve_types.get(ctype, 0) + 1
                if not _edge_pcurve_is_iso(edge, face, uv_bounds):
                    all_pcurves_iso = False
            ex.Next()

        record["wires"].append({
            "outer": bool(outer_wire is not None and wire.IsSame(outer_wire)),
            "n_edges": n_edges,
            "n_degenerated": n_degenerated,
            "curve_types": curve_types,
            "all_pcurves_iso": all_pcurves_iso,
        })
        wx.Next()

    return record


def face_supported(record: dict) -> Tuple[bool, Optional[str]]:
    """Evaluates one classify_face record against the current C++ support matrix."""
    surf_type = record["type"]
    if surf_type not in _SUPPORTED_SURFACE_TYPES:
        recognized = record.get("recognized_type")
        if recognized is not None:
            record["trim_kind"] = "recognized"
            return True, None
        return False, f"unsupported surface type '{surf_type}'"

    curve_types = set()
    for w in record["wires"]:
        curve_types.update(w["curve_types"].keys())

    if surf_type == "plane":
        bad = curve_types - _SUPPORTED_PLANAR_CURVES
        if bad:
            return False, f"plane with unsupported boundary curves: {sorted(bad)}"
        record["trim_kind"] = "wires"
        return True, None

    # quadrics: a parametric-rectangle trim, or a general line/arc/B-spline trim in the (phi, v)
    # domain (holes allowed). Only the boundary-curve type limits eligibility now; a genuinely
    # phi-wrapping (> full turn) trim is still caught at extraction time.
    bad = curve_types - _SUPPORTED_QUADRIC_CURVES
    if bad:
        record["trim_kind"] = "general"
        return False, f"{surf_type} with unsupported trim curves: {sorted(bad)}"
    is_rectangle = len(record["wires"]) == 1 and all(w["all_pcurves_iso"] for w in record["wires"])
    record["trim_kind"] = "parametric-rectangle" if is_rectangle else "general"
    return True, None


def build_surface_report(step_path: str, scale_to_cm: float, recognize_surfaces: bool = True) -> dict:
    """Builds the JSON-serializable exact-conversion eligibility report over def_shapes.

    `recognize_surfaces`: also runs the canonical-form recognition pre-pass on every face whose
    stored type has no direct extractor, and tallies a `recognized_surface_counts` /
    `recognized_stored_type_counts` breakdown in the summary (surface-only, see `classify_face`).
    """
    volumes = {}
    n_eligible = 0
    face_type_counts: Dict[str, int] = {}
    curve_type_counts: Dict[str, int] = {}
    fallback_reasons: Dict[str, int] = {}
    recognized_surface_counts: Dict[str, int] = {}
    recognized_stored_type_counts: Dict[str, int] = {}

    for lid, shape in def_shapes.items():
        faces = []
        for face in TopologyExplorer(shape).faces():
            rec = classify_face(face, scale_to_cm, recognize_surfaces=recognize_surfaces)
            ok, reason = face_supported(rec)
            rec["supported"] = ok
            if reason:
                rec["reason"] = reason
                fallback_reasons[reason] = fallback_reasons.get(reason, 0) + 1
            faces.append(rec)

            face_type_counts[rec["type"]] = face_type_counts.get(rec["type"], 0) + 1
            for w in rec["wires"]:
                for ctype, n in w["curve_types"].items():
                    curve_type_counts[ctype] = curve_type_counts.get(ctype, 0) + n
            recognized_kind = rec.get("recognized_type")
            if recognized_kind is not None:
                recognized_surface_counts[recognized_kind] = recognized_surface_counts.get(recognized_kind, 0) + 1
                recognized_stored_type_counts[rec["type"]] = recognized_stored_type_counts.get(rec["type"], 0) + 1

        eligible = bool(faces) and all(f["supported"] for f in faces)
        if eligible:
            n_eligible += 1
        volumes[lid] = {
            "name": def_names.get(lid, ""),
            "n_faces": len(faces),
            "eligible": eligible,
            "faces": faces,
        }

    return {
        "report_version": 1,
        "step_file": step_path,
        "scale_to_cm": scale_to_cm,
        "summary": {
            "n_volumes": len(volumes),
            "n_eligible": n_eligible,
            "face_type_counts": face_type_counts,
            "curve_type_counts": curve_type_counts,
            "fallback_reasons": fallback_reasons,
            "recognized_surface_counts": recognized_surface_counts,
            "recognized_stored_type_counts": recognized_stored_type_counts,
        },
        "volumes": volumes,
    }


# -------------------------------
# Surface sidecar binary IO
# -------------------------------
# Versioned binary sidecar for exact surfaces (surfaces_*.bin). The format is documented
# in scripts/geometry/BVHSurfaceSolid.md ("Surface sidecar format") and read by the C++
# loader o2::base::LoadSurfaceSolid (DetectorsBase/O2SurfaceSolidIO.h).

SURFACE_SIDECAR_MAGIC = b"O2SS"
# Version 2 appends a float64 model tolerance (cm) to the fixed header.
# Version 3 appends a uint32 edge-table size to the header and, per surface, the ordered list of
# (edgeId, flags) of the face's boundary edges -- the topological identity that turns closure from
# a proximity query into a count. The kernel reads all three.
SURFACE_SIDECAR_VERSION = 3
SURFACE_TYPE_ENUM = {"plane": 1, "cylinder": 2, "cone": 3, "sphere": 4, "torus": 5}
CURVE_TYPE_ENUM = {"line": 0, "arc": 1, "bspline": 2}
SURFACE_FLAG_INNER_WALL = 1 << 0

# Per-boundary-edge flag bits, version 3.
EDGE_FLAG_REVERSED = 1 << 0    # the face traverses the edge against the edge's own direction
EDGE_FLAG_DEGENERATE = 1 << 1  # BRep_Tool.Degenerated: a cone apex / sphere pole, no 3D curve
EDGE_FLAG_ANCHORED = 1 << 2    # entry i is trim curve i of this surface, in flattened wire order


def build_edge_table(shape):
    """Index every TopoDS_Edge of \\a shape once, and return (map, edge_id).

    `edge_id(edge)` is a 0-based `uint32` stable for the whole solid: it is what makes two faces'
    trim curves *the same edge* rather than two curves that happen to run close together. The
    measurement that licenses keying on it is in `Stream_F_EdgeIdentity.md`: over both corpora
    (9 fixtures, 12 Bagger parts, 546 edges) every edge is visited exactly twice by
    `BRepTools_WireExplorer` -- the very walk the trim extractors use -- with opposite orientation,
    a seam edge being the two visits of one face and every other edge one visit of each of two.
    """
    edge_map = TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape, TopAbs_EDGE, edge_map)

    def edge_id(edge) -> int:
        return edge_map.FindIndex(edge) - 1  # FindIndex is 1-based; 0 means "not in the map"

    return edge_map, edge_id


def face_boundary_edge_refs(face, edge_id, anchored: bool) -> List[Tuple[int, int]]:
    """The face's boundary edges as ordered (edgeId, flags) pairs.

    The order is `_face_wire_edges(face)` order -- wire by wire, edge by edge -- which is exactly
    the order `extract_planar_face`, `_quadric_trim_wire` and `_recognized_quadric_wire_block`
    emit their sidecar trim curves in. `anchored` says whether the record actually carries that
    wire block: a quadric whose trim is the plain parametric rectangle carries its trim in the
    scalar params instead, so its edge ids are still a complete statement of *which* edges the
    face owns (which is all the closure verdict needs) but there is no trim curve to hang the
    geometric deviation measurement on.
    """
    refs: List[Tuple[int, int]] = []
    base_flags = EDGE_FLAG_ANCHORED if anchored else 0
    for _wire, _is_outer, edges in _face_wire_edges(face):
        for edge, _start_vertex in edges:
            flags = base_flags
            if edge.Orientation() == TopAbs_REVERSED:
                flags |= EDGE_FLAG_REVERSED
            if BRep_Tool.Degenerated(edge):
                flags |= EDGE_FLAG_DEGENERATE
            refs.append((edge_id(edge), flags))
    return refs


def shape_model_tolerance(shape, scale_to_cm: float) -> float:
    """The model's own statement about how well its boundary is defined, in cm.

    This is the largest BRep tolerance over the shape's faces, edges and vertices, converted to
    cm. It is what the kernel needs and cannot compute: without it there is no way to know what
    epsilon two faces of the same imported solid should be expected to agree to, so every closure
    or adjacency decision falls back on a constant nobody chose for this geometry. The oracle
    reports the same quantity (occtOracle.shape_tolerance) and uses it as its gate band.
    """
    worst = 0.0
    for shape_type, getter in ((TopAbs_FACE, lambda sub: BRep_Tool.Tolerance(topods.Face(sub))),
                               (TopAbs_EDGE, lambda sub: BRep_Tool.Tolerance(topods.Edge(sub))),
                               (TopAbs_VERTEX, lambda sub: BRep_Tool.Tolerance(topods.Vertex(sub)))):
        explorer = TopExp_Explorer(shape, shape_type)
        while explorer.More():
            worst = max(worst, getter(explorer.Current()))
            explorer.Next()
    return worst * scale_to_cm


def write_surfaces_bin(path: _Path, surfaces: List[dict], model_tolerance_cm: float = 0.0,
                       n_model_edges: int = 0):
    """
    Writes the surface sidecar. `surfaces` is a list of records:
      {"type": "plane"|"cylinder"|"cone"|"sphere"|"torus",
       "inner_wall": bool,              # quadrics only, default False
       "params": [float, ...],          # fixed per-type layout, see BVHSurfaceSolid.md
       "wires": [                       # planes always; quadrics/torus when the trim is not the
                                        # plain parametric rectangle (else carried in params)
          {"role": "outer"|"inner",     # general line/arc loops (polygon, disk, rounded rect)
           "edges": [{"curve": "line"|"arc", "params": [float, ...]}, ...]},
       ],
       "edge_refs": [(edgeId, flags), ...]}   # version 3; see face_boundary_edge_refs

    `model_tolerance_cm` is the source model's declared tolerance (shape_model_tolerance), written
    into the version-2 header. Zero means "not stated"; the reader then falls back on its own
    documented constant and says so.

    `n_model_edges` is the size of the solid's edge table (`build_edge_table`), written into the
    version-3 header so the reader can size its incidence count without scanning first. Zero, or
    a record with no `edge_refs`, means this sidecar states no edge identity and the reader falls
    back on the geometric rim measurement -- which is exactly what a version-2 file does.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(SURFACE_SIDECAR_MAGIC)
        f.write(struct.pack("<III", SURFACE_SIDECAR_VERSION, len(surfaces), 0))
        f.write(struct.pack("<d", float(model_tolerance_cm)))
        f.write(struct.pack("<I", int(n_model_edges)))
        for srec in surfaces:
            stype = SURFACE_TYPE_ENUM[srec["type"]]
            flags = SURFACE_FLAG_INNER_WALL if srec.get("inner_wall") else 0
            params = [float(p) for p in srec.get("params", [])]
            f.write(struct.pack("<III", stype, flags, len(params)))
            if params:
                f.write(struct.pack(f"<{len(params)}d", *params))
            wires = srec.get("wires", [])
            f.write(struct.pack("<I", len(wires)))
            for w in wires:
                role = 0 if w.get("role", "outer") == "outer" else 1
                edges = w.get("edges", [])
                f.write(struct.pack("<II", role, len(edges)))
                for e in edges:
                    ctype = CURVE_TYPE_ENUM[e["curve"]]
                    cparams = [float(x) for x in e["params"]]
                    f.write(struct.pack("<II", ctype, len(cparams)))
                    if cparams:
                        f.write(struct.pack(f"<{len(cparams)}d", *cparams))
            edge_refs = srec.get("edge_refs", [])
            f.write(struct.pack("<I", len(edge_refs)))
            for edge_id, edge_flags in edge_refs:
                f.write(struct.pack("<IB", int(edge_id), int(edge_flags)))


def write_brep_cm(path: _Path, shape, scale_to_cm: float):
    """
    Write `shape` as an OCCT BREP file, scaled from the STEP file's own length unit to cm.

    The exact-surface extractors and the facet writer both scale their coordinates by
    `scale_to_cm` while the TopoDS shape itself stays in the STEP's units, so a BREP dumped
    as-is would silently be in mm for a typical CAD model. Scaling here keeps the BREP, the
    surfaces_*.bin sidecar and the facets_*.bin mesh in one common unit (cm), which is what the
    OCCT reference oracle needs in order to answer queries about the very same solid.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if scale_to_cm != 1.0:
        trsf = gp_Trsf()
        trsf.SetScale(gp_Pnt(0.0, 0.0, 0.0), scale_to_cm)
        shape = BRepBuilderAPI_Transform(shape, trsf, True).Shape()
    if not breptools.Write(shape, str(path)):
        raise RuntimeError(f"failed to write BREP file {path}")


# -------------------------------
# Exact-surface extraction (--exact-surfaces)
# -------------------------------
# Turn OpenCascade faces into the sidecar surface records written by write_surfaces_bin.
# This is the exact counterpart of the tessellation path: a leaf solid is emitted as an
# O2BVHSurfaceSolid only when *every* one of its faces can be extracted exactly; otherwise
# the caller keeps the tessellated fallback (auto mode) or fails (required mode).
#
# Milestone status (scripts/geometry/BVHSurfaceSolid.md): planar faces carry general line/arc/
# B-spline boundary wires (polygon, disk/annulus, rounded rectangle, slot, spline-bounded plate,
# ...) and the three quadrics (cylinder, cone, sphere) plus the torus are extracted here. A
# quadric/torus face with a single iso-parametric wire stays on the scalar parametric-rectangle
# path (byte-identical); any other supported trim (a hole/window, or a non-iso line/arc/B-spline
# boundary) is emitted as a general trim wire in the surface's (phi, height/theta) or
# (phiRing, phiTube) domain. Curved pcurves (circle/ellipse/Bezier/B-spline) are converted to a
# B-spline whose poles are pushed through the *affine* (u, v) -> (phi, height/theta) map, which is
# exact (a B-spline is closed under affine maps; an anisotropic map merely turns a circle into an
# exactly-represented ellipse). Free-form B-spline/other *surface* types still force the
# tessellated fallback.

_EXTRACT_TOL = 1.e-7


def _v_dot(a: List[float], b: List[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _v_cross(a: List[float], b: List[float]) -> List[float]:
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _planar_frame(face, scale_to_cm: float):
    """Return (origin_cm, axis_u, axis_v) for a planar face, with axisU x axisV pointing
    along the *outward* face normal (OCC surface normal for a FORWARD face, its opposite for
    a REVERSED one). Robust to a left-handed OCC ax3 where XDirection x YDirection = -normal."""
    adaptor = BRepAdaptor_Surface(face)
    pln = adaptor.Plane()
    ax3 = pln.Position()
    s = scale_to_cm
    origin_cm = _xyz(ax3.Location(), s)
    axis_u = _xyz(ax3.XDirection())
    ydir = _xyz(ax3.YDirection())
    normal = _xyz(pln.Axis().Direction())
    outward_sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
    if outward_sign * _v_dot(_v_cross(axis_u, ydir), normal) > 0.0:
        axis_v = ydir
    else:
        axis_v = [-c for c in ydir]
    return origin_cm, axis_u, axis_v


def _face_wire_edges(face):
    """Yield (wire, is_outer, [edges-in-connected-order]) for every wire of the face."""
    try:
        outer_wire = breptools.OuterWire(face)
    except Exception:
        outer_wire = None
    wx = TopExp_Explorer(face, TopAbs_WIRE)
    while wx.More():
        wire = topods.Wire(wx.Current())
        edges = []
        we = BRepTools_WireExplorer(wire, face)
        while we.More():
            edges.append((we.Current(), we.CurrentVertex()))
            we.Next()
        is_outer = outer_wire is not None and wire.IsSame(outer_wire)
        yield wire, is_outer, edges
        wx.Next()


def _planar_projector(origin_cm, axis_u, axis_v, s):
    """Return project(gp_Pnt) -> (u, v): the point's plane-local coordinates in cm."""
    def project(pnt) -> Tuple[float, float]:
        rel = [pnt.X() * s - origin_cm[0], pnt.Y() * s - origin_cm[1], pnt.Z() * s - origin_cm[2]]
        return _v_dot(rel, axis_u), _v_dot(rel, axis_v)
    return project


def _arc_edge_params(edge, curve, project, s) -> Tuple[Optional[List[float]], Optional[str]]:
    """Build [cu, cv, radius, startAngle, phiSweep] for a circular boundary edge.

    The signed sweep is recovered by sampling the 3D edge in *wire-traversal* order (the edge
    is walked backwards when its orientation is REVERSED relative to the underlying curve),
    projecting each sample into the plane frame and unwrapping the polar angle. This is robust
    to full circles (single periodic edge -> +/-2pi), arcs wider than pi, and either winding.
    """
    circ = curve.Circle()
    cu, cv = project(circ.Location())
    radius = circ.Radius() * s
    first, last = curve.FirstParameter(), curve.LastParameter()
    reversed_edge = edge.Orientation() == TopAbs_REVERSED
    angles: List[float] = []
    for tau in (0.0, 0.25, 0.5, 0.75, 1.0):
        t = (1.0 - tau) if reversed_edge else tau
        u, v = project(curve.Value(first + t * (last - first)))
        angles.append(math.atan2(v - cv, u - cu))
    unwrapped = [angles[0]]
    for a in angles[1:]:
        d = a - unwrapped[-1]
        d -= 2.0 * math.pi * math.floor((d + math.pi) / (2.0 * math.pi))  # wrap into (-pi, pi]
        unwrapped.append(unwrapped[-1] + d)
    sweep = unwrapped[-1] - unwrapped[0]
    if abs(sweep) < _EXTRACT_TOL:
        return None, "planar arc edge has a degenerate sweep"
    return [cu, cv, radius, unwrapped[0], sweep], None


def _bspline_flat_params(first: float, last: float, reversed_edge: bool, pole_xform, to_bspline):
    """Flat sidecar B-spline record [degree, nPoles, poles(2*nPoles), weights(nPoles),
    knots(nPoles+degree+1)] for a curve segment [first, last].

    `to_bspline(lo, hi)` trims the source curve to [lo, hi] and returns a clamped (Geom or Geom2d)
    BSplineCurve; `pole_xform(pole)` maps one control point to its output (u, v). The curve is
    trimmed *before* conversion so the parametrisation matches the edge; a periodic result is made
    non-periodic. Poles/weights/knots are reversed when the edge runs opposite the curve."""
    lo, hi = (first, last) if first <= last else (last, first)
    bs = to_bspline(lo, hi)
    if bs is None:
        return None
    if bs.IsPeriodic():
        bs.SetNotPeriodic()
    degree = bs.Degree()
    nb = bs.NbPoles()
    if degree < 1 or nb < degree + 1:
        return None
    poles = []
    weights = []
    for i in range(1, nb + 1):
        u, v = pole_xform(bs.Pole(i))
        poles.append((u, v))
        weights.append(bs.Weight(i))
    flat = []
    for i in range(1, bs.NbKnots() + 1):
        flat.extend([bs.Knot(i)] * bs.Multiplicity(i))
    if len(flat) != nb + degree + 1:
        return None
    if reversed_edge:
        poles.reverse()
        weights.reverse()
        span = flat[0] + flat[-1]
        flat = [span - k for k in reversed(flat)]
    params = [float(degree), float(nb)]
    for u, v in poles:
        params.extend([float(u), float(v)])
    params.extend(float(w) for w in weights)
    params.extend(float(k) for k in flat)
    return params


# Relative residual below which a sampled trim curve is accepted as EXACTLY a line or a circle.
# Deliberately at machine precision and relative to the curve's own extent, mirroring the
# discipline of the surface-side canonical recognition: a curve that is merely *almost* a circle
# must stay a B-spline, because silently replacing geometry is worse than a slow trim.
_CANONICAL_CURVE_TOL = 1.e-9


def _recognize_canonical_curve(samples, poles=None):
    """Recognize a sampled 2D trim curve as an exact line or circle in its output domain.

    CAD kernels routinely write an exact straight line or an exact circle as a B-spline; that is
    the same curve in a heavier representation, so recognizing it is exact recognition and not
    fitting. The payoff is double (see scripts/geometry/ExactTrimTopology.md, item 3): both faces
    meeting along such a seam derive it analytically and therefore agree by construction, and the
    kernel's point-in-trim test drops to its cheap line/arc path instead of flattening a B-spline
    to a polyline.

    `samples` are points already mapped into the *output* (u, v) / (phi, v) domain, in edge
    direction. `poles`, when given, are the B-spline control points in the same domain: collinear
    poles *prove* the curve is a straight segment (it lies in their convex hull), which is a
    stronger statement than agreeing with a line at the sample points.

    Returns ("line", [u0, v0, u1, v1]), ("arc", [cu, cv, r, a0, sweep]) or (None, None).
    """
    points = np.asarray(samples, dtype=float)
    if len(points) < 3:
        return None, None
    extent = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0)))
    if extent < _EXTRACT_TOL:
        return None, None

    # --- straight line
    chord = points[-1] - points[0]
    chord_length = float(np.linalg.norm(chord))
    if chord_length > _EXTRACT_TOL:
        unit = chord / chord_length
        def off_axis(candidates):
            rel = candidates - points[0]
            return float(np.abs(rel[:, 0] * unit[1] - rel[:, 1] * unit[0]).max() / extent)
        straight = off_axis(points) < _CANONICAL_CURVE_TOL
        if straight and poles is not None and len(poles) >= 2:
            straight = off_axis(np.asarray(poles, dtype=float)) < _CANONICAL_CURVE_TOL
        if straight:
            # reject a curve that doubles back along its own chord: geometrically it is not the
            # segment from the first point to the last one, however collinear the samples are
            along = (points - points[0]) @ unit
            if np.all(np.diff(along) >= -_CANONICAL_CURVE_TOL * extent):
                return "line", [float(points[0][0]), float(points[0][1]),
                                float(points[-1][0]), float(points[-1][1])]

    # --- circle: |P - C|^2 = R^2 linearized as 2 P.C + (R^2 - |C|^2) = |P|^2, one least-squares
    # solve with no initial guess. A closed loop (zero chord) lands here as well as an open arc.
    matrix = np.column_stack([2.0 * points, np.ones(len(points))])
    solution, *_ = np.linalg.lstsq(matrix, np.einsum('ij,ij->i', points, points), rcond=None)
    centre = solution[:2]
    radius_sq = solution[2] + float(centre @ centre)
    if radius_sq <= 0.0:
        return None, None
    radius = math.sqrt(radius_sq)
    if float(np.abs(np.linalg.norm(points - centre, axis=1) - radius).max() / extent) >= _CANONICAL_CURVE_TOL:
        return None, None
    # sweep by accumulating signed angle steps, so a full turn and the traversal sense survive
    angles = np.arctan2(points[:, 1] - centre[1], points[:, 0] - centre[0])
    steps = np.diff(angles)
    steps = (steps + math.pi) % (2.0 * math.pi) - math.pi
    sweep = float(steps.sum())
    if abs(sweep) < _EXTRACT_TOL:
        return None, None
    return "arc", [float(centre[0]), float(centre[1]), radius, float(angles[0]), sweep]


def _sample_curve_in_domain(curve, first, last, reversed_edge, point_map, n=64):
    """Sample an OCC curve over [first, last] and map each point into the output domain, ordered
    along the edge. `point_map(p)` takes the curve's own point type to an output (u, v)."""
    lo, hi = (first, last) if first <= last else (last, first)
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi - lo <= 0.0:
        return None
    try:
        samples = [point_map(curve.Value(float(t))) for t in np.linspace(lo, hi, n)]
    except Exception:
        return None
    if reversed_edge:
        samples.reverse()
    return samples


def _planar_bspline_edge_params(edge, project) -> Optional[List[float]]:
    """Sidecar B-spline record for a planar face's B-spline / Bezier boundary edge.

    The 3D boundary curve lies in the plane, so projecting its control poles into the plane frame
    (an affine map) yields the exact 2D B-spline. Returns None on failure (caller falls back)."""
    try:
        curve3d, first, last = BRep_Tool.Curve(edge)
        if curve3d is None:
            return None
        reversed_edge = edge.Orientation() == TopAbs_REVERSED

        def to_bspline(lo, hi):
            trimmed = Geom_TrimmedCurve(curve3d, lo, hi)
            return geomconvert.CurveToBSplineCurve(trimmed, Convert_TgtThetaOver2)

        return _bspline_flat_params(first, last, reversed_edge, project, to_bspline)
    except Exception:
        return None


def _planar_canonical_edge(edge, project, params):
    """Recognize a planar face's B-spline boundary edge as an exact line or arc in the plane frame.

    `params` is the already-extracted flat B-spline record, whose poles are reused as the convex
    hull evidence for straightness. Returns ("line"|"arc", canonical_params) or (None, None)."""
    try:
        curve3d, first, last = BRep_Tool.Curve(edge)
        if curve3d is None:
            return None, None
        reversed_edge = edge.Orientation() == TopAbs_REVERSED
        samples = _sample_curve_in_domain(curve3d, first, last, reversed_edge, project)
        if not samples:
            return None, None
        n_poles = int(params[1])
        poles = [(params[2 + 2 * i], params[3 + 2 * i]) for i in range(n_poles)]
        return _recognize_canonical_curve(samples, poles)
    except Exception:
        return None, None


def extract_planar_face(face, scale_to_cm: float, frame_override=None) -> Tuple[Optional[dict], Optional[str]]:
    """Convert a planar TopoDS face into a sidecar 'plane' surface record with general
    line/arc/B-spline boundary wires (polygon, disk/annulus, rounded rectangle, slot, ...).

    Each wire is walked in connected order; straight edges become 'line' segments and circular
    edges become 'arc' segments in the plane's local (u, v) frame. Boundary edges that are
    neither straight lines nor circular arcs (ellipses, splines, ...) force a fallback.

    `frame_override`, when given, is an (origin_cm, axis_u, axis_v) triple that replaces the
    face's own OCC plane frame (used by the canonical-recognition pre-pass to extract a face
    whose *stored* surface type is not `GeomAbs_Plane` but was recognized as flat): all boundary
    edges are still read from the face's actual 3D curves, so this only bypasses the
    `adaptor.GetType() != GeomAbs_Plane` guard and the OCC-derived frame.
    """
    if frame_override is not None:
        origin_cm, axis_u, axis_v = frame_override
    else:
        adaptor = BRepAdaptor_Surface(face)
        if adaptor.GetType() != GeomAbs_Plane:
            return None, f"not a plane ({_SURFACE_TYPE_NAMES.get(adaptor.GetType(), 'unknown')})"
        origin_cm, axis_u, axis_v = _planar_frame(face, scale_to_cm)
    s = scale_to_cm
    project = _planar_projector(origin_cm, axis_u, axis_v, s)

    wires_out: List[dict] = []
    for wire, is_outer, edges in _face_wire_edges(face):
        classified = []  # (edge, curve, geom_type, projected start (u, v))
        for edge, start_vertex in edges:
            if BRep_Tool.Degenerated(edge):
                return None, "planar face has a degenerated boundary edge"
            try:
                curve = BRepAdaptor_Curve(edge)
                gt = curve.GetType()
            except Exception:
                gt = None
            if gt not in (GeomAbs_Line, GeomAbs_Circle, GeomAbs_BSplineCurve, GeomAbs_BezierCurve):
                name = _CURVE_TYPE_NAMES.get(gt, "unknown")
                return None, f"planar boundary edge is a {name} curve (only line/circle/bspline supported)"
            classified.append((edge, curve, gt, project(BRep_Tool.Pnt(start_vertex))))

        n = len(classified)
        if n == 0:
            return None, "planar face has an empty wire"

        # Canonical-form pre-pass, as in _quadric_trim_wire. The plane projection is an isometry,
        # so unlike the quadric domains a circle written as a B-spline stays a circle here and
        # recovers the cheap exact arc path. Resolved before the polygon check below, so a wire
        # whose B-splines all turn out to be straight is validated as the polygon it now is.
        resolved = []  # per edge: ("line", None) | ("arc", params) | ("bspline", params)
        for edge, curve, gt, _start_uv in classified:
            if gt == GeomAbs_Line:
                resolved.append(("line", None))
            elif gt == GeomAbs_Circle:
                params, reason = _arc_edge_params(edge, curve, project, s)
                if params is None:
                    return None, reason
                resolved.append(("arc", params))
            else:  # B-spline / Bezier: project the 3D poles into the plane frame
                params = _planar_bspline_edge_params(edge, project)
                if params is None:
                    return None, "planar B-spline boundary edge extraction failed"
                canonical_kind, canonical = _planar_canonical_edge(edge, project, params)
                if canonical_kind == "line":
                    resolved.append(("line", None))
                elif canonical_kind == "arc":
                    resolved.append(("arc", canonical))
                else:
                    resolved.append(("bspline", params))

        n_curved = sum(1 for kind, _ in resolved if kind != "line")
        if n_curved == 0 and n < 3:
            return None, "planar polygon wire has fewer than 3 edges"

        seg_edges = []
        for i, (kind, params) in enumerate(resolved):
            if kind == "line":
                u0, v0 = classified[i][3]
                u1, v1 = classified[(i + 1) % n][3]
                seg_edges.append({"curve": "line", "params": [u0, v0, u1, v1]})
            else:
                seg_edges.append({"curve": kind, "params": params})
        wires_out.append({"role": "outer" if is_outer else "inner", "edges": seg_edges})

    if not wires_out:
        return None, "planar face has no wires"
    n_outer = sum(1 for w in wires_out if w["role"] == "outer")
    if n_outer != 1:
        return None, f"planar face has {n_outer} outer wires (expected exactly 1)"

    return {"type": "plane", "params": list(origin_cm) + list(axis_u) + list(axis_v), "wires": wires_out}, None


def _quadric_phi_range(ax3, umin: float, umax: float) -> Tuple[float, float]:
    """Map an OCC angular U-range [umin, umax] to the C++ (phiStart, phiSweep).

    The C++ bounded quadrics measure phi in a right-handed frame with YDir = axis x refU.
    OCC's stored YDirection equals that only for a *direct* (right-handed) gp_Ax3; otherwise
    it is negated, so a point at OCC parameter u sits at C++ phi = -u and the range mirrors.
    Returns a positive sweep clamped into (0, 2pi].
    """
    sweep = umax - umin
    two_pi = 2.0 * math.pi
    if sweep <= 0.0:
        sweep += two_pi
    sweep = min(sweep, two_pi)
    phi_start = umin if ax3.Direct() else -umax
    return phi_start, sweep


def _quadric_trim_wire(face, map_uv) -> Tuple[Optional[List[dict]], Optional[str]]:
    """Build general line/arc/B-spline trim wires in a quadric face's parametric (phi, v) domain.

    \a map_uv(u, v) is the *affine* map from the OCC face parameters to the C++ (phi, height/theta)
    domain (phi = +/-u + offset; height/theta linear in v). Straight pcurve edges stay lines; every
    curved pcurve (circle/ellipse/Bezier/B-spline) is converted to a B-spline and its poles are
    transformed by \a map_uv. Because the map is affine and a B-spline is closed under affine maps,
    the transformed poles describe the trim edge *exactly* in (phi, v) - even for a circle, which an
    anisotropic map turns into an ellipse (represented exactly as the rational B-spline). Line edges
    keep the shared-vertex chaining (each line ends at the next edge's start) for a tight closure.
    Returns (wires, None) with exactly one outer wire plus inner holes, or (None, reason).
    """
    wires_out: List[dict] = []
    for _wire, is_outer, edges in _face_wire_edges(face):
        parsed = []  # per edge: {"kind": "line", "start": (phi, v)} or {"kind": "bspline", ...}
        for edge, _start_vertex in edges:
            curve2d, first, last = BRep_Tool.CurveOnSurface(edge, face)
            if curve2d is None:
                return None, "quadric boundary edge has no 2D pcurve"
            reversed_edge = edge.Orientation() == TopAbs_REVERSED
            ctype = Geom2dAdaptor_Curve(curve2d).GetType()
            if ctype == GeomAbs_Line:
                param = last if reversed_edge else first
                p = curve2d.Value(param)
                parsed.append({"kind": "line", "start": map_uv(p.X(), p.Y())})
            elif ctype in (GeomAbs_Circle, GeomAbs_Ellipse, GeomAbs_BSplineCurve, GeomAbs_BezierCurve):
                def to_bspline(lo, hi, c2=curve2d):
                    trimmed = Geom2d_TrimmedCurve(c2, lo, hi)
                    return geom2dconvert.CurveToBSplineCurve(trimmed, Convert_TgtThetaOver2)

                params = _bspline_flat_params(first, last, reversed_edge,
                                              lambda p: map_uv(p.X(), p.Y()), to_bspline)
                if params is None:
                    return None, "quadric B-spline pcurve extraction failed"
                # Canonical-form pre-pass: a pcurve stored as a B-spline is very often exactly a
                # straight line in the output (phi, v) domain -- the seam of a periodic surface and
                # every iso-parametric trim edge are, and `map_uv` is affine so straightness carries
                # over. Storing it as a line is exact, and it is what makes both faces of that seam
                # agree analytically instead of each flattening its own polyline.
                samples = _sample_curve_in_domain(curve2d, first, last, reversed_edge,
                                                  lambda p: map_uv(p.X(), p.Y()))
                n_poles = int(params[1])
                poles = [(params[2 + 2 * i], params[3 + 2 * i]) for i in range(n_poles)]
                kind, canonical = _recognize_canonical_curve(samples, poles) if samples else (None, None)
                if kind == "line":
                    parsed.append({"kind": "line", "start": (canonical[0], canonical[1])})
                elif kind == "arc":
                    parsed.append({"kind": "arc", "params": canonical,
                                   "start": (canonical[0] + canonical[2] * math.cos(canonical[3]),
                                             canonical[1] + canonical[2] * math.sin(canonical[3]))})
                else:
                    parsed.append({"kind": "bspline", "params": params, "start": (params[2], params[3])})
            else:
                name = _CURVE_TYPE_NAMES.get(ctype, "unknown")
                return None, f"quadric boundary pcurve is a {name} curve (unsupported)"
        n = len(parsed)
        if n == 0:
            return None, "quadric trim wire has no edges"
        if all(p["kind"] == "line" for p in parsed) and n < 3:
            return None, "quadric line trim wire has fewer than 3 edges"
        seg_edges = []
        for i, p in enumerate(parsed):
            if p["kind"] == "line":
                u0, v0 = p["start"]
                u1, v1 = parsed[(i + 1) % n]["start"]
                seg_edges.append({"curve": "line", "params": [u0, v0, u1, v1]})
            elif p["kind"] == "arc":
                seg_edges.append({"curve": "arc", "params": p["params"]})
            else:
                seg_edges.append({"curve": "bspline", "params": p["params"]})
        wires_out.append({"role": "outer" if is_outer else "inner", "edges": seg_edges})
    if not wires_out:
        return None, "quadric face has no wires"
    n_outer = sum(1 for w in wires_out if w["role"] == "outer")
    if n_outer != 1:
        return None, f"quadric face has {n_outer} outer trim wires (expected exactly 1)"
    return wires_out, None


def _quadric_trim_fills_uv_box(face, uv_bounds) -> bool:
    """True when a quadric face's trim is exactly its parametric-rectangle UV bounding box, so the
    scalar (phiStart, phiSweep, ...) parameters describe it faithfully and no wire block is needed
    (keeping the common cylinder/tube/cone/sphere output byte-identical).

    A single line-bounded wire whose (u, v) polygon area equals the UV box area is the full
    rectangle. This is stricter than every-edge-iso: an L-shaped or notched face has all-iso edges
    too but a smaller area, and a face with a hole has more than one wire - both must take the
    general wire path instead of being silently filled in."""
    umin, umax, vmin, vmax = uv_bounds
    box_area = abs((umax - umin) * (vmax - vmin))
    if box_area <= _EXTRACT_TOL:
        return False
    wires = list(_face_wire_edges(face))
    if len(wires) != 1:
        return False
    _wire, _is_outer, edges = wires[0]
    points = []
    for edge, _start_vertex in edges:
        curve2d, first, last = BRep_Tool.CurveOnSurface(edge, face)
        if curve2d is None or Geom2dAdaptor_Curve(curve2d).GetType() != GeomAbs_Line:
            return False
        param = last if edge.Orientation() == TopAbs_REVERSED else first
        p = curve2d.Value(param)
        points.append((p.X(), p.Y()))
    area = 0.0
    n = len(points)
    for i in range(n):
        u0, v0 = points[i]
        u1, v1 = points[(i + 1) % n]
        area += u0 * v1 - u1 * v0
    return abs(0.5 * area - box_area) <= 1e-6 * box_area


def extract_cylindrical_face(face, scale_to_cm: float) -> Tuple[Optional[dict], Optional[str]]:
    """Cylindrical face -> a 'cylinder' surface record.

    OCC parametrizes the cylinder as U = azimuth, V = height along the axis. A single
    iso-parametric wire is emitted as the scalar parametric rectangle (byte-identical to before);
    any other supported trim (a hole/window, or a non-iso straight-line boundary) is emitted as a
    line trim wire in the (phi[rad], h[cm]) domain. Non-line pcurve boundaries force a fallback."""
    adaptor = BRepAdaptor_Surface(face)
    if adaptor.GetType() != GeomAbs_Cylinder:
        return None, "not a cylinder"
    umin, umax, vmin, vmax = breptools.UVBounds(face)
    cyl = adaptor.Cylinder()
    ax3 = cyl.Position()
    s = scale_to_cm
    center = _xyz(ax3.Location(), s)
    axis = _xyz(cyl.Axis().Direction())
    ref_u = _xyz(ax3.XDirection())
    radius = cyl.Radius() * s
    height_min, height_max = vmin * s, vmax * s
    if height_max - height_min <= _EXTRACT_TOL:
        return None, "cylindrical face has a degenerate height range"
    phi_start, phi_sweep = _quadric_phi_range(ax3, umin, umax)
    inner_wall = face.Orientation() == TopAbs_REVERSED
    params = list(center) + list(axis) + list(ref_u) + [radius, height_min, height_max, phi_start, phi_sweep]
    record = {"type": "cylinder", "inner_wall": inner_wall, "params": params}
    if _quadric_trim_fills_uv_box(face, (umin, umax, vmin, vmax)):
        return record, None  # trim is exactly the parametric rectangle: the scalar params suffice
    phi_of_u = (lambda u: u) if ax3.Direct() else (lambda u: -u)
    # affine (u, v) -> (phi[rad], h[cm]); OCC V is the height along the axis
    wires, reason = _quadric_trim_wire(face, lambda u, v: (phi_of_u(u), v * s))
    if wires is None:
        return None, reason
    record["wires"] = wires
    return record, None


def extract_conical_face(face, scale_to_cm: float) -> Tuple[Optional[dict], Optional[str]]:
    """Conical face trimmed to a parametric rectangle -> a 'cone' surface record.

    OCC parametrizes the cone as U = azimuth, V = distance along the ruling line, with
    radius r(v) = RefRadius + v sin(alpha) and axial position h(v) = v cos(alpha). The two
    V bounds give the (height, radius) endpoints; they are ordered so heightMin < heightMax."""
    adaptor = BRepAdaptor_Surface(face)
    if adaptor.GetType() != GeomAbs_Cone:
        return None, "not a cone"
    umin, umax, vmin, vmax = breptools.UVBounds(face)
    cone = adaptor.Cone()
    ax3 = cone.Position()
    s = scale_to_cm
    half = cone.SemiAngle()
    ref_radius = cone.RefRadius()
    cos_a, sin_a = math.cos(half), math.sin(half)
    h_lo, r_lo = vmin * cos_a, ref_radius + vmin * sin_a
    h_hi, r_hi = vmax * cos_a, ref_radius + vmax * sin_a
    if h_lo > h_hi:
        h_lo, h_hi, r_lo, r_hi = h_hi, h_lo, r_hi, r_lo
    if r_lo < -_EXTRACT_TOL or r_hi < -_EXTRACT_TOL:
        return None, "conical trim produces a negative radius"
    r_lo, r_hi = max(0.0, r_lo), max(0.0, r_hi)
    if max(r_lo, r_hi) <= _EXTRACT_TOL:
        return None, "conical face has degenerate radii"
    if (h_hi - h_lo) * s <= _EXTRACT_TOL:
        return None, "conical face has a degenerate height range"
    center = _xyz(ax3.Location(), s)
    axis = _xyz(cone.Axis().Direction())
    ref_u = _xyz(ax3.XDirection())
    phi_start, phi_sweep = _quadric_phi_range(ax3, umin, umax)
    inner_wall = face.Orientation() == TopAbs_REVERSED
    params = (list(center) + list(axis) + list(ref_u) +
              [r_lo * s, r_hi * s, h_lo * s, h_hi * s, phi_start, phi_sweep])
    record = {"type": "cone", "inner_wall": inner_wall, "params": params}
    if _quadric_trim_fills_uv_box(face, (umin, umax, vmin, vmax)):
        return record, None  # trim is exactly the parametric rectangle: the scalar params suffice
    # OCC V (ruling-line distance) maps to the C++ axial height h = v cos(alpha), in cm
    phi_of_u = (lambda u: u) if ax3.Direct() else (lambda u: -u)
    wires, reason = _quadric_trim_wire(face, lambda u, v: (phi_of_u(u), v * cos_a * s))
    if wires is None:
        return None, reason
    record["wires"] = wires
    return record, None


def extract_spherical_face(face, scale_to_cm: float) -> Tuple[Optional[dict], Optional[str]]:
    """Spherical face trimmed to a parametric rectangle -> a 'sphere' surface record.

    OCC parametrizes the sphere as U = azimuth, V = latitude in [-pi/2, pi/2]; the C++ polar
    angle measured from +polarAxis is theta = pi/2 - v, so the V range maps to [thetaMin,
    thetaMax] = [pi/2 - vmax, pi/2 - vmin]."""
    adaptor = BRepAdaptor_Surface(face)
    if adaptor.GetType() != GeomAbs_Sphere:
        return None, "not a sphere"
    umin, umax, vmin, vmax = breptools.UVBounds(face)
    sph = adaptor.Sphere()
    ax3 = sph.Position()
    s = scale_to_cm
    center = _xyz(ax3.Location(), s)
    polar_axis = _xyz(ax3.Direction())
    ref_u = _xyz(ax3.XDirection())
    radius = sph.Radius() * s
    theta_min = 0.5 * math.pi - vmax
    theta_max = 0.5 * math.pi - vmin
    if theta_max - theta_min <= _EXTRACT_TOL:
        return None, "spherical face has a degenerate polar range"
    phi_start, phi_sweep = _quadric_phi_range(ax3, umin, umax)
    params = list(center) + list(polar_axis) + list(ref_u) + [radius, theta_min, theta_max, phi_start, phi_sweep]
    inner_wall = face.Orientation() == TopAbs_REVERSED
    record = {"type": "sphere", "inner_wall": inner_wall, "params": params}
    if _quadric_trim_fills_uv_box(face, (umin, umax, vmin, vmax)):
        return record, None  # trim is exactly the parametric rectangle: the scalar params suffice
    # OCC V (latitude) maps to the C++ polar angle theta = pi/2 - v (rad)
    phi_of_u = (lambda u: u) if ax3.Direct() else (lambda u: -u)
    wires, reason = _quadric_trim_wire(face, lambda u, v: (phi_of_u(u), 0.5 * math.pi - v))
    if wires is None:
        return None, reason
    record["wires"] = wires
    return record, None


def extract_toroidal_face(face, scale_to_cm: float) -> Tuple[Optional[dict], Optional[str]]:
    """Toroidal face trimmed to a parametric rectangle -> a 'torus' surface record.

    OCC parametrizes the torus as U = phiRing (around the main axis) and V = phiTube (around the
    tube, measured from the outer equator towards the +axis pole), both in [0, 2pi] and with the
    same standard formula the C++ TorusBoundedSurface uses, so U maps to the ring phi (mirrored
    for a left-handed ax3, as for the other quadrics) and V maps directly to the tube phi. A
    single iso-parametric wire (incl. the full-torus seam rectangle) stays on the scalar path;
    any other supported trim becomes a general trim wire in the (phiRing[rad], phiTube[rad])
    domain."""
    adaptor = BRepAdaptor_Surface(face)
    if adaptor.GetType() != GeomAbs_Torus:
        return None, "not a torus"
    umin, umax, vmin, vmax = breptools.UVBounds(face)
    tor = adaptor.Torus()
    ax3 = tor.Position()
    s = scale_to_cm
    major_radius = tor.MajorRadius() * s
    minor_radius = tor.MinorRadius() * s
    if minor_radius <= _EXTRACT_TOL or major_radius <= _EXTRACT_TOL:
        return None, "toroidal face has degenerate radii"
    center = _xyz(ax3.Location(), s)
    axis = _xyz(ax3.Direction())
    ref_u = _xyz(ax3.XDirection())
    phi_start, phi_sweep = _quadric_phi_range(ax3, umin, umax)
    two_pi = 2.0 * math.pi
    tube_sweep = vmax - vmin
    if tube_sweep <= 0.0:
        tube_sweep += two_pi
    tube_sweep = min(tube_sweep, two_pi)
    tube_start = vmin
    inner_wall = face.Orientation() == TopAbs_REVERSED
    params = (list(center) + list(axis) + list(ref_u) +
              [major_radius, minor_radius, phi_start, phi_sweep, tube_start, tube_sweep])
    record = {"type": "torus", "inner_wall": inner_wall, "params": params}
    if _quadric_trim_fills_uv_box(face, (umin, umax, vmin, vmax)):
        return record, None  # trim is exactly the parametric rectangle: the scalar params suffice
    # affine (u, v) -> (phiRing[rad], phiTube[rad]); OCC V is the tube angle, unchanged by the frame
    phi_of_u = (lambda u: u) if ax3.Direct() else (lambda u: -u)
    wires, reason = _quadric_trim_wire(face, lambda u, v: (phi_of_u(u), v))
    if wires is None:
        return None, reason
    record["wires"] = wires
    return record, None


# -------------------------------
# Canonical-form recognition: recover the exact analytic model behind a stored NURBS
# -------------------------------
# The stored STEP surface type describes the exporter, not the geometry: CAD kernels routinely
# write an exact cylinder/cone/sphere/plane as a *rational* B-spline/Bezier patch (exact, not an
# approximation - a rational quadratic/cubic reproduces conics exactly). This is a numeric,
# differential-geometry recognizer (normal-field based, no initial guess/iteration needed) ported
# from the validated standalone prototype `analyze_surface_geometry.py` (see that file and
# BVHSurfaceSolid.md, milestone "Canonical-form recognition"). It is model *selection*, not
# fitting: only a fit at machine precision is accepted, so an "almost cylinder" that is really
# free-form stays free-form. Used as a pre-pass fallback in `extract_surfaces_for_shape` for faces
# whose stored surface type has no direct extractor (bspline/bezier/revolution/extrusion/...).

_RECOGNIZE_TOL_EXACT = 1.e-9


def _sample_surface_for_recognition(adaptor, umin: float, umax: float, vmin: float, vmax: float, n: int = 9):
    """Sample an (n x n) grid over the face's actual trimmed (u, v) box (from `breptools.UVBounds`,
    not the underlying surface's full natural domain). Returns (points, unit normals) in *native*
    (unscaled) CAD length units, or (None, None) if unsampleable."""
    if not all(math.isfinite(x) for x in (umin, umax, vmin, vmax)):
        return None, None
    points, normals = [], []
    for i in range(n):
        u = umin + (umax - umin) * i / (n - 1.0)
        for j in range(n):
            v = vmin + (vmax - vmin) * j / (n - 1.0)
            p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
            try:
                adaptor.D1(u, v, p, du, dv)
            except Exception:
                return None, None
            nrm = _v_cross([du.X(), du.Y(), du.Z()], [dv.X(), dv.Y(), dv.Z()])
            length = math.sqrt(_v_dot(nrm, nrm))
            if length < 1e-14:  # parametric degeneracy (pole/seam): skip this sample
                continue
            points.append([p.X(), p.Y(), p.Z()])
            normals.append([c / length for c in nrm])
    if len(points) < 3 * n:
        return None, None
    return np.array(points), np.array(normals)


def _recognize_analytic_surface(adaptor, uv_bounds) -> Optional[dict]:
    """Model-selection recognizer: plane (3 params) < sphere (4) < cylinder (5) < cone (6) <
    free-form, trying each candidate and keeping the smallest relative residual. Accepts only a
    machine-precision fit (< `_RECOGNIZE_TOL_EXACT`); returns None otherwise (freeform / degenerate
    / unsampleable). All lengths in the returned frame are in *native* (unscaled) CAD units - the
    caller applies `scale_to_cm`.
    """
    umin, umax, vmin, vmax = uv_bounds
    P, N = _sample_surface_for_recognition(adaptor, umin, umax, vmin, vmax)
    if P is None:
        return None
    scale = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0)))
    if scale < 1e-12:
        return None

    # --- plane (3): all unit normals parallel
    dev = float(np.abs(np.abs(N @ N[0]) - 1.0).max())
    if dev < 1e-12:
        normal = N[0] / np.linalg.norm(N[0])
        return {"kind": "plane", "residual": dev, "normal": normal, "point": P[0], "P": P, "N": N}

    best = ("freeform", float("inf"), {})

    # --- sphere (4): normal lines concurrent, P_i = C + r*N_i
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

    # --- cylinder (5): normals coplanar; axis = smallest right singular vector of the normal field
    _, _, Vt = np.linalg.svd(N, full_matrices=False)
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
                origin = cx * e1 + cy * e2  # a point on the axis (axial component is free)
                best = ("cylinder", res, {"axis": axis, "refu": e1, "origin": origin, "radius": R})

    # --- cone (6): N_i . (P_i - A) = 0 is linear in the apex A
    apex, *_ = np.linalg.lstsq(N, np.einsum('ij,ij->i', N, P), rcond=None)
    d = P - apex
    dn = np.linalg.norm(d, axis=1)
    ok = dn > 1e-12
    if ok.sum() > 10:
        u = d[ok] / dn[ok, None]
        res = float(np.abs(np.einsum('ij,ij->i', u, N[ok])).max())
        mean_dir = u.mean(axis=0)
        _, _, Vt2 = np.linalg.svd(u - mean_dir, full_matrices=False)
        ax2 = np.cross(Vt2[0], Vt2[1])
        n2 = np.linalg.norm(ax2)
        if n2 > 1e-12:  # constant half-angle about the ruling axis
            ax2 = ax2 / n2
            if np.dot(mean_dir, ax2) < 0.0:
                ax2 = -ax2
            res = max(res, float(np.abs(u @ ax2).std()))
            if res < best[1]:
                ref = u[0] - np.dot(u[0], ax2) * ax2
                refn = np.linalg.norm(ref)
                if refn > 1e-9:
                    half_angle = float(np.arccos(np.clip(np.abs(u @ ax2), -1.0, 1.0)).mean())
                    best = ("cone", res, {"axis": ax2, "apex": apex, "refu": ref / refn, "half_angle": half_angle})

    kind, res, extra = best
    if res >= _RECOGNIZE_TOL_EXACT:
        return None
    out = {"kind": kind, "residual": res, "P": P, "N": N}
    out.update(extra)
    return out


def _recognized_inner_wall(face, rec) -> Optional[bool]:
    """Decide, by measurement, which side of a RECOGNIZED quadric is outside the solid.

    The sidecar's `inner_wall` flag says whether the patch's outward normal points *towards* the
    axis/centre rather than away from it. For a face whose *stored* surface is an OCC canonical
    quadric, `face.Orientation()` answers that on its own: OCC's canonical quadric normal already
    points away from the axis, so FORWARD means "away" and REVERSED means "towards".

    For a NURBS-encoded quadric it does not. The stored surface is a B-spline whose du x dv may
    point either way -- that is a choice the exporter made when it wrote the parametrisation --
    and FORWARD/REVERSED is relative to *that* normal, not to the axis. The orientation flag
    therefore carries no information about which side is outside, and reading it as if it did
    inverts the outward normal of an arbitrary subset of the recognized faces.

    Measured on ALICE3/`CAD_noETA.stp`: 21 of the 87 NURBS-cylinder faces of `ST0923290_013`,
    `_018` and `_019` came out with an exactly antiparallel outward normal, which
    `DistFromOutside`/`DistFromInside` select hits by and which closure and edge identity are both
    blind to (a global sign cancels in every count they make). See
    `scripts/geometry/Stream_L_ALICE3Defect.md`.

    So: take the face's own outward normal -- OCC's rule, (du x dv) flipped for a REVERSED face --
    at every recognition sample, and compare it with the recognized quadric's canonical outward
    direction there. Returns None when the samples do not decide (the caller then keeps the
    orientation-flag answer).
    """
    samples = rec.get("P")
    normals = rec.get("N")
    if samples is None or normals is None or len(samples) == 0:
        return None
    sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
    kind = rec["kind"]
    votes = 0
    for point, normal in zip(samples, normals):
        outward = np.asarray(normal, dtype=float) * sign
        if kind == "cylinder":
            axis = np.asarray(rec["axis"], dtype=float)
            axis = axis / np.linalg.norm(axis)
            radial = point - rec["origin"]
            radial = radial - np.dot(radial, axis) * axis
        elif kind == "sphere":
            radial = point - rec["centre"]
        elif kind == "cone":
            axis = np.asarray(rec["axis"], dtype=float)
            axis = axis / np.linalg.norm(axis)
            relative = point - rec["apex"]
            radial = relative - np.dot(relative, axis) * axis
            # The cone's outward normal tilts out of the radial direction by the half angle; only
            # its sign relative to the radial direction matters here, and that tilt cannot flip it.
        else:
            return None
        length = np.linalg.norm(radial)
        if length < 1e-12:
            continue  # on the axis: this sample says nothing
        votes += 1 if float(np.dot(outward, radial / length)) > 0.0 else -1
    if votes == 0:
        return None
    return votes < 0


def _arbitrary_orthonormal_frame(axis):
    """One arbitrary orthonormal in-plane vector for an axis with no natural reference direction
    (a full/partial sphere has no preferred polar reference)."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    seed = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = seed - np.dot(seed, axis) * axis
    return e1 / np.linalg.norm(e1)


def _recognized_quadric_wire_block(face, project):
    """Build general line-only trim wires in the recognized (phi, other) domain by sampling each
    3D boundary edge directly (any curve type: line/circle/bspline...) and phi-unwrapping
    *continuously* across all samples of the whole wire loop, in wire-traversal order. Sampling
    works from the edge's own 3D curve, so this does not depend on the *stored* surface's own
    (u, v) parametrization being aligned with the recognized frame at all (unlike the existing
    `_edge_pcurve_is_iso`, which only tests alignment with the stored parametrization).

    Continuous per-sample unwrapping (rather than a coarser vertex-to-vertex or grid unwrap) is
    required for correctness: a single rim/cap edge on a patch whose angular sweep approaches pi
    connects two vertices whose raw delta-phi sits right at the +/-pi unwrap branch cut, which is
    fragile; densely sampling *along* that edge keeps every consecutive step small.

    An accepted edge is exactly a straight segment in (phi, other) - a constant-`other` edge is a
    rim/cap, a constant-`phi` edge is a generator/meridian - so no bspline pole-transform is needed
    (unlike `_quadric_trim_wire`, which reparametrizes curved *stored* pcurves under an affine map;
    here the reparametrization target is always a line by construction once an edge passes the iso
    test). A genuinely non-iso boundary edge (a slanted/curved cut in the recognized frame) is
    rejected - out of scope for this pass.

    A *degenerate* edge (a cone apex or sphere pole: a single point, no 3D curve, `other ~= 0` for
    a cone apex) is common on real CAD (countersinks, chamfers) and is handled specially: `phi` is
    numerically indeterminate there (atan2 of a near-zero radial vector), so rather than sample it,
    the point's phi is taken as the *incoming* wire phi (the last unwrapped sample before it) -
    correctly modelling the point as the transition between the two adjacent generator/meridian
    edges' phi values, with no spurious jump.

    Returns (wires, (phiStart, phiSweep, otherLo, otherHi), None) on success, or
    (None, None, reason) on failure. The window is the outer wire's own sample extent - a
    conservative bound by construction, mirroring the scalar-parameter "conservative window" role
    documented in the sidecar format for a wire-trimmed quadric.
    """
    wires_edges = list(_face_wire_edges(face))
    if not wires_edges:
        return None, None, "recognized quadric face has no wires"

    n_samples = 9
    per_wire = []  # (is_outer, [ [ (phi_raw, other), ... ] or None (degenerate) per edge ], [ other_at_degenerate_vertex or None ])
    all_other = []
    for _wire, is_outer, edges in wires_edges:
        if len(edges) < 3:
            return None, None, "recognized quadric trim wire has fewer than 3 edges"
        edge_samples = []
        degenerate_other = []
        for edge, start_vertex in edges:
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
                return None, None, "recognized quadric boundary edge has no 3D curve"
            reversed_edge = edge.Orientation() == TopAbs_REVERSED
            samples = []
            for k in range(n_samples):
                tau = k / (n_samples - 1.0)
                t = (1.0 - tau) if reversed_edge else tau
                phi, other = project(curve3d.Value(first + t * (last - first)))
                samples.append((phi, other))
                all_other.append(other)
            edge_samples.append(samples)
        per_wire.append((is_outer, edge_samples, degenerate_other))
    tol_other = 1e-6 * max(1.0, max(all_other) - min(all_other))
    tol_phi = 1e-7

    wires_out: List[dict] = []
    outer_window = None
    for is_outer, edge_samples, degenerate_other in per_wire:
        n = len(edge_samples)
        unwrapped_edges = []
        prev_phi = None
        for samples, deg_other in zip(edge_samples, degenerate_other):
            if samples is None:
                # degenerate point: carry the running phi through unchanged (see docstring)
                if prev_phi is None:
                    prev_phi = 0.0
                unwrapped_edges.append([prev_phi] * n_samples)
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

        starts = []
        all_phi_u, all_other_w = [], []
        for i, samples in enumerate(edge_samples):
            phis_u = unwrapped_edges[i]
            if samples is None:
                others = [degenerate_other[i]] * n_samples
            else:
                others = [o for _p, o in samples]
            all_phi_u.extend(phis_u)
            all_other_w.extend(others)
            is_iso_other = (max(others) - min(others)) <= tol_other
            is_iso_phi = (max(phis_u) - min(phis_u)) <= tol_phi
            if not (is_iso_other or is_iso_phi):
                return None, None, "recognized quadric boundary edge is not axis-aligned in (phi, h/theta)"
            starts.append((phis_u[0], others[0]))
        if is_outer:
            outer_window = (min(all_phi_u), max(all_phi_u) - min(all_phi_u), min(all_other_w), max(all_other_w))

        seg_edges = []
        for i in range(n):
            u0, v0 = starts[i]
            u1, v1 = starts[(i + 1) % n]
            seg_edges.append({"curve": "line", "params": [u0, v0, u1, v1]})
        wires_out.append({"role": "outer" if is_outer else "inner", "edges": seg_edges})

    n_outer = sum(1 for w in wires_out if w["role"] == "outer")
    if n_outer != 1:
        return None, None, f"recognized quadric face has {n_outer} outer trim wires (expected exactly 1)"
    return wires_out, outer_window, None


def recognize_and_extract_face(face, scale_to_cm: float) -> Tuple[Optional[dict], Optional[str]]:
    """Canonical-form recognition pre-pass: for a face whose *stored* surface has no direct
    extractor, try to recognize the exact plane/sphere/cylinder/cone hiding behind it and extract
    the same sidecar record a natively-analytic face of that kind would produce. Returns
    (None, None) when the face is genuinely not recognizable (freeform / unsampleable) - the
    caller should keep its original stored-type fallback reason in that case, not this one's.
    """
    adaptor = BRepAdaptor_Surface(face)
    try:
        uv_bounds = breptools.UVBounds(face)
    except Exception:
        return None, None
    rec = _recognize_analytic_surface(adaptor, uv_bounds)
    if rec is None:
        return None, None
    kind = rec["kind"]
    s = scale_to_cm
    # Which side of the recognized quadric is outside. The orientation flag alone answers this
    # only when the *stored* surface is an OCC canonical quadric; on a NURBS-encoded one it is
    # relative to the exporter's parametrisation and says nothing, so it is measured against the
    # face's own outward normal instead. `_recognized_inner_wall` documents why, and returns None
    # when the samples do not decide -- the orientation flag is then all there is.
    #
    # The plane branch below does not use `inner_wall` for its normal: it builds the frame from
    # the sampled normal itself, so it already carries the measurement.
    inner_wall = face.Orientation() == TopAbs_REVERSED
    measured_inner_wall = _recognized_inner_wall(face, rec)
    if measured_inner_wall is not None:
        inner_wall = measured_inner_wall

    if kind == "plane":
        normal = rec["normal"]
        e1 = _arbitrary_orthonormal_frame(normal)
        outward_sign = -1.0 if inner_wall else 1.0
        e2 = np.cross(normal, e1) * outward_sign  # axisU x axisV must equal the outward normal
        origin_cm = (rec["point"] * s).tolist()
        record, reason = extract_planar_face(face, s, frame_override=(origin_cm, e1.tolist(), e2.tolist()))
        if record is None:
            return None, f"recognized as plane but {reason}"
        record["recognized"] = {"kind": "plane", "residual": rec["residual"]}
        return record, None

    if kind == "cylinder":
        axis = rec["axis"] / np.linalg.norm(rec["axis"])
        refu = rec["refu"] - np.dot(rec["refu"], axis) * axis
        refu = refu / np.linalg.norm(refu)
        e2 = np.cross(axis, refu)
        origin_native = rec["origin"]

        def project(pnt):
            rel = np.array([pnt.X(), pnt.Y(), pnt.Z()]) - origin_native
            phi = math.atan2(np.dot(rel, e2), np.dot(rel, refu))
            return phi, float(np.dot(rel, axis)) * s

        wires, window, reason = _recognized_quadric_wire_block(face, project)
        if wires is None:
            return None, f"recognized as cylinder but {reason}"
        phi_start, phi_sweep, h_lo, h_hi = window
        if phi_sweep <= 0.0 or phi_sweep > 2.0 * math.pi + 1e-9:
            return None, "recognized cylinder trim wraps more than a full turn in phi"
        params = ((origin_native * s).tolist() + axis.tolist() + refu.tolist() +
                  [rec["radius"] * s, h_lo, h_hi, phi_start, phi_sweep])
        record = {"type": "cylinder", "inner_wall": inner_wall, "params": params, "wires": wires,
                  "recognized": {"kind": "cylinder", "residual": rec["residual"]}}
        return record, None

    if kind == "cone":
        axis = rec["axis"] / np.linalg.norm(rec["axis"])
        refu = rec["refu"] - np.dot(rec["refu"], axis) * axis
        refu = refu / np.linalg.norm(refu)
        e2 = np.cross(axis, refu)
        apex_native = rec["apex"]
        tan_half = math.tan(rec["half_angle"])

        def project(pnt):
            rel = np.array([pnt.X(), pnt.Y(), pnt.Z()]) - apex_native
            phi = math.atan2(np.dot(rel, e2), np.dot(rel, refu))
            return phi, float(np.dot(rel, axis)) * s

        wires, window, reason = _recognized_quadric_wire_block(face, project)
        if wires is None:
            return None, f"recognized as cone but {reason}"
        phi_start, phi_sweep, h_lo, h_hi = window
        if phi_sweep <= 0.0 or phi_sweep > 2.0 * math.pi + 1e-9:
            return None, "recognized cone trim wraps more than a full turn in phi"
        h_lo = max(0.0, h_lo)
        h_hi = max(h_lo, h_hi)
        params = ((apex_native * s).tolist() + axis.tolist() + refu.tolist() +
                  [h_lo * tan_half, h_hi * tan_half, h_lo, h_hi, phi_start, phi_sweep])
        record = {"type": "cone", "inner_wall": inner_wall, "params": params, "wires": wires,
                  "recognized": {"kind": "cone", "residual": rec["residual"]}}
        return record, None

    if kind == "sphere":
        centre_native = rec["centre"]
        # A sphere has no natural polar axis; any orthonormal frame is a valid (self-consistent)
        # (phi, theta) parametrization for this face.
        polar_axis = np.array([0.0, 0.0, 1.0])
        refu = _arbitrary_orthonormal_frame(polar_axis)
        e2 = np.cross(polar_axis, refu)

        def project(pnt):
            rel = (np.array([pnt.X(), pnt.Y(), pnt.Z()]) - centre_native) / rec["radius"]
            theta = math.acos(max(-1.0, min(1.0, float(np.dot(rel, polar_axis)))))
            phi = math.atan2(float(np.dot(rel, e2)), float(np.dot(rel, refu)))
            return phi, theta

        wires, window, reason = _recognized_quadric_wire_block(face, project)
        if wires is None:
            return None, f"recognized as sphere but {reason}"
        phi_start, phi_sweep, theta_lo, theta_hi = window
        if phi_sweep <= 0.0 or phi_sweep > 2.0 * math.pi + 1e-9:
            return None, "recognized sphere trim wraps more than a full turn in phi"
        params = ((centre_native * s).tolist() + polar_axis.tolist() + refu.tolist() +
                  [rec["radius"] * s, theta_lo, theta_hi, phi_start, phi_sweep])
        record = {"type": "sphere", "inner_wall": inner_wall, "params": params, "wires": wires,
                  "recognized": {"kind": "sphere", "residual": rec["residual"]}}
        return record, None

    return None, None


# Face extractors dispatched by analytic surface type. Planar faces cover both straight-edged
# polygons and circular disks/annuli; the quadrics and torus carry their trim in the surface
# parameters (or a general wire block for non-rectangular trims).
_FACE_EXTRACTORS = {
    "plane": extract_planar_face,
    "cylinder": extract_cylindrical_face,
    "cone": extract_conical_face,
    "sphere": extract_spherical_face,
    "torus": extract_toroidal_face,
}


def extract_surfaces_for_shape(shape, scale_to_cm: float,
                               recognize_surfaces: bool = True) -> Tuple[Optional[List[dict]], List[str], int]:
    """Attempt to extract every face of a leaf solid into exact sidecar surface records.

    Returns (surfaces, [], nModelEdges) when all faces are supported, or (None, reasons, 0)
    listing why the solid cannot be represented exactly (one reason per unsupported face). A shape
    with no faces is treated as unsupported.

    Every accepted record carries `edge_refs`: the identity of the face's boundary edges within
    this solid (see `face_boundary_edge_refs`). Because a single unsupported face rejects the
    whole solid, an emitted sidecar always describes *all* of the shape's faces, which is what
    makes the edge incidence count on the reader's side a complete statement rather than a
    statement about the subset that converted.

    `recognize_surfaces`: when a face's *stored* surface type has no direct extractor (typically
    bspline/bezier/revolution/extrusion), try the canonical-form recognition pre-pass
    (`recognize_and_extract_face`) before giving up on it - see "Canonical-form recognition"
    above and BVHSurfaceSolid.md.
    """
    surfaces: List[dict] = []
    reasons: List[str] = []
    n_faces = 0
    edge_map, edge_id = build_edge_table(shape)
    for face in TopologyExplorer(shape).faces():
        n_faces += 1
        adaptor = BRepAdaptor_Surface(face)
        surf_type = _SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")
        extractor = _FACE_EXTRACTORS.get(surf_type)
        if extractor is None:
            record, reason = None, f"{surf_type} face extraction not implemented yet"
        else:
            record, reason = extractor(face, scale_to_cm)
        if record is None and recognize_surfaces:
            rec_record, rec_reason = recognize_and_extract_face(face, scale_to_cm)
            if rec_record is not None:
                record, reason = rec_record, None
            elif rec_reason is not None:
                reason = f"{reason}; recognition attempted: {rec_reason}"
        if record is None:
            reasons.append(reason or f"{surf_type} face not supported")
        else:
            record["edge_refs"] = face_boundary_edge_refs(face, edge_id,
                                                         anchored=bool(record.get("wires")))
            surfaces.append(record)
    if n_faces == 0:
        return None, ["shape has no faces"], 0
    if reasons:
        return None, reasons, 0
    return surfaces, [], edge_map.Size()


# -------------------------------
# BOM / material mapping
# -------------------------------

@dataclass(frozen=True)
class BomEntry:
    part_number: str
    revision: str
    name: str
    mass_value: float  # as in CSV
    material: str

    @property
    def part_number_key(self) -> str:
        return (self.part_number or "").strip()

    @property
    def name_key(self) -> str:
        return (self.name or "").strip()


def _to_float(s: str) -> Optional[float]:
    try:
        if s is None:
            return None
        s = str(s).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


def read_bom_csv(csv_path: str) -> List[BomEntry]:
    """
    Reads a BOM CSV in the format provided by design team.

    We look for rows whose first column is 'CAD' and second is 'Mechanical/Part'.
    Columns (0-based):
      0 CAD
      1 type
      2 part number
      3 revision
      4 name/description
      5 mass
      6 material
    """
    entries: List[BomEntry] = []
    with open(csv_path, newline="", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if len(row) < 7:
                continue
            if row[0].strip() != "CAD":
                continue
            if row[1].strip() != "Mechanical/Part":
                continue

            part_no = (row[2] or "").strip()
            rev = (row[3] or "").strip()
            name = (row[4] or "").strip()
            mass = _to_float(row[5])
            mat = (row[6] or "").strip()

            if not (part_no or name):
                continue
            if mass is None:
                mass = float("nan")
            if not mat:
                mat = "Default"

            entries.append(BomEntry(part_no, rev, name, float(mass), mat))
    return entries



def normalize_material_name(mat: str) -> str:
    """
    Normalizes a BOM material string for matching / caching.

    Note: We keep the *original* string for ROOT object names; this is only used
    internally for robust matching and dictionary keys.
    """
    mat = (mat or "Default").strip()
    mat = re.sub(r"\s+", " ", mat)
    return mat


def _norm_tokens(s: str) -> List[str]:
    s = (s or "").lower()
    # common grade/format noise
    s = re.sub(r"\(.*?\)", " ", s)
    s = s.replace("en aw", " ")
    s = s.replace("en-aw", " ")
    s = s.replace("en", " ")
    s = s.replace("aw", " ")
    s = s.replace("_", " ").replace("-", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    if not s:
        return []
    toks = s.split(" ")

    # small synonym normalization
    syn = {
        "alu": "al",
        "aluminium": "aluminum",
        "silicium": "silicon",
        "inox": "stainless",
        "ss": "stainless",
        "cu": "copper",
        "fe": "iron",
        "ptfe": "teflon",
        "ti": "titanium",
        "be": "beryllium",
    }

    # Expand common element symbols to names and vice-versa so that e.g. "G4_Si" can match "silicon".
    elem_alias = {
        "h": "hydrogen", "he": "helium", "c": "carbon", "n": "nitrogen", "o": "oxygen",
        "al": "aluminum", "si": "silicon", "fe": "iron", "cu": "copper", "be": "beryllium",
        "mg": "magnesium", "mn": "manganese", "cr": "chromium", "ni": "nickel", "zn": "zinc",
        "ti": "titanium", "w": "tungsten", "pb": "lead", "sn": "tin",
    }
    name_to_sym = {v: k for k, v in elem_alias.items()}

    out: List[str] = []
    for t in toks:
        t2 = syn.get(t, t)
        out.append(t2)
        if t2 in elem_alias:
            out.append(elem_alias[t2])
        if t2 in name_to_sym:
            out.append(name_to_sym[t2])

    # de-dup while preserving order
    seen = set()
    out2: List[str] = []
    for t in out:
        if t and t not in seen:
            seen.add(t)
            out2.append(t)
    return out2


def _density_score(rho_part: Optional[float], rho_ref: Optional[float]) -> float:
    if rho_part is None or rho_ref is None or not (rho_part > 0.0) or not (rho_ref > 0.0):
        return 0.0
    # symmetric score in log-space; 1.0 is perfect match
    d = abs(math.log(rho_ref / rho_part))
    return 1.0 / (1.0 + d)


def _token_score(tokens_a: List[str], tokens_b: List[str]) -> float:
    if not tokens_a or not tokens_b:
        return 0.0
    sa = set(tokens_a)
    sb = set(tokens_b)
    inter = len(sa & sb)
    union = len(sa | sb)
    if union == 0:
        return 0.0
    return inter / union


def load_g4_nist_db(json_path: str) -> Dict[str, dict]:
    """
    Loads a JSON dump created by the 'nist_export_all' tool.
    Returns a dict: nist_name -> material record.
    """
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    mats = data.get("materials", {})
    if not isinstance(mats, dict) or not mats:
        raise RuntimeError(f"G4 NIST DB JSON seems empty or malformed: {json_path}")
    return mats

# Minimal periodic table for parsing custom alloys not present in NIST.
# Values: Z (atomic number), A (g/mol)
_ELEMENT_TABLE = {
    "H": (1, 1.00794),
    "C": (6, 12.0107),
    "N": (7, 14.0067),
    "O": (8, 15.9994),
    "Al": (13, 26.9815385),
    "Si": (14, 28.0855),
    "Fe": (26, 55.845),
    "Cu": (29, 63.546),
    "Be": (4, 9.0121831),
    "Mg": (12, 24.305),
    "Mn": (25, 54.938044),
    "Cr": (24, 51.9961),
    "Ni": (28, 58.6934),
    "Zn": (30, 65.38),
    "Ti": (22, 47.867),
    "W": (74, 183.84),
    "Pb": (82, 207.2),
    "Sn": (50, 118.71),
}


@dataclass
class ResolvedMaterial:
    bom_name: str
    nist_name: Optional[str]          # e.g. "G4_Al"
    score: float
    rho_used_g_cm3: Optional[float]   # density used in ROOT definition
    radlen_cm: Optional[float]
    intlen_cm: Optional[float]
    elements: Optional[List[dict]]    # list of {symbol,Z,A_g_mol,mass_fraction}
    note: str                         # for comments in geom.C (warnings/FIXME)

@dataclass
class MatMatchConfig:
    # Minimum combined score to accept a match.
    min_score: float = 0.35
    # If (best - second_best) < ambiguity_delta, treat as ambiguous/unresolved.
    ambiguity_delta: float = 0.05
    # Weights for the combined score = w_token * token_score + w_density * density_score
    w_token: float = 0.75
    w_density: float = 0.25
    # Optional hard filter on density proximity (in log-space). If <=0, disabled.
    # Example: max_log_density_diff=0.8 means accept within exp(0.8)~2.2x in either direction.
    max_log_density_diff: float = 0.0
    # Penalize compound matches (oxide/dioxide/carbide/...) when BOM doesn't mention those tokens.
    compound_penalty: float = 0.25


def resolve_bom_material(
    bom_material: str,
    rho_part_g_cm3: Optional[float],
    g4db: Optional[Dict[str, dict]],
    cfg: MatMatchConfig,
) -> ResolvedMaterial:
    """
    Resolves an arbitrary BOM material string to a Geant4 NIST material name using:
      - exact key match (BOM already uses e.g. "G4_Al")
      - token overlap scoring on names
      - density proximity scoring (if rho_part_g_cm3 available)

    If unresolved/ambiguous, tries to parse element symbols from the BOM string (e.g. "Cu Be")
    and emits a placeholder mixture (equal mass fractions) annotated with FIXME.
    """
    raw_bom_material = (bom_material or "").strip()
    bom_material = normalize_material_name(bom_material)

    if not g4db:
        return ResolvedMaterial(
            bom_name=bom_material,
            nist_name=None,
            score=0.0,
            rho_used_g_cm3=rho_part_g_cm3,
            radlen_cm=None,
            intlen_cm=None,
            elements=None,
            note="FIXME: No Geant4 NIST DB provided; using dummy material.",
        )

    # Trivial: BOM already provides an exact Geant4 material key
    if bom_material in g4db:
        rec = g4db[bom_material]
        rho_ref = rec.get("density_g_cm3")
        # Use NIST density for emission; CAD-derived density is used only for matching.
        rho_used = rho_ref

        rad = rec.get("radlen_cm")
        itl = rec.get("intlen_cm")

        return ResolvedMaterial(
            bom_name=bom_material,
            nist_name=bom_material,
            score=1.0,
            rho_used_g_cm3=rho_used,
            radlen_cm=rad,
            intlen_cm=itl,
            elements=rec.get("elements", []),
            note="Resolved by exact Geant4 NIST name from BOM.",
        )

    bom_toks = _norm_tokens(bom_material)
    if not bom_toks:
        return ResolvedMaterial(
            bom_name=bom_material,
            nist_name=None,
            score=0.0,
            rho_used_g_cm3=rho_part_g_cm3,
            radlen_cm=None,
            intlen_cm=None,
            elements=None,
            note="FIXME: Empty/unknown BOM material string; using dummy material.",
        )

    def _build_custom_from_elements(note_prefix: str) -> Optional[ResolvedMaterial]:
        s = raw_bom_material
        if not s:
            return None

        symbols = set(re.findall(r"\b([A-Z][a-z]?)\b", s))
        name_to_symbol = {
            "aluminum": "Al", "aluminium": "Al", "silicon": "Si", "iron": "Fe", "copper": "Cu",
            "beryllium": "Be", "magnesium": "Mg", "manganese": "Mn", "chromium": "Cr", "nickel": "Ni",
            "zinc": "Zn", "titanium": "Ti", "tungsten": "W", "lead": "Pb", "tin": "Sn",
        }
        for t in bom_toks:
            if t in name_to_symbol:
                symbols.add(name_to_symbol[t])

        symbols = [sym for sym in sorted(symbols) if sym in _ELEMENT_TABLE]
        if not symbols:
            return None

        frac = 1.0 / float(len(symbols))
        elems: List[dict] = []
        for sym in symbols:
            Z, A = _ELEMENT_TABLE[sym]
            elems.append({"symbol": sym, "Z": Z, "A_g_mol": A, "mass_fraction": frac})

        return ResolvedMaterial(
            bom_name=bom_material,
            nist_name=None,
            score=0.0,
            rho_used_g_cm3=rho_part_g_cm3,
            radlen_cm=None,
            intlen_cm=None,
            elements=elems,
            note=f"FIXME: {note_prefix} No suitable Geant4 NIST material. Emitting placeholder mixture from parsed elements {symbols} with equal mass fractions; please adjust fractions/material.",
        )

    best = (None, -1.0, 0.0, 0.0)   # (nist_name, score, dens_score, token_score)
    second = (None, -1.0, 0.0, 0.0)

    bom_has_compound = any(t in bom_toks for t in (
        "oxide", "dioxide", "carbide", "nitride", "fluoride", "chloride",
        "sulfate", "phosphate", "glass", "dioxyde"
    ))

    for nist_name, rec in g4db.items():
        nist_toks = _norm_tokens(nist_name)
        ts = _token_score(bom_toks, nist_toks)
        if ts <= 0.0:
            continue

        ds = _density_score(rho_part_g_cm3, rec.get("density_g_cm3"))

        # Optional hard density filter
        if cfg.max_log_density_diff and cfg.max_log_density_diff > 0.0 and rho_part_g_cm3 and rec.get("density_g_cm3"):
            try:
                if abs(math.log(float(rec.get("density_g_cm3")) / float(rho_part_g_cm3))) > cfg.max_log_density_diff:
                    continue
            except Exception:
                pass

        nist_has_compound = any(t in nist_toks for t in (
            "oxide", "dioxide", "carbide", "nitride", "fluoride", "chloride",
            "sulfate", "phosphate", "glass", "dioxyde"
        ))
        compound_pen = cfg.compound_penalty if (nist_has_compound and not bom_has_compound) else 0.0

        score = cfg.w_token * ts + cfg.w_density * ds - compound_pen

        if score > best[1]:
            second = best
            best = (nist_name, score, ds, ts)
        elif score > second[1]:
            second = (nist_name, score, ds, ts)

    nist_best, score_best, ds_best, ts_best = best
    nist_second, score_second, _, _ = second

    if nist_best is None or score_best < cfg.min_score:
        custom = _build_custom_from_elements("Could not resolve with enough confidence.")
        if custom is not None:
            return custom
        return ResolvedMaterial(
            bom_name=bom_material,
            nist_name=None,
            score=float(score_best if score_best > 0 else 0.0),
            rho_used_g_cm3=rho_part_g_cm3,
            radlen_cm=None,
            intlen_cm=None,
            elements=None,
            note="FIXME: Could not resolve BOM material to a Geant4 NIST material with enough confidence; using dummy material.",
        )

    if score_second > 0 and (score_best - score_second) < cfg.ambiguity_delta:
        custom = _build_custom_from_elements(
            f"Ambiguous material match (best '{nist_best}' score={score_best:.3f}, second '{nist_second}' score={score_second:.3f})."
        )
        if custom is not None:
            return custom
        return ResolvedMaterial(
            bom_name=bom_material,
            nist_name=None,
            score=float(score_best),
            rho_used_g_cm3=rho_part_g_cm3,
            radlen_cm=None,
            intlen_cm=None,
            elements=None,
            note=f"FIXME: Ambiguous material match (best '{nist_best}' score={score_best:.3f}, second '{nist_second}' score={score_second:.3f}); using dummy material.",
        )

    rec = g4db[nist_best]
    rho_ref = rec.get("density_g_cm3")
    # Use NIST density for emission; CAD-derived density is used only for matching.
    rho_used = rho_ref

    rad = rec.get("radlen_cm")
    itl = rec.get("intlen_cm")

    return ResolvedMaterial(
        bom_name=bom_material,
        nist_name=nist_best,
        score=float(score_best),
        rho_used_g_cm3=rho_used,
        radlen_cm=rad,
        intlen_cm=itl,
        elements=rec.get("elements", []),
        note=f"Resolved to '{nist_best}' (token={ts_best:.3f}, density={ds_best:.3f}, score={score_best:.3f}).",
    )


def build_volume_to_material_map(
    bom_entries: List[BomEntry],
    def_names: Dict[str, str],
) -> Dict[str, BomEntry]:
    """
    Builds a mapping def_lid -> BomEntry by matching the XCAF display name to:
      - exact part_number match
      - exact description/name match
      - substring match on part_number within the XCAF name

    This is heuristic; if nothing matches we keep no assignment for that volume.
    """
    # lookup tables
    by_part: Dict[str, BomEntry] = {}
    by_name: Dict[str, BomEntry] = {}
    for e in bom_entries:
        if e.part_number_key:
            by_part[e.part_number_key] = e
        if e.name_key and e.name_key not in by_name:
            by_name[e.name_key] = e

    out: Dict[str, BomEntry] = {}
    for lid, disp in def_names.items():
        key = (disp or "").strip()
        if not key:
            continue

        # 1) exact part number
        if key in by_part:
            out[lid] = by_part[key]
            continue
        # 2) exact name/description
        if key in by_name:
            out[lid] = by_name[key]
            continue
        # 3) substring match on any part number
        for pn, e in by_part.items():
            if pn and pn in key:
                out[lid] = e
                break
    return out


# -------------------------------
# C++ emission helpers
# -------------------------------

def trsf_to_tgeo(trsf: gp_Trsf, name: str, scale_to_cm: float) -> str:
    m = trsf.GetRotation().GetMatrix()
    t = trsf.TranslationPart()
    return f"""
  Double_t {name}_m[9] = {{
    {m.Value(1,1)}, {m.Value(1,2)}, {m.Value(1,3)},
    {m.Value(2,1)}, {m.Value(2,2)}, {m.Value(2,3)},
    {m.Value(3,1)}, {m.Value(3,2)}, {m.Value(3,3)}
  }};
  TGeoRotation *{name}_rot = new TGeoRotation();
  {name}_rot->SetMatrix({name}_m);
  TGeoCombiTrans *{name} = new TGeoCombiTrans({t.X()*scale_to_cm}, {t.Y()*scale_to_cm}, {t.Z()*scale_to_cm}, {name}_rot);
"""


def emit_cpp_prelude(exact_surfaces: bool = False, csg_shapes: bool = False) -> str:
    prelude = """#include <TGeoManager.h>
#include <TFile.h>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>

static void LoadFacets(const std::string& file, TGeoTessellated* solid, bool check=false)
{
  std::ifstream in(file, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open facet file: " + file);

  uint32_t nTri = 0;
  in.read(reinterpret_cast<char*>(&nTri), sizeof(nTri));
  if (!in) throw std::runtime_error("Bad facet header in: " + file);

  for (uint32_t i=0;i<nTri;i++) {
    float v[9];
    in.read(reinterpret_cast<char*>(v), sizeof(v));
    if (!in) throw std::runtime_error("Unexpected EOF in: " + file);

    solid->AddFacet(TGeoTessellated::Vertex_t(v[0],v[1],v[2]),
                    TGeoTessellated::Vertex_t(v[3],v[4],v[5]),
                    TGeoTessellated::Vertex_t(v[6],v[7],v[8]));
  }
  solid->CloseShape(check, true);
}
"""
    if csg_shapes:
        prelude += import_csg_hook().CPP_LOADER
    if not exact_surfaces:
        return prelude

    # Exact-surface solids need the O2 DetectorsBase library. The macro stays loadable in
    # ROOT interpreted mode: O2BVHSurfaceSolid.h is part of the ROOT dictionary module so
    # it can be included textually, while O2SurfaceSolidIO.h is not -- its single free
    # function is declared by prototype instead (the symbol resolves from
    # libO2DetectorsBase). Do not export ROOT_INCLUDE_PATH for this; R__ADD_INCLUDE_PATH
    # keeps ROOT's C++ modules intact.
    prelude += """
// --- exact-surface solid support (requires the ALICE O2 environment) ---
R__ADD_INCLUDE_PATH($O2_ROOT/include)
R__LOAD_LIBRARY(libO2DetectorsBase)
#include "DetectorsBase/O2BVHSurfaceSolid.h"
// O2SurfaceSolidIO.h is not part of the ROOT dictionary module; declare the loader
// prototype directly (the symbol resolves from libO2DetectorsBase).
namespace o2
{
namespace base
{
bool LoadSurfaceSolid(const std::string& file, O2BVHSurfaceSolid& solid);
} // namespace base
} // namespace o2

static void LoadSurfaces(const std::string& file, o2::base::O2BVHSurfaceSolid* solid, bool check=false)
{
  if (!o2::base::LoadSurfaceSolid(file, *solid)) {
    throw std::runtime_error("Cannot load surface sidecar: " + file);
  }
  solid->CloseShape(check);
  if (check && (!solid->IsClosed() || !solid->IsOrientationConsistent())) {
    throw std::runtime_error("Surface solid not closed/orientation-consistent: " + file);
  }
}
"""
    return prelude


def emit_materials_cpp(
    used_materials: Dict[str, ResolvedMaterial],
    # key: BOM material string as used in CSV after normalization
) -> Tuple[str, Dict[str, str]]:
    """
    Emits C++ code defining TGeoMaterial/TGeoMixture + TGeoMedium for all used materials.

    - If a material resolved to a Geant4 NIST entry, emit a physically correct mixture
      (element mass fractions) and set RadLen/IntLen (from Geant4) when available.
    - If unresolved/ambiguous, emit a dummy material and annotate with FIXME comments.
    """
    cpp: List[str] = []
    cpp.append("  // Default material/medium (placeholder; can be replaced later)")
    cpp.append("  TGeoMaterial *mat_Default = new TGeoMaterial(\"Default\", 0., 0., 0.);")
    cpp.append("  TGeoMedium   *med_Default = new TGeoMedium(\"Default\", 1, mat_Default);")
    cpp.append("")

    emitted_el: Dict[str, str] = {}

    def _emit_element(el: dict) -> str:
        sym = el.get("symbol", "X")
        Z = int(el.get("Z", 0))
        A = float(el.get("A_g_mol", 0.0))
        if sym in emitted_el:
            return emitted_el[sym]
        safe = sanitize_cpp_name(sym)
        var = f"el_{safe}"
        cpp.append(f"  TGeoElement *{var} = new TGeoElement(\"{sym}\", \"{sym}\", {Z}, {A:.10g});")
        emitted_el[sym] = var
        return var

    medium_var: Dict[str, str] = {"Default": "med_Default"}
    next_id = 2

    for bom_mat in sorted(used_materials.keys(), key=lambda s: s.lower()):
        rm = used_materials[bom_mat]
        safe = sanitize_cpp_name(bom_mat)
        base = safe
        k = 2
        while f"med_{safe}" in medium_var.values():
            safe = f"{base}_{k}"
            k += 1

        rho = rm.rho_used_g_cm3 if (rm.rho_used_g_cm3 and rm.rho_used_g_cm3 > 0.0) else 0.0

        cpp.append(f"  // BOM material: {rm.bom_name}")
        cpp.append(f"  // {rm.note}")

        if rm.elements:
            elems = rm.elements
            if len(elems) == 1 and abs(float(elems[0].get('mass_fraction', 1.0)) - 1.0) < 1e-6:
                el = elems[0]
                A = float(el.get("A_g_mol", 0.0))
                Z = float(el.get("Z", 0))
                cpp.append(f"  TGeoMaterial *mat_{safe} = new TGeoMaterial(\"{bom_mat}\", {A:.10g}, {Z:.10g}, {rho:.10g});")
            else:
                cpp.append(f"  TGeoMixture  *mat_{safe} = new TGeoMixture(\"{bom_mat}\", {len(elems)}, {rho:.10g});")
                for el in elems:
                    elvar = _emit_element(el)
                    w = float(el.get("mass_fraction", 0.0))
                    cpp.append(f"  mat_{safe}->AddElement({elvar}, {w:.10g});")

            if rm.radlen_cm is not None and rm.intlen_cm is not None:
                cpp.append(f"  mat_{safe}->SetRadLen({float(rm.radlen_cm):.10g}, {float(rm.intlen_cm):.10g});")
            elif rm.radlen_cm is not None:
                cpp.append(f"  mat_{safe}->SetRadLen({float(rm.radlen_cm):.10g});")
        else:
            cpp.append("  // FIXME: Unresolved material. Replace with a proper TGeoMaterial/TGeoMixture.")
            cpp.append(f"  TGeoMaterial *mat_{safe} = new TGeoMaterial(\"{bom_mat}\", 0., 0., {rho:.10g});")

        cpp.append(f"  TGeoMedium   *med_{safe} = new TGeoMedium(\"{bom_mat}\", {next_id}, mat_{safe});")
        cpp.append("")
        medium_var[bom_mat] = f"med_{safe}"
        next_id += 1

    return "\n".join(cpp), medium_var




def emit_tessellated_cpp(lid: str, vol_display_name: str, facet_abspath: str, ntriangles: int, medium_var: str) -> str:
    safe = sanitize_cpp_name(lid)
    shape_name = vol_display_name if vol_display_name else lid

    if ntriangles <= 0:
        out = []
        out.append(f'  TGeoBBox *solid_{safe} = new TGeoBBox("{shape_name}", 0.001, 0.001, 0.001);')
        out.append(f'  TGeoVolume *vol_{safe} = new TGeoVolume("{shape_name}", solid_{safe}, {medium_var});')
        return "\n".join(out)

    out = []
    out.append(f'  TGeoTessellated *solid_{safe} = new TGeoTessellated("{shape_name}", {ntriangles});')
    out.append(f'  LoadFacets("{facet_abspath}", solid_{safe}, check);')
    out.append(f'  TGeoVolume *vol_{safe} = new TGeoVolume("{shape_name}", solid_{safe}, {medium_var});')
    return "\n".join(out)


def emit_surface_solid_cpp(lid: str, vol_display_name: str, surface_abspath: str, medium_var: str) -> str:
    """Exact-surface counterpart of emit_tessellated_cpp: the volume gets an
    O2BVHSurfaceSolid filled from a surface sidecar file (see BVHSurfaceSolid.md,
    "Surface sidecar format"). Requires emit_cpp_prelude(exact_surfaces=True)."""
    safe = sanitize_cpp_name(lid)
    shape_name = vol_display_name if vol_display_name else lid

    out = []
    out.append(f'  auto *solid_{safe} = new o2::base::O2BVHSurfaceSolid("{shape_name}");')
    out.append(f'  LoadSurfaces("{surface_abspath}", solid_{safe}, check);')
    out.append(f'  TGeoVolume *vol_{safe} = new TGeoVolume("{shape_name}", solid_{safe}, {medium_var});')
    return "\n".join(out)


def emit_assembly_cpp(lid: str, asm_display_name: str) -> str:
    safe = sanitize_cpp_name(lid)
    name = asm_display_name if asm_display_name else lid
    return f'  TGeoVolumeAssembly *asm_{safe} = new TGeoVolumeAssembly("{name}");'


# -------------------------------
# CAD clipping helpers
# -------------------------------

def make_clip_box_shape(clip_box: ClipBox):
    return BRepPrimAPI_MakeBox(
        gp_Pnt(clip_box.xmin, clip_box.ymin, clip_box.zmin),
        gp_Pnt(clip_box.xmax, clip_box.ymax, clip_box.zmax),
    ).Shape()


def _compose_trsf(parent_to_world: gp_Trsf, local_to_parent: gp_Trsf) -> gp_Trsf:
    return parent_to_world.Multiplied(local_to_parent)


def _shape_is_empty(shape) -> bool:
    if shape is None:
        return True
    try:
        if shape.IsNull():
            return True
    except Exception:
        pass
    try:
        for _ in TopologyExplorer(shape).faces():
            return False
        return True
    except Exception:
        return False


def _transformed_bbox(shape, trsf: gp_Trsf) -> Optional[Tuple[float, float, float, float, float, float]]:
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    try:
        xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    except Exception:
        return None

    points = []
    for x in (xmin, xmax):
        for y in (ymin, ymax):
            for z in (zmin, zmax):
                p = gp_Pnt(x, y, z)
                p.Transform(trsf)
                points.append((p.X(), p.Y(), p.Z()))

    return (
        min(p[0] for p in points),
        min(p[1] for p in points),
        min(p[2] for p in points),
        max(p[0] for p in points),
        max(p[1] for p in points),
        max(p[2] for p in points),
    )


def _bbox_outside_clip_box(bbox: Tuple[float, float, float, float, float, float], clip_box: ClipBox) -> bool:
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    return (
        xmax < clip_box.xmin or xmin > clip_box.xmax or
        ymax < clip_box.ymin or ymin > clip_box.ymax or
        zmax < clip_box.zmin or zmin > clip_box.zmax
    )


def _bbox_inside_clip_box(bbox: Tuple[float, float, float, float, float, float], clip_box: ClipBox) -> bool:
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    return (
        xmin >= clip_box.xmin and xmax <= clip_box.xmax and
        ymin >= clip_box.ymin and ymax <= clip_box.ymax and
        zmin >= clip_box.zmin and zmax <= clip_box.zmax
    )


def _classify_shape_against_clip_box(shape, clip_box: ClipBox, local_to_world: gp_Trsf) -> Optional[str]:
    world_bbox = _transformed_bbox(shape, local_to_world)
    if world_bbox is None:
        return None
    if _bbox_outside_clip_box(world_bbox, clip_box):
        return "outside"
    if _bbox_inside_clip_box(world_bbox, clip_box):
        return "inside"
    return "overlap"


def clip_shape_to_box(shape, clip_box: ClipBox, clip_box_shape, local_to_world: gp_Trsf, lid: str):
    clip_state = _classify_shape_against_clip_box(shape, clip_box, local_to_world)
    if clip_state is None:
        return None
    if clip_state == "outside":
        return None
    if clip_state == "inside":
        return shape

    local_clip = BRepBuilderAPI_Transform(clip_box_shape, local_to_world.Inverted(), True).Shape()
    common = BRepAlgoAPI_Common(shape, local_clip)
    common.Build()
    if not common.IsDone():
        raise RuntimeError(f"Failed to clip CAD shape {lid} against --clip-box")

    clipped = common.Shape()
    if _shape_is_empty(clipped):
        return None
    return clipped


# -------------------------------
# Definition graph extraction
# -------------------------------

logical_volumes: Dict[str, list] = {}     # def_lid -> triangles
def_names: Dict[str, str] = {}           # def_lid -> human display name (may be "")
def_volumes_cm3: Dict[str, float] = {}   # def_lid -> volume in cm^3 (leaf only)
def_shapes: Dict[str, object] = {}       # def_lid -> (possibly clipped) TopoDS shape (leaf only)
assemblies = set()                       # def_lid
placements = []                          # (parent_def_lid, child_def_lid, gp_Trsf local)
top_defs = set()                         # top definition lids
visited_defs = set()                     # expanded defs


def cpp_var_for_def(lid: str) -> str:
    safe = sanitize_cpp_name(lid)
    return f"asm_{safe}" if lid in assemblies else f"vol_{safe}"


def expand_definition(
    def_label: TDF_Label,
    shape_tool,
    meshparam=None,
    scale_to_cm: float = 1.0,
    clip_box: Optional[ClipBox] = None,
    clip_box_shape=None,
    clip_deduplicate: str = "intact",
    name_filter: Optional[NameFilter] = None,
    include_subtree: bool = False,
    world_trsf: Optional[gp_Trsf] = None,
    occ_path: str = "r1",
) -> Optional[str]:
    clip_enabled = clip_box_shape is not None
    if world_trsf is None:
        world_trsf = gp_Trsf()

    def_lid = label_id(def_label)
    nm = label_name(def_label)

    subtree_included = include_subtree
    if name_filter is not None:
        if name_filter.matches_exclude(def_lid, nm):
            return None
        if name_filter.has_include and name_filter.matches_include(def_lid, nm):
            subtree_included = True

    if clip_enabled and clip_box is not None:
        try:
            shape_for_clip = shape_tool.GetShape(def_label)
        except Exception:
            shape_for_clip = None
        if shape_for_clip is not None:
            clip_state = _classify_shape_against_clip_box(shape_for_clip, clip_box, world_trsf)
            if clip_state == "outside":
                return None
            if clip_state == "inside" and clip_deduplicate == "intact":
                return expand_definition(
                    def_label,
                    shape_tool,
                    meshparam=meshparam,
                    scale_to_cm=scale_to_cm,
                    clip_box=None,
                    clip_box_shape=None,
                    clip_deduplicate=clip_deduplicate,
                    name_filter=name_filter,
                    include_subtree=subtree_included,
                )

    def_key = f"{def_lid}@{occ_path}" if clip_enabled else def_lid
    if not clip_enabled and def_lid in visited_defs:
        return def_lid
    if not clip_enabled:
        visited_defs.add(def_lid)

    if nm and def_key not in def_names:
        def_names[def_key] = nm
    elif def_key not in def_names:
        def_names[def_key] = ""

    children = TDF_LabelSequence()
    shape_tool.GetComponents(def_label, children)
    has_children = children.Length() > 0

    if has_children or shape_tool.IsAssembly(def_label):
        assemblies.add(def_key)
        kept_children = 0

        for i in range(children.Length()):
            child = children.Value(i + 1)
            child_occ_path = f"{occ_path}_{i + 1}"
            if shape_tool.IsReference(child):
                referred = TDF_Label()
                shape_tool.GetReferredShape(child, referred)

                loc = shape_tool.GetLocation(child)
                trsf = loc.Transformation()
                if clip_enabled:
                    child_key = expand_definition(
                        referred,
                        shape_tool,
                        meshparam=meshparam,
                        scale_to_cm=scale_to_cm,
                        clip_box=clip_box,
                        clip_box_shape=clip_box_shape,
                        clip_deduplicate=clip_deduplicate,
                        name_filter=name_filter,
                        include_subtree=subtree_included,
                        world_trsf=_compose_trsf(world_trsf, trsf),
                        occ_path=child_occ_path,
                    )
                    if child_key is None:
                        continue
                    placements.append((def_key, child_key, trsf))
                else:
                    child_key = expand_definition(
                        referred,
                        shape_tool,
                        meshparam=meshparam,
                        scale_to_cm=scale_to_cm,
                        clip_deduplicate=clip_deduplicate,
                        name_filter=name_filter,
                        include_subtree=subtree_included,
                    )
                    if child_key is None:
                        continue
                    placements.append((def_key, child_key, trsf))
            else:
                trsf = gp_Trsf()
                if clip_enabled:
                    child_key = expand_definition(
                        child,
                        shape_tool,
                        meshparam=meshparam,
                        scale_to_cm=scale_to_cm,
                        clip_box=clip_box,
                        clip_box_shape=clip_box_shape,
                        clip_deduplicate=clip_deduplicate,
                        name_filter=name_filter,
                        include_subtree=subtree_included,
                        world_trsf=world_trsf,
                        occ_path=child_occ_path,
                    )
                    if child_key is None:
                        continue
                    placements.append((def_key, child_key, trsf))
                else:
                    child_key = expand_definition(
                        child,
                        shape_tool,
                        meshparam=meshparam,
                        scale_to_cm=scale_to_cm,
                        clip_deduplicate=clip_deduplicate,
                        name_filter=name_filter,
                        include_subtree=subtree_included,
                    )
                    if child_key is None:
                        continue
                    placements.append((def_key, child_key, trsf))
            kept_children += 1

        if (clip_enabled or (name_filter is not None and name_filter.has_include)) and kept_children == 0:
            assemblies.discard(def_key)
            return None
        return def_key

    if shape_tool.IsSimpleShape(def_label):
        if name_filter is not None and name_filter.has_include and not subtree_included:
            return None

        if def_key not in logical_volumes:
            shape = shape_tool.GetShape(def_label)

            # store volume (for density estimation)
            try:
                volume_cm3 = volume_cm3_of_shape(shape, scale_to_cm=scale_to_cm)
            except Exception:
                volume_cm3 = 0.0

            if clip_enabled:
                shape = clip_shape_to_box(shape, clip_box, clip_box_shape, world_trsf, def_lid)
                if shape is None:
                    return None

            def_volumes_cm3[def_key] = volume_cm3
            def_shapes[def_key] = shape

            do_meshing = (meshparam is not None) and meshparam.get("do_meshing", None) is True
            logical_volumes[def_key] = triangulate_CAD_solid(shape, meshparam=meshparam, scale_to_cm=scale_to_cm) if do_meshing else triangulate_asbbox(shape, scale_to_cm=scale_to_cm)
        return def_key

    assemblies.add(def_key)
    return def_key


def extract_graph(
    step_path: str,
    meshparam=None,
    scale_to_cm: float = 1.0,
    clip_box: Optional[ClipBox] = None,
    clip_deduplicate: str = "intact",
    name_filter: Optional[NameFilter] = None,
):
    global logical_volumes, def_names, def_volumes_cm3, def_shapes, assemblies, placements, top_defs, visited_defs
    logical_volumes = {}
    def_names = {}
    def_volumes_cm3 = {}
    def_shapes = {}
    assemblies = set()
    placements = []
    top_defs = set()
    visited_defs = set()

    doc, shape_tool = load_step_with_xcaf(step_path)
    clip_box_shape = make_clip_box_shape(clip_box) if clip_box is not None else None

    roots = TDF_LabelSequence()
    shape_tool.GetFreeShapes(roots)

    for i in range(roots.Length()):
        root = roots.Value(i + 1)
        root_occ_path = f"r{i + 1}"
        if shape_tool.IsReference(root):
            ref = TDF_Label()
            shape_tool.GetReferredShape(root, ref)
            top = expand_definition(
                ref,
                shape_tool,
                meshparam=meshparam,
                scale_to_cm=scale_to_cm,
                clip_box=clip_box,
                clip_box_shape=clip_box_shape,
                clip_deduplicate=clip_deduplicate,
                name_filter=name_filter,
                occ_path=root_occ_path,
            )
        else:
            top = expand_definition(
                root,
                shape_tool,
                meshparam=meshparam,
                scale_to_cm=scale_to_cm,
                clip_box=clip_box,
                clip_box_shape=clip_box_shape,
                clip_deduplicate=clip_deduplicate,
                name_filter=name_filter,
                occ_path=root_occ_path,
            )
        if top is not None:
            top_defs.add(top)

    return doc, shape_tool


# -------------------------------
# ROOT macro emission
# -------------------------------

def emit_placement_cpp(parent_def: str, child_def: str, trsf: gp_Trsf, copy_no: int, scale_to_cm: float) -> str:
    parent_cpp = cpp_var_for_def(parent_def)
    child_cpp = cpp_var_for_def(child_def)
    tr_name = f"tr_{sanitize_cpp_name(parent_def)}_{sanitize_cpp_name(child_def)}_{copy_no}"
    return trsf_to_tgeo(trsf, tr_name, scale_to_cm) + f"  {parent_cpp}->AddNode({child_cpp}, {copy_no}, {tr_name});\n"



def _compute_density_g_cm3(
    volume_cm3: float,
    mass_value: float,
    mass_unit: str,
) -> Tuple[Optional[float], str]:
    """
    Computes an effective part density from (mass, CAD volume).

    Returns (rho_g_cm3 or None, comment). If rho is None, caller should fall back
    to the Geant4 NIST density (if resolved) or to a dummy density.
    """
    if not volume_cm3 or volume_cm3 <= 0:
        return None, "no CAD volume available for density"

    if (mass_value is None) or (isinstance(mass_value, float) and math.isnan(mass_value)):
        return None, "no BOM mass available for density"

    mass_g = float(mass_value)
    mu = (mass_unit or "kg").lower()
    if mu == "kg":
        mass_g *= 1000.0
    elif mu == "g":
        pass
    else:
        # unknown unit: assume kg
        mass_g *= 1000.0

    rho = mass_g / float(volume_cm3)
    # Guard against obvious unit/volume issues
    if not (0.01 < rho < 50.0):
        return None, f"computed density {rho:.3g} g/cm3 rejected (unit mismatch?)"

    return rho, "density from BOM mass and CAD volume"


def emit_root_macro(
    step_path: str,
    out_folder: _Path,
    meshparam=None,
    step_unit: str = "auto",
    clip_box: Optional[ClipBox] = None,
    clip_deduplicate: str = "intact",
    name_filter: Optional[NameFilter] = None,
    materials_csv: Optional[str] = None,
    bom_mass_unit: str = "kg",
    g4_nist_json: Optional[str] = None,
    mat_cfg: Optional[MatMatchConfig] = None,
    surface_report: Optional[str] = None,
    surface_files: Optional[Dict[str, str]] = None,
    exact_surfaces: str = "off",
    recognize_surfaces: str = "exact",
    dump_brep: bool = False,
    csg: str = "off",
    csg_report: Optional[str] = None,
):
    # surface_files: def_lid -> absolute path of an exact-surface sidecar (surfaces_*.bin).
    # Volumes listed here are emitted as O2BVHSurfaceSolid via emit_surface_solid_cpp;
    # all others keep the tessellated fallback. This map is normally populated internally
    # from exact_surfaces ("off"|"auto"|"required"); an explicit surface_files argument (used
    # by tests) is honoured as-is and only when exact_surfaces == "off".
    #
    # exact_surfaces mode:
    #   off      : tessellated output only (default; leaves generated output unchanged).
    #   auto     : emit O2BVHSurfaceSolid for every leaf solid whose faces all extract
    #              exactly, tessellated fallback otherwise.
    #   required : like auto, but abort if any leaf solid cannot be represented exactly.
    #
    # dump_brep: with exact_surfaces auto|required, also write the OCCT BREP of every leaf solid
    # that extracts exactly, as brep_<VOLNAME>_<LID>.brep next to its surfaces_*.bin sidecar
    # (same name suffix, so the two pair up). The shape is scaled to cm first, so the BREP is in
    # the same units as the sidecar and the facet mesh. See write_brep_cm().
    if (step_unit or "auto").lower() == "auto":
        detected = detect_step_length_unit(step_path)
        scale_to_cm = step_unit_scale_to_cm(detected)
        print(f"Detected STEP length unit: {detected} (scale to cm = {scale_to_cm})")
    else:
        scale_to_cm = step_unit_scale_to_cm(step_unit)
        print(f"Using overridden STEP length unit: {step_unit} (scale to cm = {scale_to_cm})")

    if clip_box is not None:
        print(f"Clipping CAD geometry to STEP-coordinate bounding box: {clip_box.as_tuple()}")
        print(f"Clip deduplication mode: {clip_deduplicate}")

    if name_filter is not None and name_filter.active:
        print(f"CAD name filters: {len(name_filter.include)} include regex(es), {len(name_filter.exclude)} exclude regex(es)")

    extract_graph(
        step_path,
        meshparam=meshparam,
        scale_to_cm=scale_to_cm,
        clip_box=clip_box,
        clip_deduplicate=clip_deduplicate,
        name_filter=name_filter,
    )

    out_folder = out_folder.expanduser().resolve()
    out_folder.mkdir(parents=True, exist_ok=True)

    recognize_mode = (recognize_surfaces or "exact").lower()
    if recognize_mode not in ("exact", "off"):
        raise ValueError(f"recognize_surfaces must be exact|off, got {recognize_surfaces!r}")
    recognize_flag = recognize_mode == "exact"

    # --- optional exact-surface eligibility report (does not modify the emitted geometry) ---
    if surface_report:
        report = build_surface_report(step_path, scale_to_cm, recognize_surfaces=recognize_flag)
        report_path = _Path(surface_report).expanduser().resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=1))
        summ = report["summary"]
        print(f"Surface report: {summ['n_eligible']}/{summ['n_volumes']} logical volumes eligible "
              f"for exact O2BVHSurfaceSolid conversion")
        print(f"  face types: {summ['face_type_counts']}")
        if summ["recognized_surface_counts"]:
            print(f"  recognized (stored type is not the geometry): {summ['recognized_surface_counts']}"
                  f" recovered from stored {summ['recognized_stored_type_counts']}")
        if summ["fallback_reasons"]:
            top = sorted(summ["fallback_reasons"].items(), key=lambda kv: -kv[1])[:5]
            for reason, count in top:
                print(f"  fallback ({count}x): {reason}")
        print(f"Wrote surface report: {report_path}")

    # --- exact-surface extraction (auto/required modes) ---
    exact_mode = (exact_surfaces or "off").lower()
    if exact_mode not in ("off", "auto", "required"):
        raise ValueError(f"exact_surfaces must be off|auto|required, got {exact_surfaces!r}")
    if exact_mode == "off":
        surface_files = surface_files or {}
    else:
        if surface_files:
            raise ValueError("surface_files cannot be combined with --exact-surfaces auto|required")
        surface_files = {}
        brep_files: Dict[str, str] = {}  # def_lid -> absolute path of brep_*.brep (--dump-brep)
        failures: Dict[str, List[str]] = {}  # def_lid -> unsupported-face reasons
        for lid, shape in def_shapes.items():
            surfaces, reasons, n_model_edges = extract_surfaces_for_shape(
                shape, scale_to_cm, recognize_surfaces=recognize_flag)
            if surfaces is None:
                failures[lid] = reasons
                continue
            disp = def_names.get(lid, "")
            volname = sanitize_filename(disp) if disp else "vol"
            name_suffix = f"{volname}_{sanitize_filename(lid)}"
            fpath = (out_folder / f"surfaces_{name_suffix}.bin").resolve()
            write_surfaces_bin(fpath, surfaces, shape_model_tolerance(shape, scale_to_cm), n_model_edges)
            surface_files[lid] = str(fpath)
            if dump_brep:
                bpath = (out_folder / f"brep_{name_suffix}.brep").resolve()
                write_brep_cm(bpath, shape, scale_to_cm)
                brep_files[lid] = str(bpath)
        if dump_brep:
            print(f"Wrote {len(brep_files)} reference BREP file(s) (brep_*.brep, scaled to cm)")
        n_leaf = len(def_shapes)
        print(f"Exact-surface extraction ({exact_mode}): {len(surface_files)}/{n_leaf} leaf solids "
              f"represented exactly, {len(failures)} fall back to tessellation")
        if failures:
            # Aggregate reasons for a compact, useful report.
            reason_counts: Dict[str, int] = {}
            for reasons in failures.values():
                for r in reasons:
                    reason_counts[r] = reason_counts.get(r, 0) + 1
            for reason, count in sorted(reason_counts.items(), key=lambda kv: -kv[1]):
                print(f"  fallback ({count} face(s)): {reason}")
            if exact_mode == "required":
                lines = [f"--exact-surfaces required: {len(failures)}/{n_leaf} leaf solid(s) cannot be "
                         f"represented exactly:"]
                for lid in sorted(failures):
                    name = def_names.get(lid, "") or lid
                    uniq = sorted(set(failures[lid]))
                    lines.append(f"  {name} [{lid}]: {'; '.join(uniq)}")
                raise ValueError("\n".join(lines))

    # --- CSG recognition (--csg auto|required) -- the one CSG hook ------------------------
    #
    # The cascade the project asked for is CSG -> exact surfaces -> tessellated, and this is
    # where its first tier is decided. Recognition proposes a placed primitive (or a union of
    # two) from the leaf solid's carrier structure and OCCT's symmetric-difference volume then
    # accepts or rejects it exactly; only an accepted part reaches `csg_files` and only a part
    # in `csg_files` is emitted as a native ROOT shape below. The other representations are
    # still written -- the gate scores all three side by side and this must not remove parts
    # from that comparison. See scripts/geometry/csg/hook.py and Stream_H_CSGEmitter.md.
    csg_mode = (csg or "off").lower()
    if csg_mode not in ("off", "auto", "required"):
        raise ValueError(f"csg must be off|auto|required, got {csg!r}")
    csg_files: Dict[str, str] = {}
    if csg_mode != "off":
        hook = import_csg_hook()
        csg_files, csg_records = hook.recognise_and_emit(
            def_shapes, def_names, scale_to_cm, out_folder, sanitize_filename, mode=csg_mode)
        report_path = _Path(csg_report) if csg_report else (out_folder / "csg_report.json")
        report = hook.write_report(csg_records, report_path, set(surface_files or {}),
                                   set(logical_volumes))
        hook.print_tier_table(report)
        print(f"Wrote CSG report: {report_path}")

    # --- Geant4 NIST material DB (optional but recommended) ---
    g4db: Optional[Dict[str, dict]] = None
    if g4_nist_json:
        g4db = load_g4_nist_db(g4_nist_json)
        print(f"Loaded Geant4 NIST DB with {len(g4db)} materials from: {g4_nist_json}")
    else:
        print("No --g4-nist-json provided: unresolved materials will fall back to dummy ROOT materials.")
    mat_cfg = mat_cfg or MatMatchConfig()


    # --- BOM: map volumes to materials (heuristic) ---
    lid_to_bom: Dict[str, BomEntry] = {}
    if materials_csv:
        bom_entries = read_bom_csv(materials_csv)
        lid_to_bom = build_volume_to_material_map(bom_entries, def_names)
        print(f"Loaded {len(bom_entries)} BOM entries from: {materials_csv}")
        print(f"Matched {len(lid_to_bom)} CAD logical volumes to BOM entries (by name/part-number heuristics)")
    else:
        print("No --materials-csv provided: emitting Default medium for all logical volumes")

    # --- facet files ---
    facet_files = {}  # def_lid -> absolute path string
    for lid, tris in logical_volumes.items():
        disp = def_names.get(lid, "")
        volname = sanitize_filename(disp) if disp else "vol"
        lidname = sanitize_filename(lid)
        fname = f"facets_{volname}_{lidname}.bin"
        fpath = (out_folder / fname).resolve()
        write_facets_bin(fpath, tris)
        facet_files[lid] = str(fpath).replace("\\", "\\\\")  # C++ string literal safety

    # --- which materials do we need to emit? ---
    
    # --- materials: collect unique BOM material strings actually used by leaf volumes ---
    # We resolve each unique BOM string to a Geant4 NIST material using string + density scoring.
    used_materials: Dict[str, ResolvedMaterial] = {}

    # Precompute one representative part density per BOM material (first good value wins)
    mat_to_rho: Dict[str, Optional[float]] = {}
    mat_to_rho_note: Dict[str, str] = {}

    for lid in logical_volumes.keys():
        if lid not in lid_to_bom:
            continue
        bom = lid_to_bom[lid]
        mat_name = normalize_material_name(bom.material)

        if mat_name not in mat_to_rho:
            rho_part, rho_note = _compute_density_g_cm3(
                def_volumes_cm3.get(lid, 0.0),
                bom.mass_value,
                bom_mass_unit,
            )
            mat_to_rho[mat_name] = rho_part
            mat_to_rho_note[mat_name] = rho_note

    for mat_name in sorted(mat_to_rho.keys(), key=lambda s: s.lower()):
        rho_part = mat_to_rho.get(mat_name)
        rm = resolve_bom_material(mat_name, rho_part, g4db, mat_cfg)

        # Fold density provenance into the note for geom.C comments
        rm.note = f"{rm.note} (density: {mat_to_rho_note.get(mat_name, 'n/a')})"

        if rm.nist_name is None:
            print(f"WARNING: Unresolved/ambiguous material '{mat_name}'. See FIXME in generated geom.C.")

        used_materials[mat_name] = rm

    materials_cpp, medium_var_map = emit_materials_cpp(used_materials)

    # --- emit C++ macro ---
    surface_files = surface_files or {}
    unknown_surface_lids = set(surface_files) - set(logical_volumes)
    if unknown_surface_lids:
        raise ValueError(f"surface_files references unknown logical volumes: {sorted(unknown_surface_lids)}")
    if surface_files:
        print(f"Emitting {len(surface_files)}/{len(logical_volumes)} logical volumes as exact O2BVHSurfaceSolid "
              f"(macro requires the ALICE O2 environment)")

    cpp: List[str] = []
    cpp.append(emit_cpp_prelude(exact_surfaces=bool(surface_files), csg_shapes=bool(csg_files)))

    cpp.append("TGeoVolume* build(bool check=true) {")
    cpp.append('  if (!gGeoManager) { throw std::runtime_error("gGeoManager is null. Call build_and_export() or create a TGeoManager first."); }')
    cpp.append(materials_cpp)

    for lid in logical_volumes.keys():
        ntriangles = len(logical_volumes[lid])

        # choose medium for this volume
        med = "med_Default"
        if lid in lid_to_bom:
            mat_name = normalize_material_name(lid_to_bom[lid].material)
            med = medium_var_map.get(mat_name, "med_Default")

        # The cascade, in one place: CSG, else exact surfaces, else the tessellated fallback.
        if lid in csg_files:
            shape_path = str(_Path(csg_files[lid]).expanduser().resolve()).replace("\\", "\\\\")
            cpp.append(import_csg_hook().emit_csg_shape_cpp(
                lid, def_names.get(lid, ""), shape_path, med, sanitize_cpp_name))
        elif lid in surface_files:
            sidecar = str(_Path(surface_files[lid]).expanduser().resolve()).replace("\\", "\\\\")
            cpp.append(emit_surface_solid_cpp(lid, def_names.get(lid, ""), sidecar, med))
        else:
            cpp.append(emit_tessellated_cpp(lid, def_names.get(lid, ""), facet_files[lid], ntriangles, med))

    for lid in sorted(assemblies):
        cpp.append(emit_assembly_cpp(lid, def_names.get(lid, "")))

    for idx, (parent, child, trsf) in enumerate(placements, start=1):
        cpp.append(emit_placement_cpp(parent, child, trsf, idx, scale_to_cm))

    if len(top_defs) == 1:
        top = next(iter(top_defs))
        cpp.append(f"  return {cpp_var_for_def(top)};")
    else:
        cpp.append('  TGeoVolumeAssembly *asm_WORLD = new TGeoVolumeAssembly("WORLD");')
        for i, node in enumerate(sorted(top_defs), start=1):
            cpp.append(f"  asm_WORLD->AddNode({cpp_var_for_def(node)}, {i});")
        cpp.append("  return asm_WORLD;")

    cpp.append("}")

    # exports a function allowing to export the geometry to TGeo file
    cpp.append('void build_and_export(const char* out_root = "geom.root", bool check=true) {')
    cpp.append('  if (!gGeoManager) { new TGeoManager("geom","geom"); }')
    cpp.append('  TGeoVolume* top = build(check);')
    cpp.append('  gGeoManager->SetTopVolume(top);')
    cpp.append('  gGeoManager->CloseGeometry();')
    cpp.append('  gGeoManager->CheckOverlaps();')
    cpp.append('  gGeoManager->Export(out_root);')
    cpp.append('}')

    # exports a function to get get hold of the builder function in ALICE O2
    cpp.append('std::function<TGeoVolume*()> get_builder_hook_checked() {')
    cpp.append('  return []() { return build(true); };')
    cpp.append('}')
    # exports a function to get get hold of the builder function in ALICE O2
    cpp.append('std::function<TGeoVolume*()> get_builder_hook_unchecked() {')
    cpp.append('  return []() { return build(false); };')
    cpp.append('}')

    return "\n".join(cpp)


# -------------------------------
# Geometry Tree printing (debug)
# -------------------------------

def label_entry(label):
    s = TCollection_AsciiString()
    TDF_Tool.Entry(label, s)
    return s.ToCString()


def traverse_print(label, shape_tool, depth=0):
    indent = "  " * depth
    name = label.GetLabelName()
    entry = label_entry(label)
    print(f"{indent}- {name}  =>[{entry}]") 

    if shape_tool.IsReference(label):
        ref_label = TDF_Label()
        shape_tool.GetReferredShape(label, ref_label)
        traverse_print(ref_label, shape_tool, depth + 1)
        return

    children = TDF_LabelSequence()
    shape_tool.GetComponents(label, children)
    if children.Length() > 0 or shape_tool.IsAssembly(label):
        for i in range(children.Length()):
            traverse_print(children.Value(i + 1), shape_tool, depth + 1)
        return

    if shape_tool.IsSimpleShape(label):
        shape = shape_tool.GetShape(label)
        print(f"{indent}  [LogicalShape id={id(shape)}]")


def print_geom(step_file):
    print(f"Printing GEOM hierarchy for {step_file}")
    doc, shape_tool = load_step_with_xcaf(step_file)
    roots = TDF_LabelSequence()
    shape_tool.GetFreeShapes(roots)
    for i in range(roots.Length()):
        traverse_print(roots.Value(i + 1), shape_tool)


# -------------------------------
# CLI
# -------------------------------

def main():
    ap = argparse.ArgumentParser(description="Convert STEP/XCAF to ROOT TGeo macro, facets in per-volume binary files.")
    ap.add_argument("step", help="Input STEP file")
    ap.add_argument("-o", "--out", default="geom.C", help="Output ROOT macro file name (default: geom.C)")
    ap.add_argument("--output-folder", default="./", help="Output folder for macro + facet files")
    ap.add_argument("--out-path", default=None, help="(deprecated) Alias for --output-folder")
    ap.add_argument("--mesh", action="store_true", help="Use full BRepMesh triangulation instead of bounding boxes")
    ap.add_argument("--print-tree", action="store_true", help="Just prints the geometry tree")
    ap.add_argument("--mesh-prec", default=0.1, help="meshing precision. lower --> slower")
    ap.add_argument("--step-unit", default="auto", choices=["auto", "mm", "cm", "m", "in", "ft"], help="STEP length unit override (default: auto-detect); TGeo expects cm")
    ap.add_argument("--clip-box", nargs=6, type=float, metavar=("XMIN", "YMIN", "ZMIN", "XMAX", "YMAX", "ZMAX"), default=None, help="Clip CAD geometry to this axis-aligned bounding box before meshing (coordinates in STEP file units, before conversion to cm)")
    ap.add_argument("--clip-deduplicate", default="intact", choices=["none", "intact"], help="When clipping, reuse original logical definitions for subtrees fully inside the clip box (default: intact); use 'none' for one volume per surviving occurrence")
    ap.add_argument("--include-name", action="append", default=[], help="Only convert CAD labels whose XCAF name or label entry matches this regex; may be repeated. Matching an assembly includes its subtree.")
    ap.add_argument("--exclude-name", action="append", default=[], help="Skip CAD labels/subtrees whose XCAF name or label entry matches this regex; may be repeated.")
    ap.add_argument("--name-filter-case-sensitive", action="store_true", help="Make --include-name/--exclude-name matching case-sensitive (default: case-insensitive)")
    ap.add_argument("--surface-report", default=None, metavar="PATH", help="Write a JSON report classifying each face by analytic surface type and each logical volume by exact O2BVHSurfaceSolid conversion eligibility. Does not change the generated geometry output.")
    ap.add_argument("--exact-surfaces", default="off", choices=["off", "auto", "required"], help="Emit exact O2BVHSurfaceSolid shapes instead of TGeoTessellated where possible. 'off' (default): tessellated only. 'auto': exact for each leaf solid whose faces all extract exactly, tessellated fallback otherwise. 'required': fail with a report if any leaf solid cannot be represented exactly. Writes a surfaces_*.bin sidecar per exact volume.")
    ap.add_argument("--dump-brep", action="store_true", help="With --exact-surfaces auto|required, also write brep_<VOLNAME>_<LID>.brep (OCCT BREP of the leaf solid, scaled to cm like the sidecar and the mesh) next to each surfaces_*.bin. Input for the OCCT reference oracle; changes nothing else in the output.")
    ap.add_argument("--csg", default="off", choices=["off", "auto", "required"], help="Recognise leaf solids as native ROOT CSG shapes (TGeoBBox/TGeoTube/TGeoTubeSeg/TGeoCone/TGeoSphere, and the two-cluster TGeoTube-union of a barrel and a lug) and emit each accepted one as shape_<VOLNAME>_<LID>.root. 'off' (default): unchanged behaviour. 'auto': makes the per-part cascade CSG -> exact surfaces -> tessellated. 'required': fail with a report if any leaf solid is not CSG. A part is only converted this way when OCCT's symmetric-difference volume against the CAD solid is inside the model tolerance; the description and that evidence are written to csg_<VOLNAME>_<LID>.json and csg_report.json either way.")
    ap.add_argument("--csg-report", default=None, metavar="PATH", help="Where to write the per-part CSG cascade report (default: csg_report.json in the output folder).")
    ap.add_argument("--recognize-surfaces", default="exact", choices=["exact", "off"], help="Canonical-form recognition pre-pass: recover an exact plane/sphere/cylinder/cone hiding behind a stored bspline/bezier/revolution/extrusion face (the stored STEP surface type describes the exporter, not the geometry). 'exact' (default): only accept a fit at machine precision. 'off': disable, keeping such faces on the tessellated fallback. Applies to both --surface-report and --exact-surfaces auto|required.")

    # NEW: BOM / material support
    ap.add_argument("--materials-csv", default=None, help="BOM CSV file providing material + mass per part (optional)")
    ap.add_argument("--bom-mass-unit", default="kg", choices=["kg", "g"], help="Unit of the BOM mass column (default: kg)")
    ap.add_argument("--g4-nist-json", default=None, help="Path to Geant4 NIST DB JSON dump (from nist_export_all). Enables TGeoMixture emission + RadLen/IntLen.")


    # Material matching scoring knobs (only used if --g4-nist-json is provided)
    ap.add_argument("--mat-min-score", type=float, default=0.35, help="Minimum combined score to accept a G4 NIST material match (default: 0.35)")
    ap.add_argument("--mat-ambiguity-delta", type=float, default=0.05, help="If best-second < delta, treat match as ambiguous/unresolved (default: 0.05)")
    ap.add_argument("--mat-w-token", type=float, default=0.75, help="Weight for token/name similarity score (default: 0.75)")
    ap.add_argument("--mat-w-density", type=float, default=0.25, help="Weight for density proximity score (default: 0.25)")
    ap.add_argument("--mat-max-log-density-diff", type=float, default=0.0, help="Optional hard density filter in log-space (0 disables). Example 0.8 ~ within 2.2x (default: 0.0)")
    ap.add_argument("--mat-compound-penalty", type=float, default=0.25, help="Penalty for matching to oxides/carbides/etc. when BOM doesn't mention them (default: 0.25)")

    args = ap.parse_args()

    step_path = str(_Path(args.step).expanduser().resolve())
    if args.print_tree:
        print_geom(step_path)
        return

    out_folder = _Path(args.output_folder)
    if args.out_path is not None:
        out_folder = _Path(args.out_path)

    clip_box = None
    if args.clip_box is not None:
        try:
            clip_box = ClipBox.from_values(args.clip_box)
        except ValueError as exc:
            ap.error(str(exc))

    name_filter = None
    if args.include_name or args.exclude_name:
        try:
            name_filter = NameFilter.from_patterns(
                args.include_name,
                args.exclude_name,
                case_sensitive=args.name_filter_case_sensitive,
            )
        except re.error as exc:
            ap.error(f"Invalid CAD name filter regex: {exc}")

    meshparam = {"do_meshing": args.mesh, "lin_defl": args.mesh_prec, "ang_defl": args.mesh_prec}


    mat_cfg = MatMatchConfig(
    min_score=args.mat_min_score,
    ambiguity_delta=args.mat_ambiguity_delta,
    w_token=args.mat_w_token,
    w_density=args.mat_w_density,
    max_log_density_diff=args.mat_max_log_density_diff,
    compound_penalty=args.mat_compound_penalty,
    )

    out_folder = out_folder.expanduser().resolve()
    out_folder.mkdir(parents=True, exist_ok=True)

    out_macro = (out_folder / _Path(args.out).name).resolve()
    code = emit_root_macro(
        step_path,
        out_folder,
        meshparam=meshparam,
        step_unit=args.step_unit,
        clip_box=clip_box,
        clip_deduplicate=args.clip_deduplicate,
        name_filter=name_filter,
        materials_csv=args.materials_csv,
        bom_mass_unit=args.bom_mass_unit,
        g4_nist_json=args.g4_nist_json,
        mat_cfg=mat_cfg,
        surface_report=args.surface_report,
        exact_surfaces=args.exact_surfaces,
        recognize_surfaces=args.recognize_surfaces,
        dump_brep=args.dump_brep,
        csg=args.csg,
        csg_report=args.csg_report,
    )
    out_macro.write_text(code)

    print(f"Wrote ROOT macro: {out_macro}")
    print(f"Wrote facet files into: {out_folder}")
    print("In ROOT you can do:")
    print(f"  root -l {out_macro}")
    print('  build_and_export("geom.root");')


if __name__ == "__main__":
    main()
