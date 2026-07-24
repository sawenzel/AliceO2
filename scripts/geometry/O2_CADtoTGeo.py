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
from dataclasses import dataclass
from pathlib import Path as _Path
from typing import Dict, List, Optional, Pattern, Tuple

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
from OCC.Core.GeomAbs import (
    GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus,
    GeomAbs_BezierSurface, GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
    GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface,
    GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse, GeomAbs_Hyperbola, GeomAbs_Parabola,
    GeomAbs_BezierCurve, GeomAbs_BSplineCurve, GeomAbs_OffsetCurve, GeomAbs_OtherCurve,
)
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopAbs import TopAbs_REVERSED, TopAbs_WIRE, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Extend.TopologyUtils import TopologyExplorer

from OCC.Core.STEPCAFControl import STEPCAFControl_Reader
from OCC.Core.TDocStd import TDocStd_Document
from OCC.Core.XCAFDoc import XCAFDoc_DocumentTool
from OCC.Core.IFSelect import IFSelect_RetDone

from OCC.Core.TDF import TDF_Label, TDF_LabelSequence, TDF_Tool
from OCC.Core.TCollection import TCollection_AsciiString
from OCC.Core.gp import gp_Pnt, gp_Trsf

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
#  - planar faces: general polygon/arc wires -> boundary curves must be lines or circles
#  - cylinder/cone/sphere: only parametric-rectangle trims (Add*Surface API)
_SUPPORTED_SURFACE_TYPES = {"plane", "cylinder", "cone", "sphere"}
_SUPPORTED_PLANAR_CURVES = {"line", "circle"}


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


def classify_face(face, scale_to_cm: float) -> dict:
    """Classifies a single TopoDS face: analytic type, parameters, wires and edges."""
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

    # quadrics: only a single parametric-rectangle trim is representable today
    if not all(w["all_pcurves_iso"] for w in record["wires"]):
        record["trim_kind"] = "general"
        return False, f"{surf_type} with general (non-parametric-rectangle) trim"
    if len(record["wires"]) > 1:
        record["trim_kind"] = "general"
        return False, f"{surf_type} with multiple trim wires (holes)"
    record["trim_kind"] = "parametric-rectangle"
    return True, None


def build_surface_report(step_path: str, scale_to_cm: float) -> dict:
    """Builds the JSON-serializable exact-conversion eligibility report over def_shapes."""
    volumes = {}
    n_eligible = 0
    face_type_counts: Dict[str, int] = {}
    curve_type_counts: Dict[str, int] = {}
    fallback_reasons: Dict[str, int] = {}

    for lid, shape in def_shapes.items():
        faces = []
        for face in TopologyExplorer(shape).faces():
            rec = classify_face(face, scale_to_cm)
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
SURFACE_SIDECAR_VERSION = 1
SURFACE_TYPE_ENUM = {"plane": 1, "cylinder": 2, "cone": 3, "sphere": 4}
CURVE_TYPE_ENUM = {"line": 0, "arc": 1}
SURFACE_FLAG_INNER_WALL = 1 << 0


def write_surfaces_bin(path: _Path, surfaces: List[dict]):
    """
    Writes the surface sidecar. `surfaces` is a list of records:
      {"type": "plane"|"cylinder"|"cone"|"sphere",
       "inner_wall": bool,              # quadrics only, default False
       "params": [float, ...],          # fixed per-type layout, see BVHSurfaceSolid.md
       "wires": [                       # planes only; quadrics carry their trim in params
          {"role": "outer"|"inner",     # general line/arc loops (polygon, disk, rounded rect)
           "edges": [{"curve": "line"|"arc", "params": [float, ...]}, ...]},
       ]}
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(SURFACE_SIDECAR_MAGIC)
        f.write(struct.pack("<III", SURFACE_SIDECAR_VERSION, len(surfaces), 0))
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


# -------------------------------
# Exact-surface extraction (--exact-surfaces)
# -------------------------------
# Turn OpenCascade faces into the sidecar surface records written by write_surfaces_bin.
# This is the exact counterpart of the tessellation path: a leaf solid is emitted as an
# O2BVHSurfaceSolid only when *every* one of its faces can be extracted exactly; otherwise
# the caller keeps the tessellated fallback (auto mode) or fails (required mode).
#
# Milestone status (scripts/geometry/BVHSurfaceSolid.md): planar faces carry general line/arc
# boundary wires (polygon, disk/annulus, rounded rectangle, slot, ...) and the three quadrics
# (cylinder, cone, sphere) are extracted here. A quadric face with a single iso-parametric wire
# stays on the scalar parametric-rectangle path (byte-identical); any other supported trim (a
# hole/window, or a non-iso straight-line boundary) is emitted as a general line trim wire in
# the surface's (phi, height/theta) domain. Quadric arc/spline pcurve boundaries still force a
# fallback (a circular pcurve does not survive the non-uniform u->phi, v->length/latitude remap
# as a circle), as do torus/free-form surfaces and non-line/arc planar boundary curves.

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


def extract_planar_face(face, scale_to_cm: float) -> Tuple[Optional[dict], Optional[str]]:
    """Convert a planar TopoDS face into a sidecar 'plane' surface record with general
    line/arc boundary wires (polygon, disk/annulus, rounded rectangle, slot, ...).

    Each wire is walked in connected order; straight edges become 'line' segments and circular
    edges become 'arc' segments in the plane's local (u, v) frame. Boundary edges that are
    neither straight lines nor circular arcs (ellipses, splines, ...) force a fallback.
    """
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
            if gt not in (GeomAbs_Line, GeomAbs_Circle):
                name = _CURVE_TYPE_NAMES.get(gt, "unknown")
                return None, f"planar boundary edge is a {name} curve (only line/circle supported)"
            classified.append((edge, curve, gt, project(BRep_Tool.Pnt(start_vertex))))

        n = len(classified)
        n_arcs = sum(1 for _, _, gt, _ in classified if gt == GeomAbs_Circle)
        if n_arcs == 0 and n < 3:
            return None, "planar polygon wire has fewer than 3 edges"
        if n == 0:
            return None, "planar face has an empty wire"

        seg_edges = []
        for i, (edge, curve, gt, start_uv) in enumerate(classified):
            if gt == GeomAbs_Line:
                u0, v0 = start_uv
                u1, v1 = classified[(i + 1) % n][3]
                seg_edges.append({"curve": "line", "params": [u0, v0, u1, v1]})
            else:
                params, reason = _arc_edge_params(edge, curve, project, s)
                if params is None:
                    return None, reason
                seg_edges.append({"curve": "arc", "params": params})
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


def _quadric_line_wire(face, phi_of_u, v_of_v) -> Tuple[Optional[List[dict]], Optional[str]]:
    """Build general line-only trim wires in a quadric face's parametric (phi, v) domain.

    Every boundary edge's 2D pcurve must be a straight line: a circular pcurve arc does not
    survive the non-uniform (u -> phi, v -> height/theta) remap as a circle (u is in radians and
    v carries a length/latitude scale, so the image is an ellipse), so arc/spline pcurves force a
    fallback. \a phi_of_u maps the OCC angular u to the C++ phi (identity or negation for a
    left-handed frame); \a v_of_v maps the OCC v to the C++ height (cm) or polar angle (rad).
    Returns (wires, None) with exactly one outer wire plus inner holes, or (None, reason).
    """
    wires_out: List[dict] = []
    for _wire, is_outer, edges in _face_wire_edges(face):
        points = []  # ordered (phi, v) start points in wire-traversal order
        for edge, _start_vertex in edges:
            curve2d, first, last = BRep_Tool.CurveOnSurface(edge, face)
            if curve2d is None:
                return None, "quadric boundary edge has no 2D pcurve"
            if Geom2dAdaptor_Curve(curve2d).GetType() != GeomAbs_Line:
                return None, "quadric boundary edge pcurve is not a straight line (only line trims supported)"
            param = last if edge.Orientation() == TopAbs_REVERSED else first
            p = curve2d.Value(param)
            points.append((phi_of_u(p.X()), v_of_v(p.Y())))
        if len(points) < 3:
            return None, "quadric trim wire has fewer than 3 line edges"
        seg_edges = []
        n = len(points)
        for i in range(n):
            u0, v0 = points[i]
            u1, v1 = points[(i + 1) % n]
            seg_edges.append({"curve": "line", "params": [u0, v0, u1, v1]})
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
    wires, reason = _quadric_line_wire(face, phi_of_u, lambda v: v * s)
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
    wires, reason = _quadric_line_wire(face, phi_of_u, lambda v: v * cos_a * s)
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
    wires, reason = _quadric_line_wire(face, phi_of_u, lambda v: 0.5 * math.pi - v)
    if wires is None:
        return None, reason
    record["wires"] = wires
    return record, None


# Face extractors dispatched by analytic surface type. Planar faces cover both straight-edged
# polygons and circular disks/annuli; the quadrics carry their trim in the surface parameters.
_FACE_EXTRACTORS = {
    "plane": extract_planar_face,
    "cylinder": extract_cylindrical_face,
    "cone": extract_conical_face,
    "sphere": extract_spherical_face,
}


def extract_surfaces_for_shape(shape, scale_to_cm: float) -> Tuple[Optional[List[dict]], List[str]]:
    """Attempt to extract every face of a leaf solid into exact sidecar surface records.

    Returns (surfaces, []) when all faces are supported, or (None, reasons) listing why the
    solid cannot be represented exactly (one reason per unsupported face). A shape with no
    faces is treated as unsupported.
    """
    surfaces: List[dict] = []
    reasons: List[str] = []
    n_faces = 0
    for face in TopologyExplorer(shape).faces():
        n_faces += 1
        adaptor = BRepAdaptor_Surface(face)
        surf_type = _SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")
        extractor = _FACE_EXTRACTORS.get(surf_type)
        if extractor is None:
            reasons.append(f"{surf_type} face extraction not implemented yet")
            continue
        record, reason = extractor(face, scale_to_cm)
        if record is None:
            reasons.append(reason or f"{surf_type} face not supported")
        else:
            surfaces.append(record)
    if n_faces == 0:
        return None, ["shape has no faces"]
    if reasons:
        return None, reasons
    return surfaces, []


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


def emit_cpp_prelude(exact_surfaces: bool = False) -> str:
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

    # --- optional exact-surface eligibility report (does not modify the emitted geometry) ---
    if surface_report:
        report = build_surface_report(step_path, scale_to_cm)
        report_path = _Path(surface_report).expanduser().resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=1))
        summ = report["summary"]
        print(f"Surface report: {summ['n_eligible']}/{summ['n_volumes']} logical volumes eligible "
              f"for exact O2BVHSurfaceSolid conversion")
        print(f"  face types: {summ['face_type_counts']}")
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
        failures: Dict[str, List[str]] = {}  # def_lid -> unsupported-face reasons
        for lid, shape in def_shapes.items():
            surfaces, reasons = extract_surfaces_for_shape(shape, scale_to_cm)
            if surfaces is None:
                failures[lid] = reasons
                continue
            disp = def_names.get(lid, "")
            volname = sanitize_filename(disp) if disp else "vol"
            fpath = (out_folder / f"surfaces_{volname}_{sanitize_filename(lid)}.bin").resolve()
            write_surfaces_bin(fpath, surfaces)
            surface_files[lid] = str(fpath)
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
    cpp.append(emit_cpp_prelude(exact_surfaces=bool(surface_files)))

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

        if lid in surface_files:
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
    )
    out_macro.write_text(code)

    print(f"Wrote ROOT macro: {out_macro}")
    print(f"Wrote facet files into: {out_folder}")
    print("In ROOT you can do:")
    print(f"  root -l {out_macro}")
    print('  build_and_export("geom.root");')


if __name__ == "__main__":
    main()
