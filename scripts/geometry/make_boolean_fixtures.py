#!/usr/bin/env python3
"""
Generate the synthetic Boolean fixture ladder used as the test-driven backbone for the exact
BREP solid (O2BVHSurfaceSolid) and its CAD converter.  See CodeReview_Fable.md, Section 9,
gate "G1 — the synthetic Boolean ladder".

Why this exists
---------------
Real CAD models (Bagger, as1-oc, oTOF) are far too large to debug an exact-surface solid
against: when a part fails, there is no way to tell *which* geometric feature broke it.  This
script builds a ladder of tiny, fully understood solids with pythonOCC, in increasing order of
difficulty, so that a failure names its own cause.  Each fixture is exported as a STEP file plus
a `fixtures.json` manifest carrying, per fixture, the geometric feature it exercises, the
closed-form volume where one exists, and the volume/validity OCCT itself computes.

Units
-----
Everything is modelled and written in MILLIMETRES, which is what real CAD exports use.  This is
deliberate: it exercises the converter's mm -> cm scaling path (`scale_to_cm`), which is a
frequent source of silent factor-10 bugs.  Volumes in the manifest are always reported in cm^3
(cm^3 = mm^3 * 1e-3), i.e. in the units the converted TGeo/O2BVHSurfaceSolid geometry uses.

The ladder
----------
Fixtures 1-3 are sanity/regression cases: every trim curve on every face is a line or a circle,
so they are representable by per-face 2D trim curves in closed form.

Fixtures 4, 5 and 6 (`cyl_cross_cyl`, `cyl_inter_cyl`, `tube_window`) are the mathematically
hard ones and the whole reason the ladder exists.  Their boundaries contain the intersection
curve of two orthogonal cylinders.  That curve is a *transcendental* curve in each cylinder's
own (phi, h) parametric chart -- on the cylinder of radius R with axis z, cut by a cylinder of
radius r with axis x, the curve is h(phi) = +- sqrt(r^2 - R^2 sin^2 phi), which is not a line,
circle, conic, or polynomial in (phi, h).  Consequently *no* exact per-face 2D trim curve of the
kinds a converter can emit (segment / arc / conic / B-spline of finite degree) reproduces it
exactly.  Any representation that trims each face independently in its own chart can only
approximate them.  The way out (CodeReview_Fable.md Section 4.2) is to trim against the
*neighbouring surface* rather than against a 2D curve; these three fixtures are the minimal
cases that force that design.  Fixture 6 in particular is the minimized reproducer of the
Bagger `BoomCylinderOuter` failure, and fixture 7 (`oblique_cut_cyl`, whose cut face is an exact
ellipse) is the minimized reproducer of the Bagger `Bucket` failure.

Usage
-----
  python3 make_boolean_fixtures.py                       # generate the full ladder
  python3 make_boolean_fixtures.py --list                # print the ladder, generate nothing
  python3 make_boolean_fixtures.py --only cyl_cross_cyl,tube_window
  python3 make_boolean_fixtures.py --outdir /tmp/fixtures

Requires the pythonOCC environment (same as O2_CADtoTGeo.py); see CodeReview_Fable.md,
Section 11.  The generated .step / fixtures.json are build artifacts and are not meant to be
committed.

Author: ALICE O2 simulation geometry tooling (2026)
"""

import argparse
import json
import math
import re
from pathlib import Path as _Path

from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCone,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeSphere,
    BRepPrimAPI_MakeTorus,
)
from OCC.Core.GProp import GProp_GProps
from OCC.Core.Interface import Interface_Static
from OCC.Core.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf
from OCC.Extend.TopologyUtils import TopologyExplorer

_SCRIPT_DIR = _Path(__file__).resolve().parent
_DEFAULT_OUTDIR = _SCRIPT_DIR / "STEP_examples" / "boolean_fixtures"

MM3_TO_CM3 = 1.0e-3


# -------------------------------
# small shape helpers (all lengths in mm)
# -------------------------------

def _box(dx, dy, dz, corner=(0.0, 0.0, 0.0)):
    return BRepPrimAPI_MakeBox(gp_Pnt(*corner), dx, dy, dz).Shape()


def _cylinder(radius, height, origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0)):
    ax = gp_Ax2(gp_Pnt(*origin), gp_Dir(*direction))
    return BRepPrimAPI_MakeCylinder(ax, radius, height).Shape()


def _cone(r1, r2, height, origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0)):
    ax = gp_Ax2(gp_Pnt(*origin), gp_Dir(*direction))
    return BRepPrimAPI_MakeCone(ax, r1, r2, height).Shape()


def _sphere(radius, centre=(0.0, 0.0, 0.0)):
    return BRepPrimAPI_MakeSphere(gp_Pnt(*centre), radius).Shape()


def _torus(major_r, minor_r, origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0)):
    ax = gp_Ax2(gp_Pnt(*origin), gp_Dir(*direction))
    return BRepPrimAPI_MakeTorus(ax, major_r, minor_r).Shape()


def _rotated_translated(shape, axis_dir, angle_deg, translation):
    """Rotate `shape` about the axis through the origin along `axis_dir`, then translate."""
    rot = gp_Trsf()
    rot.SetRotation(gp_Ax1(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(*axis_dir)), math.radians(angle_deg))
    tra = gp_Trsf()
    tra.SetTranslation(gp_Pnt(0.0, 0.0, 0.0), gp_Pnt(*translation))
    return BRepBuilderAPI_Transform(shape, tra.Multiplied(rot), True).Shape()


def _boolean(op_class, a, b, what):
    op = op_class(a, b)
    op.Build()
    if not op.IsDone():
        raise RuntimeError(f"{what}: boolean operation failed")
    return op.Shape()


def _fuse(a, b):
    return _boolean(BRepAlgoAPI_Fuse, a, b, "fuse")


def _common(a, b):
    return _boolean(BRepAlgoAPI_Common, a, b, "common")


def _cut(a, b):
    return _boolean(BRepAlgoAPI_Cut, a, b, "cut")


# -------------------------------
# the fixture ladder
# -------------------------------

def build_box():
    return _box(20.0, 30.0, 40.0)


def build_box_union_box():
    # Two identical boxes side by side; they share the whole x = 20 face.
    a = _box(20.0, 30.0, 40.0)
    b = _box(20.0, 30.0, 40.0, corner=(20.0, 0.0, 0.0))
    return _fuse(a, b)


def build_box_minus_cyl():
    # 40 mm cube, axial through-hole of radius 8 mm along z.
    cube = _box(40.0, 40.0, 40.0, corner=(-20.0, -20.0, -20.0))
    drill = _cylinder(8.0, 60.0, origin=(0.0, 0.0, -30.0))
    return _cut(cube, drill)


def build_cyl_cross_cyl():
    # Two r = 10 mm, L = 60 mm cylinders, axes z and x, both centred on the origin, FUSED.
    cz = _cylinder(10.0, 60.0, origin=(0.0, 0.0, -30.0), direction=(0.0, 0.0, 1.0))
    cx = _cylinder(10.0, 60.0, origin=(-30.0, 0.0, 0.0), direction=(1.0, 0.0, 0.0))
    return _fuse(cz, cx)


def build_cyl_inter_cyl():
    # The same two cylinders, intersected: the Steinmetz solid.
    cz = _cylinder(10.0, 60.0, origin=(0.0, 0.0, -30.0), direction=(0.0, 0.0, 1.0))
    cx = _cylinder(10.0, 60.0, origin=(-30.0, 0.0, 0.0), direction=(1.0, 0.0, 0.0))
    return _common(cz, cx)


def build_tube_window():
    # Tube r = 15 mm, h = 60 mm (axis z) with a transverse r = 8 mm hole drilled along x.
    tube = _cylinder(15.0, 60.0, origin=(0.0, 0.0, -30.0), direction=(0.0, 0.0, 1.0))
    drill = _cylinder(8.0, 60.0, origin=(-30.0, 0.0, 0.0), direction=(1.0, 0.0, 0.0))
    return _cut(tube, drill)


def build_oblique_cut_cyl():
    # Cylinder r = 12 mm, h = 50 mm (z in [0, 50]) cut by a plane inclined 30 deg to its axis.
    # The cutting half-space is a large box whose bottom face (z = 0) is rotated by 60 deg about
    # x -- a plane at 60 deg to the horizontal is at 30 deg to the z axis -- and then lifted to
    # z = 25 so that the elliptical section stays clear of both caps.
    cyl = _cylinder(12.0, 50.0)
    half_space = _box(400.0, 400.0, 400.0, corner=(-200.0, -200.0, 0.0))
    knife = _rotated_translated(half_space, (1.0, 0.0, 0.0), 60.0, (0.0, 0.0, 25.0))
    return _cut(cyl, knife)


def build_cyl_plus_cone():
    # Cylinder r = 10 mm, h = 30 mm with a coaxial truncated cone (10 -> 5, h = 20 mm) on top.
    cyl = _cylinder(10.0, 30.0)
    cone = _cone(10.0, 5.0, 20.0, origin=(0.0, 0.0, 30.0))
    return _fuse(cyl, cone)


def build_sphere_minus_cyl():
    # Sphere r = 20 mm with an axial r = 6 mm hole drilled through it (a "napkin ring").
    sph = _sphere(20.0)
    drill = _cylinder(6.0, 60.0, origin=(0.0, 0.0, -30.0))
    return _cut(sph, drill)


def build_torus_union_cyl():
    # Torus R = 25 mm, r = 8 mm fused with a coaxial cylinder r = 20 mm, h = 40 mm through it.
    #
    # NOTE (deviation from the original spec, which asked for a r = 10 mm cylinder): the torus
    # tube spans the radial band [R - r, R + r] = [17, 33] mm, so a coaxial cylinder of radius
    # 10 mm never touches it -- the "fuse" would return a compound of two disjoint solids with
    # no torus/cylinder intersection curve at all, defeating the purpose of the fixture.  The
    # cylinder radius is therefore 20 mm, which lies inside the tube band, giving a genuine
    # torus-cylinder junction.  Being coaxial, the junction curves are exact circles (at
    # z = +- sqrt(r^2 - (Rcyl - R)^2) = +- sqrt(39) mm), so what this fixture stresses is the
    # toroidal (quartic) *surface* and its trimming, not a transcendental trim curve.
    tor = _torus(25.0, 8.0)
    cyl = _cylinder(20.0, 40.0, origin=(0.0, 0.0, -20.0))
    return _fuse(tor, cyl)


# Closed-form volumes, in mm^3, where one exists.
_V_BOX = 20.0 * 30.0 * 40.0
_V_STEINMETZ = 16.0 * 10.0 ** 3 / 3.0                      # two equal orthogonal cylinders, r = 10

FIXTURES = [
    {
        "name": "box",
        "build": build_box,
        "description": "20 x 30 x 40 mm box.",
        "feature": "trivial sanity case: 6 planar faces, all trim curves are straight segments",
        "volume_mm3": _V_BOX,
    },
    {
        "name": "box_union_box",
        "build": build_box_union_box,
        "description": "Two 20 x 30 x 40 mm boxes fused along a shared full face at x = 20 mm.",
        "feature": "coplanar/shared-face topology: the fuse must remove the internal face and "
                   "merge the two coplanar face pairs without leaving a seam",
        "volume_mm3": 2.0 * _V_BOX,
    },
    {
        "name": "box_minus_cyl",
        "build": build_box_minus_cyl,
        "description": "40 mm cube minus an axial through-hole cylinder of radius 8 mm.",
        "feature": "plane-cylinder trims: the hole's rim on each cap is an exact circle in the "
                   "plane's 2D chart and a full-turn iso-line in the cylinder's (phi, h) chart",
        "volume_mm3": 40.0 ** 3 - math.pi * 8.0 ** 2 * 40.0,
    },
    {
        "name": "cyl_cross_cyl",
        "build": build_cyl_cross_cyl,
        "description": "Two r = 10 mm, L = 60 mm cylinders with orthogonal axes (z and x), "
                       "fused through a common centre.",
        "feature": "THE key fixture: the union boundary contains the cylinder-cylinder "
                   "intersection curve h(phi) = +- sqrt(r^2 - R^2 sin^2 phi), which is "
                   "transcendental in each cylinder's own (phi, h) chart and therefore not "
                   "exactly representable by any per-face 2D trim curve",
        "volume_mm3": 2.0 * math.pi * 10.0 ** 2 * 60.0 - _V_STEINMETZ,
    },
    {
        "name": "cyl_inter_cyl",
        "build": build_cyl_inter_cyl,
        "description": "The same two orthogonal r = 10 mm cylinders, intersected (Steinmetz "
                       "solid).",
        "feature": "the entire boundary is the transcendental cylinder-cylinder intersection "
                   "curve: two cylindrical patches, four bi-arc edges, no planar face at all",
        "volume_mm3": _V_STEINMETZ,
    },
    {
        "name": "tube_window",
        "build": build_tube_window,
        "description": "Cylinder r = 15 mm, h = 60 mm (axis z) minus a transverse r = 8 mm "
                       "cylinder (axis x) drilled through it.",
        "feature": "minimized reproducer of the Bagger 'BoomCylinderOuter' failure: unequal-"
                   "radius orthogonal cylinder-cylinder intersection, transcendental in both "
                   "charts; the window rim closes over the tube's seam line",
        # Volume is R^2*pi*h minus the intersection of two unequal orthogonal cylinders, which
        # evaluates to a complete elliptic integral, not an elementary closed form.
        "volume_mm3": None,
    },
    {
        "name": "oblique_cut_cyl",
        "build": build_oblique_cut_cyl,
        "description": "Cylinder r = 12 mm, h = 50 mm cut by a plane inclined at 30 deg to its "
                       "axis (cut via a large rotated box).",
        "feature": "minimized reproducer of the Bagger 'Bucket' failure: the cut face is an "
                   "exact ELLIPSE (semi-axes 12 and 24 mm) -- a planar face whose trim is a "
                   "conic, not an arc; on the cylinder the same edge is a sinusoid in (phi, h)",
        # The oblique plane crosses the lateral surface only (z in [25 - 12*tan60, 25 + 12*tan60]
        # = [4.2, 45.8] mm), so the remaining volume is exactly pi r^2 times the axis height.
        "volume_mm3": math.pi * 12.0 ** 2 * 25.0,
    },
    {
        "name": "cyl_plus_cone",
        "build": build_cyl_plus_cone,
        "description": "Cylinder r = 10 mm, h = 30 mm fused with a coaxial truncated cone "
                       "(r1 = 10 mm, r2 = 5 mm, h = 20 mm) stacked on top.",
        "feature": "cylinder-cone junction across a shared circular rim: two different analytic "
                   "surface types meeting tangent-discontinuously with no intervening planar "
                   "face; the shared edge must be trimmed consistently in both charts",
        "volume_mm3": math.pi * 10.0 ** 2 * 30.0
                      + math.pi * 20.0 / 3.0 * (10.0 ** 2 + 10.0 * 5.0 + 5.0 ** 2),
    },
    {
        "name": "sphere_minus_cyl",
        "build": build_sphere_minus_cyl,
        "description": "Sphere r = 20 mm minus an axial r = 6 mm through-hole (napkin ring).",
        "feature": "sphere-cylinder intersection: on the sphere the rim is a circle of constant "
                   "latitude (an iso-line in (theta, phi)), on the cylinder a full-turn "
                   "iso-line; exercises two curved charts meeting with no planar face",
        # Napkin ring: V = 4/3 pi (R^2 - a^2)^(3/2).
        "volume_mm3": 4.0 / 3.0 * math.pi * (20.0 ** 2 - 6.0 ** 2) ** 1.5,
    },
    {
        "name": "torus_union_cyl",
        "build": build_torus_union_cyl,
        "description": "Torus (R = 25 mm, r = 8 mm) fused with a coaxial cylinder r = 20 mm, "
                       "h = 40 mm passing through its centre (see build_torus_union_cyl for why "
                       "the cylinder is 20 mm and not 10 mm).",
        "feature": "toroidal (quartic) surface trimmed by its junction with a cylinder; being "
                   "coaxial the junction curves are exact circles, so this isolates the torus "
                   "surface/ray-intersection code from transcendental trimming",
        "volume_mm3": None,
    },
]

_BY_NAME = {f["name"]: f for f in FIXTURES}


# -------------------------------
# measurement and export
# -------------------------------

def shape_volume_cm3(shape):
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    return props.Mass() * MM3_TO_CM3


def shape_counts(shape):
    topo = TopologyExplorer(shape)
    return topo.number_of_faces(), topo.number_of_edges(), topo.number_of_solids()


def write_step_mm(shape, path: _Path, product_name: str):
    """Write `shape` to `path` as a STEP file with LENGTH_UNIT = millimetre.

    `product_name` becomes the STEP PRODUCT name, i.e. the XCAF label name a reader (such as
    O2_CADtoTGeo.py) picks up, so the converted volume and every artifact derived from it
    (surfaces_*.bin, facets_*.bin, brep_*.brep) carries the fixture name instead of the
    exporter's default string.
    """
    writer = STEPControl_Writer()
    # NB: the write.step.* static parameters only exist once a STEP writer has been created.
    Interface_Static.SetCVal("write.step.unit", "MM")
    Interface_Static.SetCVal("write.step.product.name", product_name)
    writer.Transfer(shape, STEPControl_AsIs)
    status = writer.Write(str(path))
    if status != 1:  # IFSelect_RetDone
        raise RuntimeError(f"STEP write failed for {path} (status {status})")

    # OCCT appends a process-global counter to the product name ("box 1", "box 2", ...), so the
    # name -- and hence every downstream file name derived from it -- would depend on how many
    # fixtures were written before it in this process, i.e. it would change under --only. Drop
    # the counter so the fixture name is stable. The substitution only shortens a quoted string
    # literal, leaving the STEP file structurally untouched.
    text = path.read_text(encoding="latin-1")
    text = re.sub(rf"'{re.escape(product_name)} \d+'", f"'{product_name}'", text)
    path.write_text(text, encoding="latin-1")


# -------------------------------
# the --transform sweep (Stream E)
# -------------------------------
#
# Every measurement recorded on this branch was taken on a part a few centimetres across, sitting
# on the origin, against kernel constants that are *absolute*: kBVHBoxTolerance 1e-3 cm,
# kWireJoinTolerance 1e-5 cm, kRayTolerance 1e-9. A number produced under those conditions is a
# statement about geometry near the origin at that size unless someone shows otherwise, and the
# cheapest way to show otherwise is to move the ladder and re-measure it. Transforming here rather
# than anywhere downstream is what makes the experiment honest: the STEP is the single source the
# converter, the sidecar, the mesh and the oracle's .brep are all derived from, so there is exactly
# one shape in the pipeline and both the kernel and the oracle see the transformed one. See
# scripts/geometry/Stream_E_Scale.md.


def parse_transform(spec: str):
    """Parse a transform spec into (gp_Trsf, volume scale factor, canonical description).

    Accepted forms (lengths in mm, i.e. STEP model units):
      translate:<dx>,<dy>,<dz>
      scale:<factor>                     uniform scaling about the origin
      <op>;<op>;...                      composition, applied left to right

    The composition exists for the one case the sweep is really aimed at: a metre-scale model with
    sub-millimetre features, i.e. `scale:0.1;translate:0,0,4000`. Shrinking and moving separately
    each probe one absolute constant; doing both at once is the combination a real ALICE geometry
    presents, and nothing else on the ladder reproduces it.
    """
    if ";" in spec:
        total = gp_Trsf()
        volume_scale = 1.0
        descriptions = []
        for part in spec.split(";"):
            if not part.strip():
                continue
            trsf, part_scale, description = parse_transform(part)
            total = trsf.Multiplied(total)  # applied after everything parsed so far
            volume_scale *= part_scale
            descriptions.append(description)
        return total, volume_scale, ";".join(descriptions)
    kind, _, rest = spec.partition(":")
    kind = kind.strip().lower()
    trsf = gp_Trsf()
    if kind == "translate":
        parts = [float(v) for v in rest.split(",")]
        if len(parts) != 3:
            raise ValueError(f"translate needs three components, got {rest!r}")
        trsf.SetTranslation(gp_Pnt(0.0, 0.0, 0.0), gp_Pnt(*parts))
        return trsf, 1.0, f"translate:{parts[0]:g},{parts[1]:g},{parts[2]:g}"
    if kind == "scale":
        factor = float(rest)
        if factor <= 0.0:
            raise ValueError(f"scale factor must be positive, got {factor}")
        trsf.SetScale(gp_Pnt(0.0, 0.0, 0.0), factor)
        return trsf, factor ** 3, f"scale:{factor:g}"
    raise ValueError(f"unknown transform {spec!r}; expected 'translate:x,y,z' or 'scale:f'")


def generate(fixture, outdir: _Path, transform=None):
    name = fixture["name"]
    shape = fixture["build"]()
    volume_scale = 1.0
    transform_desc = None
    if transform is not None:
        trsf, volume_scale, transform_desc = transform
        # copy=True: a scaling gp_Trsf is not a gp_Trsf the topology can share geometry through,
        # and sharing it would leave the untransformed poles in place on the copies.
        shape = BRepBuilderAPI_Transform(shape, trsf, True).Shape()
    step_path = (outdir / f"{name}.step").resolve()
    write_step_mm(shape, step_path, name)

    n_faces, n_edges, n_solids = shape_counts(shape)
    valid = bool(BRepCheck_Analyzer(shape).IsValid())
    occt_volume = shape_volume_cm3(shape)
    expected = fixture["volume_mm3"]
    expected_cm3 = None if expected is None else expected * MM3_TO_CM3 * volume_scale

    entry = {
        "name": name,
        "step": str(step_path),
        "transform": transform_desc,
        "description": fixture["description"],
        "feature": fixture["feature"],
        "units": "mm (STEP) / cm^3 (volumes)",
        "expected_volume_cm3": expected_cm3,
        "occt_volume_cm3": occt_volume,
        "volume_rel_error": (None if expected_cm3 in (None, 0.0)
                             else abs(occt_volume - expected_cm3) / abs(expected_cm3)),
        "n_faces": n_faces,
        "n_edges": n_edges,
        "n_solids": n_solids,
        "valid": valid,
    }
    return entry


def print_summary_line(entry):
    exp = entry["expected_volume_cm3"]
    if exp is None:
        vol = f"V={entry['occt_volume_cm3']:11.5f} cm^3 (no closed form)"
    else:
        vol = (f"V={entry['occt_volume_cm3']:11.5f} cm^3 "
               f"(analytic {exp:11.5f}, rel.err {entry['volume_rel_error']:.2e})")
    print(f"  {entry['name']:<18s} faces={entry['n_faces']:3d} edges={entry['n_edges']:3d} "
          f"solids={entry['n_solids']:2d} {vol}  valid={str(entry['valid']).lower()}")


def print_ladder():
    print("Boolean fixture ladder (increasing difficulty):")
    for i, f in enumerate(FIXTURES, 1):
        exp = f["volume_mm3"]
        exp_s = "n/a" if exp is None else f"{exp * MM3_TO_CM3:.5f} cm^3"
        print(f"{i:3d}. {f['name']}")
        print(f"       {f['description']}")
        print(f"       exercises: {f['feature']}")
        print(f"       analytic volume: {exp_s}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=str(_DEFAULT_OUTDIR),
                    help="Directory for the generated .step files and fixtures.json "
                         "(default: %(default)s)")
    ap.add_argument("--only", default=None,
                    help="Comma-separated fixture names to (re)generate instead of the full "
                         "ladder. The manifest is rewritten from the regenerated subset merged "
                         "with any previously generated entries.")
    ap.add_argument("--list", action="store_true",
                    help="Print the fixture ladder and exit without generating anything.")
    ap.add_argument("--transform", default=None,
                    help="Apply a transform to every fixture before export, for the Stream E "
                         "position/scale sweep: 'translate:dx,dy,dz' (mm) or 'scale:f' (uniform, "
                         "about the origin). Omitted, nothing is applied and the output is "
                         "byte-identical to an untransformed run.")
    args = ap.parse_args()

    if args.list:
        print_ladder()
        return

    if args.only:
        wanted = [n.strip() for n in args.only.split(",") if n.strip()]
        unknown = [n for n in wanted if n not in _BY_NAME]
        if unknown:
            raise SystemExit(f"Unknown fixture name(s): {', '.join(unknown)}\n"
                             f"Known: {', '.join(_BY_NAME)}")
        selected = [_BY_NAME[n] for n in wanted]
    else:
        selected = list(FIXTURES)

    outdir = _Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "fixtures.json"

    previous = {}
    if manifest_path.exists():
        try:
            old = json.loads(manifest_path.read_text())
            previous = {e["name"]: e for e in old.get("fixtures", [])}
        except (ValueError, KeyError):
            previous = {}

    transform = parse_transform(args.transform) if args.transform else None
    suffix = "" if transform is None else f", transform: {transform[2]}"
    print(f"Generating {len(selected)} fixture(s) into {outdir} (STEP unit: MM{suffix})")
    entries = {}
    for fixture in selected:
        entry = generate(fixture, outdir, transform)
        entries[fixture["name"]] = entry
        print_summary_line(entry)

    # Keep entries of fixtures that were not regenerated in this run, as long as their STEP
    # file is still there.
    merged = []
    for f in FIXTURES:
        entry = entries.get(f["name"]) or previous.get(f["name"])
        if entry and _Path(entry["step"]).exists():
            merged.append(entry)

    n_invalid = sum(1 for e in merged if not e["valid"])
    manifest = {
        "version": 1,
        "generator": str(_Path(__file__).resolve()),
        "step_length_unit": "mm",
        "volume_unit": "cm^3",
        "transform": None if transform is None else transform[2],
        "outdir": str(outdir),
        "fixtures": merged,
    }
    manifest_path.write_text(json.dumps(manifest, indent=1))
    print(f"\nWrote {manifest_path} ({len(merged)} fixtures, {n_invalid} invalid)")


if __name__ == "__main__":
    main()
