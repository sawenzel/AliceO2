#!/usr/bin/env python3
"""Re-tessellate a `brep_<part>.brep` at a stated linear AND angular deflection.

Why this exists
---------------
`O2_CADtoTGeo.py --mesh-prec X` sets **both** `lin_defl` and `ang_defl` of
`BRepMesh_IncrementalMesh` to the same number X. They do not have the same units and they are not
the same knob:

  * `lin_defl` is a **length**, in the model's own units, and it bounds the chordal sagitta.
  * `ang_defl` is an **angle in radians**, and it bounds the turn between consecutive chords.

So `--mesh-prec 0.5` asks for a half-unit chord tolerance *and* an angular tolerance of 0.5 rad =
**28.6 degrees**, which on a cylinder is a 12-sided prism no matter how fine the linear term is.
The coupled flag therefore cannot separate "the mesh is too coarse because of the chord bound"
from "the mesh is too coarse because of the angle bound", and the cost/accuracy curve it produces
is a curve in one variable through a two-dimensional space.

This tool takes the two apart. It reads the `.brep` the converter already dumps (`--dump-brep`),
meshes it with an explicit `--lin` and `--ang`, and writes a `facets_<part>.bin` in exactly the
converter's format, so every downstream consumer -- `O2Tessellated`, the harness, the oracle gate,
the X-ray benchmark -- reads it without knowing where it came from.

Units, and the one thing to get right
-------------------------------------
The converter writes `brep_*.brep` **already scaled to cm** (see `write_brep_cm`), so `--lin` here
is in **cm**, whereas `--mesh-prec` is in the STEP file's own units. For a model authored in mm
those differ by a factor of 10 and the two sweeps are not directly comparable unless that is
stated. It is stated in the output of this script, on every run.

Usage
-----
    remeshFromBrep.py --brep out/brep_Bagger_Boom.brep --lin 0.05 --ang 0.2 \\
                      --out out/facets_Bagger_Boom.bin

    # a whole directory, one output directory per (lin, ang) pair
    remeshFromBrep.py --brep-dir out --lin 0.05 --ang 0.2 --out-dir sweep/lin0.05_ang0.2

Needs pythonOCC; run it under the alibuild python3.10 with LD_LIBRARY_PATH / PYTHONPATH
**prepended** (never replaced -- replacing kills PyROOT on libffi.so.6).
"""

import argparse
import struct
import sys
import time
from pathlib import Path

from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepTools import breptools
from OCC.Core.TopAbs import TopAbs_REVERSED
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Extend.TopologyUtils import TopologyExplorer


def read_brep(path: Path) -> TopoDS_Shape:
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    if not breptools.Read(shape, str(path), builder):
        raise RuntimeError(f"cannot read {path}")
    return shape


def triangulate(shape, lin_defl: float, ang_defl: float):
    """Mesh `shape` and return its triangles, in the shape's own units.

    Deliberately the same traversal the converter's `triangulate_CAD_solid` performs, including
    the `TopAbs_REVERSED` vertex swap: a mesh produced by a different winding convention would
    have every outward normal inverted and would look like a kernel defect rather than like a
    different tool.
    """
    try:
        BRepMesh_IncrementalMesh(shape, lin_defl, False, ang_defl, True)
    except TypeError:
        BRepMesh_IncrementalMesh(shape, lin_defl, False, ang_defl)
    triangles = []
    for face in TopologyExplorer(shape).faces():
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation(face, loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        reverse = face.Orientation() == TopAbs_REVERSED
        for i in range(1, tri.NbTriangles() + 1):
            n1, n2, n3 = tri.Triangle(i).Get()
            p1 = tri.Node(n1).Transformed(trsf)
            p2 = tri.Node(n2).Transformed(trsf)
            p3 = tri.Node(n3).Transformed(trsf)
            if reverse:
                p2, p3 = p3, p2
            triangles.append(((p1.X(), p1.Y(), p1.Z()),
                              (p2.X(), p2.Y(), p2.Z()),
                              (p3.X(), p3.Y(), p3.Z())))
    return triangles


def write_facets_bin(path: Path, triangles):
    """The converter's format, byte for byte: uint32 count then 9 x float32 per triangle."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(struct.pack("<I", len(triangles)))
        for (a, b, c) in triangles:
            handle.write(struct.pack("<9f", *a, *b, *c))


def part_stem(brep_path: Path) -> str:
    name = brep_path.name
    if not name.startswith("brep_") or not name.endswith(".brep"):
        raise RuntimeError(f"{brep_path}: not a brep_<part>.brep")
    return name[len("brep_"):-len(".brep")]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--brep", type=Path, default=None, help="one brep_<part>.brep")
    parser.add_argument("--brep-dir", type=Path, default=None,
                        help="a directory of brep_*.brep; meshes all of them")
    parser.add_argument("--lin", type=float, required=True,
                        help="linear deflection, in CM (the .brep is written in cm)")
    parser.add_argument("--ang", type=float, required=True,
                        help="angular deflection, in RADIANS. 0.5 rad = 28.6 degrees -- the value "
                             "--mesh-prec 0.5 silently also asks for.")
    parser.add_argument("--out", type=Path, default=None, help="output facets_<part>.bin (--brep)")
    parser.add_argument("--out-dir", type=Path, default=None, help="output directory (--brep-dir)")
    args = parser.parse_args()

    if (args.brep is None) == (args.brep_dir is None):
        parser.error("give exactly one of --brep or --brep-dir")

    breps = [args.brep] if args.brep else sorted(args.brep_dir.glob("brep_*.brep"))
    if not breps:
        print("no brep_*.brep found", file=sys.stderr)
        return 1

    print(f"lin_defl = {args.lin} cm, ang_defl = {args.ang} rad "
          f"({args.ang * 180.0 / 3.141592653589793:.1f} degrees); "
          f"the .brep is already in cm, so --lin is a cm and NOT the converter's --mesh-prec")
    total_triangles = 0
    total_seconds = 0.0
    for brep_path in breps:
        stem = part_stem(brep_path)
        out = args.out if args.out else (args.out_dir / f"facets_{stem}.bin")
        started = time.time()
        triangles = triangulate(read_brep(brep_path), args.lin, args.ang)
        elapsed = time.time() - started
        write_facets_bin(out, triangles)
        total_triangles += len(triangles)
        total_seconds += elapsed
        print(f"  {stem:<40} {len(triangles):>9} triangles  {elapsed:7.2f} s  -> {out}")
    print(f"total: {len(breps)} part(s), {total_triangles} triangles, {total_seconds:.2f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
