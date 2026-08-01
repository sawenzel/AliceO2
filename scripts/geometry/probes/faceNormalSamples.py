#!/usr/bin/env python3
"""Stream L, stage 1: sample interior points of every face of a solid and record OCCT's ORIENTED
outward normal there.

The face index is `TopologyExplorer(shape).faces()` order -- the exact iteration
`extract_surfaces_for_shape` uses -- so face i is sidecar surface record i.

Output feeds `faceNormals.cxx`, which asks O2BVHSurfaceSolid::ComputeNormal(point, nullptr, n) at
the same points. A face whose two normals are antiparallel has its outward direction inverted in
the sidecar, which is invisible to every closure and edge-identity check (both are insensitive to
a global sign) and fatal to DistFromOutside/DistFromInside (both select hits by that sign).
"""
import argparse
import json
import math
import sys
from pathlib import Path

from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepTools import breptools
from OCC.Core.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.TopAbs import TopAbs_IN, TopAbs_REVERSED
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.gp import gp_Pnt2d
from OCC.Extend.TopologyUtils import TopologyExplorer

_TYPE = {0: "plane", 1: "cylinder", 2: "cone", 3: "sphere", 4: "torus", 5: "bezier",
         6: "bspline", 7: "revolution", 8: "extrusion", 9: "offset", 10: "other"}
_CURVE = {0: "line", 1: "circle", 2: "ellipse", 3: "hyperbola", 4: "parabola", 5: "bezier",
          6: "bspline", 7: "other"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--brep", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--grid", type=int, default=7)
    ap.add_argument("--samples-per-face", type=int, default=4)
    args = ap.parse_args()

    shape = TopoDS_Shape()
    breptools.Read(shape, str(args.brep), BRep_Builder())
    faces = list(TopologyExplorer(shape).faces())

    doc = {"brep": str(args.brep), "nFaces": len(faces), "faces": []}
    for idx, face in enumerate(faces):
        ad = BRepAdaptor_Surface(face)
        kind = _TYPE.get(ad.GetType(), "?")
        rev = face.Orientation() == TopAbs_REVERSED
        umin, umax, vmin, vmax = breptools.UVBounds(face)
        classifier = BRepTopAdaptor_FClass2d(face, 1.0e-9)
        surf = BRep_Tool.Surface(face)
        pts = []
        n = args.grid
        for i in range(1, n):
            for j in range(1, n):
                u = umin + (umax - umin) * i / n
                v = vmin + (vmax - vmin) * j / n
                if classifier.Perform(gp_Pnt2d(u, v)) != TopAbs_IN:
                    continue
                props = GeomLProp_SLProps(surf, u, v, 1, 1.0e-9)
                if not props.IsNormalDefined():
                    continue
                nn = props.Normal()
                vec = [nn.X(), nn.Y(), nn.Z()]
                if rev:
                    vec = [-c for c in vec]
                p = props.Value()
                # cm: the converter scales STEP mm -> cm.  The .brep the gate dumps is already cm.
                pts.append({"uv": [u, v], "p": [p.X(), p.Y(), p.Z()], "n": vec})
                if len(pts) >= args.samples_per_face:
                    break
            if len(pts) >= args.samples_per_face:
                break
        edges = []
        for e in TopologyExplorer(face).edges():
            from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
            edges.append({"kind": _CURVE.get(BRepAdaptor_Curve(e).GetType(), "?"),
                          "tol": BRep_Tool.Tolerance(e),
                          "degenerate": bool(BRep_Tool.Degenerated(e))})
        doc["faces"].append({"index": idx, "type": kind, "reversed": rev,
                             "faceTol": BRep_Tool.Tolerance(face),
                             "nEdges": len(edges), "edges": edges, "samples": pts})
    args.out.write_text(json.dumps(doc))
    print(f"{args.brep.name}: {len(faces)} faces, "
          f"{sum(len(f['samples']) for f in doc['faces'])} samples -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
