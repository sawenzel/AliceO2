#!/usr/bin/env python3
"""Find where a module hangs itself in the ALICE world, and with what matrix.

o2-sim always builds the experiment hall, and only the hall: `cave`, `barrel`
(at y = -30 in cave) and `caveRB24` come from o2::passive::Cave whatever module
list is requested.  Everything else under them is the module's own contribution.

For the closure test that distinction decides the anchoring.  Converting a
module's geometry from `cave` downward would bring three air volumes along that
the CAD run builds natively anyway, so each module is converted from its own
subtree root instead and anchored where the source put it -- same volume, same
matrix.  This script reports those roots, and is the placement-fidelity walk
that has to run BEFORE any transport: a module that does not land where the
original detector sat invalidates everything downstream, silently.

Usage:
  module_anchors.py o2sim_geometry.root [--json anchors.json]
"""

import argparse
import json
import sys

HALL = ("cave", "barrel", "caveRB24")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("geometry")
    p.add_argument("--json")
    args = p.parse_args()

    import ROOT
    ROOT.gROOT.SetBatch(True)
    mgr = ROOT.TGeoManager.Import(args.geometry)
    if not mgr:
        raise SystemExit(f"cannot read {args.geometry}")

    roots = []
    for hall in HALL:
        vol = mgr.GetVolume(hall)
        if not vol:
            continue
        for i in range(vol.GetNdaughters()):
            node = vol.GetNode(i)
            child = node.GetVolume()
            if str(child.GetName()) in HALL:
                continue
            m = node.GetMatrix()
            t = [m.GetTranslation()[k] for k in range(3)]
            r = [m.GetRotationMatrix()[k] for k in range(9)]
            box = child.GetShape()
            roots.append({
                "anchor": hall,
                "volume": str(child.GetName()),
                "node": str(node.GetName()),
                "copy": int(node.GetNumber()),
                "shape": str(box.ClassName()),
                "isAssembly": bool(child.IsAssembly()),
                "nDaughters": int(child.GetNdaughters()),
                "translation": t,
                "rotation": r,
                "isIdentity": bool(m.IsIdentity()),
            })

    print(f"{args.geometry}: {mgr.GetListOfVolumes().GetEntries()} volumes, "
          f"{len(roots)} subtree root(s) outside the hall")
    for r in roots:
        rot = "identity" if r["isIdentity"] else f"rotation {['%.6g' % x for x in r['rotation']]}"
        print(f"  {r['volume']:24s} in {r['anchor']:9s} copy {r['copy']:<4d} "
              f"{r['shape']:22s} nd={r['nDaughters']:<5d} "
              f"t={['%.6g' % x for x in r['translation']]} {rot}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"geometry": args.geometry, "roots": roots}, fh, indent=2)
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
