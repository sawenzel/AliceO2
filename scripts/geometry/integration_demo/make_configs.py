#!/usr/bin/env python3
"""Write the o2-sim external-geometry and detector-list JSON for one representation.

Both CAD models are anchored to `barrel`, which Cave.cxx places at cave (0,-30,0); a module
that is to land on the ALICE origin therefore carries translation y = +30.

Rotation: measured, not assumed. Both models are authored with their long axis along the CAD
*y* axis (IRIS: dy = 838 cm against dx = 22, dz = 43; its beam-pipe volume SOLID is centred on
CAD x = z = 0 and spans y in [-390,+390]). TGeoCombiTrans::RotateX(+90) maps local +y -> master
+z and local +z -> master -y, so rotation_deg [90,0,0] puts the CAD long axis on the beam and
the CAD beam axis on the ALICE beam axis.

Usage: make_configs.py <conv-root> <exact|tess> <output-dir>
"""
import json
import os
import sys

# ALICE-frame placement of each model, expressed as the barrel-frame translation that a
# rotation_deg of [90,0,0] needs. IRIS: the CAD origin is the interaction point, so the
# translation is just the barrel offset. Bagger: additionally moved so that its bounding-box
# centre lands at ALICE (100, 0, 0) -- clear of IRIS (which reaches |x| < 13 cm), well inside
# the barrel (r = 790.5 cm), and close enough to the origin that a box generator reaches it.
PLACEMENT = {
    "IRIS": {"translation": [0.0, 30.0, 0.0], "rotation_deg": [90.0, 0.0, 0.0]},
    "BAGR": {"translation": [120.928, 102.064, 20.327], "rotation_deg": [90.0, 0.0, 0.0]},
}


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    conv_root, rep, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    # "exact" and "tess" are the two representations of record; any other value is allowed so
    # that control variants (a coarser mesh, say) can be run through the same pipeline, as long
    # as conv/iris_<rep> and conv/bagger_<rep> exist.
    os.makedirs(outdir, exist_ok=True)

    iris_macro = os.path.abspath(os.path.join(conv_root, "conv", f"iris_{rep}", "geom.C"))
    bagr_macro = os.path.abspath(os.path.join(conv_root, "conv", f"bagger_{rep}", "geom.C"))
    for m in (iris_macro, bagr_macro):
        if not os.path.exists(m):
            print(f"missing macro {m}")
            return 1

    ext = {
        "externalDetectors": [
            {
                "name": "IRIS",
                "title": f"ALICE 3 IRIS vertex-detector CAD ({rep})",
                "macro": iris_macro,
                "anchor": "barrel",
                "detID": "TST",
                # The model carries no silicon: its BOM is mechanical. The single carbon-fibre
                # volume ST2487721_01 (CENTRAL PIPE-IRIS-2ND VACUUM) is the innermost part and
                # is one of the 20 that convert to an exact O2BVHSurfaceSolid, so it exercises
                # the exact kernel in the sensitive path.
                "sensitiveMedia": ["Carbon Fiber"],
                "placement": PLACEMENT["IRIS"],
            },
            {
                "name": "BAGR",
                "title": f"Excavator, Bucket sensitive ({rep})",
                "macro": bagr_macro,
                "anchor": "barrel",
                "detID": "FOC",
                # Substring match: this selects Bucket, BucketLink1, BucketLink2,
                # BucketCylinderInner and BucketCylinderOuter -- the whole bucket group.
                "sensitiveVolumes": ["Bucket"],
                "placement": PLACEMENT["BAGR"],
            },
        ]
    }
    detlist = {"EXTCAD": ["IRIS", "BAGR"]}

    with open(os.path.join(outdir, "externalGeometry.json"), "w") as f:
        json.dump(ext, f, indent=2)
        f.write("\n")
    with open(os.path.join(outdir, "detectorlist.json"), "w") as f:
        json.dump(detlist, f, indent=2)
        f.write("\n")
    print(f"wrote {outdir}/externalGeometry.json and {outdir}/detectorlist.json ({rep})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
