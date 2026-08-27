#!/usr/bin/env python3
"""Write the o2-sim configuration that runs the round-tripped geometry.

The CAD side of the closure test replaces PIPE, ITS, TPC and MAG with their
round-tripped selves.  Three of them are passive and go in as
o2::passive::ExternalModule; ITS is the hit detector and goes in as
o2::ext::ExternalDetector on the ITS DetID slot, so its hits are written to
o2sim_HitsITS.root like any built-in detector's.

Anchoring.  Each module was converted from its source geometry's top volume with
the experiment hall hollowed out (--hollow-volume cave barrel caveRB24), so the
converted macro reproduces the hall's *structure* -- and therefore every
subtree's transform relative to it -- while contributing none of its material.
Anchoring at `cave` with no placement is then exact by construction: nothing has
to be re-derived, and there is no second copy of the hall for the navigator to
resolve against the real one.

The external names are deliberately NOT the real module names: build_geometry.C
activates a module when its name is in the active list, so calling the external
one `PIPE` would build the native beam pipe as well.

Usage:
  make_configs.py <studydir> [--name CADCLOSURE]
"""

import argparse
import json
import os
import sys

# module -> (external name, kind).  ITS is the only sensitive one: TPC's
# sensitivity is per-electron ionisation, which a generic entrance/exit action
# cannot stand in for, so it goes in as passive material exactly as the charge
# asks.
MODULES = [
    ("PIPE", "CPIPE", "passive", "CAD round-tripped beam pipe"),
    ("TPC",  "CTPC",  "passive", "CAD round-tripped TPC (material only)"),
    ("MAG",  "CMAG",  "passive", "CAD round-tripped L3 magnet"),
    ("ITS",  "CITS",  "sensitive", "CAD round-tripped ITS"),
]

SENSITIVE_VOLUMES = {"CITS": ["ITSUSensor"]}   # substring match: ITSUSensor0..6
DET_ID = {"CITS": "ITS"}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("studydir")
    p.add_argument("--name", default="CADCLOSURE",
                   help="the detector-list key o2-sim is pointed at")
    args = p.parse_args()

    study = os.path.abspath(args.studydir)
    modules, detectors, names, missing = [], [], [], []

    for mod, name, kind, title in MODULES:
        macro = os.path.join(study, "cad", mod, "conv", "geom.C")
        if not os.path.exists(macro):
            missing.append(macro)
            continue
        entry = {"name": name, "title": title, "macro": macro, "anchor": "cave"}
        if kind == "sensitive":
            entry["detID"] = DET_ID[name]
            entry["sensitiveVolumes"] = SENSITIVE_VOLUMES[name]
            detectors.append(entry)
        else:
            modules.append(entry)
        names.append(name)

    if missing:
        raise SystemExit("no converted macro for:\n  " + "\n  ".join(missing))

    ext_path = os.path.join(study, "externalDetectors.json")
    det_path = os.path.join(study, "detectorlist.json")
    with open(ext_path, "w") as fh:
        json.dump({"externalModules": modules, "externalDetectors": detectors},
                  fh, indent=2)
    with open(det_path, "w") as fh:
        json.dump({args.name: names}, fh, indent=2)

    print(f"wrote {ext_path}")
    print(f"  {len(modules)} passive external module(s): "
          f"{', '.join(e['name'] for e in modules)}")
    print(f"  {len(detectors)} sensitive external detector(s): "
          + ", ".join(f"{e['name']} on DetID {e['detID']} "
                      f"(sensitive: {', '.join(e['sensitiveVolumes'])})"
                      for e in detectors))
    print(f"wrote {det_path}: {args.name} = {names}")
    print()
    print("run it with:")
    print(f"  o2-sim-serial -n <N> -g boxgen \\")
    print(f"      --detectorList {args.name}:{det_path} \\")
    print(f"      --extGeomFile {ext_path} \\")
    print(f"      --seed <SEED> --configKeyValues 'SimCutParams.trackSeed=true'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
