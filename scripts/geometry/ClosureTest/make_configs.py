#!/usr/bin/env python3
"""Write the o2-sim configuration that runs the round-tripped geometry.

The CAD side of the closure test replaces PIPE, ITS, TPC and MAG with their
round-tripped selves.  Three of them are passive and go in as
o2::passive::ExternalModule; ITS is the hit detector and goes in as
o2::ext::ExternalDetector on the ITS DetID slot, so its hits are written to
o2sim_HitsITS.root like any built-in detector's.

Anchoring.  A module can need more than one placement, and roundtrip_module.py
has already worked out which: everything the source hung under `barrel` is one
conversion anchored back into the real `barrel`, and a subtree hung off `cave` or
`caveRB24` is converted from its own root and placed with that root's own matrix.
Anchoring the whole module at `cave` instead does NOT work -- its detectors would
sit geometrically inside the native `barrel` while being daughters of `cave`, and
the navigator, having descended into `barrel`, never finds them.

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
    p.add_argument("--variant", default="csg", choices=("csg", "mesh"),
                   help="which back-conversion to configure: the shipped cascade "
                        "(csg) or tessellated-only (mesh), the fallback every other "
                        "CAD pipeline uses and the benchmark for the exact path")
    p.add_argument("--out-prefix", default="",
                   help="prefix for the written JSON file names, so two variants "
                        "can live side by side in one study directory")
    args = p.parse_args()

    study = os.path.abspath(args.studydir)
    modules, detectors, names, missing = [], [], [], []

    for mod, name, kind, title in MODULES:
        frag = os.path.join(study, "cad", mod, "module_entries.json")
        if not os.path.exists(frag):
            missing.append(frag)
            continue
        # One module can need several placements: everything under `barrel` goes
        # in as one piece, and a subtree the source hung off `cave` or `caveRB24`
        # goes in separately, at its own matrix. Each becomes its own external
        # module or detector, named <name> or <name>_<tag>.
        frag_entries = json.load(open(frag))["entries"]
        if isinstance(frag_entries, dict):          # variants
            if args.variant not in frag_entries:
                missing.append(f"{frag} (no '{args.variant}' variant)")
                continue
            frag_entries = frag_entries[args.variant]
        for e in frag_entries:
            suffix = "" if e["tag"] == "barrel" else "_" + e["tag"][:8].upper()
            ename = (name + suffix)[:15]
            entry = {"name": ename, "title": f"{title} [{e['tag']}]",
                     "macro": e["macro"], "anchor": e["anchor"]}
            if e.get("placement"):
                entry["placement"] = e["placement"]
            if kind == "sensitive":
                entry["detID"] = DET_ID[name]
                entry["sensitiveVolumes"] = SENSITIVE_VOLUMES[name]
                detectors.append(entry)
            else:
                modules.append(entry)
            names.append(ename)

    if missing:
        raise SystemExit("no module_entries.json for:\n  " + "\n  ".join(missing)
                         + "\n(run roundtrip_module.py for each module first)")

    pre = args.out_prefix
    ext_path = os.path.join(study, f"{pre}externalDetectors.json")
    det_path = os.path.join(study, f"{pre}detectorlist.json")
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
