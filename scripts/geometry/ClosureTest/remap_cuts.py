#!/usr/bin/env python3
"""Carry the baseline's Geant cuts and processes over to the CAD run, by name.

MaterialManager::loadCutsAndProcessesFromJSON resolves each entry with
getMediumID(module, local_id) -- so a dump is only loadable by a run that
assigns the same local indices under the same module names.  The CAD run does
neither: its modules are CPIPE/CITS/CTPC/CMAG, and remapCADMedia numbers their
media in its own traversal order.  Renaming the module keys is therefore not
enough, and a mismatched entry is skipped *silently*, which would leave the two
runs on different physics with nothing in the log to say so.

What both sides do preserve exactly is the medium NAME -- the round trip carries
it verbatim in the media sidecar, and MaterialManager prefixes it with the
module name on each side.  So the mapping is built from names:

  baseline   module ITS,   local 1, medium `ITS_AIR$`
  CAD run    module CITS,  local 7, medium `CITS_ITS_AIR$`

Strip the CAD module prefix, match on what is left, and write the baseline's
cuts under the CAD run's own module and local index.  This needs the CAD run's
own dump first, which is why the closure test runs the CAD side twice: once to
learn its indices, once with the cuts loaded.

  remap_cuts.py --baseline base.json --cad-dump cad_out.json --out cad_in.json
  remap_cuts.py --compare  base.json cad_out2.json
"""

import argparse
import json
import sys


def index_baseline(doc):
    """medium name (unprefixed by module) -> its cuts/processes record."""
    out = {}
    for key, entries in doc.items():
        if not isinstance(entries, list):
            continue
        for e in entries:
            name = e.get("medium_name")
            if not name:
                continue
            # `ITS_AIR$` under module `ITS` -> key on both the full name and the
            # part after the module prefix, so either spelling matches later.
            out.setdefault(name, e)
            if name.startswith(key + "_"):
                out.setdefault(name[len(key) + 1:], e)
    return out


def strip_module(name, module):
    return name[len(module) + 1:] if name.startswith(module + "_") else name


def build(baseline, caddump):
    by_name = index_baseline(baseline)
    out, matched, unmatched = {}, [], []
    for key, entries in caddump.items():
        if not isinstance(entries, list):
            out[key] = entries          # default / enableSpecial* pass through
            continue
        rebuilt = []
        for e in entries:
            bare = strip_module(e.get("medium_name", ""), key)
            src = by_name.get(bare) or by_name.get(e.get("medium_name", ""))
            if src is None:
                unmatched.append(f"{key}/{e.get('medium_name')}")
                continue
            rebuilt.append({
                "local_id": e["local_id"],
                "global_id": e["global_id"],
                "medium_name": e["medium_name"],
                "material_name": e.get("material_name"),
                "cuts": src.get("cuts", {}),
                "processes": src.get("processes", {}),
            })
            matched.append(f"{key}/{e.get('medium_name')} <- {bare}")
        out[key] = rebuilt
    for k in ("default", "enableSpecialCuts", "enableSpecialProcesses"):
        if k in baseline:
            out[k] = baseline[k]
    return out, matched, unmatched


def compare(baseline, caddump):
    """Do the two runs give every medium the same cuts and processes?"""
    by_name = index_baseline(baseline)
    same, differ, missing = 0, [], []
    for key, entries in caddump.items():
        if not isinstance(entries, list):
            continue
        for e in entries:
            bare = strip_module(e.get("medium_name", ""), key)
            src = by_name.get(bare) or by_name.get(e.get("medium_name", ""))
            if src is None:
                missing.append(f"{key}/{e.get('medium_name')}")
                continue
            if (src.get("cuts") == e.get("cuts")
                    and src.get("processes") == e.get("processes")):
                same += 1
            else:
                bad = [k for k in set(src.get("cuts", {})) | set(e.get("cuts", {}))
                       if src.get("cuts", {}).get(k) != e.get("cuts", {}).get(k)]
                bad += [f"proc:{k}" for k in
                        set(src.get("processes", {})) | set(e.get("processes", {}))
                        if src.get("processes", {}).get(k) != e.get("processes", {}).get(k)]
                differ.append((f"{key}/{e.get('medium_name')}", bad))
    return same, differ, missing


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--baseline", required=True)
    p.add_argument("--cad-dump", required=True)
    p.add_argument("--out")
    p.add_argument("--compare", action="store_true")
    args = p.parse_args()

    baseline = json.load(open(args.baseline))
    caddump = json.load(open(args.cad_dump))

    if args.compare:
        same, differ, missing = compare(baseline, caddump)
        print(f"media compared: {same + len(differ) + len(missing)}")
        print(f"  identical cuts and processes: {same}")
        print(f"  differing: {len(differ)}")
        for name, bad in differ[:8]:
            print(f"    {name}: {bad[:6]}")
        print(f"  no baseline medium of that name: {len(missing)}"
              + (f"  e.g. {missing[:5]}" if missing else ""))
        ok = not differ and not missing
        print("VERDICT:", "both runs give every medium the same cuts and processes"
              if ok else "DIFFERENT -- do not trust a transport comparison")
        return 0 if ok else 1

    if not args.out:
        raise SystemExit("--out is required unless --compare is given")
    out, matched, unmatched = build(baseline, caddump)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {args.out}: {len(matched)} medium/media matched by name, "
          f"{len(unmatched)} unmatched")
    if unmatched:
        print(f"  [WARN] no baseline cuts for: {unmatched[:8]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
