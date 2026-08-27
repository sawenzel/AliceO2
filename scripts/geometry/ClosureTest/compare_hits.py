#!/usr/bin/env python3
"""Compare the hit positions of two o2-sim runs, hit by hit.

The closure test runs the same events through two geometries: the hand-written
C++ TGeo baseline, and the same detector round-tripped through STEP and back.
The two sides do not produce the same hit *class* -- native ITS writes
o2::itsmft::Hit with a chip id resolved through GeometryTGeo, an external
sensitive detector writes o2::ext::Hit with its own running sensor index -- so
the detector id is deliberately never used for matching.

Hits are keyed by (event, trackID) and compared in the order the track made
them, which is well defined under SimCutParams.trackSeed=true because each
track then draws from its own seeded stream.

Usage:
  compare_hits.py A/ B/ --branch-a ITSHit --branch-b ITSHit \
      --file-a o2sim_HitsITS.root --file-b o2sim.root [--json out.json]
"""

import argparse
import json
import math
import os
import sys


def load_hits(rundir, filename, branch):
    """Return {(event, trackID): [(x, y, z), ...]} in the order the hits appear."""
    import ROOT

    path = os.path.join(rundir, filename)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {path}")
    tree = f.Get("o2sim")
    if not tree:
        raise SystemExit(f"no 'o2sim' tree in {path}")
    if not tree.GetBranch(branch):
        have = [b.GetName() for b in tree.GetListOfBranches()]
        raise SystemExit(f"no branch {branch!r} in {path}; have {have}")

    hits = {}
    total = 0
    for iev in range(tree.GetEntries()):
        tree.GetEntry(iev)
        vec = getattr(tree, branch)
        for i in range(vec.size()):
            h = vec.at(i)
            hits.setdefault((iev, h.GetTrackID()), []).append(
                (h.GetX(), h.GetY(), h.GetZ()))
            total += 1
    f.Close()
    return hits, total


def compare(a, b, tol):
    """Compare two keyed hit maps. Returns a result dict."""
    keys_a, keys_b = set(a), set(b)
    common = keys_a & keys_b

    res = {
        "tracksOnlyInA": len(keys_a - keys_b),
        "tracksOnlyInB": len(keys_b - keys_a),
        "tracksCommon": len(common),
        "tracksWithDifferentHitCount": 0,
        "hitsCompared": 0,
        "hitsWithinTolerance": 0,
        "maxDx": 0.0, "maxDy": 0.0, "maxDz": 0.0, "maxDr": 0.0,
        "sumDr": 0.0,
        "worst": None,
    }

    for key in sorted(common):
        ha, hb = a[key], b[key]
        if len(ha) != len(hb):
            res["tracksWithDifferentHitCount"] += 1
        for (xa, ya, za), (xb, yb, zb) in zip(ha, hb):
            dx, dy, dz = abs(xa - xb), abs(ya - yb), abs(za - zb)
            dr = math.sqrt(dx * dx + dy * dy + dz * dz)
            res["hitsCompared"] += 1
            res["sumDr"] += dr
            if dr <= tol:
                res["hitsWithinTolerance"] += 1
            res["maxDx"] = max(res["maxDx"], dx)
            res["maxDy"] = max(res["maxDy"], dy)
            res["maxDz"] = max(res["maxDz"], dz)
            if dr > res["maxDr"]:
                res["maxDr"] = dr
                res["worst"] = {
                    "event": key[0], "trackID": key[1],
                    "a": [xa, ya, za], "b": [xb, yb, zb], "dr": dr,
                }

    n = res["hitsCompared"]
    res["meanDr"] = res["sumDr"] / n if n else 0.0
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dir_a")
    p.add_argument("dir_b")
    p.add_argument("--file-a", default="o2sim_HitsITS.root")
    p.add_argument("--file-b", default="o2sim_HitsITS.root")
    p.add_argument("--branch-a", default="ITSHit")
    p.add_argument("--branch-b", default="ITSHit")
    p.add_argument("--tol", type=float, default=0.0,
                   help="position tolerance in cm; 0 means require exact equality")
    p.add_argument("--json", help="also write the result as JSON here")
    args = p.parse_args()

    a, na = load_hits(args.dir_a, args.file_a, args.branch_a)
    b, nb = load_hits(args.dir_b, args.file_b, args.branch_b)

    res = compare(a, b, args.tol)
    res["hitsInA"] = na
    res["hitsInB"] = nb
    res["tolerance"] = args.tol

    print(f"A: {args.dir_a}/{args.file_a}:{args.branch_a}  {na} hits, {len(a)} tracks")
    print(f"B: {args.dir_b}/{args.file_b}:{args.branch_b}  {nb} hits, {len(b)} tracks")
    print(f"tracks: {res['tracksCommon']} common, "
          f"{res['tracksOnlyInA']} only in A, {res['tracksOnlyInB']} only in B, "
          f"{res['tracksWithDifferentHitCount']} with a different hit count")
    print(f"hits compared: {res['hitsCompared']}, "
          f"within {args.tol} cm: {res['hitsWithinTolerance']}")
    print(f"max |dx| {res['maxDx']:.6g}  max |dy| {res['maxDy']:.6g}  "
          f"max |dz| {res['maxDz']:.6g}  cm")
    print(f"max |dr| {res['maxDr']:.6g} cm, mean |dr| {res['meanDr']:.6g} cm")
    if res["worst"]:
        w = res["worst"]
        print(f"worst: event {w['event']} track {w['trackID']} "
              f"A={w['a']} B={w['b']}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=2)
        print(f"wrote {args.json}")

    identical = (res["hitsInA"] == res["hitsInB"]
                 and res["tracksOnlyInA"] == 0 and res["tracksOnlyInB"] == 0
                 and res["tracksWithDifferentHitCount"] == 0
                 and res["hitsCompared"] == res["hitsWithinTolerance"])
    print("VERDICT:", "identical within tolerance" if identical else "DIFFERENT")
    return 0 if identical else 1


if __name__ == "__main__":
    sys.exit(main())
