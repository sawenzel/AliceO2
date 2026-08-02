#!/usr/bin/env python3
"""Stream N reporting probe: turn `trimEdgeCensus.py --json` output into the tables in
`../Stream_O_ImplicitTrims.md`.

Separate from the census because the census is expensive (it reloads the STEP file) and the
tables were re-cut several times. Pure stdlib + numpy; needs no OCC, so it runs under any
interpreter that has numpy.

    python3 probes/trimEdgeCensusReport.py /tmp/n/alice3.json [/tmp/n/iris.json ...]
    python3 probes/trimEdgeCensusReport.py --compare /tmp/n/alice3_legacy.json /tmp/n/alice3.json

`--compare LEGACY CURRENT` prints the reconciliation against `Stream_K_Tier0.md` §2's
1891 / 834 / 4 / 1053 table: that table was measured with the pre-fix recogniser and over a
denominator that excludes degenerate edges, and this reproduces both readings side by side.
"""

import argparse
import collections
import json
import sys

import numpy as np


def load(paths):
    out = []
    for p in paths:
        for res in json.load(open(p)):
            res["_path"] = p
            out.append(res)
    return out


def pct(values, qs=(0, 50, 90, 99, 100)):
    if not values:
        return "-"
    v = np.percentile(np.asarray(values, float), qs)
    return "  ".join(f"{q}%={x:.3g}" for q, x in zip(qs, v))


def edges_of(res):
    return [e for s in res["solids"] for e in s["edges"]]


def report_model(res):
    label = res["label"]
    solids = res["solids"]
    edges = edges_of(res)
    rej = [e for e in edges if not e["accepted"]]
    acc = [e for e in edges if e["accepted"]]
    faces = [f for s in solids for f in s["faces"]]
    print(f"\n=== {label}  ({res['model']})")
    print(f"leaf solids {len(solids)}   emitted {sum(1 for s in solids if s['emitted'])}   "
          f"faces reaching the wire block {len(faces)} "
          f"{dict(collections.Counter(f['recognizedKind'] for f in faces))}")
    print(f"boundary edges on those faces {len(edges)}   iso {len(acc)}   non-iso {len(rej)}   "
          f"verdict mismatches {sum(s['verdictMismatches'] for s in solids)}")
    print(f"buckets {dict(collections.Counter(e['bucket'] for e in rej))}   "
          f"(B threshold {res['summary']['bFactor']}x the declared BRep tolerance)")

    print("\n-- bucket detail")
    for (b, d), n in sorted(collections.Counter((e["bucket"], e["bucketDetail"])
                                                for e in rej).items()):
        print(f"   {b}  {n:5d}  {d}")

    print("\n-- surface-pair frequency, bucket B")
    bb = [e for e in rej if e["bucket"] == "B"]
    for pair, n in collections.Counter(e["pair"] for e in bb).most_common():
        sub = [e for e in bb if e["pair"] == pair]
        print(f"   {pair:20s} {n:5d}   max dev {max(e['devCm'] for e in sub):.3g} cm   "
              f"max dev/tol {max(e['devOverTolRef'] for e in sub):.4g}")

    print("\n-- deviation distribution (cm unless stated)")
    for name, key, src in (("B  dev", "devCm", bb),
                           ("B  dev/edgeLength", "devOverEdgeLength", bb),
                           ("B  dev/patchDiagonal", "devOverPatchDiagonal", bb),
                           ("B  dev/declaredTol", "devOverTolRef", bb),
                           ("ctl dev/declaredTol", "devOverTolRef",
                            [e for e in acc if "devOverTolRef" in e])):
        print(f"   {name:22s} n={len(src):5d}  {pct([e[key] for e in src])}")
    print("   ('ctl' = the edges the shipped test ACCEPTS, measured with the same instrument:")
    print("    rims and generators, which are exactly quadric-meets-cap-plane intersections.)")

    tail = sorted(bb, key=lambda e: -e["devOverTolRef"])[:10]
    print("\n-- the tail of bucket B (worst 10 by dev / declared tolerance)")
    for e in tail:
        print(f"   {e['solid']:16s} f#{e['face']:<4d} {e['pair']:18s} dev {e['devCm']:.3g} cm  "
              f"tol {e['tolRefCm']:.3g} cm  ratio {e['devOverTolRef']:.4g}  "
              f"dev/len {e['devOverEdgeLength']:.3g}")

    print("\n-- per solid (only solids with faces reaching the wire block)")
    print(f"   {'solid':18s} {'nF':>5s} {'emit':>5s} {'recF':>5s} {'badF':>5s} {'edges':>6s} "
          f"{'nonIso':>6s} {'A':>4s} {'B':>4s} {'C':>4s} {'D':>4s} {'unsupF':>6s}")
    for s in solids:
        if not s["faces"]:
            continue
        c = collections.Counter(e["bucket"] for e in s["edges"] if not e["accepted"])
        print(f"   {s['name']:18s} {s['nFaces']:5d} {str(s['emitted']):>5s} {len(s['faces']):5d} "
              f"{sum(1 for f in s['faces'] if f['nRejected'] > 0):5d} {len(s['edges']):6d} "
              f"{sum(c.values()):6d} {c['A']:4d} {c['B']:4d} {c['C']:4d} {c['D']:4d} "
              f"{s['nUnsupportedFaces']:6d}")

    ru = res["summary"]["rollup"]
    print(f"\n-- SOLID-LEVEL ROLLUP: {res['summary']['nEligibleNotEmitted']} solids are "
          f"surface-eligible but emit no sidecar")
    print(f"   fully covered by A+B: {len(ru['fullyCovered'])}   still fails: "
          f"{len(ru['stillFails'])}")
    for e in ru["fullyCovered"]:
        print(f"     covered     {e['solid']:16s} {e['nRejected']:4d} rejected {e['buckets']}")
    for e in ru["stillFails"]:
        print(f"     still fails {e['solid']:16s} uncovered {e['nUncovered']} of "
              f"{e['nRejected']}  {e['buckets']}  {e['otherWireBlockFailures']}")

    pe = [e for s in solids for e in s["planarEdges"]]
    if pe:
        print("\n-- second mechanism: planar faces declined for their trim-curve vocabulary")
        for e in pe:
            print(f"   {e['solid']:16s} f#{e['face']:<4d} {e['curveType']:10s} "
                  f"nbr={e.get('neighbourKind')}({e.get('neighbourSource')})  bucket {e['bucket']}"
                  f"  dev {e.get('devCm', float('nan')):.3g} cm  "
                  f"tol {e.get('tolRefCm', float('nan')):.3g} cm  "
                  f"ratio {e.get('devOverTolRef', float('nan')):.4g}  "
                  f"dev/len {e.get('devOverEdgeLength', float('nan')):.3g}")

    print("\n-- sphere faces: does ONE alternative polar axis rescue the face?")
    sp = [f for f in faces if f.get("sphereAlternativeAxis")]
    print(f"   sphere faces with rejected edges {len(sp)}   rescued by a single axis "
          f"{sum(1 for f in sp if f['sphereAlternativeAxis']['ok'])}")


def report_compare(legacy_path, current_path):
    print("\n=== reconciliation with Stream_K_Tier0.md §2 "
          "(1891 edges / 834 iso / 4 arc / 1053 free-form)")
    for tag, path in (("pre-fix recogniser (git 237be7f81a^)", legacy_path),
                      ("shipping recogniser", current_path)):
        # a census file may hold several models; the reconciliation is about the one with the
        # recognised faces in it, so take the richest rather than the first.
        res = max(load([path]), key=lambda r: len(edges_of(r)))
        per = {}
        for s in res["solids"]:
            for f in s["faces"]:
                per[(s["name"], f["index"])] = {"nRej": f["nRejected"], "edges": []}
            for e in s["edges"]:
                per[(s["name"], e["face"])]["edges"].append(e)
        faces = [f for s in res["solids"] for f in s["faces"]]
        bad = [v for v in per.values() if v["nRej"] > 0]
        tot = sum(len(v["edges"]) for v in bad)
        deg = sum(1 for v in bad for e in v["edges"] if e["degenerate"])
        rej = sum(v["nRej"] for v in bad)
        allrej = [e for e in edges_of(res) if not e["accepted"]]
        print(f"\n  {tag}")
        print(f"    recognised faces reaching the wire block   {len(faces)} "
              f"{dict(collections.Counter(f['recognizedKind'] for f in faces))}")
        print(f"    faces carrying at least one non-iso edge   {len(bad)}")
        print(f"    edges on those faces                        {tot}  "
              f"(degenerate {deg}, non-degenerate {tot - deg})")
        print(f"    of those: iso {tot - deg - rej}   non-iso {rej}")
        print(f"    buckets over ALL non-iso edges              "
              f"{dict(collections.Counter(e['bucket'] for e in allrej))}")
        print(f"    C by solid                                  "
              f"{dict(collections.Counter(e['solid'] for e in allrej if e['bucket'] == 'C'))}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json", nargs="*")
    ap.add_argument("--compare", nargs=2, metavar=("LEGACY", "CURRENT"))
    args = ap.parse_args()
    for res in load(args.json):
        report_model(res)
    if args.compare:
        report_compare(*args.compare)
    if not args.json and not args.compare:
        ap.error("give at least one census JSON, or --compare")
    return 0


if __name__ == "__main__":
    sys.exit(main())
