#!/usr/bin/env python3
"""Turn `implicitTrimValidate.py`'s JSON into the tables of `Stream_R_CoSurfaceTrims.md`.

stdlib + numpy only, so it runs anywhere; no OCC, no build products.

    python3 probes/implicitTrimValidateReport.py /tmp/r/*.json
    python3 probes/implicitTrimValidateReport.py --per-solid /tmp/r/alice3.json
    python3 probes/implicitTrimValidateReport.py --faces ST0923290_002 /tmp/r/alice3.json
"""

import argparse
import json
import sys
from collections import Counter, defaultdict

import numpy as np


def load(paths):
    faces = []
    for p in paths:
        d = json.load(open(p))
        for m in d["models"]:
            for s in m["solids"]:
                for f in s["faces"]:
                    f["_model"] = m["label"]
                    f["_solid"] = s["name"]
                    f["_file"] = p
                    faces.append(f)
    return faces


def q(values, qs=(0.5, 0.9, 0.99, 1.0)):
    if not values:
        return [float("nan")] * len(qs)
    a = np.asarray(values, dtype=float)
    return [float(np.quantile(a, x)) for x in qs]


def headline(faces, title):
    ok = [f for f in faces if f.get("status") == "ok"]
    n = len(ok)
    if not n:
        print(f"\n## {title}: no measured face")
        return
    samples = sum(f["nSamples"] for f in ok)
    nin = sum(f["nIn"] for f in ok)
    nout = sum(f["nOut"] for f in ok)
    fp = sum(f["nFalsePositive"] for f in ok)
    fn = sum(f["nFalseNegative"] for f in ok)
    leak = sum(f["nLeakSamples"] for f in ok)
    single = sum(1 for f in ok if f["cellsIn"] == 1)
    singlemaj = sum(1 for f in ok if f["cellsInMajor"] <= 1)
    exact = sum(1 for f in ok if f["cellsIn"] == 1 and f["nLeakSamples"] == 0
                and f["nFalsePositive"] == 0 and f["nFalseNegative"] == 0)
    noleak = sum(1 for f in ok if f["nLeakSamples"] == 0)
    print(f"\n## {title}")
    print(f"  faces measured                     {n}")
    print(f"  surface samples                    {samples}   (IN {nin} / OUT {nout})")
    print(f"  faces with NO outside sample       "
          f"{sum(1 for f in ok if f['nOut'] == 0)}  (the near field cannot over-accept there)")
    print(f"  --- the naive conjunction ---")
    print(f"  faces exact (0 FP, 0 FN)           {exact} / {n}")
    print(f"  false positives (rule IN, truth OUT) {fp}")
    print(f"  false negatives (rule OUT, truth IN) {fn}")
    print(f"  --- the arrangement ---")
    print(f"  faces whose region is ONE cell     {single} / {n}")
    print(f"  ... one *major* cell (>=1% lobes)  {singlemaj} / {n}")
    print(f"  faces with zero leak               {noleak} / {n}")
    print(f"  leaking OUT samples                {leak} / {nout}")


def distributions(faces, title):
    ok = [f for f in faces if f.get("status") == "ok"]
    if not ok:
        return
    print(f"\n## {title} -- distributions (p50 / p90 / p99 / max)")
    rows = [
        ("trimming surfaces per face K", [f["nTrimSurfaces"] for f in ok]),
        ("boundary edges with a neighbour", [f["nBoundaryEdgesWithNeighbour"] for f in ok]),
        ("distinct cells in the region", [f["cellsIn"] for f in ok]),
        ("major cells (>=1% of interior)", [f["cellsInMajor"] for f in ok]),
        ("wires per face", [f["nWires"] for f in ok]),
    ]
    for name, vals in rows:
        a, b, c, d = q(vals)
        print(f"  {name:34s} {a:8.2f} {b:8.2f} {c:8.2f} {d:8.0f}")
    fp = [f["fpMaxDepthCm"] for f in ok if f["nFalsePositive"]]
    fn = [f["fnMaxDepthCm"] for f in ok if f["nFalseNegative"]]
    lk = [f["leakMaxBoundaryDistCm"] for f in ok if f["nLeakSamples"]]
    for name, vals in (("FP depth on the wrong side (cm)", fp),
                       ("FN depth on the wrong side (cm)", fn),
                       ("leak distance from boundary (cm)", lk)):
        if vals:
            a, b, c, d = q(vals)
            print(f"  {name:34s} {a:8.2e} {b:8.2e} {c:8.2e} {d:8.2e}   n={len(vals)}")
        else:
            print(f"  {name:34s} {'--':>8s}                              n=0")


def failure_modes(faces, title):
    ok = [f for f in faces if f.get("status") == "ok"]
    allf = faces
    n = len(ok)
    if not n:
        return
    print(f"\n## {title} -- the failure modes the design must answer to")
    trans = [f.get("minTransversalitySin") for f in ok
             if f.get("minTransversalitySin") is not None]
    rows = [
        ("declined outright (see below)", sum(1 for f in allf if f.get("status") != "ok")),
        ("a neighbour lies on the face's own surface",
         sum(1 for f in ok if f["nCoincidentTrimSurfaces"])),
        ("a neighbour has no evaluable surface",
         sum(1 for f in ok if f["nNonEvaluableNeighbours"])),
        ("a seam edge (the face is its own neighbour)",
         sum(1 for f in ok if f["nSeamEdges"])),
        ("a degenerate edge (apex / pole)", sum(1 for f in ok if f["nDegenerateEdges"])),
        ("an edge shared by more than 2 faces",
         sum(1 for f in ok if f["nMultiNeighbourEdges"])),
        ("one surface bounds the face on >1 edge",
         sum(1 for f in ok if f["nSurfacesBoundingTwice"])),
        ("one surface needs BOTH senses (cuts the interior)",
         sum(1 for f in ok if f["nAmbiguousSenseSurfaces"])),
        ("an inner wire (a hole)", sum(1 for f in ok if f["nInnerWires"])),
        ("a tangent neighbour (min |sin| < 1e-3)",
         sum(1 for t in trans if t < 1e-3)),
        ("a near-tangent neighbour (min |sin| < 1e-2)",
         sum(1 for t in trans if t < 1e-2)),
    ]
    for name, v in rows:
        print(f"  {name:52s} {v:6d} / {n}")
    st = Counter(f.get("status") for f in allf if f.get("status") != "ok")
    for k, v in st.most_common():
        print(f"    declined: {k:44s} {v}")


def per_solid(faces):
    by = defaultdict(list)
    for f in faces:
        by[(f["_model"], f["_solid"])].append(f)
    print("\n## per solid")
    print(f"  {'model':10s} {'solid':24s} {'faces':>6s} {'exact':>6s} {'1cell':>6s} "
          f"{'FP':>6s} {'FN':>7s} {'leak':>6s} {'maxK':>5s} {'maxCells':>8s} {'worstDepth':>11s}")
    for (mod, sol), fs in sorted(by.items()):
        ok = [f for f in fs if f.get("status") == "ok"]
        if not ok:
            continue
        exact = sum(1 for f in ok if f["cellsIn"] == 1 and f["nLeakSamples"] == 0
                    and not f["nFalsePositive"] and not f["nFalseNegative"])
        one = sum(1 for f in ok if f["cellsIn"] == 1)
        fp = sum(f["nFalsePositive"] for f in ok)
        fn = sum(f["nFalseNegative"] for f in ok)
        lk = sum(f["nLeakSamples"] for f in ok)
        mk = max(f["nTrimSurfaces"] for f in ok)
        mc = max(f["cellsIn"] for f in ok)
        wd = max(max(f["fpMaxDepthCm"], f["fnMaxDepthCm"]) for f in ok)
        print(f"  {mod[:10]:10s} {sol[:24]:24s} {len(ok):6d} {exact:6d} {one:6d} "
              f"{fp:6d} {fn:7d} {lk:6d} {mk:5d} {mc:8d} {wd:11.3e}")


def worst(faces, n=25, key="depth"):
    ok = [f for f in faces if f.get("status") == "ok"]
    if key == "leak":
        sel = [f for f in ok if f["nLeakSamples"]]
        sel.sort(key=lambda f: -f["leakMaxBoundaryDistCm"])
        print("\n## worst faces by LEAK distance from the boundary "
              "(no DNF over these surfaces can fix these)")
        for f in sel[:n]:
            print(f"  {f['_solid']:22s} f#{f['face']:<4d} {f['selfKind']:9s} K={f['nTrimSurfaces']:<2d} "
                  f"leak {f['nLeakSamples']:5d}/{f['nOut']:<5d}  maxdist {f['leakMaxBoundaryDistCm']:.3e} cm  "
                  f"med {f['leakMedBoundaryDistCm']:.3e}  cellsIn {f['cellsIn']}")
        return
    sel = [f for f in ok if f["nFalsePositive"] or f["nFalseNegative"]]
    sel.sort(key=lambda f: -max(f["fpMaxDepthCm"], f["fnMaxDepthCm"]))
    print("\n## worst faces for the naive conjunction, by depth on the wrong side")
    for f in sel[:n]:
        print(f"  {f['_solid']:22s} f#{f['face']:<4d} {f['selfKind']:9s} K={f['nTrimSurfaces']:<2d} "
              f"cells {f['cellsIn']:<3d} FP {f['nFalsePositive']:<5d} FN {f['nFalseNegative']:<5d} "
              f"depth {max(f['fpMaxDepthCm'], f['fnMaxDepthCm']):.3e} cm  "
              f"fnDist {f['fnMedBoundaryDistCm']:.2e}  {f['trimSurfaceKinds']}")


def pairs(faces):
    c = Counter()
    for f in faces:
        if f.get("status") != "ok":
            continue
        for k, v in f.get("trimSurfaceKinds", {}).items():
            c[f["selfKind"] + " trimmed by " + k] += v
    print("\n## trimming-surface pair frequency (distinct surfaces, not edges)")
    for k, v in c.most_common():
        print(f"  {k:40s} {v}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json", nargs="+")
    ap.add_argument("--per-solid", action="store_true")
    ap.add_argument("--worst", type=int, default=20)
    ap.add_argument("--split-population", action="store_true",
                    help="report the recognised (candidate) and stored (control) faces apart")
    args = ap.parse_args()

    faces = load(args.json)
    complete = [f for f in faces if f.get("implicitComplete")]
    partial = [f for f in faces if f.get("status") == "ok" and not f.get("implicitComplete")]
    headline(faces, "all corpora, every face with an analytic surface")
    headline(complete, "all corpora, faces where EVERY boundary edge has an evaluable, "
                       "non-coincident neighbour")
    headline(partial, "all corpora, faces with at least one unrepresentable boundary edge "
                      "(the idea does not apply)")
    distributions(complete, "implicit-complete faces")
    failure_modes(faces, "all corpora")
    if args.split_population:
        for pop in ("candidate", "control"):
            sub = [f for f in complete if f.get("population") == pop]
            headline(sub, f"implicit-complete, population = {pop}")
            distributions(sub, f"implicit-complete, population = {pop}")
    by_model = defaultdict(list)
    for f in complete:
        by_model[f["_model"]].append(f)
    for mod, fs in sorted(by_model.items()):
        headline(fs, f"model {mod}, implicit-complete faces")
    if args.per_solid:
        per_solid(complete)
    pairs(complete)
    worst(complete, args.worst)
    worst(complete, args.worst, key="leak")


if __name__ == "__main__":
    main()
