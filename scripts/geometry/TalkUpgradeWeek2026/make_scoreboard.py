"""Re-render the R6 corpus report into slide-sized tables and a traceable numbers.json.

Input: the corpus JSON written by roundTripReport.py (default ~/roundtrip_report/alice_roundtrip.json).
Output: tables/*.md, tables/*.csv and numbers.json under this directory.
"""
import argparse
import collections
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def pct(a, b):
    return 100.0 * a / b if b else 0.0


def write_md(path, header, rows):
    with open(path, "w") as f:
        f.write("| " + " | ".join(header) + " |\n")
        f.write("| " + " | ".join("---" for _ in header) + " |\n")
        for r in rows:
            f.write("| " + " | ".join(str(c) for c in r) + " |\n")


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def emit(name, header, rows, outdir):
    write_md(os.path.join(outdir, name + ".md"), header, rows)
    write_csv(os.path.join(outdir, name + ".csv"), header, rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", default=os.path.expanduser("~/roundtrip_report/alice_roundtrip.json"))
    ap.add_argument("--outdir", default=os.path.join(HERE, "tables"))
    args = ap.parse_args()

    d = json.load(open(args.corpus))
    parts, mods = d["parts"], d["modules"]
    os.makedirs(args.outdir, exist_ok=True)
    numbers = {"_source": args.corpus, "_note": "every value below is recomputed from that file"}

    # ---- T1: module summary -------------------------------------------------
    rows = []
    tot = collections.Counter()
    for m in sorted(mods, key=lambda m: -m["leafSolids"]):
        t = m["tiers"]
        flat = m.get("flat", 0)
        # 'csg' in tiers counts native CSG trees plus flat parts; separate them
        native = t.get("csg", 0) - flat
        scored, failed = m.get("knownSourceScored", 0), m.get("knownSourceFailed", 0)
        rows.append([m["module"], m["leafSolids"], native, flat, t.get("surface", 0),
                     t.get("mesh", 0), f"{pct(native + flat, m['leafSolids']):.1f}",
                     f"{scored - failed}/{scored}", m.get("writerDeclined", 0)])
        tot["leaf"] += m["leafSolids"]; tot["native"] += native; tot["flat"] += flat
        tot["surface"] += t.get("surface", 0); tot["mesh"] += t.get("mesh", 0)
        tot["scored"] += scored; tot["failed"] += failed
        tot["writerDeclined"] += m.get("writerDeclined", 0)
    rows.append(["**all 17**", tot["leaf"], tot["native"], tot["flat"], tot["surface"], tot["mesh"],
                 f"{pct(tot['native'] + tot['flat'], tot['leaf']):.2f}",
                 f"{tot['scored'] - tot['failed']}/{tot['scored']}", tot["writerDeclined"]])
    emit("T1_module_summary",
         ["module", "leaf solids", "native CSG", "flat CSG", "exact surfaces", "tessellated",
          "% CSG", "agrees with source", "writer declines"], rows, args.outdir)

    numbers["corpus"] = {
        "modules": len(mods), "leafSolids": tot["leaf"], "nativeCSG": tot["native"],
        "flatCSG": tot["flat"], "exactSurfaces": tot["surface"], "tessellated": tot["mesh"],
        "csgPercent": round(pct(tot["native"] + tot["flat"], tot["leaf"]), 2),
        "knownSourceScored": tot["scored"], "knownSourceFailed": tot["failed"],
        "writerDeclinedVolumes": tot["writerDeclined"],
    }

    # ---- T2: feature matrix, source shape class vs what shipped -------------
    by_class = collections.defaultdict(collections.Counter)
    for p in parts:
        by_class[p["sourceClass"]][p["ships"]] += 1
    rows = []
    for cls, c in sorted(by_class.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(c.values())
        csg = c["csg"] + c["flatcsg"]
        rows.append([cls, n, csg, c["surface"], c["mesh"], f"{pct(csg, n):.1f}"])
    emit("T2_feature_matrix",
         ["source shape class", "parts", "CSG", "exact surfaces", "tessellated", "% CSG"],
         rows, args.outdir)
    numbers["bySourceClass"] = {cls: dict(c) for cls, c in by_class.items()}

    # ---- composite-sourced recognition, the 87% headline -------------------
    comp = [p for p in parts if p["sourceClass"] == "TGeoCompositeShape"]
    comp_ships = collections.Counter(p["ships"] for p in comp)
    numbers["compositeSourced"] = {
        "parts": len(comp), "csg": comp_ships["csg"], "flatcsg": comp_ships["flatcsg"],
        "surface": comp_ships["surface"], "mesh": comp_ships["mesh"],
        "recognisedPercent": round(pct(comp_ships["csg"] + comp_ships["flatcsg"], len(comp)), 2),
        "distinctSourceVolumes": len(set((p["module"], p["sourceVolume"]) for p in comp))}

    # ---- T3a: converter declines, why a composite did not become CSG -------
    def bucket(reason):
        r = reason or "unstated"
        # order matters: the specific reasons first, the scope buckets last
        if "free-form" in r:
            return "free-form faces (outside plane/cylinder/cone/sphere/torus)"
        if "symmetric difference" in r:
            return "symmetric difference exceeds the acceptance band"
        if "_Loop twin" in r:
            return "the emitted shape disagrees with its own loop twin"
        if "is not a whole torus" in r or "not a wedge through its axis" in r:
            return "a partial torus is not a whole torus"
        if "sphere with additional faces" in r:
            return "a sphere with extra faces is out of scope"
        if "planar face is neither a cap nor a wedge" in r:
            return "a planar face is neither a cap nor a wedge of any axis cluster"
        if "axis clusters" in r:
            return "too many axis clusters for the cell matcher"
        if "planar faces" in r:
            return "too many planar faces for the cell matcher"
        return r.split(":")[0].strip()[:60]

    why = collections.Counter()
    why_raw = collections.Counter()
    for p in comp:
        if p["ships"] in ("csg", "flatcsg"):
            continue
        why[bucket(p["whyNotCSG"])] += 1
        why_raw[(p["whyNotCSG"] or "unstated").strip()[:80]] += 1
    rows = [[w, n, f"{pct(n, len(comp)):.2f}"] for w, n in why.most_common()]
    rows.append(["**total**", sum(why.values()), f"{pct(sum(why.values()), len(comp)):.2f}"])
    emit("T3a_converter_declines",
         ["why the composite did not become CSG", "parts", "% of composite-sourced"],
         rows, args.outdir)
    emit("T3a_converter_declines_raw", ["reason as recorded", "parts"],
         [[w, n] for w, n in why_raw.most_common()], args.outdir)
    numbers["converterDeclines"] = dict(why)

    # ---- T3b: writer declines, volumes that never reached the converter ----
    wrows = collections.Counter()
    wcls = collections.Counter()
    for m in mods:
        for row in m.get("writerDeclinedRows", []) or []:
            wcls[row.get("shapeClass", "?")] += 1
            wrows[(row.get("shapeClass", "?"), (row.get("reason") or "unstated")[:60])] += 1
    rows = [[c, r, n] for (c, r), n in wrows.most_common()]
    rows.append(["**total**", "", sum(wrows.values())])
    emit("T3b_writer_declines", ["shape class", "reason", "volumes"], rows, args.outdir)
    numbers["writerDeclines"] = {"byClass": dict(wcls),
                                 "total": sum(wcls.values())}

    # ---- exact tessellation census -----------------------------------------
    te = sum(1 for p in parts if p.get("tessellationExact"))
    numbers["tessellationExact"] = {"parts": te, "percent": round(pct(te, len(parts)), 1)}

    with open(os.path.join(HERE, "numbers.json"), "w") as f:
        json.dump(numbers, f, indent=2, sort_keys=True)
    print(json.dumps(numbers["corpus"], indent=2))
    print("compositeSourced:", json.dumps(numbers["compositeSourced"]))
    print("tessellationExact:", json.dumps(numbers["tessellationExact"]))
    print("writerDeclines:", json.dumps(numbers["writerDeclines"]))


main()
