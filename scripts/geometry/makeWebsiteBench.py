#!/usr/bin/env python3
"""
Build the presentation website's per-part benchmark JSON (Track 2, `Plan_Presentation.md`) from
the JSON that the existing instruments already produce: `runOracleGate.py` (accuracy against the
OpenCascade oracle), `runXRayBench.py` (transport robustness against the same oracle) and
`o2-bench-detectorsbase-xray --perf` (per-call timing and memory for every representation).

This script does not run any of those tools itself -- it only reads their JSON output and reshapes
it into the schema `scripts/geometry/website_data/<part>.json` documents. That split is
deliberate: the three instruments need three different environments (the O2 build, the pythonOCC
converter env, and -- for timing -- a quiet single-threaded machine), and gluing them into one
subprocess pipeline would hide exactly the thing this project's ground rules insist on stating
explicitly -- which sampling produced which number. `scripts/geometry/website_data/README.md`
gives the exact command line that produces every JSON file this script consumes.

Usage
-----
    python3 makeWebsiteBench.py \\
        --gate-fixtures  <workdir>/gate_fixtures/gate.json \\
        --xray-fixtures  <workdir>/xray_fixtures/xray.json \\
        --perf-fixtures  <workdir>/perf_fixtures/perf.json \\
        --gate-bagger    <workdir>/gate_bagger/gate.json \\
        --xray-bagger    <workdir>/xray_bagger/xray.json \\
        --perf-bagger    <workdir>/perf_bagger/perf.json \\
        --gate-alice3    <workdir>/alice3_st002/oracle/ref_answers_run.json \\
        --xray-alice3    <workdir>/alice3_st002/xray/xray.json \\
        --perf-alice3    <workdir>/alice3_st002/perf/perf.json \\
        --load-avg-fixtures 1.2 --load-avg-bagger 1.2 --load-avg-alice3 1.2 \\
        --out-dir scripts/geometry/website_data

Every `--gate-*`/`--xray-*`/`--perf-*` file is the unmodified JSON one of the three instruments
already writes (a list of per-part records, keyed by `id`). This script only picks the hero
part(s) out of each file, by the `id` string the converter/fixture-generator assigned, and copies
the fields the website schema names -- it invents nothing and never overwrites a missing number
with a guess (see `_get`).
"""

import argparse
import json
from datetime import date
from pathlib import Path

GENERATED_DATE = "2026-08-22"
MACHINE = "aarch64 10-core shared (other tenants running o2-sim concurrently)"

# ------------------------------------------------------------------------------------------
# The hero parts (Plan_Presentation.md, Track 2), and which corpus/id they come from.
# ------------------------------------------------------------------------------------------
HERO_PARTS = [
    {
        "outfile": "BoomCylinderOuter",
        "part": "Bagger/BoomCylinderOuter_0_1_1_9",
        "model": "STEP_examples/Bagger.step",
        "corpus": "bagger",
        "story": "the Bagger ram: CSG tube-union, expect the shipped shape representation far "
                 "faster than the exact surface solid",
    },
    {
        "outfile": "Bucket",
        "part": "Bagger/Bucket_0_1_1_6",
        "model": "STEP_examples/Bagger.step",
        "corpus": "bagger",
        "story": "97 faces (69 plane, 22 cylinder, 4 sphere, 2 torus); the last part on this "
                 "corpus to ship tessellated, made exact by the ellipse-trim fix "
                 "(Stream_Q_EllipseTrim.md)",
    },
    {
        "outfile": "cyl_inter_cyl",
        "part": "cyl_inter_cyl/cyl_inter_cyl_0_1_1_1",
        "model": "synthetic fixture (make_boolean_fixtures.py)",
        "corpus": "fixtures",
        "story": "Steinmetz solid, two orthogonal r=10mm cylinders intersected; the minimized "
                 "reproducer of the Bagger BoomCylinderOuter failure and the watertightness case "
                 "-- the mesh's crossing list against the same OCCT oracle is the leak story",
    },
    {
        "outfile": "torus_union_cyl",
        "part": "torus_union_cyl/torus_union_cyl_0_1_1_1",
        "model": "synthetic fixture (make_boolean_fixtures.py)",
        "corpus": "fixtures",
        "story": "torus R=25mm r=8mm fused with a coaxial r=20mm cylinder; isolates the toroidal "
                 "(quartic) surface and its trimming from any transcendental trim curve",
    },
    {
        "outfile": "ST1829909_002",
        "part": "adhoc",
        "model": "ALICE_3_example/CAD_noETA.stp (leaf solid ST1829909_002)",
        "corpus": "alice3",
        "story": "ALICE3 scale: 965 analytic patches (535 plane, 308 cylinder, 102 cone, 20 "
                 "torus), the sub-patch BVH part (Stream_X_SubPatchBVH.md); exact surface only, "
                 "see README for why no CSG shape exists and how the mesh was kept safe to build",
        # This part was scored with the harness's own `--ref-answers` mode, not through
        # runOracleGate.py's driver (which always reconverts the whole 55-solid ALICE3 model --
        # unsafe at the converter's default mesh precision, see the README). That JSON therefore
        # has no `verdict` block (only the gate driver writes one); the shipped representation is
        # read here from the converter's own `csg_report.json` for this part instead (declined CSG:
        # "toroidal face (out of the recogniser's scope)"), recorded once and not re-derived.
        "shipsOverride": "surface",
    },
]

REP_ORDER = ["surface", "mesh", "shape"]
# Website schema key -> the perf JSON's own key (o2-bench-detectorsbase-xray --perf writes
# the full names "distFromOutside"/"distFromInside", NOT the "distout"/"distin" abbreviation
# runOracleGate.py's "oracle" columns use -- the two source JSONs use different conventions
# and accuracy_block() below hardcodes the gate's short names separately from this list).
FUNCS = [("contains", "contains"), ("distFromOutside", "distFromOutside"),
         ("distFromInside", "distFromInside"), ("safety", "safety")]


def _get(d, *path, default=None):
    """Chase a dotted path through nested dicts/lists; return `default` (never invent a value)."""
    cur = d
    for key in path:
        if cur is None:
            return default
        cur = cur.get(key) if isinstance(cur, dict) else None
    return default if cur is None else cur


def load_json(path):
    if path is None:
        return None
    return json.loads(Path(path).read_text())


def find_part(records, part_id):
    if records is None:
        return None
    for rec in records:
        if rec.get("id") == part_id:
            return rec
    return None


def accuracy_block(gate_rep):
    """One representation's `oracle` sub-block (against OpenCascade), reshaped to the website
    schema's `disagreements` (outside the model's own declared tolerance, gate convention) and
    `unexplained` (the strict subset with no missed-surface excuse) columns."""
    oracle = gate_rep.get("oracle") if gate_rep else None
    if oracle is None:
        return None
    disagreements = {}
    unexplained = {}
    for key in ("contains", "distout", "distin", "safety"):
        col = oracle.get(key)
        if col is None:
            continue
        disagreements[key] = col.get("nMismatchUnexplained", 0) + col.get("nMismatchMissedSurface", 0)
        unexplained[key] = col.get("nMismatchUnexplained", 0)
    cap_dev = oracle.get("capacityRelativeDeviation") if gate_rep.get("capacityComparable") else None
    return {
        "capacityRelativeDeviation": cap_dev,
        "capacityComparableNote": None if gate_rep.get("capacityComparable")
        else f"TGeoShape::Capacity() is {gate_rep.get('capacityMethod', 'not comparable')} for this "
             "representation; no capacity deviation is reported (project convention, see README)",
        "disagreements": disagreements,
        "unexplained": unexplained,
        "oracleTolerance": oracle.get("tolerance"),
        "oracleValid": oracle.get("valid"),
    }


def functions_block(perf_rep):
    if perf_rep is None:
        return None
    out = {}
    for website_key, perf_key in FUNCS:
        t = perf_rep.get(perf_key)
        if t is None:
            continue
        out[website_key] = {
            "nsPerCall": t.get("nsPerCallMedian"),
            "nsPerCallMin": t.get("nsPerCallMin"),
            "nsPerCallMax": t.get("nsPerCallMax"),
        }
    return out


def xray_block(xray_rep):
    if xray_rep is None:
        return None
    mode_a = xray_rep.get("modeA", {})
    vs_oracle = mode_a.get("vsOracle", {})
    return {
        "lost": vs_oracle.get("missingCrossings"),
        "extra": vs_oracle.get("extraCrossings"),
        "displaced": vs_oracle.get("displacedCrossings"),
        "unterminated": mode_a.get("unterminated"),
        "parity": mode_a.get("parityMismatchIntervals"),
        "parityNearBoundary": mode_a.get("parityMismatchNearBoundary"),
        "oddCrossingLists": mode_a.get("oddCrossingLists"),
        "raysTotal": mode_a.get("rays"),
        "raysIdentical": vs_oracle.get("raysIdentical"),
        "worstDeltaTCm": vs_oracle.get("worstDeltaT"),
    }


def build_part_json(cfg, gate_records, xray_records, perf_records, samples_meta, load_avg):
    gate_part = find_part(gate_records, cfg["part"])
    xray_part = find_part(xray_records, cfg["part"])
    perf_part = find_part(perf_records, cfg["part"])

    representations = []
    gate_reps = {r["name"]: r for r in (gate_part or {}).get("representations", [])}
    xray_reps = {r["name"]: r for r in (xray_part or {}).get("representations", [])}
    perf_reps = {r["name"]: r for r in (perf_part or {}).get("representations", [])}
    names = [n for n in REP_ORDER if n in gate_reps or n in perf_reps or n in xray_reps]

    for name in names:
        gate_rep = gate_reps.get(name)
        perf_rep = perf_reps.get(name)
        xray_rep = xray_reps.get(name)
        representations.append({
            "name": name,
            "class": (gate_rep or perf_rep or {}).get("class")
                     or (gate_rep or perf_rep or {}).get("shapeClass"),
            "primitives": _get(perf_rep, "primitives", default=_get(gate_rep, "primitives")),
            "primitiveKind": _get(perf_rep, "primitiveKind", default=_get(gate_rep, "primitiveKind")),
            "memoryBytes": _get(perf_rep, "structuralBytes"),
            "memoryFormula": _get(perf_rep, "structuralFormula"),
            "sidecarBytes": _get(perf_rep, "sidecarBytes"),
            "closeSeconds": _get(perf_rep, "closeSeconds"),
            "meshClosedBody": _get(gate_rep, "meshClosedBody", default=_get(perf_rep, "meshClosedBody")),
            "reliability": _get(gate_rep, "reliability"),
            "navigable": _get(gate_rep, "navigable"),
            "functions": functions_block(perf_rep) if perf_rep else None,
            "transport": {
                "nsPerRay": _get(perf_rep, "transport", "nsPerCallMedian"),
                "nsPerCrossing": _get(perf_rep, "transportNsPerCrossing"),
            } if perf_rep else None,
            "accuracy": accuracy_block(gate_rep) if gate_rep else None,
            "xray": xray_block(xray_rep) if xray_rep else None,
        })

    meta = {
        "part": cfg["outfile"],
        "partId": cfg["part"],
        "model": cfg["model"],
        "story": cfg["story"],
        "generated": GENERATED_DATE,
        "machine": MACHINE,
        "loadAverageDuringTiming": load_avg,
        "timingPreliminary": bool(load_avg is not None and load_avg > 2.0),
        "samples": samples_meta,
    }
    return {"meta": meta, "representations": representations}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    for corpus in ("fixtures", "bagger", "alice3"):
        ap.add_argument(f"--gate-{corpus}", required=True)
        ap.add_argument(f"--xray-{corpus}", required=True)
        ap.add_argument(f"--perf-{corpus}", required=True)
        ap.add_argument(f"--load-avg-{corpus}", type=float, required=True)
        ap.add_argument(f"--gate-points-{corpus}", type=int, default=2000)
        ap.add_argument(f"--gate-rays-{corpus}", type=int, default=2000)
        ap.add_argument(f"--gate-seed-{corpus}", type=int, default=1)
        ap.add_argument(f"--xray-raster-{corpus}", type=int, default=48)
        ap.add_argument(f"--xray-beams-{corpus}", type=int, default=96)
        ap.add_argument(f"--perf-points-{corpus}", type=int, default=4096)
        ap.add_argument(f"--perf-rays-{corpus}", type=int, default=4096)
        ap.add_argument(f"--perf-passes-{corpus}", type=int, default=9)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    a = vars(args)
    data = {}
    for corpus in ("fixtures", "bagger", "alice3"):
        data[corpus] = {
            "gate": load_json(a[f"gate_{corpus}"]),
            "xray": load_json(a[f"xray_{corpus}"]),
            "perf": load_json(a[f"perf_{corpus}"]),
            "load_avg": a[f"load_avg_{corpus}"],
            "samples": {
                "gatePoints": a[f"gate_points_{corpus}"],
                "gateRays": a[f"gate_rays_{corpus}"],
                "gateSeed": a[f"gate_seed_{corpus}"],
                "xrayRaster": a[f"xray_raster_{corpus}"],
                "xrayBeams": a[f"xray_beams_{corpus}"],
                "perfPoints": a[f"perf_points_{corpus}"],
                "perfRays": a[f"perf_rays_{corpus}"],
                "perfPasses": a[f"perf_passes_{corpus}"],
            },
        }

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for cfg in HERO_PARTS:
        d = data[cfg["corpus"]]
        part_json = build_part_json(cfg, d["gate"], d["xray"], d["perf"], d["samples"], d["load_avg"])
        out_path = out_dir / f"{cfg['outfile']}.json"
        out_path.write_text(json.dumps(part_json, indent=1) + "\n")
        written.append((cfg["outfile"], out_path))
        print(f"wrote {out_path}")

    # ---- summary.json: the cross-part headline table -------------------------------------
    summary = {"generated": GENERATED_DATE, "machine": MACHINE, "parts": []}
    for cfg in HERO_PARTS:
        d = data[cfg["corpus"]]
        gate_part = find_part(d["gate"], cfg["part"])
        perf_part = find_part(d["perf"], cfg["part"])
        xray_part = find_part(d["xray"], cfg["part"])
        shipped = _get(gate_part, "verdict", "shipped", "representation") or cfg.get("shipsOverride")
        perf_reps = {r["name"]: r for r in (perf_part or {}).get("representations", [])}
        xray_reps = {r["name"]: r for r in (xray_part or {}).get("representations", [])}

        def ns(rep_name, func_key):
            return _get(perf_reps.get(rep_name), func_key, "nsPerCallMedian")

        ratios = {}
        for func_key in ("contains", "distFromOutside", "distFromInside", "safety"):
            surf = ns("surface", func_key)
            fastest = None
            fastest_name = None
            for other in ("shape", "mesh"):
                v = ns(other, func_key)
                if v and (fastest is None or v < fastest):
                    fastest, fastest_name = v, other
            if surf and fastest:
                ratios[func_key] = {"vs": fastest_name, "ratio": surf / fastest}

        mesh_xray = xray_reps.get("mesh", {}).get("modeA", {})
        mesh_vs_oracle = mesh_xray.get("vsOracle", {})
        surf_xray = xray_reps.get("surface", {}).get("modeA", {})
        surf_vs_oracle = surf_xray.get("vsOracle", {})

        summary["parts"].append({
            "part": cfg["outfile"],
            "story": cfg["story"],
            "ships": shipped,
            "speedRatioSurfaceVsFastest": ratios,
            "meshLeak": {
                "lost": mesh_vs_oracle.get("missingCrossings"),
                "extra": mesh_vs_oracle.get("extraCrossings"),
                "unterminated": mesh_xray.get("unterminated"),
                "parity": mesh_xray.get("parityMismatchIntervals"),
                "raysTotal": mesh_xray.get("rays"),
            } if mesh_xray else None,
            "surfaceLeak": {
                "lost": surf_vs_oracle.get("missingCrossings"),
                "extra": surf_vs_oracle.get("extraCrossings"),
                "unterminated": surf_xray.get("unterminated"),
                "parity": surf_xray.get("parityMismatchIntervals"),
                "raysTotal": surf_xray.get("rays"),
            } if surf_xray else None,
        })
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {summary_path}")

    # ---- index.json: the file list scripts/geometry/website/js/data.js reads to discover which
    # per-part JSON files exist in this directory (it never lists a directory itself -- the
    # Artifact/static-hosting model this site is built for has no directory listing at all).
    index = {"generated": GENERATED_DATE, "files": [f"{name}.json" for name, _ in written]}
    index_path = out_dir / "index.json"
    index_path.write_text(json.dumps(index, indent=1) + "\n")
    print(f"wrote {index_path}")


if __name__ == "__main__":
    main()
