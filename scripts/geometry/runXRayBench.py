#!/usr/bin/env python3
"""Drive the X-ray / geantino transport benchmark end to end, and print the tables.

  CAD -> converter -> raster of rays -> OpenCascade crossing lists -> stepping, two ways -> score

This is the transport-loop counterpart of `runOracleGate.py`. That gate measures SINGLE-SHOT
queries: from a sampled point, how far to the surface. This one measures a TRANSPORT: step, land
on the boundary, step again from there -- which is where geometry navigation actually fails, and
where none of the failure modes (zero-length steps, ping-ponging on a face, a particle that enters
and never exits, a crossing found twice, a step that overshoots) can be expressed as a
disagreement on `distout`.

It reuses `runOracleGate.py`'s environment handling and part-database builder rather than
reimplementing them, and it never writes into that gate's outputs.

Three tables come out, and all three are the point:

  1. CROSSING LISTS vs OpenCascade, per part, per representation, per stepping mode. Lists, not
     aggregates: a missing or extra crossing localises a defect to a face, which a disagreement
     rate cannot.
  2. ROBUSTNESS -- zero-length steps, non-advancing steps, unterminated transports, odd-length
     crossing lists, and mode (a) vs mode (b) disagreements. Nothing else in this project measures
     any of it.
  3. VOLUME BY CHORD INTEGRATION, with the raster's own achieved precision stated as a MEASURED
     number at the density used. Read the scoping note in the header of the volume table before
     quoting it: this is an instrument for gross errors and for composites (where ROOT's
     `TGeoCompositeShape::Capacity()` is Monte-Carlo and nothing else independent exists), not for
     the 1e-06 capacity residuals -- the divergence-theorem capacity is 1e-11 on exact parts and
     five orders sharper than anything a raster can do.

Usage
-----
  # the ladder fixtures, converted fresh
  runXRayBench.py --workdir /tmp/xray --fixtures

  # Bagger
  runXRayBench.py --workdir /tmp/xray_bag --model scripts/geometry/STEP_examples/Bagger.step

  # reuse a finished oracle-gate workdir's part DB (no reconversion)
  runXRayBench.py --workdir /tmp/xray_bag --reuse-db /tmp/gate_bag/db

  # the quartic witness: the same ladder at one tenth the size
  runXRayBench.py --workdir /tmp/xray_x01 --fixtures --transform scale:0.1
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from runOracleGate import (_BUILD, _OCC_PYTHON, harness_env, occ_env, rebase_manifest, run,
                           sanitize_part_id)


def find_benchmark() -> Path:
    candidate = _BUILD / "stage/bin/o2-bench-detectorsbase-xray"
    if candidate.exists():
        return candidate
    found = shutil.which("o2-bench-detectorsbase-xray")
    if found:
        print(f"  [warn] using {found} from PATH; it may be stale relative to {_BUILD}")
        return Path(found)
    raise RuntimeError(f"benchmark not found at {candidate} and not on PATH")


def build_part_db(models, workdir: Path, csg_mode: str, reuse_db: Path = None) -> dict:
    db_dir = workdir / "db"
    manifest_path = db_dir / "manifest.json"
    if reuse_db is not None:
        # A finished gate workdir is a perfectly good part database and reconverting it would only
        # burn minutes and risk a different one. `manifest.json` stores absolute paths, so it is
        # read in place rather than copied -- copying is the documented trap.
        manifest_path = reuse_db / "manifest.json"
        if not manifest_path.exists():
            raise RuntimeError(f"--reuse-db given but {manifest_path} does not exist")
        print(f"[1/4] reusing part DB {reuse_db}")
        return rebase_manifest(json.loads(manifest_path.read_text()), reuse_db, manifest_path)
    print(f"[1/4] converting {len(models)} model(s) into {db_dir} (--csg {csg_mode})")
    run([_OCC_PYTHON, _HERE / "makeTestPartDB.py", "--output", db_dir, "--force",
         "--csg", csg_mode, "--models", *models], env=occ_env())
    return rebase_manifest(json.loads(manifest_path.read_text()), db_dir, manifest_path)


def run_oracle(parts, ray_dir: Path):
    print(f"[3/4] answering {len(parts)} part(s) with OpenCascade "
          "(one call per ray returns every crossing)")
    answered = []
    for part in parts:
        brep = part.get("brep")
        if not brep or not Path(brep).exists():
            print(f"  [skip] {part['id']}: no .brep (re-run the converter with --dump-brep)")
            continue
        stem = sanitize_part_id(part["id"])
        rays = ray_dir / f"xrays_{stem}.json"
        if not rays.exists():
            print(f"  [skip] {part['id']}: no ray file {rays.name}")
            continue
        run([_OCC_PYTHON, _HERE / "xrayOracle.py", "--brep", brep, "--rays", rays,
             "--out", ray_dir / f"crossings_{stem}.json"], env=occ_env())
        answered.append(part["id"])
    return answered


def score(benchmark: Path, db_dir: Path, ray_dir: Path, json_out: Path, extra=()):
    print("[4/4] stepping both modes and scoring the crossing lists")
    run([benchmark, "--db", db_dir, "--ref-crossings", ray_dir, "--json", json_out, *extra],
        env=harness_env())
    return json.loads(json_out.read_text())


# ------------------------------------------------------------------------------------------
# The three tables
# ------------------------------------------------------------------------------------------

_MODES = (("modeA", "(a) shape"), ("modeB", "(b) nav"))


def short(part_id: str) -> str:
    name = part_id.split("/")[-1]
    return name[:34]


def print_crossing_table(report):
    print("\n=== CROSSING LISTS vs OPENCASCADE ===")
    print("  Lists, not aggregates, and LOST is kept apart from DISPLACED. `LOST` is a crossing OCCT "
          "found and\n  the candidate did not -- a wall a track walks through. `extra` is the "
          "reverse. `displaced` is a\n  crossing found in the right order but more than the "
          "match tolerance away: a wrong step length,\n  not a lost wall. `identical` counts rays "
          "whose whole ordered list matched.")
    header = (f"  {'part':<36} {'repr':<8} {'mode':<10} {'identical/rays':>18} {'LOST':>6} "
              f"{'extra':>6} {'displaced':>10} {'kind':>5} {'worst dt (cm)':>14}")
    print(header)
    totals = {}
    for part in report:
        for rep in part.get("representations", []):
            for key, label in _MODES:
                block = rep.get(key, {})
                comparison = block.get("vsOracle")
                if not comparison:
                    continue
                bucket = totals.setdefault((rep["name"], key), dict(
                    identical=0, rays=0, missing=0, extra=0, displaced=0, kind=0, worst=0.0,
                    clean=0, parts=0))
                bucket["identical"] += comparison["raysIdentical"]
                bucket["rays"] += comparison["rays"]
                bucket["missing"] += comparison["missingCrossings"]
                bucket["extra"] += comparison["extraCrossings"]
                bucket["displaced"] += comparison.get("displacedCrossings", 0)
                bucket["kind"] += comparison["kindMismatch"]
                bucket["worst"] = max(bucket["worst"], comparison["worstDeltaT"])
                bucket["parts"] += 1
                bucket["clean"] += (comparison["raysIdentical"] == comparison["rays"])
                print(f"  {short(part['id']):<36} {rep['name']:<8} {label:<10} "
                      f"{comparison['raysIdentical']:>8}/{comparison['rays']:<9} "
                      f"{comparison['missingCrossings']:>6} {comparison['extraCrossings']:>6} "
                      f"{comparison.get('displacedCrossings', 0):>10} "
                      f"{comparison['kindMismatch']:>5} {comparison['worstDeltaT']:>14.3e}")
    print("\n  totals (gate-style: the count and the denominator, never one without the other)")
    for (name, key), bucket in sorted(totals.items()):
        label = dict(_MODES)[key]
        print(f"    {name:<8} {label:<10} {bucket['identical']}/{bucket['rays']} rays identical, "
              f"LOST={bucket['missing']} extra={bucket['extra']} "
              f"displaced={bucket['displaced']} kind={bucket['kind']} "
              f"worst dt={bucket['worst']:.3e} cm   "
              f"({bucket['clean']}/{bucket['parts']} part(s) fully clean)")


_ROBUST_COLUMNS = (
    ("zeroLengthSteps", "zeroStep"),
    ("nonAdvancingSteps", "noAdv"),
    ("unstickPushes", "unstick"),
    ("iterationCapHits", "capHit"),
    ("unterminated", "unterm"),
    ("oddCrossingLists", "oddList"),
    ("nonAlternating", "nonAlt"),
    ("duplicateCrossings", "dupXing"),
    ("parityMismatchIntervals", "parity"),
    ("parityMismatchNearBoundary", "parityNB"),
    ("boundaryWithoutTransition", "noTrans"),
    ("originOutsideWorld", "outWorld"),
    ("originInside", "orgIn"),
)


def print_robustness_table(report):
    print("\n=== ROBUSTNESS (the part nothing else measures) ===")
    print("  zeroStep  a step at or below 1e-9 cm            unterm   the ray ended INSIDE the solid")
    print("  noAdv     the accumulated distance did not grow  oddList  odd-length crossing list")
    print("  unstick   a stalled step repaired with a nudge   nonAlt   two crossings of the same sense")
    print("  capHit    the iteration cap was reached          dupXing  two crossings within tolerance")
    print("  parity    Contains() at an interval midpoint contradicts the crossing list "
          "(the one check\n            independent of the stepping -- both modes alternate by "
          "construction). parityNB\n            is the same event excused because the midpoint is "
          "within the match tolerance of the\n            boundary, where neither side has a "
          "defined answer.")
    print("  noTrans   mode (b): a boundary was crossed but the volume did not change")
    print("  outWorld  mode (b): the ray origin was not in the navigator world -- a "
          "MISCONFIGURATION of\n            this benchmark, never a geometry defect. Any non-zero "
          "value invalidates the row.")
    head = (f"  {'part':<30} {'repr':<8} {'mode':<10} {'steps':>9} " +
            " ".join(f"{label:>8}" for _, label in _ROBUST_COLUMNS) + f" {'a-vs-b':>9}")
    print(head)
    totals = {}
    for part in report:
        for rep in part.get("representations", []):
            for key, label in _MODES:
                block = rep.get(key)
                if not block:
                    continue
                cells = [block.get(field, 0) for field, _ in _ROBUST_COLUMNS]
                bucket = totals.setdefault((rep["name"], key),
                                           dict(steps=0, cells=[0] * len(cells), avb=0, avbrays=0))
                bucket["steps"] += block.get("steps", 0)
                for i, value in enumerate(cells):
                    bucket["cells"][i] += value
                avb = ""
                if key == "modeB" and rep.get("modeAvsB"):
                    disagree = rep["modeAvsB"]["rays"] - rep["modeAvsB"]["raysIdentical"]
                    avb = str(disagree)
                    bucket["avb"] += disagree
                    bucket["avbrays"] += rep["modeAvsB"]["rays"]
                print(f"  {short(part['id']):<30} {rep['name']:<8} {label:<10} "
                      f"{block.get('steps', 0):>9} " +
                      " ".join(f"{value:>8}" for value in cells) + f" {avb:>9}")
    print("\n  totals")
    for (name, key), bucket in sorted(totals.items()):
        label = dict(_MODES)[key]
        summary = "  ".join(f"{lbl}={value}"
                            for (_, lbl), value in zip(_ROBUST_COLUMNS, bucket["cells"]))
        print(f"    {name:<8} {label:<10} steps={bucket['steps']}  {summary}")
        if key == "modeB":
            print(f"             mode (a) vs mode (b): {bucket['avb']} of {bucket['avbrays']} rays "
                  "disagree")


def print_volume_table(report):
    print("\n=== VOLUME BY CHORD INTEGRATION ===")
    print("  SCOPE, stated before the numbers. The raster's own achieved precision is the "
          "`raster` column\n  -- OCCT's chord integral over these same rays against OCCT's exact "
          "volume. It is a 1e-4 to 1e-5\n  instrument at the densities below, which is FOUR TO "
          "FIVE ORDERS coarser than the divergence-\n  theorem capacity already reported by the "
          "oracle gate (1e-11 on exact parts). It cannot resolve\n  the 1.3e-06 capacity "
          "residuals and must not be quoted as if it could. What it is for: gross\n  errors, and "
          "composites -- `TGeoCompositeShape::Capacity()` is Monte-Carlo in ROOT (~1e-2) and this\n"
          "  is the only independent volume those parts have.\n")
    print(f"  {'part':<30} {'repr':<8} {'N':>4} {'chord V (cm^3)':>16} {'vs OCCT chord':>14} "
          f"{'raster vs exact':>16} {'Capacity vs exact':>18}")
    for part in report:
        oracle = part.get("oracle")
        if not oracle:
            continue
        raster_n = part.get("raster", {}).get("n", 0)
        exact = oracle["capacity"]
        for rep in part.get("representations", []):
            block = rep.get("modeA", {})
            volume = block.get("volumeChordCm3")
            if volume is None:
                continue
            vs_chord = (volume - oracle["volumeChordCm3"]) / oracle["volumeChordCm3"] \
                if oracle["volumeChordCm3"] else 0.0
            capacity_dev = (rep.get("capacity", 0.0) - exact) / exact if exact else 0.0
            print(f"  {short(part['id']):<30} {rep['name']:<8} {raster_n:>4} {volume:>16.8g} "
                  f"{vs_chord:>14.3e} {oracle.get('chordVsExactRelative', 0.0):>16.3e} "
                  f"{capacity_dev:>18.3e}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--model", action="append", default=[])
    parser.add_argument("--fixtures", action="store_true")
    parser.add_argument("--reuse-db", type=Path, default=None,
                        help="use an existing part DB (e.g. a finished oracle-gate workdir's db) "
                             "instead of converting")
    parser.add_argument("--transform", default=None,
                        help="applied to every --fixtures shape before conversion "
                             "('scale:0.1', 'translate:dx,dy,dz' in mm); the STEP is the only "
                             "shape in the pipeline, so the oracle moves with it")
    parser.add_argument("--csg", default="auto", choices=["off", "auto", "required"])
    parser.add_argument("--raster", type=int, default=48,
                        help="N x N rays per beam axis (default %(default)s)")
    parser.add_argument("--axes", default="xyz")
    parser.add_argument("--beams", type=int, default=0,
                        help="fire N Fibonacci-spiral beam directions instead of the axis beams. "
                             "A parallel beam is DIRECTION-POOR -- three axes are three directions "
                             "however many rays are fired -- and the torus quartic defect is "
                             "invisible to them and visible to a fan. See Stream_J_XRay.md.")
    parser.add_argument("--tilt", type=float, default=0.0,
                        help="rotate every beam off its coordinate axis by this many degrees "
                             "(default %(default)s). An axis-aligned beam is a very special family "
                             "of ray/surface configurations; a tilted one is generic. Measured: "
                             "the known torus quartic defect is INVISIBLE at tilt 0 and visible at "
                             "tilt 12 -- see Stream_J_XRay.md.")
    parser.add_argument("--representations", default=None,
                        help="comma-separated subset of surface,mesh,shape")
    parser.add_argument("--margin", type=float, default=1.0e-3,
                        help="transverse padding of the raster window over the bounding box, cm")
    parser.add_argument("--parts", default=None, help="substring filter")
    parser.add_argument("--skip-oracle", action="store_true",
                        help="reuse the crossings_*.json already in <workdir>/xray")
    args = parser.parse_args()

    args.workdir.mkdir(parents=True, exist_ok=True)
    models = list(args.model)
    if args.fixtures:
        fixture_dir = args.workdir / "fixtures"
        print(f"[0/4] generating the Boolean fixture ladder into {fixture_dir}")
        cmd = [_OCC_PYTHON, _HERE / "make_boolean_fixtures.py", "--outdir", fixture_dir]
        if args.transform:
            cmd += ["--transform", args.transform]
        run(cmd, env=occ_env())
        models += sorted(str(p) for p in fixture_dir.glob("*.step"))
    if not models and args.reuse_db is None:
        parser.error("give --model and/or --fixtures, or --reuse-db")

    benchmark = find_benchmark()
    manifest = build_part_db(models, args.workdir, args.csg, args.reuse_db)
    db_dir = args.reuse_db if args.reuse_db is not None else args.workdir / "db"
    ray_dir = args.workdir / "xray"

    parts = manifest.get("parts", [])
    if args.parts:
        parts = [p for p in parts if args.parts in p["id"] or args.parts in p.get("model", "")]

    extra = ["--parts", args.parts] if args.parts else []
    if args.representations:
        extra += ["--representations", args.representations]
    if not args.skip_oracle:
        ray_dir.mkdir(parents=True, exist_ok=True)
        run([benchmark, "--db", db_dir, "--dump-rays", ray_dir, "--raster", args.raster,
             "--axes", args.axes, "--tilt", args.tilt, "--beams", args.beams,
             "--margin", args.margin, *extra],
            env=harness_env())
        run_oracle(parts, ray_dir)
    else:
        print(f"[2-3/4] reusing the crossing lists already in {ray_dir}")

    report = score(benchmark, db_dir, ray_dir, args.workdir / "xray.json", extra)

    print_crossing_table(report)
    print_robustness_table(report)
    print_volume_table(report)
    print(f"\nFull report: {args.workdir / 'xray.json'}")

    # Exit non-zero when any representation lost or invented a crossing against OpenCascade, or
    # when the two stepping modes disagree. As with the oracle gate the code is a crude signal --
    # read the counts.
    bad = 0
    for part in report:
        for rep in part.get("representations", []):
            for key, _ in _MODES:
                comparison = rep.get(key, {}).get("vsOracle")
                if comparison:
                    bad += comparison["missingCrossings"] + comparison["extraCrossings"] \
                        + comparison["kindMismatch"]
            if rep.get("modeAvsB"):
                bad += rep["modeAvsB"]["rays"] - rep["modeAvsB"]["raysIdentical"]
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
