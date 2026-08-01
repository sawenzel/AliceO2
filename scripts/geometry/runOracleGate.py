#!/usr/bin/env python3
"""Run the exact-surface acceptance gate: CAD -> converter -> solid -> OpenCascade oracle.

This is the "gate G1" driver from scripts/geometry/CodeReview_Fable.md, Section 9. It chains the
four pieces that already exist into one reproducible command and prints a per-part verdict:

  1. `makeTestPartDB.py`   converts CAD models, emitting per-part surface sidecars, meshes and
                           (with --dump-brep) the exact BREP each sidecar was extracted from
  2. `o2-bench-...-solid-harness --dump-samples`
                           writes the seeded sample sets, which nothing outside the harness can
                           regenerate
  3. `occtOracle.py`       answers those samples from the BREP, in OpenCascade
  4. `o2-bench-...-solid-harness --ref-answers`
                           scores the exact solid against those answers

Why a gate rather than a report
-------------------------------
The previous acceptance criterion compared against a tessellated mesh "within a chording band",
which cannot certify exactness and demonstrably hid whole classes of defect. A part passes here
only if it is a closed navigable manifold *and* it agrees with ground truth outside the model's
own declared tolerance. Those are pass/fail, so progress is countable and regressions are loud.

Usage
-----
  # generate the synthetic Boolean ladder, convert it, and gate it
  runOracleGate.py --fixtures --workdir /tmp/gate

  # gate an existing CAD model
  runOracleGate.py --model STEP_examples/Bagger.step --workdir /tmp/gate

  # re-score without reconverting (fast iteration on the C++ side)
  runOracleGate.py --workdir /tmp/gate --skip-convert
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_SW = Path(os.environ.get("ALIBUILD_ARCH_ROOT", Path.home() / "alisw/sw/ubuntu2404_aarch64"))
_BUILD = Path(os.environ.get("O2_BUILD_DIR", Path.home() / "alisw/sw/BUILD/O2-latest/O2"))

# pythonOCC is built against the alibuild Python 3.10; the system interpreter cannot import it.
_OCC_PYTHON = _SW / "Python/latest/bin/python3.10"
_OCC_ENV = {
    "PYTHONPATH": f"{_SW}/pythonOCC/latest/lib/python3.10/site-packages:"
                  f"{_SW}/Python-modules/latest/lib/python3.10/site-packages",
    "LD_LIBRARY_PATH": f"{_SW}/OCCT/latest/lib:{_SW}/Python/latest/lib",
}


def occ_env():
    env = dict(os.environ)
    env.update(_OCC_ENV)
    return env


def harness_env():
    """The freshly built harness must resolve the freshly built libraries.

    The binary installed in the alibuild prefix is stale after an incremental build, and so are
    its libraries; putting the build's stage directories first is what makes a run measure the
    code that was just compiled rather than the last installed one.
    """
    env = dict(os.environ)
    stage_libs = f"{_BUILD}/stage/lib:{_BUILD}/stage/lib64"
    env["LD_LIBRARY_PATH"] = stage_libs + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def run(cmd, **kwargs):
    printable = " ".join(str(c) for c in cmd)
    print(f"  $ {printable}", flush=True)
    result = subprocess.run([str(c) for c in cmd], **kwargs)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {printable}")
    return result


def sanitize_part_id(part_id: str) -> str:
    """Must match sanitizePartId() in Detectors/Base/test/runSolidHarness.cxx."""
    return "".join(c if (c.isalnum() or c in "-.") else "_" for c in part_id)


def find_harness() -> Path:
    candidate = _BUILD / "stage/bin/o2-bench-detectorsbase-solid-harness"
    if candidate.exists():
        return candidate
    found = shutil.which("o2-bench-detectorsbase-solid-harness")
    if found:
        print(f"  [warn] using {found} from PATH; it may be stale relative to {_BUILD}")
        return Path(found)
    raise RuntimeError(f"harness not found at {candidate} and not on PATH")


def build_part_db(models, workdir: Path, skip_convert: bool) -> dict:
    db_dir = workdir / "db"
    manifest_path = db_dir / "manifest.json"
    if skip_convert:
        if not manifest_path.exists():
            raise RuntimeError(f"--skip-convert given but {manifest_path} does not exist")
        print(f"[1/4] reusing part DB {db_dir}")
        return json.loads(manifest_path.read_text())
    print(f"[1/4] converting {len(models)} model(s) into {db_dir}")
    run([_OCC_PYTHON, _HERE / "makeTestPartDB.py", "--output", db_dir, "--force",
         "--models", *models], env=occ_env())
    return json.loads(manifest_path.read_text())


def dump_samples(harness: Path, db_dir: Path, sample_dir: Path, points: int, rays: int, seed: int,
                 load_samples: Path = None):
    print(f"[2/4] dumping sample sets into {sample_dir}")
    sample_dir.mkdir(parents=True, exist_ok=True)
    cmd = [harness, "--db", db_dir, "--dump-samples", sample_dir, "--points", points,
           "--rays", rays, "--seed", seed, "--only", "contains"]
    if load_samples is not None:
        # Round-trips the frozen set through the harness so the oracle and the scoring run below
        # read the same file the scoring run will use, rather than two copies that could drift.
        cmd += ["--load-samples", load_samples]
    run(cmd, env=harness_env(), stdout=subprocess.DEVNULL)


def run_oracle(parts, sample_dir: Path, distance_limit: int):
    print(f"[3/4] answering {len(parts)} part(s) with OpenCascade")
    answered = []
    for part in parts:
        brep = part.get("brep")
        if not brep or not Path(brep).exists():
            print(f"  [skip] {part['id']}: no .brep "
                  f"(re-run the converter with --dump-brep)")
            continue
        stem = sanitize_part_id(part["id"])
        samples = sample_dir / f"samples_{stem}.json"
        if not samples.exists():
            print(f"  [skip] {part['id']}: no sample file {samples.name}")
            continue
        answers = sample_dir / f"answers_{stem}.json"
        run([_OCC_PYTHON, _HERE / "occtOracle.py", "--brep", brep, "--samples", samples,
             "--out", answers, "--distance-limit", distance_limit, "--quiet"], env=occ_env())
        answered.append(part["id"])
    return answered


def score(harness: Path, db_dir: Path, sample_dir: Path, points: int, rays: int, seed: int,
          json_out: Path, load_samples: Path = None):
    print(f"[4/4] scoring against the oracle")
    cmd = [harness, "--db", db_dir, "--ref-answers", sample_dir, "--points", points,
           "--rays", rays, "--seed", seed, "--loop-crosscheck", "--edge-identity",
           "--json", json_out]
    if load_samples is not None:
        cmd += ["--load-samples", load_samples]
    run(cmd, env=harness_env())
    return json.loads(json_out.read_text())


def verdict(part_report: dict):
    """Pass/fail per the gate definition, with the reason when it fails.

    Navigability is a precondition, not a score: on an open or non-manifold surface set the
    parity answers are undefined, so agreeing with the oracle there would be luck rather than
    correctness.
    """
    reasons = []
    navigation = part_report.get("navigation", {})
    if not navigation.get("navigable", False):
        reasons.append(f"not navigable ({navigation.get('reliability', '?')}, "
                       f"{navigation.get('boundaryEdges', 0)} boundary edges)")
    oracle = part_report.get("oracle")
    if oracle is None:
        reasons.append("no oracle answers")
        return False, reasons
    if not oracle.get("valid", False):
        reasons.append("reference BREP is not BRepCheck-valid")
    for key in ("contains", "distout", "distin", "safety"):
        column = oracle.get(key)
        if column is None:
            continue
        bad = column.get("nMismatchUnexplained", 0) + column.get("nMismatchMissedSurface", 0)
        if bad:
            reasons.append(f"{key}: {bad} disagreement(s) outside tolerance "
                           f"(missed={column.get('nMismatchMissedSurface', 0)})")
    relative_capacity = abs(oracle.get("capacityRelativeDeviation", 0.))
    if relative_capacity > 1.e-6:
        reasons.append(f"capacity off by {relative_capacity:.3g} relative")
    return not reasons, reasons


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True,
                        help="scratch directory for the part DB, samples and answers")
    parser.add_argument("--model", action="append", default=[],
                        help="CAD model to gate (repeatable)")
    parser.add_argument("--fixtures", action="store_true",
                        help="generate and gate the synthetic Boolean ladder")
    parser.add_argument("--skip-convert", action="store_true",
                        help="reuse an existing part DB in <workdir>/db")
    parser.add_argument("--points", type=int, default=2000)
    parser.add_argument("--rays", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--distance-limit", type=int, default=1000,
                        help="points per category the oracle computes exact distances for")
    parser.add_argument("--transform", default=None,
                        help="transform applied to every --fixtures shape before conversion, for "
                             "the position/scale sweep: 'translate:dx,dy,dz' (mm) or 'scale:f'. "
                             "The STEP is the only shape in the pipeline, so the converter, the "
                             "sidecar, the mesh and the oracle's .brep all move with it.")
    parser.add_argument("--load-samples", type=Path, default=None,
                        help="reuse a frozen sample set from another run's <workdir>/oracle "
                             "instead of generating one. Required to compare a transformed run "
                             "with its baseline: the generator rejection-samples through the "
                             "tessellated reference, so a differently-meshed shape otherwise gets "
                             "different points and the columns are not comparable. Transform the "
                             "frozen set by the same map first (transformSamples.py).")
    args = parser.parse_args()

    args.workdir.mkdir(parents=True, exist_ok=True)
    models = list(args.model)
    if args.fixtures:
        fixture_dir = args.workdir / "fixtures"
        print(f"[0/4] generating the Boolean fixture ladder into {fixture_dir}")
        fixture_cmd = [_OCC_PYTHON, _HERE / "make_boolean_fixtures.py", "--outdir", fixture_dir]
        if args.transform:
            fixture_cmd += ["--transform", args.transform]
        run(fixture_cmd, env=occ_env())
        models += sorted(str(p) for p in fixture_dir.glob("*.step"))
    if not models and not args.skip_convert:
        parser.error("give --model and/or --fixtures, or --skip-convert to reuse a DB")

    harness = find_harness()
    manifest = build_part_db(models, args.workdir, args.skip_convert)
    db_dir = args.workdir / "db"
    sample_dir = args.workdir / "oracle"

    dump_samples(harness, db_dir, sample_dir, args.points, args.rays, args.seed, args.load_samples)
    run_oracle(manifest.get("parts", []), sample_dir, args.distance_limit)
    report = score(harness, db_dir, sample_dir, args.points, args.rays, args.seed,
                   args.workdir / "gate.json", args.load_samples)

    print("\n=== GATE SUMMARY ===")
    passed = 0
    for part_report in report:
        ok, reasons = verdict(part_report)
        passed += ok
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {part_report['id']}")
        for reason in reasons:
            print(f"           {reason}")
    total = len(report)
    print(f"\n{passed}/{total} part(s) pass the oracle gate")
    print(f"Full report: {args.workdir / 'gate.json'}")
    return 0 if passed == total and total > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
