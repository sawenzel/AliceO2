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
only if it agrees with ground truth outside the model's own declared tolerance -- and, where the
representation has the concept, only if it is a closed navigable manifold. Those are pass/fail,
so progress is countable and regressions are loud.

What the verdict is computed on
-------------------------------
The converter emits a **cascade**: CSG, else exact surfaces, else a tessellated mesh. The gate
scores every representation a part has side by side (`Stream_G_AnyShape.md`) but used to decide
pass or fail on the exact-surface one alone -- so a Bagger ram shipping as an exact
`TGeoCompositeShape` was reported as failing on a capacity deviation of a representation it no
longer uses.

The verdict is now computed on **the representation the part ships in**, and that representation
is read from **the converter's own cascade decision** (`csg_report.json`, carried per part into
`manifest.json` as a `shipped` block) -- never from whichever representation scores best. That
distinction is a correctness constraint, not a style one: choosing the best-scoring column would
let a part pass because some representation it does not use happens to be clean, and would make
this gate structurally incapable of ever reporting a bad conversion choice.

Three things follow, and all three are visible in the output:

  * the **historical surface-representation verdict is still computed and still printed**, so the
    series stays comparable and this is auditable as a change of interpretation rather than a
    relabelling;
  * the **other representations keep their full disagreement counts** -- the mesh's numbers are
    the input to the not-yet-written `auto`-mode fallback policy and are not to be hidden;
  * the **volume criterion is per representation**: `dV_sym` (the OCCT symmetric difference the
    emitter already computed) for a CSG part, the 1e-6 capacity band where capacity is a real
    measurement, and nothing at all where `TGeoShape::Capacity()` is Monte-Carlo sampled.

`--self-test` exercises that rule, with every positive case paired with the negative one that
must fail. See `Stream_I_Verdict.md`.

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
_OCC_PREFIX = {
    "PYTHONPATH": f"{_SW}/pythonOCC/latest/lib/python3.10/site-packages:"
                  f"{_SW}/Python-modules/latest/lib/python3.10/site-packages",
    "LD_LIBRARY_PATH": f"{_SW}/OCCT/latest/lib:{_SW}/Python/latest/lib",
}


def occ_env():
    """OCC-first, but *prepended* to the inherited environment rather than replacing it.

    This used to overwrite `PYTHONPATH`/`LD_LIBRARY_PATH` outright. That is the single reason the
    CSG hook had to defer writing `shape_<part>.root`: a converter launched from here had
    pythonOCC but `import ROOT` died on `libffi.so.6`, so the cascade could recognise a part as
    CSG and then not be able to emit it -- and `csg_report.json` then recorded the part one tier
    down, because the converter never dispatches to a file it did not write. Prepending gives the
    subprocess both (the alibuild Python 3.10 pythonOCC is built against *is* the O2 interpreter,
    `Stream_H_CSGEmitter.md`), so one conversion pass now produces the shapes and a cascade record
    that describes what `geom.C` really builds. Verified: `import OCC` and `import ROOT` both
    succeed in the resulting environment.
    """
    env = dict(os.environ)
    for key, prefix in _OCC_PREFIX.items():
        inherited = env.get(key, "")
        env[key] = f"{prefix}:{inherited}" if inherited else prefix
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


def rebase_manifest(manifest: dict, db_dir: Path, manifest_path: Path) -> dict:
    """Re-root a manifest's absolute paths on the directory it was actually found in.

    `makeTestPartDB.py` writes absolute paths. Copying a finished workdir and re-scoring it with
    `--skip-convert` therefore silently reads the *original* directory: no error, just columns
    that quietly describe the wrong files. That has cost at least two sessions a wrong conclusion.
    Detecting it is one string comparison, so there is no reason for it to stay a trap.
    """
    recorded = manifest.get("output_dir")
    actual = str(db_dir.resolve())
    if not recorded or recorded == actual:
        return manifest
    print(f"  [note] manifest.json was written for {recorded} but lives in {actual}; "
          "re-rooting its absolute paths onto this copy")
    old = recorded.rstrip("/")

    def rebase(value):
        if isinstance(value, str) and (value == old or value.startswith(old + "/")):
            return actual + value[len(old):]
        return value

    def walk(node):
        if isinstance(node, dict):
            return {k: walk(v) for k, v in node.items()}
        if isinstance(node, list):
            return [walk(v) for v in node]
        return rebase(node)

    manifest = walk(manifest)
    manifest["output_dir"] = actual
    manifest["rebased_from"] = recorded
    missing = [p["id"] for p in manifest.get("parts", []) if not Path(p["surfaces"]).exists()]
    if missing:
        raise RuntimeError(f"after re-rooting, {len(missing)} part(s) still have no sidecar "
                           f"(first: {missing[0]}); the DB copy is incomplete")
    manifest_path.write_text(json.dumps(manifest, indent=1))
    return manifest


def build_part_db(models, workdir: Path, skip_convert: bool, csg_mode: str = "auto") -> dict:
    db_dir = workdir / "db"
    manifest_path = db_dir / "manifest.json"
    if skip_convert:
        if not manifest_path.exists():
            raise RuntimeError(f"--skip-convert given but {manifest_path} does not exist")
        print(f"[1/4] reusing part DB {db_dir}")
        return rebase_manifest(json.loads(manifest_path.read_text()), db_dir, manifest_path)
    print(f"[1/4] converting {len(models)} model(s) into {db_dir} (--csg {csg_mode})")
    run([_OCC_PYTHON, _HERE / "makeTestPartDB.py", "--output", db_dir, "--force",
         "--csg", csg_mode, "--models", *models], env=occ_env())
    return rebase_manifest(json.loads(manifest_path.read_text()), db_dir, manifest_path)


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


def surface_verdict(part_report: dict):
    """The historical gate verdict: the exact-surface representation, and only that.

    Kept verbatim, and still reported for every part, so the series this project has been quoting
    since the branch point stays comparable and a change of *interpretation* cannot be mistaken
    for a change of measurement. The verdict that now decides the exit code is
    `representation_verdict` below, computed on the representation the part actually ships in.

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


# ------------------------------------------------------------------------------------------
# The representation-aware verdict
# ------------------------------------------------------------------------------------------
#
# The gate scored three representations of every part side by side (Stream G) but decided pass or
# fail on the exact-surface one alone. Once the converter's cascade became real -- CSG 7, exact
# surfaces 5, tessellated 1 on Bagger -- that made the verdict describe, for seven of thirteen
# parts, a representation they do not ship in.
#
# The rule below fixes that, and the one thing about it that is a correctness constraint rather
# than a matter of taste is *where the shipped representation is read from*: the converter's own
# cascade decision (`csg_report.json`, carried into `manifest.json` as each part's `shipped`
# block), and never the best-scoring column. Picking the best score would be marking one's own
# homework -- a part would pass because some representation it does not use happens to be clean,
# and the gate would become structurally incapable of ever reporting a bad conversion choice.
# It is also the reason the surface verdict above is kept and printed: this must be auditable as a
# change of interpretation, not a relabelling that quietly turns red green.

_VOLUME_BAND_RELATIVE = 1.e-6


def find_representation(part_report: dict, name: str):
    for rep in part_report.get("representations") or []:
        if rep.get("name") == name:
            return rep
    return None


def volume_criterion(rep: dict, shipped: dict):
    """The volume test that is meaningful *for this representation*, and its name.

    Three cases, in priority order:

      * `dV_sym` -- the OCCT symmetric-difference volume the CSG emitter already computed against
        the CAD solid, in cm^3, against the model's own tolerance band. This is the criterion for
        a CSG part and it is strictly stronger than any capacity comparison, because it measures
        where the two solids differ rather than how much of each there is.
      * `capacity` -- the 1e-6 relative band the gate has always applied, used wherever the
        representation's capacity is a real measurement (`exact-divergence` for the surface solid,
        `mesh-divergence` for O2Tessellated).
      * nothing -- `TGeoCompositeShape::Capacity()` is Monte-Carlo sampled in ROOT (~1e-2 relative
        error), so it is reported and never gated. Stream G measured a geometrically exact
        composite whose ROOT capacity was 3.3e-04 off: a capacity gate would have failed a correct
        shape by 330x. `capacityComparable` is already false there, so this branch is the second
        of two independent guards against that.
    """
    evidence = (shipped or {}).get("evidence") or {}
    dv = evidence.get("symmetricDifferenceCm3")
    band = evidence.get("bandCm3")
    if dv is not None and band is not None:
        ok = abs(dv) <= band
        return ("dV_sym", ok, abs(dv), band,
                f"dV_sym = {abs(dv):.3g} cm^3 against band {band:.3g} cm^3")
    oracle = rep.get("oracle") or {}
    if rep.get("capacityComparable", False):
        rel = abs(oracle.get("capacityRelativeDeviation", 0.))
        return ("capacity", rel <= _VOLUME_BAND_RELATIVE, rel, _VOLUME_BAND_RELATIVE,
                f"capacity off by {rel:.3g} relative")
    return ("none", True, None, None,
            f"no gateable volume measurement ({rep.get('capacityMethod', '?')} is not comparable; "
            "reported only)")


def representation_verdict(part_report: dict, name: str, shipped: dict):
    """Pass/fail for one named representation of one part.

    Same oracle answers, same four columns and the same "outside tolerance" definition as the
    historical verdict; what differs is which candidate was scored, which volume criterion applies
    and whether navigability is even a question. Closure, rims and NavigationReliability are
    `O2BVHSurfaceSolid` concepts -- a `TGeoCompositeShape` has no rims -- so navigability is
    required exactly where `closureApplicable` says it means something, and is not silently
    assumed anywhere else.
    """
    result = {"representation": name, "pass": False, "reasons": []}
    rep = find_representation(part_report, name)
    if rep is None:
        result["reasons"].append(
            f"the part ships as '{name}' but has no '{name}' representation in the scorecard")
        return result
    result["shapeClass"] = rep.get("shapeClass")
    result["source"] = rep.get("source")
    result["bboxDeviationFromOracle"] = rep.get("bboxDeviationFromOracle")
    reasons = []

    oracle = rep.get("oracle")
    if oracle is None:
        result["reasons"].append("no oracle answers")
        return result
    if not oracle.get("valid", False):
        reasons.append("reference BREP is not BRepCheck-valid")

    if rep.get("closureApplicable", False):
        if not rep.get("navigable", False):
            navigation = part_report.get("navigation", {})
            reasons.append(f"not navigable ({navigation.get('reliability', '?')}, "
                           f"{navigation.get('boundaryEdges', 0)} boundary edges)")
        result["navigable"] = rep.get("navigable")
        result["reliability"] = rep.get("reliability")
    else:
        result["navigable"] = None
        result["closureApplicable"] = False

    columns = {}
    for key in ("contains", "distout", "distin", "safety"):
        column = oracle.get(key)
        if column is None:
            continue
        bad = column.get("nMismatchUnexplained", 0) + column.get("nMismatchMissedSurface", 0)
        columns[key] = bad
        if bad:
            reasons.append(f"{key}: {bad} disagreement(s) outside tolerance "
                           f"(missed={column.get('nMismatchMissedSurface', 0)})")
    result["disagreements"] = columns

    criterion, ok, value, band, text = volume_criterion(rep, shipped)
    result["volumeCriterion"] = criterion
    result["volumeValue"] = value
    result["volumeBand"] = band
    result["volumeText"] = text
    if not ok:
        reasons.append(text)

    result["pass"] = not reasons
    result["reasons"] = reasons
    return result


def shipped_block(part_report: dict, manifest_index: dict):
    """Where this part's shipped representation is stated, taken as given.

    `manifest.json` carries the converter's cascade decision per part (`makeTestPartDB.py`, which
    copies it out of `csg_report.json`). A database built before that existed has no such block;
    the honest fallback is the cascade that database was built with -- exact surfaces, since every
    part in it has a sidecar -- and it says so in `decidedBy` rather than guessing silently.
    """
    entry = manifest_index.get(part_report.get("id"))
    shipped = (entry or {}).get("shipped")
    if shipped:
        return dict(shipped)
    return {"representation": "surface", "tier": "surface",
            "decidedBy": "default (this part DB predates the cascade record)",
            "source": None, "evidence": {}}


# ------------------------------------------------------------------------------------------
# The `shape_<part>.root` sidecar: hand-written fixtures for the any-TGeoShape path
# ------------------------------------------------------------------------------------------
#
# Two ladder fixtures are, by construction, *exactly* a ROOT shape: `box` is a TGeoBBox and
# `box_minus_cyl` is a TGeoCompositeShape (box - tube). Emitting them in the sidecar convention
# and gating them proves the whole any-shape path end to end -- convention, loader, frame,
# scoring -- before the CSG emitter exists, and gives that emitter's author a smoke test on day
# one. The numbers are read straight out of make_boolean_fixtures.py (mm there, cm here).
#
# Each entry is (part id, TGeoShape builder). The builder is only called when --fixture-shapes is
# given, so nothing here needs ROOT unless it is asked for.
def _build_box_shape():
    """`box`: a 20 x 30 x 40 mm box with its corner at the origin, i.e. 2 x 3 x 4 cm.

    TGeoBBox is centred on `fOrigin`, and the fixture is not centred, so the offset is the whole
    point of this fixture: it is the cheapest way for the frame convention to be wrong and be
    caught. The oracle's own bbox is [0,0,0] to [2,3,4] cm.
    """
    import ROOT
    from array import array
    return ROOT.TGeoBBox("shape", 1.0, 1.5, 2.0, array("d", [1.0, 1.5, 2.0]))


def _build_box_minus_cyl_shape():
    """`box_minus_cyl`: a 40 mm cube centred on the origin, minus an r = 8 mm axial through-hole.

    In cm: TGeoBBox(2,2,2) - TGeoTube(0, 0.8, 2.5), the tube deliberately longer than the cube so
    the hole is a through-hole rather than a blind one with two coincident cap faces.
    """
    import ROOT
    box = ROOT.TGeoBBox("cube", 2.0, 2.0, 2.0)
    drill = ROOT.TGeoTube("drill", 0.0, 0.8, 2.5)
    # The boolean node takes ownership of both operands; without this PyROOT frees them first.
    ROOT.SetOwnership(box, False)
    ROOT.SetOwnership(drill, False)
    node = ROOT.TGeoSubtraction(box, drill, ROOT.nullptr, ROOT.nullptr)
    ROOT.SetOwnership(node, False)
    return ROOT.TGeoCompositeShape("shape", node)


_FIXTURE_SHAPES = {
    "box/box_0_1_1_1": _build_box_shape,
    "box_minus_cyl/box_minus_cyl_0_1_1_1": _build_box_minus_cyl_shape,
}


def write_fixture_shapes(manifest: dict):
    """Write `shape_<VOL>_<LID>.root` next to the sidecar for every fixture that has a builder.

    The file format is the convention documented in DetectorsBase/O2SolidHarness.h: one object
    inheriting from TGeoShape, under the key "shape", in cm. These two fixtures are already stated
    in the part's own frame, so no `placement` key is written -- and its absence *is* the identity
    (`Stream_N_PlacedPrimitives.md`), so nothing has to be added here to keep meaning what it
    meant. This is the same `WriteTObject(shape, "shape")` that `harness::saveShapeToRootFile`
    performs; the C++ side is the authority and the unit test round-trips through it.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    written = []
    for part in manifest.get("parts", []):
        builder = _FIXTURE_SHAPES.get(part["id"])
        if builder is None:
            continue
        surfaces = Path(part["surfaces"])
        target = surfaces.parent / surfaces.name.replace("surfaces_", "shape_").replace(".bin", ".root")
        shape = builder()
        out = ROOT.TFile.Open(str(target), "RECREATE")
        out.WriteTObject(shape, "shape")
        out.Close()
        print(f"  wrote {target} ({shape.ClassName()}, capacity {shape.Capacity():.6g} cm^3)")
        written.append(part["id"])
    if not written:
        print("  [warn] --fixture-shapes given but no part in the DB has a builder "
              f"(known: {', '.join(sorted(_FIXTURE_SHAPES))})")
    return written


def column_disagreements(oracle: dict, key: str):
    column = oracle.get(key)
    if column is None:
        return None
    return column.get("nMismatchUnexplained", 0) + column.get("nMismatchMissedSurface", 0)


def print_representation_scorecard(report: list, shipped_by_id: dict = None):
    """One row per (part, representation): the tiered scorecard.

    This is a *report*, not a verdict, and it deliberately does not feed the exit code. The gate
    above scores the exact-surface representation, which is what the project has always gated;
    the columns here say what every other representation of the same part would have scored
    against the same oracle answers, which is the number the converter's future fallback policy
    has to be set from.

    Three things are reported for what they are rather than folded into a pass/fail:

      * `closure` is blank wherever it is meaningless. Closure, rims and NavigationReliability are
        O2BVHSurfaceSolid concepts; a TGeoCompositeShape has no rims, and a triangle mesh has a
        different notion (`meshClosedBody`) that is not the same claim. The harness omits the keys
        rather than reporting a default, so this prints "-".
      * `capacity` is "n/a" wherever TGeoShape::Capacity() is Monte-Carlo sampled -- which it is
        for every TGeoCompositeShape (ROOT throws 10000 accepted points into the bbox, i.e. ~1e-2
        relative error against a 1e-6 band). The intended acceptance for CSG parts is OCCT
        symmetric-difference volume, which is a later step.
      * `bboxDev` is the frame check: the max deviation, in cm, between the representation's own
        bounding box and the oracle's. A large number here means the two are not answering
        questions about the same object and every other column in the row is meaningless.
    """
    shipped_by_id = shipped_by_id or {}
    rows = [p for p in report if p.get("representations")]
    if not rows:
        return
    print("\n=== REPRESENTATION SCORECARD ===")
    print("  (`*` marks the representation the converter's cascade actually ships the part in; "
          "that is the one the gate verdict is computed on. The others are measured and reported, "
          "not gated -- the mesh columns in particular are the input to the auto-mode fallback "
          "policy and are deliberately shown in full.)")
    print(f"  {' ':<1}{'part':<44} {'repr':<8} {'class':<20} {'contains':>9} {'distout':>9} "
          f"{'distin':>9} {'safety':>9} {'capacity':>10} {'bboxDev':>9}  closure")
    for part_report in rows:
        ships = (shipped_by_id.get(part_report["id"]) or {}).get("representation")
        for rep in part_report["representations"]:
            oracle = rep.get("oracle", {})
            cells = []
            for key in ("contains", "distout", "distin", "safety"):
                bad = column_disagreements(oracle, key)
                cells.append("-" if bad is None else str(bad))
            if rep.get("capacityComparable", False):
                capacity = f"{abs(oracle.get('capacityRelativeDeviation', 0.)):.2e}"
            else:
                capacity = "n/a"
            bbox = rep.get("bboxDeviationFromOracle", -1.)
            bbox_text = "-" if bbox is None or bbox < 0. else f"{bbox:.2e}"
            if rep.get("closureApplicable", False):
                closure = f"{rep.get('reliability', '?')}" + ("" if rep.get("navigable") else " (NOT navigable)")
            elif "meshClosedBody" in rep:
                closure = f"meshClosedBody={rep['meshClosedBody']}"
            else:
                closure = "-  (not applicable to this representation)"
            mark = "*" if rep["name"] == ships else " "
            print(f"  {mark}{part_report['id']:<44} {rep['name']:<8} "
                  f"{rep.get('shapeClass', '?'):<20} "
                  f"{cells[0]:>9} {cells[1]:>9} {cells[2]:>9} {cells[3]:>9} {capacity:>10} "
                  f"{bbox_text:>9}  {closure}")

    # The totals, because the invariant this project defends is a count of disagreements and it
    # must never be quoted from a table the reader has to add up by hand.
    print("\n  totals per representation (disagreements outside tolerance, summed over parts):")
    names = []
    for part_report in rows:
        for rep in part_report["representations"]:
            if rep["name"] not in names:
                names.append(rep["name"])
    for name in names:
        totals = {}
        parts_with = 0
        clean = 0
        for part_report in rows:
            for rep in part_report["representations"]:
                if rep["name"] != name:
                    continue
                parts_with += 1
                bad_here = 0
                for key in ("contains", "distout", "distin", "safety"):
                    bad = column_disagreements(rep.get("oracle", {}), key)
                    if bad is not None:
                        totals[key] = totals.get(key, 0) + bad
                        bad_here += bad
                clean += (bad_here == 0)
        summary = "  ".join(f"{key}={totals.get(key, 0)}"
                            for key in ("contains", "distout", "distin", "safety"))
        print(f"    {name:<8} {summary}   ({clean}/{parts_with} part(s) with zero disagreements)")


# ------------------------------------------------------------------------------------------
# Self-test: every positive case paired with the negative one that must fail
# ------------------------------------------------------------------------------------------
#
# A scoring rule that cannot fail is not a gate. These are hand-built part reports in exactly the
# shape the harness writes, so the rule is exercised without a build, a model or an oracle, and
# every "this passes" is paired with the minimal mutation that must turn it red.

def _fake_column(bad=0):
    return {"nCompared": 100, "nMismatchUnexplained": bad, "nMismatchWithinBand": 0,
            "nMismatchMissedSurface": 0}


def _fake_oracle(bad=0, capacity_dev=0.0, valid=True):
    return {"valid": valid, "capacityRelativeDeviation": capacity_dev,
            "contains": _fake_column(bad), "distout": _fake_column(bad),
            "distin": _fake_column(bad), "safety": _fake_column(bad)}


def _fake_part(surface_capacity_dev=0.0, mesh_bad=0, shape_bad=0, navigable=True):
    """A part with all three representations: an exact surface solid, a wrong mesh, a CSG shape."""
    return {
        "id": "fake/Part_0_1_1_1",
        "navigation": {"navigable": navigable, "reliability": "reliable", "boundaryEdges": 0},
        "oracle": _fake_oracle(capacity_dev=surface_capacity_dev),
        "representations": [
            {"name": "surface", "shapeClass": "o2::base::O2BVHSurfaceSolid",
             "capacityMethod": "exact-divergence", "capacityComparable": True,
             "closureApplicable": True, "navigable": navigable, "reliability": "reliable",
             "bboxDeviationFromOracle": 1e-7,
             "oracle": _fake_oracle(capacity_dev=surface_capacity_dev)},
            {"name": "mesh", "shapeClass": "o2::base::O2Tessellated",
             "capacityMethod": "mesh-divergence", "capacityComparable": True,
             "closureApplicable": False, "meshClosedBody": True,
             "bboxDeviationFromOracle": 1e-4,
             "oracle": _fake_oracle(bad=mesh_bad, capacity_dev=3.0e-4)},
            {"name": "shape", "shapeClass": "TGeoCompositeShape",
             "capacityMethod": "root-montecarlo", "capacityComparable": False,
             "closureApplicable": False, "bboxDeviationFromOracle": 1e-7,
             "oracle": _fake_oracle(bad=shape_bad, capacity_dev=3.3e-4)},
        ],
    }


_CSG_CLEAN = {"representation": "shape", "tier": "csg", "decidedBy": "test",
              "evidence": {"symmetricDifferenceCm3": 0.0, "bandCm3": 1.0e-7}}
_CSG_DIRTY = {"representation": "shape", "tier": "csg", "decidedBy": "test",
              "evidence": {"symmetricDifferenceCm3": 1.0e-3, "bandCm3": 1.0e-7}}
_SURFACE = {"representation": "surface", "tier": "surface", "decidedBy": "test", "evidence": {}}
_MESH = {"representation": "mesh", "tier": "mesh", "decidedBy": "test", "evidence": {}}


def self_test(verbose=True):
    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))

    # 1. The case this whole change exists for: a part that is exact as CSG but whose *surface*
    #    representation misses the capacity band. Shipped verdict passes; surface verdict fails;
    #    and the pair must disagree, or the change would be invisible.
    part = _fake_part(surface_capacity_dev=1.39e-6)
    shipped = representation_verdict(part, "shape", _CSG_CLEAN)
    old_ok, old_reasons = surface_verdict(part)
    check("CSG part, clean dV_sym: shipped verdict passes", shipped["pass"], str(shipped["reasons"]))
    check("CSG part: surface verdict still fails on capacity", not old_ok, str(old_reasons))
    check("CSG part: the two verdicts disagree (the change is observable)", shipped["pass"] != old_ok)
    check("CSG part: gated on dV_sym, never on Capacity()",
          shipped["volumeCriterion"] == "dV_sym", shipped["volumeCriterion"])

    # 2. NEGATIVE CONTROL for the volume criterion: same part, same clean oracle columns, but a
    #    symmetric difference outside the band. Must fail.
    dirty = representation_verdict(part, "shape", _CSG_DIRTY)
    check("CSG part with dV_sym outside the band FAILS", not dirty["pass"], str(dirty["reasons"]))

    # 3. NEGATIVE CONTROL for the columns: a composite whose Capacity() is Monte-Carlo (3.3e-4 off,
    #    i.e. 330x the old band) but which is geometrically exact must still pass -- and the same
    #    shape with real oracle disagreements must still fail. Both halves matter: the first says
    #    the Monte-Carlo capacity is not gated, the second says something still is.
    exact_composite = representation_verdict(_fake_part(), "shape", _CSG_CLEAN)
    check("exact composite passes despite a 3.3e-4 Monte-Carlo capacity", exact_composite["pass"],
          str(exact_composite["reasons"]))
    broken = representation_verdict(_fake_part(shape_bad=17), "shape", _CSG_CLEAN)
    check("composite with 17 disagreements per column FAILS", not broken["pass"],
          str(broken["reasons"]))

    # 4. NEGATIVE CONTROL for the shipped tier: force a part to ship tessellated. The mesh is wrong
    #    on 0/12 Bagger parts, so a part shipped there must fail even though its surface and CSG
    #    representations are clean.
    mesh_part = _fake_part(mesh_bad=411)
    mesh_shipped = representation_verdict(mesh_part, "mesh", _MESH)
    mesh_surface_ok, _ = surface_verdict(mesh_part)
    check("part forced to ship tessellated FAILS on the mesh columns", not mesh_shipped["pass"],
          str(mesh_shipped["reasons"]))
    check("...while its surface representation is clean (so the failure is the mesh's)",
          mesh_surface_ok)
    check("mesh volume criterion is its own deterministic capacity, not dV_sym",
          mesh_shipped["volumeCriterion"] == "capacity", mesh_shipped["volumeCriterion"])

    # 5. Equivalence, where it must hold: for a part that ships as `surface` the new rule and the
    #    historical one must agree, part for part. That is what makes `--csg off` inert.
    for dev in (0.0, 1.39e-6):
        p = _fake_part(surface_capacity_dev=dev)
        new = representation_verdict(p, "surface", _SURFACE)
        old, _ = surface_verdict(p)
        check(f"surface-shipped part (capacity dev {dev:g}): new rule == historical rule",
              new["pass"] == old, f"new={new['pass']} old={old}")

    # 6. Navigability stays a precondition where it means something, and is not invented where it
    #    does not: an open surface solid fails, an open surface solid shipped as CSG does not
    #    inherit that failure (the CSG shape has no rims), and closure is never assumed true.
    open_part = _fake_part(navigable=False)
    check("open surface solid shipped as surface FAILS",
          not representation_verdict(open_part, "surface", _SURFACE)["pass"])
    check("the same part shipped as CSG is not failed for the surface solid's rims",
          representation_verdict(open_part, "shape", _CSG_CLEAN)["pass"])

    # 7. A part whose cascade says CSG but which has no `shape` column must fail loudly rather than
    #    fall back to something that happens to be clean.
    no_shape = _fake_part()
    no_shape["representations"] = [r for r in no_shape["representations"] if r["name"] != "shape"]
    missing = representation_verdict(no_shape, "shape", _CSG_CLEAN)
    check("cascade says CSG but no shape column: FAILS, does not fall back",
          not missing["pass"], str(missing["reasons"]))

    # 8. The provenance rule itself: the shipped representation comes from the manifest's cascade
    #    block, and a manifest that says `mesh` yields `mesh` even when `shape` scores perfectly.
    index = {"fake/Part_0_1_1_1": {"shipped": dict(_MESH)}}
    check("shipped representation is read from the manifest, not chosen",
          shipped_block(_fake_part(), index)["representation"] == "mesh")
    check("a DB with no cascade record falls back to `surface` and says so",
          shipped_block(_fake_part(), {})["decidedBy"].startswith("default"))

    failed = [c for c in checks if not c[1]]
    if verbose:
        for name, ok, detail in checks:
            print(f"  [{'ok  ' if ok else 'FAIL'}] {name}" + (f"   {detail}" if not ok else ""))
        print(f"\n{len(checks) - len(failed)}/{len(checks)} verdict self-checks passed")
    return not failed


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, default=None,
                        help="scratch directory for the part DB, samples and answers "
                             "(required unless --self-test)")
    parser.add_argument("--model", action="append", default=[],
                        help="CAD model to gate (repeatable)")
    parser.add_argument("--fixtures", action="store_true",
                        help="generate and gate the synthetic Boolean ladder")
    parser.add_argument("--skip-convert", action="store_true",
                        help="reuse an existing part DB in <workdir>/db")
    parser.add_argument("--fixture-shapes", action="store_true",
                        help="write the hand-built shape_<part>.root sidecars for the ladder "
                             "fixtures that are exactly a ROOT shape (box -> TGeoBBox, "
                             "box_minus_cyl -> TGeoCompositeShape) before scoring, so the "
                             "any-TGeoShape path is exercised end to end without an emitter. "
                             "Needs PyROOT; combine with --fixtures or --skip-convert.")
    parser.add_argument("--csg", default="auto", choices=["off", "auto", "required"],
                        help="converter CSG mode for the conversion step (default: %(default)s). "
                             "'auto' runs the production cascade CSG -> exact surfaces -> "
                             "tessellated, which is what makes the shipped-representation verdict "
                             "mean anything; 'off' reproduces the pre-cascade database, in which "
                             "every scored part ships as `surface` and the two verdicts coincide "
                             "by construction.")
    parser.add_argument("--self-test", action="store_true",
                        help="run the verdict rule's own positive/negative checks and exit; needs "
                             "no build, no model and no oracle")
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

    if args.self_test:
        return 0 if self_test() else 1
    if args.workdir is None:
        parser.error("--workdir is required (unless --self-test)")

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
    manifest = build_part_db(models, args.workdir, args.skip_convert, args.csg)
    db_dir = args.workdir / "db"
    sample_dir = args.workdir / "oracle"

    if args.fixture_shapes:
        print("[1b/4] writing hand-built TGeoShape sidecars for the exactly-representable fixtures")
        write_fixture_shapes(manifest)

    dump_samples(harness, db_dir, sample_dir, args.points, args.rays, args.seed, args.load_samples)
    run_oracle(manifest.get("parts", []), sample_dir, args.distance_limit)
    report = score(harness, db_dir, sample_dir, args.points, args.rays, args.seed,
                   args.workdir / "gate.json", args.load_samples)

    manifest_index = {p["id"]: p for p in manifest.get("parts", [])}
    unscored = manifest.get("unscoredParts", [])
    shipped_by_id = {}

    print("\n=== GATE SUMMARY ===")
    print("  verdict on the representation the part SHIPS in (the converter's cascade decision, "
          "read from csg_report.json);")
    print("  the historical surface-representation verdict is printed beside it so the two series "
          "stay comparable.")
    passed = 0
    surface_passed = 0
    changed = []
    for part_report in report:
        shipped = shipped_block(part_report, manifest_index)
        shipped_by_id[part_report["id"]] = shipped
        result = representation_verdict(part_report, shipped["representation"], shipped)
        old_ok, old_reasons = surface_verdict(part_report)
        passed += result["pass"]
        surface_passed += old_ok
        status = "PASS" if result["pass"] else "FAIL"
        old_status = "PASS" if old_ok else "FAIL"
        if result["pass"] != old_ok:
            changed.append((part_report["id"], shipped, old_status, status, old_reasons,
                            result["reasons"]))
        print(f"  [{status}] {part_report['id']}   ships: {shipped['representation']} "
              f"(tier {shipped.get('tier', '?')}, {result.get('shapeClass', '?')})"
              f"   [surface verdict: {old_status}]")
        print(f"           volume criterion: {result.get('volumeText', 'n/a')}")
        for reason in result["reasons"]:
            print(f"           {reason}")
        # Only where the two verdicts describe different objects is the surface reason extra
        # information; when the part ships as `surface` the two are the same statement.
        if not old_ok and shipped["representation"] != "surface":
            for reason in old_reasons:
                print(f"           (surface representation, reported not gated: {reason})")
        # The verdict, both verdicts and their provenance, written back into gate.json so that a
        # reader of the artifact sees exactly what a reader of this summary sees.
        part_report["verdict"] = {
            "shipped": shipped,
            "shippedVerdict": result,
            "surfaceVerdict": {"pass": old_ok, "reasons": old_reasons,
                               "note": "the historical gate verdict, kept for series continuity"},
        }
    total = len(report)
    print(f"\n{passed}/{total} scored part(s) pass on the representation they ship in")
    print(f"{surface_passed}/{total} scored part(s) pass on the surface representation "
          "(the historical number, unchanged in definition)")

    # Requirement 5: a leaf solid with no exact sidecar never enters the part DB and used to be
    # simply absent. It is not scoreable, and pretending otherwise in either direction would be
    # worse than saying so.
    if unscored:
        print(f"\n{len(unscored)} further leaf solid(s) ship in a representation this harness "
              "cannot score, and are therefore NOT counted above:")
        for entry in unscored:
            ship = entry.get("shipped", {})
            print(f"  [UNSCORED] {entry['id']}   ships: {ship.get('representation', '?')} "
                  f"(tier {ship.get('tier', '?')}, {entry.get('nFaces', '?')} faces)")
            print(f"             {entry.get('reason', '')}")
        n_leaf = total + len(unscored)
        print(f"  => {total} of {n_leaf} leaf solid(s) in the model(s) are scored by this gate.")

    if changed:
        print("\n  verdicts that changed because the representation changed, not because a "
              "measurement did:")
        for part_id, shipped, old_status, new_status, old_reasons, new_reasons in changed:
            print(f"    {part_id}: {old_status} (surface) -> {new_status} "
                  f"({shipped['representation']}, decided by {shipped.get('decidedBy', '?')})")
            for reason in old_reasons:
                print(f"        surface said: {reason}")
            for reason in new_reasons:
                print(f"        shipped says: {reason}")

    # The gate total above and the disagreement counts below are separate numbers and neither is
    # ever to be quoted without the other; printing them adjacently is the cheapest way to keep
    # that habit. The per-representation totals in the scorecard are the invariant this project
    # actually defends -- zero unexplained disagreements, per representation.
    unexplained = {key: 0 for key in ("contains", "distout", "distin", "safety")}
    for part_report in report:
        oracle = part_report.get("oracle") or {}
        for key in unexplained:
            bad = column_disagreements(oracle, key)
            if bad is not None:
                unexplained[key] += bad
    print("oracle disagreements outside tolerance (surface representation): "
          + "  ".join(f"{key}={value}" for key, value in unexplained.items()))

    print_representation_scorecard(report, shipped_by_id)

    gate_path = args.workdir / "gate.json"
    gate_path.write_text(json.dumps(report, indent=1))
    print(f"\nFull report: {gate_path}")
    if unscored:
        (args.workdir / "unscored.json").write_text(json.dumps(unscored, indent=1))
    # An unscoreable leaf solid is not a pass. The exit code stays the crude signal it has always
    # been -- read the counts, not the code -- but it must not go green on a model the gate only
    # partly covers.
    return 0 if passed == total and total > 0 and not unscored else 1


if __name__ == "__main__":
    sys.exit(main())
