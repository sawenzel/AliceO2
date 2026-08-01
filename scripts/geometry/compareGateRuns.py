#!/usr/bin/env python3
"""Diff two gate.json reports column by column, with the scale law of each column applied.

Why not `diff`
--------------
Two rules of this project meet here. "Diff columns, not verdicts" -- a bit-identical *total* has
twice accompanied a real change, so the comparison has to be per part and per field. And "a
refuted experiment is not a refuted hypothesis" -- a comparison that cannot fail is not evidence,
so this tool is written to be run against a deliberately broken input first (`--self-test`).

What makes a gate.json comparison non-trivial is that under the position/scale sweep the two runs
are *not* meant to be identical. A count is invariant; a length scales linearly; a volume scales
as the cube; a relative deviation is invariant. Diffing raw would flag every distance column on a
x10 run and prove nothing. So every field carries a declared exponent in `_FIELD_EXPONENT` below,
the expectation is scaled by `factor ** exponent`, and what is reported is the residual after
that. Getting an exponent wrong is the obvious way to fool this tool, which is why each one is a
line of code someone can read rather than a heuristic.

What it can detect
------------------
  * any change in an integer column (surface/rim/edge counts, agreement counts, mismatch counts,
    navigability, reliability) -- exactly, since those are compared for equality;
  * any change in a real column beyond a stated relative band, after removing the scale law.

What it cannot detect
---------------------
  * a change confined to fields not present in gate.json;
  * a defect that moves a length by less than the band. The band is not a tolerance on the
    geometry, it is a floor on double arithmetic across two runs that ran different code paths; it
    is reported alongside every result so a reader can see how tight the statement is;
  * anything about points the harness never sampled.

Usage
-----
  compareGateRuns.py --baseline base/gate.json --candidate z400/gate.json --label "z+400 cm"
  compareGateRuns.py --baseline base/gate.json --candidate x10/gate.json  --scale 10
  compareGateRuns.py --baseline base/gate.json --self-test
"""

import argparse
import copy
import json
import sys
from pathlib import Path

# Timing is not a measurement here (NEXT.md: the gate's timing column moved 3x between runs of
# identical code) and checksums of point-derived quantities legitimately differ under a transform.
# `Seconds` is matched as a substring, not a suffix, because the two fields that matter most --
# closeShapeSecondsMesh and closeShapeSecondsSurface -- do not end in it.
_IGNORED_SUBSTRINGS = ("Seconds", "nsPerCall", "checksum")
_IGNORED_KEYS = {"id", "model", "worstOffenders", "rimDetail",
                 "timingCandidate", "timingReference", "timingCandidateLoop", "timingPruned",
                 "timingUnpruned"}

# The mesh is not the subject. It is built at an absolute chordal deflection, so a scaled shape
# tessellates differently on purpose; the columns that compare the exact solid against the mesh,
# and the triangle count itself, therefore move under a scaling for a reason that says nothing
# about the kernel. They are excluded by --gate-columns-only, which is how a scale comparison
# should be read; without the flag they are shown, because hiding them by default would be the
# kind of quiet exclusion this project keeps having to undo.
_MESH_DERIVED_PREFIXES = ("contains.", "distout.", "distin.", "safety.")
_MESH_DERIVED_KEYS = {"nTriangles"}

# The exponent of the length scale factor each real-valued field carries. Anything not listed is
# treated as dimensionless (exponent 0) -- counts, fractions, relative deviations, booleans.
_FIELD_EXPONENT = {
    # lengths, cm
    "maxRimIsolation": 1,
    "rimChordResolution": 1,
    "rimMatchTolerance": 1,
    "totalRimLength": 1,
    "unmatchedRimLength": 1,
    "maxSharedEdgeDeviation": 1,
    "worstDeviation": 1,
    "tolerance": 1,
    # volumes, cm^3
    "capacity": 3,
    "capacityCandidate": 3,
}

# Fields whose value is a physical measurement and must agree only to within double arithmetic
# across two independently converted shapes; everything else is required to be equal.
_REAL_BAND = 1.0e-9


def flatten(node, prefix=""):
    """Depth-first flatten of a part's report into {dotted path: scalar}."""
    flat = {}
    if isinstance(node, dict):
        for key, value in node.items():
            if key in _IGNORED_KEYS or any(s in key for s in _IGNORED_SUBSTRINGS):
                continue
            flat.update(flatten(value, f"{prefix}{key}."))
    elif isinstance(node, list):
        for i, value in enumerate(node):
            flat.update(flatten(value, f"{prefix}{i}."))
    else:
        flat[prefix.rstrip(".")] = node
    return flat


def exponent_of(path: str) -> int:
    return _FIELD_EXPONENT.get(path.rsplit(".", 1)[-1], 0)


def is_mesh_derived(path: str) -> bool:
    return (path.startswith(_MESH_DERIVED_PREFIXES) or
            path.rsplit(".", 1)[-1] in _MESH_DERIVED_KEYS)


def compare_part(base: dict, cand: dict, factor: float, gate_only: bool = False):
    """Return (list of differences, number of fields compared)."""
    flat_base = flatten(base)
    flat_cand = flatten(cand)
    if gate_only:
        flat_base = {k: v for k, v in flat_base.items() if not is_mesh_derived(k)}
        flat_cand = {k: v for k, v in flat_cand.items() if not is_mesh_derived(k)}
    differences = []
    for path in sorted(set(flat_base) | set(flat_cand)):
        if path not in flat_base:
            differences.append((path, "<absent>", flat_cand[path], "field only in candidate"))
            continue
        if path not in flat_cand:
            differences.append((path, flat_base[path], "<absent>", "field only in baseline"))
            continue
        b, c = flat_base[path], flat_cand[path]
        if isinstance(b, bool) or isinstance(c, bool) or isinstance(b, str) or isinstance(c, str):
            if b != c:
                differences.append((path, b, c, "differs"))
            continue
        if isinstance(b, int) and isinstance(c, int):
            if b != c:
                differences.append((path, b, c, f"integer differs by {c - b:+d}"))
            continue
        expected = b * (factor ** exponent_of(path))
        if expected == c:
            continue
        scale = max(abs(expected), abs(c))
        residual = abs(c - expected) / scale if scale else abs(c - expected)
        if residual > _REAL_BAND:
            differences.append((path, expected, c,
                                f"relative residual {residual:.3g} > {_REAL_BAND:g} "
                                f"(scale law: factor^{exponent_of(path)})"))
    return differences, len(set(flat_base) | set(flat_cand))


def key_reports(baseline, candidate):
    """Pair the two reports' parts up, and say how.

    The obvious key is the full part id, and it is the right one whenever the converter assigned
    the same logical-volume numbers on both sides -- which it does for a transformed fixture ladder.
    The fallback is the id's leading component, which for the fixture corpus is the fixture name
    and survives renumbering.

    It must *not* be the leading component unconditionally: every part of a single CAD model shares
    it (`Bagger/BasePin...`, `Bagger/Base...`), so keying on it silently collapses a 12-part report
    to one part and compares only whichever survived the dict. That is a comparison that cannot
    fail, and this file exists to not be one of those.
    """
    by_id = ({p["id"]: p for p in baseline}, {p["id"]: p for p in candidate})
    if set(by_id[0]) == set(by_id[1]):
        return by_id[0], by_id[1], "part id"
    by_stem = ({p["id"].split("/", 1)[0]: p for p in baseline},
               {p["id"].split("/", 1)[0]: p for p in candidate})
    collapsed = len(by_stem[0]) < len(baseline) or len(by_stem[1]) < len(candidate)
    if collapsed:
        return by_id[0], by_id[1], "part id (ids differ and the leading component is not unique)"
    return by_stem[0], by_stem[1], "leading id component"


def compare(baseline, candidate, factor, label, gate_only=False):
    base_by_key, cand_by_key, keying = key_reports(baseline, candidate)
    print(f"(paired by {keying})")
    print(f"=== {label} : {len(cand_by_key)} part(s) vs baseline's {len(base_by_key)}, "
          f"length factor {factor:g} ===")
    missing = sorted(set(base_by_key) - set(cand_by_key))
    extra = sorted(set(cand_by_key) - set(base_by_key))
    total_differences = 0
    for key in missing:
        print(f"  [MISSING]  {key}: in baseline, absent from candidate")
        total_differences += 1
    for key in extra:
        print(f"  [EXTRA]    {key}: in candidate, absent from baseline")
        total_differences += 1
    for key in sorted(set(base_by_key) & set(cand_by_key)):
        differences, n_fields = compare_part(base_by_key[key], cand_by_key[key], factor, gate_only)
        total_differences += len(differences)
        if not differences:
            print(f"  [same]     {key}: {n_fields} field(s) identical after the scale law")
            continue
        print(f"  [DIFFERS]  {key}: {len(differences)} of {n_fields} field(s)")
        for path, b, c, why in differences:
            print(f"               {path}: baseline {b!r} -> candidate {c!r}  ({why})")
    print(f"\n{total_differences} difference(s)")
    return total_differences


def _nudge(path):
    """A defect injector that multiplies a real field by 1 + 1e-8, on the first part where that
    field is not exactly zero.

    The "not exactly zero" is the point. A relative perturbation of 0.0 is 0.0, so injecting into a
    column that happens to be zero on the first fixture proves nothing and reports a miss -- which
    is how the first version of this self-test failed, and it is exactly the failure mode the
    project keeps hitting: an experiment structurally incapable of saying "yes". This picks a part
    where the field can move, and says which one.
    """
    def apply(report):
        for index, part in enumerate(report):
            node = part
            for key in path[:-1]:
                node = node[key]
            if node.get(path[-1]):
                node[path[-1]] = node[path[-1]] * (1. + 1.e-8)
                return index
        return None
    return apply


def self_test(baseline):
    """Prove the comparison is capable of saying "yes" before any run is allowed to say "no".

    Four injected defects, one per code path: an integer column, two real columns with different
    scale exponents, and a verdict string. A run that does not catch all four is a broken
    instrument and its green results mean nothing.
    """
    print("=== self-test: can this comparison detect a violation? ===")

    def bump_int(report):
        column = report[0]["oracle"]["contains"]
        column["nMismatchUnexplained"] = column["nMismatchUnexplained"] + 1
        return 0

    def downgrade(report):
        report[0]["navigation"]["reliability"] = "openBoundary"
        return 0

    cases = [
        ("integer column (one extra unexplained oracle disagreement)", bump_int),
        ("length column, exponent 1 (maxRimIsolation +1e-8 relative)",
         _nudge(["navigation", "maxRimIsolation"])),
        ("volume column, exponent 3 (capacityCandidate +1e-8 relative)",
         _nudge(["oracle", "capacityCandidate"])),
        ("verdict string (navigation reliability downgraded)", downgrade),
    ]
    caught = 0
    for what, break_it in cases:
        broken = copy.deepcopy(baseline)
        try:
            index = break_it(broken)
        except (KeyError, IndexError) as exc:
            print(f"  [SKIP] {what}: not present in this report ({exc})")
            continue
        if index is None:
            print(f"  [MISSED] {what}: no part carries a non-zero value, nothing was injected")
            continue
        differences, _ = compare_part(baseline[index], broken[index], 1.0)
        ok = bool(differences)
        caught += ok
        print(f"  [{'caught' if ok else 'MISSED'}] {what}  [part {baseline[index]['id']}]")
        for path, b, c, why in differences:
            print(f"             {path}: {b!r} -> {c!r} ({why})")
    print(f"\n{caught}/{len(cases)} injected defect(s) caught")
    return caught == len(cases)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True, type=Path)
    ap.add_argument("--candidate", type=Path)
    ap.add_argument("--scale", type=float, default=1.0,
                    help="uniform length factor applied to the candidate's geometry (1 for a pure "
                         "translation). Every real column is compared against baseline * "
                         "factor**exponent, with the exponent declared per field in this file.")
    ap.add_argument("--label", default=None)
    ap.add_argument("--gate-columns-only", action="store_true",
                    help="drop the columns that compare against the tessellated mesh, and the "
                         "triangle count. Under a scaling the mesh is deliberately not the same "
                         "mesh, so those columns move for a reason that is not the kernel's.")
    ap.add_argument("--self-test", action="store_true",
                    help="inject known defects into the baseline and report whether they are "
                         "caught; run this before believing any green comparison")
    args = ap.parse_args()

    baseline = json.loads(args.baseline.read_text())
    if args.self_test:
        return 0 if self_test(baseline) else 1
    if args.candidate is None:
        ap.error("--candidate is required unless --self-test is given")
    candidate = json.loads(args.candidate.read_text())
    label = args.label or f"{args.baseline} vs {args.candidate}"
    return 0 if compare(baseline, candidate, args.scale, label,
                        args.gate_columns_only) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
