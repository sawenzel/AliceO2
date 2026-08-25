#!/usr/bin/env python3
"""Acceptance test 3: score a converted part against the `TGeoShape` it was made from.

The other two acceptance tests are blind to the source. `csg/accept.py` measures the candidate
against the CAD B-rep, and `runOracleGate.py` measures the emitted shape against an independent
OCCT oracle; neither of them ever sees the `TGeoVolume` the CAD was written from. So a defect that
is already in the STEP file -- a writer bug, a unit slip, a mirrored placement -- passes both.

This closes that loop for the one corpus where the source *is* known: a geometry that
`O2_TGeoToCAD.py` wrote out and `O2_CADtoTGeo.py` read back. For every part the cascade carries as
CSG it finds the original volume through the writer report and asks three questions:

  * **class** -- and, where both sides are a `TGeoPcon`, the profile itself, sampled rather than
    compared array against array so that a redundant z plane is not reported as a difference in
    the shape;
  * **capacity** -- relative agreement, where both are analytic. `TGeoCompositeShape::Capacity()`
    is a Monte-Carlo sample, so a composite on either side is reported as *not comparable* rather
    than compared against a number that moves between runs;
  * **containment** -- a seeded point set in the original's bounding box, classified by both
    shapes. This is the sharp one: it is the only test here that a wrong *frame* cannot pass, and
    it composes the `shapePlacement` from `csg_report.json` exactly the way `geom.C` does.

Reading the verdict
-------------------
A **failure** is a disagreement about the shape: a wrong class, a profile off by more than the
recogniser's own declared resolution (`--profile-tolerance`, relative to the part's diagonal, 1e-6
by default because that is `csg/recognise.REL_TOL`), or a single containment disagreement. A
**flag** is a capacity that agrees to less than `--capacity-tolerance` (1e-9). The two are kept
apart on purpose: a part can be right to the pipeline's declared resolution and still not be exact
to 1e-9, and the flag is where that shows up. `--strict` makes flags fatal too.

Points nearer the boundary than `--skin` are not scored, because neither side claims to decide
them; the count of skipped points is reported so the number cannot quietly become the whole set.

One thing about the *converter's* bookkeeping it has to respect as well. An XCAF leaf label whose
shape is a compound of several disjoint solid bodies is split by `O2_CADtoTGeo.py` into one
logical volume per body (`#b1`, `#b2`, ..., which reach the report as parts named `..._b1`,
`..._b2`); TRD's eighteen BTOF barrels carry two bodies each. The *source* volume is still the
whole label, so for such a part only one direction of the containment test is a statement about
it: every point inside the emitted body must be inside the source, while a point inside the source
and outside this body belongs to a sibling body. Those parts are therefore scored one-way, their
capacity is reported as not comparable, and both facts are flagged so the weaker test is never
mistaken for the full one.

Two things about the writer's own bookkeeping that this has to respect: a geometry may hold
several distinct volumes under one name (the writer disambiguates them as `name#2`, `name#3`), and
a volume placed by a reflecting matrix is written as a Z-mirrored prototype named
`name__mirrored`. The first is resolved through the writer report's own record of the shape it
emitted, the second by reflecting z before comparing.

Usage
-----
  checkKnownSource.py --original o2sim_geometry.root --writer-report PIPE_writer_report.json \\
                      --converted /path/to/converter/output [--points 20000] [--json out.json]
  checkKnownSource.py --self-test

Exit status is non-zero if any part fails.
"""

import argparse
import json
import math
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from csg.primitives import placement_to_local  # noqa: E402

# Capacity is compared as a relative deviation; 1e-9 is far above the 1e-16 an identical analytic
# shape achieves. It is a flag rather than a failure -- see the module docstring.
CAPACITY_TOLERANCE = 1.0e-9
# The profile is compared against the resolution the recogniser itself declares (`REL_TOL`),
# relative to the part's own bounding-box diagonal.
PROFILE_TOLERANCE = 1.0e-6
# A point this close to either boundary is not scored. The two shapes are separate objects rebuilt
# through a STEP file, so their boundaries agree to a few ulp and not to zero.
DEFAULT_SKIN_CM = 1.0e-9
DEFAULT_POINTS = 20000
DEFAULT_SEED = 20260823

# `Capacity()` is a Monte-Carlo estimate for these, so it is not a number two objects can be
# compared on.
_SAMPLED_CAPACITY_CLASSES = ("TGeoCompositeShape", "TGeoUnion", "TGeoIntersection",
                             "TGeoSubtraction", "TGeoHalfSpace")


# ------------------------------------------------------------------------------------------
# profile comparison for the polycone family
# ------------------------------------------------------------------------------------------

def pcon_sections(shape, mirrored=False):
    """[(z, rmin, rmax)] of a polycone, Z-mirrored if the writer emitted the mirrored prototype."""
    rows = [(shape.GetZ(i), shape.GetRmin(i), shape.GetRmax(i)) for i in range(shape.GetNz())]
    if mirrored:
        rows = [(-z, rmin, rmax) for z, rmin, rmax in reversed(rows)]
    return rows


def _radii_at(sections, z):
    """(rmin, rmax) at `z`, clamped to the profile's own ends.

    Clamped rather than refused, because the two profiles have been through a STEP file and their
    axial extents agree to about 1e-12 cm, not exactly; the extent difference is measured
    separately in `pcon_profile_deviation`, so nothing is hidden by the clamp.
    """
    z = min(max(z, sections[0][0]), sections[-1][0])
    for i in range(len(sections) - 1):
        z0, rmin0, rmax0 = sections[i]
        z1, rmin1, rmax1 = sections[i + 1]
        if z1 <= z0:
            continue
        if z0 <= z <= z1:
            f = (z - z0) / (z1 - z0)
            return rmin0 + f * (rmin1 - rmin0), rmax0 + f * (rmax1 - rmax0)
    return sections[-1][1], sections[-1][2]


def pcon_profile_deviation(sa, sb, merge_tolerance=0.0):
    """The largest radial or axial disagreement, in cm, between two polycone profiles.

    Sampled inside every section of either shape rather than compared array against array, so a
    profile that carries a redundant plane -- legal, and what the recogniser emits when the CAD no
    longer shows one was there -- is not reported as a difference when the *shape* is the same.
    """
    levels = []
    for z in sorted({z for z, _r0, _r1 in sa} | {z for z, _r0, _r1 in sb}):
        if not levels or z - levels[-1] > merge_tolerance:
            levels.append(z)
    worst = max(abs(sa[0][0] - sb[0][0]), abs(sa[-1][0] - sb[-1][0]))
    for i in range(len(levels) - 1):
        z0, z1 = levels[i], levels[i + 1]
        if z1 <= z0:
            continue
        for f in (1.0e-9, 0.25, 0.5, 0.75, 1.0 - 1.0e-9):
            z = z0 + f * (z1 - z0)
            ra, rb = _radii_at(sa, z), _radii_at(sb, z)
            worst = max(worst, abs(ra[0] - rb[0]), abs(ra[1] - rb[1]))
    return worst


def _phi_deviation(a, b):
    """Degrees. A mirror in z leaves phi alone, so this needs no mirrored variant."""
    return max(abs(a.GetPhi1() - b.GetPhi1()), abs(a.GetDphi() - b.GetDphi()))


def shape_scale(shape):
    """The shape's bounding-box diagonal in cm, the length every relative tolerance is against."""
    return math.sqrt(shape.GetDX() ** 2 + shape.GetDY() ** 2 + shape.GetDZ() ** 2)


# ------------------------------------------------------------------------------------------
# the per-part comparison
# ------------------------------------------------------------------------------------------

def placement_is_identity(placement):
    if placement is None:
        return True
    for r in range(3):
        for c in range(4):
            want = 1.0 if r == c else 0.0
            if abs(placement[r][c] - want) > 1.0e-12:
                return False
    return True


def _bbox_of(shape):
    origin = [shape.GetOrigin()[i] for i in range(3)]
    half = [shape.GetDX(), shape.GetDY(), shape.GetDZ()]
    return origin, half


_ONE_BODY_OF_MANY = re.compile(r"_b\d+$")


def part_is_one_body_of_many(part):
    """True when the part is one body of a CAD label whose shape carried several.

    `O2_CADtoTGeo.py` gives every solid body of a multi-body XCAF leaf its own logical volume,
    keyed `#b1`, `#b2`, ...; the source volume is still the whole label. Measured on TRD: each of
    the eighteen BTOF barrels is two bodies, and `_b2` is exactly half the source's capacity.
    """
    return bool(_ONE_BODY_OF_MANY.search(part.get("part") or ""))


def contains_crosscheck(source, emitted, placement, n_points, seed, skin, max_report,
                        mirrored=False, one_way=False):
    """Classify a seeded point set against both shapes; every disagreement is reported.

    `mirrored` says the emitted shape is the source's Z-mirrored prototype, so the point is
    reflected before it enters the emitted shape's frame -- the same composition `geom.C` performs
    through the reflecting placement.

    `one_way` says the emitted shape is one body of a multi-body source, so only "inside the
    emitted shape implies inside the source" is a statement about it. The count of points the
    emitted shape *does* enclose is reported either way, so a one-way comparison cannot pass by
    enclosing nothing.
    """
    from array import array
    origin, half = _bbox_of(source)
    rng = random.Random(seed)
    scored = 0
    skipped = 0
    n_mismatches = 0
    n_inside_emitted = 0
    examples = []
    local = array("d", [0.0, 0.0, 0.0])
    probe = array("d", [0.0, 0.0, 0.0])
    for _ in range(n_points):
        point = tuple(origin[i] + rng.uniform(-half[i], half[i]) for i in range(3))
        probe[0], probe[1], probe[2] = point
        inside_source = bool(source.Contains(probe))
        if source.Safety(probe, inside_source) < skin:
            skipped += 1
            continue
        reflected = (point[0], point[1], -point[2]) if mirrored else point
        moved = placement_to_local(placement, reflected)
        local[0], local[1], local[2] = moved
        inside_emitted = bool(emitted.Contains(local))
        if emitted.Safety(local, inside_emitted) < skin:
            skipped += 1
            continue
        scored += 1
        if inside_emitted:
            n_inside_emitted += 1
        if inside_source != inside_emitted:
            if one_way and inside_source and not inside_emitted:
                continue                       # a sibling body of the same label carries it
            # Counted in full; only the first `max_report` are kept for printing, so a part with
            # thousands of disagreements does not report five.
            n_mismatches += 1
            if len(examples) < max_report:
                examples.append({"point": [float(c) for c in point],
                                 "local": [float(c) for c in moved],
                                 "source": inside_source, "emitted": inside_emitted})
    return {"points": scored, "skipped": skipped, "mismatches": n_mismatches,
            "insideEmitted": n_inside_emitted, "oneWay": bool(one_way), "examples": examples}


def reclose_flat_csg(shape):
    """Rebuild an `o2::base::O2FlatCSG`'s sub-cell boxes after it comes off a file.

    Its boxes and its BVH are transient by design (`Design_FlatCSGSolid.md` section 7 -- a stored
    BVH is a second thing that can disagree with the data it describes), so a streamed shape
    arrives un-closed and answers every query through its `_Loop` twins. Those are correct, and
    that is deliberately not good enough here: `geom.C` closes the shape it loads, so this test
    has to score the accelerated path the simulation actually runs. True for every other shape
    class, which needs nothing.
    """
    if shape.ClassName() != "o2::base::O2FlatCSG":
        return True
    if not shape.IsClosed():
        shape.CloseShape()
    return bool(shape.IsClosed())


def check_part(part, row, source_shape, emitted_shape, placement, n_points, seed, skin,
               capacity_tolerance, profile_tolerance, max_report):
    """Compare one converted part against its source shape. Returns a record."""
    mirrored = bool(row.get("mirrored"))
    one_body = part_is_one_body_of_many(part)
    source_class = row.get("shapeClass") or source_shape.ClassName()
    scale = max(shape_scale(source_shape), 1.0)
    record = {"part": part.get("part"), "volume": part.get("volume"),
              "source": row.get("name"), "mirrored": mirrored, "oneBodyOfMany": one_body,
              "sourceClass": source_class, "emittedClass": emitted_shape.ClassName(),
              "placementIsIdentity": placement_is_identity(placement),
              "classComparable": False, "classMatches": None,
              "profileDeviationCm": None, "profileToleranceCm": profile_tolerance * scale,
              "phiDeviationDeg": None,
              "capacityComparable": False, "capacitySource": None, "capacityEmitted": None,
              "capacityRelativeDeviation": None,
              "contains": None, "failures": [], "flags": []}

    # The class is compared where the two sides are *directly* comparable: same class, so the
    # parameters mean the same thing. A different class is not by itself a defect -- a TGeoCone
    # whose inner radius happens to be constant arrives as one cylinder and one cone and comes
    # back as a two-section TGeoPcon, which is the same solid stated more generally -- so it is
    # flagged and the metric checks below carry the verdict.
    same_class = source_class == emitted_shape.ClassName()
    record["classComparable"] = not one_body
    record["classMatches"] = None if one_body else same_class
    if one_body:
        record["flags"].append(
            "one body of a multi-body CAD label: the source is the whole label, so the class "
            "and the capacity are not comparable and containment is scored one-way")
    elif not same_class:
        record["flags"].append(
            f"class {emitted_shape.ClassName()} is not the source's {source_class}")

    # Under a non-identity placement the emitted profile is stated in the shape's own frame and
    # its z array is legitimately a different set of numbers; capacity and containment carry the
    # verdict there, and containment is the sharper of the two anyway.
    if (same_class and not one_body and source_class == "TGeoPcon"
            and record["placementIsIdentity"]):
        record["phiDeviationDeg"] = _phi_deviation(source_shape, emitted_shape)
        deviation = pcon_profile_deviation(pcon_sections(source_shape, mirrored),
                                           pcon_sections(emitted_shape),
                                           merge_tolerance=profile_tolerance * scale)
        record["profileDeviationCm"] = deviation
        if deviation > record["profileToleranceCm"]:
            record["failures"].append(
                f"the polycone profile is {deviation:.6g} cm off the source's, over the "
                f"{record['profileToleranceCm']:.3g} cm the recogniser claims")
        if record["phiDeviationDeg"] > 1.0e-9:
            record["failures"].append(
                f"phi differs from the source by {record['phiDeviationDeg']:.6g} deg")

    # The writer's own record of the volume it emitted is preferred over a second reading of the
    # shape: it is unambiguous even where two volumes share a name.
    capacity_source = row.get("capacity_cm3")
    if capacity_source is None and source_class not in _SAMPLED_CAPACITY_CLASSES:
        capacity_source = float(source_shape.Capacity())
    record["capacitySource"] = capacity_source
    record["capacityEmitted"] = float(emitted_shape.Capacity())
    comparable = (capacity_source is not None and capacity_source > 0.0 and not one_body
                  and source_class not in _SAMPLED_CAPACITY_CLASSES
                  and emitted_shape.ClassName() not in _SAMPLED_CAPACITY_CLASSES)
    if comparable:
        record["capacityComparable"] = True
        rel = abs(record["capacityEmitted"] - capacity_source) / capacity_source
        record["capacityRelativeDeviation"] = rel
        if rel > capacity_tolerance:
            record["flags"].append(
                f"capacity {record['capacityEmitted']:.9g} cm^3 differs from the source's "
                f"{capacity_source:.9g} cm^3 by {rel:.3g} relative")

    record["contains"] = contains_crosscheck(source_shape, emitted_shape, placement, n_points,
                                             seed, skin, max_report, mirrored, one_body)
    if record["contains"]["mismatches"]:
        record["failures"].append(
            f"{record['contains']['mismatches']} containment disagreement(s) over "
            f"{record['contains']['points']} scored point(s)")
    if record["contains"]["points"] == 0:
        record["failures"].append("no point was scored: the comparison is empty")
    if one_body and record["contains"]["insideEmitted"] == 0:
        # Without this a one-way comparison would be passed by a body that encloses nothing.
        record["failures"].append(
            f"the emitted body encloses none of the {record['contains']['points']} scored "
            "point(s): the one-way comparison is empty")
    return record


# ------------------------------------------------------------------------------------------
# driving a converter output directory
# ------------------------------------------------------------------------------------------

def _writer_index(writer_report):
    """emittedName -> the writer's row, which carries the source volume's real `name`.

    A volume that has daughters is written as an assembly plus a separate solid for the mother's
    own body, named `<name>__body`. That body is not a row of its own: the parent's row names it
    in `bodyComponent`, and the parent's `shapeClass` and `capacity_cm3` *are* the body's, because
    the body is exactly the mother's own solid. Indexing it here is what lets 152 of the 714 rows
    on these three corpora be scored at all.
    """
    index = {}
    for row in writer_report.get("volumes", []):
        emitted = row.get("emittedName") or row.get("name")
        if emitted:
            index[emitted] = row
        body = row.get("bodyComponent")
        if body and body != emitted:
            index[body] = row
    return index


def _placed_box(placement, shape):
    """A shape's axis-aligned box in the part frame: `(origin, half)`."""
    origin = [shape.GetOrigin()[i] for i in range(3)]
    half = [shape.GetDX(), shape.GetDY(), shape.GetDZ()]
    if placement is None:
        return origin, half
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    for sx in (-1.0, 1.0):
        for sy in (-1.0, 1.0):
            for sz in (-1.0, 1.0):
                local = (origin[0] + sx * half[0], origin[1] + sy * half[1],
                         origin[2] + sz * half[2])
                for i in range(3):
                    v = sum(placement[i][c] * local[c] for c in range(3)) + placement[i][3]
                    lo[i] = min(lo[i], v)
                    hi[i] = max(hi[i], v)
    return [0.5 * (lo[i] + hi[i]) for i in range(3)], [0.5 * (hi[i] - lo[i]) for i in range(3)]


def resolve_source_volume(candidates, row, emitted_shape=None, placement=None):
    """Which of several volumes sharing one name the writer's row refers to.

    A ROOT geometry may hold distinct volumes under the same name -- ITS has 21 of them, and PIPE
    has two `bellows2LowerPlie` -- and the writer disambiguates the *emitted* names (`name#2`)
    without recording an address, so the volume has to be identified by what it looks like.

    The **bounding box** decides, not the capacity. `TGeoCompositeShape::Capacity()` is a
    Monte-Carlo estimate, and PIPE's two `bellows2LowerPlie` are two halves of one bellows whose
    capacities differ by 0.09 % -- well inside that sampling noise, and the writer's own two rows
    for them differ from the geometry's by more than they differ from each other. Ranking on that
    number picked the wrong half, and the containment check then reported 601 disagreements with
    zero points inside the emitted shape, which is what a comparison against the mirror-image
    sibling looks like. A bounding box is exact for every `TGeoShape`, composites included, so it
    separates the two cleanly; capacity is kept only as a tie-break, and only where it is analytic.
    """
    if len(candidates) == 1:
        return candidates[0]
    wanted_class = row.get("shapeClass")
    wanted_capacity = row.get("capacity_cm3")
    sampled = wanted_class in _SAMPLED_CAPACITY_CLASSES
    want_box = (_placed_box(placement, emitted_shape)
                if emitted_shape is not None else None)
    best, best_key = None, None
    for volume in candidates:
        shape = volume.GetShape()
        if wanted_class and shape.ClassName() != wanted_class:
            continue
        box_score = 0.0
        if want_box is not None:
            here = _placed_box(None, shape)
            box_score = max(max(abs(here[0][i] - want_box[0][i]) for i in range(3)),
                            max(abs(here[1][i] - want_box[1][i]) for i in range(3)))
        capacity_score = (0.0 if (wanted_capacity is None or sampled)
                          else abs(shape.Capacity() - wanted_capacity)
                          / max(abs(wanted_capacity), 1.0e-30))
        key = (box_score, capacity_score)
        if best_key is None or key < best_key:
            best, best_key = volume, key
    return best


def check_run(original, writer_report_path, converted, n_points=DEFAULT_POINTS,
              seed=DEFAULT_SEED, skin=DEFAULT_SKIN_CM, capacity_tolerance=CAPACITY_TOLERANCE,
              profile_tolerance=PROFILE_TOLERANCE, max_report=5, verbose=True):
    """Compare every CSG-carried part of a converter output against its source volume."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    converted = Path(converted)
    csg_report_path = converted / "csg_report.json"
    if not csg_report_path.exists():
        raise SystemExit(f"{csg_report_path} does not exist (convert with --csg auto)")
    csg_report = json.loads(csg_report_path.read_text())
    writer_report = json.loads(Path(writer_report_path).read_text())
    index = _writer_index(writer_report)

    manager = ROOT.TGeoManager.Import(str(original))
    if manager is None:
        raise SystemExit(f"could not read a TGeoManager from {original}")
    by_name = {}
    for volume in manager.GetListOfVolumes():
        by_name.setdefault(volume.GetName(), []).append(volume)

    records = []
    open_files = []
    for part in csg_report.get("parts", []):
        if part.get("representation") != "csg":
            continue
        emitted_name = part.get("volume")
        stub = {"part": part.get("part"), "volume": emitted_name, "failures": [], "flags": []}
        row = index.get(emitted_name)
        if row is None:
            stub["failures"].append(f"no writer-report row for emittedName {emitted_name!r}")
            records.append(stub)
            continue
        shape_file = part.get("shapeFile")
        if not shape_file or not Path(shape_file).exists():
            stub["failures"].append(f"shapeFile {shape_file!r} does not exist")
            records.append(stub)
            continue
        handle = ROOT.TFile.Open(str(shape_file))
        open_files.append(handle)
        emitted_shape = handle.Get("shape")
        if not emitted_shape:
            stub["failures"].append(f"{shape_file} carries no object under the key \"shape\"")
            records.append(stub)
            continue
        if not reclose_flat_csg(emitted_shape):
            stub["failures"].append(
                f"{shape_file}: O2FlatCSG::CloseShape refused the shape after reading it, so "
                "its sub-cell boxes could not be rebuilt")
            records.append(stub)
            continue
        # The emitted shape is read *before* the source volume is resolved, because where several
        # volumes share a name its bounding box is what tells them apart.
        candidates = by_name.get(row.get("name")) or []
        source_volume = (resolve_source_volume(candidates, row, emitted_shape,
                                               part.get("shapePlacement"))
                         if candidates else None)
        if source_volume is None:
            stub["failures"].append(
                f"the original geometry has no volume named {row.get('name')!r} whose shape "
                "matches the writer's record")
            records.append(stub)
            continue
        record = check_part(part, row, source_volume.GetShape(), emitted_shape,
                            part.get("shapePlacement"), n_points, seed, skin,
                            capacity_tolerance, profile_tolerance, max_report)
        records.append(record)
        if verbose:
            print_record(record)

    n_fail = sum(1 for r in records if r["failures"])
    n_flag = sum(1 for r in records if r.get("flags"))
    if verbose:
        worst_capacity = max([r["capacityRelativeDeviation"] for r in records
                              if r.get("capacityRelativeDeviation") is not None] or [0.0])
        worst_profile = max([r["profileDeviationCm"] for r in records
                             if r.get("profileDeviationCm") is not None] or [0.0])
        print(f"\n{len(records) - n_fail}/{len(records)} CSG part(s) agree with their source "
              f"TGeoShape ({n_fail} failure(s), {n_flag} flag(s))")
        print(f"worst capacity deviation {worst_capacity:.3g} relative, worst polycone profile "
              f"deviation {worst_profile:.3g} cm")
    for handle in open_files:
        handle.Close()
    return records, n_fail, n_flag


def print_record(record):
    if record["failures"]:
        print(f"  [FAIL] {record['volume']}: " + "; ".join(record["failures"]))
        for example in (record.get("contains") or {}).get("examples", []):
            print(f"           at {example['point']} (local {example['local']}): "
                  f"source {'in' if example['source'] else 'out'}, "
                  f"emitted {'in' if example['emitted'] else 'out'}")
        return
    bits = [record["emittedClass"]]
    if record.get("mirrored"):
        bits.append("mirrored prototype")
    if record.get("oneBodyOfMany"):
        bits.append("one body of many, scored one-way")
    if record.get("classComparable"):
        bits.append("class matches" if record["classMatches"] else "class differs")
    if record.get("profileDeviationCm") is not None:
        bits.append(f"profile {record['profileDeviationCm']:.3g} cm")
    if record.get("capacityComparable"):
        bits.append(f"capacity rel {record['capacityRelativeDeviation']:.3g}")
    else:
        bits.append("capacity not comparable")
    contains = record.get("contains") or {}
    bits.append(f"Contains {contains.get('mismatches')}/{contains.get('points')} "
                f"({contains.get('skipped')} on the skin)")
    marker = "flag" if record.get("flags") else "ok  "
    print(f"  [{marker}] {record['volume']}: " + ", ".join(bits))
    for flag in record.get("flags", []):
        print(f"           flag: {flag}")


# ------------------------------------------------------------------------------------------
# self-test
# ------------------------------------------------------------------------------------------
#
# The fixtures are built in a subprocess and not here. ROOT keeps one global `gGeoManager`, and a
# process that creates a geometry and then imports another one does not survive interpreter
# shutdown -- so the checker itself only ever imports, and the self-test's source geometry is made
# next door.

_FIXTURE_BUILDER = r"""
import json, sys
from pathlib import Path
import ROOT
ROOT.gROOT.SetBatch(True)

folder = Path(sys.argv[1])
sections = [(-5.0, 1.0, 3.0), (0.0, 1.0, 3.0), (0.0, 2.0, 4.0), (5.0, 2.0, 4.0)]
names = ("GOOD", "PLACED", "MIRRORED", "TWIN", "BAD", "WRONGCLASS", "TWOBODY")


def halves_shape(zlow):
    # One half of a tube, as a composite, so its Capacity() is a Monte-Carlo estimate.
    tag = "lo" if zlow < 0.0 else "hi"
    tube = ROOT.TGeoTube("halves_t_" + tag, 0.0, 2.0, 1.0)
    slab = ROOT.TGeoBBox("halves_b_" + tag, 3.0, 3.0, 0.5)
    shift = ROOT.TGeoTranslation("halves_m_" + tag, 0.0, 0.0, zlow + 0.5)
    for obj in (tube, slab, shift):
        ROOT.SetOwnership(obj, False)
    node = ROOT.TGeoIntersection(tube, slab, ROOT.nullptr, shift)
    ROOT.SetOwnership(node, False)
    comp = ROOT.TGeoCompositeShape("halves_c_" + tag, node)
    ROOT.SetOwnership(comp, False)
    return comp


def make_pcon(name, rows, phi1=0.0, dphi=360.0):
    shape = ROOT.TGeoPcon(name, phi1, dphi, len(rows))
    for i, (z, rmin, rmax) in enumerate(rows):
        shape.DefineSection(i, z, rmin, rmax)
    ROOT.SetOwnership(shape, False)
    return shape


geometry = ROOT.TGeoManager("knownsource", "known-source self-test")
material = ROOT.TGeoMaterial("Vacuum", 0, 0, 0)
medium = ROOT.TGeoMedium("Vacuum", 1, material)
top = geometry.MakeBox("TOP", medium, 50.0, 50.0, 50.0)
geometry.SetTopVolume(top)
for name in names:
    volume = ROOT.TGeoVolume(name, make_pcon(name + "_sh", sections), medium)
    ROOT.SetOwnership(volume, False)
    top.AddNode(volume, 1)
# A second, different volume under a name that is already taken: the writer emits it as `TWIN#2`
# and the checker has to find *this* one, not the first.
twin = [(z, rmin, rmax + 1.0) for z, rmin, rmax in sections]
second = ROOT.TGeoVolume("TWIN", make_pcon("twin2_sh", twin), medium)
ROOT.SetOwnership(second, False)
top.AddNode(second, 1)
# Two COMPOSITES of one name, the two halves of a tube. Their Capacity() is Monte-Carlo and the
# two draws land within noise of each other, so only a bounding box can tell them apart.
for zlow in (-1.0, 0.0):
    half = ROOT.TGeoVolume("HALVES", halves_shape(zlow), medium)
    ROOT.SetOwnership(half, False)
    top.AddNode(half, 1)
geometry.CloseGeometry()
geometry.Export(str(folder / "source_geometry.root"))


def capacity(rows):
    return make_pcon("cap_probe", rows).Capacity()


rows = [{"name": n, "emittedName": n, "shapeClass": "TGeoPcon", "mirrored": n == "MIRRORED",
         "capacity_cm3": capacity(sections)} for n in names]
rows.append({"name": "TWIN", "emittedName": "TWIN#2", "shapeClass": "TGeoPcon",
             "mirrored": False, "capacity_cm3": capacity(twin)})
# The writer records a capacity for a composite too, and it is a different Monte-Carlo draw from
# the geometry's own -- which is exactly why the checker must not rank on it.
rows.append({"name": "HALVES", "emittedName": "HALVES", "shapeClass": "TGeoCompositeShape",
             "mirrored": False, "capacity_cm3": halves_shape(0.0).Capacity()})
(folder / "writer_report.json").write_text(json.dumps({"volumes": rows}))


def write_shape(name, shape, placement=None):
    target = folder / ("shape_%s.root" % name.replace("#", "_"))
    out = ROOT.TFile.Open(str(target), "RECREATE")
    out.WriteTObject(shape, "shape")
    out.Close()
    return {"part": name, "volume": name, "representation": "csg",
            "shapeFile": str(target), "shapePlacement": placement}


parts = [write_shape("GOOD", make_pcon("good_sh", sections))]
# A placement that is load-bearing: the profile sits 7 cm up in the shape's own frame and the
# placement brings it back, so only a consumer that composes the two agrees with the source.
shifted = [(z + 7.0, rmin, rmax) for z, rmin, rmax in sections]
parts.append(write_shape("PLACED", make_pcon("placed_sh", shifted),
                         [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, -7.0]]))
# The Z-mirrored prototype the writer emits for a volume placed by a reflecting matrix.
mirrored = [(-z, rmin, rmax) for z, rmin, rmax in reversed(sections)]
parts.append(write_shape("MIRRORED", make_pcon("mirrored_sh", mirrored)))
parts.append(write_shape("TWIN", make_pcon("twin_sh", sections)))
parts.append(write_shape("TWIN#2", make_pcon("twin2_out_sh", twin)))
# The emitted body is the LOWER half; the checker must resolve to the lower source volume.
parts.append(write_shape("HALVES", halves_shape(-1.0)))
wrong = [(z, rmin, rmax + (0.05 if i == 3 else 0.0))
         for i, (z, rmin, rmax) in enumerate(sections)]
parts.append(write_shape("BAD", make_pcon("bad_sh", wrong)))
tube = ROOT.TGeoTube("wrongclass_sh", 1.0, 4.0, 5.0)
ROOT.SetOwnership(tube, False)
parts.append(write_shape("WRONGCLASS", tube))
# One body of a CAD label that carried two. The source is the whole label; this part is the upper
# half of its profile, which is inside the source everywhere and covers only half of it.
upper = [(0.5, 2.0, 4.0), (5.0, 2.0, 4.0)]
body = write_shape("TWOBODY_b2", make_pcon("twobody_sh", upper))
body["volume"] = "TWOBODY"
parts.append(body)
# ... and the negative control for the one-way rule: a body that sticks OUT of its own label must
# still fail, or the rule would pass anything.
outside = [(0.5, 2.0, 5.0), (5.0, 2.0, 5.0)]
spill = write_shape("TWOBODYBAD_b2", make_pcon("twobodybad_sh", outside))
spill["volume"] = "TWOBODY"
parts.append(spill)
(folder / "csg_report.json").write_text(json.dumps({"parts": parts}))

# Every file is on disk. Skip the interpreter's teardown: ROOT's global geometry does not survive
# it once a closed TGeoManager has had shapes built against it, and the crash would be reported as
# a failure of a fixture that in fact succeeded.
import os
sys.stdout.flush()
os._exit(0)
"""


def self_test(verbose=True, workdir=None):
    """A geometry, a writer report and a converter output built here, with a known verdict.

    The negative controls are the point of it: a checker that cannot report a disagreement passes
    a positive-only suite exactly as well as one that works. There are three -- a displaced
    radius, a shape of the wrong class, and the placed control re-run with its placement thrown
    away, which is the mistake a consumer that forgets to compose the transform makes.

    The working folder is left behind under the system temporary directory; it is small, and it is
    the evidence for a failure.
    """
    import subprocess
    import tempfile

    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))
        if verbose:
            print(f"  [{'ok ' if condition else 'FAIL'}] {name}"
                  + (f"  {detail}" if detail else ""))

    folder = Path(tempfile.mkdtemp(prefix="knownsource_")) if workdir is None else Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    builder = folder / "_build_fixtures.py"
    builder.write_text(_FIXTURE_BUILDER)
    subprocess.run([sys.executable, str(builder), str(folder)], check=True,
                   stdout=subprocess.DEVNULL)

    records, n_fail, _n_flag = check_run(folder / "source_geometry.root",
                                         folder / "writer_report.json", folder,
                                         n_points=20000, verbose=False)
    by_name = {r["volume"]: r for r in records}

    good = by_name.get("GOOD", {})
    check("the positive control passes every comparable check",
          not good.get("failures") and not good.get("flags")
          and good.get("classMatches") is True and good.get("capacityComparable") is True
          and good.get("contains", {}).get("mismatches") == 0,
          f"failures {good.get('failures')}, flags {good.get('flags')}, capacity rel "
          f"{good.get('capacityRelativeDeviation')}, Contains "
          f"{good.get('contains', {}).get('mismatches')}/"
          f"{good.get('contains', {}).get('points')}")
    check("the positive control actually scored a useful point set",
          good.get("contains", {}).get("points", 0) > 1000,
          f"{good.get('contains', {}).get('points')} point(s) scored, "
          f"{good.get('contains', {}).get('skipped')} on the skin")

    placed = by_name.get("PLACED", {})
    check("a shape whose placement is composed correctly passes",
          not placed.get("failures") and placed.get("contains", {}).get("mismatches") == 0,
          f"failures {placed.get('failures')}")
    ignored = _placed_without_its_placement(folder)
    check("the placement is load-bearing: ignoring it must fail",
          ignored is not None and ignored["mismatches"] > 0,
          f"{ignored['mismatches'] if ignored else 'not run'} disagreement(s) with a null "
          "placement")

    mirrored = by_name.get("MIRRORED", {})
    check("a Z-mirrored prototype is compared through the mirror and passes",
          not mirrored.get("failures") and mirrored.get("mirrored") is True
          and mirrored.get("contains", {}).get("mismatches") == 0,
          f"failures {mirrored.get('failures')}")

    twin2 = by_name.get("TWIN#2", {})
    check("a name shared by two volumes resolves to the right one",
          not twin2.get("failures") and twin2.get("capacityComparable") is True
          and twin2.get("capacityRelativeDeviation") is not None
          and twin2["capacityRelativeDeviation"] < 1.0e-12,
          f"failures {twin2.get('failures')}, capacity rel "
          f"{twin2.get('capacityRelativeDeviation')}")

    halves = by_name.get("HALVES", {})
    check("two same-named composites are told apart by their box, not by a sampled capacity",
          not halves.get("failures")
          and halves.get("contains", {}).get("mismatches") == 0,
          f"failures {halves.get('failures')}, Contains "
          f"{halves.get('contains', {}).get('mismatches')}/"
          f"{halves.get('contains', {}).get('points')}")
    check("and their capacities really could not have decided it",
          _halves_capacities_are_indistinguishable(folder),
          "the two halves' Capacity() draws are within Monte-Carlo noise of each other")

    bad = by_name.get("BAD", {})
    check("the negative control is caught", bool(bad.get("failures")),
          "; ".join(bad.get("failures", [])) or "NOT CAUGHT")
    check("the negative control is caught by containment, not only by the profile",
          bad.get("contains", {}).get("mismatches", 0) > 0,
          f"{bad.get('contains', {}).get('mismatches')} disagreement(s)")
    check("the negative control's profile deviation is the displacement",
          bad.get("profileDeviationCm") is not None
          and abs(bad["profileDeviationCm"] - 0.05) < 1.0e-9,
          f"{bad.get('profileDeviationCm')}")

    wrongclass = by_name.get("WRONGCLASS", {})
    check("a shape of the wrong class is caught by the metrics, not only by its class",
          bool(wrongclass.get("failures"))
          and wrongclass.get("contains", {}).get("mismatches", 0) > 0,
          "; ".join(wrongclass.get("failures", [])) or "NOT CAUGHT")
    check("a class that differs without a geometric difference is a flag, not a failure",
          wrongclass.get("classMatches") is False and any(
              "is not the source's" in f for f in wrongclass.get("flags", [])),
          f"flags {wrongclass.get('flags')}")

    two_body = next((r for r in records if r["part"] == "TWOBODY_b2"), {})
    check("one body of a multi-body label passes on the one-way containment test",
          not two_body.get("failures") and two_body.get("oneBodyOfMany") is True
          and two_body.get("contains", {}).get("oneWay") is True
          and two_body.get("contains", {}).get("mismatches") == 0
          and two_body.get("contains", {}).get("insideEmitted", 0) > 100,
          f"failures {two_body.get('failures')}, "
          f"{two_body.get('contains', {}).get('insideEmitted')} point(s) inside the body")
    check("a multi-body part's class and capacity are reported as not comparable",
          two_body.get("capacityComparable") is False
          and two_body.get("classComparable") is False
          and any("multi-body" in f for f in two_body.get("flags", [])),
          f"flags {two_body.get('flags')}")
    spilled = next((r for r in records if r["part"] == "TWOBODYBAD_b2"), {})
    check("the one-way rule still catches a body that sticks out of its own label",
          bool(spilled.get("failures"))
          and spilled.get("contains", {}).get("mismatches", 0) > 0,
          "; ".join(spilled.get("failures", [])) or "NOT CAUGHT")

    check("the run reports exactly the three deliberately wrong parts as failures",
          n_fail == 3, f"{n_fail} failure(s) over {len(records)} part(s)")

    n_ok = sum(1 for _n, ok, _d in checks if ok)
    if verbose:
        print(f"  {n_ok}/{len(checks)} known-source self-checks passed (fixtures in {folder})")
    return n_ok, len(checks)


def _halves_capacities_are_indistinguishable(folder):
    """Are the two same-named composites' capacities within Monte-Carlo noise of each other?

    Without this the check above would prove nothing: if the two halves' capacities were far
    apart, ranking on them would have worked and the bounding-box rule would be untested. On
    PIPE's two `bellows2LowerPlie` the gap is 0.09 %, which is what this asserts in miniature.
    """
    import ROOT
    manager = ROOT.gGeoManager
    if not manager:
        return False
    capacities = [volume.GetShape().Capacity() for volume in manager.GetListOfVolumes()
                  if volume.GetName() == "HALVES"]
    if len(capacities) != 2:
        return False
    return abs(capacities[0] - capacities[1]) / max(capacities) < 0.02


def _placed_without_its_placement(folder):
    """Re-run the placed positive control with a null placement; it must then disagree.

    Uses the geometry `check_run` already imported rather than importing a second one, because
    two `TGeoManager`s in one process do not survive interpreter shutdown.
    """
    import ROOT
    manager = ROOT.gGeoManager
    if not manager:
        return None
    source = None
    for volume in manager.GetListOfVolumes():
        if volume.GetName() == "PLACED":
            source = volume.GetShape()
            break
    if source is None:
        return None
    handle = ROOT.TFile.Open(str(Path(folder) / "shape_PLACED.root"))
    emitted = handle.Get("shape")
    result = contains_crosscheck(source, emitted, None, 5000, DEFAULT_SEED, DEFAULT_SKIN_CM, 1)
    handle.Close()
    return result


# ------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original", type=Path,
                    help="the o2sim_geometry.root the STEP was written from")
    ap.add_argument("--writer-report", type=Path, dest="writer_report",
                    help="O2_TGeoToCAD.py's --report JSON, which maps emittedName -> name")
    ap.add_argument("--converted", type=Path,
                    help="the converter output folder (csg_report.json and shape_*.root)")
    ap.add_argument("--points", type=int, default=DEFAULT_POINTS,
                    help="containment samples per part; default %(default)s")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help="the fixed seed for those samples; default %(default)s")
    ap.add_argument("--skin", type=float, default=DEFAULT_SKIN_CM,
                    help="do not score points nearer than this to either boundary, in cm; "
                         "default %(default)s")
    ap.add_argument("--capacity-tolerance", type=float, default=CAPACITY_TOLERANCE,
                    dest="capacity_tolerance",
                    help="relative capacity agreement below which a part is flagged; "
                         "default %(default)s")
    ap.add_argument("--profile-tolerance", type=float, default=PROFILE_TOLERANCE,
                    dest="profile_tolerance",
                    help="polycone profile agreement demanded, relative to the part's diagonal; "
                         "default %(default)s, which is csg/recognise.REL_TOL")
    ap.add_argument("--strict", action="store_true",
                    help="treat capacity flags as failures too")
    ap.add_argument("--max-report", type=int, default=5, dest="max_report",
                    help="how many disagreeing points to print per part; default %(default)s")
    ap.add_argument("--json", type=Path, help="write the per-part records here")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        n_ok, n = self_test()
        print(f"\n{n_ok}/{n} known-source self-checks passed")
        return 0 if n_ok == n else 1

    if not (args.original and args.writer_report and args.converted):
        ap.error("give --original, --writer-report and --converted, or --self-test")

    records, n_fail, n_flag = check_run(args.original, args.writer_report, args.converted,
                                        n_points=args.points, seed=args.seed, skin=args.skin,
                                        capacity_tolerance=args.capacity_tolerance,
                                        profile_tolerance=args.profile_tolerance,
                                        max_report=args.max_report)
    if args.json:
        args.json.write_text(json.dumps(records, indent=1))
        print(f"Wrote {args.json}")
    return 1 if (n_fail or (args.strict and n_flag)) else 0


if __name__ == "__main__":
    sys.exit(main())
