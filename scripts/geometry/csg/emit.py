#!/usr/bin/env python3
"""Stream H: recognise CAD leaf solids as CSG, prove it, and emit `shape_<VOL>_<LID>.root`.

This is steps 2-4 of the MVP in `Tutorial.md` §6: the recogniser (`recognise.py`), the emitter
(here), and the two-test per-part acceptance.

    recognise  ->  OCCT symmetric-difference volume  ->  shape_<part>.root  ->  the oracle gate

A part is converted as CSG only if **both** tests pass, and they are independent by construction:
the symmetric difference is measured on an OCCT realisation of the description and the gate scores
a ROOT realisation of the same description (`primitives.py`). A builder bug in one and not the
other shows up as a disagreement between the two tests, which is a finding, not a pass.

Two ways in, one code path
--------------------------
  * `--db <gate workdir>/db` walks the `brep_*.brep` files an existing gate run already wrote,
    which are in cm and in the part's own local frame -- the very shapes the oracle answers
    about. Emitting next to them and re-scoring with `runOracleGate.py --skip-convert` is the
    loop Stream G §1 designed for.
  * `O2_CADtoTGeo.py --csg auto` calls `csg.hook`, which calls the same functions on the leaf
    solid it already has in memory.

The ROOT/pythonOCC question, answered
-------------------------------------
The alibuild Python 3.10 that pythonOCC is built against **is** the interpreter the O2
environment provides, so with `csg/occ_env.py` prepending pythonOCC's site-packages a single
process imports both `OCC` and `ROOT`. Measured, not assumed. The one place that does not hold is
inside `runOracleGate.py`, which replaces `PYTHONPATH`/`LD_LIBRARY_PATH` with OCC-only values for
the converter subprocess; a converter running there writes the description as
`csg_<VOL>_<LID>.json` and `--from-json` turns those into `.root` files afterwards.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from csg.occ_env import ensure_occ  # noqa: E402

ensure_occ()

from csg import accept, primitives as prim, recognise  # noqa: E402


# ------------------------------------------------------------------------------------------
# per-solid pipeline
# ------------------------------------------------------------------------------------------

def model_tolerance_cm(shape):
    """The shape's own statement about how well its boundary is defined, in cm.

    Same quantity `O2_CADtoTGeo.shape_model_tolerance()` and `occtOracle.shape_tolerance()`
    compute, and the same one the gate uses as its band. The shape here is already in cm, so
    there is no scale factor.
    """
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
    worst = 0.0
    for kind, getter in ((TopAbs_FACE, lambda s: BRep_Tool.Tolerance(topods.Face(s))),
                         (TopAbs_EDGE, lambda s: BRep_Tool.Tolerance(topods.Edge(s))),
                         (TopAbs_VERTEX, lambda s: BRep_Tool.Tolerance(topods.Vertex(s)))):
        exp = TopExp_Explorer(shape, kind)
        while exp.More():
            worst = max(worst, getter(exp.Current()))
            exp.Next()
    return worst


def _build_and_accept(solid, cand, tol, band_factor):
    """`(acceptance|None, reason|None)` for one candidate. Never raises on a bad candidate."""
    try:
        occ_shape = prim.build_occ(cand)
    except Exception as exc:                                     # noqa: BLE001
        return None, f"candidate failed to build in OCCT: {exc}"
    result = accept.symmetric_difference(solid, occ_shape, tol, band_factor)
    return result, (None if result.get("accepted") else result.get("reason"))


def process_solid(solid, name, tolerance=None, band_factor=1.0):
    """recognise -> build -> accept. Returns a record; `record['candidate']` is None if declined.

    A candidate the acceptance test **rejects** is retried once as a revolved profile
    (`recognise.recognise_revolved`), because the whole-part matchers answer before the revolved
    one and an all-cone stack looks like a single `TGeoCone` to them: they propose, the volume
    refuses, and without a retry the part ships one tier down. The retry runs only after a
    rejection, so a part whose first candidate is accepted never reaches it and its record and
    its candidate are what they were.
    """
    record = {"part": name, "recognised": False, "accepted": False, "candidate": None,
              "reason": None, "acceptance": None, "recogniser": None, "description": None}
    cand, reason = recognise.recognise(solid)
    if cand is None:
        record["reason"] = reason
        return record
    record["recognised"] = True
    record["recogniser"] = cand["recogniser"]
    record["description"] = prim.describe(cand)
    tol = model_tolerance_cm(solid) if tolerance is None else tolerance
    result, why_not = _build_and_accept(solid, cand, tol, band_factor)
    if result is not None:
        record["acceptance"] = result
        record["accepted"] = bool(result.get("accepted"))
    if record["accepted"]:
        record["candidate"] = cand
        return record
    record["reason"] = why_not

    alternative, alt_declined = recognise.recognise_revolved(solid)
    if alternative is None or alternative["recogniser"] == cand["recogniser"]:
        # Either the solid is not a revolved profile at all, or the revolved matcher is what
        # produced the candidate that was just refused; there is nothing else to try.
        if alternative is not None:
            return record
        record["reason"] = f"{why_not}; as a revolved profile: {alt_declined}"
        return record
    alt_result, alt_why_not = _build_and_accept(solid, alternative, tol, band_factor)
    if alt_result is not None and alt_result.get("accepted"):
        record["retriedAfter"] = {"recogniser": cand["recogniser"],
                                  "description": record["description"], "reason": why_not}
        record["recogniser"] = alternative["recogniser"]
        record["description"] = prim.describe(alternative)
        record["acceptance"] = alt_result
        record["accepted"] = True
        record["candidate"] = alternative
        record["reason"] = None
        return record
    record["reason"] = (f"{why_not}; retried as {alternative['recogniser']} "
                        f"({prim.describe(alternative)}): {alt_why_not}")
    return record


def write_shape_root(cand, path):
    """Write the description as `shape_<part>.root`, per the Stream G convention.

    One object inheriting from `TGeoShape`, under the key `shape`, in cm; and -- since
    `Stream_N_PlacedPrimitives.md` -- optionally a `TGeoHMatrix` under the key `placement`, the
    rigid transform that takes the shape from its own canonical frame into the part's local
    frame. **A file with no `placement` key means the identity**, which is what every file written
    before that change means, so nothing older has to be rewritten and nothing older changes
    meaning.

    This is the same pair `o2::base::harness::saveShapeToRootFile` writes -- the C++ side is the
    authority and the unit test round-trips through it.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    shape, placement = prim.build_root(cand, "shape")
    out = ROOT.TFile.Open(str(path), "RECREATE")
    out.WriteTObject(shape, "shape")
    matrix = prim.root_placement_matrix(placement, "placement")
    if matrix is not None:
        out.WriteTObject(matrix, "placement")
    out.Close()
    return shape


def _occ_bbox(shape):
    from OCC.Core.Bnd import Bnd_Box
    from OCC.Core.BRepBndLib import brepbndlib
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    # OCCT reports a bounding box enlarged by the shape's own tolerance ("gap"); ROOT's is tight.
    # Comparing them without removing the gap leaves a constant 1e-07 cm floor -- which is
    # exactly the number Stream G §1 reports as the `shape` representation's typical
    # `bboxDeviationFromOracle`, since the oracle's box carries the same gap.
    box.SetGap(0.0)
    return box.Get()


def crosscheck_bbox(cand, occ_shape=None):
    """Max deviation, in cm, between the ROOT realisation's bounding box and the OCCT one.

    **Read this number knowing what it can and cannot say.** For an unplaced primitive it is exact
    and a frame error moves it by the size of the error. Where a rigid transform is involved it is
    not: the ROOT box has to be carried into the part frame by transforming its eight corners and
    taking the axis-aligned hull of those, which for a rotated body is strictly larger than the
    body. The same is true of a `TGeoCompositeShape`, whose `TGeoBoolNode::ComputeBBox` does
    exactly that internally -- on the six Bagger rams the inflation is 0.13-0.62 cm and is a
    property of ROOT's bounding box, not of the geometry. Use `crosscheck_contains` for a sharp
    check.
    """
    occ_shape = occ_shape if occ_shape is not None else prim.build_occ(cand)
    xmin, ymin, zmin, xmax, ymax, zmax = _occ_bbox(occ_shape)
    shape, placement = prim.build_root(cand, "bboxprobe")
    origin = [shape.GetOrigin()[i] for i in range(3)]
    half = [shape.GetDX(), shape.GetDY(), shape.GetDZ()]
    lo_root = [origin[i] - half[i] for i in range(3)]
    hi_root = [origin[i] + half[i] for i in range(3)]
    if placement is not None:
        lo_root, hi_root = _placed_box(placement, lo_root, hi_root)
    worst = 0.0
    for i, (lo, hi) in enumerate(((xmin, xmax), (ymin, ymax), (zmin, zmax))):
        worst = max(worst, abs(lo_root[i] - lo), abs(hi_root[i] - hi))
    return worst


def _placed_box(placement, lo, hi):
    """The axis-aligned hull, in the part frame, of a local box under a rigid placement."""
    out_lo = [float("inf")] * 3
    out_hi = [float("-inf")] * 3
    for ix in (lo[0], hi[0]):
        for iy in (lo[1], hi[1]):
            for iz in (lo[2], hi[2]):
                for i in range(3):
                    v = (placement[i][0] * ix + placement[i][1] * iy + placement[i][2] * iz
                         + placement[i][3])
                    out_lo[i] = min(out_lo[i], v)
                    out_hi[i] = max(out_hi[i], v)
    return out_lo, out_hi


def crosscheck_contains(cand, original, n_points=4000, seed=1234):
    """Classify random points against the **original CAD solid** and against the emitted shape.

    This is the sharp form of the frame check, and the one that catches a transposed rotation
    that the symmetric difference cannot see: the symmetric difference measures the *OCCT*
    realisation of the description, this measures the *ROOT* one, and it measures it against the
    CAD body rather than against the other realisation. Points within one model tolerance of the
    boundary are skipped, since neither side claims to decide those.
    """
    import random
    from array import array
    from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier
    from OCC.Core.TopAbs import TopAbs_IN, TopAbs_ON
    from OCC.Core.gp import gp_Pnt

    xmin, ymin, zmin, xmax, ymax, zmax = _occ_bbox(original)
    pad = 0.05 * max(xmax - xmin, ymax - ymin, zmax - zmin)
    shape, placement = prim.build_root(cand, "containsprobe")
    tol = max(model_tolerance_cm(original), 1.0e-9)
    classifier = BRepClass3d_SolidClassifier(original)
    rng = random.Random(seed)
    disagreements = 0
    scored = 0
    for _ in range(n_points):
        p = (rng.uniform(xmin - pad, xmax + pad), rng.uniform(ymin - pad, ymax + pad),
             rng.uniform(zmin - pad, zmax + pad))
        classifier.Perform(gp_Pnt(*p), tol)
        state = classifier.State()
        if state == TopAbs_ON:
            continue
        scored += 1
        # The point is in the part frame; the shape answers in its own. Composing the placement
        # here is the same composition every other consumer performs, so a wrong placement is a
        # disagreement against the CAD body rather than a silent pass.
        local = prim.placement_to_local(placement, p)
        if bool(shape.Contains(array("d", list(local)))) != (state == TopAbs_IN):
            disagreements += 1
    return {"points": scored, "disagreements": disagreements}


# ------------------------------------------------------------------------------------------
# driving a converter output directory
# ------------------------------------------------------------------------------------------

def load_brep(path):
    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.BRepTools import breptools
    from OCC.Core.TopAbs import TopAbs_SOLID
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import TopoDS_Shape, topods
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    if not breptools.Read(shape, str(path), builder):
        raise RuntimeError(f"failed to read {path}")
    solids = []
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    while exp.More():
        solids.append(topods.Solid(exp.Current()))
        exp.Next()
    if len(solids) != 1:
        return shape, len(solids)
    return solids[0], 1


def run_db(db_dir, write_root=True, band_factor=1.0, quiet=False):
    db_dir = Path(db_dir)
    breps = sorted(db_dir.glob("*/brep_*.brep")) or sorted(db_dir.glob("brep_*.brep"))
    if not breps:
        raise SystemExit(f"no brep_*.brep under {db_dir}")
    records = []
    for brep in breps:
        suffix = brep.name[len("brep_"):-len(".brep")]
        part = f"{brep.parent.name}/{suffix}"
        solid, n_solids = load_brep(brep)
        record = process_solid(solid, part, band_factor=band_factor)
        record["brep"] = str(brep)
        record["nSolids"] = n_solids
        if record["accepted"] and write_root:
            target = brep.parent / f"shape_{suffix}.root"
            write_shape_root(record["candidate"], target)
            record["shape"] = str(target)
            record["bboxRootVsOcctCm"] = crosscheck_bbox(record["candidate"])
            record["containsCrosscheck"] = crosscheck_contains(record["candidate"], solid)
        json_target = brep.parent / f"csg_{suffix}.json"
        json_target.write_text(json.dumps(
            {"part": part, "candidate": record["candidate"], "acceptance": record["acceptance"],
             "recogniser": record["recogniser"]}, indent=1))
        records.append(record)
        if not quiet:
            _print_record(record)
    return records


def from_json(folder, quiet=False):
    """Turn every accepted `csg_<part>.json` in a folder into its `shape_<part>.root`.

    This is the second half of the split the environment forces: a converter running under
    `runOracleGate.py` has pythonOCC but not PyROOT, so it writes the description and this
    completes it. Nothing is re-recognised and nothing is re-accepted -- the description and its
    evidence are taken as given, which is the point of having a description at all.
    """
    folder = Path(folder)
    files = sorted(folder.glob("csg_*.json")) or sorted(folder.glob("*/csg_*.json"))
    written = []
    for path in files:
        payload = json.loads(path.read_text())
        if not payload.get("candidate"):
            continue
        suffix = path.name[len("csg_"):-len(".json")]
        target = path.parent / f"shape_{suffix}.root"
        shape = write_shape_root(payload["candidate"], target)
        written.append(target)
        if not quiet:
            print(f"  wrote {target} ({shape.ClassName()})")
    if not quiet:
        print(f"{len(written)} shape file(s) written from {len(files)} description(s)")
    return written


def _print_record(record):
    if record["accepted"]:
        acc = record["acceptance"]
        extra = ""
        if record.get("containsCrosscheck") is not None:
            cc = record["containsCrosscheck"]
            extra = (f", ROOT-vs-CAD Contains {cc['disagreements']}/{cc['points']}"
                     f", bbox(ROOT vs OCCT) {record['bboxRootVsOcctCm']:.2e} cm")
        print(f"  [CSG ] {record['part']}: {record['description']}  "
              f"[{record['recogniser']}]  dV_sym={acc['symmetricDifference']:.3g} cm^3 "
              f"(band {acc['band']:.3g}, rel {acc['relativeToVolume']:.2e}){extra}")
    elif record["recognised"]:
        print(f"  [rej ] {record['part']}: {record['description']} rejected -- {record['reason']}")
    else:
        print(f"  [decl] {record['part']}: {record['reason']}")


def summarise(records):
    n_csg = sum(1 for r in records if r["accepted"])
    n_rej = sum(1 for r in records if r["recognised"] and not r["accepted"])
    n_dec = sum(1 for r in records if not r["recognised"])
    print(f"\n{n_csg}/{len(records)} part(s) accepted as CSG "
          f"({n_rej} recognised but rejected by the symmetric difference, {n_dec} declined "
          f"by the recogniser)")
    return n_csg, n_rej, n_dec


# ------------------------------------------------------------------------------------------
# self-test
# ------------------------------------------------------------------------------------------

# SHA-256 of `json.dumps(candidate, sort_keys=True)` for every whole-part fixture below, taken
# from the tree before the revolved matcher and the acceptance retry existed. The floor this
# project actually has to hold is not "the suite is green" but "a part that converts today
# converts to the same bytes tomorrow", and only a recorded digest can say that. A digest here
# may be updated only together with a measured statement about which artefacts moved.
_CANDIDATE_DIGESTS_BEFORE_THE_REVOLVED_MATCHER = {
    "box":
        "99ad3ce28270c7aaef6aa08ea6d5d06a5d593384aca04f1404465801a82c968f",
    "solid cylinder":
        "7f3cd91ca1c95acb077e83c43cc978c11bdde93c2f9953617819d56d99eabb3d",
    "tube":
        "00902589f306881e06554ad039dc425451d620a0f7cc11b54eda0b94ca336cd6",
    "tube segment":
        "909afffc42726d5aa6f7f71cc01005c8a44d8f2c5fa8ac561ce8b2abae444263",
    "cone":
        "01392dc4dfbcd331a8e24a3278f1f3eee793fda0aea29e05da0771df18b279e6",
    "sphere":
        "2c5a7be0c2ee7c6126db54376e627451dccad780b1436504f2792cf99e12dc42",
    "placed tube":
        "47190bb1666e38d64f500f7e46350d62a57d2c8acf659184b986a5edf3e4e441",
    "rod-and-eye (two-cluster union)":
        "2fd3760396fb15edb0b0ad10cb99df6070678388460a7813f2501f344fd713c9",
}


# The same floor for rung 2's own emissions: SHA-256 of `json.dumps(candidate, sort_keys=True)`
# for every prism-family fixture, recorded when the matcher landed. A row here may be updated only
# together with a measured statement about which artefacts moved.
_PRISM_CANDIDATE_DIGESTS = {
    "L-shaped plate":
        "ef0c21358c50346fc6f1d19b801303791d63dbe4053538ffceff721bd7f411ba",
    "hollow 8-edge polygon (TGeoPgon)":
        "8b1b3ec817347271a435447ff843af5e42c14d42294c947f2b44784ac82a29da",
    "hollow 48-edge polygon (TGeoPgon)":
        "7dda853698a50f5d4ec9b27a62b26f80e6ebe1a84b8385563838944367e0e34a",
    "Trd1 (slanted x faces)":
        "b431f79ac92bfe13a7feed4274c300b82b7875b9fe168cfa3becf9c3001e1b3b",
    "Trd1 (taper reversed)":
        "68aa7c84b337248cddb582549445c4f7bed450b60edd2f872101166b54a79d2d",
    "Trd1 (TPC_IRB1's 0.5 % slant)":
        "26236a3814890b3441b243512c76ae135440fd2e658bea506223963b6ba86bf9",
    "Trd2 (both half-widths vary)":
        "e5b7757fda5785561119641cb371f786ec1ec13d4e9f38bd4a0725da9da33071",
    "Trd2 (isotropic taper, also a legal Xtru)":
        "c572219aa3ae096a62eeadc483cd966d840fa9cba04b3276aea65f7323c8fa3e",
    "Arb8 (parallelepiped)":
        "3532421c9f89a5ecd56d890ad321843608e7d25e7316ca39cc5eae7356dda322",
    "Arb8 (TPC_IHSTR's trapezoidal prism)":
        "3403de111d8ee4f6b7ff05fade7512f5208e92fb68807b7c0a6762eaa6f0bf27",
    "Arb8 (sheared in x only)":
        "cb46f4bad49306bb733634d2cf66e695ddd7e1aa7cc779562f33493a416758ac",
    "Arb8 (a TGeoTrap's eight corners)":
        "880bda4e01fc55380db99d46fbd139b455fbdd3d31e165f553bf5fd0620ca4ec",
    "Xtru (non-convex L section)":
        "3d3e1a41938efbef89e7bc0417dc312dc96c96c31866e6f95b0eb9f4b71f27f8",
    "Xtru (ITS ConeARibVol0's eight-corner section)":
        "991006456f492b5fe4f65b3da99b8cecd369c74c02484b6420705aa6c722a667",
    "Xtru (a triangular section)":
        "ce0549763513bc6a210183ac65614906e65d8b5e9565d36022934c7740fa7cb0",
    "Xtru (three sections, offset and scaled)":
        "0903971bfb0a916e10502b17aed38fb6e0be24ae1179096ec7ac2a31b527b37f",
    "Pgon (solid hexagonal prism)":
        "618c29e16c6396b50de9706150cc546fb83460c1085ed1dabf15075fc98ceab5",
    "Pgon (tapered eight-edge prism)":
        "782ba06f5261b5b9580ef34f7c390595e1a6ceb4fd560943821035e087e63d43",
    "Pgon (hollow 8-edge prism)":
        "e3db6ffd7d881c2e106c51cc0bfb65c6579399c066ea44a65b248262d5db0f22",
    "Pgon (hollow 48-edge prism)":
        "170bbe1bad0f51b45a67d34dde81e963063f5e233950f909518bd73dd707151d",
    "Pgon (TPC_Strip's thin 18-edge shell)":
        "878dac70124d695802949cef6e4326856d97cdee20d6f53d88d24dcbde03ce53",
    "Pgon (three hollow sections)":
        "b88d2382aea8dce7a35b6a76aaa1ac1f149d52242cea029c9ffc6ec3a5820a10",
    "Pgon (a 90 deg wedge closing on the axis)":
        "992f31e97db7c23aaff4b41d4bc1e1a7c5d728b10c15fc316ed133af4cdcfb13",
    "Pgon (a wedge across phi = 0)":
        "c3f918d0cc48576351c6ef5e63ff8724235181eb99630b79667ce7939bcbfbb2",
    "placed Trd1":
        "9b73e1a0c2a8d6cda0ed4f4790335dd7d2b0618d467d49e39a42cc7240a743d6",
    "placed Xtru (non-convex L section)":
        "3b69eaed1fa0e63b8aaa5ab2d5c5d4d2af9078535173d7bb8ea27d61e1a8e821",
}


def self_test(verbose=True, with_root=True):  # noqa: C901
    """Synthetic solids whose recognition and emission are known in closed form.

    Every positive case is paired with a negative one, because a recogniser that accepts
    everything and an acceptance test that rejects nothing both pass a positive-only suite. The
    ROOT half is checked by querying the emitted `TGeoShape` against the closed-form answer for
    the same point set, which is independent of both OCCT and the description's own builders.
    """
    import hashlib
    import math
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                         BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeSolid,
                                         BRepBuilderAPI_Sewing, BRepBuilderAPI_Transform)
    from OCC.Core.BRepFill import brepfill
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                                      BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism,
                                      BRepPrimAPI_MakeRevol, BRepPrimAPI_MakeSphere)
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.TopoDS import topods
    from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))
        if verbose:
            print(f"  [{'ok ' if condition else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))

    seen_digests = {}

    def expect(name, solid, want_recogniser, want_leaves=1):
        record = process_solid(solid, name)
        if record["accepted"]:
            seen_digests[name] = hashlib.sha256(
                json.dumps(record["candidate"], sort_keys=True).encode()).hexdigest()
        ok = record["accepted"] and record["recogniser"] == want_recogniser and \
            len(record["candidate"]["leaves"]) == want_leaves
        detail = (f"{record['recogniser']}: {record['description']}"
                  if record["recognised"] else f"declined: {record['reason']}")
        if record["recognised"] and not record["accepted"]:
            detail += f" -- rejected: {record['reason']}"
        check(f"{name} recognised as {want_recogniser} and accepted", ok, detail)
        return record

    def expect_declined(name, solid, needle=""):
        record = process_solid(solid, name)
        ok = not record["accepted"] and (needle in (record["reason"] or ""))
        check(f"{name} is not converted as CSG", ok, f"reason: {record['reason']}")
        return record

    ax = gp_Ax2(gp_Pnt(0, 0, -5), gp_Dir(0, 0, 1))

    # --- Tier 1, one per primitive the brief scopes ---
    expect("box", BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0).Shape(), "tier1-box")
    cyl = BRepPrimAPI_MakeCylinder(ax, 2.0, 10.0).Shape()
    expect("solid cylinder", cyl, "tier1-tube")
    bore = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -6), gp_Dir(0, 0, 1)), 1.0, 12.0).Shape()
    tube = BRepAlgoAPI_Cut(cyl, bore).Shape()
    expect("tube", tube, "tier1-tube")
    wedge = BRepPrimAPI_MakeCylinder(ax, 2.0, 10.0, math.radians(75.0)).Shape()
    seg = BRepAlgoAPI_Cut(wedge, bore).Shape()
    expect("tube segment", seg, "tier1-tubeseg")
    expect("cone", BRepPrimAPI_MakeCone(ax, 3.0, 1.0, 10.0).Shape(), "tier1-cone")
    expect("sphere", BRepPrimAPI_MakeSphere(gp_Pnt(1, 2, 3), 2.5).Shape(), "tier1-sphere")

    # A rotated, translated tube: the frame machinery, end to end.
    trsf = gp_Trsf()
    trsf.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), 0.7)
    shift = gp_Trsf()
    shift.SetTranslation(gp_Vec(3.0, -4.0, 5.0))
    moved = BRepBuilderAPI_Transform(tube, shift.Multiplied(trsf), True).Shape()
    moved_record = expect("placed tube", moved, "tier1-tube")

    # --- Tier 2, the Bagger ram in miniature: a rod through the wall of an eye ---
    eye = BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-0.75, 0, 0), gp_Dir(1, 0, 0)), 1.2, 1.5).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-1.0, 0, 0), gp_Dir(1, 0, 0)), 0.7, 2.0).Shape()
    ).Shape()
    rod_full = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 0.6, 8.0).Shape()
    rod = BRepAlgoAPI_Cut(rod_full, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(-0.75, 0, 0), gp_Dir(1, 0, 0)), 1.2, 1.5).Shape()).Shape()
    ram = BRepAlgoAPI_Fuse(eye, rod).Shape()
    ram_record = expect("rod-and-eye (two-cluster union)", ram, "tier2-tube-union", want_leaves=2)

    # --- the revolved profile: the shapes O2_TGeoToCAD.conv_pcon writes, read back ---
    # The fixture states the (r, z) ring itself rather than calling `primitives.pcon_profile_rz`,
    # so the profile the CAD is built from is not the profile the candidate is built from.
    def revolved(z, rmin, rmax, phi1=0.0, dphi=360.0):
        nz = len(z)
        ring = [(rmax[i], z[i]) for i in range(nz)]
        if all(r <= 0.0 for r in rmin):
            ring += [(0.0, z[nz - 1]), (0.0, z[0])]
        else:
            ring += [(rmin[i], z[i]) for i in range(nz - 1, -1, -1)]
        deduped = []
        for pt in ring:
            if deduped and abs(pt[0] - deduped[-1][0]) < 1e-12 \
                    and abs(pt[1] - deduped[-1][1]) < 1e-12:
                continue
            deduped.append(pt)
        poly = BRepBuilderAPI_MakePolygon()
        for (r, zz) in deduped:
            poly.Add(gp_Pnt(float(r), 0.0, float(zz)))
        poly.Close()
        rev = BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(poly.Wire()).Face(),
                                    gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                    math.radians(dphi))
        rev.Build()
        shape = rev.Shape()
        if abs(phi1) > 1e-12:
            spin = gp_Trsf()
            spin.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), math.radians(phi1))
            shape = BRepBuilderAPI_Transform(shape, spin, True).Shape()
        return shape

    def expect_pcon(name, solid, z, rmin, rmax, phi1=0.0, dphi=360.0):
        record = expect(name, solid, "revolved-pcon")
        if not record["accepted"]:
            return record
        p = record["candidate"]["leaves"][0]["params"]
        worst = max([abs(a - b) for a, b in zip(p["z"], z)]
                    + [abs(a - b) for a, b in zip(p["rmin"], rmin)]
                    + [abs(a - b) for a, b in zip(p["rmax"], rmax)]
                    + [abs(p["phi1"] - phi1), abs(p["dphi"] - dphi)]) \
            if len(p["z"]) == len(z) else float("inf")
        check(f"{name} reconstructs the source TGeoPcon parameters",
              len(p["z"]) == len(z) and worst < 1.0e-9,
              f"nz {len(p['z'])} vs {len(z)}, worst parameter deviation {worst:.3g}")
        return record

    # z-steps: duplicate z planes on both rmin and rmax, which the writer emits as cap annuli.
    step_z, step_rmin, step_rmax = [-5, 0, 0, 5], [1, 1, 2, 2], [3, 3, 4, 4]
    stepped = revolved(step_z, step_rmin, step_rmax)
    expect_pcon("stepped polycone (duplicate z planes)", stepped, step_z, step_rmin, step_rmax)
    # mixed cone and cylinder laterals on one axis -- the IBCYSSCone case, which the whole-part
    # matcher declines with "mixed lateral surface kinds".
    expect_pcon("cone and cylinder laterals on one axis",
                revolved([-5, 0, 5], [1, 1, 2], [2, 3, 3]), [-5, 0, 5], [1, 1, 2], [2, 3, 3])
    # rmin stepping through 0: the inner lateral is a cone that reaches the axis.
    expect_pcon("polycone whose rmin steps through 0",
                revolved([0, 5, 10], [0, 0, 2], [4, 4, 4]), [0, 5, 10], [0, 0, 2], [4, 4, 4])
    # a half turn, and a partial-phi wedge stated in absolute phi on an identity frame.
    expect_pcon("half-turn polycone", revolved([-5, 0, 5], [1, 1, 2], [2, 3, 3], 0.0, 180.0),
                [-5, 0, 5], [1, 1, 2], [2, 3, 3], 0.0, 180.0)
    expect_pcon("partial-phi stepped polycone",
                revolved(step_z, step_rmin, step_rmax, 10.0, 120.0),
                step_z, step_rmin, step_rmax, 10.0, 120.0)
    # a rotated, translated polycone: the frame machinery on a multi-section leaf.
    pcon_trsf = gp_Trsf()
    pcon_trsf.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), 0.7)
    pcon_shift = gp_Trsf()
    pcon_shift.SetTranslation(gp_Vec(3.0, -4.0, 5.0))
    pcon_place = pcon_shift.Multiplied(pcon_trsf)
    moved_pcon = BRepBuilderAPI_Transform(stepped, pcon_place, True).Shape()
    moved_pcon_record = expect("placed stepped polycone", moved_pcon, "revolved-pcon")
    check("a placed polycone travels as one leaf plus a rigid placement",
          moved_pcon_record["accepted"]
          and prim.placement_for_candidate(moved_pcon_record["candidate"]) is not None,
          "placement present" if moved_pcon_record["accepted"] else "not accepted")

    # --- negative controls: each must decline or be rejected ---
    # 1. a blind bore. This *is* a polycone -- z = [-5, 3, 3, 5], rmin = [1, 1, 0, 0] -- and the
    #    revolved matcher now converts it. Before this matcher existed it declined on its third
    #    cap plane, which is why the expectation moved rather than the shape.
    blind = BRepAlgoAPI_Cut(cyl, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, -6), gp_Dir(0, 0, 1)), 1.0, 9.0).Shape()).Shape()
    expect_pcon("cylinder with a blind bore", blind, [-5, 3, 3, 5], [1, 1, 0, 0], [2, 2, 2, 2])
    # 2. an L-shape: eight planes. This *is* a right prism on a six-corner section, and rung 2's
    #    prism matcher now converts it as a TGeoXtru. Before that matcher existed it declined on
    #    the plane count, which is why the expectation moved rather than the shape.
    ell = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 4.0, 4.0, 1.0).Shape(),
                          BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -1), 4.0, 4.0, 3.0).Shape()).Shape()
    expect("L-shaped plate", ell, "rung2-xtru")
    # 3. a torus: in scope for the surface solid, out of scope here, and it must say so.
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
    expect_declined("torus", BRepPrimAPI_MakeTorus(5.0, 1.0).Shape(), "toroidal")
    # 4. a cylinder with a flat milled off it: one cluster, but a plane that is neither cap nor
    #    wedge. The recogniser must decline rather than emit the round tube.
    flatted = BRepAlgoAPI_Cut(cyl, BRepPrimAPI_MakeBox(
        gp_Pnt(1.5, -3, -6), 3.0, 6.0, 12.0).Shape()).Shape()
    expect_declined("cylinder with a milled flat", flatted)

    # --- negative controls for the revolved matcher ---
    # 5. a TGeoPgon. `conv_pgon` writes its laterals as PLANES at the apothem radius, so a
    #    polygon is a pseudo-revolution: it looks like a polycone from a distance and is not one.
    #    Emitting it as a TGeoPcon would inflate every part by the polygon's own sagitta, so it
    #    must decline here; TGeoPgon emission is a separate recogniser.
    def prism_ring(apothem, nedges, phi1=0.0, dphi=360.0):
        dseg = math.radians(dphi) / nedges
        radius = apothem / math.cos(dseg / 2.0)
        n = nedges if abs(dphi - 360.0) < 1e-9 else nedges + 1
        return [(radius * math.cos(math.radians(phi1) + k * dseg),
                 radius * math.sin(math.radians(phi1) + k * dseg)) for k in range(n)]

    def swept_polygon(apothem, nedges, z0, z1):
        poly = BRepBuilderAPI_MakePolygon()
        for (x, y) in prism_ring(apothem, nedges):
            poly.Add(gp_Pnt(x, y, z0))
        poly.Close()
        pr = BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(poly.Wire()).Face(),
                                   gp_Vec(0, 0, z1 - z0))
        pr.Build()
        return pr.Shape()

    #    Rung 2 gave that recogniser: they are now converted as TGeoPgon, and the assertion that
    #    matters is that they are NOT converted as a polycone -- the class is the whole point.
    for nedges in (8, 48):
        pgon = BRepAlgoAPI_Cut(swept_polygon(3.0, nedges, -5.0, 5.0),
                               swept_polygon(1.5, nedges, -6.0, 6.0)).Shape()
        expect(f"hollow {nedges}-edge polygon (TGeoPgon)", pgon, "rung2-pgon")
    # 6. the same trap where it is hardest to see: polygonal laterals sharing an axis with a real
    #    cylinder, so there *is* an axis cluster and the planes reach the cap/wedge split.
    hybrid = BRepAlgoAPI_Fuse(
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -5), gp_Dir(0, 0, 1)), 3.0, 5.0).Shape(),
        swept_polygon(3.0, 6, 0.0, 5.0)).Shape()
    expect_declined("cylinder with a coaxial hexagonal section", hybrid,
                    "neither a cap nor a wedge")
    # 7. the near-miss the acceptance exists for: the recogniser cannot see a bore displaced by
    #    ten model tolerances off the axis and proposes the coaxial polycone anyway; the
    #    symmetric difference refuses it.
    for displacement in (1.0e-6, 1.0e-5):
        off = BRepAlgoAPI_Cut(
            revolved(step_z, [0, 0, 0, 0], step_rmax),
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(displacement, 0, -6), gp_Dir(0, 0, 1)),
                                     1.0, 12.0).Shape()).Shape()
        expect_declined(f"stepped polycone with the bore {displacement:g} cm off axis", off)
    # 8. a cap plane tilted off perpendicular: not a revolution any more, and it must say so
    #    rather than emit the polycone it nearly is.
    tilt = gp_Trsf()
    tilt.SetRotation(gp_Ax1(gp_Pnt(0, 0, 5), gp_Dir(1, 0, 0)), 1.0e-4)
    knife = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 4.9), gp_Dir(0, 0, 1)),
                                     10.0, 5.0).Shape()
    expect_declined("stepped polycone with a tilted top cap",
                    BRepAlgoAPI_Cut(stepped,
                                    BRepBuilderAPI_Transform(knife, tilt, True).Shape()).Shape(),
                    "neither a cap nor a wedge")

    # --- the instrument that scores the revolved candidate must be able to say "no" ---
    # `Stream_K_Tier0.md` §3: one measured quantity decides, so it is worth measuring that the
    # quantity moves. A profile with one radius displaced must be reported at that displacement.
    true_profile = prim.pcon_profile_rz({"z": [float(v) for v in step_z],
                                         "rmin": [float(v) for v in step_rmin],
                                         "rmax": [float(v) for v in step_rmax]})
    samples = [(3.0, -2.5), (4.0, 2.5), (1.0, -5.0), (2.0, 5.0), (3.5, 0.0)]
    check("the profile gap is zero on the profile's own boundary",
          recognise._profile_gap(true_profile, samples) < 1.0e-12,
          f"gap {recognise._profile_gap(true_profile, samples):.3g} cm")
    nudged_profile = [(r + (1.0e-6 if abs(r - 3.0) < 1e-12 else 0.0), z)
                      for (r, z) in true_profile]
    nudged_gap = recognise._profile_gap(nudged_profile, samples)
    check("the profile gap reports a radius displaced by ten model tolerances",
          abs(nudged_gap - 1.0e-6) < 1.0e-12, f"gap {nudged_gap:.3g} cm, expected 1e-06 cm")

    # --- the description must refuse an illegal TGeoPcon before either builder sees it ---
    for name, params in (
            ("unequal array lengths",
             {"phi1": 0.0, "dphi": 360.0, "z": [0.0, 1.0], "rmin": [0.0], "rmax": [1.0, 1.0]}),
            ("rmin above rmax",
             {"phi1": 0.0, "dphi": 360.0, "z": [0.0, 1.0], "rmin": [2.0, 2.0],
              "rmax": [1.0, 1.0]}),
            ("a single section",
             {"phi1": 0.0, "dphi": 360.0, "z": [0.0], "rmin": [0.0], "rmax": [1.0]}),
            ("z running backwards",
             {"phi1": 0.0, "dphi": 360.0, "z": [1.0, 0.0], "rmin": [0.0, 0.0],
              "rmax": [1.0, 1.0]})):
        try:
            prim.leaf("TGeoPcon", params, prim.identity_frame())
            refused = False
        except ValueError:
            refused = True
        check(f"a TGeoPcon description with {name} is refused", refused)

    # --- fix round 1: an all-cone stack, retried after the acceptance test refuses tier 1 ---
    # Two caps and two cone faces read as one TGeoCone to the whole-part matcher, which proposes a
    # cone that is not the solid; the volume refuses it and `process_solid` retries the revolved
    # matcher. Without the retry these ship one tier down.
    stack_record = expect_pcon("all-cone stack (two cones and two caps)",
                               revolved([-3, 0, 3], [0, 0, 0], [2, 3, 1]),
                               [-3, 0, 3], [0, 0, 0], [2, 3, 1])
    check("the all-cone stack was retried after tier 1 was rejected, not merely declined",
          (stack_record.get("retriedAfter") or {}).get("recogniser") == "tier1-cone",
          f"retried after {(stack_record.get('retriedAfter') or {}).get('recogniser')}: "
          f"{(stack_record.get('retriedAfter') or {}).get('reason')}")
    # An hourglass pinches to the axis in the middle, which is a legal polycone (rmax = 0 there)
    # and the case that made `BRepPrimAPI_MakeCone` raise on two identical radii.
    expect("hourglass (two cones meeting on the axis)",
           revolved([-5, 0, 5], [0, 0, 0], [2, 0, 2]), "revolved-pcon")

    # --- fix round 1: a two-section full-turn profile is said in its native class ---
    def expect_native(name, solid, want_recogniser, want_type, want_params):
        record = expect(name, solid, want_recogniser)
        if not record["accepted"]:
            return record
        lf = record["candidate"]["leaves"][0]
        worst = max(abs(lf["params"][k] - v) for k, v in want_params.items()) \
            if lf["type"] == want_type else float("inf")
        check(f"{name} emits a native {want_type} with the source's parameters",
              lf["type"] == want_type and worst < 1.0e-9,
              f"{lf['type']}, worst parameter deviation {worst:.3g}")
        return record

    # The ABSO pair: a TGeoCone with one radius constant arrives as one cylinder and one cone, so
    # the whole-part matcher declines it with "mixed lateral surface kinds" and the revolved
    # matcher must give the class back rather than a two-section polycone.
    expect_native("cone with a cylindrical bore (constant rmin)",
                  revolved([-25, 25], [4.5, 4.5], [16.22, 25.04]), "revolved-cone", "TGeoCone",
                  {"dz": 25.0, "rmin1": 4.5, "rmax1": 16.22, "rmin2": 4.5, "rmax2": 25.04})
    expect_native("cylinder with a conical bore (constant rmax)",
                  revolved([-3, 3], [6.99, 7.374], [26.02, 26.02]), "revolved-cone", "TGeoCone",
                  {"dz": 3.0, "rmin1": 6.99, "rmax1": 26.02, "rmin2": 7.374, "rmax2": 26.02})
    # ... and the two shapes that must NOT be canonicalised, because two sections do not describe
    # them: a wedge needs phi and a step needs the duplicated plane.
    expect_pcon("two-section wedge stays a polycone", revolved([-5, 5], [1, 1], [2, 3], 0.0,
                                                               120.0),
                [-5, 5], [1, 1], [2, 3], 0.0, 120.0)
    expect_pcon("a stepped profile stays a polycone", stepped, step_z, step_rmin, step_rmax)
    # The TGeoTube branch of the canonicalisation is not reachable from CAD -- a full-turn
    # two-section profile with both radii constant is a tube, and the whole-part matcher gets
    # there first -- so it is exercised on the description directly.
    tube_leaf, tube_tag = recognise._canonical_revolved_leaf(
        prim.leaf("TGeoPcon", {"phi1": 0.0, "dphi": 360.0, "z": [-4.0, 6.0],
                               "rmin": [1.0, 1.0], "rmax": [2.0, 2.0]}, prim.identity_frame()),
        (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 1.0e-9)
    check("a two-section profile with constant radii canonicalises to a TGeoTube",
          tube_tag == "revolved-tube" and tube_leaf["type"] == "TGeoTube"
          and abs(tube_leaf["params"]["dz"] - 5.0) < 1e-12
          and abs(tube_leaf["frame"]["origin"][2] - 1.0) < 1e-12,
          f"{tube_tag}, {tube_leaf['type']}, dz {tube_leaf['params']['dz']}, origin "
          f"{tube_leaf['frame']['origin']}")

    # --- rung 2: the prism family, the shapes `_prism_from_rings` writes, read back ---
    # The fixture sews its own explicit faces from ring coordinates it states itself, so the CAD
    # a candidate is measured against is not built by the code that builds the candidate.
    def prism(rings, inner=None):
        stacks = [[[tuple(float(c) for c in q) for q in ring] for ring in rings]]
        if inner is not None:
            stacks.append([[tuple(float(c) for c in q) for q in ring] for ring in inner])
        faces = []
        for stack in stacks:
            nv = len(stack[0])
            for k in range(len(stack) - 1):
                lo, hi = stack[k], stack[k + 1]
                for i in range(nv):
                    j = (i + 1) % nv
                    poly = BRepBuilderAPI_MakePolygon()
                    for q in (lo[i], lo[j], hi[j], hi[i]):
                        poly.Add(gp_Pnt(*q))
                    poly.Close()
                    made = BRepBuilderAPI_MakeFace(poly.Wire())
                    if made.IsDone():
                        faces.append(made.Face())
        for idx in (0, -1):
            poly = BRepBuilderAPI_MakePolygon()
            for q in stacks[0][idx]:
                poly.Add(gp_Pnt(*q))
            poly.Close()
            made = BRepBuilderAPI_MakeFace(poly.Wire())
            if len(stacks) == 2:
                hole = BRepBuilderAPI_MakePolygon()
                for q in stacks[1][idx]:
                    hole.Add(gp_Pnt(*q))
                hole.Close()
                made.Add(topods.Wire(hole.Wire().Reversed()))
            faces.append(made.Face())
        extent = max(abs(c) for stack in stacks for r in stack for q in r for c in q) or 1.0
        sew = BRepBuilderAPI_Sewing(1.0e-7 * extent)
        for face in faces:
            sew.Add(face)
        sew.Perform()
        ms = BRepBuilderAPI_MakeSolid(topods.Shell(sew.SewedShape()))
        ms.Build()
        solid = ms.Solid()
        props = GProp_GProps()
        brepgprop.VolumeProperties(solid, props)
        if props.Mass() < 0.0:
            solid = topods.Solid(solid.Reversed())
        return solid

    def polygon_ring(corners, z):
        return [(x, y, z) for (x, y) in corners]

    def regular_ring(apothem, nedges, z, phi1=0.0, dphi=360.0):
        dseg = math.radians(dphi) / nedges
        radius = apothem / math.cos(dseg / 2.0)
        n = nedges if abs(dphi - 360.0) < 1e-9 else nedges + 1
        return [(radius * math.cos(math.radians(phi1) + k * dseg),
                 radius * math.sin(math.radians(phi1) + k * dseg), z) for k in range(n)]

    def expect_prism(name, solid, want_recogniser, want_type, want_params):
        record = expect(name, solid, want_recogniser)
        if not record["accepted"]:
            return record
        lf = record["candidate"]["leaves"][0]
        worst = 0.0
        if lf["type"] != want_type:
            worst = float("inf")
        else:
            for key, want in want_params.items():
                got = lf["params"][key]
                if isinstance(want, (list, tuple)):
                    worst = (float("inf") if len(got) != len(want)
                             else max([worst] + [abs(a - b) for a, b in zip(got, want)]))
                else:
                    worst = max(worst, abs(got - want))
        check(f"{name} emits a native {want_type} with the source's parameters",
              lf["type"] == want_type and worst < 1.0e-9,
              f"{lf['type']}, worst parameter deviation {worst:.3g}")
        return record

    def trd_rings(dx1, dx2, dy1, dy2, dz):
        return [[(-dx1, -dy1, -dz), (dx1, -dy1, -dz), (dx1, dy1, -dz), (-dx1, dy1, -dz)],
                [(-dx2, -dy2, dz), (dx2, -dy2, dz), (dx2, dy2, dz), (-dx2, dy2, dz)]]

    # TGeoTrd1: the slanted prism behind TPC's 44 "a box face has no opposite partner" declines.
    expect_prism("Trd1 (slanted x faces)", prism(trd_rings(3, 1, 2, 2, 5)), "rung2-trd1",
                 "TGeoTrd1", {"dx1": 3.0, "dx2": 1.0, "dy": 2.0, "dz": 5.0})
    # The taper the other way round, and the barely-slanted TPC_IRB1, whose two half-widths differ
    # by 0.076 cm on 14.2 -- the case a per-class angular criterion would have called a box.
    expect_prism("Trd1 (taper reversed)", prism(trd_rings(1, 3, 2, 2, 4)), "rung2-trd1",
                 "TGeoTrd1", {"dx1": 1.0, "dx2": 3.0, "dy": 2.0, "dz": 4.0})
    expect_prism("Trd1 (TPC_IRB1's 0.5 % slant)",
                 prism(trd_rings(14.205637404580152, 14.281551908396947, 2.06, 2.06, 0.2)),
                 "rung2-trd1", "TGeoTrd1",
                 {"dx1": 14.205637404580152, "dx2": 14.281551908396947, "dy": 2.06, "dz": 0.2})
    # TGeoTrd2: both half-widths vary. The isotropic one is also a legal TGeoXtru, and must not
    # be said as one -- the more specific class wins.
    expect_prism("Trd2 (both half-widths vary)", prism(trd_rings(3, 1, 2, 4, 5)), "rung2-trd2",
                 "TGeoTrd2", {"dx1": 3.0, "dx2": 1.0, "dy1": 2.0, "dy2": 4.0, "dz": 5.0})
    expect_prism("Trd2 (isotropic taper, also a legal Xtru)", prism(trd_rings(3, 1.5, 2, 1, 5)),
                 "rung2-trd2", "TGeoTrd2",
                 {"dx1": 3.0, "dx2": 1.5, "dy1": 2.0, "dy2": 1.0, "dz": 5.0})

    # TGeoArb8: a sheared hexahedron, and TPC_IHSTR's trapezoidal prism stated corner for corner.
    para = prism([[(-2, -2, -3), (2, -2, -3), (2, 2, -3), (-2, 2, -3)],
                  [(-1, -1.5, 3), (3, -1.5, 3), (3, 2.5, 3), (-1, 2.5, 3)]])
    expect_prism("Arb8 (parallelepiped)", para, "rung2-arb8", "TGeoArb8",
                 {"dz": 3.0, "vertices": [-2, -2, 2, -2, 2, 2, -2, 2,
                                          -1, -1.5, 3, -1.5, 3, 2.5, -1, 2.5]})
    ihstr = [(0.0, 0.0), (0.0, 1.08), (2.3, 1.08), (3.38, 0.0)]
    expect_prism("Arb8 (TPC_IHSTR's trapezoidal prism)",
                 prism([polygon_ring(ihstr, -0.6), polygon_ring(ihstr, 0.6)]),
                 "rung2-arb8", "TGeoArb8",
                 {"dz": 0.6, "vertices": [0, 0, 3.38, 0, 2.3, 1.08, 0, 1.08,
                                          0, 0, 3.38, 0, 2.3, 1.08, 0, 1.08]})
    # A hexahedron sheared in x only: neither a Trd (the section centres do not share a line) nor
    # an Xtru (the scale is not isotropic).
    expect("Arb8 (sheared in x only)",
           prism([[(-2, -1, -2), (2, -1, -2), (2, 1, -2), (-2, 1, -2)],
                  [(-2, -1, 2), (4, -1, 2), (4, 1, 2), (-2, 1, 2)]]), "rung2-arb8")
    # ABSO's two TGeoTrap volumes arrive as a planar hexahedron and are said as one. The corners
    # are ROOT's own, from `TGeoTrap(5, 10, 20, 2, 3, 4, 5, 2, 3, 4, 5).GetVertices()`.
    trap_bottom = [(-4.003443140137866, -2.301536896070458),
                   (-4.653488486034171, 1.698463103929542),
                   (3.3465115139658295, 1.698463103929542),
                   (1.996556859862133, -2.301536896070458)]
    trap_top = [(-2.346511513965829, -1.698463103929542),
                (-2.996556859862133, 2.301536896070458),
                (5.003443140137866, 2.301536896070458),
                (3.653488486034171, -1.698463103929542)]
    expect("Arb8 (a TGeoTrap's eight corners)",
           prism([polygon_ring(trap_bottom, -5.0), polygon_ring(trap_top, 5.0)]), "rung2-arb8")

    # TGeoXtru: ITS's 23 Xtru volumes are all right prisms on a general, often non-convex polygon.
    ell_poly = [(0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)]
    expect_prism("Xtru (non-convex L section)",
                 prism([polygon_ring(ell_poly, -2), polygon_ring(ell_poly, 2)]),
                 "rung2-xtru", "TGeoXtru",
                 {"x": [0, 3, 3, 1, 1, 0], "y": [0, 0, 1, 1, 3, 3], "z": [-2, 2],
                  "xoff": [0, 0], "yoff": [0, 0], "scale": [1, 1]})
    rib = [(0, 0), (4.2, 0), (4.2, 0.1), (5.05, 0.1), (9.803, 1.83), (5.9, 1.83), (5.0, 2.73),
           (0, 2.73)]
    expect("Xtru (ITS ConeARibVol0's eight-corner section)",
           prism([polygon_ring(rib, -0.045), polygon_ring(rib, 0.045)]), "rung2-xtru")
    expect("Xtru (a triangular section)",
           prism([polygon_ring([(0, 0), (0.05, 0), (0, 0.074)], -14.5),
                  polygon_ring([(0, 0), (0.05, 0), (0, 0.074)], 14.5)]), "rung2-xtru")
    # Three sections with a per-section offset and an isotropic scale, which is the full
    # vocabulary ROOT gives a TGeoXtru and the only part of it OCC can build exactly.
    scaled_poly = [(0, 0), (2, 0), (2, 1), (1, 2), (0, 2)]
    scaled = prism([[(0.0 + 1.0 * x, 0.0 + 1.0 * y, -3.0) for x, y in scaled_poly],
                    [(0.5 + 1.4 * x, -0.25 + 1.4 * y, 0.0) for x, y in scaled_poly],
                    [(1.0 + 0.6 * x, 0.0 + 0.6 * y, 3.0) for x, y in scaled_poly]])
    expect_prism("Xtru (three sections, offset and scaled)", scaled, "rung2-xtru", "TGeoXtru",
                 {"z": [-3, 0, 3], "xoff": [0, 0.5, 1.0], "yoff": [0, -0.25, 0],
                  "scale": [1.0, 1.4, 0.6]})

    # TGeoPgon: the laterals are planes at the APOTHEM radius, so the corners sit at
    # `r / cos(dseg/2)` and that is what the recogniser inverts. The 48-edge case is the one that
    # looks like a polycone from a distance: it must be said as a polygon, not revolved.
    expect_prism("Pgon (solid hexagonal prism)",
                 prism([regular_ring(3, 6, -5), regular_ring(3, 6, 5)]), "rung2-pgon",
                 "TGeoPgon", {"nedges": 6, "phi1": 0.0, "dphi": 360.0, "z": [-5, 5],
                              "rmin": [0, 0], "rmax": [3, 3]})
    expect_prism("Pgon (tapered eight-edge prism)",
                 prism([regular_ring(3, 8, -5), regular_ring(1.5, 8, 5)]), "rung2-pgon",
                 "TGeoPgon", {"nedges": 8, "phi1": 0.0, "dphi": 360.0, "z": [-5, 5],
                              "rmin": [0, 0], "rmax": [3, 1.5]})
    for nedges in (8, 48):
        hollow = prism([regular_ring(3, nedges, -5), regular_ring(3, nedges, 5)],
                       inner=[regular_ring(1.5, nedges, -5), regular_ring(1.5, nedges, 5)])
        expect_prism(f"Pgon (hollow {nedges}-edge prism)", hollow, "rung2-pgon", "TGeoPgon",
                     {"nedges": nedges, "phi1": 0.0, "dphi": 360.0, "z": [-5, 5],
                      "rmin": [1.5, 1.5], "rmax": [3, 3]})
    # TPC_Strip: 18 edges, a 1 mm wall on an 85 cm radius, 250 cm long -- the annular section no
    # single wire can express and the one case only TGeoPgon carries.
    expect_prism("Pgon (TPC_Strip's thin 18-edge shell)",
                 prism([regular_ring(85.235, 18, -124.8), regular_ring(85.235, 18, 124.8)],
                       inner=[regular_ring(85.225, 18, -124.8), regular_ring(85.225, 18, 124.8)]),
                 "rung2-pgon", "TGeoPgon",
                 {"nedges": 18, "phi1": 0.0, "dphi": 360.0, "z": [-124.8, 124.8],
                  "rmin": [85.225, 85.225], "rmax": [85.235, 85.235]})
    # Three sections, hollow, with the radii stepping: the profile a TGeoPcon would carry, on a
    # polygon.
    expect_prism("Pgon (three hollow sections)",
                 prism([regular_ring(3, 6, -5), regular_ring(3, 6, 0), regular_ring(4, 6, 5)],
                       inner=[regular_ring(1, 6, -5), regular_ring(1, 6, 0),
                              regular_ring(2, 6, 5)]),
                 "rung2-pgon", "TGeoPgon",
                 {"nedges": 6, "phi1": 0.0, "dphi": 360.0, "z": [-5, 0, 5],
                  "rmin": [1, 1, 2], "rmax": [3, 3, 4]})
    # A phi wedge, closing on its own axis, and one that crosses phi = 0 -- where the angular
    # reading has to find the gap rather than sort.
    expect_prism("Pgon (a 90 deg wedge closing on the axis)",
                 prism([regular_ring(4, 3, -2, 10.0, 90.0) + [(0.0, 0.0, -2.0)],
                        regular_ring(4, 3, 2, 10.0, 90.0) + [(0.0, 0.0, 2.0)]]),
                 "rung2-pgon", "TGeoPgon",
                 {"nedges": 3, "phi1": 10.0, "dphi": 90.0, "z": [-2, 2], "rmin": [0, 0],
                  "rmax": [4, 4]})
    expect_prism("Pgon (a wedge across phi = 0)",
                 prism([regular_ring(4, 2, -2, 350.0, 20.0) + [(0.0, 0.0, -2.0)],
                        regular_ring(4, 2, 2, 350.0, 20.0) + [(0.0, 0.0, 2.0)]]),
                 "rung2-pgon", "TGeoPgon",
                 {"nedges": 2, "phi1": 350.0, "dphi": 20.0, "z": [-2, 2], "rmin": [0, 0],
                  "rmax": [4, 4]})

    # A placed prism: the frame machinery on a leaf that has no origin of its own.
    prism_trsf = gp_Trsf()
    prism_trsf.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), 0.7)
    prism_shift = gp_Trsf()
    prism_shift.SetTranslation(gp_Vec(3.0, -4.0, 5.0))
    prism_place = prism_shift.Multiplied(prism_trsf)
    moved_trd = BRepBuilderAPI_Transform(prism(trd_rings(3, 1, 2, 2, 5)), prism_place,
                                         True).Shape()
    moved_trd_record = expect("placed Trd1", moved_trd, "rung2-trd1")
    check("a placed Trd1 travels as one leaf plus a rigid placement",
          moved_trd_record["accepted"]
          and prim.placement_for_candidate(moved_trd_record["candidate"]) is not None,
          "placement present" if moved_trd_record["accepted"] else "not accepted")
    moved_xtru = BRepBuilderAPI_Transform(
        prism([polygon_ring(ell_poly, -2), polygon_ring(ell_poly, 2)]), prism_place,
        True).Shape()
    expect("placed Xtru (non-convex L section)", moved_xtru, "rung2-xtru")

    # --- rung 2 negative controls ---
    # 1. A genuinely twisted TGeoArb8. `_quad_face` writes its four lateral patches as ruled
    #    brepfill surfaces, which are B-splines and not planes, so the solid is declined before
    #    any matcher sees it. That is this rung's scope ruling and it is asserted, not assumed.
    twisted_faces = []
    twist_bottom = [(-2, -2, -2), (-2, 2, -2), (2, 2, -2), (2, -2, -2)]
    twist_top = [(-1.41, -2.73, 2), (-2.73, 1.41, 2), (1.41, 2.73, 2), (2.73, -1.41, 2)]
    for i in range(4):
        j = (i + 1) % 4
        e1 = BRepBuilderAPI_MakeEdge(gp_Pnt(*twist_bottom[i]), gp_Pnt(*twist_bottom[j])).Edge()
        e2 = BRepBuilderAPI_MakeEdge(gp_Pnt(*twist_top[i]), gp_Pnt(*twist_top[j])).Edge()
        twisted_faces.append(brepfill.Face(e1, e2))
    for ring in (twist_bottom, twist_top):
        poly = BRepBuilderAPI_MakePolygon()
        for q in ring:
            poly.Add(gp_Pnt(*q))
        poly.Close()
        twisted_faces.append(BRepBuilderAPI_MakeFace(poly.Wire()).Face())
    sew_twist = BRepBuilderAPI_Sewing(1.0e-6)
    for face in twisted_faces:
        sew_twist.Add(face)
    sew_twist.Perform()
    twist_solid = BRepBuilderAPI_MakeSolid(topods.Shell(sew_twist.SewedShape()))
    twist_solid.Build()
    expect_declined("twisted hexahedron (a ruled TGeoArb8 side)", twist_solid.Solid(),
                    "free-form faces")

    # 2. The near-miss, at two displacements, and they fail on different legs. A rectangular
    #    three-section prism whose middle section is stretched in y only is not a TGeoXtru (that
    #    scale is isotropic), and every lateral of it is still a plane, so the structure alone
    #    cannot say no.
    #      1e-06 cm, ten model tolerances: the recogniser cannot see it and proposes the Xtru;
    #                the symmetric difference refuses it at 1.28e-05 cm^3 against a 6.4e-06 band.
    #      1e-05 cm: the recogniser's own measured gap refuses it, at 8.94e-06 cm on a 6 cm
    #                diagonal against the 1e-06 relative bound.
    for displacement in (1.0e-6, 1.0e-5):
        near = prism([[(-2, -1, -2), (2, -1, -2), (2, 1, -2), (-2, 1, -2)],
                      [(-2, -1 - displacement, 0), (2, -1 - displacement, 0),
                       (2, 1 + displacement, 0), (-2, 1 + displacement, 0)],
                      [(-2, -1, 2), (2, -1, 2), (2, 1, 2), (-2, 1, 2)]])
        expect_declined(f"prism with one section {displacement:g} cm out of similarity", near)

    # 3. A polycone must not be taken by a prism template, and a many-edged polygon must not be
    #    taken by the revolved one. The pair is the point: the two matchers key on the same
    #    structure from opposite sides, and only `_match_box` declining first separates them.
    check("a polycone reaches the revolved matcher, not the prism one",
          process_solid(stepped, "pcon-vs-prism")["recogniser"] == "revolved-pcon",
          f"{process_solid(stepped, 'pcon-vs-prism')['recogniser']}")

    # --- the instrument that scores a prism candidate must be able to say "no" ---
    exact_ring = [(-2.0, -1.0, -2.0), (2.0, -1.0, -2.0), (2.0, 1.0, -2.0), (-2.0, 1.0, -2.0)]
    nudged = [(x + (1.0e-6 if i == 0 else 0.0), y, z)
              for i, (x, y, z) in enumerate(exact_ring)]
    check("the point-set gap is zero on the point set itself",
          recognise._point_set_gap(exact_ring, exact_ring) == 0.0,
          f"gap {recognise._point_set_gap(exact_ring, exact_ring):.3g} cm")
    nudged_gap = recognise._point_set_gap(exact_ring, nudged)
    check("the point-set gap reports a corner displaced by ten model tolerances",
          abs(nudged_gap - 1.0e-6) < 1.0e-15, f"gap {nudged_gap:.3g} cm, expected 1e-06 cm")
    # ... and it must report a hexahedron read out in the wrong corner order, which has exactly
    # the same eight corners. This is why `prism_samples` carries the edge midpoints.
    good_arb8 = prim.leaf("TGeoArb8", {"dz": 3.0,
                                       "vertices": [-2, -2, 2, -2, 2, 2, -2, 2,
                                                    -1, -1.5, 3, -1.5, 3, 2.5, -1, 2.5]},
                          prim.identity_frame())
    swapped = list(good_arb8["params"]["vertices"])
    swapped[2:4], swapped[4:6] = swapped[4:6], swapped[2:4]
    bad_arb8 = prim.leaf("TGeoArb8", {"dz": 3.0, "vertices": swapped}, prim.identity_frame())
    corner_only_gap = recognise._point_set_gap(
        [tuple(q) for q in prim.prism_samples(good_arb8)[0::3]],
        [tuple(q) for q in prim.prism_samples(bad_arb8)[0::3]])
    order_gap = recognise._point_set_gap(prim.prism_samples(good_arb8),
                                         prim.prism_samples(bad_arb8))
    check("the edge midpoints are what catch a hexahedron read in the wrong corner order",
          corner_only_gap == 0.0 and order_gap > 0.1,
          f"corners alone {corner_only_gap:.3g} cm, corners and edge midpoints "
          f"{order_gap:.3g} cm")

    # --- the description must refuse an illegal prism before either builder sees it ---
    for name, kind, params in (
            ("a TGeoXtru whose z runs backwards", "TGeoXtru",
             {"x": [0, 1, 0], "y": [0, 0, 1], "z": [1.0, 0.0], "xoff": [0, 0], "yoff": [0, 0],
              "scale": [1, 1]}),
            ("a TGeoXtru with a repeated corner", "TGeoXtru",
             {"x": [0, 1, 1], "y": [0, 0, 0], "z": [0.0, 1.0], "xoff": [0, 0], "yoff": [0, 0],
              "scale": [1, 1]}),
            ("a TGeoXtru with two corners", "TGeoXtru",
             {"x": [0, 1], "y": [0, 0], "z": [0.0, 1.0], "xoff": [0, 0], "yoff": [0, 0],
              "scale": [1, 1]}),
            ("a TGeoXtru with a zero scale", "TGeoXtru",
             {"x": [0, 1, 0], "y": [0, 0, 1], "z": [0.0, 1.0], "xoff": [0, 0], "yoff": [0, 0],
              "scale": [1, 0]}),
            ("a TGeoArb8 with fifteen coordinates", "TGeoArb8",
             {"dz": 1.0, "vertices": [0.0] * 15}),
            ("a TGeoArb8 with a collapsed face", "TGeoArb8",
             {"dz": 1.0, "vertices": [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1]}),
            ("a TGeoTrd1 with both half-widths zero", "TGeoTrd1",
             {"dx1": 0.0, "dx2": 0.0, "dy": 1.0, "dz": 1.0}),
            ("a TGeoTrd2 with a negative half-width", "TGeoTrd2",
             {"dx1": -1.0, "dx2": 1.0, "dy1": 1.0, "dy2": 1.0, "dz": 1.0}),
            ("a TGeoPgon with no edges", "TGeoPgon",
             {"phi1": 0.0, "dphi": 360.0, "nedges": 0, "z": [0.0, 1.0], "rmin": [0.0, 0.0],
              "rmax": [1.0, 1.0]})):
        try:
            prim.leaf(kind, params, prim.identity_frame())
            refused = False
        except ValueError:
            refused = True
        check(f"{name} is refused", refused)

    # The two independent array lengths a TGeoXtru needs -- a polygon count and a section count --
    # go through the same declaration Pcon's single length does, so a mismatch inside one group
    # is still refused.
    try:
        prim.leaf("TGeoXtru", {"x": [0, 1, 0], "y": [0, 0], "z": [0.0, 1.0], "xoff": [0, 0],
                               "yoff": [0, 0], "scale": [1, 1]}, prim.identity_frame())
        refused = False
    except ValueError:
        refused = True
    check("a TGeoXtru whose x and y differ in length is refused", refused)
    xtru_two_lengths = prim.leaf(
        "TGeoXtru", {"x": [0, 2, 2, 0], "y": [0, 0, 1, 1], "z": [-1.0, 0.0, 1.0],
                     "xoff": [0, 0, 0], "yoff": [0, 0, 0], "scale": [1, 1, 1]},
        prim.identity_frame())
    check("a TGeoXtru carries four corners and three sections in one description",
          len(xtru_two_lengths["params"]["x"]) == 4 and len(xtru_two_lengths["params"]["z"]) == 3,
          f"{len(xtru_two_lengths['params']['x'])} corners, "
          f"{len(xtru_two_lengths['params']['z'])} sections")

    # --- the floor for rung 2's own emissions ---
    check("every prism-family candidate is byte-identical to its recorded digest",
          all(seen_digests.get(name) == digest for name, digest
              in _PRISM_CANDIDATE_DIGESTS.items()),
          "; ".join(f"{name}: {seen_digests.get(name)} != {digest}" for name, digest
                    in _PRISM_CANDIDATE_DIGESTS.items()
                    if seen_digests.get(name) != digest)
          or f"{len(_PRISM_CANDIDATE_DIGESTS)} candidates unchanged")

    # --- the floor: nothing that converted before this matcher existed converts differently ---
    check("every whole-part candidate is byte-identical to before the revolved matcher",
          all(seen_digests.get(name) == digest for name, digest
              in _CANDIDATE_DIGESTS_BEFORE_THE_REVOLVED_MATCHER.items()),
          "; ".join(f"{name}: {seen_digests.get(name)} != {digest}" for name, digest
                    in _CANDIDATE_DIGESTS_BEFORE_THE_REVOLVED_MATCHER.items()
                    if seen_digests.get(name) != digest)
          or f"{len(_CANDIDATE_DIGESTS_BEFORE_THE_REVOLVED_MATCHER)} candidates unchanged")

    # --- the ROOT half: the emitted TGeoShape must answer like the closed form ---
    if with_root:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        from array import array
        import random
        shape, placement = prim.build_root(moved_record["candidate"], "probe_moved")
        # `Stream_N_PlacedPrimitives.md`: a placed primitive is now the bare primitive plus a
        # transform, not the self-union composite. Assert the class, because everything below
        # would also pass on a composite and the class is the whole point of the change.
        check("a rotated, translated tube emits a bare TGeoTube, not a TGeoCompositeShape",
              shape.ClassName() == "TGeoTube" and placement is not None,
              f"{shape.ClassName()}, placement {'present' if placement else 'absent'}")
        # closed form for the placed tube: transform the point into the tube's frame and test
        # 1 <= r <= 2, |z| <= 5. The frame is stated here from the *transform that built the
        # OCCT solid*, not from the description, so this is an independent statement.
        frame = moved_record["candidate"]["leaves"][0]["frame"]
        bad = 0
        random.seed(11)
        for _ in range(20000):
            p = (random.uniform(-2, 8), random.uniform(-9, 1), random.uniform(0, 10))
            rel = prim._sub(p, tuple(frame["origin"]))
            zc = prim._dot(rel, tuple(frame["z"]))
            rc = math.sqrt(max(prim._dot(rel, rel) - zc * zc, 0.0))
            want = (1.0 <= rc <= 2.0) and abs(zc) <= 5.0
            got = bool(shape.Contains(array("d", list(prim.placement_to_local(placement, p)))))
            if want != got and min(abs(rc - 1.0), abs(rc - 2.0), abs(abs(zc) - 5.0)) > 1e-9:
                bad += 1
        check("the emitted placed tube answers Contains like the closed form",
              bad == 0, f"{bad} disagreement(s) over 20000 points")
        # An analytic Capacity() is what the composite cost us. pi (rmax^2 - rmin^2) 2 dz for
        # rmin=1, rmax=2, dz=5, and it is invariant under the placement.
        want_capacity = math.pi * (2.0 ** 2 - 1.0 ** 2) * 10.0
        rel_capacity = abs(shape.Capacity() - want_capacity) / want_capacity
        check("the placed tube's Capacity() is analytic", rel_capacity < 1.0e-14,
              f"{shape.Capacity():.12f} vs {want_capacity:.12f}, rel {rel_capacity:.2e}")
        # negative control on that check itself
        wrong, wrong_pl = prim.build_root(prim.candidate("primitive", [prim.leaf(
            "TGeoTube", {"rmin": 1.0, "rmax": 2.05, "dz": 5.0}, frame)], "probe"), "probe_wrong")
        bad_wrong = 0
        random.seed(11)
        for _ in range(20000):
            p = (random.uniform(-2, 8), random.uniform(-9, 1), random.uniform(0, 10))
            rel = prim._sub(p, tuple(frame["origin"]))
            zc = prim._dot(rel, tuple(frame["z"]))
            rc = math.sqrt(max(prim._dot(rel, rel) - zc * zc, 0.0))
            want = (1.0 <= rc <= 2.0) and abs(zc) <= 5.0
            if want != bool(wrong.Contains(array("d", list(prim.placement_to_local(wrong_pl, p))))):
                bad_wrong += 1
        check("the same check does report a wrong radius", bad_wrong > 0,
              f"{bad_wrong} disagreement(s) with rmax 2.05")
        # ... and the placement itself must be load-bearing: transposing its rotation has to move
        # the count. Without this the check above would pass on a shape whose placement is
        # ignored, which is precisely the mistake a composition-order bug makes.
        transposed = [[placement[r][c] for r in range(3)] + [placement[c][3]] for c in range(3)]
        bad_transposed = 0
        random.seed(11)
        for _ in range(20000):
            p = (random.uniform(-2, 8), random.uniform(-9, 1), random.uniform(0, 10))
            rel = prim._sub(p, tuple(frame["origin"]))
            zc = prim._dot(rel, tuple(frame["z"]))
            rc = math.sqrt(max(prim._dot(rel, rel) - zc * zc, 0.0))
            want = (1.0 <= rc <= 2.0) and abs(zc) <= 5.0
            got = bool(shape.Contains(array("d", list(prim.placement_to_local(transposed, p)))))
            if want != got:
                bad_transposed += 1
        check("a transposed placement rotation does move the count", bad_transposed > 0,
              f"{bad_transposed} disagreement(s) with R^T")
        # the round trip through the artefact: placement written, placement read back
        placed_target = Path("/tmp/csg_selftest_placed.root")
        write_shape_root(moved_record["candidate"], placed_target)
        fp = ROOT.TFile.Open(str(placed_target))
        back_shape = fp.Get("shape")
        back_matrix = fp.Get("placement")
        back_placement = prim.placement_from_root_matrix(back_matrix) if back_matrix else None
        worst_pl = (max(abs(back_placement[r][c] - placement[r][c])
                        for r in range(3) for c in range(4))
                    if back_placement is not None else float("inf"))
        check("shape_<part>.root round-trips the placement under the key \"placement\"",
              back_shape is not None and back_shape.ClassName() == "TGeoTube"
              and worst_pl < 1.0e-15,
              f"read {back_shape.ClassName() if back_shape else 'nothing'}, worst placement "
              f"element deviation {worst_pl:.3g}")
        fp.Close()
        # the two-leaf union must round-trip through a file and keep its class
        target = Path("/tmp/csg_selftest_shape.root")
        written = write_shape_root(ram_record["candidate"], target)
        f = ROOT.TFile.Open(str(target))
        back = f.Get("shape")
        check("a two-leaf union round-trips through shape_<part>.root",
              back and back.InheritsFrom("TGeoShape"),
              f"wrote {written.ClassName()}, read {back.ClassName() if back else 'nothing'}")
        f.Close()
        dev = crosscheck_bbox(ram_record["candidate"])
        check("the OCCT and ROOT realisations agree on the bounding box", dev < 1.0e-9,
              f"max deviation {dev:.3g} cm")

        # An axis-aligned box must come out as a *bare* TGeoBBox carrying its own origin, not as
        # the self-union a rotated primitive needs. This is the only path on which ROOT reports
        # an analytic Capacity(), and it is the common case in mechanical CAD (the census counts
        # 62560 placed six-plane boxes in oTOF), so it is worth asserting rather than assuming.
        box_record = process_solid(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0).Shape(),
                                   "box-emission")
        box_shape, box_placement = prim.build_root(box_record["candidate"], "boxprobe")
        origin = [box_shape.GetOrigin()[i] for i in range(3)]
        check("an axis-aligned box emits a bare TGeoBBox with its own origin",
              box_shape.ClassName() == "TGeoBBox" and box_placement is None
              and max(abs(origin[0] - 1.0), abs(origin[1] - 1.5), abs(origin[2] - 2.0)) < 1e-12
              and abs(box_shape.Capacity() - 24.0) < 1e-12,
              f"{box_shape.ClassName()}, origin {origin}, capacity {box_shape.Capacity():.6f}, "
              f"placement {'present' if box_placement else 'absent'}")

        # A genuine multi-leaf boolean is out of scope for this change and must stay a composite,
        # unplaced. Asserted so that "the composite is gone" can never quietly become true of the
        # cases that legitimately need one.
        ram_shape, ram_placement = prim.build_root(ram_record["candidate"], "ramprobe")
        check("a genuine two-leaf union is still an unplaced TGeoCompositeShape",
              ram_shape.ClassName() == "TGeoCompositeShape" and ram_placement is None,
              f"{ram_shape.ClassName()}, placement "
              f"{'present' if ram_placement else 'absent'}")

        # --- the ROOT half of the revolved matcher ---
        stepped_record = process_solid(stepped, "pcon-emission")
        pcon_shape, pcon_placement = prim.build_root(stepped_record["candidate"], "pconprobe")
        # 100 pi: pi (3^2 - 1^2) 5 below z = 0 and pi (4^2 - 2^2) 5 above it. Stated here from the
        # fixture's own numbers, so it is independent of both builders.
        want_capacity = math.pi * ((3.0 ** 2 - 1.0 ** 2) * 5.0 + (4.0 ** 2 - 2.0 ** 2) * 5.0)
        rel_capacity = abs(pcon_shape.Capacity() - want_capacity) / want_capacity
        check("an axis-aligned polycone emits a bare TGeoPcon with an analytic Capacity()",
              pcon_shape.ClassName() == "TGeoPcon" and pcon_placement is None
              and rel_capacity < 1.0e-14,
              f"{pcon_shape.ClassName()}, capacity {pcon_shape.Capacity():.9f} vs "
              f"{want_capacity:.9f} (rel {rel_capacity:.2e}), placement "
              f"{'present' if pcon_placement else 'absent'}")

        placed_pcon_shape, placed_pcon_placement = prim.build_root(
            moved_pcon_record["candidate"], "placedpconprobe")
        check("a placed polycone is a bare TGeoPcon plus a placement",
              placed_pcon_shape.ClassName() == "TGeoPcon" and placed_pcon_placement is not None,
              f"{placed_pcon_shape.ClassName()}, placement "
              f"{'present' if placed_pcon_placement else 'absent'}")
        # The closed form is stated through the inverse of the transform that *built* the OCCT
        # solid, never through the description's own frame, so a placement that is wrong in the
        # same way in both builders is still caught.
        pcon_inverse = pcon_place.Inverted()
        bad_pcon = 0
        scored_pcon = 0
        random.seed(23)
        for _ in range(20000):
            p3 = (random.uniform(-3, 9), random.uniform(-10, 2), random.uniform(-1, 11))
            probe = gp_Pnt(*p3)
            probe.Transform(pcon_inverse)
            zc, rc = probe.Z(), math.hypot(probe.X(), probe.Y())
            if min(abs(zc + 5.0), abs(zc), abs(zc - 5.0), abs(rc - 1.0), abs(rc - 2.0),
                   abs(rc - 3.0), abs(rc - 4.0)) < 1.0e-6:
                continue
            scored_pcon += 1
            want = (1.0 <= rc <= 3.0) if -5.0 <= zc <= 0.0 else (
                (2.0 <= rc <= 4.0) if 0.0 < zc <= 5.0 else False)
            got = bool(placed_pcon_shape.Contains(
                array("d", list(prim.placement_to_local(placed_pcon_placement, p3)))))
            if want != got:
                bad_pcon += 1
        check("the emitted placed polycone answers Contains like the closed form",
              bad_pcon == 0, f"{bad_pcon} disagreement(s) over {scored_pcon} points")
        cc = crosscheck_contains(moved_pcon_record["candidate"], moved_pcon)
        check("the ROOT polycone and the CAD solid agree on Contains",
              cc["disagreements"] == 0,
              f"{cc['disagreements']} disagreement(s) over {cc['points']} points")

        pcon_target = Path("/tmp/csg_selftest_pcon.root")
        write_shape_root(stepped_record["candidate"], pcon_target)
        fpcon = ROOT.TFile.Open(str(pcon_target))
        back_pcon = fpcon.Get("shape")
        sections_ok = (back_pcon is not None and back_pcon.ClassName() == "TGeoPcon"
                       and back_pcon.GetNz() == 4
                       and max(abs(back_pcon.GetZ(i) - step_z[i]) for i in range(4)) < 1e-15
                       and max(abs(back_pcon.GetRmin(i) - step_rmin[i]) for i in range(4)) < 1e-15
                       and max(abs(back_pcon.GetRmax(i) - step_rmax[i]) for i in range(4)) < 1e-15)
        check("shape_<part>.root round-trips a TGeoPcon with all its sections", sections_ok,
              f"read {back_pcon.ClassName() if back_pcon else 'nothing'}, nz "
              f"{back_pcon.GetNz() if back_pcon else 0}")
        fpcon.Close()

        # --- the ROOT half of the prism family ---
        # Each class must come out as itself, not as a composite and not as the more general
        # class next door, and each must report an analytic Capacity() -- which is the whole
        # reason for emitting the specific class rather than a general one.
        for name, solid, want_class, want_capacity in (
                ("Trd1", prism(trd_rings(3, 1, 2, 2, 5)), "TGeoTrd1",
                 4.0 * 2.0 * (3.0 + 1.0) * 5.0),
                ("Trd2", prism(trd_rings(3, 1, 2, 4, 5)), "TGeoTrd2", None),
                ("Arb8", para, "TGeoArb8", None),
                ("Xtru", prism([polygon_ring(ell_poly, -2), polygon_ring(ell_poly, 2)]),
                 "TGeoXtru", 5.0 * 4.0),
                ("Pgon", prism([regular_ring(3, 6, -5), regular_ring(3, 6, 5)]), "TGeoPgon",
                 6.0 * 9.0 * math.tan(math.pi / 6.0) * 10.0)):
            record = process_solid(solid, f"{name}-emission")
            if not record["accepted"]:
                check(f"an axis-aligned {want_class} emits a bare {want_class}", False,
                      f"not accepted: {record['reason']}")
                continue
            shape, placed = prim.build_root(record["candidate"], f"{name}probe")
            ok = shape.ClassName() == want_class and placed is None
            detail = (f"{shape.ClassName()}, capacity {shape.Capacity():.9f}, placement "
                      f"{'present' if placed else 'absent'}")
            if want_capacity is not None:
                rel = abs(shape.Capacity() - want_capacity) / want_capacity
                ok = ok and rel < 1.0e-12
                detail += f", closed form {want_capacity:.9f} (rel {rel:.2e})"
            check(f"an axis-aligned {want_class} emits a bare {want_class} with an analytic "
                  "Capacity()", ok, detail)

        # A placed Trd1: the closed form is stated through the inverse of the transform that
        # *built* the OCCT solid, never through the description's own frame, so a placement that
        # is wrong in the same way in both builders is still caught.
        trd_shape, trd_placement = prim.build_root(moved_trd_record["candidate"], "movedtrdprobe")
        check("a placed Trd1 is a bare TGeoTrd1 plus a placement",
              trd_shape.ClassName() == "TGeoTrd1" and trd_placement is not None,
              f"{trd_shape.ClassName()}, placement "
              f"{'present' if trd_placement else 'absent'}")
        trd_inverse = prism_place.Inverted()
        bad_trd = 0
        scored_trd = 0
        random.seed(37)
        for _ in range(20000):
            p3 = (random.uniform(-3, 9), random.uniform(-10, 2), random.uniform(-2, 12))
            probe = gp_Pnt(*p3)
            probe.Transform(trd_inverse)
            xc, yc, zc = probe.X(), probe.Y(), probe.Z()
            half = 2.0 - 0.2 * zc                      # dx1 = 3, dx2 = 1, dz = 5
            if min(abs(abs(zc) - 5.0), abs(abs(yc) - 2.0), abs(abs(xc) - half)) < 1.0e-6:
                continue
            scored_trd += 1
            want = abs(zc) <= 5.0 and abs(yc) <= 2.0 and abs(xc) <= half
            got = bool(trd_shape.Contains(
                array("d", list(prim.placement_to_local(trd_placement, p3)))))
            if want != got:
                bad_trd += 1
        check("the emitted placed Trd1 answers Contains like the closed form",
              bad_trd == 0, f"{bad_trd} disagreement(s) over {scored_trd} points")
        cc_prism = crosscheck_contains(moved_trd_record["candidate"], moved_trd)
        check("the ROOT Trd1 and the CAD solid agree on Contains",
              cc_prism["disagreements"] == 0,
              f"{cc_prism['disagreements']} disagreement(s) over {cc_prism['points']} points")

        # The artefact must carry a TGeoXtru's polygon and its sections, which is the one class
        # here whose parameters do not fit in a constructor call.
        xtru_record = process_solid(scaled, "xtru-emission")
        xtru_target = Path("/tmp/csg_selftest_xtru.root")
        write_shape_root(xtru_record["candidate"], xtru_target)
        fxtru = ROOT.TFile.Open(str(xtru_target))
        back_xtru = fxtru.Get("shape")
        xtru_ok = (back_xtru is not None and back_xtru.ClassName() == "TGeoXtru"
                   and back_xtru.GetNvert() == 5 and back_xtru.GetNz() == 3
                   and max(abs(back_xtru.GetZ(k) - z) for k, z in enumerate((-3.0, 0.0, 3.0)))
                   < 1e-12
                   and max(abs(back_xtru.GetScale(k) - v)
                           for k, v in enumerate((1.0, 1.4, 0.6))) < 1e-12)
        check("shape_<part>.root round-trips a TGeoXtru with its polygon and its sections",
              xtru_ok, f"read {back_xtru.ClassName() if back_xtru else 'nothing'}, "
                       f"nvert {back_xtru.GetNvert() if back_xtru else 0}, "
                       f"nz {back_xtru.GetNz() if back_xtru else 0}")
        fxtru.Close()

    n_ok = sum(1 for _n, ok, _d in checks if ok)
    if verbose:
        print(f"  {n_ok}/{len(checks)} recognise/emit self-checks passed")
    return n_ok, len(checks)


# ------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, help="a gate workdir's db/ directory (walks brep_*.brep)")
    ap.add_argument("--brep", type=Path, help="a single .brep file, in cm")
    ap.add_argument("--from-json", type=Path, dest="from_json",
                    help="write shape_*.root for every accepted csg_*.json in this folder "
                         "(needs PyROOT only; nothing is re-recognised)")
    ap.add_argument("--report", type=Path, help="write the per-part record as JSON")
    ap.add_argument("--no-root", action="store_true",
                    help="recognise and accept, but do not write shape_*.root (no PyROOT needed)")
    ap.add_argument("--band-factor", type=float, default=1.0,
                    help="multiplier on the acceptance band (model tolerance x area); "
                         "default %(default)s")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--no-self-test", action="store_true",
                    help="skip the self-test that otherwise runs before any emission")
    args = ap.parse_args()

    if args.self_test:
        ok_a, n_a = accept.self_test()
        ok_e, n_e = self_test(with_root=not args.no_root)
        print(f"\n{ok_a + ok_e}/{n_a + n_e} self-checks passed")
        return 0 if (ok_a == n_a and ok_e == n_e) else 1

    if args.from_json:
        from_json(args.from_json)
        return 0

    if not args.db and not args.brep:
        ap.error("give --db, --brep, --from-json or --self-test")

    if not args.no_self_test:
        ok_a, n_a = accept.self_test(verbose=False)
        ok_e, n_e = self_test(verbose=False, with_root=not args.no_root)
        if ok_a != n_a or ok_e != n_e:
            raise SystemExit(f"self-test failed ({ok_a}/{n_a} acceptance, {ok_e}/{n_e} "
                             "recognise/emit); refusing to emit")
        print(f"[self-test] {ok_a + ok_e}/{n_a + n_e} checks passed")

    if args.brep:
        solid, _n = load_brep(args.brep)
        suffix = args.brep.name[len("brep_"):-len(".brep")]
        record = process_solid(solid, suffix, band_factor=args.band_factor)
        if record["accepted"] and not args.no_root:
            target = args.brep.parent / f"shape_{suffix}.root"
            write_shape_root(record["candidate"], target)
            record["shape"] = str(target)
            record["bboxRootVsOcctCm"] = crosscheck_bbox(record["candidate"])
            record["containsCrosscheck"] = crosscheck_contains(record["candidate"], solid)
        _print_record(record)
        records = [record]
    else:
        records = run_db(args.db, write_root=not args.no_root, band_factor=args.band_factor)
    summarise(records)
    if args.report:
        args.report.write_text(json.dumps(records, indent=1))
        print(f"Wrote {args.report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
