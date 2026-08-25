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

# The shape-tolerance helper lives in `csg/accept.py`, because `csg/recognise.py` needs it too
# and must not import this module. Re-exported here under its long-standing name.
model_tolerance_cm = accept.model_tolerance_cm


# What `process_solid` tries, in order, when the cascade's own candidate is REJECTED by the
# acceptance test rather than declined by a matcher. Each entry answers about the whole solid and
# skips every earlier matcher, so the order here is the order of increasing generality.
_RETRIES = (("a revolved profile", recognise.recognise_revolved),
            ("a single cell", recognise.recognise_single_cell),
            ("a union of cells", recognise.recognise_union_of_cells),
            # Last, and only after the union of cells has been tried and refused, for the same
            # reason the cascade orders them that way: no part that converts today changes
            # representation.
            ("flat cells", recognise.recognise_flat_cells))


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

    notes = []
    for label, propose in _RETRIES:
        alternative, alt_declined = propose(solid)
        if alternative is None:
            notes.append(f"as {label}: {alt_declined}")
            continue
        if alternative["recogniser"] == cand["recogniser"]:
            # This matcher is what produced the candidate that was just refused; retrying it
            # would refuse it again.
            notes.append(f"as {label}: the same proposal that was just rejected")
            continue
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
        notes.append(f"retried as {alternative['recogniser']} "
                     f"({prim.describe(alternative)}): {alt_why_not}")
    record["reason"] = "; ".join([why_not] + notes)
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
        if payload["candidate"].get("op") == "flatCells":
            # The macro loads `flatcsg_<part>.bin`, not the `.root` file, so completing a
            # deferred flat part means writing the sidecar too. `csg/hook.py` already writes it
            # in the same folder without needing ROOT; rewriting it here is idempotent and keeps
            # this entry point usable on a folder it did not produce.
            from csg import flat as flat_writer
            blocks, cells = prim.flat_sidecar_records(payload["candidate"])
            flat_writer.write_sidecar(path.parent / f"flatcsg_{suffix}.bin", blocks, cells)
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


# The flat-CSG programme's R1 and R2 emissions, frozen the same way: the torus carrier and the
# elliptic cylinder. `torus shell (a bellows ply)` and `half a bellows ply` are the two shapes
# PIPE's sixteen plies reduce to, so what they emit is worth being unable to change by accident.
_TORUS_ELTU_CANDIDATE_DIGESTS = {
    "solid torus":
        "d72f862e7d29c9775d0ba20669f6495d9b2670fe863b80cb4a855ec0339b737e",
    "torus shell (a bellows ply)":
        "9bf469b23b0f8de4db091924473463e893c29662cf28703259e903934648ff78",
    "hollow torus wedge":
        "ce6dbc08dbb1df8613edf26c12f1f9ae6a67b7ba4b484d404ea218f981c63355",
    "half a bellows ply":
        "4a4a8485968902fcf1d1c99672d46090d140bf2c4118da4b6aa500b0a8e33471",
    "elliptic cylinder, a > b":
        "ca335cfe52684bd5d7e4680d0ebb50783cef954de2c7fd3e884d4a757f32a36a",
    "elliptic cylinder, a < b":
        "dfae936208a9a94b23286472b24be3c0b91ef9021b6f0ab83069db2eb50fbab1",
    "elliptic cylinder with equal semi-axes":
        "cbe32ffd7120066ee6b01803e29966ddc46294bdc8b7a5a9dfa2239a171bcd48",
}

# The single cell's candidates, frozen the same way. `tube_window` and `cyl_inter_cyl` are the
# two parts `Stream_AA_FlatCSG.md` §5 step 2 exists to retire, so what they emit is worth being
# unable to change by accident.
_CELL_CANDIDATE_DIGESTS = {
    "Steinmetz solid (two cylinders intersected)":
        "9f87b85a40ae93b2fe3e94c582974225d0954fa4c39337358fdd9bf7fee8db57",
    "tube with a transverse window":
        "e85097015a53ee1ff72b99a842a1f97a319a44e27e407beaa11c974f83c2660c",
    "cylinder cut by an oblique plane":
        "77388f8d14e7e706e3fce4efb3b7ef2b7347c185278fc97c209c679d43f5b6ec",
    "cube with an axial through-hole":
        "506948171bc7ae669699ccc49013126000b4cb6f4b5327f00293959de4b7e56f",
    "cylinder with a milled flat":
        "21b1e7406b7f17dca37349d160742135c17af3b6865673225fdb07d61459018b",
}

# Rung 4's own emissions: the two-level DNF. What is frozen is the whole `unionOfCells`
# description -- the cells, their order, and every leaf in them -- because a decomposition that
# silently reorders or re-folds its cells tomorrow is exactly the drift a digest exists to catch.
_UNION_OF_CELLS_CANDIDATE_DIGESTS = {
    "a cylinder with a hexagonal collar":
        "056503ec6d01870b075b59bef50c0702e251a543da8869141fbb3c33d40e37f8",
    "two rods sharing no edge":
        "450a9f7857fd1114ba1e8115abcd02fec08121ee5c349a6a7f976d201a29c35a",
    "a torus with a cylinder through it":
        "ad897f27f24c758cc0f853025fbb339956f7270ff92843bba6172464a5de706e",
    "three disjoint boxes":
        "d5db83ca1fd3baf02430ad34d29937a673415ebd7221050710de545fccaafb65",
}


# Rung 3's own emissions: the whole-part fixtures whose every carrier arrives Tier-0 canonicalised
# from a stored B-spline. What is worth freezing here is not just that they convert but that they
# convert to the SAME description a natively-analytic twin does -- so these digests are asserted
# against the corresponding rows above as well as against themselves.
_TIER0_CANDIDATE_DIGESTS = {
    "NURBS-encoded box":
        "3ce4f94e64c58317f51cdeee12835f804335d3eafdf0c72125c80f0fd180a806",
    "NURBS-encoded solid cylinder":
        "8fb22fddb5e8480b8f7a3fae32754a94dea1b826dec21358856885071844673e",
    "NURBS-encoded tube segment":
        "3a18fceec54f3665cef60fa94575d0441e82c027a60881da74fd316dd3a99f95",
    "NURBS-encoded cone":
        "b1b1168b57f1f4e1315b695170f70ddd5cfed936482293f8603452737274818a",
    "NURBS-encoded sphere":
        "6e9a7e981bb74bcd688d129a9bbe556618ca52e8482e5a18e9b104058cd1ca8d",
    "NURBS-encoded solid torus":
        "510e17434832263effc284e3d9689aba4ba9436cef23d8d083ca30bedc30dd60",
    "NURBS-encoded hollow torus wedge":
        "5f0a993c08cb82afd05e0b35b3b7e3338b723b45d36f8d57b925a42457d6ffed",
    "NURBS-encoded cube with an axial through-hole":
        "04193f9ddb417d5d5cd01e42c67bce2c8ce55cbcfad290fc950083d5d9046e48",
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
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                         BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeSolid,
                                         BRepBuilderAPI_MakeWire, BRepBuilderAPI_Sewing,
                                         BRepBuilderAPI_Transform)
    from OCC.Core.BRepFill import brepfill
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                                      BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism,
                                      BRepPrimAPI_MakeRevol, BRepPrimAPI_MakeSphere,
                                      BRepPrimAPI_MakeTorus)
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.TopoDS import topods
    from OCC.Core.GeomAPI import GeomAPI_Interpolate
    from OCC.Core.TColgp import TColgp_HArray1OfPnt
    from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Elips, gp_Pnt, gp_Trsf, gp_Vec

    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))
        if verbose:
            print(f"  [{'ok ' if condition else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))

    seen_digests = {}
    seen_recognisers = {}

    def expect(name, solid, want_recogniser, want_leaves=1):
        record = process_solid(solid, name)
        seen_recognisers[name] = record["recogniser"] if record["accepted"] else None
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

    def expect_single_cell_declined(name, solid, needle):
        """The single-cell verdict, asserted at the level the control was written about.

        Rung 4 converts genuinely multi-cell bodies, so several controls that used to end at
        "not converted as CSG" now end one rung further down. What they were written to protect
        is the ONE-CELL read's verdict -- that a body with a trusted concave edge is not one
        cell, and that it says so -- and that is asserted here, unchanged, against the matcher
        that makes it. What the part then converts to is a separate assertion beside this one.
        """
        _cand, why = recognise.recognise_single_cell(solid)
        ok = _cand is None and needle in (why or "")
        check(f"{name} is refused by the one-cell read", ok, f"reason: {why}")
        return why

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
    # 3. a torus. It used to end the cascade at `_face_records`; the torus carrier landed with
    #    the flat-CSG programme's R1, so the expectation moved and the shape did not.
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
    expect("torus", BRepPrimAPI_MakeTorus(5.0, 1.0).Shape(), "tier1-torus")
    # 4. a cylinder with a flat milled off it: one cluster, but a plane that is neither cap nor
    #    wedge. The recogniser must decline rather than emit the round tube.
    flatted = BRepAlgoAPI_Cut(cyl, BRepPrimAPI_MakeBox(
        gp_Pnt(1.5, -3, -6), 3.0, 6.0, 12.0).Shape()).Shape()
    # A milled flat is a cylinder, its two caps and one plane: four halfspaces, no concave edge,
    # one cell. The single-cell emitter converts it exactly, so the expectation moved and the
    # shape did not -- the same flip the blind bore made when the revolved matcher landed.
    flat_record = expect("cylinder with a milled flat", flatted, "cell-intersection",
                         want_leaves=2)

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
    # It is not a polycone and not a prism, and the whole-part matchers must say so. Since
    # rung 4 it *is* converted -- as the two cells it genuinely is, asserted below beside the
    # rest of the decomposition controls -- so what is checked here is the class it must never
    # be given, which is what this ladder has always been about.
    hybrid_single = recognise.recognise_single_cell(hybrid)[1]
    check("a cylinder with a coaxial hexagonal section is no whole-part primitive",
          recognise.recognise(hybrid)[0]["recogniser"] == "cells-union"
          and "neither a cap nor a wedge" in (recognise.recognise_revolved(hybrid)[1] or ""),
          f"one-cell read: {(hybrid_single or '')[:90]}")
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

    # --- rung 3: the single cell, the three fixtures Stream_AA §5 step 2 names ---
    # `BRepPrimAPI` builds these the way `make_boolean_fixtures.py` does, in cm rather than mm,
    # so the ladder's own solids and these are the same constructions.
    def cyl_along(radius, length, origin, direction):
        return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*origin), gp_Dir(*direction)),
                                        radius, length).Shape()

    # Two orthogonal r = 1 cylinders intersected: the Steinmetz solid, whose whole boundary is
    # the transcendental cylinder-cylinder curve and which has no planar face at all.
    steinmetz = BRepAlgoAPI_Common(cyl_along(1.0, 6.0, (0, 0, -3), (0, 0, 1)),
                                   cyl_along(1.0, 6.0, (-3, 0, 0), (1, 0, 0))).Shape()
    steinmetz_record = expect("Steinmetz solid (two cylinders intersected)", steinmetz,
                              "cell-intersection", want_leaves=2)
    check("the Steinmetz solid reaches the cell emitter only after a rejection",
          (steinmetz_record.get("retriedAfter") or {}).get("recogniser") == "tier2-tube-union",
          f"retried after {(steinmetz_record.get('retriedAfter') or {}).get('recogniser')}")
    # A tube with a transverse hole: four halfspaces, one of them the hole wall, whose material
    # is OUTSIDE its own carrier and which therefore enters as a subtraction.
    window = BRepAlgoAPI_Cut(cyl_along(1.5, 6.0, (0, 0, -3), (0, 0, 1)),
                             cyl_along(0.8, 6.0, (-3, 0, 0), (1, 0, 0))).Shape()
    window_record = expect("tube with a transverse window", window, "cell-intersection",
                           want_leaves=2)
    check("the window's hole wall is a complemented leaf and its barrel is not",
          window_record["accepted"]
          and not window_record["candidate"]["leaves"][0].get("outside")
          and window_record["candidate"]["leaves"][1].get("outside") is True,
          window_record["description"])
    check("the barrel and its two caps folded into one TGeoTube",
          window_record["accepted"]
          and window_record["candidate"]["leaves"][0]["type"] == "TGeoTube"
          and abs(window_record["candidate"]["leaves"][0]["params"]["dz"] - 3.0) < 1e-12
          and window_record["candidate"]["notes"]["nCarriers"] == 4,
          f"{window_record['candidate']['notes'] if window_record['accepted'] else 'n/a'}")
    # A cylinder cut by an oblique plane: the cut face is an exact ellipse, and the remaining
    # cap plane folds into the tube while the oblique one stays a halfspace.
    oblique_knife = BRepPrimAPI_MakeBox(gp_Pnt(-20, -20, 0), 40.0, 40.0, 40.0).Shape()
    oblique_spin = gp_Trsf()
    oblique_spin.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), math.radians(60.0))
    oblique_lift = gp_Trsf()
    oblique_lift.SetTranslation(gp_Vec(0.0, 0.0, 2.5))
    oblique = BRepAlgoAPI_Cut(
        cyl_along(1.2, 5.0, (0, 0, 0), (0, 0, 1)),
        BRepBuilderAPI_Transform(oblique_knife, oblique_lift.Multiplied(oblique_spin),
                                 True).Shape()).Shape()
    expect("cylinder cut by an oblique plane", oblique, "cell-intersection", want_leaves=2)
    # A cube with an axial through-hole: six planes that are a TGeoBBox, and a hole wall.
    drilled = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(-2, -2, -2), 4.0, 4.0, 4.0).Shape(),
                              cyl_along(0.8, 6.0, (0, 0, -3), (0, 0, 1))).Shape()
    drilled_record = expect("cube with an axial through-hole", drilled, "cell-intersection",
                            want_leaves=2)
    check("the cube's six plane carriers folded into one TGeoBBox",
          drilled_record["accepted"]
          and drilled_record["candidate"]["leaves"][0]["type"] == "TGeoBBox"
          and drilled_record["candidate"]["notes"]["nCarriers"] == 7,
          drilled_record["description"])

    # --- rung 3 negative controls: the ladder of near-misses, and which test catches each ---
    # A shallow V notch cut into a cylinder is genuinely TWO cells -- the solid is the union of
    # two halfspaces' complements, not their intersection -- and the notch angle decides which of
    # the three tests notices. Every rung of this ladder must refuse; the last must accept,
    # because at that angle the part IS the cell to the resolution the pipeline declares.
    def notched_cylinder(angle):
        def knife(sign):
            slab = BRepPrimAPI_MakeBox(gp_Pnt(1.5, -10.0, -10.0), 20.0, 20.0, 20.0).Shape()
            spin = gp_Trsf()
            spin.SetRotation(gp_Ax1(gp_Pnt(1.5, 0.0, 0.0), gp_Dir(0, 0, 1)), sign * angle)
            return BRepBuilderAPI_Transform(slab, spin, True).Shape()
        return BRepAlgoAPI_Cut(cyl, BRepAlgoAPI_Common(knife(1.0), knife(-1.0)).Shape()).Shape()

    notch_trusted = expect_single_cell_declined(
        "a cylinder with a 2e-03 rad notch (a trusted concave edge)",
        notched_cylinder(2.0e-3), "trusted concave edge")
    check("the concave decline names how many edges it counted",
          "1 trusted concave edge(s) of 9" in (notch_trusted or ""),
          (notch_trusted or "")[:120])
    # And since rung 4 the notch above the trust filter is not merely refused as one cell, it is
    # CONVERTED as the two cells it is. That is the ladder's whole point read forwards: the
    # filter's job was always to say which side of the blind band a reflex is on, and on the
    # trusted side there is now something to do about it.
    notch_converted = process_solid(notched_cylinder(2.0e-3), "notched cylinder (trusted)")
    check("the notch above the trust filter converts as two cells, exactly",
          notch_converted["accepted"] and notch_converted["recogniser"] == "cells-union"
          and notch_converted["candidate"]["notes"]["nCells"] == 2
          and notch_converted["acceptance"]["symmetricDifference"] == 0.0,
          f"{notch_converted['recogniser']}: {notch_converted['description']}, "
          f"dV_sym={notch_converted['acceptance']['symmetricDifference'] if notch_converted['accepted'] else 'n/a'}")
    notch_gap = expect_declined("cylinder with a 1e-05 rad notch (below the trust filter)",
                                notched_cylinder(1.0e-5), "the cell's boundary is")
    check("the gap is what refuses the notch the trust filter let through",
          "the cell's boundary is" in (notch_gap["reason"] or "")
          and notch_gap["recogniser"] is None,
          (notch_gap["reason"] or "")[-140:])
    notch_volume = expect_declined("cylinder with a 1e-06 rad notch (ten model tolerances deep)",
                                   notched_cylinder(1.0e-6), "symmetric difference")
    check("the volume is what refuses a notch too shallow for the gap to see",
          notch_volume["recogniser"] == "cell-intersection"
          and not notch_volume["accepted"],
          (notch_volume["reason"] or "")[:140])
    # ... and the other edge of the same knife: one model tolerance deep is inside every
    # declared resolution and must be accepted, or the ladder above would prove nothing.
    expect("cylinder with a 1e-07 rad notch (one model tolerance deep)",
           notched_cylinder(1.0e-7), "cell-intersection", want_leaves=2)

    # A genuine two-cell body, asked of the cell emitter directly, because the cascade answers it
    # correctly as a union long before the emitter would see it.
    crossed = BRepAlgoAPI_Fuse(cyl_along(1.0, 6.0, (0, 0, -3), (0, 0, 1)),
                               cyl_along(1.0, 6.0, (-3, 0, 0), (1, 0, 0))).Shape()
    crossed_cand, crossed_why = recognise.recognise_single_cell(crossed)
    check("two fused cylinders are refused by the cell emitter, naming the concave edges",
          crossed_cand is None and "trusted concave edge(s)" in (crossed_why or ""),
          (crossed_why or "")[:120])
    # An all-planar body is the prism family's. Asked with a chamfered box -- convex, seven
    # planes, no concave edge at all, so the cell test itself would say yes.
    chamfer = BRepPrimAPI_MakeBox(gp_Pnt(1.2, -9.0, -9.0), 20.0, 20.0, 20.0).Shape()
    chamfer_spin = gp_Trsf()
    chamfer_spin.SetRotation(gp_Ax1(gp_Pnt(1.2, 0.0, 0.0), gp_Dir(0, 1, 0)),
                             math.radians(35.0))
    chamfered = BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeBox(gp_Pnt(-2, -2, -2), 4.0, 4.0, 4.0).Shape(),
        BRepBuilderAPI_Transform(chamfer, chamfer_spin, True).Shape()).Shape()
    planar_cand, planar_why = recognise.recognise_single_cell(chamfered)
    check("an all-planar solid is handed to the prism family, not read as halfspaces",
          planar_cand is None and "belongs to the prism family" in (planar_why or ""),
          (planar_why or "")[:120])
    # The L-plate has a concave edge, so the cell emitter refuses it on that count instead.
    ell_cand, ell_why = recognise.recognise_single_cell(ell)
    check("the L-plate is refused by the cell emitter on its concave edge",
          ell_cand is None and "trusted concave edge(s)" in (ell_why or ""),
          (ell_why or "")[:120])

    # --- the floor: the parts the earlier rungs own are not intercepted ---
    for name, want in (("L-shaped plate", "rung2-xtru"),
                       ("placed Xtru (non-convex L section)", "rung2-xtru"),
                       ("stepped polycone (duplicate z planes)", "revolved-pcon"),
                       ("box", "tier1-box"), ("tube", "tier1-tube"),
                       ("rod-and-eye (two-cluster union)", "tier2-tube-union")):
        check(f"{name} is still recognised as {want}", seen_recognisers.get(name) == want,
              f"{seen_recognisers.get(name)}")

    check("every single-cell candidate is byte-identical to its recorded digest",
          all(seen_digests.get(name) == digest for name, digest
              in _CELL_CANDIDATE_DIGESTS.items()),
          "; ".join(f"{name}: {seen_digests.get(name)} != {digest}" for name, digest
                    in _CELL_CANDIDATE_DIGESTS.items()
                    if seen_digests.get(name) != digest)
          or f"{len(_CELL_CANDIDATE_DIGESTS)} candidates unchanged")

    # --- the description must refuse an ill-formed intersection ---
    unit_box = prim.leaf("TGeoBBox", {"dx": 1.0, "dy": 1.0, "dz": 1.0}, prim.identity_frame())
    hole = prim.leaf("TGeoTube", {"rmin": 0.0, "rmax": 0.5, "dz": 2.0},
                     prim.identity_frame(), True)
    for label, op, leaves in (("a single leaf", "intersection", [unit_box]),
                              ("a complement first", "intersection", [hole, unit_box]),
                              ("a complement in a union", "union", [unit_box, hole])):
        try:
            prim.candidate(op, leaves, "self-test")
            refused = False
        except ValueError:
            refused = True
        check(f"a candidate with {label} is refused", refused)

    # --- flat-CSG R1: the torus carrier ---
    def torus_at(major, minor, angle=None, origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0),
                 ref=(1.0, 0.0, 0.0)):
        axis = gp_Ax2(gp_Pnt(*origin), gp_Dir(*direction), gp_Dir(*ref))
        maker = (BRepPrimAPI_MakeTorus(axis, major, minor) if angle is None
                 else BRepPrimAPI_MakeTorus(axis, major, minor, angle))
        maker.Build()
        return maker.Shape()

    def expect_torus(name, solid, r, rmin, rmax, phi1=0.0, dphi=360.0):
        record = expect(name, solid, "tier1-torus")
        if not record["accepted"]:
            return record
        p = record["candidate"]["leaves"][0]["params"]
        want = {"r": r, "rmin": rmin, "rmax": rmax, "phi1": phi1, "dphi": dphi}
        worst = max(abs(p[k] - v) for k, v in want.items())
        check(f"{name} reconstructs the source TGeoTorus parameters", worst < 1.0e-9,
              f"worst parameter deviation {worst:.3g}")
        return record

    solid_torus = torus_at(4.0, 1.0)
    solid_torus_record = expect_torus("solid torus", solid_torus, 4.0, 0.0, 1.0)
    # A shell: two concentric tori of the same major radius, which is a bellows ply's section.
    ply = BRepAlgoAPI_Cut(torus_at(5.0, 0.30), torus_at(5.0, 0.28)).Shape()
    ply_record = expect_torus("torus shell (a bellows ply)", ply, 5.0, 0.28, 0.30)
    # A phi wedge, hollow, whose two cut planes pass through the axis.
    wedge_torus = BRepAlgoAPI_Cut(
        torus_at(4.0, 1.0, math.radians(120.0), ref=(math.cos(math.radians(20.0)),
                                                     math.sin(math.radians(20.0)), 0.0)),
        torus_at(4.0, 0.8, math.radians(120.0) + 1.0e-4,
                 ref=(math.cos(math.radians(20.0)), math.sin(math.radians(20.0)), 0.0))).Shape()
    wedge_torus_record = expect_torus("hollow torus wedge", wedge_torus,
                                      4.0, 0.8, 1.0, 20.0, 120.0)
    # Placed, so the frame machinery is exercised on a torus too.
    torus_spin = gp_Trsf()
    torus_spin.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), 0.7)
    torus_shift = gp_Trsf()
    torus_shift.SetTranslation(gp_Vec(3.0, -4.0, 5.0))
    placed_torus = BRepBuilderAPI_Transform(solid_torus,
                                            torus_shift.Multiplied(torus_spin), True).Shape()
    placed_torus_record = expect("placed torus", placed_torus, "tier1-torus")
    check("a placed torus travels as one leaf plus a rigid placement",
          placed_torus_record["accepted"]
          and prim.placement_for_candidate(placed_torus_record["candidate"]) is not None,
          "placement present" if placed_torus_record["accepted"] else "not accepted")

    # The torus as a cell-emitter carrier: a ply cut by a plane is a cell of two toroidal
    # halfspaces, the bore's one complemented, and one box.
    half_ply = BRepAlgoAPI_Common(
        ply, BRepPrimAPI_MakeBox(gp_Pnt(-10, -10, 0), 20.0, 20.0, 20.0).Shape()).Shape()
    half_ply_record = expect("half a bellows ply", half_ply, "cell-intersection", want_leaves=3)
    check("the ply's bore enters the cell as a complemented TGeoTorus",
          half_ply_record["accepted"]
          and sum(1 for lf in half_ply_record["candidate"]["leaves"]
                  if lf["type"] == "TGeoTorus") == 2
          and any(lf.get("outside") and lf["type"] == "TGeoTorus"
                  for lf in half_ply_record["candidate"]["leaves"]),
          half_ply_record["description"])

    # --- R1 negative controls ---
    # A torus fused with a cylinder through it: two trusted concave edges, genuinely two cells,
    # and the fixture ladder's `torus_union_cyl`. It must decline, and say which of the two it
    # failed on.
    # `torus_union_cyl` from the fixture ladder, in cm: the cylinder radius sits inside the
    # tube band so the junction really exists, and it is concave on both circles. Two cells,
    # and the emitter must say so rather than propose the torus it partly is.
    torus_cyl = BRepAlgoAPI_Fuse(
        torus_at(2.5, 0.8),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -2.0), gp_Dir(0, 0, 1)),
                                 2.0, 4.0).Shape()).Shape()
    torus_cyl_why = expect_single_cell_declined("a torus fused with a coaxial cylinder through it",
                                               torus_cyl, "trusted concave edge")
    # The torus template's own verdict, asked of the template rather than read out of a decline
    # chain: since rung 4 the part converts, so there is no chain to read it out of, but what the
    # control was written to protect -- that the template says what it found instead of being
    # silently skipped -- is asserted here directly.
    torus_cyl_records, _why = recognise._face_records(torus_cyl)
    torus_cyl_diag = recognise._bbox_diagonal(torus_cyl)
    try:
        recognise._match_torus(torus_cyl, torus_cyl_records,
                               recognise.REL_TOL * max(torus_cyl_diag, 1.0), torus_cyl_diag)
        torus_template_why = "the template accepted it"
    except recognise.Declined as declined:
        torus_template_why = str(declined)
    check("the torus template says what it found before the cell test refuses it",
          "is not a whole torus" in torus_template_why, torus_template_why[:120])
    # A shell whose bore is displaced off the barrel's axis, as a ladder, because the answer
    # changes with the displacement and only a ladder says where. `tol` here is the recogniser's
    # declared resolution, REL_TOL times the part's 16.2 cm diagonal, i.e. 1.6e-05 cm.
    #
    #   below tol : the template merges the two carriers and proposes ONE coaxial torus. That is
    #               safe rather than lucky -- a rigid shift of the bore is volume-preserving, and
    #               the measured symmetric difference of the merged proposal is 7.1e-15 cm^3
    #               against a band of 1.1e-05, so the two are the same solid to every resolution
    #               the pipeline declares. Refusing here would be the wrong answer.
    #   above tol : the carriers separate, the template refuses on concentricity, and the cell
    #               emitter represents the part exactly as two toroidal halfspaces.
    #
    # So this shape has no misrepresenting displacement, and that is the claim being asserted --
    # not "a near-miss declines", which would be false.
    for displacement, want in ((1.0e-6, "tier1-torus"), (1.0e-5, "tier1-torus"),
                               (3.0e-5, "cell-intersection"), (1.0e-3, "cell-intersection")):
        skewed = BRepAlgoAPI_Cut(
            torus_at(5.0, 0.30),
            torus_at(5.0, 0.28, origin=(displacement, 0.0, 0.0))).Shape()
        skewed_record = process_solid(skewed, f"shell, bore {displacement:g} cm off axis")
        acceptance = skewed_record.get("acceptance") or {}
        check(f"a shell whose bore is {displacement:g} cm off the axis converts as {want}, "
              "within the band",
              skewed_record["accepted"] and skewed_record["recogniser"] == want
              and acceptance.get("symmetricDifference", 1.0) <= acceptance.get("band", 0.0),
              f"{skewed_record['recogniser']}: dV="
              f"{acceptance.get('symmetricDifference')} band={acceptance.get('band')}")
        if want == "cell-intersection":
            # And it is the concentricity test that hands it over, not an accident further on.
            records, _reason = recognise._face_records(skewed)
            skewed_diag = recognise._bbox_diagonal(skewed)
            try:
                recognise._match_torus(skewed, records,
                                       recognise.REL_TOL * max(skewed_diag, 1.0), skewed_diag)
                refused = None
            except recognise.Declined as declined:
                refused = str(declined)
            check(f"and the torus template is what refuses it at {displacement:g} cm",
                  refused is not None and "concentric" in refused,
                  refused or "IT PROPOSED ONE")

    # --- the fillet blend that used to kill the conversion ---
    # ALICE3 carries fillet blends whose torus has a minor radius LARGER than its major one --
    # a self-intersecting torus, which `TGeoTorus` cannot state and the validator rightly
    # refuses. `_cell_leaf` called `prim.leaf` without catching that, so the refusal escaped
    # `recognise()` (which catches only `Declined`) and killed the whole ALICE3 conversion on
    # the first blend it met. The numbers below are the ones from that crash.
    blend_lobe = BRepAlgoAPI_Common(
        torus_at(0.0428825434729, 0.1),
        BRepPrimAPI_MakeBox(gp_Pnt(0.06, -1.0, -1.0), 2.0, 2.0, 2.0).Shape()).Shape()
    blend_record = expect_declined("a lobe of a self-intersecting fillet torus", blend_lobe,
                                   "self-intersecting torus")
    check("the fillet blend reaches the cell path and declines there, naming the blend",
          "as a single cell: TGeoTorus: rmax" in (blend_record["reason"] or "")
          and "fillet blend" in (blend_record["reason"] or ""),
          (blend_record["reason"] or "")[:150])
    # ... and the refusal it declines on is a real one: the description layer does refuse those
    # numbers, so the control is not passing because nothing was checked.
    try:
        prim.leaf("TGeoTorus", {"r": 0.0428825434729, "rmin": 0.0, "rmax": 0.1,
                                "phi1": 0.0, "dphi": 360.0}, prim.identity_frame())
        refused_kind = None
    except prim.InvalidDescription:
        refused_kind = "InvalidDescription"
    except ValueError:
        refused_kind = "ValueError"
    check("the description layer refuses those numbers as an illegal solid",
          refused_kind == "InvalidDescription", f"raised {refused_kind}")
    # The altitude matters: an illegal SOLID declines, a bug in the matcher still raises. Without
    # this half the wrapper could swallow a missing parameter and nobody would know.
    for label, kind, params in (("a missing parameter", "TGeoTorus", {"r": 1.0}),
                                ("an unknown leaf type", "TGeoNotAShape", {})):
        try:
            recognise._leaf(kind, params, prim.identity_frame())
            outcome = "returned a leaf"
        except recognise.Declined:
            outcome = "declined"
        except ValueError:
            outcome = "raised"
        check(f"{label} still raises rather than declining", outcome == "raised", outcome)

    # --- flat-CSG R2: the elliptic cylinder ---
    def elliptic_cylinder(a, b, dz, ref=None):
        axis = gp_Ax2(gp_Pnt(0, 0, -dz), gp_Dir(0, 0, 1),
                      gp_Dir(*(ref if ref is not None else (1.0, 0.0, 0.0))))
        major, minor = max(a, b), min(a, b)
        edge = BRepBuilderAPI_MakeEdge(gp_Elips(axis, major, minor)).Edge()
        face = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(edge).Wire()).Face()
        prism = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, 2 * dz))
        prism.Build()
        return prism.Shape()

    def expect_eltu(name, solid, a, b, dz):
        record = expect(name, solid, "tier1-eltu")
        if not record["accepted"]:
            return record
        p = record["candidate"]["leaves"][0]["params"]
        worst = max(abs(p["a"] - a), abs(p["b"] - b), abs(p["dz"] - dz))
        check(f"{name} reconstructs the source TGeoEltu parameters", worst < 1.0e-9,
              f"a={p['a']:.6g} b={p['b']:.6g} dz={p['dz']:.6g}, worst {worst:.3g}")
        return record

    # `conv_eltu` puts the major axis on y when the source had a < b, so both orders are built
    # the way the writer builds them and both must come back as the source's own two numbers.
    eltu_solid = elliptic_cylinder(3.0, 1.5, 5.0)
    eltu_record = expect_eltu("elliptic cylinder, a > b", eltu_solid, 3.0, 1.5, 5.0)
    expect_eltu("elliptic cylinder, a < b",
                elliptic_cylinder(1.5, 3.0, 5.0, ref=(0.0, 1.0, 0.0)), 1.5, 3.0, 5.0)
    # a == b is a circle, and it is still a TGeoEltu: the carrier is an extrusion, never a
    # cylinder, so nothing can confuse the two. Asserted rather than left to chance.
    circle_eltu = expect_eltu("elliptic cylinder with equal semi-axes",
                              elliptic_cylinder(2.0, 2.0, 5.0), 2.0, 2.0, 5.0)
    check("an ellipse with equal semi-axes stays a TGeoEltu and is not read as a tube",
          circle_eltu["accepted"]
          and circle_eltu["candidate"]["leaves"][0]["type"] == "TGeoEltu",
          circle_eltu["description"])
    placed_eltu = BRepBuilderAPI_Transform(elliptic_cylinder(3.0, 1.5, 5.0),
                                           torus_shift.Multiplied(torus_spin), True).Shape()
    placed_eltu_record = expect("placed elliptic cylinder", placed_eltu, "tier1-eltu")
    check("a placed elliptic cylinder travels as one leaf plus a rigid placement",
          placed_eltu_record["accepted"]
          and prim.placement_for_candidate(placed_eltu_record["candidate"]) is not None,
          "placement present" if placed_eltu_record["accepted"] else "not accepted")

    # --- R2 negative controls ---
    # An extruded oval that is NOT an ellipse: a B-spline racetrack. Its basis curve is not a
    # GeomAbs_Ellipse, so it stays free-form and says so.
    racetrack = []
    for i in range(24):
        ang = 2.0 * math.pi * i / 24.0
        racetrack.append(gp_Pnt(3.0 * math.cos(ang),
                                1.5 * math.sin(ang) * (1.0 + 0.15 * math.cos(2 * ang)), -5.0))
    spline_pts = TColgp_HArray1OfPnt(1, len(racetrack))
    for i, pnt in enumerate(racetrack, start=1):
        spline_pts.SetValue(i, pnt)
    interp = GeomAPI_Interpolate(spline_pts, True, 1.0e-7)
    interp.Perform()
    oval_edge = BRepBuilderAPI_MakeEdge(interp.Curve()).Edge()
    oval_face = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(oval_edge).Wire()).Face()
    oval_prism = BRepPrimAPI_MakePrism(oval_face, gp_Vec(0, 0, 10.0))
    oval_prism.Build()
    expect_declined("extruded B-spline racetrack (not an ellipse)", oval_prism.Shape(),
                    "free-form faces")

    # --- both new templates' instruments must be able to say "no" ---
    # Ten model tolerances on one semi-axis, reported at their true size by the same measurement
    # that gates the proposal.
    true_eltu = prim.candidate("primitive", [prim.leaf(
        "TGeoEltu", {"a": 3.0, "b": 1.5, "dz": 5.0}, prim.identity_frame())], "self-test")
    nudged_eltu = prim.candidate("primitive", [prim.leaf(
        "TGeoEltu", {"a": 3.0 + 1.0e-6, "b": 1.5, "dz": 5.0},
        prim.identity_frame())], "self-test")
    eltu_gap = recognise._boundary_gap(prim.build_occ(true_eltu), prim.build_occ(nudged_eltu))
    check("the gap reports a semi-axis displaced by ten model tolerances",
          abs(eltu_gap - 1.0e-6) < 1.0e-9, f"gap {eltu_gap:.3g} cm, expected 1e-06 cm")
    true_torus = prim.candidate("primitive", [prim.leaf(
        "TGeoTorus", {"r": 4.0, "rmin": 0.0, "rmax": 1.0, "phi1": 0.0, "dphi": 360.0},
        prim.identity_frame())], "self-test")
    nudged_torus = prim.candidate("primitive", [prim.leaf(
        "TGeoTorus", {"r": 4.0, "rmin": 0.0, "rmax": 1.0 + 1.0e-6, "phi1": 0.0, "dphi": 360.0},
        prim.identity_frame())], "self-test")
    torus_gap = recognise._boundary_gap(prim.build_occ(true_torus), prim.build_occ(nudged_torus))
    check("the gap reports a tube radius displaced by ten model tolerances",
          abs(torus_gap - 1.0e-6) < 1.0e-9, f"gap {torus_gap:.3g} cm, expected 1e-06 cm")

    # --- the descriptions must refuse illegal parameters ---
    for label, kind, params in (
            ("a torus whose tube is wider than its major radius", "TGeoTorus",
             {"r": 1.0, "rmin": 0.0, "rmax": 2.0, "phi1": 0.0, "dphi": 360.0}),
            ("a torus with rmin above rmax", "TGeoTorus",
             {"r": 4.0, "rmin": 1.0, "rmax": 0.5, "phi1": 0.0, "dphi": 360.0}),
            ("a torus with dphi zero", "TGeoTorus",
             {"r": 4.0, "rmin": 0.0, "rmax": 1.0, "phi1": 0.0, "dphi": 0.0}),
            ("an elliptic cylinder with a zero semi-axis", "TGeoEltu",
             {"a": 0.0, "b": 1.5, "dz": 5.0})):
        try:
            prim.leaf(kind, params, prim.identity_frame())
            refused = False
        except ValueError:
            refused = True
        check(f"a description of {label} is refused", refused)

    check("every torus and elliptic-cylinder candidate is byte-identical to its recorded digest",
          all(seen_digests.get(name) == digest for name, digest
              in _TORUS_ELTU_CANDIDATE_DIGESTS.items()),
          "; ".join(f"{name}: {seen_digests.get(name)} != {digest}" for name, digest
                    in _TORUS_ELTU_CANDIDATE_DIGESTS.items()
                    if seen_digests.get(name) != digest)
          or f"{len(_TORUS_ELTU_CANDIDATE_DIGESTS)} candidates unchanged")

    # --- Tier 0: the quadric a stored B-spline face already is ---
    #
    # `Stream_K_Tier0.md` §3's lesson, applied to this rung: the instrument is checked against a
    # known displacement BEFORE any face is. Every control below is built in-process from OCC
    # primitives, `BRepBuilderAPI_NurbsConvert` being the exporter artefact the whole service
    # exists for, rebuilt.
    from csg import tier0
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform, BRepBuilderAPI_NurbsConvert
    from OCC.Core.BRepTools import breptools
    from OCC.Core.GeomAbs import GeomAbs_Cylinder
    from OCC.Core.TopAbs import TopAbs_FACE
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.gp import gp_Ax3, gp_Cylinder, gp_GTrsf, gp_Mat

    check("the canonicaliser's band is the cascade's own",
          tier0.REL_TOL == recognise.REL_TOL,
          f"tier0 {tier0.REL_TOL:.0e} vs recognise {recognise.REL_TOL:.0e}")

    def nurbs(shape):
        return BRepBuilderAPI_NurbsConvert(shape, True).Shape()

    def faces_of(shape):
        found = []
        walk = TopExp_Explorer(shape, TopAbs_FACE)
        while walk.More():
            found.append(topods.Face(walk.Current()))
            walk.Next()
        return found

    def samples_of(face, n):
        adaptor = BRepAdaptor_Surface(face, True)
        import O2_CADtoTGeo as converter
        return converter._sample_surface_for_recognition(adaptor, *breptools.UVBounds(face), n=n)

    # (c) the instrument: a model displaced by a known amount must be reported at that size.
    probe_cylinder = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                              5.0, 8.0).Shape()
    probe_points, _probe_normals = samples_of(
        [f for f in faces_of(probe_cylinder)
         if BRepAdaptor_Surface(f, True).GetType() == GeomAbs_Cylinder][0], 17)
    probe_sphere_points, _ = samples_of(faces_of(BRepPrimAPI_MakeSphere(5.0).Shape())[0], 17)
    probe_torus_points, _ = samples_of(faces_of(BRepPrimAPI_MakeTorus(6.0, 1.5).Shape())[0], 17)
    for displacement in (1.0e-3, 1.0e-6, 1.0e-9):
        for label, kind, model, points in (
                ("cylinder radius", "cylinder",
                 {"axis": [0.0, 0.0, 1.0], "origin": [0.0, 0.0, 0.0],
                  "radius": 5.0 + displacement}, probe_points),
                ("sphere radius", "sphere",
                 {"centre": [0.0, 0.0, 0.0], "radius": 5.0 + displacement},
                 probe_sphere_points),
                ("torus tube radius", "torus",
                 {"axis": [0.0, 0.0, 1.0], "centre": [0.0, 0.0, 0.0], "major": 6.0,
                  "minor": 1.5 + displacement}, probe_torus_points)):
            measured = tier0.surface_gap(kind, model, points)
            check(f"the gap reports a {label} displaced by {displacement:.0e} cm at its true size",
                  abs(measured - displacement) <= 1.0e-9 * displacement + 1.0e-13,
                  f"measured {measured:.6g} cm, displaced {displacement:.0e} cm")

    # Positive controls: the same solid, written as NURBS, must convert to the SAME BODY. This
    # is the sharpest available statement that canonicalisation is exact -- not merely that the
    # part converts, but that it converts to what its natively-analytic twin converts to.
    tier0_pairs = (
        ("box", BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0).Shape(), "tier1-box", 1),
        ("solid cylinder", BRepPrimAPI_MakeCylinder(ax, 2.0, 10.0).Shape(), "tier1-tube", 1),
        ("tube segment", BRepPrimAPI_MakeCylinder(ax, 2.0, 10.0, math.radians(72.0)).Shape(),
         "tier1-tubeseg", 1),
        ("cone", BRepPrimAPI_MakeCone(ax, 3.0, 1.0, 6.0).Shape(), "tier1-cone", 1),
        ("sphere", BRepPrimAPI_MakeSphere(3.0).Shape(), "tier1-sphere", 1),
        ("solid torus", BRepPrimAPI_MakeTorus(6.0, 1.5).Shape(), "tier1-torus", 1),
        ("hollow torus wedge",
         BRepAlgoAPI_Cut(BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                               6.0, 1.5, math.radians(140.0)).Shape(),
                         BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                               6.0, 0.7, math.radians(140.0)).Shape()).Shape(),
         "tier1-torus", 1),
        ("cube with an axial through-hole",
         BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(-3, -3, -3), 6.0, 6.0, 6.0).Shape(),
                         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -5), gp_Dir(0, 0, 1)),
                                                  1.5, 10.0).Shape()).Shape(),
         "cell-intersection", 2),
    )

    def realisation_gap(one, other):
        """The largest distance between the two candidates' realised boundaries, in cm.

        Compared as solids rather than as parameter lists, because the same solid has several
        equally correct descriptions: the fit is free to choose the sense of a cylinder's axis,
        which swaps a cone's two radii, and free to choose the direction phi is measured from,
        which shifts a wedge's `phi1` against its frame's rotation by the same angle. A distance
        between the two bodies sees through all of that and would still see a wrong radius.
        """
        if one is None or other is None:
            return float("inf")
        if (one["op"], one["recogniser"], len(one["leaves"])) != \
                (other["op"], other["recogniser"], len(other["leaves"])):
            return float("inf")
        if [lf["type"] for lf in one["leaves"]] != [lf["type"] for lf in other["leaves"]]:
            return float("inf")
        return recognise._boundary_gap(prim.build_occ(one), prim.build_occ(other))

    for label, solid, want_recogniser, want_leaves in tier0_pairs:
        native_record = process_solid(solid, f"tier0 native {label}")
        encoded = expect(f"NURBS-encoded {label}", nurbs(solid), want_recogniser, want_leaves)
        deviation = realisation_gap(native_record["candidate"], encoded["candidate"])
        notes = (encoded["candidate"] or {}).get("notes", {})
        check(f"the NURBS-encoded {label} realises the analytic one's solid",
              deviation <= 1.0e-9,
              f"{notes.get('tier0Faces', 0)} canonicalised carrier(s) at a worst gap of "
              f"{notes.get('tier0WorstGapRelative', float('nan')):.3g} of the part; the two "
              f"realisations are {deviation:.3g} cm apart")

    check("every Tier-0 candidate is byte-identical to its recorded digest",
          all(seen_digests.get(name) == digest for name, digest
              in _TIER0_CANDIDATE_DIGESTS.items()),
          "; ".join(f"{name}: {seen_digests.get(name)} != {digest}" for name, digest
                    in _TIER0_CANDIDATE_DIGESTS.items()
                    if seen_digests.get(name) != digest)
          or f"{len(_TIER0_CANDIDATE_DIGESTS)} candidates unchanged")

    # (a) the negative control that matters most: a genuinely free-form face must NOT
    # canonicalise, and its decline must carry a number rather than a bare "free-form". That
    # number is the smallest gap any PROPOSED model reaches -- not the distance to the nearest
    # quadric in the world, which nothing here computes. Said plainly because it is the sort of
    # claim this project has over-stated before.
    import O2_CADtoTGeo as converter
    for label, face in (
            ("free-form saddle", converter._self_test_bezier_patch(
                lambda s, t: (10 * s - 5, 10 * t - 5, (10 * s - 5) * (10 * t - 5) / 10.0), 6, 6)),
            ("narrow free-form ridge", converter._self_test_bezier_patch(
                lambda s, t: (20 * s - 10, 0.5 * t,
                              0.02 * (20 * s - 10) ** 2 + 0.3 * (20 * s - 10) * t), 6, 6)),
            ("swept non-circular profile (bulge 1e-2)",
             converter._self_test_tapered_near_circle(1.0e-2, 1.0e-4))):
        adaptor = BRepAdaptor_Surface(face, True)
        carrier, gap = tier0.canonicalise(face, adaptor, 20.0)
        check(f"a {label} is not canonicalised, and the gap says how far off it is",
              carrier is None and gap is not None and gap > 10.0 * tier0.REL_TOL * 20.0,
              f"{'declined' if carrier is None else 'ACCEPTED as ' + carrier['kind']}, best "
              f"proposal {gap:.4g} cm away, {gap / 20.0:.3g} of the part against "
              f"{tier0.REL_TOL:.0e}")

    # (b) a disguised cylinder displaced by ten model tolerances must be refused BY THE GAP --
    # and the same construction at a tenth of one tolerance must be accepted, because a criterion
    # that only ever says no is not a criterion either. An exact cylinder is squashed along x by
    # a `gp_GTrsf`, which scales the poles of the rational patch exactly, so what comes back is
    # an elliptic cylinder stored as a B-spline whose distance to the best circle is known in
    # closed form: `(a - b) / 2 = R eps / 2`.
    squash_radius, squash_scale = 5.0, 20.0
    squash_base = nurbs(BRepBuilderAPI_MakeFace(
        gp_Cylinder(gp_Ax3(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)),
                    squash_radius), 0.0, 2.0 * math.pi, 1.0, 9.0).Shape())
    squash_measured = {}
    for multiple in (0.1, 1.0, 10.0):
        intended = multiple * tier0.REL_TOL * squash_scale
        transform = gp_GTrsf()
        transform.SetVectorialPart(gp_Mat(1.0 + 2.0 * intended / squash_radius, 0.0, 0.0,
                                          0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
        squashed = faces_of(BRepBuilderAPI_GTransform(squash_base, transform, True).Shape())[0]
        carrier, gap = tier0.canonicalise(squashed, BRepAdaptor_Surface(squashed, True),
                                          squash_scale)
        squash_measured[multiple] = gap
        want_accepted = multiple < 1.0
        check(f"a disguised cylinder displaced by {multiple:g} model tolerance(s) is "
              f"{'accepted' if want_accepted else 'refused by the gap'}",
              (carrier is not None) == want_accepted,
              f"{'accepted as ' + carrier['kind'] if carrier else 'declined'}, "
              f"measured gap {gap:.4g} cm, {gap / squash_scale:.3g} of the part against "
              f"{tier0.REL_TOL:.0e}")
    ratios = [squash_measured[m] / (m * tier0.REL_TOL * squash_scale) for m in (0.1, 1.0, 10.0)]
    check("the measured gap is proportional to the displacement that caused it",
          max(ratios) - min(ratios) <= 1.0e-3 * max(ratios),
          f"gap / displacement = {', '.join(f'{r:.4f}' for r in ratios)}")

    # Canonicalisation carried ALICE3's `ST0923290_01#b12` as far as the cell emitter, where its
    # four pins on four parallel axes fold to an intersection of nothing. The gap then reported
    # every sample as unmeasurable and blamed OCCT; an empty proposal has to say so itself.
    empty_common = BRepAlgoAPI_Common(
        BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1.0, 1.0, 1.0).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(9, 9, 9), 1.0, 1.0, 1.0).Shape()).Shape()
    try:
        recognise._boundary_gap(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 1.0, 1.0, 1.0).Shape(),
                                empty_common)
        empty_reason = "no decline"
    except recognise.Declined as declined:
        empty_reason = str(declined)
    check("an empty proposal declines as empty, not as an OCCT measurement failure",
          "the proposal is empty" in empty_reason, empty_reason)

    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.TopoDS import TopoDS_Compound

    # --- Rung 4: the union of cells ---
    #
    # `Stream_AA_FlatCSG.md` §5 step 3. Every body here is one whose cell count is known in
    # closed form, which is what makes "2 cells" or "3 cells" an assertion rather than an
    # observation; the probe's own self-test asserts the same counts on the same shapes through
    # the same loop (`probes/cellCountProbe.py --self-test`, 7/7).
    from csg import decompose as decomp

    def expect_cells(name, solid, want_cells, want_leaves=None, **kwargs):
        record = process_solid(solid, name, **kwargs)
        seen_recognisers[name] = record["recogniser"] if record["accepted"] else None
        if record["accepted"]:
            seen_digests[name] = hashlib.sha256(
                json.dumps(record["candidate"], sort_keys=True).encode()).hexdigest()
        notes = (record["candidate"] or {}).get("notes", {})
        ok = (record["accepted"] and record["recogniser"] == "cells-union"
              and notes.get("nCells") == want_cells
              and (want_leaves is None or notes.get("nLeaves") == want_leaves))
        detail = (f"{notes.get('nCells')} cell(s) of {notes.get('cellLeaves')} leaves, "
                  f"{notes.get('nSplits')} split(s), volume drift "
                  f"{notes.get('volumeDriftRelative', float('nan')):.3g}, gap "
                  f"{notes.get('cellGapCm', float('nan')):.3g} cm"
                  if record["accepted"] else f"declined: {record['reason']}")
        check(f"{name} converts as {want_cells} cells", ok, detail)
        return record

    l_plate = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 4.0, 4.0, 1.0).Shape(),
                              BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -1), 4.0, 4.0, 3.0).Shape()).Shape()
    grooved = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 6.0, 4.0, 3.0).Shape(),
                              BRepPrimAPI_MakeBox(gp_Pnt(2, -1, 1), 2.0, 6.0, 3.0).Shape()).Shape()
    # An L-plate and a grooved block are the prism family's before they are the decomposition's,
    # and that ordering is the point: the multi-cell path runs only on what everything above it
    # declines. They are driven through `recognise_union_of_cells` directly so that the cell
    # counts are still asserted on the shapes whose counts are known.
    for label, solid, want in (("an L-plate", l_plate, 2), ("a grooved block", grooved, 3)):
        cand, why = recognise.recognise_union_of_cells(solid)
        gap = (None if cand is None else
               recognise._boundary_gap(prim.build_occ(cand), solid))
        check(f"{label} decomposes into {want} cells and realises the solid",
              cand is not None and cand["notes"]["nCells"] == want and gap <= 1.0e-9,
              (f"{cand['notes']['nCells']} cells of {cand['notes']['cellLeaves']} leaves, "
               f"{gap:.3g} cm from the part" if cand else f"declined: {why}"))

    # A hexagonal collar on a cylinder: genuinely two cells, and one of them is eight halfspaces
    # wide, so it exercises a cell that is neither a box nor a tube.
    hex_collar = BRepAlgoAPI_Fuse(
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -5), gp_Dir(0, 0, 1)), 3.0, 5.0).Shape(),
        swept_polygon(3.0, 6, 0.0, 5.0)).Shape()
    expect_cells("a cylinder with a hexagonal collar", hex_collar, 2, want_leaves=9)

    # Two rods on parallel axes, sharing no edge at all: R3 §6.1's lesson as a control. It has
    # ZERO trusted concave edges, so nothing but the connectivity split can find its two cells,
    # and before that split existed the cell emitter read it as an empty intersection.
    disjoint = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(disjoint)
    builder.Add(disjoint, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.0, 5.0).Shape())
    builder.Add(disjoint, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(6, 0, 0), gp_Dir(0, 0, 1)), 1.0, 5.0).Shape())
    disjoint_record = expect_cells("two rods sharing no edge", disjoint, 2, want_leaves=2)
    check("the disjoint pair is found by connectivity and needs no split at all",
          disjoint_record["accepted"]
          and disjoint_record["candidate"]["notes"]["nComponents"] == 2
          and disjoint_record["candidate"]["notes"]["nSplits"] == 0
          and decomp.count_trusted_concave(disjoint) == 0,
          f"{decomp.count_trusted_concave(disjoint)} trusted concave edge(s), "
          f"{(disjoint_record['candidate'] or {}).get('notes', {}).get('nSplits')} split(s)")

    # A torus with a cylinder through it -- the fixture ladder's `torus_union_cyl`, and the one
    # multi-cell body in the suite whose cells are not all planar.
    torus_through = BRepAlgoAPI_Fuse(
        torus_at(2.5, 0.8),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -2.0), gp_Dir(0, 0, 1)),
                                 2.0, 4.0).Shape()).Shape()
    expect_cells("a torus with a cylinder through it", torus_through, 2, want_leaves=3)

    # (a) THE VOLUME GUARD, against a real corruption rather than a turned knob: a decomposition
    # that genuinely loses a piece. Three disjoint boxes, with the component walk made to hand
    # back only two of them -- so the pieces really do sum to less than the part, by the third
    # box's volume -- must be refused, and the decline must say by how much.
    three_boxes = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(three_boxes)
    for x in (0.0, 4.0, 8.0):
        builder.Add(three_boxes, BRepPrimAPI_MakeBox(gp_Pnt(x, 0, 0), 2.0, 2.0, 2.0).Shape())
    expect_cells("three disjoint boxes", three_boxes, 3, want_leaves=3)
    intact_components = decomp.solid_components
    try:
        decomp.solid_components = lambda shape: intact_components(shape)[:-1]
        _lost, lost_why = recognise.recognise_union_of_cells(three_boxes)
    finally:
        decomp.solid_components = intact_components
    check("a decomposition that loses a cell is refused by the volume guard",
          _lost is None and "volume" in (lost_why or ""), lost_why or "ACCEPTED")
    check("the volume guard reports the drift it measured, at its true size",
          _lost is None and "0.333" in (lost_why or ""),
          f"one box of three is 1/3 of the part; the decline says: "
          f"{(lost_why or '')[:120]}")

    # (b) the budgets, each declining by name. Passed as arguments rather than by reaching into
    # the module, so the control says what it is testing.
    _over_cells, cells_why = recognise.recognise_union_of_cells(grooved, max_cells=2)
    check("a part over the cell budget declines naming the bound",
          _over_cells is None and "cell budget of 2" in (cells_why or ""), cells_why or "ACCEPTED")
    _over_leaves, leaves_why = recognise.recognise_union_of_cells(grooved, max_leaves=2)
    check("a part over the leaf budget declines naming the bound",
          _over_leaves is None and "part budget of 2" in (leaves_why or ""),
          leaves_why or "ACCEPTED")

    # (c) the DNF is two levels and the emitter refuses a third.
    flat_cell = prim.cell("primitive", [prim.leaf("TGeoBBox", {"dx": 1.0, "dy": 1.0, "dz": 1.0},
                                                 prim.identity_frame())])
    for label, cells_in in (
            ("a cell that is itself a union",
             [flat_cell, {"op": "union", "leaves": [flat_cell["leaves"][0]] * 2}]),
            ("a cell carrying a recogniser of its own",
             [flat_cell, {"op": "primitive", "leaves": flat_cell["leaves"],
                          "recogniser": "nested"}]),
            ("a single cell called a union", [flat_cell])):
        try:
            prim.union_of_cells(cells_in, "self-test")
            refused = False
        except (ValueError, prim.InvalidDescription):
            refused = True
        check(f"a description with {label} is refused", refused)

    # The balanced tree, asserted as a depth rather than as a hope: N cells must give a
    # ceil(log2 N) tree, which is the whole of `Stream_P` §3 reading 4's free win.
    if with_root:
        import ROOT as _ROOT

        def union_depth(shape):
            if shape.ClassName() != "TGeoCompositeShape":
                return 0
            node = shape.GetBoolNode()
            return 1 + max(union_depth(node.GetLeftShape()), union_depth(node.GetRightShape()))

        ladder = []
        for n_cells in (2, 3, 5, 8):
            comp = TopoDS_Compound()
            builder = BRep_Builder()
            builder.MakeCompound(comp)
            for i in range(n_cells):
                builder.Add(comp, BRepPrimAPI_MakeBox(gp_Pnt(4.0 * i, 0, 0),
                                                      2.0, 2.0, 2.0).Shape())
            cand, why = recognise.recognise_union_of_cells(comp)
            shape, _placement = prim.build_root(cand, f"balanced{n_cells}") if cand else (None, None)
            want = math.ceil(math.log2(n_cells))
            got = union_depth(shape) if shape is not None else -1
            ladder.append((n_cells, got, want))
            check(f"{n_cells} cells emit a balanced union tree of depth {want}", got == want,
                  f"depth {got}" if cand else f"declined: {why}")
        check("the union tree's depth is logarithmic in the cell count, not linear",
              all(got == want for _n, got, want in ladder),
              ", ".join(f"{n}->{got}" for n, got, _w in ladder))

    # The description has to survive the round trip through `csg_<part>.json`, because that is
    # how a two-level candidate reaches ROOT at all under `runOracleGate.py`: the converter
    # there has pythonOCC and no PyROOT, writes the JSON, and `--from-json` completes it in a
    # second interpreter. A `unionOfCells` is the first description with a nested level to make
    # that trip.
    if with_root:
        round_trip = json.loads(json.dumps(disjoint_record["candidate"]))
        rebuilt, rebuilt_placement = prim.build_root(round_trip, "roundtrip")
        direct, _direct_placement = prim.build_root(disjoint_record["candidate"], "direct")
        check("a two-level description survives the JSON round trip byte for byte",
              json.dumps(round_trip, sort_keys=True)
              == json.dumps(disjoint_record["candidate"], sort_keys=True)
              and rebuilt.ClassName() == direct.ClassName() and rebuilt_placement is None,
              f"{rebuilt.ClassName()}, placement "
              f"{'present' if rebuilt_placement else 'absent'}")
        gap = recognise._boundary_gap(prim.build_occ(round_trip),
                                      prim.build_occ(disjoint_record["candidate"]))
        check("the round-tripped description realises the same solid", gap <= 1.0e-12,
              f"{gap:.3g} cm apart")

    check("every union-of-cells candidate is byte-identical to its recorded digest",
          all(seen_digests.get(name) == digest for name, digest
              in _UNION_OF_CELLS_CANDIDATE_DIGESTS.items()),
          "; ".join(f"{name}: {seen_digests.get(name)} != {digest}" for name, digest
                    in _UNION_OF_CELLS_CANDIDATE_DIGESTS.items()
                    if seen_digests.get(name) != digest)
          or f"{len(_UNION_OF_CELLS_CANDIDATE_DIGESTS)} candidates unchanged")

    # --- the flat emitter's sign convention, measured against the shipped cell emitter --------
    # `Design_FlatCSGSolid.md` section 3.3: an inverted halfspace is still a solid, so the sign
    # composition -- the carrier's own orientation against the `side` field -- is the one thing
    # here that can be silently wrong. `recognise._cell_leaf` has shipped 1775 known-source-clean
    # parts and bounds each carrier into a padded native primitive whose conjunction IS the cell
    # over the part's neighbourhood, so it is the oracle the flat halfspaces are measured against.
    import struct
    from csg import decompose, flat as flatmod

    def _flat_gradient(block, point):
        """`|grad f|` at a point, for turning a quadric value into a first-order distance."""
        c = block["c"]
        x, y, z = point
        if block["kind"] == "torus":
            return 1.0                    # the torus block already IS a signed distance
        gx = 2.0 * (c[0] * x + c[1] * y + c[2] * z + c[6])
        gy = 2.0 * (c[1] * x + c[3] * y + c[4] * z + c[7])
        gz = 2.0 * (c[2] * x + c[4] * y + c[5] * z + c[8])
        return math.sqrt(gx * gx + gy * gy + gz * gz)

    def _flat_oracle(name, solid, seed=20260824, samples=4000):
        """The flat blocks of a one-cell solid, and `_cell_leaf`'s verdict on sampled points.

        Sampling is done once and the verdicts kept, so the sign-inversion controls below re-use
        them instead of rebuilding the padded OCCT solid per flip.

        Points within `REL_TOL x max(diag, 1)` of a carrier surface are not scored, and neither
        are points the classifier itself answers `TopAbs_ON` for. `|f| / |grad f|` is the
        first-order distance to that block's surface for every block type here -- with `grad`
        forced to 1 for the torus, whose block already is a signed distance -- so the band is a
        band around the boundary and nothing else. It cannot hide the defect this test exists to
        find: an inverted sign is a volumetric error, and a volumetric error has interior points.
        """
        import random
        from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier
        from OCC.Core.TopAbs import TopAbs_IN, TopAbs_ON
        from OCC.Core.gp import gp_Pnt
        diag = decompose.bbox_diagonal(solid)
        tol = recognise.REL_TOL * max(diag, 1.0)
        carriers = recognise._halfspace_carriers(solid, tol)
        box = recognise._CellBox(solid, diag)
        blocks = flatmod.blocks_from_carriers(carriers)
        leaves = [recognise._cell_leaf(c, box) for c in carriers]
        cand = prim.cell("intersection" if len(leaves) > 1 else "primitive", leaves)
        classifier = BRepClass3d_SolidClassifier(prim.build_occ(cand))
        rng = random.Random(seed)
        (xlo, ylo, zlo, xhi, yhi, zhi) = recognise._bbox_of(solid)
        scored = []
        for _ in range(samples):
            point = (rng.uniform(xlo, xhi), rng.uniform(ylo, yhi), rng.uniform(zlo, zhi))
            near = min(abs(flatmod.eval_block(b, point))
                       / max(_flat_gradient(b, point), 1.0e-300) for b in blocks)
            if near <= tol:
                continue
            classifier.Perform(gp_Pnt(*point), tol)
            state = classifier.State()
            if state == TopAbs_ON:
                continue
            scored.append((point, state == TopAbs_IN))
        worst_plane = max((flatmod.plane_scaling_error(b) for b in blocks
                           if flatmod.plane_scaling_error(b) is not None), default=None)
        return {"name": name, "solid": solid, "carriers": carriers, "blocks": blocks,
                "points": scored, "kinds": sorted({c["kind"] for c in carriers}),
                "sides": sorted({c["side"] for c in carriers}), "worstPlane": worst_plane}

    def _flat_disagreements(blocks, points):
        return sum(1 for point, occ_inside in points
                   if flatmod.flat_contains(blocks, point) != occ_inside)

    flat_axis = gp_Ax2(gp_Pnt(0, 0, -5), gp_Dir(0, 0, 1))
    flat_tube = BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeCylinder(flat_axis, 2.0, 10.0).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -6), gp_Dir(0, 0, 1)), 1.0, 12.0).Shape()
    ).Shape()
    # a box with a spherical scoop taken out of one corner: six planes and an EXTERIOR sphere
    flat_scooped = BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeBox(gp_Pnt(-3, -3, -3), 6.0, 6.0, 6.0).Shape(),
        BRepPrimAPI_MakeSphere(gp_Pnt(3, 3, 3), 2.5).Shape()).Shape()
    # A cylinder about (1, 1, 1): the ONLY fixture with off-diagonal quadric coefficients. Every
    # axis-aligned carrier leaves a01, a02 and a12 at zero, so a transposition in `_quadric`'s
    # slots -- or a mismatch against `EvalHalfspace`'s 2*(c1 xy + c2 xz + c4 yz) -- would pass a
    # suite built only from them and ship a wrong solid. A = I - dd^T here is -1/3 off-diagonal.
    flat_tilted = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(-1, -1, -1), gp_Dir(1, 1, 1)), 1.5, 6.0).Shape()
    flat_cases = (
        ("a box", BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0).Shape()),
        ("a tube, whose bore is an exterior cylinder", flat_tube),
        ("a cone frustum", BRepPrimAPI_MakeCone(flat_axis, 3.0, 1.0, 10.0).Shape()),
        ("a hemisphere", BRepAlgoAPI_Common(
            BRepPrimAPI_MakeSphere(gp_Pnt(1, 2, 3), 2.5).Shape(),
            BRepPrimAPI_MakeBox(gp_Pnt(-3, -1, 3), 9.0, 9.0, 9.0).Shape()).Shape()),
        ("a box with a spherical scoop, an exterior sphere", flat_scooped),
        ("a torus ply", BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                              4.0, 1.0).Shape()),
        ("a cylinder tilted about (1,1,1), whose quadric is dense", flat_tilted),
    )
    flat_results = [_flat_oracle(name, solid) for name, solid in flat_cases]
    for result in flat_results:
        bad = _flat_disagreements(result["blocks"], result["points"])
        check(f"the flat halfspaces of {result['name']} classify exactly as _cell_leaf's "
              "primitives",
              bad == 0 and len(result["points"]) > 0.5 * 4000,
              f"{bad} of {len(result['points'])} scored points disagree; carriers "
              f"{'+'.join(result['kinds'])} ({'+'.join(result['sides'])})")

    # All five carrier kinds must be covered, or the oracle comparison above proves less than it
    # reads: a kind nobody exercises is a kind whose sign nothing measured.
    flat_kinds_seen = sorted({k for r in flat_results for k in r["kinds"]})
    check("the flat oracle comparison covers all five carrier kinds",
          flat_kinds_seen == ["cone", "cylinder", "plane", "sphere", "torus"],
          f"covered {flat_kinds_seen}")
    check("the flat oracle comparison exercises a complemented (exterior) carrier",
          any("exterior" in r["sides"] for r in flat_results),
          "; ".join(f"{r['name']}: {'+'.join(r['sides'])}" for r in flat_results))
    # and a quadric with genuinely non-zero off-diagonal terms, per the note on `flat_tilted`
    flat_dense = [r["name"] for r in flat_results
                  if any(b["kind"] == "quadric" and max(abs(b["c"][1]), abs(b["c"][2]),
                                                        abs(b["c"][4])) > 1.0e-3
                         for b in r["blocks"])]
    check("the flat oracle comparison exercises off-diagonal quadric coefficients",
          bool(flat_dense), f"dense-quadric fixtures: {flat_dense}")

    # The negative control on the comparison itself. Zero disagreements would read the same way
    # if the sampler scored nothing that discriminates, so invert one halfspace at a time -- the
    # exact defect section 3.3 says is silent -- across EVERY fixture, and require every
    # inversion to be caught. This is what makes the checks above evidence rather than assertion.
    flat_missed = []
    for result in flat_results:
        for index in range(len(result["blocks"])):
            flipped = [dict(b, sign=-b["sign"]) if i == index else b
                       for i, b in enumerate(result["blocks"])]
            if _flat_disagreements(flipped, result["points"]) == 0:
                flat_missed.append(f"{result['name']}[{index}]")
    flat_flips = sum(len(r["blocks"]) for r in flat_results)
    check("inverting any one halfspace of any fixture is caught by the same comparison",
          flat_flips > 0 and not flat_missed,
          f"{flat_flips} inversion(s) over {len(flat_results)} fixtures, missed {flat_missed}")

    # --- the cone's mirror nappe, which the oracle sampler structurally cannot see -------------
    # `sign*Q <= 0` for a cone is `rho <= |r + k u|`: the DOUBLE cone. The mirror nappe lies
    # beyond the apex, which for an interior cone is outside the part bbox -- so neither the
    # sampler above nor `accept.contains_disagreements` (part bbox + 5%) ever looks at it. It is
    # not harmless there: `Contains_Loop` tests halfspaces with no cell-box clip while `Contains`
    # walks boxes built inside the declared box, so a nappe leaving that box makes O2FlatCSG's
    # own twins disagree. The assertion therefore has to be direct.
    flat_cone_result = next(r for r in flat_results if r["name"] == "a cone frustum")
    flat_cone_carrier = next(c for c in flat_cone_result["carriers"] if c["kind"] == "cone")
    flat_apex, flat_k = flatmod.cone_apex(flat_cone_carrier)
    flat_axis_d = flat_cone_carrier["d"]
    flat_ref = flat_cone_carrier["x"]

    def _flat_along_apex(steps, radial=0.0):
        """A point `steps` along the axis from the apex, positive being the material side.

        `rho <= r + k u` holds on the side where `k * (x - apex).d > 0`, so the material side is
        `+sign(k) d` and the mirror nappe is the other one, whatever the sign of `k`.
        """
        walk = math.copysign(1.0, flat_k) * steps
        return tuple(flat_apex[i] + walk * flat_axis_d[i] + radial * flat_ref[i]
                     for i in range(3))

    # a point strictly inside the mirror nappe: |r + k u| = |k| * 5 there, and the radius is half
    flat_mirror = _flat_along_apex(-5.0, radial=0.5 * abs(flat_k) * 5.0)
    flat_real = _flat_along_apex(5.0, radial=0.5 * abs(flat_k) * 5.0)
    # the emitter's contract for ONE cone carrier, isolated from the fixture's caps
    flat_cone_blocks = flatmod.blocks_from_carriers([flat_cone_carrier])
    flat_cone_quadric = [b for b in flat_cone_blocks
                         if not (b["kind"] == "quadric" and all(b["c"][i] == 0.0
                                                                for i in range(6)))]
    check("the cone quadric alone would admit a point on the mirror nappe",
          len(flat_cone_quadric) == 1 and flatmod.flat_contains(flat_cone_quadric, flat_mirror),
          f"the point {tuple(round(v, 6) for v in flat_mirror)} beyond the apex "
          f"{tuple(round(v, 6) for v in flat_apex)}")
    check("the emitted interior cone excludes the mirror nappe beyond its apex",
          len(flat_cone_blocks) == 2
          and not flatmod.flat_contains(flat_cone_blocks, flat_mirror),
          f"{len(flat_cone_blocks)} block(s) for one carrier, apex plane included")
    check("the apex plane cuts nothing on the cone's real nappe",
          flatmod.flat_contains(flat_cone_blocks, flat_real),
          f"the mirrored point {tuple(round(v, 6) for v in flat_real)} is still material")
    flat_apex_plane = flatmod.cone_apex_plane(flat_cone_carrier)
    check("the apex plane obeys the 2b = n convention like any other plane",
          flatmod.plane_scaling_error(flat_apex_plane) < 1.0e-15,
          f"residual {flatmod.plane_scaling_error(flat_apex_plane)}")

    # An EXTERIOR cone cannot be repaired that way -- the complement of one nappe is a union, so
    # a plane there would cut real material -- and must be DECLINED wherever the cell box reaches
    # past the apex. That needs the box, so it lives in `check_cell_box`, which the CONVERTER
    # owes per cell with the same lo/hi it passes to `SetCellBBox`.
    flat_exterior_cone = dict(flat_cone_carrier, side="exterior")
    check("an exterior cone is not silently given an apex plane",
          flatmod.cone_apex_plane(flat_exterior_cone) is None,
          "cone_apex_plane declines to repair a complemented cone")

    def _flat_cube_at(centre, half=0.5):
        return ([centre[i] - half for i in range(3)], [centre[i] + half for i in range(3)])

    flat_past_lo, flat_past_hi = _flat_cube_at(_flat_along_apex(-5.0))
    try:
        flatmod.check_cell_box([flat_exterior_cone], flat_past_lo, flat_past_hi)
        flat_box_reason = ""
    except recognise.Declined as why:
        flat_box_reason = str(why)
    check("an exterior cone whose cell box reaches past its apex is declined",
          "mirror nappe" in flat_box_reason, f"reason: {flat_box_reason or 'nothing raised'}")
    flat_short_lo, flat_short_hi = _flat_cube_at(_flat_along_apex(5.0))
    try:
        flatmod.check_cell_box([flat_exterior_cone], flat_short_lo, flat_short_hi)
        flat_stay_ok = True
    except recognise.Declined:
        flat_stay_ok = False
    check("an exterior cone whose cell box stays short of its apex is not declined",
          flat_stay_ok, f"box {tuple(round(v, 3) for v in flat_short_lo)} .. "
                        f"{tuple(round(v, 3) for v in flat_short_hi)}")
    # and an INTERIOR cone is never refused by that check, since its apex plane already fixed it
    try:
        flatmod.check_cell_box([flat_cone_carrier], flat_past_lo, flat_past_hi)
        flat_interior_ok = True
    except recognise.Declined:
        flat_interior_ok = False
    check("an interior cone is not refused for reaching past its apex",
          flat_interior_ok, "the apex plane already removed the mirror nappe")

    # The plane convention of design section 3.1, asserted where it is CREATED. `s.Q <= 0` names
    # the same halfspace under any positive rescaling, so no geometry test above would notice a
    # plane stored with |2b| != 1 -- what it costs is bit identity between O2FlatCSG's
    # accelerated queries and their _Loop twins. The C++ side cannot make this check, because it
    # cannot tell an emitter-produced halfspace from a hand-built one.
    flat_plane_worst = max((r["worstPlane"] for r in flat_results
                            if r["worstPlane"] is not None), default=None)
    check("every emitted plane block stores 2b = n for a unit normal",
          flat_plane_worst is not None and flat_plane_worst < 1.0e-15,
          f"worst | |2b| - 1 | over the fixtures: "
          f"{'no plane blocks' if flat_plane_worst is None else f'{flat_plane_worst:.3g}'}")
    # negative control on that check itself: a plane rescaled by 3 must be caught
    flat_tripled = {"kind": "quadric", "sign": 1.0,
                    "c": [0.0] * 6 + [1.5, 0.0, 0.0, -3.0, 0.0]}
    check("a plane block rescaled by three is refused by the convention check",
          abs(flatmod.plane_scaling_error(flat_tripled) - 2.0) < 1.0e-15,
          f"residual {flatmod.plane_scaling_error(flat_tripled)}")

    # A carrier kind with no quadric form declines rather than emitting a wrong halfspace.
    try:
        flatmod.quadric_from_carrier({"kind": "torus", "side": "interior"})
        flat_declined = ""
    except recognise.Declined as why:
        flat_declined = str(why)
    check("a carrier with no quadric form is declined, not guessed at",
          "no quadric form" in flat_declined, f"reason: {flat_declined or 'nothing raised'}")

    # A torus axis is normalised on the way into a block, because `AddTorus` normalises on load
    # and two evaluators reading one file must not differ.
    flat_long_axis = flatmod.blocks_from_carriers(
        [{"kind": "torus", "side": "interior", "p": (0.0, 0.0, 0.0), "d": (0.0, 0.0, 3.0),
          "r": 4.0, "rt": 1.0}])[0]
    check("a torus block's axis is a unit vector whatever the carrier carried",
          abs(math.sqrt(sum(flat_long_axis["c"][3 + i] ** 2 for i in range(3))) - 1.0) < 1.0e-15,
          f"axis {tuple(flat_long_axis['c'][3:6])}")

    # The sidecar writer's record sizes, against the format `LoadFlatCSG` checks before reading:
    # 20-byte header, 100-byte halfspace, 64-byte cell, little-endian, no padding.
    flat_sidecar = Path("/tmp/csg_selftest_flatcsg.bin")
    flat_probe_blocks = flat_results[1]["blocks"]
    flat_probe_cells = [{"first": 0, "count": len(flat_probe_blocks), "volume": 1.5,
                         "lo": [-2.0, -2.0, -5.0], "hi": [2.0, 2.0, 5.0]}]
    flatmod.write_sidecar(flat_sidecar, flat_probe_blocks, flat_probe_cells)
    flat_bytes = flat_sidecar.read_bytes()
    check("the sidecar is magic + version + two counts + fixed-length records",
          len(flat_bytes) == 20 + 100 * len(flat_probe_blocks) + 64 * len(flat_probe_cells)
          and flat_bytes[:8] == flatmod.SIDECAR_MAGIC
          and struct.unpack("<III", flat_bytes[8:20]) == (flatmod.SIDECAR_VERSION,
                                                          len(flat_probe_blocks),
                                                          len(flat_probe_cells)),
          f"{len(flat_bytes)} bytes for {len(flat_probe_blocks)} halfspace(s) and "
          f"{len(flat_probe_cells)} cell(s)")
    # and it round-trips through the format's own reader, field by field
    flat_read_back = []
    for index in range(len(flat_probe_blocks)):
        at = 20 + 100 * index
        kind = struct.unpack("<i", flat_bytes[at:at + 4])[0]
        sign = struct.unpack("<d", flat_bytes[at + 4:at + 12])[0]
        coeff = struct.unpack("<11d", flat_bytes[at + 12:at + 100])
        flat_read_back.append((kind, sign, coeff))
    check("every halfspace block reads back from the sidecar bit for bit",
          all(rb[0] == (1 if b["kind"] == "torus" else 0) and rb[1] == b["sign"]
              and list(rb[2]) == list(b["c"]) + [0.0] * (11 - len(b["c"]))
              for rb, b in zip(flat_read_back, flat_probe_blocks)),
          f"{len(flat_read_back)} block(s)")
    flat_cell_at = 20 + 100 * len(flat_probe_blocks)
    flat_cell_read = struct.unpack("<iid3d3d", flat_bytes[flat_cell_at:flat_cell_at + 64])
    check("the cell record reads back from the sidecar bit for bit",
          flat_cell_read == (0, len(flat_probe_blocks), 1.5,
                             -2.0, -2.0, -5.0, 2.0, 2.0, 5.0),
          f"{flat_cell_read}")

    # --- R5: the flat path takes what the TREE budget refuses, and nothing else --------------
    #
    # `Design_FlatCSGSolid.md` section 8. The two paths read the same decomposition and the same
    # cells; what separates them is the budget, and the budget is the whole point -- a 66-leaf
    # boolean tree costs 66 virtual calls per `Contains` and `O2FlatCSG` is not a tree.
    #
    # The brief for this task names `make_boolean_fixtures.many_celled_solid()`, which does not
    # exist (that module exposes `build_*` fixture entries, and none of them is over any budget).
    # `hex_collar` is used instead: it is the self-test's own two-cell, nine-leaf body from the
    # rung-4 block above, so the budgets below are stated against a shape whose cell and leaf
    # counts are already asserted there rather than against a new one.
    import itertools
    import random as flat_random
    tree_declined, tree_why = recognise.recognise_union_of_cells(hex_collar, max_leaves=8)
    check("a part over the tree's leaf budget still declines on the tree path",
          tree_declined is None and "part budget of 8" in (tree_why or ""),
          tree_why or "ACCEPTED")
    flat_record, flat_why = recognise.recognise_flat_cells(hex_collar)
    check("the same part is accepted on the flat path",
          flat_record is not None and flat_record["recogniser"] == "flat-cells", flat_why)
    check("the flat record carries a bounding box per cell",
          flat_record is not None
          and all(len(c["lo"]) == 3 and len(c["hi"]) == 3 for c in flat_record["cells"]),
          "a cell is missing its box")
    check("the flat description carries no leaves and no op on a cell",
          flat_record is not None and "leaves" not in flat_record
          and all(set(c) == set(prim.FLAT_CELL_KEYS) for c in flat_record["cells"]),
          f"{sorted(flat_record) if flat_record else None}")

    over_cells, over_why = recognise.recognise_flat_cells(hex_collar, max_cells=1)
    check("a part over the FLAT cell budget declines naming the bound",
          over_cells is None and "flat part budget of 1 cells" in (over_why or ""),
          over_why or "ACCEPTED")
    over_halfspaces, over_hs_why = recognise.recognise_flat_cells(hex_collar, max_halfspaces=3)
    check("a part over the FLAT halfspace budget declines naming the bound",
          over_halfspaces is None and "flat part budget of 3 halfspaces" in (over_hs_why or ""),
          over_hs_why or "ACCEPTED")

    # the false-accept guard must run here too: OCCT Cut can report IsDone with zero solids both
    # ways, which is how ST0923290_01#b19 got through the symmetric difference in rung 4
    check("the flat path runs the containment corroboration",
          "containsScored" in (flat_record or {}).get("notes", {})
          and (flat_record or {})["notes"].get("containsScored", 0) > 0,
          "no containment corroboration on the flat record")

    # ROUTING: the flat path runs only on what the tree path declines. `hex_collar` converts as
    # a union of cells today and must go on doing so -- if the flat branch were reached before
    # the union one, this is the check that would say so.
    routed, _routed_why = recognise.recognise(hex_collar)
    check("a part the tree path accepts is NOT intercepted by the flat path",
          routed is not None and routed["recogniser"] == "cells-union",
          routed["recogniser"] if routed else "declined")

    # --- R5: the multi-cell exterior cone, which is where the CELL box and the PART box differ --
    #
    # Task 8 left `flat.check_cell_box` an obligation on this task's call site and tested it only
    # on a single cell equal to the whole part, where the two boxes are the same box. This is the
    # interaction that was untested: an L whose lower arm carries a conical through-hole, and
    # whose cone apex sits ABOVE that arm -- outside the arm's own cell box, inside the part's.
    # Handing `check_cell_box` the part's box would refuse a part that is perfectly sound.
    l_cone_solid = BRepAlgoAPI_Cut(
        BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 10.0, 4.0, 4.0).Shape(),
                         BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 4), 4.0, 4.0, 6.0).Shape()).Shape(),
        BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(7, 2, 0), gp_Dir(0, 0, 1)),
                             1.5, 0.0, 8.0).Shape()).Shape()
    cone_record, cone_why = recognise.recognise_flat_cells(l_cone_solid)
    check("a multi-cell part with an exterior cone converts on the flat path",
          cone_record is not None and cone_record["notes"]["nCells"] == 2
          and cone_record["notes"]["cellGapCm"] <= 1.0e-9,
          cone_why or f"{cone_record['notes']['nCells']} cells, gap "
                      f"{cone_record['notes']['cellGapCm']:.3g} cm")

    # the discrimination itself, on the same carriers: the cell's own box accepts, the part's
    # box refuses. Without it the check above would pass on a part where the distinction is moot.
    l_cone_diag = recognise._bbox_diagonal(l_cone_solid)
    l_cone_tol = recognise.REL_TOL * max(l_cone_diag, 1.0)
    l_cone_report = decomp.split_into_cells(l_cone_solid, scale=max(l_cone_diag, 1.0))
    l_cone_part_box = recognise._bbox_of(l_cone_solid)
    cone_own, cone_part = [], []
    for piece in l_cone_report["pieces"]:
        _lv, piece_carriers, _out = recognise._cell_leaves(
            piece, l_cone_tol, decomp.bbox_diagonal(piece), whole_part=False)
        if not any(c["kind"] == "cone" and c["side"] == "exterior" for c in piece_carriers):
            continue
        piece_lo, piece_hi = recognise._flat_cell_box(
            piece, recognise._FLAT_BOX_MARGIN * max(l_cone_diag, 1.0))
        for label, box_lo, box_hi in (("own", piece_lo, piece_hi),
                                      ("part", list(l_cone_part_box[:3]),
                                       list(l_cone_part_box[3:]))):
            try:
                flatmod.check_cell_box(piece_carriers, box_lo, box_hi)
                (cone_own if label == "own" else cone_part).append("accepted")
            except recognise.Declined:
                (cone_own if label == "own" else cone_part).append("declined")
    check("the exterior cone is judged against its CELL's box, which the PART's box would fail",
          cone_own == ["accepted"] and cone_part == ["declined"],
          f"own box {cone_own}, part box {cone_part}")

    # and the call site really does hand `check_cell_box` the boxes it writes, per cell
    seen_boxes = []
    intact_check = flatmod.check_cell_box
    try:
        def _recording_check(carriers, lo, hi):
            seen_boxes.append(([float(v) for v in lo], [float(v) for v in hi]))
            return intact_check(carriers, lo, hi)
        flatmod.check_cell_box = _recording_check
        boxed_record, _boxed_why = recognise.recognise_flat_cells(l_cone_solid)
    finally:
        flatmod.check_cell_box = intact_check
    check("check_cell_box is called once per cell with exactly the box the sidecar carries",
          boxed_record is not None
          and len(seen_boxes) == len(boxed_record["cells"])
          and all(seen == ([float(v) for v in c["lo"]], [float(v) for v in c["hi"]])
                  for seen, c in zip(seen_boxes, boxed_record["cells"])),
          f"{len(seen_boxes)} call(s) for "
          f"{len(boxed_record['cells']) if boxed_record else '?'} cell(s)")

    # a Declined out of `check_cell_box` must reach the converter as a decline reason naming the
    # cell, not as an exception that kills the conversion
    try:
        def _refusing_check(carriers, lo, hi):
            raise recognise.Declined("a self-test refusal from check_cell_box")
        flatmod.check_cell_box = _refusing_check
        refused, refused_why = recognise.recognise_flat_cells(l_cone_solid)
    finally:
        flatmod.check_cell_box = intact_check
    check("a check_cell_box refusal becomes a decline naming the cell it came from",
          refused is None and "a self-test refusal from check_cell_box" in (refused_why or "")
          and "cell 1 of 2" in (refused_why or ""), refused_why or "ACCEPTED")

    # --- R5: the cell bounding box is an OUTER bound, checked rather than assumed --------------
    # `Design_FlatCSGSolid.md` section 4.2 makes this a correctness obligation on the converter:
    # O2FlatCSG builds no sub-cell box outside the declared one, so a cell reaching past its box
    # is material the accelerated `Contains` cannot find while `Contains_Loop` still can.
    box_escapes = []
    for record_label, record_cand in (("hex collar", flat_record), ("L with a cone", cone_record)):
        for index, c in enumerate(record_cand["cells"]):
            span = [c["hi"][i] - c["lo"][i] for i in range(3)]
            rng = flat_random.Random(90210 + index)
            for _ in range(3000):
                point = tuple(c["lo"][i] - span[i] + rng.random() * 3.0 * span[i]
                              for i in range(3))
                inside_box = all(c["lo"][i] <= point[i] <= c["hi"][i] for i in range(3))
                if not inside_box and flatmod.flat_contains(c["blocks"], point):
                    box_escapes.append(f"{record_label} cell {index}")
                    break
    check("no cell reaches outside the bounding box its record declares",
          not box_escapes, "; ".join(box_escapes) or "2 records, 4 cells, 12000 points sampled")
    escaped = None
    try:
        # one plane, `x <= 0`: an unbounded cell, and the box cannot hold it
        recognise._flat_box_holds_cell(
            [{"kind": "quadric", "sign": 1.0, "c": [0.0] * 6 + [0.5, 0.0, 0.0, 0.0] + [0.0]}],
            [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    except recognise.Declined as declined:
        escaped = str(declined)
    check("the outward probe catches a cell that is not closed up by its own halfspaces",
          escaped is not None and "past its own bounding box" in escaped,
          escaped or "ACCEPTED an unbounded cell")

    # --- R5: the sidecar the macro loads, and the shape the gate scores, are one solid ---------
    flat_blocks, flat_sidecar_cells = prim.flat_sidecar_records(cone_record)
    check("the sidecar's cell table indexes its concatenated halfspace blocks",
          len(flat_blocks) == cone_record["notes"]["nHalfspaces"]
          and [c["count"] for c in flat_sidecar_cells]
              == [len(c["blocks"]) for c in cone_record["cells"]]
          and [c["first"] for c in flat_sidecar_cells]
              == list(itertools.accumulate([0] + [len(c["blocks"])
                                                  for c in cone_record["cells"]][:-1]))
          and all(c["volume"] > 0.0 for c in flat_sidecar_cells),
          f"{len(flat_blocks)} block(s), {len(flat_sidecar_cells)} cell(s)")

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

        # --- the ROOT half of the single cell ---
        window_shape, window_placement = prim.build_root(window_record["candidate"],
                                                         "cellwindowprobe")
        node = window_shape.GetBoolNode()
        check("a single cell emits an unplaced TGeoCompositeShape over a TGeoSubtraction node",
              window_shape.ClassName() == "TGeoCompositeShape" and window_placement is None
              and node.ClassName() == "TGeoSubtraction"
              and node.GetLeftShape().ClassName() == "TGeoTube"
              and node.GetRightShape().ClassName() == "TGeoTube",
              f"{window_shape.ClassName()} over {node.ClassName()}"
              f"({node.GetLeftShape().ClassName()}, {node.GetRightShape().ClassName()}), "
              f"placement {'present' if window_placement else 'absent'}")
        steinmetz_shape, _pl = prim.build_root(steinmetz_record["candidate"], "cellsteinprobe")
        check("an intersection cell emits a TGeoIntersection node",
              steinmetz_shape.GetBoolNode().ClassName() == "TGeoIntersection",
              steinmetz_shape.GetBoolNode().ClassName())
        for label, record, solid_of in (("the window", window_record, window),
                                        ("the Steinmetz solid", steinmetz_record, steinmetz),
                                        ("the drilled cube", drilled_record, drilled)):
            cc = crosscheck_contains(record["candidate"], solid_of, n_points=20000)
            check(f"the ROOT cell and the CAD solid agree on Contains for {label}",
                  cc["disagreements"] == 0,
                  f"{cc['disagreements']} disagreement(s) over {cc['points']} points")
            dev = crosscheck_bbox(record["candidate"])
            check(f"the OCCT and ROOT realisations agree on the bounding box for {label}",
                  dev < 1.0e-9, f"max deviation {dev:.3g} cm")
        # 16/3 r^3 is the Steinmetz volume in closed form. `TGeoCompositeShape::Capacity()` is a
        # Monte-Carlo estimate, which is exactly why `checkKnownSource.py` marks a composite's
        # capacity not comparable -- the band here is that sampling noise, not a tolerance.
        want_steinmetz = 16.0 / 3.0
        rel_steinmetz = abs(steinmetz_shape.Capacity() - want_steinmetz) / want_steinmetz
        check("the emitted Steinmetz composite has the closed-form volume, to sampling noise",
              rel_steinmetz < 0.02,
              f"{steinmetz_shape.Capacity():.6f} vs {want_steinmetz:.6f} "
              f"(rel {rel_steinmetz:.2e}, Monte-Carlo)")
        cell_target = Path("/tmp/csg_selftest_cell.root")
        write_shape_root(window_record["candidate"], cell_target)
        fcell = ROOT.TFile.Open(str(cell_target))
        back_cell = fcell.Get("shape")
        check("shape_<part>.root round-trips a single cell as a TGeoCompositeShape",
              back_cell is not None and back_cell.ClassName() == "TGeoCompositeShape"
              and back_cell.GetBoolNode().ClassName() == "TGeoSubtraction",
              f"read {back_cell.ClassName() if back_cell else 'nothing'}")
        fcell.Close()

        # --- the ROOT half of the two new leaf types ---
        # The bounding box is asserted against the CLOSED FORM, not against OCCT's, for the two
        # tori: OCCT's `Bnd_Box` of a toroidal face is conservative -- it reports +-5.412 where
        # the body ends at +-5.000 -- while `TGeoTorus::ComputeBBox` is exact at (R + rmax,
        # R + rmax, rmax). Comparing the two would report 0.41 cm of OCCT's slack as if it were
        # this recogniser's error. The elliptic cylinder has no such slack and is compared both
        # ways.
        for label, record, solid_of, want_class, want_capacity, want_half in (
                ("the solid torus", solid_torus_record, solid_torus, "TGeoTorus",
                 2.0 * math.pi ** 2 * 4.0 * 1.0 ** 2, (5.0, 5.0, 1.0)),
                ("the torus shell", ply_record, ply, "TGeoTorus",
                 2.0 * math.pi ** 2 * 5.0 * (0.30 ** 2 - 0.28 ** 2), (5.3, 5.3, 0.3)),
                ("the elliptic cylinder", eltu_record, eltu_solid, "TGeoEltu",
                 math.pi * 3.0 * 1.5 * 10.0, (3.0, 1.5, 5.0))):
            shape, placement = prim.build_root(record["candidate"], f"probe_{want_class}")
            rel = abs(shape.Capacity() - want_capacity) / want_capacity
            check(f"{label} emits a bare {want_class} with the closed-form Capacity()",
                  shape.ClassName() == want_class and placement is None and rel < 1.0e-12,
                  f"{shape.ClassName()}, capacity {shape.Capacity():.9f} vs "
                  f"{want_capacity:.9f} (rel {rel:.2e}), placement "
                  f"{'present' if placement else 'absent'}")
            cc = crosscheck_contains(record["candidate"], solid_of, n_points=20000)
            check(f"the ROOT {want_class} and the CAD solid agree on Contains for {label}",
                  cc["disagreements"] == 0,
                  f"{cc['disagreements']} disagreement(s) over {cc['points']} points")
            half = (shape.GetDX(), shape.GetDY(), shape.GetDZ())
            worst = max(abs(h - w) for h, w in zip(half, want_half))
            check(f"the emitted {want_class}'s bounding box is the closed form for {label}",
                  worst < 1.0e-12,
                  f"{tuple(round(h, 9) for h in half)} vs {want_half}, worst {worst:.3g} cm")
        # A phi wedge is where a torus's own frame convention could be mirrored without any
        # volume noticing, so it gets the placement-composing check of its own.
        wedge_shape, wedge_placement = prim.build_root(wedge_torus_record["candidate"],
                                                       "probe_toruswedge")
        cc = crosscheck_contains(wedge_torus_record["candidate"], wedge_torus, n_points=20000)
        check("the ROOT TGeoTorus and the CAD solid agree on Contains for the hollow wedge",
              wedge_shape.ClassName() == "TGeoTorus" and cc["disagreements"] == 0,
              f"{wedge_shape.ClassName()}, {cc['disagreements']} disagreement(s) over "
              f"{cc['points']} points")
        placed_torus_shape, placed_torus_placement = prim.build_root(
            placed_torus_record["candidate"], "probe_placedtorus")
        cc = crosscheck_contains(placed_torus_record["candidate"], placed_torus, n_points=20000)
        check("a placed torus is a bare TGeoTorus plus a placement that composes correctly",
              placed_torus_shape.ClassName() == "TGeoTorus"
              and placed_torus_placement is not None and cc["disagreements"] == 0,
              f"{placed_torus_shape.ClassName()}, {cc['disagreements']} disagreement(s) over "
              f"{cc['points']} points")
        placed_eltu_shape, placed_eltu_placement = prim.build_root(
            placed_eltu_record["candidate"], "probe_placedeltu")
        cc = crosscheck_contains(placed_eltu_record["candidate"], placed_eltu, n_points=20000)
        check("a placed elliptic cylinder is a bare TGeoEltu plus a placement",
              placed_eltu_shape.ClassName() == "TGeoEltu"
              and placed_eltu_placement is not None and cc["disagreements"] == 0,
              f"{placed_eltu_shape.ClassName()}, {cc['disagreements']} disagreement(s) over "
              f"{cc['points']} points")
        for label, record in (("a TGeoTorus", solid_torus_record),
                              ("a TGeoEltu", eltu_record)):
            target = Path(f"/tmp/csg_selftest_{record['candidate']['leaves'][0]['type']}.root")
            write_shape_root(record["candidate"], target)
            handle = ROOT.TFile.Open(str(target))
            back = handle.Get("shape")
            check(f"shape_<part>.root round-trips {label}",
                  back is not None
                  and back.ClassName() == record["candidate"]["leaves"][0]["type"],
                  f"read {back.ClassName() if back else 'nothing'}")
            handle.Close()

        # --- the loaded shape reads the slots the writer meant -----------------------------
        # Byte parity with `WriteFlatCSG` says the two writers agree; it says nothing about
        # whether `LoadFlatCSG` and `O2FlatCSG::EvalHalfspace` INTERPRET slot 1 as a01 and slot 4
        # as a12. This closes that, and closes the duplication between `flat.eval_block` and
        # `EvalHalfspace` at the same time: a sidecar written by `flat.write_sidecar`, loaded and
        # closed in C++, must answer `Contains` exactly as `flat.flat_contains` does. The tilted
        # cylinder is in the loop on purpose -- it is the only fixture with off-diagonal terms,
        # so a transposition in either implementation shows up here.
        ROOT.gInterpreter.AddIncludePath(f"{ROOT.gSystem.Getenv('O2_ROOT')}/include")
        ROOT.gSystem.Load("libO2DetectorsBase")
        ROOT.gInterpreter.Declare(
            '#include "DetectorsBase/O2FlatCSG.h"\n'
            'namespace o2 { namespace base {\n'
            'bool LoadFlatCSG(const std::string& file, O2FlatCSG& solid);\n'
            '} }')
        flat_rt_bad = flat_rt_scored = 0
        flat_rt_failed = []
        for flat_index, result in enumerate(flat_results):
            blocks = result["blocks"]
            xlo, ylo, zlo, xhi, yhi, zhi = recognise._bbox_of(result["solid"])
            # the part bbox is an outer bound of this cell, because the cell IS the part here
            cells = [{"first": 0, "count": len(blocks), "volume": 1.0,
                      "lo": [xlo, ylo, zlo], "hi": [xhi, yhi, zhi]}]
            sidecar = Path(f"/tmp/csg_selftest_flatrt_{flat_index}.bin")
            flatmod.write_sidecar(sidecar, blocks, cells)
            loaded = ROOT.o2.base.O2FlatCSG(f"probe_flat_{flat_index}")
            if not ROOT.o2.base.LoadFlatCSG(str(sidecar), loaded):
                flat_rt_failed.append(f"{result['name']}: LoadFlatCSG refused the sidecar")
                continue
            loaded.CloseShape()
            if loaded.GetNhalfspaces() != len(blocks) or loaded.GetNcells() != 1:
                flat_rt_failed.append(f"{result['name']}: loaded "
                                      f"{loaded.GetNhalfspaces()}/{loaded.GetNcells()}")
                continue
            for point, _occ in result["points"]:
                flat_rt_scored += 1
                if bool(loaded.Contains(array("d", list(point)))) != \
                        flatmod.flat_contains(blocks, point):
                    flat_rt_bad += 1
        check("a sidecar written in Python and loaded in C++ answers Contains identically",
              not flat_rt_failed and flat_rt_bad == 0 and flat_rt_scored > 0,
              f"{flat_rt_bad} of {flat_rt_scored} points disagree over {len(flat_results)} "
              f"fixtures" + ("; " + "; ".join(flat_rt_failed) if flat_rt_failed else ""))

        # --- R5: the shipped shape of a multi-cell flat candidate ---------------------------
        # `prim.build_root` realises a `flatCells` description by writing its sidecar and loading
        # it back through `LoadFlatCSG`, so this exercises the exact artifact `geom.C` loads --
        # on a MULTI-CELL part, which the fixture round-trip above never was (every one of its
        # fixtures is a single cell equal to the whole part).
        flat_shape, flat_placement = prim.build_root(cone_record, "probe_flat_cells")
        check("a multi-cell flat candidate builds an O2FlatCSG through its own sidecar",
              flat_shape.ClassName() == "o2::base::O2FlatCSG" and flat_shape.IsClosed()
              and flat_shape.GetNcells() == len(cone_record["cells"])
              and flat_shape.GetNhalfspaces() == cone_record["notes"]["nHalfspaces"]
              and flat_placement is None,
              f"{flat_shape.GetNcells()} cell(s), {flat_shape.GetNhalfspaces()} halfspace(s), "
              f"{flat_shape.GetNboxes()} sub-cell box(es)")
        # the accelerated queries against the twin that defines them, and both against the
        # Python side that wrote the file: three implementations, one answer
        flat_rng = flat_random.Random(5150)
        blo = [min(c["lo"][i] for c in cone_record["cells"]) for i in range(3)]
        bhi = [max(c["hi"][i] for c in cone_record["cells"]) for i in range(3)]
        twin_bad = python_bad = 0
        for _ in range(20000):
            point = [blo[i] + flat_rng.random() * (bhi[i] - blo[i]) for i in range(3)]
            probe = array("d", point)
            accelerated = bool(flat_shape.Contains(probe))
            if accelerated != bool(flat_shape.Contains_Loop(probe)):
                twin_bad += 1
            if accelerated != any(flatmod.flat_contains(c["blocks"], tuple(point))
                                  for c in cone_record["cells"]):
                python_bad += 1
        check("the shipped flat shape agrees with its own _Loop twin and with csg/flat.py",
              twin_bad == 0 and python_bad == 0,
              f"{twin_bad} twin and {python_bad} emitter disagreement(s) over 20000 points")
        check("the flat shape's Capacity is the sum of its cells' own volumes",
              abs(flat_shape.Capacity()
                  - sum(c["volume"] for c in cone_record["cells"])) <= 1.0e-9,
              f"{flat_shape.Capacity():.9g} vs "
              f"{sum(c['volume'] for c in cone_record['cells']):.9g} cm^3")


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
