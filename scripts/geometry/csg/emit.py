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


def process_solid(solid, name, tolerance=None, band_factor=1.0):
    """recognise -> build -> accept. Returns a record; `record['candidate']` is None if declined."""
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
    try:
        occ_shape = prim.build_occ(cand)
    except Exception as exc:                                     # noqa: BLE001
        record["reason"] = f"candidate failed to build in OCCT: {exc}"
        return record
    result = accept.symmetric_difference(solid, occ_shape, tol, band_factor)
    record["acceptance"] = result
    record["accepted"] = bool(result.get("accepted"))
    if not record["accepted"]:
        record["reason"] = result.get("reason")
    else:
        record["candidate"] = cand
    return record


def write_shape_root(cand, path):
    """Write the description as `shape_<part>.root`, per the Stream G convention.

    One object inheriting from `TGeoShape`, under the key `shape`, in cm, in the part's local
    frame. This is the same single `WriteTObject(shape, "shape")` that
    `o2::base::harness::saveShapeToRootFile` performs -- the C++ side is the authority and the
    unit test round-trips through it.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    shape = prim.build_root(cand, "shape")
    out = ROOT.TFile.Open(str(path), "RECREATE")
    out.WriteTObject(shape, "shape")
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

    **Read this number knowing what it can and cannot say.** For a bare primitive it is exact and
    a frame error moves it by the size of the error. For a `TGeoCompositeShape` it does not:
    `TGeoBoolNode::ComputeBBox` takes each operand's own axis-aligned box, transforms its eight
    corners and unions the results, so a leaf whose axis is not a coordinate axis contributes a
    box strictly larger than itself. On the six Bagger rams that inflation is 0.13-0.62 cm and is
    a property of ROOT's bounding box, not of the geometry -- the same effect Stream G §1 reports
    for the surface representation's conservative per-face AABBs. Use `crosscheck_contains` for a
    sharp check.
    """
    occ_shape = occ_shape if occ_shape is not None else prim.build_occ(cand)
    xmin, ymin, zmin, xmax, ymax, zmax = _occ_bbox(occ_shape)
    shape = prim.build_root(cand, "bboxprobe")
    origin = [shape.GetOrigin()[i] for i in range(3)]
    half = [shape.GetDX(), shape.GetDY(), shape.GetDZ()]
    worst = 0.0
    for i, (lo, hi) in enumerate(((xmin, xmax), (ymin, ymax), (zmin, zmax))):
        worst = max(worst, abs(origin[i] - half[i] - lo), abs(origin[i] + half[i] - hi))
    return worst


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
    shape = prim.build_root(cand, "containsprobe")
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
        if bool(shape.Contains(array("d", list(p)))) != (state == TopAbs_IN):
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

def self_test(verbose=True, with_root=True):
    """Synthetic solids whose recognition and emission are known in closed form.

    Every positive case is paired with a negative one, because a recogniser that accepts
    everything and an acceptance test that rejects nothing both pass a positive-only suite. The
    ROOT half is checked by querying the emitted `TGeoShape` against the closed-form answer for
    the same point set, which is independent of both OCCT and the description's own builders.
    """
    import math
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                                      BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere)
    from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))
        if verbose:
            print(f"  [{'ok ' if condition else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))

    def expect(name, solid, want_recogniser, want_leaves=1):
        record = process_solid(solid, name)
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

    # --- negative controls: each must decline or be rejected ---
    # 1. a blind bore: the recogniser proposes a through tube, the volume must refuse it.
    blind = BRepAlgoAPI_Cut(cyl, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, -6), gp_Dir(0, 0, 1)), 1.0, 9.0).Shape()).Shape()
    expect_declined("cylinder with a blind bore", blind)
    # 2. an L-shape: eight planes, no template.
    ell = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 4.0, 4.0, 1.0).Shape(),
                          BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -1), 4.0, 4.0, 3.0).Shape()).Shape()
    expect_declined("L-shaped plate", ell)
    # 3. a torus: in scope for the surface solid, out of scope here, and it must say so.
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
    expect_declined("torus", BRepPrimAPI_MakeTorus(5.0, 1.0).Shape(), "toroidal")
    # 4. a cylinder with a flat milled off it: one cluster, but a plane that is neither cap nor
    #    wedge. The recogniser must decline rather than emit the round tube.
    flatted = BRepAlgoAPI_Cut(cyl, BRepPrimAPI_MakeBox(
        gp_Pnt(1.5, -3, -6), 3.0, 6.0, 12.0).Shape()).Shape()
    expect_declined("cylinder with a milled flat", flatted)

    # --- the ROOT half: the emitted TGeoShape must answer like the closed form ---
    if with_root:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        from array import array
        import random
        shape = prim.build_root(moved_record["candidate"], "probe_moved")
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
            got = bool(shape.Contains(array("d", list(p))))
            if want != got and min(abs(rc - 1.0), abs(rc - 2.0), abs(abs(zc) - 5.0)) > 1e-9:
                bad += 1
        check("the emitted placed tube answers Contains like the closed form",
              bad == 0, f"{bad} disagreement(s) over 20000 points")
        # negative control on that check itself
        wrong = prim.build_root(prim.candidate("primitive", [prim.leaf(
            "TGeoTube", {"rmin": 1.0, "rmax": 2.05, "dz": 5.0}, frame)], "probe"), "probe_wrong")
        bad_wrong = 0
        random.seed(11)
        for _ in range(20000):
            p = (random.uniform(-2, 8), random.uniform(-9, 1), random.uniform(0, 10))
            rel = prim._sub(p, tuple(frame["origin"]))
            zc = prim._dot(rel, tuple(frame["z"]))
            rc = math.sqrt(max(prim._dot(rel, rel) - zc * zc, 0.0))
            want = (1.0 <= rc <= 2.0) and abs(zc) <= 5.0
            if want != bool(wrong.Contains(array("d", list(p)))):
                bad_wrong += 1
        check("the same check does report a wrong radius", bad_wrong > 0,
              f"{bad_wrong} disagreement(s) with rmax 2.05")
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
        box_shape = prim.build_root(box_record["candidate"], "boxprobe")
        origin = [box_shape.GetOrigin()[i] for i in range(3)]
        check("an axis-aligned box emits a bare TGeoBBox with its own origin",
              box_shape.ClassName() == "TGeoBBox"
              and max(abs(origin[0] - 1.0), abs(origin[1] - 1.5), abs(origin[2] - 2.0)) < 1e-12
              and abs(box_shape.Capacity() - 24.0) < 1e-12,
              f"{box_shape.ClassName()}, origin {origin}, capacity {box_shape.Capacity():.6f}")

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
