"""The converter's single CSG integration point.

`O2_CADtoTGeo.py --csg auto|required` calls `recognise_and_emit()` once, with the leaf solids it
already has in memory, and gets back the map of logical volumes that are to be emitted as a native
ROOT shape plus a per-part record of what happened and on what evidence. Everything else lives in
this package; the converter gains one flag, one call and one emission branch.

The cascade
-----------
With `--csg auto` the converter's per-part choice becomes

    CSG  ->  exact surfaces (O2BVHSurfaceSolid)  ->  tessellated

and the reason for each choice is written to `csg_report.json` next to the other artifacts.

**The other representations are still written.** A part carried as CSG keeps its
`surfaces_*.bin`, `facets_*.bin` and `brep_*.brep` if those were asked for, because the gate
scores every representation of a part side by side (Stream G) and dropping them would silently
remove parts from the corpus that defends this project's zero-disagreement invariant. The cascade
decides what `geom.C` *builds*, which is the decision that matters in production.

Why the ROOT file may be written later
--------------------------------------
Recognition and acceptance need pythonOCC; writing the `TGeoShape` needs PyROOT. One interpreter
can have both (`csg/emit.py`), but `runOracleGate.py` launches the converter with
`PYTHONPATH`/`LD_LIBRARY_PATH` replaced by OCC-only values, and there `import ROOT` fails on
`libffi.so.6`. So this hook always writes the description as `csg_<VOL>_<LID>.json` and writes
`shape_<VOL>_<LID>.root` only when ROOT is importable; `csg/emit.py --from-json <dir>` completes
the job afterwards. A part whose `.root` file does not exist is *not* dispatched to CSG in
`geom.C` -- the converter never emits a reference to a file it did not write.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from csg import emit, primitives as prim  # noqa: E402


def have_root():
    try:
        import ROOT  # noqa: F401
        return True
    except Exception:                                            # noqa: BLE001
        return False


def _scaled_to_cm(shape, scale_to_cm):
    """Same transform `O2_CADtoTGeo.write_brep_cm()` applies, so the recogniser sees the frame
    and the units the sidecar, the mesh, the `.brep` and the oracle all use."""
    if scale_to_cm == 1.0:
        return shape
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.gp import gp_Pnt, gp_Trsf
    trsf = gp_Trsf()
    trsf.SetScale(gp_Pnt(0.0, 0.0, 0.0), scale_to_cm)
    return BRepBuilderAPI_Transform(shape, trsf, True).Shape()


def recognise_and_emit(def_shapes, def_names, scale_to_cm, out_folder, sanitize_filename,
                       mode="auto", band_factor=1.0, verbose=True):
    """Recognise every leaf solid; emit what both acceptance tests admit.

    Returns `(csg_files, records)`: `csg_files` maps a logical-volume id to the absolute path of
    its `shape_*.root`, and `records` is the per-part evidence for `csg_report.json`.
    """
    out_folder = Path(out_folder)
    root_available = have_root()
    csg_files = {}
    records = []
    for lid, shape in def_shapes.items():
        display = def_names.get(lid, "")
        volname = sanitize_filename(display) if display else "vol"
        suffix = f"{volname}_{sanitize_filename(lid)}"
        solid = _scaled_to_cm(shape, scale_to_cm)
        record = emit.process_solid(solid, suffix, band_factor=band_factor)
        record["lid"] = lid
        record["volume"] = display
        # The placement is derived from the description alone (no ROOT needed), so the deferred
        # `--from-json` path and this one cannot disagree about it. None means identity.
        record["placement"] = (prim.placement_for_candidate(record["candidate"])
                               if record["candidate"] else None)
        (out_folder / f"csg_{suffix}.json").write_text(json.dumps(
            {"part": suffix, "lid": lid, "candidate": record["candidate"],
             "acceptance": record["acceptance"], "recogniser": record["recogniser"],
             "placement": record["placement"]}, indent=1))
        if record["accepted"]:
            if root_available:
                target = (out_folder / f"shape_{suffix}.root").resolve()
                emit.write_shape_root(record["candidate"], target)
                record["shape"] = str(target)
                record["bboxRootVsOcctCm"] = emit.crosscheck_bbox(record["candidate"])
                record["containsCrosscheck"] = emit.crosscheck_contains(record["candidate"], solid)
                csg_files[lid] = str(target)
            else:
                record["shape"] = None
                record["shapeDeferred"] = True
        records.append(record)
        if verbose:
            emit._print_record(record)

    n_csg = sum(1 for r in records if r["accepted"])
    if verbose:
        print(f"CSG recognition ({mode}): {n_csg}/{len(records)} leaf solid(s) accepted as native "
              f"ROOT shapes ({len(csg_files)} written)")
        if n_csg and not root_available:
            print("  [warn] PyROOT is not importable in this interpreter, so no shape_*.root was "
                  "written; run `csg/emit.py --from-json <output folder>` under the O2 "
                  "environment to complete them. Those parts fall back one tier in geom.C.")
    if mode == "required":
        failed = [r for r in records if not r["accepted"]]
        if failed:
            lines = [f"--csg required: {len(failed)}/{len(records)} leaf solid(s) are not CSG:"]
            for r in failed:
                lines.append(f"  {r['volume'] or r['lid']}: {r['reason']}")
            raise ValueError("\n".join(lines))
    return csg_files, records


def write_report(records, path, surface_lids, facet_lids):
    """The tiered scorecard the project asked for: which representation carries each part, and why.

    A single coverage fraction cannot describe a system with three representations, so this
    reports one row per part with the representation that accepted it and the evidence that
    admitted it -- the symmetric-difference volume for CSG, the exact-surface extractor's own
    verdict for the surface solid, and nothing at all for the mesh, which is the point.
    """
    rows = []
    tiers = {"csg": 0, "surface": 0, "mesh": 0}
    for record in records:
        lid = record["lid"]
        if record["accepted"] and record.get("shape"):
            tier = "csg"
            evidence = {
                "recogniser": record["recogniser"],
                "description": record["description"],
                "symmetricDifferenceCm3": record["acceptance"]["symmetricDifference"],
                "bandCm3": record["acceptance"]["band"],
                "relativeToVolume": record["acceptance"]["relativeToVolume"],
                "rootVsCadContains": record.get("containsCrosscheck"),
            }
        elif lid in surface_lids:
            tier = "surface"
            evidence = {"declinedCsgBecause": record["reason"]}
        else:
            tier = "mesh"
            evidence = {"declinedCsgBecause": record["reason"]}
        tiers[tier] += 1
        # `part` is the artifact stem -- `<VOLNAME>_<LID>`, the same suffix that names
        # surfaces_/facets_/brep_/shape_ -- and it is written here so that a consumer can join
        # this row to `manifest.json` and to `gate.json` without re-implementing the converter's
        # filename sanitiser. `runOracleGate.py` reads the cascade decision from exactly this
        # field to decide which representation a part's verdict is computed on
        # (scripts/geometry/Stream_I_Verdict.md).
        rows.append({"lid": lid, "part": record.get("part"), "volume": record["volume"],
                     "representation": tier, "shapeFile": record.get("shape"),
                     "shapeDeferred": bool(record.get("shapeDeferred", False)),
                     # 3x4 row-major [R | t] taking the shape's own frame into the part frame,
                     # or null for identity (Stream_N_PlacedPrimitives.md). Mirrors the
                     # TGeoHMatrix under the `placement` key of shape_<part>.root.
                     "shapePlacement": record.get("placement"),
                     "evidence": evidence})
    report = {"tiers": tiers, "nLeafSolids": len(records), "parts": rows}
    Path(path).write_text(json.dumps(report, indent=1))
    return report


def print_tier_table(report):
    print("\n=== REPRESENTATION CASCADE (per leaf solid) ===")
    print(f"  {'volume':<28} {'carried by':<10} evidence")
    for row in report["parts"]:
        ev = row["evidence"]
        if row["representation"] == "csg":
            detail = (f"{ev['description']} [{ev['recogniser']}], dV_sym="
                      f"{ev['symmetricDifferenceCm3']:.3g} cm^3 (band {ev['bandCm3']:.3g})")
        else:
            detail = f"declined CSG: {ev['declinedCsgBecause']}"
        print(f"  {(row['volume'] or row['lid'])[:28]:<28} {row['representation']:<10} {detail}")
    tiers = report["tiers"]
    print(f"  tiers: CSG {tiers['csg']}, exact surfaces {tiers['surface']}, "
          f"tessellated {tiers['mesh']}  (of {report['nLeafSolids']} leaf solids)")


def csg_placement_var(lid, sanitize_cpp_name):
    """The macro variable holding a CSG part's shape placement. One namer, two call sites."""
    return f"shapePlace_{sanitize_cpp_name(lid)}"


def emit_csg_shape_cpp(lid, vol_display_name, shape_abspath, medium_var, sanitize_cpp_name):
    """geom.C branch for a CSG part: load the TGeoShape and its placement from its own file.

    Symmetrical with the sidecar and facet branches -- the macro reads its geometry from files
    next to itself rather than inlining thousands of numbers -- and it is the same file, in the
    same convention, that the oracle gate scored.

    Since `Stream_N_PlacedPrimitives.md` the shape may be expressed in its **own canonical
    frame** rather than the part's, with the rigid transform between them stored beside it. The transform is loaded here
    and composed at every `AddNode` of this volume as `partPlacement * shapePlacement`; it is
    never applied to the shape, because a `TGeoShape` cannot hold it -- that is the whole reason
    it travels separately. `LoadShapePlacement` returns the identity when the file records none,
    so a file written before Stream N (placed primitives) composes to exactly what it did before.
    """
    safe = sanitize_cpp_name(lid)
    shape_name = vol_display_name if vol_display_name else lid
    return "\n".join([
        f'  TGeoShape *solid_{safe} = LoadShape("{shape_abspath}", "{shape_name}");',
        f'  TGeoVolume *vol_{safe} = new TGeoVolume("{shape_name}", solid_{safe}, {medium_var});',
        f'  TGeoHMatrix *{csg_placement_var(lid, sanitize_cpp_name)} = '
        f'LoadShapePlacement("{shape_abspath}");',
    ])


def emit_csg_composed_placement_cpp(matrix_var, placement_var, composed_var):
    """`composed = partPlacement * shapePlacement`, in that order.

    Order, stated once because getting it wrong is silent: `matrix_var` maps the part frame into
    the parent's, `placement_var` maps the shape's own frame into the part's, and a point is
    carried shape -> part -> parent. `TGeoHMatrix::Multiply(right)` computes `this = this * right`,
    so the part placement is the one that is copied and the shape placement is the right operand.
    Reversing them, or transposing either rotation, moves a point measurably -- which is what
    `O2_CADtoTGeo.py --self-test` asserts.
    """
    return "\n".join([
        f"  TGeoHMatrix *{composed_var} = new TGeoHMatrix(*{matrix_var});",
        f"  {composed_var}->Multiply({placement_var});",
    ])


CPP_LOADER = r'''
// --- CSG parts: one ROOT-serialised TGeoShape per part, written by scripts/geometry/csg ---
// The file holds exactly one object inheriting from TGeoShape under the key "shape", in cm; and
// optionally a TGeoHMatrix under the key "placement", the rigid transform from the shape's own
// canonical frame into the part's local frame. No "placement" key means the identity, which is
// what every file written before that change means (scripts/geometry/Stream_G_AnyShape.md,
// scripts/geometry/Stream_N_PlacedPrimitives.md, and O2SolidHarness.h next to the C++ loader that
// reads the same convention).
TGeoHMatrix* LoadShapePlacement(const char* path) {
  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) {
    throw std::runtime_error(std::string("cannot open CSG shape file: ") + path);
  }
  auto* stored = dynamic_cast<TGeoHMatrix*>(f->Get("placement"));
  // Identity when the file records none. Returning a matrix rather than a null pointer keeps the
  // composition below unconditional, so the placed and unplaced cases go down one code path.
  auto* placement = stored ? new TGeoHMatrix(*stored) : new TGeoHMatrix("identity");
  f->Close();
  delete f;
  return placement;
}

TGeoShape* LoadShape(const char* path, const char* name) {
  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) {
    throw std::runtime_error(std::string("cannot open CSG shape file: ") + path);
  }
  auto* shape = dynamic_cast<TGeoShape*>(f->Get("shape"));
  if (!shape) {
    delete f;
    throw std::runtime_error(std::string("no TGeoShape under key \"shape\" in ") + path);
  }
  // The shape registers itself with gGeoManager on construction and is owned by it; the file can
  // go away.
  shape->SetName(name);
  f->Close();
  delete f;
  return shape;
}
'''
