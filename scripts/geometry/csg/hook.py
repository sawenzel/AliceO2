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

from csg import emit, planar, primitives as prim  # noqa: E402


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

    Returns `(csg_files, flat_files, records)`: `csg_files` maps a logical-volume id to the
    absolute path of its `shape_*.root`, `flat_files` maps the ids that ship as
    `o2::base::O2FlatCSG` to their `flatcsg_*.bin` sidecar, and `records` is the per-part
    evidence for `csg_report.json`.

    A part is in exactly one of the two maps. Both artifacts are written for a flat part -- the
    sidecar is what `geom.C` loads (`Design_FlatCSGSolid.md` section 7) and the `shape_*.root` is
    what the oracle gate and `checkKnownSource.py` score -- but the macro has one emission branch
    per part, and the flat branch is the one that gets the `::Fatal` on a sidecar that does not
    load. The `.root` file is built by loading that same sidecar, so the two cannot describe
    different solids.
    """
    out_folder = Path(out_folder)
    root_available = have_root()
    csg_files = {}
    flat_files = {}
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
            is_flat = record["candidate"]["op"] == "flatCells"
            if root_available:
                # Built ONCE, measured, and only then written. Every check below runs before any
                # file exists, so a part the twin-parity gate refuses leaves neither a
                # `shape_<part>.root` nor a `flatcsg_<part>.bin` behind for a later reader to
                # find and wonder about.
                built = prim.build_root(record["candidate"], "shape")
                record["twinParity"] = emit.twin_parity(built[0])
                record["bboxRootVsOcctCm"] = emit.crosscheck_bbox(record["candidate"], built=built)
                record["containsCrosscheck"] = emit.crosscheck_contains(
                    record["candidate"], solid, built=built)
                # Two independent samplings of the same property, and either one refuses the
                # part. `twinParity` walks the shape's own declared cell boxes doubled, so 7/8 of
                # its points are OUTSIDE every box -- which is where a cell that reaches past its
                # box lives. `containsCrosscheck`'s twin column walks the CAD solid's bounding
                # box plus 5 %, which is where the part is. A disagreement in either means a cell
                # reaches past the bounding box the converter declared for it, the one defect
                # `SetCellBBox` cannot detect for itself (design section 4.2), so the part falls
                # a tier rather than ship a solid that is not the same solid as its own reference.
                parity = record["twinParity"]
                cross_twin = (record["containsCrosscheck"] or {}).get("twinDisagreements")
                if parity is not None and parity["disagreements"]:
                    record["accepted"] = False
                    record["reason"] = emit.twin_decline_reason(parity)
                elif cross_twin:
                    record["accepted"] = False
                    record["reason"] = emit.twin_decline_reason(
                        {"disagreements": cross_twin,
                         "points": record["containsCrosscheck"]["points"]})
                if not record["accepted"]:
                    record["shape"] = None
                    record["flatSidecar"] = None
                else:
                    if is_flat:
                        # Written inside this branch on purpose: a deferred part is not
                        # dispatched to CSG in `geom.C`, and a report row advertising a sidecar
                        # for a part that ships as a mesh would be a lie about what the macro
                        # builds. `emit.py --from-json` writes it when it completes the deferred
                        # shape, behind the same gate.
                        record["flatSidecar"] = write_flat_sidecar(
                            record["candidate"], out_folder, suffix)
                    target = (out_folder / f"shape_{suffix}.root").resolve()
                    emit.write_shape_object(built[0], built[1], target)
                    record["shape"] = str(target)
                    if is_flat:
                        flat_files[lid] = record["flatSidecar"]
                    else:
                        csg_files[lid] = str(target)
            else:
                record["shape"] = None
                record["shapeDeferred"] = True
                # NEXT.md item 2: this used to leave `reason` empty, so `csg_report.json` said
                # "declined CSG: None" for a part whose candidate was *accepted*. The reason must
                # name the real cause -- the environment, not the geometry.
                record["reason"] = ("csg deferred: ROOT unavailable in this interpreter; the "
                                    f"accepted candidate is in csg_{suffix}.json -- run "
                                    "csg/emit.py --from-json on the output folder to complete it")
        records.append(record)
        if verbose:
            emit._print_record(record)
            if record.get("shapeDeferred"):
                print(f"  [WARN] {display or lid}: accepted as CSG but NOT emitted -- "
                      "ROOT unavailable; geom.C will dispatch this part one tier down")

    n_csg = sum(1 for r in records if r["accepted"])
    if verbose:
        print(f"CSG recognition ({mode}): {n_csg}/{len(records)} leaf solid(s) accepted as native "
              f"ROOT shapes ({len(csg_files) + len(flat_files)} written, of which "
              f"{len(flat_files)} as flat halfspace solids)")
        if n_csg and not root_available:
            n_deferred = sum(1 for r in records if r.get("shapeDeferred"))
            print(f"  [WARN] PyROOT is not importable in this interpreter: {n_deferred} accepted "
                  "CSG part(s) were NOT emitted and geom.C dispatches them one tier down. "
                  "csg_report.json records each as 'csg deferred: ROOT unavailable'. Run "
                  "`csg/emit.py --from-json <output folder>` under the O2 environment, then "
                  "reconvert (or re-run the gate), to ship them as CSG.")
    if mode == "required":
        failed = [r for r in records if not r["accepted"]]
        if failed:
            lines = [f"--csg required: {len(failed)}/{len(records)} leaf solid(s) are not CSG:"]
            for r in failed:
                lines.append(f"  {r['volume'] or r['lid']}: {r['reason']}")
            raise ValueError("\n".join(lines))
    return csg_files, flat_files, records


def write_flat_sidecar(cand, out_folder, suffix):
    """Write `flatcsg_<part>.bin` for a `flatCells` candidate; returns its absolute path.

    The version-1 flat-CSG sidecar of `BVHSurfaceSolid.md` ("Flat-CSG sidecar format"), written
    by `csg/flat.py` and read by `o2::base::LoadFlatCSG`. Byte-compatible with `WriteFlatCSG`,
    which `csg/emit.py --self-test` pins.
    """
    from csg import flat
    target = (Path(out_folder) / f"flatcsg_{suffix}.bin").resolve()
    blocks, cells = prim.flat_sidecar_records(cand)
    flat.write_sidecar(target, blocks, cells)
    return str(target)


def write_report(records, path, surface_lids, facet_lids):
    """The tiered scorecard the project asked for: which representation carries each part, and why.

    A single coverage fraction cannot describe a system with three representations, so this
    reports one row per part with the representation that accepted it and the evidence that
    admitted it -- the symmetric-difference volume for CSG, the exact-surface extractor's own
    verdict for the surface solid, and nothing at all for the mesh, which is the point.

    Each row also carries `tessellationExact`: whether this part's *mesh* is the same solid as its
    exact surfaces rather than an approximation of them (`csg/planar.py`). That is a property of
    the part, not a choice between representations, and it is recorded rather than acted on --
    see the module docstring for what the measurements say the policy could be.

    `surface_lids` may be a set of lids or the lid -> sidecar-path mapping the converter builds;
    only the mapping lets the exactness be computed, and a bare set simply leaves it null.
    """
    surface_paths = surface_lids if isinstance(surface_lids, dict) else {}
    rows = []
    tiers = {"csg": 0, "surface": 0, "mesh": 0}
    exactness = {"exact": 0, "approximate": 0, "unknown": 0}
    for record in records:
        lid = record["lid"]
        if record["accepted"] and record.get("shape"):
            tier = "csg"
            why_not_csg = None
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
            why_not_csg = record["reason"]
            evidence = {"declinedCsgBecause": record["reason"]}
        else:
            tier = "mesh"
            why_not_csg = record["reason"]
            evidence = {"declinedCsgBecause": record["reason"]}
        tiers[tier] += 1
        sidecar = surface_paths.get(lid)
        if sidecar:
            mesh_exact, mesh_reason, mesh_census = planar.tessellation_is_exact(sidecar)
        else:
            mesh_exact, mesh_reason, mesh_census = None, "no exact sidecar for this part", None
        exactness["exact" if mesh_exact else
                  ("approximate" if mesh_exact is False else "unknown")] += 1
        # `part` is the artifact stem -- `<VOLNAME>_<LID>`, the same suffix that names
        # surfaces_/facets_/brep_/shape_ -- and it is written here so that a consumer can join
        # this row to `manifest.json` and to `gate.json` without re-implementing the converter's
        # filename sanitiser. `runOracleGate.py` reads the cascade decision from exactly this
        # field to decide which representation a part's verdict is computed on
        # (scripts/geometry/Stream_I_Verdict.md).
        rows.append({"lid": lid, "part": record.get("part"), "volume": record["volume"],
                     "representation": tier, "shapeFile": record.get("shape"),
                     # The `flatcsg_*.bin` a part carried by `o2::base::O2FlatCSG` ships with,
                     # null for every other part. `geom.C` loads this file; `shapeFile` above is
                     # the same solid, serialised, for the gate and the known-source check.
                     "flatSidecar": record.get("flatSidecar"),
                     # The brief, machine-readable decline reason (None for a part that ships as
                     # CSG). Same string as evidence.declinedCsgBecause, promoted to a top-level
                     # field so consumers (decline_catalogue.py, the website) need not know the
                     # evidence layout. A deferred part carries "csg deferred: ROOT unavailable".
                     "whyNotCSG": why_not_csg,
                     "shapeDeferred": bool(record.get("shapeDeferred", False)),
                     # 3x4 row-major [R | t] taking the shape's own frame into the part frame,
                     # or null for identity (Stream_N_PlacedPrimitives.md). Mirrors the
                     # TGeoHMatrix under the `placement` key of shape_<part>.root.
                     "shapePlacement": record.get("placement"),
                     # Is this part's TESSELLATION the same solid as its exact surfaces? True only
                     # when every face is a planar polygon; null when there is no sidecar to ask.
                     # Recorded, not acted on -- csg/planar.py says why it is worth having.
                     "tessellationExact": mesh_exact,
                     "tessellationExactWhy": mesh_reason,
                     "surfaceCensus": mesh_census,
                     "evidence": evidence})
    report = {"tiers": tiers, "tessellationExactness": exactness,
              "nLeafSolids": len(records), "parts": rows}
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
    exact = report.get("tessellationExactness") or {}
    if exact.get("exact"):
        total = sum(exact.values()) or 1
        print(f"  tessellation is EXACT (every face a planar polygon) for {exact['exact']} of "
              f"{total} part(s) -- {100.0 * exact['exact'] / total:.1f} %; for those the mesh is "
              f"not an approximation of the part, it is the part")
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


def emit_flat_csg_shape_cpp(lid, vol_display_name, sidecar_abspath, medium_var,
                            sanitize_cpp_name):
    """geom.C branch for a part carried by `o2::base::O2FlatCSG`: construct, load, close.

    Deliberately NOT the `LoadShape` branch above. The flat solid's sub-cell boxes and its BVH are
    transient (`Design_FlatCSGSolid.md` section 7), so a streamed shape has to be closed again by
    whoever reads it, and the sidecar is the form the design chose for the macro: `geom.C` stays
    small, Cling never sees a 59-cell part inlined, and the loader is the one the surface path has
    already proven.

    **A sidecar that fails to load is fatal.** `NEXT.md` item 2 is a live defect of exactly the
    opposite shape -- a JIT namespace bug let `geom.C` continue with a module silently absent and
    a simulation ran without it. A geometry that cannot be built must stop the job, not shrink.
    """
    safe = sanitize_cpp_name(lid)
    shape_name = vol_display_name if vol_display_name else lid
    return "\n".join([
        f'  auto *solid_{safe} = new o2::base::O2FlatCSG("{shape_name}");',
        f'  if (!o2::base::LoadFlatCSG("{sidecar_abspath}", *solid_{safe})) {{',
        f'    ::Fatal("geom", "flat-CSG sidecar for {shape_name} failed to load: '
        f'{sidecar_abspath}");',
        '  }',
        f'  solid_{safe}->CloseShape();',
        f'  if (!solid_{safe}->IsClosed()) {{',
        f'    ::Fatal("geom", "flat-CSG shape {shape_name} refused to close; see the Error above");',
        '  }',
        f'  TGeoVolume *vol_{safe} = new TGeoVolume("{shape_name}", solid_{safe}, {medium_var});',
    ])


FLAT_CPP_PRELUDE = r'''
// --- flat-CSG parts: o2::base::O2FlatCSG filled from a flatcsg_*.bin sidecar ---
// Both headers are included, never declared by prototype: loadCADGeometryHook JITs this
// macro inside a unique namespace and hoists only '#' lines to global scope, so a
// `namespace o2 { namespace base {` block here becomes `<wrapper>::o2::base` and shadows
// the real one -- every later o2::base:: name then fails to resolve and the module
// silently does not load. O2SurfaceSolidIO.h declares LoadFlatCSG and LoadSurfaceSolid both.
R__ADD_INCLUDE_PATH($O2_ROOT/include)
R__LOAD_LIBRARY(libO2DetectorsBase)
#include "DetectorsBase/O2FlatCSG.h"
#include "DetectorsBase/O2SurfaceSolidIO.h"
#include <TError.h>
'''


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
