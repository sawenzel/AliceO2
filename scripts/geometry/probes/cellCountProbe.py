#!/usr/bin/env python3
"""The cell-count table Tier 3 was gated on: how many flat CSG cells a decomposition needs.

`CSG_Pipeline.md` §8 step 7 gates the flat `O2CSGSolid`/`TGeoBVHCSG` decision on "the cell-count
and success-rate table over the eligible ALICE3 and Bagger solids", and `Stream_A_CSG.md` §3
tightened the gate to *trusted* concave edges. This probe produces that table by actually running
the decomposition, not by estimating it:

    while a piece has a trusted concave (or mixed) edge:
        extend the carrier of one of the edge's two faces to a full surface;
        split the piece with BRepAlgoAPI_Splitter;
        recurse on the pieces.
    a piece with no trusted concave edge is one CSG cell
    (the intersection of its own faces' oriented halfspaces, Stream_A §2.1).

Per solid it reports: faces, distinct carrier halfspaces, axis clusters, trusted concave edges,
the number of cells reached, the number of splits and split failures, the per-cell halfspace
counts (the AND-lengths of the flat DNF), and whether the piece volumes still sum to the
original (a conservation guard on the boolean kernel). A budget (cells, splits, wall clock)
turns a blow-up into a reported lower bound instead of a hang.

This measures the *size* of a flat decomposition, not its exactness: the terminal pieces are
OCCT solids, and turning one into a sign-vector cell still requires the halfspace orientation
walk. The count is the number the cost/benefit argument needs.

Usage
-----
  cellCountProbe.py --self-test
  cellCountProbe.py --db <gate workdir>/db [--include NAME ...]      # brep_*.brep, in cm
  cellCountProbe.py --model .../CAD_noETA.stp [--include NAME ...]   # census prototypes
  ... --json out.json --markdown
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from csg.occ_env import ensure_occ  # noqa: E402

ensure_occ()

from csg import census  # noqa: E402
from csg.census import (NEAR_TANGENTIAL_SIN, FaceEdgeOrientations,  # noqa: E402
                        edge_dihedral, _carriers, carrier_clusters, distinct_carriers,
                        volume_of, solid_faces)


def bbox_diag(shape):
    box = census.bounding_box(shape)
    if box is None:
        return 1.0
    xmin, ymin, zmin, xmax, ymax, zmax = box
    return math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)

from OCC.Core.BRep import BRep_Tool  # noqa: E402
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Splitter  # noqa: E402
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace  # noqa: E402
from OCC.Core.Geom import (Geom_ConicalSurface, Geom_CylindricalSurface, Geom_Plane,  # noqa: E402
                           Geom_SphericalSurface, Geom_ToroidalSurface)
from OCC.Core.GeomAbs import (GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane,  # noqa: E402
                              GeomAbs_Sphere, GeomAbs_Torus)
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID  # noqa: E402
from OCC.Core.TopExp import TopExp_Explorer, topexp  # noqa: E402
from OCC.Core.TopTools import (TopTools_IndexedDataMapOfShapeListOfShape,  # noqa: E402
                               TopTools_ListOfShape)
from OCC.Core.TopoDS import topods  # noqa: E402


# ------------------------------------------------------------------------------------------
# finding the split witness: the first trusted concave/mixed edge, with its two faces
# ------------------------------------------------------------------------------------------

def first_trusted_concave_edge(solid):
    """(edge, face1, face2) of a trusted concave or mixed dihedral, or None.

    "Trusted" is `Stream_A_CSG.md` §1.4's filter: a concave/mixed verdict whose |n1 x n2| stays
    below NEAR_TANGENTIAL_SIN is a blend seam at the noise floor and no decomposition should
    split on it. Preference order: the sharpest (largest max-sin) witness, so the split happens
    where the reflex is best conditioned.
    """
    amap = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, amap)
    orients = FaceEdgeOrientations(solid)
    best = None
    for i in range(1, amap.Size() + 1):
        edge = topods.Edge(amap.FindKey(i))
        if BRep_Tool.Degenerated(edge):
            continue
        faces = list(amap.FindFromIndex(i))
        distinct = []
        for f in faces:
            if not any(f.IsSame(g) for g in distinct):
                distinct.append(f)
        if len(distinct) != 2:
            continue
        f1, f2 = topods.Face(distinct[0]), topods.Face(distinct[1])
        verdict, max_sin = edge_dihedral(edge, f1, f2, orients)
        if verdict in ("concave", "mixed") and max_sin >= NEAR_TANGENTIAL_SIN:
            if best is None or max_sin > best[3]:
                best = (edge, f1, f2, max_sin)
    return None if best is None else best[:3]


def count_trusted_concave(solid):
    counts = census.edge_census(solid)
    return (counts["concave"] + counts["mixed"]
            - counts["concaveNearTangential"] - counts["mixedNearTangential"])


# ------------------------------------------------------------------------------------------
# extending a face's carrier into a splitting tool
# ------------------------------------------------------------------------------------------

def carrier_tool_face(face, extent):
    """A face covering the whole carrier of `face`, big enough to cut anything within `extent`.

    For a plane and the lateral quadrics the parametric window is widened to `extent` (the
    caller passes a few bounding-box diagonals); a sphere and a torus are already closed.
    Returns None when the carrier type is not one of the five families.
    """
    ad = BRepAdaptor_Surface(face, True)
    t = ad.GetType()
    two_pi = 2.0 * math.pi
    try:
        if t == GeomAbs_Plane:
            surf = Geom_Plane(ad.Plane())
            return BRepBuilderAPI_MakeFace(surf, -extent, extent, -extent, extent,
                                           1.0e-7).Face()
        if t == GeomAbs_Cylinder:
            surf = Geom_CylindricalSurface(ad.Cylinder())
            return BRepBuilderAPI_MakeFace(surf, 0.0, two_pi, -extent, extent, 1.0e-7).Face()
        if t == GeomAbs_Cone:
            surf = Geom_ConicalSurface(ad.Cone())
            return BRepBuilderAPI_MakeFace(surf, 0.0, two_pi, -extent, extent, 1.0e-7).Face()
        if t == GeomAbs_Sphere:
            surf = Geom_SphericalSurface(ad.Sphere())
            return BRepBuilderAPI_MakeFace(surf, 0.0, two_pi, -0.5 * math.pi, 0.5 * math.pi,
                                           1.0e-7).Face()
        if t == GeomAbs_Torus:
            surf = Geom_ToroidalSurface(ad.Torus())
            return BRepBuilderAPI_MakeFace(surf, 0.0, two_pi, 0.0, two_pi, 1.0e-7).Face()
    except Exception:
        return None
    return None


def split_solid(piece, tool):
    """Split `piece` by `tool`; returns the list of solids, or None on failure."""
    splitter = BRepAlgoAPI_Splitter()
    args = TopTools_ListOfShape()
    args.Append(piece)
    tools = TopTools_ListOfShape()
    tools.Append(tool)
    splitter.SetArguments(args)
    splitter.SetTools(tools)
    try:
        splitter.Build()
    except Exception:
        return None
    if not splitter.IsDone():
        return None
    out = []
    exp = TopExp_Explorer(splitter.Shape(), TopAbs_SOLID)
    while exp.More():
        out.append(topods.Solid(exp.Current()))
        exp.Next()
    return out or None


# ------------------------------------------------------------------------------------------
# the decomposition loop
# ------------------------------------------------------------------------------------------

def decompose(solid, max_cells=128, max_splits=256, timeout_s=300.0, verbose=False):
    """Split at trusted concave edges until every piece is one cell, or the budget runs out."""
    t0 = time.time()
    diag = bbox_diag(solid)
    extent = 4.0 * max(diag, 1.0)
    vol0 = volume_of(solid)
    pending = [solid]
    cells = []
    unresolved = []          # pieces whose witness could not be split (both tools failed)
    n_splits = 0
    n_split_failures = 0
    budget = None
    while pending:
        if len(cells) + len(pending) + len(unresolved) > max_cells:
            budget = f"cell budget {max_cells} exceeded"
            break
        if n_splits >= max_splits:
            budget = f"split budget {max_splits} exceeded"
            break
        if time.time() - t0 > timeout_s:
            budget = f"timeout {timeout_s:.0f} s exceeded"
            break
        piece = pending.pop()
        witness = first_trusted_concave_edge(piece)
        if witness is None:
            cells.append(piece)
            continue
        _edge, f1, f2 = witness
        # Prefer the planar carrier as the knife: OCCT's boolean core is most robust on planes.
        faces = sorted((f1, f2), key=lambda f: BRepAdaptor_Surface(f, True).GetType()
                       != GeomAbs_Plane)
        parts = None
        for face in faces:
            tool = carrier_tool_face(face, extent)
            if tool is None:
                continue
            parts = split_solid(piece, tool)
            if parts is not None and len(parts) > 1:
                break
            parts = None
        if parts is None:
            n_split_failures += 1
            unresolved.append(piece)
            continue
        n_splits += 1
        pending.extend(parts)
        if verbose:
            print(f"    split {n_splits}: {len(parts)} piece(s), "
                  f"{len(cells)} cell(s) so far, {len(pending)} pending")

    per_cell_halfspaces = []
    residual_untrusted = 0
    for cell in cells:
        faces = solid_faces(cell)
        carriers = _carriers(faces)
        n_h = distinct_carriers(carriers, max(diag, 1.0)) if carriers else len(faces)
        per_cell_halfspaces.append(n_h if n_h is not None else len(faces))
        counts = census.edge_census(cell)
        residual_untrusted += (counts["concaveNearTangential"] + counts["mixedNearTangential"])

    vol_cells = sum(volume_of(c) for c in cells) + sum(volume_of(u) for u in unresolved)
    conserved = (abs(vol_cells - vol0) <= 1.0e-6 * max(abs(vol0), 1.0)) if budget is None else None
    complete = budget is None and not unresolved
    return {
        "cells": len(cells) + len(unresolved) + len(pending),
        "cellsClean": len(cells),
        "unresolvedPieces": len(unresolved),
        "pendingAtBudget": len(pending),
        "splits": n_splits,
        "splitFailures": n_split_failures,
        "budgetStop": budget,
        "complete": complete,
        "volumeConserved": conserved,
        "volumeOriginal": vol0,
        "volumePieces": vol_cells,
        "perCellHalfspaces": sorted(per_cell_halfspaces, reverse=True),
        "sumCellHalfspaces": sum(per_cell_halfspaces),
        "residualUntrustedConcave": residual_untrusted,
        "seconds": round(time.time() - t0, 2),
    }


def probe_solid(solid, name, **kwargs):
    faces = solid_faces(solid)
    carriers = _carriers(faces)
    diag = bbox_diag(solid)
    clusters = carrier_clusters(carriers, max(diag, 1.0)) if carriers else None
    row = {
        "name": name,
        "faces": len(faces),
        "quadricOnly": carriers is not None,
        "halfspaces": distinct_carriers(carriers, max(diag, 1.0)) if carriers else None,
        "axisClusters": clusters["nAxisClusters"] if clusters else None,
        "planeDirections": clusters["nPlaneDirections"] if clusters else None,
        "concaveTrusted": count_trusted_concave(solid),
    }
    row.update(decompose(solid, **kwargs))
    return row


# ------------------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------------------

def solids_from_db(db_dir, include):
    from csg.emit import load_brep
    breps = sorted(Path(db_dir).glob("*/brep_*.brep")) or sorted(Path(db_dir).glob("brep_*.brep"))
    for brep in breps:
        name = brep.name[len("brep_"):-len(".brep")]
        if include and not any(pat in name for pat in include):
            continue
        solid, _n = load_brep(brep)
        yield name, solid


def solids_from_model(path, include):
    """The model's prototype solids (one per IsPartner class), as the census counts them."""
    protos = []
    for name, solid in census.load_step_solids(str(path)):
        if any(solid.IsPartner(p) for _n, p in protos):
            continue
        protos.append((name, solid))
    for name, solid in protos:
        if include and not any(pat in name for pat in include):
            continue
        yield name, solid


# ------------------------------------------------------------------------------------------
# self-test: solids whose cell count is known in closed form
# ------------------------------------------------------------------------------------------

def self_test(verbose=True):
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt

    checks = []

    def check(name, cond, detail=""):
        checks.append((name, bool(cond)))
        if verbose:
            print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))

    def cyl(r, h, o=(0.0, 0.0, 0.0), d=(0.0, 0.0, 1.0)):
        return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*o), gp_Dir(*d)), r, h).Shape()

    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 4.0, 3.0, 2.0).Shape()
    r = probe_solid(box, "box")
    check("a box is one cell, zero splits", r["cells"] == 1 and r["splits"] == 0
          and r["complete"] and r["volumeConserved"], json.dumps(r["perCellHalfspaces"]))

    plate = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 6.0, 6.0, 1.0).Shape(),
                            cyl(1.0, 3.0, (3.0, 3.0, -1.0))).Shape()
    r = probe_solid(plate, "through-hole plate")
    check("a through-hole plate is ONE cell (a hole adds no concave edge)",
          r["cells"] == 1 and r["splits"] == 0 and r["complete"], f"cells={r['cells']}")

    ell = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 4.0, 4.0, 1.0).Shape(),
                          BRepPrimAPI_MakeBox(gp_Pnt(2, 2, -1), 4.0, 4.0, 3.0).Shape()).Shape()
    r = probe_solid(ell, "L-plate")
    check("an L-plate needs exactly one split into two cells",
          r["cells"] == 2 and r["splits"] == 1 and r["complete"] and r["volumeConserved"],
          f"cells={r['cells']} splits={r['splits']} conserved={r['volumeConserved']}")

    blind = BRepAlgoAPI_Cut(cyl(2.0, 10.0), cyl(1.0, 9.0, (0.0, 0.0, 1.0))).Shape()
    r = probe_solid(blind, "blind bore")
    check("a blind bore splits at the bore bottom into two cells",
          r["cells"] == 2 and r["complete"] and r["volumeConserved"], f"cells={r['cells']}")

    groove = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 6.0, 4.0, 3.0).Shape(),
                             BRepPrimAPI_MakeBox(gp_Pnt(2, -1, 1), 2.0, 6.0, 3.0).Shape()).Shape()
    r = probe_solid(groove, "grooved block")
    check("a grooved block decomposes into three cells",
          r["cells"] == 3 and r["complete"] and r["volumeConserved"], f"cells={r['cells']}")

    steinmetz = BRepAlgoAPI_Common(cyl(1.0, 6.0, (0.0, 0.0, -3.0)),
                                   cyl(1.0, 6.0, (0.0, -3.0, 0.0), (0.0, 1.0, 0.0))).Shape()
    r = probe_solid(steinmetz, "Steinmetz")
    check("the Steinmetz intersection is one cell of two halfspaces",
          r["cells"] == 1 and r["complete"] and r["perCellHalfspaces"] == [2],
          f"cells={r['cells']} halfspaces={r['perCellHalfspaces']}")

    tiny = probe_solid(ell, "budget", max_splits=0)
    check("an exhausted budget is reported as a lower bound, not an answer",
          tiny["budgetStop"] is not None and not tiny["complete"], tiny["budgetStop"] or "")

    n_ok = sum(1 for _n, ok in checks if ok)
    print(f"  {n_ok}/{len(checks)} cell-count probe self-checks passed")
    return n_ok, len(checks)


# ------------------------------------------------------------------------------------------

def markdown(rows):
    out = ["| part | faces | halfsp | clusters | concave (trusted) | cells | splits | fail | "
           "sum cell-halfsp | max cell | complete | vol ok | s |",
           "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: "
           "| :---: | ---: |"]
    for r in rows:
        cells = str(r["cells"]) if r["complete"] else f">={r['cells']}"
        top = r["perCellHalfspaces"][0] if r["perCellHalfspaces"] else "-"
        out.append(f"| `{r['name']}` | {r['faces']} | {r['halfspaces'] or '-'} "
                   f"| {r['axisClusters'] if r['axisClusters'] is not None else '-'} "
                   f"| {r['concaveTrusted']} | **{cells}** | {r['splits']} "
                   f"| {r['splitFailures']} | {r['sumCellHalfspaces'] or '-'} | {top} "
                   f"| {'Y' if r['complete'] else r['budgetStop'] or 'n'} "
                   f"| {'Y' if r['volumeConserved'] else '.'} | {r['seconds']} |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, help="gate workdir db/ (walks brep_*.brep, in cm)")
    ap.add_argument("--model", type=Path, help="a STEP model (census prototypes)")
    ap.add_argument("--include", nargs="*", default=[],
                    help="substring filter(s) on part names")
    ap.add_argument("--max-cells", type=int, default=128)
    ap.add_argument("--max-splits", type=int, default=256)
    ap.add_argument("--timeout", type=float, default=300.0, help="per-solid seconds")
    ap.add_argument("--json", type=Path, help="write the table as JSON")
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        ok, n = self_test()
        return 0 if ok == n else 1

    if not args.db and not args.model:
        ap.error("give --db or --model (or --self-test)")

    source = (solids_from_db(args.db, args.include) if args.db
              else solids_from_model(args.model, args.include))
    rows = []
    for name, solid in source:
        print(f"[probe] {name} ...", flush=True)
        row = probe_solid(solid, name, max_cells=args.max_cells, max_splits=args.max_splits,
                          timeout_s=args.timeout, verbose=args.verbose)
        rows.append(row)
        cells = row["cells"] if row["complete"] else f">={row['cells']} ({row['budgetStop']})"
        print(f"        faces={row['faces']} halfspaces={row['halfspaces']} "
              f"concaveTrusted={row['concaveTrusted']} -> cells={cells} "
              f"splits={row['splits']} failures={row['splitFailures']} "
              f"in {row['seconds']} s", flush=True)

    if args.json:
        args.json.write_text(json.dumps(rows, indent=1))
        print(f"Wrote {args.json}")
    if args.markdown:
        print(markdown(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
