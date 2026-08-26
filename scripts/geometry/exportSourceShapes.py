#!/usr/bin/env python3
"""Export the `TGeoShape` a round-tripped part was made FROM, as its own traceable artefact.

The converter's output holds what came *out* of the round trip -- the recognised composite, the
flat halfspace solid, the exact surface sidecar, the tessellation. It does not hold what went
*in*, and that is the one representation a reader actually starts from: the shape the detector
geometry ships today. Without it every comparison on the presentation website is a comparison
between four reconstructions, and the question a reader really asks -- is any of this better than
what we already have -- cannot be put on the page at all.

This writes `original_<stem>.root` next to the converted artefacts, one per part, carrying the
source volume's own shape under the key "shape", so `website/fetch_testdata.sh` picks it up like
any other artefact and `tgeoRayService.py` loads it like any other `.root`.

Finding the source volume is not a name lookup
----------------------------------------------
It is exactly the problem `checkKnownSource.py` already solves, and this script calls that module
rather than re-deriving it: the writer disambiguates duplicate names (`name#2`), a volume with
daughters is written as an assembly plus a `<name>__body` solid that is not a row of its own, and
several volumes can share one name and be told apart only by their bounding box. Every one of
those has already been paid for once.

What it refuses, and why
------------------------
A part is exported only when its source shape is in the SAME frame as the converted artefacts,
which means `shapePlacement` is the identity and the part is not a mirrored prototype. The website
traces every artefact in its own file's frame and has no per-artefact transform; shipping a shape
that needs one would put a correctly converted part next to a source shape drawn somewhere else
and call the difference a defect. A refused part is reported with its reason and simply has no
`original.root` -- a mirrored prototype loses nothing by it, because the volume it mirrors is
itself in the corpus.

Usage
-----
    exportSourceShapes.py --original   <corpus>/o2sim_geometry.root \\
                          --writer-report <corpus>/ITS_writer_report.json \\
                          --converted  <convdir> \\
                          [--parts IBCYSSFlangeC,BREF1] [--json original_report.json]

`--converted` is a directory holding the converter's `csg_report.json`; the files are written into
it. With no `--parts` every CSG-carried part of the report is exported.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import checkKnownSource as cks  # noqa: E402  (the path insert above is what makes it importable)


def boolean_shape_size(shape):
    """`(depth, leaves, {class: count})` of a TGeo boolean tree; a primitive is depth 0, 1 leaf.

    The depth is the recursion depth of `TGeoCompositeShape`'s own node tree, which is what a
    navigator descends on every query, and it is the number that says how a hand-written CSG
    expression is actually shaped -- a chain of 33 subtractions is depth 35, not depth 6.
    """
    from collections import Counter
    if not shape.InheritsFrom("TGeoCompositeShape"):
        return 0, 1, Counter([shape.ClassName()])
    node = shape.GetBoolNode()
    dl, nl, cl = boolean_shape_size(node.GetLeftShape())
    dr, nr, cr = boolean_shape_size(node.GetRightShape())
    return max(dl, dr) + 1, nl + nr, cl + cr


def describe(shape):
    depth, leaves, classes = boolean_shape_size(shape)
    return {
        "class": shape.ClassName(),
        "booleanDepth": depth,
        "leaves": leaves,
        "leafClasses": dict(sorted(classes.items())),
    }


def export_run(original, writer_report_path, converted, parts=None, verbose=True, write=True,
               tiers=("csg",)):
    """Write one `original_<stem>.root` per exportable part. Returns the per-part records.

    With `write=False` nothing is written and every part is still described -- the class, the
    boolean depth, the leaf count and the leaf classes of the shape it was made from. That is what
    a corpus-wide report wants: seventeen modules' worth of source shapes is a lot of ROOT files
    nobody asked for, and the description is the part of it a report reads.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libO2DetectorsBase")   # the emitted shape may be an O2 class

    converted = Path(converted)
    report_path = converted / "csg_report.json"
    if not report_path.exists():
        raise SystemExit(f"{report_path} does not exist (convert with --csg auto)")
    report = json.loads(report_path.read_text())
    writer_report = json.loads(Path(writer_report_path).read_text())
    index = cks._writer_index(writer_report)

    manager = ROOT.TGeoManager.Import(str(original))
    if manager is None:
        raise SystemExit(f"could not read a TGeoManager from {original}")
    by_name = {}
    for volume in manager.GetListOfVolumes():
        by_name.setdefault(volume.GetName(), []).append(volume)

    wanted = set(parts) if parts else None
    records = []
    for part in report.get("parts", []):
        # Writing is only meaningful for a part the cascade carried as CSG, because only then is
        # there an emitted shape to disambiguate against. A *report* wants every tier: the source
        # shape of a part that fell through to the mesh is exactly the interesting one.
        if part.get("representation") not in tiers:
            continue
        emitted_name = part.get("volume")
        stem = part.get("part")
        if wanted is not None and emitted_name not in wanted and stem not in wanted:
            continue
        record = {"part": stem, "volume": emitted_name, "written": None, "refused": None}
        row = index.get(emitted_name)
        if row is None:
            record["refused"] = f"no writer-report row for emittedName {emitted_name!r}"
            records.append(record)
            continue

        placement = part.get("shapePlacement")
        mirrored = emitted_name.endswith("__mirrored")
        # Both refusals below are about the FRAME, and a frame only matters when a file is being
        # written for something to trace. A description -- class, boolean depth, leaf classes --
        # is a property of the source shape and is the same whichever frame the part sits in, so
        # describe-only mode records the reason and carries on rather than dropping the part.
        frame_problem = None
        if mirrored:
            frame_problem = (
                f"a Z-mirrored prototype of {row.get('name')!r}: its source shape needs a "
                "reflection to reach this part's frame, and the volume it mirrors is in the "
                "corpus in its own right")
        elif not cks.placement_is_identity(placement):
            frame_problem = ("the emitted shape carries a non-identity placement, so the "
                             "source shape is not in the same frame as the other artefacts")
        if frame_problem:
            record["refused"] = frame_problem
            if write:
                records.append(record)
                continue

        candidates = by_name.get(row.get("name")) or []
        # The emitted shape tells apart several volumes that share a name, by its bounding box
        # (checkKnownSource.resolve_source_volume) -- and it is needed for nothing else. So the
        # file is opened only when there IS an ambiguity, and closed the moment resolution is
        # done. Holding one TFile per part open until the end costs a file descriptor each, and
        # MFT has 19 487 of them against a default limit of 1024.
        handle, emitted_shape = None, None
        if len(candidates) > 1:
            shape_file = part.get("shapeFile")
            if shape_file and Path(shape_file).exists():
                handle = ROOT.TFile.Open(str(shape_file))
                emitted_shape = handle.Get("shape") if handle else None
                if emitted_shape:
                    cks.reclose_flat_csg(emitted_shape)
        volume = (cks.resolve_source_volume(candidates, row, emitted_shape, placement)
                  if candidates else None)
        if handle:
            handle.Close()          # resolution is done; the emitted shape is not read again
        if volume is None:
            record["refused"] = (f"the original geometry has no volume named {row.get('name')!r} "
                                 "whose shape matches the writer's record")
            records.append(record)
            continue

        shape = volume.GetShape()
        record["sourceVolume"] = row.get("name")
        record.update(describe(shape))
        if write and not frame_problem:
            target = converted / f"original_{stem}.root"
            out = ROOT.TFile.Open(str(target), "RECREATE")
            out.WriteTObject(shape, "shape")
            out.Close()
            record["written"] = str(target)
        else:
            target = None
        records.append(record)
        if verbose and write:
            print(f"  {emitted_name:32s} <- {row.get('name'):28s} {record['class']:22s} "
                  f"depth {record['booleanDepth']:>3} / {record['leaves']:>3} leaves  -> "
                  f"{target.name}")

    if verbose:
        written = sum(1 for r in records if r["written"])
        print(f"{written}/{len(records)} part(s) exported"
              if write else f"{len(records)} part(s) described (nothing written)")
        for r in records:
            if r["refused"]:
                print(f"  refused {r['volume']}: {r['refused']}")
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original", required=True, help="the source o2sim_geometry.root")
    ap.add_argument("--writer-report", required=True, help="<MOD>_writer_report.json")
    ap.add_argument("--converted", required=True, help="the converter output directory")
    ap.add_argument("--parts", help="comma-separated part or volume names (default: all CSG parts)")
    ap.add_argument("--json", help="write the per-part records here")
    ap.add_argument("--no-write", action="store_true",
                    help="describe every part but write no original_*.root (for a corpus report)")
    ap.add_argument("--tiers", default="csg",
                    help="which cascade tiers to cover: csg,surface,mesh (default: csg)")
    args = ap.parse_args()

    parts = [p for p in (args.parts or "").split(",") if p] or None
    records = export_run(args.original, args.writer_report, args.converted, parts,
                         write=not args.no_write,
                         tiers=tuple(t for t in args.tiers.split(",") if t))
    if args.json:
        Path(args.json).write_text(json.dumps({"parts": records}, indent=1))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
