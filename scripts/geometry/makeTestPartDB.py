#!/usr/bin/env python3
"""
Build a test-part database for the solid-navigation harness (see
scripts/geometry/SolidNavigationHarness.md, Step 1).

For each input CAD model, runs O2_CADtoTGeo.py with `--exact-surfaces auto --mesh
--surface-report --dump-brep`, which writes (for every leaf logical volume that converts exactly)
a `surfaces_<VOL>_<LID>.bin` sidecar, a `facets_<VOL>_<LID>.bin` mesh and a
`brep_<VOL>_<LID>.brep` OCCT reference solid (all in cm) into the same output directory. This
script does not compute anything geometric itself: it only orchestrates the converter runs and
indexes the paired artifacts (a part enters the database only when both the sidecar and the mesh
exist; the BREP is indexed as `"brep"` when present and omitted otherwise, so a database built
before --dump-brep existed still loads) into a single `manifest.json` that the harness reads.

A fourth, optional artifact is indexed the same way: `shape_<VOL>_<LID>.root`, one ROOT-serialised
TGeoShape under the key `"shape"`, in cm and in the part's local frame. Nothing writes it today --
it is the hand-over format for the planned CSG emitter -- but the harness scores it side by side
with the other two representations whenever it is present.

Usage:
  python3 makeTestPartDB.py --output <db-dir>
  python3 makeTestPartDB.py --models Bagger.step as1-oc-214.stp --output <db-dir> --force

Requires the O2 + pythonOCC environment (same as O2_CADtoTGeo.py itself); see
SolidNavigationHarness.md, "Build and environment notes".
"""

import argparse
import datetime
import json
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_CONVERTER = _SCRIPT_DIR / "O2_CADtoTGeo.py"
_DEFAULT_MODEL_DIR = _SCRIPT_DIR / "STEP_examples"

# The MVP model set (see SolidNavigationHarness.md, Step 1): three small/fast models that
# together exercise plane/cylinder/sphere/torus analytic surfaces, recognized-NURBS-as-cylinder,
# and an all-planar part large enough to matter for scaling.
_DEFAULT_MODELS = [
    "Bagger.step",
    "as1-oc-214.stp",
    "oTOF System V3-R92cm.step",
]


def _sanitize_filename(s: str) -> str:
    """Mirror of O2_CADtoTGeo.py's sanitize_filename(); keep in sync with that copy."""
    safe = re.sub(r"[^0-9a-zA-Z]", "_", s)
    return safe or "x"


def _slugify_model(model_path: Path) -> str:
    return _sanitize_filename(model_path.stem)


def _resolve_model(model_arg: str) -> Path:
    p = Path(model_arg)
    if not p.is_absolute():
        candidate = _DEFAULT_MODEL_DIR / model_arg
        if candidate.exists():
            p = candidate
        else:
            p = Path(model_arg).expanduser().resolve()
    return p.resolve()


def _read_facets_summary(path: Path):
    """Return (nTriangles, bboxMin, bboxMax) by scanning a facets_*.bin file."""
    with open(path, "rb") as f:
        header = f.read(4)
        (n_tri,) = struct.unpack("<I", header)
        bbox_min = [float("inf")] * 3
        bbox_max = [float("-inf")] * 3
        rec_size = 9 * 4
        for _ in range(n_tri):
            data = f.read(rec_size)
            if len(data) != rec_size:
                raise ValueError(f"{path}: truncated facet record (expected {n_tri} triangles)")
            v = struct.unpack("<9f", data)
            for vi in range(3):
                x, y, z = v[3 * vi], v[3 * vi + 1], v[3 * vi + 2]
                bbox_min[0] = min(bbox_min[0], x)
                bbox_min[1] = min(bbox_min[1], y)
                bbox_min[2] = min(bbox_min[2], z)
                bbox_max[0] = max(bbox_max[0], x)
                bbox_max[1] = max(bbox_max[1], y)
                bbox_max[2] = max(bbox_max[2], z)
    if n_tri == 0:
        bbox_min = [0.0, 0.0, 0.0]
        bbox_max = [0.0, 0.0, 0.0]
    return n_tri, bbox_min, bbox_max


def _convert_model(model_path: Path, out_dir: Path, skip_existing: bool, force: bool):
    report_path = out_dir / "surface_report.json"

    if out_dir.exists() and any(out_dir.iterdir()):
        if force:
            shutil.rmtree(out_dir)
        elif skip_existing and report_path.exists():
            print(f"  [skip] {out_dir} already converted (--skip-existing)")
            return json.loads(report_path.read_text()), None
        else:
            raise RuntimeError(
                f"{out_dir} already exists and is non-empty. Pass --skip-existing to reuse it "
                "or --force to regenerate it."
            )

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(_CONVERTER),
        str(model_path),
        "--exact-surfaces", "auto",
        "--mesh",
        "--dump-brep",
        "--surface-report", str(report_path),
        "--output-folder", str(out_dir),
        "-o", "geom.C",
    ]
    print(f"  running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    if not report_path.exists():
        raise RuntimeError(f"{model_path}: converter did not produce {report_path}")
    return json.loads(report_path.read_text()), cmd


def _index_parts(model_name: str, slug: str, out_dir: Path, report: dict):
    """Pair surfaces_*.bin / facets_*.bin by <VOL>_<LID> suffix and read the report's
    (raw lid -> volume name) map to recover the manifest's `volume`/`lid` fields."""
    suffix_to_lid = {}
    for lid, info in report.get("volumes", {}).items():
        name = info.get("name") or ""
        volname = _sanitize_filename(name) if name else "vol"
        lidname = _sanitize_filename(lid)
        suffix_to_lid[f"{volname}_{lidname}"] = (lid, name)

    parts = []
    warnings = []
    for surf_path in sorted(out_dir.glob("surfaces_*.bin")):
        suffix = surf_path.name[len("surfaces_"):-len(".bin")]
        facet_path = out_dir / f"facets_{suffix}.bin"
        if not facet_path.exists():
            warnings.append(f"{surf_path.name}: no matching facets_{suffix}.bin, skipped")
            continue
        lid, volname = suffix_to_lid.get(suffix, (None, None))
        if lid is None:
            warnings.append(f"{surf_path.name}: suffix not found in surface_report.json volumes")
        n_tri, bbox_min, bbox_max = _read_facets_summary(facet_path)
        part = {
            "id": f"{slug}/{suffix}",
            "model": model_name,
            "volume": volname,
            "lid": lid,
            "surfaces": str(surf_path),
            "facets": str(facet_path),
            "nTriangles": n_tri,
            "bbox": {"min": bbox_min, "max": bbox_max},
        }
        # OCCT reference solid in cm (converter --dump-brep); absent for older databases.
        brep_path = out_dir / f"brep_{suffix}.brep"
        if brep_path.exists():
            part["brep"] = str(brep_path)
        # A third representation of the same part: one ROOT-serialised TGeoShape, under the key
        # "shape", in cm and in the part's own local frame -- the same frame as the sidecar, the
        # mesh and the .brep. This is the hand-over format for the CSG emitter (TGeoTube /
        # TGeoBBox / TGeoCompositeShape trees); nothing writes it today, so it is optional and a
        # part without one is indexed exactly as before. The harness also derives the same path
        # from the sidecar's, so dropping the file in after the fact and re-scoring with
        # --skip-convert works without re-indexing. See scripts/geometry/Stream_G_AnyShape.md.
        shape_path = out_dir / f"shape_{suffix}.root"
        if shape_path.exists():
            part["shape"] = str(shape_path)
        parts.append(part)
    return parts, warnings


def build_db(models, output: Path, skip_existing: bool, force: bool):
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        "version": 1,
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "output_dir": str(output.resolve()),
        "models": [],
        "parts": [],
    }

    for model_arg in models:
        model_path = _resolve_model(model_arg)
        if not model_path.exists():
            raise RuntimeError(f"Model not found: {model_arg} (resolved to {model_path})")
        slug = _slugify_model(model_path)
        out_dir = output / slug
        print(f"[{slug}] {model_path}")

        report, cmd = _convert_model(model_path, out_dir, skip_existing, force)
        parts, warnings = _index_parts(model_path.name, slug, out_dir, report)
        for w in warnings:
            print(f"  [warn] {w}")

        summary = report.get("summary", {})
        manifest["models"].append({
            "model": model_path.name,
            "model_path": str(model_path),
            "slug": slug,
            "output_dir": str(out_dir.resolve()),
            "command": cmd,
            "surface_report": str((out_dir / "surface_report.json").resolve()),
            "n_volumes": summary.get("n_volumes"),
            "n_eligible": summary.get("n_eligible"),
            "n_paired": len(parts),
            "warnings": warnings,
        })
        manifest["parts"].extend(parts)
        print(f"  -> {len(parts)} parts paired (of {summary.get('n_eligible')} exact-eligible / "
              f"{summary.get('n_volumes')} total volumes)")

    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=1))
    print(f"\nWrote {manifest_path} ({len(manifest['parts'])} parts across {len(manifest['models'])} models)")
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", nargs="+", default=_DEFAULT_MODELS,
                    help="CAD model files (relative names are resolved against "
                         f"{_DEFAULT_MODEL_DIR}). Default: the three-model MVP set.")
    ap.add_argument("--output", default=str(_SCRIPT_DIR / "test_part_db"),
                    help="Database output directory (default: %(default)s)")
    ap.add_argument("--skip-existing", action="store_true",
                    help="Reuse a model's output directory if already converted, re-indexing only.")
    ap.add_argument("--force", action="store_true",
                    help="Delete and regenerate a model's output directory even if it exists.")
    args = ap.parse_args()

    build_db(args.models, Path(args.output), args.skip_existing, args.force)


if __name__ == "__main__":
    main()
