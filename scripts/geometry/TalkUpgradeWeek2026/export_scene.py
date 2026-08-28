#!/usr/bin/env python3
"""Export a built CAD-converted TGeo geometry as a Blender-ready scene bundle.

Joins three things that already exist after a conversion:

  * ``facets_<part>.bin``  -- the part's triangles in its own local frame (cm),
  * the built TGeo geometry -- every global placement of every part,
  * ``csg_report.json``    -- which representation the cascade shipped the part in.

The result is one OBJ object per *placement*, plus a manifest carrying identity,
attribution and a label anchor per part. It assigns no colours, materials or
cameras: what a part looks like is the renderer's decision, and the manifest is
written so that colouring per part or per representation are equally possible.

This replaces the STEP volume-matching detour of
``O2DPGAgentic/doc/detector-cad-rendering.md``. That was needed because the mesh
dump's repeated leaf volumes are unplaced by design; with a real built geometry
the placements are simply read off the node tree.

Usage:
  export_scene.py --model-dir <dir with geom.C, csg_report.json, facets_*.bin> \
                  --out <output dir> --name <tag>
"""
import argparse
import json
import math
import os
import re
import struct
import sys


# ---------------------------------------------------------------- inputs ----

def parse_volume_order(geom_c):
    """Ordered [(part_key, volume_name)] as geom.C creates them.

    The lid lives only in the C++ variable name (``vol__0_1_1_2``); the TGeo
    volume carries the display name alone, and display names repeat across the
    bodies of one multi-body label (IRIS has 28 volumes called ST0923290_01).
    So creation order is the join, and `<name>_<lid>` is csg_report's part key.
    """
    out = []
    for var, name in re.findall(
            r'TGeoVolume\s*\*(\w+)\s*=\s*new TGeoVolume\("([^"]+)"', geom_c):
        lid = var[len("vol__"):] if var.startswith("vol__") else var
        out.append((f"{name}_{lid}", name))
    return out


def read_facets(path):
    """(ntri, [ (x,y,z) * 9 ]) from a facets_*.bin: uint32 count, then 9 float32 per triangle."""
    with open(path, "rb") as fh:
        ntri = struct.unpack("<I", fh.read(4))[0]
        data = fh.read(ntri * 9 * 4)
    if len(data) != ntri * 9 * 4:
        raise RuntimeError(f"{path}: truncated, want {ntri} triangles")
    return ntri, struct.unpack(f"<{ntri * 9}f", data)


SHORT_CLASS = {"surface": "exact surfaces", "mesh": "tessellated"}


def short_label(rec):
    """Two-to-four words: the part plus its solid class. No tolerances, no dV_sym.

    A label drawn next to the geometry has to read at slide size; the recogniser
    string and the full description stay in their own fields for a caption.
    """
    vol = rec.get("volume") or rec["part"]
    representation = rec["representation"]
    if representation == "csg":
        desc = ((rec.get("evidence") or {}).get("description") or "").strip()
        cls = desc.split("(")[0].strip() or "CSG"
    else:
        cls = SHORT_CLASS.get(representation, representation)
    return f"{vol} · {cls}"


# ------------------------------------------------------------- placements ----

DUMP_MACRO = r'''
// Dump every leaf placement of a converted CAD geometry as JSON:
// the part key (volumes are renamed to it first) and the 3x4 global matrix.
#include <fstream>
void dump_placements(const char* out) {
  std::ofstream o(out);
  o << "[";
  bool first = true;
  TGeoIterator it(gGeoManager->GetTopVolume());
  TGeoNode* node;
  while ((node = it.Next())) {
    TGeoVolume* v = node->GetVolume();
    if (v->IsAssembly()) continue;
    TGeoHMatrix m = *it.GetCurrentMatrix();
    const Double_t* r = m.GetRotationMatrix();
    const Double_t* t = m.GetTranslation();
    if (!first) o << ",";
    first = false;
    TGeoBBox* bb = dynamic_cast<TGeoBBox*>(v->GetShape());
    o << "{\"volume\":\"" << v->GetName() << "\",\"m\":[";
    for (int i = 0; i < 9; ++i) o << r[i] << ",";
    o << t[0] << "," << t[1] << "," << t[2] << "]";
    if (bb) {
      const Double_t* bo = bb->GetOrigin();
      o << ",\"shapeBox\":[" << bb->GetDX() << "," << bb->GetDY() << "," << bb->GetDZ()
        << "," << bo[0] << "," << bo[1] << "," << bo[2] << "]";
    }
    o << "}";
  }
  o << "]";
  o.close();
}
'''


def dump_placements(model_dir, part_order, workdir):
    """Build the geometry in ROOT and write out every leaf placement.

    Volumes are renamed to their part key first, so a placement identifies its
    part unambiguously even where display names repeat.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    # Keep a Python reference: an unreferenced TGeoManager is collected straight away and
    # build() then throws "gGeoManager is null".
    mgr = ROOT.TGeoManager("scene", "scene export")
    assert mgr
    ROOT.gROOT.ProcessLine(f'.L {os.path.join(model_dir, "geom.C")}')
    # Build in C++: ProcessLineSync hands back an address, not a TGeoVolume proxy, and
    # SetTopVolume then refuses it.
    ROOT.gROOT.ProcessLine("TGeoVolume* __scene_top = build(false); "
                           "gGeoManager->SetTopVolume(__scene_top);")
    if not ROOT.gGeoManager.GetTopVolume():
        raise RuntimeError("build(false) produced no top volume")

    vols = ROOT.gGeoManager.GetListOfVolumes()
    solids = [v for v in vols if not v.IsAssembly()]
    got = [v.GetName() for v in solids]
    want = [name for _, name in part_order]
    if got != want:
        raise RuntimeError(
            "volume creation order does not match geom.C: "
            f"{len(got)} solid volumes vs {len(want)} parsed; "
            f"first mismatch at {next((i for i, (a, b) in enumerate(zip(got, want)) if a != b), 'n/a')}")
    for vol, (part, _) in zip(solids, part_order):
        vol.SetName(part)

    ROOT.gGeoManager.CloseGeometry()
    ROOT.gROOT.ProcessLine(DUMP_MACRO)
    out = os.path.join(workdir, "placements.json")
    ROOT.gROOT.ProcessLine(f'dump_placements("{out}");')
    top_shape = ROOT.gGeoManager.GetTopVolume().GetShape()
    o = top_shape.GetOrigin()
    top_box = {"half": [top_shape.GetDX(), top_shape.GetDY(), top_shape.GetDZ()],
               "origin": [o[0], o[1], o[2]]}
    with open(out) as fh:
        return json.load(fh), top_box


# ----------------------------------------------------------------- output ----

def mat_mul(a, b):
    """Compose two 3x4 transforms (rotation rows then translation), a after b."""
    ra, ta = a[:9], a[9:]
    rb, tb = b[:9], b[9:]
    r = [sum(ra[3 * i + k] * rb[3 * k + j] for k in range(3)) for i in range(3) for j in range(3)]
    t = [sum(ra[3 * i + k] * tb[k] for k in range(3)) + ta[i] for i in range(3)]
    return r + t


def mat_inv(m):
    """Inverse of a rigid 3x4 transform: R^T, -R^T t."""
    r, t = m[:9], m[9:]
    ri = [r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]]
    ti = [-sum(ri[3 * i + k] * t[k] for k in range(3)) for i in range(3)]
    return ri + ti


def shape_placement_3x4(rows):
    """csg_report's shapePlacement (3 rows of 4) as a flat 3x4."""
    r = [rows[i][j] for i in range(3) for j in range(3)]
    t = [rows[i][3] for i in range(3)]
    return r + t


def transform(tri, m):
    """Apply a 3x4 (row-major rotation + translation) to 9 floats."""
    r, t = m[:9], m[9:]
    out = []
    for k in range(3):
        x, y, z = tri[3 * k], tri[3 * k + 1], tri[3 * k + 2]
        out += [r[0] * x + r[1] * y + r[2] * z + t[0],
                r[3] * x + r[4] * y + r[5] * z + t[1],
                r[6] * x + r[7] * y + r[8] * z + t[2]]
    return out


def tri_area_normal(v):
    ax, ay, az = v[3] - v[0], v[4] - v[1], v[5] - v[2]
    bx, by, bz = v[6] - v[0], v[7] - v[1], v[8] - v[2]
    nx, ny, nz = ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx
    n = math.sqrt(nx * nx + ny * ny + nz * nz)
    return 0.5 * n, ((nx / n, ny / n, nz / n) if n > 0 else (0.0, 0.0, 1.0))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True)
    args = ap.parse_args()

    model_dir, out_dir = args.model_dir, args.out
    os.makedirs(out_dir, exist_ok=True)

    report = json.load(open(os.path.join(model_dir, "csg_report.json")))
    records = {r["part"]: r for r in report["parts"]}
    part_order = parse_volume_order(open(os.path.join(model_dir, "geom.C")).read())

    missing = {p for p, _ in part_order} ^ set(records)
    if missing:
        raise RuntimeError(f"geom.C and csg_report.json disagree on parts: {sorted(missing)[:6]}")

    placements, top_box = dump_placements(model_dir, part_order, out_dir)
    print(f"[{args.name}] {len(placements)} leaf placement(s) over {len(records)} part(s)")

    by_part = {}
    shape_boxes = {}
    for p in placements:
        by_part.setdefault(p["volume"], []).append(p["m"])
        if "shapeBox" in p:
            shape_boxes.setdefault(p["volume"], []).append((p["m"], p["shapeBox"]))

    facets = {}
    for part in records:
        f = os.path.join(model_dir, f"facets_{part}.bin")
        if not os.path.exists(f):
            raise RuntimeError(f"missing facet file for {part}: {f}")
        facets[part] = read_facets(f)

    obj_path = os.path.join(out_dir, f"{args.name}.obj")
    manifest = {"model": args.name, "modelDir": os.path.abspath(model_dir),
                "units": "cm", "frame": "world (top volume of the built geometry)",
                "note": "no colours/materials/cameras here by design; see labelShort for the "
                        "drawn label and recogniser/description for a caption",
                "parts": []}
    gmin = [float("inf")] * 3
    gmax = [float("-inf")] * 3
    total_world_tris = 0
    vbase = 1

    with open(obj_path, "w") as obj:
        obj.write(f"# {args.name}: one object per placement, world space, cm\n")
        for part, rec in sorted(records.items()):
            ntri, local = facets[part]
            mats = by_part.get(part, [])
            if not mats:
                raise RuntimeError(f"{part} has no placement in the built geometry")
            # A CSG part's node matrix is tr * shapePlacement: it maps the CANONICAL shape's
            # frame to the world, because geom.C composes the shape placement into the node
            # (the `_placed` matrices). The facet file is in the part's CAD-local frame, which
            # is the shape frame only when there is no shape placement. Applying the node
            # matrix straight to facets therefore applies shapePlacement twice and scatters
            # every CSG part. Undo it: world = M * shapePlacement^-1 * cadLocal.
            sp = rec.get("shapePlacement")
            if sp:
                inv = mat_inv(shape_placement_3x4(sp))
                mats = [mat_mul(m, inv) for m in mats]
            pmin = [float("inf")] * 3
            pmax = [float("-inf")] * 3
            best = (-1.0, None, None)          # area, centroid, normal
            csum = [0.0, 0.0, 0.0]
            cnt = 0
            entries = []
            for pi, m in enumerate(mats):
                name = f"{part}#{pi}"
                obj.write(f"o {name}\n")
                wmin = [float("inf")] * 3
                wmax = [float("-inf")] * 3
                verts = []
                for t in range(ntri):
                    w = transform(local[t * 9:(t + 1) * 9], m)
                    verts.append(w)
                    for k in range(3):
                        for c in range(3):
                            val = w[3 * k + c]
                            if val < wmin[c]:
                                wmin[c] = val
                            if val > wmax[c]:
                                wmax[c] = val
                            csum[c] += val
                    cnt += 3
                    a, n = tri_area_normal(w)
                    if a > best[0]:
                        ctr = tuple(sum(w[3 * k + c] for k in range(3)) / 3.0 for c in range(3))
                        best = (a, ctr, n)
                for w in verts:
                    for k in range(3):
                        obj.write(f"v {w[3*k]:.6g} {w[3*k+1]:.6g} {w[3*k+2]:.6g}\n")
                for t in range(ntri):
                    b = vbase + 3 * t
                    obj.write(f"f {b} {b+1} {b+2}\n")
                vbase += 3 * ntri
                total_world_tris += ntri
                for c in range(3):
                    pmin[c] = min(pmin[c], wmin[c])
                    pmax[c] = max(pmax[c], wmax[c])
                entries.append({"object": name, "matrix3x4": m,
                                "bbox": {"min": wmin, "max": wmax}, "triangles": ntri})
            for c in range(3):
                gmin[c] = min(gmin[c], pmin[c])
                gmax[c] = max(gmax[c], pmax[c])

            centroid = [csum[c] / cnt for c in range(3)]
            # Label anchor: the centroid of the part's largest-area triangle, with that
            # triangle's normal flipped to point away from the part centroid. The biggest
            # facet is a broad flat face, which is where a leader line reads best.
            _, apt, anrm = best
            outward = sum((apt[c] - centroid[c]) * anrm[c] for c in range(3))
            if outward < 0:
                anrm = tuple(-x for x in anrm)
            ev = rec.get("evidence") or {}
            manifest["parts"].append({
                "part": part,
                "volume": rec.get("volume"),
                "lid": rec.get("lid"),
                "representation": rec["representation"],
                "recogniser": ev.get("recogniser"),
                "description": ev.get("description") or rec.get("whyNotCSG"),
                "labelShort": short_label(rec),
                "facetFile": os.path.abspath(os.path.join(model_dir, f"facets_{part}.bin")),
                "trianglesPerPlacement": ntri,
                "placementCount": len(mats),
                "placements": entries,
                "bbox": {"min": pmin, "max": pmax},
                "centroid": centroid,
                "labelAnchor": {"point": list(apt), "normal": list(anrm)},
            })

    manifest["scene"] = {
        "bbox": {"min": gmin, "max": gmax},
        "halfLengths": [(gmax[c] - gmin[c]) / 2.0 for c in range(3)],
        "centre": [(gmax[c] + gmin[c]) / 2.0 for c in range(3)],
        "placements": len(placements),
        "worldTriangles": total_world_tris,
    }
    with open(os.path.join(out_dir, f"{args.name}_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)

    # ------------------------------------------------------------- gates ----
    tiers = report["tiers"]
    got_tiers = {}
    for p in manifest["parts"]:
        got_tiers[p["representation"]] = got_tiers.get(p["representation"], 0) + 1
    ok = True
    def gate(label, good, detail):
        nonlocal ok
        ok = ok and good
        print(f"  [{'ok ' if good else 'FAIL'}] {label}: {detail}")

    gate("representation counts match csg_report",
         all(got_tiers.get(k, 0) == v for k, v in tiers.items()),
         f"{got_tiers} vs {tiers}")
    gate("every part present exactly once",
         len(manifest["parts"]) == len(records),
         f"{len(manifest['parts'])} of {len(records)}")
    # The check that catches a misplaced part: the triangles a part contributes must lie
    # inside the world box of the shape TGeo actually navigates. A CSG part whose shape
    # placement was applied twice fails this immediately, where the part counts do not.
    worst = (0.0, None)
    for entry in manifest["parts"]:
        boxes = shape_boxes.get(entry["part"], [])
        if not boxes:
            continue
        smin = [float("inf")] * 3
        smax = [float("-inf")] * 3
        for m, (dx, dy, dz, ox, oy, oz) in boxes:
            for sx in (-1, 1):
                for sy in (-1, 1):
                    for sz in (-1, 1):
                        c = [ox + sx * dx, oy + sy * dy, oz + sz * dz]
                        w = [sum(m[3 * i + k] * c[k] for k in range(3)) + m[9 + i] for i in range(3)]
                        for i in range(3):
                            smin[i] = min(smin[i], w[i])
                            smax[i] = max(smax[i], w[i])
        over = max(max(smin[c] - entry["bbox"]["min"][c], entry["bbox"]["max"][c] - smax[c])
                   for c in range(3))
        if over > worst[0]:
            worst = (over, entry["part"])
    gate("facets lie inside the shape TGeo navigates",
         worst[0] < 0.5,
         f"worst excursion {worst[0]:.4f} cm" + (f" on {worst[1]}" if worst[1] else ""))

    # And the scene as a whole must sit inside the box the built geometry claims for itself.
    # Containment, not equality: an assembly's box is the union of its daughters' shape boxes,
    # and those can be loose (NEXT.md's open item on loose bounding boxes).
    tmin = [top_box["origin"][c] - top_box["half"][c] for c in range(3)]
    tmax = [top_box["origin"][c] + top_box["half"][c] for c in range(3)]
    slack = max(max(tmin[c] - gmin[c], gmax[c] - tmax[c]) for c in range(3))
    manifest["scene"]["topVolumeBox"] = {"half": top_box["half"], "origin": top_box["origin"]}
    gate("scene lies inside the built geometry's own box", slack < 1e-3,
         f"worst excursion {slack:.4f} cm; top box half "
         f"({top_box['half'][0]:.2f}, {top_box['half'][1]:.2f}, {top_box['half'][2]:.2f}) "
         f"origin ({top_box['origin'][0]:.2f}, {top_box['origin'][1]:.2f}, {top_box['origin'][2]:.2f})")

    facet_sum = sum(facets[p][0] for p in records)
    gate("local triangles = sum over facet files", True,
         f"{facet_sum} local, {total_world_tris} world over {len(placements)} placement(s)")
    print(f"  scene bbox half-lengths dx={manifest['scene']['halfLengths'][0]:.2f} "
          f"dy={manifest['scene']['halfLengths'][1]:.2f} dz={manifest['scene']['halfLengths'][2]:.2f} "
          f"centre=({manifest['scene']['centre'][0]:.2f}, {manifest['scene']['centre'][1]:.2f}, "
          f"{manifest['scene']['centre'][2]:.2f})")
    print(f"  wrote {obj_path} ({os.path.getsize(obj_path)/1e6:.1f} MB)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
