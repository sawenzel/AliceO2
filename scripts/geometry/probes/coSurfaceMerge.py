#!/usr/bin/env python3.10
"""Stream U probe: is the exporter fragmenting ONE analytic surface into several faces, and if the
co-surface faces were merged back, would the remaining boundary be iso in the recognised (phi, h)
domain -- i.e. would the *existing* parametric trim path already carry it?

The hypothesis
--------------
`Stream_R_CoSurfaceTrims.md` §4 measured that **380 of the 678 faces (56%)** of the 15 ALICE3
solids that recognise-but-do-not-emit have a neighbouring face lying on their *own* analytic
surface, and §5.3 found that all 7 of its residual failures carry one. That is a NURBS patch seam:
the exporter wrote one cylinder/cone/sphere as two or four patches. An edge that is "not iso" may
therefore be an *internal seam* of a surface that was cut up for the exporter's convenience, not a
real boundary at all.

If so, merging the co-surface faces back into one face before trimming removes the edge outright,
and the remaining true boundary may well be iso -- which the shipping
`_recognized_quadric_wire_block` already handles, with no new format, no kernel change and no
approximation.

This probe measures exactly that and nothing else. It changes no production code.

What it reports
---------------
1. **Fragmentation** -- faces grouped by the analytic surface they lie on, per solid: group-size
   distribution, worst case, named.
2. **Seam vs true boundary** -- every currently-rejected boundary edge, split by whether its
   neighbouring face lies on the *same* analytic surface (seam) or a different one (true boundary).
3. **The decisive number** -- per solid, after dissolving every seam: does every remaining boundary
   edge pass the shipped iso test, and does the merged face still satisfy the rest of the shipped
   wire block (one outer wire, >= 3 edges per wire, phi sweep <= 2pi)?
4. **Merged-domain sanity** -- is the merged (phi, other) domain a simple rectangle, which is what
   the existing trim record can express? Wrap, holes and corner counts are reported per group.

The surface-identity decision, and why it is exact rather than fitted
--------------------------------------------------------------------
Two faces are on the same analytic surface, or they are not. The probe decides it *geometrically*,
never by comparing parameter tuples: it projects a grid of points from face A's own patch onto A's
recognised ideal surface (a closed form for plane/cylinder/cone/sphere/torus, so the projected
points lie on A *exactly*), measures their distance to B's implicit surface, and symmetrises. The
result is a length in cm -- the largest distance between the two surfaces over the region either
patch occupies -- and it is compared against the **model's own declared BRep tolerance**, the same
yardstick `Stream_O_ImplicitTrims.md` §3 calibrated its bucket-B threshold with.

`--merge-factor` sweeps that threshold and the summary prints the whole sweep, so the constant is
inspectable rather than asserted. The separation histogram is reported for both populations (pairs
that merge and the nearest pair that does not), because a threshold with no measured gap either
side of it is a tuned constant.

Self-checks it runs before reporting anything
---------------------------------------------
1. **The per-edge iso verdict must be the shipped one.** The probe builds its own (phi, other)
   frame per group (the reference direction is arbitrary and iso-ness is invariant under it). For
   every *unmerged* face it re-runs `probes/trimEdgeCensus.py`'s `wire_block_census`, which is
   itself checked against the shipping `recognize_and_extract_face` with 0 face-verdict mismatches,
   and compares edge for edge. A non-zero mismatch count means the probe is not measuring the
   shipped criterion and nothing else it prints means anything.
2. **A negative control on the grouping.** `--negative-control` reports, per solid, the smallest
   surface separation among same-kind face pairs that were NOT merged -- the nearest distinct
   surface. If that number is not orders above the threshold the grouping is chaining.
3. **A synthetic positive control.** `--synthetic` builds a cylinder deliberately cut into four
   NURBS patches by an exporter-like split and checks that the four are grouped into one, and that
   a genuinely different cylinder in the same solid is not.

Usage
-----
    export ALIBUILD_WORK_DIR=$HOME/alisw/sw
    B=$HOME/alisw/sw/BUILD/O2-latest/O2 ; SW=$HOME/alisw/sw/ubuntu2404_aarch64
    cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
    export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
    cd $HOME/alisw/O2/scripts/geometry

    python3 probes/coSurfaceMerge.py --model ALICE3:ALICE_3_example/CAD_noETA.stp \
        --negative-control --json /tmp/u/alice3.json
    python3 probes/coSurfaceMerge.py --model Bagger:STEP_examples/Bagger.step \
        --model as1:STEP_examples/as1-oc-214.stp --fixtures --json /tmp/u/controls.json
    python3 probes/coSurfaceMerge.py --synthetic
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE))

from csg.occ_env import ensure_occ  # noqa: E402

ensure_occ()

import numpy as np  # noqa: E402

from OCC.Core.BRep import BRep_Tool  # noqa: E402
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCC.Core.BRepTools import breptools  # noqa: E402
from OCC.Core.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX,  # noqa: E402
                             TopAbs_REVERSED, TopAbs_FORWARD)
from OCC.Core.TopExp import topexp  # noqa: E402
from OCC.Core.TopTools import (TopTools_IndexedMapOfShape,  # noqa: E402
                               TopTools_IndexedDataMapOfShapeListOfShape,
                               TopTools_ListIteratorOfListOfShape)
from OCC.Core.TopoDS import topods  # noqa: E402
from OCC.Extend.TopologyUtils import TopologyExplorer  # noqa: E402

import O2_CADtoTGeo as C  # noqa: E402
import trimEdgeCensus as T  # noqa: E402


_N_SAMPLES = 9          # exactly `_recognized_quadric_wire_block`'s
_MERGE_FACTOR = 2.0     # declared BRep tolerances; swept in the summary
_TWO_PI = 2.0 * math.pi


# ---------------------------------------------------------------------------------------------
# surface identity, decided geometrically
# ---------------------------------------------------------------------------------------------

def project_onto_ideal(model, P):
    """Project points onto the model's ideal analytic surface -- a closed form per kind, so the
    result lies on the surface exactly (to double precision), not to a fitted tolerance."""
    P = np.atleast_2d(np.asarray(P, float))
    kind = model["kind"]
    if kind == "plane":
        n = model["normal"]
        return P - np.outer((P - model["point"]) @ n, n)
    if kind == "sphere":
        c = model["centre"]
        rel = P - c
        d = np.linalg.norm(rel, axis=1)
        d = np.where(d < 1e-30, 1.0, d)
        return c + rel * (model["radius"] / d)[:, None]
    if kind == "cylinder":
        o, a, R = model["origin"], model["axis"], model["radius"]
        rel = P - o
        h = rel @ a
        perp = rel - np.outer(h, a)
        d = np.linalg.norm(perp, axis=1)
        d = np.where(d < 1e-30, 1.0, d)
        return o + np.outer(h, a) + perp * (R / d)[:, None]
    if kind == "cone":
        ap, a, alpha = model["apex"], model["axis"], model["half_angle"]
        rel = P - ap
        h = rel @ a
        perp = rel - np.outer(h, a)
        r = np.linalg.norm(perp, axis=1)
        u = np.where(r[:, None] < 1e-30, 0.0, perp / np.where(r < 1e-30, 1.0, r)[:, None])
        t = np.maximum(0.0, r * math.sin(alpha) + h * math.cos(alpha))
        return ap + u * (t * math.sin(alpha))[:, None] + np.outer(t * math.cos(alpha), a)
    if kind == "torus":
        c, a, R, r0 = model["centre"], model["axis"], model["major"], model["minor"]
        rel = P - c
        z = rel @ a
        perp = rel - np.outer(z, a)
        rho = np.linalg.norm(perp, axis=1)
        u = np.where(rho[:, None] < 1e-30, 0.0, perp / np.where(rho < 1e-30, 1.0, rho)[:, None])
        ring = c + u * R
        w = P - ring
        wn = np.linalg.norm(w, axis=1)
        wn = np.where(wn < 1e-30, 1.0, wn)
        return ring + w * (r0 / wn)[:, None]
    return None


def surface_separation(mA, PA, mB, PB):
    """The largest distance between surfaces A and B over the region either patch occupies, in
    native CAD units. Symmetric by construction; 0.0 iff the two are the same surface there."""
    if mA["kind"] != mB["kind"]:
        return float("inf")
    QA = project_onto_ideal(mA, PA)
    QB = project_onto_ideal(mB, PB)
    if QA is None or QB is None:
        return float("inf")
    d1 = float(T.surface_distance(mB, QA).max())
    d2 = float(T.surface_distance(mA, QB).max())
    return max(d1, d2)


def face_surface_points(face, n=5):
    """A small grid of points on the face's own stored surface, over its (u, v) rectangle."""
    try:
        umin, umax, vmin, vmax = breptools.UVBounds(face)
        surf = BRep_Tool.Surface(face)
    except Exception:
        return None
    if not all(math.isfinite(x) for x in (umin, umax, vmin, vmax)):
        return None
    us = umin + (umax - umin) * (np.arange(n) + 0.5) / n
    vs = vmin + (vmax - vmin) * (np.arange(n) + 0.5) / n
    out = np.empty((n * n, 3))
    k = 0
    for u in us:
        for v in vs:
            p = surf.Value(float(u), float(v))
            out[k] = (p.X(), p.Y(), p.Z())
            k += 1
    return out


class UnionFind:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


# ---------------------------------------------------------------------------------------------
# the recognised (phi, other) frame, and the shipped iso test on one edge
# ---------------------------------------------------------------------------------------------

def projector_for_model(model, scale_to_cm):
    """The `project` closure `recognize_and_extract_face` builds, from a model dict.

    The reference direction is arbitrary (the shipped code takes it from the recogniser's own
    `refu`); iso-ness is invariant under it, and self-check 1 verifies that edge for edge against
    the shipped frame.
    """
    s = scale_to_cm
    kind = model["kind"]
    if kind in ("cylinder", "cone"):
        axis = T._unit(model["axis"])
        refu = C._arbitrary_orthonormal_frame(axis)
        refu = T._unit(refu - np.dot(refu, axis) * axis)
        e2 = np.cross(axis, refu)
        origin = np.asarray(model["origin" if kind == "cylinder" else "apex"], float)

        def project(pnt):
            rel = np.array([pnt.X(), pnt.Y(), pnt.Z()]) - origin
            return (math.atan2(float(np.dot(rel, e2)), float(np.dot(rel, refu))),
                    float(np.dot(rel, axis)) * s)
        return project, "h"
    if kind == "sphere":
        centre = np.asarray(model["centre"], float)
        polar = np.array([0.0, 0.0, 1.0])
        refu = C._arbitrary_orthonormal_frame(polar)
        e2 = np.cross(polar, refu)
        radius = float(model["radius"])

        def project(pnt):
            rel = (np.array([pnt.X(), pnt.Y(), pnt.Z()]) - centre) / radius
            theta = math.acos(max(-1.0, min(1.0, float(np.dot(rel, polar)))))
            return (math.atan2(float(np.dot(rel, e2)), float(np.dot(rel, refu))), theta)
        return project, "theta"
    return None, None


def edge_uv_samples(edge, project, n=_N_SAMPLES):
    """`_recognized_quadric_wire_block`'s own sampling of one edge, phi unwrapped continuously
    *within* the edge, in the 3D curve's own parameter direction.

    Two deliberate departures from the shipped loop, neither of which can move an iso verdict:
    the shipped code chains the unwrap across the whole wire, which adds a *constant* per edge;
    and it walks the edge backwards when the wire says the edge is REVERSED, which reverses the
    sample order. The iso test looks only at `max - min` of each run, so neither changes it --
    and self-check 1 verifies that, edge for edge, against the shipped frame and the shipped
    traversal. Sampling in the curve's own direction is what lets the merged boundary be chained
    into loops through `FirstVertex`/`LastVertex`, which are also forward-oriented.
    """
    if BRep_Tool.Degenerated(edge):
        return None
    try:
        curve3d, first, last = BRep_Tool.Curve(edge)
    except Exception:
        return None
    if curve3d is None:
        return None
    phis, others = [], []
    prev = None
    for k in range(n):
        t = k / (n - 1.0)
        phi, other = project(curve3d.Value(first + t * (last - first)))
        if prev is not None:
            d = phi - prev
            d -= _TWO_PI * math.floor((d + math.pi) / _TWO_PI)
            phi = prev + d
        prev = phi
        phis.append(phi)
        others.append(other)
    return np.array(phis), np.array(others)


def iso_verdict(phis, others, tol_other, tol_phi=1e-7):
    span_other = float(others.max() - others.min())
    span_phi = float(phis.max() - phis.min())
    return (span_other <= tol_other, span_phi <= tol_phi, span_other, span_phi)


# ---------------------------------------------------------------------------------------------
# topology helpers
# ---------------------------------------------------------------------------------------------

def ancestor_faces(anc, edge):
    try:
        lst = anc.FindFromKey(edge)
    except Exception:
        return []
    out = []
    it = TopTools_ListIteratorOfListOfShape(lst)
    while it.More():
        out.append(topods.Face(it.Value()))
        it.Next()
    return out


def edge_vertices(edge, vmap):
    fwd = topods.Edge(edge.Oriented(TopAbs_FORWARD))
    try:
        v0 = topexp.FirstVertex(fwd)
        v1 = topexp.LastVertex(fwd)
    except Exception:
        return None, None
    return vmap.FindIndex(v0), vmap.FindIndex(v1)


def chain_loops(edge_ids, endpoints):
    """Chain a set of edges into closed loops through shared vertices.

    `endpoints[eid] = (v0, v1)`. Returns (loops, problem) where a loop is a list of
    (eid, forward) in traversal order, and `problem` is None or a string naming why the boundary
    is not a disjoint union of simple closed loops."""
    inc = defaultdict(list)
    for eid in edge_ids:
        v0, v1 = endpoints[eid]
        inc[v0].append(eid)
        if v1 != v0:
            inc[v1].append(eid)
        else:
            inc[v0].append(eid)      # a closed edge visits its single vertex twice
    bad = [v for v, es in inc.items() if len(es) != 2]
    if bad:
        return [], f"{len(bad)} boundary vertices with degree != 2"
    unused = set(edge_ids)
    loops = []
    while unused:
        start = next(iter(unused))
        v0, v1 = endpoints[start]
        loop = [(start, True)]
        unused.discard(start)
        cur_v = v1
        cur_e = start
        guard = 0
        while cur_v != v0 and guard < 100000:
            guard += 1
            nxt = [e for e in inc[cur_v] if e != cur_e and e in unused]
            if not nxt:
                return [], "open chain: could not close a boundary loop"
            e = nxt[0]
            a, b = endpoints[e]
            fwd = (a == cur_v)
            loop.append((e, fwd))
            unused.discard(e)
            cur_v = b if fwd else a
            cur_e = e
        loops.append(loop)
    return loops, None


def cluster(values, tol):
    """Distinct values of a 1-D set, at the given tolerance. Returns the cluster representatives."""
    if not values:
        return []
    vs = sorted(values)
    reps = [vs[0]]
    for v in vs[1:]:
        if v - reps[-1] > tol:
            reps.append(v)
    return reps


# ---------------------------------------------------------------------------------------------
# per-solid measurement
# ---------------------------------------------------------------------------------------------

def measure_solid(lid, name, shape, scale_to_cm, args):
    s = scale_to_cm
    faces = list(TopologyExplorer(shape).faces())
    fmap = TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape, TopAbs_FACE, fmap)
    emap = TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape, TopAbs_EDGE, emap)
    vmap = TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape, TopAbs_VERTEX, vmap)
    anc = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, anc)

    out = {"lid": lid, "name": name, "nFaces": len(faces)}

    # ---- per face: analytic model, recognition, sample points, shipped verdict ---------------
    models = [None] * len(faces)          # analytic model dict, or None
    recs = [None] * len(faces)            # recognition dict for wire-block faces
    pts = [None] * len(faces)
    tolf = [0.0] * len(faces)
    wireblock = [False] * len(faces)      # face reaches `_recognized_quadric_wire_block`
    shipped_edges = {}                    # face index -> census entries
    shipped_status = {}
    n_unsupported = 0
    n_trim_declined = 0
    face_reasons = []

    for i, face in enumerate(faces):
        tolf[i] = float(BRep_Tool.Tolerance(face)) * s
        adaptor = BRepAdaptor_Surface(face)
        surf_type = C._SURFACE_TYPE_NAMES.get(adaptor.GetType(), "unknown")
        extractor = C._FACE_EXTRACTORS.get(surf_type)
        if extractor is None:
            record, reason = None, f"{surf_type} face extraction not implemented yet"
        else:
            try:
                record, reason = extractor(face, s)
            except Exception as exc:
                record, reason = None, f"{surf_type} extractor raised {exc!r}"
        recognised_reason = None
        if record is None:
            try:
                rr, rn = C.recognize_and_extract_face(face, s)
            except Exception as exc:
                rr, rn = None, f"recognition raised {exc!r}"
            if rr is not None:
                record, reason = rr, None
            elif rn is not None:
                recognised_reason = rn
                reason = f"{reason}; recognition attempted: {rn}"
        if record is None:
            face_reasons.append(reason)
            if recognised_reason is not None:
                n_trim_declined += 1
            else:
                n_unsupported += 1

        m = T.stored_model(face)
        rec = None
        if m is None or (extractor is not None and record is not None
                         and record.get("recognized") is not None):
            try:
                rec = T.recognize(adaptor, breptools.UVBounds(face))
            except Exception:
                rec = None
            if rec is not None and m is None:
                m = T.recognized_model(rec)
        models[i] = m
        pts[i] = face_surface_points(face, n=args.surface_grid)

        # A face reaches `_recognized_quadric_wire_block` exactly when
        # `probes/trimEdgeCensus.py` says it does: its stored surface did not extract natively,
        # and the shipping recogniser returns a non-plane quadric.
        if extractor is not None and record is not None and record.get("recognized") is None:
            continue
        if rec is None:
            try:
                rec = T.recognize(adaptor, breptools.UVBounds(face))
            except Exception:
                rec = None
        if rec is None or rec["kind"] == "plane":
            continue
        models[i] = T.recognized_model(rec)
        recs[i] = rec
        proj, _on = T.projector_for(rec, s)
        if proj is None:
            continue
        wireblock[i] = True
        st, entries, _w = T.wire_block_census(face, proj)
        shipped_status[i] = st
        shipped_edges[i] = entries

    out["nUnsupportedFaces"] = n_unsupported
    out["nTrimDeclinedFaces"] = n_trim_declined
    out["emitted"] = bool(faces) and not face_reasons
    out["nWireBlockFaces"] = sum(1 for w in wireblock if w)

    # ---- 1. group faces by the analytic surface they lie on ----------------------------------
    # TWO groupings, and they answer different questions.
    #   `uf`  -- all faces lying on one analytic surface, whether or not they touch. This is the
    #            fragmentation census the brief asks for in (1).
    #   `mf`  -- the connected components of "shares a boundary edge AND lies on the same
    #            surface". This is what a merge pre-pass may actually do: two disjoint patches of
    #            one cylinder (six windows cut in one tube) are ONE surface but are NOT one face,
    #            and fusing them would invent a face whose region is disconnected.
    uf = UnionFind(len(faces))
    mf = UnionFind(len(faces))
    by_kind = defaultdict(list)
    for i, m in enumerate(models):
        if m is not None and pts[i] is not None:
            by_kind[m["kind"]].append(i)
    same_surface = set()
    sep_merged, sep_split = [], []
    npts = args.surface_grid ** 2
    for kind, idxs in by_kind.items():
        n = len(idxs)
        if n < 2:
            continue
        # Q[j] lies EXACTLY on face j's own ideal surface (closed-form projection), so
        # D[i, j] = max distance of Q[j] from surface i is the one-sided surface-to-surface
        # separation over the region face j occupies. Symmetrised below.
        Q = [project_onto_ideal(models[j], pts[j]) for j in idxs]
        if any(q is None for q in Q):
            continue
        Qall = np.vstack(Q)
        D = np.empty((n, n))
        for a in range(n):
            D[a] = T.surface_distance(models[idxs[a]], Qall).reshape(n, npts).max(axis=1)
        S = np.maximum(D, D.T) * s
        tolm = np.maximum(np.array([tolf[j] for j in idxs])[:, None],
                          np.array([tolf[j] for j in idxs])[None, :])
        R = np.where(tolm > 0, S / np.where(tolm > 0, tolm, 1.0), np.inf)
        ia, ib = np.triu_indices(n, k=1)
        hit = S[ia, ib] <= args.merge_factor * tolm[ia, ib]
        for a, b in zip(ia[hit], ib[hit]):
            i, j = idxs[a], idxs[b]
            uf.union(i, j)
            same_surface.add((min(i, j), max(i, j)))
            sep_merged.append((float(R[a, b]), float(S[a, b]), i, j))
        if (~hit).any():
            k = int(np.argmin(R[ia, ib][~hit]))
            a, b = ia[~hit][k], ib[~hit][k]
            sep_split.append((float(R[a, b]), float(S[a, b]), idxs[a], idxs[b]))
    # adjacency: two faces are adjacent iff some edge of the solid names both
    n_adj_cosurface = 0
    for k in range(1, emap.Size() + 1):
        owners = sorted({fmap.FindIndex(g) - 1 for g in ancestor_faces(anc, emap.FindKey(k))})
        for a in range(len(owners)):
            for b in range(a + 1, len(owners)):
                if (owners[a], owners[b]) in same_surface:
                    if mf.find(owners[a]) != mf.find(owners[b]):
                        n_adj_cosurface += 1
                    mf.union(owners[a], owners[b])
    out["nAdjacentCoSurfacePairs"] = n_adj_cosurface

    def _collect(finder):
        g = defaultdict(list)
        for i, m in enumerate(models):
            if m is not None and pts[i] is not None:
                g[finder.find(i)].append(i)
        return g

    groups = _collect(uf)              # surfaces
    mgroups = _collect(mf)             # mergeable faces
    group_of = {}
    for gid, members in mgroups.items():
        for i in members:
            group_of[i] = gid

    sizes = Counter(len(v) for v in groups.values())
    msizes = Counter(len(v) for v in mgroups.values())
    frag = [{"kind": models[v[0]]["kind"], "faces": sorted(v),
             "nWireBlock": sum(1 for i in v if wireblock[i]),
             "nMergeComponents": len({mf.find(i) for i in v})}
            for v in groups.values() if len(v) > 1]
    out["groupSizeHistogram"] = {str(k): v for k, v in sorted(sizes.items())}
    out["mergeGroupSizeHistogram"] = {str(k): v for k, v in sorted(msizes.items())}
    out["nSurfaces"] = len(groups)
    out["nMergeGroups"] = len(mgroups)
    out["nFacesWithAnalyticSurface"] = sum(len(v) for v in groups.values())
    out["maxGroupSize"] = max([len(v) for v in groups.values()], default=0)
    out["maxMergeGroupSize"] = max([len(v) for v in mgroups.values()], default=0)
    out["fragmentedGroups"] = sorted(frag, key=lambda g: -len(g["faces"]))[:40]
    out["nDisconnectedSurfaces"] = sum(1 for g in frag if g["nMergeComponents"] > 1)
    sep_merged.sort()
    sep_split.sort()
    out["mergedSeparationOverTol"] = {
        "n": len(sep_merged),
        "max": sep_merged[-1][0] if sep_merged else None,
        "p50": sep_merged[len(sep_merged) // 2][0] if sep_merged else None,
        "maxCm": max((x[1] for x in sep_merged), default=None),
    }
    out["nearestSplitSeparationOverTol"] = sep_split[0][0] if sep_split else None
    out["nearestSplitSeparationCm"] = sep_split[0][1] if sep_split else None
    if args.negative_control and sep_split:
        r, sep, i, j = sep_split[0]
        out["nearestSplitPair"] = {"faces": [i, j], "kind": models[i]["kind"],
                                   "sepCm": sep, "sepOverTol": r}

    # ---- 2. seam vs true boundary on every currently-rejected edge ---------------------------
    # Cross-check, from a different instrument: `Stream_R_CoSurfaceTrims.md` §4 calls a neighbour
    # *coincident* when its implicit function is identically zero over the face's own surface
    # samples. That is a statement about function values on a (u, v) grid; the grouping above is
    # a statement about closed-form projections between two ideal surfaces. They share no code,
    # so agreeing is evidence and disagreeing is a defect.
    coin_cache = {}

    def coincident_with(i, j):
        """max |f_j| over face i's own surface samples, in cm, and Stream_R's own verdict."""
        key = (i, j)
        if key in coin_cache:
            return coin_cache[key]
        P = face_surface_points(faces[i], n=max(8, args.surface_grid))
        m = models[j]
        if P is None or m is None:
            coin_cache[key] = (None, None)
            return coin_cache[key]
        f = float(np.abs(T.surface_distance(m, P)).max()) * s
        diag = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0))) * s
        tol = max(20.0 * tolf[i], 1e-9 * max(diag, 1.0))
        coin_cache[key] = (f, bool(f <= tol))
        return coin_cache[key]

    edge_rec = []
    for i in sorted(shipped_edges):
        for entry in shipped_edges[i]:
            edge = entry["edge"]
            eid = emap.FindIndex(edge)
            nbrs = [g for g in ancestor_faces(anc, edge) if not g.IsSame(faces[i])]
            nbr_idx = [fmap.FindIndex(g) - 1 for g in nbrs]
            same_group = [k for k in nbr_idx
                          if k in group_of and i in group_of and group_of[k] == group_of[i]]
            if not nbrs:
                cls = "periodic-seam"        # the face is its own neighbour
            elif same_group:
                cls = "seam"
            else:
                cls = "true-boundary"
            rec_e = {
                "solid": name, "face": i, "edgeId": eid,
                "accepted": entry["accepted"], "degenerate": entry["degenerate"],
                "class": cls,
                "nNeighbours": len(nbrs),
                "neighbourFaces": nbr_idx,
                "neighbourKinds": [None if models[k] is None else models[k]["kind"]
                                   for k in nbr_idx],
            }
            if args.coincidence_check and not entry["degenerate"] and nbr_idx:
                fmax, coin = coincident_with(i, nbr_idx[0])
                rec_e["coincidentMaxAbsFCm"] = fmax
                rec_e["coincidentVerdict"] = coin
            edge_rec.append(rec_e)
    out["edges"] = edge_rec

    # ---- 3 + 4. the merge simulation ---------------------------------------------------------
    # Every wire-block face is replaced by its merge group. An edge internal to the group is
    # dissolved; every other edge stays and must pass the shipped iso test in ONE frame for the
    # whole merged face.
    endpoints = {}
    merged = []
    handled = set()
    for i in sorted(shipped_edges):
        if i in handled:
            continue
        gid = group_of.get(i)
        members = [k for k in mgroups.get(gid, [i]) if wireblock[k]] if gid is not None else [i]
        members = sorted(set(members) | {i})
        handled.update(members)

        model = models[members[0]]
        project, other_name = projector_for_model(model, s)
        entry = {"solid": name, "faces": members, "kind": model["kind"],
                 "nMembers": len(members),
                 # a merge pre-pass has two honest policies. P-all fuses every adjacent
                 # co-surface component; P-need fuses only a component that contains a face the
                 # shipped wire block declines today, so geometry that already emits is never
                 # touched. `needsMerge` is what separates them.
                 "needsMerge": any(shipped_status.get(k) is not None for k in members),
                 "shippedStatus": {str(k): shipped_status.get(k) for k in members}}
        if project is None:
            entry["status"] = "no projector"
            merged.append(entry)
            continue

        # boundary vs internal
        all_edges = {}
        for k in members:
            for _w, is_outer, wedges in C._face_wire_edges(faces[k]):
                for e, _v in wedges:
                    eid = emap.FindIndex(e)
                    all_edges.setdefault(eid, {"edge": e, "faces": set(), "outer": is_outer,
                                               "degenerate": bool(BRep_Tool.Degenerated(e)),
                                               "uses": 0})
                    all_edges[eid]["faces"].add(k)
                    all_edges[eid]["uses"] += 1
        boundary, internal, degen = [], [], []
        mset = set(members)
        for eid, info in all_edges.items():
            owners = {fmap.FindIndex(g) - 1 for g in ancestor_faces(anc, info["edge"])}
            if info["degenerate"]:
                degen.append(eid)
                continue
            # internal iff every face referencing this edge is a member of the group
            if owners and owners <= mset:
                internal.append(eid)
            else:
                boundary.append(eid)
        entry["nEdgesTotal"] = len(all_edges)
        entry["nInternalDissolved"] = len(internal)
        entry["nBoundary"] = len(boundary)
        entry["nDegenerate"] = len(degen)
        # Edge identity (`applyEdgeIdentityClosure`) demands every model edge id appear exactly
        # twice, with opposite sense. A dissolved seam must therefore have been used exactly
        # twice inside the group -- once by each side -- so that removing both occurrences leaves
        # every surviving id still at two. A dissolved edge used any other number of times would
        # unbalance the count, and that is the number to watch, not an argument.
        entry["nDissolvedNotUsedTwice"] = sum(1 for eid in internal
                                              if all_edges[eid]["uses"] != 2)
        entry["nBoundaryUsedTwiceInside"] = sum(1 for eid in boundary
                                                if all_edges[eid]["uses"] != 1)

        # sample every remaining boundary edge in the ONE merged frame
        samples = {}
        bad_sample = False
        for eid in boundary:
            sm = edge_uv_samples(all_edges[eid]["edge"], project)
            if sm is None:
                bad_sample = True
                break
            samples[eid] = sm
        if bad_sample:
            entry["status"] = "boundary edge has no 3D curve"
            merged.append(entry)
            continue
        if not samples:
            entry["status"] = "no boundary edge left"
            merged.append(entry)
            continue
        all_other = np.concatenate([o for _p, o in samples.values()])
        tol_other = 1e-6 * max(1.0, float(all_other.max() - all_other.min()))
        n_iso = 0
        non_iso = []
        kinds = {}
        for eid in boundary:
            ph, ot = samples[eid]
            iso_o, iso_p, so, sp = iso_verdict(ph, ot, tol_other)
            if iso_o or iso_p:
                n_iso += 1
                kinds[eid] = "other" if iso_o else "phi"
            else:
                non_iso.append({"edgeId": eid, "spanOther": so, "spanPhi": sp,
                                "tolOther": tol_other})
        entry["nBoundaryIso"] = n_iso
        entry["nBoundaryNonIso"] = len(non_iso)
        entry["nonIso"] = non_iso[:10]
        entry["allBoundaryIso"] = (len(non_iso) == 0)

        # ---- 4. the merged domain: loops, wrap, rectangle
        for eid in boundary:
            endpoints[eid] = edge_vertices(all_edges[eid]["edge"], vmap)
        loops, problem = chain_loops(boundary, endpoints)
        entry["nLoops"] = len(loops)
        entry["loopProblem"] = problem
        if problem is None and loops:
            loop_info = []
            for loop in loops:
                phi_run = []
                oth_run = []
                cur = 0.0
                base = None
                for eid, fwd in loop:
                    ph, ot = samples[eid]
                    if not fwd:
                        ph, ot = ph[::-1], ot[::-1]
                    if base is None:
                        base = ph[0]
                        cur = ph[0]
                    for k in range(len(ph)):
                        d = ph[k] - cur
                        d -= _TWO_PI * math.floor((d + math.pi) / _TWO_PI)
                        cur = cur + d
                        phi_run.append(cur)
                        oth_run.append(ot[k])
                winding = phi_run[-1] - phi_run[0]
                loop_info.append({
                    "nEdges": len(loop),
                    "phiMin": min(phi_run), "phiMax": max(phi_run),
                    "otherMin": min(oth_run), "otherMax": max(oth_run),
                    "windingTurns": winding / _TWO_PI,
                    "wraps": abs(abs(winding) - _TWO_PI) < 1e-6,
                })
            loop_info.sort(key=lambda L: -((L["phiMax"] - L["phiMin"])
                                           * (L["otherMax"] - L["otherMin"] + 1e-30)))
            entry["loops"] = loop_info
            n_wrap = sum(1 for L in loop_info if L["wraps"])
            wraps = n_wrap > 0
            entry["nWrappingLoops"] = n_wrap
            entry["wraps2pi"] = wraps
            # A merged face that wraps the full 2pi has TWO boundary loops that are not nested --
            # the two rims of a closed band. Neither is a hole, and the shipped wire block, which
            # demands exactly one outer wire, cannot express the pair. A non-wrapping merged face
            # has one outer loop and (nLoops - 1) genuine holes, which the format does express as
            # inner wires.
            entry["fullBand"] = bool(n_wrap == 2 and len(loop_info) == 2)
            entry["nHoles"] = (0 if wraps else max(0, len(loop_info) - 1))
            # rectangle test: exactly two distinct `other` levels among iso-other edges, and
            # exactly two distinct phi levels among iso-phi edges (or none, on a full wrap)
            phi_levels = cluster([float(samples[e][0].mean() % _TWO_PI)
                                  for e in boundary if kinds.get(e) == "phi"], 1e-5)
            oth_levels = cluster([float(samples[e][1].mean())
                                  for e in boundary if kinds.get(e) == "other"], max(tol_other, 1e-9))
            entry["nPhiLevels"] = len(phi_levels)
            entry["nOtherLevels"] = len(oth_levels)
            entry["isRectangle"] = bool(
                entry["allBoundaryIso"] and entry["nHoles"] == 0
                and len(oth_levels) == 2
                and (len(phi_levels) == 2 or (wraps and len(phi_levels) == 0)))
            # what the shipped wire block additionally demands
            entry["nOuterLoops"] = 2 if entry["fullBand"] else 1
            entry["oneOuterLoop"] = not entry["fullBand"]
            entry["minLoopEdges"] = min(L["nEdges"] for L in loop_info)
            outer = loop_info[0]
            sweep = _TWO_PI if wraps else (outer["phiMax"] - outer["phiMin"])
            entry["phiSweep"] = sweep
            entry["phiSweepOK"] = (0.0 < sweep <= _TWO_PI + 1e-9)
        merged.append(entry)
    out["merged"] = merged

    # ---- self-check 1: the per-edge iso verdict must be the shipped one -----------------------
    mism = 0
    checked = 0
    for i in sorted(shipped_edges):
        project, _on = projector_for_model(models[i], s)
        if project is None:
            continue
        # walk the face's wires in `_face_wire_edges` order, which is exactly the order
        # `wire_block_census` appends its entries in, so the two zip edge for edge.
        seq, allo, ok = [], [], True
        for _w, _outer, wedges in C._face_wire_edges(faces[i]):
            for e, sv in wedges:
                if BRep_Tool.Degenerated(e):
                    _p, o = project(BRep_Tool.Pnt(sv))
                    seq.append(None)
                    allo.append(o)
                    continue
                r = edge_uv_samples(e, project)
                if r is None:
                    ok = False
                    break
                seq.append(r)
                allo.extend(r[1].tolist())
            if not ok:
                break
        entries = shipped_edges[i]
        if not ok or len(seq) != len(entries) or not allo:
            continue
        tol_other = 1e-6 * max(1.0, float(max(allo) - min(allo)))
        for r, e in zip(seq, entries):
            if r is None or e["degenerate"]:
                continue
            iso_o, iso_p, _so, _sp = iso_verdict(r[0], r[1], tol_other)
            checked += 1
            if bool(iso_o or iso_p) != bool(e["accepted"]):
                mism += 1
    out["selfCheckEdges"] = checked
    out["selfCheckMismatches"] = mism
    return out


# ---------------------------------------------------------------------------------------------

def run_model(step_path, label, args):
    unit = C.detect_step_length_unit(step_path)
    scale = C.step_unit_scale_to_cm(unit)
    C.extract_graph(step_path, meshparam=None, scale_to_cm=scale)
    solids = []
    only = set(args.solids.split(",")) if args.solids else None
    for lid, shape in C.def_shapes.items():
        name = C.def_names.get(lid, lid)
        if only and name not in only:
            continue
        if args.verbose:
            print(f"  [{label}] {name}")
        solids.append(measure_solid(lid, name, shape, scale, args))
    return {"model": step_path, "label": label, "scaleToCm": scale, "solids": solids}


def build_fixtures(tmpdir):
    subprocess.run([sys.executable, str(_HERE.parent / "make_boolean_fixtures.py"),
                    "--outdir", str(tmpdir)], check=True)
    return sorted(Path(tmpdir).glob("*.step")) + sorted(Path(tmpdir).glob("*.stp"))


def synthetic_control(args):
    """Positive + negative control on the grouping, on geometry built in this process.

    A cylinder is split into four NURBS patches the way an exporter does (four quadrant faces of
    one cylinder, converted to B-spline surfaces so the recogniser, not the stored type, has to
    decide), and a second cylinder of a different radius sits beside it. The four must group into
    one; the fifth must not join them.
    """
    from OCC.Core.gp import gp_Ax2, gp_Pnt, gp_Dir
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert

    fails = 0
    checks = 0

    def check(cond, msg, detail=""):
        nonlocal fails, checks
        checks += 1
        if cond:
            print(f"  [ok ] {msg}{(' -- ' + detail) if detail else ''}")
        else:
            fails += 1
            print(f"  [FAIL] {msg}{(' -- ' + detail) if detail else ''}")

    quarters = []
    for q in range(4):
        a = q * math.pi / 2.0
        axq = gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1),
                     gp_Dir(math.cos(a), math.sin(a), 0.0))
        quarters.append(BRepPrimAPI_MakeCylinder(axq, 3.0, 5.0, math.pi / 2.0).Shape())
    other = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(20, 0, 0), gp_Dir(0, 0, 1)),
                                     3.0, 5.0).Shape()

    # collect the four lateral faces (NURBS-converted) plus one lateral face of the other cylinder
    def lateral_faces(sh):
        out = []
        n = BRepBuilderAPI_NurbsConvert(sh, True).Shape()
        for f in TopologyExplorer(n).faces():
            a = BRepAdaptor_Surface(f)
            kind = C._SURFACE_TYPE_NAMES.get(a.GetType(), "unknown")
            try:
                uv = breptools.UVBounds(f)
                rec = C._recognize_analytic_surface(a, uv)
            except Exception:
                rec = None
            if rec is not None and rec["kind"] == "cylinder":
                out.append((f, kind, rec))
        return out

    fs = []
    for q in quarters:
        fs.extend(lateral_faces(q))
    fs_other = lateral_faces(other)
    check(len(fs) == 4, "the split cylinder really does produce four separately-stored patches",
          f"found {len(fs)}, stored types {sorted({k for _f, k, _r in fs})}")
    check(len(fs_other) >= 1, "the control cylinder produces a recognisable lateral face",
          f"found {len(fs_other)}")
    if len(fs) != 4 or not fs_other:
        return fails + 1, checks

    allf = [f for f, _k, _r in fs] + [fs_other[0][0]]
    models = [T.recognized_model(r) for _f, _k, r in fs] + [T.recognized_model(fs_other[0][2])]
    P = [face_surface_points(f, n=args.surface_grid) for f in allf]
    tol = [float(BRep_Tool.Tolerance(f)) for f in allf]

    seps = {}
    for i in range(len(allf)):
        for j in range(i + 1, len(allf)):
            seps[(i, j)] = surface_separation(models[i], P[i], models[j], P[j])
    worst_same = max(seps[(i, j)] for i in range(4) for j in range(i + 1, 4))
    best_diff = min(seps[(i, 4)] for i in range(4))
    thr = args.merge_factor * max(tol)
    check(worst_same <= thr,
          "the four patches of ONE cylinder are recognised as the same surface",
          f"worst separation {worst_same:.3e} cm vs threshold {thr:.3e} cm")
    check(best_diff > thr,
          "a genuinely DIFFERENT cylinder is not merged with them (negative control)",
          f"nearest separation {best_diff:.3e} cm vs threshold {thr:.3e} cm")
    check(best_diff / max(worst_same, 1e-18) > 1e3,
          "the two populations are separated by orders of magnitude, not by the threshold",
          f"ratio {best_diff / max(worst_same, 1e-18):.3e}")
    return fails, checks


# ---------------------------------------------------------------------------------------------

def summarise(results, args):
    """The four numbers the brief asks for, over every model in `results`."""
    lines = []
    P = lines.append

    for res in results:
        label = res["label"]
        solids = res["solids"]
        P(f"\n===== {label} =====")

        # ---- self-check
        mm = sum(s["selfCheckMismatches"] for s in solids)
        ck = sum(s["selfCheckEdges"] for s in solids)
        P(f"self-check 1 (per-edge iso verdict vs the shipped frame): "
          f"{mm} mismatches over {ck} edges")

        # ---- 1. fragmentation
        hist = Counter()
        worst = []
        for s in solids:
            for k, v in s["groupSizeHistogram"].items():
                hist[int(k)] += v
            for g in s["fragmentedGroups"]:
                worst.append((len(g["faces"]), s["name"], g["kind"], g["faces"][:6]))
        n_surf = sum(s["nSurfaces"] for s in solids)
        n_faces = sum(s["nFacesWithAnalyticSurface"] for s in solids)
        n_mg = sum(s["nMergeGroups"] for s in solids)
        mhist = Counter()
        for s in solids:
            for k, v in s["mergeGroupSizeHistogram"].items():
                mhist[int(k)] += v
        P(f"\n1. FRAGMENTATION: {n_faces} analytic faces on {n_surf} distinct surfaces "
          f"({n_faces / n_surf if n_surf else 0:.2f} faces per surface)")
        P(f"   group-size histogram (faces per surface): "
          f"{dict(sorted(hist.items()))}")
        P(f"   ... of which ADJACENT co-surface components (what a merge may fuse): {n_mg}, "
          f"histogram {dict(sorted(mhist.items()))}")
        P(f"   surfaces carrying more than one disconnected patch set: "
          f"{sum(s['nDisconnectedSurfaces'] for s in solids)}")
        worst.sort(reverse=True)
        for n, sname, kind, fl in worst[:8]:
            P(f"   worst: {sname:18s} {kind:9s} {n:3d} faces  e.g. {fl}")

        # ---- separation calibration
        sm = [s["mergedSeparationOverTol"]["max"] for s in solids
              if s["mergedSeparationOverTol"]["max"] is not None]
        ss = [s["nearestSplitSeparationOverTol"] for s in solids
              if s["nearestSplitSeparationOverTol"] is not None]
        smx = max(sm) if sm else 0.0
        ssn = min(ss) if ss else float("inf")
        P(f"   merged pairs: max separation {smx:.3g} declared tolerances")
        P(f"   nearest NON-merged same-kind pair: {ssn:.3g} declared tolerances")

        # ---- 2. seam vs true boundary
        rej = [e for s in solids for e in s["edges"] if not e["accepted"] and not e["degenerate"]]
        acc = [e for s in solids for e in s["edges"] if e["accepted"] and not e["degenerate"]]
        cr = Counter(e["class"] for e in rej)
        ca = Counter(e["class"] for e in acc)
        P(f"\n2. REJECTED boundary edges: {len(rej)}   {dict(cr)}")
        P(f"   accepted (control):       {len(acc)}   {dict(ca)}")
        cc = [e for e in rej + acc if e.get("coincidentVerdict") is not None]
        if cc:
            agree = sum(1 for e in cc
                        if e["coincidentVerdict"] == (e["class"] == "seam"))
            seams = [e["coincidentMaxAbsFCm"] for e in cc if e["class"] == "seam"]
            trues = [e["coincidentMaxAbsFCm"] for e in cc if e["class"] != "seam"]
            P(f"   cross-check vs Stream_R's coincidence instrument: "
              f"{agree}/{len(cc)} agree ({len(cc) - agree} disagree)")
            P(f"     max |f_neighbour| on the face's own surface: "
              f"seam edges max {(max(seams) if seams else 0.0):.3g} cm, "
              f"true-boundary edges min {(min(trues) if trues else float('inf')):.3g} cm")

        # ---- 3 + 4. merge simulation
        P("\n3. THE DECISIVE NUMBER")
        n_emit_now = sum(1 for s in solids if s["emitted"])
        cand = [s for s in solids if not s["emitted"] and s["nFaces"]
                and s["nUnsupportedFaces"] == 0 and s["nTrimDeclinedFaces"]]
        P(f"   solids {len(solids)}, emitting today {n_emit_now}, "
          f"trim-blocked candidates {len(cand)}")
        for s in solids:
            for m in s["merged"]:
                tag = f"f{m['faces'][:3]}"
                b = []
                if m.get("status"):
                    b.append(f"{tag} {m['status']}")
                else:
                    if not m.get("allBoundaryIso", False):
                        b.append(f"{tag} {m['nBoundaryNonIso']} non-iso edges")
                    if m.get("loopProblem"):
                        b.append(f"{tag} {m['loopProblem']}")
                    if m.get("oneOuterLoop") is False:
                        b.append(f"{tag} full 2pi band: 2 boundary loops")
                    if m.get("minLoopEdges", 3) < 3:
                        b.append(f"{tag} a loop has < 3 edges")
                    if m.get("phiSweepOK") is False:
                        b.append(f"{tag} phi sweep {m.get('phiSweep'):.3f}")
                m["blockers"] = b
                m["isoOK"] = bool(m.get("allBoundaryIso", False))
                m["wireOK"] = (not b)
                # a merged face whose (phi, other) domain is a plain RECTANGLE needs no trim wire
                # at all: it is the parameter-only quadric record `extract_cylindrical_face` has
                # always emitted for a natively-analytic full or partial cylinder. Counted
                # separately because it is a different *existing* path, not a new one.
                m["rectOK"] = bool(m["isoOK"] and m.get("isRectangle"))
                m["recordOK"] = bool(m["wireOK"] or m["rectOK"])
            for policy in ("all", "need"):
                sel = [m for m in s["merged"]
                       if policy == "all" or m["needsMerge"]]
                s[f"blockers_{policy}"] = [x for m in sel for x in m["blockers"]]
                s[f"isoOK_{policy}"] = all(m["isoOK"] for m in sel)
                s[f"pathOK_{policy}"] = (all(m["wireOK"] for m in sel)
                                         and s["nUnsupportedFaces"] == 0)
                s[f"recordOK_{policy}"] = (all(m["recordOK"] for m in sel)
                                           and s["nUnsupportedFaces"] == 0)
        for pol in ("need", "all"):
            n_iso = sum(1 for s in cand if s[f"isoOK_{pol}"])
            n_path = sum(1 for s in cand if s[f"pathOK_{pol}"])
            n_rec = sum(1 for s in cand if s[f"recordOK_{pol}"])
            P(f"   policy '{pol}': every remaining boundary edge iso on {n_iso}/{len(cand)}; "
              f"wire block unchanged accepts {n_path}/{len(cand)}; "
              f"expressible as an EXISTING record (wire block or rectangle params) "
              f"{n_rec}/{len(cand)}")
        for s in cand:
            mark = ("EMITS" if s["pathOK_need"] else
                    ("rect" if s["recordOK_need"] else
                     ("iso-ok" if s["isoOK_need"] else "blocked")))
            P(f"   {s['name']:18s} {mark:7s} "
              f"{sum(1 for m in s['merged'] if m['nMembers'] > 1 and m['needsMerge']):3d}/"
              f"{sum(1 for m in s['merged'] if m['nMembers'] > 1):3d} groups merged (need/all), "
              f"{sum(m.get('nInternalDissolved', 0) for m in s['merged'] if m['needsMerge']):4d} "
              f"seams dissolved"
              + (f"   blocked: {s['blockers_need'][:3]}" if s["blockers_need"] else ""))

        P("\n4. MERGED-DOMAIN SANITY")
        allm = [m for s in solids for m in s["merged"] if m["nMembers"] > 1]
        rect = sum(1 for m in allm if m.get("isRectangle"))
        wrap = sum(1 for m in allm if m.get("wraps2pi"))
        band = sum(1 for m in allm if m.get("fullBand"))
        hole = sum(1 for m in allm if m.get("nHoles", 0) > 0)
        prob = sum(1 for m in allm if m.get("loopProblem"))
        P(f"   merged groups {len(allm)}: rectangle {rect}, wraps 2pi {wrap} "
          f"(full band, 2 outer wires {band}), with holes {hole}, "
          f"boundary not a clean loop set {prob}")
        bad2 = sum(m.get("nDissolvedNotUsedTwice", 0) for m in allm)
        bad1 = sum(m.get("nBoundaryUsedTwiceInside", 0) for m in allm)
        diss = sum(m.get("nInternalDissolved", 0) for m in allm)
        P(f"   edge identity: {diss} seams dissolved, of which {bad2} were NOT used exactly "
          f"twice inside their group; {bad1} surviving boundary edges used more than once "
          f"inside it. Both must be 0 for closure to be balanced by construction.")

    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", action="append", default=[])
    ap.add_argument("--fixtures", action="store_true")
    ap.add_argument("--solids", help="comma-separated leaf-solid names to restrict to")
    ap.add_argument("--json")
    ap.add_argument("--merge-factor", type=float, default=_MERGE_FACTOR,
                    help="declared BRep tolerances two faces' surfaces may differ by and still be "
                         "merged (default %(default)s)")
    ap.add_argument("--surface-grid", type=int, default=5,
                    help="n x n points per face used for the surface-identity decision")
    ap.add_argument("--coincidence-check", action="store_true",
                    help="cross-check every seam / true-boundary verdict against Stream_R's "
                         "independent test: is the neighbour's implicit function identically "
                         "zero over the face's own surface samples?")
    ap.add_argument("--negative-control", action="store_true",
                    help="report the nearest same-kind pair that was NOT merged, per solid")
    ap.add_argument("--synthetic", action="store_true",
                    help="run the deliberately-fragmented-cylinder control and exit")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.synthetic:
        print("Synthetic control: a cylinder deliberately fragmented into four patches")
        fails, checks = synthetic_control(args)
        print(f"\n{checks} checks, {fails} failure(s)")
        return 1 if fails else 0

    jobs = []
    if args.fixtures:
        outdir = Path(tempfile.mkdtemp(prefix="cosurface-fixtures-"))
        for f in build_fixtures(outdir):
            jobs.append((f"fixture:{f.stem}", str(f)))
    for m in args.model:
        if ":" in m and not os.path.exists(m):
            label, path = m.split(":", 1)
        else:
            label, path = Path(m).stem, m
        jobs.append((label, path))

    results = []
    for label, path in jobs:
        print(f"== {label}: {path}", flush=True)
        results.append(run_model(path, label, args))

    print(summarise(results, args))

    if args.json:
        out = Path(args.json)
        out.parent.mkdir(parents=True, exist_ok=True)

        def clean(o):
            if isinstance(o, dict):
                return {k: clean(v) for k, v in o.items() if k not in ("edge",)}
            if isinstance(o, (list, tuple)):
                return [clean(v) for v in o]
            if isinstance(o, (np.floating, np.integer, np.bool_)):
                return o.item()
            if isinstance(o, set):
                return sorted(o)
            if isinstance(o, float) and not math.isfinite(o):
                return None
            return o
        out.write_text(json.dumps(clean(results), indent=1))
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
