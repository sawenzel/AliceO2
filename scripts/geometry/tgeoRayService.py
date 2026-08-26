#!/usr/bin/env python3
"""A local ray-tracing bridge to the real O2 geometry kernel.

Serves the website's RemoteEngine: rays in, kernel answers out, so the browser view is the
actual O2BVHSurfaceSolid (or any TGeoShape) rather than a port of it.

  POST /load   JSON {"path": "<surfaces_*.bin, flatcsg_*.bin, facets_*.bin or shape_*.root>"}
               -> {"ok": true, "kind": "surface|flatcsg|mesh|shape", "bbox": [...], ...}

               All four representations one converted part can ship as, so the same part can be
               traced and benchmarked in each of them by the real kernel. A .root file is a
               TGeoShape; a .bin file is dispatched on its own magic ("O2SS" -> the exact surface
               sidecar, "O2FLTCSG" -> the flat-CSG sidecar, neither -> the facet mesh, which has
               no magic and begins with its triangle count).
  POST /bench  JSON {"samples": N, "bbox": [x0,y0,z0,x1,y1,z1], "seed": S}
               -> per-call timings for Contains / DistFromOutside / DistFromInside / Safety.
                  `bbox` and `seed` fix the sample points, so the SAME points can be put to every
                  representation of one part; without them the loaded shape's own box is used.

  POST /trace  raw Float32Array, 6 floats per ray (origin, unit direction, cm)
               -> raw Float32Array, 5 floats per ray: (t, nx, ny, nz, startedInside);
                  t < 0 means no hit

Run inside the O2 environment (needs pyROOT + libO2DetectorsBase):

  source <env with O2/latest-swenzel-bvhsurfacesolid-o2>
  python3 scripts/geometry/tgeoRayService.py --port 8077

The tracer loop is JIT-compiled C++, split over threads; the per-ray cost is the kernel's own
(measured microseconds for the surface solid, nanoseconds for composites).
"""
import argparse
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np
import ROOT

CPP = r"""
#include "DetectorsBase/O2BVHSurfaceSolid.h"
#include "DetectorsBase/O2FlatCSG.h"
#include "DetectorsBase/O2Tessellated.h"
#include "DetectorsBase/O2SurfaceSolidIO.h"
#include "TGeoShape.h"
#include "TFile.h"
#include <thread>
#include <chrono>
#include <vector>

namespace raysvc {

o2::base::O2BVHSurfaceSolid* loadSurface(const char* path)
{
  auto* solid = new o2::base::O2BVHSurfaceSolid("raysvc");
  if (!o2::base::LoadSurfaceSolid(path, *solid)) {
    delete solid;
    return nullptr;
  }
  solid->CloseShape();
  return solid;
}

o2::base::O2FlatCSG* loadFlatCSG(const char* path)
{
  auto* solid = new o2::base::O2FlatCSG("raysvcflat");
  if (!o2::base::LoadFlatCSG(path, *solid)) {
    delete solid;
    return nullptr;
  }
  solid->CloseShape();
  return solid;
}

o2::base::O2Tessellated* loadFacets(const char* path)
{
  auto* solid = new o2::base::O2Tessellated("raysvcmesh");
  if (!o2::base::LoadFacetSolid(path, *solid)) {
    delete solid;
    return nullptr;
  }
  // The mesh is a converter artefact, not a hand-written one: check and fix flipped facets, but
  // do not print a report per load.
  solid->CloseShape(true, true, false);
  return solid;
}

TGeoShape* loadShapeFile(const char* path)
{
  TFile* file = TFile::Open(path, "READ");
  if (file == nullptr || file->IsZombie()) {
    return nullptr;
  }
  auto* shape = file->Get<TGeoShape>("shape");
  return shape; // the file stays open on purpose; the shape lives in it
}

void traceChunk(const TGeoShape* shape, const float* rays, int begin, int end, float* out)
{
  for (int i = begin; i < end; ++i) {
    const float* r = rays + 6 * i;
    float* o = out + 5 * i;
    const double p[3] = {r[0], r[1], r[2]};
    const double d[3] = {r[3], r[4], r[5]};
    const bool inside = const_cast<TGeoShape*>(shape)->Contains(p);
    const double t = inside ? const_cast<TGeoShape*>(shape)->DistFromInside(p, d, 3, 1e30, nullptr)
                            : const_cast<TGeoShape*>(shape)->DistFromOutside(p, d, 3, 1e30, nullptr);
    if (!(t < 1e29)) {
      o[0] = -1.f;
      o[1] = o[2] = o[3] = 0.f;
      o[4] = inside ? 1.f : 0.f;
      continue;
    }
    const double hp[3] = {p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2]};
    double n[3] = {0., 0., 0.};
    const_cast<TGeoShape*>(shape)->ComputeNormal(hp, const_cast<double*>(d), n);
    o[0] = static_cast<float>(t);
    o[1] = static_cast<float>(n[0]);
    o[2] = static_cast<float>(n[1]);
    o[3] = static_cast<float>(n[2]);
    o[4] = inside ? 1.f : 0.f;
  }
}

// Is this shape safe to query from several threads at once?
//
// A TGeoCompositeShape is NOT. Every TGeoBoolNode keeps per-thread scratch state -- which operand
// the last DistFrom* selected, and cached operand points -- behind TGeoBoolNode::GetThreadData(),
// which indexes an array by TGeoManager::ThreadId(). ThreadId() returns 0 for EVERY thread unless
// the manager is multi-threaded, and a manager only becomes multi-threaded through
// SetMaxThreads(), which refuses before the geometry is closed. This service has no geometry to
// close -- it navigates bare shapes -- so every worker here would be thread 0 and they would all
// write over one another's scratch state.
//
// That is not a theoretical race. Measured on ITS BREF1 with 8 threads, the composite MISSED 13033
// of 35420 hits its own exact solid found, reported phantom hits, and returned surface normals up
// to 90 degrees wrong on ~19% of pixels -- all of it gone, exactly, at one thread. Registering the
// threads instead is not an option here: ThreadId() hands out ids from a counter that is never
// recycled, this pool creates fresh std::threads per band, and TGeoBoolNode::GetThreadData() has
// its bounds check commented out, so the ids would run past the array and read out of bounds.
//
// So a boolean shape is traced on one thread. It costs nothing that matters: composites are the
// fast subjects, and the slow ones -- the surface solid, the flat solid, the mesh -- hold no
// mutable per-query state and stay parallel.
bool isThreadSafe(const TGeoShape* shape)
{
  return shape != nullptr && !shape->InheritsFrom("TGeoCompositeShape");
}

void trace(const TGeoShape* shape, const float* rays, int n, float* out, int nThreads)
{
  if (!isThreadSafe(shape)) {
    nThreads = 1;
  }
  if (nThreads <= 1 || n < 512) {
    traceChunk(shape, rays, 0, n, out);
    return;
  }
  std::vector<std::thread> pool;
  const int chunk = (n + nThreads - 1) / nThreads;
  for (int t = 0; t < nThreads; ++t) {
    const int begin = t * chunk;
    const int end = std::min(n, begin + chunk);
    if (begin < end) {
      pool.emplace_back(traceChunk, shape, rays, begin, end, out);
    }
  }
  for (auto& th : pool) {
    th.join();
  }
}

double benchFunction(const TGeoShape* shape, int which, const float* origins, const float* dirs,
                     int n, int repeats)
{
  // which: 0 Contains, 1 DistFromOutside, 2 DistFromInside, 3 Safety. Single-threaded on
  // purpose: the number is ns per call, not throughput.
  auto* sh = const_cast<TGeoShape*>(shape);
  volatile double sink = 0.;
  const auto start = std::chrono::steady_clock::now();
  for (int r = 0; r < repeats; ++r) {
    for (int i = 0; i < n; ++i) {
      const double p[3] = {origins[3 * i], origins[3 * i + 1], origins[3 * i + 2]};
      const double d[3] = {dirs[3 * i], dirs[3 * i + 1], dirs[3 * i + 2]};
      switch (which) {
        case 0: sink = sink + (sh->Contains(p) ? 1. : 0.); break;
        case 1: sink = sink + sh->DistFromOutside(p, d, 3, 1e30, nullptr); break;
        case 2: sink = sink + sh->DistFromInside(p, d, 3, 1e30, nullptr); break;
        default: sink = sink + sh->Safety(p, false); break;
      }
    }
  }
  const std::chrono::duration<double, std::nano> elapsed = std::chrono::steady_clock::now() - start;
  return elapsed.count() / (double(n) * repeats);
}

} // namespace raysvc
"""


def jit_setup(o2_src: str):
    ROOT.gSystem.Load("libO2DetectorsBase")
    ROOT.gInterpreter.AddIncludePath(os.path.join(o2_src, "Detectors/Base/include"))
    if not ROOT.gInterpreter.Declare(CPP):
        raise RuntimeError("JIT compilation of the tracer failed")
    # Release the GIL around the long C++ calls: a handler thread that holds it through a
    # multi-second trace starves every other request, and a worker thread that needs the
    # interpreter for lazy symbol resolution while the caller holds the GIL deadlocks the
    # whole server (observed on the first threaded trace of a surface solid).
    for fn in (ROOT.raysvc.trace, ROOT.raysvc.loadSurface, ROOT.raysvc.loadFlatCSG,
               ROOT.raysvc.loadFacets, ROOT.raysvc.loadShapeFile, ROOT.raysvc.benchFunction):
        try:
            fn.__release_gil__ = True
        except AttributeError:
            pass


# The four representations one converted part can ship as. A .root file is a streamed TGeoShape;
# a .bin sidecar says which one it is in its own first bytes, so a caller never has to encode the
# kind in the path it sends -- the website simply posts whatever artefact the manifest names.
SIDECAR_MAGIC = (("O2FLTCSG", "flatcsg"), ("O2SS", "surface"))


def sidecar_kind(path):
    """Which representation `path` holds, from its extension and its magic."""
    if path.endswith(".root"):
        return "shape"
    with open(path, "rb") as handle:
        head = handle.read(8)
    for magic, kind in SIDECAR_MAGIC:
        if head.startswith(magic.encode()):
            return kind
    # facets_*.bin carries no magic at all: it opens with its uint32 triangle count.
    return "mesh"


def load_shape(kind, path):
    """The kernel loader for one kind, called on the JIT-compiled C++ side."""
    return {"shape": ROOT.raysvc.loadShapeFile, "surface": ROOT.raysvc.loadSurface,
            "flatcsg": ROOT.raysvc.loadFlatCSG, "mesh": ROOT.raysvc.loadFacets}[kind](path)


def boolean_shape_size(shape):
    """`(depth, leaves)` of a TGeo boolean tree; a shape that is not one is `(0, 1)`."""
    if not shape.InheritsFrom("TGeoCompositeShape"):
        return 0, 1
    node = shape.GetBoolNode()
    left = boolean_shape_size(node.GetLeftShape())
    right = boolean_shape_size(node.GetRightShape())
    return max(left[0], right[0]) + 1, left[1] + right[1]


def shape_info(kind, shape):
    """What this representation can say about itself, beyond its bounding box: the counts the
    part card prints. Every number is read off the loaded shape, never off the file."""
    if kind == "surface":
        return {"nSurfaces": shape.GetNsurfaces(),
                "reliability": shape.GetNavigationReliabilityName(shape.GetNavigationReliability())}
    if kind == "flatcsg":
        return {"nCells": shape.GetNcells(), "nHalfspaces": shape.GetNhalfspaces(),
                "nBoxes": shape.GetNboxes(), "bvhBytes": int(shape.GetBVHMemory())}
    if kind == "mesh":
        return {"nFacets": shape.GetNfacets(), "nVertices": shape.GetNvertices()}
    info = {"className": shape.ClassName()}
    # A .root may hold a boolean tree, and how that tree is SHAPED is the thing a reader wants to
    # know about it -- ITS's IBCYSSFlangeC is 36 leaves at recursion depth 35, because it is a
    # chain of 33 subtractions and not a balanced tree, and a navigator descends that depth on
    # every query. Reporting the class alone hides exactly that.
    depth, leaves = boolean_shape_size(shape)
    if depth:
        info["booleanDepth"] = depth
        info["leaves"] = leaves
        # See raysvc::isThreadSafe: a TGeoBoolNode's scratch state is shared between this
        # service's workers, so a boolean shape is traced on one thread and the page says so
        # rather than quietly rendering a corrupted picture.
        info["singleThreaded"] = True
    return info


class State:
    shape = None
    kind = None
    root = "."
    cache = {}        # path -> (shape, kind): shapes register themselves in the TGeoManager and
                      # are NEVER replaced or deleted; a reload is a cache hit and a pointer swap
    lock = None       # serialises /load; traces read state.shape once and keep their pointer


def make_handler(state: State, threads: int):
    class Handler(BaseHTTPRequestHandler):
        def _headers(self, code, ctype, length):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(length))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()

        def do_OPTIONS(self):  # CORS preflight for binary POSTs from another localhost port
            self._headers(204, "text/plain", 0)

        def _reply_json(self, obj, code=200):
            body = json.dumps(obj).encode()
            self._headers(code, "application/json", len(body))
            self.wfile.write(body)

        def do_POST(self):
            body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
            if self.path == "/load":
                try:
                    path = json.loads(body)["path"]
                except Exception:
                    return self._reply_json({"ok": False, "error": "bad request"}, 400)
                # the website sends its own relative testdata paths; resolve them against --root
                if not os.path.isabs(path):
                    path = os.path.join(state.root, path)
                path = os.path.realpath(path)
                if not os.path.exists(path):
                    return self._reply_json({"ok": False, "error": f"no such file: {path}"}, 422)
                with state.lock:
                    if path in state.cache:
                        shape, kind = state.cache[path]
                    else:
                        kind = sidecar_kind(path)
                        shape = load_shape(kind, path)
                        if not shape:
                            return self._reply_json({"ok": False, "error": f"cannot load {path}"}, 422)
                        state.cache[path] = (shape, kind)
                    state.shape, state.kind = shape, kind
                # Warm-up: one single-threaded ray, so every lazy-JIT symbol of this shape's
                # navigation path is resolved on this thread before any std::thread runs it.
                warm_rays = np.zeros(6, dtype=np.float32); warm_rays[3] = 1.0
                warm_out = np.empty(5, dtype=np.float32)
                ROOT.raysvc.trace(state.shape, warm_rays, 1, warm_out, 1)
                origin = [shape.GetOrigin()[i] for i in range(3)]
                half = [shape.GetDX(), shape.GetDY(), shape.GetDZ()]
                info = {"ok": True, "kind": kind,
                        "bbox": [origin[0] - half[0], origin[1] - half[1], origin[2] - half[2],
                                 origin[0] + half[0], origin[1] + half[1], origin[2] + half[2]]}
                info.update(shape_info(kind, state.shape))
                return self._reply_json(info)
            if self.path == "/trace":
                if state.shape is None:
                    return self._reply_json({"ok": False, "error": "no shape loaded"}, 409)
                rays = np.frombuffer(body, dtype=np.float32)
                if rays.size % 6 != 0:
                    return self._reply_json({"ok": False, "error": "ray buffer not n*6 floats"}, 400)
                n = rays.size // 6
                out = np.empty(n * 5, dtype=np.float32)
                ROOT.raysvc.trace(state.shape, rays, n, out, threads)
                payload = out.tobytes()
                self._headers(200, "application/octet-stream", len(payload))
                self.wfile.write(payload)
                return
            if self.path == "/bench":
                if state.shape is None:
                    return self._reply_json({"ok": False, "error": "no shape loaded"}, 409)
                try:
                    cfg = json.loads(body) if body else {}
                except Exception:
                    cfg = {}
                n = int(cfg.get("samples", 4000))
                n = max(100, min(n, 50000))
                sh = state.shape
                # The sample box may be given, and for a comparison it MUST be: /bench is called
                # once per representation of the same part, and the representations do not share a
                # bounding box -- a TGeoCompositeShape's is the union of its padded leaf boxes and
                # is wider than the flat solid's or the sidecar's. Drawing each one's points from
                # its own box would time four kernels on four different point sets and call the
                # result a comparison. The caller therefore passes one box for all of them;
                # falling back to the loaded shape's own box keeps a single /bench honest.
                box = cfg.get("bbox")
                if isinstance(box, (list, tuple)) and len(box) == 6 and all(
                        isinstance(v, (int, float)) for v in box):
                    o = [0.5 * (box[i] + box[i + 3]) for i in range(3)]
                    h = [max(1e-9, 0.5 * abs(box[i + 3] - box[i])) for i in range(3)]
                    box_source = "caller"
                else:
                    o = [sh.GetOrigin()[i] for i in range(3)]
                    h = [sh.GetDX(), sh.GetDY(), sh.GetDZ()]
                    box_source = "loaded shape"
                rng = np.random.default_rng(int(cfg.get("seed", 20260822)))
                # points spread over 1.5x the bbox: a mix of inside and outside
                pts = np.stack([rng.uniform(o[i] - 1.5 * h[i], o[i] + 1.5 * h[i], n)
                                for i in range(3)], axis=1).astype(np.float32)
                dirs = rng.normal(size=(n, 3))
                dirs /= np.linalg.norm(dirs, axis=1)[:, None]
                dirs = dirs.astype(np.float32)
                # classify once through /trace machinery to split inside/outside honestly
                probe = np.empty(n * 5, dtype=np.float32)
                rays = np.concatenate([pts, dirs], axis=1).astype(np.float32).ravel()
                with state.lock:
                    ROOT.raysvc.trace(sh, rays, n, probe, 1)
                    answers = probe.reshape(-1, 5)
                    started_inside = answers[:, 4] > 0.5
                    inside_pts = [pts[started_inside]]
                    inside_dirs = [dirs[started_inside]]
                    harvested = 0
                    # A uniform draw over the bounding box of a thin part -- a frame, a rail, a
                    # shell -- lands inside it essentially never, and then DistFromInside is not
                    # measured at all. The probe already knows where every outside ray ENTERS the
                    # solid, so a hair past that entry point is an interior point that cost nothing
                    # extra. Contains is asked about each one and only the confirmed ones are
                    # kept, so a grazing tangent cannot smuggle an outside point in.
                    entered = (~started_inside) & (answers[:, 0] >= 0)
                    if entered.any():
                        span = 2.0 * float(np.hypot(h[0], np.hypot(h[1], h[2])))
                        step = (answers[entered, 0] + 1e-6 * span)[:, None]
                        candidates = (pts[entered] + step * dirs[entered]).astype(np.float32)
                        entered_dirs = dirs[entered]
                        keep = np.fromiter(
                            (bool(sh.Contains(np.asarray(pt, dtype=np.float64))) for pt in candidates),
                            dtype=bool, count=len(candidates))
                        harvested = int(keep.sum())
                        inside_pts.append(candidates[keep])
                        inside_dirs.append(entered_dirs[keep])
                    inside = np.ascontiguousarray(np.concatenate(inside_pts))
                    din = np.ascontiguousarray(np.concatenate(inside_dirs))
                    outside = np.ascontiguousarray(pts[~started_inside])
                    dout = np.ascontiguousarray(dirs[~started_inside])
                    repeats = max(1, 20000 // n)
                    result = {"ok": True, "samples": n, "repeats": repeats,
                              "insideSamples": int(len(inside)),
                              "insideFromEntry": harvested,
                              "bboxSource": box_source,
                              "loadAverage": os.getloadavg()[0],
                              "functions": {}}
                    result["functions"]["contains"] = {
                        "nsPerCall": ROOT.raysvc.benchFunction(sh, 0, pts.ravel(), dirs.ravel(), n, repeats)}
                    if len(outside):
                        result["functions"]["distFromOutside"] = {
                            "nsPerCall": ROOT.raysvc.benchFunction(sh, 1, outside.ravel(), dout.ravel(), len(outside), repeats)}
                    if len(inside):
                        result["functions"]["distFromInside"] = {
                            "nsPerCall": ROOT.raysvc.benchFunction(sh, 2, inside.ravel(), din.ravel(), len(inside), repeats)}
                    result["functions"]["safety"] = {
                        "nsPerCall": ROOT.raysvc.benchFunction(sh, 3, pts.ravel(), dirs.ravel(), n, repeats)}
                return self._reply_json(result)
            self._reply_json({"ok": False, "error": "unknown endpoint"}, 404)

        def log_message(self, fmt, *args):  # quiet
            pass

    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8077)
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--o2-src", default=os.path.expanduser("~/alisw/O2"))
    parser.add_argument("--load", help="optionally load a shape at startup")
    parser.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "website"),
                        help="directory that relative /load paths resolve against (default: the website dir)")
    args = parser.parse_args()

    jit_setup(args.o2_src)
    import threading
    ROOT.TGeoManager("raysvc", "ray service geometry")  # one manager, created once, up front
    state = State()
    state.lock = threading.Lock()
    state.root = args.root
    if args.load:
        state.kind = sidecar_kind(args.load)
        state.shape = load_shape(state.kind, args.load)
        if not state.shape:
            print(f"cannot load {args.load}", file=sys.stderr)
            return 1
        print(f"loaded {args.load}")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(state, args.threads))
    print(f"tgeoRayService on 127.0.0.1:{args.port}  (threads={args.threads})")
    server.serve_forever()


if __name__ == "__main__":
    sys.exit(main() or 0)
