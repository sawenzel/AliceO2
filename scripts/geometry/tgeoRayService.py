#!/usr/bin/env python3
"""A local ray-tracing bridge to the real O2 geometry kernel.

Serves the website's RemoteEngine: rays in, kernel answers out, so the browser view is the
actual O2BVHSurfaceSolid (or any TGeoShape) rather than a port of it.

  POST /load   JSON {"path": "<surfaces_*.bin or shape_*.root>"}
               -> {"ok": true, "kind": "surface|shape", "nSurfaces": N, "bbox": [...]}
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
#include "DetectorsBase/O2SurfaceSolidIO.h"
#include "TGeoShape.h"
#include "TFile.h"
#include <thread>
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

void trace(const TGeoShape* shape, const float* rays, int n, float* out, int nThreads)
{
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
    for fn in (ROOT.raysvc.trace, ROOT.raysvc.loadSurface, ROOT.raysvc.loadShapeFile):
        try:
            fn.__release_gil__ = True
        except AttributeError:
            pass


class State:
    shape = None
    kind = None
    root = "."


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
                if not os.path.exists(path):
                    return self._reply_json({"ok": False, "error": f"no such file: {path}"}, 422)
                if path.endswith(".root"):
                    shape, kind = ROOT.raysvc.loadShapeFile(path), "shape"
                else:
                    shape, kind = ROOT.raysvc.loadSurface(path), "surface"
                if not shape:
                    return self._reply_json({"ok": False, "error": f"cannot load {path}"}, 422)
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
                if kind == "surface":
                    info["nSurfaces"] = state.shape.GetNsurfaces()
                    info["reliability"] = state.shape.GetNavigationReliabilityName(
                        state.shape.GetNavigationReliability())
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
    state = State()
    state.root = args.root
    if args.load:
        loader = ROOT.raysvc.loadShapeFile if args.load.endswith(".root") else ROOT.raysvc.loadSurface
        state.shape = loader(args.load)
        state.kind = "shape" if args.load.endswith(".root") else "surface"
        if not state.shape:
            print(f"cannot load {args.load}", file=sys.stderr)
            return 1
        print(f"loaded {args.load}")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(state, args.threads))
    print(f"tgeoRayService on 127.0.0.1:{args.port}  (threads={args.threads})")
    server.serve_forever()


if __name__ == "__main__":
    sys.exit(main() or 0)
