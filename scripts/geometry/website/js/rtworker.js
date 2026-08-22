// The LocalEngine's worker: it parses the sidecar itself and answers batches of rays.
//
// The worker is handed the raw sidecar bytes rather than a built solid, because a SurfaceSolid is
// a graph of classes and cannot be structured-cloned. Parsing costs a millisecond or two per part
// and keeps one code path for building the geometry.

import { parseSidecar, parseFacets } from './sidecar.js';
import { SurfaceSolid, K_RAY_TOLERANCE, BIG } from './solid.js';
import { MeshSolid } from './meshtrace.js';
import { generateEvents } from './gun.js';

let solid = null;
let mesh = null;

function traceExact(rays) {
  const n = rays.length / 6;
  const out = new Float32Array(n * 5);
  for (let i = 0; i < n; ++i) {
    const o = i * 6;
    const hit = solid ? solid.firstHit(rays[o], rays[o + 1], rays[o + 2], rays[o + 3], rays[o + 4], rays[o + 5]) : null;
    const q = i * 5;
    if (!hit) { out[q] = -1; continue; }
    const len = Math.hypot(hit.nx, hit.ny, hit.nz) || 1;
    out[q] = hit.t;
    out[q + 1] = hit.nx / len;
    out[q + 2] = hit.ny / len;
    out[q + 3] = hit.nz / len;
    out[q + 4] = hit.boundary ? 1 : 0;
  }
  return out;
}

function traceMesh(rays) {
  const n = rays.length / 6;
  const out = new Float32Array(n * 5);
  for (let i = 0; i < n; ++i) {
    const o = i * 6;
    const hit = mesh ? mesh.firstHit(rays[o], rays[o + 1], rays[o + 2], rays[o + 3], rays[o + 4], rays[o + 5]) : null;
    const q = i * 5;
    if (!hit) { out[q] = -1; continue; }
    const len = Math.hypot(hit.nx, hit.ny, hit.nz) || 1;
    // A mesh facet's winding need not face the camera, so orient the normal against the ray: the
    // shading is about the surface, not about the facet's bookkeeping.
    const sign = (hit.nx * rays[o + 3] + hit.ny * rays[o + 4] + hit.nz * rays[o + 5]) > 0 ? -1 : 1;
    out[q] = hit.t;
    out[q + 1] = sign * hit.nx / len;
    out[q + 2] = sign * hit.ny / len;
    out[q + 3] = sign * hit.nz / len;
    out[q + 4] = 0;
  }
  return out;
}

function parity(rays, which) {
  const n = rays.length / 6;
  const out = new Int32Array(n);
  for (let i = 0; i < n; ++i) {
    const o = i * 6;
    if (which === 'mesh') {
      out[i] = mesh ? mesh.crossingCount(rays[o], rays[o + 1], rays[o + 2], rays[o + 3], rays[o + 4], rays[o + 5]) : 0;
    } else {
      out[i] = solid ? solid.crossingCount(rays[o], rays[o + 1], rays[o + 2], rays[o + 3], rays[o + 4], rays[o + 5]) : 0;
    }
  }
  return out;
}

self.onmessage = (event) => {
  const message = event.data;
  try {
    if (message.type === 'load') {
      solid = null;
      mesh = null;
      if (message.sidecar) {
        solid = new SurfaceSolid(parseSidecar(message.sidecar, message.label || 'sidecar'), message.label);
      }
      if (message.facets) {
        mesh = new MeshSolid(parseFacets(message.facets, message.label || 'facets').positions, message.label);
      }
      self.postMessage({
        type: 'loaded', id: message.id,
        nSurfaces: solid ? solid.nSurfaces : 0,
        nTriangles: mesh ? mesh.nTriangles : 0,
        aabb: solid ? solid.aabb : (mesh ? mesh.aabb : null),
      });
      return;
    }
    if (message.type === 'trace' || message.type === 'traceMesh') {
      const out = message.type === 'trace' ? traceExact(message.rays) : traceMesh(message.rays);
      self.postMessage({ type: 'traced', id: message.id, band: message.band, results: out }, [out.buffer]);
      return;
    }
    if (message.type === 'generate') {
      if (!solid) { throw new Error('no exact solid loaded'); }
      self.postMessage({ type: 'generated', id: message.id, doc: generateEvents(solid, message.options || {}) });
      return;
    }
    if (message.type === 'parity') {
      const out = parity(message.rays, message.which);
      self.postMessage({ type: 'parity', id: message.id, band: message.band, counts: out }, [out.buffer]);
      return;
    }
  } catch (e) {
    self.postMessage({ type: 'error', id: message.id, band: message.band, message: e.message });
  }
};
