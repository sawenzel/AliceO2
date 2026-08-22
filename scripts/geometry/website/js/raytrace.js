// The raytracer tab: a progressive CPU raytrace of the EXACT surfaces, the mesh, and the
// differences between them.
//
// Nothing here knows how a pixel was produced -- it asks an Engine for (t, normal) per ray and
// shades the result. That is what makes the engine-difference view possible: the same rays go to
// the JS port and to the bridge into the real O2 kernel, and the two answers are subtracted.

import { LocalEngine, RemoteEngine } from './engine.js';

const BAND_ROWS = 8;

export const VIEWS = [
  { key: 'exact', label: 'exact surfaces (local JS port)' },
  { key: 'exactBridge', label: 'exact surfaces (bridge: the real kernel)' },
  { key: 'mesh', label: 'tessellation' },
  { key: 'diff', label: 'difference: exact vs mesh' },
  { key: 'engine', label: 'difference: local vs bridge' },
  { key: 'parityExact', label: 'watertightness: exact' },
  { key: 'parityMesh', label: 'watertightness: mesh' },
];

/// The two views the engine timing pane compares: the SAME picture -- the exact surfaces -- asked
/// of each engine in turn. A mesh frame or a parity frame is a different question and is not put
/// next to a bridge frame, however single-engine it is.
export const PERF_VIEWS = { exact: 'local', exactBridge: 'remote' };

/// Which engine answers the rays of each view. A view served by exactly one engine is a view
/// whose frame time IS that engine's frame time, which is what makes the two comparable.
export const VIEW_ENGINES = {
  exact: ['local'],
  exactBridge: ['remote'],
  mesh: ['local'],
  diff: ['local'],
  engine: ['local', 'remote'],
  parityExact: ['local'],
  parityMesh: ['local'],
};

function normalize(v) {
  const n = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / n, v[1] / n, v[2] / n];
}
function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }

export class Camera {
  constructor() {
    this.target = [0, 0, 0];
    this.origin = [3, 3, 3];
    this.fovY = 45 * Math.PI / 180;
    this.upHint = [0, 1, 0];
  }
  frameBox(box, direction = [1, 0.75, 1]) {
    this.target = [(box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2];
    const radius = 0.5 * Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
    const distance = radius / Math.sin(this.fovY / 2) * 1.15;
    const d = normalize(direction);
    this.origin = [this.target[0] + d[0] * distance, this.target[1] + d[1] * distance, this.target[2] + d[2] * distance];
  }
  basis() {
    const forward = normalize(sub(this.target, this.origin));
    let right = cross(forward, this.upHint);
    if (Math.hypot(right[0], right[1], right[2]) < 1e-9) { right = [1, 0, 0]; }
    right = normalize(right);
    const up = normalize(cross(right, forward));
    return { forward, right, up };
  }
  orbit(dx, dy) {
    const offset = sub(this.origin, this.target);
    const radius = Math.hypot(offset[0], offset[1], offset[2]);
    let theta = Math.atan2(offset[0], offset[2]);
    let phi = Math.acos(Math.max(-1, Math.min(1, offset[1] / radius)));
    theta -= dx * 0.006;
    phi = Math.max(1e-3, Math.min(Math.PI - 1e-3, phi - dy * 0.006));
    this.origin = [
      this.target[0] + radius * Math.sin(phi) * Math.sin(theta),
      this.target[1] + radius * Math.cos(phi),
      this.target[2] + radius * Math.sin(phi) * Math.cos(theta),
    ];
  }
  dolly(factor) {
    const offset = sub(this.origin, this.target);
    const radius = Math.max(1e-4, Math.hypot(offset[0], offset[1], offset[2]) * factor);
    const d = normalize(offset);
    this.origin = [this.target[0] + d[0] * radius, this.target[1] + d[1] * radius, this.target[2] + d[2] * radius];
  }
}

// ------------------------------------------------------------------------------------------
// shading
// ------------------------------------------------------------------------------------------

const BACKGROUND_TOP = [16, 20, 27];
const BACKGROUND_BOTTOM = [10, 12, 16];
const GRID_COLOUR = [40, 54, 70];

/// The background: a vertical gradient, with a ground grid where the ray meets the reference plane
/// so the eye keeps its bearings while orbiting.
function background(out, o, dirX, dirY, dirZ, ox, oy, oz, groundY, gridStep) {
  const skyMix = Math.max(0, Math.min(1, 0.5 + 0.5 * dirY));
  let r = BACKGROUND_BOTTOM[0] + (BACKGROUND_TOP[0] - BACKGROUND_BOTTOM[0]) * skyMix;
  let g = BACKGROUND_BOTTOM[1] + (BACKGROUND_TOP[1] - BACKGROUND_BOTTOM[1]) * skyMix;
  let b = BACKGROUND_BOTTOM[2] + (BACKGROUND_TOP[2] - BACKGROUND_BOTTOM[2]) * skyMix;
  if (dirY < -1e-6) {
    const t = (groundY - oy) / dirY;
    if (t > 0) {
      const hx = ox + dirX * t, hz = oz + dirZ * t;
      const fx = Math.abs(hx / gridStep - Math.round(hx / gridStep));
      const fz = Math.abs(hz / gridStep - Math.round(hz / gridStep));
      const line = Math.min(fx, fz);
      // Widen the line with distance and fade it out: a fixed-width line in cell units aliases
      // into speckle where the plane is seen almost edge-on.
      const cells = t / Math.max(gridStep, 1e-9);
      const halfWidth = Math.min(0.2, 0.02 * (1 + cells / 25));
      const fade = cells > 120 ? 0 : Math.exp(-cells / 45);
      if (line < halfWidth && fade > 0.02) {
        const w = (1 - line / halfWidth) * fade;
        r += (GRID_COLOUR[0] - r) * w; g += (GRID_COLOUR[1] - g) * w; b += (GRID_COLOUR[2] - b) * w;
      }
    }
  }
  out[0] = r; out[1] = g; out[2] = b;
}

const MATERIAL = { exact: [150, 178, 214], mesh: [214, 150, 110] };

const SPECULAR_EXPONENT = 48;

/// Headlight plus one cool fill, from the analytic normal, with a Blinn-Phong specular on each.
/// This is where the exact surfaces pay off visually: a cylinder shades smoothly, and its highlight
/// runs as one unbroken band, because its normal is the cylinder's and not a facet's.
function shade(out, nx, ny, nz, dirX, dirY, dirZ, right, up, base) {
  // headlight along the view direction, fill from up-left in camera space
  const headlight = Math.max(0, -(nx * dirX + ny * dirY + nz * dirZ));
  const fillDir = normalize([-right[0] * 0.6 + up[0] * 0.7 - dirX * 0.4,
                             -right[1] * 0.6 + up[1] * 0.7 - dirY * 0.4,
                             -right[2] * 0.6 + up[2] * 0.7 - dirZ * 0.4]);
  const fill = Math.max(0, nx * fillDir[0] + ny * fillDir[1] + nz * fillDir[2]);
  const ambient = 0.16;
  const intensity = ambient + 0.74 * headlight + 0.28 * fill;

  // Blinn-Phong: the half vector between the light and the eye. The eye direction is -ray, so for
  // the headlight the half vector is the eye direction itself and its highlight is the headlight
  // term sharpened; the fill's highlight is the one that travels across a curved face.
  const halfFill = normalize([fillDir[0] - dirX, fillDir[1] - dirY, fillDir[2] - dirZ]);
  const specFill = Math.pow(Math.max(0, nx * halfFill[0] + ny * halfFill[1] + nz * halfFill[2]), SPECULAR_EXPONENT);
  const specHead = Math.pow(headlight, SPECULAR_EXPONENT);
  const specular = 150 * specFill + 55 * specHead;

  // a touch of rim light so a silhouette against the background stays readable
  const rim = Math.pow(1 - headlight, 3) * 0.25;
  out[0] = Math.min(255, base[0] * intensity + 90 * rim + specular);
  out[1] = Math.min(255, base[1] * intensity + 100 * rim + specular * 1.02);
  out[2] = Math.min(255, base[2] * intensity + 120 * rim + specular * 1.06);
}

/// Blue -> yellow -> red, for the depth-difference heatmap.
function heat(out, x) {
  const v = Math.max(0, Math.min(1, x));
  if (v < 0.5) {
    const u = v / 0.5;
    out[0] = 40 + 200 * u; out[1] = 90 + 140 * u; out[2] = 220 - 60 * u;
  } else {
    const u = (v - 0.5) / 0.5;
    out[0] = 240; out[1] = 230 - 190 * u; out[2] = 160 - 150 * u;
  }
}

export class Raytracer {
  constructor(canvas) {
    this.canvas = canvas;
    this.context = canvas.getContext('2d');
    this.camera = new Camera();
    this.view = 'exact';
    this.width = 480;
    this.height = 320;
    this.groundY = 0;
    this.gridStep = 1;
    this.local = new LocalEngine();
    this.remote = null;
    this.generation = 0;
    this.onProgress = null;
    this.onDone = null;
    this.counters = {};
    this.rendering = false;
    // The last frame each engine rendered ALONE: same camera, same resolution, one engine, so the
    // two numbers can be put side by side without apportioning a shared frame between them.
    this.perf = { local: null, remote: null };
    // Only the pixels the part's AABB can possibly cover are traced; the rest is background, which
    // costs a gradient and a grid line and no ray at all.
    this.useScissor = true;
    this.scissor = null;
    // One extra bounce: the mirror direction off every hit, answered by the SAME engine, so the
    // reflection is as exact (or as faceted) as the representation being shown.
    this.reflect = false;
    this.reflectivity = 0.85;   // the grazing-angle limit; Schlick weights it per pixel
  }

  /// Which views take a second, reflected ray batch, and which engine answers it.
  get reflectView() {
    return this.reflect && (this.view === 'exact' || this.view === 'exactBridge' || this.view === 'mesh');
  }

  /// The part's world AABB projected to this frame's pixel rectangle, padded by a couple of pixels
  /// so a silhouette is never clipped by the projection's own rounding. Returns the whole frame
  /// when the box straddles the camera plane, where a perspective projection has no finite rect.
  computeScissor() {
    const full = { x0: 0, y0: 0, x1: this.width, y1: this.height, full: true };
    if (!this.useScissor || !this.box) { return full; }
    const { forward, right, up } = this.camera.basis();
    const tanHalf = Math.tan(this.camera.fovY / 2);
    const aspect = this.width / this.height;
    const o = this.camera.origin;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let corner = 0; corner < 8; ++corner) {
      const px = this.box[(corner & 1) ? 3 : 0];
      const py = this.box[(corner & 2) ? 4 : 1];
      const pz = this.box[(corner & 4) ? 5 : 2];
      const vx = px - o[0], vy = py - o[1], vz = pz - o[2];
      const z = vx * forward[0] + vy * forward[1] + vz * forward[2];
      if (z <= 1e-9) { return full; }        // a corner at or behind the eye: no finite rectangle
      const a = (vx * right[0] + vy * right[1] + vz * right[2]) / z;
      const b = (vx * up[0] + vy * up[1] + vz * up[2]) / z;
      const sx = (a / (tanHalf * aspect) + 1) * this.width / 2 - 0.5;
      const sy = (1 - b / tanHalf) * this.height / 2 - 0.5;
      if (sx < minX) { minX = sx; } if (sx > maxX) { maxX = sx; }
      if (sy < minY) { minY = sy; } if (sy > maxY) { maxY = sy; }
    }
    const pad = 2;
    const x0 = Math.max(0, Math.floor(minX) - pad);
    const y0 = Math.max(0, Math.floor(minY) - pad);
    const x1 = Math.min(this.width, Math.ceil(maxX) + 1 + pad);
    const y1 = Math.min(this.height, Math.ceil(maxY) + 1 + pad);
    if (x1 <= x0 || y1 <= y0) { return { x0: 0, y0: 0, x1: 0, y1: 0, full: false }; }
    return { x0, y0, x1, y1, full: x0 === 0 && y0 === 0 && x1 === this.width && y1 === this.height };
  }

  /// Paint a rectangle of pure background -- no ray is cast for any pixel in it.
  _fillBackground(x0, y0, x1, y1) {
    if (x1 <= x0 || y1 <= y0) { return; }
    const data = this.image.data;
    const rgb = [0, 0, 0];
    const { forward, right, up } = this.camera.basis();
    const aspect = this.width / this.height;
    const tanHalf = Math.tan(this.camera.fovY / 2);
    const o = this.camera.origin;
    for (let y = y0; y < y1; ++y) {
      const ndcY = (1 - 2 * (y + 0.5) / this.height) * tanHalf;
      for (let x = x0; x < x1; ++x) {
        const ndcX = (2 * (x + 0.5) / this.width - 1) * tanHalf * aspect;
        const dx = forward[0] + right[0] * ndcX + up[0] * ndcY;
        const dy = forward[1] + right[1] * ndcX + up[1] * ndcY;
        const dz = forward[2] + right[2] * ndcX + up[2] * ndcY;
        const n = Math.hypot(dx, dy, dz) || 1;
        background(rgb, null, dx / n, dy / n, dz / n, o[0], o[1], o[2], this.groundY, this.gridStep);
        const p = (y * this.width + x) * 4;
        data[p] = rgb[0] | 0; data[p + 1] = rgb[1] | 0; data[p + 2] = rgb[2] | 0; data[p + 3] = 255;
      }
    }
  }

  setSize(width, height) {
    this.width = Math.max(16, Math.round(width));
    this.height = Math.max(16, Math.round(height));
    this.canvas.width = this.width;
    this.canvas.height = this.height;
    this.image = this.context.createImageData(this.width, this.height);
    this.image.data.fill(255);
  }

  async load(part) {
    await this.local.load({ sidecar: part.sidecar, facets: part.facets, label: part.label });
    const box = this.local.stats.aabb || [-1, -1, -1, 1, 1, 1];
    this.box = box;
    this.groundY = box[1] - 0.08 * Math.max(1e-6, box[4] - box[1]);
    const span = Math.max(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
    this.gridStep = Math.pow(10, Math.round(Math.log10(span / 6)));
    this.camera.frameBox(box);
    return this.local.stats;
  }

  /// Connect (or reconnect) the bridge. Never throws: an offline bridge is reported, not fatal.
  async connectBridge(port, path) {
    this.remote = new RemoteEngine(port);
    try {
      const info = await this.remote.load({ path });
      return { ok: true, info };
    } catch (e) {
      return { ok: false, error: e.message };
    }
  }

  get bridgeReady() { return !!(this.remote && this.remote.ready); }

  /// A signature of the current view, so two engine timings can say whether they are comparable.
  cameraKey() {
    return this.camera.origin.map(v => v.toFixed(4)).join(',') + '|' +
           this.camera.target.map(v => v.toFixed(4)).join(',') + '|' + this.camera.fovY.toFixed(5);
  }

  cancel() { this.generation += 1; }

  /// Render the current view progressively: a quarter-resolution pass first so the frame appears
  /// at once, then the full-resolution pass band by band.
  async render() {
    this.cancel();
    const generation = this.generation;
    this.rendering = true;
    const started = performance.now();
    const scissor = this.computeScissor();
    this.scissor = scissor;
    const scissorPixels = Math.max(0, (scissor.x1 - scissor.x0) * (scissor.y1 - scissor.y0));
    this.counters = {
      hit: 0, total: this.width * this.height, mismatch: 0, differing: 0, maxDeltaT: 0,
      parityBreaks: 0, raysTraced: 0,
      scissorPixels, scissorRect: [scissor.x0, scissor.y0, scissor.x1, scissor.y1],
      scissorSaving: 1 - scissorPixels / Math.max(1, this.width * this.height),
    };
    // The frame outside the box is background, painted once, with no ray cast for any of it.
    this._fillBackground(0, 0, this.width, scissor.y0);
    this._fillBackground(0, scissor.y1, this.width, this.height);
    this._fillBackground(0, scissor.y0, scissor.x0, scissor.y1);
    this._fillBackground(scissor.x1, scissor.y0, this.width, scissor.y1);
    this.context.putImageData(this.image, 0, 0);
    try {
      await this._pass(generation, 4);
      if (generation !== this.generation) { return; }
      await this._pass(generation, 1);
    } finally {
      if (generation === this.generation) {
        this.rendering = false;
        const ms = Math.round(performance.now() - started);
        this.counters.ms = ms;
        const slot = PERF_VIEWS[this.view];
        if (slot && !this.counters.error && this.counters.raysTraced > 0) {
          this.perf[slot] = {
            ms, rays: this.counters.raysTraced, width: this.width, height: this.height,
            view: this.view, camera: this.cameraKey(),
            engine: slot === 'local' ? this.local.name : (this.remote ? this.remote.name : 'bridge'),
          };
        }
        if (this.onDone) { this.onDone(this.counters); }
      }
    }
  }

  /// True when the current view sends its primary rays to the bridge.
  _usesBridge() {
    return this.bridgeReady && (this.view === 'exactBridge' || this.view === 'engine');
  }

  async _pass(generation, step) {
    const { y0: top, y1: bottom } = this.scissor;
    // The bridge pays an HTTP round trip per band and threads inside one request, so it wants
    // few, fat batches; the local worker pool wants many thin ones so every worker stays busy.
    const bandRows = this._usesBridge() ? BAND_ROWS * 8 : BAND_ROWS;
    const bands = [];
    for (let y = top; y < bottom; y += bandRows * step) {
      bands.push([y, Math.min(bottom, y + bandRows * step)]);
    }
    // Reset the counters on the full-resolution pass so they describe the finished frame.
    if (step === 1) {
      this.counters.hit = 0; this.counters.mismatch = 0; this.counters.differing = 0;
      this.counters.maxDeltaT = 0; this.counters.parityBreaks = 0;
    }
    const inFlight = [];
    for (const [y0, y1] of bands) {
      if (generation !== this.generation) { return; }
      inFlight.push(this._band(generation, y0, y1, step));
      // keep a few bands in flight so every worker has work, but not the whole image at once;
      // the bridge threads inside each request, so more than a couple in flight just oversubscribes
      const maxInFlight = this._usesBridge() ? 2 : this.local.workers.length * 2;
      if (inFlight.length >= maxInFlight) { await inFlight.shift(); }
    }
    await Promise.all(inFlight);
  }

  async _band(generation, y0, y1, step) {
    const left = this.scissor.x0, rightEdge = this.scissor.x1;
    const rowsRendered = [];
    for (let y = y0; y < y1; y += step) { rowsRendered.push(y); }
    if (!rowsRendered.length) { return; }

    // On a coarse pass only every `step`-th row and column is traced; the result is replicated
    // into the rows and columns between, so the frame appears whole immediately.
    const sampleRows = rowsRendered.length;
    const sampleCols = Math.ceil((rightEdge - left) / step);
    if (sampleCols <= 0) { return; }
    const rays = new Float32Array(sampleRows * sampleCols * 6);
    const { forward, right, up } = this.camera.basis();
    const aspect = this.width / this.height;
    const tanHalf = Math.tan(this.camera.fovY / 2);
    let k = 0;
    for (const y of rowsRendered) {
      const ndcY = (1 - 2 * (y + 0.5) / this.height) * tanHalf;
      for (let sx = 0; sx < sampleCols; ++sx) {
        const x = Math.min(rightEdge - 1, left + sx * step);
        const ndcX = (2 * (x + 0.5) / this.width - 1) * tanHalf * aspect;
        const dx = forward[0] + right[0] * ndcX + up[0] * ndcY;
        const dy = forward[1] + right[1] * ndcX + up[1] * ndcY;
        const dz = forward[2] + right[2] * ndcX + up[2] * ndcY;
        const n = Math.hypot(dx, dy, dz) || 1;
        rays[k++] = this.camera.origin[0]; rays[k++] = this.camera.origin[1]; rays[k++] = this.camera.origin[2];
        rays[k++] = dx / n; rays[k++] = dy / n; rays[k++] = dz / n;
      }
    }
    const raysCopy = rays.slice(0);   // the engine transfers its input; keep our own for shading
    this.counters.raysTraced += sampleRows * sampleCols;

    let primary = null, secondary = null, counts = null;
    try {
      if (this.view === 'exact') {
        primary = await this.local.traceRays(rays);
      } else if (this.view === 'exactBridge') {
        if (!this.bridgeReady) { throw new Error('bridge not connected'); }
        primary = await this.remote.traceRays(rays);
      } else if (this.view === 'mesh') {
        primary = await this.local.traceRaysMesh(rays);
      } else if (this.view === 'diff') {
        primary = await this.local.traceRays(rays);
        secondary = await this.local.traceRaysMesh(raysCopy.slice(0));
      } else if (this.view === 'engine') {
        primary = await this.local.traceRays(rays);
        secondary = this.bridgeReady ? await this.remote.traceRays(raysCopy.slice(0)) : null;
      } else if (this.view === 'parityExact' || this.view === 'parityMesh') {
        counts = await this.local.parity(rays, this.view === 'parityMesh' ? 'mesh' : 'exact');
        primary = await this.local.traceRays(raysCopy.slice(0));
      }
    } catch (e) {
      if (generation === this.generation) { this.counters.error = e.message; }
      return;
    }
    if (generation !== this.generation) { return; }

    let bounce = null;
    if (this.reflectView && primary) {
      try { bounce = await this._bounce(raysCopy, primary); } catch (e) {
        if (generation === this.generation) { this.counters.error = e.message; }
      }
      if (generation !== this.generation) { return; }
    }

    this._paint(rowsRendered, sampleCols, step, raysCopy, primary, secondary, counts, right, up, left, rightEdge, bounce);
    this.context.putImageData(this.image, 0, 0);
    if (this.onProgress) { this.onProgress(this.counters); }
  }

  /// The second ray batch: the mirror direction off every pixel that hit something, sent to the
  /// same engine that answered the first. Only the hits are sent, so the batch is as small as the
  /// silhouette. Returns { slot, rays, results } with slot[sample] the index into the batch, or -1.
  async _bounce(rays, primary) {
    const n = rays.length / 6;
    const slot = new Int32Array(n).fill(-1);
    let m = 0;
    for (let i = 0; i < n; ++i) { if (primary[i * 5] >= 0) { slot[i] = m++; } }
    if (!m) { return { slot, rays: new Float32Array(0), results: new Float32Array(0) }; }

    const span = Math.max(1e-9, Math.hypot(this.box[3] - this.box[0], this.box[4] - this.box[1], this.box[5] - this.box[2]));
    const offset = 1e-5 * span;    // step off the surface, or the bounce re-hits the face it left
    const out = new Float32Array(m * 6);
    for (let i = 0; i < n; ++i) {
      const j = slot[i];
      if (j < 0) { continue; }
      const o = i * 6, q = i * 5;
      const t = primary[q];
      const nx = primary[q + 1], ny = primary[q + 2], nz = primary[q + 3];
      const dx = rays[o + 3], dy = rays[o + 4], dz = rays[o + 5];
      const dot = dx * nx + dy * ny + dz * nz;
      const rx = dx - 2 * dot * nx, ry = dy - 2 * dot * ny, rz = dz - 2 * dot * nz;
      const len = Math.hypot(rx, ry, rz) || 1;
      const k = j * 6;
      out[k] = rays[o] + dx * t + nx * offset;
      out[k + 1] = rays[o + 1] + dy * t + ny * offset;
      out[k + 2] = rays[o + 2] + dz * t + nz * offset;
      out[k + 3] = rx / len; out[k + 4] = ry / len; out[k + 5] = rz / len;
    }
    const copy = out.slice(0);
    this.counters.raysTraced += m;
    let results;
    if (this.view === 'mesh') { results = await this.local.traceRaysMesh(out); }
    else if (this.view === 'exactBridge') { results = await this.remote.traceRays(out); }
    else { results = await this.local.traceRays(out); }
    return { slot, rays: copy, results };
  }

  _paint(rowsRendered, sampleCols, step, rays, primary, secondary, counts, right, up, left, rightEdge, bounce) {
    const width = this.width, height = this.height;
    const data = this.image.data;
    const rgb = [0, 0, 0];
    const mirror = [0, 0, 0];
    const scale = Math.max(1e-9, Math.hypot(this.box[3] - this.box[0], this.box[4] - this.box[1], this.box[5] - this.box[2]));
    const full = step === 1;

    for (let r = 0; r < rowsRendered.length; ++r) {
      const y = rowsRendered[r];
      for (let sx = 0; sx < sampleCols; ++sx) {
        const s = r * sampleCols + sx;
        const o = s * 6;
        const dirX = rays[o + 3], dirY = rays[o + 4], dirZ = rays[o + 5];
        const q = s * 5;
        const t = primary ? primary[q] : -1;
        const hit = t >= 0;

        if (this.view === 'diff' || this.view === 'engine') {
          const t2 = secondary ? secondary[q] : -1;
          const hit2 = t2 >= 0;
          if (!hit && !hit2) {
            background(rgb, null, dirX, dirY, dirZ, rays[o], rays[o + 1], rays[o + 2], this.groundY, this.gridStep);
          } else if (hit !== hit2) {
            // one representation sees a surface here and the other does not
            rgb[0] = 235; rgb[1] = 60; rgb[2] = 70;
            if (full) { this.counters.mismatch += 1; }
          } else {
            const delta = Math.abs(t - t2);
            if (full) {
              this.counters.hit += 1;
              this.counters.maxDeltaT = Math.max(this.counters.maxDeltaT, delta);
              if (delta > 1e-6 * scale) { this.counters.differing += 1; }
            }
            // the heatmap saturates at one part in a thousand of the model size
            heat(rgb, delta / (1e-3 * scale));
          }
        } else if (this.view === 'parityExact' || this.view === 'parityMesh') {
          const crossings = counts ? counts[s] : 0;
          if ((crossings & 1) !== 0) {
            rgb[0] = 255; rgb[1] = 0; rgb[2] = 220;      // entered and never came out
            if (full) { this.counters.parityBreaks += 1; }
          } else if (hit) {
            shade(rgb, primary[q + 1], primary[q + 2], primary[q + 3], dirX, dirY, dirZ, right, up,
                  crossings > 0 ? [90, 110, 130] : [70, 80, 95]);
            if (full) { this.counters.hit += 1; }
          } else {
            background(rgb, null, dirX, dirY, dirZ, rays[o], rays[o + 1], rays[o + 2], this.groundY, this.gridStep);
          }
        } else if (hit) {
          const base = this.view === 'mesh' ? MATERIAL.mesh : MATERIAL.exact;
          shade(rgb, primary[q + 1], primary[q + 2], primary[q + 3], dirX, dirY, dirZ, right, up, base);
          if (bounce) {
            const j = bounce.slot[s];
            if (j >= 0) {
              const bo = j * 6, bq = j * 5;
              const bx = bounce.rays[bo + 3], by = bounce.rays[bo + 4], bz = bounce.rays[bo + 5];
              if (bounce.results[bq] >= 0) {
                shade(mirror, bounce.results[bq + 1], bounce.results[bq + 2], bounce.results[bq + 3],
                      bx, by, bz, right, up, base);
              } else {
                background(mirror, null, bx, by, bz, bounce.rays[bo], bounce.rays[bo + 1], bounce.rays[bo + 2],
                           this.groundY, this.gridStep);
              }
              // Schlick: a face seen edge-on is a mirror, a face seen straight on is barely one.
              // Without this the whole part just gets darker, because the sky it reflects is dark.
              const cosine = Math.max(0, -(primary[q + 1] * dirX + primary[q + 2] * dirY + primary[q + 3] * dirZ));
              const k = this.reflectivity * (0.05 + 0.95 * Math.pow(1 - cosine, 5));
              rgb[0] += (mirror[0] - rgb[0]) * k;
              rgb[1] += (mirror[1] - rgb[1]) * k;
              rgb[2] += (mirror[2] - rgb[2]) * k;
            }
          }
          if (full) { this.counters.hit += 1; }
        } else {
          background(rgb, null, dirX, dirY, dirZ, rays[o], rays[o + 1], rays[o + 2], this.groundY, this.gridStep);
        }

        const r8 = rgb[0] | 0, g8 = rgb[1] | 0, b8 = rgb[2] | 0;
        const x0 = Math.min(rightEdge - 1, left + sx * step);
        for (let yy = y; yy < Math.min(this.scissor.y1, y + step); ++yy) {
          for (let xx = x0; xx < Math.min(rightEdge, x0 + step); ++xx) {
            const p = (yy * width + xx) * 4;
            data[p] = r8; data[p + 1] = g8; data[p + 2] = b8; data[p + 3] = 255;
          }
        }
      }
    }
  }

  dispose() { this.local.dispose(); }
}
