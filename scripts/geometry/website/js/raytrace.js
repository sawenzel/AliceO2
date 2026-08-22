// The raytracer tab: a progressive CPU raytrace of the EXACT surfaces, the mesh, and the
// differences between them.
//
// Nothing here knows how a pixel was produced -- it asks an Engine for (t, normal) per ray and
// shades the result. That is what makes the engine-difference view possible: the same rays go to
// the JS port and to the bridge into the real O2 kernel, and the two answers are subtracted.

import { LocalEngine, RemoteEngine } from './engine.js';

const BAND_ROWS = 8;

export const VIEWS = [
  { key: 'exact', label: 'exact surfaces' },
  { key: 'mesh', label: 'tessellation' },
  { key: 'diff', label: 'difference: exact vs mesh' },
  { key: 'engine', label: 'difference: local vs bridge' },
  { key: 'parityExact', label: 'watertightness: exact' },
  { key: 'parityMesh', label: 'watertightness: mesh' },
];

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
  distance() { const o = sub(this.origin, this.target); return Math.hypot(o[0], o[1], o[2]); }
}

/// Rays for one horizontal band, as the flat n*6 buffer the engine contract asks for.
export function bandRays(camera, width, height, y0, y1) {
  const { forward, right, up } = camera.basis();
  const aspect = width / height;
  const tanHalf = Math.tan(camera.fovY / 2);
  const rows = y1 - y0;
  const rays = new Float32Array(rows * width * 6);
  let k = 0;
  for (let y = y0; y < y1; ++y) {
    const ndcY = (1 - 2 * (y + 0.5) / height) * tanHalf;
    for (let x = 0; x < width; ++x) {
      const ndcX = (2 * (x + 0.5) / width - 1) * tanHalf * aspect;
      const dx = forward[0] + right[0] * ndcX + up[0] * ndcY;
      const dy = forward[1] + right[1] * ndcX + up[1] * ndcY;
      const dz = forward[2] + right[2] * ndcX + up[2] * ndcY;
      const n = Math.hypot(dx, dy, dz) || 1;
      rays[k++] = camera.origin[0]; rays[k++] = camera.origin[1]; rays[k++] = camera.origin[2];
      rays[k++] = dx / n; rays[k++] = dy / n; rays[k++] = dz / n;
    }
  }
  return rays;
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

/// Headlight plus one cool fill, from the analytic normal. This is where the exact surfaces pay
/// off visually: a cylinder shades smoothly because its normal is the cylinder's, not a facet's.
function shade(out, nx, ny, nz, dirX, dirY, dirZ, right, up, base) {
  // headlight along the view direction, fill from up-left in camera space
  const headlight = Math.max(0, -(nx * dirX + ny * dirY + nz * dirZ));
  const fillDir = normalize([-right[0] * 0.6 + up[0] * 0.7 - dirX * 0.4,
                             -right[1] * 0.6 + up[1] * 0.7 - dirY * 0.4,
                             -right[2] * 0.6 + up[2] * 0.7 - dirZ * 0.4]);
  const fill = Math.max(0, nx * fillDir[0] + ny * fillDir[1] + nz * fillDir[2]);
  const ambient = 0.16;
  const intensity = ambient + 0.78 * headlight + 0.3 * fill;
  // a touch of rim light so a silhouette against the background stays readable
  const rim = Math.pow(1 - headlight, 3) * 0.25;
  out[0] = Math.min(255, base[0] * intensity + 90 * rim);
  out[1] = Math.min(255, base[1] * intensity + 100 * rim);
  out[2] = Math.min(255, base[2] * intensity + 120 * rim);
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

  cancel() { this.generation += 1; }

  /// Render the current view progressively: a quarter-resolution pass first so the frame appears
  /// at once, then the full-resolution pass band by band.
  async render() {
    this.cancel();
    const generation = this.generation;
    this.rendering = true;
    const started = performance.now();
    this.counters = { hit: 0, total: this.width * this.height, mismatch: 0, differing: 0, maxDeltaT: 0, parityBreaks: 0 };
    try {
      await this._pass(generation, 4);
      if (generation !== this.generation) { return; }
      await this._pass(generation, 1);
    } finally {
      if (generation === this.generation) {
        this.rendering = false;
        this.counters.ms = Math.round(performance.now() - started);
        if (this.onDone) { this.onDone(this.counters); }
      }
    }
  }

  async _pass(generation, step) {
    const width = this.width, height = this.height;
    const bands = [];
    for (let y = 0; y < height; y += BAND_ROWS * step) {
      bands.push([y, Math.min(height, y + BAND_ROWS * step)]);
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
      // keep a few bands in flight so every worker has work, but not the whole image at once
      if (inFlight.length >= this.local.workers.length * 2) { await inFlight.shift(); }
    }
    await Promise.all(inFlight);
  }

  async _band(generation, y0, y1, step) {
    const width = this.width;
    const rowsRendered = [];
    for (let y = y0; y < y1; y += step) { rowsRendered.push(y); }
    if (!rowsRendered.length) { return; }

    // On a coarse pass only every `step`-th row and column is traced; the result is replicated
    // into the rows and columns between, so the frame appears whole immediately.
    const sampleRows = rowsRendered.length;
    const sampleCols = Math.ceil(width / step);
    const rays = new Float32Array(sampleRows * sampleCols * 6);
    const { forward, right, up } = this.camera.basis();
    const aspect = width / this.height;
    const tanHalf = Math.tan(this.camera.fovY / 2);
    let k = 0;
    for (const y of rowsRendered) {
      const ndcY = (1 - 2 * (y + 0.5) / this.height) * tanHalf;
      for (let sx = 0; sx < sampleCols; ++sx) {
        const x = Math.min(width - 1, sx * step);
        const ndcX = (2 * (x + 0.5) / width - 1) * tanHalf * aspect;
        const dx = forward[0] + right[0] * ndcX + up[0] * ndcY;
        const dy = forward[1] + right[1] * ndcX + up[1] * ndcY;
        const dz = forward[2] + right[2] * ndcX + up[2] * ndcY;
        const n = Math.hypot(dx, dy, dz) || 1;
        rays[k++] = this.camera.origin[0]; rays[k++] = this.camera.origin[1]; rays[k++] = this.camera.origin[2];
        rays[k++] = dx / n; rays[k++] = dy / n; rays[k++] = dz / n;
      }
    }
    const raysCopy = rays.slice(0);   // the engine transfers its input; keep our own for shading

    let primary = null, secondary = null, counts = null;
    try {
      if (this.view === 'exact') {
        primary = await this.local.traceRays(rays);
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

    this._paint(rowsRendered, sampleCols, step, raysCopy, primary, secondary, counts, right, up);
    this.context.putImageData(this.image, 0, 0);
    if (this.onProgress) { this.onProgress(this.counters); }
  }

  _paint(rowsRendered, sampleCols, step, rays, primary, secondary, counts, right, up) {
    const width = this.width, height = this.height;
    const data = this.image.data;
    const rgb = [0, 0, 0];
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
          shade(rgb, primary[q + 1], primary[q + 2], primary[q + 3], dirX, dirY, dirZ, right, up,
                this.view === 'mesh' ? MATERIAL.mesh : MATERIAL.exact);
          if (full) { this.counters.hit += 1; }
        } else {
          background(rgb, null, dirX, dirY, dirZ, rays[o], rays[o + 1], rays[o + 2], this.groundY, this.gridStep);
        }

        const r8 = rgb[0] | 0, g8 = rgb[1] | 0, b8 = rgb[2] | 0;
        const x0 = Math.min(width - 1, sx * step);
        for (let yy = y; yy < Math.min(height, y + step); ++yy) {
          for (let xx = x0; xx < Math.min(width, x0 + step); ++xx) {
            const p = (yy * width + xx) * 4;
            data[p] = r8; data[p + 1] = g8; data[p + 2] = b8; data[p + 3] = 255;
          }
        }
      }
    }
  }

  dispose() { this.local.dispose(); }
}
