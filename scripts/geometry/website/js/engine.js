// Where pixels come from.
//
// One interface, two implementations:
//
//   LocalEngine   the JS port of the O2 kernel, in a pool of workers. It is the default, it works
//                 offline, and it is the only engine a published (CSP-locked) copy of this page
//                 can use.
//   RemoteEngine  a bridge on localhost that traces the same rays through the REAL O2 kernel, so
//                 the two can be subtracted pixel by pixel. When it is not there, the page says
//                 so and keeps using LocalEngine -- an offline bridge is never an error.
//
// The contract is deliberately tiny:
//
//   async traceRays(Float32Array rays)  -> Float32Array results
//     rays:    n * 6 floats (ox, oy, oz, dx, dy, dz), cm, directions unit
//     results: n * 5 floats (t, nx, ny, nz, flags); t < 0 means no hit, flags bit 0 reserved
//
// Everything above the interface -- shading, the difference views, the parity overlay -- works on
// those two arrays and does not know which engine produced them.

export class Engine {
  get name() { return 'engine'; }
  get available() { return true; }
  async load() {}
  async traceRays(_rays) { throw new Error('not implemented'); }
  dispose() {}
}

// ------------------------------------------------------------------------------------------
// LocalEngine
// ------------------------------------------------------------------------------------------

export class LocalEngine extends Engine {
  constructor(workerCount = Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 4) - 1))) {
    super();
    this.workers = [];
    this.pending = new Map();
    this.nextId = 1;
    this.ready = false;
    this.stats = { nSurfaces: 0, nTriangles: 0, aabb: null };
    for (let i = 0; i < workerCount; ++i) {
      const worker = new Worker(new URL('./rtworker.js', import.meta.url), { type: 'module' });
      worker.onmessage = (event) => this._onMessage(worker, event.data);
      worker.busy = false;
      this.workers.push(worker);
    }
    this.queue = [];
  }

  get name() { return `local (${this.workers.length} workers)`; }

  _onMessage(worker, message) {
    const entry = this.pending.get(message.id);
    this.pending.delete(message.id);
    worker.busy = false;
    if (entry) {
      if (message.type === 'error') { entry.reject(new Error(message.message)); } else { entry.resolve(message); }
    }
    this._drain();
  }

  _drain() {
    for (const worker of this.workers) {
      if (worker.busy || !this.queue.length) { continue; }
      const job = this.queue.shift();
      worker.busy = true;
      this.pending.set(job.message.id, job);
      worker.postMessage(job.message, job.transfer);
    }
  }

  _post(message, transfer = []) {
    message.id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.queue.push({ message, transfer, resolve, reject });
      this._drain();
    });
  }

  /// Hand every worker the same bytes. ArrayBuffers cannot be transferred to several workers, so
  /// each gets its own copy -- a sidecar is kilobytes and a mesh a few hundred.
  async load({ sidecar, facets, label }) {
    this.ready = false;
    const results = await Promise.all(this.workers.map((worker) => new Promise((resolve, reject) => {
      const id = this.nextId++;
      this.pending.set(id, { resolve, reject });
      worker.busy = true;
      worker.postMessage({
        type: 'load', id, label,
        sidecar: sidecar ? sidecar.slice(0) : null,
        facets: facets ? facets.slice(0) : null,
      });
    })));
    this.stats = { nSurfaces: results[0].nSurfaces, nTriangles: results[0].nTriangles, aabb: results[0].aabb };
    this.ready = true;
    return this.stats;
  }

  async traceRays(rays) {
    const message = { type: 'trace', rays };
    const reply = await this._post(message, [rays.buffer]);
    return reply.results;
  }

  /// Local-only: the same first-hit question asked of the triangle mesh.
  async traceRaysMesh(rays) {
    const reply = await this._post({ type: 'traceMesh', rays }, [rays.buffer]);
    return reply.results;
  }

  /// Local-only: crossings along the whole ray, for the watertightness overlay.
  async parity(rays, which = 'exact') {
    const reply = await this._post({ type: 'parity', rays, which }, [rays.buffer]);
    return reply.counts;
  }

  dispose() { for (const worker of this.workers) { worker.terminate(); } this.workers = []; }
}

// ------------------------------------------------------------------------------------------
// RemoteEngine
// ------------------------------------------------------------------------------------------

/// The bridge to scripts/geometry/tgeoRayService.py, which traces through the real O2 kernel.
///
///   POST /load  {"path": "<abs path to surfaces_*.bin or shape_*.root>"}
///        -> {"ok": true, "kind": "surface"|"shape", "nSurfaces": N, "bbox": [...]}
///   POST /trace raw Float32Array, n*6  ->  raw Float32Array, n*5
export class RemoteEngine extends Engine {
  constructor(port = 8077, host = '127.0.0.1') {
    super();
    this.port = port;
    this.host = host;
    this.ready = false;
    this.info = null;
    this.lastError = null;
  }

  get base() { return `http://${this.host}:${this.port}`; }
  get name() { return `bridge ${this.host}:${this.port}`; }
  get available() { return this.ready; }

  async load({ path }) {
    this.ready = false;
    this.lastError = null;
    try {
      const response = await fetch(`${this.base}/load`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ path }),
      });
      if (!response.ok) { throw new Error(`/load returned HTTP ${response.status}`); }
      const info = await response.json();
      if (!info || info.ok !== true) { throw new Error(info && info.error ? info.error : '/load did not report ok'); }
      this.info = info;
      this.ready = true;
      return info;
    } catch (e) {
      // An unreachable bridge is the normal case, not a failure of the page. A TypeError from
      // fetch means the request never completed: the service is not running, or it is running but
      // the browser refused the cross-origin answer.
      this.lastError = e instanceof TypeError
        ? `no answer from ${this.base} (service not running, or it does not send ` +
          'Access-Control-Allow-Origin and answer the OPTIONS preflight)'
        : e.message;
      this.ready = false;
      throw new Error(this.lastError);
    }
  }

  async traceRays(rays) {
    const response = await fetch(`${this.base}/trace`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/octet-stream' },
      body: rays.buffer,
    });
    if (!response.ok) { throw new Error(`/trace returned HTTP ${response.status}`); }
    const buffer = await response.arrayBuffer();
    const expected = (rays.length / 6) * 5 * 4;
    if (buffer.byteLength !== expected) {
      throw new Error(`/trace returned ${buffer.byteLength} bytes, expected ${expected}`);
    }
    return new Float32Array(buffer);
  }
}
