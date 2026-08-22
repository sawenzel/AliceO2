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
// A /load JIT-compiles the navigation path of a shape it has not seen, which is slow once and
// fast thereafter; a /trace of a full frame is well under a second. Both get a ceiling, because a
// fetch with no timeout turns a service that has stopped answering into a page that hangs for
// ever with no message -- and this service can stop answering.
const LOAD_TIMEOUT_MS = 180000;
const TRACE_TIMEOUT_MS = 60000;

function timeoutSignal(ms) {
  return (typeof AbortSignal !== 'undefined' && AbortSignal.timeout) ? AbortSignal.timeout(ms) : undefined;
}

export class RemoteEngine extends Engine {
  constructor(port = 8077, host = '127.0.0.1') {
    super();
    this.port = port;
    this.host = host;
    this.ready = false;
    this.info = null;
    this.lastError = null;
    // The service holds ONE shape, so a /load must not overlap a /trace in EITHER direction: a
    // load waits for the traces already in the air, and a trace issued while a load runs waits
    // for it. A cancelled frame's bands are still on the wire when the next view asks for the
    // other file, and a trace that reaches the service while it is swapping shapes wedges it.
    this.pending = new Set();
    this._loading = null;
  }

  get base() { return `http://${this.host}:${this.port}`; }
  get name() { return `bridge ${this.host}:${this.port}`; }
  get available() { return this.ready; }

  async load({ path }) {
    let release;
    const gate = new Promise((resolve) => { release = resolve; });
    const previous = this._loading;
    this._loading = gate;              // set before the first await: a trace started now waits
    try {
      if (previous) { await previous; }
      await this.waitIdle();
      return await this._load(path);
    } finally {
      if (this._loading === gate) { this._loading = null; }
      release();
    }
  }

  async _load(path) {
    this.ready = false;
    this.lastError = null;
    try {
      const response = await fetch(`${this.base}/load`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ path }),
        signal: timeoutSignal(LOAD_TIMEOUT_MS),
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
        : (e && e.name === 'TimeoutError'
          ? `${this.base} did not answer /load within ${LOAD_TIMEOUT_MS / 1000} s`
          : e.message);
      this.ready = false;
      throw new Error(this.lastError);
    }
  }

  /// Resolves once every trace this engine has sent has answered (or failed): a cancelled
  /// frame's bands are awaited, not abandoned.
  async waitIdle() {
    while (this.pending.size) { await Promise.allSettled([...this.pending]); }
  }

  async traceRays(rays) {
    while (this._loading) { await this._loading; }
    const job = this._traceRays(rays);
    this.pending.add(job);
    try { return await job; } finally { this.pending.delete(job); }
  }

  async _traceRays(rays) {
    let response;
    try {
      response = await fetch(`${this.base}/trace`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/octet-stream' },
        body: rays.buffer,
        signal: timeoutSignal(TRACE_TIMEOUT_MS),
      });
    } catch (e) {
      // A timeout here means the service accepted the rays and never answered. Say so, and let
      // the page fall back to the local engine rather than waiting for ever.
      this.ready = false;
      throw new Error(e && e.name === 'TimeoutError'
        ? `${this.base} accepted the rays and did not answer within ${TRACE_TIMEOUT_MS / 1000} s ` +
          '(the service is wedged; restart it and press "connect bridge")'
        : `/trace: ${e.message}`);
    }
    if (!response.ok) { throw new Error(`/trace returned HTTP ${response.status}`); }
    const buffer = await response.arrayBuffer();
    const expected = (rays.length / 6) * 5 * 4;
    if (buffer.byteLength !== expected) {
      throw new Error(`/trace returned ${buffer.byteLength} bytes, expected ${expected}`);
    }
    return new Float32Array(buffer);
  }
}
