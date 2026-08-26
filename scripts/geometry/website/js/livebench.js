// Timing measured here and now, on this machine, through the bridge.
//
// The records in website_data/ are one machine's numbers taken once, and a part that has no
// record has no bars at all. This module asks the bridge to measure the part the reader is
// looking at: POST /load a representation, POST /bench it, and hand the answer to the chart
// labelled as what it is -- a live measurement, with the load average it was taken under.
//
// Every representation the bridge has a kernel loader for is measurable, and since the loaders
// were generalised that is all four of them: the CSG composite (shape.root), the flat halfspace
// solid (flatcsg.bin), the exact sidecar (surfaces.bin) and the facet mesh (facets.bin). The
// service dispatches on the file's own magic, so this module only has to name the artefacts.
//
// /bench answers timing and nothing else. Accuracy, capacity and the X-ray counters come from the
// Track-2 records and are not touched here.

import { REPRESENTATIONS } from './representations.js';

const STORAGE_KEY = 'o2surfaces.bridge';        // the same settings the raytracer tab writes
const PROBE_TIMEOUT_MS = 4000;
const LOAD_TIMEOUT_MS = 180000;                 // a shape the service has not seen JIT-compiles once
const BENCH_TIMEOUT_MS = 120000;

/// The sample count /bench is asked for. The service picks its own repeat count from it and
/// reports both back, so the caption can say what the number rests on.
export const SAMPLES = 4000;

/// part name -> the last live record measured for it, for as long as the page is open.
export const liveResults = new Map();

/// The subjects the bridge can load, in the order they are measured -- the baseline first,
/// because that is what a reader compares against, then the cascade's own order.
export const MEASURABLE = REPRESENTATIONS.map(rep => ({ key: rep.key, field: rep.field,
                                                       role: rep.role }));

/// What a measurement covers unless it is asked for more. A `comparison` subject is left out: the
/// cells-tree exists to reproduce one specific claim about the flat solid, and putting it in the
/// default chart asks every reader to work out what it is before they can read the four bars they
/// came for -- the shape the part was made from, and the three the converter can ship.
export function defaultSubjects() { return MEASURABLE.filter(m => m.role !== 'comparison'); }

/// Whose bounding box fixes the sample points for ALL of them. The exact surface solid's is the
/// tightest and is the CAD extent itself; a TGeoCompositeShape's is the union of its padded leaf
/// boxes and is wider, sometimes much wider, so using one of those would spread the points over a
/// volume that no representation actually fills. Falls back to whatever the part has.
const BOX_PREFERENCE = ['surface', 'flatcsg', 'mesh', 'shape', 'cellstree', 'original'];

function settings() {
  try { return JSON.parse(localStorage.getItem(STORAGE_KEY)) || {}; } catch (e) { return {}; }
}

export function defaultPort() { return Number(settings().port) || 8077; }

function base(port) { return `http://127.0.0.1:${port}`; }

function timeoutSignal(ms) {
  return (typeof AbortSignal !== 'undefined' && AbortSignal.timeout) ? AbortSignal.timeout(ms) : undefined;
}

/// One JSON POST to the bridge. A service that is not running fails as a TypeError from fetch;
/// that is the normal offline case and gets a sentence rather than a stack.
async function post(port, endpoint, body, ms) {
  let response;
  try {
    response = await fetch(`${base(port)}${endpoint}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body || {}),
      signal: timeoutSignal(ms),
    });
  } catch (e) {
    throw new Error(e instanceof TypeError
      ? `no answer from ${base(port)} (scripts/geometry/tgeoRayService.py is not running there)`
      : (e && e.name === 'TimeoutError' ? `${base(port)}${endpoint} did not answer within ${ms / 1000} s` : e.message));
  }
  let payload = null;
  try { payload = await response.json(); } catch (e) { payload = null; }
  if (!response.ok || !payload || payload.ok !== true) {
    throw new Error((payload && payload.error) || `${endpoint} returned HTTP ${response.status}`);
  }
  return payload;
}

/// Is the bridge there at all? An unknown endpoint answers 404 with a JSON body and the CORS
/// headers, which is proof of life that disturbs nothing -- a /load or a /bench would change the
/// shape the service holds. The console records that 404: it is this probe, and it is expected.
/// (A bare OPTIONS cannot be used instead: its own preflight asks for Access-Control-Allow-Methods,
/// which the service does not send, so the browser refuses it and a live bridge looks dead.)
export async function probeBridge(port = defaultPort()) {
  try {
    const response = await fetch(`${base(port)}/ping`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: '{}',
      signal: timeoutSignal(PROBE_TIMEOUT_MS),
    });
    await response.text();
    return { ok: true, port };
  } catch (e) {
    return {
      ok: false,
      port,
      error: e instanceof TypeError
        ? `no bridge on ${base(port)} -- start scripts/geometry/tgeoRayService.py --port ${port}`
        : `bridge on ${base(port)}: ${e.message}`,
    };
  }
}

/// What a live measurement of this part would cover. `comparisons` adds the opt-in subjects; what
/// is left out because the part does not carry it is reported separately from what is left out
/// because it was not asked for, so neither can be mistaken for the other.
export function measurablePlan(entry, { comparisons = false } = {}) {
  const carried = MEASURABLE.filter(m => entry && entry[m.field]);
  const can = carried.filter(m => comparisons || m.role !== 'comparison');
  return { can, skipped: [], optional: carried.filter(m => !can.includes(m)) };
}

/// Put the bridge back to the file the raytracer tab expects to be there. The tracer tracks its
/// own loaded path, so leaving it holding a benchmark's shape would make the next frame silently
/// trace the wrong solid. The service caches loads, so this costs a pointer swap.
async function restoreRaytracer() {
  try {
    const { activeRaytracer } = await import('./app.js');
    const tracer = activeRaytracer();
    if (!tracer || !tracer.bridgeReady) { return null; }
    return await tracer.reloadBridgeShape();
  } catch (e) {
    return null;
  }
}

/// Measure every representation of this part the bridge can load, in turn.
///
/// Returns { part, when, port, samples, repeats, loadAverage, reps: [{ key, path, functions }] }
/// and remembers it in `liveResults`.
export async function measurePart(entry, { port = defaultPort(), samples = SAMPLES,
                                           comparisons = false, onStep = () => {} } = {}) {
  const { can } = measurablePlan(entry, { comparisons });
  if (!can.length) {
    throw new Error(`${entry.name} carries none of the four representation artefacts: the bridge has nothing to load for it`);
  }
  // The service resolves a relative path against its own --root, which is not this page's URL
  // space. The raytracer tab remembers the prefix that worked; try it first, then the documented
  // two, and keep whichever answers for the rest of this run.
  const stored = settings();
  const prefixes = [];
  for (const p of [stored.prefix, '', 'scripts/geometry/website/']) {
    if (typeof p === 'string' && !prefixes.includes(p)) { prefixes.push(p); }
  }
  let prefix = null;
  const reps = [];
  const boxKey = BOX_PREFERENCE.find(key => can.some(m => m.key === key)) || can[0].key;
  const boxOrder = [MEASURABLE.find(m => m.key === boxKey), ...can.filter(m => m.key !== boxKey)];
  // The sample points must be the same for every subject or the bars are not a comparison: the
  // artefacts of one part do NOT share a bounding box (a TGeoCompositeShape's is the union of its
  // padded leaf boxes), and /bench draws from the loaded shape's box unless it is told otherwise.
  // One subject is loaded first purely to read its box, and that box is then passed to every
  // /bench. The load is cached on the service, so the subject's own turn costs nothing extra.
  let sharedBox = null;
  let boxFrom = null;
  for (const target of boxOrder) {
    const file = `testdata/${entry[target.field]}`;
    onStep(`loading ${file} on the bridge ...`);
    let info = null, failure = null;
    for (const candidate of (prefix === null ? prefixes : [prefix])) {
      try { info = await post(port, '/load', { path: candidate + file }, LOAD_TIMEOUT_MS); prefix = candidate; break; } catch (e) { failure = e; }
    }
    if (!info) { throw failure || new Error(`the bridge could not load ${file}`); }
    if (target.key === boxKey && Array.isArray(info.bbox) && info.bbox.length === 6) {
      sharedBox = info.bbox;
      boxFrom = boxKey;
    }
    onStep(`benchmarking ${target.key} (${samples} samples) ...`);
    const measured = await post(port, '/bench', sharedBox ? { samples, bbox: sharedBox } : { samples },
                                BENCH_TIMEOUT_MS);
    reps.push({
      key: target.key,
      path: prefix + file,
      kind: info.kind,
      info,
      functions: measured.functions || {},
      samples: measured.samples,
      repeats: measured.repeats,
      insideSamples: measured.insideSamples,
      insideFromEntry: measured.insideFromEntry,
      bboxSource: measured.bboxSource,
      loadAverage: measured.loadAverage,
    });
  }
  const record = {
    part: entry.name,
    when: new Date(),
    port,
    sharedBox,
    boxFrom,
    samples: reps[0].samples,
    repeats: reps[0].repeats,
    insideSamples: reps[0].insideSamples,
    loadAverage: Math.max(...reps.map(r => (typeof r.loadAverage === 'number' ? r.loadAverage : 0))),
    reps,
  };
  record.reps.sort((a, b) => MEASURABLE.findIndex(m => m.key === a.key) -
                             MEASURABLE.findIndex(m => m.key === b.key));
  liveResults.set(entry.name, record);
  onStep('restoring the shape the raytracer view needs ...');
  record.restored = await restoreRaytracer();
  return record;
}
