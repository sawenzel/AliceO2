// The one place the page reads data from.
//
// Everything else in the site asks this module, never fetch() directly, so that a later
// single-file build only has to define window.__INLINE_DATA__ (base64 for binaries, objects for
// JSON) before the modules load and nothing else changes. The Artifact CSP blocks every external
// request, and a locally served copy fetches the same relative URLs.

const INLINE = (typeof window !== 'undefined' && window.__INLINE_DATA__) || null;

function base64ToArrayBuffer(b64) {
  const binary = atob(b64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; ++i) { bytes[i] = binary.charCodeAt(i); }
  return bytes.buffer;
}

/// Fetch a binary asset by site-relative path, e.g. "testdata/box/surfaces.bin".
export async function loadBinary(path) {
  if (INLINE && INLINE.binary && INLINE.binary[path]) { return base64ToArrayBuffer(INLINE.binary[path]); }
  const response = await fetch(path);
  if (!response.ok) { throw new Error(`${path}: HTTP ${response.status}`); }
  return response.arrayBuffer();
}

/// Fetch a JSON asset by site-relative path. Returns null when it is simply not there, so a page
/// served without optional data still works.
export async function loadJSON(path, { optional = false } = {}) {
  if (INLINE && INLINE.json && Object.prototype.hasOwnProperty.call(INLINE.json, path)) { return INLINE.json[path]; }
  try {
    const response = await fetch(path);
    if (!response.ok) {
      if (optional) { return null; }
      throw new Error(`${path}: HTTP ${response.status}`);
    }
    return await response.json();
  } catch (e) {
    if (optional) { return null; }
    throw e;
  }
}

/// The parts available to the viewer, from testdata/manifest.json (written by fetch_testdata.sh).
/// Returns [] with a reason when the test data has not been fetched into this checkout.
export async function listParts() {
  const manifest = await loadJSON('testdata/manifest.json', { optional: true });
  if (!manifest || !Array.isArray(manifest.parts)) {
    return { parts: [], reason: 'testdata/manifest.json is missing -- run ./fetch_testdata.sh <gate-workdir> to populate testdata/.' };
  }
  return { parts: manifest.parts, reason: null, generated: manifest.generated, sources: manifest.sources };
}

/// The Track-2 benchmark JSON. Read from website_data/ first (produced by the benchmark run,
/// one directory up), then from the committed sample_data/ fallback.
export async function loadBenchmarks() {
  const index = await loadJSON('../website_data/index.json', { optional: true });
  if (index && Array.isArray(index.files) && index.files.length) {
    const out = [];
    for (const file of index.files) {
      const doc = await loadJSON(`../website_data/${file}`, { optional: true });
      if (doc) { out.push({ source: `website_data/${file}`, doc }); }
    }
    if (out.length) { return { origin: 'website_data', benchmarks: out }; }
  }
  // No index: try the sample set, which follows the same schema and is committed.
  const sampleIndex = await loadJSON('sample_data/index.json', { optional: true });
  if (sampleIndex && Array.isArray(sampleIndex.files)) {
    const out = [];
    for (const file of sampleIndex.files) {
      const doc = await loadJSON(`sample_data/${file}`, { optional: true });
      if (doc) { out.push({ source: `sample_data/${file}`, doc }); }
    }
    return { origin: 'sample_data', benchmarks: out };
  }
  return { origin: 'none', benchmarks: [] };
}

/// The event-display replay. Real o2-sim output when it exists, the synthetic sample otherwise.
export async function loadEvents() {
  const real = await loadJSON('../website_data/events.json', { optional: true });
  if (real) { return { origin: 'website_data/events.json', doc: real }; }
  const sample = await loadJSON('sample_data/events_sample.json', { optional: true });
  if (sample) { return { origin: 'sample_data/events_sample.json', doc: sample }; }
  return { origin: 'none', doc: null };
}
