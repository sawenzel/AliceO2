// Page wiring: the part selector, the tabs, and the state every tab shares.

import { listParts, loadBinary, loadBenchmarks } from './data.js';
import { parseSidecar, parseFacets } from './sidecar.js';
import { SurfaceSolid } from './solid.js';
import { Viewer3D } from './viewer3d.js';
import { renderBenchmarks } from './charts.js';
import { PartSelector } from './partselect.js';

export const state = {
  parts: [],           // every manifest entry
  part: null,          // the manifest entry
  parsed: null,        // parseSidecar result
  solid: null,         // SurfaceSolid
  sidecarBuffer: null, // the raw bytes, so a worker can re-parse them itself
  facets: null,        // { nTriangles, positions }
  facetsBuffer: null,  // the raw facets_*.bin bytes, for the worker
  listeners: [],
};

export function onPartChanged(fn) { state.listeners.push(fn); }

const statusEl = document.getElementById('load-status');
function setStatus(text, isError = false) {
  statusEl.textContent = text;
  statusEl.className = isError ? 'status error' : 'status';
}

// --- tabs --------------------------------------------------------------------------------------

const tabInitialisers = {};
const tabInitialised = {};
export function registerTab(name, initialiser) { tabInitialisers[name] = initialiser; }

let activeTab = 'mesh';
function selectTab(name) {
  activeTab = name;
  for (const button of document.querySelectorAll('nav.tabs button')) {
    button.setAttribute('aria-selected', String(button.dataset.tab === name));
  }
  for (const section of document.querySelectorAll('section.tab')) {
    section.classList.toggle('active', section.id === `tab-${name}`);
  }
  const init = tabInitialisers[name];
  if (init && !tabInitialised[name]) { tabInitialised[name] = true; init(); }
  window.dispatchEvent(new CustomEvent('tabshown', { detail: { tab: name } }));
}
document.getElementById('tabs').addEventListener('click', (e) => {
  const button = e.target.closest('button[data-tab]');
  if (button) { selectTab(button.dataset.tab); }
});

// --- the mesh viewer -----------------------------------------------------------------------------

let viewer = null;
function meshViewer() {
  if (!viewer) {
    viewer = new Viewer3D(document.getElementById('mesh-viewport'));
    document.getElementById('opt-wireframe').addEventListener('change', (e) => viewer.setWireframe(e.target.checked));
    document.getElementById('opt-mesh').addEventListener('change', (e) => viewer.setMeshVisible(e.target.checked));
    document.getElementById('opt-edges').addEventListener('change', (e) => viewer.setEdgesVisible(e.target.checked));
    document.getElementById('btn-frame').addEventListener('click', () => viewer.frame(state.solid ? state.solid.aabb : [-1, -1, -1, 1, 1, 1]));
  }
  return viewer;
}
export function sharedViewer() { return meshViewer(); }

function renderMeshTab() {
  const v = meshViewer();
  const { solid, facets } = state;
  if (facets) {
    v.setMesh(facets.positions);
    v.setWireframe(document.getElementById('opt-wireframe').checked);
  } else {
    v.clearMesh();
  }
  const polylines = [];
  for (const surface of solid.surfaces) {
    try { for (const loop of surface.patchOutline()) { polylines.push(loop); } } catch (e) { /* a face whose outline cannot be sampled is simply not drawn */ }
  }
  v.setExactEdges(polylines);
  v.setGrid(solid.aabb);
  v.setMeshVisible(document.getElementById('opt-mesh').checked);
  v.setEdgesVisible(document.getElementById('opt-edges').checked);
  v.frame(solid.aabb);

  const box = solid.aabb;
  document.getElementById('mesh-hud').textContent =
    `${state.part.name}\n${solid.nSurfaces} exact faces` +
    (facets ? ` / ${facets.nTriangles} triangles` : ' / no mesh') +
    `\nbbox ${(box[3] - box[0]).toFixed(2)} x ${(box[4] - box[1]).toFixed(2)} x ${(box[5] - box[2]).toFixed(2)} cm`;
}

// --- loading a part ------------------------------------------------------------------------------

async function loadPart(entry) {
  setStatus(`loading ${entry.name}...`);
  const sidecarBuffer = await loadBinary(`testdata/${entry.surfaces}`);
  const parsed = parseSidecar(sidecarBuffer, entry.surfaces);
  const solid = new SurfaceSolid(parsed, entry.name);
  let facets = null;
  let facetsBuffer = null;
  if (entry.facets) {
    try {
      facetsBuffer = await loadBinary(`testdata/${entry.facets}`);
      facets = parseFacets(facetsBuffer, entry.facets);
    } catch (e) { facets = null; facetsBuffer = null; }
  }
  state.part = entry;
  state.parsed = parsed;
  state.solid = solid;
  state.sidecarBuffer = sidecarBuffer;
  state.facets = facets;
  state.facetsBuffer = facetsBuffer;

  const notes = [];
  if (parsed.warnings.length) { notes.push(parsed.warnings.length + ' warning(s)'); }
  if (solid.failed.length) { notes.push(`${solid.failed.length} record(s) rejected`); }
  if (solid.unsupported.length) { notes.push(`${solid.unsupported.length} unsupported record(s)`); }
  setStatus(`${entry.name}: ${solid.nSurfaces} faces` + (facets ? `, ${facets.nTriangles} triangles` : '') +
            (notes.length ? ` (${notes.join('; ')})` : ''), solid.failed.length > 0);

  renderMeshTab();
  if (tabInitialised.bench) { renderBenchTab(); }
  for (const fn of state.listeners) {
    try { fn(state); } catch (e) { console.error(e); }
  }
}

// --- boot -------------------------------------------------------------------------------------

registerTab('raytracer', async () => {
  const { initRaytracerTab } = await import('./rtui.js');
  initRaytracerTab();
});

registerTab('check', () => {
  const button = document.getElementById('btn-selfcheck');
  const status = document.getElementById('selfcheck-status');
  const log = document.getElementById('selfcheck-log');
  button.addEventListener('click', async () => {
    button.disabled = true;
    log.textContent = '';
    status.textContent = 'running...';
    try {
      const [{ runSelfCheck }, { listParts }] = await Promise.all([import('./selfcheck.js'), import('./data.js')]);
      const listed = await listParts();
      if (!listed.parts.length) { status.textContent = listed.reason; button.disabled = false; return; }
      const load = async (name) => {
        const entry = listed.parts.find(p => p.name === name);
        return {
          sidecar: await loadBinary(`testdata/${entry.surfaces}`),
          facets: entry.facets ? await loadBinary(`testdata/${entry.facets}`) : null,
        };
      };
      const started = performance.now();
      const report = await runSelfCheck(listed.parts.map(p => p.name), load, (line, ok) => {
        const span = document.createElement('span');
        if (ok === true) { span.className = 'pass'; } else if (ok === false) { span.className = 'fail'; }
        span.textContent = line + '\n';
        log.appendChild(span);
      });
      status.textContent = `${report.pass} passed, ${report.fail} failed in ${Math.round(performance.now() - started)} ms`;
      status.className = report.fail ? 'status error' : 'status';
    } catch (e) {
      status.textContent = `self-check failed to run: ${e.message}`;
      status.className = 'status error';
    }
    button.disabled = false;
  });
});

registerTab('events', async () => {
  const { initEventsTab } = await import('./events.js');
  initEventsTab();
});

let benchmarksLoaded = null;
async function renderBenchTab() {
  const container = document.getElementById('bench-body');
  if (!benchmarksLoaded) {
    container.textContent = 'loading...';
    try { benchmarksLoaded = await loadBenchmarks(); } catch (e) { container.textContent = `benchmarks: ${e.message}`; return; }
  }
  renderBenchmarks(container, benchmarksLoaded, state);
}
registerTab('bench', renderBenchTab);

async function boot() {
  const { parts, reason } = await listParts();
  if (!parts.length) {
    setStatus(reason, true);
    document.getElementById('mesh-hud').textContent = 'no test data';
    return;
  }
  state.parts = parts;
  const selector = new PartSelector(document.getElementById('part-select'), {
    onChange: (entry) => loadPart(entry).catch(e => setStatus(`${entry.name}: ${e.message}`, true)),
  });
  selector.setParts(parts);
  const first = parts.find(p => p.name === 'Bucket') || parts[0];
  selector.select(first.name);
  await loadPart(first).catch(e => setStatus(`${first.name}: ${e.message}`, true));
}

boot();
