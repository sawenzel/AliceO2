// Page wiring: the part selector, the tabs, and the state every tab shares.

import { listParts, loadBinary, loadJSON, loadBenchmarks, loadDeclineReasons } from './data.js';
import { parseSidecar, parseFacets } from './sidecar.js';
import { SurfaceSolid } from './solid.js';
import { Viewer3D } from './viewer3d.js';
import { renderBenchmarks, csgStructure, shipsKeys, SHIPS_LABEL } from './charts.js';
import { PartSelector } from './partselect.js';

export const state = {
  parts: [],           // every manifest entry
  part: null,          // the manifest entry
  parsed: null,        // parseSidecar result
  solid: null,         // SurfaceSolid
  sidecarBuffer: null, // the raw bytes, so a worker can re-parse them itself
  facets: null,        // { nTriangles, positions }
  facetsBuffer: null,  // the raw facets_*.bin bytes, for the worker
  csg: null,           // testdata/<part>/csg.json: what the CSG recogniser found, and its verdict
  decline: null,       // the matching website_data/decline_reasons.json entry, when there is one
  aabb: null,          // the part's extent, from the exact solid or, failing that, the mesh
  listeners: [],
};

/// The extent of a triangle soup, for a part that has no exact solid to ask.
export function facetsBox(facets) {
  if (!facets) { return [-1, -1, -1, 1, 1, 1]; }
  const p = facets.positions;
  const box = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  for (let i = 0; i < p.length; i += 3) {
    for (let k = 0; k < 3; ++k) {
      if (p[i + k] < box[k]) { box[k] = p[i + k]; }
      if (p[i + k] > box[k + 3]) { box[k + 3] = p[i + k]; }
    }
  }
  return Number.isFinite(box[0]) ? box : [-1, -1, -1, 1, 1, 1];
}

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
    document.getElementById('btn-frame').addEventListener('click', () => viewer.frame(state.aabb || [-1, -1, -1, 1, 1, 1]));
  }
  return viewer;
}
export function sharedViewer() { return meshViewer(); }

// The raytracer tab owns the bridge connection. The benchmarks tab borrows the bridge for a live
// measurement and has to hand it back, so this is where the two find each other.
let raytracer = null;
export function registerRaytracer(tracer) { raytracer = tracer; }
export function activeRaytracer() { return raytracer; }

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
  if (solid) {
    for (const surface of solid.surfaces) {
      try { for (const loop of surface.patchOutline()) { polylines.push(loop); } } catch (e) { /* a face whose outline cannot be sampled is simply not drawn */ }
    }
  }
  v.setExactEdges(polylines);
  const box = state.aabb;
  v.setGrid(box);
  v.setMeshVisible(document.getElementById('opt-mesh').checked);
  v.setEdgesVisible(document.getElementById('opt-edges').checked);
  v.frame(box);

  // The compact per-part line: what it is made of, how big it is, and -- for a part the CSG
  // recogniser accepted -- the composite it ships as, named in TGeo's own shape classes.
  const structure = csgStructure(state.csg);
  // Every representation this part has at full quality, the cascade's own pick first.
  const ships = shipsKeys(state.part.ships, state.part, state.csg)
    .map(key => SHIPS_LABEL[key] || key.toUpperCase()).join(' + ');
  document.getElementById('mesh-hud').textContent =
    `${state.part.name}\n` + (solid ? `${solid.nSurfaces} exact faces` : 'tessellated only -- no exact sidecar') +
    (facets ? ` / ${facets.nTriangles} triangles` : ' / no mesh') +
    `\nbbox ${(box[3] - box[0]).toFixed(2)} x ${(box[4] - box[1]).toFixed(2)} x ${(box[5] - box[2]).toFixed(2)} cm` +
    `\nships ${ships}` + (structure ? `: ${structure}` : '');
}

// --- loading a part ------------------------------------------------------------------------------

async function loadPart(entry) {
  setStatus(`loading ${entry.name}...`);
  // A part the converter declined for exact extraction has no sidecar at all. It is still shown --
  // that is the coverage story -- with every exact view turned off and said so.
  let sidecarBuffer = null, parsed = null, solid = null;
  if (entry.surfaces) {
    sidecarBuffer = await loadBinary(`testdata/${entry.surfaces}`);
    parsed = parseSidecar(sidecarBuffer, entry.surfaces);
    solid = new SurfaceSolid(parsed, entry.name);
  }
  // Optional, and small: the recogniser's record for this part, and its reason for declining.
  const csg = entry.csg ? await loadJSON(`testdata/${entry.csg}`, { optional: true }) : null;
  const reasons = await loadDeclineReasons();
  const decline = reasons ? (reasons.get(entry.name) || null) : null;
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
  state.csg = csg;
  state.decline = decline;

  state.aabb = solid ? solid.aabb : facetsBox(facets);

  const notes = [];
  if (parsed && parsed.warnings.length) { notes.push(parsed.warnings.length + ' warning(s)'); }
  if (solid && solid.failed.length) { notes.push(`${solid.failed.length} record(s) rejected`); }
  if (solid && solid.unsupported.length) { notes.push(`${solid.unsupported.length} unsupported record(s)`); }
  setStatus(`${entry.name}: ` + (solid ? `${solid.nSurfaces} faces` : 'tessellated only -- no exact sidecar') +
            (facets ? `, ${facets.nTriangles} triangles` : '') +
            (notes.length ? ` (${notes.join('; ')})` : ''), !!(solid && solid.failed.length));

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
      // A tessellated-only part has nothing for these assertions to assert against.
      const exactParts = listed.parts.filter(p => p.surfaces);
      const load = async (name) => {
        const entry = listed.parts.find(p => p.name === name);
        return {
          sidecar: await loadBinary(`testdata/${entry.surfaces}`),
          facets: entry.facets ? await loadBinary(`testdata/${entry.facets}`) : null,
        };
      };
      const started = performance.now();
      const report = await runSelfCheck(exactParts.map(p => p.name), load, (line, ok) => {
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
