// The raytracer tab's controls, counters and camera interaction.

import { Raytracer, VIEWS, CSG_VIEWS, DIFF_VIEWS } from './raytrace.js';
import { state, onPartChanged } from './app.js';
import { Viewer3D } from './viewer3d.js';

// How long the camera must stand still after a drag before the frame is traced.
const SETTLE_MS = 300;

// Raster width; the height follows the viewport's aspect. "native" instead traces one ray per
// physical device pixel of the viewport, which is the only setting whose picture is not upscaled.
const PRESETS = [
  { label: '360p', width: 480 },
  { label: '480p', width: 640 },
  { label: '720p', width: 960 },
  { label: '1080p', width: 1440 },
  { label: '1440p', width: 1920 },
  { label: 'native (device pixels)', native: true },
];

const STORAGE_KEY = 'o2surfaces.bridge';

function loadBridgeSettings() {
  try { return JSON.parse(localStorage.getItem(STORAGE_KEY)) || {}; } catch (e) { return {}; }
}
function saveBridgeSettings(settings) {
  try { localStorage.setItem(STORAGE_KEY, JSON.stringify(settings)); } catch (e) { /* private mode: not worth a message */ }
}

export function initRaytracerTab() {
  const section = document.getElementById('tab-raytracer');
  section.innerHTML = `
    <div class="split">
      <div>
        <div class="viewport" id="rt-viewport">
          <div class="rt-proxy" id="rt-proxy"></div>
          <canvas id="rt-canvas" class="rtcanvas rt-overlay"></canvas>
          <div class="hud" id="rt-hud"></div>
          <div class="rt-badge" id="rt-badge" hidden>orbiting the WebGL proxy &mdash; release to trace this view</div>
        </div>
        <div id="rt-scale"></div>
        <p class="muted small" id="rt-footnote"></p>
      </div>
      <div>
        <div class="pane">
          <h3>View</h3>
          <div class="row">
            <select id="rt-view" class="grow"></select>
          </div>
          <p class="small" id="rt-coverage" hidden></p>
          <div class="row">
            <label>resolution <select id="rt-res"></select></label>
            <button id="rt-render" class="primary">render</button>
          </div>
          <div class="row"><span class="status" id="rt-resnote"></span></div>
          <div class="row">
            <label><input type="checkbox" id="rt-scissor" checked> bounding-box scissor</label>
            <label><input type="checkbox" id="rt-reflect"> mirror reflection</label>
          </div>
          <p class="muted small">The scissor projects the part's world AABB to a screen rectangle and casts
            rays only inside it. Everything outside is background &mdash; a gradient and a grid line, no ray.</p>
          <p class="muted small">The mirror sends one more ray off every hit, to the <em>same</em> engine, so a
            reflection in the exact view is exact and a reflection in the tessellation is faceted. It costs one
            extra ray per hit pixel, and it applies to the three shaded views only, not to the difference or
            parity overlays.</p>
          <div class="row">
            <button id="rt-frame">re-frame</button>
            <button id="rt-fromview">camera from the mesh tab</button>
          </div>
          <p class="muted small">Drag to orbit, wheel to zoom, shift-drag to pan. Orbiting moves a WebGL
            proxy of the part &mdash; the mesh with the gold exact boundaries over it &mdash; and the ray
            trace of that fixed view starts 300&nbsp;ms after the camera stands still. Nothing is raytraced
            while you are dragging.</p>
          <p class="muted small" id="rt-viewnote"></p>
        </div>
        <div class="pane">
          <h3>Engine</h3>
          <p class="muted small">The picture is produced by an engine that answers rays. The local engine is
            this page's JS port of the O2 kernel. The bridge sends the same rays to
            <code>scripts/geometry/tgeoRayService.py</code>, which traces them through the real kernel, so the
            two can be subtracted. The bridge holds <strong>one</strong> shape at a time, so a switch between an
            exact view and a CSG view re-loads the part's other file &mdash; the sidecar or the
            <code>shape.root</code> composite &mdash; before the first band goes out.</p>
          <div class="row">
            <label>port <input type="number" id="rt-port" value="8077" min="1" max="65535" style="width:6.5em"></label>
            <button id="rt-connect">connect bridge</button>
          </div>
          <div class="row">
            <label class="grow">path on the bridge host
              <input type="text" id="rt-path" class="grow" style="width:100%">
            </label>
          </div>
          <div class="row"><span class="status" id="rt-bridge-status">bridge not connected</span></div>
        </div>
        <div class="pane">
          <h3>ms per frame, per engine</h3>
          <p class="muted small">The same picture &mdash; the exact surfaces &mdash; asked of each engine in
            turn, and then the same part's <em>CSG composite</em> asked of the same kernel. Each row is the
            last such frame, and the whole of it was that engine's work, so they are directly comparable when
            the camera and the resolution agree. A mesh or parity frame is a different question and does not
            appear here.</p>
          <table id="rt-perf"></table>
          <div class="row" style="margin-top:8px">
            <button id="rt-timeboth">time every engine</button>
          </div>
        </div>
        <div class="pane">
          <h3>Counters</h3>
          <dl class="facts" id="rt-counters"></dl>
        </div>
      </div>
    </div>`;

  const canvas = document.getElementById('rt-canvas');
  const tracer = new Raytracer(canvas);
  // Orbiting a CPU raytrace is unusable, so the drag happens on a WebGL proxy of the same part --
  // the mesh with the gold exact boundaries over it, the mesh tab's own look -- and the traced
  // frame is produced once the camera stands still.
  const proxy = new Viewer3D(document.getElementById('rt-proxy'));
  const badge = document.getElementById('rt-badge');
  const viewSelect = document.getElementById('rt-view');
  const resSelect = document.getElementById('rt-res');
  const hud = document.getElementById('rt-hud');
  const counters = document.getElementById('rt-counters');
  const bridgeStatus = document.getElementById('rt-bridge-status');
  const pathInput = document.getElementById('rt-path');
  const portInput = document.getElementById('rt-port');
  const viewNote = document.getElementById('rt-viewnote');
  const coverage = document.getElementById('rt-coverage');
  const footnote = document.getElementById('rt-footnote');

  // Every view but these two asks the exact solid a question; a tessellated-only part cannot answer.
  const MESH_ONLY_VIEWS = new Set(['mesh', 'parityMesh']);
  for (const view of VIEWS) {
    const option = document.createElement('option');
    option.value = view.key;
    option.textContent = view.label;
    viewSelect.appendChild(option);
  }

  /// Why this part has no CSG composite to trace, in the recogniser's own words where they exist.
  function whyNotCSG() {
    const reason = state.decline && state.decline.whyNotCSG;
    if (reason) { return reason; }
    const acceptance = state.csg && state.csg.acceptance;
    if (acceptance && acceptance.reason) { return acceptance.reason; }
    return 'no shape_*.root was written for this part, so there is no CSG composite to trace';
  }

  /// Turn the exact views off for a part that has no sidecar, and the CSG views off for a part
  /// with no shape_*.root -- in both cases saying why rather than failing.
  function applyCoverage() {
    const exact = !!state.sidecarBuffer;
    const csg = !!(state.part && state.part.shape);
    const reason = csg ? '' : whyNotCSG();
    for (const option of viewSelect.options) {
      const needsCSG = CSG_VIEWS.has(option.value);
      option.disabled = (!exact && !MESH_ONLY_VIEWS.has(option.value)) || (needsCSG && !csg);
      if (needsCSG) { option.title = csg ? `traces testdata/${state.part.shape} on the bridge` : reason; }
    }
    if (viewSelect.selectedOptions[0] && viewSelect.selectedOptions[0].disabled) {
      const fallback = exact ? 'exact' : 'mesh';
      viewSelect.value = fallback;
      tracer.view = fallback;
    }
    // The mirror still works on a mesh-only part: its bounce goes to the triangles.
    for (const id of ['rt-connect', 'rt-timeboth']) {
      const el = document.getElementById(id);
      if (el) { el.disabled = !exact; }
    }
    const notes = [];
    if (!exact) {
      notes.push(`${state.part ? state.part.name : 'This part'} is tessellated only -- no exact ` +
        'sidecar. The converter declined it for exact extraction, so there is nothing here to intersect a ray ' +
        'with but triangles: the exact, bridge, difference and exact-parity views are turned off.');
    }
    if (!csg) {
      notes.push(`No CSG composite for this part: ${reason}. The two CSG views are turned off.`);
    }
    coverage.hidden = !notes.length;
    coverage.textContent = notes.join(' ');
  }
  PRESETS.forEach((preset, index) => {
    const option = document.createElement('option');
    option.value = String(index);
    option.textContent = preset.label;
    resSelect.appendChild(option);
  });
  resSelect.value = '0';

  const settings = loadBridgeSettings();
  if (settings.port) { portInput.value = settings.port; }

  function defaultPath() {
    const part = state.part;
    const prefix = settings.prefix || '';
    return part && part.surfaces ? `${prefix}testdata/${part.surfaces}` : '';
  }

  /// The CSG composite for this part, on the bridge host. The service keeps one shape loaded, so
  /// this file and the sidecar take turns; the tracer swaps them when the view changes.
  function shapePath() {
    const part = state.part;
    const prefix = settings.prefix || '';
    return part && part.shape ? `${prefix}testdata/${part.shape}` : null;
  }

  /// What the bridge is, and which of the part's two files it is holding right now.
  function describeBridge() {
    if (!tracer.bridgeReady) { return; }
    const info = tracer.bridgeInfo || {};
    const holding = tracer.bridgeLoaded ? ` holding ${tracer.bridgeLoaded}` : '';
    bridgeStatus.textContent = `bridge ready: ${info.kind || '?'}` +
      (info.nSurfaces === undefined ? '' : `, ${info.nSurfaces} surfaces`) + holding;
    bridgeStatus.className = 'status';
  }

  const resNote = document.getElementById('rt-resnote');

  function applySize() {
    const preset = PRESETS[Number(resSelect.value)];
    const viewport = document.getElementById('rt-viewport');
    const viewWidth = Math.max(1, viewport.clientWidth), viewHeight = Math.max(1, viewport.clientHeight);
    if (preset.native) {
      const ratio = window.devicePixelRatio || 1;
      tracer.setSize(Math.round(viewWidth * ratio), Math.round(viewHeight * ratio));
    } else {
      const aspect = Math.max(0.4, Math.min(2.5, viewWidth / viewHeight));
      tracer.setSize(preset.width, Math.round(preset.width / aspect));
    }
    canvas.style.width = '100%';
    canvas.style.height = 'auto';
    describeSize();
  }

  function describeSize() {
    const rays = tracer.width * tracer.height;
    const preset = PRESETS[Number(resSelect.value)];
    const ratio = window.devicePixelRatio || 1;
    resNote.textContent = `${tracer.width} x ${tracer.height} = ` +
      (rays >= 1e6 ? `${(rays / 1e6).toFixed(2)} Mrays` : `${Math.round(rays / 1e3)} krays`) + ' per full pass' +
      (preset.native ? ` (devicePixelRatio ${ratio})` : '');
  }

  const scaleBar = document.getElementById('rt-scale');
  function renderScale() {
    scaleBar.innerHTML = '';
    if (!DIFF_VIEWS.has(tracer.view)) { return; }
    const span = tracer.box
      ? Math.hypot(tracer.box[3] - tracer.box[0], tracer.box[4] - tracer.box[1], tracer.box[5] - tracer.box[2]) : 1;
    const wrapper = document.createElement('div');
    wrapper.className = 'legend';
    const gradient = document.createElement('span');
    gradient.style.cssText = 'display:inline-block;width:170px;height:11px;border-radius:3px;' +
      'background:linear-gradient(90deg,#285ADC,#F0E6A0,#F0280A)';
    const lo = document.createElement('span'); lo.textContent = '0';
    const hi = document.createElement('span');
    hi.textContent = `${(1e-3 * span).toExponential(1)} cm depth difference`;
    const red = document.createElement('span');
    red.className = 'legend-item';
    const swatch = document.createElement('span');
    swatch.className = 'swatch'; swatch.style.background = '#eb3c46';
    red.appendChild(swatch);
    red.appendChild(document.createTextNode('one hit, the other not'));
    wrapper.append(lo, gradient, hi, red);
    scaleBar.appendChild(wrapper);
  }

  function describeView() {
    renderScale();
    const key = tracer.view;
    const notes = {
      exactBridge: 'The same picture, every pixel answered by scripts/geometry/tgeoRayService.py driving the REAL O2 kernel through TGeo. Identical shading, so what you see is the kernel\'s own hit points and its own analytic normals -- not a difference image of them.',
      exact: 'Every pixel is a ray intersected with the real trimmed patches: quadratics for the plane, cylinder, cone and sphere, the quartic for the torus, and 2D winding against the trim wires. The shading normal is analytic, so a cylinder is smooth.',
      mesh: 'The same camera, against the triangles of facets_*.bin, with flat facet normals. The silhouette is a polygon.',
      diff: 'Red where one representation is hit and the other is not; otherwise a heatmap of the depth difference, saturating at one part in a thousand of the model size.',
      csgBridge: 'The part\'s CSG composite -- the TGeoShape the recogniser produced and the converter ships -- traced by the real kernel through the bridge, shaded exactly like the two exact views. Same pixels, a different representation: the ms/frame row is where the measured Contains ratio becomes something you can watch.',
      engine: 'The same rays sent to the JS port and to the bridge. This is the port\'s validation instrument: any structure in this image is a disagreement with the real kernel.',
      csgDiff: 'The exact surface solid against the CSG composite, on the same rays. This is the acceptance test made visible: a part whose dV_sym came out 0 should show no red at all and a flat blue heatmap, because the two representations are the same solid.',
      parityExact: 'Crossings are counted along the whole ray. Magenta means an odd count -- the ray entered the solid and never came out, which is a hole in the surface set.',
      parityMesh: 'The same count against the tessellation. A mesh that loses rays lights up here; the exact solid should not.',
    };
    viewNote.textContent = (notes[key] || '') +
      (tracer.reflectView ? ' A second, mirrored ray batch is sent off every hit: the reflection is answered by ' +
        'the same engine, so it carries the same normals the surface itself does.' : '');
    const bits = [];
    if (state.solid && state.solid.bsplineTrimFaces) {
      bits.push(`${state.solid.bsplineTrimFaces} face(s) carry B-spline trim curves: their boundaries are the adaptively flattened polyline at 1e-5, exactly as the kernel navigates them, not the rational curve itself.`);
    }
    if (state.solid && state.solid.unsupported.length) {
      bits.push(`${state.solid.unsupported.length} record(s) of an unsupported kind are not drawn at all and are counted in the panel.`);
    }
    if (state.solid && state.solid.failed.length) {
      bits.push(`${state.solid.failed.length} record(s) were rejected while building and are not drawn.`);
    }
    footnote.textContent = bits.join(' ');
  }

  const perfTable = document.getElementById('rt-perf');

  function rate(entry) {
    if (!entry || !entry.ms) { return 'n/a'; }
    const raysPerSecond = entry.rays / (entry.ms / 1000);
    return raysPerSecond >= 1e6 ? `${(raysPerSecond / 1e6).toFixed(2)} Mray/s` : `${Math.round(raysPerSecond / 1e3)} kray/s`;
  }

  /// The two engines side by side. A row appears only once that engine has rendered a frame of
  /// its own; a frame in which both ran is not apportioned between them and is not shown here.
  function showPerf() {
    const rows = [['local', tracer.perf.local], ['bridge', tracer.perf.remote], ['bridge CSG', tracer.perf.csg]];
    perfTable.innerHTML = '';
    const head = document.createElement('tr');
    for (const label of ['engine', 'ms/frame', 'rays', 'rate']) {
      const th = document.createElement('th'); th.textContent = label; head.appendChild(th);
    }
    perfTable.appendChild(head);
    for (const [key, entry] of rows) {
      const tr = document.createElement('tr');
      const cells = entry
        ? [`${entry.engine}${key === 'bridge CSG' ? ', CSG composite' : ''}`, `${entry.ms}`, `${entry.rays}`, rate(entry)]
        : [key, 'not rendered yet', '-', '-'];
      cells.forEach((text, i) => {
        const td = document.createElement('td');
        td.textContent = text;
        if (!entry) { td.className = 'na'; }
        if (i === 0) { td.style.textAlign = 'left'; }
        tr.appendChild(td);
      });
      perfTable.appendChild(tr);
    }
    const note = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = 4;
    td.style.textAlign = 'left';
    td.className = 'na';
    const a = tracer.perf.local, b = tracer.perf.remote, c = tracer.perf.csg;
    const matches = (x, y) => x && y && x.width === y.width && x.height === y.height && x.camera === y.camera;
    if (a && b) {
      const bits = [];
      if (matches(a, b)) { bits.push(`bridge / local ${(b.ms / a.ms).toFixed(2)}x`); }
      if (matches(b, c)) { bits.push(`bridge exact / bridge CSG ${(b.ms / c.ms).toFixed(2)}x`); }
      td.textContent = bits.length
        ? `same camera at ${a.width}x${a.height}: ${bits.join(', ')}`
        : 'the frames used a different camera or resolution -- press "time every engine" for a matched set';
    } else {
      td.textContent = 'render the "exact surfaces" view and the "exact surfaces (bridge)" view to fill both rows';
    }
    note.appendChild(td);
    perfTable.appendChild(note);
  }

  function showCounters(c) {
    counters.innerHTML = '';
    const add = (key, value) => {
      const dt = document.createElement('dt'); dt.textContent = key;
      const dd = document.createElement('dd'); dd.textContent = value;
      counters.appendChild(dt); counters.appendChild(dd);
    };
    if ((tracer.view === 'engine' || tracer.view === 'csgDiff') && tracer.bridgeReady) {
      add('engine', `local vs ${tracer.remote.name}`);
    } else if (tracer.view === 'exactBridge' || tracer.view === 'csgBridge') {
      add('engine', tracer.remote ? tracer.remote.name : 'bridge (not connected)');
    } else { add('engine', tracer.local.name); }
    if (CSG_VIEWS.has(tracer.view)) { add('representation', 'CSG composite (shape.root) on the bridge'); }
    add('pixels', `${c.total}`);
    add('hit', c.hit === undefined ? 'n/a' : `${c.hit} (${(100 * c.hit / Math.max(1, c.total)).toFixed(1)}%)`);
    if (DIFF_VIEWS.has(tracer.view)) {
      add('hit / no-hit disagree', `${c.mismatch}`);
      add('depth differs', `${c.differing}`);
      add('max |dt|', `${c.maxDeltaT.toExponential(2)} cm`);
    }
    if (tracer.view.startsWith('parity')) { add('parity breaks', `${c.parityBreaks}`); }
    if (c.scissorPixels !== undefined) {
      add('traced pixels', `${c.scissorPixels} (${(100 * (1 - c.scissorSaving)).toFixed(1)}% of the frame)`);
      add('scissor saved', `${(100 * c.scissorSaving).toFixed(1)}% of the pixels`);
      add('scissor rect', c.scissorRect.join(', '));
    }
    if (c.bridgeLoadMs !== undefined) { add('bridge /load', `${c.bridgeLoadMs} ms (the shape swap, not counted in the frame)`); }
    if (c.ms !== undefined) { add('time', `${c.ms} ms`); }
    if (c.error) { add('error', c.error); }
  }

  tracer.onProgress = showCounters;
  tracer.onDone = (c) => {
    showCounters(c);
    showPerf();
    describeBridge();
    const a = tracer.perf.local, b = tracer.perf.remote, csg = tracer.perf.csg;
    const same = (x, y) => x && y && x.width === y.width && x.height === y.height && x.camera === y.camera;
    const both = same(a, b)
      ? `\nlocal ${a.ms} ms  vs  bridge ${b.ms} ms  (${(b.ms / a.ms).toFixed(2)}x)` +
        (same(b, csg) ? `  vs  bridge CSG ${csg.ms} ms  (${(b.ms / csg.ms).toFixed(2)}x)` : '') : '';
    const scissor = c.scissorSaving > 0.001
      ? `\nscissor ${c.scissorPixels} px traced, ${(100 * c.scissorSaving).toFixed(1)}% saved` : '';
    hud.textContent = `${state.part ? state.part.name : ''} - ${VIEWS.find(v => v.key === tracer.view).label}\n` +
      `${tracer.width}x${tracer.height} in ${c.ms} ms` + scissor + both;
  };

  let renderTimer = null;
  function scheduleRender(delay = 30) {
    clearTimeout(renderTimer);
    renderTimer = setTimeout(() => { showTraced(); tracer.render(); }, delay);
  }

  /// Show the raytraced frame; the proxy stays underneath it and is what an orbit reveals.
  function showTraced() { canvas.classList.remove('stale'); badge.hidden = true; }
  function showProxy() { canvas.classList.add('stale'); badge.hidden = false; }

  function syncCameraFromProxy() {
    const spec = proxy.cameraSpec();
    tracer.camera.origin = spec.origin.slice();
    tracer.camera.target = [proxy.controls.target.x, proxy.controls.target.y, proxy.controls.target.z];
    tracer.camera.fovY = spec.fovY;
  }

  let dragging = false;
  let settleTimer = null;
  function settleSoon() {
    clearTimeout(settleTimer);
    settleTimer = setTimeout(() => {
      if (dragging) { settleSoon(); return; }   // still on the mouse: wait for the release
      syncCameraFromProxy();
      showTraced();
      tracer.render();
    }, SETTLE_MS);
  }

  /// Any camera move abandons the frame in flight and hands the view back to the proxy.
  proxy.controls.onChange = () => {
    proxy.requestRender();
    tracer.cancel();
    clearTimeout(renderTimer);
    showProxy();
    settleSoon();
  };
  const proxyDom = proxy.renderer.domElement;
  proxyDom.addEventListener('pointerdown', () => { dragging = true; });
  for (const type of ['pointerup', 'pointercancel']) {
    proxyDom.addEventListener(type, () => { dragging = false; settleSoon(); });
  }

  /// Give the proxy the same picture the mesh tab shows: the tessellation, and the exact trimmed
  /// boundaries in gold over it.
  function loadProxy() {
    if (state.facets) { proxy.setMesh(state.facets.positions); } else { proxy.clearMesh(); }
    const polylines = [];
    if (state.solid) {
      for (const surface of state.solid.surfaces) {
        try { for (const loop of surface.patchOutline()) { polylines.push(loop); } } catch (e) { /* not drawable */ }
      }
    }
    proxy.setExactEdges(polylines);
    const box = tracer.box || [-1, -1, -1, 1, 1, 1];
    proxy.setGrid(box);
    proxy.frame(box);
    syncCameraFromProxy();
  }

  viewSelect.addEventListener('change', () => {
    tracer.view = viewSelect.value;
    const needsBridge = tracer.view === 'engine' || tracer.view === 'exactBridge' || CSG_VIEWS.has(tracer.view);
    if (needsBridge && !tracer.bridgeReady) {
      bridgeStatus.textContent = 'bridge offline: connect it, or this view has nothing to render';
      bridgeStatus.className = 'status error';
    }
    describeView();
    scheduleRender();
  });
  resSelect.addEventListener('change', () => { applySize(); scheduleRender(); });
  document.getElementById('rt-scissor').addEventListener('change', (e) => {
    tracer.useScissor = e.target.checked;
    scheduleRender();
  });
  document.getElementById('rt-reflect').addEventListener('change', (e) => {
    tracer.reflect = e.target.checked;
    describeView();
    scheduleRender();
  });
  document.getElementById('rt-render').addEventListener('click', () => { syncCameraFromProxy(); showTraced(); tracer.render(); });
  document.getElementById('rt-frame').addEventListener('click', () => {
    if (tracer.box) { proxy.frame(tracer.box); syncCameraFromProxy(); scheduleRender(); }
  });
  document.getElementById('rt-fromview').addEventListener('click', async () => {
    const { sharedViewer } = await import('./app.js');
    const viewer = sharedViewer();
    const spec = viewer.cameraSpec();
    proxy.applyCameraSpec(spec, [viewer.controls.target.x, viewer.controls.target.y, viewer.controls.target.z]);
    syncCameraFromProxy();
    scheduleRender();
  });

  /// Try the bridge without being asked. The service resolves a relative path against its own
  /// working directory, which is not this page's, so a couple of documented candidates are tried
  /// and whichever answers is remembered.
  async function autoConnect() {
    if (!state.part || !state.part.surfaces) { return false; }
    const suffix = `testdata/${state.part.surfaces}`;
    const candidates = [];
    if (settings.prefix) { candidates.push(settings.prefix); }
    candidates.push('', 'scripts/geometry/website/');
    const port = Number(portInput.value) || 8077;
    for (const prefix of candidates) {
      const path = prefix + suffix;
      const shape = state.part.shape ? `${prefix}testdata/${state.part.shape}` : null;
      const result = await tracer.connectBridge(port, path, { shapePath: shape });
      if (result.ok) {
        settings.prefix = prefix;
        settings.port = port;
        saveBridgeSettings(settings);
        pathInput.value = path;
        describeBridge();
        return true;
      }
      // A dead service fails the same way for every candidate; do not hammer it.
      if (/no answer from/.test(result.error || '')) { break; }
    }
    bridgeStatus.textContent = 'bridge offline; the local engine answers every view that does not need it';
    bridgeStatus.className = 'status';
    return false;
  }

  document.getElementById('rt-timeboth').addEventListener('click', async () => {
    const button = document.getElementById('rt-timeboth');
    if (!tracer.bridgeReady) {
      bridgeStatus.textContent = 'connect the bridge first: with only one engine there is nothing to compare';
      bridgeStatus.className = 'status error';
      return;
    }
    button.disabled = true;
    const restore = viewSelect.value;
    // The CSG composite joins the matched set when the part has one: same camera, same rays, the
    // same kernel -- only the representation differs.
    const views = tracer.csgReady ? ['exact', 'exactBridge', 'csgBridge'] : ['exact', 'exactBridge'];
    for (const view of views) {
      tracer.view = view;
      viewSelect.value = view;
      describeView();
      await tracer.render();
    }
    viewSelect.value = restore;
    tracer.view = restore;
    describeView();
    button.disabled = false;
    await tracer.render();
  });

  document.getElementById('rt-connect').addEventListener('click', async () => {
    const port = Number(portInput.value) || 8077;
    const path = pathInput.value.trim();
    settings.port = port;
    // remember whatever prefix the user put in front of testdata/, so the next part keeps it
    const suffix = state.part ? `testdata/${state.part.surfaces}` : '';
    settings.prefix = suffix && path.endsWith(suffix) ? path.slice(0, path.length - suffix.length) : '';
    saveBridgeSettings(settings);
    bridgeStatus.textContent = `connecting to 127.0.0.1:${port} ...`;
    bridgeStatus.className = 'status';
    const result = await tracer.connectBridge(port, path, { shapePath: shapePath() });
    if (result.ok) {
      describeBridge();
      if (tracer.view === 'engine' || CSG_VIEWS.has(tracer.view)) { scheduleRender(); }
    } else {
      bridgeStatus.textContent = `bridge offline (${result.error}); staying on the local engine`;
      bridgeStatus.className = 'status error';
    }
  });

  async function loadCurrentPart() {
    if (!state.part) { return; }
    hud.textContent = 'loading...';
    applySize();
    await tracer.load({
      sidecar: state.sidecarBuffer,
      facets: state.facetsBuffer ? state.facetsBuffer.slice(0) : null,
      label: state.part.name,
    });
    loadProxy();
    applyCoverage();
    pathInput.value = defaultPath();
    describeView();
    showPerf();
    showTraced();
    tracer.render();
    if (state.sidecarBuffer) { await autoConnect(); }
  }

  onPartChanged(() => { if (initialised) { loadCurrentPart(); } });
  window.addEventListener('tabshown', (e) => {
    if (e.detail.tab === 'raytracer') { proxy.resize(); applySize(); if (!tracer.rendering) { scheduleRender(); } }
  });

  let initialised = false;
  loadCurrentPart().then(() => { initialised = true; });
}
