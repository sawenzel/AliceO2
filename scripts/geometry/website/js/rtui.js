// The raytracer tab's controls, counters and camera interaction.

import { Raytracer, VIEWS } from './raytrace.js';
import { state, onPartChanged } from './app.js';

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
        <div class="viewport" id="rt-viewport" style="display:flex;align-items:center;justify-content:center;">
          <canvas id="rt-canvas" class="rtcanvas"></canvas>
          <div class="hud" id="rt-hud"></div>
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
          <div class="row">
            <label>resolution <select id="rt-res"></select></label>
            <button id="rt-render" class="primary">render</button>
          </div>
          <div class="row"><span class="status" id="rt-resnote"></span></div>
          <div class="row">
            <label><input type="checkbox" id="rt-scissor" checked> bounding-box scissor</label>
          </div>
          <p class="muted small">The scissor projects the part's world AABB to a screen rectangle and casts
            rays only inside it. Everything outside is background &mdash; a gradient and a grid line, no ray.</p>
          <div class="row">
            <button id="rt-frame">re-frame</button>
            <button id="rt-fromview">camera from 3D view</button>
          </div>
          <p class="muted small" id="rt-viewnote"></p>
        </div>
        <div class="pane">
          <h3>Engine</h3>
          <p class="muted small">The picture is produced by an engine that answers rays. The local engine is
            this page's JS port of the O2 kernel. The bridge sends the same rays to
            <code>scripts/geometry/tgeoRayService.py</code>, which traces them through the real kernel, so the
            two can be subtracted.</p>
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
          <p class="muted small">Each row is the last frame that engine rendered <em>on its own</em>: the whole
            frame was its work, so the two are directly comparable when the camera and the resolution agree.</p>
          <table id="rt-perf"></table>
          <div class="row" style="margin-top:8px">
            <button id="rt-timeboth">time both engines</button>
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
  const viewSelect = document.getElementById('rt-view');
  const resSelect = document.getElementById('rt-res');
  const hud = document.getElementById('rt-hud');
  const counters = document.getElementById('rt-counters');
  const bridgeStatus = document.getElementById('rt-bridge-status');
  const pathInput = document.getElementById('rt-path');
  const portInput = document.getElementById('rt-port');
  const viewNote = document.getElementById('rt-viewnote');
  const footnote = document.getElementById('rt-footnote');

  for (const view of VIEWS) {
    const option = document.createElement('option');
    option.value = view.key;
    option.textContent = view.label;
    viewSelect.appendChild(option);
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
    return part ? `${prefix}testdata/${part.surfaces}` : '';
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
    if (tracer.view !== 'diff' && tracer.view !== 'engine') { return; }
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
      engine: 'The same rays sent to the JS port and to the bridge. This is the port\'s validation instrument: any structure in this image is a disagreement with the real kernel.',
      parityExact: 'Crossings are counted along the whole ray. Magenta means an odd count -- the ray entered the solid and never came out, which is a hole in the surface set.',
      parityMesh: 'The same count against the tessellation. A mesh that loses rays lights up here; the exact solid should not.',
    };
    viewNote.textContent = notes[key] || '';
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
    const rows = [['local', tracer.perf.local], ['bridge', tracer.perf.remote]];
    perfTable.innerHTML = '';
    const head = document.createElement('tr');
    for (const label of ['engine', 'ms/frame', 'rays', 'rate']) {
      const th = document.createElement('th'); th.textContent = label; head.appendChild(th);
    }
    perfTable.appendChild(head);
    for (const [key, entry] of rows) {
      const tr = document.createElement('tr');
      const cells = entry
        ? [entry.engine, `${entry.ms}`, `${entry.rays}`, rate(entry)]
        : [key === 'local' ? 'local' : 'bridge', 'not rendered yet', '-', '-'];
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
    const a = tracer.perf.local, b = tracer.perf.remote;
    if (a && b) {
      const comparable = a.width === b.width && a.height === b.height && a.camera === b.camera;
      td.textContent = comparable
        ? `same camera at ${a.width}x${a.height}: bridge / local ${(b.ms / a.ms).toFixed(2)}x`
        : 'the two frames used a different camera or resolution -- press "time both engines" for a matched pair';
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
    if (tracer.view === 'engine' && tracer.bridgeReady) { add('engine', `local vs ${tracer.remote.name}`); }
    else if (tracer.view === 'exactBridge') { add('engine', tracer.remote ? tracer.remote.name : 'bridge (not connected)'); }
    else { add('engine', tracer.local.name); }
    add('pixels', `${c.total}`);
    add('hit', c.hit === undefined ? 'n/a' : `${c.hit} (${(100 * c.hit / Math.max(1, c.total)).toFixed(1)}%)`);
    if (tracer.view === 'diff' || tracer.view === 'engine') {
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
    if (c.ms !== undefined) { add('time', `${c.ms} ms`); }
    if (c.error) { add('error', c.error); }
  }

  tracer.onProgress = showCounters;
  tracer.onDone = (c) => {
    showCounters(c);
    showPerf();
    const a = tracer.perf.local, b = tracer.perf.remote;
    const both = a && b && a.width === b.width && a.height === b.height && a.camera === b.camera
      ? `\nlocal ${a.ms} ms  vs  bridge ${b.ms} ms  (${(b.ms / a.ms).toFixed(2)}x)` : '';
    const scissor = c.scissorSaving > 0.001
      ? `\nscissor ${c.scissorPixels} px traced, ${(100 * c.scissorSaving).toFixed(1)}% saved` : '';
    hud.textContent = `${state.part ? state.part.name : ''} - ${VIEWS.find(v => v.key === tracer.view).label}\n` +
      `${tracer.width}x${tracer.height} in ${c.ms} ms` + scissor + both;
  };

  let renderTimer = null;
  function scheduleRender(delay = 30) {
    clearTimeout(renderTimer);
    renderTimer = setTimeout(() => { tracer.render(); }, delay);
  }

  viewSelect.addEventListener('change', () => {
    tracer.view = viewSelect.value;
    if ((tracer.view === 'engine' || tracer.view === 'exactBridge') && !tracer.bridgeReady) {
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
  document.getElementById('rt-render').addEventListener('click', () => tracer.render());
  document.getElementById('rt-frame').addEventListener('click', () => {
    if (tracer.box) { tracer.camera.frameBox(tracer.box); scheduleRender(); }
  });
  document.getElementById('rt-fromview').addEventListener('click', async () => {
    const { sharedViewer } = await import('./app.js');
    const viewer = sharedViewer();
    const spec = viewer.cameraSpec();
    tracer.camera.origin = spec.origin.slice();
    tracer.camera.target = [viewer.controls.target.x, viewer.controls.target.y, viewer.controls.target.z];
    tracer.camera.fovY = spec.fovY;
    scheduleRender();
  });

  // camera interaction directly on the raytraced canvas
  let dragging = false, last = { x: 0, y: 0 };
  canvas.addEventListener('pointerdown', (e) => {
    dragging = true; last = { x: e.clientX, y: e.clientY };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener('pointermove', (e) => {
    if (!dragging) { return; }
    const scale = canvas.clientWidth / tracer.width;
    tracer.camera.orbit((e.clientX - last.x) / scale, (e.clientY - last.y) / scale);
    last = { x: e.clientX, y: e.clientY };
    scheduleRender(90);
  });
  const endDrag = (e) => {
    if (!dragging) { return; }
    dragging = false;
    try { canvas.releasePointerCapture(e.pointerId); } catch (_) { /* already released */ }
    scheduleRender(0);
  };
  canvas.addEventListener('pointerup', endDrag);
  canvas.addEventListener('pointercancel', endDrag);
  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    tracer.camera.dolly(Math.exp(e.deltaY * 0.0015));
    scheduleRender(90);
  }, { passive: false });

  /// Try the bridge without being asked. The service resolves a relative path against its own
  /// working directory, which is not this page's, so a couple of documented candidates are tried
  /// and whichever answers is remembered.
  async function autoConnect() {
    if (!state.part) { return false; }
    const suffix = `testdata/${state.part.surfaces}`;
    const candidates = [];
    if (settings.prefix) { candidates.push(settings.prefix); }
    candidates.push('', 'scripts/geometry/website/');
    const port = Number(portInput.value) || 8077;
    for (const prefix of candidates) {
      const path = prefix + suffix;
      const result = await tracer.connectBridge(port, path);
      if (result.ok) {
        settings.prefix = prefix;
        settings.port = port;
        saveBridgeSettings(settings);
        pathInput.value = path;
        bridgeStatus.textContent = `bridge ready: ${result.info.kind || '?'}, ${result.info.nSurfaces ?? '?'} surfaces`;
        bridgeStatus.className = 'status';
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
    for (const view of ['exact', 'exactBridge']) {
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
    const result = await tracer.connectBridge(port, path);
    if (result.ok) {
      const info = result.info;
      bridgeStatus.textContent = `bridge ready: ${info.kind || '?'}, ${info.nSurfaces ?? '?'} surfaces`;
      bridgeStatus.className = 'status';
      if (tracer.view === 'engine') { scheduleRender(); }
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
    pathInput.value = defaultPath();
    describeView();
    showPerf();
    tracer.render();
    await autoConnect();
  }

  onPartChanged(() => { if (initialised) { loadCurrentPart(); } });
  window.addEventListener('tabshown', (e) => {
    if (e.detail.tab === 'raytracer') { applySize(); if (!tracer.rendering) { scheduleRender(); } }
  });

  let initialised = false;
  loadCurrentPart().then(() => { initialised = true; });
}
