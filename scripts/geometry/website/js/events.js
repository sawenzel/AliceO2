// The event-display tab: track polylines through the part, the step point cloud, and the
// radiograph accumulating on a virtual screen behind it.
//
// This is Track 3b layer 1 -- batch replay. It reads events.json (whatever produced it: the
// synthetic sample generator today, an MCStepLogger reduction of a real o2-sim run tomorrow) and
// plays it back over the geometry. The screen is post-processing: it is not in the geometry, so
// it can be moved without re-running anything.

import * as THREE from '../vendor/three.module.min.js';
import { Viewer3D } from './viewer3d.js';
import { loadEvents, loadBinary } from './data.js';
import { parseSidecar, parseFacets } from './sidecar.js';
import { SurfaceSolid } from './solid.js';
import { state, onPartChanged } from './app.js';

// Colour by species, falling back to charge. Identity is never colour alone here either: the
// legend names every species that appears in the loaded file.
const SPECIES_COLOUR = {
  22: { colour: 0xeda100, label: 'photon' },
  11: { colour: 0x3987e5, label: 'electron' },
  '-11': { colour: 0x86b6ef, label: 'positron' },
  211: { colour: 0xd95926, label: 'pi+' },
  '-211': { colour: 0xe87ba4, label: 'pi-' },
  2212: { colour: 0x199e70, label: 'proton' },
  2112: { colour: 0x8a94a6, label: 'neutron' },
  13: { colour: 0x9085e9, label: 'mu-' },
  '-13': { colour: 0xc0b8ff, label: 'mu+' },
  0: { colour: 0x7c8899, label: 'geantino' },
};
function speciesOf(track) {
  const byPdg = SPECIES_COLOUR[String(track.pdg)];
  if (byPdg) { return byPdg; }
  if (track.charge > 0) { return { colour: 0xd95926, label: `pdg ${track.pdg} (+)` }; }
  if (track.charge < 0) { return { colour: 0x3987e5, label: `pdg ${track.pdg} (-)` }; }
  return { colour: 0x8a94a6, label: `pdg ${track.pdg} (neutral)` };
}

function normalize3(v) {
  const n = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / n, v[1] / n, v[2] / n];
}
function cross3(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }

/// The virtual screen: a plane with an in-plane frame and a pixel grid. Crossings are accumulated
/// into it, which is what makes the radiograph.
export class Screen {
  constructor(spec) {
    this.origin = spec.origin.slice();
    this.normal = normalize3(spec.normal);
    const upRaw = spec.up || [0, 0, 1];
    this.right = normalize3(cross3(upRaw, this.normal));
    this.up = normalize3(cross3(this.normal, this.right));
    this.halfWidth = spec.halfWidth;
    this.halfHeight = spec.halfHeight;
    // The file states a pixel grid; the display starts coarser, because a few hundred crossings
    // spread over 220x220 is a dot field rather than an image. The control resets it.
    this.pixels = [60, 60];
    this.filePixels = spec.pixels ? spec.pixels.slice() : [200, 200];
    this.offset = 0; // how far the screen has been pushed along its own normal, cm
    this.reset();
  }
  reset() {
    this.counts = new Float32Array(this.pixels[0] * this.pixels[1]);
    this.maximum = 0;
    this.crossings = 0;
  }
  setResolution(n) {
    this.pixels = [n, n];
    this.reset();
  }
  planeOrigin() {
    return [this.origin[0] + this.normal[0] * this.offset,
            this.origin[1] + this.normal[1] * this.offset,
            this.origin[2] + this.normal[2] * this.offset];
  }
  /// Where a segment crosses the plane, in screen pixels, or null.
  crossingPixel(a, b) {
    const o = this.planeOrigin();
    const da = (a[0] - o[0]) * this.normal[0] + (a[1] - o[1]) * this.normal[1] + (a[2] - o[2]) * this.normal[2];
    const db = (b[0] - o[0]) * this.normal[0] + (b[1] - o[1]) * this.normal[1] + (b[2] - o[2]) * this.normal[2];
    if ((da > 0) === (db > 0) || da === db) { return null; }
    const t = da / (da - db);
    const p = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
    const rx = (p[0] - o[0]) * this.right[0] + (p[1] - o[1]) * this.right[1] + (p[2] - o[2]) * this.right[2];
    const ry = (p[0] - o[0]) * this.up[0] + (p[1] - o[1]) * this.up[1] + (p[2] - o[2]) * this.up[2];
    if (Math.abs(rx) > this.halfWidth || Math.abs(ry) > this.halfHeight) { return null; }
    const px = Math.min(this.pixels[0] - 1, Math.max(0, Math.floor((rx + this.halfWidth) / (2 * this.halfWidth) * this.pixels[0])));
    const py = Math.min(this.pixels[1] - 1, Math.max(0, Math.floor((this.halfHeight - ry) / (2 * this.halfHeight) * this.pixels[1])));
    return [px, py, p];
  }
  add(a, b) {
    const hit = this.crossingPixel(a, b);
    if (!hit) { return null; }
    const index = hit[1] * this.pixels[0] + hit[0];
    this.counts[index] += 1;
    this.maximum = Math.max(this.maximum, this.counts[index]);
    this.crossings += 1;
    return hit[2];
  }
  /// The four corners of the screen rectangle, for drawing it in the 3D view.
  corners() {
    const o = this.planeOrigin();
    const out = [];
    for (const [sx, sy] of [[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]]) {
      out.push([
        o[0] + this.right[0] * sx * this.halfWidth + this.up[0] * sy * this.halfHeight,
        o[1] + this.right[1] * sx * this.halfWidth + this.up[1] * sy * this.halfHeight,
        o[2] + this.right[2] * sx * this.halfWidth + this.up[2] * sy * this.halfHeight,
      ]);
    }
    return out;
  }
}

/// Separable box blur, used only for display.
function boxBlur(source, nx, ny, radius) {
  const temporary = new Float32Array(nx * ny);
  const output = new Float32Array(nx * ny);
  const width = 2 * radius + 1;
  for (let y = 0; y < ny; ++y) {
    for (let x = 0; x < nx; ++x) {
      let sum = 0;
      for (let k = -radius; k <= radius; ++k) {
        const xx = Math.min(nx - 1, Math.max(0, x + k));
        sum += source[y * nx + xx];
      }
      temporary[y * nx + x] = sum / width;
    }
  }
  for (let y = 0; y < ny; ++y) {
    for (let x = 0; x < nx; ++x) {
      let sum = 0;
      for (let k = -radius; k <= radius; ++k) {
        const yy = Math.min(ny - 1, Math.max(0, y + k));
        sum += temporary[yy * nx + x];
      }
      output[y * nx + x] = sum / width;
    }
  }
  return output;
}

export function initEventsTab() {
  const section = document.getElementById('tab-events');
  section.innerHTML = `
    <div class="split">
      <div>
        <div class="viewport" id="ev-viewport"><div class="hud" id="ev-hud"></div></div>
        <p class="muted small" id="ev-note"></p>
      </div>
      <div>
        <div class="pane">
          <h3>Playback</h3>
          <div class="row">
            <button id="ev-play" class="primary">play</button>
            <button id="ev-step">next event</button>
            <button id="ev-reset">reset</button>
          </div>
          <div class="row">
            <label>event <input type="range" id="ev-index" min="0" max="0" value="0" style="width:9em"></label>
            <span class="status" id="ev-index-label">0</span>
          </div>
          <div class="row">
            <label>speed <input type="range" id="ev-speed" min="1" max="40" value="14" style="width:8em"></label>
            <label><input type="checkbox" id="ev-accumulate" checked> accumulate radiograph</label>
          </div>
          <div class="row">
            <label><input type="checkbox" id="ev-points" checked> step points</label>
            <label><input type="checkbox" id="ev-geometry" checked> geometry</label>
          </div>
        </div>
        <div class="pane">
          <h3>Generate here</h3>
          <p class="muted small">The committed replay is small on purpose. A sharper radiograph is a matter of
            statistics, so the same synthetic gun can be run here, in a worker, against the exact solid of the
            part currently selected.</p>
          <div class="row">
            <label>events <input type="number" id="ev-gen-events" value="5" min="1" max="50" style="width:5em"></label>
            <label>tracks each <input type="number" id="ev-gen-tracks" value="3000" min="10" max="40000" style="width:6.5em"></label>
          </div>
          <div class="row">
            <button id="ev-generate">generate</button>
            <span class="status" id="ev-gen-status"></span>
          </div>
        </div>
        <div class="pane">
          <h3>Radiograph</h3>
          <canvas id="ev-screen" class="rtcanvas" style="width:100%;image-rendering:pixelated"></canvas>
          <div class="row" style="margin-top:8px">
            <label>screen offset <input type="range" id="ev-offset" min="-100" max="100" value="0" style="width:8em"></label>
            <span class="status" id="ev-offset-label">0 cm</span>
          </div>
          <div class="row">
            <label>pixels <select id="ev-pixels">
              <option value="40">40</option><option value="60" selected>60</option><option value="100">100</option>
              <option value="160">160</option>
            </select></label>
            <span class="status" id="ev-crossings"></span>
          </div>
          <p class="muted small">The screen is post-processing, not geometry: moving it changes nothing that was
            simulated. Each track segment that crosses it deposits one count.</p>
        </div>
        <div class="pane">
          <h3>Species</h3>
          <div class="legend" id="ev-legend"></div>
          <dl class="facts" id="ev-facts"></dl>
        </div>
      </div>
    </div>`;

  const viewer = new Viewer3D(document.getElementById('ev-viewport'));
  const hud = document.getElementById('ev-hud');
  const note = document.getElementById('ev-note');
  const screenCanvas = document.getElementById('ev-screen');
  const screenContext = screenCanvas.getContext('2d');

  let doc = null, origin = 'none', screen = null;
  let generatorWorker = null, generatorPart = null;
  let eventIndex = 0, front = 0, playing = false, lastFrame = 0;
  let trackObjects = [], pointObject = null, screenObject = null;
  let maxSegments = 0;

  function clearOverlay() {
    viewer.clearOverlay();
    trackObjects = [];
    pointObject = null;
    screenObject = null;
  }

  /// Build the geometry for one event: one Line per track, plus one Points cloud of every step.
  function buildEvent(index) {
    clearOverlay();
    const event = doc.events[index];
    if (!event) { return; }
    maxSegments = 0;
    const allPoints = [];
    for (const track of event.tracks) {
      const spec = speciesOf(track);
      const points = track.points || [];
      maxSegments = Math.max(maxSegments, points.length);
      const positions = new Float32Array(points.length * 3);
      for (let i = 0; i < points.length; ++i) {
        positions[3 * i] = points[i][0]; positions[3 * i + 1] = points[i][1]; positions[3 * i + 2] = points[i][2];
        allPoints.push(points[i]);
      }
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      geometry.setDrawRange(0, 0);
      const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: spec.colour }));
      line.userData.track = track;
      viewer.overlayGroup.add(line);
      trackObjects.push(line);
    }
    const cloud = new Float32Array(allPoints.length * 3);
    for (let i = 0; i < allPoints.length; ++i) {
      cloud[3 * i] = allPoints[i][0]; cloud[3 * i + 1] = allPoints[i][1]; cloud[3 * i + 2] = allPoints[i][2];
    }
    const pointGeometry = new THREE.BufferGeometry();
    pointGeometry.setAttribute('position', new THREE.BufferAttribute(cloud, 3));
    const span = viewer.controls.target.length() || 1;
    pointObject = new THREE.Points(pointGeometry, new THREE.PointsMaterial({
      color: 0xbfd4ea, size: Math.max(0.01, 0.004 * Math.max(1, span)), sizeAttenuation: true, opacity: 0.55, transparent: true,
    }));
    pointObject.visible = document.getElementById('ev-points').checked;
    viewer.overlayGroup.add(pointObject);
    drawScreenRectangle();
    // Selecting an event shows it whole; playback is what reveals it progressively.
    front = playing ? 0 : maxSegments;
    applyFront();
  }

  function drawScreenRectangle() {
    if (screenObject) { viewer.overlayGroup.remove(screenObject); }
    if (!screen) { return; }
    const corners = screen.corners();
    const positions = new Float32Array(corners.length * 3);
    for (let i = 0; i < corners.length; ++i) {
      positions[3 * i] = corners[i][0]; positions[3 * i + 1] = corners[i][1]; positions[3 * i + 2] = corners[i][2];
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    screenObject = new THREE.Line(geometry, new THREE.LineDashedMaterial({ color: 0x4a5b70 }));
    viewer.overlayGroup.add(screenObject);
    viewer.requestRender();
  }

  /// Reveal the tracks up to the wavefront, and deposit whatever crossed the screen on the way.
  function applyFront() {
    for (const line of trackObjects) {
      const count = line.geometry.getAttribute('position').count;
      line.geometry.setDrawRange(0, Math.min(count, front));
    }
    if (pointObject) {
      pointObject.visible = document.getElementById('ev-points').checked;
    }
    viewer.requestRender();
  }

  function depositEvent(index) {
    const event = doc.events[index];
    if (!event || !screen) { return; }
    for (const track of event.tracks) {
      const points = track.points || [];
      for (let i = 0; i + 1 < points.length; ++i) { screen.add(points[i], points[i + 1]); }
    }
    drawRadiograph();
  }

  function drawRadiograph() {
    if (!screen) { return; }
    const [nx, ny] = screen.pixels;
    screenCanvas.width = nx;
    screenCanvas.height = ny;
    const image = screenContext.createImageData(nx, ny);
    // A few hundred crossings over thousands of pixels is a dot field; one box-blur pass turns it
    // into the density it stands for. The counts themselves are untouched.
    const radius = Math.max(1, Math.round(nx / 40));
    const blurred = boxBlur(screen.counts, nx, ny, radius);
    let maximum = 0;
    for (const v of blurred) { maximum = Math.max(maximum, v); }
    maximum = Math.max(1e-9, maximum);
    for (let i = 0; i < nx * ny; ++i) {
      // sqrt so a thin part of the shadow is still visible next to a hot pixel
      const v = Math.sqrt(blurred[i] / maximum);
      const p = i * 4;
      image.data[p] = 20 + 235 * v;
      image.data[p + 1] = 24 + 210 * v;
      image.data[p + 2] = 32 + 150 * v;
      image.data[p + 3] = 255;
    }
    screenContext.putImageData(image, 0, 0);
    document.getElementById('ev-crossings').textContent = `${screen.crossings} crossings, peak ${screen.maximum}`;
  }

  function showEvent(index, { deposit = true } = {}) {
    eventIndex = Math.max(0, Math.min(doc.events.length - 1, index));
    document.getElementById('ev-index').value = String(eventIndex);
    document.getElementById('ev-index-label').textContent =
      `${eventIndex + 1} / ${doc.events.length} (${doc.events[eventIndex].tracks.length} tracks)`;
    buildEvent(eventIndex);
    if (deposit) { depositEvent(eventIndex); }
  }

  function tick(now) {
    if (!playing || !doc) { return; }
    const speed = Number(document.getElementById('ev-speed').value);
    if (now - lastFrame > 1000 / 60) {
      lastFrame = now;
      front += speed;
      if (front > maxSegments) {
        if (eventIndex + 1 < doc.events.length) {
          showEvent(eventIndex + 1);
        } else if (document.getElementById('ev-accumulate').checked) {
          playing = false;
          document.getElementById('ev-play').textContent = 'play';
        } else {
          showEvent(0);
        }
      } else {
        applyFront();
      }
    }
    requestAnimationFrame(tick);
  }

  async function showGeometryFor(partName) {
    // The events file names the part it was made from; show that geometry when the manifest has it.
    let entry = null;
    try {
      const { listParts } = await import('./data.js');
      const listed = await listParts();
      entry = listed.parts.find(p => p.name === partName) || null;
    } catch (e) { entry = null; }
    if (!entry && state.part) { entry = state.part; }
    if (!entry) { return null; }
    const sidecar = await loadBinary(`testdata/${entry.surfaces}`);
    const solid = new SurfaceSolid(parseSidecar(sidecar, entry.name), entry.name);
    if (entry.facets) {
      try {
        const facets = parseFacets(await loadBinary(`testdata/${entry.facets}`), entry.facets);
        viewer.setMesh(facets.positions, { color: 0x6d7d92, opacity: 0.42 });
      } catch (e) { viewer.clearMesh(); }
    }
    const polylines = [];
    for (const surface of solid.surfaces) {
      try { for (const loop of surface.patchOutline()) { polylines.push(loop); } } catch (e) { /* not drawable */ }
    }
    viewer.setExactEdges(polylines);
    viewer.setGrid(solid.aabb);
    return { entry, solid };
  }

  /// Run the synthetic gun in a worker against the part currently selected, and replay the result.
  async function generateHere() {
    const status = document.getElementById('ev-gen-status');
    if (!state.part || !state.sidecarBuffer) { status.textContent = 'no part loaded'; return; }
    const events = Number(document.getElementById('ev-gen-events').value) || 5;
    const tracks = Number(document.getElementById('ev-gen-tracks').value) || 1000;
    status.textContent = `generating ${events} x ${tracks} tracks ...`;
    const started = performance.now();

    if (!generatorWorker) {
      generatorWorker = new Worker(new URL('./rtworker.js', import.meta.url), { type: 'module' });
    }
    const call = (message, transfer = []) => new Promise((resolve, reject) => {
      const onMessage = (event) => {
        generatorWorker.removeEventListener('message', onMessage);
        if (event.data.type === 'error') { reject(new Error(event.data.message)); } else { resolve(event.data); }
      };
      generatorWorker.addEventListener('message', onMessage);
      generatorWorker.postMessage(message, transfer);
    });

    try {
      if (generatorPart !== state.part.name) {
        await call({ type: 'load', id: 1, label: state.part.name, sidecar: state.sidecarBuffer.slice(0), facets: null });
        generatorPart = state.part.name;
      }
      const reply = await call({ type: 'generate', id: 2, options: { events, tracksPerEvent: tracks, seed: 20260822 } });
      const ms = Math.round(performance.now() - started);
      await adopt(reply.doc, `generated here (${ms} ms)`);
      status.textContent = `${events} x ${tracks} tracks in ${ms} ms`;
    } catch (e) {
      status.textContent = `generation failed: ${e.message}`;
    }
  }

  /// Take a freshly loaded or generated document and set the whole tab up for it.
  async function adopt(newDoc, newOrigin) {
    doc = newDoc;
    origin = newOrigin;
    screen = new Screen(doc.screen);
    const geometry = await showGeometryFor(doc.meta && doc.meta.part);
    const box = geometry ? geometry.solid.aabb : [-1, -1, -1, 1, 1, 1];
    const wide = box.slice();
    for (const corner of screen.corners()) {
      for (let i = 0; i < 3; ++i) {
        wide[i] = Math.min(wide[i], corner[i]);
        wide[i + 3] = Math.max(wide[i + 3], corner[i]);
      }
    }
    viewer.frame(wide);

    const slider = document.getElementById('ev-index');
    slider.max = String(doc.events.length - 1);
    const span = Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
    const offsetSlider = document.getElementById('ev-offset');
    offsetSlider.min = String(-Math.round(span));
    offsetSlider.max = String(Math.round(span));
    offsetSlider.step = String(span / 100);
    offsetSlider.value = '0';
    document.getElementById('ev-offset-label').textContent = '0.00 cm';

    describeDocument();
    showEvent(0);
  }

  function describeDocument() {
    const meta = doc.meta || {};
    const facts = document.getElementById('ev-facts');
    facts.innerHTML = '';
    const fact = (k, v) => {
      const dt = document.createElement('dt'); dt.textContent = k;
      const dd = document.createElement('dd'); dd.textContent = v;
      facts.appendChild(dt); facts.appendChild(dd);
    };
    fact('source', origin);
    fact('part', meta.part || '?');
    fact('generator', meta.generator || '?');
    fact('seed', meta.seed === undefined ? 'n/a' : String(meta.seed));
    fact('representation', meta.representation || 'n/a');
    fact('events', String(doc.events.length));
    fact('tracks', String(doc.events.reduce((n, ev) => n + ev.tracks.length, 0)));
    if (meta.absorbed !== undefined) { fact('absorbed', String(meta.absorbed)); }

    note.textContent = meta.synthetic
      ? 'These tracks are SYNTHETIC: a fixed-seed gun with straight or helical steps and crude scattering and ' +
        'absorption. Only the material is real -- a step counts as inside when the exact solid says so. Real ' +
        'Geant4 transport arrives as the same file, from an MCStepLogger reduction.'
      : `Replay of ${origin}.`;

    const legend = document.getElementById('ev-legend');
    legend.innerHTML = '';
    const seen = new Map();
    for (const event of doc.events) {
      for (const track of event.tracks) {
        const spec = speciesOf(track);
        if (!seen.has(spec.label)) { seen.set(spec.label, spec.colour); }
      }
    }
    for (const [label, colour] of seen) {
      const item = document.createElement('span');
      item.className = 'legend-item';
      const swatch = document.createElement('span');
      swatch.className = 'swatch';
      swatch.style.background = `#${colour.toString(16).padStart(6, '0')}`;
      item.appendChild(swatch);
      item.appendChild(document.createTextNode(label));
      legend.appendChild(item);
    }
    hud.textContent = `${meta.part || ''}\n${doc.events.length} events`;
  }

  async function boot() {
    const loaded = await loadEvents();
    if (!loaded.doc) {
      note.textContent = 'No events.json found. Track 3b writes scripts/geometry/website_data/events.json; ' +
        'tools/make_sample_data.mjs writes the synthetic sample_data/events_sample.json, and the ' +
        '"generate here" panel makes one in the browser.';
      hud.textContent = 'no events';
      return;
    }
    await adopt(loaded.doc, loaded.origin);
  }

  document.getElementById('ev-generate').addEventListener('click', () => generateHere());

  document.getElementById('ev-play').addEventListener('click', (e) => {
    playing = !playing;
    e.target.textContent = playing ? 'pause' : 'play';
    if (playing) { lastFrame = performance.now(); requestAnimationFrame(tick); }
  });
  document.getElementById('ev-step').addEventListener('click', () => {
    if (doc) { showEvent((eventIndex + 1) % doc.events.length); }
  });
  document.getElementById('ev-reset').addEventListener('click', () => {
    if (!doc) { return; }
    if (screen) { screen.reset(); }
    showEvent(0);
    drawRadiograph();
  });
  document.getElementById('ev-index').addEventListener('input', (e) => {
    if (!doc) { return; }
    if (!document.getElementById('ev-accumulate').checked && screen) { screen.reset(); }
    showEvent(Number(e.target.value));
  });
  document.getElementById('ev-points').addEventListener('change', applyFront);
  document.getElementById('ev-geometry').addEventListener('change', (e) => {
    viewer.setMeshVisible(e.target.checked);
    viewer.setEdgesVisible(e.target.checked);
  });
  document.getElementById('ev-offset').addEventListener('input', (e) => {
    if (!screen) { return; }
    screen.offset = Number(e.target.value);
    document.getElementById('ev-offset-label').textContent = `${screen.offset.toFixed(2)} cm`;
    screen.reset();
    drawScreenRectangle();
    depositEvent(eventIndex);
  });
  document.getElementById('ev-pixels').addEventListener('change', (e) => {
    if (!screen) { return; }
    screen.setResolution(Number(e.target.value));
    depositEvent(eventIndex);
  });

  window.addEventListener('tabshown', (e) => { if (e.detail.tab === 'events') { viewer.resize(); } });

  boot().catch((e) => { note.textContent = `event display: ${e.message}`; });
}
