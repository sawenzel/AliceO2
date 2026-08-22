// The Assembly tab's controls: load the placement table, put the barrel on screen, and drive it.

import { loadAssemblyData, AssemblyView, statsLine, barrelLine, grouped, AXIS_NAME } from './assembly.js';

let view = null;
let model = null;
let loading = null;

const el = (id) => document.getElementById(id);

/// Where the slice sliders sit, mapped onto the axis extent. The width uses a squared response so
/// the narrow end -- one ring out of twenty-eight -- is reachable with the slider, not by luck.
function sliceWindow() {
  const lo0 = model.box[model.axis], hi0 = model.box[model.axis + 3];
  const extent = hi0 - lo0;
  let centre = lo0 + extent * (Number(el('asm-centre').value) / 1000);
  const t = Number(el('asm-width').value) / 1000;
  const width = model.step + (extent - model.step) * t * t;
  // A window no wider than the ring spacing snaps to the nearest station, so dragging it walks
  // the barrel ring by ring instead of landing between two of them and showing nothing.
  if (width < 1.5 * model.step && model.stations.length) {
    centre = model.stations.reduce((best, s) => (Math.abs(s - centre) < Math.abs(best - centre) ? s : best), model.stations[0]);
  }
  return { lo: centre - width / 2, hi: centre + width / 2, width, extent, centre };
}

/// The extent of what the slice leaves on screen: the assembly box, cut to the window.
function sliceBox() {
  const { lo, hi } = sliceWindow();
  const box = model.box.slice();
  box[model.axis] = Math.max(box[model.axis], lo);
  box[model.axis + 3] = Math.min(box[model.axis + 3], hi);
  if (box[model.axis + 3] <= box[model.axis]) { return model.box; }
  return box;
}

function applySlice() {
  const { lo, hi, width, extent } = sliceWindow();
  view.setSlice(lo, hi);
  const whole = width >= extent - 1e-6;
  el('asm-slice-readout').textContent = whole
    ? `whole barrel (${AXIS_NAME[model.axis]} ${model.box[model.axis].toFixed(0)} to ${model.box[model.axis + 3].toFixed(0)} cm)`
    : `${AXIS_NAME[model.axis]} ${lo.toFixed(1)} to ${hi.toFixed(1)} cm (${width.toFixed(1)} cm wide)`;
  renderHud();
}

function renderHud() {
  const lines = [
    model.index.label || model.index.name,
    statsLine(model),
    barrelLine(model),
    `${grouped(view.visibleSolids)} solids drawn`,
  ];
  if (view.spinning && view.fps) {
    const rate = view.fps >= 10 ? view.fps.toFixed(0) : view.fps.toFixed(1);
    lines.push(`${rate} fps (median frame)`);
  }
  lines.push('WebGL only -- no raytracing, no bridge');
  el('assembly-hud').textContent = lines.join('\n');
}

function renderLegend() {
  const box = el('asm-legend');
  box.innerHTML = '';
  for (const entry of view.legend) {
    const item = document.createElement('div');
    item.className = 'legend-item';
    const swatch = document.createElement('span');
    swatch.className = 'swatch';
    swatch.style.background = entry.colour;
    const text = document.createElement('span');
    text.className = 'mono';
    text.textContent = `${entry.name} x${grouped(entry.instances)}`;
    item.append(swatch, text);
    item.title = `${entry.triangles} triangles per instance`;
    box.appendChild(item);
  }
}

export function initAssemblyTab() {
  if (loading) { return loading; }
  const status = el('asm-status');
  status.textContent = 'loading the placement table...';
  loading = loadAssemblyData(el('tab-assembly').dataset.root || 'testdata/otof_assembly').then((built) => {
    model = built;
    view = new AssemblyView(el('assembly-viewport'), model);
    const totals = model.totals;
    status.textContent = `${grouped(totals.triangles || 0)} triangles in ${grouped(view.meshes.length)} draw calls`;

    el('asm-wireframe').addEventListener('change', (e) => view.setWireframe(e.target.checked));
    el('asm-spin').addEventListener('change', (e) => {
      view.setSpin(e.target.checked, () => renderHud());
      renderHud();
    });
    el('asm-frame').addEventListener('click', () => view.frame(sliceBox()));
    el('asm-centre').addEventListener('input', applySlice);
    el('asm-width').addEventListener('input', applySlice);
    el('asm-whole').addEventListener('click', () => {
      el('asm-width').value = '1000'; el('asm-centre').value = '500'; applySlice();
    });
    el('asm-ring').addEventListener('click', () => { el('asm-width').value = '0'; applySlice(); });

    renderLegend();
    applySlice();
    return view;
  }).catch((e) => {
    loading = null;
    status.textContent = e.message;
    status.className = 'status error';
    throw e;
  });
  return loading;
}

/// Re-frame when the tab becomes visible: the viewport has no size while it is display:none, so a
/// view built off-screen would otherwise come up with a stale aspect ratio.
export function refreshAssemblyTab() {
  if (view) { view.viewer.resize(); view.viewer.requestRender(); }
}
