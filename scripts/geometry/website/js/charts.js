// Benchmark charts, drawn as plain SVG -- no chart library, so nothing is fetched and the whole
// thing inlines into one file.
//
// The measure is nanoseconds per call, per navigation function, per representation. The
// representations differ by one to two orders of magnitude, so the axis is logarithmic; a linear
// axis would render every fast bar as a flat line against the slow one.

const SERIES = [
  { key: 'surface', label: 'surface (exact)', color: '#3987e5' },
  { key: 'mesh', label: 'mesh (tessellated)', color: '#d95926' },
  { key: 'shape', label: 'shape (CSG)', color: '#199e70' },
];
const SERIES_BY_KEY = Object.fromEntries(SERIES.map(s => [s.key, s]));

const FUNCTIONS = [
  { key: 'contains', label: 'Contains' },
  { key: 'distFromOutside', label: 'DistFromOutside' },
  { key: 'distFromInside', label: 'DistFromInside' },
  { key: 'safety', label: 'Safety' },
];

const svgNS = 'http://www.w3.org/2000/svg';
function el(name, attrs = {}, text = null) {
  const node = document.createElementNS(svgNS, name);
  for (const [k, v] of Object.entries(attrs)) { node.setAttribute(k, v); }
  if (text !== null) { node.textContent = text; }
  return node;
}

/// null / undefined means "not measured" and is rendered as n/a -- never as zero.
function fmt(value, digits = 3) {
  if (value === null || value === undefined || Number.isNaN(value)) { return 'n/a'; }
  if (value === 0) { return '0'; }
  const a = Math.abs(value);
  if (a >= 1e4 || a < 1e-3) { return value.toExponential(2); }
  return Number(value.toPrecision(digits)).toString();
}

function nsPerCall(representation, fnKey) {
  const fns = representation && representation.functions;
  if (!fns || !fns[fnKey]) { return null; }
  const v = fns[fnKey].nsPerCall;
  return (typeof v === 'number' && Number.isFinite(v)) ? v : null;
}

/// Grouped bars: one group per navigation function, one bar per representation present.
export function barChart(doc, { width = 720, height = 320 } = {}) {
  const reps = (doc.representations || []).filter(r => r && r.name);
  const values = [];
  for (const r of reps) { for (const f of FUNCTIONS) { const v = nsPerCall(r, f.key); if (v && v > 0) { values.push(v); } } }
  const figure = document.createElement('figure');
  figure.className = 'chart';
  if (!values.length) {
    const p = document.createElement('p');
    p.className = 'muted';
    p.textContent = 'No per-function timings in this record (nsPerCall not measured).';
    figure.appendChild(p);
    return figure;
  }

  const margin = { top: 28, right: 16, bottom: 52, left: 62 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;
  const dataMin = Math.min(...values), dataMax = Math.max(...values);
  const lo = Math.pow(10, Math.floor(Math.log10(dataMin)));
  const hi = Math.pow(10, Math.ceil(Math.log10(dataMax)));
  const yOf = (v) => plotH - (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)) * plotH;

  const svg = el('svg', { viewBox: `0 0 ${width} ${height}`, width: '100%', role: 'img' });
  svg.setAttribute('aria-label', 'nanoseconds per call by navigation function and representation, log scale');
  const plot = el('g', { transform: `translate(${margin.left},${margin.top})` });
  svg.appendChild(plot);

  // recessive decade grid
  for (let d = Math.log10(lo); d <= Math.log10(hi) + 1e-9; d += 1) {
    const v = Math.pow(10, d);
    const y = yOf(v);
    plot.appendChild(el('line', { x1: 0, y1: y, x2: plotW, y2: y, stroke: 'var(--grid)', 'stroke-width': 1 }));
    plot.appendChild(el('text', { x: -8, y: y + 4, 'text-anchor': 'end', class: 'axis-label' }, v >= 1000 ? `${v / 1000}k` : `${v}`));
  }
  plot.appendChild(el('text', {
    x: -margin.left + 12, y: -12, class: 'axis-title',
  }, 'ns / call (log)'));

  const groupW = plotW / FUNCTIONS.length;
  const present = reps.filter(r => FUNCTIONS.some(f => nsPerCall(r, f.key) !== null));
  const barW = Math.min(34, (groupW - 18) / Math.max(1, present.length));

  FUNCTIONS.forEach((fn, gi) => {
    const gx = gi * groupW;
    const cx = gx + groupW / 2;
    plot.appendChild(el('text', { x: cx, y: plotH + 20, 'text-anchor': 'middle', class: 'axis-label' }, fn.label));

    present.forEach((rep, si) => {
      const spec = SERIES_BY_KEY[rep.name] || { color: '#8a94a6', label: rep.name };
      const v = nsPerCall(rep, fn.key);
      const x = cx - (present.length * barW) / 2 + si * barW;
      if (v === null || !(v > 0)) {
        plot.appendChild(el('text', { x: x + barW / 2, y: plotH - 6, 'text-anchor': 'middle', class: 'na-label' }, 'n/a'));
        return;
      }
      const y = yOf(v);
      // 2px surface gap between adjacent bars, 4px rounded data end
      const bar = el('rect', {
        x: x + 1, y, width: Math.max(1, barW - 2), height: Math.max(1, plotH - y),
        rx: 4, fill: spec.color,
      });
      bar.appendChild(el('title', {}, `${spec.label} - ${fn.label}: ${fmt(v)} ns/call`));
      plot.appendChild(bar);
      plot.appendChild(el('text', { x: x + barW / 2, y: y - 5, 'text-anchor': 'middle', class: 'value-label' }, fmt(v, 3)));
    });

    // the headline of this whole story: how many times slower the mesh is than the exact surface
    const surfaceValue = nsPerCall(reps.find(r => r.name === 'surface'), fn.key);
    const meshValue = nsPerCall(reps.find(r => r.name === 'mesh'), fn.key);
    if (surfaceValue && meshValue && surfaceValue > 0) {
      plot.appendChild(el('text', { x: cx, y: plotH + 38, 'text-anchor': 'middle', class: 'ratio-label' },
        `mesh / exact ${(meshValue / surfaceValue).toFixed(1)}x`));
    }
  });

  figure.appendChild(svg);

  const legend = document.createElement('div');
  legend.className = 'legend';
  for (const rep of present) {
    const spec = SERIES_BY_KEY[rep.name] || { color: '#8a94a6', label: rep.name };
    const item = document.createElement('span');
    item.className = 'legend-item';
    const swatch = document.createElement('span');
    swatch.className = 'swatch';
    swatch.style.background = spec.color;
    item.appendChild(swatch);
    item.appendChild(document.createTextNode(spec.label));
    legend.appendChild(item);
  }
  figure.appendChild(legend);
  return figure;
}

function table(rows, headers) {
  const t = document.createElement('table');
  const thead = document.createElement('thead');
  const hr = document.createElement('tr');
  for (const h of headers) { const th = document.createElement('th'); th.textContent = h; hr.appendChild(th); }
  thead.appendChild(hr);
  t.appendChild(thead);
  const tbody = document.createElement('tbody');
  for (const row of rows) {
    const tr = document.createElement('tr');
    for (const cell of row) {
      const td = document.createElement('td');
      if (cell && cell.node) { td.appendChild(cell.node); } else { td.textContent = cell === null || cell === undefined ? 'n/a' : String(cell); }
      if (cell === null || cell === undefined || cell === 'n/a') { td.className = 'na'; }
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
  t.appendChild(tbody);
  return t;
}

function countOf(map) {
  if (!map || typeof map !== 'object') { return null; }
  let total = 0;
  for (const v of Object.values(map)) { if (typeof v === 'number') { total += v; } }
  return total;
}

/// Accuracy and X-ray counters, as the table view the charts owe the reader.
export function countersTable(doc) {
  const reps = doc.representations || [];
  const rows = reps.map(rep => {
    const acc = rep.accuracy || {};
    const xray = rep.xray || {};
    return [
      (SERIES_BY_KEY[rep.name] || { label: rep.name }).label,
      acc.capacityRelativeDeviation === undefined ? null : fmt(acc.capacityRelativeDeviation),
      countOf(acc.disagreements),
      countOf(acc.unexplained),
      xray.raysTotal ?? null,
      xray.lost ?? null,
      xray.extra ?? null,
      xray.unterminated ?? null,
      xray.parity ?? null,
    ];
  });
  return table(rows, ['representation', 'capacity dev.', 'disagreements', 'unexplained', 'rays', 'lost', 'extra', 'unterminated', 'parity']);
}

export function partCard(entry) {
  const doc = entry.doc;
  const card = document.createElement('section');
  card.className = 'card';
  const meta = doc.meta || {};
  const h = document.createElement('h3');
  h.textContent = meta.part || '(unnamed part)';
  card.appendChild(h);

  const sub = document.createElement('p');
  sub.className = 'muted small';
  const bits = [];
  for (const [k, v] of Object.entries(meta)) { if (k !== 'part') { bits.push(`${k}: ${v}`); } }
  bits.push(`source: ${entry.source}`);
  sub.textContent = bits.join(' · ');
  card.appendChild(sub);

  card.appendChild(barChart(doc));
  const h2 = document.createElement('h4');
  h2.textContent = 'accuracy and X-ray counters';
  card.appendChild(h2);
  card.appendChild(countersTable(doc));
  return card;
}

export function renderBenchmarks(container, loaded) {
  container.innerHTML = '';
  const note = document.createElement('p');
  note.className = 'muted';
  if (loaded.origin === 'none') {
    note.textContent = 'No benchmark JSON found. Track 2 writes scripts/geometry/website_data/*.json; sample_data/ holds a schema-shaped stand-in.';
    container.appendChild(note);
    return;
  }
  note.innerHTML = loaded.origin === 'website_data'
    ? 'Measured data from <code>scripts/geometry/website_data/</code>.'
    : 'No measured data present yet: these are the <strong>synthetic sample records</strong> in <code>sample_data/</code>, shaped exactly like the Track-2 schema.';
  container.appendChild(note);
  for (const entry of loaded.benchmarks) { container.appendChild(partCard(entry)); }
}
