// Benchmark charts, drawn as plain SVG -- no chart library, so nothing is fetched and the whole
// thing inlines into one file.
//
// The measure is nanoseconds per call, per navigation function, per representation. The
// representations differ by one to two orders of magnitude, so the axis is logarithmic; a linear
// axis would render every fast bar as a flat line against the slow one.
//
// A bar can also be measured on the spot: livebench.js drives the bridge's /bench and hands the
// answer back here, where it is drawn next to the recorded one and labelled as live.

import { liveResults, measurablePlan, measurePart, probeBridge, defaultPort } from './livebench.js';
import { REPRESENTATIONS, PRODUCTS, repLabel, repNote, partHas, partNote } from './representations.js';

// Colour is the only thing the chart decides for itself. Every NAME comes from
// representations.js, so a bar, a badge, a raytracer view and the part card all call the same
// thing by the same word -- which is the whole reason that module exists. The baseline and the
// comparison are grey on purpose: neither is something the converter ships, and giving them a
// hue would put them on the same footing as the four that are.
const SERIES_COLOR = {
  surface: '#3987e5', mesh: '#d95926', shape: '#199e70', flatcsg: '#a163d8',
  original: '#8d94a4', cellstree: '#5c6473',
};

/// The drawing spec for one series of the chart. `name` is either a representation key, for a
/// recorded bar, or `<key>-live` for one measured on the spot.
function seriesSpec(name) {
  const live = typeof name === 'string' && name.endsWith('-live');
  const key = live ? name.slice(0, -'-live'.length) : name;
  const colour = SERIES_COLOR[key];
  if (!colour) { return { color: '#8a94a6', label: String(name) }; }
  if (!live) { return { color: colour, label: `${repLabel(key)}, recorded` }; }
  return { color: LIVE_TINT[key] || colour, stroke: colour, live: true,
           label: `${repLabel(key)}, measured live` };
}

const FUNCTIONS = [
  { key: 'contains', label: 'Contains' },
  { key: 'distFromOutside', label: 'DistFromOutside' },
  { key: 'distFromInside', label: 'DistFromInside' },
  { key: 'safety', label: 'Safety' },
];

// A live series sits beside its recorded twin in the same hue, lighter and dash-outlined, so a
// bar measured a second ago is never mistaken for the record it is being compared against.
const LIVE_TINT = {
  surface: '#7fb2ea', mesh: '#eda184', shape: '#5fc4a2', flatcsg: '#c9a4ec',
  original: '#b9bfca', cellstree: '#8a92a1',
};

const svgNS = 'http://www.w3.org/2000/svg';
function el(name, attrs = {}, text = null) {
  const node = document.createElementNS(svgNS, name);
  for (const [k, v] of Object.entries(attrs)) { node.setAttribute(k, v); }
  if (text !== null) { node.textContent = text; }
  return node;
}

/// The one ratio printed under a group of bars. A series is looked up in its recorded form first
/// and then in its live form, so a part with no Track-2 record still gets the annotation.
function ratioFor(reps, fnKey) {
  const value = (key) => nsPerCall(reps.find(r => r.name === key), fnKey)
                      ?? nsPerCall(reps.find(r => r.name === `${key}-live`), fnKey);
  // Named both ways round, because "FlatCSG / original 4.5x" does not say which of the two is the
  // fast one and a reader should not have to look at the bars to find out. The names are the short
  // forms: one group of this chart is about 160 px wide and the full labels collide.
  const SHORT = { original: 'original', shape: 'CSG', cellstree: 'cells-tree',
                  flatcsg: 'FlatCSG', surface: 'exact', mesh: 'mesh' };
  const faster = (aKey, aVal, bKey, bVal) => (aVal <= bVal
    ? `${SHORT[aKey] || aKey} ${(bVal / aVal).toFixed(1)}x faster`
    : `${SHORT[bKey] || bKey} ${(aVal / bVal).toFixed(1)}x faster`);
  const original = value('original');
  if (original && original > 0) {
    for (const key of ['flatcsg', 'shape', 'surface']) {
      const v = value(key);
      if (v && v > 0) { return faster('original', original, key, v); }
    }
  }
  const surface = value('surface'), mesh = value('mesh');
  if (surface && mesh && surface > 0) { return faster('surface', surface, 'mesh', mesh); }
  return null;
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
      const spec = seriesSpec(rep.name);
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
      // A live bar is outlined in a dash, so a screenshot says which numbers were measured here.
      if (spec.live) { bar.setAttribute('stroke', spec.stroke || spec.color); bar.setAttribute('stroke-dasharray', '3 2'); }
      bar.appendChild(el('title', {}, `${spec.label} - ${fn.label}: ${fmt(v)} ns/call`));
      plot.appendChild(bar);
      plot.appendChild(el('text', { x: x + barW / 2, y: y - 5, 'text-anchor': 'middle', class: 'value-label' }, fmt(v, 3)));
    });

    // One ratio under each group, and WHICH ratio depends on what the part carries. Where the
    // shape the part was made from is present that is the comparison worth printing -- everything
    // else on this chart is a reconstruction of it. Failing that, the mesh against the exact
    // surfaces, which was this chart's original headline.
    const ratio = ratioFor(reps, fn.key);
    if (ratio) {
      plot.appendChild(el('text', { x: cx, y: plotH + 38, 'text-anchor': 'middle', class: 'ratio-label' },
        ratio));
    }
  });

  figure.appendChild(svg);

  const legend = document.createElement('div');
  legend.className = 'legend';
  for (const rep of present) {
    const spec = seriesSpec(rep.name);
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
      repLabel(rep.name),
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

/// The benchmark record for one part, matched on the part name the selector holds.
export function benchmarkFor(loaded, partName) {
  if (!partName) { return null; }
  return (loaded.benchmarks || []).find(e => e.doc && e.doc.meta && e.doc.meta.part === partName) || null;
}

/// Which representation the converter's cascade actually emits for this part, and where that
/// statement comes from. website_data's own verdict wins; failing that, what this checkout holds.
export function shipsVerdict(partName, entry, summary) {
  const row = summary && Array.isArray(summary.parts) ? summary.parts.find(p => p.part === partName) : null;
  if (row && row.ships) { return { ships: row.ships, source: 'website_data/summary.json' }; }
  if (entry && entry.ships) { return { ships: entry.ships, source: 'the converter that wrote testdata/' }; }
  if (!entry) { return null; }
  // Failing both, the cascade's own order decides: the first artefact this part carries is the
  // furthest the cascade got before it stopped.
  const first = PRODUCTS.find(rep => partHas(entry, rep.key));
  return { ships: first ? first.key : 'mesh', source: 'inferred from what testdata/ holds' };
}

export const SHIPS_LABEL = Object.fromEntries(REPRESENTATIONS.map(rep => [rep.key, rep.badge]));

/// Every representation this part carries at full quality, the cascade's own choice first.
///
/// A part can have more than one: `BoomCylinderInner` has both an exact sidecar and a shape.root
/// the acceptance test passed, and it ships CSG *and* SURFACE. The mesh is named only when it is
/// what the part ships as -- a tessellation sitting next to an exact sidecar is not a
/// representation the part has, it is the thing the sidecar replaced.
export function shipsKeys(primary, entry, csg) {
  const keys = [];
  if (primary) { keys.push(primary); }
  const rejected = !!(csg && csg.acceptance && csg.acceptance.accepted === false);
  // Only a PRODUCT can be what a part ships: the original is the conversion's input and the
  // cells-tree is a comparison built for this page, so neither is ever a badge.
  for (const rep of PRODUCTS) {
    if (keys.includes(rep.key) || !partHas(entry, rep.key)) { continue; }
    // The mesh is named only when it is what the part ships as (see the note above).
    if (rep.key === 'mesh') { continue; }
    // A candidate the acceptance test threw out is not a representation this part has.
    if (rep.key === 'shape' && rejected) { continue; }
    keys.push(rep.key);
  }
  return keys;
}

const SHIPS_NOTE = Object.fromEntries(REPRESENTATIONS.map(rep => [rep.key, rep.note]));

// ------------------------------------------------------------------------------------------
// the CSG record
// ------------------------------------------------------------------------------------------

/// How the recogniser's boolean reads. Only the operators the cascade can emit are named; anything
/// else is printed as the JSON's own word, so a new operator shows up rather than disappearing.
const CSG_OPERATOR = { union: ' \u222a ', intersection: ' \u2229 ', subtraction: ' \u2212 ', difference: ' \u2212 ' };

/// One line naming the composite: "TGeoTube \u222a TGeoTube", or the primitive's own class.
/// Returns null when there is no accepted candidate -- a part that is not CSG says so elsewhere.
export function csgStructure(csg) {
  const candidate = csg && csg.candidate;
  if (!candidate) { return null; }
  // A flat-DNF candidate has no leaf list at all: it is cells over halfspaces, and saying so is
  // the only honest one-liner for it.
  if (candidate.op === 'flatCells') {
    const notes = candidate.notes || {};
    const cells = notes.nCells || (candidate.cells || []).length;
    return notes.nHalfspaces
      ? `${cells} intersection-cells over ${notes.nHalfspaces} halfspaces`
      : `${cells} intersection-cells`;
  }
  if (!Array.isArray(candidate.leaves) || !candidate.leaves.length) { return null; }
  const types = candidate.leaves.map(leaf => leaf.type || '?');
  if (candidate.op === 'primitive' || types.length === 1) { return types[0]; }
  return types.join(CSG_OPERATOR[candidate.op] || ` ${candidate.op} `);
}

/// A leaf's parameters in the order TGeo takes them, with the numbers rounded to something a
/// reader can hold in their head. Angles keep their degrees.
function csgLeafParams(leaf) {
  const params = (leaf && leaf.params) || {};
  const entries = Object.entries(params).filter(([, v]) => typeof v === 'number');
  if (!entries.length) { return ''; }
  return entries.map(([k, v]) => `${k} ${fmt(v, 5)}`).join(', ');
}

/// The recogniser's record for this part, as the part card shows it: what it found, what the
/// acceptance test said about it, and whether a shape_*.root sits next to it for the raytracer.
function csgSection(card, part, csg, decline) {
  const heading = document.createElement('h4');
  heading.textContent = 'CSG, as the recogniser found it';
  card.appendChild(heading);

  const candidate = csg && csg.candidate;
  const acceptance = csg && csg.acceptance;
  const whyNot = decline && decline.whyNotCSG;

  // A flat-DNF candidate is not a tree of native primitives and must not be printed as one --
  // `(candidate.leaves || []).every(...)` is vacuously true on an empty list and would claim
  // every leaf carries a frame. It gets its own section; this one says where to look.
  if (candidate && candidate.op === 'flatCells') {
    const line = document.createElement('p');
    line.className = 'muted small';
    line.textContent = 'The recogniser did not produce a native CSG tree for this part: its cells make a ' +
      'composite wider than the routing threshold, so the cascade emitted the flat halfspace solid ' +
      'instead. What it found is in the FlatCSG section below.' +
      (partHas(part, 'cellstree')
        ? ` The cells-tree beside it is ${partNote(part, 'cellstree') || 'the decomposed cells as a plain composite'}` +
          ' -- a comparison built for this page, not what the converter ships.'
        : '');
    card.appendChild(line);
    return;
  }

  if (candidate) {
    const pairs = [
      ['structure', csgStructure(csg)],
      ['recogniser', candidate.recogniser || (csg && csg.recogniser) || null],
    ];
    (candidate.leaves || []).forEach((leaf, index) => {
      pairs.push([
        (candidate.leaves.length > 1 ? `leaf ${index + 1}` : 'primitive'),
        `${leaf.type || '?'}  ${csgLeafParams(leaf)}`.trim(),
      ]);
    });
    pairs.push([(candidate.leaves || []).length > 1 ? 'leaf placement' : 'primitive placement',
      (candidate.leaves || []).every(leaf => leaf && leaf.frame)
      ? 'every leaf carries its own world frame (origin and axes)'
      : 'not every leaf carries a frame']);
    pairs.push(['solid placement', csg.placement
      ? 'one 3 x 4 world matrix on the composite itself'
      : 'none -- the leaves are already in world coordinates']);
    if (acceptance) {
      const dv = acceptance.symmetricDifference;
      pairs.push(['acceptance', acceptance.accepted
        ? `accepted: dV_sym ${fmt(dv)} cm^3 against a band of ${fmt(acceptance.band)} cm^3`
        : `REJECTED: ${acceptance.reason || 'no reason recorded'}`]);
      pairs.push(['volume', acceptance.volumeOriginal === undefined ? null
        : `CAD ${fmt(acceptance.volumeOriginal, 6)} cm^3, candidate ${fmt(acceptance.volumeCandidate, 6)} cm^3`]);
    }
    pairs.push(['traceable here', part && part.shape
      ? `testdata/${part.shape} -- the "CSG (bridge)" raytracer view traces this very file`
      : 'no shape_*.root in testdata/ for this part']);
    card.appendChild(dl(pairs));
    return;
  }

  const line = document.createElement('p');
  line.className = 'muted small';
  if (acceptance && acceptance.reason) {
    line.textContent = 'The recogniser produced a candidate and the acceptance test threw it out: ' +
      acceptance.reason + '. The part therefore ships as something else.';
  } else if (whyNot) {
    line.textContent = `Not CSG: ${whyNot}`;
  } else if (csg) {
    line.textContent = 'The recogniser found no candidate for this part, and recorded no reason beyond that.';
  } else {
    line.textContent = 'No csg_*.json was copied for this part, and website_data/decline_reasons.json ' +
      'has no entry for it, so this page cannot say why it is not CSG.';
  }
  card.appendChild(line);
  if (whyNot && acceptance && acceptance.reason && whyNot !== acceptance.reason) {
    const second = document.createElement('p');
    second.className = 'muted small';
    second.textContent = `The cascade's own summary of the same decision: ${whyNot}`;
    card.appendChild(second);
  }
}

/// Where this part came from, when it came from an existing detector geometry rather than from
/// CAD. This is the row a reader should read first: every other number on the card is about a
/// reconstruction, and this one is the thing being reconstructed.
function baselineSection(card, part) {
  if (!partHas(part, 'original')) { return; }
  const heading = document.createElement('h4');
  heading.textContent = 'the shape this part was made from';
  card.appendChild(heading);
  card.appendChild(dl([
    ['source', partNote(part, 'original') || 'the TGeoShape this part was made from'],
    ['artefact', `testdata/${part.original}, written by scripts/geometry/exportSourceShapes.py`],
    ['what it is for', 'the baseline. It is the input to the round trip, not one of its outputs, ' +
      'so it is never counted in what the part ships -- but it is what every other representation ' +
      'here has to beat, and the "Original TGeo (bridge)" view and bars are it'],
  ]));
}

/// What the flat halfspace solid IS, read off the very bytes the kernel loads: the cell table,
/// the halfspace kinds, and the boxes the converter had to supply because an intersection of
/// halfspaces does not bound itself. The recogniser's own notes fill in what the file cannot say.
function flatCsgSection(card, part, csg, flat) {
  if (!partHas(part, 'flatcsg') && !(csg && csg.candidate && csg.candidate.op === 'flatCells')) { return; }
  const heading = document.createElement('h4');
  heading.textContent = 'FlatCSG, as the halfspace solid';
  card.appendChild(heading);

  if (!flat) {
    const line = document.createElement('p');
    line.className = 'muted small';
    line.textContent = partHas(part, 'flatcsg')
      ? 'The flat-CSG sidecar is named in the manifest but this page could not read it.'
      : 'The recogniser produced a flat-DNF candidate for this part, but no flatcsg_*.bin was copied ' +
        'into testdata/, so there is nothing here to read.';
    card.appendChild(line);
    return;
  }

  const notes = (csg && csg.candidate && csg.candidate.notes) || {};
  const kinds = [];
  if (flat.kinds.quadric) { kinds.push(`${flat.kinds.quadric} quadric`); }
  if (flat.kinds.torus) { kinds.push(`${flat.kinds.torus} torus`); }
  const pairs = [
    ['structure', `a union of ${flat.nCells} intersection-cell(s) over ${flat.nHalfspaces} signed halfspaces`],
    ['halfspaces', `${kinds.join(', ')} (${flat.signs.interior} interior, ${flat.signs.exterior} exterior)`],
    ['per cell', flat.minCellHalfspaces === flat.maxCellHalfspaces
      ? `${flat.minCellHalfspaces} halfspaces in every cell`
      : `${flat.minCellHalfspaces} to ${flat.maxCellHalfspaces} halfspaces`],
    ['cell boxes', `${flat.nCells} AABB(s) supplied by the converter -- an intersection of halfspaces ` +
      'does not bound itself, so the decomposition is the only thing that can say where a cell ends'],
    ['sidecar', `version ${flat.version}, ${(flat.byteLength / 1024).toFixed(1)} kB`],
    ['cell volume', `${fmt(flat.volume, 6)} cm^3, summed over the cells`],
  ];
  if (notes.nSplits !== undefined) { pairs.push(['decomposition splits', String(notes.nSplits)]); }
  if (notes.cellGapCm !== undefined) { pairs.push(['worst cell gap', `${fmt(notes.cellGapCm)} cm`]); }
  if (notes.cellBoxMarginCm !== undefined) {
    pairs.push(['cell-box margin', `${fmt(notes.cellBoxMarginCm)} cm`]);
  }
  if (notes.containsScored !== undefined) {
    pairs.push(['Contains corroborated at', `${notes.containsScored} point(s) against the CAD solid`]);
  }
  pairs.push(['traceable here', partHas(part, 'flatcsg')
    ? `testdata/${part.flatcsg} -- the "FlatCSG (bridge)" raytracer view traces this very file, ` +
      'and "measure this part now" benchmarks it'
    : 'no flatcsg_*.bin in testdata/ for this part']);
  if (partHas(part, 'cellstree')) {
    pairs.push(['measured against', `${partNote(part, 'cellstree') || 'the decomposed cells as a plain composite'} ` +
      `(testdata/${part.cellstree}) -- the "CSG cells-tree" view and bars are that composite`]);
  }
  if (partHas(part, 'original')) {
    pairs.push(['and against the baseline', `${partNote(part, 'original') || 'the source TGeoShape'} ` +
      `(testdata/${part.original}) -- what this part looked like before the round trip`]);
  }
  card.appendChild(dl(pairs));
}

function dl(pairs) {
  const list = document.createElement('dl');
  list.className = 'facts';
  for (const [key, value] of pairs) {
    if (value === undefined) { continue; }
    const dt = document.createElement('dt'); dt.textContent = key;
    const dd = document.createElement('dd');
    if (value && value.nodeType) { dd.appendChild(value); } else { dd.textContent = value === null ? 'n/a' : String(value); }
    if (value === null) { dd.className = 'na'; }
    list.appendChild(dt); list.appendChild(dd);
  }
  return list;
}

function bytes(n) {
  if (n === null || n === undefined) { return null; }
  if (n < 1024) { return `${n} B`; }
  if (n < 1024 * 1024) { return `${(n / 1024).toFixed(1)} kB`; }
  return `${(n / (1024 * 1024)).toFixed(2)} MB`;
}

/// The key numbers per representation: what it is made of, what it costs in memory, how far its
/// volume is from the CAD model's, and whether the kernel calls it navigable.
function representationTable(doc) {
  const reps = (doc && doc.representations) || [];
  const rows = reps.map((rep) => {
    const acc = rep.accuracy || {};
    return [
      repLabel(rep.name),
      rep.primitives === undefined ? null : `${rep.primitives} ${rep.primitiveKind || ''}`.trim(),
      bytes(rep.memoryBytes),
      bytes(rep.sidecarBytes),
      acc.capacityRelativeDeviation === undefined || acc.capacityRelativeDeviation === null
        ? null : fmt(acc.capacityRelativeDeviation),
      rep.reliability || (rep.meshClosedBody === true ? 'closed body' : rep.meshClosedBody === false ? 'NOT a closed body' : null),
    ];
  });
  return table(rows, ['representation', 'made of', 'memory', 'sidecar', 'capacity dev.', 'reliability']);
}

/// Everything the page knows about the selected part: what ships, the key numbers, the sidecar's
/// own properties, and then the per-function timing bars for this part alone.
export function partCard(entry, state, summary) {
  const doc = entry ? entry.doc : null;
  const part = state && state.part ? state.part : null;
  const partName = part ? part.name : (doc && doc.meta ? doc.meta.part : '(no part)');
  const card = document.createElement('section');
  card.className = 'card';

  const heading = document.createElement('h3');
  heading.textContent = partName;
  const verdict = shipsVerdict(partName, part, summary);
  // One badge per representation the part actually has, not one per part: the cascade's choice is
  // first and bold, and anything else at full quality stands next to it.
  const keys = verdict ? shipsKeys(verdict.ships, part, state && state.csg) : [];
  keys.forEach((key, index) => {
    const badge = document.createElement('span');
    badge.className = `badge ships ships-${key}` + (index ? ' alt' : '');
    badge.textContent = `ships ${SHIPS_LABEL[key] || key.toUpperCase()}`;
    badge.title = `${SHIPS_NOTE[key] || ''} (${index === 0 ? `${verdict.source}; what the cascade picked` : 'also carried by this part, at full quality'})`;
    heading.append(' ', badge);
  });
  // Not a representation the part ships -- a property of one of them. A tessellation of a part
  // whose every face is a planar polygon is the SAME SOLID, not an approximation, and that turns
  // the mesh from a fallback into a candidate on its merits.
  const meshExact = state && state.meshExact;
  if (meshExact && meshExact.exact && part && part.facets) {
    const badge = document.createElement('span');
    badge.className = 'badge exact';
    badge.textContent = 'tessellation is EXACT';
    badge.title = meshExact.reason;
    heading.append(' ', badge);
  }
  card.appendChild(heading);

  if (verdict) {
    const line = document.createElement('p');
    line.className = 'muted small';
    line.textContent = keys.length > 1
      ? `${SHIPS_NOTE[verdict.ships] || ''} - ${verdict.source}. This part also carries ` +
        keys.slice(1).map(k => `${SHIPS_LABEL[k]} (${SHIPS_NOTE[k] || ''})`).join(', ') +
        `, so every one of the ${keys.length} representations named above is available at full quality.`
      : `${SHIPS_NOTE[verdict.ships] || ''} - ${verdict.source}.`;
    card.appendChild(line);
  }

  if (doc && doc.meta) {
    const sub = document.createElement('p');
    sub.className = 'muted small';
    const bits = [];
    for (const [k, v] of Object.entries(doc.meta)) {
      if (k === 'part' || v === null || typeof v === 'object') { continue; }
      bits.push(`${k}: ${v}`);
    }
    bits.push(`source: ${entry.source}`);
    sub.textContent = bits.join(' \u00b7 ');
    card.appendChild(sub);
    for (const [key, value] of Object.entries(doc.meta).filter(([, v]) => v && typeof v === 'object')) {
      const line = document.createElement('p');
      line.className = 'muted small';
      line.textContent = `${key}: ` + Object.entries(value).map(([k, v]) => `${k} ${v}`).join(', ');
      card.appendChild(line);
    }
  }

  // --- what this checkout actually loaded, straight off the sidecar ------------------------
  if (part) {
    const h = document.createElement('h4');
    h.textContent = 'this part, as loaded here';
    card.appendChild(h);
    const solid = state.solid, parsed = state.parsed, facets = state.facets;
    const records = document.createElement('span');
    if (!solid) {
      records.className = 'badge warn';
      records.textContent = 'tessellated only -- no exact sidecar';
    } else if (solid.failed.length || solid.unsupported.length) {
      records.className = 'badge bad';
      records.textContent = `${solid.failed.length} rejected, ${solid.unsupported.length} unsupported`;
      records.title = [...solid.failed, ...solid.unsupported].map(f => `#${f.index}: ${f.reason}`).join('\n');
    } else {
      records.className = 'badge ok';
      records.textContent = 'all records built';
    }
    const meshExactness = state.meshExact;
    card.appendChild(dl([
      ['exact faces', solid ? String(solid.nSurfaces) : null],
      ['by type', solid ? (Object.entries(solid.counts).map(([k, n]) => `${n} ${k}`).join(', ') || '-') : null],
      ['triangles', facets ? String(facets.nTriangles) : null],
      ['tessellation', !(meshExactness && facets) ? null : (meshExactness.exact
        ? `EXACT -- ${meshExactness.reason}. The mesh is not an approximation of this part, it is ` +
          'the same solid, so choosing it costs no accuracy at all'
        : `an approximation -- ${meshExactness.reason}`)],
      ['sidecar', parsed ? `version ${parsed.version}, ${(parsed.byteLength / 1024).toFixed(1)} kB` : null],
      ['wire-trimmed', solid ? `${solid.wireTrimFaces} face(s)` : null],
      ['B-spline trims', solid ? `${solid.bsplineTrimFaces} face(s)` : null],
      ['model tolerance', parsed ? (parsed.modelToleranceStated ? `${parsed.modelTolerance.toExponential(2)} cm` : 'not stated (v1)') : null],
      ['worst join gap', solid ? `${solid.worstJoinGap.toExponential(2)} cm (band ${solid.joinTolerance.toExponential(1)})` : null],
      ['records', records],
    ]));
  }

  if (part) { baselineSection(card, part); }
  if (part) { csgSection(card, part, state.csg, state.decline); }
  if (part) { flatCsgSection(card, part, state.csg, state.flatcsg); }

  // Whatever the bridge measured on this machine, this run. It is timing only, so it joins the
  // bars and nothing else; the Track-2 record stays the authority on accuracy and the X-ray counts.
  const live = liveResults.get(partName) || null;
  const liveReps = live ? live.reps.map(rep => ({
    name: `${rep.key}-live`,
    functions: rep.functions,
  })) : [];
  const staticReps = (doc && doc.representations) || [];

  if (!doc) {
    const missing = document.createElement('p');
    missing.className = 'muted small';
    missing.textContent = 'No measured Track-2 record for this part, so there is no memory model and there are no ' +
      'X-ray counters here. The numbers above come from the sidecar this page loaded' +
      (liveReps.length ? ', and the bars below were measured live on this machine through the bridge.'
        : '; press "measure now (bridge)" above to time it on this machine.');
    card.appendChild(missing);
  } else {
    const h1 = document.createElement('h4');
    h1.textContent = 'representations, as measured';
    card.appendChild(h1);
    card.appendChild(representationTable(doc));
  }

  if (staticReps.length || liveReps.length) {
    const h2 = document.createElement('h4');
    h2.textContent = 'ns per call, per navigation function';
    card.appendChild(h2);
    card.appendChild(barChart({ representations: [...staticReps, ...liveReps] }));
    if (live) {
      const caption = document.createElement('p');
      caption.className = 'muted small';
      caption.textContent = `The dashed bars were measured live on this machine, load avg ` +
        `${fmt(live.loadAverage, 2)}, at ${live.when.toLocaleTimeString()}: ${live.samples} deterministic samples ` +
        `x ${live.repeats} repeats per function, single-threaded, on the shape the bridge itself loaded ` +
        `(${live.reps.map(r => `${repLabel(r.key)} from ${r.path}`).join('; ')}). ` +
        (live.reps.length > 1
          ? 'Every subject was put the SAME sample points -- the artefacts of one part do not share a ' +
            `bounding box, so the ${repLabel(live.boxFrom || 'surface')} box was passed to each /bench. ` +
            `${live.insideSamples} of them are interior points, ` +
            `${live.reps[0].insideFromEntry || 0} of those harvested from where a ray enters the solid, which ` +
            'is what makes DistFromInside measurable on a thin part at all. '
          : '') +
        (staticReps.length
          ? 'The solid bars are the Track-2 record, taken on the machine named in website_data/summary.json. '
          : '') +
        'A live number is timing only: it says nothing about accuracy.';
      card.appendChild(caption);
    }
  }

  if (doc) {
    const h3 = document.createElement('h4');
    h3.textContent = 'accuracy and X-ray counters';
    card.appendChild(h3);
    card.appendChild(countersTable(doc));
  }
  return card;
}

// What the last measurement of a part said. The pane is rebuilt when the result lands -- the
// bars have to be redrawn -- so the sentence has to outlive the element that wrote it.
const liveStatus = new Map();

// Whether the reader has asked for the opt-in comparison subjects. Lives here rather than in the
// pane, for the same reason: the pane is thrown away and rebuilt on every measurement.
let wantComparisons = false;

/// The "measure now" pane: the button, what it will do, and what it cannot do.
function livePane(state, rerender) {
  const part = state && state.part;
  const pane = document.createElement('div');
  pane.className = 'pane';
  const heading = document.createElement('h3');
  heading.textContent = 'measure this part now, on this machine';
  pane.appendChild(heading);

  const explain = document.createElement('p');
  explain.className = 'muted small';
  explain.innerHTML = 'The bars above are one machine\'s numbers, taken once. This button asks the bridge ' +
    '(<code>scripts/geometry/tgeoRayService.py</code>) to <code>/load</code> each representation of this part ' +
    'and <code>/bench</code> it on the spot: single-threaded, deterministic sample points, the real O2 kernel. ' +
    'The result joins the chart as a dashed bar next to the recorded one. It is <strong>timing only</strong>; ' +
    'accuracy, capacity and the X-ray counters stay with the Track-2 record. The bridge is put back to the file ' +
    'the raytracer tab expects afterwards.';
  pane.appendChild(explain);

  const row = document.createElement('div');
  row.className = 'row';
  const button = document.createElement('button');
  button.className = 'primary';
  button.textContent = 'measure now (bridge)';
  button.disabled = true;
  const status = document.createElement('span');
  status.className = 'status';
  status.textContent = 'looking for the bridge ...';
  row.append(button, status);
  pane.appendChild(row);

  const plan = measurablePlan(part, { comparisons: wantComparisons });
  const coverage = document.createElement('p');
  coverage.className = 'muted small';
  coverage.textContent = plan.can.length
    ? `Measurable here: ${plan.can.map(m => repLabel(m.key)).join(', ')}. ` +
      'Each one is a separate /load of a separate file into the same kernel, benchmarked on the same ' +
      'deterministic sample points, so the bars differ only by which representation answered.'
    : 'This part carries none of the representation artefacts, so there is nothing the bridge can load for it.';
  pane.appendChild(coverage);

  // The opt-in subjects, named for what they are rather than hidden. A part with none of them
  // gets no checkbox at all rather than a dead one.
  if (plan.optional.length || wantComparisons) {
    const optional = measurablePlan(part, { comparisons: true }).can
      .filter(m => m.role === 'comparison');
    if (optional.length) {
      const row = document.createElement('div');
      row.className = 'row';
      const label = document.createElement('label');
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.checked = wantComparisons;
      box.addEventListener('change', () => { wantComparisons = box.checked; rerender(); });
      label.appendChild(box);
      label.appendChild(document.createTextNode(
        ` also measure ${optional.map(m => repLabel(m.key)).join(', ')}`));
      row.appendChild(label);
      pane.appendChild(row);
      const why = document.createElement('p');
      why.className = 'muted small';
      why.textContent = 'Off by default. The cells-tree is not a representation this part ships -- ' +
        'it is the decomposition emitted as a plain composite, and it exists so the flat solid has ' +
        'something to be measured against. Turn it on to reproduce that specific ratio; leave it off ' +
        'to compare the shape this part was made from against what the converter can ship.';
      pane.appendChild(why);
    }
  }

  if (!part || !plan.can.length) {
    status.textContent = 'nothing to measure';
    return pane;
  }

  const outcome = liveStatus.get(part.name);
  if (outcome) {
    status.textContent = outcome.text;
    status.className = outcome.error ? 'status error' : 'status';
  }

  const port = defaultPort();
  probeBridge(port).then((probe) => {
    if (!probe.ok) {
      button.disabled = true;
      button.title = probe.error;
      status.textContent = probe.error;
      status.className = 'status error';
      return;
    }
    button.disabled = false;
    if (!liveStatus.has(part.name)) {
      status.textContent = `bridge answering on 127.0.0.1:${port}`;
      status.className = 'status';
    }
  });

  button.addEventListener('click', async () => {
    button.disabled = true;
    status.className = 'status';
    liveStatus.delete(part.name);
    try {
      const record = await measurePart(part, { port, comparisons: wantComparisons,
                                              onStep: (line) => { status.textContent = line; } });
      liveStatus.set(part.name, {
        text: `measured ${record.reps.length} representation(s) at load avg ${fmt(record.loadAverage, 2)}` +
          (record.restored ? `; the bridge is back on ${record.restored}` : ''),
      });
      rerender();     // rebuilds this pane, which reads the sentence back out of liveStatus
    } catch (e) {
      liveStatus.set(part.name, { text: `measurement failed: ${e.message}`, error: true });
      status.textContent = `measurement failed: ${e.message}`;
      status.className = 'status error';
      button.disabled = false;
    }
  });
  return pane;
}

export function renderBenchmarks(container, loaded, state) {
  container.innerHTML = '';
  const partName = state && state.part ? state.part.name : null;
  const note = document.createElement('p');
  note.className = 'muted';
  if (loaded.origin === 'none') {
    note.innerHTML = 'No benchmark JSON found. Track 2 writes <code>scripts/geometry/website_data/*.json</code>; ' +
      '<code>sample_data/</code> holds a schema-shaped stand-in. The part card below is read from the sidecar.';
  } else {
    note.innerHTML = loaded.origin === 'website_data'
      ? 'Measured data from <code>scripts/geometry/website_data/</code>, for the selected part only.'
      : 'No measured data present yet: these are the <strong>synthetic sample records</strong> in ' +
        '<code>sample_data/</code>, shaped exactly like the Track-2 schema.';
  }
  container.appendChild(note);
  // The button re-renders this whole tab when it has a result, so the bars pick the live record up.
  container.appendChild(livePane(state, () => renderBenchmarks(container, loaded, state)));
  container.appendChild(partCard(benchmarkFor(loaded, partName), state, loaded.summary));
}
