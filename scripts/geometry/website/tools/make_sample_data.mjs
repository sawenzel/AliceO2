// Generate the committed stand-in data the page falls back on when no measured data is present:
//
//   sample_data/index.json          the list of benchmark records
//   sample_data/<part>.json         one benchmark record per part, in the Track-2 schema
//   sample_data/events_sample.json  a synthetic event replay through one test part
//
// The benchmark numbers are SYNTHETIC and labelled as such in their own meta block; they exist so
// the chart code can be developed and reviewed before the measured run lands, and the page says
// so wherever it draws them. The event replay is synthetic too, but its geometry is real: the
// tracks are transported through the actual exact solid, and a step is recorded where the parity
// test says the point is inside material.
//
//   node tools/make_sample_data.mjs [part-for-events]

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { parseSidecar } from '../js/sidecar.js';
import { SurfaceSolid } from '../js/solid.js';
import { generateEvents } from '../js/gun.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, '..');
const outDir = path.join(root, 'sample_data');
fs.mkdirSync(outDir, { recursive: true });

// A deterministic generator, so re-running produces the same file. mulberry32 rather than a
// textbook LCG: the LCG's multiply overflows 2^53 in double arithmetic, which silently throws
// away the low bits it lives on, and the "random" aim then lands on a coarse lattice.
let seedState = 20260822;
function rnd() {
  seedState = (seedState + 0x6D2B79F5) | 0;
  let t = seedState;
  t = Math.imul(t ^ (t >>> 15), t | 1);
  t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}

// ---------------------------------------------------------------------------------------------
// benchmark records
// ---------------------------------------------------------------------------------------------

function record(part, { faces, triangles, surfaceNs, meshFactor, shapeNs, xray, capacity }) {
  const scale = (base, factor) => Number((base * factor).toFixed(1));
  const fns = (base) => ({
    contains: { nsPerCall: scale(base, 1) },
    distFromOutside: { nsPerCall: scale(base, 2.4) },
    distFromInside: { nsPerCall: scale(base, 2.1) },
    safety: { nsPerCall: scale(base, 1.6) },
  });
  const representations = [
    {
      name: 'surface',
      functions: fns(surfaceNs),
      accuracy: { capacityRelativeDeviation: capacity, disagreements: {}, unexplained: {} },
      xray: xray.surface,
    },
    {
      name: 'mesh',
      functions: fns(surfaceNs * meshFactor),
      accuracy: { capacityRelativeDeviation: null, disagreements: {}, unexplained: {} },
      xray: xray.mesh,
    },
  ];
  if (shapeNs !== null) {
    representations.push({
      name: 'shape',
      functions: fns(shapeNs),
      accuracy: { capacityRelativeDeviation: 0, disagreements: {}, unexplained: {} },
      xray: xray.shape || null,
    });
  }
  return {
    meta: {
      part,
      synthetic: true,
      note: 'SYNTHETIC placeholder in the Track-2 schema; replace with scripts/geometry/website_data/*.json',
      faces, triangles,
    },
    representations,
  };
}

const zeroXray = (rays) => ({ lost: 0, extra: 0, unterminated: 0, parity: 0, raysTotal: rays });

const records = [
  record('box', {
    faces: 6, triangles: 12, surfaceNs: 41, meshFactor: 3.2, shapeNs: 12,
    capacity: 0, xray: { surface: zeroXray(4096), mesh: zeroXray(4096), shape: zeroXray(4096) },
  }),
  record('cyl_inter_cyl', {
    faces: 6, triangles: 604, surfaceNs: 96, meshFactor: 41, shapeNs: null,
    capacity: 2.8e-7,
    xray: {
      surface: zeroXray(4096),
      // the part whose mesh leaks: rays enter the tessellation and never come out
      mesh: { lost: 6, extra: 0, unterminated: 6, parity: 6, raysTotal: 4096 },
    },
  }),
  record('torus_union_cyl', {
    faces: 6, triangles: 23432, surfaceNs: 128, meshFactor: 145, shapeNs: 31,
    capacity: 4.1e-8,
    xray: { surface: zeroXray(4096), mesh: { lost: 0, extra: 2, unterminated: 0, parity: 0, raysTotal: 4096 }, shape: zeroXray(4096) },
  }),
  record('tube_window', {
    faces: 4, triangles: 1248, surfaceNs: 88, meshFactor: 52, shapeNs: null,
    capacity: 1.6e-7, xray: { surface: zeroXray(4096), mesh: zeroXray(4096) },
  }),
  record('Bucket', {
    faces: 97, triangles: 9556, surfaceNs: 210, meshFactor: 68, shapeNs: null,
    capacity: 3.3e-7, xray: { surface: zeroXray(4096), mesh: { lost: 0, extra: 0, unterminated: 0, parity: 0, raysTotal: 4096 } },
  }),
  record('BoomCylinderInner', {
    faces: 6, triangles: 1612, surfaceNs: 79, meshFactor: 44, shapeNs: null,
    capacity: 2.8e-7, xray: { surface: zeroXray(4096), mesh: zeroXray(4096) },
  }),
];

// Safety is not measured for every representation of every part; leave one hole so the page's
// "n/a, never zero" rule is exercised by the committed data itself.
delete records[1].representations[1].functions.safety;

const files = [];
for (const rec of records) {
  const name = `bench_${rec.meta.part}.json`;
  fs.writeFileSync(path.join(outDir, name), JSON.stringify(rec, null, 2) + '\n');
  files.push(name);
}
fs.writeFileSync(path.join(outDir, 'index.json'), JSON.stringify({ synthetic: true, files }, null, 2) + '\n');

// ---------------------------------------------------------------------------------------------
// event replay
// ---------------------------------------------------------------------------------------------

const eventPart = process.argv[2] || 'torus_union_cyl';
const sidecarPath = path.join(root, 'testdata', eventPart, 'surfaces.bin');
if (!fs.existsSync(sidecarPath)) {
  console.error(`no ${sidecarPath}; run ./fetch_testdata.sh first (benchmark samples were written)`);
  process.exit(0);
}
const buffer = fs.readFileSync(sidecarPath);
const solid = new SurfaceSolid(
  parseSidecar(buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength), eventPart), eventPart);

// The committed sample is deliberately small -- it demonstrates the mechanism, and the page can
// regenerate a high-statistics replay in a worker whenever a sharper radiograph is wanted.
const eventsDoc = generateEvents(solid, {
  events: Number(process.argv[3] || 5),
  tracksPerEvent: Number(process.argv[4] || 300),
  seed: 20260822,
});
fs.writeFileSync(path.join(outDir, 'events_sample.json'), JSON.stringify(eventsDoc) + '\n');

const trackCount = eventsDoc.events.reduce((n, ev) => n + ev.tracks.length, 0);
const pointCount = eventsDoc.events.reduce((n, ev) => n + ev.tracks.reduce((m, t) => m + t.points.length, 0), 0);
console.log(`wrote ${files.length} benchmark records and events_sample.json (${eventPart}, ` +
            `${eventsDoc.events.length} events, ${trackCount} tracks, ${pointCount} points, ` +
            `${eventsDoc.meta.absorbed} absorbed) into ${outDir}`);
