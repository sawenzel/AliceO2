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
import { SurfaceSolid, CONTAINS_TEST_DIRECTION } from '../js/solid.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, '..');
const outDir = path.join(root, 'sample_data');
fs.mkdirSync(outDir, { recursive: true });

// A deterministic generator, so re-running produces the same file.
let seedState = 20260822;
function rnd() { seedState = (seedState * 1103515245 + 12345) & 0x7fffffff; return seedState / 0x7fffffff; }

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

const box = solid.aabb;
const center = [(box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2];
const radius = 0.5 * Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]);

// The gun sits on -x and fires towards +x; the screen plane sits behind the part on +x.
const gun = [center[0] - 3 * radius, center[1], center[2]];
const screen = {
  origin: [center[0] + 2.2 * radius, center[1], center[2]],
  normal: [1, 0, 0],
  up: [0, 0, 1],
  halfWidth: 1.4 * radius,
  halfHeight: 1.4 * radius,
  pixels: [220, 220],
};

const SPECIES = [
  { pdg: 22, charge: 0, energyGeV: 1.0, weight: 0.35 },     // photon: straight
  { pdg: 11, charge: -1, energyGeV: 0.6, weight: 0.3 },     // electron: bends
  { pdg: 211, charge: 1, energyGeV: 2.0, weight: 0.25 },    // pion: bends the other way
  { pdg: 2212, charge: 1, energyGeV: 3.0, weight: 0.1 },    // proton: stiff
];

function pickSpecies() {
  const r = rnd();
  let acc = 0;
  for (const s of SPECIES) { acc += s.weight; if (r <= acc) { return s; } }
  return SPECIES[SPECIES.length - 1];
}

/// Step a track from the gun towards the part, bending it in a uniform field along +z, and record
/// the polyline. Steps inside material are shortened, which is what makes the point cloud show
/// where the material is.
function trackPoints(species) {
  // aim at a random point on the part's bounding sphere face
  const aim = [center[0], center[1] + (2 * rnd() - 1) * radius * 0.9, center[2] + (2 * rnd() - 1) * radius * 0.9];
  let dir = [aim[0] - gun[0], aim[1] - gun[1], aim[2] - gun[2]];
  const len = Math.hypot(...dir);
  dir = dir.map(c => c / len);

  const points = [gun.slice()];
  let position = gun.slice();
  const vacuumStep = radius * 0.25;
  const materialStep = radius * 0.03;
  // curvature ~ charge / momentum, in the x-y plane (field along z)
  const kappa = species.charge === 0 ? 0 : species.charge * 0.35 / species.energyGeV / Math.max(radius, 1e-6);
  const maxRange = 8 * radius;
  let travelled = 0;
  let insideBefore = false;
  while (travelled < maxRange) {
    const inside = solid.containsAlong(position[0], position[1], position[2],
      CONTAINS_TEST_DIRECTION[0], CONTAINS_TEST_DIRECTION[1], CONTAINS_TEST_DIRECTION[2]);
    const step = inside ? materialStep : vacuumStep;
    if (kappa !== 0) {
      const angle = kappa * step;
      const cos = Math.cos(angle), sin = Math.sin(angle);
      const nx = dir[0] * cos - dir[1] * sin;
      const ny = dir[0] * sin + dir[1] * cos;
      dir = [nx, ny, dir[2]];
      const n = Math.hypot(...dir);
      dir = dir.map(c => c / n);
    }
    position = [position[0] + dir[0] * step, position[1] + dir[1] * step, position[2] + dir[2] * step];
    points.push(position.slice());
    travelled += step;
    if (inside !== insideBefore) { insideBefore = inside; }
    if (position[0] > screen.origin[0] + 0.1 * radius) { break; }
  }
  return points;
}

const events = [];
for (let e = 0; e < 6; ++e) {
  const tracks = [];
  const nTracks = 6 + Math.floor(rnd() * 6);
  for (let t = 0; t < nTracks; ++t) {
    const species = pickSpecies();
    tracks.push({
      pdg: species.pdg,
      charge: species.charge,
      energyGeV: Number((species.energyGeV * (0.6 + 0.8 * rnd())).toFixed(3)),
      points: trackPoints(species).map(p => p.map(c => Number(c.toFixed(4)))),
    });
  }
  events.push({ tracks });
}

const eventsDoc = {
  meta: {
    part: eventPart,
    generator: 'website/tools/make_sample_data.mjs (synthetic gun; transport is straight/helical, not Geant4)',
    seed: 20260822,
    representation: 'exact',
    synthetic: true,
  },
  screen,
  events,
};
fs.writeFileSync(path.join(outDir, 'events_sample.json'), JSON.stringify(eventsDoc) + '\n');

console.log(`wrote ${files.length} benchmark records and events_sample.json (${eventPart}, ${events.length} events, ` +
            `${events.reduce((n, ev) => n + ev.tracks.length, 0)} tracks) into ${outDir}`);
