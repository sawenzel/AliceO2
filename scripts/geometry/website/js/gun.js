// A synthetic particle gun, so the event display has something to play before real o2-sim output
// exists -- and so the page can raise the statistics on demand instead of shipping a huge file.
//
// The transport is NOT Geant4. It is a straight or gently helical step with two crude material
// effects (multiple scattering and absorption), and the only thing that is real about it is the
// material: a step counts as "inside" when the exact solid's own parity test says the point is
// inside. That is enough for the two things the display is for -- the step cloud thickens where
// the part is, and the part throws a shadow on the screen -- and the page never claims more.

import { CONTAINS_TEST_DIRECTION } from './solid.js';

export const SPECIES = [
  { pdg: 22, charge: 0, energyGeV: 1.0, weight: 0.4 },      // photon: straight
  { pdg: 11, charge: -1, energyGeV: 0.6, weight: 0.22 },    // electron: bends one way
  { pdg: 211, charge: 1, energyGeV: 2.0, weight: 0.26 },    // pion: bends the other
  { pdg: 2212, charge: 1, energyGeV: 3.0, weight: 0.12 },   // proton: stiff
];

/// mulberry32. A textbook LCG is wrong here: its multiply overflows 2^53 in double arithmetic,
/// which throws away the low bits it lives on, and the "random" aim then lands on a coarse
/// lattice -- which is exactly what a radiograph makes visible.
export function makeRandom(seed) {
  let state = seed | 0;
  return function random() {
    state = (state + 0x6D2B79F5) | 0;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/// The default gun and screen geometry for a part: the gun on -x, the screen behind it on +x,
/// covering the part with a margin.
export function defaultSetup(box) {
  const center = [(box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2];
  const radius = 0.5 * Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
  const halfY = Math.max(1e-3, 0.5 * (box[4] - box[1]));
  const halfZ = Math.max(1e-3, 0.5 * (box[5] - box[2]));
  return {
    center, radius,
    gun: [center[0] - 3 * radius, center[1], center[2]],
    screen: {
      origin: [center[0] + 2.2 * radius, center[1], center[2]],
      normal: [1, 0, 0],
      up: [0, 0, 1],
      halfWidth: 1.3 * halfY,
      halfHeight: 1.3 * halfZ,
      pixels: [120, 120],
    },
  };
}

export function generateEvents(solid, options = {}) {
  const {
    events: eventCount = 5,
    tracksPerEvent = 300,
    seed = 20260822,
    scatterPerStep = 0.06,      // rad, rms, per step in material
    absorptionPerStep = 0.05,   // probability of stopping, per step in material
  } = options;

  const box = solid.aabb;
  const setup = defaultSetup(box);
  const { center, radius, gun, screen } = setup;
  const random = makeRandom(seed);

  const gauss = () => {
    const u1 = Math.max(1e-9, random()), u2 = random();
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  };
  const pickSpecies = () => {
    const r = random();
    let acc = 0;
    for (const s of SPECIES) { acc += s.weight; if (r <= acc) { return s; } }
    return SPECIES[SPECIES.length - 1];
  };

  const transport = (species) => {
    const aim = [screen.origin[0],
                 center[1] + (2 * random() - 1) * screen.halfWidth,
                 center[2] + (2 * random() - 1) * screen.halfHeight];
    let dir = [aim[0] - gun[0], aim[1] - gun[1], aim[2] - gun[2]];
    const length = Math.hypot(dir[0], dir[1], dir[2]);
    dir = dir.map(c => c / length);

    const points = [gun.slice()];
    let position = gun.slice();
    const vacuumStep = radius * 0.5;
    const materialStep = radius * 0.04;
    // a gentle bend: about 0.1 rad over the whole flight for a 1 GeV charged track
    const kappa = species.charge === 0 ? 0 : species.charge * 0.02 / species.energyGeV / Math.max(radius, 1e-6);
    const maxRange = 10 * radius;
    let travelled = 0;
    let lastInside = false;
    while (travelled < maxRange) {
      const inside = solid.containsAlong(position[0], position[1], position[2],
        CONTAINS_TEST_DIRECTION[0], CONTAINS_TEST_DIRECTION[1], CONTAINS_TEST_DIRECTION[2]);
      const step = inside ? materialStep : vacuumStep;
      if (kappa !== 0) {
        const angle = kappa * step;
        const cos = Math.cos(angle), sin = Math.sin(angle);
        dir = [dir[0] * cos - dir[1] * sin, dir[0] * sin + dir[1] * cos, dir[2]];
      }
      if (inside) {
        dir = [dir[0], dir[1] + scatterPerStep * gauss() * 0.35, dir[2] + scatterPerStep * gauss() * 0.35];
        if (random() < absorptionPerStep) {
          points.push(position.slice());
          return { points, absorbed: true };
        }
      }
      const n = Math.hypot(dir[0], dir[1], dir[2]);
      dir = dir.map(c => c / n);
      position = [position[0] + dir[0] * step, position[1] + dir[1] * step, position[2] + dir[2] * step];
      // a vertex at every material step and at every material boundary; a straight neutral track
      // in vacuum needs no intermediate points
      if (inside || lastInside || kappa !== 0 || points.length < 2) { points.push(position.slice()); }
      lastInside = inside;
      travelled += step;
      if (position[0] > screen.origin[0] + 0.15 * radius) { points.push(position.slice()); break; }
    }
    return { points, absorbed: false };
  };

  const eventList = [];
  let absorbed = 0;
  for (let e = 0; e < eventCount; ++e) {
    const tracks = [];
    for (let t = 0; t < tracksPerEvent; ++t) {
      const species = pickSpecies();
      const result = transport(species);
      if (result.absorbed) { absorbed += 1; }
      tracks.push({
        pdg: species.pdg,
        charge: species.charge,
        energyGeV: Number((species.energyGeV * (0.6 + 0.8 * random())).toFixed(3)),
        points: result.points.map(p => p.map(c => Number(c.toFixed(3)))),
      });
    }
    eventList.push({ tracks });
  }

  return {
    meta: {
      part: solid.label,
      generator: 'website/js/gun.js (synthetic gun: straight/helical steps with crude scattering and ' +
                 'absorption; the material is the real exact solid, the transport is not Geant4)',
      seed,
      representation: 'exact',
      synthetic: true,
      tracksPerEvent,
      absorbed,
    },
    screen,
    events: eventList,
  };
}
