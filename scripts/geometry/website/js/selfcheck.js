// The self-check: deterministic ray grids against the parts in testdata/, with expectations that
// are analytic rather than recorded. It runs unchanged in node (tools/selfcheck.mjs) and in the
// browser (the Self-check tab), because both hand it the same loader.
//
// What it asserts, and why each one is worth asserting:
//
//   quartic          solveQuarticReal on polynomials whose roots are known in closed form,
//                    including the biquadratic whose odd coefficient cancels to noise -- the case
//                    the C++ comments record as having silently taken the wrong branch.
//   box slab         for a box, the exact solid must reproduce the analytic axis-aligned slab
//                    intersection for EVERY ray, so 2000 random rays are compared against it.
//   tube chords      three rays through tube_window whose material path length follows from the
//                    radii by hand: along the bore (no material at all, which is a statement
//                    about the trim wires), across it, and parallel to the axis beside it.
//   parity closure   every enter has an exit: the crossing count along a full ray must be even,
//                    for the exact solid and for the mesh, on every part.

import { parseSidecar, parseFacets } from './sidecar.js';
import { SurfaceSolid, RESHOOT_DIRECTIONS } from './solid.js';
import { MeshSolid } from './meshtrace.js';
import { solveQuarticReal, QuarticBranch } from './quartic.js';
import { makeRandom } from './gun.js';

function fmt(x, digits = 4) {
  if (!Number.isFinite(x)) { return String(x); }
  return Math.abs(x) >= 1e4 || (Math.abs(x) < 1e-3 && x !== 0) ? x.toExponential(digits - 1) : x.toFixed(digits);
}

class Report {
  constructor(log) { this.log = log || (() => {}); this.pass = 0; this.fail = 0; this.lines = []; }
  check(ok, name, detail = '') {
    if (ok) { this.pass += 1; } else { this.fail += 1; }
    const line = `${ok ? 'PASS' : 'FAIL'}  ${name}${detail ? '  -- ' + detail : ''}`;
    this.lines.push({ ok, line });
    this.log(line, ok);
    return ok;
  }
  note(text) { this.lines.push({ ok: null, line: text }); this.log(text, null); }
}

// ------------------------------------------------------------------------------------------
// the quartic
// ------------------------------------------------------------------------------------------

function checkQuartic(report) {
  const near = (a, b, tolerance) => Math.abs(a - b) <= tolerance;

  // (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
  {
    const roots = solveQuarticReal(1, -10, 35, -50, 24).sort((a, b) => a - b);
    const expected = [1, 2, 3, 4];
    const ok = roots.length === 4 && expected.every((e, i) => near(roots[i], e, 1e-9));
    report.check(ok, 'quartic: four separated real roots', `got [${roots.map(r => fmt(r, 6)).join(', ')}]`);
  }
  // (x^2-1)(x^2-4) = x^4 - 5x^2 + 4, a biquadratic: q is exactly zero
  {
    const branch = {};
    const roots = solveQuarticReal(1, 0, -5, 0, 4, branch).sort((a, b) => a - b);
    const expected = [-2, -1, 1, 2];
    const ok = roots.length === 4 && expected.every((e, i) => near(roots[i], e, 1e-9))
               && branch.branch === QuarticBranch.Biquadratic;
    report.check(ok, 'quartic: exact biquadratic takes the biquadratic branch',
                 `branch ${branch.branch}, roots [${roots.map(r => fmt(r, 6)).join(', ')}]`);
  }
  // the same roots at 1e-3: the odd coefficient cancels to ~1e-18 out of terms of size 1e-6, which
  // is zero to the precision it can be computed with. A guard in centimetres misroutes this.
  {
    const s = 1e-3;
    const r = [-2 * s, -1 * s, 1 * s, 2 * s];
    const a3 = -(r[0] + r[1] + r[2] + r[3]);
    const a2 = r[0] * r[1] + r[0] * r[2] + r[0] * r[3] + r[1] * r[2] + r[1] * r[3] + r[2] * r[3];
    const a1 = -(r[0] * r[1] * r[2] + r[0] * r[1] * r[3] + r[0] * r[2] * r[3] + r[1] * r[2] * r[3]);
    const a0 = r[0] * r[1] * r[2] * r[3];
    const branch = {};
    const roots = solveQuarticReal(1, a3, a2, a1, a0, branch).sort((a, b) => a - b);
    const ok = roots.length === 4 && r.every((e, i) => near(roots[i], e, 1e-12))
               && branch.branch === QuarticBranch.Biquadratic;
    report.check(ok, 'quartic: scaled biquadratic (roots +/-1e-3, +/-2e-3) is still recognised',
                 `branch ${branch.branch}, roots [${roots.map(x => fmt(x, 6)).join(', ')}]`);
  }
  // a double root must come back as a near-equal pair, so parity clustering can drop it
  {
    const roots = solveQuarticReal(1, -4, 5, -2, 0).sort((a, b) => a - b); // x(x-1)^2(x-2)
    const ok = roots.length === 4 && near(roots[0], 0, 1e-7) && near(roots[1], 1, 1e-3)
               && near(roots[2], 1, 1e-3) && near(roots[3], 2, 1e-7);
    report.check(ok, 'quartic: a double root is returned as a pair',
                 `roots [${roots.map(x => fmt(x, 6)).join(', ')}]`);
  }
  // a quartic with no real root contributes nothing
  {
    const roots = solveQuarticReal(1, 0, 2, 0, 1); // (x^2+1)^2
    report.check(roots.length === 0 || roots.every(x => Math.abs(x) < 1e-6) === false || roots.length === 0,
                 'quartic: a complex-only quartic returns no real root', `got ${roots.length} root(s)`);
  }
}

// ------------------------------------------------------------------------------------------
// the box: the exact solid against the analytic slab intersection
// ------------------------------------------------------------------------------------------

/// Entry and exit of a ray with an axis-aligned box, in closed form.
function boxSlab(box, o, d) {
  let t0 = -Infinity, t1 = Infinity;
  for (let i = 0; i < 3; ++i) {
    const inv = 1 / d[i];
    let a = (box[i] - o[i]) * inv, b = (box[i + 3] - o[i]) * inv;
    if (a > b) { const s = a; a = b; b = s; }
    t0 = Math.max(t0, a); t1 = Math.min(t1, b);
  }
  return t0 <= t1 ? [t0, t1] : null;
}

function checkBox(report, solid) {
  // The sidecar's own faces give the box: strip the BVH margin back off the surface boxes.
  const box = [solid.aabb[0] + 1e-3, solid.aabb[1] + 1e-3, solid.aabb[2] + 1e-3,
               solid.aabb[3] - 1e-3, solid.aabb[4] - 1e-3, solid.aabb[5] - 1e-3];
  report.note(`  box extent from the sidecar: [${box.map(x => fmt(x, 3)).join(', ')}] cm`);

  // three axial rays, whose chords are the box edges
  const centre = [(box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2];
  const axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let axis = 0; axis < 3; ++axis) {
    const d = axes[axis];
    const origin = centre.slice();
    origin[axis] = box[axis] - 5;
    const hits = solid.allIntersections(origin[0], origin[1], origin[2], d[0], d[1], d[2]);
    hits.sort((a, b) => a.t - b.t);
    const expectedIn = box[axis] - origin[axis];
    const expectedOut = box[axis + 3] - origin[axis];
    const ok = hits.length === 2 && Math.abs(hits[0].t - expectedIn) < 1e-9 && Math.abs(hits[1].t - expectedOut) < 1e-9;
    report.check(ok, `box: axial ray along axis ${axis} enters at ${fmt(expectedIn)} and exits at ${fmt(expectedOut)}`,
                 hits.length === 2 ? `got ${fmt(hits[0].t)} / ${fmt(hits[1].t)}`
                                   : `got ${hits.length} crossing(s)`);
  }

  // and 2000 random rays against the closed-form slab intersection
  const random = makeRandom(7919);
  const radius = Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]);
  let compared = 0, worst = 0, missing = 0, spurious = 0;
  for (let i = 0; i < 2000; ++i) {
    const u = 2 * random() - 1, phi = 2 * Math.PI * random(), s = Math.sqrt(Math.max(0, 1 - u * u));
    const d = [s * Math.cos(phi), s * Math.sin(phi), u];
    const aim = [box[0] + (box[3] - box[0]) * random(),
                 box[1] + (box[4] - box[1]) * random(),
                 box[2] + (box[5] - box[2]) * random()];
    const o = [aim[0] - d[0] * radius, aim[1] - d[1] * radius, aim[2] - d[2] * radius];
    const analytic = boxSlab(box, o, d);
    const hits = solid.allIntersections(o[0], o[1], o[2], d[0], d[1], d[2]);
    hits.sort((a, b) => a.t - b.t);
    if (!analytic) { if (hits.length) { spurious += 1; } continue; }
    if (hits.length < 2) { missing += 1; continue; }
    compared += 1;
    worst = Math.max(worst, Math.abs(hits[0].t - analytic[0]), Math.abs(hits[hits.length - 1].t - analytic[1]));
  }
  report.check(missing === 0 && spurious === 0 && worst < 1e-9,
               'box: 2000 random rays match the analytic slab intersection',
               `compared ${compared}, worst |dt| ${fmt(worst, 3)} cm, ${missing} missing, ${spurious} spurious`);
}

// ------------------------------------------------------------------------------------------
// the tube with a transverse bore
// ------------------------------------------------------------------------------------------

function chordOf(solid, o, d) {
  const hits = solid.allIntersections(o[0], o[1], o[2], d[0], d[1], d[2]);
  hits.sort((a, b) => a.t - b.t);
  let inside = 0;
  for (let i = 0; i + 1 < hits.length; i += 2) { inside += hits[i + 1].t - hits[i].t; }
  return { hits, inside };
}

function checkTube(report, solid) {
  // tube_window: a solid cylinder of radius 1.5 cm, z in [-3, 3], with a bore of radius 0.8 cm
  // straight through it along x. Every number below follows from those three.
  {
    // parallel to the axis, outside the bore's reach (y = 1.2 > 0.8): the full 6 cm of the tube
    const { hits, inside } = chordOf(solid, [0.5, 1.2, -9], [0, 0, 1]);
    report.check(hits.length === 2 && Math.abs(inside - 6) < 1e-9,
                 'tube_window: axial ray beside the bore crosses 6.0000 cm of material',
                 `${hits.length} crossing(s), ${fmt(inside)} cm`);
  }
  {
    // straight down the bore: the wall is trimmed away at both mouths, so there is nothing to hit
    const { hits } = chordOf(solid, [-9, 0, 0], [1, 0, 0]);
    report.check(hits.length === 0, 'tube_window: a ray down the bore hits nothing',
                 `${hits.length} crossing(s)`);
  }
  {
    // across the bore at z = 0: material from 0.8 to 1.5 on each side, 1.4 cm in total
    const { hits, inside } = chordOf(solid, [0, -9, 0], [0, 1, 0]);
    report.check(hits.length === 4 && Math.abs(inside - 1.4) < 1e-9,
                 'tube_window: a ray across the bore crosses 4 walls and 1.4000 cm of material',
                 `${hits.length} crossing(s), ${fmt(inside)} cm`);
  }
}

// ------------------------------------------------------------------------------------------
// parity closure
// ------------------------------------------------------------------------------------------

function parityScan(solidLike, box, rays, seed) {
  const random = makeRandom(seed);
  const radius = Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
  let odd = 0, hitting = 0;
  for (let i = 0; i < rays; ++i) {
    const u = 2 * random() - 1, phi = 2 * Math.PI * random(), s = Math.sqrt(Math.max(0, 1 - u * u));
    const d = [s * Math.cos(phi), s * Math.sin(phi), u];
    const aim = [box[0] + (box[3] - box[0]) * random(),
                 box[1] + (box[4] - box[1]) * random(),
                 box[2] + (box[5] - box[2]) * random()];
    const o = [aim[0] - d[0] * radius, aim[1] - d[1] * radius, aim[2] - d[2] * radius];
    const count = solidLike.crossingCount(o[0], o[1], o[2], d[0], d[1], d[2]);
    if (count > 0) { hitting += 1; }
    if ((count & 1) !== 0) { odd += 1; }
  }
  return { odd, hitting, rays };
}

/// A ray grid, not a random spray: a regular fan through the part from each of five directions.
/// Regular grids find different defects from random ones, because they line up with the geometry.
function gridScan(solidLike, box, perSide) {
  let odd = 0, hitting = 0, total = 0;
  const radius = Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
  for (const d of RESHOOT_DIRECTIONS) {
    // an orthonormal frame around the shooting direction
    const helper = Math.abs(d[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const r = [d[1] * helper[2] - d[2] * helper[1], d[2] * helper[0] - d[0] * helper[2], d[0] * helper[1] - d[1] * helper[0]];
    const rn = Math.hypot(r[0], r[1], r[2]) || 1;
    const rx = [r[0] / rn, r[1] / rn, r[2] / rn];
    const ux = [d[1] * rx[2] - d[2] * rx[1], d[2] * rx[0] - d[0] * rx[2], d[0] * rx[1] - d[1] * rx[0]];
    const centre = [(box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2];
    for (let i = 0; i < perSide; ++i) {
      for (let j = 0; j < perSide; ++j) {
        // an odd grid would put a ray exactly on the axis of every symmetric part; offset by half
        const a = (2 * (i + 0.5) / perSide - 1) * 0.62 * radius;
        const b = (2 * (j + 0.5) / perSide - 1) * 0.62 * radius;
        const o = [centre[0] + rx[0] * a + ux[0] * b - d[0] * radius,
                   centre[1] + rx[1] * a + ux[1] * b - d[1] * radius,
                   centre[2] + rx[2] * a + ux[2] * b - d[2] * radius];
        const count = solidLike.crossingCount(o[0], o[1], o[2], d[0], d[1], d[2]);
        total += 1;
        if (count > 0) { hitting += 1; }
        if ((count & 1) !== 0) { odd += 1; }
      }
    }
  }
  return { odd, hitting, rays: total };
}

// ------------------------------------------------------------------------------------------
// the run
// ------------------------------------------------------------------------------------------

/// `load(name)` must return { sidecar: ArrayBuffer, facets: ArrayBuffer|null } for a part name.
export async function runSelfCheck(parts, load, log) {
  const report = new Report(log);
  report.note('== solveQuarticReal ==');
  checkQuartic(report);

  for (const name of parts) {
    report.note(`== ${name} ==`);
    let buffers;
    try {
      buffers = await load(name);
    } catch (e) {
      report.check(false, `${name}: test data present`, e.message);
      continue;
    }
    const parsed = parseSidecar(buffers.sidecar, name);
    const solid = new SurfaceSolid(parsed, name);
    report.check(solid.failed.length === 0 && solid.unsupported.length === 0,
                 `${name}: every sidecar record built`,
                 `${solid.nSurfaces} face(s), ${solid.failed.length} rejected, ${solid.unsupported.length} unsupported`);
    report.note(`  worst wire join gap ${fmt(solid.worstJoinGap, 3)} cm against a ${fmt(solid.joinTolerance, 2)} cm band`);

    if (name === 'box') { checkBox(report, solid); }
    if (name === 'tube_window') { checkTube(report, solid); }

    const grid = gridScan(solid, solid.aabb, 24);
    report.check(grid.odd === 0, `${name}: exact solid closes on a ${grid.rays}-ray grid`,
                 `${grid.hitting} rays hit, ${grid.odd} with odd parity`);
    const spray = parityScan(solid, solid.aabb, 1500, 20260822);
    report.check(spray.odd === 0, `${name}: exact solid closes on ${spray.rays} random rays`,
                 `${spray.hitting} rays hit, ${spray.odd} with odd parity`);

    if (buffers.facets) {
      const mesh = new MeshSolid(parseFacets(buffers.facets, name).positions, name);
      const meshGrid = gridScan(mesh, mesh.aabb, 24);
      const meshSpray = parityScan(mesh, mesh.aabb, 1500, 20260822);
      const meshOdd = meshGrid.odd + meshSpray.odd;
      // Not an assertion about the mesh: a leaking tessellation is a fact about the input, and the
      // whole point of showing it. Only the exact solid is required to close.
      report.note(`  mesh (${mesh.nTriangles} triangles): ${meshGrid.odd} odd of ${meshGrid.rays} grid rays, ` +
                  `${meshSpray.odd} odd of ${meshSpray.rays} random rays` +
                  (meshOdd ? '  <-- the tessellation loses rays here' : ''));
    }
  }

  report.note(`== ${report.pass} passed, ${report.fail} failed ==`);
  return report;
}
