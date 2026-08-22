// The exact solid: a set of trimmed analytic surfaces, and ray queries over them.
// Ported from the parity / crossing-cluster machinery of Detectors/Base/src/O2BVHSurfaceSolid.cxx.
//
// A cluster of near-equal hits is one shared edge seen by several faces. A cluster whose members
// all enter, or all exit, is one genuine crossing; a cluster carrying both is a graze -- the ray
// touched an edge and came back out the side it came from, so it crossed nothing. Keeping that
// distinction is what makes the parity count trustworthy, and it is the reason the watertightness
// overlay means something.

import { buildSurface, sameIntersection } from './surfaces.js';
import { wireJoinToleranceFor } from './sidecar.js';
import { K_TOLERANCE } from './curve2d.js';

export const K_RAY_TOLERANCE = 1e-9;
export const BIG = 1e30;

export const Sense = { Entering: 1, Exiting: -1, Tangential: 0 };

/// Which side of the surface a single hit is on: the outward normal opposes the direction when
/// entering; within kTolerance of tangency the surface is not actually crossed.
export function crossingSense(hit, dx, dy, dz) {
  const alignment = hit.nx * dx + hit.ny * dy + hit.nz * dz;
  if (alignment < -K_TOLERANCE) { return Sense.Entering; }
  if (alignment > K_TOLERANCE) { return Sense.Exiting; }
  return Sense.Tangential;
}

/// Sort `hits` and walk their clusters in increasing distance, calling
/// visitor(firstIndex, endIndex, sense); return false from the visitor to stop.
/// Membership is compared against the cluster's FIRST member, never its predecessor: chaining
/// neighbour to neighbour is transitive and lets N hits merge into one cluster N windows wide.
export function forEachCrossingCluster(hits, dx, dy, dz, visitor) {
  hits.sort((a, b) => a.t - b.t);
  let i = 0;
  while (i < hits.length) {
    let entering = false, exiting = false;
    let end = i;
    while (end < hits.length && (end === i || sameIntersection(hits[end].t, hits[i].t))) {
      const s = crossingSense(hits[end], dx, dy, dz);
      if (s === Sense.Entering) { entering = true; } else if (s === Sense.Exiting) { exiting = true; }
      ++end;
    }
    const sense = entering === exiting ? Sense.Tangential : (entering ? Sense.Entering : Sense.Exiting);
    if (!visitor(i, end, sense)) { return; }
    i = end;
  }
}

export class SurfaceSolid {
  /// `parsed` is the object parseSidecar returns.
  constructor(parsed, label = 'solid') {
    this.label = label;
    this.version = parsed.version;
    this.modelTolerance = parsed.modelTolerance;
    this.modelToleranceStated = parsed.modelToleranceStated;
    this.joinTolerance = wireJoinToleranceFor(parsed.modelTolerance);
    this.warnings = parsed.warnings.slice();
    this.surfaces = [];
    this.unsupported = [];      // records this port cannot render, kept and counted
    this.failed = [];           // records that threw while building, kept and counted
    this.counts = {};           // per surface-type-name count of built faces
    this.bsplineTrimFaces = 0;
    this.wireTrimFaces = 0;
    this.worstJoinGap = 0;

    for (const record of parsed.surfaces) {
      try {
        const s = buildSurface(record, this.joinTolerance);
        if (s.reason) {
          this.unsupported.push({ index: record.index, reason: s.reason });
        } else {
          this.counts[s.typeName] = (this.counts[s.typeName] || 0) + 1;
          if (s.hasWireTrim) { this.wireTrimFaces += 1; }
          if (s.hasBSplineTrim) { this.bsplineTrimFaces += 1; }
          this.worstJoinGap = Math.max(this.worstJoinGap, s.worstJoinGap || 0);
          this.surfaces.push(s);
        }
      } catch (e) {
        this.failed.push({ index: record.index, type: record.type, reason: e.message });
      }
    }

    this.aabb = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
    for (const s of this.surfaces) {
      for (let i = 0; i < 3; ++i) {
        this.aabb[i] = Math.min(this.aabb[i], s.aabb[i]);
        this.aabb[i + 3] = Math.max(this.aabb[i + 3], s.aabb[i + 3]);
      }
    }
    if (!Number.isFinite(this.aabb[0])) { this.aabb = [-1, -1, -1, 1, 1, 1]; }
  }

  get nSurfaces() { return this.surfaces.length; }

  /// Every crossing of the trimmed patches along the ray, unsorted.
  allIntersections(ox, oy, oz, dx, dy, dz, tMin = K_RAY_TOLERANCE, tMax = BIG) {
    const hits = [];
    const invx = 1 / dx, invy = 1 / dy, invz = 1 / dz;
    for (const s of this.surfaces) {
      const box = s.aabb;
      // slab test; the boxes carry a kBVHBoxTolerance margin so a grazing hit is never pruned
      let t0 = tMin, t1 = tMax;
      let a = (box[0] - ox) * invx, b = (box[3] - ox) * invx;
      if (a > b) { const tmp = a; a = b; b = tmp; }
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      a = (box[1] - oy) * invy; b = (box[4] - oy) * invy;
      if (a > b) { const tmp = a; a = b; b = tmp; }
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      a = (box[2] - oz) * invz; b = (box[5] - oz) * invz;
      if (a > b) { const tmp = a; a = b; b = tmp; }
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      if (!(t0 <= t1)) { continue; }
      s.appendIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax, hits);
    }
    return hits;
  }

  /// Nearest genuine crossing: { t, nx, ny, nz, sense, surface, boundary } or null.
  firstHit(ox, oy, oz, dx, dy, dz, tMin = K_RAY_TOLERANCE, tMax = BIG) {
    const hits = this.allIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax);
    let result = null;
    forEachCrossingCluster(hits, dx, dy, dz, (first, end, sense) => {
      if (sense === Sense.Tangential) { return true; }
      // Representative of the cluster: the member whose sense is the cluster's.
      let pick = hits[first];
      for (let i = first; i < end; ++i) {
        if (crossingSense(hits[i], dx, dy, dz) === sense) { pick = hits[i]; break; }
      }
      result = { t: pick.t, nx: pick.nx, ny: pick.ny, nz: pick.nz, sense, surface: pick.surface, boundary: pick.boundary };
      return false;
    });
    return result;
  }

  /// Genuine crossings along the whole ray, and whether their count is odd. An odd count means
  /// the ray entered the solid and never came out: a hole in the surface set.
  crossingCount(ox, oy, oz, dx, dy, dz) {
    const hits = this.allIntersections(ox, oy, oz, dx, dy, dz, K_RAY_TOLERANCE, BIG);
    let crossings = 0;
    forEachCrossingCluster(hits, dx, dy, dz, (first, end, sense) => {
      if (sense !== Sense.Tangential) { ++crossings; }
      return true;
    });
    return crossings;
  }

  /// Parity-based containment along one direction (the C++ parityAlong).
  containsAlong(x, y, z, dx, dy, dz) { return (this.crossingCount(x, y, z, dx, dy, dz) & 1) !== 0; }

  /// Containment by majority vote over five quasi-uniform directions, as the kernel's
  /// containsByVote does: one unlucky aim through an edge no longer decides the answer.
  contains(x, y, z) {
    let inside = 0, outside = 0;
    for (const d of RESHOOT_DIRECTIONS) {
      if (this.containsAlong(x, y, z, d[0], d[1], d[2])) { ++inside; } else { ++outside; }
      if (inside >= 3 || outside >= 3) { break; }
    }
    return inside > outside;
  }
}

/// The arbitrary skew test direction the kernel uses for single-shot parity: it probes all normals
/// and avoids the evident symmetries.
export const CONTAINS_TEST_DIRECTION = (() => {
  const v = [1, 1.41421356237, 1.73205080757];
  const n = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / n, v[1] / n, v[2] / n];
})();

/// Five directions on a golden-angle spiral: quasi-uniform on the sphere, so none aligns with a
/// coordinate axis or a 45-degree symmetry plane and no two are correlated.
export const RESHOOT_DIRECTIONS = (() => {
  const golden = Math.PI * (3 - Math.sqrt(5));
  const out = [];
  for (let i = 0; i < 5; ++i) {
    const cosTheta = 1 - 2 * (i + 0.5) / 5;
    const sinTheta = Math.sqrt(Math.max(0, 1 - cosTheta * cosTheta));
    const phi = golden * i;
    out.push([sinTheta * Math.cos(phi), sinTheta * Math.sin(phi), cosTheta]);
  }
  return out;
})();
