// The six analytic bounded-surface kernels, built from parsed sidecar records.
// Ported from Detectors/Base/src/BoundedSurface.h (PlanarBoundedSurface, CurvedPlanarBoundedSurface,
// CylindricalBoundedSurface, SphericalBoundedSurface, ConicalBoundedSurface, TorusBoundedSurface).
//
// Each surface answers appendIntersections(origin, direction, tMin, tMax, hits), pushing every
// crossing of the *trimmed* patch with its outward-oriented analytic normal. Parity depends on
// this being every crossing, not the nearest one, so tangential grazes are suppressed rather than
// reported once.

import { Curve2D, v2, K_TOLERANCE, K_TOLERANCE_SQ, TWO_PI } from './curve2d.js';
import {
  CurveWire, PolygonWire, ParametricMetric, trimContains, trimWindow, unwrapAngleInto,
  planeParametricMetric, cylinderParametricMetric, coneParametricMetric, sphereParametricMetric,
  torusParametricMetric,
} from './wire.js';
import { solveQuarticReal } from './quartic.js';
import { SURFACE_TYPE, CURVE_TYPE, parseBSplineEdge, surfaceTypeName } from './sidecar.js';

export const K_INTERSECTION_TOLERANCE = 1e-7;
export const K_BVH_BOX_TOLERANCE = 1e-3;

export function sameIntersection(a, b) {
  return Math.abs(a - b) <= K_INTERSECTION_TOLERANCE * Math.max(1, Math.max(Math.abs(a), Math.abs(b)));
}

export function angularTolerance(radius) { return K_TOLERANCE / Math.max(radius, K_TOLERANCE); }

export function angleInSweepRange(angle, start, sweep, tolerance) {
  if (sweep >= TWO_PI - K_TOLERANCE) { return true; }
  let delta = angle - start;
  delta -= TWO_PI * Math.floor(delta / TWO_PI);
  return delta <= sweep + tolerance || delta >= TWO_PI - tolerance;
}

const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm3 = (a) => Math.sqrt(dot3(a, a));
const scale3 = (a, s) => [a[0] * s, a[1] * s, a[2] * s];

/// Orthonormal frame (U, V, W) with W along `axis` and U the projection of `referenceAxisU`
/// perpendicular to W, so that U x V = W (CylindricalBoundedSurface::makeFrame).
export function makeFrame(axis, referenceAxisU) {
  const wLen = norm3(axis);
  if (wLen <= K_TOLERANCE) { throw new Error('surface axis is degenerate'); }
  const W = scale3(axis, 1 / wLen);
  const proj = [referenceAxisU[0] - W[0] * dot3(referenceAxisU, W),
                referenceAxisU[1] - W[1] * dot3(referenceAxisU, W),
                referenceAxisU[2] - W[2] * dot3(referenceAxisU, W)];
  const uLen = norm3(proj);
  if (uLen <= K_TOLERANCE) { throw new Error('surface reference axis is parallel to the main axis'); }
  const U = scale3(proj, 1 / uLen);
  const V = cross3(W, U);
  return { U, V, W };
}

/// Range of a cos(t) + b sin(t) over [t0, t1]: the endpoint values widened to +/- the amplitude
/// wherever the crest or trough falls inside the interval. Exact, so a cover box carries no slack.
export function sinusoidRange(a, b, t0, t1) {
  const f = (t) => a * Math.cos(t) + b * Math.sin(t);
  let lo = Math.min(f(t0), f(t1)), hi = Math.max(f(t0), f(t1));
  const amplitude = Math.hypot(a, b);
  if (amplitude <= 0) { return [lo, hi]; }
  const crest = Math.atan2(b, a);
  for (let k = -2; k <= 2; ++k) {
    const tCrest = crest + k * TWO_PI;
    if (tCrest > t0 && tCrest < t1) { hi = amplitude; }
    const tTrough = crest + Math.PI + k * TWO_PI;
    if (tTrough > t0 && tTrough < t1) { lo = -amplitude; }
  }
  return [lo, hi];
}

// ------------------------------------------------------------------------------------------
// Trim construction from a sidecar record's wire block
// ------------------------------------------------------------------------------------------

function edgeToCurve(edge) {
  const p = edge.params;
  if (edge.curveType === CURVE_TYPE.LINE) {
    if (p.length < 4) { throw new Error('malformed line edge'); }
    return { curve: Curve2D.makeLine(v2(p[0], p[1]), v2(p[2], p[3])), curved: false };
  }
  if (edge.curveType === CURVE_TYPE.ARC) {
    if (p.length < 5) { throw new Error('malformed arc edge'); }
    return { curve: Curve2D.makeArc(v2(p[0], p[1]), p[2], p[3], p[3] + p[4]), curved: true };
  }
  if (edge.curveType === CURVE_TYPE.BSPLINE) {
    const b = parseBSplineEdge(p);
    if (!b) { throw new Error('malformed bspline edge'); }
    return { curve: Curve2D.makeBSpline(b.degree, b.poles, b.weights, b.knots), curved: true };
  }
  throw new Error(`unsupported wire edge curve type ${edge.curveType}`);
}

/// Convert a record's wire block into { outer, inners, anyCurved }, with each loop a list of
/// Curve2D. Mirrors wireToCurves + collectQuadricTrim, minus the join-gap rejection: a sidecar
/// that the kernel accepts has already passed that test, and a viewer that refuses to draw a face
/// teaches nothing. Join gaps are measured and reported instead (see buildSurface).
function collectTrim(record, metric, joinTolerance) {
  let outer = null;
  const inners = [];
  let anyCurved = false;
  let worstJoinGap = 0;
  for (const wire of record.wires) {
    const curves = [];
    for (const edge of wire.edges) {
      const { curve, curved } = edgeToCurve(edge);
      anyCurved = anyCurved || curved;
      curves.push(curve);
    }
    for (let e = 0; e < curves.length; ++e) {
      const end = curves[e].endPoint();
      const nextStart = curves[(e + 1) % curves.length].startPoint();
      worstJoinGap = Math.max(worstJoinGap, Math.sqrt(metric.distanceSq(end, nextStart)));
    }
    if (wire.role === 0) {
      if (outer) { throw new Error('surface has more than one outer trim wire'); }
      outer = curves;
    } else {
      inners.push(curves);
    }
  }
  if (!outer) { throw new Error('trim block has no outer wire'); }
  return { outer, inners, anyCurved, worstJoinGap, joinTolerance };
}

// ------------------------------------------------------------------------------------------
// Surfaces
// ------------------------------------------------------------------------------------------

class BaseSurface {
  constructor(index, type, innerWall) {
    this.index = index;
    this.type = type;
    this.typeName = surfaceTypeName(type);
    this.innerWall = innerWall;
    this.aabb = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
    this.hasWireTrim = false;
    this.hasBSplineTrim = false;
    this.worstJoinGap = 0;
  }
  _expandAabb(margin = K_BVH_BOX_TOLERANCE) {
    for (let i = 0; i < 3; ++i) { this.aabb[i] -= margin; this.aabb[i + 3] += margin; }
  }

  /// The face's trimmed boundary as 3D polylines, for the exact-edge overlay drawn on top of the
  /// tessellation. Curved trim edges are sampled here purely for display; the intersection kernels
  /// never see these points.
  patchOutline(segmentsPerTurn = 96) {
    const out = [];
    if (this.outerWire) {
      for (const wire of [this.outerWire, ...(this.innerWires || [])]) {
        const loop = wire.displayPolyline(segmentsPerTurn);
        if (loop.length > 1) { out.push(loop.map(p => this.paramToGlobal(p.u, p.v))); }
      }
      return out;
    }
    const w = this.paramWindow ? this.paramWindow() : null;
    if (!w) { return out; }
    // The parametric rectangle: the two u = const edges stay single segments, the two v = const
    // edges are chorded, because in 3D they are curves.
    const steps = Math.max(2, Math.round(segmentsPerTurn * (w.uMax - w.uMin) / TWO_PI));
    const loop = [];
    for (let i = 0; i <= steps; ++i) { loop.push(this.paramToGlobal(w.uMin + (w.uMax - w.uMin) * i / steps, w.vMin)); }
    for (let i = steps; i >= 0; --i) { loop.push(this.paramToGlobal(w.uMin + (w.uMax - w.uMin) * i / steps, w.vMax)); }
    loop.push(loop[0]);
    out.push(loop);
    return out;
  }
}

export class PlaneSurface extends BaseSurface {
  constructor(record, joinTolerance) {
    super(record.index, record.type, record.innerWall);
    const p = record.params;
    this.origin = [p[0], p[1], p[2]];
    this.axisU = [p[3], p[4], p[5]];
    this.axisV = [p[6], p[7], p[8]];
    const n = cross3(this.axisU, this.axisV);
    this.areaScale = norm3(n);
    if (this.areaScale <= K_TOLERANCE) { throw new Error('surface frame axes are degenerate'); }
    this.normal = scale3(n, 1 / this.areaScale);

    const g = planeParametricMetric(this.axisU, this.axisV);
    this.gUU = g[0]; this.gUV = g[1]; this.gVV = g[2];
    const det = this.gUU * this.gVV - this.gUV * this.gUV;
    if (Math.abs(det) <= K_TOLERANCE_SQ) { throw new Error('surface frame metric is singular'); }
    this.invDet = 1 / det;
    this.metric = new ParametricMetric(() => g);

    const trim = collectTrim(record, this.metric, joinTolerance);
    this.worstJoinGap = trim.worstJoinGap;
    this.hasWireTrim = true;
    // A pure line-segment loop keeps the polygon path; any arc or bspline edge routes to the
    // curved path -- the same split the loader makes between AddPlanarSurface and
    // AddCurvedPlanarSurface.
    if (trim.anyCurved) {
      this.outerWire = new CurveWire(trim.outer, 0, this.metric);
      this.innerWires = trim.inners.map(l => new CurveWire(l, 1, this.metric));
      this.hasBSplineTrim = trim.outer.concat(...trim.inners).some(c => c.isBSpline);
    } else {
      this.outerWire = new PolygonWire(trim.outer.map(c => c.lineStart), 0, this.metric);
      this.innerWires = trim.inners.map(l => new PolygonWire(l.map(c => c.lineStart), 1, this.metric));
    }

    const lower = v2(Infinity, Infinity), upper = v2(-Infinity, -Infinity);
    this.outerWire.parametricBounds(lower, upper);
    // The affine image of the parametric AABB contains the patch; its corners bound the 3D AABB.
    for (const cu of [lower.u, upper.u]) {
      for (const cv of [lower.v, upper.v]) {
        const g3 = this.toGlobal(cu, cv);
        for (let i = 0; i < 3; ++i) {
          this.aabb[i] = Math.min(this.aabb[i], g3[i]);
          this.aabb[i + 3] = Math.max(this.aabb[i + 3], g3[i]);
        }
      }
    }
    this._expandAabb();
  }

  toGlobal(u, v) {
    return [this.origin[0] + this.axisU[0] * u + this.axisV[0] * v,
            this.origin[1] + this.axisU[1] * u + this.axisV[1] * v,
            this.origin[2] + this.axisU[2] * u + this.axisV[2] * v];
  }
  toLocal(x, y, z) {
    const rx = x - this.origin[0], ry = y - this.origin[1], rz = z - this.origin[2];
    const pu = rx * this.axisU[0] + ry * this.axisU[1] + rz * this.axisU[2];
    const pv = rx * this.axisV[0] + ry * this.axisV[1] + rz * this.axisV[2];
    return v2((pu * this.gVV - pv * this.gUV) * this.invDet, (pv * this.gUU - pu * this.gUV) * this.invDet);
  }

  containsLocal(uv) { return trimContains(this.outerWire, this.innerWires, uv); }

  paramToGlobal(u, v) { return this.toGlobal(u, v); }

  appendIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax, hits) {
    const n = this.normal;
    const denominator = n[0] * dx + n[1] * dy + n[2] * dz;
    if (Math.abs(denominator) <= K_TOLERANCE) { return; }
    const t = ((this.origin[0] - ox) * n[0] + (this.origin[1] - oy) * n[1] + (this.origin[2] - oz) * n[2]) / denominator;
    if (t < tMin || t > tMax) { return; }
    const uv = this.toLocal(ox + t * dx, oy + t * dy, oz + t * dz);
    const c = this.containsLocal(uv);
    if (!c.contained) { return; }
    hits.push({ t, nx: n[0], ny: n[1], nz: n[2], boundary: c.boundary, surface: this.index });
  }

  normalAt() { return this.normal; }
}

export class CylinderSurface extends BaseSurface {
  constructor(record, joinTolerance) {
    super(record.index, record.type, record.innerWall);
    const p = record.params;
    this.center = [p[0], p[1], p[2]];
    const frame = makeFrame([p[3], p[4], p[5]], [p[6], p[7], p[8]]);
    this.U = frame.U; this.V = frame.V; this.W = frame.W;
    this.radius = p[9];
    this.heightMin = p[10]; this.heightMax = p[11];
    this.phiStart = p[12]; this.phiSweep = Math.min(p[13], TWO_PI);
    if (this.radius <= K_TOLERANCE) { throw new Error('cylindrical surface needs a positive radius'); }
    this.normalSign = this.innerWall ? -1 : 1;
    const g = cylinderParametricMetric(this.radius);
    this.metric = new ParametricMetric(() => g);

    if (record.wires.length) {
      const trim = collectTrim(record, this.metric, joinTolerance);
      this.worstJoinGap = trim.worstJoinGap;
      this.outerWire = new CurveWire(trim.outer, 0, this.metric);
      this.innerWires = trim.inners.map(l => new CurveWire(l, 1, this.metric));
      const w = trimWindow(this.outerWire);
      if (!w) { throw new Error('quadric trim wire spans more than a full turn in phi'); }
      this.phiStart = w.lower.u;
      this.phiSweep = Math.min(TWO_PI, w.upper.u - w.lower.u);
      this.heightMin = w.lower.v; this.heightMax = w.upper.v;
      this.hasWireTrim = true;
      this.hasBSplineTrim = trim.outer.concat(...trim.inners).some(c => c.isBSpline);
    }
    this._bounds();
  }

  _bounds() {
    const [lo, hi] = [this.phiStart, this.phiStart + this.phiSweep];
    for (let i = 0; i < 3; ++i) {
      const [sLo, sHi] = sinusoidRange(this.U[i], this.V[i], lo, hi);
      const hLo = Math.min(this.W[i] * this.heightMin, this.W[i] * this.heightMax);
      const hHi = Math.max(this.W[i] * this.heightMin, this.W[i] * this.heightMax);
      this.aabb[i] = this.center[i] + this.radius * sLo + hLo;
      this.aabb[i + 3] = this.center[i] + this.radius * sHi + hHi;
    }
    this._expandAabb();
  }

  pointAt(phi, h) {
    const c = Math.cos(phi), s = Math.sin(phi);
    return [this.center[0] + this.W[0] * h + (this.U[0] * c + this.V[0] * s) * this.radius,
            this.center[1] + this.W[1] * h + (this.U[1] * c + this.V[1] * s) * this.radius,
            this.center[2] + this.W[2] * h + (this.U[2] * c + this.V[2] * s) * this.radius];
  }

  pointInTrim(phi, height) {
    const u = unwrapAngleInto(phi, this.phiStart, this.phiStart + this.phiSweep);
    return trimContains(this.outerWire, this.innerWires, v2(u, height));
  }

  paramToGlobal(phi, h) { return this.pointAt(phi, h); }
  paramWindow() { return { uMin: this.phiStart, uMax: this.phiStart + this.phiSweep, vMin: this.heightMin, vMax: this.heightMax }; }

  appendIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax, hits) {
    const rx = ox - this.center[0], ry = oy - this.center[1], rz = oz - this.center[2];
    const lox = rx * this.U[0] + ry * this.U[1] + rz * this.U[2];
    const loy = rx * this.V[0] + ry * this.V[1] + rz * this.V[2];
    const loz = rx * this.W[0] + ry * this.W[1] + rz * this.W[2];
    const ldx = dx * this.U[0] + dy * this.U[1] + dz * this.U[2];
    const ldy = dx * this.V[0] + dy * this.V[1] + dz * this.V[2];
    const ldz = dx * this.W[0] + dy * this.W[1] + dz * this.W[2];

    const a = ldx * ldx + ldy * ldy;
    if (a <= K_TOLERANCE_SQ) { return; } // parallel to the axis: no transversal crossing
    const b = 2 * (lox * ldx + loy * ldy);
    const c = lox * lox + loy * loy - this.radius * this.radius;
    const disc = b * b - 4 * a * c;
    if (disc <= 0) { return; }
    const sq = Math.sqrt(disc);
    const t0 = (-b - sq) / (2 * a), t1 = (-b + sq) / (2 * a);
    if (sameIntersection(t0, t1)) { return; } // tangential graze: parity stays even

    for (const t of [t0, t1]) {
      if (t < tMin || t > tMax) { continue; }
      const hu = lox + t * ldx, hv = loy + t * ldy, hh = loz + t * ldz;
      const phi = Math.atan2(hv, hu);
      let boundary = false;
      if (this.hasWireTrim) {
        const r = this.pointInTrim(phi, hh);
        if (!r.contained) { continue; }
        boundary = r.boundary;
      } else {
        if (hh < this.heightMin - K_TOLERANCE || hh > this.heightMax + K_TOLERANCE) { continue; }
        if (!angleInSweepRange(phi, this.phiStart, this.phiSweep, angularTolerance(this.radius))) { continue; }
      }
      const rho = Math.hypot(hu, hv);
      const s = this.normalSign;
      hits.push({
        t,
        nx: (this.U[0] * (hu / rho) + this.V[0] * (hv / rho)) * s,
        ny: (this.U[1] * (hu / rho) + this.V[1] * (hv / rho)) * s,
        nz: (this.U[2] * (hu / rho) + this.V[2] * (hv / rho)) * s,
        boundary, surface: this.index,
      });
    }
  }
}

export class ConeSurface extends BaseSurface {
  constructor(record, joinTolerance) {
    super(record.index, record.type, record.innerWall);
    const p = record.params;
    this.center = [p[0], p[1], p[2]];
    const frame = makeFrame([p[3], p[4], p[5]], [p[6], p[7], p[8]]);
    this.U = frame.U; this.V = frame.V; this.W = frame.W;
    const radiusAtMin = p[9], radiusAtMax = p[10];
    this.heightMin = p[11]; this.heightMax = p[12];
    this.phiStart = p[13]; this.phiSweep = Math.min(p[14], TWO_PI);
    if (this.heightMax - this.heightMin <= K_TOLERANCE) { throw new Error('conical surface needs a positive height range'); }
    this.slope = (radiusAtMax - radiusAtMin) / (this.heightMax - this.heightMin);
    this.radius0 = radiusAtMin - this.slope * this.heightMin;
    this.normalSign = this.innerWall ? -1 : 1;
    this.metric = new ParametricMetric((u, v) => coneParametricMetric(this.radiusAt(v), this.slope));

    if (record.wires.length) {
      const trim = collectTrim(record, this.metric, joinTolerance);
      this.worstJoinGap = trim.worstJoinGap;
      this.outerWire = new CurveWire(trim.outer, 0, this.metric);
      this.innerWires = trim.inners.map(l => new CurveWire(l, 1, this.metric));
      const w = trimWindow(this.outerWire);
      if (!w) { throw new Error('quadric trim wire spans more than a full turn in phi'); }
      this.phiStart = w.lower.u;
      this.phiSweep = Math.min(TWO_PI, w.upper.u - w.lower.u);
      this.heightMin = w.lower.v; this.heightMax = w.upper.v;
      this.hasWireTrim = true;
      this.hasBSplineTrim = trim.outer.concat(...trim.inners).some(c => c.isBSpline);
    }
    this._bounds();
  }

  radiusAt(h) { return this.radius0 + this.slope * h; }
  meanRadius() { return 0.5 * (this.radiusAt(this.heightMin) + this.radiusAt(this.heightMax)); }

  _bounds() {
    const [lo, hi] = [this.phiStart, this.phiStart + this.phiSweep];
    const rA = Math.max(0, this.radiusAt(this.heightMin)), rB = Math.max(0, this.radiusAt(this.heightMax));
    const rLo = Math.min(rA, rB), rHi = Math.max(rA, rB);
    for (let i = 0; i < 3; ++i) {
      const [sLo, sHi] = sinusoidRange(this.U[i], this.V[i], lo, hi);
      const products = [rLo * sLo, rLo * sHi, rHi * sLo, rHi * sHi];
      const hLo = Math.min(this.W[i] * this.heightMin, this.W[i] * this.heightMax);
      const hHi = Math.max(this.W[i] * this.heightMin, this.W[i] * this.heightMax);
      this.aabb[i] = this.center[i] + Math.min(...products) + hLo;
      this.aabb[i + 3] = this.center[i] + Math.max(...products) + hHi;
    }
    this._expandAabb();
  }

  pointInTrim(phi, height) {
    const u = unwrapAngleInto(phi, this.phiStart, this.phiStart + this.phiSweep);
    return trimContains(this.outerWire, this.innerWires, v2(u, height));
  }

  paramToGlobal(phi, h) {
    const c = Math.cos(phi), s = Math.sin(phi), r = this.radiusAt(h);
    return [this.center[0] + this.W[0] * h + (this.U[0] * c + this.V[0] * s) * r,
            this.center[1] + this.W[1] * h + (this.U[1] * c + this.V[1] * s) * r,
            this.center[2] + this.W[2] * h + (this.U[2] * c + this.V[2] * s) * r];
  }
  paramWindow() { return { uMin: this.phiStart, uMax: this.phiStart + this.phiSweep, vMin: this.heightMin, vMax: this.heightMax }; }

  appendIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax, hits) {
    const rx = ox - this.center[0], ry = oy - this.center[1], rz = oz - this.center[2];
    const lox = rx * this.U[0] + ry * this.U[1] + rz * this.U[2];
    const loy = rx * this.V[0] + ry * this.V[1] + rz * this.V[2];
    const loz = rx * this.W[0] + ry * this.W[1] + rz * this.W[2];
    const ldx = dx * this.U[0] + dy * this.U[1] + dz * this.U[2];
    const ldy = dx * this.V[0] + dy * this.V[1] + dz * this.V[2];
    const ldz = dx * this.W[0] + dy * this.W[1] + dz * this.W[2];

    const rAtOrigin = this.radius0 + this.slope * loz;
    const a = ldx * ldx + ldy * ldy - this.slope * this.slope * ldz * ldz;
    const b = 2 * (lox * ldx + loy * ldy - this.slope * ldz * rAtOrigin);
    const c = lox * lox + loy * loy - rAtOrigin * rAtOrigin;

    const candidates = [];
    if (Math.abs(a) <= K_TOLERANCE_SQ) {
      if (Math.abs(b) <= K_TOLERANCE_SQ) { return; } // along the cone surface or its asymptote
      candidates.push(-c / b);
    } else {
      const disc = b * b - 4 * a * c;
      if (disc <= 0) { return; }
      const sq = Math.sqrt(disc);
      const t0 = (-b - sq) / (2 * a), t1 = (-b + sq) / (2 * a);
      if (sameIntersection(t0, t1)) { return; } // tangential graze (also rays through the apex)
      candidates.push(Math.min(t0, t1), Math.max(t0, t1));
    }

    for (const t of candidates) {
      if (t < tMin || t > tMax) { continue; }
      const hh = loz + t * ldz;
      if (this.radiusAt(hh) < -K_TOLERANCE) { continue; } // mirror nappe of the infinite cone
      const hu = lox + t * ldx, hv = loy + t * ldy;
      const rho = Math.hypot(hu, hv);
      if (rho <= K_TOLERANCE) { continue; } // apex hit: the normal is undefined there
      const phi = Math.atan2(hv, hu);
      let boundary = false;
      if (this.hasWireTrim) {
        const r = this.pointInTrim(phi, hh);
        if (!r.contained) { continue; }
        boundary = r.boundary;
      } else {
        if (hh < this.heightMin - K_TOLERANCE || hh > this.heightMax + K_TOLERANCE) { continue; }
        if (!angleInSweepRange(phi, this.phiStart, this.phiSweep, angularTolerance(this.meanRadius()))) { continue; }
      }
      const s = this.normalSign / Math.sqrt(1 + this.slope * this.slope);
      hits.push({
        t,
        nx: (this.U[0] * (hu / rho) + this.V[0] * (hv / rho) - this.W[0] * this.slope) * s,
        ny: (this.U[1] * (hu / rho) + this.V[1] * (hv / rho) - this.W[1] * this.slope) * s,
        nz: (this.U[2] * (hu / rho) + this.V[2] * (hv / rho) - this.W[2] * this.slope) * s,
        boundary, surface: this.index,
      });
    }
  }
}

export class SphereSurface extends BaseSurface {
  constructor(record, joinTolerance) {
    super(record.index, record.type, record.innerWall);
    const p = record.params;
    this.center = [p[0], p[1], p[2]];
    const frame = makeFrame([p[3], p[4], p[5]], [p[6], p[7], p[8]]);
    this.U = frame.U; this.V = frame.V; this.W = frame.W;
    this.radius = p[9];
    this.thetaMin = p[10]; this.thetaMax = p[11];
    this.phiStart = p[12]; this.phiSweep = Math.min(p[13], TWO_PI);
    if (this.radius <= K_TOLERANCE) { throw new Error('spherical surface needs a positive radius'); }
    this.normalSign = this.innerWall ? -1 : 1;
    this.metric = new ParametricMetric((u, v) => sphereParametricMetric(this.radius, v));

    if (record.wires.length) {
      const trim = collectTrim(record, this.metric, joinTolerance);
      this.worstJoinGap = trim.worstJoinGap;
      this.outerWire = new CurveWire(trim.outer, 0, this.metric);
      this.innerWires = trim.inners.map(l => new CurveWire(l, 1, this.metric));
      const w = trimWindow(this.outerWire);
      if (!w) { throw new Error('quadric trim wire spans more than a full turn in phi'); }
      this.phiStart = w.lower.u;
      this.phiSweep = Math.min(TWO_PI, w.upper.u - w.lower.u);
      this.thetaMin = Math.max(0, w.lower.v); this.thetaMax = Math.min(Math.PI, w.upper.v);
      this.hasWireTrim = true;
      this.hasBSplineTrim = trim.outer.concat(...trim.inners).some(c => c.isBSpline);
    }
    // Conservative: the full sphere's box. A trimmed cap is a subset of it.
    for (let i = 0; i < 3; ++i) {
      this.aabb[i] = this.center[i] - this.radius;
      this.aabb[i + 3] = this.center[i] + this.radius;
    }
    this._expandAabb();
  }

  pointInTrim(phi, theta) {
    const u = unwrapAngleInto(phi, this.phiStart, this.phiStart + this.phiSweep);
    return trimContains(this.outerWire, this.innerWires, v2(u, theta));
  }

  paramToGlobal(phi, theta) {
    const st = Math.sin(theta), ct = Math.cos(theta), cp = Math.cos(phi), sp = Math.sin(phi);
    const R = this.radius;
    return [this.center[0] + (this.U[0] * st * cp + this.V[0] * st * sp + this.W[0] * ct) * R,
            this.center[1] + (this.U[1] * st * cp + this.V[1] * st * sp + this.W[1] * ct) * R,
            this.center[2] + (this.U[2] * st * cp + this.V[2] * st * sp + this.W[2] * ct) * R];
  }
  paramWindow() { return { uMin: this.phiStart, uMax: this.phiStart + this.phiSweep, vMin: this.thetaMin, vMax: this.thetaMax }; }

  /// directionInTrim, on the local point.
  _directionInTrim(lx, ly, lz) {
    const r = Math.hypot(lx, ly, lz);
    if (r <= K_TOLERANCE) { return { contained: true, boundary: false }; }
    const thetaTol = angularTolerance(this.radius);
    const theta = Math.acos(Math.max(-1, Math.min(1, lz / r)));
    const transverse = Math.hypot(lx, ly);
    if (this.hasWireTrim) {
      if (transverse <= K_TOLERANCE) {
        const ok = theta >= this.thetaMin - thetaTol && theta <= this.thetaMax + thetaTol;
        return { contained: ok, boundary: false };
      }
      return this.pointInTrim(Math.atan2(ly, lx), theta);
    }
    if (theta < this.thetaMin - thetaTol || theta > this.thetaMax + thetaTol) { return { contained: false, boundary: false }; }
    if (transverse <= K_TOLERANCE) { return { contained: true, boundary: false }; }
    return { contained: angleInSweepRange(Math.atan2(ly, lx), this.phiStart, this.phiSweep, thetaTol), boundary: false };
  }

  appendIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax, hits) {
    const rx = ox - this.center[0], ry = oy - this.center[1], rz = oz - this.center[2];
    const a = dx * dx + dy * dy + dz * dz;
    if (a <= K_TOLERANCE_SQ) { return; }
    const b = 2 * (rx * dx + ry * dy + rz * dz);
    const c = rx * rx + ry * ry + rz * rz - this.radius * this.radius;
    const disc = b * b - 4 * a * c;
    if (disc <= 0) { return; }
    const sq = Math.sqrt(disc);
    const t0 = (-b - sq) / (2 * a), t1 = (-b + sq) / (2 * a);
    if (sameIntersection(t0, t1)) { return; }

    for (const t of [t0, t1]) {
      if (t < tMin || t > tMax) { continue; }
      const gx = ox + t * dx - this.center[0], gy = oy + t * dy - this.center[1], gz = oz + t * dz - this.center[2];
      const lx = gx * this.U[0] + gy * this.U[1] + gz * this.U[2];
      const ly = gx * this.V[0] + gy * this.V[1] + gz * this.V[2];
      const lz = gx * this.W[0] + gy * this.W[1] + gz * this.W[2];
      const r = this._directionInTrim(lx, ly, lz);
      if (!r.contained) { continue; }
      const s = this.normalSign / this.radius;
      hits.push({
        t,
        nx: (this.U[0] * lx + this.V[0] * ly + this.W[0] * lz) * s,
        ny: (this.U[1] * lx + this.V[1] * ly + this.W[1] * lz) * s,
        nz: (this.U[2] * lx + this.V[2] * ly + this.W[2] * lz) * s,
        boundary: r.boundary, surface: this.index,
      });
    }
  }
}

export class TorusSurface extends BaseSurface {
  constructor(record, joinTolerance) {
    super(record.index, record.type, record.innerWall);
    const p = record.params;
    this.center = [p[0], p[1], p[2]];
    const frame = makeFrame([p[3], p[4], p[5]], [p[6], p[7], p[8]]);
    this.U = frame.U; this.V = frame.V; this.W = frame.W;
    this.majorRadius = p[9]; this.minorRadius = p[10];
    this.phiStart = p[11]; this.phiSweep = Math.min(p[12], TWO_PI);
    this.tubeStart = p[13]; this.tubeSweep = Math.min(p[14], TWO_PI);
    if (this.majorRadius <= K_TOLERANCE || this.minorRadius <= K_TOLERANCE) {
      throw new Error('toroidal surface needs positive major and minor radii');
    }
    this.normalSign = this.innerWall ? -1 : 1;
    this.metric = new ParametricMetric((u, v) => torusParametricMetric(this.majorRadius, this.minorRadius, v));

    if (record.wires.length) {
      const trim = collectTrim(record, this.metric, joinTolerance);
      this.worstJoinGap = trim.worstJoinGap;
      this.outerWire = new CurveWire(trim.outer, 0, this.metric);
      this.innerWires = trim.inners.map(l => new CurveWire(l, 1, this.metric));
      const w = trimWindow(this.outerWire);
      if (!w) { throw new Error('toroidal trim wire spans more than a full turn in the ring angle'); }
      this.phiStart = w.lower.u;
      this.phiSweep = Math.min(TWO_PI, w.upper.u - w.lower.u);
      this.tubeStart = w.lower.v;
      this.tubeSweep = Math.min(TWO_PI, w.upper.v - w.lower.v);
      this.hasWireTrim = true;
      this.hasBSplineTrim = trim.outer.concat(...trim.inners).some(c => c.isBSpline);
    }
    // Conservative: the full torus's box.
    const R = this.majorRadius + this.minorRadius, r = this.minorRadius;
    for (let i = 0; i < 3; ++i) {
      const transverse = Math.hypot(this.U[i], this.V[i]);
      const extent = R * transverse + r * Math.abs(this.W[i]);
      this.aabb[i] = this.center[i] - extent;
      this.aabb[i + 3] = this.center[i] + extent;
    }
    this._expandAabb();
  }

  pointInTrim(phiRing, phiTube) {
    const u = unwrapAngleInto(phiRing, this.phiStart, this.phiStart + this.phiSweep);
    const v = unwrapAngleInto(phiTube, this.tubeStart, this.tubeStart + this.tubeSweep);
    return trimContains(this.outerWire, this.innerWires, v2(u, v));
  }

  paramToGlobal(phiRing, phiTube) {
    const ringRadius = this.majorRadius + this.minorRadius * Math.cos(phiTube);
    const c = Math.cos(phiRing), s = Math.sin(phiRing), z = this.minorRadius * Math.sin(phiTube);
    return [this.center[0] + (this.U[0] * c + this.V[0] * s) * ringRadius + this.W[0] * z,
            this.center[1] + (this.U[1] * c + this.V[1] * s) * ringRadius + this.W[1] * z,
            this.center[2] + (this.U[2] * c + this.V[2] * s) * ringRadius + this.W[2] * z];
  }
  paramWindow() { return { uMin: this.phiStart, uMax: this.phiStart + this.phiSweep, vMin: this.tubeStart, vMax: this.tubeStart + this.tubeSweep }; }

  localNormal(lx, ly, lz) {
    const rho = Math.hypot(lx, ly);
    if (rho <= K_TOLERANCE) {
      const s = lz >= 0 ? this.normalSign : -this.normalSign;
      return [this.W[0] * s, this.W[1] * s, this.W[2] * s];
    }
    const radialFactor = (rho - this.majorRadius) / rho;
    const nx = radialFactor * lx, ny = radialFactor * ly, nz = lz;
    const len = Math.hypot(nx, ny, nz);
    if (len <= K_TOLERANCE) { return scale3(this.U, this.normalSign); }
    const s = this.normalSign / len;
    return [(this.U[0] * nx + this.V[0] * ny + this.W[0] * nz) * s,
            (this.U[1] * nx + this.V[1] * ny + this.W[1] * nz) * s,
            (this.U[2] * nx + this.V[2] * ny + this.W[2] * nz) * s];
  }

  appendIntersections(ox, oy, oz, dx, dy, dz, tMin, tMax, hits) {
    const rx = ox - this.center[0], ry = oy - this.center[1], rz = oz - this.center[2];
    const lox = rx * this.U[0] + ry * this.U[1] + rz * this.U[2];
    const loy = rx * this.V[0] + ry * this.V[1] + rz * this.V[2];
    const loz = rx * this.W[0] + ry * this.W[1] + rz * this.W[2];
    const ldx = dx * this.U[0] + dy * this.U[1] + dz * this.U[2];
    const ldy = dx * this.V[0] + dy * this.V[1] + dz * this.V[2];
    const ldz = dx * this.W[0] + dy * this.W[1] + dz * this.W[2];

    // (|X|^2 + R^2 - r^2)^2 = 4 R^2 (x^2 + y^2), substituting X = O + t D.
    const dirDotDir = ldx * ldx + ldy * ldy + ldz * ldz;
    if (dirDotDir <= K_TOLERANCE_SQ) { return; }
    const originDotDir = lox * ldx + loy * ldy + loz * ldz;
    const originDotOrigin = lox * lox + loy * loy + loz * loz;
    const constantK = this.majorRadius * this.majorRadius - this.minorRadius * this.minorRadius;
    const transverseE = ldx * ldx + ldy * ldy;
    const transverseF = lox * ldx + loy * ldy;
    const transverseG = lox * lox + loy * loy;
    const fourRSquared = 4 * this.majorRadius * this.majorRadius;

    const coeff4 = dirDotDir * dirDotDir;
    const coeff3 = 4 * dirDotDir * originDotDir;
    const coeff2 = 4 * originDotDir * originDotDir + 2 * dirDotDir * (originDotOrigin + constantK) - fourRSquared * transverseE;
    const coeff1 = 4 * originDotDir * (originDotOrigin + constantK) - 2 * fourRSquared * transverseF;
    const coeff0 = (originDotOrigin + constantK) * (originDotOrigin + constantK) - fourRSquared * transverseG;

    const candidates = solveQuarticReal(coeff4, coeff3, coeff2, coeff1, coeff0);
    if (!candidates.length) { return; }
    candidates.sort((a, b) => a - b);

    // Cluster near-equal roots: an even-sized cluster is a tangential (double) root that must not
    // be reported so crossing parity stays consistent; an odd-sized cluster is a genuine crossing
    // reported once at its mean.
    let i = 0;
    while (i < candidates.length) {
      let end = i + 1;
      let sum = candidates[i];
      while (end < candidates.length && sameIntersection(candidates[end], candidates[end - 1])) {
        sum += candidates[end];
        ++end;
      }
      const size = end - i;
      i = end;
      if ((size & 1) === 0) { continue; }
      const t = sum / size;
      if (t < tMin || t > tMax) { continue; }
      const gx = ox + t * dx - this.center[0], gy = oy + t * dy - this.center[1], gz = oz + t * dz - this.center[2];
      const lx = gx * this.U[0] + gy * this.U[1] + gz * this.U[2];
      const ly = gx * this.V[0] + gy * this.V[1] + gz * this.V[2];
      const lz = gx * this.W[0] + gy * this.W[1] + gz * this.W[2];
      const rho = Math.hypot(lx, ly);
      if (rho <= K_TOLERANCE) { continue; }
      const phiTube = Math.atan2(lz, rho - this.majorRadius);
      const phiRing = Math.atan2(ly, lx);
      let boundary = false;
      if (this.hasWireTrim) {
        const r = this.pointInTrim(phiRing, phiTube);
        if (!r.contained) { continue; }
        boundary = r.boundary;
      } else {
        if (!angleInSweepRange(phiRing, this.phiStart, this.phiSweep, angularTolerance(this.majorRadius))) { continue; }
        if (!angleInSweepRange(phiTube, this.tubeStart, this.tubeSweep, angularTolerance(this.minorRadius))) { continue; }
      }
      const n = this.localNormal(lx, ly, lz);
      hits.push({ t, nx: n[0], ny: n[1], nz: n[2], boundary, surface: this.index });
    }
  }
}

/// A record kind this port does not support: it occupies its slot, is counted, and is drawn in a
/// distinct colour rather than silently dropped.
export class UnsupportedSurface extends BaseSurface {
  constructor(record, reason) {
    super(record.index, record.type, false);
    this.reason = reason;
    this.aabb = [0, 0, 0, 0, 0, 0];
  }
  appendIntersections() {}
  patchOutline() { return []; }
}

export function buildSurface(record, joinTolerance) {
  if (record.unsupported) { return new UnsupportedSurface(record, record.unsupported); }
  switch (record.type) {
    case SURFACE_TYPE.PLANE: return new PlaneSurface(record, joinTolerance);
    case SURFACE_TYPE.CYLINDER: return new CylinderSurface(record, joinTolerance);
    case SURFACE_TYPE.CONE: return new ConeSurface(record, joinTolerance);
    case SURFACE_TYPE.SPHERE: return new SphereSurface(record, joinTolerance);
    case SURFACE_TYPE.TORUS: return new TorusSurface(record, joinTolerance);
    default: return new UnsupportedSurface(record, `unknown surface type ${record.type}`);
  }
}
