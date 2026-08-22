// Closed trim loops in a surface's parametric domain, and point classification against them.
// Ported from SurfaceWire / CurveWire / curveTrimContains in Detectors/Base/src/BoundedSurface.h.
//
// Containment is winding by half-open scanline crossings -- a topological count that needs no
// metric at all. The metric enters only to size the on-boundary band, because kTolerance is a
// length in cm and these domains mix radians with centimetres.
//
// One deliberate omission against the C++: the wire's winding is NOT normalised (outer CCW, inner
// CW). Orientation changes no crossing parity and nothing here integrates a signed area, so the
// loop is used in the order the file gives it, and the outer/inner role is taken from the file.

import { pointSegmentDistanceSq, v2, K_TOLERANCE, TWO_PI } from './curve2d.js';

export const WireClass = { Inside: 'Inside', Outside: 'Outside', Boundary: 'Boundary' };

/// A surface's first fundamental form, as a function of (u, v). `g` returns [gUU, gUV, gVV].
export class ParametricMetric {
  constructor(g) { this.g = g || null; }
  lengthSq(uv, du, dv) {
    if (!this.g) { return du * du + dv * dv; }
    const [gUU, gUV, gVV] = this.g(uv.u, uv.v);
    return gUU * du * du + 2 * gUV * du * dv + gVV * dv * dv;
  }
  distanceSq(from, to) { return this.lengthSq(from, to.u - from.u, to.v - from.v); }
  /// The largest 3D length a unit parametric displacement can span at (u, v): the square root of
  /// the larger eigenvalue of the form.
  maxScale(uv) {
    if (!this.g) { return 1; }
    const [gUU, gUV, gVV] = this.g(uv.u, uv.v);
    const trace = gUU + gVV;
    const det = gUU * gVV - gUV * gUV;
    const disc = Math.max(0, trace * trace - 4 * det);
    return Math.sqrt(Math.max(0, 0.5 * (trace + Math.sqrt(disc))));
  }
}

export const IDENTITY_METRIC = new ParametricMetric(null);

// --- the first fundamental forms, by family (the closed forms the kernel and the loader share) --
export function planeParametricMetric(axisU, axisV) {
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return [dot(axisU, axisU), dot(axisU, axisV), dot(axisV, axisV)];
}
export function cylinderParametricMetric(radius) { return [radius * radius, 0, 1]; }
export function coneParametricMetric(radiusAtHeight, slope) {
  return [radiusAtHeight * radiusAtHeight, 0, 1 + slope * slope];
}
export function sphereParametricMetric(radius, theta) {
  const s = radius * Math.sin(theta);
  return [s * s, 0, radius * radius];
}
export function torusParametricMetric(majorRadius, minorRadius, phiTube) {
  const ring = majorRadius + minorRadius * Math.cos(phiTube);
  return [ring * ring, 0, minorRadius * minorRadius];
}

/// One closed loop of Curve2D segments. Role 0 = outer, 1 = inner (hole).
export class CurveWire {
  constructor(curves, role, metric) {
    this.curves = curves;
    this.role = role;
    this.metric = metric || IDENTITY_METRIC;
    // Fix one vertex value per seam and give it to both curves that meet there, so winding and
    // point-to-curve distance measure against one polyline rather than two that differ by the
    // join drift.
    const n = curves.length;
    for (let i = 0; i < n; ++i) {
      curves[i].setCanonicalEndpoints(curves[i].startPoint(), curves[(i + 1) % n].startPoint());
    }
    this._repTolerance = 0;
    for (const c of curves) { this._repTolerance = Math.max(this._repTolerance, c.representationTolerance()); }
  }

  representationTolerance() { return this._repTolerance; }

  /// The half-width of the on-boundary band at `point`, in parametric units: kTolerance converted
  /// through the metric's largest scale, floored by what the representation actually is.
  boundaryBand(point) {
    const scale = this.metric.maxScale(point);
    const lengthFloor = scale > K_TOLERANCE ? K_TOLERANCE / scale : 0;
    return Math.max(lengthFloor, this._repTolerance);
  }

  classify(point) {
    const band = this.boundaryBand(point);
    const bandSq = band * band;
    for (const curve of this.curves) {
      // withinBandSq prunes by the curve's own box first and then, for a B-spline, walks only the
      // chords in this v's bucket -- the band is tiny, so that is a handful of them.
      if (curve.withinBandSq(point, bandSq)) { return WireClass.Boundary; }
    }
    let crossings = 0;
    for (const curve of this.curves) {
      if (curve.scanlineMisses(point)) { continue; }
      crossings += curve.rightwardCrossings(point, curve.loopStart(), curve.loopEnd());
    }
    return (crossings % 2 === 1) ? WireClass.Inside : WireClass.Outside;
  }

  parametricBounds(lower, upper) { for (const c of this.curves) { c.extendBounds(lower, upper); } }
  tightParametricBounds(lower, upper) { for (const c of this.curves) { c.extendTightBounds(lower, upper); } }

  /// Ordered, closed display polyline of the loop.
  displayPolyline(segmentsPerTurn = 48) {
    const out = [];
    for (const c of this.curves) { for (const p of c.displaySamples(segmentsPerTurn)) { out.push(p); } }
    if (out.length) { out.push(out[0]); }
    return out;
  }
}

/// A polygon wire (a plane whose trim is pure line segments), stored as its vertex ring.
export class PolygonWire {
  constructor(vertices, role, metric) {
    this.vertices = vertices;
    this.role = role;
    this.metric = metric || IDENTITY_METRIC;
  }
  representationTolerance() { return 0; }
  classify(point) {
    const scale = this.metric.maxScale(point);
    const band = scale > K_TOLERANCE ? K_TOLERANCE / scale : 0;
    const bandSq = band * band;
    let inside = false;
    const n = this.vertices.length;
    for (let i = 0; i < n; ++i) {
      const a = this.vertices[i], b = this.vertices[(i + 1) % n];
      if (pointSegmentDistanceSq(point, a, b) <= bandSq) { return WireClass.Boundary; }
      if ((a.v > point.v) !== (b.v > point.v)) {
        const intersectU = a.u + (point.v - a.v) * (b.u - a.u) / (b.v - a.v);
        if (point.u < intersectU) { inside = !inside; }
      }
    }
    return inside ? WireClass.Inside : WireClass.Outside;
  }
  parametricBounds(lower, upper) {
    for (const p of this.vertices) {
      lower.u = Math.min(lower.u, p.u); lower.v = Math.min(lower.v, p.v);
      upper.u = Math.max(upper.u, p.u); upper.v = Math.max(upper.v, p.v);
    }
  }
  tightParametricBounds(lower, upper) { this.parametricBounds(lower, upper); }
  displayPolyline() { return this.vertices.concat(this.vertices.length ? [this.vertices[0]] : []); }
}

/// Inside or on the boundary of the outer loop, and not strictly inside any hole.
/// Returns { contained, boundary }.
export function trimContains(outerWire, innerWires, point) {
  const outer = outerWire.classify(point);
  if (outer === WireClass.Outside) { return { contained: false, boundary: false }; }
  if (outer === WireClass.Boundary) { return { contained: true, boundary: true }; }
  for (const inner of innerWires) {
    const cls = inner.classify(point);
    if (cls === WireClass.Boundary) { return { contained: true, boundary: true }; }
    if (cls === WireClass.Inside) { return { contained: false, boundary: false }; }
  }
  return { contained: true, boundary: false };
}

/// Shift `angle` by whole turns to lie as close as possible to the window [uMin, uMax], so a raw
/// atan2 result can be classified against a trim wire whose phi range straddles the branch cut.
export function unwrapAngleInto(angle, uMin, uMax) {
  const windowCenter = 0.5 * (uMin + uMax);
  return angle - TWO_PI * Math.round((angle - windowCenter) / TWO_PI);
}

/// Build the trim window from an outer loop, mirroring buildCurveTrim: the conservative pole-hull
/// box, re-measured on the curves themselves when it exceeds a full turn in u.
export function trimWindow(outerWire) {
  const lower = v2(Infinity, Infinity), upper = v2(-Infinity, -Infinity);
  outerWire.parametricBounds(lower, upper);
  if (!Number.isFinite(lower.u) || !Number.isFinite(upper.u)) { return null; }
  if (upper.u - lower.u > TWO_PI + K_TOLERANCE) {
    const tl = v2(Infinity, Infinity), tu = v2(-Infinity, -Infinity);
    outerWire.tightParametricBounds(tl, tu);
    if (!Number.isFinite(tl.u) || tu.u - tl.u > TWO_PI + K_TOLERANCE) { return null; }
    return { lower: tl, upper: tu };
  }
  return { lower, upper };
}
