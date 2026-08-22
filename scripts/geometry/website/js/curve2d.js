// Curve2D: one trimmed boundary curve in a surface's parametric (u, v) domain -- a straight line,
// a circular arc, or a clamped (rational) B-spline. Ported from the Curve2D value type in
// Detectors/Base/src/BoundedSurface.h.
//
// The three details that matter, and that a casual implementation gets wrong (they each closed a
// real defect in the C++):
//  - endpoints of a B-spline are EVALUATED, never read off the first/last pole: an unclamped knot
//    vector does not interpolate its poles (K1);
//  - the adaptive flattening never calls an interval flat across an interior knot, probes at
//    1/4, 1/2 and 3/4 rather than at the midpoint alone, and switches to distance-to-point when
//    the chord is degenerate (a closed curve has a zero-length chord) (K4);
//  - one seam vertex per join is fixed by the wire and substituted into both neighbours, so
//    winding and point-to-curve distance measure against the same polyline (K5).

export const K_TOLERANCE = 1e-9;
export const K_TOLERANCE_SQ = K_TOLERANCE * K_TOLERANCE;
export const K_BSPLINE_FLATNESS = 1e-5;
export const K_BSPLINE_FLATNESS_SQ = K_BSPLINE_FLATNESS * K_BSPLINE_FLATNESS;
export const TWO_PI = 2 * Math.PI;
export const HALF_PI = 0.5 * Math.PI;

export const CurveKind = { Line: 'Line', Arc: 'Arc', BSpline: 'BSpline' };

export function v2(u, v) { return { u, v }; }
export function distanceSq2(a, b) { const du = a.u - b.u, dv = a.v - b.v; return du * du + dv * dv; }

export function pointSegmentDistanceSq(point, segStart, segEnd) {
  const su = segEnd.u - segStart.u, sv = segEnd.v - segStart.v;
  const lengthSq = su * su + sv * sv;
  if (lengthSq <= K_TOLERANCE_SQ) { return distanceSq2(point, segStart); }
  let t = ((point.u - segStart.u) * su + (point.v - segStart.v) * sv) / lengthSq;
  t = Math.max(0, Math.min(1, t));
  const du = point.u - (segStart.u + t * su), dv = point.v - (segStart.v + t * sv);
  return du * du + dv * dv;
}

export class Curve2D {
  constructor() {
    this.kind = CurveKind.Line;
    this.lineStart = v2(0, 0);
    this.lineEnd = v2(0, 0);
    this.center = v2(0, 0);
    this.radius = 0;
    this.startAngle = 0;
    this.endAngle = 0;
    this.degree = 0;
    this.poles = [];
    this.weights = [];
    this.knots = [];
    this._samples = null;
    this._box = null;
    this._index = null;
    this.canonicalStart = null;
    this.canonicalEnd = null;
  }

  static makeLine(start, end) {
    const c = new Curve2D();
    c.kind = CurveKind.Line; c.lineStart = { ...start }; c.lineEnd = { ...end };
    return c;
  }
  static makeArc(center, radius, startAngle, endAngle) {
    const c = new Curve2D();
    c.kind = CurveKind.Arc; c.center = { ...center }; c.radius = radius;
    c.startAngle = startAngle; c.endAngle = endAngle;
    return c;
  }
  static makeBSpline(degree, poles, weights, knots) {
    const c = new Curve2D();
    c.kind = CurveKind.BSpline; c.degree = degree;
    c.poles = poles.map(p => (Array.isArray(p) ? v2(p[0], p[1]) : { ...p }));
    c.weights = weights ? Array.from(weights) : [];
    c.knots = Array.from(knots);
    return c;
  }

  get isLine() { return this.kind === CurveKind.Line; }
  get isArc() { return this.kind === CurveKind.Arc; }
  get isBSpline() { return this.kind === CurveKind.BSpline; }
  sweep() { return this.endAngle - this.startAngle; }

  setCanonicalEndpoints(start, end) {
    this.canonicalStart = { ...start };
    this.canonicalEnd = { ...end };
    this._box = null;
    this._index = null;
    if (this.kind === CurveKind.BSpline) { this._samples = null; }
  }
  get hasCanonical() { return this.canonicalStart !== null; }
  loopStart() { return this.hasCanonical ? this.canonicalStart : this.startPoint(); }
  loopEnd() { return this.hasCanonical ? this.canonicalEnd : this.endPoint(); }

  // --- B-spline machinery (Piegl & Tiller) -------------------------------------------------
  bsplineT0() { return this.knots[this.degree]; }
  bsplineT1() { return this.knots[this.poles.length]; }

  bsplineSpan(t) {
    const lastPole = this.poles.length - 1;
    if (t >= this.knots[lastPole + 1]) { return lastPole; }
    if (t <= this.knots[this.degree]) { return this.degree; }
    let low = this.degree, high = lastPole + 1, mid = (low + high) >> 1;
    while (t < this.knots[mid] || t >= this.knots[mid + 1]) {
      if (t < this.knots[mid]) { high = mid; } else { low = mid; }
      mid = (low + high) >> 1;
    }
    return mid;
  }

  /// Non-zero degree-p basis functions and their first derivatives (DersBasisFuns, k = 1 only).
  bsplineBasis(span, t, basis, basisDeriv) {
    const p = this.degree, knots = this.knots;
    const ndu = [];
    for (let i = 0; i <= p; ++i) { ndu.push(new Float64Array(p + 1)); }
    const left = new Float64Array(p + 1), right = new Float64Array(p + 1);
    ndu[0][0] = 1;
    for (let j = 1; j <= p; ++j) {
      left[j] = t - knots[span + 1 - j];
      right[j] = knots[span + j] - t;
      let saved = 0;
      for (let r = 0; r < j; ++r) {
        ndu[j][r] = right[r + 1] + left[j - r];
        const temp = ndu[r][j - 1] / ndu[j][r];
        ndu[r][j] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    for (let j = 0; j <= p; ++j) { basis[j] = ndu[j][p]; }
    for (let r = 0; r <= p; ++r) {
      let d = 0;
      const pk = p - 1;
      if (r >= 1) { d += (1 / ndu[pk + 1][r - 1]) * ndu[r - 1][pk]; }
      if (r <= pk) { d += (-1 / ndu[pk + 1][r]) * ndu[r][pk]; }
      basisDeriv[r] = d * p;
    }
  }

  /// (Rational) B-spline point and knot-parameter derivative at knot value t.
  bsplineEval(t) {
    const p = this.degree;
    const span = this.bsplineSpan(t);
    const basis = new Float64Array(p + 1), basisDeriv = new Float64Array(p + 1);
    this.bsplineBasis(span, t, basis, basisDeriv);
    let su = 0, sv = 0, du = 0, dv = 0, wTotal = 0, wDeriv = 0;
    for (let j = 0; j <= p; ++j) {
      const idx = span - p + j;
      const w = this.weights.length ? this.weights[idx] : 1;
      const pole = this.poles[idx];
      su += basis[j] * w * pole.u; sv += basis[j] * w * pole.v;
      wTotal += basis[j] * w;
      du += basisDeriv[j] * w * pole.u; dv += basisDeriv[j] * w * pole.v;
      wDeriv += basisDeriv[j] * w;
    }
    const invW = Math.abs(wTotal) > K_TOLERANCE ? 1 / wTotal : 0;
    return {
      point: v2(su * invW, sv * invW),
      deriv: v2((du * wTotal - su * wDeriv) * invW * invW, (dv * wTotal - sv * wDeriv) * invW * invW),
    };
  }

  bsplinePointAt(parameter) {
    const t = this.bsplineT0() + parameter * (this.bsplineT1() - this.bsplineT0());
    return this.bsplineEval(t).point;
  }

  /// Does the knot interval (lowT, highT) contain a knot strictly inside it? An interval that
  /// still does is never called flat, whatever the probes say.
  spansInteriorKnot(lowT, highT) {
    const firstInterior = this.degree + 1;
    const endInterior = Math.min(this.poles.length, this.knots.length);
    for (let i = firstInterior; i < endInterior; ++i) {
      const k = this.knots[i];
      if (k > lowT && k < highT) { return true; }
    }
    return false;
  }

  bsplineSampleInto(samples, flatnessSq = K_BSPLINE_FLATNESS_SQ, maxDepth = 16) {
    const t0 = this.bsplineT0(), t1 = this.bsplineT1();
    const p0 = this.bsplineEval(t0).point;
    const p1 = this.bsplineEval(t1).point;
    samples.push(p0);
    this._sampleRecursive(t0, t1, p0, p1, flatnessSq, maxDepth, samples);
  }

  _sampleRecursive(t0, t1, p0, p1, flatnessSq, depth, samples) {
    const tMid = 0.5 * (t0 + t1);
    const midPoint = this.bsplineEval(tMid).point;
    // A short chord alone must NOT end the recursion: a closed curve has p0 == p1 exactly, so its
    // chord is zero however much curve lies between them.
    const degenerateChord = distanceSq2(p0, p1) <= flatnessSq;
    const deviationSq = (pt) => (degenerateChord ? distanceSq2(pt, p0) : pointSegmentDistanceSq(pt, p0, p1));
    // Three interior probes, not one: a curve symmetric about its parameter midpoint defeats a
    // single midpoint probe exactly.
    let flatness = deviationSq(midPoint);
    for (const fraction of [0.25, 0.75]) {
      flatness = Math.max(flatness, deviationSq(this.bsplineEval(t0 + (t1 - t0) * fraction).point));
    }
    if (depth <= 0 || (flatness <= flatnessSq && !this.spansInteriorKnot(t0, t1))) {
      samples.push(p1);
      return;
    }
    this._sampleRecursive(t0, tMid, p0, midPoint, flatnessSq, depth - 1, samples);
    this._sampleRecursive(tMid, t1, midPoint, p1, flatnessSq, depth - 1, samples);
  }

  /// The cached flattened polyline: one polyline, and it is the canonical one, so winding, closest
  /// point and any visualization all see the identical boundary.
  bsplineSamples() {
    if (this._samples === null) {
      const s = [];
      this.bsplineSampleInto(s);
      if (this.hasCanonical && s.length >= 2) {
        s[0] = this.canonicalStart;
        s[s.length - 1] = this.canonicalEnd;
      }
      this._samples = s;
    }
    return this._samples;
  }

  // --- geometry ------------------------------------------------------------------------------
  pointAtAngle(angle) {
    return v2(this.center.u + this.radius * Math.cos(angle), this.center.v + this.radius * Math.sin(angle));
  }

  pointAt(parameter) {
    if (this.kind === CurveKind.Line) {
      return v2(this.lineStart.u + parameter * (this.lineEnd.u - this.lineStart.u),
                this.lineStart.v + parameter * (this.lineEnd.v - this.lineStart.v));
    }
    if (this.kind === CurveKind.BSpline) { return this.bsplinePointAt(parameter); }
    return this.pointAtAngle(this.startAngle + parameter * this.sweep());
  }

  startPoint() {
    if (this.kind === CurveKind.Line) { return this.lineStart; }
    if (this.kind === CurveKind.BSpline) { return this.bsplineEval(this.bsplineT0()).point; }
    return this.pointAtAngle(this.startAngle);
  }
  endPoint() {
    if (this.kind === CurveKind.Line) { return this.lineEnd; }
    if (this.kind === CurveKind.BSpline) { return this.bsplineEval(this.bsplineT1()).point; }
    return this.pointAtAngle(this.endAngle);
  }

  angleInSweep(angle) {
    const total = this.sweep();
    const magnitude = Math.abs(total);
    if (magnitude >= TWO_PI - K_TOLERANCE) { return true; }
    let delta = total >= 0 ? angle - this.startAngle : this.startAngle - angle;
    delta -= TWO_PI * Math.floor(delta / TWO_PI);
    return delta <= magnitude + K_TOLERANCE;
  }

  angleParameter(angle) {
    const total = this.sweep();
    if (Math.abs(total) <= K_TOLERANCE) { return 0; }
    let delta = total >= 0 ? angle - this.startAngle : this.startAngle - angle;
    delta -= TWO_PI * Math.floor(delta / TWO_PI);
    return Math.max(0, Math.min(1, delta / Math.abs(total)));
  }

  /// Endpoints plus, for an arc, the axis-extreme points inside the sweep.
  _includeAnalyticExtremes(include) {
    if (this.kind === CurveKind.Line) { include(this.lineStart); include(this.lineEnd); return; }
    if (this.kind === CurveKind.BSpline) { include(this.startPoint()); include(this.endPoint()); return; }
    include(this.startPoint()); include(this.endPoint());
    for (let k = -4; k <= 4; ++k) {
      for (const base of [0, HALF_PI, Math.PI, 1.5 * Math.PI]) {
        const angle = base + k * TWO_PI;
        if (this.angleInSweep(angle)) { include(this.pointAtAngle(angle)); }
      }
    }
  }

  /// Conservative extent: a B-spline contributes its control-pole hull.
  extendBounds(lower, upper) {
    const include = (p) => {
      lower.u = Math.min(lower.u, p.u); lower.v = Math.min(lower.v, p.v);
      upper.u = Math.max(upper.u, p.u); upper.v = Math.max(upper.v, p.v);
    };
    if (this.kind === CurveKind.BSpline) { for (const p of this.poles) { include(p); } return; }
    this._includeAnalyticExtremes(include);
  }

  /// Tight extent: a B-spline contributes its sampled polyline instead of its pole hull.
  extendTightBounds(lower, upper) {
    const include = (p) => {
      lower.u = Math.min(lower.u, p.u); lower.v = Math.min(lower.v, p.v);
      upper.u = Math.max(upper.u, p.u); upper.v = Math.max(upper.v, p.v);
    };
    if (this.kind === CurveKind.BSpline) { for (const p of this.bsplineSamples()) { include(p); } return; }
    this._includeAnalyticExtremes(include);
  }

  closestPoint(point) {
    if (this.kind === CurveKind.BSpline) {
      const poly = this.bsplineSamples();
      if (poly.length < 2) { return this.startPoint(); }
      let best = poly[0], bestSq = Infinity;
      for (let i = 0; i + 1 < poly.length; ++i) {
        const a = poly[i], b = poly[i + 1];
        const su = b.u - a.u, sv = b.v - a.v;
        const lenSq = su * su + sv * sv;
        let t = 0;
        if (lenSq > K_TOLERANCE_SQ) {
          t = Math.max(0, Math.min(1, ((point.u - a.u) * su + (point.v - a.v) * sv) / lenSq));
        }
        const cand = v2(a.u + t * su, a.v + t * sv);
        const dsq = distanceSq2(point, cand);
        if (dsq < bestSq) { bestSq = dsq; best = cand; }
      }
      return best;
    }
    if (this.kind === CurveKind.Line) {
      const su = this.lineEnd.u - this.lineStart.u, sv = this.lineEnd.v - this.lineStart.v;
      const lenSq = su * su + sv * sv;
      if (lenSq <= K_TOLERANCE_SQ) { return this.lineStart; }
      const t = Math.max(0, Math.min(1, ((point.u - this.lineStart.u) * su + (point.v - this.lineStart.v) * sv) / lenSq));
      return v2(this.lineStart.u + t * su, this.lineStart.v + t * sv);
    }
    // arc: project radially onto the circle, then clamp the angle to the sweep
    const du = point.u - this.center.u, dv = point.v - this.center.v;
    if (du * du + dv * dv <= K_TOLERANCE_SQ) { return this.startPoint(); }
    const angle = Math.atan2(dv, du);
    if (this.angleInSweep(angle)) { return this.pointAtAngle(angle); }
    const sp = this.startPoint(), ep = this.endPoint();
    return distanceSq2(point, sp) <= distanceSq2(point, ep) ? sp : ep;
  }

  distanceSq(point) { return distanceSq2(point, this.closestPoint(point)); }

  /// A scanline index over the flattened polyline: the v range is cut into buckets and each
  /// segment is listed in every bucket its v span touches. Both per-point queries -- "does a +u
  /// scanline at this v cross the curve" and "is this point within the on-boundary band" -- only
  /// ever look at one v, so this turns a walk over ~400 chords into a walk over a handful. It is
  /// an index over the same polyline, so it changes no answer.
  _scanIndex() {
    if (this._index) { return this._index; }
    const poly = this.bsplineSamples();
    const n = poly.length;
    const box = this.tightBox();
    const vMin = box.lower.v, vMax = box.upper.v;
    const span = vMax - vMin;
    const buckets = Math.max(1, Math.min(256, Math.ceil(n / 4)));
    const lists = [];
    for (let b = 0; b < buckets; ++b) { lists.push([]); }
    const bucketOf = (v) => {
      if (!(span > 0)) { return 0; }
      return Math.max(0, Math.min(buckets - 1, Math.floor((v - vMin) / span * buckets)));
    };
    for (let i = 0; i + 1 < n; ++i) {
      const lo = bucketOf(Math.min(poly[i].v, poly[i + 1].v));
      const hi = bucketOf(Math.max(poly[i].v, poly[i + 1].v));
      for (let b = lo; b <= hi; ++b) { lists[b].push(i); }
    }
    this._index = { poly, buckets, vMin, vMax, span, lists, bucketOf };
    return this._index;
  }

  /// True when the point lies within sqrt(bandSq) of the curve. Exact for lines and arcs; for a
  /// B-spline it is the distance to the flattened polyline, which is what the boundary band is
  /// sized from anyway.
  withinBandSq(point, bandSq) {
    if (this.boxDistanceSq(point) > bandSq) { return false; }
    if (this.kind !== CurveKind.BSpline) { return this.distanceSq(point) <= bandSq; }
    const idx = this._scanIndex();
    const band = Math.sqrt(bandSq);
    const loBucket = idx.bucketOf(point.v - band), hiBucket = idx.bucketOf(point.v + band);
    for (let b = loBucket; b <= hiBucket; ++b) {
      for (const i of idx.lists[b]) {
        if (pointSegmentDistanceSq(point, idx.poly[i], idx.poly[i + 1]) <= bandSq) { return true; }
      }
    }
    return false;
  }

  /// Tight (u, v) box of the curve itself, cached. Used only to prune the two per-point queries
  /// -- a point far outside the box is neither within the on-boundary band nor able to have its
  /// +u scanline cross the curve -- so it must be tight and correct, never conservative-large.
  tightBox() {
    if (!this._box) {
      const lower = v2(Infinity, Infinity), upper = v2(-Infinity, -Infinity);
      this.extendTightBounds(lower, upper);
      // The loop-canonical seam vertices are part of the boundary this curve stands for, and they
      // can sit a join tolerance outside its own geometry; a box that excluded them would prune
      // away a crossing the unpruned walk counts.
      for (const p of [this.loopStart(), this.loopEnd()]) {
        lower.u = Math.min(lower.u, p.u); lower.v = Math.min(lower.v, p.v);
        upper.u = Math.max(upper.u, p.u); upper.v = Math.max(upper.v, p.v);
      }
      this._box = { lower, upper };
    }
    return this._box;
  }

  /// Squared distance to the box, 0 when inside: a lower bound on distanceSq that costs nothing.
  boxDistanceSq(point) {
    const b = this.tightBox();
    const du = Math.max(b.lower.u - point.u, 0, point.u - b.upper.u);
    const dv = Math.max(b.lower.v - point.v, 0, point.v - b.upper.v);
    return du * du + dv * dv;
  }

  /// True when no +u scanline through `point` can cross this curve.
  scanlineMisses(point) {
    const b = this.tightBox();
    return point.v < b.lower.v || point.v > b.upper.v || point.u >= b.upper.u;
  }

  /// How far this curve's representation can sit from the curve it stands for, in parametric
  /// units. Lines and arcs are exact; a B-spline is only as good as kBSplineFlatness.
  representationTolerance() { return this.kind === CurveKind.BSpline ? K_BSPLINE_FLATNESS : 0; }

  /// Number of times a horizontal ray from `point` towards +u crosses this curve, half-open so
  /// that shared endpoints and tangent extrema are not double-counted.
  rightwardCrossings(point, canonicalStart, canonicalEnd) {
    const cs = canonicalStart || this.loopStart();
    const ce = canonicalEnd || this.loopEnd();

    if (this.kind === CurveKind.Line) {
      const firstAbove = cs.v > point.v, secondAbove = ce.v > point.v;
      if (firstAbove === secondAbove) { return 0; }
      const intersectU = cs.u + (point.v - cs.v) * (ce.u - cs.u) / (ce.v - cs.v);
      return point.u < intersectU ? 1 : 0;
    }

    if (this.kind === CurveKind.BSpline) {
      // The flattened polyline already ends on the loop-canonical vertices, so the same half-open
      // segment rule applies with no substitution here. Only the segments in this v's bucket can
      // straddle the scanline, so only they are walked.
      const poly = this.bsplineSamples();
      if (poly.length < 2) { return 0; }
      const idx = this._scanIndex();
      let crossings = 0;
      for (const i of idx.lists[idx.bucketOf(point.v)]) {
        const a = poly[i], b = poly[i + 1];
        const firstAbove = a.v > point.v, secondAbove = b.v > point.v;
        if (firstAbove === secondAbove) { continue; }
        const intersectU = a.u + (point.v - a.v) * (b.u - a.u) / (b.v - a.v);
        if (point.u < intersectU) { ++crossings; }
      }
      return crossings;
    }

    // Split the arc into v-monotonic sub-arcs at its top/bottom extreme angles (theta = pi/2 + k pi).
    // On each sub-arc cos(theta) keeps a constant sign, so the exact u at v = point.v is
    // centre.u +/- r*sqrt(1 - sin^2), with the sign taken from the sub-arc.
    const totalSweep = this.sweep();
    if (Math.abs(totalSweep) <= K_TOLERANCE || this.radius <= K_TOLERANCE) { return 0; }
    const breaks = [0];
    const lowAngle = Math.min(this.startAngle, this.endAngle);
    const highAngle = Math.max(this.startAngle, this.endAngle);
    const firstK = Math.floor((lowAngle - HALF_PI) / Math.PI) - 1;
    const lastK = Math.ceil((highAngle - HALF_PI) / Math.PI) + 1;
    for (let k = firstK; k <= lastK && breaks.length < 7; ++k) {
      const extremeAngle = HALF_PI + k * Math.PI;
      if (extremeAngle <= lowAngle + K_TOLERANCE || extremeAngle >= highAngle - K_TOLERANCE) { continue; }
      const p = (extremeAngle - this.startAngle) / totalSweep;
      if (p > K_TOLERANCE && p < 1 - K_TOLERANCE) { breaks.push(p); }
    }
    breaks.push(1);
    breaks.sort((a, b) => a - b);

    let ratio = (point.v - this.center.v) / this.radius;
    ratio = Math.max(-1, Math.min(1, ratio));
    const cosMagnitude = Math.sqrt(Math.max(0, 1 - ratio * ratio));

    let crossings = 0;
    for (let i = 0; i + 1 < breaks.length; ++i) {
      const subStart = i === 0 ? cs : this.pointAt(breaks[i]);
      const subEnd = i + 2 === breaks.length ? ce : this.pointAt(breaks[i + 1]);
      const midAngle = this.startAngle + 0.5 * (breaks[i] + breaks[i + 1]) * totalSweep;
      const cosSign = Math.cos(midAngle) >= 0 ? 1 : -1;
      const intersectU = this.center.u + cosSign * this.radius * cosMagnitude;
      const firstAbove = subStart.v > point.v, secondAbove = subEnd.v > point.v;
      if (firstAbove !== secondAbove && point.u < intersectU) { ++crossings; }
    }
    return crossings;
  }

  /// Ordered polyline of this curve for display: lines give their endpoints, arcs are chorded,
  /// B-splines give the adaptive flattening. The closing point is NOT appended (the next curve's
  /// start reproduces it).
  displaySamples(segmentsPerTurn = 48) {
    if (this.kind === CurveKind.Line) { return [this.loopStart()]; }
    if (this.kind === CurveKind.BSpline) {
      const s = this.bsplineSamples();
      return s.slice(0, s.length - 1);
    }
    const steps = Math.max(1, Math.round(segmentsPerTurn * Math.abs(this.sweep()) / TWO_PI));
    const out = [];
    for (let i = 0; i < steps; ++i) { out.push(this.pointAt(i / steps)); }
    return out;
  }
}
