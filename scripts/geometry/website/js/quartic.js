// Real roots of a cubic and a quartic, ported from Detectors/Base/src/BoundedSurface.h
// (solveDepressedCubic / solveQuarticReal). The port is deliberately literal: the Cauchy-bound
// power-of-two normalisation and the dimensionless 32*eps guards are the whole point of that
// code, and a naive Ferrari drops grazing torus rays without them.

export const EPS = 2.220446049250313e-16;
/// How many double-precision roundings a quantity may be away from zero and still be called zero,
/// relative to the terms that produced it (kQuarticEpsilon).
export const K_QUARTIC_EPSILON = 32 * EPS;
export const TWO_PI = 2 * Math.PI;

export const QuarticBranch = { NotAQuartic: 'NotAQuartic', Biquadratic: 'Biquadratic', Resolvent: 'Resolvent' };

/// Real roots of the depressed cubic w^3 + P w + Q = 0, appended to `roots`.
/// The two branches are separated by an exact structural condition (P >= 0 forces a non-negative
/// discriminant), not by a tolerance.
export function solveDepressedCubic(coeffP, coeffQ, roots) {
  const discriminant = coeffQ * coeffQ / 4 + coeffP * coeffP * coeffP / 27;
  if (!(coeffP < 0) || discriminant > 0) {
    const sqrtDiscriminant = Math.sqrt(Math.max(0, discriminant));
    roots.push(Math.cbrt(-0.5 * coeffQ + sqrtDiscriminant) + Math.cbrt(-0.5 * coeffQ - sqrtDiscriminant));
    return;
  }
  // three real roots: coeffP < 0 here, so the trigonometric form is well defined
  const magnitude = 2 * Math.sqrt(-coeffP / 3);
  const cosineArgument = Math.max(-1, Math.min(1, 3 * coeffQ / (coeffP * magnitude)));
  const baseAngle = Math.acos(cosineArgument);
  for (let branch = 0; branch < 3; ++branch) {
    roots.push(magnitude * Math.cos((baseAngle - TWO_PI * branch) / 3));
  }
}

// frexp: the binary exponent of x, such that x = m * 2^e with m in [0.5, 1). frexp(0) has e = 0.
function binaryExponent(x) {
  if (x === 0 || !Number.isFinite(x)) { return 0; }
  const a = Math.abs(x);
  let e = Math.ceil(Math.log2(a));
  // log2 rounding can be off by one at exact powers of two; fix it by testing the mantissa
  while (a / Math.pow(2, e) >= 1) { e += 1; }
  while (a / Math.pow(2, e - 1) < 1) { e -= 1; }
  return e;
}

/// Real roots of a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0 (requires a4 != 0), via Ferrari's method
/// with a resolvent cubic and two Newton polishing steps. A tangential (double) root is returned
/// as a near-equal pair so crossing-parity callers can cluster and drop it.
///
/// The root variable is rescaled by the Cauchy bound rounded UP to a power of two, so every
/// coefficient lands in [-1, 1] and each branch test ("did q cancel to zero", "is the resolvent
/// resolved") compares against 32*eps -- a statement about arithmetic, not about centimetres.
/// The power of two makes the substitution exact, so normalisation cannot change an answer the
/// unnormalised code already got right; only the guards can.
export function solveQuarticReal(a4, a3, a2, a1, a0, branchOut) {
  const note = (b) => { if (branchOut) { branchOut.branch = b; } };
  note(QuarticBranch.NotAQuartic);
  const roots = [];
  if (!(Math.abs(a4) > 0)) { return roots; }

  let coeffB = a3 / a4, coeffC = a2 / a4, coeffD = a1 / a4, coeffE = a0 / a4;
  if (!Number.isFinite(coeffB) || !Number.isFinite(coeffC) || !Number.isFinite(coeffD) || !Number.isFinite(coeffE)) {
    return roots;
  }
  const rootBound = Math.max(Math.abs(coeffB), Math.sqrt(Math.abs(coeffC)),
                             Math.cbrt(Math.abs(coeffD)), Math.sqrt(Math.sqrt(Math.abs(coeffE))));
  const boundExponent = binaryExponent(rootBound);
  const scale = Math.pow(2, boundExponent);
  coeffB /= scale;
  coeffC /= scale * scale;
  coeffD /= scale * scale * scale;
  coeffE /= scale * scale * scale * scale;

  // depress with y = z - b/4: z^4 + p z^2 + q z + r
  const termP = coeffC - 3 * coeffB * coeffB / 8;
  const termQ = coeffD - coeffB * coeffC / 2 + coeffB * coeffB * coeffB / 8;
  const termR = coeffE - coeffB * coeffD / 4 + coeffB * coeffB * coeffC / 16 - 3 * Math.pow(coeffB, 4) / 256;
  const shift = -coeffB / 4;

  const addQuadraticRoots = (quadB, quadC) => {
    const discriminant = quadB * quadB - 4 * quadC;
    if (discriminant < 0) { return; }
    const sq = Math.sqrt(discriminant);
    roots.push(shift + 0.5 * (-quadB - sq));
    roots.push(shift + 0.5 * (-quadB + sq));
  };

  const addBiquadraticRoots = () => {
    const discriminant = termP * termP - 4 * termR;
    if (discriminant < 0) { return; }
    const sq = Math.sqrt(discriminant);
    for (const zSquared of [0.5 * (-termP + sq), 0.5 * (-termP - sq)]) {
      if (zSquared >= 0) {
        const z = Math.sqrt(zSquared);
        roots.push(shift + z);
        roots.push(shift - z);
      }
    }
  };

  // q is zero exactly when it is zero to the precision with which it can be computed; after the
  // substitution, 1 is the unit of this quartic.
  let biquadratic = Math.abs(termQ) <= K_QUARTIC_EPSILON;
  if (!biquadratic) {
    note(QuarticBranch.Resolvent);
    const cubicA2 = termP;
    const cubicA1 = termP * termP / 4 - termR;
    const cubicA0 = -termQ * termQ / 8;
    const cubicP = cubicA1 - cubicA2 * cubicA2 / 3;
    const cubicQ = 2 * cubicA2 * cubicA2 * cubicA2 / 27 - cubicA2 * cubicA1 / 3 + cubicA0;
    const cubicRoots = [];
    solveDepressedCubic(cubicP, cubicQ, cubicRoots);
    let resolvent = 0;
    for (const r of cubicRoots) { resolvent = Math.max(resolvent, r - cubicA2 / 3); }
    // Below a few ulps of the resolvent cubic's own Cauchy bound the computed resolvent is noise;
    // q is then at the same noise level and the biquadratic branch is the better-conditioned
    // answer. Returning nothing turns an ill-conditioned crossing into a missing wall.
    const resolventScale = Math.max(Math.abs(cubicA2), Math.sqrt(Math.abs(cubicA1)), Math.cbrt(Math.abs(cubicA0)));
    if (resolvent > K_QUARTIC_EPSILON * resolventScale) {
      const sqrtTwoResolvent = Math.sqrt(2 * resolvent);
      const linearTerm = sqrtTwoResolvent * termQ / (4 * resolvent);
      addQuadraticRoots(-sqrtTwoResolvent, termP / 2 + resolvent + linearTerm);
      addQuadraticRoots(sqrtTwoResolvent, termP / 2 + resolvent - linearTerm);
    } else {
      biquadratic = true;
    }
  }
  if (biquadratic) {
    note(QuarticBranch.Biquadratic);
    addBiquadraticRoots();
  }

  // Newton polishing against the monic quartic. Every root of the normalised quartic satisfies
  // |z| <= 2 by the Cauchy bound, so a longer step is meaningless whatever produced it -- which
  // also rejects the non-finite step a zero derivative gives.
  const quartic = (x) => (((x + coeffB) * x + coeffC) * x + coeffD) * x + coeffE;
  const quarticDerivative = (x) => ((4 * x + 3 * coeffB) * x + 2 * coeffC) * x + coeffD;
  for (let i = 0; i < roots.length; ++i) {
    for (let it = 0; it < 2; ++it) {
      const step = quartic(roots[i]) / quarticDerivative(roots[i]);
      if (Number.isFinite(step) && Math.abs(step) <= 2) { roots[i] -= step; }
    }
  }
  for (let i = 0; i < roots.length; ++i) { roots[i] *= scale; } // exact: scale is a power of two
  return roots;
}
