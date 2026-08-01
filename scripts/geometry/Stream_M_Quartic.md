# Stream M — the quartic solver's dimensionally-inconsistent guards

Date: 2026-08-02. Branch `swenzel/bvhsurfacesolid`. Companion to
[`Stream_E_Scale.md`](Stream_E_Scale.md) §5 (which found the defect and scoped the fix),
[`Stream_L_ALICE3Defect.md`](Stream_L_ALICE3Defect.md) §3 (which reduced it to five coefficients on
real, unscaled geometry), [`Stream_J_XRay.md`](Stream_J_XRay.md) §6.2 (the instrument that detects
it, and the warning about `--beams`) and [`TolerancePolicy.md`](TolerancePolicy.md) §14 (the
register entry for the new criterion).

**In one paragraph.** `solveQuarticReal` decided all three of its branches by comparing
`kTolerance` — 1e-9 **cm, a length** — against quantities that are not lengths: `termQ` scales as
L³, the resolvent as L², the Newton derivative as L³. The repair is not a better constant. It is to
**remove the dimensions**: the monic coefficients give a Cauchy bound on the roots, and rounding
that *up to a power of two* and substituting `x = scale · y` leaves every coefficient in [−1, 1], so
`p`, `q`, `r`, the resolvent and the derivative become dimensionless O(1) numbers and a threshold in
units of the machine epsilon is a statement about arithmetic rather than about centimetres. The
power of two is what makes the substitution free — every scaling is then exact, so on any input
where the old code was already right it performs the same operations with the same roundings and
only the exponents shifted, and **the guards are the only thing that can change**. Measured: the
ALICE3 production-scale reproducer goes from **0 roots to 2 correct ones**; the ×0.1 `torus_union_cyl`
X-ray goes **LOST 4 → 0, parity 2 → 0**; **ALICE3 goes LOST 14 → 0, unterminated 2 → 0, mode-(a)-vs-(b)
6 → 0, and 17/18 parts clean → 18/18**; the ×1 ladder parity audit goes **2 → 0**. Both oracle gates
are unchanged on every verdict and every disagreement count, and their `gate.json` reports differ in
**9 leaf fields of 8487 on the fixtures and 0 of 15140 on Bagger** — all nine on the one toroidal
part, all in the mesh-reference column, all in the last two or three digits.

---

## 1. The defect, restated as three lines of code

```cpp
if (std::abs(termQ) <= kTolerance)      // termQ  ~ L^3 ; selects the biquadratic branch
if (resolvent > kTolerance)             // resolvent ~ L^2 ; if it fails, NO ROOTS ARE PUSHED
if (std::abs(derivative) > kTolerance)  // derivative ~ L^3 ; licenses a Newton polishing step
```

Two failure modes, both confirmed before anything was changed:

* the resolvent test fails ⇒ the function returns the **empty root set**, so a ray silently misses a
  torus it does cross;
* the `termQ` test misroutes an asymmetric quartic into the biquadratic branch, which *assumes*
  `q = 0` and therefore forces its roots symmetric about `−b/4` ⇒ it returns **confidently wrong
  roots**. That is worse than a miss, because a miss at least leaves a visible gap.

**The trigger is the ratio of the ray's lever arm to the feature it hits, not the model's scale.**
Stream E measured it at ×0.1 on a fixture and Stream L reduced it to a single call on **unscaled**
ALICE3 geometry — a ray 375 cm from a 0.1 cm tube.

### 1.1 The red test, and its recorded failure

Five cases went in first, in a delimited `// --- Stream M: ... ---` block in
`Detectors/Base/test/testBVHSurfaceSolid.cxx`, and were run against the unchanged solver. Each is
the first failure of its case; the full run reported **5 of 91 cases failed, 68 of 404 assertions
failed**:

```
StreamM_QuarticFindsTheALICE3ProductionScaleRoots
  fatal error: critical check roots.size() == 2u has failed [0 != 2]

StreamM_QuarticIsScaleInvariantOnAnAsymmetricQuartic
  error: check std::abs(roots[i] - expected[i]) < 1.e-9 * std::abs(expected[i]) has failed
         [1.7093415414235272e-08 >= 1.0000000000000002e-12]
  error: check quarticBackwardError(c, root) < 1.e-12 has failed [1.068390966516601e-06 >= 1e-12]

StreamM_QuarticStillSelectsTheBiquadraticBranch
  fatal error: critical check roots.size() == 4u has failed [0 != 4]

StreamM_QuarticStillDeclinesDegenerateConfigurations
  error: check branch == surf::QuarticBranch::Resolvent has failed

StreamM_QuarticHasNoCliffAsAQuarticApproachesBiquadratic
  error: check std::abs(roots[i] - expected[i]) < 1.e-8 * std::abs(expected[i]) has failed
         [2.564730126064762e-06 >= 2e-12]
  error: check quarticBackwardError(c, root) < 1.e-12 has failed [0.0075646312136552077 >= 1e-12]
```

Read the fourth one: `(x² + k²)((x + k)² + 4k²)` has no real roots at any `k`, so the *count* is 0
either way and only the **branch** assertion catches that the shipped code routes it into the
biquadratic branch at small `k`. That is why the branch is observable at all — see §2.4.

---

## 2. The fix

### 2.1 Normalise the root variable, exactly

```cpp
const double rootBound = std::max({|b|, sqrt(|c|), cbrt(|d|), sqrt(sqrt(|e|))});  // Cauchy: |x| <= 2*bound
int boundExponent = 0; std::frexp(rootBound, &boundExponent);
const double scale = std::ldexp(1., boundExponent);           // round UP to a power of two
```

After dividing `b, c, d, e` by `scale¹⁻⁴`, every coefficient lies in [−1, 1] with at least one of
the four bound-defining quantities in (0.5, 1]. So `|p| ≤ 1.375`, `|q| ≤ 1.625`, `|r| ≤ 1.32`, the
resolvent cubic's coefficients are O(1), the roots satisfy `|z| ≤ 2`, and **1 is the unit of the
problem**. A dimensionless epsilon now means something.

**Why the power of two matters more than the bound does.** Scaling by 2^k is exact in binary
floating point, and so are all the derived scalings: `sqrt` of an even power, `cbrt` of a multiple
of three (the resolvent cubic's coefficients scale as `scale³`, `scale⁴`, `scale⁶`, all handled).
So every intermediate quantity in the algorithm is exactly the old one with its exponent shifted,
every rounding is the same rounding, and scaling the roots back is exact too. **The normalisation
cannot change an answer. Only the guards can.** That is the whole inertness argument, and §5
measures it rather than trusting it.

A plain Cauchy bound would not have this property: dividing by 1501.728 rounds, and it cost 1.7e-11
relative on a true biquadratic in an intermediate version of this work before the `frexp`/`ldexp`
rounding was added.

### 2.2 The `termQ` criterion

```cpp
bool biquadratic = std::abs(termQ) <= kQuarticEpsilon;      // 32 * DBL_EPSILON
```

`q = d − bc/2 + b³/8`, and normalisation has put `|b|, |c|, |d|` in [0, 1], so the three terms sum
to at most 1.625 and the running-error bound of their cancellation is a small multiple of the
machine epsilon with no further factor. The constant is an allowance for that multiple, not a
tuned value.

**A bound taken over q's own terms does not work, and this is a measured correction to the first
version of this fix.** When `b = 0` the terms reduce to `|d|` alone, and `d` is then itself the
residual of a cancellation the *caller* performed, so a floor proportional to `|d|` is proportional
to the very noise it has to dominate. Measured: the quartic with roots {−2, −1, 1, 2}·1e-3 has
`|d| = 2.9e-18` in the normalised variable — a true biquadratic that a term-relative floor would
never recognise, and did not.

### 2.3 The resolvent criterion, and the fall-through

```cpp
const double resolventScale = std::max({|A2|, sqrt(|A1|), cbrt(|A0|)});   // the cubic's own Cauchy bound
if (resolvent > kQuarticEpsilon * resolventScale) { ...Ferrari... } else { biquadratic = true; }
```

Two separate points, and the second is the one that removes the missing walls.

**The test is about resolution, not admissibility.** The resolvent cubic is `−q²/8 < 0` at `m = 0`
and grows, so whenever `q ≠ 0` its largest real root is *mathematically* strictly positive. There is
nothing to admit. What the test has to decide is whether the computed root is distinguishable from
zero, and the scale at which that cubic's roots are resolved is its own Cauchy bound — hence the
factor.

**When it fails, take the biquadratic branch — never return nothing.** A resolvent at noise level
means `q` is at noise level, and a quartic whose `q` is at noise level *is* biquadratic to the
precision available. Returning the empty set, as the absolute guard did, converts an
ill-conditioned crossing into a missing wall, which is the single worst thing a ray/surface solver
can do. The two branches are the two sides of one limit and the code now treats them that way; §4.3
measures that there is no longer a cliff between them.

### 2.4 The Newton polish, with no tolerance at all

```cpp
const double step = quartic(root) / quarticDerivative(root);
if (std::isfinite(step) && std::abs(step) <= 2.) { root -= step; }
```

Every root of the normalised quartic satisfies `|z| ≤ 2` by the Cauchy bound, so a step longer than
that is meaningless whatever produced it, and a zero derivative makes the step non-finite and is
rejected by the same test. The old guard existed only to avoid dividing by ~0; this is the same
protection expressed structurally, and it is strictly safer — the old form allowed an arbitrarily
long step from a merely *small* derivative.

### 2.5 `solveDepressedCubic` loses its tolerance too

Its opener was `|P| <= kTolerance && |Q| <= kTolerance -> push 0`, whose only real job was to keep
`P = Q = 0` out of the `0/0` in the trigonometric branch. The branches now split on an **exact
structural condition**: `P >= 0` forces the discriminant `Q²/4 + P³/27 >= 0`, so Cardano always
applies there; `P < 0` makes the trigonometric magnitude `2√(−P/3)` strictly positive, so that form
is always well defined. No gap, no threshold, and `P = Q = 0` lands in Cardano and returns the same
0 it used to.

### 2.6 The branch is now observable

`solveQuarticReal` takes an optional `QuarticBranch*` out-parameter recording which branch ran. The
kernel never reads it and every existing call site is unchanged. It exists because "it returned the
right roots" and "it took the right branch" are different statements on disjoint inputs, and the
two-directional positive control of §3 needs the second one directly rather than inferred. It went
in with the *red* commit, before any behaviour changed, so the pre-fix branch decisions are on the
record too.

---

## 3. The two-directional positive control

**A guard that always passes is not a fix.** The biquadratic branch is genuinely correct when `q`
really is zero, and a non-positive resolvent genuinely means something. Both directions are
asserted, and each catches the mutation the other cannot.

| direction | what is asserted | which mutation it catches |
| --- | --- | --- |
| **still selects** | `x⁴ − 5x² + 4` and the shifted/scaled family `{c−2, c−1, c+1, c+2}·k` for `c ∈ {0, 1, 100}`, `k ∈ {1, 1e-3, 1e3}` report `Biquadratic` **and** return the four exact roots | a `termQ` floor narrowed away: everything goes down the resolvent branch, whose resolvent is then ~0 |
| **still declines** | `a4 = 0` → 0 roots, `NotAQuartic`; `(x²+k²)(x²+4k²)` → 0 roots, `Biquadratic`; `(x²+k²)((x+k)²+4k²)` → 0 roots, `Resolvent`; all at `k ∈ {1e4, 1, 1e-4, 1e-8}` | a guard widened into irrelevance, and a branch that stops being able to say "no real roots" |
| **and the reverse of "still selects"** | `{1, 2, 3, 7}·k` reports `Resolvent` at every `k` from 1e6 to 1e-8 | a `termQ` floor widened: an asymmetric quartic misrouted into the branch that assumes symmetry |

Measured with the shipped code, all green:

```
[3] x^4 - 5x^2 + 4                     n=4 : -2 -1 1 2            residual 0
    biquadratic about x=100            n=4 : 98 99 101 102        residual 0
[4] x^4 + 1                            n=0
    x^4 + x^2 + 1                      n=0
    (x^2+1)(x^2+2x+5)                  n=0
    (x-1)^4                            n=4 : 1 1 1 1              residual 0
    (x^2-1)^2                          n=4 : -1 -1 1 1            residual 0
```

The tangency convention survives: `(x−1)(x−1−ε)(x²+1)` returns a **pair** at every ε from 0 to
1e-3, so the parity clustering still sees an even-sized cluster and drops it.

---

## 4. What the solver now does, measured

`residual` throughout is the *relative backward error* `|p(x)| / Σ|aᵢxⁱ|` — scale-free, so it means
the same thing for a torus 400 cm away and one 0.1 cm across.

### 4.1 The ALICE3 production-scale reproducer

```
solveQuarticReal(1.0, -1501.7280000044018, 845808.25396968238,
                 -211752288.545858, 19882619385.616932)
```

| | roots | residual |
| --- | --- | --- |
| before | **none** | — |
| after | 375.339229707, 375.524770358 | **1.2e-17** |

The reference roots are Newton run to convergence on the *exact binary values* of those five
coefficients in 60-digit decimal arithmetic: **375.33922957799471** and **375.52477042409094**. The
returned roots are wrong by 1.3e-07 and 6.6e-08 cm — **below one ulp of the input**: `p'` at the
first root is −14.47, so one ulp of `a0` (3.8e-06 at 1.99e+10) moves it by 2.6e-07 cm. The
coefficients do not determine the roots better than that, and the test's 1e-06 cm band is a few
times the conditioning limit and four orders below the 0.1 cm tube whose crossing this is.

*(`Tutorial.md` §5.5 and `Stream_L_ALICE3Defect.md` §3.3 quote the true roots as 375.3392298 and
375.5247703. Those are within 2.2e-07 of the exact values above, i.e. also inside the conditioning
limit — not an error, but they should not be treated as more precise than they are.)*

### 4.2 The asymmetric quartic {1, 2, 3, 7}, swept over 14 decades

| k | before | after |
| --- | --- | --- |
| 1e6 … 1e-2 | 4 roots, residual ≤ 7.4e-17 | 4 roots, residual ≤ 7.2e-17 |
| **1e-3** | 4 roots, **residual 1.07e-06** (roots wrong in the 5th digit) | 4 roots, residual 5.4e-17 |
| **1e-4 … 1e-8** | **2 roots** — 6.5093e-04·(k/1e-4) and −9.3257e-07·(k/1e-4), residual **1** | 4 roots, residual ≤ 7.9e-17 |

The two wrong roots at 1e-4 are exactly the "6.5093 and −0.0093" of the brief, scaled. They are
symmetric about `−b/4` by construction, which is the signature of the misrouting.

### 4.3 There is no cliff any more

The defect is a *cliff*: an absolute threshold crossed by a quantity that carries units. The repair
has to be continuous instead. `{−2, −1, 1, 2 + δ}` walks a quartic away from exactly biquadratic:

| δ, at scale 1e-4 | before, worst residual | after |
| --- | --- | --- |
| 1e-1 | **1.56e-02** | 9.4e-17 |
| 1e-3 | **1.50e-04** | 4.9e-17 |
| 1e-6 | **1.50e-07** | 4.9e-17 |
| 1e-9 | **1.50e-10** | 9.9e-17 |
| 1e-12 | 1.5e-13 | 2.5e-17 |
| 0 (exactly biquadratic) | 2.0e-16 | 4.9e-17 |

At scale 1 both are exact at every δ; the whole column is a scale effect and it is gone.

### 4.4 The Stream E fixture ray

| global scale | before | after |
| --- | --- | --- |
| 1, 0.15 | 2 roots | 2 roots, identical |
| **0.12, 0.1, 0.05, 0.01** | **0 roots** | 2 roots, exactly the ×1 roots times the scale |

`StreamE_TorusQuarticLosesEveryRootBelowTheResolventGuard` asserted 2 / 0 / 0 and carried a comment
saying that failing after the repair was the intended behaviour and that the case must be updated to
say 2 / 2 / 2. It is renamed `...KeepsEveryRootBelowTheOldResolventGuard`, updated, and given a
scale-covariance check over four scales.

---

## 5. Inertness — measured, because this touches every torus intersection in the project

### 5.1 Both gates, totals and disagreement counts together

| | before | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 86 cases, green | **91 cases, green** (5 appended; 1 existing case rewritten per its own instruction) |
| fixtures gate, shipped | 9/9 | **9/9** |
| fixtures gate, surface (historical) | 9/9 | **9/9** |
| fixtures scored | 9 of 10 leaf solids | **9 of 10** |
| Bagger gate, shipped | 12/12 | **12/12** |
| Bagger gate, surface (historical) | 9/12 | **9/12** |
| Bagger scored | 12 of 13 leaf solids | **12 of 13** |
| unexplained oracle disagreements, **surface**, fixtures | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| unexplained oracle disagreements, **surface**, Bagger | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| disagreements, **shape**, fixtures (2 parts) / Bagger (7 parts) | 0/0/0/0 | **0/0/0/0** |
| disagreements, **mesh**, fixtures | 283 / 6936 / 5504 / 5561 (2/9 clean) | **identical** |
| disagreements, **mesh**, Bagger | 418 / 6721 / 7973 / 10299 (0/12 clean) | **identical** |
| `runOracleGate.py --self-test` | 17/17 | **17/17** |
| `o2-bench-detectorsbase-xray --self-test` | 17 checks, 0 failures | **17 checks, 0 failures** |

### 5.2 The column diff — which fields moved, and which did not

`gate.json` compared leaf by leaf with `timing*`, `*Seconds`, `nsPerCall`, `checksum` and absolute
workdir paths stripped:

| corpus | leaf fields | removed | **changed** | added |
| --- | ---: | ---: | ---: | ---: |
| fixtures (9 parts) | 8487 | 0 | **9** | 0 |
| Bagger (12 parts) | 15140 | 0 | **0** | 0 |

**Bagger is bit-identical**, including its 2 toroidal faces of 288 — no ray of its 5000-point sample
set is conditioned badly enough for a branch to move.

All nine fixture changes are on **one part, `torus_union_cyl`**, and all nine are inside its
**mesh-reference** block `distout.validation`:

```
~ [7].distout.validation.worstDeviation           0.018052601964193116 -> 0.01805260196419267
~ [7].distout.validation.worstOffenders[0..9].candidateValue / .deviation   (4 offenders)
```

They are the last two or three digits — 2.5e-14 relative. **Not one count, classification or verdict
moved**, and every `oracle` column on that part reads 0 disagreements with `worstDeviation` exactly
`0.0` before *and* after. The residual movement is the expected footprint of the change: on four
rays the branch decision itself moved (a `q` now recognised as noise), so Ferrari's differently but
equally validly rounded root comes out.

*A comparison that cannot fail is not evidence.* The differ was given exactly one nudged leaf
(`nMismatchMissedSurface` 0 → 1) and reported exactly that one, on 8487 fields.

> **`compareGateRuns.py` could not be used and is broken on this branch.** It raises
> `TypeError: unsupported operand type(s) for *: 'NoneType' and 'float'` at line 135 on this
> `gate.json` — a null leaf (the CSG parts' `capacityRelativeDeviation`) reaches
> `b * factor ** exponent`. It fails its own `--self-test` the same way, so the failure is in the
> tool, not in the inputs. The gate scripts are not this stream's to edit; the diff above was done
> with a throwaway differ carrying its own positive control. **Reported, not fixed.**

### 5.3 The §4.2 direction-independence sweep

Root finding decides crossing parity, so the standing rule applies. 11000 bbox-spread points per
part, 13 golden-spiral directions, `Contains` against `ContainsAlongDirection`, every `Reliable`
part of both gated corpora — and, as the positive control, the ×0.1 fixture ladder where the defect
is known to be live.

| corpus | points × directions | before | after |
| --- | --- | ---: | ---: |
| fixtures + Bagger, 21 parts, all `Reliable` | 231000 × 13 | **1** split | **1** split — *the same point, same part* |
| ×0.1 fixtures, 9 parts, all `Reliable` | 99000 × 13 | **32** splits | **29** splits |
| — of which `torus_union_cyl` | | **3** | **0** |
| — `cyl_cross_cyl` / `cyl_inter_cyl` / `tube_window` | | 10 / 14 / 5 | 10 / 14 / 5 — unchanged |

This is the sweep's self-check as much as its result: an instrument that reported "no change" here
without being able to see the defect would prove nothing, and it *does* see it — 3 direction-split
points on the ×0.1 torus fixture before, 0 after, with every other part's count identical to the
unit. Nothing moved onto or off `Contains`'s single-shot fast path: all 21 parts are `Reliable`
before and after.

The one split on the gated corpora is **identical before and after** and is not this stream's:
`cyl_inter_cyl`, point `(-0.49825084855313678, 0.43950780674062478, -0.4362455074381224)`,
`Contains` = inside, `Safety` = 0.336 cm. It is deep inside the solid rather than near a surface, so
it is not the trim-sliver family of `TolerancePolicy.md` §11; see §7.

---

## 6. End to end — the three places the defect was observed

### 6.1 `torus_union_cyl` at ×0.1, `--beams 96 --raster 24`

| | before | after |
| --- | --- | --- |
| rays identical, mode (a) | 55294 / 55296 | **55296 / 55296** |
| **LOST** | **4** (in pairs: enter *and* exit) | **0** |
| **parity** (`parityNB` 0 both) | **2** | **0** |
| extra / displaced / kind | 0 / 0 / 0 | 0 / 0 / 0 |
| zero / non-advancing / stalled / capped steps | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| mode (a) vs mode (b) | 0 of 55296 | 0 of 55296 |

### 6.2 ALICE3, `--raster 16 --representations surface`

| | before | after |
| --- | --- | --- |
| `surface` (a) shape | 13814 / 13822 rays identical | **13822 / 13822** |
| `surface` (b) nav | 13816 / 13822 | **13822 / 13822** |
| **parts fully clean** | **17 / 18** | **18 / 18** |
| **LOST** | **14** | **0** |
| **unterminated** | **2** | **0** |
| extra / displaced / kind / parity | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| mode (a) vs mode (b) | **6** of 13824 rays | **0** of 13824 |
| zero / non-advancing / stalled / capped steps | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |

`Stream_L_ALICE3Defect.md` §2.5 left exactly this residual after fixing mechanism 1 — "14 LOST, 2
unterminated, 6 mode-(a)-vs-(b) disagreements, all on `ST2487462_01`" — and it is gone. **ALICE3 is
now 18 of 18 transport-clean at this density.** The `--raster 16` caveat of Stream J §7 still
applies: three axis beams is direction-poor, so read this as "no defect at this density".

### 6.3 The shipped ×1 ladder parity audit — and a trap that nearly cost the measurement

| | before | after |
| --- | --- | --- |
| `--parts torus --beams 256 --raster 32`, **all representations** | 262144/262144 rays, LOST 0, **parity 2** | 262144/262144, LOST 0, **parity 0** |
| the same with **`--representations surface`** | LOST 0, **parity 0** | LOST 0, parity 0 |

**Read the second row before the first.** Adding `--representations surface` makes this check
structurally incapable of reproducing the event it exists to detect. The raster window is derived
from the tightest bounding box available, and with no mesh converted it comes from the surface
solid's conservative box instead — a different window is a different ray lattice (563688 steps
versus 563732), and a 2-in-262144 event that one lattice samples the other does not. The first run
of this measurement reported parity 0 *before* the fix and would have been written up as "the ×1
parity mismatch is not reproducible on this branch". It is reproducible; the experiment was not
capable. This is `Tutorial.md` §4.5 again, in the same benchmark that taught it (Stream J §6.1), one
flag further down.

---

## 7. Found, not fixed

1. **`compareGateRuns.py` is broken** — `TypeError` on a null leaf, including in its own
   `--self-test`. §5.2. It is the project's declared instrument for "diff columns, not verdicts",
   so this matters more than its size suggests.
2. **The depression step is the real conditioning ceiling, and it is unchanged by this work.**
   `p, q, r` are formed from `b, c, d, e` by cancelling numbers of size `centreᵏ` to leave numbers
   of size `spreadᵏ`, so a quartic whose roots agree to N significant figures has lost ~2N digits
   before any branch is chosen. Measured on `{c−2, c−1, c+1, c+2}` with rounded coefficients,
   **identically before and after**:

   | relative root spread | worst relative root error |
   | --- | --- |
   | 2e-01 | 3.9e-14 |
   | 2e-02 | 6.2e-11 |
   | 2e-03 | 9.3e-08 |
   | **2e-04** | **no roots returned at all** |

   This is why the §3 "still selects" family stops at centre 100. Removing it means abandoning
   Ferrari for a companion-matrix or Jenkins–Traub solve, which is a different piece of work with a
   different cost. For a torus it bounds how close two crossings may be relative to the ray's lever
   arm before the solve stops resolving them.

3. **A double-double tangency can lose both pairs.** `((x−k)²(x−3k)²)` at `k = 1e-6` returns 0 roots
   before *and* after, because the biquadratic discriminant `p² − 4r` is exactly 0 in exact
   arithmetic and the input coefficients' own rounding takes it negative. Benign for parity — an
   even cluster would have been dropped as a tangency anyway — but it is a real zero-vs-four
   discontinuity and it is not this stream's.
4. **One direction-split point survives on the gated corpora**, identical before and after:
   `cyl_inter_cyl`, `Safety` = 0.336 cm, `Contains` = inside, 1 of 13 directions. It is *not* near a
   surface, so the trim-sliver diagnosis of `TolerancePolicy.md` §11 does not obviously apply, and it
   is a different part and a different point from the `cyl_cross_cyl` offender §10.1 recorded. The
   earlier sweeps used a different point generator, so 1-in-231000 here and 0-in-154000 there are
   different samples of the same rate, not a before/after.
5. **`sameIntersection` is still `|t1 − t2| ≤ 1e-7 · max(1, |t1|, |t2|)`.** The `max(1., …)` stops it
   being relative below 1 cm and makes it an absolute 1e-7 cm band — the same class of defect as
   this one, on the same code path, flagged by `Stream_E_Scale.md` §5.3 and still not implicated by
   any measurement. Out of scope here; it should be looked at by whoever next touches parity.
6. **The ×0.1 ladder still has 29 direction-split points** on `cyl_cross_cyl`, `cyl_inter_cyl` and
   `tube_window`, unchanged by this work. They are the K5 trim-sliver family: `kBSplineFlatness` is
   an absolute 1e-5, so at ×0.1 the sliver it licenses is ten times larger relative to the part.
   That is the *next* dimensionally-inconsistent constant in the same file, and the ×0.1 sweep is
   its reproducer.

---

## 8. Reproducing all of it

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH   # stage FIRST
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid \
      o2-bench-detectorsbase-solid-harness o2-bench-detectorsbase-xray
ctest -R BVHSurfaceSolid                                     # 91 cases

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate  --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step

# the three end-to-end measurements of section 6 -- note the THIRD one takes no
# --representations flag, on purpose (section 6.3)
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xt --fixtures \
    --transform scale:0.1 --parts torus --beams 96 --raster 24 --representations surface
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xa \
    --model scripts/geometry/ALICE_3_example/CAD_noETA.stp --raster 16 --representations surface
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xp --fixtures \
    --parts torus --beams 256 --raster 32
```

The fastest loop by far, and how every number in §4 was produced, is a standalone `.cxx` calling
`solveQuarticReal` directly — no CAD, no ROOT geometry, two minutes to write
(`TolerancePolicy.md` §6):

```bash
g++ -std=c++20 -O2 -w -o /tmp/probe /tmp/probe.cxx \
    -I$B/stage/include -I$HOME/alisw/O2/Detectors/Base/include \
    -I$HOME/alisw/O2/Detectors/Base/src $(root-config --cflags --libs) \
    -L$B/stage/lib -lO2DetectorsBase
```

The §5.3 sweep is the same recipe plus `-lGeom`, walking every `surfaces_*.bin` under a gate
workdir's `db/`, calling `LoadSurfaceSolid` + `CloseShape(false)` and comparing `Contains` against
`ContainsAlongDirection` over 13 golden-spiral directions at 11000 deterministic bbox-spread points.
Point it at a ×0.1 conversion first: if it does not report 3 splits on `torus_union_cyl` against the
pre-fix library, it is not measuring what it claims to.
