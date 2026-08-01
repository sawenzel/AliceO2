# Stream C (wave 0) — kernel hygiene: the capacity residual, and two latent defects

Date: 2026-08-01. Branch `swenzel/bvhsurfacesolid`, on top of `2c25a1a20a`.
Brief: [`Workstreams.md`](Workstreams.md) Stream C. Findings: [`CodeReview_Fable_v2.md`](CodeReview_Fable_v2.md) §3.

**Headline: the capacity column is solved.** The residual that survived five sessions and had been
declared "unexplained" is a quadrature defect after all, and it was two defects wearing each other
as camouflage.

| | before | after | OCCT |
| --- | --- | --- | --- |
| `BoomCylinderInner` | 4.88e-04 | **2.8e-07** | — |
| `StickCylinderInner` | 6.87e-04 | **passes** | — |
| `BucketCylinderOuter` | 1.22e-04 | **passes** | — |
| `BoomCylinderOuter` | 1.08e-04 | 1.39e-06 | — |
| `BucketCylinderInner` | 3.04e-04 | 1.31e-06 | — |
| `StickCylinderOuter` | 3.78e-04 | 1.39e-06 | — |
| `tube_window` (fixture) | 3.02e-05 | **passes** | — |

Three of the six Bagger cylinders and the one failing fixture now pass their capacity column
outright; the other three sit at 1.3e-6 against a 1e-6 band, a 100-500x improvement, and that
residual is **not** quadrature (see §4).

## 1. How it was found — the instrument, not the hypothesis

`Capacity()` is one number for a whole solid, so it can say *that* a volume is wrong and never
*where*. Every previous attempt reasoned about mechanisms; this one built the localiser first.

- `O2BVHSurfaceSolid::GetSurfaceCapacityContributions()` — each face's own divergence-theorem term,
  in `GetSurfaceRecords()` order. New public diagnostic.
- OpenCascade answers the same question per face with `BRepGProp::VolumeProperties` on a single
  face. **Caveat found by self-check and worth recording: it returns exactly 0 for a planar face.**
  On `BasePin` the curved term is 20.943951 and the volume 31.415927, and the missing 10.471976 is
  precisely the two end caps' analytic contribution — so the identity holds and the comparison is
  meaningful on quadrics, with planes checked by difference against the total.

On `BoomCylinderInner` that localised the entire 1.1e-2 cm³ error to **one face** — the fat tube's
wall, the only face carrying a hole — in a single run. Faces 0 and 5 matched OpenCascade to 7e-6
and to every printed digit respectively.

## 2. The defect: Gauss-Legendre across knot spans

`contourIntegralAlongCurve` spent one 20-node Gauss-Legendre rule across a B-spline's *entire* knot
domain. Its comment justified the blanket order with "the antiderivatives the quadric surfaces
supply are entire in u" — true, and about the wrong variable. The integration variable is `t`, and
`F(u(t), v(t)) · v'(t)` on a B-spline is analytic only *within* a knot span. A 179-pole cubic has
~176 spans and was given 8 intervals. `Curve2D::signedAreaContribution`, 1000 lines above in the
same header, has integrated the same curve family span by span since it was written.

Fix: subdivide at the curve's interior knots (`Curve2D::appendInteriorKnots`) before applying the
interval cap.

## 3. Why it hid: a second defect made the first untestable

The obvious experiment — shrink `kContourMaxSpanU` and watch the error fall — was run first, at
`π/4 → π/256`, a 64x finer split. **Nothing moved beyond the last printed digit**, and on that
evidence the quadrature hypothesis was written off (`CodeReview_Fable_v2.md` §3/N1 records the
refutation as it stood).

The sweep was blind by construction. `spanU` was the difference between the curve's *endpoints*,
and the face carrying the whole error is trimmed by a **closed** loop, whose endpoints coincide.
So `spanU = 0`, and `max(1, ceil(0 / x))` is `1` **at every value of x** — the one curve that
needed splitting was the one curve the knob could not reach.

Fix: `Curve2D::uVariation`, the travel of `u` along the curve — exact for a line, exact for an arc
via its turning points at angle 0 and π, and bounded by the per-span pole hull for a B-spline.

**Attributed, because two fixes landed together and only one of them is the cause**
(`BoomCylinderInner`, OCCT 22.989411):

| configuration | `Capacity()` |
| --- | --- |
| baseline | 23.000623 |
| `uVariation` only | 23.000579 |
| knot subdivision only | (test fails, see below) |
| both | **22.989418** |

So the knot subdivision is the fix and `uVariation` is worth 4e-5 cm³ here. `uVariation` stays
because it is correct and because it is what made the first experiment lie — but it is a
correctness fix, not the cause, and this document says so rather than letting the commit imply
otherwise.

## 4. What the remaining 1.3e-6 is not

Doubling `kContourQuadratureOrder` from 20 to 40 leaves the three residuals bit-identical
(1.39e-06, 1.31e-06, 1.39e-06). **The quadrature is converged**; what is left is the trim data —
the pcurve as a statement about where the BREP's face actually ends. That is the same object
Stream F (edge identity) and the CSG route address, and it is a reasonable floor to stop at:
250-500x better than the start, and below the threshold at which anything else on the board can
be measured.

## 5. The other two findings

- **N2** — the loader read a B-spline edge's endpoints off its first and last poles, which are the
  endpoints only for a *clamped* knot vector, while the kernel has evaluated them since K1. On an
  unclamped input the two layers measured the same wire join between different points and could
  disagree in both directions. `edgeEndpoints` now evaluates. Latent on this corpus — the
  converter's canonical recognition converts such curves away before the loader sees them — and
  pinned by a test instead.
- **N4** — `containsByVote` counted a shot whose parity rested on a trim-boundary tie-break exactly
  like one that did not, so a point could be decided three-to-two by three coin flips and report
  nothing. It now decides among the shots that stand on the geometry whenever those are not
  themselves tied. Inert on this corpus (every gate column bit-identical), which is the expected
  outcome given §16's measured fire rate of one shot in 104000.
- **N6** — the 38 untracked entries in `scripts/geometry/` turned out to contain **nothing
  tracked**: there was no `git rm` to do, contrary to the brief. A `.gitignore` removes the editor
  litter (`*~`, `#..#`, `._*`, `a.out`, `*.o`, `__pycache__`) from `git status` without deleting
  anything. The stale `O2_TGeoToCAD*.py` variants are deliberately left visible: Stream E is about
  to revive one for the TGeo → STEP round trip, and choosing which needs a human.

## 6. Acceptance

- `ctest -R BVHSurfaceSolid` green, **64 cases** (from 61).
- Oracle gate: fixtures **8/9**, Bagger **6/12** — totals unchanged, because every part that still
  fails does so on navigability, which is Stream F's subject and which this stream does not touch.
- **Every gate column bit-identical except `oracle`**, and that only on the 6 Bagger parts and 3
  fixtures with B-spline trims: checksums, `contains`/`distout`/`distin`/`safety` validation,
  `navigation.rimDetail` and the offender lists all match the baseline exactly. Diffed with the
  timing fields stripped, per `NEXT.md`.
- The capacity regression test was verified to **fail** (by 8e-4 absolute) with the knot
  subdivision removed.

## 7. Two things for whoever reads the numbers next

- **`BRepGProp::VolumeProperties` on a single planar face returns 0.** Any per-face comparison must
  either restrict itself to curved faces or account for the planes by difference. This cost a
  confused half hour and is not documented anywhere in OCCT's own text.
- **A refuted experiment is not a refuted hypothesis.** The `kContourMaxSpanU` sweep was a sound
  test of an unsound assumption — that the knob reached the curve. When an experiment says "no",
  check that it was capable of saying "yes"; here it never could have.
