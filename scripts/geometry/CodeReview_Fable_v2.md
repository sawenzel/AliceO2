# BVHSurfaceSolid — second deep review: code, feasibility, scale-up, and the CSG alternative

Date: 2026-08-01. Reviewer: Claude (Fable 5), on branch `swenzel/bvhsurfacesolid` at `2c25a1a20a`.
Method: full read of [`CodeReview_Fable.md`](CodeReview_Fable.md) (v1, sections 1-20),
[`NEXT.md`](NEXT.md), the section index of [`TolerancePolicy.md`](TolerancePolicy.md), the five
commits since v1 was opened, and targeted reading of the code they touched — the rim machinery and
the Green's-theorem capacity in `BoundedSurface.h`, the re-shoot policy in `O2BVHSurfaceSolid.cxx`,
the sidecar writer in `O2_CADtoTGeo.py`, the loader in `O2SurfaceSolidIO.cxx`, and the gate scripts.

**The review was written from reading; §3/N1 was then measured and acted on**, in a second pass on
the same day, after this document's own first diagnosis was put to the test — baseline gates, a
sidecar census, a 64× quadrature sweep, a chart-contour dump, a Monte-Carlo volume witness, and
finally a per-face localiser that found it. That work is Stream C of
[`Workstreams.md`](Workstreams.md) and its record is [`Stream_C_Hygiene.md`](Stream_C_Hygiene.md);
§3/N1 below keeps the intermediate refutation in place, because how it was wrong transfers and the
conclusion alone does not. Everything else here is unmeasured: where a claim is a prediction it says
so and gives the experiment that would falsify it. This document supersedes
v1's *plan* (sections 9-10); v1's *findings* remain the register of record and are referenced by
their labels (K*, S*, C*).

This document answers the questions posed for this session:

1. What is the state of the code after the last five commits, and what is wrong with the new material?
2. Are the goals still feasible, and are they still the right goals?
3. What are the next steps toward large-scale testing (the ALICE beampipe)?
4. The orthogonal proposal: a simplification pipeline that recognises and emits **boolean solids**
   — could tube-tube be handled precisely and more cheaply that way? And extrusions?

---

## 1. Verdict

**The navigation problem is solved and the representation problem is not, and the project should
now stop trying to solve the second one inside a single trimmed-surface solid.**

That is a stronger statement than v1 could make, and the branch earned it: over both corpora — 9
converted fixtures and 12 Bagger parts, 5000 points and 5000 rays each — **every disagreement with
the OpenCascade oracle outside tolerance, on every column, is zero.** `Contains`,
`DistFromOutside`, `DistFromInside` and `Safety` agree with a real CAD kernel on real CAD, including
on the seven parts the solid itself declares unreliable. That was 35 disagreements two sessions ago
and 4588 "unexplained" points a week ago. The kernel, the BVH, the parity rule, the tie-break
policy and the oracle gate are, as an assembly, doing their job.

What is left over is entirely bookkeeping about the *boundary representation*: a residual open rim
of 0.3-0.9 cm per part on six Bagger cylinders, and a capacity error of 1e-4 to 7e-4 relative on
those same parts. Neither is a wrong answer to a navigation query. Both come from the same root: two
faces meeting at a curve that **neither of their charts can express**, so each carries its own
approximation of it and the two are compared with a band that has to be inflated by the sampling
error of the comparison itself. `NEXT.md` reaches the honest end of that road: "the sagitta bounds
each polyline against *its own* curve and says nothing about the other face's… do not paper over
that with a fudge factor."

Three things follow, and they are this review's recommendations:

- **The capacity residual is solved** (§3/N1, and [`Stream_C_Hygiene.md`](Stream_C_Hygiene.md)).
  It fails on exactly those parts carrying a trim B-spline with more than one interior knot span —
  21 of 21 parts, both corpora — and the cause is Gauss-Legendre spanning knots, hidden behind a
  second defect that made the obvious experiment incapable of detecting it. Fixed: 4.9e-4 → 2.8e-7
  on the worst part, three of six Bagger cylinders and the failing fixture now pass their capacity
  column outright, and all three of `NEXT.md` thread 2's candidates are dead. §3/N1 keeps the
  refuted intermediate diagnosis, because the way it was wrong is the transferable part.
- **Closure should be decided by topology, not by proximity.** The converter knows, per trim curve,
  which `TopoDS_Edge` it came from and therefore which other face shares it. It throws that away.
  Carrying an edge ID in the sidecar turns "is this rim matched" from a geometric search with a
  self-referential tolerance into a structural check that cannot be wrong, and demotes the geometric
  disagreement from a *verdict* to a *reported number*. This is ~40 lines each side and it dissolves
  the band circularity permanently (Section 3.3).
- **The user's orthogonal thought is right, and it is stronger than they put it.** Emitting boolean
  solids is not merely a cheaper special case — for solids bounded by planes, cylinders, cones,
  spheres and tori it is a *general, exact, watertight-by-construction* alternative to the entire
  trim/pcurve/closure edifice, with two decades of production precedent in the fusion-neutronics
  community (McCAD, GEOUNED). It subsumes v1's Phase 2 (adjacency trims) for every case Phase 2 was
  designed for, it is a better fit for the AOT-codegen direction already chosen, and it fails
  gracefully — non-convertible parts keep the existing solid. Section 5 evaluates it properly,
  including where it does *not* work and what it costs. **My recommendation is to do it instead of
  Phase 2, and to keep BVHSurfaceSolid as the fallback leaf.**

On the beampipe: the corpus for it does not exist in this repository and acquiring it is on the
critical path. Section 6 proposes getting a ground-truth one for free by round-tripping O2's *own*
`Pipe.cxx`/`PipeRun4.cxx` TGeo beampipe out to STEP and back — which gives something neither Bagger
nor ALICE3 can, an oracle that is exact rather than tolerant, and a scoring answer for the
recognition pipeline that is known in advance.

---

## 2. What the last five commits did, and how they hold up

| commit | subject | verdict |
| --- | --- | --- |
| `d60584b` | oracle categorises the gate's own rays | **Sound, and important.** The rule never consults the candidate, which is exactly what makes a gate change legitimate. It removed 78 phantom disagreements on `BucketLink2` alone. The `nRelabelled` counter on the harness line is the right hygiene. |
| `01f64b4` | capacity by Green's theorem on wire-trimmed quadrics | **Right idea, correct reduction, defective quadrature — confirmed by measurement and fixed in Stream C.** The four antiderivatives are right and the `deltaV != 0` bridge skip is exactly right (`∮F dv` over a `dv = 0` segment is identically zero). But the integrator spanned knots, and its interval count came from the curve's endpoints — which is zero for every closed hole loop, so the cap could not reach the curve that carried the error. Together those were the whole of the residual capacity column. §3/N1. |
| `bcc220d` | scope item 4, correct the premise | **Model behaviour.** Three measurements killed the premise the work was scoped on (no face is missing; the two trims *are* the same curve; it is not a matching tolerance) and the session reported that instead of proceeding. `maxGap` being independent of `rimEpsilon` "by construction" is the kind of tell that is easy to look past. |
| `80edd15` | name the rims, find the curve that was never sampled | **The best commit on the branch.** Instrumenting first (per-rim records) and reading the answer rather than inferring it found K4 — a 179-pole junction spline replaced by a straight line by a single-midpoint flatness probe — which was v1's oldest finding with no reproducer. The `RimState`↔`NavigationReliability` identity as a self-check is a good pattern. |
| `2c25a1a` | record items 1-3, hand over | Accurate. `NEXT.md`'s "six things that will not be obvious" is the most useful page in the repository. |

Two process observations worth keeping:

- **"A bit-identical numeric column can accompany a real change, and did"** (`NEXT.md`) is the
  single most valuable lesson of the last session and it generalises: this project's gate totals are
  a lagging, low-resolution indicator, and three of the last four real advances were invisible in
  them. Keep diffing columns.
- The K4 fix arrived through an instrument built for a different purpose. That is the third time on
  this branch (S6's crossing dump, the oracle itself, now the rim records). The pattern —
  *build the instrument that names the object, then read it* — is worth making explicit policy.

---

## 3. Findings — new material

Labels continue v1's series. **N** = new this review.

### N1. The capacity residual — a perfect predictor, four dead hypotheses, and one instrument left

**This section records a diagnosis of mine that the measurement refuted.** It is kept in full,
because the refutation is the useful part and because the eliminated hypotheses are the ones the
next session would otherwise spend a week on. Everything below was measured on 2026-08-01 against
the branch state, baseline gates re-run first (fixtures 8/9, Bagger 6/12, matching `NEXT.md`).

#### N1.1 The predictor is exact: a trim B-spline with more than one interior knot span

Parsing the sidecars directly (no kernel, no OCCT) for the degree, pole count and distinct interior
knots of every B-spline trim curve:

| corpus | part | interior knot spans (max) | capacity relDev | gate |
| --- | --- | --- | --- | --- |
| fixtures | `box`, `box_union_box`, `box_minus_cyl`, `cyl_plus_cone`, `sphere_minus_cyl`, `torus_union_cyl` | no B-spline trims | ≤1e-13 | PASS |
| fixtures | `cyl_cross_cyl` | **1** (degree 8, 16 poles) | 1.8e-12 | PASS |
| fixtures | `cyl_inter_cyl` | **1** (degree 8, 16 poles) | 1.0e-11 | PASS |
| fixtures | `tube_window` | **21** (degree 7, 134 poles) | 3.02e-05 | FAIL |
| Bagger | `Base`, `BasePin`, `Boom`, `Stick`, `BucketLink1`, `BucketLink2` | no B-spline trims | ≤1e-11 | PASS |
| Bagger | `BoomCylinderInner`, `BoomCylinderOuter`, `StickCylinderInner`, `StickCylinderOuter` | **175** (degree 3, 179 poles) | 1.1e-04 … 6.9e-04 | FAIL |
| Bagger | `BucketCylinderInner` | **153** | 3.04e-04 | FAIL |
| Bagger | `BucketCylinderOuter` | **197** | 1.22e-04 | FAIL |

**21 of 21 parts, no exceptions.** Every part whose trim wires carry only lines, arcs, or a
single-span (Bézier) spline is exact to machine precision; every part carrying a many-span spline
deviates, and only on capacity. This is a much sharper statement than "the six cylinders", and it
holds across two corpora built by different means.

#### N1.2 The obvious mechanism — and it is wrong

`contourIntegralAlongCurve` (`BoundedSurface.h:2032`) splits a trim curve into
`ceil(spanU / (π/4))` intervals and puts 20 Gauss-Legendre nodes on each. Its comment justifies the
blanket order by "the antiderivatives the quadric surfaces supply are entire in u". They are — but
the integration variable is `t`, and `F(u(t), v(t)) · v'(t)` on a B-spline is only *piecewise*
smooth: analytic within a knot span, `C^{p-m}` across a knot. Gauss-Legendre's geometric
convergence needs analyticity on the whole interval it covers. A 179-pole cubic has ~176 spans and
gets 8 quadrature intervals. And `Curve2D::signedAreaContribution`, 1000 lines above in the same
header, integrates the very same curve family **span by span** for exactly this reason.

The argument is sound, the asymmetry between the two integrators is real, and the predictor in
N1.1 fits it perfectly. It is still not the cause.

**Measured.** `kContourMaxSpanU` set from `π/4` to `π/256` — a 64× finer split, 512 intervals over
a full turn, comfortably more than the 175 spans — rebuilt, Bagger gate re-run with
`--skip-convert`:

| part | `π/4` | `π/256` |
| --- | --- | --- |
| `BoomCylinderInner` | 4.88e-04 | 4.88e-04 |
| `BoomCylinderOuter` | 1.08e-04 | 1.08e-04 |
| `BucketCylinderInner` | 3.04e-04 | 3.04e-04 |
| `BucketCylinderOuter` | 1.22e-04 | 1.23e-04 |
| `StickCylinderInner` | 6.87e-04 | 6.87e-04 |
| `StickCylinderOuter` | 3.78e-04 | 3.77e-04 |

Nothing moves beyond the last printed digit. **The quadrature is converged**; it is computing the
contour integral of the emitted contour correctly. The constant was reverted; the tree is clean.

(The two secondary observations in that function stand and are still worth fixing on their own
merits: `spanU` is taken from the endpoints only, so it reports **zero** for a closed curve and
gives one interval for a whole loop; and the bridge `Curve2D` is constructed per seam per call.
Neither is the residual — the sweep would have exposed both.)

#### N1.3 The chart-wrapping candidate is dead too, and so is the seam-bridge one

`NEXT.md` thread 2 lists three candidates. All three are now closed.

- **Seam bridges skipped when `deltaV == 0`** — closed on paper, no run needed. The contour
  integral is `∮ F dv`; a bridge with `dv = 0` contributes identically zero whatever its `u` extent.
  The skip is exact.
- **A wire wrapping the periodic chart without its cut segments.** This was the strongest remaining
  candidate and it is refuted by the data. `F` for a cylinder contains the non-periodic term
  `(s r/3)·r·φ`, so a wire that wraps `u` needs the chart cut explicitly or the integral is wrong by
  `(s·2π r²/3)·Δv`. Dumping `BoomCylinderInner` face 0's wire as a chart contour shows the cuts are
  **present as ordinary edges**:

  ```
  edge 0: line                     ( +0.00000, +0.00000) -> ( +0.00000,-16.33000)
  edge 1: bspline(deg=3,poles=179) ( +0.00000,-16.33000) -> ( +6.28319,-16.33000)
  edge 2: line                     ( +6.28319,-16.33000) -> ( +6.28319, +0.00000)
  edge 3: line                     ( +6.28319, +0.00000) -> ( +0.00000, +0.00000)
  LOOP: total du = 0, total dv = 0   (curve sums close to 7e-12, bridges cancel them)
  ```

  A closed, counter-clockwise loop in the unrolled chart with both vertical cuts at `u = 0` and
  `u = 2π` carrying `dv = ∓16.33`. Green's theorem applies as written. Every other face of that part
  is either a full parametric rectangle or a plane with arc wires, all of which check out the same
  way.
- **Contour vs. the retained midpoint integrator** — not run, and now low value: with the quadrature
  shown converged, the midpoint rule's own O(1/N) staircase (recorded oscillating between 22.93 and
  23.03 on this part) cannot resolve a 1.1e-2 cm³ difference.

#### N1.4 The `Contains` witness works but is three times too coarse

The measurement of §13's H2, re-run: Monte-Carlo the solid's volume using `Contains` — which agrees
with OpenCascade on 5000/5000 points of every one of these parts — and compare against `Capacity()`.
6M points each:

| part | MC volume (`Contains`) | `Capacity()` | OCCT | verdict |
| --- | --- | --- | --- | --- |
| `BucketLink2` (control, passes) | 17.06899 ± 0.04230 | 17.07926 | 17.07926 | probe validated |
| `BoomCylinderInner` | 23.00327 ± 0.03352 | 23.00062 | 22.98941 | **inconclusive** |
| `StickCylinderInner` | 23.00681 ± 0.03344 | 23.00520 | 22.98941 | **inconclusive** |
| `BucketCylinderOuter` | 16.00726 ± 0.01730 | 16.01370 | 16.01173 | **inconclusive** |

The discrepancy to be resolved is 1.1e-2 cm³ and the MC uncertainty is 3.4e-2 cm³ — the witness
cannot separate `Capacity()` from OCCT at this precision. Reaching 3σ needs ~500M points per part
(the bbox acceptance is only 7%), which is possible but crude. Recorded so nobody re-runs it
expecting an answer.

#### N1.5 Resolved the same day — the instrument found it, and the refutation was itself wrong

The instrument named in this section was built (`GetSurfaceCapacityContributions()` plus
OpenCascade's per-face `BRepGProp::VolumeProperties`) and it localised the whole 1.1e-2 cm³ error
of `BoomCylinderInner` to **one face** in a single run: the fat tube's wall, the only face carrying
a hole. Faces 0 and 5 matched OpenCascade to 7e-6 and to every printed digit.

**N1.2's mechanism was right and N1.2's experiment was incapable of showing it.** The quadrature
*is* span-blind across knots. The `kContourMaxSpanU` sweep found nothing because `spanU` was the
gap between the curve's *endpoints*, and the face carrying the error is trimmed by a **closed**
loop whose endpoints coincide: `spanU = 0`, and `max(1, ceil(0 / x))` is `1` at every `x`. The one
curve that needed splitting was the one curve the knob could not reach.

Both were fixed — subdivision at interior knots, and `Curve2D::uVariation` for the travel of `u` —
and attributed separately, because two fixes landing together must not be allowed to share credit:
on `BoomCylinderInner` (OCCT 22.989411) the baseline is 23.000623, `uVariation` alone gives
23.000579, and both together give **22.989418**. The knot subdivision is the cause; `uVariation` is
a correctness fix worth 4e-5 cm³ here, and it is what made the first experiment lie.

Doubling the quadrature order to 40 then leaves the three remaining residuals bit-identical, so
what is left (1.3e-6 on three parts) is the trim data, not the integration.

Full record, measurements and acceptance: [`Stream_C_Hygiene.md`](Stream_C_Hygiene.md).

**The lesson worth keeping: a refuted experiment is not a refuted hypothesis.** When an experiment
says "no", check that it was capable of saying "yes".

### N2. The loader still reads B-spline endpoints off the poles — K1's IO twin is open

v1 recorded K1 fixed: "endpoints of an unclamped B-spline are now evaluated instead of read off the
poles." That is true of the kernel. `O2SurfaceSolidIO.cxx:149-157` still does the old thing:

```cpp
if (edge.curveType == kBSpline2D) {
  // for a clamped B-spline the first/last poles are the curve endpoints
  ...
  start = curve.poles.front();
  end = curve.poles.back();
```

and nothing in `parseBSplineEdge` or `Curve2D::valid()` validates clamping — `valid()` checks
`knots.size() == nPoles + degree + 1` and monotonicity, which an unclamped vector satisfies.

Consequence: on an unclamped or periodic input the loader's wire-join check measures a gap between
two points **that are not on the curve**, while the kernel measures it between two points that are.
They can disagree in both directions: a sound wire rejected at load (this is the shape of the
long-standing `ST1829909_01` loader rejection), or a genuinely open wire accepted. It is latent on
the current corpus for the same reason K1 was — the converter's canonical recognition converts most
such curves away before they are written — which is exactly why it will surface on the first model
that is not Bagger.

Fix: three lines, `curve.pointAt(0.)` / `curve.pointAt(1.)`. And add clamping to `valid()` as a
*reported* property rather than a rejection, so the two layers can never diverge silently again.

### N3. Rim matching is geometric where the data is topological, and that is why the band cannot be chosen

`measureRimClosure` decides whether a chord is shared by searching a uniform grid for the nearest
chord of a *different face* and comparing against `epsilon + own sagitta + partner sagitta`. The
implementation is careful — the two-band split (sampling-aware for "shared", declared-tolerance-only
for "non-manifold") is correct and the reasoning for keeping them apart is right. But the whole
construction is answering a topological question with a proximity query, and `NEXT.md` has already
walked into the wall that guarantees:

> the sagitta bounds each polyline against *its own* curve and says nothing about the other face's,
> which is the term that dominates. It works at the shipped tolerance because the two happen to be
> within an order of magnitude — a coincidence of scale, not a criterion.

That is correct and it is not fixable by choosing a better band, because the quantity being
thresholded (how far apart two independent approximations of one curve are) is not the quantity
being decided (whether the two faces share an edge). Tightening the flattening improves the first
and makes the verdict *worse* — §13.8 measured exactly that, isolation 4.8e-5 → 2.6e-6 while the
open length went 0.325 → 0.984 cm.

**The information that answers it is thrown away in the converter.** `_quadric_trim_wire` and
`extract_planar_face` walk `TopoDS_Edge`s; `TopExp.MapShapesAndAncestors` gives the two faces of
each edge in one call (v1 §4.2 already noted it is "on the shelf"); the sidecar writer
(`write_surfaces_bin`, format v2) emits `{curve type, params}` per edge and no identity.

**Recommendation — sidecar v3 carries an edge ID per trim curve** (a `uint32` index into a
per-model edge table, plus a sense bit). Then:

- Closure becomes: every edge ID appears exactly twice, with opposite sense. Zero tolerance, zero
  sampling, zero band. `boundaryRims`, `reversedRims` and `nonManifoldRims` become exact counts.
- The geometric disagreement between the two faces' realisations of one edge becomes a **reported
  measurement** (`maxSharedEdgeDeviation`, in cm, per edge and per solid) rather than a verdict — a
  number that Phase 2 or the CSG route can be scored against, which today does not exist.
- `CloseShape` stops depending on `kBSplineFlatness`, so tightening the flattening is once again a
  pure improvement.
- It costs nothing at run time and shrinks nothing: 5 bytes per trim curve.

This is worth doing **whatever is decided about Section 5**, because it is what makes any claim of
"watertight" checkable rather than tuned, and because the same edge table is the input Phase 2 would
need anyway.

### N4. The vote does not check its own shots

`containsByParity` (`O2BVHSurfaceSolid.cxx:1621-1632`) escalates to `containsByVote` when the single
shot's parity rests on a hit flagged `onTrimBoundary`. Good — and §16's measurement (24 wrong → 0,
0 broken, fires on 1 shot in 104000) is convincing. But `containsByVote` calls
`parityAlong(point, direction, useBVH)` with the ambiguity pointer **omitted** (line 633), so the
five voting shots are counted equally whether or not each of them also rests on a tie-break. A point
where three of five directions are ambiguous is decided by a majority of coin flips, and reports
nothing.

Cheap and strictly better: collect `(answer, ambiguous)` per direction; decide the majority among
unambiguous shots if there is one; only fall back to all five if every shot is ambiguous — and in
that case the point is genuinely undecidable by this solid and something should say so. Expected
incidence on the current corpus: near zero, which is the point — it should be *measured* to be near
zero rather than assumed. Same battery as §16.

### N5. Two silent tie-breaks remain, and they are now the only ones

v1 §16 states the scope honestly and it has not changed: `DistFromOutside`, `DistFromInside` and
`ComputeNormal` consume an `onTrimBoundary` hit without noticing. Since parity is now self-checking
and the oracle columns are at zero, these are the last places where the trim-accuracy sliver
converts an undecidable input into a confident wrong answer. Their failure mode is worse than
parity's: a wrong distance is a step into a wall, and `nearestCrossing` has no equivalent of the
vote to escalate to.

This is not urgent on the current corpus (0 disagreements on both distance columns), and it becomes
*moot* on any part converted by Section 5. Recommendation: leave it, and mark it explicitly as
"deliberately unhandled, mooted by the CSG route where that applies, to be revisited if the residual
corpus keeps it" — rather than letting it drift as an unlabelled known gap.

### N6. Directory hygiene is now a real cost

v1 §7 mentioned it; it has not moved and the folder has grown. `scripts/geometry/` currently holds
five stale `O2_TGeoToCAD*` variants, four `O2_CADtoTGeo_with_*` variants, editor backups (`*.py~`,
`#...#`), macOS resource forks (`._*`), `a.out`, `tgeo_to_cad.o`, and a `__pycache__`. The
TGeo→CAD scripts in particular are about to matter again (Section 6.1) and it is currently not
possible to tell which of the five is the live one. One commit, `git rm` plus a `.gitignore`, before
the beampipe work starts.

---

## 4. Feasibility, re-evaluated

v1 §9's gates, honestly re-scored, and re-judged for whether they are still the right target.

| gate | v1 target | today | still the right goal? |
| --- | --- | --- | --- |
| **G0** kernel integrity | ctest green + adversarial cases | **61 cases green**, BVH == Loop everywhere, K1/K2/K4/K5/K8 closed with regression tests | **Yes.** Keep as the permanent floor. |
| **G1** synthetic Boolean ladder | 10/10 | **8/9 converted** (`oblique_cut_cyl` does not convert at all — planar ellipse trim, C-series) | **Yes**, and it should grow — see §6.3. |
| **G2** as1-oc-214 | 5/5 exact + oracle-clean | not re-measured under the oracle gate | **Downgrade.** as1 is subsumed by the ladder; measure it once, do not track it. |
| **G3** Bagger 13/13 | 13/13 exact, Reliable, oracle-clean | **6/12 pass; 12/12 oracle-clean** | **Restate.** "Oracle-clean" is achieved and is the property that matters. The 6/12 is a *closure and capacity* score. Report the two separately, permanently. |
| **G4** oTOF + ALICE3 coverage | honest coverage, no silent unreliable shipping | 20/55 ALICE3 extracted; the `auto`-mode fallback policy is **still not implemented** | **Yes, and it is now the most overdue item.** Shipping an unreliable exact part silently is the one live production risk on this branch. |
| **G5** performance | Safety BVH; within ~2x of tessellated | untouched; gate timing column is "unable to resolve anything below a few per cent" | **Yes but defer.** Section 5 changes what the performance target even is. |

**What has genuinely become more feasible:** exact navigation on real CAD. It is done, measured
against an independent kernel, on 21 parts. Nobody should re-litigate that.

**What has become less feasible, or at least differently shaped:** "one trimmed-surface solid that
is watertight on imported CAD". Every session since v1 has narrowed the residual and the residual has
not gone to zero; the last two narrowed it to a single irreducible term — *two charts cannot agree
about a curve neither can express* — which v1 §4.2 predicted, correctly, three weeks ago, and which
Phase 2 was scoped to remove at a cost of 1-2 weeks of kernel work plus a sidecar format change plus
a converter change, leaving a free-form remainder (E3) still to do afterwards.

That estimate should now be compared honestly against the alternative, which did not exist as a
concrete option when v1 was written.

**On `oblique_cut_cyl`**: v1 listed the planar-ellipse trim as a 3-line converter fallback
(C-series). It is still blocking a G1 fixture and a Bagger part (`Bucket`). It is the cheapest
un-taken item on the board and it should just be done, independently of everything else.

---

## 5. The orthogonal proposal: recognise and emit boolean solids

The user's thought, stated back: instead of representing a part as a set of trimmed surfaces, detect
when it *is* a boolean combination of primitives and emit that — union, intersection, subtraction of
tubes, cones, boxes, tori — which would handle tube-tube exactly and might be cheaper. And bring
extrusions in as well.

**This is correct, it is stronger than stated, and it has substantial production precedent.** It is
also, in its strongest form, a direct replacement for v1's Phase 2 rather than a complement to it.

### 5.1 What it actually is, and what is known about it

The general problem is **B-rep → CSG conversion**, and it has a real literature and a real answer
for exactly our surface families:

- **Shapiro & Vossler** (*CAD* 23(1) 1991; *ACM TOG* 12(1) 1993, "Separation for boundary to CSG
  conversion") proved that a solid bounded by *natural quadrics* (planes, cylinders, cones, spheres)
  admits a CSG representation over the halfspaces of its own faces, provided that halfspace set is
  "separating", and gave a construction — including how to add auxiliary separating halfspaces when
  it is not.
- **Buchele & Crawford** (2004) give a practical 3D halfspace CSG-tree construction from implicit
  B-reps.
- The **fusion-neutronics community has been doing this in production for fifteen years**, because
  MCNP accepts only CSG. **McCAD** (KIT) and, more relevantly, **GEOUNED** (UNED/CERN, open source,
  2023-) convert STEP straight to MCNP/OpenMC/Serpent CSG cells. GEOUNED is built on **OCCT via
  FreeCAD** — the same kernel this project already has installed and already drives from Python. It
  handles planes, cylinders, cones, spheres and tori; it is used on ITER-scale models.

The core algorithm they converged on is not pattern-matching and is not heuristic:

> **Decompose by extending the carriers.** Take the solid's faces' carrier surfaces. Repeatedly cut
> the solid by an extended carrier that resolves a non-convexity, until every piece is convex. Each
> convex piece is then simply the **intersection of the oriented halfspaces of its own faces**, and
> the solid is the **union of the pieces**:
>
>     solid = ∪_i ( ∩_j H_ij )
>
> — a two-level disjunctive normal form over halfspaces, with no trims, no pcurves, no wires, no
> edges and no watertightness question of any kind.

Everything this project has fought for a month — pcurves that disagree, rims that will not pair,
bands that cannot be chosen, closure checks that structurally cannot pass — **does not exist in that
representation.** Membership is `OR` over cells of `AND` over sign tests. Two faces cannot disagree
about a shared curve because there are no faces.

The user's instinct that "tube-tube intersections can be treated precisely" is exactly this: a tube
crossing a tube is `H_cyl1 ∩ H_cyl2 ∩ H_plane…` or a subtraction of the two — and the transcendental
intersection curve, the object that has blocked this project since day one, **is never
represented at all**. It is implied by two sign tests.

### 5.2 Three tiers, in increasing power and increasing cost

They should be built in this order, and each is independently shippable.

**Tier 1 — primitive recognition (cheapest, highest immediate value, no new kernel).**
Recognise that a part's whole face set is one ROOT primitive: `TGeoTube`/`TGeoTubeSeg`,
`TGeoCone`/`TGeoConeSeg`, `TGeoSphere`, `TGeoTorus`, `TGeoBBox`, `TGeoPcon`, `TGeoPgon`,
`TGeoArb8`. Detection is a small closed-form matcher over the recognised carriers the converter
*already produces* (v1 §2 credits canonical recognition as sound and machine-precision). Emission is
a one-line TGeo constructor. No new C++ at all.

Value: this is not a marginal case. `BVHSurfaceSolid.md` records ALICE3 as 1419 cylinder + 1398 cone
+ 72 sphere faces recognised at machine precision; a mechanical assembly is full of plain washers,
pins, collars, shims and spacers. Every one that matches becomes an exact, fast, ROOT-native,
Geant4-native, visualisable shape with no sidecar, no BVH, no closure check and no risk.

**Tier 2 — recognised boolean of primitives (the user's proposal literally).**
Detect the common two- and three-primitive patterns: `tube − tube` (the Bagger junction), `tube ∪
tube`, `box − cyl` (bolt holes), `primitive − N holes`. Emit `TGeoCompositeShape` over
`TGeoBoolNode`. Verification is Section 5.4. Still no new C++ kernel.

Value: this covers the six Bagger cylinder parts — the exact set that is failing today — and `tube
− cyl` covers a very large share of bracket/plate geometry. Note that `TGeoBoolNode` accepts *any*
`TGeoShape`, so an `O2BVHSurfaceSolid` may be a leaf: a part that is "one awkward trimmed body minus
four bolt holes" can be emitted as a composite with the awkward part kept in the existing
representation. That hybrid is a real and underrated option.

**Tier 3 — general convex decomposition (the GEOUNED/Shapiro construction).**
For any part whose faces are all analytic: split by extended carriers until convex, emit the DNF.
This is the general answer and it is where the research risk lives (Section 5.6).

For Tier 3 I would **not** emit a deep `TGeoCompositeShape`. A DNF with 20 cells is 20 nested
boolean nodes, each of which re-transforms the point and re-dispatches virtually; ROOT's boolean
navigation is the wrong shape for that. Emit instead a flat `O2CSGSolid`: a `std::vector` of cells,
each a small contiguous run of oriented quadric coefficient blocks. Then

- `Contains` = `OR_i AND_j (sign of a quadratic form)` — a few dozen FLOPs, no branches worth
  speaking of, no memory indirection, no BVH.
- `DistFromOutside`/`DistFromInside` = per cell, intersect the ray's parameter interval with each
  halfspace's inside-interval (a quadratic root pair, or a half-line for a plane) — a Cyrus-Beck
  clip generalised to quadrics — then take the min/max over cells. All the root solvers this needs
  **already exist and are already validated** in `BoundedSurface.h`.
- `Safety` = min over cells of the min over halfspaces of the exact distance-to-quadric already
  implemented.
- `Capacity` = no quadrature at all for Tier 1/2; for Tier 3, the divergence theorem over the cells,
  or simply keep the existing exact primitives' closed forms.

And this is the representation that **fits the AOT-codegen decision recorded in v1 §10 Phase 4**
better than anything the trimmed-surface path can offer: a cell list of fixed length with fixed
coefficients is a straight-line generated function with no loops over heterogeneous virtual
surfaces. That direction was chosen for the current solid; it is much more natural here.

### 5.3 Cost — is it actually cheaper?

The user asked. My expectation, stated as a prediction to be measured, not a claim:

| | `O2BVHSurfaceSolid` (Bagger cylinder part, ~23 faces) | Tier-1/2 primitive or composite | Tier-3 flat CSG (~5 cells × ~5 halfspaces) |
| --- | --- | --- | --- |
| `Contains` | BVH descent + N patch intersections, each with a point-in-wire query walking B-spline polylines; **plus** a possible 5-direction re-shoot | one or two exact primitive `Contains` | ~25 quadratic-form evaluations, no allocation, no traversal |
| `DistFromOutside` | ray-BVH traversal, per-patch quartic/quadratic roots, cluster classification, possible unpruned redo | ROOT primitive distances, composite bookkeeping | per-cell interval clipping, min over cells |
| build cost | sidecar parse + BVH build + `CloseShape` (rim matching over all chords, capacity quadrature) | nothing | nothing |
| memory | surfaces + wires + cached polylines + BVH | a few dozen bytes | a few hundred bytes |
| correctness risk | trims, bands, closure | none beyond the primitive | none beyond acceptance |

I expect Tier 1 and Tier 3 to be **substantially** faster than the current solid and Tier 2
(`TGeoCompositeShape`) to be the uncertain one — ROOT's boolean navigation has a reputation for
being the slow path, and a composite's `Safety` in particular is weak. That uncertainty is cheap to
remove: **the harness already times any `TGeoShape`.** Building `cyl_inter_cyl` by hand as a
`TGeoCompositeShape` and running it through `runSolidHarness` against the existing solid and
`O2Tessellated` is an afternoon and it decides Tier 2's emission target empirically. Do that before
committing to `TGeoCompositeShape` as the Tier-2 output.

The honest caveat on all of the above: the gate's timing column "cannot resolve anything below a few
per cent" and moved 3× between identical runs (`NEXT.md`). Any of these numbers must come from a
dedicated fixed-point-set loop, not the gate.

### 5.4 Verification — this is why it is feasible *now* and was not before

A recognition pipeline is only as good as its acceptance test, and the branch has spent two sessions
building exactly the right one. Three independent checks are available and all three are cheap:

1. **OCCT symmetric difference — the decisive one.** Build the candidate CSG in OCCT
   (`BRepPrimAPI_MakeCylinder` etc. + `BRepAlgoAPI_Cut`/`Fuse`/`Common`), then compute
   `volume(A − B) + volume(B − A)` against the original shape. **Zero (to model tolerance) is a
   proof of equality**, not a sample-based argument. This is a dozen lines of pythonOCC against a
   kernel already installed and already driven by `occtOracle.py`. Nothing in the trimmed-surface
   path has ever had an acceptance test this strong.
2. **The existing oracle gate**, unchanged. A recognised part is still a `TGeoShape`; `runOracleGate.py`
   scores it on `contains`/`distout`/`distin`/`safety`/`capacity` against OpenCascade exactly as
   today. The gate does not care how the shape is implemented — which is the payoff for having built
   it against an external reference rather than against ourselves.
3. **`BRepGProp` volume and bounding box**, as a fast pre-filter before the expensive boolean.

So the pipeline can be as heuristic, incomplete and greedy as we like on the *proposal* side,
because the *acceptance* side is exact. That is the architecture: **propose cheaply, verify
exactly, fall back silently.** A part that fails acceptance keeps `O2BVHSurfaceSolid`; a part that
fails that keeps `O2Tessellated`. Three tiers of graceful degradation, each strictly better than the
next, each independently verified.

### 5.5 Extrusions and revolutions — and why the beampipe makes this urgent

The user is right that these become important, and the beampipe is the reason.

**Revolutions.** A part all of whose carriers are surfaces of revolution about one common axis —
planes ⊥ to it, cylinders, cones, spheres and tori coaxial with it — *is* a solid of revolution, and
its `(r, z)` profile is a closed loop of line and **arc** segments. `TGeoPcon` can express the line
part only; the arcs (every fillet, every bellows corrugation, every torus knuckle) must be
polygonised, which is precisely the exactness loss this project exists to avoid.

An exact `O2RevolvedSolid` — a `(r, z)` profile of lines and arcs, swept through a φ range — is a
**small, well-scoped new solid that reuses machinery already written and validated**: `Curve2D`
already carries lines and arcs, `CurveWire` already does winding and Green's-theorem areas, and the
ray intersection reduces to a 2D problem in `(r, z)` per segment (a line segment gives a cone/plane
root, an arc gives a torus/sphere root — all four solvers exist). Capacity is Pappus/Green in closed
form. `Safety` is a 2D distance to the profile.

**The ALICE beampipe is a surface of revolution almost everywhere.** `Detectors/Passive/src/Pipe.cxx`
builds it today from `TGeoTube`, `TGeoPcon`, `TGeoTorus` and `TGeoCompositeShape` — i.e. someone has
already hand-written the CSG. A CAD beampipe converted through an exact revolved-profile solid would
be *one shape per section*, exact, with ~10 parameters, instead of dozens of trimmed patches with
B-spline seams. This is the single highest-value new solid for the stated target.

**Extrusions.** Symmetrically: all carriers either ⊥ to a common direction `ê` (two cap planes) or
∥ to it (planes and cylinders) ⇒ a prism over a 2D profile of lines and arcs. Same kernel reuse,
same 2D reduction (ray → 2D ray-vs-profile plus a slab clip). `TGeoXtru` is **not** an adequate
target because its profile is polygonal only — arcs are exactly what mechanical profiles are full
of. So this too wants a small exact solid, `O2ExtrudedSolid`, sharing the `Curve2D` profile code
with the revolved one. Between them they cover a very large fraction of mechanical CAD: brackets,
plates, flanges, sheet metal, collars, shims.

Both are also trivially AOT-codegen-able and both are trivially watertight — a profile loop that
closes in 2D closes in 3D, by construction.

Note the pleasing consequence: `O2RevolvedSolid` and `O2ExtrudedSolid` make **`tube_window`, the one
remaining failing G1 fixture, and most of the beampipe, exact without any CSG at all.**

### 5.6 Where it does not work, and what it costs

Stated plainly, because the case for this is strong enough not to need overselling:

- **Free-form faces are out.** Any part with a genuine B-spline *surface* cannot be a halfspace CSG.
  ALICE3 is 5775 free-form faces by the recorded classification — an upper bound, but large. Those
  parts keep `O2BVHSurfaceSolid` (if their trims are analytic) or the tessellated fallback. **The
  CSG route does not retire the existing solid; it removes its hardest inputs.**
- **Tier 3 can blow up.** The decomposition's cell count is worst-case exponential in the number of
  faces, and the boolean operations that perform the splits (`BOPAlgo`) are the least robust part of
  OCCT, notoriously so on tangencies and near-coincident faces. GEOUNED's practical experience is
  that the pathologies are real and are handled by falling back per part. Mitigation is structural:
  gate Tier 3 by face count (say ≤ 40) and by a cell-count cap, and fall back on failure. The
  1505-face oTOF plate is *not* a CSG candidate and should never be attempted as one — it is exactly
  what the BVH solid is good at. **The two representations are complementary in precisely the right
  way: CSG scales badly in face count and beautifully in boolean complexity; the BVH solid is the
  reverse.**
- **Fillets and blends.** A constant-radius fillet is a cylinder/torus halfspace and is fine. A
  variable-radius or setback blend is free-form and is not.
- **Tolerance.** A CSG built from carriers recognised at machine precision is exact; the acceptance
  test is against OCCT's *tolerant* model, so "equal" means equal to model tolerance, same bar as
  today. Nothing gets worse.
- **Effort.** Tier 1: ~3-4 days including the acceptance harness (Python only). Tier 2: ~3 days on
  top. `O2RevolvedSolid` + `O2ExtrudedSolid`: ~1.5-2 weeks together, C++ kernel plus converter plus
  tests, and they share most of their code. Tier 3: ~2-3 weeks with real risk of overrun; **do not
  start it until Tiers 1-2 and the profile solids have been measured**, because they may well leave
  little for it to do.

Compare against v1 Phase 2 (adjacency trims): 1-2 weeks, kernel + sidecar v3 + converter, benefiting
only the analytic-analytic seams, leaving E3 (shared-edge sampling for free-form) still to build
afterwards, and leaving the closure/band machinery in place forever because free-form parts still
need it.

### 5.7 Recommendation

**Do the simplification pipeline, and do it instead of Phase 2, in the order Tier 1 → profile solids
→ Tier 2 → (only if still needed) Tier 3.**

Reasoning, compressed:

- It removes the failure class rather than approximating around it. Phase 2 makes two charts agree;
  CSG removes the charts.
- Every case Phase 2 was scoped for (v1 §4.2's table: all of Bagger, all of as1, 373 ALICE3 faces) is
  a case the CSG route covers, and it covers them with *less* new kernel code, not more.
- It has an exact acceptance test (symmetric difference) where Phase 2 has only the same gate.
- It produces shapes that are portable downstream — ROOT visualisation, ROOT's overlap checker,
  Geant4 native solids, and any future VecGeom path — where a custom `TGeoShape` is a permanent
  integration liability. (O2 today has no VecGeom dependency, so nothing breaks either way; but
  emitting `TGeoTube` is strictly more portable than emitting anything of ours.)
- It is the right fit for the AOT-codegen direction already chosen.
- It fails gracefully into work that is already done and already validated.

**Keep from the current plan, regardless:** N3 (edge IDs — they are what makes any watertightness
claim checkable, and they are also the input Tier 3's decomposition wants), N1 (the capacity fix,
because `Capacity` must be right for the parts that stay on the BVH solid), the `auto`-mode fallback
policy (G4), and the ellipse trim (C-series).

**Drop or defer:** Phase 2 as specified; Phase 3 (shared-edge sampling) — its remaining
justification is free-form parts, and those are a coverage question, not a correctness one; and any
further tuning of the rim match band, which N3 makes moot.

---

## 6. Scaling up: the ALICE beampipe

### 6.1 The corpus does not exist here, and there is a better one available for free

There is no beampipe CAD in this repository. `scripts/geometry/STEP_examples/` holds `as1-oc-214`,
`Bagger.step` and `oTOF System V3-R92cm.step`; `ALICE_3_example/` holds the 32 MB ALICE3 assembly.
Acquiring a real beampipe STEP from the ALICE CAD repository is on the critical path and should be
started immediately, in parallel with everything else, because it has a lead time that no amount of
local work shortens.

Meanwhile — and this is worth doing **even after** the real STEP arrives — there is a ground-truth
beampipe already in this repository, written in TGeo:

`Detectors/Passive/src/Pipe.cxx` (3176 lines) and `PipeRun4.cxx` build the Run 3 beampipe from
`TGeoTube`, `TGeoPcon`, `TGeoTorus`, `TGeoCone` and `TGeoCompositeShape` — the beryllium section,
the bellows, the adaptors, the carbon support with its screw holes and sight holes, the vespel
rings. Round-trip it:

    TGeo (exact, known)  →[O2_TGeoToCAD]→  STEP  →[O2_CADtoTGeo]→  shape  →  compare against the original TGeo

This gives four things nothing else does:

1. **An exact oracle.** Not OCCT's tolerant classification — the original `TGeoPcon`. Every
   disagreement is a real defect with no tolerance excuse. The harness already compares two
   `TGeoShape`s; this is a reference-shape swap, not new machinery.
2. **A known answer for the recognition pipeline.** We know the source was a tube, a pcon, a torus,
   a composite — so Tier 1/2 recognition can be *scored*, not just accepted. "Did we recover the
   shape the engineer drew" is a question no CAD corpus can ask and this one answers exactly.
3. **A torture test for the exact families we claim to support**: tori (bellows), pcons (tapers),
   composites with holes (the carbon support), all at real detector coordinates.
4. **Scale-dependence coverage, at last.** v1 S12 flagged that nothing has ever been tested away
   from the origin, and that "ALICE geometries live at metres". The beampipe spans z ≈ ±4 m with
   features of ~0.08 cm (the beryllium wall). That is a 10⁴ dynamic range against absolute constants
   like `kBVHBoxTolerance = 1e-3` cm and `kWireJoinTolerance = 1e-5` cm. **I expect this to break
   things**, and finding out is most of the value of the exercise.

Cost: the `O2_TGeoToCAD*.py` scripts exist (five stale variants — see N6). Budget ~2 days to revive
one, plus the round-trip harness. Do it before the real CAD arrives, so the machinery is ready.

### 6.2 What will break at scale, and the measurement for each

Ranked by my estimate of the probability of biting, each with what to run:

| # | risk | why | measurement |
| --- | --- | --- | --- |
| 1 | **Absolute tolerances at metre coordinates** (S12, still open) | `kBVHBoxTolerance` 1e-3 cm, `kWireJoinTolerance` 1e-5 cm, `kRayTolerance` 1e-9, `sameIntersection` at 1e-7·\|t\| — a ray from z = −400 cm carries ~1e-5 cm of relative-merge width by the time it reaches z = +400 | run the whole G1 ladder translated by (0, 0, 400) and scaled ×10 and ×0.1; **every column must be unchanged** |
| 2 | **`CloseShape` cost per part × part count** | rim matching is a grid search over all chords of all faces; capacity was 0.15 s/part before the Green's fix and 0.012 s after — but ALICE3 is 55 volumes and a real assembly is thousands | time `CloseShape` vs face count on oTOF's 1505-face part; establish the scaling law before assuming it |
| 3 | **Memory** | sidecar + wires + cached flattened polylines + BVH, per part, all resident | measure resident set for the full ALICE3 conversion; a 179-pole spline per seam is not small |
| 4 | **The gate does not scale** | `runOracleGate.py` is per-part and serial; the oracle is ~0.3 s per 2000 samples | parallelise per part; add a `--sample-budget` so a 1000-part model is scored in bounded time |
| 5 | **Silent unreliable shipping** (G4) | `auto` mode ships an exact part that is not `Reliable` and only warns | implement the fallback policy; it is a ~20-line change and it is the only production-risk item on the branch |
| 6 | **Conversion wall-clock** | ALICE3 is 32 MB of STEP; recognition adds work per face | profile the converter once; it has never been profiled |

Item 1 deserves emphasis: **it is the one that could invalidate results already recorded.** Every
measurement on this branch was taken on parts a few centimetres across near the origin. A single
translated-ladder run is cheap and either confirms the whole corpus or reframes it.

### 6.3 What the harness and fixtures need

- **Translate/scale sweep** over the existing ladder (item 1 above) — cheapest high-value addition.
- **Ladder extensions** for the new representations: `pcon_with_fillets` (arcs in the profile — the
  bellows knuckle), `extruded_bracket_with_holes` (lines + arcs profile + through holes),
  `tube_with_bolt_circle` (Tier-2 subtraction of N cylinders), `oblique_cut_cyl` (already there,
  still unconverted), and one deliberately **non**-CSG part with a free-form face, to pin the
  fallback path.
- **A recognition scorecard** in the gate output: per part, which tier accepted it, what the
  symmetric-difference volume was, and what it fell back to. Coverage becomes a measured, tiered
  number instead of a single fraction.
- **Per-part wall-clock and memory** in `gate.json`, so §6.2's scaling laws accumulate for free.

---

## 7. Recommended plan, ordered

Each item is small enough to finish and verify in one session, and each states its gate.

**Now — three things that are cheap, independent, and unblock reading everything else.**

1. **N1: fix the contour quadrature.** Run the two experiments first (`kContourMaxSpanU` sweep;
   midpoint integrator at high N), then make `contourIntegralAlongCurve` span-aware. *Gate:* the six
   Bagger cylinders and `tube_window` reach ≤1e-11 relative capacity; the two 1e-9 wire-trim capacity
   tests stay green; every other column bit-identical. This closes the last failing numeric column
   on the board.
2. **N2 + N4 + N6.** Loader B-spline endpoints (3 lines + a clamping check reported not enforced);
   ambiguity-aware voting (measure the incidence, expect ~0); directory cleanup. *Gate:* ctest green,
   both gates bit-identical.
3. **G4: the `auto`-mode fallback policy.** An exact part that is not `Reliable` must not ship
   silently. *Gate:* a test that a deliberately-unreliable part converts to tessellated under `auto`
   and throws under `required`.

**Next — decide the architecture on evidence, in about a week.**

4. **N3: sidecar v3 carries edge IDs.** Converter emits a per-model edge table and a
   `(edgeId, sense)` per trim curve; the loader carries them; `validateClosure` decides on identity
   and *reports* `maxSharedEdgeDeviation` in cm. *Gate:* on every currently-`Reliable` part, the
   verdict is unchanged and the new deviation number is below the old match band; on the six
   cylinders, the verdict becomes closed-by-topology with a stated deviation. This is also the first
   time the project can quote a defensible "how far apart are the two faces" number.
5. **Tier 1 recognition, in the converter, Python only.** Match whole face sets to ROOT primitives;
   accept by OCCT symmetric difference; emit the primitive. *Gate:* every G1 fixture that is a
   primitive is recognised; symmetric-difference volume ≤ model tolerance; the oracle gate passes on
   the emitted shape; recognition coverage reported per model.
6. **The `TGeoCompositeShape` timing probe.** One afternoon, decides Tier 2's emission target.
7. **Tier 2: `tube − tube` and `primitive − N holes`.** *Gate:* the six Bagger cylinder parts and
   `tube_window` convert as booleans, pass the oracle gate, and pass symmetric difference. **This is
   the point at which G3 becomes 12/12 and the tube-tube problem is finished.**

**Then — the beampipe, and the solids it needs.**

8. **Beampipe round-trip harness** (§6.1) plus the translate/scale sweep (§6.2 item 1). *Gate:* the
   G1 ladder is column-identical at (0,0,400) and at ×10/×0.1 scale — or the constants that fail are
   named and fixed.
9. **`O2RevolvedSolid` and `O2ExtrudedSolid`** — exact `(r,z)` and `(x,y)` profiles of lines and
   arcs, sharing `Curve2D`. *Gate:* new ladder fixtures; exact capacity in closed form; oracle-clean;
   `Pipe.cxx`'s bellows section reproduced exactly through the round trip.
10. **Real beampipe STEP**, when it arrives, through everything above, with the tiered coverage
    scorecard.

**Deferred, explicitly:** Phase 2 (adjacency trims) — revisit only if Tiers 1-2 and the profile
solids leave an analytic-analytic residual that matters. Phase 3 (shared-edge sampling) — a coverage
item, not a correctness one. Tier 3 (general convex decomposition) — only on measured evidence that
Tiers 1-2 leave enough on the table to justify its risk. `Safety` BVH and the rest of G5 — after
correctness, and after Section 5 has changed what is being optimised.

---

## 8. What would falsify this review

Stated so it can be checked rather than believed:

- **N1 is wrong** if the `kContourMaxSpanU` sweep does not move the six parts' capacity. Then the
  defect is in the trim *region*, the midpoint comparison will show it, and `NEXT.md` thread 2's
  remaining candidate — a wire that wraps the periodic chart without its cut segments — becomes the
  live one. (I checked that one from the source and believe it is handled: the converter appears to
  emit the vertical seam runs as ordinary curve-to-curve joins, which the bridge integrates because
  they have `deltaV != 0`. But it is the natural second suspect and it is testable by printing, per
  wire, whether the contour closes in the chart.)
- **The Green's-theorem seam-bridge skip is correct** and I claim this without needing a run: a
  bridge with `deltaV == 0` contributes `∫ F dv = 0` identically, whatever its `u` extent. `NEXT.md`
  thread 2's first candidate can be closed on paper.
- **Section 5 is wrong for this project** if Tier 1 recognition, run over Bagger + as1 + oTOF +
  ALICE3, accepts only a small share of parts. That is one Python session to find out and it should
  be the *first* thing done in step 5 — measure the recognition rate before building the emission.
  If the rate is low, the CSG route is a Bagger-specific optimisation and Phase 2 comes back.
- **The beampipe round-trip is worthless** if `O2_TGeoToCAD` cannot produce a faithful STEP of a
  `TGeoCompositeShape`. Check that on one composite before budgeting the two days.
- **Item 1 of §6.2 could invalidate recorded results.** If the translated ladder disagrees, then
  every "0 disagreements" number on this branch is a statement about geometry near the origin only,
  and must be re-qualified in v1 §20, `NEXT.md` and `TolerancePolicy.md`.

---

## 9. Note on the documents

`CodeReview_Fable.md` (v1) is now 1078 lines and interleaves findings, plans, corrections and
executed-session records; three of its plan sections are superseded and two of its findings (K7,
K4-as-open) are corrected in later sections of the same file. It is still the register of record for
the K/S/C findings and should stay, but the *plan* in §9-10 should carry a one-line pointer to this
document so nobody executes it by accident. `NEXT.md` is the right size and the right idea and
should keep being rewritten every session. `TolerancePolicy.md` (1026 lines) has become the
measurement archive, which is a good role for it — but §3.3, §9.1 and §12.2 all contain predictions
that were later refuted in the same file, and each should carry a superseded-by marker at the point
of the claim, not only at the point of the correction.
