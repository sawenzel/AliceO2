# Stream E — position, scale, and what the branch's numbers actually mean

Date: 2026-08-01. Companion to [`Workstreams.md`](Workstreams.md) (Stream E brief),
[`TolerancePolicy.md`](TolerancePolicy.md) (the constants) and
[`Stream_F_EdgeIdentity.md`](Stream_F_EdgeIdentity.md) (whose result this stream ended up
measuring from the outside).

Written to be folded into `NEXT.md` by the integrating session. Nothing here edits a document this
stream does not own.

---

## 0. The verdict, first

**No. This branch's recorded results are position-independent but *not* scale-independent.**

- **Position:** at (0, 0, 400) cm the whole G1 ladder is unchanged on every gate-decisive column.
  Zero unexplained oracle disagreements, zero missed surfaces, every verdict, every count, every
  reliability string identical. 9/9, exactly as at the origin. The only movements are diagnostic
  reals whose *noise floor* moves with the coordinate magnitude — which is the correct and expected
  behaviour of double arithmetic, not of a tolerance.
- **Scale:** at ×10 the ladder is likewise clean, and on a *relatively ten times tighter* band than
  the baseline (§3). At **×0.1 it is not**: `torus_union_cyl` loses one of 2000 `distout` rays as a
  missed surface, and the gate goes 9/9 → 8/9. The invariant this project defends —
  *zero disagreements with the oracle outside tolerance on every column* — **is broken by shrinking
  the geometry by a factor of ten, and by nothing else.**

The cause is named, localised to one line, and reproduced from three independent directions:

> `solveQuarticReal` (`Detectors/Base/src/BoundedSurface.h:944`) guards Ferrari's second stage with
> ```cpp
> if (resolvent > kTolerance) {
> ```
> `kTolerance` is **1e-9 cm**; `resolvent` is the largest root of the resolvent cubic and has units
> of **cm²**. The comparison is dimensionally inconsistent, and when it fails the function returns
> **the empty root set** — not a merged pair, not a tangency, no roots at all — for a ray that
> transversally enters and exits the torus at 69° incidence.

**What has to be re-qualified.** Every "0 disagreements" statement on this branch is a statement
about geometry *at least a few centimetres across*. It is safe at metre coordinates. It is safe at
ten times the current size. It is **not** safe at a tenth, and §5 shows the defect is already
firing at the shipped size at a rate of ~2e-4 per torus-hitting ray — the ladder missed it by
about half an expected event, not by a margin.

`torus_union_cyl` is the only part on either corpus affected, because it is the only part with a
toroidal surface. Nothing else in the sweep moved a verdict at any scale between ×0.1 and ×10.

---

## 1. How the experiment was made like-for-like

The sweep is only evidence if the transformed run asks *the same questions* of a shape that
differs from the baseline's *only* by the transform. Three things had to be arranged for that, and
each was checked rather than assumed.

**The shape.** `make_boolean_fixtures.py --transform` applies a `gp_Trsf` (translation,
uniform scaling about the origin, or a `;`-separated composition) to each fixture **before the
STEP is written**. The STEP is the single source from which the converter, the `surfaces_*.bin`
sidecar the kernel loads, the `facets_*.bin` mesh and the `brep_*.brep` the oracle answers from are
all derived. So there is exactly one shape in the pipeline and **the oracle is transformed by
construction**, not by a second arrangement that could drift from the first.

**The samples — and why this was the hard part.** `generateSamples()` rejection-samples through the
*tessellated* reference, and BRepMesh does **not** produce the transformed triangulation of a
transformed shape. Measured, on `torus_union_cyl` translated by 400 cm: same 23432 triangles, but
only **1284 of them** survive an undo-the-shift-and-sort comparison. On `sphere_minus_cyl`,
5606 of 7048. Under scaling it is worse still, because the chordal deflection is absolute: ×10
nearly doubles the ladder's triangle count (36256 → 67138) and ×0.001 collapses it to 954.

So the sample sets would have differed for a reason that has nothing to do with the kernel. The fix
is `--load-samples` (new): the harness reads a frozen sample set instead of generating one, and
`transformSamples.py` pushes the baseline's set through the **same `gp_Trsf`** — using OCCT's own
transform, and `gp_Pnt` for ray origins against `gp_Dir` for ray directions, because transforming a
unit direction as a point would denormalise it under a scaling and silently change the question.
Point *i* of the transformed run is then exactly *T*(point *i*) of the baseline, in the same
category, in the same order, and the oracle answers those same points on the transformed `.brep`.

*Self-check:* feeding a generated sample set straight back in through `--load-samples` reproduces
**112 of 112 fields identical on all 9 parts**, including every oracle column. The JSON round trip
is bit-exact, so `--load-samples` adds nothing of its own.

**The comparison.** `compareGateRuns.py` diffs two `gate.json` reports field by field, with a
declared scale exponent per field (`_FIELD_EXPONENT`): counts and relative deviations are
invariant, lengths carry exponent 1, capacities exponent 3. Timing (`timing*`, anything containing
`Seconds`, `nsPerCall`, `checksum`) is stripped, per the contract. `--gate-columns-only` drops the
columns that compare against the mesh, since the mesh is deliberately not the same mesh.

**What the comparison can and cannot detect.** It compares integers for equality and reals against
`baseline × factor^exponent` within 1e-9 relative. It therefore detects any change of verdict,
count, or reliability string exactly, and any change of a real column above that floor. It cannot
detect a change in a field `gate.json` does not carry, a change below 1e-9 relative, a relative
perturbation of a field that is exactly zero, or anything about points the harness never sampled.

---

## 2. The positive control

Required before any green result may be reported, because *a refuted experiment is not a refuted
hypothesis* and the last session lost a day to a sweep that was structurally incapable of saying
"yes".

**(a) The comparator, against injected defects.** `compareGateRuns.py --self-test` breaks the
baseline four ways — one integer column, one exponent-1 length, one exponent-3 volume, one verdict
string — and reports whether each is caught. **4/4.**

This earned its place immediately. Its *first* version scored **2/3**, missing the length-column
injection, because it perturbed `maxRimIsolation` on the first fixture by 1e-8 *relative* and that
field is exactly `0.0` on `box`. That is the same failure mode as `max(1, ceil(0/x))`: an
experiment that cannot move the number it is testing. The injector now searches for a part where
the field is non-zero and names it, and reports "no part carries a non-zero value, nothing was
injected" rather than silently scoring a miss.

**(b) The sweep, against real violations.** The sweep found one without being asked to, at ×0.1
(§3), which is the strongest positive control available. Pushed further:

| transform | parts | unexplained | missed surfaces | gate |
| --- | --- | --- | --- | --- |
| baseline | 9 | 0 | 0 | 9/9 |
| translate (0, 0, 400) cm | 9 | 0 | 0 | 9/9 |
| scale ×10 | 9 | 0 | 0 | 9/9 |
| **scale ×0.1** | 9 | 0 | **1** | **8/9** |
| **×0.1 then +400 cm** | 9 | 0 | **1** | **8/9** |
| **scale ×0.001** | 9 | **2060** | **525** | **5/9** |

(`scale:1000` was also attempted and abandoned — see §8; it is a tessellation problem, not a kernel
one.)

The ×0.001 row is the deliberate break: the ladder is then 20–40 µm across, against
`kBVHBoxTolerance` 1e-3 cm = 10 µm, and four parts fail on three columns. The sweep says "yes"
loudly when there is something to say.

Note the fourth and fifth rows are identical. **Translation adds nothing to the failure** — the
defect is purely a function of size.

**(c) The comparator itself had a can-never-fail bug, and was caught by being used.** `compareGateRuns.py`
originally keyed parts by the leading component of the part id, which is the fixture name on G1 and
therefore unique — but on a single CAD model *every* part shares it (`Bagger/BasePin…`,
`Bagger/Base…`). The first Bagger comparison it was pointed at reported "1 part(s) vs baseline's 1"
and compared exactly one of twelve. Two things follow, and both are the point of writing this down:
the fixture results above were never affected (their keys are genuinely unique, and the pairing is
now printed on every run), and **a comparison tool has to be run against a corpus it was not
designed on before it is believed**. It now keys on the full part id, falls back to the leading
component only when the ids differ *and* that component is unique, and states which it used.
The corrected Bagger run: 12 parts, 8 differing fields each, all of them the new edge-identity
fields — 96 differences, nothing pre-existing moved.

---

## 3. The sweep tables

Invariant chosen, per column:

| column | invariant | why |
| --- | --- | --- |
| `oracle.<c>.nMismatchUnexplained` + `nMismatchMissedSurface` | **identical integer** | this is the invariant the project defends; it is dimensionless |
| `navigation.reliability`, `navigable`, rim/edge counts, `nSurfaces` | **identical** | topology |
| `oracle.capacityRelativeDeviation` | **identical to within its own noise** | a ratio; the volumes scale as *f*³ and cancel |
| lengths (`maxRimIsolation`, `maxSharedEdgeDeviation`, `totalRimLength`, `worstDeviation`) | **scale as *f*¹** | cm |
| `oracle.tolerance`, `navigation.rimMatchTolerance` | **do not scale** — see §4 | absolute constants of the oracle and the closure check |

### 3.1 Oracle disagreements outside tolerance, per part (the invariant)

| part | baseline | +400 cm | ×10 | ×0.1 | ×0.1 +400 cm |
|---|---|---|---|---|---|
| `box` | 0 | 0 | 0 | 0 | 0 |
| `box_minus_cyl` | 0 | 0 | 0 | 0 | 0 |
| `box_union_box` | 0 | 0 | 0 | 0 | 0 |
| `cyl_cross_cyl` | 0 | 0 | 0 | 0 | 0 |
| `cyl_inter_cyl` | 0 | 0 | 0 | 0 | 0 |
| `cyl_plus_cone` | 0 | 0 | 0 | 0 | 0 |
| `sphere_minus_cyl` | 0 | 0 | 0 | 0 | 0 |
| **`torus_union_cyl`** | **0** | **0** | **0** | **1** | **1** |
| `tube_window` | 0 | 0 | 0 | 0 | 0 |

Every part is `reliable` and `navigable` in every run, including ×0.001.

### 3.2 Capacity, relative to OpenCascade

| part | baseline | +400 cm | ×10 | ×0.1 |
|---|---|---|---|---|
| `box` | 1.5e-16 | 1.5e-16 | 0 | 2.9e-16 |
| `box_minus_cyl` | 2.4e-15 | 3.9e-15 | 2.5e-15 | 2.5e-15 |
| `box_union_box` | 0 | 0 | 0 | 2.9e-16 |
| `cyl_cross_cyl` | 1.8e-12 | −7.3e-12 | 1.8e-12 | 1.8e-12 |
| `cyl_inter_cyl` | −1.0e-11 | 7.6e-11 | −1.0e-11 | −1.0e-11 |
| `cyl_plus_cone` | 6.5e-14 | 8.7e-12 | 8.5e-14 | 7.2e-14 |
| `sphere_minus_cyl` | −7.0e-14 | 6.5e-13 | −9.2e-14 | −1.2e-13 |
| `torus_union_cyl` | 9.3e-14 | 8.5e-14 | 7.6e-14 | 6.8e-14 |
| `tube_window` | −9.44e-08 | −9.39e-08 | −9.44e-08 | −9.44e-08 |

Capacity is **scale-exact to the last digits** — the ×10 and ×0.1 columns reproduce the baseline to
3–4 significant figures on every part, and `tube_window`'s residual (the known trim-data floor) is
identical to 3 digits at all three sizes. Under translation the exact-arithmetic parts lose about
three digits (1e-14 → 1e-11): at *z* ≈ 400 cm the contour integrand carries coordinates 400× larger
and cancels them, which is ordinary conditioning, not a tolerance. The combined
"×0.1 then +400 cm" case loses the same three digits (worst 8.8e-10) and is still four orders
inside the 1e-6 gate.

### 3.3 The closure numbers — and an unplanned result about Stream F

| part | metric | baseline | +400 cm | ×10 | ×0.1 |
|---|---|---|---|---|---|
| `cyl_cross_cyl` | `boundaryEdges` (proximity) | 9 | 9 | **1175** | **11** |
| `cyl_cross_cyl` | `maxSharedEdgeDeviation` (identity) | 3.812e-09 | 3.812e-09 | **3.812e-08** | **3.812e-10** |
| `tube_window` | `boundaryEdges` (proximity) | 1269 | 1269 | **3523** | **692** |
| `tube_window` | `maxSharedEdgeDeviation` (identity) | 4.712e-09 | 4.712e-09 | **4.712e-08** | **4.712e-10** |

Read those two rows against each other. The proximity-based boundary-edge count swings by **two
orders of magnitude** over a 100× size range on the *same shape*; `maxSharedEdgeDeviation` scales
**exactly linearly, to four significant digits**, in both directions.

This was not what the sweep was aimed at, and it is the cleanest external confirmation Stream F
could have asked for. It also says something sharper: at ×10, `cyl_cross_cyl` and `cyl_inter_cyl`
report **1175 boundary edges and are still `reliable` and `navigable`**, because the verdict now
comes from edge identity. Under the proximity rule those two parts would have gone non-navigable
at ten times their current size. **Had this sweep been run one wave earlier, it would have
reported a scale-dependent navigability verdict.** Stream F removed that before it was measured.

`maxRimIsolation` tells the same story from the other side: its *absolute* value barely moves over
a 100× size range (6.7e-6 → 4.2e-6 at ×10 → 1.07e-5 cm at ×0.1), i.e. it is not a length of the
geometry at all, so relative to the part it degrades by ~10× per decade of shrinking. That is the
"coincidence of scale, not a criterion" the previous hand-over named, now measured over two
decades.

---

## 4. What is absolute and stays absolute (correctly or not)

Three quantities did **not** scale with the geometry, and each is worth stating because a reader of
the tables will otherwise think they are failures.

1. **`oracle.tolerance` is 1e-7 cm at every scale.** It is OCCT's `BRep_Tool::Tolerance`, an
   absolute number the modeller attaches, and it is the gate's own comparison band
   (`oracleOpt.meshBand`). Consequence: **the ×10 run passed on a relatively ten times tighter band
   than the baseline, and the ×0.1 run on a ten times looser one.** The ×10 result is therefore
   *stronger* than the baseline, and the ×0.1 failure is a failure despite a more forgiving band.
2. **`navigation.rimMatchTolerance` is 1e-8 cm** (or the model's declared tolerance) at every
   scale. Absolute by design — see `kRimMatchTolerance`'s comment. Now largely moot for the
   verdict, per §3.3.
3. **The oracle's own `_RAY_EPS` = 1e-9 cm** classifier tolerance. Visible as `nNoVerdict`, the
   count of sample points OCCT calls `ON` the boundary: 18002 at baseline, 18001 at ×10, 18011 at
   ×0.1, 19452 at ×0.001. It moves in exactly the direction an absolute band must, and it is the
   *oracle's* constant, not the kernel's.

The kernel constants the brief listed behaved as follows. `kBVHBoxTolerance` (1e-3 cm) was never
implicated: the BVH handed the surface kernel 2 candidate patches on the offending ray at *every*
scale, and `--loop-crosscheck` reported BVH == `_Loop` on 2000/2000 rays in every run. Down at
×0.001, where the whole ladder is 20–40 µm and 1e-3 cm is 10 µm, it presumably is implicated, but
that regime is not a claim anyone needs. `kWireJoinTolerance` (1e-6 cm) and the rim machinery never
moved a verdict at any scale. **The one constant that broke something is `kTolerance`, and not in a
tolerance role at all — in a dimensionally-inconsistent guard inside a root solver.**

---

## 5. The defect, in detail

### 5.1 The measurement chain

The offender is `x0.1/gate.json`, `torus_union_cyl`, `oracle.distout.worstOffenders[0]`: origin
(0.20943, 0.32925, 0.19348) cm, direction (−0.75473, −0.031547, −0.65528), oracle distance
0.24384686727698293 cm, candidate `1e30` (miss).

**Step 1 — it is the same ray.** The sweep's sample sets are exact images, so this ray is
literally *T*(a baseline ray) that the kernel answers correctly. Asking each sidecar the
correspondingly scaled question (`/tmp` probe against `-lO2DetectorsBase`, per
`TolerancePolicy.md` §6):

| sidecar | expected | `DistFromOutside` | `DistFromOutside_Loop` | BVH candidates |
| --- | --- | --- | --- | --- |
| ×0.1 | 0.2438468673 | **1e30 (MISS)** | **1e30 (MISS)** | 2 |
| baseline | 2.438468673 | 2.438468673 | 2.438468673 | 2 |
| ×10 | 24.38468673 | 24.38468673 | 24.38468673 | 2 |

Relative deviation where it hits: 1.6e-15 and 4.7e-15. So the kernel is exact at both larger sizes
and returns *no intersection at all* at the smaller one, for the identical configuration.

**Step 2 — it is not the BVH and not the trim.** BVH and `_Loop` agree (both miss), the BVH hands
over 2 candidate patches at every scale, and the crossing list along that ray is **empty**
(`DescribeContainsCrossings`: BVH 0, Loop 0) where at ×0.15 it is `t=0.36577 ENTER, t=0.51997 EXIT`.
Meanwhile, at the very same point, the *other* two kernels are correct:

```
scale 0.12   t=0.99x expected: Contains=0  Safety=0.0010678
             t=1.001x expected: Contains=1  Safety=0.000104662
```

`Contains` (which shoots in a different direction) and `Safety` both know the surface is there. Only
this ray's intersection is lost.

**Step 3 — the incidence is ordinary.** The hit point lies at exactly the minor radius from the
tube circle, and dot(direction, surface normal) = **−0.358**, i.e. 69° incidence. This is not a
tangency being lost to conditioning; it is a transversal crossing being discarded.

**Step 4 — the line.** Replaying the torus quartic (`TorusBoundedSurface::appendIntersections`)
for that ray at a series of scales and printing what `solveQuarticReal` branches on:

| global scale | R [cm] | r [cm] | `termP` | `termQ` | `resolvent` | roots |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.5 | 0.8 | 17.463 | −6.36e-03 | 6.250e-08 | 2 |
| 0.5 | 1.25 | 0.4 | 4.366 | −7.95e-04 | 1.563e-08 | 2 |
| 0.3 | 0.75 | 0.24 | 1.572 | −1.72e-04 | 5.625e-09 | 2 |
| 0.2 | 0.5 | 0.16 | 0.699 | −5.09e-05 | 2.500e-09 | 2 |
| 0.15 | 0.375 | 0.12 | 0.393 | −2.15e-05 | 1.406e-09 | 2 |
| **0.12** | 0.3 | 0.096 | 0.251 | −1.10e-05 | **9.001e-10** | **0** |
| **0.1** | 0.25 | 0.08 | 0.175 | −6.36e-06 | **6.250e-10** | **0** |
| **0.05** | 0.125 | 0.04 | 0.044 | −7.95e-07 | **1.563e-10** | **0** |

`resolvent ∝ L²` exactly (6.250e-08 · s²). The guard is `resolvent > kTolerance = 1e-9`, so the
predicted threshold is s = √(1e-9 / 6.2504e-8) = **0.1265**. The measured threshold, by converting
the fixture at nine scales and probing each, is **between 0.12 and 0.15**. The model predicts the
bisection.

### 5.2 How close is the shipped geometry?

`resolvent` is configuration-dependent as well as scale-dependent, so "the fixture is 62× above the
guard on this one ray" is not the number that matters. Firing 40 000 isotropic rays at a torus of
each size and counting the ones where the quartic changes sign (ground truth, no tolerance
involved) but `solveQuarticReal` returns nothing:

| torus | R [cm] | r [cm] | rays that hit | reported as misses | rate |
| --- | --- | --- | --- | --- | --- |
| ×1000 | 2500 | 800 | 5162 | 0 | 0 |
| ×100 | 250 | 80 | 5162 | 0 | 0 |
| ×10 | 25 | 8 | 5162 | 0 | 0 |
| **×1 (as shipped)** | 2.5 | 0.8 | 5162 | **1** | **1.9e-04** |
| ×0.5 | 1.25 | 0.4 | 5162 | 1 | 1.9e-04 |
| ×0.2 | 0.5 | 0.16 | 5162 | 1 | 1.9e-04 |
| ×0.15 | 0.375 | 0.12 | 5162 | 2 | 3.9e-04 |
| ×0.12 | 0.3 | 0.096 | 5162 | 3 | 5.8e-04 |
| ×0.1 | 0.25 | 0.08 | 5162 | 4 | 7.8e-04 |
| ×0.01 | 0.025 | 0.008 | 5162 | 11 | 2.1e-03 |

**The defect is live at the shipped fixture size**, at about 2e-4 of torus-hitting rays. The G1
gate fires 2000 outside rays at `torus_union_cyl`, of which a fraction hit the torus, so the
expected count of lost rays at baseline is well under one — which is precisely why five sessions of
green gates never saw it. Shrinking by ten raises the rate fourfold and the ladder caught it.
Above ~25 cm major radius it does not fire at all in this sample.

For ALICE geometry the practical reading is: **toroidal surfaces smaller than a few centimetres are
at risk today, at the origin, at production scale.** Beampipe bellows knuckles and small fillets are
exactly that size class.

### 5.3 What the fix is not

Not this stream's commit, deliberately (the brief: report with evidence, fix as a separate
decision). But three things are worth recording so the fix is not scoped wrong.

- The mathematical condition for Ferrari's non-biquadratic branch is `resolvent > 0`. Any numerical
  guard has to be relative to the coefficient magnitudes — something of the order
  `eps · max(1, |termP|, √|termR|)` — not a fixed length.
- **The same file has two more dimensionally-inconsistent guards in the same function.**
  `if (std::abs(termQ) <= kTolerance)` selects the biquadratic branch, and `termQ ∝ L³`
  (−6.4e-03 at ×1, −6.4e-09 at ×0.01): below about ×0.01 the solver silently takes the biquadratic
  branch on a quartic that is not biquadratic. And the Newton polishing guard
  `std::abs(derivative) > kTolerance` compares an *L³* quantity against 1e-9 too. Fixing only the
  resolvent guard leaves two.
- `sameIntersection` is `|t1 − t2| ≤ 1e-7 · max(1, |t1|, |t2|)`. The `max(1., …)` means it stops
  being relative below 1 cm and becomes an absolute 1e-7 cm band. It was **not** implicated here
  (the roots are 0.1 cm apart), but it is the same class of defect and is on the same code path.

---

## 6. Step 2 — `--edge-identity`

Requested by Stream F, which correctly did not add it because this stream owns the file.
`HasEdgeIdentity()`, the six source-edge counts (total, shared, boundary, non-manifold, reversed,
degenerate) and `GetMaxSharedEdgeDeviation()` now go into `gate.json`'s `navigation` block
unconditionally, and to stdout under `--edge-identity`, which `runOracleGate.py` always passes.

What it reports on the two corpora today: **`hasEdgeIdentity` is true on all 9 fixtures and all 12
Bagger parts**; every part has `sharedSourceEdges == sourceEdges` with zero boundary, non-manifold,
reversed and degenerate source edges. `maxSharedEdgeDeviation` ranges from exactly `0.0`
(`torus_union_cyl`) through 4.1e-13 cm (`cyl_plus_cone`) to 4.7e-09 cm (`tube_window`), the largest
on the fixture ladder. As §3.3 shows, it is the only closure number in `gate.json` that is exactly
scale-covariant.

---

## 7. Instruments added, and where they live

| file | what it is |
| --- | --- |
| `runSolidHarness.cxx --load-samples <dir>` | read a frozen sample set instead of generating one; the prerequisite for any cross-shape comparison |
| `runSolidHarness.cxx --edge-identity` | Stream F's numbers on stdout; always in `--json` |
| `make_boolean_fixtures.py --transform` | `translate:x,y,z` (mm), `scale:f`, or `a;b` composition, applied to the shape before the STEP is written |
| `transformSamples.py` | pushes a frozen sample set through the same `gp_Trsf`; `gp_Pnt` for origins, `gp_Dir` for directions |
| `runOracleGate.py --transform / --load-samples` | the two above, wired into the gate |
| `compareGateRuns.py` | per-field diff with a declared scale exponent per field, plus `--self-test` |
| `testBVHSurfaceSolid.cxx`, Stream E block | scale covariance of the torus quartic where it works, and a *characterisation* test of the resolvent guard that will fail loudly when the guard is fixed |

Reproducing the whole sweep:

```bash
# baseline
runOracleGate.py --workdir /tmp/base --fixtures
# transformed run, comparable point by point
transformSamples.py --in /tmp/base/oracle --out /tmp/s_x01 --transform scale:0.1
runOracleGate.py --workdir /tmp/x01 --fixtures --transform scale:0.1 --load-samples /tmp/s_x01
# read it
compareGateRuns.py --baseline /tmp/base/gate.json --candidate /tmp/x01/gate.json \
                   --scale 0.1 --gate-columns-only
compareGateRuns.py --baseline /tmp/base/gate.json --self-test
```

`transformSamples.py` takes the spec in **millimetres**, the same string as
`make_boolean_fixtures.py`, and converts to cm itself — the two cannot disagree about units.

---

## 8. Two smaller findings, recorded because they will be needed

**The tessellated fallback is stored in `float32`.** `facets_*.bin` is `<9f` per triangle while the
exact sidecar is `<d`. Half a float32 ulp at *z* = 400 cm is 1.5e-05 cm, and that is exactly the
size of the movement seen in every mesh-reference column of the translated run (worst 2.7e-05 cm on
`cyl_inter_cyl`'s `distin`). At ALICE's outer radii the tessellated fallback's vertex quantisation
is ~3e-05 cm; at 40 m it is ~3e-03 cm. This is a property of the fallback representation, not of
the exact solid, but `auto` mode ships it (Stream E step 3, still open) and the number belongs on
the record before that policy is decided.

**Meshing cost is driven by an absolute deflection, and it is a wall, not a slope.** `--mesh-prec`
defaults to 0.1 in *model units* and is passed as **both** the linear and the angular deflection —
the second of which is in radians and must not be scaled with the first. ×10 nearly doubles the
ladder's triangle count (36256 → 67138); ×0.001 collapses it to 954 triangles for nine parts, i.e.
an unusable reference. **×1000 never completed:** converting a single fixture,
`sphere_minus_cyl` at a 2 m radius, reached **22.9 GB resident and 9 minutes of CPU** before it was
killed to save the machine, on a 31 GB box. One sphere.

That is a hard result for Stream E step 2 (the ALICE3 corpus gate, 55 solids and 9266 faces at
detector scale) and it lands before that work starts: **the tessellation step, not `CloseShape`, is
the thing that will not scale**, and `makeTestPartDB.py` currently offers no way to set the mesh
precision per model. The `--sample-budget` in the brief is necessary but not sufficient; a
`--mesh-prec` passthrough with a *relative* linear deflection and a fixed angular one is the
prerequisite. The ×1000 sweep point is therefore not reported above — the interesting direction was
downwards and ×0.001 covers it.

---

## 9. Before and after

| | ctest | fixtures gate | Bagger gate | unexplained disagreements |
| --- | --- | --- | --- | --- |
| before (branch state at hand-over) | green, 73 cases | 9/9 | 9/12 | 0 / 0 (both corpora, all four columns) |
| after | green, **75** cases | 9/9 | 9/12 | 0 / 0 (both corpora, all four columns) |

The two new cases are the Stream E block. Nothing this stream committed changes an answer: the
new-harness fixture baseline differs from the old one in exactly 72 fields (8 new edge-identity
fields × 9 parts) and the Bagger run in exactly 96 (× 12 parts), and in nothing else. The three
Bagger failures remain the
capacity floor (`BoomCylinderOuter` 1.39e-06, `BucketCylinderInner` 1.31e-06,
`StickCylinderOuter` 1.39e-06 against a 1e-06 band), which is wave 0's trim data.

## 10. What is left of the Stream E brief

- **Step 1 — done**, and it changed a recorded result. §0 is the re-qualification.
- **Step 2 — done.** §6.
- Step 3 (`--sample-budget`, per-part wall clock and RSS in `gate.json`, the ALICE3 corpus gate)
  and Step 4 (G4, the `auto`-mode fallback policy) are **not started**. On Stream F's question of
  whether `NavigationReliability` should gate on `maxSharedEdgeDeviation`: §3.3 is the evidence
  that it should — it is the only closure number that is a length of the geometry rather than of
  the sampling — but the threshold has to be *relative* to the part, or it reintroduces exactly the
  defect §5 documents. That is a decision for whoever does step 4, with this table in front of
  them.
- The beampipe round-trip and the ladder extensions are untouched.
