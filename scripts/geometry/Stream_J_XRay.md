# Stream J — the X-ray / geantino transport benchmark

Date: 2026-08-01. Companion to [`Tutorial.md`](Tutorial.md) (the map),
[`SolidNavigationHarness.md`](SolidNavigationHarness.md) (the single-shot harness),
[`Stream_G_AnyShape.md`](Stream_G_AnyShape.md) (the representation abstraction this reuses) and
[`Stream_E_Scale.md`](Stream_E_Scale.md) §3/§5 (the quartic defect used here as a live positive
control). Written by the session that built it; the integrating session folds it.

**In one paragraph.** Everything this project measures is a **single-shot** query: from a sampled
point, how far to the surface. A transport loop is different in kind — step, land *on* the
boundary, step again from there — and that is where geometry navigation actually fails: zero-length
steps, ping-ponging on a face, a particle that enters and never exits, a crossing found twice, a
step that overshoots into the next volume. **None of those can be expressed as a disagreement on
`distout` from an interior sample, so the existing gate is structurally blind to all of them.**
This stream builds the instrument that can express them: a structured parallel-beam raster fired
through a part, producing the **ordered crossing list** per ray by *stepping*, two independent ways
— a direct shape-API loop and the real `TGeoNavigator` — compared against OpenCascade **as lists,
not as aggregates**. The result on the shipped representations is clean: **62208/62208 fixture rays
and 82944/82944 Bagger rays with bit-for-bit identical crossing lists on the surface solid, 0
mode (a) vs mode (b) disagreements anywhere, and every robustness counter zero.** The mesh is a
different story (§3, §4), and the known torus quartic defect *is* detected — but only after the
instrument's first version failed to detect it and had to be shown capable of doing so, which is
§6 and is the most useful thing in this document. It also **runs on ALICE3** — 18 parts, full
oracle round trip, 82 s and 361 MB, no mesh needed — where it finds four parts that edge-identity
closure certifies as watertight and that a transported particle nonetheless walks through (§7).

---

## 1. What it measures, and how to run it

Three stages, mirroring the oracle gate, and one driver that chains them.

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH   # stage FIRST

cd $HOME/alisw/O2
# convert + raster + oracle + score, in one command
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xray --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xray_bag \
    --model scripts/geometry/STEP_examples/Bagger.step

# or reuse a finished oracle-gate workdir's part database, with no reconversion
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xray_bag \
    --reuse-db /tmp/gate_bag/db --raster 48
```

Both self-tests need no database, no model and no oracle:

```bash
$B/stage/bin/o2-bench-detectorsbase-xray --self-test        # 19 checks
$SW/Python/latest/bin/python3.10 scripts/geometry/xrayOracle.py --self-test   # 11 checks
```

| piece | what it is |
| --- | --- |
| `Detectors/Base/test/XRayTransport.h` | the algorithms: the stepping loop, the per-ray audit, the list comparator, the raster. Header-only and **not** in `DetectorsBase`, so no object file of the gate path is rebuilt differently. Included by both the benchmark and the unit tests, so the tests exercise the code that runs. |
| `Detectors/Base/test/runXRayBenchmark.cxx` → `o2-bench-detectorsbase-xray` | the front end: representations, both stepping modes, the navigator world, JSON |
| `scripts/geometry/xrayOracle.py` | ground truth: the ordered crossing list per ray, from the part's `.brep` |
| `scripts/geometry/runXRayBench.py` | the driver, and the three tables |

### The two stepping modes, and why both

| mode | what it is |
| --- | --- |
| **(a) `shape`** | `Contains()` to establish the starting state, then alternating `DistFromOutside()` / `DistFromInside()`, advancing the point along the ray until the accumulated distance leaves the raster window. **Depends on nothing but the shape.** |
| **(b) `nav`** | the part placed in a `TGeoVolume` inside a minimal world, transported with `TGeoNavigator::FindNextBoundaryAndStep()`. **The production path.** |

Build (a) first and report both. If they disagree, that isolates *the shape* from *the navigator*
immediately; with only (b) you cannot tell which one lied. This is not hypothetical — it is exactly
how the undersized-world defect in §6.2 was caught and correctly attributed to the benchmark rather
than to ROOT.

Two details of mode (a) that are load-bearing:

* **the push.** After landing on a boundary the point must be advanced *past* it or the next query
  re-finds the same crossing at zero. The default is `--push 1e-9` cm = the kernel's own
  `kRayTolerance`, so the recorded crossing distances carry a known bias of at most (k−1)·push over
  a k-crossing ray, i.e. below 1e-8 cm — two orders under the 1e-6 cm comparison band. A *stalled*
  step is repaired with a larger `--unstick-push` (1e-6 cm) and **every use is counted**; a repair
  that is not counted is a lie.
* **`stepmax` is deliberately not used to bound the query.** Its semantics differ between shape
  implementations (some return the crossing, some return `stepmax`, some return `Big`), and this
  loop must measure the crossing list rather than that convention. The window is enforced on the
  returned crossing distance instead.

### The oracle, and why it is not just the intersector

`IntCurvesFace_ShapeIntersector` returns every crossing along a ray in **one call**, which is what
makes the oracle side cheap (§7). But its hits are *face* hits, not *solid* transitions: a ray
through a shared edge hits two faces at the same parameter, a ray tangent to a cylinder hits it
twice and never goes inside. Taking those at face value would put spurious crossings into the
ground truth, and a ground truth with spurious crossings fails a correct kernel.

So the intersections are only **candidate positions**. The list is decided by classifying the
**midpoint of every interval** with `BRepClass3d_SolidClassifier` and keeping only the positions
where the classification changes. That is a different algorithm from the one under test, it is
immune to tangency and edge double-hits, and it produces an alternating enter/exit list *by
construction* — which is what licenses treating a non-alternating candidate list as a defect. A
midpoint OCCT calls `ON` flags the ray `amb`; the benchmark excludes it rather than scoring it,
exactly as the sample gate's `nNoVerdict` points are excluded. **Across both corpora, 0 rays were
ambiguous.**

The oracle's self-test pins this: a ray tangent to the outer wall of a hollow cylinder must produce
**zero** crossings, not two. It does.

### The comparison keeps LOST apart from DISPLACED

This is the whole localising value of comparing lists rather than aggregates:

| column | meaning |
| --- | --- |
| **`LOST`** | a crossing OpenCascade found and the candidate did not — **a wall a track walks through** |
| `extra` | the reverse: a crossing the candidate invented |
| `displaced` | a crossing found in the right order and sense but more than the match tolerance away — **a wrong step length, not a lost wall** |
| `kind` | right position, wrong sense |
| `identical` | rays whose *whole ordered list* matched. This is the number a transport actually depends on. |

The first version merged the first three and reported `tube_window`'s mesh as *2952 missing + 2952
extra* when the truth is **0 lost, 2952 displaced by at most 1.29e-03 cm** — the chording sagitta.
The two statements have completely different consequences and the first one is wrong.

---

## 2. Crossing lists against OpenCascade

Raster 48 × 48 per beam, three axis beams = **6912 rays per part**, match tolerance 1e-6 cm.
Every representation a part has, both stepping modes, against the same oracle answers.

### Fixture ladder — 9 parts

| representation | mode | rays identical | LOST | extra | displaced | worst \|Δt\| (cm) | parts clean |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `surface` | (a) shape | **62208 / 62208** | **0** | **0** | **0** | 4.93e-14 | **9/9** |
| `surface` | (b) nav | **62208 / 62208** | **0** | **0** | **0** | 4.93e-14 | **9/9** |
| `shape` (CSG) | (a) shape | **13824 / 13824** | **0** | **0** | **0** | 6.22e-15 | **2/2** |
| `shape` (CSG) | (b) nav | **13824 / 13824** | **0** | **0** | **0** | 6.22e-15 | **2/2** |
| `mesh` | (a) shape | 28711 / 62208 | 40 | 6 | 70750 | 2.70e-02 | 2/9 |
| `mesh` | (b) nav | 28711 / 62208 | 40 | 6 | 70750 | 2.70e-02 | 2/9 |

### Bagger — 12 parts

| representation | mode | rays identical | LOST | extra | displaced | worst \|Δt\| (cm) | parts clean |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `surface` | (a) shape | **82944 / 82944** | **0** | **0** | **0** | 3.47e-11 | **12/12** |
| `surface` | (b) nav | **82944 / 82944** | **0** | **0** | **0** | 3.47e-11 | **12/12** |
| `shape` (CSG) | (a) shape | **48384 / 48384** | **0** | **0** | **0** | 2.48e-11 | **7/7** |
| `shape` (CSG) | (b) nav | **48384 / 48384** | **0** | **0** | **0** | 2.47e-11 | **7/7** |
| `mesh` | (a) shape | 46088 / 82944 | 1430 | 0 | 92081 | 2.93e-02 | 0/12 |
| `mesh` | (b) nav | 46088 / 82944 | 1430 | 0 | 92081 | 2.93e-02 | 0/12 |

**The two shipped representations reproduce OpenCascade's crossing list ray for ray, at both
corpora, in both stepping modes, to 5e-14 cm on the ladder and 3.5e-11 cm on Bagger.** That is a
stronger statement than the gate's, because it is about a *sequence* of queries each starting from
where the previous one landed, not about independent samples.

### Where the mesh loses walls

All **1430** Bagger LOST crossings are on **`Bagger/BucketLink2`** — the one part in either corpus
whose mesh is not a closed body (`meshClosedBody=false`, Stream G §4 item 5). The gate reported
that part as 366 `contains` disagreements with a 2.2 cm worst deviation; the X-ray benchmark says
what that costs a *track*: 1430 boundary crossings that a transported particle would walk straight
through, on 2605 of 6912 rays.

On the ladder the LOST crossings are on the two curved-intersection fixtures (`cyl_cross_cyl` 16,
`cyl_inter_cyl` 24), and both also produce **unterminated transports** — see §4.

---

## 3. The robustness table — the part nothing else measures

Same runs. Every counter is a count of events, never a rate.

| corpus | repr | mode | steps | zeroStep | noAdv | unstick | capHit | **unterm** | oddList | nonAlt | dupXing | **parity** | noTrans | outWorld | (a) vs (b) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| fixtures | surface | (a) | 183560 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | – | – | — |
| fixtures | surface | (b) | 183560 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | 0 | 0 | **0 / 62208** |
| fixtures | shape | (a) | 39424 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | – | – | — |
| fixtures | shape | (b) | 39424 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | 0 | 0 | **0 / 13824** |
| fixtures | mesh | (a) | 183526 | 0 | 0 | 0 | 0 | **6** | 6 | 0 | 0 | **15** | – | – | — |
| fixtures | mesh | (b) | 183526 | 0 | 0 | 0 | 0 | **6** | 6 | 0 | 0 | **15** | 0 | 0 | **0 / 62208** |
| Bagger | surface | (a) | 212758 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | – | – | — |
| Bagger | surface | (b) | 212758 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | 0 | 0 | **0 / 82944** |
| Bagger | shape | (a) | 123072 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | – | – | — |
| Bagger | shape | (b) | 123072 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **0** | 0 | 0 | **0 / 48384** |
| Bagger | mesh | (a) | 211328 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **534** | – | – | — |
| Bagger | mesh | (b) | 211328 | 0 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | **534** | 0 | 0 | **0 / 82944** |

Column meanings:

| column | meaning |
| --- | --- |
| `zeroStep` | a step at or below 1e-9 cm |
| `noAdv` | the accumulated distance did not increase |
| `unstick` | a stalled step that had to be nudged to continue |
| `capHit` | the 512-iteration cap reached without leaving the window |
| **`unterm`** | **the ray ended INSIDE the solid: entered and never left** |
| `oddList` | odd-length crossing list (the same event, counted as the brief names it) |
| `nonAlt` | two consecutive crossings of the same sense |
| `dupXing` | two crossings closer together than the match tolerance |
| **`parity`** | **`Contains()` at an interval midpoint contradicts the crossing list** |
| `parityNB` | the same, excused because the midpoint is within the match tolerance of the boundary (0 everywhere on both corpora) |
| `noTrans` | mode (b): a boundary was crossed but the volume did not change |
| `outWorld` | mode (b): the ray origin was not in the navigator world — **a misconfiguration of this benchmark, never a geometry defect** |

**Four things to read off this table.**

1. **Zero zero-length steps, zero non-advancing steps, zero unstick pushes, zero iteration-cap
   hits — on every representation, both modes, 1.15 million steps.** The transport loop never
   stalls on any part of either corpus. That was the single most likely failure mode going in and
   it does not happen. A 1e-9 cm push off the boundary is enough for every one of these shapes.
2. **`unterm` is the headline new failure mode, and only the mesh has it.** Six rays enter
   `cyl_inter_cyl`'s tessellation and never leave it. A single-shot `distin` query cannot express
   this: it reports one wrong distance and the sample is scored. In a transport it is a lost track.
3. **`parity` is the sharpest column, and it is the only one independent of the stepping.** Both
   modes produce an alternating list *by construction*, so `nonAlt` can never fire on them; asking
   the shape's own `Contains()` at the midpoint of every interval is the only way this instrument
   can contradict itself. It fires **534 times on Bagger's mesh** (all on `BucketLink2`) and **15
   times on the ladder's mesh** — i.e. `Contains()` and `DistFrom*` on the *same* `O2Tessellated`
   disagree about which side of the surface a point is on. **Zero on the surface solid and zero on
   the CSG shapes.** `parityNB` is 0 everywhere, so none of these is excused by near-tangency.
4. **Mode (a) and mode (b) never disagree — 0 rays out of 62208 + 13824 + 82944 + 48384.**
   `TGeoNavigator` transports every part of both corpora exactly as the bare shape API does. That
   is a real result: it says the navigator adds no error of its own on any representation here, and
   it is why every number above can be attributed to the shape.

### A robustness property that fell out by accident, and is worth keeping

The fixtures corpus was run twice with ray positions differing by at most **1.8e-15 cm** (≈ 8
double ulp; the two runs computed the same lattice through two algebraically identical but
differently-rounded expressions). The `surface` and `shape` representations were **bit-identical**
across the two runs — 0 LOST both times. The **mesh moved from 14 to 40 LOST crossings**, all on
`cyl_cross_cyl` and `cyl_inter_cyl`. Bagger's mesh did not move (1430 both times).

So a tessellated part's crossing list is not reproducible under a last-ulp change of ray position,
where the exact representations are. That is what a chorded silhouette does — a ray grazing it can
fall on either side of a triangle for a perturbation far below any physical tolerance. It is one
more measured input to the `auto`-mode fallback policy, and it is the reason mesh counts in this
document are quoted with the run that produced them.

---

## 4. Volume by chord integration — and what it can and cannot resolve

**Read this before the numbers.** An earlier session built a Monte-Carlo `Contains` witness aimed
at the 1.1e-2 cm³ capacity residual, validated it, and found it **3× too coarse** to separate what
it was aimed at. A structured raster is better than random sampling, but boundary cells still
dominate. The measurement below says the achieved precision at N = 48 is **2.9e-05 to 4.9e-03
relative**, against a divergence-theorem capacity that already sits at **1e-11** on exact parts.
**It is four to eight orders coarser and it cannot resolve the 1.3e-06 capacity residuals. It must
not be quoted as if it could.**

What it *is* for: **gross errors, and composites** — `TGeoCompositeShape::Capacity()` is
Monte-Carlo in ROOT (~1e-2 relative) and nothing else independent exists for a CSG part.

### Three volume numbers, answering three different questions

| number | what it measures |
| --- | --- |
| candidate chord volume **vs the oracle's chord volume over the same rays** | the part's volume agreement. **Immune to the raster's own error**, which cancels exactly. |
| oracle chord volume **vs OCCT's exact volume** | **the raster's achieved precision**, measured, at the stated density |
| `Capacity()` vs OCCT's exact volume | the number the sample gate already scores |

### The raster's achieved precision, measured

Two closed-form self-checks fix it without any CAD at all:

* **A box's chord integral is EXACT at every raster density** when the window is its own bounding
  box — no convergence argument and no tolerance. `BOOST_CHECK` at N = 5, 16, 41 and in the
  benchmark's self-test at N = 7, 32.
* **A sphere of radius 1**, single beam, window = bounding box:

| N | chord volume (cm³) | exact | relative |
| ---: | ---: | ---: | ---: |
| 24 | 4.18792146 | 4.18879020 | 2.07e-04 |
| 48 | 4.18853263 | | 6.15e-05 |
| 96 | 4.18915612 | | 8.74e-05 |
| 192 | 4.18897443 | | 4.40e-05 |

**The convergence is not monotone in N**, because the silhouette cells realign with the lattice at
every density. The honest statement is therefore an *envelope at a stated density*, never an
extrapolation, and the unit test asserts the measured envelope (2e-3 over N = 24…192) rather than a
rate.

### The measured table, fixture ladder, N = 48, three beams

`raster vs exact` is OCCT's own chord integral against OCCT's exact volume — the instrument's
error, identical for every row of a part. `vs OCCT chord` is the candidate against that same
integral, and is the column that means something about the part.

| part | repr | chord V (cm³) | vs OCCT chord | raster vs exact | `Capacity()` vs exact |
| --- | --- | ---: | ---: | ---: | ---: |
| `box` | surface / mesh / shape | 24.034679 | **0** / **0** / **0** | 1.45e-03 | 1.5e-16 |
| `box_minus_cyl` | surface | 56.057979 | **1.2e-14** | 1.80e-03 | 2.4e-15 |
| | mesh | 56.060830 | 5.1e-05 | | 6.0e-05 |
| `box_union_box` | surface | 48.053348 | **−3.0e-16** | 1.11e-03 | 0 |
| | mesh | 48.053348 | −3.0e-16 | | 0 |
| `cyl_cross_cyl` | surface | 32.524976 | **2.6e-14** | 4.92e-03 | 1.8e-12 |
| | shape (CSG) | 32.524976 | **2.6e-14** | | **−4.2e-03** |
| | mesh | 32.457042 | −2.1e-03 | | −3.0e-04 |
| `cyl_inter_cyl` | surface | 5.3334864 | **−8.0e-15** | 2.87e-05 | −1.0e-11 |
| | mesh | 5.3205914 | −2.4e-03 | | −3.7e-04 |
| `cyl_plus_cone` | surface | 13.094830 | **−6.6e-15** | 3.71e-04 | 6.5e-14 |
| | mesh | 13.087162 | −5.9e-04 | | −7.1e-04 |
| `sphere_minus_cyl` | surface | 29.095832 | **5.7e-15** | 2.08e-04 | −7.0e-14 |
| | mesh | 29.038402 | −2.0e-03 | | −1.9e-03 |
| `torus_union_cyl` | surface | 78.694491 | **9.0e-16** | −9.46e-04 | 9.3e-14 |
| | mesh | 78.652308 | −5.4e-04 | | −5.5e-04 |
| `tube_window` | surface | 36.638195 | **0** | 9.75e-04 | −9.4e-08 |
| | mesh | 36.629954 | −2.2e-04 | | −3.4e-04 |

**The result worth having is the CSG column.** On the same rays, the exact surface solid and the
`TGeoCompositeShape` produce chord volumes agreeing to **2.6e-14 relative** — while ROOT's own
`TGeoCompositeShape::Capacity()` for that shape is **4.2e-03** off. On Bagger the same comparison
holds on all seven CSG-carried rams, with `Capacity()` deviations of 9.3e-04 to 1.45e-02:

| Bagger part | surface chord V | shape chord V | shape vs surface | ROOT `Capacity()` vs OCCT |
| --- | ---: | ---: | ---: | ---: |
| `BasePin` | (identical to 1e-13) | | ≤ 1e-13 | −1.7e-03 |
| `BoomCylinderInner` | 22.535358 | 22.535358 | −7.5e-14 | −1.3e-03 |
| `BoomCylinderOuter` | 48.822893 | 48.822893 | 1.1e-14 | **1.13e-02** |
| `BucketCylinderInner` | 9.0573779 | 9.0573779 | 1.6e-12 | 3.8e-03 |
| `BucketCylinderOuter` | 16.046199 | 16.046199 | −5.0e-13 | **−1.05e-02** |
| `StickCylinderInner` | 22.859457 | 22.859457 | 2.2e-13 | **1.45e-02** |
| `StickCylinderOuter` | 48.490352 | 48.490352 | −1.1e-14 | 9.3e-04 |

**This is the first volume check a CSG part has ever had that is not Monte-Carlo**, and it is nine
orders sharper than the one ROOT provides. It is not a substitute for `dV_sym` (which is stronger
still, and exact), but it is a cheap independent third opinion computed from a completely different
functional, and `Stream_I_Verdict.md` §10 explicitly asked for one.

**Two caveats that must travel with the table.**

1. `raster vs exact` on a box-like part is dominated by the **window excess**, not by quadrature.
   The window is the bounding box plus `--margin` (1e-3 cm by default, to cover the fact that a
   tessellated bounding box is *inscribed*), and for `box` the predicted mean excess over the three
   beams is 1.4450e-03 against a measured 1.4450e-03 — it is the entire error. Set `--margin 0`
   when the bounding box is known exactly and the volume becomes exact.
2. On thin, **tilted** parts the raster is much worse: Bagger's rams reach 2.0e-02, because a
   diagonal tube of small cross-section puts most of its cells on the silhouette. The per-beam
   spread is reported and is the honest error bar.

---

## 5. Was the instrument perturbed? Both gates and ctest, before and after

Gate totals and disagreement counts, quoted together as the contract requires.

| | before | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 80 cases, green | **86 cases, green** (6 appended, §8) |
| fixtures gate, shipped | 9/9 | **9/9** |
| fixtures gate, surface (historical) | 9/9 | **9/9** |
| fixtures scored | 9 of 10 leaf solids | **9 of 10** |
| Bagger gate, shipped | 12/12 | **12/12** |
| Bagger gate, surface (historical) | 9/12, the same three capacity values | **9/12, the same three, the same values** |
| Bagger scored | 12 of 13 leaf solids | **12 of 13** |
| unexplained disagreements, **surface**, fixtures | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| unexplained disagreements, **surface**, Bagger | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| disagreements, **shape**, fixtures / Bagger | 0/0/0/0 (2 parts) / 0/0/0/0 (7 parts) | **identical** |
| disagreements, **mesh**, Bagger | contains 418, distout 6721, distin 7973, safety 10299 (0/12 clean) | **identical** |
| `runOracleGate.py --self-test` | 17/17 | **17/17** |

`gate.json` was compared field by field with `timing*`, `*Seconds`, `nsPerCall` and `checksum`
stripped, on both corpora, before and after everything in this stream:

| corpus | leaf fields | removed | **changed** | added |
| --- | ---: | ---: | ---: | ---: |
| fixtures (9 parts) | 8525 | 0 | **0** | 0 |
| Bagger (12 parts) | 15195 | 0 | **0** | 0 |

**No leaf field added, removed or changed.** That is the expected outcome — nothing in the gate
path was touched — and *a comparison that cannot fail is not evidence*, so the differ was pointed
at a baseline with exactly two fields nudged (`nMismatchUnexplained` 0 → 1 and a
`capacityRelativeDeviation` at +1e-6 relative) and reported **exactly those two**.

The only build-system change is one new `o2_add_executable(xray …)` target;
`Detectors/Base/CMakeLists.txt` reconfigured cleanly and the documented incremental-`ninja`
rootcling failure did not recur.

---

## 6. The positive control, and the three defects the instrument found in itself

This is the most useful section of this document, and it is a failure report.

### 6.1 The first version reported the known quartic defect as clean — and the experiment was incapable

`Stream_E_Scale.md` §5 documents a real, unfixed defect: `solveQuarticReal` guards its branches
with `kTolerance` (1e-9 cm, a *length*) against quantities of dimension L² and L³, and at ×0.1 the
`torus_union_cyl` fixture **silently loses a ray it does cross**. The brief is explicit: if the
benchmark runs on that fixture and reports nothing wrong, the benchmark is broken.

It reported nothing wrong. 6912 rays, three axis beams, ×0.1: **6912/6912 identical, 0 LOST.**

*A refuted experiment is not a refuted hypothesis.* Before concluding anything, the instrument was
fed the **exact ray `Stream_E_Scale.md` §5.1 recorded** — origin (0.20943, 0.32925, 0.19348) cm,
direction (−0.75473, −0.031547, −0.65528) — plus a deterministic 4000-ray fan around it, against
the same ×0.1 `.brep`:

```
(a) shape API : 4001 rays, 5834 crossings | zero=0 stall=0 nonAdv=0 cap=0 unterm=0
    vs OCCT : 3994/4001 rays identical, LOST=14 extra=0 displaced=0 kind=0
    worst   : MISSING crossing at o=(0.208337, 0.318952, 0.18302) d=(-0.7662, -0.02475, -0.6421)
```

**14 lost crossings.** The instrument was capable; the raster had simply never sampled a dangerous
configuration.

### 6.2 Why: a parallel-beam raster is DIRECTION-POOR

Three axis beams are **three directions**, however many rays are fired. Impact-parameter density
and direction density are different resolutions, and this defect's trigger — the resolvent of
Ferrari's cubic falling below 1e-9 — depends on the *direction*. Tilting the beams off the
coordinate axes (`--tilt`) does not help either: it still gives three directions.

`--beams N` now fires **N Fibonacci-spiral directions**. With that, at 96 directions × 24 × 24:

| fixture scale | rays | rays identical | **LOST** | parity | verdict |
| --- | ---: | ---: | ---: | ---: | --- |
| ×1 (as shipped) | 55296 | 55296 / 55296 | **0** | 0 | clean |
| **×0.1** | 55296 | 55294 / 55296 | **4** | **2** | **defect detected** |

and at 256 directions × 32 × 32:

| fixture scale | rays | rays identical | **LOST** | parity |
| --- | ---: | ---: | ---: | ---: |
| ×1 (as shipped) | 262144 | 262144 / 262144 | **0** | **2** |
| **×0.1** | 262144 | 262140 / 262144 | **8** | **12** |

**The measured deviation, and it matches the documented mechanism exactly.** The lost crossings come
in **pairs** — 4 = 2 rays × (enter + exit), 8 = 4 rays × (enter + exit) — which is precisely
`Stream_E_Scale.md`'s "the ray silently misses a torus it does cross": the whole intersection is
discarded, not one root of it. The rate is 2 rays in 55296 (3.6e-05) and 4 in 262144 (1.5e-05) of
*all* rays, of which only a small fraction hit the torus at all; that is consistent with Stream E's
7.8e-04 per torus-hitting ray at ×0.1.

For contrast, the same fixture with three axis beams, at three tilts, 96 × 96 per beam = 27648 rays:

| tilt | ×1 | ×0.1 |
| --- | --- | --- |
| 0° | 0 LOST | **0 LOST** |
| 12° | 0 LOST | **0 LOST** |
| 27° | 0 LOST | **0 LOST** |

**Direction diversity, not ray count, is what finds it.** 27648 rays in 3 directions see nothing;
55296 rays in 96 directions see it. Anyone hunting a numerical defect with this benchmark must use
`--beams`.

**One further reading, at the shipped size.** At 262144 rays the ×1 fixture has **0 LOST but 2
parity mismatches** — `Contains()` at an interval midpoint contradicting a crossing list that
itself matches OpenCascade exactly, with `parityNB` = 0 so neither is excused by near-tangency. The
parity audit is therefore *more sensitive* than the crossing-list comparison on this defect, and it
fires at production scale. It is not localised further here; it is on the only toroidal fixture in
the corpus and it is reported, not fixed. **The quartic is not touched by this stream** — that is a
separate, deliberate decision, and `Stream_E_Scale.md` §5.3 already scopes what the fix has to be.

### 6.3 Two more defects, both in the benchmark, both caught by its own controls

* **The navigator world was too small for tilted beams.** A rotated lattice reaches outside the
  part's axis-aligned bounding box, so ray origins fell outside the world and `InitTrack` refused
  them — and the first version folded that into `iterationCapHits`, so it surfaced as **5358 LOST
  crossings at a 27° tilt in mode (b) only**. A misconfiguration of the benchmark wearing the
  costume of a geometry defect. It was caught *because* mode (a) was clean on the same rays, which
  is the entire argument for building both modes. The world is now derived from every ray of the
  raster, start to end, and `originOutsideWorld` is its own loud counter. It is **0 everywhere** in
  every table above.
* **The parity audit excused everything.** It asked `Safety()` with the state the crossing list
  *expects* rather than the state the shape *reports*, and `TGeoBBox::Safety(p, in=true)` goes
  negative at an outside point — so every real contradiction looked like near-tangency and the
  counter read 0 on a deliberately truncated list. Caught by its own negative control, which is why
  every positive check in both self-tests is paired with the mutation that must turn it red.

### 6.4 The synthetic controls

`o2-bench-detectorsbase-xray --self-test` — 19 checks, no database, no oracle:

* a box's two crossings and a hollow tube's **four along a diameter**, at analytic distances;
* the box chord integral exact at N = 7 and 32; the sphere envelope over N = 24…192;
* **the comparator against four injected defects** — a crossing moved by 1e-3 cm (must be
  `displaced` = 1 and `missing` = 0, *not* a lost/extra pair), a dropped crossing (`missing` = 1),
  an inserted one (`extra` = 1), a crossing with the wrong sense (`kindMismatch` = 1);
* **the parity audit against a truncated list and against an invented one**, both of which must
  fire.

`xrayOracle.py --self-test` — 11 checks: a box's chord, a hollow cylinder's four crossings, **a
tangent ray that must produce zero crossings** (the case that makes the midpoint classification
necessary rather than decorative), and the sphere volume at two densities.

---

## 7. Cost, and ALICE3

### The benchmark is far cheaper per part than the sample gate

Whole-corpus wall clock, single-threaded, one machine, 6912 rays per part, all representations,
both stepping modes, including the OCCT oracle:

| corpus | parts | rays scored | total wall | peak RSS |
| --- | ---: | ---: | ---: | ---: |
| fixture ladder | 9 | 62208 × (2 or 3 reps) × 2 modes | **10.7 s** | 316 MB |
| Bagger | 12 | 82944 × (2 or 3 reps) × 2 modes | **63.8 s** | 316 MB |

The OCCT side is the reason: **one `IntCurvesFace_ShapeIntersector` call yields every crossing along
a ray**, where `occtOracle.py` needs one call per sample. Measured: 10691 rays/s on the ladder and
1484 rays/s on Bagger, for 5.8 s and 55.9 s of oracle time respectively.

Kernel throughput, mode (a), `O2BVHSurfaceSolid`:

| | rays/s |
| --- | ---: |
| ladder, per part | 24 425 – 1 086 578 (median ~430 000) |
| Bagger, per part | 27 585 – 1 171 642 (median ~160 000) |
| mode (b) `TGeoNavigator` relative to mode (a) | **0.81 – 0.93 ×** (the navigator costs 8–20 % more) |

### ALICE3 — it runs, and it is the first instrument here that does

`scripts/geometry/ALICE_3_example/CAD_noETA.stp`, 55 solids / 9266 faces. **The X-ray benchmark
needs no tessellated mesh** — the raster is structured, so unlike `generateSamples()` nothing
rejection-samples through `O2Tessellated`. That removes the wall Stream E §8 documented (one 2 m
sphere reached 22.9 GB resident before being killed).

Converting **without `--mesh`**:

```bash
O2_CADtoTGeo.py scripts/geometry/ALICE_3_example/CAD_noETA.stp \
    --exact-surfaces auto --dump-brep --output-folder /tmp/J_alice3 -o geom.C
```

| | |
| --- | --- |
| wall clock | **24.2 s** |
| peak RSS | **581 MB** |
| exact sidecars written | **20** of 55 leaf solids (the rest have free-form faces — the known coverage ceiling) |
| `.brep` files written | 20 |
| output size | 22 MB |

**The tessellation wall is gone.** 24 s and 0.6 GB against "22.9 GB for one sphere". The X-ray
benchmark's precondition is the exact sidecar and the `.brep`, and neither needs a mesh.

Two of the 20 sidecars do not load (`ST1829909_004`, `ST1829909_01` — trim wires with 4.0e-06 and
5.4e-06 cm gaps against a 1e-06 cm join tolerance), the converter-side defect
`SolidNavigationHarness.md` already records. They are skipped and named; **18 parts are scored.**

#### The cost curve, measured

Mode (a) only, 18 parts, all three axis beams:

| raster N | rays | crossings | stepping | rays/s | load + `CloseShape` | total wall | peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 3 456 | 6 744 | 0.31 s | 11 327 | 14.1 s | **28.8 s** | 370 MB |
| 16 | 13 824 | 27 580 | 0.69 s | 20 152 | 14.0 s | **29.2 s** | 376 MB |
| 32 | 55 296 | 107 874 | 1.79 s | 30 857 | 14.1 s | **30.4 s** | 376 MB |
| 48 | 124 416 | 243 380 | 3.57 s | 34 831 | 13.9 s | **31.9 s** | 371 MB |

**The stepping is not the cost; loading is.** 14 s of `LoadSurfaceSolid` + `CloseShape` is a fixed
tax per run, of which 9.7 s is the single 965-patch `ST1829909_002`, and it dominates every
density up to N = 48. Ray throughput *improves* with N (fixed setup amortised), so raster density
on ALICE3 is nearly free: N = 48 costs 11 % more wall clock than N = 8 for 36 × the rays.

Full round trip including the OpenCascade oracle, 18 parts at N = 16 (13 824 rays), both stepping
modes: **82.3 s wall, 361 MB RSS.** For comparison, the sample-based oracle gate cannot run on this
model at all (no sample budget, no parallelism, and it requires a mesh).

**So yes: this is the first instrument in this project that runs on ALICE3.** Not a projection —
the table below is its output.

#### And it finds something the existing gate cannot

18 parts, N = 16, three axis beams, 768 rays per part, against the OCCT crossing-list oracle:

| | mode (a) shape | mode (b) nav |
| --- | ---: | ---: |
| rays identical | **13578 / 13822** | 13580 / 13822 |
| parts fully clean | **14 / 18** | 14 / 18 |
| LOST | **418** | 416 |
| extra | 0 | 0 |
| displaced | 0 | 0 |
| kind mismatch | **236** | 236 |
| worst \|Δt\| | 5.03e-08 cm | 1.76e-07 cm |
| unterminated | **70** | 68 |
| parity mismatch | **336** (parityNB 0) | 336 |
| zero / non-advancing / stalled / capped steps | **0 / 0 / 0 / 0** | 0 / 0 / 0 / 0 |
| **mode (a) vs mode (b)** | — | **6 of 13824 rays disagree** |

The four parts that are not clean:

| part | patches | rays identical | LOST | kind | unterm | parity | closure verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `ST0923290_013` | 20 | 700 / 768 | 68 | 68 | **68** | 132 | **reliable, navigable** |
| `ST0923290_018` | 48 | 696 / 768 | 144 | 72 | 0 | 96 | **reliable, navigable** |
| `ST0923290_019` | 44 | 672 / 768 | 192 | 96 | 0 | 108 | **reliable, navigable** |
| `ST2487462_01` | 80 | 760 / 768 | 14 | 0 | 2 | 0 | **reliable, navigable** |

**Every one of the 18 loaded ALICE3 parts is `reliable` and `navigable`** — checked directly
(`GetNavigationReliability`, `IsNavigable`, edge identity, zero non-manifold edges on all 18). So
these four are not the open-surface-set population: they are parts that **closure certifies as
watertight and that a transported particle nonetheless walks through**. That is precisely the
caveat `Tutorial.md` §5.4 attaches to edge identity — *"identity certifies that the source B-rep's
topology survived conversion. It does **not** certify each face's geometry — a mis-trimmed face
still reads closed"* — and this is the first measurement of it. The `kind` column (a crossing at
the right position with the **wrong sense**) points at face orientation or trim rather than at the
ray solver, and `ST0923290_013`'s 68 unterminated transports say a ray enters and never leaves on
8.9 % of its rays.

**Reported, not diagnosed, and not fixed here.** Three of the four are the same `ST0923290` family,
which is the natural place to start. Note also that the other 14 parts — including the 965-patch
`ST1829909_002` and the 236-patch `ST1829909_003` — are clean to 2.5e-13 and 6.9e-13 cm, so this is
a per-part defect and not a scale effect.

Two caveats on the ALICE3 numbers:

* **the raster window comes from the surface solid's conservative bounding box** (there is no mesh
  and no CSG shape to take a tight one from). Measured cross-section excess is 1.9e-04 to 3.4e-02,
  so **the volume column is not trustworthy on ALICE3** and is reported for completeness only. The
  crossing-list and robustness columns are unaffected: a larger window costs raster density, never
  correctness.
* 768 rays per part is a thin raster, and only three directions. The clean parts should be read as
  "no defect at this density", and §6.2 is the standing warning about direction poverty.

---

## 8. Unit tests

Six cases appended to `Detectors/Base/test/testBVHSurfaceSolid.cxx` in a delimited block ending
`// --- Stream J: X-ray / geantino transport benchmark ---`. `ctest -R BVHSurfaceSolid`
**80 → 86 cases, green.**

They include `XRayTransport.h` — **the same header the benchmark steps with**. A test written
against a second implementation of the same idea tests neither.

| case | what it pins |
| --- | --- |
| `XRayCrossingListsMatchClosedFormOnPrimitives` | a box's two crossings and a hollow tube's four, at analytic distances and senses; and zero stalls |
| `XRayCrossingListsAgreeBetweenBVHAndLoopOnAllFixtures` | **the transport-level BVH == `_Loop` guard**: whole ordered crossing lists, bit-identical, over an 11-direction fan on all seven kernel fixtures. Harder than the existing single-query version because every query after the first starts where the previous one landed. |
| `XRayCrossingComparatorCatchesInjectedDefects` | the comparator's positive and four negative controls, including that a *displaced* crossing is not reported as a lost one and that a sub-tolerance nudge is reported as nothing |
| `XRayParityAuditContradictsATruncatedList` | the parity audit against a truncated list and an invented one |
| `XRayChordIntegralIsExactForABoxAndConvergesForASphere` | the box exactness at three densities and the sphere's measured envelope |
| `XRayRasterRaysStartOutsideAndFansAreDirectionDiverse` | every ray starts and ends strictly outside the solid, the beam frames are orthonormal, and direction diversity is asserted as a number |

---

## 9. What this leaves open

* **Assembly-level transport is out of scope for this pass, deliberately**, and it is the obvious
  next step. What it would take: the world currently holds **one** `TGeoVolume`, so the crossing
  list is a two-state (in-part / out-of-part) sequence. An assembly needs (i) the converter's
  placement matrices carried into the manifest — they exist in `geom.C` but nothing machine-readable
  indexes them per part; (ii) the crossing list generalised from a boolean to a *volume id* per
  segment, which changes `Crossing`, the comparator and the oracle's answer format; (iii) an oracle
  that answers about the assembly rather than one leaf solid, i.e. `IntCurvesFace_ShapeIntersector`
  loaded with a compound and the classifier run per solid; and (iv) a new failure mode to count —
  **leaking**, where a track exits volume A and the navigator does not report it entering B. That
  last one is the whole point of doing it, and none of the three counters here can express it.
  Estimated: the oracle change is the small half; the volume-id crossing list is the real work.
* **The quartic is not fixed.** §6.2 is evidence, with a reproduction and a rate. The fix is scoped
  in `Stream_E_Scale.md` §5.3 and is a separate decision.
* **The parity mismatches at ×1 (2 in 262144 rays) are reported and not localised.** They are on the
  only toroidal fixture and they are the sharpest signal this instrument produced at production
  scale.
* **`O2Tessellated`'s `Contains()` disagrees with its own `DistFrom*`** — 534 intervals on
  Bagger/`BucketLink2`, 15 on the ladder. That is an internal inconsistency of one representation,
  measured here for the first time, and it belongs to whoever writes the `auto`-mode fallback
  policy.
* **The benchmark is single-threaded**, like the gate. Its cost per part is low enough that this has
  not mattered yet; on the full ALICE3 model with a dense raster it would.
* **No verdict, deliberately.** This is a report and an exit code, not a gate. Turning the crossing
  lists into an acceptance criterion means choosing a band for `displaced` per representation, and
  the mesh numbers say that band cannot be one number.
