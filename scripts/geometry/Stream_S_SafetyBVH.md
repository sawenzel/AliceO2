# Stream S — `Safety()` and `ComputeNormal()` learn about the BVH

Date: 2026-08-02. Follows [`Stream_P_RepresentationBench.md`](Stream_P_RepresentationBench.md) §2.3,
which found the defect and measured it. Companion to [`Tutorial.md`](Tutorial.md) (the three
representations) and [`Stream_J_XRay.md`](Stream_J_XRay.md) (the transport instrument).

**In one paragraph.** `O2BVHSurfaceSolid::Safety()` and `ComputeNormal()` ask the same question —
*which trimmed patch is nearest to this point* — and both answered it with a bare linear scan over
every patch, on a solid that has had a BVH since `CloseShape()` existed. On ALICE3's 965-patch
`ST1829909_002` that scan is **835 µs per call**; measured again here, it is the largest per-call
number anywhere in this project and it is nearly two orders of magnitude above the same solid's
`DistFromOutside`. Both kernels now walk the tree, pruning a node whose bounding-box distance
already exceeds the best patch distance found so far. **On `ST1829909_002` that is 835 305 ns →
22 376 ns, a 37.3× speedup, visiting 8.63 patches of 965.** On Bagger it is 1.7× to 5.4×. The
answers are **bit-identical** to the scan they replaced — 56 320 points on real parts from both
corpora, 18 009 in the unit tests, zero disagreements on either kernel — and a deliberately broken
pruning bound is shown to break that agreement, in the dangerous direction, on every fixture that
has anything to prune. Nothing else moved: both gates, both X-ray corpora at 96 beams and both
stepping modes reproduce their previous counters exactly, including the one pre-existing defect.
**The remaining cost is not the traversal, it is `distanceSqToPatch` itself** — a visited candidate
costs 3× to 4.6× the corpus-average patch, because the near patch is always the one that takes the
expensive trim-wire branch. That is the same function pair Stream P §2.2 named, and it is now the
whole of what is left.

Everything below is measured on `swenzel/bvhsurfacesolid` at `c716b25489`, single-threaded, on the
same 10-core aarch64 box as Stream P, with other tenants. Read §1 before any number.

---

## 1. What was wrong, and what the fix is

### 1.1 The defect

```cpp
double bestDistanceSq = std::numeric_limits<double>::infinity();
for (const auto& surface : fImpl->surfaces) {
  bestDistanceSq = std::min(bestDistanceSq, surface->distanceSqToPatch(testPoint));
}
return std::nextafter(std::sqrt(bestDistanceSq), 0.);
```

No BVH, no bound, no early exit; `ComputeNormal()` immediately below had the same body plus a
`closestSurface` pointer. Θ(patches) with a constant that is the full per-patch trim cost, on a
query `TGeoNavigator` and Geant4 both issue on every step.

### 1.2 The fix

An **ordered depth-first traversal with a running best**. Each stack entry carries the bound it was
pushed with and is re-tested when popped, so a node queued before a better patch turned up is still
dropped; children are visited nearer-box-first, so the best tightens as early as possible. No
priority queue and no heap allocation on the hot path — the stack is a reused `thread_local`
vector, the same convention the ray kernels' hit buffers already use.

`O2Tessellated::SafetyKernel` does the same job with a best-first priority queue and float
arithmetic, and it takes a shortcut this one deliberately does not: when the point is outside the
root box it returns the distance to that box. That is a *valid* safety and a much cheaper one, but
it is not the same number the scan returns, and this round's whole claim is that the number does
not change.

### 1.3 Which direction each bound errs, and why that is the entire safety argument

A valid safety may **never exceed** the true distance to the boundary: too little costs a navigator
one boundary computation it could have skipped, too much lets it step through a wall. The existing
`std::nextafter(std::sqrt(...), 0.)` already treats the two directions asymmetrically and is kept
verbatim. Four things could have pushed the answer up, and each errs down instead:

| quantity | errs | why |
| --- | --- | --- |
| the **patch distance** `distanceSqToPatch` | **down** | pre-existing and documented: on a wire-trimmed quadric it returns the distance to the *untrimmed* window, which contains the patch. Unchanged by this round |
| the **leaf box** | **down** (nearer) | `buildBVH` inflates each surface's `conservativeBounds` by `kBVHBoxTolerance = 1e-3` and rounds outward in float. A larger box is closer to the point |
| the **box-distance arithmetic** | forced **down** | computed in double from the float corners (which convert exactly), then scaled by `(1 - 1e-12)`. Three subtractions, three squarings and two adds carry ≲ 4 ulps ≈ 9e-16 relative, so the guard is three orders of magnitude of headroom. Scaling a lower bound down keeps it a lower bound |
| the **prune test** | conservative | `>=` for `Safety` (a node whose bound ties the best cannot improve a minimum); **strict `>`** for `ComputeNormal`, see §1.5 |

The load-bearing claim is the one in the middle rows: **the distance from the point to a node's box
is a lower bound on the distance to every patch in that subtree.** It holds because every
`distanceSqToPatch` in `BoundedSurface.h` is realised on the patch's *untrimmed* window — the wire
itself for the two planar families, the full rim band for a cylinder or cone, the full sphere, the
full torus — and each family's `conservativeBounds()` encloses exactly that window. Internal node
boxes in bvh2 are float min/max unions of their children's, which is exact, so the property lifts
from the leaves to the root.

Reading the code says that. Because a future surface family that broke it would make `Safety()`
silently answer *too much* and nothing else in the tree would notice,
`StreamS_PatchDistanceIsNeverBelowTheDistanceToItsOwnBoundingBox` now measures it directly: seven
concrete surfaces (polygon, disk, a B-spline-trimmed planar face, partial-sweep cylinder, cone,
sphere, torus), 4 000 points each in four distance regimes out to 1e6 cm, **28 000 checks, 0
violations**, asserted with the same 1e-12 relative guard the traversal applies and no tolerance
beyond it.

### 1.4 What did *not* change

Trim semantics, the sidecar format, every distance kernel, `Contains`, the ray tmax pruning, and
`Safety`'s returned value. The `Bool_t in` argument stays ignored, as before. Before `CloseShape()`
there is no tree, and both kernels fall back to the scan exactly as the distance kernels do.

### 1.5 The tie-break, which is not decoration

The scan uses a strict `<`, so the **lowest-indexed** patch wins an exact tie. Exact ties are
routine, not exotic: the centre of a box is exactly equidistant from all six faces, and so is the
axis of a hollow tube from its inner wall and both caps. A traversal in BVH order would legitimately
pick a different patch there and `ComputeNormal` would return a different normal — an answer that
depends on how the tree happened to be built.

There are callers that would see it. `scripts/geometry/probes/faceNormals.cxx` compares the kernel's
outward normal against a per-face reference, and `O2SolidHarness.cxx` asks for a normal at sampled
points; both read the *identity* of the chosen face out of the normal. So the tie-break is
preserved explicitly: the traversal carries the index, prefers the lower one at an exact tie, and
therefore may **not** prune a node whose bound merely *equals* the current best. That is the only
place where `ComputeNormal` prunes less than `Safety`, and it is measured — on the 144-patch ring
fixture the two visit the *same* 3.16 candidates per call, and across all 16 real parts in §4
`ComputeNormal` sits within 4 % of `Safety` and scatters both ways, i.e. the concession is below the
noise floor rather than merely small.

---

## 2. Correctness: exact equality against the scan it replaced

The scan is kept, as `Safety_Loop()` and `ComputeNormal_Loop()`, following the `_Loop` convention
of `Contains_Loop` / `DistFromOutside_Loop`. (Those are public, so these are too; the brief asked
for a private test-only reference, and matching the existing convention won.) The contract is
**exact equality**, not agreement to a tolerance: both minimise the same `distanceSqToPatch` over
the same patches under the same tie-break, so any difference at all is a traversal or pruning bug.

| where | points compared | `Safety` disagreements | `ComputeNormal` disagreements |
| --- | ---: | ---: | ---: |
| **Bagger, all 13 parts** (`--perf`, 4096 points each) | **53 248** | **0** | **0** |
| **ALICE3**, `ST1829909_002` / `_003` / `ST0923290_019` (1024 each) | **3 072** | **0** | **0** |
| unit tests, 9 synthetic fixtures × 2 001 points (each compared under 2 `Safety` and 2 `ComputeNormal` calls) | **18 009** | **0** | **0** |
| **total** | **74 329** | **0** | **0** |

`Safety` is compared with both `in` labels and `ComputeNormal` both with and without a direction
(the direction flips the chosen patch's normal, and a wrong patch can hide behind that flip). The
worst signed gap `accelerated − brute force` over the unit-test sample is asserted to be
**identically 0.**, so this is the strong claim and not the weak "accelerated ≤ brute force" one.

The unit-test sample is deliberately loaded rather than uniform, because pruning is trivially
correct where one patch is far nearer than the rest and only bites where several compete: one fifth
of the points are walked **onto** the surface along their own normal (distance ≈ 0), one fifth to
1e-9 cm off it, one fifth well outside the bounding box, one fifth at 1e7 cm where the node bound is
large and its rounding is worst, one fifth generic — plus the exact centre of each solid, which is
the tie. The nine fixtures are the seven of the existing `_Loop` cross-validation, a 12-box ring
(72 patches, heavily overlapping boxes), and a cylinder carrying four **B-spline** trim windows,
which is the trim family Stream P priced at 2–6 µs per candidate patch.

### 2.1 The negative control, and it fires

A cross-check that cannot fail has not passed. `SetSafetyBoundUnsoundForTest(true)` replaces the
node bound with the distance to the box **centre**, which is larger than the box distance for any
box with extent and therefore bounds nothing. It must break the comparison, and it must break it in
a stated direction.

```
sabotaged bound caught on safetyBox:         995 disagreements over 401 points
sabotaged bound caught on safetyTube:       1587 disagreements over 401 points
sabotaged bound caught on safetyHollowTube: 1345 disagreements over 401 points
sabotaged bound caught on safetyCone:       1696 disagreements over 401 points
sabotaged bound NOT caught on safetySphere  (1 patch)
sabotaged bound NOT caught on safetyTorus   (1 patch)
sabotaged bound caught on safetyCapsule:     529 disagreements over 401 points
sabotaged bound caught on safetyRing:       1915 disagreements over 401 points
sabotaged bound caught on safetyWireTrim:    510 disagreements over 401 points
```

**7 of 7 fixtures that have anything to prune are caught.** The two that are not are the sphere and
the torus, which are *single-patch* solids: their BVH is one leaf, every bound sound or not is
applied to a node that must be entered anyway, and there is nothing for a bad bound to discard. That
is not a hole in the control, it is the statement that the control tests pruning — so the test
excludes them by the explicit criterion `GetNsurfaces() > 1` and asserts that the uncaught ones are
exactly the one-patch solids, rather than fudging the expected count.

And the failure direction is asserted: over the whole sabotaged sample, **every** disagreeing
`Safety` came out **larger** than the reference (`safetyTooLarge > 0`, `safetyTooSmall == 0`). Too
large is the failure mode that lets a navigator walk through a wall, and it is the one the sabotage
reproduces — so the healthy path is watched by a test that has been shown to see the thing it
exists to see. Turning the sabotage off restores exact agreement on the same sample, so the
disagreements belong to the bound and not to the fixtures.

### 2.2 The independent oracle also still agrees

The gate scores `Safety()` against OpenCascade, which is a reference the `_Loop` twin cannot supply:

```
fixtures gate: safety(cand) scored=5000  agree=100.00%  mismatch(band=0,missed=0,unexplained=0)  worstDev=0
               safety(ref)  scored=5000  agree=100.00%  mismatch(band=0,missed=0,unexplained=0)  worstDev=0
```

---

## 3. Nothing moved anywhere else

Re-run in full. Totals **and** disagreement counts, together.

### 3.1 Gates

| gate | scored / pass | `surface` contains/distout/distin/safety | `shape` | exit |
| --- | --- | --- | --- | ---: |
| fixtures | **10/10 scored parts pass** | **0 / 0 / 0 / 0** on **10/10 parts** | **0 / 0 / 0 / 0** on **2/2** | **0** |
| Bagger | **13/13 scored parts pass** | **0 / 0 / 0 / 0** on **13/13 parts** | **0 / 0 / 0 / 0** on **7/7** | **0** |

Bagger's separate `10/13 pass on the surface representation` line is unchanged and unrelated: it is
the **capacity** criterion on three rams that ship as `shape` (`BoomCylinderOuter`,
`BucketCylinderInner`, `StickCylinderOuter`, each 1.3–1.4e-06 relative), reported and not gated.
Bagger ships no tessellated part: CSG 7, exact surfaces 6, of 13.

### 3.2 X-ray, `--beams 96`, both stepping modes, all robustness counters

`--raster 16 --beams 96` = 24 576 rays per part; 10 fixture parts and 13 Bagger parts.

| corpus | repr | mode | rays identical | LOST | extra | displaced | worst \|Δt\| cm | parts clean |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| fixtures | `surface` | (a) shape | 245 758 / 245 760 | **0** | **2** | **0** | 5.764e-13 | 9/10 |
| fixtures | `surface` | (b) nav | 245 758 / 245 760 | **0** | **2** | **0** | 5.764e-13 | 9/10 |
| fixtures | `shape` | (a) / (b) | 49 152 / 49 152 | **0** | **0** | **0** | 1.146e-13 | 2/2 |
| Bagger | `surface` | (a) shape | **319 488 / 319 488** | **0** | **0** | **0** | 4.571e-10 | **13/13** |
| Bagger | `surface` | (b) nav | **319 488 / 319 488** | **0** | **0** | **0** | 4.573e-10 | **13/13** |
| Bagger | `shape` | (a) / (b) | 172 032 / 172 032 | **0** | **0** | **0** | 8.436e-11 | 7/7 |

Robustness, every counter:

```
fixtures surface (a) steps=599870  zeroStep=0 noAdv=0 unstick=0 capHit=0 unterm=2 oddList=2
                                   nonAlt=0 dupXing=0 parity=2 parityNB=0 noTrans=0 outWorld=0 orgIn=0
fixtures surface (b) steps=599870  ... identical ...
Bagger   surface (a) steps=612616  zeroStep=0 noAdv=0 unstick=0 capHit=0 unterm=0 oddList=0
                                   nonAlt=0 dupXing=0 parity=0 parityNB=0 noTrans=0 outWorld=0 orgIn=0
Bagger   surface (b) steps=612616  ... identical ...
Bagger   shape   (a)/(b) steps=285178  every counter 0
```

The fixtures' `2 extra / 2 unterminated / 2 parity` is **`cyl_inter_cyl`, and it is the
pre-existing defect Stream P §4 reported at exactly these counts** — a `surface`-representation
disagreement that appears at 96 Fibonacci directions and not at 3 axis beams. It reproduces
unchanged, on the same part, at the same count. It is not touched by this round and is not
diagnosed here.

### 3.3 The instrument's own controls

`o2-bench-detectorsbase-xray --self-test`: **34 checks, 0 failures** (unchanged).
`ctest -R BVHSurfaceSolid`: **102 → 108 cases, 291 919 assertions, green.**

---

## 4. The speedup, with the instrument that found the defect

`--perf` writes `safetyBVHNs` / `safetyLoopNs` / `normalBVHNs` / `normalLoopNs` /
`bvhCandidatesPerSafetyCall` / `safetyDisagreements` / `normalDisagreements` into the surface
solid's `localise` block. **Before and after are the same binary, the same sample set and the same
pass structure**, because the `_Loop` twin *is* the code Stream P measured, verbatim — which is a
stronger comparison than two runs on a shared machine, and it is why §4.1's `_Loop` column
reproduces Stream P's `Safety` column to a few per cent.

### 4.1 Bagger, 13 parts — `--perf --raster 16 --beams 8`, 4096 points, 9 passes

| part | patches | `Safety` **before** (`_Loop`) | `Safety` **after** | **speedup** | candidates/call | `ComputeNormal` before → after |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `BasePin` | 3 | 323.9 | **103.6** (2.1 %) | **3.1×** | 1.15 / 3 | 333.6 → 106.2 |
| `Base` | 44 | 4 886.7 | **901.8** (7.9 %) | **5.4×** | 3.38 / 44 | 5 032.1 → 895.2 |
| `BoomCylinderInner` | 6 | 581.3 | **167.0** (22 %) | **3.5×** | 1.34 / 6 | 585.9 → 167.8 |
| `BoomCylinderOuter` | 8 | 669.6 | **223.6** (3.8 %) | **3.0×** | 2.29 / 8 | 676.5 → 221.2 |
| `Boom` | 31 | 3 695.8 | **2 098.3** (138 % ‡) | **1.8×** | 5.00 / 31 | 3 752.5 → 2 013.8 |
| `BucketCylinderInner` | 6 | 586.2 | **174.2** (6.0 %) | **3.4×** | 1.37 / 6 | 593.9 → 175.5 |
| `BucketCylinderOuter` | 10 | 920.5 | **228.2** (7.6 %) | **4.0×** | 2.26 / 10 | 939.1 → 234.1 |
| `BucketLink1` | 22 | 2 618.4 | **849.0** (9.0 %) | **3.1×** | 3.04 / 22 | 2 655.7 → 838.9 |
| `BucketLink2` | 23 | 3 557.4 | **1 060.7** (3.7 %) | **3.3×** | 3.86 / 23 | 3 610.3 → 1 065.2 |
| `Bucket` | **97** | 9 890.6 | **2 403.9** (3.7 %) | **4.2×** | 5.14 / 97 | 9 939.9 → 2 359.0 |
| `StickCylinderInner` | 6 | 600.1 | **160.4** (10 %) | **3.8×** | 1.33 / 6 | 598.4 → 164.2 |
| `StickCylinderOuter` | 8 | 675.6 | **236.7** (10 %) | **2.9×** | 2.34 / 8 | 676.3 → 235.2 |
| `Stick` | 24 | 2 783.3 | **1 397.9** (2.8 %) | **2.0×** | 3.95 / 24 | 2 832.9 → 1 385.0 |

ns per call, median of 9 passes; the parenthesised figure is the pass-to-pass spread of the *after*
column. ‡ `Boom`'s 138 % is one bad pass on a shared machine, not a property of the kernel: the
median is 2 098 ns and the independent `localise` timing of the same call in the same run is
2 049 ns.

### 4.2 ALICE3 — where the 812 µs was

1024 points, 5 passes, sidecars from a fresh `--exact-surfaces auto` conversion without `--mesh`.

| part | patches | `Safety` before | `Safety` after | speedup | candidates/call | `ComputeNormal` before → after |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| **`ST1829909_002`** | **965** | **835 305** (Stream P: 811 870) | **22 376** (1.2 %) | **37.3×** | **8.63 / 965** | 837 217 → **22 286** |
| `ST1829909_003` | 236 | 251 883 (Stream P: 246 376) | **40 819** (1.8 %) | **6.2×** | 4.43 / 236 | 250 636 → 40 713 |
| `ST0923290_019` | 44 | 31 348 (Stream P: 31 885) | **6 608** (6.0 %) | **4.8×** | 3.37 / 44 | 31 287 → 6 671 |

**The headline number of Stream P §2.3 is gone.** `Safety()` on `ST1829909_002` was 94× that solid's
`DistFromOutside`; it is now **2.5×**, and it is no longer the largest per-call figure on the part —
transport is (20 972 ns/ray).

### 4.3 What did not move: the other three kernels and transport

Not touched, and measured to confirm it. Against Stream P §4's Bagger table, on the same four
kernels the benchmark also times (`--perf` passes `safe = nullptr`, so the distance kernels do not
call `Safety` here):

| part | `Contains` P → now | `dOut` P → now | `dIn` P → now | ns/ray P → now |
| --- | --- | --- | --- | --- |
| `BasePin` | 129 → 135 | 211 → 217 | 142 → 146 | 404 → 411 |
| `Base` | 476 → 489 | 642 → 652 | 626 → 645 | 1 559 → 1 585 |
| `BoomCylinderOuter` | 2 633 → 2 710 | 6 625 → 6 795 | 8 468 → 8 698 | 10 105 → 10 403 |
| `Bucket` | 639 → 655 | 1 139 → 1 128 | 785 → 780 | 2 135 → 2 127 |
| `StickCylinderOuter` | 2 676 → 2 756 | 6 693 → 6 851 | 8 444 → 8 562 | 9 834 → 10 107 |
| `Stick` | 701 → 710 | 623 → 636 | 886 → 903 | 1 723 → 1 763 |

**Everything is up by 0 – 3 %, uniformly, across all four columns and all thirteen parts** — the
signature of ambient load on a shared box between two runs a day apart, not of a code change. Two
things make that reading solid rather than convenient: no kernel moves *more* than any other, and
the `Safety_Loop` column of §4.1, which is the *old* `Safety` code running in the *new* binary,
lands within a few per cent of Stream P's `Safety` column on every part. If the build had got
slower, that column would have moved with the rest; it did not.

**Reported as a finding, not measured here:** `DistFromOutside`/`DistFromInside` *do* call `Safety()`
when a navigator passes a non-null `safe`, which is how `TGeoNavigator` uses them. The benchmark
passes `nullptr`, so that path is invisible to every number above. Arithmetically, on
`ST1829909_002` a `DistFromOutside` with safety was 8 912 + 835 305 = **844 µs** and is now
8 912 + 22 376 = **31 µs**, i.e. **27×** — but that is composition of two measured numbers, not a
measurement, and nothing in this project times that call as a navigator makes it.

### 4.4 Pruning statistics, and where the remaining cost is

Patches handed to `distanceSqToPatch` per `Safety` call. **Before** is `GetNsurfaces()` by
construction — the scan visits every patch, always.

| part | patches (before) | after | patches avoided | speedup | **ns per patch, scan** | **ns per candidate, BVH** | ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Bagger/Boom` | 31 | 5.00 | 6.2× | 1.8× | 119 | 420 | 3.5× |
| `Bagger/Stick` | 24 | 3.95 | 6.1× | 2.0× | 116 | 354 | 3.1× |
| `Bagger/Bucket` | 97 | 5.14 | 19× | 4.2× | 102 | 468 | 4.6× |
| `Bagger/Base` | 44 | 3.38 | 13× | 5.4× | 111 | 267 | 2.4× |
| `ALICE3/ST1829909_002` | 965 | 8.63 | **112×** | 37× | 866 | 2 593 | 3.0× |
| `ALICE3/ST1829909_003` | 236 | 4.43 | 53× | 6.2× | 1 067 | 9 214 | 8.6× |

**The speedup is always well below the patch-count reduction, and the last three columns say why:
the patches the traversal keeps are the expensive ones.** A visited candidate costs 2.4× to 8.6× the
corpus-average patch, because the nearest patch is precisely the one whose point falls outside its
trim and takes the branch that walks the flattened trim wire, while a far patch takes the cheap
in-domain branch. Pruning removes patches at the average cost and leaves the tail.

That bounds what any further work on the traversal can buy — at 4.43 candidates of 236 on
`ST1829909_003` there is at most another 4× available from pruning *ever*, and the eight remaining
candidates on the 965-patch part are 9 µs of genuine geometry. **The next lever is
`distanceSqToPatch` itself**, i.e. `Curve2D::closestPoint` — the same function Stream P §2.2 found
holding 61 % of cycles for a different reason. Two things about it are visible from here and were
not before: it is entered on the *near* patch, where its cost is worst, and it is now called 5–9
times per `Safety` instead of 100–1000, so its per-call cost is the whole number rather than a term
in it.

---

## 5. What this invalidates elsewhere

Recorded here rather than edited into the central documents, per the working agreement. **No table
in `Stream_P_RepresentationBench.md` is edited**; it is the baseline this stream measures against.

| document | claim | now |
| --- | --- | --- |
| `Stream_P_RepresentationBench.md` §2.3 (**do not edit — this is the baseline**) | *"`Safety()` has no acceleration structure at all … 811 870 ns on `ST1829909_002` … This is the single largest per-call number anywhere in this document and it is a structural property of one function"* | fixed. 835 305 → **22 376 ns**; 94× `DistFromOutside` → 2.5×. The §2.3 table, and the `Safety` column of the §4 and §6.2 tables, are superseded by §4 here. Everything else in that document stands and was re-measured unchanged (§4.3) |
| `Stream_P_RepresentationBench.md` §7, ALICE3 row | *"the exact path stays, and its 9.75 s `CloseShape` and 812 µs `Safety()` become the thing to fix"* | half done. The 812 µs is fixed; **`CloseShape` is untouched** and is now the largest remaining ALICE3 number by three orders of magnitude |
| `Stream_P_RepresentationBench.md` §6.2, ALICE3 medians (`surface`: `Safety` 11 994 ns over 18 parts) | | superseded for the `Safety` column only; the three parts re-measured here fall 4.8× to 37× |
| `Tutorial.md` §5.6 / `NEXT.md` | nothing in either mentions `Safety`'s cost | **no claim invalidated.** Listed so the reconciliation does not have to go looking |
| `NEXT.md`, "Where the branch stands" | `ctest -R BVHSurfaceSolid` **97 cases** | **108** (Stream P took it to 102) |

---

## 6. What a reviewer should challenge

In the order I would attack this.

1. **The lower-bound argument rests on a property of `distanceSqToPatch` that no interface
   enforces.** §1.3's claim — every patch distance is realised inside the surface's own
   `conservativeBounds` — is true of all six families today, is now measured on all six, and is
   written down in the traversal's doc comment. It is still not a *contract*: `BoundedSurface`
   declares `conservativeBounds` and `distanceSqToPatch` as independent virtuals with no stated
   relation between them. A seventh family could satisfy both docstrings and break `Safety`. The
   right fix is a sentence in the base class saying the two are coupled, and I did not write it
   because it is a change to a header this round promised not to touch beyond its own additions.
2. **The negative control proves the bound is load-bearing; it does not prove *this* bound is
   sound.** Breaking the bound is caught. That establishes the test can see a wrong bound — it does
   not, on its own, establish that no *other* wrong bound would slip through. What carries that is
   §1.3's argument plus the 28 000-check invariant case plus 74 329 points of exact agreement. A
   reviewer who wants more should ask for the sabotage to be parameterised (e.g. a bound scaled by
   1 + ε for shrinking ε) to find where the sample stops noticing.
3. **`ComputeNormal`'s tie-break is preserved on the strength of two callers I found by grep.** I
   argued in §1.5 that `faceNormals.cxx` and `O2SolidHarness.cxx` read the identity of the chosen
   face out of the normal. Neither has a test that would fail if the choice changed, so my evidence
   that the tie-break matters is a reading of intent, not a red test. Preserving it is free (measured
   at < 1 % on 16 parts), so the decision is cheap either way — but the *reason* is softer than the
   rest of this document.
4. **The before/after comparison uses the `_Loop` twin, not the previous commit.** §4's "before"
   column is the old code running inside the new binary. That is deliberate and I argue in §4 it is
   stronger. It is also not the same thing as checking out `fe2339bdfd` and re-running, which I did
   not do. If the two paths shared a cache effect or an inlining decision, the twin would
   under-report the win. The check against Stream P's independently-measured column (§4.3) is what
   bounds that, and it is a few per cent.
5. **Every number is warm-cache, one machine, one run, other tenants.** Same caveat as Stream P §1,
   and one column (`Boom`, 138 % spread) visibly took a hit from it. The medians are over 9 passes
   (5 on ALICE3) and the spreads travel with them, but nothing here is reproduced on a second
   machine and the 0–3 % uniform drift in §4.3 is exactly the size of the effect.
6. **The unit-test sample is loaded towards the surface by construction, using `Safety_Loop` and
   `ComputeNormal_Loop` to place the points.** That is circular in a narrow sense: the oracle is
   used to build the sample the oracle is then checked against. It is not circular in the sense that
   matters — the *comparison* is still accelerated against brute force at a point, whatever produced
   the point — but a sample built from the accelerated path would test a different set of points,
   and I did not build one.
7. **ALICE3 is measured on 3 of 20 parts.** The three that Stream P §6.2 named as the worst rows.
   The other 15 are covered for *correctness* (§2's zero-disagreement claim covers only 3 of them
   too) but not for cost. A full 20-part `--perf` on ALICE3 costs ~15 s of `CloseShape` per part and
   I judged three worst cases plus 13 Bagger parts enough; a reviewer may not.
8. **`kBVHBoxTolerance = 1e-3` is a large inflation for a bound.** Every leaf box is 1e-3 cm bigger
   than its patch in each direction, so the bound is that much looser than it needs to be, and on
   sub-millimetre features that is a real fraction of the box. It costs pruning, never correctness,
   and it is the existing constant the ray path already uses — but nobody has asked whether the
   *distance* queries want the same slack the *ray* queries do.

---

## 7. Reproducing it

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid o2-bench-detectorsbase-xray
ctest -R BVHSurfaceSolid                       # 108 cases
$B/stage/tests/o2-test-detectorsbase-BVHSurfaceSolid --run_test='StreamS_*' --log_level=message
$B/stage/bin/o2-bench-detectorsbase-xray --self-test                      # 34 checks

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/S_fix --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/S_bag \
    --model scripts/geometry/STEP_examples/Bagger.step
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/S_xbag \
    --reuse-db /tmp/S_bag/db --beams 96 --raster 16

# the speedup: the `localise` block of every surface row now carries safetyBVHNs / safetyLoopNs /
# safetySpeedup / bvhCandidatesPerSafetyCall / safetyDisagreements / normalDisagreements
$B/stage/bin/o2-bench-detectorsbase-xray --db /tmp/S_bag/db --perf --raster 16 --beams 8 \
    --json /tmp/perf_bagger.json

# ALICE3, one sidecar at a time (no db needed)
SW=$HOME/alisw/sw/ubuntu2404_aarch64    # convert ALICE3 WITHOUT --mesh; see NEXT.md
$B/stage/bin/o2-bench-detectorsbase-xray --surfaces /tmp/S_a3/surfaces_ST1829909_002_0_1_1_36.bin \
    --perf --raster 6 --beams 3 --perf-points 1024 --perf-rays 1024 --perf-passes 5
```
