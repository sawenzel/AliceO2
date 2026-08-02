# Stream P — what each representation costs, per part, per call

Date: 2026-08-02. Companion to [`Tutorial.md`](Tutorial.md) (the three representations),
[`Stream_J_XRay.md`](Stream_J_XRay.md) (the transport instrument this extends) and
[`MeshHealing.md`](MeshHealing.md) (why a good accuracy number is not a safety statement).

**In one paragraph.** The question was: per part, side by side, what do CSG, `O2BVHSurfaceSolid`
and `O2Tessellated` cost to *ask* and to *hold*, how wrong is the tessellation, and where is the
knee of the mesh-precision curve. The answer inverts the assumption the project has been carrying.
**`TGeoCompositeShape` is not slow — it is by far the fastest thing here**, 43× to 82× faster than
the exact surface solid on the six Bagger rams, and the preliminary 77× observation on
`StickCylinderOuter` reproduces at **72× under warm-up, repetition and load exclusion**. The
expensive representation is the one we thought was the cheap one: **`O2BVHSurfaceSolid` costs 130 ns
to 13 µs per query, and the spread is not patch count but trim-wire complexity.** 76 % of its cycles
are in one function pair — `Curve2D::closestPoint` (61 %) and `Curve2D::rightwardCrossings` (15 %),
both called from `curveTrimContains` — and `Safety()` has no acceleration structure at all, which on
ALICE3's 965-patch solid is **812 µs per call**. On ALICE3 the coarse mesh beats the exact solid on
cost on every axis at once — **14× faster in transport, 1.7× smaller, 5× faster to load** — and then
**loses 160 walls**, so the exact path stays and its cost becomes the thing to fix. The mesh-precision
sweep says the user's guessed 0.5 sits on the wrong side of the knee by ~60× in accuracy, that
`--mesh-prec` is an angular knob in disguise, and that refining past 0.1 makes the one broken mesh in
the corpus **worse**, not better.

Everything below is measured on `swenzel/bvhsurfacesolid` at `298a08ac79`, single-threaded, on a
10-core aarch64 box with other tenants. **Read §1 before any number.**

---

## 1. How to read these numbers, and what they are not

Six statements that qualify every table in this document. They are here rather than in a footnote
because each of them changes what a ratio means.

1. **Warm cache, and that is a ceiling not a prediction.** Every kernel is warmed and then timed
   over 9 complete passes over the same sample set; the reported figure is the **median** and the
   `[min .. max]` spread travels with it. Every solid measured here fits in cache — the largest is
   13 MB and most are under 1 MB. A real simulation holds thousands of solids and misses to DRAM,
   and the representation with the larger working set (the mesh: 2.5 MB where the surface solid is
   0.2 MB on ALICE3 *per part*, but 5.4 MB vs 0.2 MB summed over Bagger at `--mesh-prec 0.05`)
   loses more there than these numbers show. **The mesh's advantage is an upper bound. The
   composite's is not** — a `TGeoCompositeShape` of two tubes is a few hundred bytes and will still
   be resident when everything else has been evicted.
2. **The same questions from the same sample sets.** Per part, one set of 4096 points and
   2 × 4096 rays is built once, from a *named* reference representation's own `Contains()`, and
   handed unchanged to all three. The name is in every JSON record (`partitionedBy`) and it is
   `surface` for every part of both corpora. Letting each representation partition its own
   inside/outside set would time `DistFromInside` on three different sets and call the ratio a
   speed comparison.
3. **Load is excluded from the per-call cost and reported separately**, split into the file read and
   `CloseShape()`. That split is not cosmetic: on ALICE3's `ST1829909_002`, loading is 0.012 s and
   `CloseShape` is **9.75 s**, so lumping them would price an acceleration structure as a file
   format. `Stream_J_XRay.md` §7's "loading dominates" is confirmed and now localised: it is not
   loading, it is closure.
4. **`Safety()` is asked with a fixed in/out label** from the reference partition, not with each
   shape's own `Contains()`. Otherwise a shape that disagreed about a point would be timed on a
   different branch, and the column would price two kernels as one.
5. **Every distance timing carries its hit fraction.** A `DistFromOutside` column measured on rays
   that all miss prices the early-out. Outside rays are aimed into the bounding box; hit fractions
   run 40 %–100 % and are in the JSON.
6. **Single-threaded, deliberately**, like the gate and like `Stream_J_XRay.md`, so the numbers are
   comparable with everything else in this project.

### The instrument's own controls

`NEXT.md`: *a refuted experiment is not a refuted hypothesis.* A timing harness that cannot tell a
slowed shape from a fast one is not measuring what it claims, so `o2-bench-detectorsbase-xray
--self-test` went **17 → 37 checks** and every measurement family got a negative control:

| control | what it does | what it caught |
| --- | --- | --- |
| `BallastShape` | a `TGeoBBox` that answers identically and burns cycles; the ratio must exceed 2 **on all four kernels** | nothing — it passes, which is the point |
| checksum + pass count | the timed loop must not be elided | nothing |
| 64 MB touched allocation | both memory columns must see it; the heap one must see the release | **`mallinfo2().uordblks` moved by exactly zero.** glibc services anything over the mmap threshold into `hblkhd`. The heap column was blind to 64 MB until `hblkhd` was added |
| ladder structure | a fixture claiming 32 leaves must hold 32, and chain/balanced must differ in depth | nothing |
| ladder scaling | K=32 must cost measurably more per `Contains()` than K=2 | nothing |
| sample partition | the labels must agree with the reference that made them, and rays must hit | nothing |

`ctest -R BVHSurfaceSolid`: **97 → 102 cases, green.** The five new ones exercise
`RepresentationBench.h`, the same header the benchmark times with.

---

## 2. The headline: composites are the fast ones

### 2.1 The six Bagger rams, all three representations, same rays

Transport is the composite number a simulation actually pays. `--raster 16 --beams 8`, 2048 rays,
mode (a).

| part | surface ns/ray | mesh ns/ray | **shape ns/ray** | surface / shape | mesh / shape |
| --- | ---: | ---: | ---: | ---: | ---: |
| `BasePin` (1 `TGeoTube`) | 404 | 448 | **68** | 5.9× | 6.6× |
| `BoomCylinderInner` | 4 630 | 796 | **107** | **43×** | 7.4× |
| `BoomCylinderOuter` | 10 105 | 1 624 | **134** | **75×** | 12× |
| `BucketCylinderInner` | 4 159 | 669 | **106** | **39×** | 6.3× |
| `BucketCylinderOuter` | 9 909 | 1 634 | **143** | **69×** | 11× |
| `StickCylinderInner` | 4 634 | 796 | **108** | **43×** | 7.4× |
| `StickCylinderOuter` | 9 834 | 1 649 | **137** | **72×** | 12× |

**The 77× observation holds.** On `StickCylinderOuter` it reproduces at **71.7×**, with warm-up, load
excluded, 9 repeated passes and a spread of 1.6 % (surface) and 8.9 % (shape) — and the two
representations produce the **identical 1596 crossings** over the same 2048 rays, so this is not one
of them doing less work. The single-shot kernels say the same thing and say it harder:

| `StickCylinderOuter` | `Contains` | `Safety` | `DistFromOutside` | `DistFromInside` |
| --- | ---: | ---: | ---: | ---: |
| `surface` (8 patches) | 2 676 ns | 678 ns | 6 693 ns | 8 444 ns |
| `mesh` (2244 tri) | 704 ns | 5 902 ns | 1 685 ns | 1 269 ns |
| **`shape`** (2-leaf union) | **19 ns** | **23 ns** | **82 ns** | **62 ns** |
| surface / shape | **141×** | 29× | **82×** | **136×** |

A recognised *plain* primitive is cheaper still: `BasePin` as a `TGeoTube` is **1.9 ns** for
`Contains` and **68 ns/ray** in transport. The `box` fixture as a `TGeoBBox` is 2.5 ns and 64 ns/ray.

### 2.2 Why — and it is not the BVH

The BVH is doing its job. On these parts it hands **1.0 to 2.7** patches to the leaf callback per
`DistFromOutside` call out of 6–44 patches, and switching its ray-`tmax` pruning off costs 1.1× to
2.0×. On ALICE3's 965-patch solid, the non-BVH `_Loop` twin is **442 µs** against the BVH's
**9.0 µs — a 49× win**. The acceleration structure is not the problem.

The problem is what one candidate patch costs. Under `perf record` on `StickCylinderOuter`, with the
surface representation only:

| symbol | share of cycles |
| --- | ---: |
| `surface::Curve2D::closestPoint` | **61.1 %** |
| `surface::Curve2D::rightwardCrossings` | **15.2 %** |
| `surface::measureRimClosure` (closure, i.e. `CloseShape`, not a query) | 9.1 % |

Both of the first two are called from **`curveTrimContains`** — the 2D point-in-trim-wire test.
Reading the code confirms the shape of it: `CurveWire::classify` first asks
`curve.distanceSq(point)` of **every curve of the wire, unconditionally**, to decide whether the
point is inside a boundary band, and only then counts crossings; and for a B-spline curve
`Curve2D::closestPoint` walks the **entire cached flattened polyline**. So the cost of one patch is
Θ(total flattened trim-wire length), it is paid on every query whether or not the point is anywhere
near the boundary, and it is paid **twice over** (once for the band, once for the crossings).

The measurement matches that prediction across both corpora. Grouping parts by their sidecar bytes
per patch — a direct proxy for flattened trim length — gives a clean separation with nothing in
between:

| trim data per patch | parts | ns per candidate patch, `DistFromOutside` |
| --- | --- | ---: |
| 160–460 B (line/arc trims only) | `box`, `box_minus_cyl`, `box_union_box`, `cyl_plus_cone`, `sphere_minus_cyl`, `torus_union_cyl`, `BasePin`, `Base`, `Boom`, `Stick`, `BucketLink1`, `BucketLink2`, `Bucket` | **135 – 337** |
| 700–4700 B (B-spline / ellipse trims from curved-surface intersections) | `oblique_cut_cyl`, `cyl_inter_cyl`, `cyl_cross_cyl`, `tube_window`, the six rams | **2 253 – 5 908** |

That is a **10× to 40× step**, and it is the whole story of the surface solid's cost. The six rams
are `cylinder ∩ cylinder`, so their trims are exactly the expensive kind.

### 2.3 `Safety()` has no acceleration structure — and it shows at scale

`O2BVHSurfaceSolid::Safety()` loops over every surface with no BVH, no early-out and no bound on the
current best. Measured, on ALICE3:

| part | patches | `Safety` | `DistFromOutside` (BVH) | ratio |
| --- | ---: | ---: | ---: | ---: |
| `ST1829909_002` | 965 | **811 870 ns** | 8 658 ns | 94× |
| `ST1829909_003` | 236 | 246 376 ns | 7 668 ns | 32× |
| `ST0923290_019` | 44 | 31 885 ns | 4 473 ns | 7.1× |
| `Bagger/Bucket` | 97 | 10 767 ns | 1 139 ns | 9.5× |

811 µs / 965 patches = **841 ns per patch**, i.e. exactly the per-patch trim cost with no pruning at
all. `TGeoNavigator` calls `Safety()` on every step. This is the single largest per-call number
anywhere in this document and it is a structural property of one function, not of the geometry.

**Nothing here was optimised.** Both findings are reported, not fixed; the numbers in this document
would be invalidated by a change in the same commit that produced them.

---

## 3. The leaf-count scaling: what a BVH-over-primitives CSG solid would buy

Every genuine boolean in the corpus is a **2-leaf** union of two `TGeoTube`s
(`Stream_N_PlacedPrimitives.md`), so the corpus cannot answer this and no amount of running it
harder will make it. `--ladder 1,2,4,8,16,32` builds the missing fixture: K overlapping tubes
(pitch 0.8 against radius 0.5, so every leaf genuinely shares interior with its neighbour), unioned
two ways — a left-deep **chain**, the natural output of a fold over a list, and a **balanced** binary
tree, which is what a BVH over primitives gives you for free minus the box rejection.

| K | tree | depth | `Contains` | `Safety` | `DistFromOutside` | `DistFromInside` |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | (a bare `TGeoTube`) | 1 | **1.8** | 4.3 | 22.3 | 13.5 |
| 2 | chain | 2 | 12.6 | 20.8 | 72.4 | 57.2 |
| 2 | balanced | 2 | 12.8 | 21.1 | 73.2 | 57.2 |
| 4 | chain | 4 | 26.6 | 61.3 | 145.1 | 115.7 |
| 4 | balanced | 3 | 27.0 | 59.2 | 144.8 | 108.1 |
| 8 | chain | 8 | 49.8 | 206.7 | 256.8 | 247.9 |
| 8 | balanced | 4 | 50.4 | 162.0 | 267.5 | 228.4 |
| 16 | chain | 16 | 92.9 | **868.2** | 438.8 | **724.8** |
| 16 | balanced | 5 | 96.8 | 414.9 | 440.4 | 341.9 |
| 32 | chain | 32 | 184.3 | **3 342.9** | 812.4 | **2 374.1** |
| 32 | balanced | 6 | 182.3 | **1 058.7** | 761.5 | **641.9** |

Four readings, and the fourth is the answer to the question that prompted this stream.

1. **`TGeoCompositeShape` has no spatial pruning whatsoever.** `Contains` is 5.7 ns per leaf,
   dead linear from K=2 to K=32 (12.6 → 184.3 ns, i.e. 14.6× for 16× the leaves).
   `DistFromOutside` is 24–36 ns per leaf, also linear. Every query visits every leaf.
2. **Entering the boolean machinery costs 7× on its own.** K=1 (a bare `TGeoTube`) is 1.8 ns;
   K=2 is 12.6 ns. That is not two tubes' worth of work, it is the `TGeoBoolNode` dispatch, the
   matrix and the virtual calls. It is also why the corpus's own 2-leaf rams are 19 ns rather than
   4 ns for `Contains`.
3. **Tree shape does not matter for `Contains` or `DistFromOutside`** (within 5 %), because the
   traversal is exhaustive either way. It matters a great deal for **`Safety`** (3.2× at K=32) and
   **`DistFromInside`** (3.7× at K=32), where a left-deep chain is *super*-linear: `Safety` goes
   20.8 → 3343 ns for 16× the leaves, i.e. **161×**.
4. **Therefore: a BVH over primitives would pay, and only above about 8 leaves.** Turning Θ(K) into
   Θ(log K) on `DistFromOutside` would save nothing at K=2 (the entire corpus today), ~1.5× at K=8,
   and ~5× at K=32 — while an emitter that simply builds a *balanced* tree instead of a fold
   already recovers 3.2–3.7× of the `Safety`/`DistFromInside` cost at K=32 **at zero engineering
   cost**. And the whole ladder is dwarfed by the surface solid: a **32-leaf** chain union is
   812 ns for `DistFromOutside`, still **8× faster** than `StickCylinderOuter`'s 8-patch exact
   surface solid at 6 693 ns.

**So the honest recommendation on "should we build a high-performance general CSG solid" is: not
for speed, not yet.** Nothing in the corpus is above 2 leaves; at 2 leaves the composite is already
the fastest representation by a factor of 40–80; and the measured cost of the thing we *do* have is
80× the thing we were worried about. If a general decomposition (Tier 3) later produces 8- to
32-leaf bodies, build the balanced tree first and measure again — that is a one-line change to the
emitter and it recovers most of what a BVH would.

---

## 4. Per-part cost and memory, fixture ladder and Bagger

`--perf --raster 16 --beams 8`, 4096 points, 4096+4096 rays, 9 passes. `struct` is the exact
structural byte count (for the mesh: `nVert × 24 + nFacets × 20 + nFacets × 24`; for the surface
solid: the sidecar on disk, because the in-memory trim arrays are not introspectable from outside;
for a composite: (leaves + nodes) × ~200 B of ROOT object). `heap` is the measured `mallinfo2`
delta across load + `CloseShape`. LOST / worst |Δt| are against the OpenCascade crossing-list
oracle over 24 576 rays (96 Fibonacci beams × 16²).

### Fixture ladder — 10 parts

| part | repr | prim | Cont | Safety | dOut | dIn | ns/ray | struct B | heap B | LOST | worst \|Δt\| cm |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `box` | surface | 6 | 185 | 163 | 244 | 197 | 622 | 1 708 | 6 320 | 0 | 6.2e-15 |
| | mesh | 12 | 49 | 112 | 106 | 75 | 204 | 720 | 480 | 0 | 6.2e-15 |
| | **shape** | 1 | **2** | **3** | **16** | **17** | **64** | **200** | 1 888 | 0 | 6.2e-15 |
| `box_minus_cyl` | surface | 7 | 295 | 445 | 356 | 310 | 1 032 | 2 150 | 11 968 | 0 | 2.0e-14 |
| | mesh | 520 | 146 | 540 | 377 | 204 | 697 | 29 120 | 72 208 | 0 | 5.4e-03 |
| `box_union_box` | surface | 10 | 205 | 340 | 248 | 236 | 672 | 2 828 | 11 696 | 0 | 9.8e-15 |
| | mesh | 20 | 55 | 137 | 120 | 86 | 219 | 1 168 | 480 | 0 | 9.8e-15 |
| `cyl_cross_cyl` | surface | 8 | 2 012 | 726 | 4 770 | 4 263 | 7 153 | 7 634 | 179 056 | 0 | 7.3e-14 |
| | mesh | 1 612 | 136 | 1 297 | 429 | 306 | 593 | 90 320 | 170 496 | 0 | 1.1e-02 |
| | **shape** | 2 | **14** | **21** | **78** | **51** | **155** | **600** | 2 336 | 0 | 1.1e-13 |
| `cyl_inter_cyl` | surface | 6 | 6 960 | 367 | **13 515** | 8 057 | **22 207** | 6 942 | 166 064 | 0 (2 extra) | 2.8e-14 |
| | mesh | 604 | 135 | 1 038 | 355 | 229 | 538 | 33 872 | 79 744 | 0 | 2.1e-02 |
| `cyl_plus_cone` | surface | 4 | 168 | 351 | 274 | 200 | 462 | 814 | 11 776 | 0 | 1.2e-13 |
| | mesh | 1 760 | 161 | 1 483 | 482 | 254 | 643 | 98 608 | 175 936 | **6** | 6.4e-02 |
| `oblique_cut_cyl` | surface | 3 | 4 086 | 4 835 | 6 054 | 4 072 | 12 534 | 2 118 | 151 392 | 0 | 4.3e-14 |
| | mesh | 566 | 121 | 893 | 304 | 252 | 515 | 31 744 | 77 344 | **2** | 3.9e-02 |
| `sphere_minus_cyl` | surface | 2 | 145 | **40** | 232 | 180 | 417 | 500 | 26 832 | 0 | 4.8e-14 |
| | mesh | 7 048 | 331 | 2 599 | 880 | 567 | 1 212 | 394 688 | 721 984 | **10** | 8.9e-02 |
| `torus_union_cyl` | surface | 6 | 432 | 394 | 563 | 632 | 1 130 | 958 | 42 160 | 0 | 5.8e-13 |
| | mesh | 23 432 | 296 | 3 565 | 770 | 498 | 1 087 | 1 312 240 | 2 694 128 | **12** | 3.6e-02 |
| `tube_window` | surface | 4 | 4 587 | 352 | 9 826 | 6 496 | 17 179 | 18 820 | 190 672 | 0 | 4.2e-14 |
| | mesh | 1 248 | 189 | 935 | 500 | 280 | 791 | 69 888 | 156 864 | **12** | 6.0e-02 |

**One new defect, on the exact representation.** `cyl_inter_cyl`'s `surface` produces **2 extra
crossings, 2 unterminated transports and 2 parity mismatches** in 24 576 rays under a 96-direction
fan. `Stream_J_XRay.md` §2 reports this fixture clean; it was clean under three axis beams. This is
`NEXT.md`'s standing warning firing again — direction diversity, not ray count — and it is the only
non-zero surface-solid counter in either corpus. Reported, not diagnosed.

### Bagger — 13 parts (all 13 now carry an exact sidecar)

| part | repr | prim | Cont | Safety | dOut | dIn | ns/ray | struct B | heap B | LOST | worst \|Δt\| cm |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BasePin` | surface | 3 | 129 | 315 | 211 | 142 | 404 | 486 | 7 040 | 0 | 1.8e-13 |
| | mesh | 500 | 104 | 903 | 303 | 197 | 448 | 28 048 | 47 984 | 0 | 1.9e-02 |
| | **shape** `TGeoTube` | 1 | **2** | **4** | **26** | **6** | **68** | **200** | 1 920 | 0 | 1.8e-13 |
| `Base` | surface | 44 | 476 | 4 708 | 642 | 626 | 1 559 | 14 720 | 134 672 | 0 | 3.9e-13 |
| | mesh | 4 364 | 189 | 967 | 457 | 410 | 928 | 244 096 | 589 824 | 1 | 7.6e-02 |
| `BoomCylinderInner` | surface | 6 | 1 019 | 587 | 6 101 | 4 997 | 4 630 | 13 110 | 104 288 | 0 | 2.7e-12 |
| | mesh | 1 612 | 401 | 3 420 | 819 | 774 | 796 | 90 272 | 163 616 | 0 | 5.1e-03 |
| | **shape** | 2 | **18** | **21** | **80** | **53** | **107** | **600** | 2 832 | 0 | 2.7e-12 |
| `BoomCylinderOuter` | surface | 8 | 2 633 | 653 | 6 625 | 8 468 | 10 105 | 23 756 | 348 928 | 0 | 6.4e-12 |
| | mesh | 2 244 | 698 | 5 774 | 1 666 | 1 157 | 1 624 | 125 664 | 288 320 | 2 | 8.5e-03 |
| | **shape** | 2 | **19** | **23** | **82** | **62** | **134** | **600** | 2 656 | 0 | 6.4e-12 |
| `Boom` | surface | 31 | 768 | 3 555 | 838 | 1 054 | 1 812 | 10 678 | 98 752 | 0 | 7.4e-12 |
| | mesh | 4 520 | 204 | 1 650 | 390 | 472 | 754 | 252 832 | 590 848 | 2 | 1.4e-02 |
| `BucketCylinderInner` | surface | 6 | 1 133 | 571 | 4 462 | 4 508 | 4 159 | 11 702 | 101 024 | 0 | 5.8e-11 |
| | mesh | 1 636 | 384 | 3 520 | 873 | 758 | 669 | 91 616 | 162 720 | 0 | 9.2e-03 |
| | **shape** | 2 | **20** | **24** | **87** | **54** | **106** | **600** | 2 912 | 0 | 5.8e-11 |
| `BucketCylinderOuter` | surface | 10 | 2 328 | 895 | 5 142 | 6 697 | 9 909 | 22 144 | 347 552 | 0 | 5.6e-11 |
| | mesh | 2 178 | 708 | 5 329 | 1 647 | 1 198 | 1 634 | 121 968 | 286 512 | 0 | 8.6e-03 |
| | **shape** | 2 | **22** | **26** | **82** | **59** | **143** | **600** | 3 072 | 0 | 5.6e-11 |
| `BucketLink1` | surface | 22 | 426 | 2 570 | 647 | 635 | 1 401 | 6 800 | 85 600 | 0 | 4.6e-10 |
| | mesh | 4 316 | 335 | 1 409 | 714 | 665 | 1 211 | 241 600 | 593 408 | 0 | 2.9e-02 |
| `BucketLink2` | surface | 23 | 446 | 3 545 | 741 | 775 | 1 446 | 7 746 | 98 544 | 0 | 2.9e-10 |
| | mesh **(open!)** | 4 628 | 111 | 1 163 | 253 | 346 | 421 | 263 632 | 609 312 | **6 591** | 7.2e-02 |
| `Bucket` | surface | 97 | 639 | 10 767 | 1 139 | 785 | 2 135 | 44 302 | 481 696 | 0 | 1.3e-10 |
| | mesh | 9 556 | 144 | 764 | 315 | 377 | 787 | 534 992 | 1 215 888 | **0** | 1.4e-02 |
| `StickCylinderInner` | surface | 6 | 1 031 | 603 | 6 092 | 4 965 | 4 634 | 13 110 | 104 128 | 0 | 8.4e-11 |
| | mesh | 1 612 | 406 | 3 410 | 817 | 766 | 796 | 90 272 | 162 304 | 4 | 2.1e-02 |
| | **shape** | 2 | **18** | **21** | **77** | **53** | **108** | **600** | 2 912 | 0 | 8.4e-11 |
| `StickCylinderOuter` | surface | 8 | 2 676 | 678 | 6 693 | 8 444 | 9 834 | 23 756 | 350 224 | 0 | 2.2e-11 |
| | mesh | 2 244 | 704 | 5 902 | 1 685 | 1 269 | 1 649 | 125 664 | 288 880 | 0 | 2.4e-02 |
| | **shape** | 2 | **19** | **23** | **82** | **62** | **137** | **600** | 3 392 | 0 | 2.2e-11 |
| `Stick` | surface | 24 | 701 | 2 761 | 623 | 886 | 1 723 | 8 312 | 89 072 | 0 | 1.7e-10 |
| | mesh | 4 412 | 149 | 747 | 336 | 370 | 548 | 246 784 | 583 424 | 0 | 2.0e-02 |

**`Bagger/Bucket` — the open question in `MeshHealing.md` is answered: its mesh is fine.**
`meshClosedBody = true`, **0 LOST crossings** in 24 576 rays at both `--mesh-prec 0.1` and `0.05`,
0 unterminated, 0 parity, worst deviation 1.4e-02 cm on a 19.4 cm diagonal (7.0e-04 relative). The
one part that used to ship as a mesh does **not** have the shared-edge problem. (Note that a
concurrent change has since given `Bucket` an exact sidecar, so Bagger ships no mesh part at all;
the measurement stands as a corpus measurement.)

### The two memory numbers, and where they disagree

They disagree by a factor of **1.2 to 3.5**, always in the same direction, and the reason is
allocator granularity rather than anything geometric:

| representation | structural | measured heap | ratio | why |
| --- | ---: | ---: | ---: | --- |
| `Bagger/Bucket` mesh | 535 kB | 1 216 kB | 2.3× | `std::vector` growth doubling, plus the BVH built in `CloseShape` (which the structural formula does not count) |
| `Bagger/Base` mesh | 244 kB | 590 kB | 2.4× | same |
| `Bagger/StickCylinderOuter` surface | 24 kB (sidecar) | 350 kB | 15× | the sidecar is a *compressed* description; the flattened B-spline polylines are built at load and are not in it |
| any `shape` | ~0.2–0.6 kB | 1.9–3.4 kB | ~6× | ROOT object overhead, `TNamed`, the `TGeoBoolNode` |

**The structural number is the exact one and it is the one to extrapolate with**, except for the
surface solid, where the sidecar is a lower bound only and the measured heap is the number that
matters. `residentBytes` (`/proc/self/statm`) is the noisiest of the three — it moves in 64 kB steps
and reads 0 for most parts because the pages were already mapped. It is reported and should be
ignored below ~1 MB.

---

## 5. The mesh-precision sweep — and the knee is not where the guess was

`--mesh-prec X` sets **both** `lin_defl` and `ang_defl` of `BRepMesh_IncrementalMesh` to X. They are
not the same knob: `lin_defl` is a **length in the model's units** bounding the chordal sagitta;
`ang_defl` is an **angle in radians** bounding the turn between chords. `--mesh-prec 0.5` therefore
asks for a half-unit chord bound *and* **0.5 rad = 28.6°**, which makes a cylinder a 12-sided prism
however fine the linear term is.

So the sweep was run twice. **Sweep A** reproduces the coupled flag. **Sweep B** decouples them,
holding `ang_defl = 0.1 rad` and varying `lin_defl` alone, via `remeshFromBrep.py` — a new tool that
re-tessellates the `brep_<part>.brep` the converter already dumps and writes `facets_<part>.bin` in
the converter's exact format. Its equivalence to the converter was checked first: at the converter's
own settings it reproduces **12 of 13** Bagger parts to the triangle, the 13th (`BucketLink2`, the
one with the broken mesh) differing by 8 triangles in 4 628 — the converter meshes the solid in
STEP units, this tool meshes the cm-scaled `.brep`, and OCCT's internal absolute tolerances make
that a marginal difference exactly on the degenerate part. **The `.brep` is written in cm, so
`--lin` is a cm and the converter's `--mesh-prec` is not.**

Bagger, 13 parts, 319 488 rays against the OCCT crossing-list oracle. `cost` columns are the median
over the 13 parts.

| sweep | lin | ang | triangles | struct MB | `Cont` | `Safety` | `dOut` | `dIn` | ns/ray | identical/rays | LOST | worst \|Δt\| cm | unterm | parity |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **A** `--mesh-prec 1.0` | 0.1 cm | 1.0 rad | 4 704 | 0.26 | 117 | 540 | 266 | 287 | 395 | 235 075/319 488 | 9 137 | **5.16** | 676 | 1 365 |
| **A** `--mesh-prec 0.5` | 0.05 cm | 0.5 rad | 8 710 | 0.49 | 131 | 642 | 316 | 327 | 537 | 235 242/319 488 | 7 196 | **4.39** | 682 | 1 375 |
| **A** `--mesh-prec 0.1` | 0.01 cm | 0.1 rad | 43 830 | 2.46 | 292 | 1 586 | 695 | 643 | 815 | 235 331/319 488 | 6 600 | **7.56e-02** | 673 | 1 359 |
| **A** `--mesh-prec 0.05` | 0.005 cm | 0.05 rad | 96 579 | 5.42 | 414 | 2 439 | 1 003 | 881 | 1 081 | 234 204/319 488 | **10 697** | **2.07e-02** | **1 843** | **4 537** |
| **B** decoupled | 0.1 cm | 0.1 rad | 42 020 | 2.36 | 300 | 1 592 | 665 | 608 | 738 | 235 331/319 488 | 6 600 | 7.56e-02 | 673 | 1 359 |
| **B** decoupled | 0.05 cm | 0.1 rad | 43 734 | 2.45 | 293 | 1 592 | 681 | 649 | 809 | 235 332/319 488 | 6 602 | 7.56e-02 | 673 | 1 359 |
| **B** decoupled | 0.01 cm | 0.1 rad | 43 830 | 2.46 | 312 | 1 602 | 736 | 684 | 809 | 235 331/319 488 | 6 600 | 7.56e-02 | 673 | 1 359 |
| **B** decoupled | 0.005 cm | 0.1 rad | 43 824 | 2.46 | 290 | 1 620 | 674 | 673 | 809 | 235 332/319 488 | 6 602 | 7.56e-02 | 673 | 1 359 |
| — `surface`, same rays | | | 288 patches | 0.20 | 768 | 895 | 1 139 | 1 054 | 2 135 | **319 488/319 488** | **0** | **4.6e-10** | **0** | **0** |
| — `shape` (7 parts) | | | 13 leaves | 0.003 | 19 | 23 | 82 | 54 | 108 | **172 032/172 032** | **0** | **8.4e-11** | **0** | **0** |

**Four results, and three of them are not what was expected.**

1. **The knee is between `--mesh-prec 0.5` and `0.1`, and it is a cliff, not a knee.** The worst
   deviation goes 4.39 cm → 7.56e-02 cm — a factor of **58** — for 5× the triangles. The user's
   guessed 0.5 sits on the wrong side of it. Below 0.1 the returns are ordinary: 0.05 buys another
   3.7× in deviation for 2.2× the triangles and 1.3× the per-call cost.
2. **`--mesh-prec` is an angular knob in disguise.** Sweep B changes `lin_defl` by **20×** (0.1 →
   0.005 cm) at fixed `ang_defl`, and the triangle count moves 4 %, the worst deviation not at all
   (7.558e-02 cm in every row), LOST by 2 in 6 600. On Bagger — cylinders and planes throughout —
   **the linear deflection does nothing**; the accuracy is set entirely by `ang_defl`. Row A
   `--mesh-prec 0.1` and row B `lin 0.01 / ang 0.1` are bit-identical, which is the positive control
   for that claim. Anyone tuning `--mesh-prec` on quadric CAD is tuning the angle and paying for a
   chord bound they are not using.
3. **Mesh validity does not degrade monotonically with precision, and the FINEST mesh is the
   WORST.** Refining from `--mesh-prec 0.1` to `0.05` takes LOST from 6 600 to **10 697**,
   unterminated from 673 to **1 843**, parity from 1 359 to **4 537** — while the chordal deviation
   improves 3.7×. Every one of those events is on **`BucketLink2`**, the one part with
   `meshClosedBody = false`: a finer tessellation of two independently-meshed curved faces puts
   *more* triangles along the seam they do not share, so there are more gaps, not fewer. **A
   chordal-deviation number is not a safety number, and the sweep that improves one degrades the
   other.** `MeshHealing.md` says a mesh can be invalid rather than inaccurate; this is that
   sentence with a slope on it.
4. **The mesh is 2 to 5× cheaper per call than the exact surface solid on Bagger, at every
   precision tested**, and at `--mesh-prec 1.0` it is 5.4× cheaper in transport — with a 5 cm worst
   deviation, which is useless. The usable trade is `--mesh-prec 0.1`: 2.6× cheaper transport than
   `surface`, 12× more memory, 7.6e-04 relative accuracy, and 6 600 LOST crossings that are entirely
   one broken part's.

### Per-part chordal error, normalised

At `--mesh-prec 0.1`, |Δt| divided by the part's bounding-box diagonal:

| part | tri | worst \|Δt\| cm | diag cm | relative | LOST |
| --- | ---: | ---: | ---: | ---: | ---: |
| `BoomCylinderInner` | 1 612 | 5.13e-03 | 20.6 | **2.5e-04** | 0 |
| `Boom` | 4 520 | 1.36e-02 | 48.8 | 2.8e-04 | 2 |
| `BoomCylinderOuter` | 2 244 | 8.53e-03 | 21.3 | 4.0e-04 | 2 |
| `BucketCylinderOuter` | 2 178 | 8.56e-03 | 18.8 | 4.6e-04 | 0 |
| `BucketCylinderInner` | 1 636 | 9.18e-03 | 18.9 | 4.9e-04 | 0 |
| `Stick` | 4 412 | 1.99e-02 | 40.5 | 4.9e-04 | 0 |
| `Bucket` | 9 556 | 1.35e-02 | 19.4 | 7.0e-04 | 0 |
| `StickCylinderInner` | 1 612 | 2.10e-02 | 20.6 | 1.0e-03 | 4 |
| `StickCylinderOuter` | 2 244 | 2.37e-02 | 21.3 | 1.1e-03 | 0 |
| `BasePin` | 500 | 1.92e-02 | 10.4 | 1.8e-03 | 0 |
| `BucketLink1` | 4 316 | 2.88e-02 | 7.5 | 3.8e-03 | 0 |
| `Base` | 4 364 | 7.56e-02 | 18.2 | 4.1e-03 | 1 |
| **`BucketLink2`** | 4 636 | 7.18e-02 | 9.9 | **7.3e-03** | **6 591** |

The relative deviation is **2.5e-04 to 7.3e-03** — i.e. the tessellation error is a few parts in a
thousand of the part's own size, and it is bounded by `ang_defl` × the local radius, exactly as the
theory says. That is the answer to "quantify the numerical error due to tessellation": **at
`--mesh-prec 0.1`, a few 1e-4 to a few 1e-3 of part size, and 1e-2 to 8e-2 cm in absolute terms.**
Against `surface`'s 4.6e-10 cm and `shape`'s 8.4e-11 cm, that is **eight orders of magnitude.**

---

## 6. ALICE3 — and the tessellation wall is a `--mesh-prec` wall, not an ALICE3 wall

`scripts/geometry/ALICE_3_example/CAD_noETA.stp`, 55 leaf solids / 9 266 faces.

### 6.1 It meshes. The documented 22.9 GB is reproducible and it is two precision steps away

`Tutorial.md` §5.6 and `NEXT.md` record: *"converting a single 2 m sphere reached 22.9 GB resident
before being killed"*, and instruct converting ALICE3 **without `--mesh`**. That is true at the
converter's **default** `--mesh-prec 0.1` and **false two steps coarser.** Measured, whole model,
under an 8 GB `ulimit -v`:

| `--mesh-prec` | outcome | wall | peak RSS | triangles | parts meshed |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1.0 | **completed** | 14.6 s | **679 MB** | 138 923 | 55 / 55 |
| 0.5 | **completed** | 16.7 s | **841 MB** | 345 517 | 55 / 55 |
| 0.1 (the default) | **killed at the 8 GB cap** | 2 m 08 s | 7.34 GB | — | 0 / 55 |

So **ALICE3 can be tessellated** — at 0.5 the whole model is 345 k triangles and 841 MB, which is
nothing. The wall is real, it is the default, and it is a *linear-deflection* wall: at 0.1 mm on a
2 m feature the triangle count is a ratio of 20 000. This invalidates the flat instruction in
`NEXT.md`'s "Traps in the environment" and in `Tutorial.md` §5.6; the correct statement is *convert
ALICE3 with `--mesh-prec 0.5` or coarser, or without `--mesh`*.

The same is true from the `.brep` side: `remeshFromBrep.py` at `lin 1.0 cm / ang 0.5 rad` meshed all
**20** exact-sidecar solids in **0.7 s** total for 48 072 triangles, with a 4 GB cap and a 180 s
timeout per part, and nothing came close to either.

### 6.2 Three-way cost on ALICE3, 20 parts

Conversion without `--mesh`: **25.1 s, 918 MB, 20 exact sidecars, 2 CSG shapes** (both
`TGeoTubeSeg`) of 55 leaf solids. Meshes from `--mesh-prec 0.5`. 1024 points, 1024+1024 rays, 5
passes, raster 6 × 3 axis beams. Two of the 20 sidecars still do not load (`ST1829909_004`,
`ST1829909_01` — the known 1e-06 cm wire-join tolerance defect), so 18 surface rows.

Medians over parts, and totals for the memory and load columns:

| repr | parts | primitives | `Cont` | `Safety` | `dOut` | `dIn` | **ns/ray** | heap total | load+close total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `surface` | 18 | 1 936 patches | 2 970 | 11 994 | 3 431 | 2 575 | **7 865** | **22.3 MB** | **14.63 s** |
| `mesh` (prec 0.5) | 20 | 112 930 tri | **115** | **669** | **398** | **290** | **570** | **12.8 MB** | **2.78 s** |
| `shape` | 2 | 2 leaves | **6** | **7** | **31** | **22** | **69** | 0.003 MB | 0.00 s |

**On ALICE3 the coarse mesh beats the exact surface solid on every axis simultaneously: 14× faster
in transport, 26× faster on `Contains`, 18× on `Safety`, 1.7× smaller in memory and 5× faster to
load.** That is the case the user anticipated — *"we might stay with tessellated in some cases"* —
and ALICE3 is that case, on cost. Whether it is that case on *correctness* is §6.3.

The worst individual rows are worth naming, because they are where a production run would hurt:

| part | patches | `Safety` | `dOut` | ns/ray | `CloseShape` | heap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ST1829909_002` | 965 | **811 870 ns** | 8 658 | 20 914 | **9.75 s** | 13.2 MB |
| `ST1829909_003` | 236 | 246 376 ns | 11 174 | 31 305 | 1.11 s | 3.2 MB |
| `ST2487462_01` | 80 | 29 799 ns | 12 738 | 32 848 | 2.50 s | 1.6 MB |
| the same three, as `mesh` at prec 0.5 | 19 304 / 12 598 / 2 196 tri | 1 062 / 1 386 / 511 | 380 / 744 / 322 | 570 / 801 / 368 | 0.36 / 0.15 / 0.01 s | 2.5 / 1.3 / 0.3 MB |

`CloseShape` is **9.75 s for one part**, and it is the whole of `Stream_J_XRay.md` §7's "14 s of
load dominates every density". It is closure, not loading, and it is per part per process.

### 6.3 And the correctness verdict goes the other way

Full oracle round trip on the same 20 parts, three axis beams, 768 rays per part, both stepping
modes, against OpenCascade's crossing lists — i.e. exactly the comparison §6.2's cost table cannot
make:

| repr | rays identical | LOST | extra | displaced | worst \|Δt\| cm | parts clean | unterm | parity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `surface` (a) / (b) | **13 824 / 13 824** | **0** | **0** | **0** | 1.14e-07 / 4.85e-07 | **18 / 18** | 0 | 0 |
| `shape` (a) / (b) | **1 536 / 1 536** | **0** | **0** | **0** | 7.1e-15 / 2.1e-15 | **2 / 2** | 0 | 0 |
| `mesh` @ `--mesh-prec 0.5` | 8 575 / 15 360 | **160** | **109** | 21 682 | **1.300** | **0 / 20** | 1 | 1 |

**A coarse ALICE3 mesh loses 160 walls and displaces crossings by up to 1.3 cm.** Worst offenders:
`ST2487721_01` and `ST1A38494_003` at 40 LOST each, `ST1782525_01` at 24, and — note — the two
extremes of triangle count are both bad (`ST1A38494_003` has **68** triangles and loses 40
crossings; `ST1829909_003` has **12 598** and loses 14). Zero of twenty parts are clean.

So the ALICE3 story is: **the mesh is 14× cheaper and it is not correct enough at the precision that
makes it cheap.** At `--mesh-prec 0.1` — the precision that would fix the accuracy — the model does
not mesh at all (§6.1). The cost measurement therefore does not become a recommendation; what it
becomes is a strong argument for spending effort on `curveTrimContains` and on `Safety()`, because
the exact path is the only one that is *right* on ALICE3 and it is the one that is slow.

The 18 loaded surface parts are **18/18 clean with every robustness counter zero over 46 318 steps**,
which reproduces `Stream_L_ALICE3Defect.md` / `Stream_M_Quartic.md`'s post-fix state on a
differently-built raster.

---

## 7. Which parts are legitimately better off tessellated

Reading §4, §5 and §6 together, and taking correctness first:

| case | representation | why |
| --- | --- | --- |
| anything a CSG recogniser accepts (7 of 13 Bagger parts, 2 of 55 ALICE3 solids) | **`shape`** | exact (`dV_sym = 0`), 40–140× faster than every alternative, ~600 B, no closure cost. There is no argument against it |
| quadric parts with **line/arc trims only** — `Base`, `Boom`, `Stick`, `BucketLink1`, `Bucket`, `box*`, `sphere_minus_cyl`, `torus_union_cyl` | **`surface`** | 135–337 ns per candidate patch, 0 LOST, 1e-10 cm; 1.5–2.5× the mesh's per-call cost but 10–50× smaller and exact. The mesh saves little here and gives up eight orders of accuracy |
| quadric parts whose trims are **curved-surface intersections** — the six rams, `cyl_cross_cyl`, `cyl_inter_cyl`, `oblique_cut_cyl`, `tube_window` | **`shape` where recognised; otherwise a genuine question** | the surface solid is 2.2–5.9 µs per candidate patch here. Where CSG recognises the body (the six rams) it wins outright. Where it does not (`cyl_inter_cyl` at 22 µs/ray, `tube_window` at 17 µs/ray) the mesh is **41×** and **22×** faster with a 2.1e-02 / 6.0e-02 cm error and 0 / 12 LOST. **That is the trade worth having a policy about**, and today nothing computes it |
| **ALICE3-scale free-form-trimmed solids** — the 18 loaded ALICE3 parts, especially `ST1829909_002/003`, `ST2487462_01` | **`surface`, despite the cost** | the mesh is 14× faster, 1.7× smaller and 5× faster to load — **and loses 160 walls with a 1.3 cm worst displacement, 0/20 parts clean** (§6.3). The precision that would fix that (`--mesh-prec 0.1`) does not fit in 8 GB. So the exact path stays, and its 9.75 s `CloseShape` and 812 µs `Safety()` become the thing to fix |
| **any part whose mesh reports `meshClosedBody = false`** | **never the mesh** | `BucketLink2`: 6 591 LOST crossings, and refining the mesh makes it **worse**. There is still no policy attached to that flag |

**The one row that is a genuine open question** is the third: quadric parts with curved-intersection
trims that CSG does *not* recognise. `cyl_inter_cyl` costs 22 µs/ray exact and 0.54 µs/ray meshed —
**41×** — for 2.1e-02 cm of error and **zero** LOST crossings at 96 directions. That is a real trade
and today nothing computes it. Extending the CSG recogniser to cover those bodies would remove the
question entirely, which is a better answer than choosing a side.

---

## 8. What a reviewer should challenge

In the order I would attack this document.

1. **Warm cache is doing work in the mesh's favour and I have not bounded it.** Every solid here is
   resident and most fit in L2. The mesh's per-call advantage on ALICE3 (14× in transport) is
   measured on a working set of 2.5 MB against the surface solid's 0.2 MB. In a real simulation the
   mesh misses more. I did not measure a cold or interleaved case, and §1 says so. It happens not to
   change the ALICE3 verdict — §6.3 settles that on correctness, not on cost — but every mesh-vs-
   surface ratio in §4 and §5 is exposed to it.
2. **Bagger's meshes are compared at four precisions; its surface and CSG representations at one.**
   There is no precision knob on those, but there *is* a tolerance knob (`kBSplineFlatness`, the
   flattening the whole cost of §2.2 hangs on) and I did not sweep it. A coarser flattening would
   make the surface solid cheaper and less accurate — i.e. the surface solid has a cost/accuracy
   curve too, and this document only measures the mesh's.
3. **The trim-complexity explanation is a correlation plus a profile, not a controlled experiment.**
   The clean 10–40× step between the two groups in §2.2, and the 76 % of cycles in
   `curveTrimContains`, are strong. But I did not build the decisive fixture: the *same* solid with
   the *same* patch count and two different trim representations. That fixture would take an hour
   and would settle it.
4. **`--ladder` uses overlapping tubes on a line.** That is one spatial arrangement. A ladder of
   *disjoint* leaves would let a short-circuiting evaluator stop early and might scale differently;
   a ladder of *coaxial* leaves (which is what the corpus actually has) might too. The Θ(K) result
   is robust — `Contains` is dead linear over four octaves — but the constant is arrangement-
   dependent and the "a BVH would pay above 8 leaves" threshold rests on the constant.
5. **The 2 extra crossings on `cyl_inter_cyl`'s surface solid** (§4) are new, unexplained, and on
   the representation this project trusts. They appear at 96 directions and not at 3. Someone should
   decide whether that is a real defect or an artefact of my raster.
6. **`ns/crossing` is a derived number and it flatters parts that find fewer crossings.**
   `BucketLink2`'s mesh looks fast partly because it loses 358 crossings out of 1 728. I report
   crossings beside it, but the column invites the mistake.
7. **The structural memory column for `surface` is the sidecar on disk**, which is a compressed
   description and 15× smaller than the measured heap. It is honest (the formula says so) but it is
   not comparable with the mesh's structural column, which is exact. Any memory ratio that mixes the
   two is wrong.
8. **One machine, one run of each configuration, other tenants present.** The spreads are reported
   (mostly 1–10 %, occasionally 50 % on the sub-10 ns rows where the clock resolution bites) and the
   median is over 9 passes, but nothing here is reproduced on a second machine.

---

## 9. What this invalidates elsewhere

Recorded here rather than edited into the central documents, per the working agreement.

| document | claim | measurement |
| --- | --- | --- |
| `NEXT.md`, "Traps": *"Convert ALICE3 **without `--mesh`** (24 s / 581 MB). With meshing, one 2 m sphere reached **22.9 GB** resident and was killed."* | true at the default `--mesh-prec 0.1` (reproduced: killed at an 8 GB cap, 7.34 GB, 2 m 08 s) and **false at 0.5 or 1.0**, where the whole 55-solid model meshes in 15–17 s and 679–841 MB. The trap is the default, not the model | §6.1 |
| `Tutorial.md` §5.6: *"Tessellation, not the exact path, is the scaling problem."* | **half inverted.** On *cost* it is wrong at ALICE3 scale: the exact path spends 9.75 s in `CloseShape` on one solid and 811 µs per `Safety()` call on it, while the coarse mesh of the same 20 parts loads 5× faster and queries 14× faster. On *correctness* it stands: that mesh loses 160 walls | §6.2, §6.3 |
| `Tutorial.md` §5.2: *"A dedicated `O2CSGSolid` (flat DNF) is only needed if native composite trees turn out too slow."* | native composite trees are **not** too slow — they are the fastest representation measured, by 40–140×. They are Θ(leaves) with no pruning, so the case reopens above ~8 leaves, which nothing in the corpus reaches | §2, §3 |
| `MeshHealing.md`: *"Nobody has specifically checked whether [`Bagger/Bucket`'s] curved faces produce the same shared-edge problem."* | checked: **they do not.** `meshClosedBody = true`, 0 LOST in 24 576 rays at two precisions | §4 |
| `Stream_J_XRay.md` §7: *"loading dominates … 14 s of `LoadSurfaceSolid` + `CloseShape`"* | localised: it is **`CloseShape`**, not loading. On `ST1829909_002`, 0.012 s to read the sidecar and 9.75 s to close it | §6.2 |
| `Stream_J_XRay.md` §2: fixtures, `surface`, 0 extra crossings | at **96** Fibonacci beams `cyl_inter_cyl` shows 2 extra crossings, 2 unterminated and 2 parity. §2 was measured at 3 axis beams | §4 |

---

## 10. Reproducing it

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH

# the instrument's own controls -- no database, no oracle, 37 checks
$B/stage/bin/o2-bench-detectorsbase-xray --self-test

# the synthetic boolean ladder -- no database either
$B/stage/bin/o2-bench-detectorsbase-xray --ladder 1,2,4,8,16,32 --json ladder.json

# a corpus, cost and memory, all representations from one sample set
cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/P_bag \
    --model scripts/geometry/STEP_examples/Bagger.step --beams 96 --raster 16
$B/stage/bin/o2-bench-detectorsbase-xray --db /tmp/P_bag/db --perf --raster 16 --beams 8 \
    --json /tmp/P_bag/perf.json

# the mesh sweep: re-tessellate the .brep at a stated lin AND ang, then re-score the SAME rays
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH   # PREPEND
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$PYTHONPATH
$SW/Python/latest/bin/python3.10 scripts/geometry/remeshFromBrep.py \
    --brep-dir /tmp/P_bag/db/Bagger --lin 0.01 --ang 0.1 --out-dir /tmp/P_sweep
# then point a copy of manifest.json's `facets` at /tmp/P_sweep and re-run with
#   --ref-crossings /tmp/P_bag/xray --representations mesh
```

`--perf` writes everything above as JSON: per representation, `contains` / `safety` /
`distFromOutside` / `distFromInside` / `transport` each with `nsPerCallMedian`, `nsPerCallMin`,
`nsPerCallMax`, `spread`, `checksum` and (for ray kernels) `hitFraction`; plus `structuralBytes`
with its `structuralFormula`, `sidecarBytes`, `heapBytesLoad`, `heapBytesClose`,
`residentBytesTotal`, `meshClosedBody`, and for the surface solid a `localise` block with BVH
candidates per call, pruned/unpruned/`_Loop` timings and ns per candidate patch.
