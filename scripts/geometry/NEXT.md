# NEXT — session-start instruction for the CAD → TGeo work

This file is the current hand-over. Whoever finishes a session should **rewrite it**.
Last updated 2026-08-09, after the sub-patch BVH session; the bulk of the file below is still the
2026-08-03 hand-over and its numbers predate Stream X.

Branch `swenzel/bvhsurfacesolid`. Everything below is committed.

**(2026-08-09) The sub-patch BVH, the safety anchor seed, and approximate safety** — one BVH leaf
per *cover box* instead of per surface (a full cylinder is 8 phi-chunk boxes, a sphere a 4×8 grid
of the whole ball — the sphere/torus `distanceSqToPatch` realises on the whole surface, so their
covers must too), per-query surface dedup keeps every answer bit-identical to the `_Loop` twins.
ALICE3 `ST1829909_002`: Safety candidates 8.63 → **3.94**/call, distout 4.2 → **1.7**, transport
−8 %; both gates re-run clean (fixtures 10/10, Bagger 13/13, exit 0). `Safety()` is seeded from
24 on-patch anchor points (VecGeom master's tessellated trick) and gained an opt-in approximate
mode `SetSafetySlack(s)` — a guaranteed underestimate within (1−s) of exact; a far point answers
from the root box in **127 ns instead of 39 µs, zero patches evaluated**. Unit suite now
**112 cases**. `Stream_X_SubPatchBVH.md` has the measurements and the soundness argument.
(Also noted there: `surfaces_ST1829909_01…bin`, the largest ALICE3 sidecar, fails validation at
load with a 5.4e-6 cm wire-join gap — a converter-side item, visible only now that someone tried
to load it.)

## Read this first, and possibly only this

**[`Tutorial.md`](Tutorial.md)** is the orientation document: the three representations and the
cascade between them, the pipeline, the artefacts, the four commands, the seven principles, and a
per-capability status with its explicit gaps. It was written because five investigation records
and two reviews described this branch and none of them answered *"what exists, how do I drive it,
what is missing"*. Start there. This file is only the hand-over on top of it.

## Where the branch stands

| | |
| --- | --- |
| `ctest -R BVHSurfaceSolid` | **108 cases**, green |
| ladder fixtures | **10/10 scored** of 10 leaf solids — `oblique_cut_cyl` emits since the ellipse trim |
| `Bagger.step` | **13/13 scored, 13/13 pass, gate exits 0** |
| unexplained oracle disagreements | **0/0/0/0** on `surface`, both corpora; **0/0/0/0** on `shape` |
| `runOracleGate.py --self-test` | 17/17 |
| `O2_CADtoTGeo.py --self-test` | **36 checks**, 0 failures (18 recognition + 8 placement + 10 planar trim) |
| `csg/emit.py --self-test` | 33/33 |
| `o2-bench-detectorsbase-xray --self-test` | **34 checks**, 0 failures |

**Gate totals and disagreement counts are separate numbers. Never quote one without the other.**

**Bagger is now represented entirely exactly: CSG 7, exact surfaces 6, tessellated 0**, of 13 leaf
solids. Every CSG part is exact (`dV_sym = 0`) and every part is oracle-clean. The gate exit code is
`0` for the first time. The MVP was a *mixture* including one mesh; this is past it.

`Bucket` was the last mesh part, and it was never about its spheres and tori — it has 69 planes,
22 cylinders, 4 spheres, 2 tori and exactly **two** unsupported faces, both planes bounded by an
**ellipse**. A conic *is* a rational quadratic B-spline, exactly, and the sidecar's planar B-spline
record is rational, so widening one set literal was the whole fix. See `Stream_Q_EllipseTrim.md`.

## What landed today, in order

- **Capacity** — the residual that had survived five sessions was Gauss-Legendre spanning a
  B-spline's knot spans, hidden behind an interval count taken from a *closed* curve's endpoints
  (so `max(1, ceil(0/x))` was 1 at every x). 4.9e-04 → 2.8e-07. `Stream_C_Hygiene.md`.
- **Edge identity (sidecar v3)** — closure stopped deciding a topological question with a
  proximity query. Parts failing on navigability went **7/21 → 0/21**. `Stream_F_EdgeIdentity.md`.
- **Position/scale sweep** — position-independent yes, scale-independent **no**; found the quartic
  guard defect below. `Stream_E_Scale.md`.
- **CSG census** — 1004 of ALICE3's 2377 "B-spline" surfaces are quadrics in disguise.
  `Stream_A_CSG.md`.
- **Any-`TGeoShape` gating**, then the **CSG recogniser/emitter**, then the
  **representation-aware verdict**. `Stream_G_AnyShape.md`, `Stream_H_CSGEmitter.md`,
  `Stream_I_Verdict.md`.
- **X-ray transport benchmark** — the only instrument here that *steps* rather than asking
  single-shot questions, and the only one that runs on ALICE3. `Stream_J_XRay.md`.
- **(2026-08-02) Placed primitives** — the self-union composite is gone: a recognised primitive is
  emitted at the origin with its rigid transform stored beside it, and `Bagger/BasePin` is a
  `TGeoTube` again with an analytic `Capacity()` (2.045e-04 Monte-Carlo → 6.785e-16).
  `Stream_N_PlacedPrimitives.md`.
- **(2026-08-02) The ellipse trim** — Bagger 12→13 and the fixtures 9→10, at machine precision
  (7.9e-15 cm), by widening one set literal. **Bagger now ships no tessellated part.**
  `Stream_Q_EllipseTrim.md`.
- **(2026-08-02) The representation benchmark** — per-call ns, memory and accuracy for all three
  representations. **`TGeoCompositeShape` is the *fastest*, by 40–145×**, not the slow one.
  `Stream_P_RepresentationBench.md`.
- **(2026-08-02) `Safety`/`ComputeNormal` use the BVH** — 835 µs → 22.4 µs on ALICE3's 965-patch
  solid, bit-for-bit identical to the loop over 74329 points. `Stream_S_SafetyBVH.md`.
- **(2026-08-03) The assembly oracle and the overlap census** — and the finding that outranks
  everything else here: **the models are not legal for Geant4.** `Stream_T_AssemblyOracle.md`.
- **(2026-08-03) Two exact-coverage routes measured and closed.** Co-surface half-space trims carry
  **1 of 15** target solids, not 15 (`Stream_R_CoSurfaceTrims.md`); merging the exporter's
  co-surface patches is worth **1 of 16** (`Stream_U_CoSurfaceMerge.md`).

## Open, in the order I would take them

0a. **(2026-08-09) `surfaces_ST1829909_01…bin` — ALICE3's largest exact sidecar (1.46 MB) — fails
   `LoadSurfaceSolid` validation**: `surface 1006: wire edge 1 end does not join the next edge
   start (gap 5.41e-06 cm, tolerance 1e-06 cm)`. Nobody had loaded this sidecar standalone
   before, so the failure predates Stream X and is not caused by it. Two candidate causes, in
   the order to check: (1) the loader judges join gaps against the fixed fallback
   `kWireJoinTolerance = 1e-6 cm` (`O2SurfaceSolidIO.cxx:279`) even when the sidecar (v2+)
   declares the model's own tolerance — the constant's own doc block calls it "a fallback, not a
   measurement of the model", and `GetRimMatchTolerance` already prefers the declared tolerance
   for the same kind of decision; (2) the converter emitted genuinely non-joining edges for this
   face and should share endpoints exactly. Decide on evidence (read the sidecar's declared
   tolerance, measure the gap's effect on containment), not by widening a constant.

0. **THE MODELS ARE NOT LEGAL FOR GEANT4, AND THIS BLOCKS THE STATED NEXT GOAL.** Neither
   `Bagger.step` nor ALICE3 composes into a world TGeo or Geant4 will accept: both contain placed
   solids with **positive shared volume**. Bagger has **3 of 78** pairs — `Base`∩`BoomCylinderOuter`
   8.66e-02 cm³, `Stick`∩`Bucket` 7.27e-02, `BucketLink2`∩`BucketLink1` 3.56e-02 — all at joints,
   tenths of a millimetre, where the modeller plainly meant an exact touch (11 pairs *are* exact
   touches, which are legal) or a designed clearance and got neither. ALICE3 is worse: **7 of 378**
   pairs in the `ST0923290` sub-assembly alone, **2.3 mm** deep, one part with **32%** of its volume
   inside a neighbour — and they are the same parts (`_013`, `_018`, `_019`) that `Stream_J_XRay.md`
   §7 named as ALICE3's only non-clean single-solid transports. `as1-oc-214.stp` is clean, 0/153.

   **The *scale* on ALICE3 is unknown.** The full 206-instance / 1699-pair census was killed at 59
   minutes without writing output; only the `ST0923290` sub-assembly is complete. The verdict does
   not depend on it — 7 > 0 already — but nobody knows how many of the other 1641 pairs are also
   bad. `Stream_T_AssemblyOracle.md` §5 calls this the highest-value hour left on that document.
   Run it **detached**, writing to `--out`, not stdout.

   **Every "0 disagreements" number on this branch is a statement about solids in isolation.** They
   remain true and they do not add up to a valid world. Fixing the source geometry is the *model
   owner's* decision, not a conversion bug — do not silently repair it.

   Two corrections travel with this. The converter's "55 solids" is **definitions**; ALICE3's world
   has **206 placed instances**, i.e. 21115 pairs. And `TGeoManager::CheckOverlaps` **does not work
   on `O2BVHSurfaceSolid`** — starved by default (`GetPointsOnSegments` returns without filling its
   buffer) and wrong when fed (it called a provably 0.09 cm gap a 0.41 cm overlap, and its
   `BasePin|Base` entry is a 24-gon's chord error to five significant figures). So there is
   currently **no working overlap check on the representation we ship**, which is ours to fix and is
   a second blocking item.

1. **Assembly-level transport, the C++ half.** The oracle now exists — `assemblyOracle.py` gives an
   ordered per-ray crossing list whose interval answer is a **set of occupants**, not a boolean, so
   touching, gaps, nesting and interpenetration all fall out of one design (28 analytic self-checks,
   0 failures). What is missing is the navigator side and the **leaking** counter. Note that on a
   geometry with item 0's overlaps, "leaking" and "overlapping" are not separable — fix or exclude
   the overlapping pairs first, or the counter measures the model rather than the code.

2. **The exact-coverage ladder for ALICE3 has one rung left, and it is expensive.** ALICE3 emits
   **20 of 55**. Everything cheap has now been measured and closed:
   - *fitted trim curves* — refused: approximate, and an independently-fitted shared edge breaks
     edge identity;
   - *co-surface half-space trims* — the naive conjunction carries **1 of 15** target solids and
     gets Bagger **1 of 13** (it ships 13/13 today); the L-shaped `l-bracket` face misses by 4.94 cm
     and needs 3 arrangement cells. `Stream_R_CoSurfaceTrims.md`;
   - *merging the exporter's co-surface patches* — **1 of 16**. The fragmentation premise was right
     (92% of wire-block faces have a neighbour on their own surface) but seams are 37.6% of all
     boundary edges and only 31.5% of the *rejected* ones: it is everywhere, and almost all of it is
     on edges that already work. `Stream_U_CoSurfaceMerge.md`.

   The only remaining exact route is **exact arrangement-cell enumeration** (`Stream_R` §9.3) — a
   union of cells rather than one conjunction, which reaches 13/15 but *cannot* be derived by
   sampling (cell counts move with grid density on 3.6% of faces; the leak verdict flips on 2.2%).
   That is research-grade work. **Absent it, the honest answer for those solids is tessellated** —
   which is the standing bargain, now evidence-backed rather than assumed.

   Two bounded exact wins do remain: `ST0923290_014` via a per-solid merge pre-pass (20 → 21), and
   `ST0923290_021`, which needs more than the three-edge guard.

3. **The tessellated fallback is not automatically sound** — see `MeshHealing.md`. A mesh can be
   *invalid*, not merely inaccurate, and chordal accuracy does not detect it. This matters more now
   that item 2 makes tessellation the answer for 35 ALICE3 solids. Mesh validity is also
   **non-monotone in precision**: refining `BucketLink2` from 0.1 to 0.05 takes it from 6600 to
   10697 LOST crossings. CGAL is already an O2 dependency, so the repair options cost a call rather
   than an adoption.

4. **`Curve2D::closestPoint` is now the hottest thing in the kernel** — 61% of `Contains`, and the
   reason a visited candidate patch costs 2.4–8.6× the corpus average. The surface solid is 40–145×
   slower than a `TGeoCompositeShape` doing identical work, and this is where that lives.
   `Stream_P_RepresentationBench.md` §2.2, `Stream_S_SafetyBVH.md`.

5. **A face-geometry gate column.** The criterion `Stream_F_EdgeIdentity.md` says edge identity
   lacks: **no face's outward normal may be antiparallel to the source face's**. It caught a defect
   closure, edge identity and the sampling gate all missed — three parts read `reliable` with nine
   inward-pointing faces. **Closure and edge identity are sign-blind.** Still not a harness flag or
   a gate column; it should be both. *(Flagged in three consecutive hand-overs now.)*

6. **The oTOF corpus is unreachable from the converter, and it is the one that would pay best.**
   `O2_CADtoTGeo.py` sees **3** leaf solids in `oTOF System V3-R92cm.step`, not the 20 prototypes in
   62628 placements `Stream_A_CSG.md`'s census reports, and without `--mesh` it does not even
   complete — `triangulate_asbbox()` dies with `Standard_ConstructionError: Bnd_Box is void`. The
   census walks the STEP assembly differently from the converter's own leaf-solid extraction, and
   only the converter's view can produce artefacts. ~62560 *placed boxes*, i.e. exactly the case
   `Stream_N_PlacedPrimitives.md` made cheap. Measured 2026-08-02.

7. **Free-form surfaces** — the genuine remainder: **19 of 55 solids, 1373 faces**, not the 36/2377
   the older brief states. Largest effort. Must report its own achieved tolerance honestly rather
   than claiming the exactness the analytic path has.

8. **Small, named, and unexplained.** ALICE3's 9.75 s `CloseShape` is now the largest number on that
   corpus by three orders of magnitude. `cyl_inter_cyl`'s *surface* shows 2 extra crossings / 2
   unterminated / 2 parity at 96 beams and is clean at 3 axis beams. `Bagger/Bucket` contributes 1
   extra crossing in 221184 rays. Each needs a per-face X-ray localiser that does not exist.

   *(Done and no longer here: `compareGateRuns.py` (`11ba928968`); the ALICE3 transport defect and
   the quartic guards; the ellipse trim; the placed-primitive composite; `Safety`'s missing BVH.)*


1. **A face-geometry gate column.** The ALICE3 diagnosis produced the criterion
   `Stream_F_EdgeIdentity.md` says edge identity lacks: **no face's outward normal may be
   antiparallel to the source face's**. It caught a defect that closure, edge identity and the
   sampling gate all missed — three parts read `reliable` with nine inward-pointing faces.
   **Closure and edge identity are sign-blind.** It is not yet a harness flag or a gate column;
   it should be both.
2. **The oTOF corpus is unreachable from the converter, and it is the one that would pay best.**
   `O2_CADtoTGeo.py` sees **3** leaf solids in `oTOF System V3-R92cm.step`, not the 20 prototypes
   in 62628 placements `Stream_A_CSG.md`'s census reports, and without `--mesh` it does not even
   complete — `triangulate_asbbox()` dies with `Standard_ConstructionError: Bnd_Box is void`. The
   census walks the STEP assembly differently from the converter's own leaf-solid extraction, and
   only the converter's view can produce artefacts. That corpus is ~62560 *placed boxes*, i.e.
   exactly the case `Stream_N_PlacedPrimitives.md` just made cheap. Measured 2026-08-02.

   *(`compareGateRuns.py` was fixed in `11ba928968` and is working: a Bagger pre/post diff ran
   clean and accounted for all 115 moved fields. The item that used to sit here is done.)*

   *(The ALICE3 transport defect and the quartic guards are both fixed and verified — ALICE3 is
   **13822/13822 rays identical to OpenCascade, 18/18 parts clean**, every robustness counter zero
   in both stepping modes. `Stream_L_ALICE3Defect.md`, `Stream_M_Quartic.md`.)*
3. **The trim generalisation — this is what "Tier 0" actually turned out to be.** Recognition
   already worked; ALICE3 emits **20** sidecars against **36** eligible solids and the 16 missing
   ones fail in `_recognized_quadric_wire_block`, on boundary edges that are not iso in the
   recognised (φ, h) domain. **Biggest converter-side coverage item on the board — and a decision,
   not just work.**

   **The "fitted curve" half of this item is superseded — read
   [`Stream_O_ImplicitTrims.md`](Stream_O_ImplicitTrims.md) before acting on it.** The claim that
   1053 of 1891 edges are "genuinely free-form" and that "no exact representation exists"
   reproduces only against the *pre-fix* recogniser; against the shipping converter it is 763 of
   1303. More importantly, **691 of those 763 are exactly the intersection of two analytic
   surfaces we already recognise**, and a per-edge census puts **15 of the 16 solids** entirely in
   that bucket — 443 edges, none needing a fit. The route is an *implicit / co-surface* trim, not
   a fitted B-spline. Smallest first step: an ellipse boundary on a **planar** face
   (`plane ∩ cylinder`), which is exact at machine precision and takes Bagger 12→13 and the ladder
   fixtures 9→10, i.e. it moves the two corpora the gate can actually score.
4. **Free-form surfaces** — the genuine remainder after Tier 0: **19 of 55 solids, 1373 faces**,
   not the 36/2377 the older brief states. Largest effort. Must report its own achieved tolerance
   honestly rather than claiming the exactness the analytic path has.
5. **Assembly-level transport** — several parts in one world, where *leaking between volumes*
   would show. Scoped in `Stream_J_XRay.md` §9; needs a **leaking** counter none of the current
   metrics can express.

## Things that will not be obvious

**A refuted experiment is not a refuted hypothesis.** The most expensive mistake in this project's
history, and it recurred *today* inside the X-ray benchmark itself: its first version reported a
known real defect as clean, because a **parallel beam is direction-poor** — three axis beams are
three directions however many rays you fire. Use `--beams 96`. When an experiment fails to move a
number, check that it was *capable* of moving it.

**An aggregate says *that*, never *where*.** The capacity residual survived five sessions of
correct reasoning and fell in one run once a per-face localiser existed. Build the instrument that
names the object before theorising.

**Prepend, never replace, `PYTHONPATH`/`LD_LIBRARY_PATH`.** The alibuild python3.10 that pythonOCC
needs *is* the O2 interpreter, so one process can import both `OCC` and `ROOT` — but replacing them
kills PyROOT on `libffi.so.6`. That bug was live in `runOracleGate.py` and silently recorded every
part one CSG tier down.

**`manifest.json` stores absolute paths.** Copying a finished workdir and re-running with
`--skip-convert` silently scores the *original* directory — no error, just a missing column that
looks like a failed emitter. `rebase_manifest()` exists now; use it.

**No `TGeoShape` in ROOT 6.36 carries a rigid transform** — `TGeoBBox` has `fOrigin`, nothing else
has anything, only `TGeoCompositeShape` holds a matrix. That used to force a placed primitive to be
emitted as its union with an identical copy of itself. **Since Stream N (placed primitives) it is not**: the shape is
written in its own canonical frame and the transform travels beside it, as a `TGeoHMatrix` under the
key `placement` in `shape_<part>.root` and as a 3x4 array in `csg_*.json` / `manifest.json` /
`gate.json`. **No placement recorded means the identity**, so every older artefact still loads
unchanged. Consumers compose it: the harness and the X-ray benchmark transform points and rays into
the shape's frame; `geom.C` places the volume with `partPlacement * shapePlacement`, in that order.
`Stream_N_PlacedPrimitives.md`.

**`TGeoCompositeShape::Capacity()` is Monte-Carlo in ROOT.** A geometrically exact composite came
back 3.3e-04 off; a capacity gate would have failed a correct shape by 330×. Composites are
accepted on OCCT symmetric-difference volume, and cross-checked by X-ray chord integration.

**`BRepGProp::VolumeProperties` on a single planar face returns exactly 0.** Documented nowhere.
Any per-face comparison must restrict itself to curved faces and account for planes by difference.

**The X-ray's chord-integral volume is a 1e-4 to 1e-5 instrument and is non-monotone in N.** It
cannot resolve the 1.3e-06 capacity residuals and must never be quoted as if it could. Its value
is gross errors and **composites**.

**Do not read `reliable` as evidence about a face's geometry.** Three ALICE3 parts reported
`reliable`, `navigable`, zero non-manifold edges — with nine faces whose outward normal pointed
*inward*. **Closure and edge identity are sign-blind.** The normal-direction check is the fix for
that blindness and is not yet in the gate.

**A recognised NURBS quadric's `face.Orientation()` says nothing about its axis.** It describes the
exporter's parametrisation. Measure the outward side; never read it off the flag. This is what
mechanism 1 of `Stream_L_ALICE3Defect.md` was.

## Traps in the environment

- The `eval` and the `export` must be in the **same shell command** as whatever runs next, or
  anything linking O2 dies on `libarrow_acero.so`. Applies even to `--help`.
- Put `$B/stage/lib:$B/stage/lib64` **first** on `LD_LIBRARY_PATH`, or stale installed libs win.
- Test binary `$B/stage/tests/o2-test-detectorsbase-BVHSurfaceSolid`; harness and X-ray in
  `$B/stage/bin/`.
- `--skip-convert` reuses `<workdir>/db` — valid only when the **converter** did not change.
- The gate exits non-zero whenever any part fails. Normal state. Read the counts.
- **Never** write generated artifacts into `STEP_examples/`. Use a scratch folder.
- Convert ALICE3 **without `--mesh`** (24 s / 581 MB). With meshing, one 2 m sphere reached
  **22.9 GB** resident and was killed. Tessellation, not the exact path, is the scaling problem.
- A worktree created by the agent tooling branches from **`dev`**, not from this branch. Reset it
  onto `swenzel/bvhsurfacesolid` before working.
- Quick kernel probe: standalone `.cxx` against `$B/stage/lib -lO2DetectorsBase -lGeom` plus
  `-I$HOME/alisw/O2/Detectors/Base/src`. Most of this project's key measurements were written that
  way in minutes each — usually faster than adding a flag.
- **A long run must write to `--out`, never rely on stdout, and use a uniquely-named path.** When
  several sessions work this tree at once their scratch directory is **shared, not private**: one
  session created a directory where another had a log file, the log vanished underneath a live
  process, and a 59-minute census produced nothing. The process itself was never affected — only
  its output. Detach anything over a few minutes.
- **Only one `ninja` at a time.** The build directory is shared and two concurrent builds collide,
  which is why parallel work here is split into one build-owning stream plus Python-only streams.

## Commands

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ninja O2lib-DetectorsBase o2-test-detectorsbase-BVHSurfaceSolid \
      o2-bench-detectorsbase-solid-harness o2-bench-detectorsbase-xray
ctest -R BVHSurfaceSolid

cd $HOME/alisw/O2
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate2 \
    --model scripts/geometry/STEP_examples/Bagger.step
O2_BUILD_DIR=$B python3 scripts/geometry/runXRayBench.py --workdir /tmp/xray --fixtures --beams 96
```

Name a part's rims, edge identity and open loops:

```bash
$B/stage/bin/o2-bench-detectorsbase-solid-harness --db <workdir>/db --parts <part> \
    --points 1 --rays 1 --rims --edge-identity
```

First thing to run on any new model: compare `nFaces` in
`<workdir>/oracle/answers_<model>_<part>.json` against `nSurfaces` in `<workdir>/gate.json`.
