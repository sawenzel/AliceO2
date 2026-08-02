# NEXT — session-start instruction for the CAD → TGeo work

This file is the current hand-over. Whoever finishes a session should **rewrite it**.
Last updated 2026-08-01, end of a long day: wave 0, wave 1, the MVP, and the X-ray benchmark.

Branch `swenzel/bvhsurfacesolid`. Everything below is committed.

## Read this first, and possibly only this

**[`Tutorial.md`](Tutorial.md)** is the orientation document: the three representations and the
cascade between them, the pipeline, the artefacts, the four commands, the seven principles, and a
per-capability status with its explicit gaps. It was written because five investigation records
and two reviews described this branch and none of them answered *"what exists, how do I drive it,
what is missing"*. Start there. This file is only the hand-over on top of it.

## Where the branch stands

| | |
| --- | --- |
| `ctest -R BVHSurfaceSolid` | **97 cases**, green |
| ladder fixtures | **9/9 scored** — of **10** leaf solids; `oblique_cut_cyl` has no sidecar and never has |
| `Bagger.step` | **12/12 shipped**, 9/12 surface, 1 unscored (`Bucket`, ships as mesh) |
| unexplained oracle disagreements | **0/0/0/0** on `surface`, both corpora; **0/0/0/0** on `shape` |
| `runOracleGate.py --self-test` | 17/17 |
| `O2_CADtoTGeo.py --self-test` | **26 checks**, 0 failures (18 recognition + 8 placement) |
| `csg/emit.py --self-test` | 33/33 |
| `o2-bench-detectorsbase-xray --self-test` | 17 checks, 0 failures |

**Gate totals and disagreement counts are separate numbers. Never quote one without the other.**

Bagger now converts as a **mixture**: CSG 7, exact surfaces 5, tessellated 1, of 13 leaf solids —
every CSG part exact (`dV_sym = 0`) and oracle-clean. That was the MVP and it is done.

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

## Open, in the order I would take them

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
