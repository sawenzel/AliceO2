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
| `ctest -R BVHSurfaceSolid` | **86 cases**, green |
| ladder fixtures | **9/9 scored** — of **10** leaf solids; `oblique_cut_cyl` has no sidecar and never has |
| `Bagger.step` | **12/12 shipped**, 9/12 surface, 1 unscored (`Bucket`, ships as mesh) |
| unexplained oracle disagreements | **0/0/0/0** on `surface`, both corpora; **0/0/0/0** on `shape` |
| `runOracleGate.py --self-test` | 17/17 |
| `o2-bench-detectorsbase-xray --self-test` | 19/19 |

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

## Open, in the order I would take them

1. **The ALICE3 defect** *(in progress at the time of writing)*. The X-ray found that **4 of 18
   converting ALICE3 parts lose 418 crossings, invert the sense of 236, leave 70 transports
   unterminated and contradict `Contains()` on 336 intervals — while all 18 report `reliable` and
   `navigable` with zero non-manifold edges.** Three of the four are the `ST0923290` family. This
   is the first concrete instance of the caveat `Stream_F_EdgeIdentity.md` states in writing:
   identity certifies that the source *topology* survived conversion, **not** each face's
   *geometry*. Diagnosis goes in `Stream_L_ALICE3Defect.md`.
2. **The quartic guards.** `solveQuarticReal` compares quantities of dimension **L²** and **L³**
   against `kTolerance`, a **length** (1e-9 cm). Two live consequences: the resolvent guard fails
   and a ray **silently misses a torus it does cross**; the `termQ` guard misroutes an asymmetric
   quartic into the biquadratic branch and returns **two confidently wrong roots instead of four
   right ones**. Not a small-geometry curiosity: 1 ray in 5162 is already lost at the shipped
   fixture size. Bagger has 2 toroidal faces of 288, ALICE3 has 350.
3. **Tier 0** — decode the NURBS-encoded quadrics. Takes ALICE3 quadric-only solids **15/55 →
   36/55**, `as1-oc-214` **0/5 → 5/5**. Helps the exact path and the CSG path equally, because it
   changes what a solid *is* before either looks at it. Highest-value coverage item.
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
has anything, only `TGeoCompositeShape` holds a matrix. A *placed* primitive is emitted as its
union with an identical copy, so a recognised plain tube reads as a composite in `gate.json`.

**`TGeoCompositeShape::Capacity()` is Monte-Carlo in ROOT.** A geometrically exact composite came
back 3.3e-04 off; a capacity gate would have failed a correct shape by 330×. Composites are
accepted on OCCT symmetric-difference volume, and cross-checked by X-ray chord integration.

**`BRepGProp::VolumeProperties` on a single planar face returns exactly 0.** Documented nowhere.
Any per-face comparison must restrict itself to curved faces and account for planes by difference.

**The X-ray's chord-integral volume is a 1e-4 to 1e-5 instrument and is non-monotone in N.** It
cannot resolve the 1.3e-06 capacity residuals and must never be quoted as if it could. Its value
is gross errors and **composites**.

**Do not read ALICE3's `reliable` as evidence.** Item 1 above is exactly the counterexample.

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
