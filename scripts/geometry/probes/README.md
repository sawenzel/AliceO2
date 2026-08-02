# Diagnostic probes

Standalone instruments, none of them part of the build, each written because the measurement it
makes is not available anywhere else. Read the stream document first — it is the result; these are
the instruments.

| group | probes | result document |
| --- | --- | --- |
| Stream L — ALICE3 transport | `faceNormalSamples.py`, `faceNormals.cxx`, `faceAttrib.cxx`, `torusQuartic.cxx` | [`../Stream_L_ALICE3Defect.md`](../Stream_L_ALICE3Defect.md), [`../Stream_M_Quartic.md`](../Stream_M_Quartic.md) |
| Stream N — implicit trims | `trimEdgeCensus.py`, `trimEdgeCensusReport.py` | [`../Stream_O_ImplicitTrims.md`](../Stream_O_ImplicitTrims.md) |
| Stream R — co-surface trims | `implicitTrimValidate.py`, `implicitTrimValidateReport.py` | [`../Stream_R_CoSurfaceTrims.md`](../Stream_R_CoSurfaceTrims.md) |

---

## Stream L — transport-defect attribution

Nothing here is part of the build. Each `.cxx` is a standalone probe compiled the way
`Workstreams.md` §1 documents:

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH   # stage FIRST

g++ -std=c++20 -O2 -w -o faceNormals faceNormals.cxx \
    -I$HOME/alisw/O2/Detectors/Base/include -I$HOME/alisw/O2/Detectors/Base/src \
    -I$HOME/alisw/O2/Detectors/Base/test -I$(root-config --incdir) \
    -I$ALIBUILD_WORK_DIR/ubuntu2404_aarch64/nlohmann_json/latest/include \
    -L$B/stage/lib -lO2DetectorsBase $(root-config --libs) -lGeom
```

The python side needs the alibuild python3.10 with pythonOCC; `../csg/occ_env.py` is the clean
reference for that environment (**prepend**, never replace, `PYTHONPATH`/`LD_LIBRARY_PATH`).

| probe | what it measures |
| --- | --- |
| `faceNormalSamples.py` | OCCT's **oriented outward normal** at interior samples of every face of a `.brep`, in `TopologyExplorer(shape).faces()` order — the converter's own face iteration, so face *i* is sidecar record *i*. |
| `faceNormals.cxx` | the kernel's outward normal (`ComputeNormal(p, nullptr, n)`) at those same points, and the count of faces where the two are **antiparallel**. This is the check that closure and edge identity cannot make: both are blind to a global sign. |
| `faceAttrib.cxx` | the per-ray transport diff against the OpenCascade crossing list, stepping with `XRayTransport.h`'s own loop, keeping **which** reference crossing was lost or sense-inverted and what the kernel's hit list says at that distance. Its summary line must reproduce the benchmark's per-part counts. |
| `torusQuartic.cxx` | one torus record and one ray reduced to the coefficients of a single `solveQuarticReal` call, swept over the ray's start point, with the failing guard named. |

Typical run, on a workdir produced by `runXRayBench.py`:

```bash
python3.10 faceNormalSamples.py --brep  <db>/brep_<part>.brep --out /tmp/samples.json
./faceNormals --surfaces <db>/surfaces_<part>.bin --samples /tmp/samples.json
./faceAttrib  --surfaces <db>/surfaces_<part>.bin \
              --crossings <workdir>/xray/crossings_<part>.json --out /tmp/diff.json
```

**Self-check before trusting any of them**, which is how they were used: `faceNormals` must report
0 inverted faces on a part that transports cleanly, and `faceAttrib` must reproduce the benchmark's
own LOST / kind / unterminated counts for the part exactly.

---

## Stream N — why a recognised quadric's boundary edge is rejected

Converter-side, pure python, no build products involved. Read
[`../Stream_O_ImplicitTrims.md`](../Stream_O_ImplicitTrims.md) §2 for the full method.

| probe | what it measures |
| --- | --- |
| `trimEdgeCensus.py` | for every boundary edge that `_recognized_quadric_wire_block` declines as non-iso: whether it is iso after all (**A**), whether it is exactly the intersection of two analytic surfaces (**B**, with the measured max distance from *both* implicit surfaces, absolutely and normalised by edge length, patch diagonal and the edge's own BRep tolerance), whether the neighbour is free-form (**C**), or something else, named (**D**). Rolls up to "which leaf solids would emit if implicit trims existed". Also censuses the *second* mechanism — planar faces declined for an elliptical boundary edge. |
| `trimEdgeCensusReport.py` | turns that JSON into the document's tables; stdlib + numpy only, so it runs anywhere. `--compare LEGACY CURRENT` reproduces `Stream_K_Tier0.md` §2's 1891/834/1057 against the pre-fix recogniser. |

```bash
# env as above, PLUS pythonOCC prepended (see ../csg/occ_env.py)
python3 trimEdgeCensus.py --model ../ALICE_3_example/CAD_noETA.stp --fixtures \
        --json /tmp/n/all.json          # ~25 s for ALICE3
python3 trimEdgeCensusReport.py /tmp/n/all.json
```

**Two self-checks it runs itself, and both must hold before any of its output means anything.**
(i) It re-implements the shipped per-edge iso test so it can continue past the first bad edge; its
predicted *face* verdict is compared against the shipping `recognize_and_extract_face` on every
face and the mismatch count must be **0**. (ii) Its deviation instrument is calibrated on the edges
the shipped test already *accepts* — rims and generators, which are known-exact
quadric-meets-cap-plane intersections — because an instrument that only ever sees small numbers on
the population it is arguing about has not been controlled. `--b-factor` exposes the resulting
threshold and the summary prints the whole sweep.

`--legacy-recognizer` runs the census against a verbatim copy of the pre-fix
`_recognize_analytic_surface` (`git 237be7f81a^`). It exists solely to reproduce numbers measured
before `Stream_K_Tier0.md` §5 landed; the production function is never modified and the default
path calls the shipping one.

---

## Stream R — is the trimmed face the intersection of its neighbours' half-spaces?

Converter-side, pure python. Read [`../Stream_R_CoSurfaceTrims.md`](../Stream_R_CoSurfaceTrims.md)
§2 for the full method; `Stream_N`'s census says the *edges* are co-surface intersections, this asks
whether the resulting *containment rule* is the same set as the face.

| probe | what it measures |
| --- | --- |
| `implicitTrimValidate.py` | per face: samples the face's own surface on an `N x N` grid over its `(u, v)` rectangle, classifies every sample with **`BRepTopAdaptor_FClass2d`** (ground truth) and with the half-space conjunction (the rule), and reports false positives / false negatives with the 3D distance to the boundary and the depth on the wrong side. Also reports the **arrangement-cell** structure — `cellsIn` (how many DNF terms an exact description would need) and `leak` (points no DNF over those surfaces can reject) — both independent of any sense convention. Censuses the failure modes: a neighbour on the face's own surface, a free-form neighbour, a seam, a hole, one surface bounding twice, one surface needing both senses, tangency. |
| `implicitTrimValidateReport.py` | the document's tables; stdlib + numpy only. `--per-solid`, `--solids`, `--split-population`. |

```bash
# env as for Stream N
python3 implicitTrimValidate.py --model ALICE3:../ALICE_3_example/CAD_noETA.stp \
        --grid 32 --far-grid 16 --json /tmp/r/alice3.json        # ~5 min for ALICE3
python3 implicitTrimValidateReport.py /tmp/r/alice3.json --per-solid
```

**Three self-checks, and the numbers mean nothing without them.** (i) `--flip-sense k` inverts one
stored sense per face and `--perturb-radius R` displaces every trimming surface by R cm; both must
move the disagreement counts, and a **10 µm** perturbation does. (ii) `--verify N` re-checks false
positives against `BRepExtrema_DistShapeShape`, which shares no code with `FClass2d`; 268 of 268
checked were confirmed genuinely off the face. (iii) `--solid-crosscheck N` compares the ground truth
against `BRepClass3d_SolidClassifier`. Sweep `--grid` before believing any *positive* result: the
instrument can prove a face wrong but not prove one right (§5.5).
