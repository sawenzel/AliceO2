# Stream L diagnostic probes

Built to attribute the ALICE3 transport losses to individual source faces; kept because the
measurement they make is not available anywhere else. Read
[`../Stream_L_ALICE3Defect.md`](../Stream_L_ALICE3Defect.md) first — it is the result; these are
the instruments.

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
