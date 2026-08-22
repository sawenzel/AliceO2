# The presentation website: exact CAD surfaces in the browser

This directory is a self-contained page that loads a converted CAD part's `surfaces_*.bin` sidecar
and renders the **exact trimmed analytic faces** — the same records `o2::base::LoadSurfaceSolid`
feeds to `O2BVHSurfaceSolid` — next to the triangle mesh they replace. Nothing is fetched from a
network: three.js is vendored under `vendor/`, everything else is written here, and the ray
tracing is JavaScript in a worker.

## Serving it

```
cd scripts/geometry/website
./fetch_testdata.sh /path/to/gate_fixtures2 /path/to/gate_bagger2   # once, to populate testdata/
python3 -m http.server 8231
```

then open <http://127.0.0.1:8231/>. The page needs a server rather than `file://` because it
fetches `.bin` and `.json` by relative URL and runs its workers as ES modules.

`testdata/` is **not committed** (see `.gitignore`) — it holds copies of gate output.
`fetch_testdata.sh` takes one or more gate workdirs (a directory with a `db/` in it), copies a
representative set of parts, and writes `testdata/manifest.json`, which is what the page reads.
The set is: `box`, `cyl_inter_cyl`, `torus_union_cyl`, `tube_window`, `BoomCylinderInner` and
`Bucket`.

All data loading goes through one small module, `js/data.js`. A later single-file build only has to
define `window.__INLINE_DATA__` (base64 for the binaries, objects for the JSON) before the modules
load; nothing else changes. That matters because a published Claude Artifact runs under a CSP that
blocks every external request.

## The tabs

**Mesh viewer.** The tessellation from `facets_*.bin`, flat-shaded so the faceting is visible, with
a wireframe toggle — and over it, in gold, the **exact** trimmed face boundaries evaluated from the
sidecar. Where the grey mesh cuts a corner off a gold loop, that is the tessellation error at this
part's mesh precision.

**Exact raytracer.** A progressive CPU raytrace, in a pool of workers, with six views:

| view | what it shows |
| --- | --- |
| exact surfaces | rays intersected with the real trimmed patches; the shading normal is analytic, so a cylinder shades smoothly |
| tessellation | the same camera against the triangles, with flat facet normals |
| difference: exact vs mesh | red where one representation is hit and the other is not; otherwise a heatmap of the depth difference, saturating at 1e-3 of the model size |
| difference: local vs bridge | the same rays sent to the JS port and to the real O2 kernel (see below) |
| watertightness: exact | crossings counted along the whole ray; magenta where the count is odd, i.e. the ray entered and never came out |
| watertightness: mesh | the same count against the tessellation |

**Benchmarks.** Per-part grouped bars of `nsPerCall` per navigation function per representation on
a log axis (the ratios are one to two orders of magnitude, so a linear axis would be useless), plus
the accuracy and X-ray counters as a table. It reads `scripts/geometry/website_data/*.json` first
(Track 2's output, via an `index.json` listing the files) and falls back to the **synthetic**
records in `sample_data/`, saying so on the page. A field that was not measured renders as `n/a`,
never as `0`.

**Event display.** Track 3b layer 1: batch replay of an `events.json` over the geometry — track
polylines coloured by pdg, the step points as a cloud, and a virtual screen behind the part
accumulating a radiograph. The screen is post-processing, not geometry: it can be moved along its
normal and re-binned without re-running anything. The committed replay is synthetic and labelled
so; the same gun can be re-run in the browser at higher statistics (5 × 3000 tracks takes about
0.5 s), which is what makes the part's shadow legible on the screen.

**Self-check.** Runs the assertions below in the page and prints PASS/FAIL. The same code runs from
the command line as `node tools/selfcheck.mjs`.

## The bridge to the real kernel

The raytracer's pixels come from an **engine**, and there are two:

- `LocalEngine` — the JS port in this directory, in workers. Default, offline, and the only engine
  a published copy can use.
- `RemoteEngine` — `scripts/geometry/tgeoRayService.py` on localhost, which traces the same rays
  through the real O2 kernel.

The contract is one call: `async traceRays(Float32Array rays) -> Float32Array results`, with
`n × 6` floats in (`ox, oy, oz, dx, dy, dz`, cm, unit directions) and `n × 5` floats out
(`t, nx, ny, nz, flags`; `t < 0` means no hit). Over HTTP:

```
POST /load   {"path": "<path to surfaces_*.bin or shape_*.root>"}
             -> {"ok": true, "kind": "surface"|"shape", "nSurfaces": N, "bbox": [...]}
POST /trace  raw Float32Array n*6  ->  raw Float32Array n*5
```

Because the page and the service are different origins, the service must send
`Access-Control-Allow-Origin` and answer the `OPTIONS` preflight (both content types used here are
outside the CORS safelist). The service in this repository does. If it is not running, the page
says "bridge offline" and stays on the local engine — an absent bridge is never an error.

### What the engine difference measured

Against the running service, at 480 × 294 (141 120 rays per part), every test part came back
**identical**:

| part | faces | pixels hit | hit/no-hit disagreements | max abs depth difference |
| --- | --- | --- | --- | --- |
| box | 6 | 29 498 | 0 | 0 cm |
| cyl_inter_cyl | 6 | 20 104 | 0 | 0 cm |
| torus_union_cyl | 6 | 15 312 | 0 | 0 cm |
| tube_window | 4 | 19 853 | 0 | 0 cm |
| BoomCylinderInner | 6 | 3 368 | 0 | 0 cm |
| Bucket | 97 | 16 243 | 0 | 0 cm |

The comparison is in float32 on both sides, which is the precision the contract carries; it is a
statement that the port takes the same branches and lands on the same surfaces, not that both
sides agree to the last double bit.

## The sidecar subset that is supported

The parser (`js/sidecar.js`) mirrors `Detectors/Base/src/O2SurfaceSolidIO.cxx` record for record
and accepts **versions 1, 2 and 3**, including the version-3 edge identities (read, range-checked
against the model edge table, and otherwise unused here — closure is not a question this page
asks). A version-1 file is loaded with the 1e-6 cm extractor-precision fallback tolerance, and the
page says so.

Supported surface records: **plane** (polygon and curved-boundary), **cylinder**, **cone**,
**sphere**, **torus** — all five, with and without a wire trim, with the `innerWall` flag honoured
so a bore's normal points the right way. Supported trim curves: **line**, **circular arc**,
**clamped rational B-spline**.

Any record kind outside that set is **kept, counted and reported** — in the part panel, in the
raytracer's footnote, and as a failure in the self-check — rather than silently dropped. On the six
test parts nothing is unsupported and nothing is rejected.

The kernels are ports of `Detectors/Base/src/BoundedSurface.h`:

- the stable quadratics for cylinder, cone and sphere, with `sameIntersection` graze suppression so
  crossing parity survives a tangency;
- `solveQuarticReal` for the torus, **including** the Cauchy-bound power-of-two normalisation and
  the dimensionless `32 * eps` guards, and including the rule that an unresolved resolvent routes
  to the biquadratic branch rather than to "no roots";
- trim containment by half-open scanline winding, exact for lines and for arcs (split at their
  v-extreme angles so each sub-arc is v-monotone), and against the flattened polyline for
  B-splines;
- `bsplineSampleInto`'s three interior probes, its refusal to call an interval flat across an
  interior knot, and its degenerate-chord handling;
- one canonical seam vertex per join, substituted into both neighbours, so winding and
  point-to-curve distance measure against the same polyline;
- the crossing-cluster walk (anchored on the cluster's first member) that decides enter / exit /
  graze, which is what makes the parity counter mean something.

## Known limitations

- **B-spline trims are polylines.** A rational B-spline trim curve enters this page exactly as it
  enters the kernel: as the adaptively flattened polyline at `kBSplineFlatness = 1e-5`. The
  *surfaces* stay exact; the trim boundary is accurate to that flatness, and the on-boundary band
  is widened to match, as the kernel does. The raytracer footnote names how many faces of the
  loaded part carry such a trim (10 of 97 on `Bucket`, 6 of 6 on `cyl_inter_cyl`).
- **No safety, no distance-to-in/out queries.** This page asks two questions only — nearest
  crossing and crossing count — because that is what a picture and a parity counter need.
- **The wire winding is not normalised.** The kernel re-orients a loop so an outer wire winds
  counter-clockwise; nothing here integrates a signed area, and crossing parity does not depend on
  orientation, so the loops are used in file order. The outer/inner role comes from the file.
- **The loader does not reject a wire whose join gap exceeds the band.** The kernel does, and must.
  A viewer that refuses to draw a face teaches nothing, so the gap is measured and reported instead
  (the part panel shows the worst gap next to the band it would be judged against).
- **Per-surface AABBs, no BVH over the faces.** The parts here have at most 97 faces and a slab
  test each is enough. A part with thousands of faces would want the real thing.
- **The mesh view is only as good as the mesh.** `facets_*.bin` carries no normals, so the shading
  normal is the facet's own geometric normal, oriented against the ray.

## Self-check results

`node tools/selfcheck.mjs` (and the Self-check tab, which runs the same code) — **30 passed,
0 failed**, in 0.33 s on the command line and 296 ms in the browser. The assertions are analytic,
not recorded:

```
== solveQuarticReal ==
PASS  four separated real roots                       (x-1)(x-2)(x-3)(x-4)
PASS  exact biquadratic takes the biquadratic branch
PASS  scaled biquadratic (roots +/-1e-3, +/-2e-3) is still recognised
PASS  a double root is returned as a pair             so parity clustering can drop it
PASS  a complex-only quartic returns no real root

== box (6 faces, extent 2 x 3 x 4 cm from the sidecar) ==
PASS  axial ray along x enters at 5.0000 and exits at 7.0000
PASS  axial ray along y enters at 5.0000 and exits at 8.0000
PASS  axial ray along z enters at 5.0000 and exits at 9.0000
PASS  2000 random rays match the analytic slab intersection
      worst |dt| 1.78e-15 cm, 0 missing, 0 spurious

== tube_window (cylinder r = 1.5, z in [-3, 3], bore r = 0.8 through it along x) ==
PASS  axial ray beside the bore crosses 6.0000 cm of material   (2 crossings)
PASS  a ray down the bore hits nothing                          (0 crossings)
PASS  a ray across the bore crosses 4 walls and 1.4000 cm       (4 crossings)

== parity closure, every part ==
PASS  box                2880 grid rays (800 hit)  + 1500 random rays (1500 hit):  0 odd
PASS  cyl_inter_cyl      2880 grid rays (686 hit)  + 1500 random rays (1300 hit):  0 odd
PASS  torus_union_cyl    2880 grid rays (476 hit)  + 1500 random rays (1164 hit):  0 odd
PASS  tube_window        2880 grid rays (632 hit)  + 1500 random rays (1384 hit):  0 odd
PASS  BoomCylinderInner  2880 grid rays (84 hit)   + 1500 random rays (353 hit):   0 odd
PASS  Bucket             2880 grid rays (572 hit)  + 1500 random rays (1192 hit):  0 odd
```

The worst wire join gap found while building, against the 1e-6 cm band each model would be judged
at: `box` 0, `cyl_inter_cyl` 2.4e-16, `torus_union_cyl` 9.8e-16, `tube_window` 6.2e-13,
`BoomCylinderInner` 4.2e-12, `Bucket` 3.7e-11 cm.

### Where this page does *not* reproduce a recorded number

`cyl_inter_cyl` is on record as the part whose **mesh leaks**: six rays enter the tessellation and
never leave, in the X-ray transport benchmark. **This page does not see that.** Its mesh
watertightness overlay and its self-check both report zero odd-parity rays for that mesh, and a
separate 60 000-ray scan run while writing this found zero as well (0 odd of 60 000 for
`cyl_inter_cyl`, `torus_union_cyl` and `Bucket` alike).

The two measurements are not the same question. This page counts **triangle crossings** along an
infinite ray — a statement about whether the triangle set is closed. The X-ray benchmark **steps**
through the geometry with `DistFromOutside` / `DistFromInside` on `O2Tessellated`, where a ray can
be lost to a navigation tolerance at a facet seam even though the triangle set it walks is closed.
A closed mesh that loses stepping rays is exactly the kind of defect the X-ray instrument exists to
find, and the crossing count cannot see it. So the honest statement for the talk is: the
tessellation's leak is a navigation failure, not a hole in the triangle set, and the picture that
demonstrates it is the X-ray counter table, not this overlay.

## Layout

```
index.html          the page; tabs, and nothing else
css/style.css
js/data.js          the ONE place data is loaded from (the single-file build's seam)
js/sidecar.js       surfaces_*.bin / facets_*.bin parsers, mirroring O2SurfaceSolidIO.cxx
js/quartic.js       solveDepressedCubic + solveQuarticReal
js/curve2d.js       lines, arcs, clamped rational B-splines; flattening, winding, distance
js/wire.js          trim loops, the first fundamental forms, containment
js/surfaces.js      the six surface kernels
js/solid.js         the solid: crossing clusters, first hit, parity, containment by vote
js/meshtrace.js     Moller-Trumbore over a median-split BVH
js/engine.js        LocalEngine (workers) and RemoteEngine (the bridge)
js/rtworker.js      the worker: parses the sidecar itself and answers ray batches
js/raytrace.js      camera, bands, shading, the difference and parity views
js/rtui.js          the raytracer tab's controls
js/viewer3d.js      the shared three.js view
js/orbit.js         a minimal orbit control (so no bare-specifier import map is needed)
js/charts.js        the benchmark charts, as plain SVG
js/events.js        the event-display tab
js/gun.js           the synthetic gun, shared by the sample-data tool and the worker
js/selfcheck.js     the assertions, runnable in node and in the browser
js/app.js           part selection, tabs, shared state
tools/              make_sample_data.mjs, selfcheck.mjs
sample_data/        committed synthetic stand-ins, all labelled synthetic in their own meta
vendor/             three.js r185 (MIT, licence included)
testdata/           gate output, NOT committed; see fetch_testdata.sh
```
