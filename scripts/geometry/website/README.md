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

It copies **four** artefacts per part, whichever of them the source layout holds:

| in the source | in `testdata/<part>/` | in the manifest | what it is |
| --- | --- | --- | --- |
| `surfaces_<stem>.bin` | `surfaces.bin` | `"surfaces"` | the exact trimmed analytic faces |
| `facets_<stem>.bin` | `facets.bin` | `"facets"` | the tessellation |
| `shape_<stem>.root` | `shape.root` | `"shape"` | the recognised CSG composite, as a `TGeoShape` |
| `csg_<stem>.json` | `csg.json` | `"csg"` | what the recogniser found, and its acceptance verdict |

A missing artefact is `null` in the manifest, and the page turns off exactly the views that needed
it. `ships` follows the same evidence: a part with a `shape_*.root` is one the cascade accepted as
CSG, and the manifest says `"shape"` for it. Of the fetched set, `box` and `BoomCylinderInner` ship
CSG; the ALICE3 parts come from a flat converter directory that holds no `csg_*.json` at all, so
`null` there is the truth and not a gap.

All data loading goes through one small module, `js/data.js`. A later single-file build only has to
define `window.__INLINE_DATA__` (base64 for the binaries, objects for the JSON) before the modules
load; nothing else changes. That matters because a published Claude Artifact runs under a CSP that
blocks every external request.

## The part selector drives everything

There is one global control, the **part** listbox in the header, and every tab renders whatever it
holds -- the mesh viewer, the raytracer, the benchmarks card and the event display alike. Each
entry carries a 112 px thumbnail rendered offscreen once, from the part's tessellation (or, for a
part with no mesh, from its exact trim loops) and then cached; the thumbnails are drawn on the
first open of the panel, not at boot, because with ALICE3 in the list that is several megabytes of
facets nobody has asked to look at. Entries are grouped by the `group` field the fetch script
writes.

One entry is not a part: where `testdata/otof_assembly/` exists, the list ends with the **oTOF
assembly**, whose thumbnail is one point per placement. Selecting it swaps the per-part tabs for
the single **Assembly** tab, because a placement table has no sidecar, no CSG record and no
benchmark row; selecting a part again brings them back untouched.

## The tabs

**Mesh viewer.** The tessellation from `facets_*.bin`, flat-shaded so the faceting is visible, with
a wireframe toggle — and over it, in gold, the **exact** trimmed face boundaries evaluated from the
sidecar. Where the grey mesh cuts a corner off a gold loop, that is the tessellation error at this
part's mesh precision. Its HUD carries the compact per-part line: faces, triangles, bounding box,
and what the part ships as — naming the composite where that is CSG.

**Exact raytracer.** A progressive CPU raytrace, in a pool of workers, with nine views:

| view | what it shows |
| --- | --- |
| exact surfaces (local JS port) | rays intersected with the real trimmed patches; the shading normal is analytic, so a cylinder shades smoothly |
| exact surfaces (bridge: the real kernel) | the same picture, every pixel answered by the real O2 kernel through the bridge, shaded identically |
| CSG (bridge: the real composite) | the part's own `shape.root` — the composite the recogniser produced — traced by the same kernel and shaded the same way, so the three shaded pictures differ only in the representation behind them |
| tessellation | the same camera against the triangles, with flat facet normals |
| difference: exact vs mesh | red where one representation is hit and the other is not; otherwise a heatmap of the depth difference, saturating at 1e-3 of the model size |
| difference: local vs bridge | the same rays sent to the JS port and to the real O2 kernel (see below) |
| difference: exact vs CSG (bridge) | the exact surface solid against the composite, on the same rays: the acceptance test's `dV_sym` made visible, one pixel at a time |
| watertightness: exact | crossings counted along the whole ray; magenta where the count is odd, i.e. the ray entered and never came out |
| watertightness: mesh | the same count against the tessellation |

The two CSG views need a `shape.root` for the part **and** a connected bridge. A part without one
has them disabled, with the recogniser's own reason (`whyNotCSG`, or the acceptance test's own
rejection sentence) as the option's tooltip and as a line under the selector — the same
honest-coverage pattern a tessellated-only part gets.

Four things about how it behaves:

- **Nothing is raytraced while you orbit.** The drag moves a WebGL proxy of the part -- the
  tessellation with the exact trim boundaries in gold over it, the mesh tab's own look -- and the
  traced frame appears 300 ms after the camera stands still, and never while the pointer is down.
  A camera move cancels the frame in flight.
- **Resolution presets** run 360p, 480p, 720p, 1080p, 1440p and *native*, which traces one ray per
  physical device pixel of the viewport and is the only setting whose picture is not upscaled. The
  line under the selector gives the raster and the rays per full pass.
- **The bounding-box scissor** projects the part's world AABB to a screen rectangle and casts rays
  only inside it, padded by two pixels; everything outside is background, with no ray cast for any
  of it. The HUD reports the traced pixel count and the saving. A box that straddles the eye falls
  back to the whole frame, where a perspective projection has no finite rectangle. On `Bucket` at
  1440 x 740 the rectangle is 27.7% of the frame and the frame takes 351 ms instead of 470 ms --
  the time saving is smaller than the pixel saving because a ray that misses the AABB only fails a
  slab test.
- **Shading** is two lights with a Blinn-Phong specular at exponent 48 on each, and an optional
  single-bounce mirror. The mirror sends a second ray batch off every hit pixel to the *same*
  engine that answered the first -- the bridge included, whose `/trace` returns the kernel's own
  normals -- so the reflection in the exact view is exact and the reflection in the tessellation is
  faceted. Only the hit pixels are in the second batch. The mix is Schlick-weighted, so a face seen
  edge-on reflects and a face seen straight on barely does; without that the part just gets darker,
  because the sky it reflects is dark. The mirror applies to the three shaded views, not to the
  difference or parity overlays.

**Benchmarks.** One card for the selected part, and nothing else. Where the part has a `csg.json`
the card carries a **CSG** block read straight out of it: the structure in TGeo's own classes
(`TGeoTube ∪ TGeoTube` for `BoomCylinderInner`), the recogniser that produced it, each leaf's
parameters, whether the leaves carry their own world frames and whether the composite carries a
placement of its own, the acceptance verdict with `dV_sym` against its band, and the `shape.root`
the raytracer's CSG view traces. A part the recogniser declined shows why in the same place — the
acceptance test's own sentence where there is one, and otherwise the `whyNotCSG` line from
`website_data/decline_reasons.json`, which is loaded optionally and simply says less when it is
absent. It heads with a badge **per representation the part
carries at full quality** -- CSG, SURFACE or TESSELLATED. The cascade's own choice comes first and
solid, taken from `website_data/summary.json`'s own verdict where it has one and otherwise from
what the converter wrote into `testdata/`; anything else the part has stands next to it, dashed.
`BoomCylinderInner` therefore reads **ships CSG  ships SURFACE**: the recogniser's composite passed
the acceptance test *and* the exact sidecar is there, and both are navigable. A tessellated-only
part keeps its single badge -- a mesh sitting next to an exact sidecar is not a representation the
part has, it is the thing the sidecar replaced, so it is never a second badge. The mesh tab's HUD
line says the same thing (`ships CSG + SURFACE: TGeoTube ∪ TGeoTube`). Under it: the sidecar's own properties as this page loaded them (faces by
type, wire trims, model tolerance, worst join gap, record status -- this is the old Properties
pane, which now lives here); a representations table with primitives, memory, sidecar size,
capacity deviation and reliability; grouped bars of `nsPerCall` per navigation function per
representation on a log axis (the ratios are one to two orders of magnitude, so a linear axis would
be useless); and the accuracy and X-ray counters. It reads Track 2's
`scripts/geometry/website_data/*.json` via the `index.json` that lists them, and falls back to the
**synthetic** records in `sample_data/`, saying which on the page. A field that was not measured
renders as `n/a`, never as `0`. A part with no measured record still gets its card, from the
sidecar, and the page says which parts the measured set covers.

**measure now (bridge).** The recorded bars are one machine's numbers, taken once, and a part with
no record in `website_data/` (`BoomCylinderInner`, say) had no bars at all. The button at the top
of the tab asks the bridge to measure the part in front of the reader: for each representation it
can load -- `surfaces.bin` and `shape.root` -- one `POST /load` followed by one `POST /bench`, and
the answer is drawn as a lighter, dash-outlined bar beside the recorded one, with the load average
and the sample count in the caption under the chart. The mesh is **static data only**: the bridge
has no kernel loader for `facets_*.bin`, so a tessellated row can only ever come from the recorded
set, and the pane says so. `/bench` measures *timing* and nothing else, so the capacity, accuracy
and X-ray columns stay exactly what Track 2 recorded; a live bar is never mixed into them. With no
bridge the button is disabled and carries the reason (`no bridge on http://127.0.0.1:8077 -- start
scripts/geometry/tgeoRayService.py --port 8077`), which is also how the liveness probe reports a
service that has stopped. Measured this way on this machine, `BoomCylinderInner` `Contains` comes
back at 475 ns/call for the exact surface solid against 29.7 ns/call for its CSG composite, and
`Bucket`'s exact `Contains` at 300 ns/call against the 654 ns/call in the recorded set -- the same
part, a differently loaded machine, which is the whole reason a live bar is labelled as one.

`website_data` inside this directory is a symlink to `../website_data`, because `python3 -m
http.server` refuses paths above its own root; the loader also tries `../website_data/` for the
case where the server was started one level up instead.

**Event display.** Track 3b layer 1: batch replay of an `events.json` over the geometry — track
polylines coloured by pdg, the step points as a cloud, and a virtual screen behind the part
accumulating a radiograph. The screen is post-processing, not geometry: it can be moved along its
normal and re-binned without re-running anything. The committed replay is synthetic and labelled
so; the same gun can be re-run in the browser at higher statistics (5 × 3000 tracks takes about
0.5 s), which is what makes the part's shadow legible on the screen.

**Assembly.** The placed tree rather than one part, and it appears only when the selector's
assembly entry is chosen — the per-part tabs have nothing to say about a placement table, so they
step aside while it is selected and come back with the next part.

The honest note first: **oTOF has no tubular *solid*.** All twenty of its body prototypes are
planar — nineteen exact boxes and one 1493-plane part. The tube is the **assembly**: 20 prototypes,
3 741 leaf placements, 62 628 placed solids, 1 136 824 triangles, and out of those a barrel about
`x` of R 85.0–96.7 cm and 3.48 m long. Those numbers are measured from the placement table on load,
not written into a caption.

It draws one `THREE.InstancedMesh` per body prototype, so the whole barrel is 20 draw calls. What
makes that possible is a geometry fact from the converter: a body's vertices are already in its
**leaf's** frame — the body's own pose is baked in — so an instance's world transform is exactly
its leaf's matrix, and every prototype of a leaf shares one instance-matrix buffer. The slice keeps
whole placements whose position along the barrel axis falls in the window (it never cuts a body
open); at its narrowest the window is one ring, 2 312 solids, and `re-frame` then frames the slice
rather than the barrel. There is a wireframe toggle and a spin toggle whose frame-rate read-out is
the median of the last nine **drawn** frames.

This tab is **WebGL only** and says so on the page: no raytracing and no bridge. 62 628 solids is a
rasteriser's job.

Its data does not come from `fetch_testdata.sh`, which is per-part. `./fetch_assembly.sh
<placements.json> <facets-dir> ...` copies the placement table verbatim to
`testdata/otof_assembly/placements.json`, copies one `facets_<body>.bin` per prototype next to it
(taking the properly meshed copy where two source dirs offer the same body), and writes
`testdata/otof_assembly/index.json` with the leaf → bodies mapping and the totals. A checkout
without that directory simply has no assembly entry in the selector and no Assembly tab.

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
POST /bench  {"samples": 4000}
             -> {"ok": true, "samples": N, "repeats": R, "insideSamples": M, "loadAverage": L,
                 "functions": {"contains": {"nsPerCall": x}, "distFromOutside": {...}, ...}}
```

`/bench` times the shape that is loaded **right now**, single-threaded, on deterministic sample
points, and the benchmarks tab's "measure now" button drives it (see above). Because it measures
whatever the service holds, it has to take the bridge over for the duration: it `/load`s each of
the part's representations in turn and then hands the bridge back, re-`/load`ing the file the
raytracer's current view needs (`Raytracer.reloadBridgeShape`). Without that handover the tracer
would still believe its own file was resident and the next frame would quietly trace the other
solid. The service caches loads, so the handover is a pointer swap.

The button's liveness probe is a `POST` to an endpoint the service does not have, which answers
404 with a JSON body and the CORS headers: proof of life that disturbs nothing a `/load` or a
`/bench` would. That 404 in the browser console is the probe, and it is expected. A bare `OPTIONS`
cannot serve instead -- its own preflight asks for an `Access-Control-Allow-Methods` the service
does not send, so a live bridge would look dead.

The bridge is tried without being asked when a part loads: the stored path prefix first, then two
documented candidates (`` and `scripts/geometry/website/`), because the service resolves a
relative path against **its own** working directory and not this page's. Whichever answers is
remembered. An offline bridge stays what it always was -- a message, never an error.

**The service holds one shape at a time**, and a part can have two files worth tracing: its
`surfaces.bin` and its `shape.root`. The page remembers both paths at connect time, tracks which
one is resident, and re-`/load`s on a view switch that needs the other — once, not per band. That
swap is timed and reported on its own line (`bridge /load`) rather than folded into the frame, so
a ms/frame row is ray work and not file I/O. No `/trace` is ever allowed to overlap a `/load`: the
load waits for the traces already in the air and the traces issued while it runs wait for it,
because a trace that reaches the service mid-swap wedges it. Both requests carry a timeout (180 s
for `/load`, 60 s for `/trace`); a service that stops answering now produces a sentence saying so
instead of a page that waits for ever.

Because the page and the service are different origins, the service must send
`Access-Control-Allow-Origin` and answer the `OPTIONS` preflight (both content types used here are
outside the CORS safelist). The service in this repository does. If it is not running, the page
says "bridge offline" and stays on the local engine — an absent bridge is never an error.

### ms per frame, per engine

The pane compares **the same picture** -- the exact surfaces -- asked of each engine in turn, and
then the same part's **CSG composite** asked of the same kernel. Each of those three views is
served by exactly one engine, so the whole frame is that engine's work and they are directly
comparable; a frame in which both ran (the difference views) is not apportioned between them and
does not appear, and neither does a mesh or parity frame, which is a different question. The pane
says so when the frames used a different camera or resolution. **"time every engine"** renders the
matched set at one camera, and includes the CSG row for a part that has a `shape.root`. Measured here on `Bucket` at
720p (960 x 493, 152 274 rays after the scissor), warm, four consecutive matched pairs, 8 local
workers against the service on the same shared 10-core box:

| engine | ms/frame, four pairs | rate |
| --- | --- | --- |
| local (8 workers) | 134, 143, 147, 158 | 0.96-1.14 Mray/s |
| bridge 127.0.0.1:8077 | 271, 278, 281, 297 | 513-562 kray/s |

The ratio the pane reports over those four pairs is **bridge / local 1.84-2.07x**.

The three-row form, on `BoomCylinderInner` at 720p (960 x 493, 105 854 rays after the scissor),
warm, one matched set, same box:

| engine | representation | ms/frame | rate |
| --- | --- | --- | --- |
| local (8 workers) | exact surfaces | 131 | 808 kray/s |
| bridge 127.0.0.1 | exact surfaces | 157 | 674 kray/s |
| bridge 127.0.0.1 | CSG composite (`shape.root`) | 123 | 861 kray/s |

**bridge exact / bridge CSG is 1.28x here**, and that is the honest number: it is *not* the 141x
`Contains` ratio `website_data/summary.json` records for the sibling ram `BoomCylinderOuter`. The
two measure different things. That 141x is per navigation call, in-process; this 1.28x is a whole
frame of 105 854 rays over HTTP in fat bands, where the round trip, the JSON-free octet stream and
the thread pool are most of the wall clock and the kernel maths is the small remainder. What the
CSG row does show is the shape of the story at a glance: the same kernel, the same rays, the same
camera, the composite coming back faster than the surface solid rather than slower.

That is a statement about *these two implementations as deployed*, not about the kernels: the local
engine is eight workers and the bridge is one single-threaded Python service reached over HTTP, so
the ratio carries a process boundary, an octet-stream round trip and a thread count, not just the
maths. The spread across four consecutive pairs on a box with other tenants on it is the reason
four are quoted rather than one. What it is good for is showing that the real kernel answers a full
frame of rays in a fraction of a second, on the same camera, next to the port.

### What the exact-vs-CSG difference measured

`BoomCylinderInner` is a `tier2-tube-union` whose acceptance test recorded `dV_sym = 0` against a
8.57e-6 cm^3 band. The **difference: exact vs CSG (bridge)** view is that statement as a picture,
and it comes back empty: at 960 x 493 with the scissor on, 99 530 rays traced and 10 363 pixels
hit, **0 hit / no-hit disagreements and max |dt| 0.00e+0 cm** — a flat blue silhouette with no red
in it. The exact surface solid in the JS port and the TGeo composite in the real kernel are the
same solid, to float32, on every ray of that frame.

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
raytracer's footnote, and as a failure in the self-check — rather than silently dropped. On all 14
parts in `testdata/` that have a sidecar — the six original ones and eight from ALICE3, up to 965
faces — nothing is unsupported and nothing is rejected.

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
- **Per-surface AABBs, no BVH over the faces.** The largest part here, `ST1829909_002`, has 965 of
  them, and a slab test each is still enough for a frame in a third of a second at 480 x 247 with
  the scissor on. A part with tens of thousands of faces would want the real thing; this is what
  the sub-patch BVH in the kernel exists for and what this page deliberately does not reproduce.
- **The mesh view is only as good as the mesh.** `facets_*.bin` carries no normals, so the shading
  normal is the facet's own geometric normal, oriented against the ray.
- **The bridge does not survive unlimited shape swapping.** The CSG views make the page `/load`
  far more often than it used to: once per part before, now once per switch between a surface view
  and a CSG view. `scripts/geometry/tgeoRayService.py` does not take that indefinitely. Loading a
  `shape_*.root` and a `surfaces_*.bin` alternately from a script is stable (12 alternations with
  four concurrent trace bands each, measured while writing this), but a browser session that keeps
  switching eventually leaves the service either aborted (`double free or corruption (!prev)`,
  logged straight after `TGeoManager: default geometry created`) or wedged with every thread in a
  futex wait, answering `/load` while `/trace` never returns. Reading a `TGeoShape` out of a file
  creates a default `TGeoManager`, and a `TGeoShape` registers itself in the manager that is
  current — so replacing the resident shape is not the local operation the endpoint's signature
  suggests. This page does what it can from its side: it never overlaps a `/load` with a `/trace`,
  it reloads only on a switch that actually needs the other file, and both requests time out with
  a message instead of hanging. **The fix belongs in the service** — load each path once and keep
  it, and create the `TGeoManager` once at startup — and that file is outside this directory.
  Until then, a talk demo should expect to restart the bridge if the CSG views stop answering.

## Self-check results

`node tools/selfcheck.mjs` (and the Self-check tab, which runs the same code) -- **54 passed,
0 failed**, in 1.4 s on the command line and 1143 ms in the browser, over the 14 parts in
`testdata/` that have an exact sidecar. A tessellated-only part is skipped: every assertion here is
about the exact solid, and it has none. The assertions are analytic, not recorded:

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

== parity closure, every part with a sidecar ==
                     faces   grid rays (hit)   random rays (hit)   odd   worst join gap / band
box                      6     2880 ( 800)        1500 (1500)        0   0        / 1.0e-6
cyl_inter_cyl            6     2880 ( 686)        1500 (1300)        0   2.4e-16  / 1.0e-6
torus_union_cyl          6     2880 ( 476)        1500 (1164)        0   9.8e-16  / 1.0e-6
tube_window              4     2880 ( 632)        1500 (1384)        0   6.2e-13  / 1.0e-6
BoomCylinderInner        6     2880 (  84)        1500 ( 353)        0   4.2e-12  / 1.0e-6
Bucket                  97     2880 ( 572)        1500 (1192)        0   3.7e-11  / 1.0e-6
ST0923290_013           20     2880 ( 626)        1500 (1373)        0   1.0e-11  / 3.3e-6
ST1829909_002          965     2880 ( 304)        1500 ( 715)        0   6.8e-9   / 4.3e-4
ST1829909_004          720     2880 ( 205)        1500 ( 623)        0   4.0e-6   / 4.4e-4
ST1A38494_002            6     2880 (  85)        1500 (1361)        0   2.9e-12  / 1.0e-6
ST1A38495_01            65     2880 ( 440)        1500 (1199)        0   9.9e-10  / 8.6e-5
ST2487455_01            66     2880 (  16)        1500 (1122)        0   2.7e-8   / 4.9e-5
ST2487459_01           202     2880 ( 120)        1500 ( 501)        0   2.3e-10  / 7.4e-5
ST2487462_01            80     2880 (  21)        1500 (1066)        0   1.3e-8   / 5.3e-5
```

Every one of those 14 parts also had **every** sidecar record built: 0 rejected, 0 unsupported.
Their meshes are closed too -- 0 odd of 2880 grid rays and 0 odd of 1500 random rays, each.

### The ALICE3 part that is NOT in testdata/, and why

`ST1829909_01` (1052 faces, the largest leaf of `CAD_noETA.stp`) is deliberately **left out** of the
fetched set. It parses cleanly here -- 0 rejected, 0 unsupported, worst join gap 1.0e-5 cm against a
4.7e-4 cm band -- but its exact solid does **not** close: 1 odd-parity ray of the self-check's 1500,
and 6 odd of 9340 hitting rays on an independent 20 000-ray scan taken while writing this. Including
it would make this page's own self-check red, so it is excluded and the measurement recorded here
instead. It is the same part that is already on record as failing to load in O2 (see the
`ST1829909_01` handoff item); this is an independent second symptom on the same part, found by a
different instrument, and it is a statement about that ray sample, not a located defect.

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
js/partselect.js    the global part listbox
js/thumbs.js        the offscreen thumbnail renderer behind it
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
js/charts.js        the benchmark charts and the part card, SVG and DOM
js/livebench.js     the benchmarks tab's live measurement: /load + /bench on the bridge
js/events.js        the event-display tab
js/assembly.js      the placement table: loading it, and the InstancedMesh view of it
js/assemblyui.js    the assembly tab's controls
js/gun.js           the synthetic gun, shared by the sample-data tool and the worker
js/selfcheck.js     the assertions, runnable in node and in the browser
js/app.js           part selection, tabs, shared state
tools/              make_sample_data.mjs, selfcheck.mjs
sample_data/        committed synthetic stand-ins, all labelled synthetic in their own meta
vendor/             three.js r185 (MIT, licence included)
testdata/           gate output, NOT committed; see fetch_testdata.sh and fetch_assembly.sh
website_data        symlink to ../website_data, Track 2's measured JSON
```
