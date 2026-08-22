// Reader for the exact-surface sidecar format (surfaces_*.bin, versions 1 to 3).
//
// A direct port of Detectors/Base/src/O2SurfaceSolidIO.cxx, which is the authoritative reader;
// the layout is documented in scripts/geometry/BVHSurfaceSolid.md ("Surface sidecar format").
// Everything the C++ reader rejects is rejected here, and with the same message, so a file that
// loads in the browser is a file the kernel loads too.

export const SURFACE_TYPE = { PLANE: 1, CYLINDER: 2, CONE: 3, SPHERE: 4, TORUS: 5 };
export const CURVE_TYPE = { LINE: 0, ARC: 1, BSPLINE: 2 };

export const SIDECAR_VERSION_MIN = 1;
export const SIDECAR_VERSION_MAX = 3;
// What a version-1 sidecar's model tolerance is taken to be, in cm (kSidecarV1FallbackTolerance).
export const SIDECAR_V1_FALLBACK_TOLERANCE = 1e-6;
const FLAG_INNER_WALL = 1;

const TYPE_NAME = { 1: 'plane', 2: 'cylinder', 3: 'cone', 4: 'sphere', 5: 'torus' };
export function surfaceTypeName(t) { return TYPE_NAME[t] || `unknown(${t})`; }

// Expected parameter counts, per the Add*Surface signatures.
const EXPECTED_PARAMS = { 1: 9, 2: 14, 3: 15, 4: 14, 5: 15 };

class Cursor {
  constructor(buffer) {
    this.view = new DataView(buffer);
    this.bytes = new Uint8Array(buffer);
    this.offset = 0;
    this.size = buffer.byteLength;
  }
  get remaining() { return this.size - this.offset; }
  u32() {
    if (this.remaining < 4) { throw new Error('truncated (uint32)'); }
    const v = this.view.getUint32(this.offset, true); this.offset += 4; return v;
  }
  u8() {
    if (this.remaining < 1) { throw new Error('truncated (uint8)'); }
    return this.bytes[this.offset++];
  }
  f64() {
    if (this.remaining < 8) { throw new Error('truncated (float64)'); }
    const v = this.view.getFloat64(this.offset, true); this.offset += 8; return v;
  }
  doubles(n) {
    // Size every allocation against what the file can actually hold: a count of billions is a
    // truncation error, not an out-of-memory kill (the C++ reader's bytesRemaining guard).
    if (n * 8 > this.remaining) { throw new Error(`claims ${n} doubles, more than the file holds`); }
    const out = new Float64Array(n);
    for (let i = 0; i < n; ++i) { out[i] = this.view.getFloat64(this.offset + 8 * i, true); }
    this.offset += 8 * n;
    return out;
  }
}

/// Parse a bspline (curveType 2) edge record laid out as
/// [degree, nPoles, poles(2*nPoles), weights(nPoles), knots(nPoles+degree+1)].
export function parseBSplineEdge(params) {
  if (params.length < 2) { return null; }
  const degree = Math.round(params[0]);
  const nPoles = Math.round(params[1]);
  if (degree < 1 || nPoles < degree + 1) { return null; }
  const nKnots = nPoles + degree + 1;
  const expected = 2 + 2 * nPoles + nPoles + nKnots;
  if (params.length < expected) { return null; }
  let o = 2;
  const poles = new Array(nPoles);
  for (let i = 0; i < nPoles; ++i) { poles[i] = [params[o], params[o + 1]]; o += 2; }
  const weights = new Array(nPoles);
  for (let i = 0; i < nPoles; ++i) { weights[i] = params[o++]; }
  const knots = new Array(nKnots);
  for (let i = 0; i < nKnots; ++i) { knots[i] = params[o++]; }
  return { degree, poles, weights, knots };
}

/// Parse a surfaces_*.bin sidecar from an ArrayBuffer into a plain-object model:
///   { version, nSurfaces, modelTolerance, modelToleranceStated, nModelEdges, surfaces: [...] }
/// Each surface: { type, innerWall, params: Float64Array, wires: [{role, edges: [{curveType, params}]}],
///                 edgeRefs: [{id, flags}] }
export function parseSidecar(buffer, fileLabel = 'sidecar') {
  const c = new Cursor(buffer);
  if (c.remaining < 16) { throw new Error(`${fileLabel}: truncated header`); }
  const magic = String.fromCharCode(c.bytes[0], c.bytes[1], c.bytes[2], c.bytes[3]);
  c.offset = 4;
  if (magic !== 'O2SS') { throw new Error(`${fileLabel} is not a surface sidecar file (bad magic)`); }
  const version = c.u32();
  const nSurfaces = c.u32();
  const reserved = c.u32();
  if (version < SIDECAR_VERSION_MIN || version > SIDECAR_VERSION_MAX) {
    throw new Error(`${fileLabel}: unsupported sidecar version ${version} (reader supports ${SIDECAR_VERSION_MIN}..${SIDECAR_VERSION_MAX})`);
  }
  let modelTolerance = SIDECAR_V1_FALLBACK_TOLERANCE;
  let modelToleranceStated = false;
  let nModelEdges = 0;
  const warnings = [];
  if (version >= 2) {
    modelTolerance = c.f64();
    modelToleranceStated = true;
    if (version >= 3) { nModelEdges = c.u32(); }
  } else {
    warnings.push(`${fileLabel} is a version-1 sidecar and states no model tolerance; assuming ${SIDECAR_V1_FALLBACK_TOLERANCE} cm (the extractor's precision).`);
  }

  const surfaces = [];
  for (let s = 0; s < nSurfaces; ++s) {
    const type = c.u32();
    const flags = c.u32();
    const nParams = c.u32();
    const params = c.doubles(nParams);

    const nWires = c.u32();
    if (nWires * 8 > c.remaining) {
      throw new Error(`${fileLabel}: surface ${s} claims ${nWires} wires, more than the file holds`);
    }
    const wires = [];
    for (let w = 0; w < nWires; ++w) {
      const role = c.u32();
      const nEdges = c.u32();
      if (nEdges * 8 > c.remaining) {
        throw new Error(`${fileLabel}: surface ${s} claims ${nEdges} wire edges, more than the file holds`);
      }
      const edges = [];
      for (let e = 0; e < nEdges; ++e) {
        const curveType = c.u32();
        const nCurveParams = c.u32();
        edges.push({ curveType, params: c.doubles(nCurveParams) });
      }
      wires.push({ role, edges });
    }

    const edgeRefs = [];
    if (version >= 3) {
      const nEdgeRefs = c.u32();
      if (nEdgeRefs * 5 > c.remaining) {
        throw new Error(`${fileLabel}: surface ${s} claims ${nEdgeRefs} edge identities, more than the file holds`);
      }
      for (let e = 0; e < nEdgeRefs; ++e) {
        const id = c.u32();
        const eflags = c.u8();
        if (nModelEdges > 0 && id >= nModelEdges) {
          throw new Error(`${fileLabel}: surface ${s} edge identity ${e} is ${id}, outside the model's ${nModelEdges} edge(s)`);
        }
        edgeRefs.push({ id, flags: eflags });
      }
    }

    const expected = EXPECTED_PARAMS[type];
    if (expected === undefined) {
      // Unknown record kinds are carried, not dropped: the renderer counts and flags them.
      surfaces.push({ index: s, type, innerWall: false, params, wires, edgeRefs, unsupported: `unknown surface type ${type}` });
      continue;
    }
    if (nParams !== expected) {
      throw new Error(`${fileLabel}: ${surfaceTypeName(type)} surface ${s} has ${nParams} parameters, expected ${expected}`);
    }
    surfaces.push({
      index: s,
      type,
      innerWall: (flags & FLAG_INNER_WALL) !== 0,
      params,
      wires,
      edgeRefs,
      unsupported: null,
    });
  }

  return { version, nSurfaces, reserved, modelTolerance, modelToleranceStated, nModelEdges, surfaces, warnings, byteLength: buffer.byteLength };
}

/// The wire-join acceptance band for a model that declares its own tolerance (wireJoinToleranceFor).
export const K_WIRE_JOIN_TOLERANCE = 1e-6;
export function wireJoinToleranceFor(modelTolerance) {
  return Math.max(K_WIRE_JOIN_TOLERANCE, modelTolerance > 0 ? modelTolerance : 0);
}

/// Parse facets_*.bin: uint32 nTriangles, then 9 x float32 per triangle (cm).
/// Returns { nTriangles, positions: Float32Array(9*n) }.
export function parseFacets(buffer, fileLabel = 'facets') {
  if (buffer.byteLength < 4) { throw new Error(`${fileLabel}: truncated header`); }
  const view = new DataView(buffer);
  const n = view.getUint32(0, true);
  const need = 4 + n * 36;
  if (buffer.byteLength < need) {
    throw new Error(`${fileLabel}: truncated facet data (${n} triangles need ${need} bytes, file has ${buffer.byteLength})`);
  }
  // The payload is little-endian float32 and every platform this runs on is little-endian, so a
  // typed-array view over the buffer is the same bytes without a copy loop.
  const positions = new Float32Array(buffer.slice(4, need));
  return { nTriangles: n, positions };
}
