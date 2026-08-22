// The Assembly view: the placed tree, not one part.
//
// oTOF has no tubular *solid* -- all twenty of its body prototypes are planar, nineteen of them
// exact boxes. The tube is the *assembly*: 62 628 placed bodies that together are a barrel. So
// this view draws the placement table rather than a sidecar, one three.js InstancedMesh per body
// prototype, and the tube is what comes out.
//
// The geometry fact that makes it one draw call per prototype: a body's vertices are already in
// its LEAF's frame -- the converter bakes the body's own pose in -- so a body instance's world
// transform is exactly its leaf's matrix, and every prototype of a leaf shares one instance-matrix
// buffer.
//
// This is WebGL only. There is no raytracing and no bridge here; 62 628 solids is a rasteriser's
// job, and saying so on the page is part of the point.

import * as THREE from '../vendor/three.module.min.js';
import { loadBinary, loadJSON } from './data.js';
import { parseFacets } from './sidecar.js';
import { Viewer3D } from './viewer3d.js';

const AXIS_NAME = ['x', 'y', 'z'];

/// Digits grouped with a non-breaking space: 62628 -> "62 628".
export function grouped(n) { return String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, ' '); }

// --- the data ------------------------------------------------------------------------------------

const cache = new Map(); // root -> Promise of the built model

/// Load an assembly directory: index.json, its placement table, and every body prototype's mesh.
/// Memoised, because the selector's thumbnail and the view both want it.
export function loadAssemblyData(root = 'testdata/otof_assembly') {
  if (!cache.has(root)) { cache.set(root, build(root).catch((e) => { cache.delete(root); throw e; })); }
  return cache.get(root);
}

/// Read the index alone, to decide whether this checkout has an assembly at all.
export async function loadAssemblyIndex(root = 'testdata/otof_assembly') {
  const index = await loadJSON(`${root}/index.json`, { optional: true });
  return index && Array.isArray(index.leaves) ? index : null;
}

/// A 3x4 row-major [rotation | translation] as three.js wants it: 16 column-major floats.
function toColumnMajor(m, out, offset) {
  out[offset + 0] = m[0][0]; out[offset + 1] = m[1][0]; out[offset + 2] = m[2][0]; out[offset + 3] = 0;
  out[offset + 4] = m[0][1]; out[offset + 5] = m[1][1]; out[offset + 6] = m[2][1]; out[offset + 7] = 0;
  out[offset + 8] = m[0][2]; out[offset + 9] = m[1][2]; out[offset + 10] = m[2][2]; out[offset + 11] = 0;
  out[offset + 12] = m[0][3]; out[offset + 13] = m[1][3]; out[offset + 14] = m[2][3]; out[offset + 15] = 1;
}

function localBoxCorners(positions) {
  const box = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += 3) {
    for (let k = 0; k < 3; ++k) {
      const v = positions[i + k];
      if (v < box[k]) { box[k] = v; }
      if (v > box[k + 3]) { box[k + 3] = v; }
    }
  }
  if (!Number.isFinite(box[0])) { return []; }
  const corners = [];
  for (const x of [box[0], box[3]]) {
    for (const y of [box[1], box[4]]) {
      for (const z of [box[2], box[5]]) { corners.push(x, y, z); }
    }
  }
  return corners;
}

async function build(root) {
  const index = await loadAssemblyIndex(root);
  if (!index) { throw new Error(`${root}/index.json is missing -- run ./fetch_assembly.sh`); }
  const table = await loadJSON(`${root}/${index.placements || 'placements.json'}`);
  const byEntry = new Map(table.leaves.map(leaf => [leaf.entry, leaf]));

  const leaves = [];
  for (const leafIndex of index.leaves) {
    const placed = byEntry.get(leafIndex.entry);
    if (!placed || !leafIndex.bodies.length) { continue; }  // an empty node places nothing
    const bodies = [];
    for (const body of leafIndex.bodies) {
      const facets = parseFacets(await loadBinary(`${root}/${body.facets}`), body.facets);
      bodies.push({ name: body.name, facets, corners: localBoxCorners(facets.positions) });
    }
    leaves.push({ entry: leafIndex.entry, name: leafIndex.name, bodies, matrices: placed.matrices });
  }
  if (!leaves.length) { throw new Error(`${root}: no leaf carries a body prototype`); }

  // The extent of the placed assembly, from every body's world-transformed local box corners.
  const box = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  const sweep = (visit) => {
    for (const leaf of leaves) {
      for (const m of leaf.matrices) {
        for (const body of leaf.bodies) {
          const c = body.corners;
          for (let i = 0; i < c.length; i += 3) {
            const x = c[i], y = c[i + 1], z = c[i + 2];
            visit(m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3],
                  m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3],
                  m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3]);
          }
        }
      }
    }
  };
  sweep((x, y, z) => {
    const p = [x, y, z];
    for (let k = 0; k < 3; ++k) {
      if (p[k] < box[k]) { box[k] = p[k]; }
      if (p[k] > box[k + 3]) { box[k + 3] = p[k]; }
    }
  });

  // The barrel axis is the long one, and its centre line is the middle of the two short spans --
  // measured, not assumed, so an assembly laid out along another axis still reads correctly.
  const span = [box[3] - box[0], box[4] - box[1], box[5] - box[2]];
  const axis = span.indexOf(Math.max(...span));
  const transverse = [0, 1, 2].filter(k => k !== axis);
  const centre = transverse.map(k => (box[k] + box[k + 3]) / 2);
  let rMin = Infinity, rMax = 0;
  sweep((x, y, z) => {
    const p = [x, y, z];
    const r = Math.hypot(p[transverse[0]] - centre[0], p[transverse[1]] - centre[1]);
    if (r < rMin) { rMin = r; }
    if (r > rMax) { rMax = r; }
  });

  // Every leaf's placements sorted along the axis, as one shared instance-matrix buffer: the slice
  // then costs a copy of a contiguous run rather than a rebuild.
  for (const leaf of leaves) {
    const n = leaf.matrices.length;
    const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => leaf.matrices[a][axis][3] - leaf.matrices[b][axis][3]);
    const matrices = new Float32Array(16 * n);
    const along = new Float64Array(n);
    for (let i = 0; i < n; ++i) {
      toColumnMajor(leaf.matrices[order[i]], matrices, 16 * i);
      along[i] = leaf.matrices[order[i]][axis][3];
    }
    leaf.sorted = matrices;
    leaf.along = along;
    leaf.count = n;
  }

  // The stations of the barrel: where along the axis the leaf that places the most bodies puts
  // them. One station is one ring, and the spacing between two of them is the finest slice worth
  // offering. Positions are clustered first, because two placements a fraction of a millimetre
  // apart are the same ring and not two of them.
  const dominant = leaves.reduce((a, b) => (a.count * a.bodies.length >= b.count * b.bodies.length ? a : b));
  const medianGap = (values) => {
    const gaps = [];
    for (let i = 1; i < values.length; ++i) { if (values[i] - values[i - 1] > 1e-6) { gaps.push(values[i] - values[i - 1]); } }
    if (!gaps.length) { return 0; }
    gaps.sort((a, b) => a - b);
    return gaps[gaps.length >> 1];
  };
  const raw = Array.from(new Set(Array.from(dominant.along))).sort((a, b) => a - b);
  const tolerance = medianGap(raw) / 4;
  const stations = [];
  for (const value of raw) {
    if (!stations.length || value - stations[stations.length - 1] > tolerance) { stations.push(value); }
  }
  let step = medianGap(stations);
  if (!(step > 0)) { step = span[axis] || 1; }

  return {
    root, index, leaves, box, axis, span, step, stations,
    rMin, rMax, centre, transverse,
    totals: index.totals || {},
  };
}

// --- colours ---------------------------------------------------------------------------------

// Muted, and one family per leaf: the modules dominate the picture, so they are one quiet slate
// range and the support frame is a warmer one that reads against it without shouting.
const FAMILY = [
  { h: 0.105, s: 0.30 },  // amber-olive: the oTOF support bodies
  { h: 0.575, s: 0.20 },  // slate blue: the modules
];

function bodyColour(leafIndex, bodyIndex, bodyCount) {
  const family = FAMILY[leafIndex % FAMILY.length];
  const t = bodyCount > 1 ? bodyIndex / (bodyCount - 1) : 0.5;
  return new THREE.Color().setHSL(family.h, family.s, 0.42 + 0.24 * t);
}

// --- the view --------------------------------------------------------------------------------

export class AssemblyView {
  constructor(container, model) {
    this.model = model;
    this.viewer = new Viewer3D(container);
    this.viewer.setGrid(model.box);
    this.meshes = [];   // {mesh, leaf}
    this.legend = [];   // {name, colour, triangles, instances}
    this.spinning = false;
    this.fps = 0;
    this.visibleSolids = 0;
    this._buildInstances();
    // Down the axis and a little off it, so the bore and the length are both in the first frame.
    this.frame();
  }

  _buildInstances() {
    const { model } = this;
    for (let li = 0; li < model.leaves.length; ++li) {
      const leaf = model.leaves[li];
      // One buffer per leaf, shared by every prototype it places: same transforms, same slice.
      leaf.attribute = new THREE.InstancedBufferAttribute(leaf.sorted.slice(), 16);
      leaf.attribute.setUsage(THREE.DynamicDrawUsage);
      for (let bi = 0; bi < leaf.bodies.length; ++bi) {
        const body = leaf.bodies[bi];
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(body.facets.positions, 3));
        geometry.computeVertexNormals();
        const colour = bodyColour(li, bi, leaf.bodies.length);
        const material = new THREE.MeshPhongMaterial({
          color: colour, flatShading: true, side: THREE.DoubleSide, shininess: 12,
        });
        const mesh = new THREE.InstancedMesh(geometry, material, leaf.count);
        mesh.instanceMatrix = leaf.attribute;
        mesh.count = leaf.count;
        mesh.frustumCulled = false;   // twenty draw calls; culling them is not worth a bounds pass
        this.viewer.overlayGroup.add(mesh);
        this.meshes.push({ mesh, leaf });
        this.legend.push({
          name: body.name, colour: `#${colour.getHexString()}`,
          triangles: body.facets.nTriangles, instances: leaf.count,
        });
      }
    }
    this.visibleSolids = this.meshes.reduce((sum, m) => sum + m.mesh.count, 0);
  }

  /// Keep only the placements whose position along the barrel axis is inside [lo, hi].
  setSlice(lo, hi) {
    for (const leaf of this.model.leaves) {
      const along = leaf.along;
      let first = 0, last = along.length;
      while (first < last && along[first] < lo) { ++first; }
      while (last > first && along[last - 1] > hi) { --last; }
      const n = last - first;
      if (n > 0) { leaf.attribute.array.set(leaf.sorted.subarray(16 * first, 16 * last), 0); }
      leaf.attribute.needsUpdate = true;
      leaf.visible = n;
    }
    for (const { mesh, leaf } of this.meshes) { mesh.count = leaf.visible; }
    this.visibleSolids = this.meshes.reduce((sum, m) => sum + m.mesh.count, 0);
    this.viewer.requestRender();
  }

  setWireframe(on) {
    for (const { mesh } of this.meshes) { mesh.material.wireframe = on; }
    this.viewer.requestRender();
  }

  /// Frame the assembly, or -- when the slice has narrowed it -- just the part still on screen.
  frame(box = this.model.box) {
    const direction = new THREE.Vector3(0, 0, 0);
    direction.setComponent(this.model.axis, 1);
    direction.setComponent(this.model.transverse[0], 0.42);
    direction.setComponent(this.model.transverse[1], 0.62);
    this.viewer.controls.frameBox(box, direction);
    this.viewer.requestRender();
  }

  /// Turn the model at a steady rate, and report the frame rate that comes back. The turn is
  /// driven from the renderer's own frame hook, not from an animation-frame callback: on a slow
  /// rasteriser a callback can fire many times per drawn frame, and a rate counted from those
  /// would flatter the view rather than measure it.
  setSpin(on, onFrame) {
    this.spinning = on;
    if (!on) { this.viewer.onFrame = null; this.fps = 0; return; }
    let previous = performance.now();
    this.viewer.onFrame = (now) => {
      if (!this.spinning) { return; }
      const dt = now - previous;
      previous = now;
      if (dt > 0 && dt < 4000) { this.fps = this.fps ? this.fps * 0.85 + (1000 / dt) * 0.15 : 1000 / dt; }
      this.viewer.controls.rotate(Math.min(6, dt * 0.09), 0);   // a steady turn, whatever the rate
      if (onFrame) { onFrame(this.fps); }
    };
    this.viewer.requestRender();
  }
}

/// The line that says what is on screen: prototypes, leaf placements, solids.
export function statsLine(model) {
  const t = model.totals;
  const prototypes = t.prototypes != null ? t.prototypes : model.leaves.reduce((n, l) => n + l.bodies.length, 0);
  const placements = t.placements != null ? t.placements : model.leaves.reduce((n, l) => n + l.count, 0);
  const solids = t.solids != null ? t.solids : model.leaves.reduce((n, l) => n + l.count * l.bodies.length, 0);
  return `${grouped(prototypes)} prototypes · ${grouped(placements)} leaf placements · ${grouped(solids)} solids`;
}

/// What the barrel measures, taken from the placed geometry rather than from a caption.
export function barrelLine(model) {
  const axis = AXIS_NAME[model.axis];
  return `barrel about ${axis} · R ${model.rMin.toFixed(1)}–${model.rMax.toFixed(1)} cm · ` +
         `${(model.span[model.axis] / 100).toFixed(2)} m long`;
}

export { AXIS_NAME };
