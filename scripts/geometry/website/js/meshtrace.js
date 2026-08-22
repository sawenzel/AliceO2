// Ray tracing against the triangle mesh, so the exact solid and the tessellation can be rendered
// from the same camera and subtracted.
//
// Moller-Trumbore over a median-split BVH. The mesh answers the same two questions the exact solid
// does: the nearest crossing (for the picture) and the number of crossings along the whole ray
// (for the watertightness overlay). The second is the interesting one -- a tessellation that loses
// a ray shows up there and nowhere else.

const EPSILON = 1e-12;

export class MeshSolid {
  /// `positions` is the flat Float32Array of 9 floats per triangle that parseFacets returns.
  constructor(positions, label = 'mesh') {
    this.label = label;
    this.positions = positions;
    this.nTriangles = positions.length / 9;
    this.aabb = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
    for (let i = 0; i < positions.length; i += 3) {
      for (let c = 0; c < 3; ++c) {
        const v = positions[i + c];
        if (v < this.aabb[c]) { this.aabb[c] = v; }
        if (v > this.aabb[c + 3]) { this.aabb[c + 3] = v; }
      }
    }
    if (!Number.isFinite(this.aabb[0])) { this.aabb = [-1, -1, -1, 1, 1, 1]; }
    this._buildBVH();
  }

  _buildBVH() {
    const n = this.nTriangles;
    const centroids = new Float32Array(n * 3);
    const boxes = new Float32Array(n * 6);
    for (let t = 0; t < n; ++t) {
      const o = t * 9;
      for (let c = 0; c < 3; ++c) {
        const a = this.positions[o + c], b = this.positions[o + 3 + c], d = this.positions[o + 6 + c];
        boxes[t * 6 + c] = Math.min(a, b, d);
        boxes[t * 6 + 3 + c] = Math.max(a, b, d);
        centroids[t * 3 + c] = (a + b + d) / 3;
      }
    }
    const order = new Int32Array(n);
    for (let i = 0; i < n; ++i) { order[i] = i; }

    // nodes: [minx,miny,minz,maxx,maxy,maxz, leftOrFirst, countOrZero]
    const nodes = [];
    const build = (begin, end) => {
      const box = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
      for (let i = begin; i < end; ++i) {
        const t = order[i];
        for (let c = 0; c < 3; ++c) {
          box[c] = Math.min(box[c], boxes[t * 6 + c]);
          box[c + 3] = Math.max(box[c + 3], boxes[t * 6 + 3 + c]);
        }
      }
      const index = nodes.length;
      nodes.push({ box, left: -1, right: -1, first: begin, count: end - begin });
      if (end - begin <= 8) { return index; }
      let axis = 0, extent = -1;
      for (let c = 0; c < 3; ++c) { const e = box[c + 3] - box[c]; if (e > extent) { extent = e; axis = c; } }
      const mid = (begin + end) >> 1;
      const slice = Array.from(order.subarray(begin, end));
      slice.sort((a, b) => centroids[a * 3 + axis] - centroids[b * 3 + axis]);
      order.set(slice, begin);
      const left = build(begin, mid);
      const right = build(mid, end);
      nodes[index].left = left;
      nodes[index].right = right;
      nodes[index].count = 0;
      return index;
    };
    if (n > 0) { build(0, n); } else { nodes.push({ box: this.aabb.slice(), left: -1, right: -1, first: 0, count: 0 }); }
    this.nodes = nodes;
    this.order = order;
  }

  _triangleHit(t, ox, oy, oz, dx, dy, dz) {
    const o = t * 9;
    const p = this.positions;
    const ax = p[o], ay = p[o + 1], az = p[o + 2];
    const e1x = p[o + 3] - ax, e1y = p[o + 4] - ay, e1z = p[o + 5] - az;
    const e2x = p[o + 6] - ax, e2y = p[o + 7] - ay, e2z = p[o + 8] - az;
    const hx = dy * e2z - dz * e2y, hy = dz * e2x - dx * e2z, hz = dx * e2y - dy * e2x;
    const a = e1x * hx + e1y * hy + e1z * hz;
    if (a > -EPSILON && a < EPSILON) { return null; }
    const f = 1 / a;
    const sx = ox - ax, sy = oy - ay, sz = oz - az;
    const u = f * (sx * hx + sy * hy + sz * hz);
    if (u < 0 || u > 1) { return null; }
    const qx = sy * e1z - sz * e1y, qy = sz * e1x - sx * e1z, qz = sx * e1y - sy * e1x;
    const v = f * (dx * qx + dy * qy + dz * qz);
    if (v < 0 || u + v > 1) { return null; }
    const distance = f * (e2x * qx + e2y * qy + e2z * qz);
    // the geometric normal, unnormalised
    return { t: distance, nx: e1y * e2z - e1z * e2y, ny: e1z * e2x - e1x * e2z, nz: e1x * e2y - e1y * e2x };
  }

  _traverse(ox, oy, oz, dx, dy, dz, tMin, tMax, visit) {
    const invx = 1 / dx, invy = 1 / dy, invz = 1 / dz;
    const stack = [0];
    while (stack.length) {
      const node = this.nodes[stack.pop()];
      const box = node.box;
      let t0 = tMin, t1 = tMax;
      let a = (box[0] - ox) * invx, b = (box[3] - ox) * invx;
      if (a > b) { const s = a; a = b; b = s; }
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      a = (box[1] - oy) * invy; b = (box[4] - oy) * invy;
      if (a > b) { const s = a; a = b; b = s; }
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      a = (box[2] - oz) * invz; b = (box[5] - oz) * invz;
      if (a > b) { const s = a; a = b; b = s; }
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      if (!(t0 <= t1)) { continue; }
      if (node.count > 0) {
        for (let i = node.first; i < node.first + node.count; ++i) {
          const hit = this._triangleHit(this.order[i], ox, oy, oz, dx, dy, dz);
          if (hit && hit.t >= tMin && hit.t <= tMax) { visit(hit); }
        }
      } else {
        stack.push(node.left, node.right);
      }
    }
  }

  firstHit(ox, oy, oz, dx, dy, dz, tMin = 1e-9, tMax = 1e30) {
    let best = null;
    this._traverse(ox, oy, oz, dx, dy, dz, tMin, tMax, (hit) => {
      if (!best || hit.t < best.t) { best = hit; }
    });
    return best;
  }

  /// How many triangles the ray crosses. A closed mesh gives an even count from outside; an odd
  /// one means the ray got in and never came out -- the leak the X-ray benchmark counts.
  crossingCount(ox, oy, oz, dx, dy, dz, tMin = 1e-9, tMax = 1e30) {
    const distances = [];
    this._traverse(ox, oy, oz, dx, dy, dz, tMin, tMax, (hit) => distances.push(hit.t));
    if (distances.length < 2) { return distances.length; }
    // Cluster hits that land on a shared edge, as the exact side does, so a ray through a seam
    // between two triangles counts once rather than twice.
    distances.sort((a, b) => a - b);
    let count = 0;
    let i = 0;
    while (i < distances.length) {
      // anchored on the cluster's first member, never on its predecessor: chaining neighbour to
      // neighbour is transitive and would swallow a thin wall at large ray parameter
      const anchor = distances[i];
      let end = i + 1;
      while (end < distances.length &&
             Math.abs(distances[end] - anchor) <= 1e-7 * Math.max(1, Math.abs(anchor))) { ++end; }
      ++count;
      i = end;
    }
    return count;
  }
}
