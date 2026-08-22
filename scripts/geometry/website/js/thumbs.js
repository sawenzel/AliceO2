// Tiny 3D thumbnails for the part selector.
//
// One offscreen WebGL renderer, reused for every part, rendering the tessellation (or, for a part
// that has no mesh, the exact trim-edge loops) once into a data URL that is then cached. The
// selector shows a picture of the thing you are choosing rather than a name you have to recognise.

import * as THREE from '../vendor/three.module.min.js';
import { loadBinary } from './data.js';
import { parseFacets, parseSidecar } from './sidecar.js';
import { SurfaceSolid } from './solid.js';

const SIZE = 112;

const cache = new Map();    // part name -> data URL (or null when nothing could be drawn)
const inFlight = new Map(); // part name -> Promise

let renderer = null, scene = null, camera = null, content = null;

function ensureRenderer() {
  if (renderer) { return true; }
  const canvas = document.createElement('canvas');
  canvas.width = SIZE; canvas.height = SIZE;
  try {
    renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true, preserveDrawingBuffer: true });
  } catch (e) {
    renderer = null;
    return false;   // no WebGL here: the selector simply shows no pictures
  }
  renderer.setPixelRatio(1);
  renderer.setSize(SIZE, SIZE, false);
  renderer.setClearColor(0x151a21, 1);
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(38, 1, 0.01, 1000);
  const key = new THREE.DirectionalLight(0xffffff, 2.2); key.position.set(1, 1, 1);
  camera.add(key);
  const fill = new THREE.DirectionalLight(0x88aaff, 0.6); fill.position.set(-1, -0.4, -0.6);
  camera.add(fill);
  scene.add(camera);
  scene.add(new THREE.AmbientLight(0xffffff, 0.3));
  content = new THREE.Group();
  scene.add(content);
  return true;
}

function clearContent() {
  while (content.children.length) {
    const child = content.children.pop();
    if (child.geometry) { child.geometry.dispose(); }
    if (child.material) { child.material.dispose(); }
  }
}

function frame(box) {
  const centre = new THREE.Vector3((box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2);
  const radius = 0.5 * Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
  const distance = radius / Math.sin((camera.fov * Math.PI / 180) / 2) * 1.05;
  const direction = new THREE.Vector3(1, 0.65, 1).normalize().multiplyScalar(distance);
  camera.position.copy(centre).add(direction);
  camera.near = Math.max(1e-4, distance - 4 * radius);
  camera.far = distance + 8 * radius;
  camera.updateProjectionMatrix();
  camera.lookAt(centre);
}

function boxOf(positions) {
  const box = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += 3) {
    for (let k = 0; k < 3; ++k) {
      const v = positions[i + k];
      if (v < box[k]) { box[k] = v; }
      if (v > box[k + 3]) { box[k + 3] = v; }
    }
  }
  return Number.isFinite(box[0]) ? box : [-1, -1, -1, 1, 1, 1];
}

async function renderThumbnail(entry) {
  if (!ensureRenderer()) { return null; }
  clearContent();
  let box = null;

  if (entry.facets) {
    try {
      const facets = parseFacets(await loadBinary(`testdata/${entry.facets}`), entry.facets);
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(facets.positions, 3));
      geometry.computeVertexNormals();
      content.add(new THREE.Mesh(geometry, new THREE.MeshPhongMaterial({
        color: entry.surfaces ? 0x9fb4cf : 0xb08a6a, flatShading: true, side: THREE.DoubleSide, shininess: 20,
      })));
      box = boxOf(facets.positions);
    } catch (e) { /* fall through to the edge-only thumbnail */ }
  }

  if (!box && entry.surfaces) {
    try {
      const solid = new SurfaceSolid(parseSidecar(await loadBinary(`testdata/${entry.surfaces}`), entry.name), entry.name);
      const material = new THREE.LineBasicMaterial({ color: 0xffc857 });
      for (const surface of solid.surfaces) {
        let loops = [];
        try { loops = surface.patchOutline(); } catch (e) { continue; }
        for (const loop of loops) {
          if (loop.length < 2) { continue; }
          const positions = new Float32Array(loop.length * 3);
          for (let i = 0; i < loop.length; ++i) {
            positions[3 * i] = loop[i][0]; positions[3 * i + 1] = loop[i][1]; positions[3 * i + 2] = loop[i][2];
          }
          const geometry = new THREE.BufferGeometry();
          geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
          content.add(new THREE.Line(geometry, material.clone()));
        }
      }
      material.dispose();
      box = solid.aabb;
    } catch (e) { /* nothing drawable */ }
  }

  if (!box) { return null; }
  frame(box);
  renderer.render(scene, camera);
  const url = renderer.domElement.toDataURL('image/png');
  clearContent();
  return url;
}

/// The thumbnail for a manifest entry, rendered once and then served from the cache. Returns null
/// when nothing could be drawn (no WebGL, or neither a mesh nor a drawable sidecar).
export function thumbnailFor(entry) {
  if (cache.has(entry.name)) { return Promise.resolve(cache.get(entry.name)); }
  if (inFlight.has(entry.name)) { return inFlight.get(entry.name); }
  const promise = renderThumbnail(entry)
    .catch(() => null)
    .then((url) => { cache.set(entry.name, url); inFlight.delete(entry.name); return url; });
  inFlight.set(entry.name, promise);
  return promise;
}
