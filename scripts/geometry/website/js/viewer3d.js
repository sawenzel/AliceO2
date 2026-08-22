// The shared three.js view: the tessellation, the exact trim-edge overlay, and whatever the event
// player adds on top of them.

import * as THREE from '../vendor/three.module.min.js';
import { OrbitControl } from './orbit.js';

export class Viewer3D {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.container.appendChild(this.renderer.domElement);
    this.renderer.domElement.style.display = 'block';
    this.renderer.domElement.style.width = '100%';
    this.renderer.domElement.style.height = '100%';
    this.renderer.domElement.style.touchAction = 'none';

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x0e1116);
    this.camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
    this.camera.position.set(3, 3, 3);

    this.controls = new OrbitControl(this.camera, this.renderer.domElement);
    this.controls.onChange = () => this.requestRender();

    const key = new THREE.DirectionalLight(0xffffff, 2.1);
    key.position.set(1, 1, 1);
    this.camera.add(key);
    const fill = new THREE.DirectionalLight(0x88aaff, 0.55);
    fill.position.set(-1, -0.4, -0.6);
    this.camera.add(fill);
    this.scene.add(this.camera);
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.28));

    this.meshGroup = new THREE.Group();
    this.edgeGroup = new THREE.Group();
    this.overlayGroup = new THREE.Group();
    this.gridGroup = new THREE.Group();
    this.scene.add(this.meshGroup, this.edgeGroup, this.overlayGroup, this.gridGroup);

    // A lost context comes back blank unless something asks for a frame; a dirty-flag render loop
    // otherwise sits idle and shows black.
    this.renderer.domElement.addEventListener('webglcontextlost', (e) => e.preventDefault());
    this.renderer.domElement.addEventListener('webglcontextrestored', () => this.requestRender());

    this._dirty = true;
    this._running = false;
    this._observer = new ResizeObserver(() => this.resize());
    this._observer.observe(this.container);
    this.resize();
    this.start();
  }

  requestRender() { this._dirty = true; }

  start() {
    if (this._running) { return; }
    this._running = true;
    const tick = () => {
      if (!this._running) { return; }
      if (this._dirty) { this._dirty = false; this.renderer.render(this.scene, this.camera); }
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  }

  resize() {
    const w = Math.max(1, this.container.clientWidth), h = Math.max(1, this.container.clientHeight);
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.requestRender();
  }

  _clear(group) {
    while (group.children.length) {
      const child = group.children.pop();
      if (child.geometry) { child.geometry.dispose(); }
      if (child.material) { child.material.dispose(); }
    }
  }

  clearMesh() { this._clear(this.meshGroup); this.requestRender(); }
  clearEdges() { this._clear(this.edgeGroup); this.requestRender(); }
  clearOverlay() { this._clear(this.overlayGroup); this.requestRender(); }

  /// Set the tessellation from a flat Float32Array of 9 floats per triangle. Flat shading, so the
  /// faceting is visible -- that is the point of this view, not a defect of it.
  setMesh(positions, { color = 0x9fb4cf, opacity = 1 } = {}) {
    this.clearMesh();
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.computeVertexNormals();
    const material = new THREE.MeshPhongMaterial({
      color, flatShading: true, side: THREE.DoubleSide, shininess: 18,
      transparent: opacity < 1, opacity,
    });
    this.solidMesh = new THREE.Mesh(geometry, material);
    this.meshGroup.add(this.solidMesh);

    const wireMaterial = new THREE.MeshBasicMaterial({ color: 0x2b3a4d, wireframe: true });
    this.wireMesh = new THREE.Mesh(geometry, wireMaterial);
    this.wireMesh.visible = false;
    this.meshGroup.add(this.wireMesh);
    this.requestRender();
    return geometry;
  }

  setWireframe(on) {
    if (this.wireMesh) { this.wireMesh.visible = on; }
    if (this.solidMesh) { this.solidMesh.material.opacity = on ? 0.55 : 1; this.solidMesh.material.transparent = on; }
    this.requestRender();
  }

  setMeshVisible(on) { this.meshGroup.visible = on; this.requestRender(); }
  setEdgesVisible(on) { this.edgeGroup.visible = on; this.requestRender(); }

  /// Draw the exact trimmed face boundaries as line loops: where the mesh is faceted, these run
  /// through the true curve.
  setExactEdges(polylines, color = 0xffc857) {
    this.clearEdges();
    const material = new THREE.LineBasicMaterial({ color });
    for (const line of polylines) {
      if (line.length < 2) { continue; }
      const positions = new Float32Array(line.length * 3);
      for (let i = 0; i < line.length; ++i) {
        positions[3 * i] = line[i][0]; positions[3 * i + 1] = line[i][1]; positions[3 * i + 2] = line[i][2];
      }
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      this.edgeGroup.add(new THREE.Line(geometry, material.clone()));
    }
    material.dispose();
    this.requestRender();
  }

  /// A ground grid sized to the model, so the eye has a scale reference (cm).
  setGrid(box) {
    this._clear(this.gridGroup);
    const size = Math.max(box[3] - box[0], box[4] - box[1], box[5] - box[2]) * 2 || 2;
    const grid = new THREE.GridHelper(size, 20, 0x394b60, 0x1d2733);
    grid.position.set((box[0] + box[3]) / 2, box[1] - 0.05 * size, (box[2] + box[5]) / 2);
    this.gridGroup.add(grid);
    this.requestRender();
  }

  frame(box) { this.controls.frameBox(box); this.requestRender(); }

  /// The camera as the raytracer needs it: origin, orthonormal basis and vertical field of view.
  cameraSpec() {
    const position = this.camera.position;
    const forward = new THREE.Vector3();
    this.camera.getWorldDirection(forward);
    const up = new THREE.Vector3(0, 1, 0).applyQuaternion(this.camera.quaternion);
    const right = new THREE.Vector3().crossVectors(forward, up).normalize();
    return {
      origin: [position.x, position.y, position.z],
      forward: [forward.x, forward.y, forward.z],
      up: [up.x, up.y, up.z],
      right: [right.x, right.y, right.z],
      fovY: this.camera.fov * Math.PI / 180,
    };
  }

  /// Point the three.js camera at exactly the spec a raytraced frame used, so switching tabs does
  /// not move the model.
  applyCameraSpec(spec, target) {
    this.camera.position.set(spec.origin[0], spec.origin[1], spec.origin[2]);
    if (target) { this.controls.target.set(target[0], target[1], target[2]); }
    this.camera.lookAt(this.controls.target);
    this.camera.fov = spec.fovY * 180 / Math.PI;
    this.camera.updateProjectionMatrix();
    this.requestRender();
  }
}
