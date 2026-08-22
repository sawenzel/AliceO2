// A minimal orbit control: drag to rotate, wheel to dolly, right-drag or shift-drag to pan.
// Written here rather than vendored so the page needs no bare-specifier import map and the
// single-file build stays a concatenation.

import * as THREE from '../vendor/three.module.min.js';

export class OrbitControl {
  constructor(camera, domElement, target = new THREE.Vector3()) {
    this.camera = camera;
    this.dom = domElement;
    this.target = target.clone();
    this.minDistance = 1e-4;
    this.maxDistance = 1e7;
    this.rotateSpeed = 0.006;
    this.zoomSpeed = 0.0015;
    this.enabled = true;
    this.onChange = null;

    this._dragging = 0; // 0 none, 1 rotate, 2 pan
    this._last = { x: 0, y: 0 };

    this._onPointerDown = (e) => {
      if (!this.enabled) { return; }
      this._dragging = (e.button === 2 || e.shiftKey) ? 2 : 1;
      this._last = { x: e.clientX, y: e.clientY };
      this.dom.setPointerCapture(e.pointerId);
      e.preventDefault();
    };
    this._onPointerMove = (e) => {
      if (!this._dragging) { return; }
      const dx = e.clientX - this._last.x, dy = e.clientY - this._last.y;
      this._last = { x: e.clientX, y: e.clientY };
      if (this._dragging === 1) { this.rotate(dx, dy); } else { this.pan(dx, dy); }
    };
    this._onPointerUp = (e) => {
      this._dragging = 0;
      try { this.dom.releasePointerCapture(e.pointerId); } catch (_) { /* already released */ }
    };
    this._onWheel = (e) => {
      if (!this.enabled) { return; }
      e.preventDefault();
      this.dolly(Math.exp(e.deltaY * this.zoomSpeed));
    };
    this._onContext = (e) => e.preventDefault();

    this.dom.addEventListener('pointerdown', this._onPointerDown);
    this.dom.addEventListener('pointermove', this._onPointerMove);
    this.dom.addEventListener('pointerup', this._onPointerUp);
    this.dom.addEventListener('pointercancel', this._onPointerUp);
    this.dom.addEventListener('wheel', this._onWheel, { passive: false });
    this.dom.addEventListener('contextmenu', this._onContext);
  }

  dispose() {
    this.dom.removeEventListener('pointerdown', this._onPointerDown);
    this.dom.removeEventListener('pointermove', this._onPointerMove);
    this.dom.removeEventListener('pointerup', this._onPointerUp);
    this.dom.removeEventListener('pointercancel', this._onPointerUp);
    this.dom.removeEventListener('wheel', this._onWheel);
    this.dom.removeEventListener('contextmenu', this._onContext);
  }

  _changed() { if (this.onChange) { this.onChange(); } }

  rotate(dx, dy) {
    const offset = this.camera.position.clone().sub(this.target);
    const radius = offset.length();
    let theta = Math.atan2(offset.x, offset.z);
    let phi = Math.acos(Math.max(-1, Math.min(1, offset.y / radius)));
    theta -= dx * this.rotateSpeed;
    phi = Math.max(1e-4, Math.min(Math.PI - 1e-4, phi - dy * this.rotateSpeed));
    offset.set(radius * Math.sin(phi) * Math.sin(theta), radius * Math.cos(phi), radius * Math.sin(phi) * Math.cos(theta));
    this.camera.position.copy(this.target).add(offset);
    this.camera.lookAt(this.target);
    this._changed();
  }

  pan(dx, dy) {
    const offset = this.camera.position.clone().sub(this.target);
    const distance = offset.length();
    const fovScale = 2 * Math.tan((this.camera.fov * Math.PI / 180) / 2) * distance;
    const height = this.dom.clientHeight || 1;
    const right = new THREE.Vector3().setFromMatrixColumn(this.camera.matrix, 0);
    const up = new THREE.Vector3().setFromMatrixColumn(this.camera.matrix, 1);
    const move = right.multiplyScalar(-dx * fovScale / height).add(up.multiplyScalar(dy * fovScale / height));
    this.camera.position.add(move);
    this.target.add(move);
    this.camera.lookAt(this.target);
    this._changed();
  }

  dolly(factor) {
    const offset = this.camera.position.clone().sub(this.target);
    const distance = Math.max(this.minDistance, Math.min(this.maxDistance, offset.length() * factor));
    offset.setLength(distance);
    this.camera.position.copy(this.target).add(offset);
    this.camera.lookAt(this.target);
    this._changed();
  }

  /// Place the camera so that a bounding box [xmin,ymin,zmin,xmax,ymax,zmax] fills the view.
  frameBox(box, direction = new THREE.Vector3(1, 0.7, 1)) {
    const center = new THREE.Vector3((box[0] + box[3]) / 2, (box[1] + box[4]) / 2, (box[2] + box[5]) / 2);
    const radius = 0.5 * Math.hypot(box[3] - box[0], box[4] - box[1], box[5] - box[2]) || 1;
    const distance = radius / Math.sin((this.camera.fov * Math.PI / 180) / 2) * 1.1;
    this.target.copy(center);
    this.camera.position.copy(center).add(direction.clone().normalize().multiplyScalar(distance));
    this.camera.near = Math.max(1e-4, distance - 4 * radius);
    this.camera.far = distance + 8 * radius;
    this.camera.updateProjectionMatrix();
    this.camera.lookAt(center);
    this._changed();
  }
}
