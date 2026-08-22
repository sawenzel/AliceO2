// The part selector: the one global control on this page.
//
// It is a listbox rather than a <select> because an <option> cannot hold a picture, and the whole
// point of this control is that you pick a part by looking at it. Every tab renders whatever this
// selector says; nothing else on the page chooses a part.

import { thumbnailFor } from './thumbs.js';

export class PartSelector {
  constructor(container, { onChange } = {}) {
    this.container = container;
    this.onChange = onChange || (() => {});
    this.parts = [];
    this.value = null;
    this.open = false;
    this._thumbsRequested = false;

    container.classList.add('partsel');
    container.innerHTML = `
      <button type="button" class="partsel-button" aria-haspopup="listbox" aria-expanded="false">
        <span class="partsel-thumb"></span>
        <span class="partsel-name">-</span>
        <span class="partsel-sub"></span>
        <span class="partsel-caret" aria-hidden="true">&#9662;</span>
      </button>
      <div class="partsel-panel" role="listbox" hidden></div>`;
    this.button = container.querySelector('.partsel-button');
    this.panel = container.querySelector('.partsel-panel');

    this.button.addEventListener('click', () => this.toggle());
    this.button.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); this.step(e.key === 'ArrowDown' ? 1 : -1); }
    });
    document.addEventListener('pointerdown', (e) => {
      if (this.open && !container.contains(e.target)) { this.setOpen(false); }
    });
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape' && this.open) { this.setOpen(false); this.button.focus(); } });
  }

  setParts(parts) {
    this.parts = parts;
    this._buildPanel();
  }

  /// Move the selection by one entry, so the whole page can be stepped through from the keyboard.
  step(delta) {
    if (!this.parts.length) { return; }
    const index = Math.max(0, this.parts.findIndex(p => p.name === this.value));
    const next = this.parts[Math.max(0, Math.min(this.parts.length - 1, index + delta))];
    if (next && next.name !== this.value) { this.select(next.name, true); }
  }

  select(name, notify = false) {
    const entry = this.parts.find(p => p.name === name);
    if (!entry) { return; }
    this.value = name;
    this.button.querySelector('.partsel-name').textContent = entry.name;
    this.button.querySelector('.partsel-sub').textContent = subtitleOf(entry);
    const slot = this.button.querySelector('.partsel-thumb');
    slot.innerHTML = '';
    thumbnailFor(entry).then((url) => {
      if (this.value !== name) { return; }
      slot.innerHTML = '';
      if (url) { const img = new Image(); img.src = url; img.alt = ''; slot.appendChild(img); }
    });
    for (const item of this.panel.querySelectorAll('.partsel-item')) {
      item.setAttribute('aria-selected', String(item.dataset.name === name));
    }
    if (notify) { this.onChange(entry); }
  }

  toggle() { this.setOpen(!this.open); }

  setOpen(open) {
    this.open = open;
    this.panel.hidden = !open;
    this.button.setAttribute('aria-expanded', String(open));
    // Thumbnails are rendered on the first open, not at boot: for a part list with ALICE3 in it
    // that is a few megabytes of facets nobody has asked to look at yet.
    if (open && !this._thumbsRequested) { this._thumbsRequested = true; this._fillThumbnails(); }
  }

  _buildPanel() {
    this.panel.innerHTML = '';
    let group = null;
    for (const entry of this.parts) {
      if ((entry.group || '') !== group) {
        group = entry.group || '';
        if (group) {
          const heading = document.createElement('div');
          heading.className = 'partsel-group';
          heading.textContent = group;
          this.panel.appendChild(heading);
        }
      }
      const item = document.createElement('button');
      item.type = 'button';
      item.className = 'partsel-item';
      item.setAttribute('role', 'option');
      item.dataset.name = entry.name;
      item.setAttribute('aria-selected', String(entry.name === this.value));
      const thumb = document.createElement('span');
      thumb.className = 'partsel-thumb';
      const text = document.createElement('span');
      text.className = 'partsel-text';
      const name = document.createElement('span');
      name.className = 'partsel-itemname';
      name.textContent = entry.name;
      const sub = document.createElement('span');
      sub.className = 'partsel-itemsub';
      sub.textContent = subtitleOf(entry);
      text.append(name, sub);
      item.append(thumb, text);
      item.addEventListener('click', () => { this.setOpen(false); this.select(entry.name, true); });
      this.panel.appendChild(item);
    }
  }

  async _fillThumbnails() {
    // Sequentially, so one shared offscreen renderer is enough and a long list does not stall the
    // page while it draws.
    for (const entry of this.parts) {
      const item = this.panel.querySelector(`.partsel-item[data-name="${cssEscape(entry.name)}"]`);
      if (!item) { continue; }
      const url = await thumbnailFor(entry);
      const slot = item.querySelector('.partsel-thumb');
      slot.innerHTML = '';
      if (url) { const img = new Image(); img.src = url; img.alt = ''; slot.appendChild(img); }
      else { slot.classList.add('empty'); }
    }
  }
}

function cssEscape(s) { return s.replace(/["\\]/g, '\\$&'); }

/// The one line under a part's name: what representations this checkout actually has for it.
function subtitleOf(entry) {
  const bits = [];
  if (entry.surfaces) { bits.push('exact'); }
  if (entry.facets) { bits.push('mesh'); }
  if (entry.shape) { bits.push('CSG'); }
  if (!bits.length) { return 'nothing to load'; }
  if (!entry.surfaces) { return entry.shape ? bits.join(' + ') : 'tessellated only'; }
  if (!entry.facets) { return `${bits.join(' + ')}, no mesh`; }
  return bits.join(' + ');
}
