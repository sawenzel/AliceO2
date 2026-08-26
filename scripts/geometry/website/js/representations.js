// The subjects one part can be traced and benchmarked as, and the page's whole vocabulary for them.
//
// A CAD part is converted through a cascade: native CSG first, then the flat halfspace DNF, then
// the exact trimmed surfaces, then the tessellation. The artefact the converter wrote for a part
// IS the cascade's answer for it. But every one of them is a real TGeoShape the O2 kernel
// navigates, so a part carrying more than one can be put on the same camera and the same sample
// points and the pictures and timings compared -- which is what this site is for.
//
// Two of the six are NOT products of the conversion, and keeping them in their own slots is the
// whole point rather than a detail. Calling three different .root files "CSG" is what makes a
// page unreadable:
//
//   original   what the part was made FROM -- the shape the detector geometry ships today, before
//              it ever went to STEP. The only baseline that answers "is any of this better than
//              what we already have". Exists only for a round-tripped part.
//   cellstree  the decomposition emitted as a plain composite. NOTHING ships this. It exists so
//              the flat solid has something to be measured against: the cells the cascade refused
//              to emit as a tree because they would make one wider than the routing threshold.
//
// `role` is what keeps them out of every statement about what a part ships:
//   'product'    the converter emits it, and it can be what the part ships
//   'baseline'   the input to the conversion, not an output
//   'comparison' built for this page, shipped by nothing
//
// The array is in the order the page prints them: the baseline first, because that is what a
// reader starts from, then the cascade's own order.

export const REPRESENTATIONS = [
  {
    key: 'original',
    role: 'baseline',
    badge: 'ORIGINAL',
    field: 'original',
    label: 'Original TGeo',
    klass: 'the source TGeoShape',
    note: 'the shape this part was made from, before the round trip -- what the detector geometry ships today',
  },
  {
    key: 'shape',
    role: 'product',
    badge: 'CSG',
    field: 'shape',
    label: 'CSG',
    klass: 'TGeoCompositeShape / TGeo primitive',
    note: 'native CSG primitives, navigated by a TGeo composite shape',
  },
  {
    key: 'cellstree',
    role: 'comparison',
    badge: 'CELLS-TREE',
    field: 'cellstree',
    label: 'CSG cells-tree',
    klass: 'TGeoCompositeShape',
    note: 'the decomposed cells emitted as a plain TGeo composite -- a comparison for the flat ' +
          'solid, not something the converter ships',
  },
  {
    key: 'flatcsg',
    role: 'product',
    badge: 'FLATCSG',
    field: 'flatcsg',
    label: 'FlatCSG',
    klass: 'o2::base::O2FlatCSG',
    note: 'a union of intersection-cells over signed implicit halfspaces, with a BVH over ' +
          'sub-boxes of those cells, navigated by O2FlatCSG',
  },
  {
    key: 'surface',
    role: 'product',
    badge: 'SURFACE',
    field: 'surfaces',
    label: 'Surface',
    klass: 'o2::base::O2BVHSurfaceSolid',
    note: 'the exact trimmed analytic faces, navigated by O2BVHSurfaceSolid',
  },
  {
    key: 'mesh',
    role: 'product',
    badge: 'TESSELLATED',
    field: 'facets',
    label: 'Tessellated',
    klass: 'o2::base::O2Tessellated',
    note: 'the triangle mesh, navigated by O2Tessellated -- the fallback',
  },
];

export const REP_BY_KEY = new Map(REPRESENTATIONS.map(rep => [rep.key, rep]));

/// The cascade's own order, for anything that reports what a part ships. The baseline and the
/// comparison are not in it, because neither is something a part can ship.
export const PRODUCTS = REPRESENTATIONS.filter(rep => rep.role === 'product');

export function repBadge(key) { const r = REP_BY_KEY.get(key); return r ? r.badge : String(key).toUpperCase(); }
export function repField(key) { const r = REP_BY_KEY.get(key); return r ? r.field : null; }
export function repLabel(key) { const r = REP_BY_KEY.get(key); return r ? r.label : String(key); }
export function repNote(key) { const r = REP_BY_KEY.get(key); return r ? r.note : ''; }
export function repClass(key) { const r = REP_BY_KEY.get(key); return r ? r.klass : ''; }
export function repRole(key) { const r = REP_BY_KEY.get(key); return r ? r.role : 'product'; }

/// Does this manifest entry carry the artefact for this subject?
export function partHas(entry, key) {
  const field = repField(key);
  return !!(entry && field && entry[field]);
}

/// The site-relative path of one subject's artefact, or null.
export function partPath(entry, key) {
  return partHas(entry, key) ? `testdata/${entry[repField(key)]}` : null;
}

/// Every subject this part carries an artefact for, in print order.
export function partSubjects(entry) {
  return REPRESENTATIONS.filter(rep => partHas(entry, rep.key)).map(rep => rep.key);
}

/// What this part's copy of a subject actually is, where `fetch_testdata.sh` recorded it -- the
/// source volume and its boolean depth for the original, the cell and leaf counts for the
/// cells-tree. Returns null when there is nothing recorded beyond the generic note.
export function partNote(entry, key) {
  return (entry && entry.notes && entry.notes[key]) || null;
}
