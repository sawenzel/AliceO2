// Run the self-check from the command line against testdata/:
//
//   node tools/selfcheck.mjs [part ...]
//
// It is the same code the Self-check tab runs; only the loader differs. Exits non-zero on failure.

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { runSelfCheck } from '../js/selfcheck.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, '..');
const testdata = path.join(root, 'testdata');

function readBuffer(file) {
  const b = fs.readFileSync(file);
  return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
}

let parts = process.argv.slice(2);
if (!parts.length) {
  const manifestPath = path.join(testdata, 'manifest.json');
  if (!fs.existsSync(manifestPath)) {
    console.error(`no ${manifestPath}; run ./fetch_testdata.sh <gate-workdir> first`);
    process.exit(2);
  }
  // A tessellated-only part has no sidecar, and every assertion here is about the exact solid.
  parts = JSON.parse(fs.readFileSync(manifestPath, 'utf8')).parts.filter(p => p.surfaces).map(p => p.name);
}

const load = async (name) => {
  const surfaces = path.join(testdata, name, 'surfaces.bin');
  const facets = path.join(testdata, name, 'facets.bin');
  if (!fs.existsSync(surfaces)) { throw new Error(`${surfaces} is missing`); }
  return { sidecar: readBuffer(surfaces), facets: fs.existsSync(facets) ? readBuffer(facets) : null };
};

const report = await runSelfCheck(parts, load, (line) => console.log(line));
process.exit(report.fail === 0 ? 0 : 1);
