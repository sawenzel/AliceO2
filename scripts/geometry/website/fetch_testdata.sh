#!/bin/bash
# Copy a small representative set of gate outputs into website/testdata/ and write its manifest.
#
# Usage: ./fetch_testdata.sh <gate-workdir> [<gate-workdir> ...]
#
# A gate workdir is what runOracleGate.py produces: it contains a db/ directory with one
# subdirectory per part, each holding surfaces_<PART>.bin and facets_<PART>.bin. Typically:
#
#   ./fetch_testdata.sh /path/to/gate_fixtures2 /path/to/gate_bagger2
#
# The binaries are deliberately not committed (see .gitignore); this script is how a checkout
# gets them back. Parts not present in the given workdirs are skipped with a note.

set -u
here="$(cd "$(dirname "$0")" && pwd)"
dest="$here/testdata"

# The demo set: a box (six planes, the simplest exact solid), the part whose *mesh leaks*, a
# torus part for the quartic, a cylinder with a spline-trimmed window and a hole, the Bagger
# Bucket (97 faces, spheres and tori) and one tube-tube seam part.
wanted=(box cyl_inter_cyl torus_union_cyl tube_window Bucket BoomCylinderInner)

if [ $# -lt 1 ]; then
  echo "usage: $0 <gate-workdir> [<gate-workdir> ...]" >&2
  exit 2
fi

mkdir -p "$dest"
entries=()
found=()

for workdir in "$@"; do
  db="$workdir/db"
  [ -d "$db" ] || db="$workdir"
  [ -d "$db" ] || { echo "skipping $workdir: no db/ directory" >&2; continue; }
  for partdir in "$db"/*/; do
    [ -d "$partdir" ] || continue
    for surfaces in "$partdir"surfaces_*.bin; do
      [ -e "$surfaces" ] || continue
      base="$(basename "$surfaces")"
      stem="${base#surfaces_}"; stem="${stem%.bin}"
      # stem is <PART>_<a>_<b>_<c>_<lid>; the part name is everything before the four trailing ids
      part="$(echo "$stem" | sed -E 's/(_[0-9]+){4}$//')"
      keep=no
      for w in "${wanted[@]}"; do [ "$w" = "$part" ] && keep=yes; done
      [ "$keep" = yes ] || continue
      facets="$partdir/facets_$stem.bin"
      mkdir -p "$dest/$part"
      cp "$surfaces" "$dest/$part/surfaces.bin"
      if [ -e "$facets" ]; then cp "$facets" "$dest/$part/facets.bin"; fi
      found+=("$part")
      hasmesh=false
      [ -e "$dest/$part/facets.bin" ] && hasmesh=true
      entries+=("    {\"name\": \"$part\", \"stem\": \"$stem\", \"surfaces\": \"$part/surfaces.bin\", \"facets\": $( [ "$hasmesh" = true ] && echo "\"$part/facets.bin\"" || echo null )}")
    done
  done
done

if [ ${#entries[@]} -eq 0 ]; then
  echo "no wanted parts found in: $*" >&2
  exit 1
fi

{
  echo '{'
  echo '  "generated": "'"$(date -u +%Y-%m-%dT%H:%M:%SZ)"'",'
  echo '  "sources": "'"$*"'",'
  echo '  "parts": ['
  printf '%s,\n' "${entries[@]}" | sed '$ s/,$//'
  echo '  ]'
  echo '}'
} > "$dest/manifest.json"

echo "copied ${#found[@]} part(s) into $dest:"
du -sh "$dest"
for p in "${found[@]}"; do echo "  $p"; done
echo "manifest: $dest/manifest.json"
