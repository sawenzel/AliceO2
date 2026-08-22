#!/bin/bash
# Copy converted parts into website/testdata/ and write the manifest the page reads.
#
# Usage: ./fetch_testdata.sh [OPTIONS] <source-dir> [<source-dir> ...]
#
#   --parts LIST    comma-separated part names to copy (default: the six demo parts below)
#   --all           copy every part the source dirs hold
#   --group NAME    the heading these parts get in the part selector
#   --append        add to the existing testdata/manifest.json instead of replacing it
#   --list          print what the source dirs hold, copy nothing, and exit
#   -h, --help      this text
#
# Three source layouts are understood, and a source dir may be any of them:
#
#   gate workdir     <dir>/db/<part>/surfaces_<stem>.bin, facets_<stem>.bin   (runOracleGate.py)
#   per-part dirs    <dir>/<part>/surfaces_<stem>.bin, facets_<stem>.bin
#   flat converter   <dir>/surfaces_<stem>.bin, <dir>/facets_<stem>.bin       (O2_CADtoTGeo.py
#                    --output-folder: every part of a whole model in one directory)
#
# A part the converter declined for exact extraction has a facets_*.bin and no surfaces_*.bin.
# Such a part is copied and marked in the manifest as tessellated-only; the page then says so and
# turns its exact views off, which is honest coverage rather than a gap.
#
# The binaries are deliberately not committed (see .gitignore); this script is how a checkout gets
# them back. Parts not present in the given source dirs are skipped with a note.
#
# Example -- the six demo parts, then a representative ALICE3 selection on top:
#
#   ./fetch_testdata.sh --group fixtures --parts box,cyl_inter_cyl,torus_union_cyl,tube_window \
#       /path/to/gate_fixtures2
#   ./fetch_testdata.sh --append --group Bagger --parts Bucket,BoomCylinderInner \
#       /path/to/gate_bagger2
#   ./fetch_testdata.sh --append --group ALICE3 --parts ST1829909_002,ST2487458_01 \
#       /path/to/alice3_conv

set -u
here="$(cd "$(dirname "$0")" && pwd)"
dest="$here/testdata"

# The default demo set: a box (six planes, the simplest exact solid), the part whose *mesh leaks*,
# a torus part for the quartic, a cylinder with a spline-trimmed window and a hole, the Bagger
# Bucket (97 faces, spheres and tori) and one tube-tube seam part.
default_parts="box,cyl_inter_cyl,torus_union_cyl,tube_window,Bucket,BoomCylinderInner"

parts="$default_parts"
group=""
append=no
listonly=no
takeall=no

while [ $# -gt 0 ]; do
  case "$1" in
    --parts) parts="${2:-}"; shift 2 ;;
    --all) takeall=yes; shift ;;
    --group) group="${2:-}"; shift 2 ;;
    --append) append=yes; shift ;;
    --list) listonly=yes; shift ;;
    -h|--help) sed -n '2,45p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --*) echo "unknown option: $1" >&2; exit 2 ;;
    *) break ;;
  esac
done

if [ $# -lt 1 ]; then
  echo "usage: $0 [--parts a,b,c | --all] [--group NAME] [--append] [--list] <source-dir> ..." >&2
  exit 2
fi

# The part name is the stem minus the four trailing placement ids the converter appends.
partname() { echo "$1" | sed -E 's/(_[0-9]+){4}$//'; }

declare -A SURF FACE
keys=()

for dir in "$@"; do
  if [ ! -d "$dir" ]; then echo "skipping $dir: not a directory" >&2; continue; fi
  found=0
  while IFS= read -r file; do
    base="$(basename "$file")"
    parent="$(dirname "$file")"
    case "$base" in
      surfaces_*) kind=s; stem="${base#surfaces_}" ;;
      facets_*)   kind=f; stem="${base#facets_}" ;;
      *) continue ;;
    esac
    stem="${stem%.bin}"
    key="$parent|$stem"
    if [ -z "${SURF[$key]+x}" ] && [ -z "${FACE[$key]+x}" ]; then keys+=("$key"); fi
    if [ "$kind" = s ]; then SURF["$key"]="$file"; else FACE["$key"]="$file"; fi
    found=1
  done < <(find "$dir" -maxdepth 3 -type f \( -name 'surfaces_*.bin' -o -name 'facets_*.bin' \) | sort)
  [ "$found" = 1 ] || echo "skipping $dir: no surfaces_*.bin or facets_*.bin under it" >&2
done

if [ ${#keys[@]} -eq 0 ]; then
  echo "no parts found in: $*" >&2
  exit 1
fi

if [ "$listonly" = yes ]; then
  echo "parts available in: $*"
  for key in "${keys[@]}"; do
    stem="${key#*|}"
    what=""
    [ -n "${SURF[$key]+x}" ] && what="exact"
    [ -n "${FACE[$key]+x}" ] && what="${what:+$what+}mesh"
    printf '  %-34s %-11s %s\n' "$(partname "$stem")" "$what" "$stem"
  done | sort -u
  exit 0
fi

mkdir -p "$dest"
entries_file="$(mktemp)"
trap 'rm -f "$entries_file"' EXIT
copied=()
requested=()
if [ "$takeall" != yes ]; then
  IFS=',' read -ra requested <<< "$parts"
fi

for key in "${keys[@]}"; do
  stem="${key#*|}"
  part="$(partname "$stem")"
  if [ "$takeall" != yes ]; then
    keep=no
    for want in "${requested[@]}"; do [ "$want" = "$part" ] && keep=yes; done
    [ "$keep" = yes ] || continue
  fi
  mkdir -p "$dest/$part"
  surfaces=null
  facets=null
  # What the converter's cascade chose for this part IS whether it wrote a sidecar for it.
  ships=mesh
  if [ -n "${SURF[$key]+x}" ]; then
    cp "${SURF[$key]}" "$dest/$part/surfaces.bin"
    surfaces="\"$part/surfaces.bin\""
    ships=surface
  fi
  if [ -n "${FACE[$key]+x}" ]; then
    cp "${FACE[$key]}" "$dest/$part/facets.bin"
    facets="\"$part/facets.bin\""
  fi
  if [ "$surfaces" = null ] && [ "$facets" = null ]; then continue; fi
  printf '{"name": "%s", "stem": "%s", "group": "%s", "surfaces": %s, "facets": %s, "ships": "%s"}\n' \
    "$part" "$stem" "$group" "$surfaces" "$facets" "$ships" >> "$entries_file"
  copied+=("$part")
done

if [ ${#copied[@]} -eq 0 ]; then
  echo "none of the requested parts were found in: $*" >&2
  echo "  requested: $parts" >&2
  echo "  try: $0 --list $*" >&2
  exit 1
fi

# The manifest is JSON, so it is merged in python3 rather than pasted together in shell: --append
# has to read what is already there, replace an entry of the same name and keep the rest in order.
APPEND="$append" SOURCES="$*" python3 - "$dest/manifest.json" "$entries_file" <<'PYTHON'
import json, os, sys, datetime

manifest_path, entries_path = sys.argv[1], sys.argv[2]
append = os.environ.get("APPEND") == "yes"
sources = os.environ.get("SOURCES", "")

parts, previous_sources = [], ""
if append and os.path.exists(manifest_path):
    with open(manifest_path) as handle:
        old = json.load(handle)
    parts = old.get("parts", [])
    previous_sources = old.get("sources", "")

with open(entries_path) as handle:
    fresh = [json.loads(line) for line in handle if line.strip()]

by_name = {entry["name"]: index for index, entry in enumerate(parts)}
for entry in fresh:
    if entry["name"] in by_name:
        parts[by_name[entry["name"]]] = entry
    else:
        by_name[entry["name"]] = len(parts)
        parts.append(entry)

document = {
    "generated": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "sources": (previous_sources + " " + sources).strip() if append else sources,
    "parts": parts,
}
with open(manifest_path, "w") as handle:
    handle.write("{\n")
    handle.write('  "generated": %s,\n' % json.dumps(document["generated"]))
    handle.write('  "sources": %s,\n' % json.dumps(document["sources"]))
    handle.write('  "parts": [\n')
    handle.write(",\n".join("    " + json.dumps(entry) for entry in parts))
    handle.write("\n  ]\n}\n")
print("manifest holds %d part(s)" % len(parts))
PYTHON

echo "copied ${#copied[@]} part(s) into $dest:"
for p in "${copied[@]}"; do
  what="exact+mesh"
  [ -e "$dest/$p/surfaces.bin" ] || what="tessellated only"
  [ -e "$dest/$p/facets.bin" ] || what="exact, no mesh"
  printf '  %-34s %s\n' "$p" "$what"
done
du -sh "$dest"
echo "manifest: $dest/manifest.json"
