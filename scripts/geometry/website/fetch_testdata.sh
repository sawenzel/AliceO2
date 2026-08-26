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
# Seven artefacts per part are picked up where the source layout has them:
#
#   surfaces_<stem>.bin   the exact trimmed analytic faces  -> <part>/surfaces.bin
#   facets_<stem>.bin     the tessellation                  -> <part>/facets.bin
#   flatcsg_<stem>.bin    the flat halfspace solid, a union of intersection-cells over signed
#                         implicit halfspaces                -> <part>/flatcsg.bin
#   shape_<stem>.root     the recognised CSG composite, as a TGeoShape the bridge can trace
#                                                           -> <part>/shape.root
#   treeshape_<stem>.root the SAME cells as a plain TGeoCompositeShape (mktree.py)
#                                                           -> <part>/cellstree.root
#   original_<stem>.root  the TGeoShape the part was made FROM, before the round trip, written by
#                         exportSourceShapes.py             -> <part>/original.root
#   csg_<stem>.json       what the recogniser found and how the acceptance test judged it
#                                                           -> <part>/csg.json
#
# A part the converter declined for exact extraction has a facets_*.bin and no surfaces_*.bin.
# Such a part is copied and marked in the manifest as tessellated-only; the page then says so and
# turns its exact views off, which is honest coverage rather than a gap. In the same way a part
# with no shape_*.root gets the CSG views turned off, and a csg_*.json with a null candidate is
# still copied: it is the recogniser saying, in its own words, that it found nothing.
#
# One rule is not a preference and is worth stating here. A part that ships the flat solid ALSO
# has a shape_<stem>.root beside it, and that file is the O2FlatCSG streamed -- the same solid, not
# a CSG composite. Copying it as <part>/shape.root would make the page claim the part ships CSG
# when it does not, so when a flatcsg_*.bin is present the shape_*.root is skipped.
#
# The three .root artefacts are kept in three separate slots BECAUSE they are three different
# things, and calling all of them "CSG" is what makes a page unreadable:
#
#   shape.root      what the converter SHIPS for this part, when the recogniser produced a tree
#   cellstree.root  the decomposition emitted as a plain composite. Nothing ships this. It exists
#                   so the flat solid has something to be measured against -- the cells the
#                   cascade refused to emit as a tree because they would make one wider than the
#                   routing threshold
#   original.root   what the part was made FROM: the shape the detector geometry ships today,
#                   before it ever went to STEP. Not a product of the conversion at all, and the
#                   only baseline that answers "is any of this better than what we already have"
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
    -h|--help) sed -n '2,66p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
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

declare -A SURF FACE SHAPE TREE ORIG FLAT CSGJ
keys=()

for dir in "$@"; do
  if [ ! -d "$dir" ]; then echo "skipping $dir: not a directory" >&2; continue; fi
  found=0
  while IFS= read -r file; do
    base="$(basename "$file")"
    parent="$(dirname "$file")"
    case "$base" in
      surfaces_*)  kind=s; stem="${base%.bin}"; stem="${stem#surfaces_}" ;;
      facets_*)    kind=f; stem="${base%.bin}"; stem="${stem#facets_}" ;;
      flatcsg_*)   kind=x; stem="${base%.bin}"; stem="${stem#flatcsg_}" ;;
      treeshape_*) kind=t; stem="${base%.root}"; stem="${stem#treeshape_}" ;;
      original_*)  kind=o; stem="${base%.root}"; stem="${stem#original_}" ;;
      shape_*)     kind=r; stem="${base%.root}"; stem="${stem#shape_}" ;;
      # csg_report.json is the whole model's report, not a part's record.
      csg_report.json) continue ;;
      csg_*)       kind=c; stem="${base%.json}"; stem="${stem#csg_}" ;;
      *) continue ;;
    esac
    key="$parent|$stem"
    if [ -z "${SURF[$key]+x}" ] && [ -z "${FACE[$key]+x}" ] && [ -z "${SHAPE[$key]+x}" ] \
       && [ -z "${TREE[$key]+x}" ] && [ -z "${ORIG[$key]+x}" ] && [ -z "${FLAT[$key]+x}" ] \
       && [ -z "${CSGJ[$key]+x}" ]; then keys+=("$key"); fi
    case "$kind" in
      s) SURF["$key"]="$file" ;;
      f) FACE["$key"]="$file" ;;
      r) SHAPE["$key"]="$file" ;;
      t) TREE["$key"]="$file" ;;
      o) ORIG["$key"]="$file" ;;
      x) FLAT["$key"]="$file" ;;
      c) CSGJ["$key"]="$file" ;;
    esac
    found=1
  done < <(find "$dir" -maxdepth 3 -type f \( -name 'surfaces_*.bin' -o -name 'facets_*.bin' \
             -o -name 'flatcsg_*.bin' -o -name 'shape_*.root' -o -name 'treeshape_*.root' \
             -o -name 'original_*.root' -o -name 'csg_*.json' \) | sort)
  [ "$found" = 1 ] || echo "skipping $dir: no surfaces_*.bin, facets_*.bin, flatcsg_*.bin, shape_*.root, treeshape_*.root, original_*.root or csg_*.json under it" >&2
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
    [ -n "${FLAT[$key]+x}" ] && what="${what:+$what+}flatcsg"
    [ -n "${SHAPE[$key]+x}" ] && what="${what:+$what+}csg"
    [ -n "${TREE[$key]+x}" ] && what="${what:+$what+}cellstree"
    [ -n "${ORIG[$key]+x}" ] && what="${what:+$what+}original"
    printf '  %-34s %-11s %s\n' "$(partname "$stem")" "$what" "$stem"
  done | sort -u
  exit 0
fi

mkdir -p "$dest"
entries_file="$(mktemp)"
trap 'rm -f "$entries_file"' EXIT
copied=()
stale=()
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
  shape=null
  cellstree=null
  original=null
  flatcsg=null
  csg=null
  # What the converter's cascade chose for this part IS which artefact it wrote for it: a
  # shape_*.root means the CSG recogniser was accepted, a flatcsg_*.bin means the flat halfspace
  # solid, a surfaces_*.bin means the exact solid. The order below is the cascade's own, so the
  # last assignment to `ships` wins the way the cascade wins.
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
  if [ -n "${FLAT[$key]+x}" ]; then
    cp "${FLAT[$key]}" "$dest/$part/flatcsg.bin"
    flatcsg="\"$part/flatcsg.bin\""
    ships=flatcsg
  fi
  # The header explains this: a shape_*.root beside a flat sidecar is the flat solid streamed,
  # not a CSG composite, so it is skipped rather than presented as one.
  if [ -z "${FLAT[$key]+x}" ] && [ -n "${SHAPE[$key]+x}" ]; then
    cp "${SHAPE[$key]}" "$dest/$part/shape.root"
    shape="\"$part/shape.root\""
    ships=shape
  fi
  if [ -n "${TREE[$key]+x}" ]; then
    cp "${TREE[$key]}" "$dest/$part/cellstree.root"
    cellstree="\"$part/cellstree.root\""
  fi
  if [ -n "${ORIG[$key]+x}" ]; then
    cp "${ORIG[$key]}" "$dest/$part/original.root"
    original="\"$part/original.root\""
  fi
  if [ -n "${CSGJ[$key]+x}" ]; then
    cp "${CSGJ[$key]}" "$dest/$part/csg.json"
    csg="\"$part/csg.json\""
  fi
  if [ "$surfaces" = null ] && [ "$facets" = null ] && [ "$flatcsg" = null ]; then continue; fi
  # An artefact file this run did NOT write, left behind by an earlier fetch with a different
  # source or a different layout, is not referenced by the manifest and so is invisible to the
  # page -- but it is very visible to anyone reading the directory, and it made this script's own
  # summary claim a part had a CSG shape when the manifest said it had none. Nothing is deleted
  # here; the stale names are collected and printed at the end so they can be dealt with.
  for artefact in surfaces.bin facets.bin flatcsg.bin shape.root cellstree.root original.root csg.json; do
    var=""
    case "$artefact" in
      surfaces.bin) var="$surfaces" ;; facets.bin) var="$facets" ;; flatcsg.bin) var="$flatcsg" ;;
      shape.root) var="$shape" ;; cellstree.root) var="$cellstree" ;; original.root) var="$original" ;;
      csg.json) var="$csg" ;;
    esac
    if [ "$var" = null ] && [ -e "$dest/$part/$artefact" ]; then
      stale+=("$part/$artefact")
    fi
  done
  printf '{"name": "%s", "stem": "%s", "group": "%s", "surfaces": %s, "facets": %s, "shape": %s, "cellstree": %s, "original": %s, "flatcsg": %s, "csg": %s, "ships": "%s"}\n' \
    "$part" "$stem" "$group" "$surfaces" "$facets" "$shape" "$cellstree" "$original" "$flatcsg" "$csg" "$ships" >> "$entries_file"
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

def cellstree_note(entry):
    """What a part's cellstree.root is, in the counts that matter.

    The number the measured crossover is expressed in is that tree's LEAF count -- not the
    halfspace count, because `_cell_leaf` folds six planes into one TGeoBBox. The candidate carries
    the OCCT realisation's cells, so the leaf count is read from there rather than guessed."""
    if not entry.get("cellstree") or not entry.get("csg"):
        return None
    try:
        payload = json.loads(open(os.path.join(os.path.dirname(manifest_path), entry["csg"])).read())
        cells = payload["candidate"]["notes"]["occCells"]
        leaves = sum(len(cell["leaves"]) for cell in cells)
    except Exception:
        return "the decomposed cells, emitted as a plain TGeoCompositeShape"
    return ("the same %d cells as a plain TGeoCompositeShape, %d leaves" % (len(cells), leaves))


def original_note(entry):
    """What a part's original.root is: the source volume, its class and how its tree is shaped.

    `exportSourceShapes.py` measured all of that when it wrote the file and left it in
    `original_report.json` beside it, so this reads that record rather than opening a ROOT file --
    which this step, running in whatever python the shell has, cannot do."""
    if not entry.get("original"):
        return None
    for source in sources.split():
        path = os.path.join(source, "original_report.json")
        if not os.path.exists(path):
            continue
        try:
            report = json.loads(open(path).read())
        except Exception:
            continue
        for row in report.get("parts", []):
            if row.get("part") != entry["stem"] or not row.get("written"):
                continue
            classes = ", ".join("%d %s" % (n, k) for k, n in
                                sorted(row.get("leafClasses", {}).items(), key=lambda kv: -kv[1]))
            return ("%s from volume %s, boolean depth %s, %s leaves (%s)" %
                    (row.get("class"), row.get("sourceVolume"), row.get("booleanDepth"),
                     row.get("leaves"), classes))
    return "the TGeoShape this part was made from, before the round trip"


by_name = {entry["name"]: index for index, entry in enumerate(parts)}
for entry in fresh:
    notes = {}
    for key, fn in (("cellstree", cellstree_note), ("original", original_note)):
        note = fn(entry)
        if note:
            notes[key] = note
    if notes:
        entry["notes"] = notes
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
# What each part got is read back out of the manifest this run just wrote, not off the
# filesystem: a stale file from an earlier fetch is not what this part carries.
python3 - "$dest/manifest.json" "${copied[@]}" <<'PYTHON'
import json, sys
manifest = json.load(open(sys.argv[1]))
wanted = set(sys.argv[2:])
for entry in manifest["parts"]:
    if entry["name"] not in wanted:
        continue
    bits = []
    bits.append("exact" if entry.get("surfaces") else "tessellated only")
    if entry.get("facets"):
        bits[0] += "+mesh"
    elif entry.get("surfaces"):
        bits[0] += ", no mesh"
    for field, label in (("flatcsg", "FlatCSG"), ("shape", "CSG shape"),
                         ("cellstree", "cells-tree"), ("original", "original")):
        if entry.get(field):
            bits.append(label)
    if entry.get("csg") and not entry.get("shape") and not entry.get("flatcsg"):
        bits.append("CSG declined")
    print("  %-34s %s" % (entry["name"], ", ".join(bits)))
PYTHON
if [ ${#stale[@]} -gt 0 ]; then
  echo
  echo "note: ${#stale[@]} artefact file(s) in $dest are NOT referenced by the manifest, left by an"
  echo "      earlier fetch. The page ignores them; remove them if you want the directory to match:"
  for f in "${stale[@]}"; do echo "        $f"; done
fi
du -sh "$dest"
echo "manifest: $dest/manifest.json"
