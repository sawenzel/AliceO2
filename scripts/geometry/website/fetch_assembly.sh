#!/bin/bash
# Copy an assembly -- one placement table plus the body prototypes it instances -- into
# website/testdata/<name>/ and write the index the Assembly view reads.
#
# Usage: ./fetch_assembly.sh [OPTIONS] <placements.json> <facets-dir> [<facets-dir> ...]
#
#   --name NAME     the directory under testdata/ and the selector entry (default: otof_assembly)
#   --label TEXT    the name shown in the part selector (default: "oTOF assembly")
#   --group NAME    the heading the selector entry gets (default: the oTOF group)
#   -h, --help      this text
#
# The placement table is the converter's own dump of the placed tree, and is copied verbatim:
#
#   {"units": "cm",
#    "leaves": [{"entry": "0:1:1:4", "name": "Module v1", "count": 3672,
#                "matrices": [[[3x4 row-major rotation | translation in cm]], ...]}]}
#
# Each leaf's bodies are the converter's per-body prototypes for that leaf, named
# <leaf name with spaces as underscores>_<entry with colons as underscores>_b<N> -- so leaf
# "Module v1" at "0:1:1:4" owns facets_Module_v1_0_1_1_4_b1.bin .. _b17.bin. A body's vertices are
# already in its LEAF's frame (its own pose is baked in by the converter), so a body instance's
# world transform is exactly its leaf's matrix and there is no per-body transform to carry.
#
# A leaf with no bodies is kept in the index with an empty body list: it is an empty node of the
# tree, it places nothing, and the placement count on the page says so rather than quietly
# dropping it.
#
# Written out (testdata/ is not committed -- see .gitignore -- so this script is how a checkout
# gets the assembly back):
#
#   testdata/<name>/placements.json          the placement table, byte-for-byte
#   testdata/<name>/facets_<body>.bin        one tessellated prototype per body
#   testdata/<name>/index.json               leaves -> bodies -> facet file, and the totals
#
# Example -- the oTOF barrel, from the converter output this page was built against:
#
#   ./fetch_assembly.sh /path/to/otof_assembly_placements.json /path/to/otof/otof_mesh

set -u
here="$(cd "$(dirname "$0")" && pwd)"

name=otof_assembly
label="oTOF assembly"
group="oTOF (STEP_examples/oTOF System V3-R92cm.step)"

while [ $# -gt 0 ]; do
  case "$1" in
    --name) name="${2:-}"; shift 2 ;;
    --label) label="${2:-}"; shift 2 ;;
    --group) group="${2:-}"; shift 2 ;;
    -h|--help) sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --*) echo "unknown option: $1" >&2; exit 2 ;;
    *) break ;;
  esac
done

if [ $# -lt 2 ]; then
  echo "usage: $0 [--name NAME] [--label TEXT] [--group NAME] <placements.json> <facets-dir> ..." >&2
  exit 2
fi

placements="$1"; shift
if [ ! -f "$placements" ]; then echo "$placements: not a file" >&2; exit 1; fi

dest="$here/testdata/$name"
mkdir -p "$dest"
cp "$placements" "$dest/placements.json"

LABEL="$label" GROUP="$group" NAME="$name" SOURCES="$placements $*" \
python3 - "$dest" "$@" <<'PYTHON'
import datetime, glob, json, os, re, shutil, struct, sys

dest, facet_dirs = sys.argv[1], sys.argv[2:]
label = os.environ["LABEL"]
group = os.environ["GROUP"]
name = os.environ["NAME"]
sources = os.environ["SOURCES"]

with open(os.path.join(dest, "placements.json")) as handle:
    table = json.load(handle)

def triangles(path):
    """The triangle count in the facets_*.bin header (uint32 nTriangles, then 9 float32 each)."""
    with open(path, "rb") as handle:
        count = struct.unpack("<I", handle.read(4))[0]
    size = os.path.getsize(path)
    if size < 4 + 36 * count:
        raise SystemExit("%s: header claims %d triangles, file holds %d bytes" % (path, count, size))
    return count

def body_index(path):
    match = re.search(r"_b(\d+)\.bin$", path)
    return int(match.group(1)) if match else 0

leaves, total_solids, total_triangles, total_placements = [], 0, 0, 0
for leaf in table["leaves"]:
    prefix = "%s_%s_b" % (leaf["name"].replace(" ", "_"), leaf["entry"].replace(":", "_"))
    found = []
    for directory in facet_dirs:
        found += glob.glob(os.path.join(directory, "facets_%s*.bin" % prefix))
    # One body may be offered by several source dirs; take the one with the most triangles, which
    # is the properly meshed copy rather than a coarse fallback.
    best = {}
    for path in sorted(found):
        key = os.path.basename(path)
        n = triangles(path)
        if key not in best or n > best[key][1]:
            best[key] = (path, n)
    bodies = []
    for key in sorted(best, key=lambda k: body_index(k)):
        path, n = best[key]
        shutil.copyfile(path, os.path.join(dest, key))
        bodies.append({"name": key[len("facets_"):-len(".bin")], "facets": key, "triangles": n})
    placed = len(leaf["matrices"])
    total_placements += placed
    total_solids += placed * len(bodies)
    total_triangles += placed * sum(body["triangles"] for body in bodies)
    leaves.append({"entry": leaf["entry"], "name": leaf["name"], "prefix": prefix,
                   "placements": placed, "bodies": bodies})
    if not bodies:
        print("  leaf %s (%s): no bodies under facets_%s*.bin -- kept as an empty node"
              % (leaf["entry"], leaf["name"], prefix))

document = {
    "generated": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "sources": sources,
    "name": name,
    "label": label,
    "group": group,
    "units": table.get("units", "cm"),
    "placements": "placements.json",
    "leaves": leaves,
    "totals": {
        "prototypes": sum(len(leaf["bodies"]) for leaf in leaves),
        "placements": total_placements,
        "solids": total_solids,
        "triangles": total_triangles,
    },
}
with open(os.path.join(dest, "index.json"), "w") as handle:
    json.dump(document, handle, indent=1)
    handle.write("\n")

for leaf in leaves:
    print("  %-9s %-14s %5d placements x %2d bodies = %7d solids"
          % (leaf["entry"], leaf["name"], leaf["placements"], len(leaf["bodies"]),
             leaf["placements"] * len(leaf["bodies"])))
print("%d prototypes, %d leaf placements, %d solids, %d triangles"
      % (document["totals"]["prototypes"], document["totals"]["placements"],
         document["totals"]["solids"], document["totals"]["triangles"]))
PYTHON

echo "assembly '$name' in $dest"
du -sh "$dest"
