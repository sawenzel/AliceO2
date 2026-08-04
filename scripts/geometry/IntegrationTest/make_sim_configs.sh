#!/bin/bash
# Stage 3 of the CAD->TGeo Geant integration test: write the o2-sim configuration --
# the detector-list JSON and one ExternalModule JSON per representation.
#
# Usage: make_sim_configs.sh <workdir>     (the <workdir> given to convert.sh)
#
# Placement values are MEASURED, not assumed (Handoff section 4 requires this), with
# cadModuleProbe on the converted models:
#
#   ALICE3 IRIS (CAD_noETA.stp): local bbox is long in Y (836 cm, beam pipe included),
#     i.e. the CAD is authored with the beam along +Y; its own origin is the IP.
#     RotateX(90) maps local +Y -> world +Z. Translation is only the barrel-anchor
#     compensation: `barrel` sits at y = -30 cm in the cave (Cave.cxx), so +30 lands the
#     CAD origin on the ALICE origin. Measured ALICE-frame extents after placement:
#     x [-12.9, 8.7], y [-7.7, 34.2], z [-398.1, 439.6] cm (the transverse asymmetry is
#     the CAD's own content, present identically in its local frame).
#
#   oTOF (oTOF System V3-R92cm.step): local bbox x [-0.15, 347.8], y [-11.3, 182.0],
#     z [-96.6, 96.6] -- the cylinder axis runs along local X through (y, z) = (85.330, 0)
#     and the CAD frame is corner-anchored, NOT IP-anchored. RotateY(90) maps local +X ->
#     world -Z; the translation then recentres the axis onto the beam line and the barrel
#     longitudinally at z = 0: (0, 30 - 85.330, +173.805). Measured ALICE-frame extents
#     after placement: x [-96.6, 96.6], y [-96.6, 96.6], z [-174.0, 174.0] cm.
#
# Both fit far inside `barrel` (R 790.5 cm, z +-714.6 cm), and the two modules do not
# intersect: oTOF's shell starts at r ~ 85 cm, ALICE3's transverse extent ends at 34 cm.

set -euo pipefail
WORK=${1:?usage: make_sim_configs.sh <workdir>}
mkdir -p "$WORK/sim"

cat > "$WORK/sim/detlist.json" <<EOF
{ "INTEG": [ "A3IRIS", "OTOF" ] }
EOF

for repr in exact mesh; do
  cat > "$WORK/sim/extgeom_$repr.json" <<EOF
{
  "externalModules": [
    {
      "name": "A3IRIS",
      "title": "ALICE 3 IRIS vertex detector ($repr)",
      "macro": "$WORK/conv/alice3_$repr/geom.C",
      "anchor": "barrel",
      "placement": { "translation": [0, 30, 0], "rotation_deg": [90, 0, 0] }
    },
    {
      "name": "OTOF",
      "title": "outer TOF R92cm, whole module POLYSTYRENE ($repr)",
      "macro": "$WORK/conv/otof_$repr/geom.C",
      "anchor": "barrel",
      "placement": { "translation": [0, -55.330, 173.805], "rotation_deg": [0, 90, 0] }
    }
  ]
}
EOF
done
echo "wrote $WORK/sim/{detlist.json,extgeom_exact.json,extgeom_mesh.json}"
