#!/usr/bin/env python3
"""One module through the round trip, anchored where the source geometry put it.

  roundtrip_module.py <studydir> <MODULE> [csg,mesh]

Why this is not one conversion per module.  o2-sim always builds the experiment
hall -- `cave`, `barrel` (at y = -30 in cave) and `caveRB24` -- whatever module
list is asked for, and those are real solids.  Converting a module from `cave`
downward and anchoring the result at `cave` puts its detectors geometrically
inside the native `barrel` while making them daughters of `cave`, and the
navigator, having descended into `barrel`, never sees them: measured, a ray from
the IP crossed `barrel` then `cave` and nothing else, and transport fell to 161
steps per event against the baseline's 36000.

So each module is converted per hall anchor instead, and each piece is anchored
where the source hung it:

  * everything under `barrel` -- which is all of ITS, TPC and MAG, and 8 of the
    beam pipe's 10 subtrees -- becomes one conversion with `--top barrel`, the
    hall volume hollowed so it contributes no second copy of the air, and placed
    back into the real `barrel` with the identity;
  * anything under `cave` or `caveRB24` is converted from its own subtree root
    and placed with that root's own matrix.

Writes <studydir>/cad/<MODULE>/ and a module_entries.json fragment that
make_configs.py assembles into the o2-sim external-geometry file.
"""

import json
import math
import os
import shlex
import subprocess
import sys

HALL = ("cave", "barrel", "caveRB24")


def sh(cmd, cwd, log):
    """Run one step in its own shell, with its output kept in a log."""
    with open(os.path.join(cwd, log), "w") as fh:
        r = subprocess.run(["bash", "-c", cmd], cwd=cwd, stdout=fh,
                           stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(os.path.join(cwd, log)).read()[-3000:])
        raise SystemExit(f"step failed ({r.returncode}): {cmd[:120]}")


def euler_deg(rot, env_o2):
    """The rotation_deg triple ExternalModule's JSON wants, verified.

    makePlacementFromJSON applies RotateX, then RotateY, then RotateZ to a fresh
    TGeoCombiTrans. Rather than trust an analytic decomposition against ROOT's
    own multiplication order, the candidate is rebuilt through ROOT and compared
    with the target -- a placement that silently differs would move a whole
    subdetector, and nothing downstream would say so.
    """
    if all(abs(rot[i] - (1.0 if i in (0, 4, 8) else 0.0)) < 1e-12 for i in range(9)):
        return None                       # identity: omit the rotation entirely
    # ROOT lives in the o2 environment, not in this driver's interpreter, so the
    # candidate is rebuilt there.
    probe = f"""
import ROOT, json, itertools, sys
target = {list(rot)!r}
for cand in itertools.product((0,90,-90,180),repeat=3):
    c = ROOT.TGeoCombiTrans()
    c.RotateX(cand[0]); c.RotateY(cand[1]); c.RotateZ(cand[2])
    m = c.GetRotationMatrix()
    if all(abs(m[i]-target[i]) < 1e-9 for i in range(9)):
        print(json.dumps(list(cand))); sys.exit(0)
sys.exit(3)
"""
    r = subprocess.run(["bash", "-c", f'{env_o2}; python3 -c {shlex.quote(probe)}'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"cannot express this rotation as rotation_deg: {rot}\n"
                         "ExternalModule's JSON carries Euler angles only; this "
                         "placement needs a full matrix and the loader does not "
                         "take one yet.")
    return json.loads(r.stdout.strip().splitlines()[-1])


def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    study, mod = os.path.abspath(sys.argv[1]), sys.argv[2]
    d = os.path.join(study, "cad", mod)
    os.makedirs(d, exist_ok=True)
    ct = os.path.dirname(os.path.abspath(__file__))
    env_o2 = f'source "{study}/env_o2.sh" >/dev/null 2>&1'
    env_cv = f'{env_o2}; source "{study}/env_converter.sh"'
    py = '"$SW/Python/latest/bin/python3.10"'

    print(f"=== {mod}: the source geometry")
    sh(f'{env_o2}; o2-sim-serial -n 0 -g boxgen -m {mod} -o o2sim', d, "geom.log")

    print(f"=== {mod}: where does it hang itself?")
    sh(f'{env_o2}; python3 "{ct}/module_anchors.py" o2sim_geometry.root '
       f'--json anchors.json', d, "anchors.log")
    roots = json.load(open(os.path.join(d, "anchors.json")))["roots"]
    in_barrel = [r for r in roots if r["anchor"] == "barrel"]
    elsewhere = [r for r in roots if r["anchor"] != "barrel"]
    print(f"    {len(in_barrel)} subtree(s) under barrel, "
          f"{len(elsewhere)} elsewhere: "
          f"{[(r['volume'], r['anchor']) for r in elsewhere]}")

    entries = {}
    variants = sys.argv[3].split(",") if len(sys.argv) > 3 else ["csg", "mesh"]

    def convert(tag, top, hollow, anchor, placement, variant="csg"):
        """One --top conversion plus its media sidecar, scored.

        `variant` picks the back-conversion: "csg" is the shipped cascade (CSG,
        else exact surfaces, else mesh) and "mesh" is tessellated-only, which is
        the fallback every other CAD pipeline uses and the natural benchmark to
        score the exact path against.
        """
        out = f"conv_{tag}" if variant == "csg" else f"conv_{variant}_{tag}"
        cascade = ("--csg auto --exact-surfaces auto --mesh" if variant == "csg"
                   else "--mesh")
        hollow_args = " ".join(f"--hollow-volume {h}" for h in hollow)
        tagarg = f'--hollow-tag {mod}' if hollow else ""
        print(f"=== {mod}: --top {top} -> anchor {anchor}")
        # No --carve-mothers. TGeo carves a mother implicitly -- a daughter takes
        # precedence over its mother's solid -- and the converter restores that
        # nesting from the sidecar, which is exact and costs no boolean. Carving
        # is also not a substitute: it subtracts daughter SOLIDS, and an assembly
        # daughter has none, so ITSUWrapVol0 (whose only daughter is an assembly)
        # came back uncarved and swallowed the whole inner barrel.
        sh(f'{env_cv}; {py} "{ct}/../O2_TGeoToCAD.py" o2sim_geometry.root {tag}.step '
           f'--top {top} --report {tag}_writer_report.json '
           f'--media-json {tag}_media.json {hollow_args} {tagarg}',
           d, f"writer_{tag}.log")
        sh(f'{env_cv}; {py} "{ct}/../O2_CADtoTGeo.py" {tag}.step -o geom.C '
           f'--output-folder {out} {cascade} '
           f'--media-json {tag}_media.json', d, f"{out}.log")
        for line in open(os.path.join(d, f"{out}.log")):
            if "tiers:" in line or "Media from sidecar" in line or "[WARN]" in line:
                print("   ", line.rstrip())
        sh(f'{env_o2}; python3 "{ct}/check_media.py" --original o2sim_geometry.root '
           f'--macro {out}/geom.C --rtol 1e-6 '
           f'--writer-report {tag}_writer_report.json --json media_{out}.json',
           d, f"media_{out}.log")
        for line in open(os.path.join(d, f"media_{out}.log")):
            if line.startswith(("VERDICT", "  media identical", "  left on")):
                print("   ", line.rstrip())
        e = {"tag": tag, "macro": os.path.join(d, out, "geom.C"), "anchor": anchor}
        if placement:
            e["placement"] = placement
        entries.setdefault(variant, []).append(e)

    # The STEP is written once per anchor; only the back-conversion differs, so a
    # variant comparison is a comparison of representations and not of two exports.
    for v in variants:
        if in_barrel:
            convert("barrel", "barrel", ["barrel"], "barrel", None, v)
        for r in elsewhere:
            rot = euler_deg(r["rotation"], env_o2)
            pl = {"translation": [float(x) for x in r["translation"]]}
            if rot:
                pl["rotation_deg"] = rot
            convert(r["volume"], r["volume"], [], r["anchor"], pl, v)

    with open(os.path.join(d, "module_entries.json"), "w") as fh:
        json.dump({"module": mod, "entries": entries}, fh, indent=2)
    print(f"=== {mod}: " + ", ".join(f"{len(v)} {k} placement(s)"
                                     for k, v in entries.items())
          + f" -> {d}/module_entries.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
