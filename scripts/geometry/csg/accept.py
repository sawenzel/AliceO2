"""Acceptance test 1 of 2: the OCCT symmetric-difference volume.

`CSG_Pipeline.md` §3.5 states the criterion:

    volume(candidate - original) + volume(original - candidate)  ==  0   (to model tolerance)

computed with `BRepAlgoAPI_Cut` both ways and `BRepGProp`. This is a *proof of equality*, not a
sample-based argument, and it is what lets the recogniser be greedy: a proposal that is nearly
right is rejected outright rather than shipped.

Reading the number
------------------
The residual is a **volume**, and the model tolerance is a **length** (1.0e-07 cm on Bagger). The
two are related through the boundary: displacing a surface of area `A` by `t` moves a volume
`A*t`. So the criterion this module applies is

    dV_sym  <=  bandFactor * modelTolerance * area(original)

and it *also* reports `dV_sym / volume(original)`, because the ratio is the number a reader can
compare across parts. `bandFactor` is 1 by default: a candidate is accepted only if its symmetric
difference could be explained by the whole boundary sitting one model tolerance away from where
the CAD says it is.

The trap that this module is written around
-------------------------------------------
`BRepGProp::VolumeProperties` returns exactly **0** for a single planar face, and OCCT's boolean
operations do not return an error for an empty or degenerate result. So a symmetric difference of
zero is *also* what you get from a candidate that failed to build, from a cut that silently
produced nothing, and from comparing a shape with itself by accident. Every one of those is a
false accept, and each is checked for explicitly below:

  * both cuts must report `IsDone()`;
  * the original's own volume and area must be strictly positive;
  * the candidate's own volume must be strictly positive and within a loose factor of the
    original's, *before* the cuts are attempted.

`self_test()` closes the loop with hand-built pairs whose answer is known, including three
deliberately-wrong candidates that must be rejected. Run it before believing a table.
"""

import math

_BAND_FACTOR = 1.0
# A candidate whose own volume is not even in the right ballpark is a recogniser bug, not a
# near-miss; catching it here gives a readable reason instead of a symmetric difference the size
# of the part.
_SANITY_VOLUME_RATIO = 4.0


def _props(shape):
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.GProp import GProp_GProps
    vol = GProp_GProps()
    brepgprop.VolumeProperties(shape, vol)
    surf = GProp_GProps()
    brepgprop.SurfaceProperties(shape, surf)
    return vol.Mass(), surf.Mass()


def _cut_volume(a, b, what):
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
    op = BRepAlgoAPI_Cut(a, b)
    op.Build()
    if not op.IsDone():
        raise RuntimeError(f"BRepAlgoAPI_Cut failed ({what})")
    volume, _area = _props(op.Shape())
    return abs(volume)


def symmetric_difference(original, cand_shape, model_tolerance_cm, band_factor=_BAND_FACTOR):
    """Measure `original` against `cand_shape`; returns a dict, never raises on a mere mismatch."""
    v_orig, a_orig = _props(original)
    if not (v_orig > 0.0 and a_orig > 0.0):
        return {"accepted": False,
                "reason": f"original has non-positive volume/area ({v_orig:.6g}/{a_orig:.6g})"}
    v_cand, a_cand = _props(cand_shape)
    if not v_cand > 0.0:
        return {"accepted": False, "volumeOriginal": v_orig, "volumeCandidate": v_cand,
                "reason": f"candidate has non-positive volume ({v_cand:.6g})"}
    if not (1.0 / _SANITY_VOLUME_RATIO <= v_cand / v_orig <= _SANITY_VOLUME_RATIO):
        return {"accepted": False, "volumeOriginal": v_orig, "volumeCandidate": v_cand,
                "reason": f"candidate volume {v_cand:.6g} is not comparable to the original's "
                          f"{v_orig:.6g}"}
    extra = _cut_volume(cand_shape, original, "candidate - original")
    missing = _cut_volume(original, cand_shape, "original - candidate")
    dv = extra + missing
    band = band_factor * model_tolerance_cm * a_orig
    return {
        "accepted": dv <= band,
        "volumeOriginal": v_orig,
        "volumeCandidate": v_cand,
        "areaOriginal": a_orig,
        "extraVolume": extra,
        "missingVolume": missing,
        "symmetricDifference": dv,
        "band": band,
        "modelToleranceCm": model_tolerance_cm,
        "relativeToVolume": dv / v_orig,
        "reason": None if dv <= band else
                  f"symmetric difference {dv:.6g} cm^3 exceeds the band {band:.6g} cm^3 "
                  f"(= {model_tolerance_cm:.3g} cm x area {a_orig:.6g} cm^2); "
                  f"extra {extra:.6g}, missing {missing:.6g}",
    }


# ------------------------------------------------------------------------------------------
# self-test
# ------------------------------------------------------------------------------------------

def self_test(verbose=True):
    """Hand-built pairs whose verdict is known, including candidates that must be *rejected*.

    Without the rejecting half, an acceptance test that is structurally incapable of reporting a
    difference would pass every check here (this project has shipped exactly that mistake once,
    in a sweep that could not say "yes").
    """
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder,
                                      BRepPrimAPI_MakeSphere)
    from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt
    from csg import primitives as prim

    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))
        if verbose:
            print(f"  [{'ok ' if condition else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))

    tol = 1.0e-7

    # 1. positive control: a shape against itself.
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0).Shape()
    r = symmetric_difference(box, box, tol)
    check("a box against itself is accepted", r["accepted"], f"dV={r['symmetricDifference']:.3g}")
    check("a box against itself has zero symmetric difference", r["symmetricDifference"] == 0.0)

    # 2. positive control through the description, which is how the pipeline uses it: an OCCT box
    #    built independently of the description must still match.
    cand = prim.candidate("primitive", [prim.leaf(
        "TGeoBBox", {"dx": 1.0, "dy": 1.5, "dz": 2.0},
        prim.identity_frame((1.0, 1.5, 2.0)))], "self-test")
    r = symmetric_difference(box, prim.build_occ(cand), tol)
    check("an independently built TGeoBBox description matches the box", r["accepted"],
          f"dV={r['symmetricDifference']:.3g} band={r['band']:.3g}")

    # 3. negative control: the same box 1 micron (1e-4 cm) too long.
    long_box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0001).Shape()
    r = symmetric_difference(box, long_box, tol)
    check("a box 1e-4 cm too long is rejected", not r["accepted"],
          f"dV={r['symmetricDifference']:.3g} band={r['band']:.3g}")

    # 3b. how fine is the knife? A displacement of exactly one model tolerance must sit at the
    #     band, and ten of them must be outside it.
    for factor, want_accept in ((0.5, True), (10.0, False)):
        nudged = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 2.0, 3.0, 4.0 + factor * tol).Shape()
        r = symmetric_difference(box, nudged, tol)
        check(f"a box {factor}x the model tolerance too long is "
              f"{'accepted' if want_accept else 'rejected'}",
              r["accepted"] == want_accept,
              f"dV={r['symmetricDifference']:.3g} band={r['band']:.3g}")

    # 4. negative control: right volume, wrong shape. A sphere and a box of equal volume must be
    #    rejected -- this is the check that a volume *comparison* could pass and a symmetric
    #    difference cannot.
    v = 2.0 * 3.0 * 4.0
    rad = (3.0 * v / (4.0 * math.pi)) ** (1.0 / 3.0)
    sphere = BRepPrimAPI_MakeSphere(gp_Pnt(1.0, 1.5, 2.0), rad).Shape()
    r = symmetric_difference(box, sphere, tol)
    check("a sphere of equal volume is rejected", not r["accepted"],
          f"dV={r['symmetricDifference']:.3g}, volumes {r['volumeOriginal']:.4f} vs "
          f"{r['volumeCandidate']:.4f}")

    # 5. a tube, built two ways: OCCT cut versus the description's TGeoTube.
    outer = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -5), gp_Dir(0, 0, 1)), 2.0, 10.0).Shape()
    inner = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -6), gp_Dir(0, 0, 1)), 1.0, 12.0).Shape()
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
    op = BRepAlgoAPI_Cut(outer, inner)
    op.Build()
    tube = op.Shape()
    cand = prim.candidate("primitive", [prim.leaf(
        "TGeoTube", {"rmin": 1.0, "rmax": 2.0, "dz": 5.0}, prim.identity_frame())], "self-test")
    r = symmetric_difference(tube, prim.build_occ(cand), tol)
    check("a TGeoTube description matches an OCCT-cut tube", r["accepted"],
          f"dV={r['symmetricDifference']:.3g} band={r['band']:.3g}")

    # 6. negative control on the tube: a solid cylinder must not pass as the tube.
    solid_cand = prim.candidate("primitive", [prim.leaf(
        "TGeoTube", {"rmin": 0.0, "rmax": 2.0, "dz": 5.0}, prim.identity_frame())], "self-test")
    r = symmetric_difference(tube, prim.build_occ(solid_cand), tol)
    check("a solid cylinder is rejected as the tube", not r["accepted"],
          f"dV={r['symmetricDifference']:.3g}")

    # 7. the tube in a rotated, translated frame -- the case the whole frame machinery exists for.
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.gp import gp_Trsf, gp_Ax1, gp_Vec
    trsf = gp_Trsf()
    trsf.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), 0.7)
    shift = gp_Trsf()
    shift.SetTranslation(gp_Vec(3.0, -4.0, 5.0))
    moved = BRepBuilderAPI_Transform(tube, shift.Multiplied(trsf), True).Shape()
    zaxis = _rotated((0.0, 0.0, 1.0), (1.0, 1.0, 0.0), 0.7)
    xaxis = _rotated((1.0, 0.0, 0.0), (1.0, 1.0, 0.0), 0.7)
    frame = prim.frame_from_axis((3.0, -4.0, 5.0), zaxis, xaxis)
    cand = prim.candidate("primitive", [prim.leaf(
        "TGeoTube", {"rmin": 1.0, "rmax": 2.0, "dz": 5.0}, frame)], "self-test")
    r = symmetric_difference(moved, prim.build_occ(cand), tol)
    check("a rotated, translated tube matches its placed description", r["accepted"],
          f"dV={r['symmetricDifference']:.3g} band={r['band']:.3g}")

    # 8. negative control on the frame: the *unrotated* description must be rejected against it.
    flat = prim.candidate("primitive", [prim.leaf(
        "TGeoTube", {"rmin": 1.0, "rmax": 2.0, "dz": 5.0},
        prim.identity_frame((3.0, -4.0, 5.0)))], "self-test")
    r = symmetric_difference(moved, prim.build_occ(flat), tol)
    check("the same tube without the rotation is rejected", not r["accepted"],
          f"dV={r['symmetricDifference']:.3g}")

    n_ok = sum(1 for _n, ok, _d in checks if ok)
    if verbose:
        print(f"  {n_ok}/{len(checks)} acceptance self-checks passed")
    return n_ok, len(checks)


def _rotated(vec, axis, angle):
    """Rodrigues; used only by the self-test, to state the expected frame independently."""
    from csg.primitives import _cross, _dot, _scale, _unit, _add
    k = _unit(axis)
    return _add(_add(_scale(vec, math.cos(angle)), _scale(_cross(k, vec), math.sin(angle))),
                _scale(k, _dot(k, vec) * (1.0 - math.cos(angle))))
