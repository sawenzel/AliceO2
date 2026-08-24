"""Tier 0: the plane, cylinder, cone, sphere or torus a stored B-spline face already is.

A CAD exporter routinely writes an exact quadric as a rational B-spline patch, and the stored
surface type then describes the exporter rather than the geometry. On ALICE3's `CAD.stp` **998 of
the 2450 faces whose stored surface is not analytic are exactly a quadric** -- 786 cylinders, 176
cones, 36 spheres -- at a worst measured gap of 3.8e-10 cm. Three instruments already agreed on
that population (`Stream_A_CSG.md` §2.3, `Stream_K_Tier0.md` §3.1, OCCT's own
`ShapeAnalysis_CanonicalRecognition`); what did not exist is those faces **participating** as
carriers. Until this module, `recognise._face_records` counted every one of them as free-form, so
a part carrying a single such face declined CSG outright whatever its structure.

This service decides, per face, whether the face IS a canonical carrier, and returns it in the
same vocabulary a natively-analytic face gets, so the cascade above it needs no second code path.

Propose cheap, accept measured
------------------------------
The proposals for plane / sphere / cylinder / cone are `O2_CADtoTGeo._analytic_surface_proposals`,
**reused rather than rewritten**: those are the shipping recogniser's own linear solves, they need
no initial guess, and its self-test holds them against closed-form controls. (Splitting that
generator out of `_recognize_analytic_surface` is the one edit this rung made to the converter,
and it is why: the converter accepts a face for sidecar emission against the *patch's* diagonal at
1e-9, this module accepts a carrier against the *part's* at 1e-6, and a service whose documented
criterion is silently overridden by a different one upstream is exactly the defect
`Stream_K_Tier0.md` §3 is about. The converter's own behaviour is unchanged, byte for byte.)
A **torus** proposal is added here, because the converter deliberately carries none -- its
self-test asserts the absence, and this module does not change that -- as one global linear solve
for the axis, polished by Gauss-Newton on the gap itself.

The acceptance is ONE measured quantity, and it is not the proposer's:

    gap = the largest distance, in cm, from a sampled face point to the candidate surface
    admissible  <=>  gap <= REL_TOL * max(part bounding-box diagonal, 1 cm)

and among the proposals that are admissible the fewest-parameter one is the answer -- the rule
`_recognize_analytic_surface` already follows for its plane, extended to the whole ladder. Nothing
but the gap is ever compared against anything.

`Stream_K_Tier0.md` §3 is why that sentence is written this way. The recogniser used to score four
candidates with four different expressions -- two distances, two angles -- against one tolerance,
and 184 ALICE3 faces passed as cones that miss their own surface by up to 79 cm. An acceptance
criterion that is not a measured criterion has burned this exact path once. So: no angles, no
per-class thresholds, one number, and `csg/emit.py --self-test` checks the instrument itself
against a known displacement at its true size before it checks any face.

The samples the acceptance measures on are **not** the samples the proposal was fitted to: the
proposal comes off the converter's 9x9 grid and the gap is measured on an independent 17x17 one.
A model that fits its own fitting points and nothing else fails here.

What the gap does and does not bound
------------------------------------
The gap bounds the candidate's error **on the patch**. A carrier is used well beyond its patch --
`recognise._cell_leaf` extends it across the whole part -- and no measurement made on the patch
bounds the extension. This module is therefore a *proposal service*: what admits a part is still
the part-level symmetric difference in `csg/accept.py` and the oracle gate, exactly as for a
natively-analytic carrier. A canonicalised face that is subtly wrong costs coverage there, not
correctness.

One more thing this module does not do: it says nothing about the face's **trim**. The canonical
`(u, v)` box it reports is the trim's parametric bounding box in the canonical chart, which is
all the cluster and wedge matchers read; the trim curves themselves stay in the B-spline's own
chart and are `Stream_K_Tier0.md` §2's separate, unbuilt work item.
"""

import math

# The same band, and the same reasoning, as `recognise.REL_TOL`: relative to the part's own
# bounding-box diagonal, because CAD arrives with ~1e-7 relative agreement between faces the
# engineer meant to be identical. `csg/emit.py --self-test` asserts the two have not drifted
# apart -- the cascade cannot have one opinion of "the same" for its carriers and another for
# the faces they came from.
REL_TOL = 1.0e-6

# The proposal grid (the converter's own) and the independent, denser acceptance grid.
_PROPOSE_N = 9
_ACCEPT_N = 17


class _Unavailable(Exception):
    """The converter module could not be imported, so nothing here can run."""


_CONVERTER = None


def _converter():
    """`O2_CADtoTGeo`, imported lazily and kept.

    It is a 4000-line module that pulls in the whole of pythonOCC, and the common case -- a part
    whose faces are all natively analytic -- never needs it. When the converter itself is the
    caller the module is already loaded and this costs a dictionary lookup.
    """
    global _CONVERTER
    if _CONVERTER is None:
        import sys
        from pathlib import Path
        root = str(Path(__file__).resolve().parent.parent)
        if root not in sys.path:
            sys.path.insert(0, root)
        try:
            import O2_CADtoTGeo
        except Exception as exc:                                 # noqa: BLE001
            raise _Unavailable(str(exc)) from None
        _CONVERTER = O2_CADtoTGeo
    return _CONVERTER


# ------------------------------------------------------------------------------------------
# the instrument
# ------------------------------------------------------------------------------------------

def surface_gap(kind, model, points):
    """The largest distance, in cm, from any of `points` to the candidate surface.

    The one quantity this module decides on. For plane / sphere / cylinder / cone it IS
    `O2_CADtoTGeo._analytic_surface_gap`, so the converter's surface path and this one cannot
    drift apart in what they mean by "how far off is it".
    """
    if kind == "torus":
        return _torus_gap(points, model)
    return _converter()._analytic_surface_gap(kind, model, points)


def _torus_residual(points, centre, axis, major, minor):
    import numpy as np
    h = (points - centre) @ axis
    rho = np.linalg.norm(points - centre - np.outer(h, axis), axis=1)
    return np.sqrt((rho - major) ** 2 + h ** 2) - minor


def _torus_gap(points, model):
    import numpy as np
    return float(np.abs(_torus_residual(points, model["centre"], model["axis"],
                                        model["major"], model["minor"])).max())


# ------------------------------------------------------------------------------------------
# the torus proposal (the one model the converter's recogniser does not carry)
# ------------------------------------------------------------------------------------------

def _torus_radii(points, centre, axis):
    """`(R, r)` by least squares once the axis is fixed, or None if the solve is not a torus.

    With `h` the axial and `rho` the radial coordinate about the axis, a torus is
    `(rho - R)^2 + h^2 = r^2`, i.e. `rho^2 + h^2 = 2R rho + (r^2 - R^2)` -- linear in the two
    unknowns `2R` and `r^2 - R^2`.
    """
    import numpy as np
    h = (points - centre) @ axis
    rho = np.linalg.norm(points - centre - np.outer(h, axis), axis=1)
    design = np.column_stack([rho, np.ones_like(rho)])
    sol, *_ = np.linalg.lstsq(design, rho ** 2 + h ** 2, rcond=None)
    major = 0.5 * float(sol[0])
    minor_sq = float(sol[1]) + major * major
    if not (major > 0.0 and minor_sq > 0.0):
        return None
    return major, math.sqrt(minor_sq)


def _propose_torus(points, normals, refinements=25):
    """`{axis, centre, major, minor}` for the torus these samples propose, or None.

    Every point of a torus lies in a meridian plane that also contains the axis direction `d` and
    the point's own normal, so `det[d, N, P - c] = 0` for every sample. Written out, with
    `g = c x d`, that is

        (N_i x P_i) . d  +  N_i . g  =  0,

    which is **linear and homogeneous in the six unknowns `(d, g)` jointly**, so one SVD answers
    it globally with no initial guess and no iteration; `c` comes back as `(d x g) / |d|^2`, its
    free component along the axis being zero by construction. Solving instead for `d` and `c`
    alternately -- the obvious reading of the same equation -- converges to the wrong axis on a
    torus wedge, because for a narrow range of azimuths the meridian normals no longer pin it;
    measured, that reading missed a 90-degree wedge of a `R = 6, r = 1.5` torus completely, and
    this one recovers both radii to 1e-14 cm from 30 degrees of sweep upwards.

    A cylinder has all its normals in one plane, which makes `(d = 0, g = axis)` an exact
    solution of the system; that degenerate answer is caught by the `|d|` test and declined, and
    everything else is left in and simply scored by the gap -- the same discipline
    `_recognize_analytic_surface` follows for its runaway cone apex.

    The Gauss-Newton polish that follows works on the gap itself with a numeric Jacobian in the
    seven parameters. It exits the moment it stops improving, so on the exact input above it
    costs one iteration and changes nothing.
    """
    import numpy as np

    design = np.column_stack([np.cross(normals, points), normals])
    _, _singular, right = np.linalg.svd(design, full_matrices=False)
    solution = right[-1]
    axis, moment = solution[:3], solution[3:]
    length = float(np.linalg.norm(axis))
    if length < 1.0e-6:
        return None                     # d = 0: coplanar normals, i.e. a cylinder, not a torus
    axis = axis / length
    centre = np.cross(axis, moment / length)

    radii = _torus_radii(points, centre, axis)
    if radii is None:
        return None
    major, minor = radii
    span = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0))) or 1.0
    step = 1.0e-7 * span
    for _ in range(refinements):
        tangent_a = np.cross(axis, [1.0, 0.0, 0.0])
        if np.linalg.norm(tangent_a) < 1e-6:
            tangent_a = np.cross(axis, [0.0, 1.0, 0.0])
        tangent_a = tangent_a / np.linalg.norm(tangent_a)
        tangent_b = np.cross(axis, tangent_a)
        base = _torus_residual(points, centre, axis, major, minor)

        def at(delta):
            tilted = axis + delta[3] * tangent_a + delta[4] * tangent_b
            tilted = tilted / np.linalg.norm(tilted)
            return _torus_residual(points, centre + delta[:3], tilted,
                                   major + delta[5], minor + delta[6])

        jacobian = np.zeros((len(points), 7))
        for column in range(7):
            probe = np.zeros(7)
            probe[column] = step
            jacobian[:, column] = (at(probe) - base) / step
        try:
            delta, *_ = np.linalg.lstsq(jacobian, -base, rcond=None)
        except np.linalg.LinAlgError:
            break
        if np.abs(at(delta)).max() >= np.abs(base).max():
            break                                    # no longer improving: keep what converged
        centre = centre + delta[:3]
        axis = axis + delta[3] * tangent_a + delta[4] * tangent_b
        axis = axis / np.linalg.norm(axis)
        major += float(delta[5])
        minor += float(delta[6])
    if not (major > 0.0 and minor > 0.0):
        return None
    return {"axis": axis, "centre": centre, "major": float(major), "minor": float(minor)}


# ------------------------------------------------------------------------------------------
# the service
# ------------------------------------------------------------------------------------------

def canonicalise(face, adaptor, scale):
    """`(carrier, gap)`: the canonical carrier this face IS, and the gap that decided.

    `carrier` is None when the face is not canonical; `gap` is then still the smallest distance
    any of the models PROPOSED for it achieves, so a decline can say how far from canonical the
    face is instead of only that it is. That is not the distance to the nearest quadric in the
    world -- nothing here computes that -- and the wording of every message built from it says
    so. Both are None where the face cannot be sampled at all.

    `scale` is `max(part bounding-box diagonal, 1 cm)` -- the acceptance denominator, so that one
    part's faces are all judged against one length and a small patch on a big part is not held to
    a tighter standard than the part it belongs to.

    The returned record speaks `recognise._face_records`' vocabulary (`kind` plus the carrier
    parameters) and adds three fields of its own: `canonicalised`, `tier0GapCm` and
    `tier0GapRelative`. `uv` is the trim's bounding box in the CANONICAL chart for the three kinds
    that read it (cylinder, cone and torus, in `_axial_extent` and `_phi_range`) and is
    `None` for a plane and a sphere, whose canonical chart nothing reads and which carry no
    reference direction of their own to state one against, so a future reader that needs one
    fails loudly instead of silently taking a B-spline's own parameters for an angle.
    """
    try:
        conv = _converter()
    except _Unavailable:
        return None, None
    from OCC.Core.BRepTools import breptools

    try:
        uv_bounds = breptools.UVBounds(face)
    except Exception:                                            # noqa: BLE001
        return None, None
    propose_points, propose_normals = conv._sample_surface_for_recognition(
        adaptor, *uv_bounds, n=_PROPOSE_N)
    if propose_points is None:
        return None, None
    accept_points, _accept_normals = conv._sample_surface_for_recognition(
        adaptor, *uv_bounds, n=_ACCEPT_N)
    if accept_points is None:
        return None, None

    # The candidates, in order of parsimony extended by one rung: plane (3 parameters) <
    # sphere (4) < cylinder (5) < cone (6) < torus (7). The first four come from the converter's
    # own solves; the torus is this module's.
    proposals = list(conv._analytic_surface_proposals(propose_points, propose_normals))
    torus = _propose_torus(propose_points, propose_normals)
    if torus is not None:
        proposals.append(("torus", torus))

    # **The one measured quantity decides whether a proposal is admissible at all**; among the
    # proposals it admits, the fewest-parameter one wins, which is the rule
    # `_recognize_analytic_surface` already follows for its plane. Taking the smallest gap
    # instead is worse and measurably so: a sphere is a torus of zero major radius, so the torus
    # proposal fits a sphere at 1e-16 -- fractionally better than the sphere itself -- and a
    # NURBS-encoded sphere then arrives as a self-intersecting `TGeoTorus` that cannot be stated.
    kind, model, gap, best_gap = None, None, None, float("inf")
    for candidate_kind, candidate in proposals:
        try:
            candidate_gap = surface_gap(candidate_kind, candidate, accept_points)
        except Exception:                                        # noqa: BLE001
            continue
        if not math.isfinite(candidate_gap):
            continue
        best_gap = min(best_gap, candidate_gap)
        if kind is None and candidate_gap <= REL_TOL * scale:
            kind, model, gap = candidate_kind, candidate, candidate_gap

    if kind is None:
        return None, (None if not math.isfinite(best_gap) else best_gap)
    record = _carrier_record(kind, model, adaptor, uv_bounds)
    if record is None:
        return None, gap
    record["canonicalised"] = True
    record["tier0GapCm"] = gap
    record["tier0GapRelative"] = gap / scale
    return record, gap


def carrier_side(face, adaptor, carrier):
    """`interior` / `exterior` for a canonicalised face, by `census`'s one rule."""
    from csg import census
    return census.halfspace_side_of(face, adaptor, carrier)


def _carrier_record(kind, model, adaptor, uv_bounds):
    """The canonical carrier as `recognise._face_records` states one."""
    import numpy as np
    if kind == "plane":
        normal = np.asarray(model["normal"], dtype=float)
        normal = normal / np.linalg.norm(normal)
        # Unflipped, i.e. the underlying surface's own normal: both callers apply the face's
        # REVERSED flag themselves, exactly as they do for a native plane.
        return {"kind": "plane", "n": tuple(float(c) for c in normal),
                "p": tuple(float(c) for c in model["point"]), "uv": None}
    if kind == "sphere":
        return {"kind": "sphere", "p": tuple(float(c) for c in model["centre"]),
                "r": float(model["radius"]), "uv": None}
    if kind == "torus":
        axis = _unit_array(model["axis"])
        # A fitted torus brings no reference direction of its own, so one is chosen here and the
        # chart below is measured against that same one -- the two cannot disagree.
        ref = np.asarray(_perpendicular_to(axis), dtype=float)
        chart = _canonical_chart(adaptor, uv_bounds, np.asarray(model["centre"], dtype=float),
                                 axis, ref, semi_angle=None, major=float(model["major"]))
        if chart is None:
            return None
        return {"kind": "torus", "d": tuple(float(c) for c in axis),
                "p": tuple(float(c) for c in model["centre"]),
                "x": tuple(float(c) for c in ref), "r": float(model["major"]),
                "rt": float(model["minor"]), "uv": chart}
    if kind == "cylinder":
        axis = _unit_array(model["axis"])
        origin = np.asarray(model["origin"], dtype=float)
        ref = _orthonormalise(np.asarray(model["refu"], dtype=float), axis)
        if ref is None:
            return None
        chart = _canonical_chart(adaptor, uv_bounds, origin, axis, ref, semi_angle=None)
        if chart is None:
            return None
        return {"kind": "cylinder", "d": tuple(float(c) for c in axis),
                "p": tuple(float(c) for c in origin), "x": tuple(float(c) for c in ref),
                "r": float(model["radius"]), "uv": chart}
    if kind == "cone":
        axis = _unit_array(model["axis"])
        apex = np.asarray(model["apex"], dtype=float)
        ref = _orthonormalise(np.asarray(model["refu"], dtype=float), axis)
        if ref is None:
            return None
        half = float(model["half_angle"])
        if not (1.0e-9 < half < 0.5 * math.pi - 1.0e-9):
            return None
        chart = _canonical_chart(adaptor, uv_bounds, apex, axis, ref, semi_angle=half)
        if chart is None:
            return None
        # Stated at the apex, where the reference radius is zero -- the chart `_axial_extent`
        # and `_cell_leaf` read is `r = RefRadius + v sin(a)`, `t = v cos(a)`, which is OCC's
        # own `gp_Cone` parametrisation.
        return {"kind": "cone", "d": tuple(float(c) for c in axis),
                "p": tuple(float(c) for c in apex), "x": tuple(float(c) for c in ref),
                "r": 0.0, "a": half, "uv": chart}
    return None


def _unit_array(vec):
    import numpy as np
    v = np.asarray(vec, dtype=float)
    return v / np.linalg.norm(v)


def _orthonormalise(vec, axis):
    import numpy as np
    ref = np.asarray(vec, dtype=float)
    ref = ref - float(ref @ axis) * axis
    length = float(np.linalg.norm(ref))
    if length < 1.0e-9:
        return None
    return ref / length


def _perpendicular_to(axis):
    import numpy as np
    seed = np.array([1.0, 0.0, 0.0]) if abs(float(axis[0])) < 0.9 else np.array([0.0, 1.0, 0.0])
    ref = _orthonormalise(seed, axis)
    return tuple(float(c) for c in ref)


_CHART_N = 33


def _canonical_chart(adaptor, uv_bounds, origin, axis, ref, semi_angle, major=None):
    """`(umin, umax, vmin, vmax)`: the trim's bounding box in the carrier's OWN chart.

    This is the one thing that cannot be taken from the face: a B-spline's `(u, v)` is some
    reparametrisation of the quadric's, so `breptools.UVBounds` of a canonicalised face says
    nothing about an angle or an axial distance, and `_axial_extent` and `_phi_range` read both
    as if it did.

    Measured along the two midlines of the face's parametric rectangle, because that is where the
    azimuth is a monotone curve and can be unwrapped: accumulating the 2pi jumps is what makes a
    patch straddling `atan2`'s branch cut, and a face that closes on itself at exactly a full
    turn, come out right, where taking the min and max of the raw angles would report a full turn
    for the first and a short arc for the second. Both midlines pass through the patch centre, so
    they are anchored to one branch there and their ranges can then simply be unioned.
    """
    import numpy as np
    umin, umax, vmin, vmax = uv_bounds
    umid, vmid = 0.5 * (umin + umax), 0.5 * (vmin + vmax)
    binormal = np.cross(axis, ref)

    def chart_of(u, v):
        try:
            point = adaptor.Value(u, v)
        except Exception:                                        # noqa: BLE001
            return None
        rel = np.array([point.X(), point.Y(), point.Z()]) - origin
        axial = float(rel @ axis)
        perp = rel - axial * axis
        if float(np.linalg.norm(perp)) < 1.0e-30:
            return None
        azimuth = math.atan2(float(perp @ binormal), float(perp @ ref))
        if major is not None:                                    # a torus: the meridian angle
            return azimuth, math.atan2(axial, float(np.linalg.norm(perp)) - major)
        return azimuth, axial if semi_angle is None else axial / math.cos(semi_angle)

    anchor = chart_of(umid, vmid)
    if anchor is None:
        return None
    phis, axials = [anchor[0]], [anchor[1]]
    for fixed, lo, hi, along_u in ((vmid, umin, umax, True), (umid, vmin, vmax, False)):
        samples, centre_index = [], None
        for k in range(_CHART_N):
            t = lo + (hi - lo) * k / (_CHART_N - 1.0)
            got = chart_of(t, fixed) if along_u else chart_of(fixed, t)
            if got is None:
                continue
            if centre_index is None and t >= 0.5 * (lo + hi):
                centre_index = len(samples)
            samples.append(got)
        if len(samples) < 2 or centre_index is None:
            continue
        unwrapped = np.unwrap(np.array([s[0] for s in samples]))
        unwrapped += 2.0 * math.pi * round((anchor[0] - unwrapped[centre_index])
                                           / (2.0 * math.pi))
        phis.extend(float(p) for p in unwrapped)
        axials.extend(s[1] for s in samples)
    return (min(phis), max(phis), min(axials), max(axials))
