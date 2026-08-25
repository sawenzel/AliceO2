"""The flat-DNF emitter: a cell's carriers as signed implicit halfspaces.

`recognise._halfspace_carriers` already produces exactly what this module needs -- an oriented
carrier per face with the material side decided -- and `recognise._cell_leaf` then *bounds* each
one into a padded native primitive, because ROOT has no usable halfspace and a
`TGeoCompositeShape` is the only thing the converter could emit. `O2FlatCSG` removes that
constraint, so this module maps a carrier to the halfspace it always was.

The sign is the whole risk
--------------------------
An inverted halfspace still produces a solid, so a sign error here is silent. Two things compose
into the single sign: the carrier's own orientation (a plane's normal is already flipped for
`TopAbs_REVERSED` upstream) and the `side` field. `csg/emit.py --self-test` measures the result
against `_cell_leaf`'s padded conjunction on sampled points, because `_cell_leaf` has shipped 1775
known-source-clean parts and is the strongest oracle available for it.

The plane scaling is an obligation, not a taste
-----------------------------------------------
A plane is stored with `2b = n` for a *unit* outward normal `n`. Any positive rescaling describes
the same halfspace, so no geometry test would notice a different one -- what it costs is bit
identity between `O2FlatCSG`'s accelerated queries and their `_Loop` twins, which is that class's
whole self-check discipline (`Design_FlatCSGSolid.md` section 3.1). This module is the emitter that
obligation is addressed to, so `emit.py --self-test` asserts `|2b| == 1` on every plane block it
produces -- at the place the convention is created, since the C++ side cannot tell an
emitter-produced halfspace from a hand-built one.
"""

import math
import struct

SIDECAR_MAGIC = b"O2FLTCSG"
SIDECAR_VERSION = 1

# The two `FlatCSGHalfspace::Kind` values, as the sidecar spells them.
KIND_QUADRIC = 0
KIND_TORUS = 1

# One quadric block is ten doubles; the record on file carries eleven, the last unused.
QUADRIC_COEFFICIENTS = 10
BLOCK_COEFFICIENTS = 11


def _outer(u, v):
    return [[u[i] * v[j] for j in range(3)] for i in range(3)]


def _quadric(a, b, c):
    """Pack A (3x3 symmetric), b (3) and c into the ten-double block the shape stores."""
    return [a[0][0], a[0][1], a[0][2], a[1][1], a[1][2], a[2][2], b[0], b[1], b[2], c]


def quadric_from_carrier(carrier):
    """`(sign, block)` for a plane, sphere, cylinder or cone carrier.

    The material side is `sign * Q(x) <= 0`. `side == "exterior"` means the material is on the
    far side of the carrier from its own inside, which is exactly a flipped sign.
    """
    from csg import recognise
    sign = -1.0 if carrier["side"] == "exterior" else 1.0
    kind = carrier["kind"]

    if kind == "plane":
        n = carrier["n"]
        p = carrier["p"]
        # Q(x) = n.(x - p); the material side of an outward normal is Q <= 0. `2b = n` for a unit
        # normal is the convention of design section 3.1 and is not free to vary.
        return sign, _quadric([[0.0] * 3 for _ in range(3)],
                              [0.5 * n[0], 0.5 * n[1], 0.5 * n[2]],
                              -(n[0] * p[0] + n[1] * p[1] + n[2] * p[2]))

    if kind == "sphere":
        p = carrier["p"]
        r = carrier["r"]
        identity = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        return sign, _quadric(identity, [-p[0], -p[1], -p[2]],
                              p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - r * r)

    if kind in ("cylinder", "cone"):
        d = carrier["d"]
        p = carrier["p"]
        r = carrier["r"]
        k = 0.0 if kind == "cylinder" else math.tan(carrier["a"])
        scale = 1.0 + k * k
        dd = _outer(d, d)
        a = [[(1.0 if i == j else 0.0) - scale * dd[i][j] for j in range(3)] for i in range(3)]
        ap = [sum(a[i][j] * p[j] for j in range(3)) for i in range(3)]
        pd = sum(p[i] * d[i] for i in range(3))
        b = [-ap[i] - r * k * d[i] for i in range(3)]
        c = sum(p[i] * ap[i] for i in range(3)) + 2.0 * r * k * pd - r * r
        return sign, _quadric(a, b, c)

    raise recognise.Declined(f"a {kind} carrier has no quadric form")


def torus_from_carrier(carrier):
    """`(sign, centre, axis, major, minor)` for a torus carrier."""
    sign = -1.0 if carrier["side"] == "exterior" else 1.0
    return sign, list(carrier["p"]), list(carrier["d"]), carrier["r"], carrier["rt"]


def blocks_from_carriers(carriers):
    """One halfspace block per carrier, in the carriers' own order."""
    blocks = []
    for carrier in carriers:
        if carrier["kind"] == "torus":
            sign, centre, axis, major, minor = torus_from_carrier(carrier)
            blocks.append({"kind": "torus", "sign": sign,
                           "c": centre + axis + [major, minor, 0.0, 0.0, 0.0]})
        else:
            sign, block = quadric_from_carrier(carrier)
            blocks.append({"kind": "quadric", "sign": sign, "c": block + [0.0]})
    return blocks


def eval_block(block, point):
    """`sign * f(point)`; the halfspace contains the point when this is `<= 0`.

    The same arithmetic as `O2FlatCSG::EvalHalfspace`, and it exists so the self-test can
    classify a point without a built shape. The C++ side is the authority; if the two ever
    disagree the test that finds it is `emit.py --self-test`.
    """
    c = block["c"]
    x, y, z = point
    if block["kind"] == "torus":
        offset = (x - c[0], y - c[1], z - c[2])
        along = offset[0] * c[3] + offset[1] * c[4] + offset[2] * c[5]
        radial = tuple(offset[i] - along * c[3 + i] for i in range(3))
        rho = math.sqrt(sum(v * v for v in radial))
        return block["sign"] * (math.hypot(rho - c[6], along) - c[7])
    quadratic = (c[0] * x * x + c[3] * y * y + c[5] * z * z +
                 2.0 * (c[1] * x * y + c[2] * x * z + c[4] * y * z))
    return block["sign"] * (quadratic + 2.0 * (c[6] * x + c[7] * y + c[8] * z) + c[9])


def flat_contains(blocks, point):
    """True when every block contains the point: one cell's membership test."""
    return all(eval_block(block, point) <= 0.0 for block in blocks)


def plane_scaling_error(block):
    """`| |2b| - 1 |` for a plane block, or None when the block is not a plane.

    A plane is the quadric whose `A` vanishes identically. Design section 3.1 obliges this module
    to store it with `2b = n` for a unit `n`, so this is the residual the self-test asserts on.
    """
    if block["kind"] != "quadric":
        return None
    c = block["c"]
    if any(c[index] != 0.0 for index in range(6)):
        return None
    two_b = math.sqrt(4.0 * (c[6] * c[6] + c[7] * c[7] + c[8] * c[8]))
    return abs(two_b - 1.0)


def write_sidecar(path, blocks, cells):
    """Write the version-1 flat-CSG sidecar. Must stay byte-compatible with `WriteFlatCSG`.

    Field by field, never as a struct: the halfspace record's natural C++ layout pads to 104
    bytes while the file packs it at 100 (`BVHSurfaceSolid.md`, "Flat-CSG sidecar format").
    """
    with open(path, "wb") as handle:
        handle.write(SIDECAR_MAGIC)
        handle.write(struct.pack("<III", SIDECAR_VERSION, len(blocks), len(cells)))
        for block in blocks:
            handle.write(struct.pack("<i", KIND_TORUS if block["kind"] == "torus"
                                     else KIND_QUADRIC))
            handle.write(struct.pack("<d", block["sign"]))
            coefficients = list(block["c"])
            if len(coefficients) > BLOCK_COEFFICIENTS:
                raise ValueError(f"a halfspace block carries {len(coefficients)} coefficients, "
                                 f"more than the {BLOCK_COEFFICIENTS} the sidecar has room for")
            coefficients += [0.0] * (BLOCK_COEFFICIENTS - len(coefficients))
            handle.write(struct.pack("<11d", *coefficients))
        for cell in cells:
            handle.write(struct.pack("<ii", cell["first"], cell["count"]))
            handle.write(struct.pack("<d", cell["volume"]))
            handle.write(struct.pack("<3d", *cell["lo"]))
            handle.write(struct.pack("<3d", *cell["hi"]))
