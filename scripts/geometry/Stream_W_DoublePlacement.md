# Stream W — ALICE3 built the whole detector twice, and now it does not

Date: 2026-08-03. Acts on the verdict in [`Stream_T_AssemblyOracle.md`](Stream_T_AssemblyOracle.md) §0,
which found the defect and left the fix open. Companion to [`Tutorial.md`](Tutorial.md) (the map).

---

## 0. The verdict

**ALICE3's converted world no longer contains a single coincident duplicate.**
`CAD_noETA.stp` declares 206 leaf solid placements holding only 103 distinct
(definition, world transform) pairs — every solid twice, at an *exactly* identical matrix. The
converter now emits **103**, multiplicity `{1: 103}`, and refuses to emit anything else.

Bagger, the control, is untouched: **13 in, 13 out**, and every gate number on it is unchanged.

---

## 1. Attribution: the STEP file declares it. This is a modelling defect, not a traversal bug.

Read out of the file's own entities rather than inferred, with a nine-line text parser over the
`NEXT_ASSEMBLY_USAGE_OCCURRENCE` (NAUO) graph:

* There is exactly **one** root product definition, `#5 = PRODUCT_DEFINITION('design','',#6,#9)`,
  product name `Unnamed`, reached from the root `#3 = SHAPE_DEFINITION_REPRESENTATION(#4,#10)`.
  It has **no** parent in any NAUO. Its shape representation is
  `#10 = SHAPE_REPRESENTATION('',(#11,#15,#19,#23,#27),#31)` — an origin plus **four** child
  placements.
* Those four children are, entity for entity:

  | NAUO | occurrence id | child | also a child of |
  | --- | ---: | --- | --- |
  | `#66` | `'120'` | `ST2487728_01` — the whole detector, 103 leaves | — |
  | `#452894` | `'121'` | `ST2487730_01`, 1 leaf | **`#39` = `ST2487728_01`** |
  | `#452899` | `'122'` | `ST1A38517_01`, 74 leaves | **`#39` = `ST2487728_01`** |
  | `#452904` | `'123'` | `ST0923290_058`, 28 leaves | **`#39` = `ST2487728_01`** |

  1 + 74 + 28 = 103. The root lists the complete detector **and, as its siblings, that detector's
  own three children.**
* All four placements are the identity: `#15`, `#19`, `#23`, `#27` are each
  `AXIS2_PLACEMENT_3D` at `(0,0,0)` with axis `(0,0,1)` and ref direction `(1,0,0)`, and each of
  the three duplicate NAUOs carries an `ITEM_DEFINED_TRANSFORMATION('','',#11,#19|#23|#27)`.
  So both copies land on top of each other, 100 % coincident.
* The three duplicate NAUOs are **`#452894`, `#452899` and `#452904` — the last three assembly
  entities in a 454 632-entity file**, carrying occurrence ids `'121'`, `'122'`, `'123'`
  immediately after the legitimate `'120'`. They were appended after the model was otherwise
  complete.

A *plain*, non-XCAF read agrees, which is the independent confirmation that the structure and not
our importer carries it — `STEPControl_Reader` + `OneShape()`, solids explored in world
coordinates, grouped by (volume, world bounding box) rounded to 1e-9:

```
CAD_noETA.stp               : 206 placed solids | 103 distinct signatures | multiplicity {2: 103}
Bagger.step                 :  13 placed solids |  13 distinct signatures | multiplicity {1: 13}
as1-oc-214.stp              :  18 placed solids |  18 distinct signatures | multiplicity {1: 18}
oTOF System V3-R92cm.step   : 62628 placed solids | 62628 distinct signatures | multiplicity {1: 62628}
```

**Therefore the converter is compensating for a modelling defect, not fixing a traversal bug**, and
the log message says exactly that:

> `WARNING: this CAD model DECLARES 103 leaf solid placement(s) that coincide exactly with another placement of the same solid.`
> `  The assembly structure in the file says so -- these are not an artefact of this traversal, which walks the STEP product structure edge for edge.`

That distinction matters downstream: the CAD should be fixed at source, and until it is, every
model from that exporter is suspect. It does **not** change what the converter must do — a
geometry with two volumes in the same place has no defined transport through it either way.

---

## 2. The rule: key on (definition, world transform), never on the definition

**The single thing that could have made this fix worse than the bug** is a rule keyed on the
*definition*. Two placements of one definition at *different* transforms are what instancing IS —
ALICE3 places `ST2487458_01` twelve times — and collapsing those would have silently deleted most
of the detector while making the count look better. Every rule below keys on the **pair**.

`deduplicate_placements()` in `O2_CADtoTGeo.py`, applied once, after the traversal:

1. **Root containment** (`root-containment`) — the structural rule, and the one that names the
   defect. A child of a top-level definition that a *sibling* child of the same top-level
   definition already places, at the same world transform, is not placed at top level a second
   time. Safe without further argument: a top-level definition occurs exactly once, so dropping
   one of its edges cannot remove a placement from some other instance of the same parent.

2. **Coincident-occurrence de-duplication** (`coincident-occurrence`) — the defensive rule, for
   whatever the structural one misses: one sub-assembly reached at the same place down two
   different assembly paths, or a third instance that coincides with the second. A placement edge
   *every* one of whose occurrences lands on a (definition, world transform) pair another edge
   already discovered contributes nothing, and is dropped. Run to a fixed point.

   **"Every one of whose occurrences" is load-bearing.** `placements` is keyed by *definition*, so
   an edge that is redundant under one instance of its parent and needed under another cannot be
   dropped — that would delete the needed one too. Such edges are **reported and left alone**;
   suppressing them would need the parent definition split per instance, which is a larger change
   than this defect justifies. No corpus in the repository has one.

**On ALICE3 it is rule 1 that fires, and only rule 1**: 3 edges, `3 by root-containment, 0 by
coincident-occurrence`, dropping 103 duplicate leaf placements. Rule 1 running first is not
cosmetic — it decides *which* of the two coincident copies survives. Rule 2 alone would, depending
on the order XCAF hands back the root's children, have been just as entitled to gut
`ST2487728_01` and leave its three children hanging off the root: the same 103 leaves at the same
103 matrices, but with the assembly structure destroyed.

Transforms are compared through `trsf_signature()` — the 12 matrix entries rounded to 1e-9. The
rounding can only ever cost a *missed* duplicate, and a missed duplicate leaves geometry exactly
where the CAD put it. The unsafe direction, deleting a real instance, is unreachable: two
genuinely distinct placements are never within 1e-9 model units of each other.

---

## 3. The permanent check

`verify_placement_invariant()` runs on **every** conversion, not as a one-off:

* leaf placements in == leaf placements out (distinct), and
* no two leaf placements share a definition *and* a world transform.

It **raises**, it does not warn. A geometry with two volumes in the same place has no defined
transport through it, so emitting one is worse than refusing to.

`report_duplicate_placements()` prints on every conversion too, including the clean case:

```
Placement check: 13 leaf placement(s), all at distinct world transforms.
```

so that silence never means "not checked".

---

## 4. Evidence

### 4.1 The tests were written first and watched to fail

Four fixtures, 12 checks, added as the fourth `--self-test` block
(`run_duplicate_placement_self_test()`). With the rule's result computed but not applied, **5 of
the 12 failed, each for the stated reason**:

```
  [FAIL] it converts to 3 leaf placements, not 6 -- 6 placed
  [FAIL] and no two of them share a definition and a world transform -- 3 distinct (definition, world transform) pair(s)
  [FAIL] a third instance that coincides with the second is the one that goes -- 3 declared -> 3 placed at 2 distinct transform(s)
  [FAIL] the same sub-assembly at the same transform down two assembly paths is placed once -- 4 declared -> 4 placed
  [FAIL] and the two paths' own, distinct parts both survive -- 3 distinct (definition, world transform) pair(s)
```

The seven that passed either way are the ones that assert the *fixture* has the shape it claims
(`6 declared, multiplicity {2: 3}`) and the negative controls — which is the correct division: a
negative control that only passes after the fix is not a control.

The fixtures are built as XCAF assemblies **in memory** and pushed through the production
`expand_free_shapes()` traversal. In memory rather than through a STEP file on purpose: the
structure under test is pathological, and building it directly is the only way to be sure the
fixture has the shape the test claims rather than whatever a STEP writer normalised it into.

1. **The ALICE3 shape** — a root whose first child contains its own three siblings, all at
   identity. 6 declared → **3** placed, `{1: 3}`, `3 by root-containment`.
2. **The negative control that matters most** — one sub-assembly instanced twice at **different**
   transforms. Still **2** placements, at 2 distinct world transforms, **nothing suppressed**.
   If the rule ever keys on the definition instead of the pair, this reports 1.
   Plus its sibling: bolt a *third* instance coincident with the second onto the same model and
   exactly one goes, 3 → 2.
3. **A diamond** — the same sub-assembly at the same transform down two assembly paths: 4 declared
   → **3** placed, and it is the *defensive* rule that fires here, not the structural one. Each
   path's own distinct part survives.
4. **The invariant on a real corpus** — `Bagger.step`, 13 placed solids in and 13 out.

### 4.2 Numbers

`O2_CADtoTGeo.py --self-test`: **48 checks, 0 failures** in four blocks (was 36 in three).
`runOracleGate.py --self-test` 17/17. `csg/emit.py --self-test` 33/33.

Leaf placements, before → after, measured on the converted definition graph:

| corpus | declared by the STEP file | before | after |
| --- | ---: | ---: | ---: |
| `CAD_noETA.stp` (ALICE3) | 206, multiplicity `{2: 103}` | 206 `{2: 103}` | **103 `{1: 103}`** |
| `Bagger.step` | 13, `{1: 13}` | 13 `{1: 13}` | 13 `{1: 13}` |
| `as1-oc-214.stp` | 18, `{1: 18}` | 18 `{1: 18}` | 18 `{1: 18}` |
| `oTOF System V3-R92cm.step` | 62628, `{1: 62628}` | does not convert | does not convert, identically |

oTOF still dies in `triangulate_asbbox()` on `OpenCASCADE Error [Standard_ConstructionError]:
Bnd_Box is void (in Bnd_Box::Get)`, before and after, byte for byte the same message. That failure
is inside the traversal, upstream of anything this stream touches, and oTOF declares no coincident
duplicates anyway. It remains open and is not this stream's.

ALICE3's exact-surface sidecar coverage is **unchanged at 20 of 55** definitions — as it must be:
the rule removes *placements*, never definitions.

### 4.3 The independent confirmation, on the converted geometry

The nine-line signature check re-run against the **output** rather than the input — `geom.C` built
in ROOT, every leaf physical node signed by (volume name, world-frame shape bounding box, rotation)
rounded to 1e-9, walked with `TGeoIterator`:

```
converted geometry: 103 leaf placements | 103 distinct signatures | multiplicity histogram {1: 103}
```

Against the same check on the STEP input, `206 | 103 | {2: 103}`. This shares no code with the
converter's own bookkeeping.

### 4.4 No regression

| gate | before | after |
| --- | --- | --- |
| fixtures, 10/10 scored | `surface` 0/0/0/0 on 10/10; `shape` 0/0/0/0 on 2/2 | identical |
| Bagger, 13/13 scored | `surface` 0/0/0/0 on 13/13; `shape` 0/0/0/0 on 7/7, exit 0 | identical |

(`mesh` totals also unchanged; they are large by construction and are the input to the unwritten
`auto`-mode fallback policy, not a verdict.)

---

## 5. What this invalidates elsewhere

Not edited here — `NEXT.md` and `Tutorial.md` are reconciled centrally. For whoever does that:

* **`Stream_T_AssemblyOracle.md` §0** — its closing paragraph offers two readings, "converter
  defect" or "the CAD is wrong", and says the first "may be a one-line fix". §1 above settles it:
  the file declares the duplication, so it is a **CAD defect that the converter now compensates
  for**, and the fix is a rule, not a line. Everything else in §0 stands, including the
  measurement.
* **Every per-instance count in `Stream_T_AssemblyOracle.md`** was, as §0 already warned, inflated
  2×. They can now simply be re-measured: the converted model has 103 placed solids and 5253
  distinct pairs.
* **`Stream_T_AssemblyOracle.md` §3.4**'s full 206-instance overlap census should be re-scoped to
  103. The 6-of-66 surviving pairs it found among the 12 `ST2487458_01` placements were *precisely*
  the 6 duplicate pairs, and those are gone — so that census should now terminate instead of
  running 39 minutes on a solid intersected with an exact copy of itself.
* **`NEXT.md`'s ALICE3 overlap item**: `CheckOverlaps` on ALICE3 was measuring a geometry in which
  every volume had a coincident twin. Any overlap number taken on ALICE3 before today is void.
* Nothing in `Stream_J_XRay.md` moves: it scores *definitions*, and no definition changed.
