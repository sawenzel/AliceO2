# Stream I — the representation-aware gate verdict

Date: 2026-08-01. Closes the gap named in [`Tutorial.md`](Tutorial.md) §5.3 and §6 and recorded in
[`Stream_H_CSGEmitter.md`](Stream_H_CSGEmitter.md) §7 and §9. Written by the session that built it;
the integrating session folds it.

**In one paragraph.** The converter ships a cascade — CSG 7, exact surfaces 5, tessellated 1 of
Bagger's 13 leaf solids — and `runOracleGate.py` scored all three representations side by side but
still computed **pass/fail on the surface representation alone**. Three Bagger rams were therefore
reported as failing on a `Capacity()` deviation of ~1.3e-06 belonging to a representation they no
longer ship in, while their shipped form has a symmetric difference against OpenCascade of exactly
zero. The verdict is now computed on the representation the part **actually ships in**, read from
the converter's own cascade record and from nothing else. Bagger goes **9/12 → 12/12 on the shipped
representation**, with the surface number **still printed, still 9/12, and still computed by the
same unchanged function**. No numeric field in `gate.json` moved on either corpus. Two things fell
out on the way that are findings rather than plumbing: the gate has always been quietly scoring
**9 of 10** fixture leaf solids and **12 of 13** Bagger ones, and `bboxDeviationFromOracle` on the
`shape` column is **0.13–0.62 cm on every CSG-carried ram**, not the ~1e-07 the Stream G table
predicted (§6).

---

## 1. The rule, and the one part of it that is a correctness constraint

A part passes if **the representation the converter dispatches it to** is oracle-clean — same four
columns, same "outside tolerance" definition, same answer set — with two per-representation
refinements:

| | criterion |
| --- | --- |
| **navigability** | required exactly where `closureApplicable` is true, i.e. on the surface solid. A `TGeoCompositeShape` has no rims and a mesh's `meshClosedBody` is a different claim; neither is silently assumed true. |
| **volume** | `dV_sym` where the CSG emitter computed one; else the 1e-6 relative capacity band where capacity is a real measurement (`exact-divergence`, `mesh-divergence`); else **nothing**, reported and not gated. |

`TGeoCompositeShape::Capacity()` is **never** gated on. It is Monte-Carlo sampled in ROOT (~1e-2
relative), and Stream G measured a geometrically exact composite whose ROOT capacity was 3.3e-04
off — a capacity gate would have failed a correct shape by 330×. There are two independent guards:
`capacityComparable` is already false for `root-montecarlo`, and `dV_sym` takes priority wherever
the emitter produced one.

**The constraint that matters.** The shipped representation is read from the converter's cascade
decision — `csg_report.json`, copied per part into `manifest.json` as a `shipped` block — and
**never** by picking whichever representation scores best. That is not a stylistic preference. A
best-of rule would let a part pass because some representation it does not use happens to be clean,
and would make the gate structurally incapable of ever reporting a bad conversion choice. It would
also make this change unfalsifiable, which is the failure mode the project distrusts most.

**Why that source cannot be gamed.** `csg_report.json` is written by `csg/hook.write_report()` from
the same `records` that decide what `geom.C` emits, in the same converter run, before any oracle
answer exists. A row is `csg` only when the emitter both accepted the candidate *and* wrote the
`.root` file — the converter never dispatches to a file it did not write — so the gate's notion of
"shipped" is the same predicate the macro uses. The gate has no code path that inspects scores
before choosing a representation; `shipped_block()` takes the manifest's statement or, when the
database predates the record, says `default (this part DB predates the cascade record)` in the
output rather than guessing quietly. And the failure mode is loud in the right direction: a part
whose cascade says `csg` but whose `shape` column is absent **fails**, it does not fall back to a
clean surface column (measured, §5.2).

---

## 2. The per-part verdict table

Bagger, `--csg auto`, 2000 points / 2000 rays, seed 1. "old" is the historical surface-only
verdict; "new" is the verdict on the shipped representation.

| part | ships as | old | new | evidence for the new verdict |
| --- | --- | --- | --- | --- |
| `BasePin` | shape (csg) | PASS | PASS | `TGeoCompositeShape`, oracle 0/0/0/0, `dV_sym = 0` cm³ (band 6.91e-06) |
| `Base` | surface | PASS | PASS | oracle 0/0/0/0, capacity 4.83e-15 |
| `BoomCylinderInner` | shape (csg) | PASS | PASS | oracle 0/0/0/0, `dV_sym = 0` (band 8.57e-06) |
| **`BoomCylinderOuter`** | **shape (csg)** | **FAIL** | **PASS** | oracle 0/0/0/0 on the composite; `dV_sym = 0` (band 2.11e-05). Old reason: surface `Capacity()` off by **1.39e-06** — a property of the surface representation's trim data, still reported |
| `Boom` | surface | PASS | PASS | oracle 0/0/0/0, capacity 1.82e-12 |
| **`BucketCylinderInner`** | **shape (csg)** | **FAIL** | **PASS** | oracle 0/0/0/0; `dV_sym = 0` (band 5.34e-06). Old reason: surface capacity **1.31e-06** |
| `BucketCylinderOuter` | shape (csg) | PASS | PASS | oracle 0/0/0/0, `dV_sym = 0` (band 1.15e-05) |
| `BucketLink1` | surface | PASS | PASS | oracle 0/0/0/0, capacity 2.94e-12 |
| `BucketLink2` | surface | PASS | PASS | oracle 0/0/0/0, capacity 1.00e-12 |
| `StickCylinderInner` | shape (csg) | PASS | PASS | oracle 0/0/0/0, `dV_sym = 0` (band 8.57e-06) |
| **`StickCylinderOuter`** | **shape (csg)** | **FAIL** | **PASS** | oracle 0/0/0/0; `dV_sym = 0` (band 2.11e-05). Old reason: surface capacity **1.39e-06** |
| `Stick` | surface | PASS | PASS | oracle 0/0/0/0, capacity 1.52e-12 |
| `Bucket` | mesh | *absent* | **UNSCORED** | 97 faces, no exact sidecar; the harness cannot score it (§4) |

**Exactly three verdicts changed**, all in the same direction and all for the same reason: the
criterion that failed them belonged to `O2BVHSurfaceSolid` and those parts ship as
`TGeoCompositeShape`. Nothing about the surface representation changed — its capacity deviations
are still 1.39/1.31/1.39e-06, still computed, still printed on the same line as the new verdict.

Fixtures, `--csg auto` — cascade **CSG 2, exact surfaces 7, tessellated 1** of 10 leaf solids:
**9/9 → 9/9**, unchanged. Two fixtures now ship as CSG (`box` →
`TGeoBBox`, `cyl_cross_cyl` → `TGeoCompositeShape`, both `dV_sym = 0`) and passed either way, so
the corpus contains a case where the rule changed the *object* being gated without changing the
*verdict* — which is the control that says the fixtures result is not an artefact of nothing having
happened.

---

## 3. Totals and disagreement counts, quoted together

Never one without the other.

| | before | after |
| --- | --- | --- |
| `ctest -R BVHSurfaceSolid` | 80 cases, green | **80 cases, green** (no C++ changed) |
| fixtures gate, shipped | 9/9 | **9/9** |
| fixtures gate, surface (historical) | 9/9 | **9/9** |
| Bagger gate, shipped | 9/12 | **12/12** |
| Bagger gate, surface (historical) | 9/12 | **9/12, the same three, the same values** |
| unexplained disagreements, **surface**, fixtures | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| unexplained disagreements, **surface**, Bagger | 0 / 0 / 0 / 0 | **0 / 0 / 0 / 0** |
| disagreements, **shape**, fixtures | 0/0/0/0 (2 parts) | **0/0/0/0 (2 parts, now CSG-emitted rather than hand-built)** |
| disagreements, **shape**, Bagger | 0/0/0/0 (7 parts) | **0/0/0/0 (7 parts)** |
| disagreements, **mesh**, Bagger | contains 418, distout 6721, distin 7973, safety 10299 — **0/12 clean** | **identical, and still printed in full** |
| leaf solids scored | 12 of 13 (Bagger), **unstated** | 12 of 13, **stated** |

The invariant is unchanged: **zero unexplained oracle disagreements outside tolerance, per
representation, on both corpora, on every column.** The gate total moved because the rule now names
a different object; the measurements did not move at all, which is the next section.

---

## 4. Diff columns, not verdicts

The comparison that matters is *same converter settings, before and after the code change*. Both
corpora were reconverted with `--csg off` — the setting under which the two verdicts are identical
by construction — and `gate.json` compared leaf field by leaf field with `timing*`, `*Seconds` and
`nsPerCall` stripped.

| corpus | leaf fields (baseline) | removed | changed | of which non-path | added |
| --- | ---: | ---: | ---: | ---: | ---: |
| fixtures (9 parts) | 8315 | **0** | 20 | **0** | 189 |
| Bagger (12 parts) | 14571 | **0** | 24 | **0** | 258 |

Every changed field is a `representations[..].source` path differing only because the two runs used
different scratch directories. **No numeric field moved on either corpus.** All added fields are
under the new per-part `verdict` key. `compareGateRuns.py`, which applies the declared scale law per
column, reports the same: the only non-added differences are `source` strings.

**Positive control, because a comparison that cannot fail is not evidence.** Two independent ones
were run before believing the table above:

* `compareGateRuns.py --self-test` injected four defects into the Bagger baseline (an extra
  unexplained disagreement, a length at +1e-8 relative, a volume at +1e-8 relative, a downgraded
  reliability string) — **4/4 caught**;
* the flatten-and-compare used for the table was pointed at a baseline with exactly two fields
  nudged (`nMismatchUnexplained` 0 → 1, `capacityRelativeDeviation` +1e-6 relative) and reported
  **exactly those two**, both classified non-path.

---

## 5. The negative control: the rule still fails what deserves to fail

A scoring rule that cannot fail is not a gate. Three controls, of increasing realism.

### 5.1 Forcing a part into the tessellated tier — end to end

`BasePin`'s row in a *copy* of a finished Bagger workdir's `csg_report.json` was edited from
`csg` to `mesh`, the database re-indexed (`makeTestPartDB.py --skip-existing`, which copies the
cascade decision into `manifest.json`), and the gate re-run with `--skip-convert`. Nothing else
changed; the mesh in question is the one the converter actually wrote.

```
[FAIL] Bagger/BasePin_0_1_1_2   ships: mesh (tier mesh, o2::base::O2Tessellated)   [surface verdict: PASS]
         volume criterion: capacity off by 0.000414 relative
         contains:  19 disagreement(s) outside tolerance (missed=0)
         distout: 1011 disagreement(s) outside tolerance (missed=0)
         distin:   922 disagreement(s) outside tolerance (missed=0)
         safety:   573 disagreement(s) outside tolerance (missed=0)

11/12 scored part(s) pass on the representation they ship in
9/12  scored part(s) pass on the surface representation
```

The part's surface and CSG representations are both still clean, and the new rule fails it anyway,
because the mesh is what it was declared to ship in. Exit code 1.

### 5.2 A cascade that claims CSG for a shape that is not there

The same copied workdir with `shape_BasePin_0_1_1_2.root` deleted **from the copy only** (the
original still on disk). The gate re-rooted the manifest onto the copy, found no `shape` column,
and failed the part:

```
[note] manifest.json was written for /tmp/I_new_bag/db but lives in /tmp/I_reb_bag/db;
       re-rooting its absolute paths onto this copy
[FAIL] Bagger/BasePin_0_1_1_2   ships: shape (tier csg, ?)   [surface verdict: PASS]
         the part ships as 'shape' but has no 'shape' representation in the scorecard
```

Two statements at once. The verdict does **not** quietly fall back to a representation that happens
to be clean. And the run scored the *copy*: pre-change it would have read the original's files
through the absolute paths in `manifest.json` and reported a perfectly healthy `shape` column,
which is exactly the trap that cost Stream H a wrong conclusion.

### 5.3 The rule's own self-checks

`runOracleGate.py --self-test` — no build, no model, no oracle, **17 checks, every positive case
paired with the negative one that must fail**: `dV_sym` inside/outside the band; an exact composite
passing despite a 3.3e-04 Monte-Carlo capacity *and* the same composite failing on real column
disagreements; a part forced to the mesh tier failing while its surface stays clean; the new rule
agreeing with the historical one for every surface-shipped part (which is what makes `--csg off`
inert); an open surface solid failing as `surface` but not inheriting that failure as CSG; a
missing `shape` column failing rather than falling back; and the provenance rule itself — the
manifest's `mesh` beats a perfectly scoring `shape`.

---

## 6. Two findings that were not the task

### 6.1 The gate has always been scoring fewer parts than the corpus has

Requirement 5 was about Bagger's `Bucket`: 13 leaf solids, 12 scored, nothing in the output saying
so. Making that visible surfaced the same hole in the **fixture ladder**, which nobody had noticed:

```
1 further leaf solid(s) ship in a representation this harness cannot score:
  [UNSCORED] oblique_cut_cyl/oblique_cut_cyl_0_1_1_1   ships: mesh (tier mesh, 3 faces)
  => 9 of 10 leaf solid(s) in the model(s) are scored by this gate.
```

`oblique_cut_cyl` is a three-face ladder fixture that has never had an exact sidecar and has
therefore never been in any gate total this project has quoted. "Fixtures 9/9" has always meant
"9 of the 10 fixtures that could be scored". The number is now printed with its denominator, both
corpora, every run, and an unscoreable leaf solid keeps the exit code non-zero.

Making these parts *scoreable* is a different job: the harness's sample generator and reference both
hang off the surface solid, so a sidecar-less part cannot enter `manifest["parts"]` without
restructuring `runSolidHarness.cxx`. They are listed in `manifest["unscoredParts"]` (ignored by the
C++ side) and in `<workdir>/unscored.json`.

### 6.2 `bboxDeviationFromOracle` on the `shape` column is not ~1e-07

`Stream_G_AnyShape.md` tabulates the frame check as "`shape` (analytic ROOT shape): 1.0e-07 cm —
`TGeoBBox`/`TGeoCompositeShape` compute a tight box". Measured on the six Bagger rams, it is
**0.13 to 0.62 cm**:

| part | `shape` | `surface` | `mesh` |
| --- | ---: | ---: | ---: |
| `BasePin` | 1.00e-07 | 9.90e-08 | 3.11e-04 |
| `BoomCylinderInner` | 4.9706e-01 | 4.9706e-01 | 1.32e-04 |
| `BoomCylinderOuter` | 6.2132e-01 | 6.2132e-01 | 2.21e-04 |
| `BucketCylinderInner` | 1.2971e-01 | 1.2971e-01 | 1.36e-04 |
| `BucketCylinderOuter` | 1.5846e-01 | 1.5846e-01 | 3.12e-04 |
| `StickCylinderInner` | 4.9473e-01 | 4.9473e-01 | 1.39e-04 |
| `StickCylinderOuter` | 6.1841e-01 | 6.1841e-01 | 2.31e-04 |

Stream G's number was measured on two **axis-aligned** hand-built fixtures. These composites are
tubes on a **tilted** axis, and `TGeoBoolNode`'s bounding box for a rotated operand is conservative
in precisely the way `O2BVHSurfaceSolid`'s per-patch AABBs are — the two agree to **1e-9 cm**, which
is the surface solid's own documented padding, so this is one mechanism seen twice and not two
coincidences. The mesh's box is tight (1e-4) and matches the oracle's, so the oracle's box is not
the loose one and **this is not a frame error**: all four oracle columns are 0 and `dV_sym` is
exactly 0 on every one of these parts.

**Why it matters for this stream.** `bboxDeviationFromOracle` was deliberately left *out* of the
verdict, and reported only. Had it been included at anything like a 1e-3 band — which reading the
Stream G table would have suggested was safe — all six rams would have failed a frame check that is
not detecting a frame error, and the 12/12 above would have been 6/12 for a reason that has nothing
to do with geometry. The Stream G table's `shape` row should be read as "tight for an axis-aligned
primitive", not "tight".

---

## 7. What changed, in files

| file | change |
| --- | --- |
| `csg/hook.py` | `write_report()` adds `part` (the `<VOLNAME>_<LID>` artifact stem), `shapeFile` and `shapeDeferred` to each row of `csg_report.json`. The cascade decision was already machine-readable; it lacked the one key needed to join it to `manifest.json` without re-implementing the converter's filename sanitiser. No tier assignment changed. |
| `makeTestPartDB.py` | `--csg {off,auto,required}`, default **auto** — the script ran only the last two tiers of the cascade, so the database never contained the representation a CSG part ships in. Copies the cascade decision into each part's `shipped` block, and lists leaf solids it cannot index under `unscoredParts`. |
| `runOracleGate.py` | the verdict (§1), both verdicts written into `gate.json` under `verdict`, the `*` marker in the scorecard, `--csg` pass-through, `--self-test`, `occ_env()` now **prepends** (§8), and `rebase_manifest()` for the absolute-path trap. |
| *nothing in C++* | no kernel, harness or test source changed; `ctest` is 80 cases before and after. The rule's tests are the 17 self-checks, which need no build. |

`O2_CADtoTGeo.py` was **not** touched: the cascade decision was already recorded there.

---

## 8. The environment change, and why it is a fix rather than a workaround

`runOracleGate.py` used to **replace** `PYTHONPATH`/`LD_LIBRARY_PATH` with OCC-only values for its
subprocesses. That single line is why the CSG hook had to defer writing `shape_<part>.root`: a
converter launched from the gate had pythonOCC but `import ROOT` died on `libffi.so.6`, so a part
could be recognised and accepted as CSG and then not emitted — and `csg_report.json` honestly
recorded it one tier down, because the converter never dispatches to a file it did not write. Under
the gate, the cascade was therefore *always* "CSG 0" no matter what the recogniser found.

`csg/occ_env.py` had already established the fix and written down why (prepending keeps OCCT's
libraries winning while leaving the rest intact; the alibuild Python 3.10 pythonOCC is built against
*is* the O2 interpreter). `runOracleGate.py` was the one place that had not adopted it. It now
prepends, verified directly (`import OCC` and `import ROOT` both succeed in the resulting
environment) and end to end:

```
CSG recognition (auto): 7/13 leaf solid(s) accepted as native ROOT shapes (7 written)
tiers: CSG 7, exact surfaces 5, tessellated 1  (of 13 leaf solids)
```

Before, the same run wrote **0**. One conversion pass now produces both the shapes and a cascade
record that describes what `geom.C` really builds, and `csg/emit.py --from-json` is no longer part
of the gate loop.

---

## 9. Reproducing everything here

```bash
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
B=$HOME/alisw/sw/BUILD/O2-latest/O2
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-dev-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
ctest -R BVHSurfaceSolid                                  # 80 cases

cd $HOME/alisw/O2
python3 scripts/geometry/runOracleGate.py --self-test      # 17 checks, no build needed

# the two gates, with the production cascade
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate_fix --fixtures
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate_bag \
    --model scripts/geometry/STEP_examples/Bagger.step

# inertness: --csg off must reproduce the historical numbers exactly
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate_off \
    --model scripts/geometry/STEP_examples/Bagger.step --csg off
python3 scripts/geometry/compareGateRuns.py --baseline <old>/gate.json --candidate /tmp/gate_off/gate.json
python3 scripts/geometry/compareGateRuns.py --baseline <old>/gate.json --self-test

# the negative control (§5.1)
cp -r /tmp/gate_bag /tmp/gate_neg
python3 - <<'PY'
import json; p='/tmp/gate_neg/db/Bagger/csg_report.json'; d=json.load(open(p))
for r in d['parts']:
    if r['volume'] == 'BasePin': r['representation'] = 'mesh'
json.dump(d, open(p,'w'), indent=1)
PY
python3 scripts/geometry/makeTestPartDB.py --output /tmp/gate_neg/db --skip-existing \
    --models scripts/geometry/STEP_examples/Bagger.step
O2_BUILD_DIR=$B python3 scripts/geometry/runOracleGate.py --workdir /tmp/gate_neg --skip-convert
```

The gate exits non-zero whenever any part fails **or any leaf solid is unscoreable**, which on
Bagger is now the normal state. Read the counts.

---

## 10. What this leaves open

* **`Bucket` and `oblique_cut_cyl` are still unscored.** Scoring a sidecar-less part means making
  the harness's sample generator and reference independent of `O2BVHSurfaceSolid` — a real change to
  `runSolidHarness.cxx`, not a policy one. Until then the tessellated tier is measured only on parts
  that also have an exact sidecar, i.e. never on the parts that actually fall to it.
* **The `auto`-mode fallback policy is still unwritten**, and it is the reason the mesh's
  disagreement counts stay reported in full. The measurement now exists on 21 parts across two
  corpora; the policy that consumes it does not.
* **`bboxDeviationFromOracle` is reported and not gated**, on any representation. §6.2 says why a
  band would have to be per-representation and per-part-shape before it could be one; the honest
  version of that check is a subdivision BVH for the surface solid and a tighter boolean-node box
  for composites, neither of which exists.
* **A part shipping as CSG is gated on `dV_sym` plus the oracle columns, but its `Capacity()` is
  never checked at all.** That is correct today — the only composite capacity available is
  Monte-Carlo — but it means the CSG tier has no volume check independent of OCCT. An exact analytic
  capacity for a two-leaf tube union would be a cheap third opinion.
* **Nothing here is exercised beyond 13 parts.** The verdict rule is O(parts) and adds no queries,
  but the gate's lack of parallelism and sample budget is unchanged and is still what stops an
  ALICE3-sized model being gateable in bounded time.
