# Index — how to read this directory

Reading order for anyone arriving fresh, then a status table of every document. Statuses:
**reference** (current, consult as needed), **record** (a finished investigation — its numbers
were true when written and later documents may supersede them), **open** (describes work not yet
done), **historical** (kept for the reasoning, no longer steering anything).

## Reading order

1. [`Tutorial.md`](Tutorial.md) — the map: the three representations, the pipeline, the
   commands, the principles. Written 2026-08-01; a few numbers postdate it (see `NEXT.md`).
2. [`Review_2026-09.md`](Review_2026-09.md) — the deep review of the mathematics and
   algorithms, read from the source, with the honest gaps.
3. [`NEXT.md`](NEXT.md) — the session-to-session hand-over: exact current state, open items,
   environment traps. Rewritten by whoever finishes a session.
4. [`Plan_Presentation.md`](Plan_Presentation.md) — the current plan of record.
5. Everything else, on demand, via the table.

## Documents

| document | status | one line |
| --- | --- | --- |
| `Tutorial.md` | reference | the orientation map |
| `Review_2026-09.md` | reference | deep review 2026-09, math and algorithms |
| `NEXT.md` | reference | the live hand-over |
| `Plan_Presentation.md` | reference | plan for the WG presentation |
| `Roadmap.md` | reference | deliberately deferred work, with reasons |
| `README.md` | reference | user-facing converter/o2-sim usage |
| `BVHSurfaceSolid.md` | reference | the solid's design log and the sidecar format |
| `TolerancePolicy.md` | reference | every tolerance constant, its value's justification |
| `SolidNavigationHarness.md` | reference | the harness: options, output, perf entry point |
| `CSG_Pipeline.md` | reference | the B-rep → CSG design and tiers |
| `Handoff_IntegrationTest.md` | open | the Geant integration test brief; basis of Track 1 |
| `MeshHealing.md` | open | mesh validity ≠ accuracy; repair options unbuilt |
| `Workstreams.md` | historical | the parallel-streams contract (waves 0–3, executed) |
| `CodeReview_Fable.md`, `_v2.md` | historical | the two defect registers that drove waves 0–1 |
| `ExactTrimTopology.md` | historical | why closure moved from proximity to topology |
| `simulating_CAD_modules.md` | historical | predates the branch; §4 build pattern still valid |
| `TODO.md` | historical | superseded by `Roadmap.md` / `NEXT.md` |

## Stream records (investigation logs, chronological by wave)

| document | status | outcome in one line |
| --- | --- | --- |
| `Stream_A_CSG.md` | record | the census: 1004 of ALICE3's 2377 "B-splines" are quadrics in disguise |
| `Stream_C_Hygiene.md` | record | the capacity residual: closed-loop `uVariation` was the mechanism |
| `Stream_E_Scale.md` | record | position-independent yes, scale-independent no → found the quartic guards |
| `Stream_F_EdgeIdentity.md` | record | closure by identity, sidecar v3; navigability 7/21 → 0/21 |
| `Stream_G_AnyShape.md` | record | the gate scores any `TGeoShape` |
| `Stream_H_CSGEmitter.md` | record | the Bagger-scope recogniser/emitter |
| `Stream_I_Verdict.md` | record | the representation-aware verdict |
| `Stream_J_XRay.md` | record | the X-ray transport benchmark; §9 scopes assembly transport |
| `Stream_K_Tier0.md` | record | recognition already existed; the vacuous cone criterion fixed |
| `Stream_L_ALICE3Defect.md` | record | the ALICE3 transport defect: three mechanisms, all fixed |
| `Stream_M_Quartic.md` | record | the dimensionless quartic; guards structural, not tuned |
| `Stream_N_PlacedPrimitives.md` | record | a placed primitive is a plain shape plus a transform |
| `Stream_O_ImplicitTrims.md` | record | 691 of 763 "free-form" trim edges are analytic∩analytic |
| `Stream_P_RepresentationBench.md` | record | composites are 40–145× faster; `closestPoint` is 61% of `Contains` |
| `Stream_Q_EllipseTrim.md` | record | an ellipse is a rational quadratic B-spline; Bagger 12→13 |
| `Stream_R_CoSurfaceTrims.md` | record | co-surface conjunction carries 1/15, not 15 — closed negative |
| `Stream_S_SafetyBVH.md` | record | Safety/ComputeNormal through the BVH, bit-identical, 835→22 µs |
| `Stream_T_AssemblyOracle.md` | record | occupancy-sequence oracle; **the models are not legal worlds** (open) |
| `Stream_U_CoSurfaceMerge.md` | record | merging exporter co-surface patches is worth 1/16 — closed negative |
| `Stream_W_DoublePlacement.md` | record | ALICE3's duplicate identity placement: 103 leaf parts, not 206 |
| `Stream_X_SubPatchBVH.md` | record | sub-patch cover boxes; safety candidates 8.6 → 3.9/call |
| `Stream_Y_SidecarJoinTolerance.md` | record | the wire-join gate honours the model's declared tolerance |

`attic/` holds parked stale script variants (pre-branch `O2_TGeoToCAD*` experiments and scratch);
nothing in it is live, and it is untracked on purpose.
