# Raw run products of the integration demo

Not committed (see `.gitignore`) — regenerate with `../run_all.sh <conv-root>` and copy from
`<conv-root>/runs/<tag>/`. Kept here so Track 3b (the event display) has a stable path.

| file | what it is |
| --- | --- |
| `MCStepLoggerOutput_<tag>.root` | MCStepLogger per-step tree, one entry per event |
| `MCStepLoggerVolMap_<tag>.dat` | the volume-name map MCStepLogger writes alongside it |
| `o2sim_geometry_<rep>.root` | the assembled TGeo world, for looking volumes and media up |
| `o2sim_<species>_<rep>.root` | the o2-sim output tree, carrying `IRISHit` and `BAGRHit` |

`<tag>` is `<species>_<rep>` with species in {geantino, electron, pion, matfan} and rep in
{exact, tess}, plus `matfan_coarse` (the mesh-precision 2.0 control).

## Reading the step tree

```
TTree StepLoggerTree
  Steps    std::vector<o2::StepInfo>   one record per Geant step
  Lookups  o2::StepLookups             volume id -> name / medium / module
  Calls    magnetic-field call records
```

`o2::StepInfo` carries `stepid, volId, copyNo, trackID, t, x, y, z, E, px, py, pz, edep, mass,
step, nsecondaries, prodprocess, entered, exited, stopped, insensitiveRegion`. The volume name
of a step is `*Lookups->volidtovolname[s.volId]`. Reading it needs
`libMCStepLoggerCore` on `LD_LIBRARY_PATH` and `MCStepLogger/StepInfo.h` on
`ROOT_INCLUDE_PATH`; `../analyse_all.sh` sets both.

Track polylines for an event display come straight from `(x, y, z)` grouped by `trackID`,
in tree-entry (event) order.

## Reading the hits

`o2-sim-serial` leaves external-detector hits in the monolithic `o2sim.root`, tree `o2sim`,
branches `IRISHit` and `BAGRHit`, both `std::vector<o2::ext::Hit>`. See `../count_hits.C`.
