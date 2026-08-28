This is plan for **what** I would like to present a main developments
in the area of CAD-Geant simulations for the ALICE upgrade week in Sardinia September 2026.

The presentation is meant to give an intro in the "why" we are doing this
which is quicker R&D, less C++ coding as in the old days of TGeo.
The client is ALICE3 detector R&D. ALICE3 is the LHC Run5 installation of the ALICE experiment at CERN.

I have given a previous presentation on the topic, available as PDF here
https://indico.cern.ch/event/1585123/contributions/7013953/attachments/3250424/5800648/SimulationsFromCADModels_TB.pdf

which focused on

- problem description
- fast tessellated solid as foundation
- CAD->TGeo pipeline
- very initial o2-sim demonstration

In the meantime, we have worked and achieved a lot. The new presentation should highlight on

- Pipeline / algorithms to make TGeo parts an exact representation

  - New BVHSurfaceSolid !!
  - New FlatCSG solid !!
  - with nice visualizations of simple cases of each, explaining the core concepts

- Example driven development for the recognition program: First do the inverse TGeo->CAD to have "truth"

- The shape recognition program
   - how does it work in a few, easily understable sentences

- The website that can inspect solids, what features they support and benchmark
  to find out best possible representation
- Generalization of external modules: make it easier to compose a detector from multiple CAD sources via json and not hard coded into build_geometry
- Support for sensitive actions via JIT / Cling ROOT macros
- A full system closure test demonstrating everything end-to-end using current ALICE geometry
- A full example with new ALICE3 beam-pipe, inner vertex detector (IRIS) + oTOF --> giving hits in IRIS and oTOF (still to be done)
- Side-effects:

  - Uncovered duplication bugs in MFT --> Should sell this as story to be careful about this when developing ... a geomtry doctor might automatically see this in the future

Things that I would like to have which would make this talk awesome:

- nicely Blender rendered visualizations of say parts of IRIS and oTOF:
  - one BEFORE: every part is a Tessellated
  - one AFTER: parts are labeled either (SurfaceSolid/exact, CSG/exact, tessellated)
  - labels are shiny and nicely readable and directly attached to a surface

- Tables showing statistics of how many parts are recognized exactly

- A benchmark showing the cost of "all tessellated" vs "most exact" and differences in some observables.

- A benchmark investigating "accurate safety" vs "fast return safety + safety caching " (to be done)


Pain points:
- media : How do we assign magnetic field properties and physics cuts on 
- We will have to prepare a complete integration test for the ALICE3 case. Models should be VTX_noETA.root as well as the models obtained from https://mattermost.web.cern.ch/alice/messages/@njacazio (There is one oTOF and one empty beampipe which might be part of VTX_noETA or not). I don't have materials for oTOF but they could be deduced from AliceO2 source code which also has a C++ version somewhere.); Nicolo also has some raytrace image in https://mattermost.web.cern.ch/alice/messages/@njacazio/r9ninodfgb8qfk963f4rpktimh but I doubt this is about oTOF.

