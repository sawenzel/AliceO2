| | TB talk (SimulationsFromCADModels_TB) | now |
| --- | --- | --- |
| representations available | tessellated only | tessellated, exact BVH surface solid, native CSG tree, `O2FlatCSG` |
| what a converted part is | an approximation at a chosen mesh precision | for 99.4 % of Run 3 leaf solids, an exact native shape |
| ground truth | none — judged by eye | the TGeo → STEP writer gives every part a known source; 22 747 of 22 756 agree |
| recognition | none | one recogniser per primitive class plus a cell matcher; 2 987 of 3 125 composite-sourced parts |
| breadth tested | single example models | all 17 Run 3 modules, 22 901 leaf solids |
| `o2-sim` | a first demonstration | sensitive external detectors writing real hits; a full-system closure test on PIPE/ITS/TPC/MAG |
| material fidelity | not measured | median per-ray x/X0 difference 2.6e-12 over 2 000 Fibonacci rays |
| composing a detector | hard-coded in `build_geometry` | `externalModules` / `externalDetectors` JSON, sensitive actions as JIT ROOT macros |
| inspection and benchmarking | none | a website that loads every representation of a part and benchmarks them on one shared sample set |
