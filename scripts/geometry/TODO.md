- implement a BVHSurface solid as an exact representation of a CAD solid
- complete geometry configurable as JSON --> even the world volume ?
- surfaces_ST1829909_01 (largest ALICE3 exact sidecar) fails load validation: wire-join gap
  5.41e-6 cm vs the fixed 1e-6 fallback tolerance; loader should probably honour the sidecar's
  declared model tolerance, else the converter must emit exactly-joining edges — NEXT.md item 0a

 