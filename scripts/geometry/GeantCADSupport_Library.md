The present CAD/STEP->TGeo tooling is ALICE specific or works for
customers in the VirtualMC space using TGeo as their backend.

Things that have been achieved:

- Fast TGeoTessellated (O2Tessellated)
- New/fast BVHSurfaceSolid
- Fast TGeoAssembly using BVH as a means to compose 
- Possibly new "BVHCompositeCSG" --> flat boolean vs recursive constructs as generalization of G4MultiUnion
- STEP->TGeo translation including primitive/CSG recognition
- shape healing for meshes (via CGAL)
- binary formats doing better than GDML
- attribution of materials


# Goal:

Ultimately, I would like to provide all of this as a library that all CERN/LHC experiments
could use. Users should be able to extend their core Geant4 setup with functionality from the library.
Extension could be favoured over making this a core component of Geant4 itself because Geant4 doesn't like external dependencies (such as Embree/OCCT/...)

The library should support TGeo as well as native Geant4 backends.
AliceO2 can get rid of the current directory and just use this library.

I would like to write a comprehensive paper about this library. Something like "Redefining Geant CAD simulations"... totally humble.

# Ideas:

The library can have multiple components:
- foundations (BVH)
- solids
- STEP converters (including primitive recognition, healing)
- materials
- persistence (GMDL or binary)

Solids would exist in TGeo and G4 variants.

The BVH could come from madman BVH open source implementation on github
or use IntelEmbree speed-of-light implementation with vectorization. This can be optional.
Either we code all slids in 2 implemenentations or we find a way of (zero-cost) abstraction.

Additional Python interfaces in the sense of pyg4ometry would be nice.

I do not want to compete with pyg4ometry but maybe add additional functionality.
