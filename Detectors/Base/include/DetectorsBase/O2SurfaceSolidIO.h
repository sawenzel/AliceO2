// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_BASE_O2SURFACESOLIDIO_
#define ALICEO2_BASE_O2SURFACESOLIDIO_

#include <string>

namespace o2
{
namespace base
{

class O2BVHSurfaceSolid;
class O2Tessellated;

/// Load an exact-surface sidecar file (surfaces_*.bin, version-1 format documented in
/// scripts/geometry/BVHSurfaceSolid.md, "Surface sidecar format") into \a solid by
/// dispatching to its Add*Surface methods. The caller is expected to call
/// solid.CloseShape() afterwards. Returns false (with an error message) on I/O or
/// format problems; the solid may then contain a partial surface list and should be
/// discarded.
bool LoadSurfaceSolid(const std::string& file, O2BVHSurfaceSolid& solid);

/// Load a facet mesh sidecar file (facets_*.bin, written by O2_CADtoTGeo.py's
/// write_facets_bin: little-endian uint32 triangle count followed by that many records of 9
/// float32 -- ax ay az bx by bz cx cy cz) into \a solid via AddFacet(). Mirrors the LoadFacets
/// helper that O2_CADtoTGeo.py's generated macro embeds (emit_cpp_prelude), so the harness does
/// not re-implement the format a second time. The caller is expected to call
/// solid.CloseShape() afterwards. Returns false (with an error message) on I/O or format
/// problems; the solid may then contain a partial facet list and should be discarded. Facets that
/// O2Tessellated rejects as degenerate are a property of the input mesh rather than of the file,
/// so they are skipped and counted in a warning instead of failing the load.
bool LoadFacetSolid(const std::string& file, O2Tessellated& solid);

} // namespace base
} // namespace o2

#endif
