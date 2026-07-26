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

/// \file O2SurfaceSolidIO.cxx
/// \brief Reader for the exact-surface sidecar format (surfaces_*.bin, version 1).
///
/// The binary layout is documented in scripts/geometry/BVHSurfaceSolid.md
/// ("Surface sidecar format") and written by write_surfaces_bin in
/// scripts/geometry/O2_CADtoTGeo.py. Keep the three places in sync.

#include "DetectorsBase/O2SurfaceSolidIO.h"
#include "DetectorsBase/O2BVHSurfaceSolid.h"
#include "DetectorsBase/O2Tessellated.h"

#include <TError.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <vector>

namespace o2
{
namespace base
{

namespace
{

constexpr uint32_t kSidecarVersion = 1;
constexpr uint32_t kFlagInnerWall = 1u << 0;

enum SurfaceType : uint32_t {
  kPlane = 1,
  kCylinder = 2,
  kCone = 3,
  kSphere = 4,
  kTorus = 5,
};

enum CurveType : uint32_t {
  kLineSegment = 0,
  kCircularArc = 1,
  kBSpline2D = 2,
};

/// Parse a bspline (curveType 2) edge record laid out as
/// [degree, nPoles, poles(2*nPoles), weights(nPoles), knots(nPoles+degree+1)] into a public
/// PlanarBoundaryCurve. Returns false on a malformed / truncated record.
bool parseBSplineEdge(const std::vector<double>& params, O2BVHSurfaceSolid::PlanarBoundaryCurve& curve)
{
  if (params.size() < 2) {
    return false;
  }
  const int degree = static_cast<int>(std::lround(params[0]));
  const int nPoles = static_cast<int>(std::lround(params[1]));
  if (degree < 1 || nPoles < degree + 1) {
    return false;
  }
  const size_t nKnots = static_cast<size_t>(nPoles) + degree + 1;
  const size_t expected = 2 + 2 * static_cast<size_t>(nPoles) + static_cast<size_t>(nPoles) + nKnots;
  if (params.size() < expected) {
    return false;
  }
  std::vector<O2BVHSurfaceSolid::Point2D> poles(nPoles);
  size_t offset = 2;
  for (int i = 0; i < nPoles; ++i) {
    poles[i] = {params[offset], params[offset + 1]};
    offset += 2;
  }
  std::vector<double> weights(nPoles);
  for (int i = 0; i < nPoles; ++i) {
    weights[i] = params[offset++];
  }
  std::vector<double> knots(nKnots);
  for (size_t i = 0; i < nKnots; ++i) {
    knots[i] = params[offset++];
  }
  curve = O2BVHSurfaceSolid::PlanarBoundaryCurve::makeBSpline(degree, std::move(poles), std::move(weights),
                                                             std::move(knots));
  return true;
}

struct SidecarEdge {
  uint32_t curveType = 0;
  std::vector<double> params;
};

struct SidecarWire {
  uint32_t role = 0; // 0 = outer, 1 = inner
  std::vector<SidecarEdge> edges;
};

bool readU32(std::ifstream& in, uint32_t& value)
{
  in.read(reinterpret_cast<char*>(&value), sizeof(value));
  return static_cast<bool>(in);
}

bool readDoubles(std::ifstream& in, std::vector<double>& values, uint32_t n)
{
  values.resize(n);
  in.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(n) * sizeof(double));
  return static_cast<bool>(in);
}

O2BVHSurfaceSolid::Point3D point3(const std::vector<double>& p, size_t offset)
{
  return {p[offset], p[offset + 1], p[offset + 2]};
}

/// Start/end (u, v) endpoints of a sidecar edge (line segment or circular arc).
bool edgeEndpoints(const SidecarEdge& edge, O2BVHSurfaceSolid::Point2D& start, O2BVHSurfaceSolid::Point2D& end)
{
  if (edge.curveType == kLineSegment && edge.params.size() >= 4) {
    start = {edge.params[0], edge.params[1]};
    end = {edge.params[2], edge.params[3]};
    return true;
  }
  if (edge.curveType == kCircularArc && edge.params.size() >= 5) {
    const double cu = edge.params[0], cv = edge.params[1], r = edge.params[2];
    const double a0 = edge.params[3], a1 = edge.params[3] + edge.params[4];
    start = {cu + r * std::cos(a0), cv + r * std::sin(a0)};
    end = {cu + r * std::cos(a1), cv + r * std::sin(a1)};
    return true;
  }
  if (edge.curveType == kBSpline2D) {
    // for a clamped B-spline the first/last poles are the curve endpoints
    O2BVHSurfaceSolid::PlanarBoundaryCurve curve;
    if (!parseBSplineEdge(edge.params, curve)) {
      return false;
    }
    start = curve.poles.front();
    end = curve.poles.back();
    return true;
  }
  return false;
}

/// Convert a sidecar wire into the public PlanarBoundaryCurve loop expected by
/// AddCurvedPlanarSurface, validating that consecutive edges join end-to-start. Sets
/// \a anyArc when the wire carries at least one circular-arc edge.
bool wireToCurves(const std::string& file, size_t surfaceIndex, const SidecarWire& wire,
                  std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& curves, bool& anyArc)
{
  using Curve = O2BVHSurfaceSolid::PlanarBoundaryCurve;
  curves.clear();
  curves.reserve(wire.edges.size());
  // The extractor projects/samples curve endpoints in the surface frame, so a line vertex and the
  // shared endpoint of a neighbouring arc or B-spline can differ by the extractor's ~1e-6 sampling
  // precision. The kernel's own CurveWire closure tolerance is looser, so accept that join gap
  // here (a stricter 1e-9 wrongly rejected mixed line/arc/bspline loops such as D-shaped caps).
  constexpr double kJoinTolerance = 1.e-5;
  for (size_t e = 0; e < wire.edges.size(); ++e) {
    const auto& edge = wire.edges[e];
    O2BVHSurfaceSolid::Point2D start{}, end{};
    if (!edgeEndpoints(edge, start, end)) {
      ::Error("LoadSurfaceSolid", "%s: surface %zu: unsupported or malformed wire edge %zu", file.c_str(),
              surfaceIndex, e);
      return false;
    }
    O2BVHSurfaceSolid::Point2D nextStart{}, nextEnd{};
    if (!edgeEndpoints(wire.edges[(e + 1) % wire.edges.size()], nextStart, nextEnd)) {
      return false;
    }
    if (std::abs(end[0] - nextStart[0]) > kJoinTolerance || std::abs(end[1] - nextStart[1]) > kJoinTolerance) {
      ::Error("LoadSurfaceSolid", "%s: surface %zu: wire edge %zu end does not join the next edge start",
              file.c_str(), surfaceIndex, e);
      return false;
    }
    if (edge.curveType == kCircularArc) {
      anyArc = true;
      curves.push_back(Curve::makeArc({edge.params[0], edge.params[1]}, edge.params[2], edge.params[3],
                                      edge.params[3] + edge.params[4]));
    } else if (edge.curveType == kBSpline2D) {
      anyArc = true; // a bspline is a curved edge, so route the plane through AddCurvedPlanarSurface
      Curve curve;
      if (!parseBSplineEdge(edge.params, curve)) {
        ::Error("LoadSurfaceSolid", "%s: surface %zu: malformed bspline wire edge %zu", file.c_str(), surfaceIndex, e);
        return false;
      }
      curves.push_back(std::move(curve));
    } else {
      curves.push_back(Curve::makeLine(start, end));
    }
  }
  return true;
}

/// Collect a quadric surface's trim wire block into one outer loop plus inner (hole) loops of
/// PlanarBoundaryCurve segments, expressed in the surface's parametric (u, v) domain
/// (u = phi[rad]; v = height[cm] for cylinder/cone, v = theta[rad] for sphere).
bool collectQuadricTrim(const std::string& file, size_t surfaceIndex, const std::vector<SidecarWire>& wires,
                        std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& outer,
                        std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>>& inners)
{
  bool haveOuter = false;
  bool anyArc = false; // quadric domains accept both line and arc trim edges
  for (const auto& wire : wires) {
    std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> curves;
    if (!wireToCurves(file, surfaceIndex, wire, curves, anyArc)) {
      return false;
    }
    if (wire.role == 0) {
      if (haveOuter) {
        ::Error("LoadSurfaceSolid", "%s: quadric surface %zu has more than one outer trim wire", file.c_str(),
                surfaceIndex);
        return false;
      }
      outer = std::move(curves);
      haveOuter = true;
    } else {
      inners.push_back(std::move(curves));
    }
  }
  if (!haveOuter) {
    ::Error("LoadSurfaceSolid", "%s: quadric surface %zu trim block has no outer wire", file.c_str(), surfaceIndex);
    return false;
  }
  return true;
}

} // namespace

bool LoadSurfaceSolid(const std::string& file, O2BVHSurfaceSolid& solid)
{
  std::ifstream in(file, std::ios::binary);
  if (!in) {
    ::Error("LoadSurfaceSolid", "Cannot open surface sidecar file %s", file.c_str());
    return false;
  }

  char magic[4];
  in.read(magic, sizeof(magic));
  if (!in || std::memcmp(magic, "O2SS", 4) != 0) {
    ::Error("LoadSurfaceSolid", "%s is not a surface sidecar file (bad magic)", file.c_str());
    return false;
  }

  uint32_t version = 0, nSurfaces = 0, reserved = 0;
  if (!readU32(in, version) || !readU32(in, nSurfaces) || !readU32(in, reserved)) {
    ::Error("LoadSurfaceSolid", "%s: truncated header", file.c_str());
    return false;
  }
  if (version != kSidecarVersion) {
    ::Error("LoadSurfaceSolid", "%s: unsupported sidecar version %u (reader supports %u)", file.c_str(), version,
            kSidecarVersion);
    return false;
  }

  for (size_t s = 0; s < nSurfaces; ++s) {
    uint32_t surfaceType = 0, flags = 0, nParams = 0;
    if (!readU32(in, surfaceType) || !readU32(in, flags) || !readU32(in, nParams)) {
      ::Error("LoadSurfaceSolid", "%s: truncated surface record %zu", file.c_str(), s);
      return false;
    }
    std::vector<double> p;
    if (!readDoubles(in, p, nParams)) {
      ::Error("LoadSurfaceSolid", "%s: truncated parameters of surface %zu", file.c_str(), s);
      return false;
    }

    // The wire block is self-describing; read it unconditionally.
    uint32_t nWires = 0;
    if (!readU32(in, nWires)) {
      ::Error("LoadSurfaceSolid", "%s: truncated wire count of surface %zu", file.c_str(), s);
      return false;
    }
    std::vector<SidecarWire> wires(nWires);
    for (auto& wire : wires) {
      uint32_t nEdges = 0;
      if (!readU32(in, wire.role) || !readU32(in, nEdges)) {
        ::Error("LoadSurfaceSolid", "%s: truncated wire header in surface %zu", file.c_str(), s);
        return false;
      }
      wire.edges.resize(nEdges);
      for (auto& edge : wire.edges) {
        uint32_t nCurveParams = 0;
        if (!readU32(in, edge.curveType) || !readU32(in, nCurveParams) || !readDoubles(in, edge.params, nCurveParams)) {
          ::Error("LoadSurfaceSolid", "%s: truncated edge record in surface %zu", file.c_str(), s);
          return false;
        }
      }
    }

    const bool innerWall = (flags & kFlagInnerWall) != 0;
    bool added = false;

    switch (surfaceType) {
      case kPlane: {
        if (nParams != 9) {
          ::Error("LoadSurfaceSolid", "%s: plane surface %zu has %u parameters, expected 9", file.c_str(), s, nParams);
          return false;
        }
        // Read every wire as a general line/arc loop. A pure line-segment loop keeps the
        // polygon path (AddPlanarSurface, general-metric); any arc routes to the curved path.
        std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> outer;
        std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>> inners;
        bool haveOuter = false;
        bool anyArc = false;
        for (const auto& wire : wires) {
          std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> curves;
          if (!wireToCurves(file, s, wire, curves, anyArc)) {
            return false;
          }
          if (wire.role == 0) {
            if (haveOuter) {
              ::Error("LoadSurfaceSolid", "%s: plane surface %zu has more than one outer wire", file.c_str(), s);
              return false;
            }
            outer = std::move(curves);
            haveOuter = true;
          } else {
            inners.push_back(std::move(curves));
          }
        }
        if (!haveOuter) {
          ::Error("LoadSurfaceSolid", "%s: plane surface %zu has no outer wire", file.c_str(), s);
          return false;
        }
        if (anyArc) {
          added = solid.AddCurvedPlanarSurface(point3(p, 0), point3(p, 3), point3(p, 6), outer, inners);
        } else {
          const auto toPolygon = [](const std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& curves) {
            std::vector<O2BVHSurfaceSolid::Point2D> polygon;
            polygon.reserve(curves.size());
            for (const auto& c : curves) {
              polygon.push_back(c.lineStart);
            }
            return polygon;
          };
          std::vector<std::vector<O2BVHSurfaceSolid::Point2D>> innerPolys;
          innerPolys.reserve(inners.size());
          for (const auto& inner : inners) {
            innerPolys.push_back(toPolygon(inner));
          }
          added = solid.AddPlanarSurface(point3(p, 0), point3(p, 3), point3(p, 6), toPolygon(outer), innerPolys);
        }
        break;
      }
      case kCylinder: {
        if (nParams != 14) {
          ::Error("LoadSurfaceSolid", "%s: cylinder surface %zu has %u parameters, expected 14", file.c_str(), s,
                  nParams);
          return false;
        }
        if (wires.empty()) {
          added = solid.AddCylindricalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12],
                                              p[13], innerWall);
        } else {
          std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> outer;
          std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>> inners;
          if (!collectQuadricTrim(file, s, wires, outer, inners)) {
            return false;
          }
          added = solid.AddCylindricalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12],
                                              p[13], innerWall, outer, inners);
        }
        break;
      }
      case kCone: {
        if (nParams != 15) {
          ::Error("LoadSurfaceSolid", "%s: cone surface %zu has %u parameters, expected 15", file.c_str(), s, nParams);
          return false;
        }
        if (wires.empty()) {
          added = solid.AddConicalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                          p[14], innerWall);
        } else {
          std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> outer;
          std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>> inners;
          if (!collectQuadricTrim(file, s, wires, outer, inners)) {
            return false;
          }
          added = solid.AddConicalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                          p[14], innerWall, outer, inners);
        }
        break;
      }
      case kSphere: {
        if (nParams != 14) {
          ::Error("LoadSurfaceSolid", "%s: sphere surface %zu has %u parameters, expected 14", file.c_str(), s,
                  nParams);
          return false;
        }
        if (wires.empty()) {
          added = solid.AddSphericalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12],
                                            p[13], innerWall);
        } else {
          std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> outer;
          std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>> inners;
          if (!collectQuadricTrim(file, s, wires, outer, inners)) {
            return false;
          }
          added = solid.AddSphericalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12],
                                            p[13], innerWall, outer, inners);
        }
        break;
      }
      case kTorus: {
        if (nParams != 15) {
          ::Error("LoadSurfaceSolid", "%s: torus surface %zu has %u parameters, expected 15", file.c_str(), s,
                  nParams);
          return false;
        }
        if (wires.empty()) {
          added = solid.AddToroidalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                           p[14], innerWall);
        } else {
          std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> outer;
          std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>> inners;
          if (!collectQuadricTrim(file, s, wires, outer, inners)) {
            return false;
          }
          added = solid.AddToroidalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                           p[14], innerWall, outer, inners);
        }
        break;
      }
      default:
        ::Error("LoadSurfaceSolid", "%s: surface %zu has unknown surface type %u", file.c_str(), s, surfaceType);
        return false;
    }

    if (!added) {
      ::Error("LoadSurfaceSolid", "%s: surface %zu was rejected by O2BVHSurfaceSolid", file.c_str(), s);
      return false;
    }
  }

  return true;
}

bool LoadFacetSolid(const std::string& file, O2Tessellated& solid)
{
  std::ifstream in(file, std::ios::binary);
  if (!in) {
    ::Error("LoadFacetSolid", "Cannot open facet sidecar file %s", file.c_str());
    return false;
  }

  uint32_t nTriangles = 0;
  if (!readU32(in, nTriangles)) {
    ::Error("LoadFacetSolid", "%s: truncated header", file.c_str());
    return false;
  }

  float v[9];
  uint32_t nDegenerate = 0;
  for (uint32_t i = 0; i < nTriangles; ++i) {
    in.read(reinterpret_cast<char*>(v), sizeof(v));
    if (!in) {
      ::Error("LoadFacetSolid", "%s: truncated facet record %u of %u", file.c_str(), i, nTriangles);
      return false;
    }
    const O2Tessellated::Vertex_t p0(v[0], v[1], v[2]);
    const O2Tessellated::Vertex_t p1(v[3], v[4], v[5]);
    const O2Tessellated::Vertex_t p2(v[6], v[7], v[8]);
    if (!solid.AddFacet(p0, p1, p2)) {
      // AddFacet rejects a facet whose vertices collapse (TGeoFacet::CompactFacet leaves < 3).
      // That is a mesh-quality property of the input, not an I/O or format error, and real CAD
      // tessellations do contain a few slivers: discarding a whole 200k-triangle reference mesh
      // over one of them would be the wrong trade. Count and carry on -- only an unreadable or
      // truncated file is fatal, as the header contract says.
      ++nDegenerate;
      continue;
    }
  }
  if (nDegenerate > 0) {
    ::Warning("LoadFacetSolid", "%s: skipped %u degenerate facet(s) of %u", file.c_str(), nDegenerate, nTriangles);
  }

  return true;
}

} // namespace base
} // namespace o2
