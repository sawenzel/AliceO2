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
  kPlanarDisk = 5,
};

enum CurveType : uint32_t {
  kLineSegment = 0,
  kCircularArc = 1,
};

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

/// Convert the line-segment edges of a sidecar wire into the closed polygon vertex list
/// expected by AddPlanarSurface. Consecutive edges must be connected end-to-start.
bool wireToPolygon(const std::string& file, size_t surfaceIndex, const SidecarWire& wire,
                   std::vector<O2BVHSurfaceSolid::Point2D>& polygon)
{
  polygon.clear();
  polygon.reserve(wire.edges.size());
  for (size_t e = 0; e < wire.edges.size(); ++e) {
    const auto& edge = wire.edges[e];
    if (edge.curveType != kLineSegment || edge.params.size() < 4) {
      ::Error("LoadSurfaceSolid", "%s: surface %zu: planar wires support only line-segment edges in version 1",
              file.c_str(), surfaceIndex);
      return false;
    }
    const auto& next = wire.edges[(e + 1) % wire.edges.size()];
    const double gapU = edge.params[2] - next.params[0];
    const double gapV = edge.params[3] - next.params[1];
    constexpr double kJoinTolerance = 1.e-9;
    if (std::abs(gapU) > kJoinTolerance || std::abs(gapV) > kJoinTolerance) {
      ::Error("LoadSurfaceSolid", "%s: surface %zu: wire edge %zu end does not join the next edge start",
              file.c_str(), surfaceIndex, e);
      return false;
    }
    polygon.push_back({edge.params[0], edge.params[1]});
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
        std::vector<O2BVHSurfaceSolid::Point2D> outer;
        std::vector<std::vector<O2BVHSurfaceSolid::Point2D>> inners;
        bool haveOuter = false;
        for (const auto& wire : wires) {
          std::vector<O2BVHSurfaceSolid::Point2D> polygon;
          if (!wireToPolygon(file, s, wire, polygon)) {
            return false;
          }
          if (wire.role == 0) {
            if (haveOuter) {
              ::Error("LoadSurfaceSolid", "%s: plane surface %zu has more than one outer wire", file.c_str(), s);
              return false;
            }
            outer = std::move(polygon);
            haveOuter = true;
          } else {
            inners.push_back(std::move(polygon));
          }
        }
        if (!haveOuter) {
          ::Error("LoadSurfaceSolid", "%s: plane surface %zu has no outer wire", file.c_str(), s);
          return false;
        }
        added = solid.AddPlanarSurface(point3(p, 0), point3(p, 3), point3(p, 6), outer, inners);
        break;
      }
      case kCylinder: {
        if (nParams != 14) {
          ::Error("LoadSurfaceSolid", "%s: cylinder surface %zu has %u parameters, expected 14", file.c_str(), s,
                  nParams);
          return false;
        }
        added = solid.AddCylindricalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                            innerWall);
        break;
      }
      case kCone: {
        if (nParams != 15) {
          ::Error("LoadSurfaceSolid", "%s: cone surface %zu has %u parameters, expected 15", file.c_str(), s, nParams);
          return false;
        }
        added = solid.AddConicalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                        p[14], innerWall);
        break;
      }
      case kSphere: {
        if (nParams != 14) {
          ::Error("LoadSurfaceSolid", "%s: sphere surface %zu has %u parameters, expected 14", file.c_str(), s,
                  nParams);
          return false;
        }
        added = solid.AddSphericalSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10], p[11], p[12], p[13],
                                          innerWall);
        break;
      }
      case kPlanarDisk: {
        if (nParams != 11) {
          ::Error("LoadSurfaceSolid", "%s: planar-disk surface %zu has %u parameters, expected 11", file.c_str(), s,
                  nParams);
          return false;
        }
        added = solid.AddPlanarDiskSurface(point3(p, 0), point3(p, 3), point3(p, 6), p[9], p[10]);
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

} // namespace base
} // namespace o2
