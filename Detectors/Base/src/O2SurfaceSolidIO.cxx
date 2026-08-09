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
/// \brief Reader for the exact-surface sidecar format (surfaces_*.bin, versions 1 to 3).
///
/// The binary layout is documented in scripts/geometry/BVHSurfaceSolid.md
/// ("Surface sidecar format") and written by write_surfaces_bin in
/// scripts/geometry/O2_CADtoTGeo.py. Keep the three places in sync.

#include "DetectorsBase/O2SurfaceSolidIO.h"
#include "DetectorsBase/O2BVHSurfaceSolid.h"
#include "DetectorsBase/O2Tessellated.h"

#include "BoundedSurface.h"

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

/// Sidecar versions this reader understands. Version 2 appends a float64 model tolerance (cm) to
/// the fixed header; everything after the header is unchanged, so a v1 file is a v2 file that
/// simply does not say what its model's tolerance is. Version 3 appends a uint32 edge-table size
/// to the header and, after each surface's wire block, that face's boundary edge identities -- so
/// a v2 file is a v3 file that states no identities, and it loads exactly as it always did and
/// gets exactly the geometric closure verdict it always got.
constexpr uint32_t kSidecarVersionMin = 1;
constexpr uint32_t kSidecarVersionMax = 3;

/// What a version-1 sidecar's model tolerance is taken to be, in cm. Version 1 carried no such
/// statement, so this is the project's measured extractor precision standing in for one -- the
/// same number kWireJoinTolerance is set from. It is a fallback, not a measurement of the model,
/// and the reader says so.
constexpr double kSidecarV1FallbackTolerance = 1.e-6;

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

bool readU8(std::ifstream& in, uint8_t& value)
{
  in.read(reinterpret_cast<char*>(&value), sizeof(value));
  return static_cast<bool>(in);
}

/// How many bytes are left to read, or 0 once the stream has gone bad.
///
/// Every count in this format is a `uint32` read from the file, and a file that is not the file it
/// claims to be -- a version-2 body behind a version-3 header, say -- yields counts of billions.
/// `resize`-ing to them is not a parse error, it is an out-of-memory kill, and the reader then
/// reports nothing at all. Sizing every allocation against what the file can actually hold turns
/// that back into the truncation error it is.
uint64_t bytesRemaining(std::ifstream& in, std::streamoff fileSize)
{
  if (!in) {
    return 0;
  }
  const std::streamoff here = in.tellg();
  return here < 0 || here > fileSize ? 0 : static_cast<uint64_t>(fileSize - here);
}

bool readDoubles(std::ifstream& in, std::vector<double>& values, uint32_t n, std::streamoff fileSize)
{
  if (static_cast<uint64_t>(n) * sizeof(double) > bytesRemaining(in, fileSize)) {
    return false;
  }
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
    O2BVHSurfaceSolid::PlanarBoundaryCurve curve;
    if (!parseBSplineEdge(edge.params, curve)) {
      return false;
    }
    // Evaluate the curve rather than read its first and last poles. Those coincide with the
    // endpoints only for a *clamped* knot vector, and nothing here or in Curve2D::valid() checks
    // clamping -- a periodic or unclamped spline (what OCC writes for a tube-tube seam before
    // SetNotPeriodic) has off-curve poles. The kernel stopped assuming this when K1 was fixed; the
    // loader did not, so the two layers measured the same wire join between different points and
    // could disagree in both directions: a sound wire rejected at load, or an open one accepted
    // (CodeReview_Fable_v2.md, N2).
    std::vector<surface::Vec2> poles;
    poles.reserve(curve.poles.size());
    for (const auto& pole : curve.poles) {
      poles.push_back({pole[0], pole[1]});
    }
    const surface::Curve2D evaluated =
      surface::Curve2D::makeBSpline(curve.degree, std::move(poles), curve.weights, curve.knots);
    const surface::Vec2 first = evaluated.startPoint();
    const surface::Vec2 last = evaluated.endPoint();
    start = {first.uCoord, first.vCoord};
    end = {last.uCoord, last.vCoord};
    return true;
  }
  return false;
}

/// The first fundamental form of the surface a sidecar record describes, so the reader can judge
/// its wires' join gaps in centimetres exactly as the kernel will. It reads the record's own
/// parameter block and defers to the kernel's closed forms; the layouts below are the argument
/// orders of the matching O2BVHSurfaceSolid::Add*Surface overload.
struct RecordMetric {
  uint32_t surfaceType = 0;
  const double* params = nullptr;

  static void evaluate(const void* context, const surface::Vec2& uv, double& gUU, double& gUV, double& gVV)
  {
    const auto& record = *static_cast<const RecordMetric*>(context);
    const double* p = record.params;
    switch (record.surfaceType) {
      case kPlane:
        surface::planeParametricMetric({p[3], p[4], p[5]}, {p[6], p[7], p[8]}, gUU, gUV, gVV);
        return;
      case kCylinder:
        surface::cylinderParametricMetric(p[9], gUU, gUV, gVV);
        return;
      case kCone: {
        // r(h) = radiusAtMin + slope * (h - heightMin), with slope from the two radii/heights
        const double slope = (p[10] - p[9]) / (p[12] - p[11]);
        surface::coneParametricMetric(p[9] + slope * (uv.vCoord - p[11]), slope, gUU, gUV, gVV);
        return;
      }
      case kSphere:
        surface::sphereParametricMetric(p[9], uv.vCoord, gUU, gUV, gVV);
        return;
      case kTorus:
        surface::torusParametricMetric(p[9], p[10], uv.vCoord, gUU, gUV, gVV);
        return;
      default:
        // an unknown type is rejected further down; the identity keeps this total meanwhile
        gUU = 1.;
        gUV = 0.;
        gVV = 1.;
        return;
    }
  }

  surface::ParametricMetric metric() const { return {&evaluate, this}; }
};

/// Convert a sidecar wire into the public PlanarBoundaryCurve loop expected by
/// AddCurvedPlanarSurface, validating that consecutive edges join end-to-start. Sets
/// \a anyArc when the wire carries at least one circular-arc edge.
///
/// \a metric is the surface's first fundamental form, built from the record's own parameters by
/// RecordMetric above. The join gap is judged as a 3D length in cm against exactly the band the
/// kernel will apply to the same wire moments later -- the loader used to apply a bare 1e-5 to a
/// mixed rad/cm separation, a different rule, which is finding S10 (and the reason the known
/// ST1829909_01 rejection happened at all: six joins drifting under 3e-5 rad on cylinder trims,
/// a negligible arc length, read as three times over tolerance).
///
/// \a joinTolerance is wireJoinToleranceFor(the sidecar's declared model tolerance): the model's
/// own statement of how closely its edges meet when it makes one, the extractor-precision
/// constant as the floor. Judging a declared-4.7e-4-cm model at the bare 1e-6 fallback re-created
/// the same class of rejection with correct units (surface 1006 again, a 5.41e-6 cm seam).
/// \a toleranceOrigin says which of the two the band is, for the diagnostic.
bool wireToCurves(const std::string& file, size_t surfaceIndex, const SidecarWire& wire,
                  std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>& curves, bool& anyArc,
                  const surface::ParametricMetric& metric, double joinTolerance, const char* toleranceOrigin)
{
  using Curve = O2BVHSurfaceSolid::PlanarBoundaryCurve;
  curves.clear();
  curves.reserve(wire.edges.size());
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
    const double joinGapSq = metric.distanceSq({end[0], end[1]}, {nextStart[0], nextStart[1]});
    if (joinGapSq > joinTolerance * joinTolerance) {
      ::Error("LoadSurfaceSolid",
              "%s: surface %zu: wire edge %zu end does not join the next edge start (gap %.3g cm, tolerance %.3g cm, "
              "%s)",
              file.c_str(), surfaceIndex, e, std::sqrt(joinGapSq), joinTolerance, toleranceOrigin);
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
                        std::vector<std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve>>& inners,
                        const surface::ParametricMetric& metric, double joinTolerance, const char* toleranceOrigin)
{
  bool haveOuter = false;
  bool anyArc = false; // quadric domains accept both line and arc trim edges
  for (const auto& wire : wires) {
    std::vector<O2BVHSurfaceSolid::PlanarBoundaryCurve> curves;
    if (!wireToCurves(file, surfaceIndex, wire, curves, anyArc, metric, joinTolerance, toleranceOrigin)) {
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

/// Permute a face's edge identities from the sidecar's wire order into the kernel's.
///
/// The sidecar lists wires in the order the CAD walk found them, each tagged outer or inner; the
/// Add*Surface API takes the outer wire as one argument and the inner wires as another, so the
/// kernel's flattened trim-curve index is "outer wire's curves, then the inner wires' in order".
/// Those two orders coincide only when the outer wire happens to come first, and an edge identity
/// that is off by one wire pairs the wrong two curves -- which is silent, because a wrong pairing
/// still produces a number. Permute rather than hope.
///
/// A face with no wire block (a quadric trimmed to its parametric rectangle) is left alone: its
/// identities are not anchored to any curve, and their order carries no meaning.
void reorderEdgeRefsToKernelOrder(const std::vector<SidecarWire>& wires, std::vector<unsigned int>& edgeIds,
                                  std::vector<unsigned char>& edgeFlags)
{
  size_t totalEdges = 0;
  for (const auto& wire : wires) {
    totalEdges += wire.edges.size();
  }
  if (wires.empty() || totalEdges != edgeIds.size()) {
    return;
  }
  // kernel offset of each sidecar wire: the outer wire first, then the inner wires in file order
  std::vector<size_t> kernelOffset(wires.size(), 0);
  size_t running = 0;
  for (size_t w = 0; w < wires.size(); ++w) {
    if (wires[w].role == 0) {
      kernelOffset[w] = 0;
      running = wires[w].edges.size();
      break;
    }
  }
  for (size_t w = 0; w < wires.size(); ++w) {
    if (wires[w].role != 0) {
      kernelOffset[w] = running;
      running += wires[w].edges.size();
    }
  }

  std::vector<unsigned int> permutedIds(edgeIds.size());
  std::vector<unsigned char> permutedFlags(edgeFlags.size());
  size_t sidecarOffset = 0;
  for (size_t w = 0; w < wires.size(); ++w) {
    for (size_t e = 0; e < wires[w].edges.size(); ++e) {
      permutedIds[kernelOffset[w] + e] = edgeIds[sidecarOffset + e];
      permutedFlags[kernelOffset[w] + e] = edgeFlags[sidecarOffset + e];
    }
    sidecarOffset += wires[w].edges.size();
  }
  edgeIds.swap(permutedIds);
  edgeFlags.swap(permutedFlags);
}

} // namespace

bool LoadSurfaceSolid(const std::string& file, O2BVHSurfaceSolid& solid)
{
  std::ifstream in(file, std::ios::binary);
  if (!in) {
    ::Error("LoadSurfaceSolid", "Cannot open surface sidecar file %s", file.c_str());
    return false;
  }

  in.seekg(0, std::ios::end);
  const std::streamoff fileSize = in.tellg();
  in.seekg(0, std::ios::beg);

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
  if (version < kSidecarVersionMin || version > kSidecarVersionMax) {
    ::Error("LoadSurfaceSolid", "%s: unsupported sidecar version %u (reader supports %u..%u)", file.c_str(), version,
            kSidecarVersionMin, kSidecarVersionMax);
    return false;
  }

  uint32_t nModelEdges = 0;
  if (version >= 2) {
    double modelTolerance = 0.;
    in.read(reinterpret_cast<char*>(&modelTolerance), sizeof(modelTolerance));
    if (!in) {
      ::Error("LoadSurfaceSolid", "%s: truncated version-2 header (no model tolerance)", file.c_str());
      return false;
    }
    solid.SetModelTolerance(modelTolerance);
    if (version >= 3 && !readU32(in, nModelEdges)) {
      ::Error("LoadSurfaceSolid", "%s: truncated version-3 header (no edge table size)", file.c_str());
      return false;
    }
  } else {
    ::Warning("LoadSurfaceSolid",
              "%s is a version-1 sidecar and states no model tolerance; assuming %g cm (the extractor's precision). "
              "Re-run the converter to record the model's own value.",
              file.c_str(), kSidecarV1FallbackTolerance);
    solid.SetModelTolerance(kSidecarV1FallbackTolerance);
  }

  // The wire-join acceptance band, fixed once from the header: the model's own declared tolerance
  // when it is looser than the extractor floor, the floor otherwise -- exactly the band the
  // kernel's Add*Surface will re-apply to the same wires (which is why the loader must not judge
  // more strictly than it, nor less).
  const double joinTolerance = surface::wireJoinToleranceFor(solid.GetModelTolerance());
  const char* toleranceOrigin = joinTolerance > surface::kWireJoinTolerance
                                  ? "declared by the model"
                                  : "the extractor-precision fallback";

  for (size_t s = 0; s < nSurfaces; ++s) {
    uint32_t surfaceType = 0, flags = 0, nParams = 0;
    if (!readU32(in, surfaceType) || !readU32(in, flags) || !readU32(in, nParams)) {
      ::Error("LoadSurfaceSolid", "%s: truncated surface record %zu", file.c_str(), s);
      return false;
    }
    std::vector<double> p;
    if (!readDoubles(in, p, nParams, fileSize)) {
      ::Error("LoadSurfaceSolid", "%s: truncated parameters of surface %zu", file.c_str(), s);
      return false;
    }

    // The wire block is self-describing; read it unconditionally.
    uint32_t nWires = 0;
    if (!readU32(in, nWires)) {
      ::Error("LoadSurfaceSolid", "%s: truncated wire count of surface %zu", file.c_str(), s);
      return false;
    }
    // 8 bytes of header per wire is the floor, so a count beyond that cannot be honest
    if (static_cast<uint64_t>(nWires) * 8u > bytesRemaining(in, fileSize)) {
      ::Error("LoadSurfaceSolid", "%s: surface %zu claims %u wires, more than the file holds", file.c_str(), s, nWires);
      return false;
    }
    std::vector<SidecarWire> wires(nWires);
    for (auto& wire : wires) {
      uint32_t nEdges = 0;
      if (!readU32(in, wire.role) || !readU32(in, nEdges)) {
        ::Error("LoadSurfaceSolid", "%s: truncated wire header in surface %zu", file.c_str(), s);
        return false;
      }
      if (static_cast<uint64_t>(nEdges) * 8u > bytesRemaining(in, fileSize)) {
        ::Error("LoadSurfaceSolid", "%s: surface %zu claims %u wire edges, more than the file holds", file.c_str(), s,
                nEdges);
        return false;
      }
      wire.edges.resize(nEdges);
      for (auto& edge : wire.edges) {
        uint32_t nCurveParams = 0;
        if (!readU32(in, edge.curveType) || !readU32(in, nCurveParams) ||
            !readDoubles(in, edge.params, nCurveParams, fileSize)) {
          ::Error("LoadSurfaceSolid", "%s: truncated edge record in surface %zu", file.c_str(), s);
          return false;
        }
      }
    }

    // Version 3: the face's boundary edge identities, in the sidecar's own wire order.
    std::vector<unsigned int> edgeIds;
    std::vector<unsigned char> edgeFlags;
    if (version >= 3) {
      uint32_t nEdgeRefs = 0;
      if (!readU32(in, nEdgeRefs)) {
        ::Error("LoadSurfaceSolid", "%s: truncated edge identity count of surface %zu", file.c_str(), s);
        return false;
      }
      if (static_cast<uint64_t>(nEdgeRefs) * 5u > bytesRemaining(in, fileSize)) {
        ::Error("LoadSurfaceSolid", "%s: surface %zu claims %u edge identities, more than the file holds",
                file.c_str(), s, nEdgeRefs);
        return false;
      }
      edgeIds.resize(nEdgeRefs);
      edgeFlags.resize(nEdgeRefs);
      for (uint32_t e = 0; e < nEdgeRefs; ++e) {
        uint32_t edgeId = 0;
        uint8_t edgeFlag = 0;
        if (!readU32(in, edgeId) || !readU8(in, edgeFlag)) {
          ::Error("LoadSurfaceSolid", "%s: truncated edge identity %u of surface %zu", file.c_str(), e, s);
          return false;
        }
        if (nModelEdges > 0 && edgeId >= nModelEdges) {
          ::Error("LoadSurfaceSolid", "%s: surface %zu edge identity %u is %u, outside the model's %u edge(s)",
                  file.c_str(), s, e, edgeId, nModelEdges);
          return false;
        }
        edgeIds[e] = edgeId;
        edgeFlags[e] = edgeFlag;
      }
    }

    const bool innerWall = (flags & kFlagInnerWall) != 0;
    const RecordMetric recordMetric{surfaceType, p.data()};
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
          if (!wireToCurves(file, s, wire, curves, anyArc, recordMetric.metric(), joinTolerance, toleranceOrigin)) {
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
          if (!collectQuadricTrim(file, s, wires, outer, inners, recordMetric.metric(), joinTolerance, toleranceOrigin)) {
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
          if (!collectQuadricTrim(file, s, wires, outer, inners, recordMetric.metric(), joinTolerance, toleranceOrigin)) {
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
          if (!collectQuadricTrim(file, s, wires, outer, inners, recordMetric.metric(), joinTolerance, toleranceOrigin)) {
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
          if (!collectQuadricTrim(file, s, wires, outer, inners, recordMetric.metric(), joinTolerance, toleranceOrigin)) {
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
    if (!edgeIds.empty()) {
      reorderEdgeRefsToKernelOrder(wires, edgeIds, edgeFlags);
      solid.SetSurfaceBoundaryEdges(static_cast<int>(s), edgeIds, edgeFlags);
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
