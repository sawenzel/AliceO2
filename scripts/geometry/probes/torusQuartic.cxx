// Stream L, stage 6: reduce the missing torus crossing to the quartic root solve, on the REAL
// record and the REAL ray, and vary only how far away the ray starts.

#include "BoundedSurface.h"
#include "DetectorsBase/O2BVHSurfaceSolid.h"
#include "DetectorsBase/O2SurfaceSolidIO.h"

#include "TGeoManager.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2::base;
using o2::base::surface::solveQuarticReal;

int main(int argc, char** argv)
{
  std::string surfaces;
  int which = -1;
  double o[3] = {0, 0, 0}, d[3] = {1, 0, 0};
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--surfaces") {
      surfaces = argv[++i];
    } else if (a == "--surface") {
      which = std::atoi(argv[++i]);
    } else if (a == "--o") {
      for (int k = 0; k < 3; ++k) {
        o[k] = std::atof(argv[++i]);
      }
    } else if (a == "--d") {
      for (int k = 0; k < 3; ++k) {
        d[k] = std::atof(argv[++i]);
      }
    }
  }
  auto* manager = new TGeoManager("xL6", "torus quartic");
  auto* solid = new O2BVHSurfaceSolid("part");
  if (!LoadSurfaceSolid(surfaces, *solid)) {
    return 2;
  }
  solid->CloseShape(true);
  (void)manager;
  const auto& rec = solid->GetSurfaceRecords().at(which);
  const double R = rec.scalars[0], r = rec.scalars[1];
  const double C[3] = {rec.origin[0], rec.origin[1], rec.origin[2]};
  const double W[3] = {rec.axisA[0], rec.axisA[1], rec.axisA[2]};
  const double U[3] = {rec.axisB[0], rec.axisB[1], rec.axisB[2]};
  const double V[3] = {W[1] * U[2] - W[2] * U[1], W[2] * U[0] - W[0] * U[2], W[0] * U[1] - W[1] * U[0]};
  std::printf("surface %d: torus R=%.10g r=%.10g centre=(%.10g, %.10g, %.10g)\n", which, R, r, C[0],
              C[1], C[2]);
  std::printf("%12s %14s  roots (absolute t)\n", "start offset", "remaining");
  for (double shift : {0., 1., 10., 100., 200., 300., 350., 370., 374., 375.}) {
    double P[3];
    for (int k = 0; k < 3; ++k) {
      P[k] = o[k] + shift * d[k];
    }
    const double rel[3] = {P[0] - C[0], P[1] - C[1], P[2] - C[2]};
    const double lo[3] = {rel[0] * U[0] + rel[1] * U[1] + rel[2] * U[2],
                          rel[0] * V[0] + rel[1] * V[1] + rel[2] * V[2],
                          rel[0] * W[0] + rel[1] * W[1] + rel[2] * W[2]};
    const double ld[3] = {d[0] * U[0] + d[1] * U[1] + d[2] * U[2],
                          d[0] * V[0] + d[1] * V[1] + d[2] * V[2],
                          d[0] * W[0] + d[1] * W[1] + d[2] * W[2]};
    const double dd = ld[0] * ld[0] + ld[1] * ld[1] + ld[2] * ld[2];
    const double od = lo[0] * ld[0] + lo[1] * ld[1] + lo[2] * ld[2];
    const double oo = lo[0] * lo[0] + lo[1] * lo[1] + lo[2] * lo[2];
    const double K = R * R - r * r;
    const double E = ld[0] * ld[0] + ld[1] * ld[1];
    const double F = lo[0] * ld[0] + lo[1] * ld[1];
    const double G = lo[0] * lo[0] + lo[1] * lo[1];
    const double R4 = 4. * R * R;
    const double a4 = dd * dd;
    const double a3 = 4. * dd * od;
    const double a2 = 4. * od * od + 2. * dd * (oo + K) - R4 * E;
    const double a1 = 4. * od * (oo + K) - 2. * R4 * F;
    const double a0 = (oo + K) * (oo + K) - R4 * G;
    auto roots = solveQuarticReal(a4, a3, a2, a1, a0);
    std::printf("%12.6g %14.6g  n=%zu :", shift, 375.33923 - shift, roots.size());
    for (double x : roots) {
      std::printf(" %.10g", x + shift);
    }
    std::printf("\n             a4=%.17g a3=%.17g a2=%.17g a1=%.17g a0=%.17g\n", a4, a3, a2, a1, a0);
    // Which guard, exactly.
    {
      const double b = a3 / a4, cc = a2 / a4, dd2 = a1 / a4, e = a0 / a4;
      const double p = cc - 3. * b * b / 8.;
      const double q = dd2 - b * cc / 2. + b * b * b / 8.;
      const double rr2 = e - b * dd2 / 4. + b * b * cc / 16. - 3. * b * b * b * b / 256.;
      const double A2 = p, A1 = p * p / 4. - rr2, A0 = -q * q / 8.;
      const double cp = A1 - A2 * A2 / 3., cq = 2. * A2 * A2 * A2 / 27. - A2 * A1 / 3. + A0;
      std::vector<double> cr;
      o2::base::surface::solveDepressedCubic(cp, cq, cr);
      double res = 0.;
      for (double x : cr) {
        res = std::max(res, x - A2 / 3.);
      }
      std::printf("             p=%.10g q=%.10g r=%.10g | cubicRoots=%zu resolvent=%.10g %s\n", p, q,
                  rr2, cr.size(), res, res > 1.e-9 ? "" : "<== RESOLVENT GUARD FAILS: no roots");
      if (res > 1.e-9) {
        const double s = std::sqrt(2. * res), lin = s * q / (4. * res);
        std::printf("             quad discriminants: %.10g  %.10g %s\n", s * s - 4. * (p / 2. + res + lin),
                    s * s - 4. * (p / 2. + res - lin),
                    (s * s - 4. * (p / 2. + res + lin) < 0. && s * s - 4. * (p / 2. + res - lin) < 0.)
                      ? "<== BOTH COMPLEX: no roots"
                      : "");
      }
    }
  }
  return 0;
}
