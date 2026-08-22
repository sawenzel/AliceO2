#!/usr/bin/env python3
"""Make an exact-surface geom.C loadable through o2-sim's external-geometry mechanism.

WORKAROUND, not a fix. `o2::base::loadCADGeometryHook` (Detectors/Base/src/CADGeometryUtils.cxx)
JITs the macro inside a unique namespace `o2_cadgeom_<TAG>_<N>`, hoisting only lines that start
with '#' to global scope. `O2_CADtoTGeo.py::emit_cpp_prelude(exact_surfaces=True)` writes a
hand-rolled forward declaration

    namespace o2 { namespace base { bool LoadSurfaceSolid(...); } }

which the wrapper turns into `o2_cadgeom_X::o2::base`. Every later `o2::base::O2BVHSurfaceSolid`
then resolves against that nested `o2` instead of the global one, and the macro fails to
compile with

    error: no type named 'O2BVHSurfaceSolid' in namespace 'o2_cadgeom_IRIS_0::o2::base'

so no exact-surface CAD geometry can be loaded into o2-sim at all.

This script replaces the forward declaration with a real #include of the header that declares
the same symbol. Being a '#' line it is hoisted to global scope, where `o2::base` is the real
one. The header is installed ($O2_ROOT/include/DetectorsBase/O2SurfaceSolidIO.h) and the macro
already adds that include path.

Usage: patch_exact_macro.py <geom.C> [...]   (idempotent)
"""
import sys

BLOCK = """// O2SurfaceSolidIO.h is not part of the ROOT dictionary module; declare the loader
// prototype directly (the symbol resolves from libO2DetectorsBase).
namespace o2
{
namespace base
{
bool LoadSurfaceSolid(const std::string& file, O2BVHSurfaceSolid& solid);
} // namespace base
} // namespace o2
"""

REPLACEMENT = """// PATCHED by integration_demo/patch_exact_macro.py: the emitted forward declaration is
// nested by the JIT namespace wrapper in CADGeometryUtils.cxx and shadows ::o2. A '#include'
// is hoisted to global scope by that wrapper, so it declares the right symbol.
#include "DetectorsBase/O2SurfaceSolidIO.h"
"""


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    rc = 0
    for path in sys.argv[1:]:
        text = open(path).read()
        if REPLACEMENT.splitlines()[-1] in text:
            print(f"{path}: already patched")
            continue
        if BLOCK not in text:
            if "O2BVHSurfaceSolid" not in text:
                print(f"{path}: no exact-surface prelude, nothing to do")
            else:
                print(f"{path}: ERROR prelude not recognised -- converter output changed")
                rc = 1
            continue
        open(path, "w").write(text.replace(BLOCK, REPLACEMENT))
        print(f"{path}: patched")
    return rc


if __name__ == "__main__":
    sys.exit(main())
