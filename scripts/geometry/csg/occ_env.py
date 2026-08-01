"""Make `import OCC` work regardless of which interpreter started us.

pythonOCC on this machine is built against the alibuild Python 3.10; the system 3.12 cannot
import it, and the failure mode (a bare ImportError) has cost time before. `occtOracle.py` and
`runOracleGate.py` solve this by *launching* the right interpreter as a subprocess; a library
that wants to be importable from either side needs the mirror image of that, so this module
re-executes the current script under the correct interpreter with the correct environment and
then never has to think about it again.

Call `ensure_occ()` as the first statement of `main()`, before any `from OCC...` import.
"""

import os
import sys
from pathlib import Path

# Same locations `runOracleGate.py` uses; kept in sync deliberately rather than imported, because
# that script is owned by Stream E and this package must not add a dependency on it.
_SW = Path(os.environ.get("ALIBUILD_ARCH_ROOT", Path.home() / "alisw/sw/ubuntu2404_aarch64"))
OCC_PYTHON = _SW / "Python/latest/bin/python3.10"
OCC_ENV = {
    "PYTHONPATH": f"{_SW}/pythonOCC/latest/lib/python3.10/site-packages:"
                  f"{_SW}/Python-modules/latest/lib/python3.10/site-packages",
    "LD_LIBRARY_PATH": f"{_SW}/OCCT/latest/lib:{_SW}/Python/latest/lib",
}

_GUARD = "O2_CSG_OCC_REEXEC"


def have_occ() -> bool:
    try:
        import OCC  # noqa: F401
        return True
    except Exception:
        return False


def ensure_occ() -> None:
    """Re-exec this process under the pythonOCC interpreter if OCC is not importable."""
    if have_occ():
        return
    if os.environ.get(_GUARD):
        raise SystemExit(
            f"cannot import OCC even under {OCC_PYTHON}; check the pythonOCC installation")
    if not OCC_PYTHON.exists():
        raise SystemExit(f"pythonOCC interpreter not found: {OCC_PYTHON}")
    env = dict(os.environ)
    env.update(OCC_ENV)
    env[_GUARD] = "1"
    script = Path(sys.argv[0]).resolve()
    os.execve(str(OCC_PYTHON), [str(OCC_PYTHON), str(script)] + sys.argv[1:], env)
