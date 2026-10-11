"""Shared fixtures for the eOn test suite.

The integration tests under tests/ shell out to the built binaries. Those
are products of a meson build rather than of the python package, so they
are only present once one has been run and installed. Resolving them at
import time turns their absence into a collection error for the whole
module; resolving them here turns it into a skip that names what is
missing.
"""

import ctypes
import os
import shutil
from pathlib import Path

import pytest

# Fork after threads. OpenBLAS's pthread server holds server_lock. This
# test process later forks eonclient through sh, which keeps only the
# calling thread, so the lock stays held and numpy.linalg.solve waits in
# blas_thread_init. One thread here, before collection imports numpy.
# Not an environment variable: the client process keeps its own pool.
def _one_openblas_thread() -> None:
    for name in ("libopenblas.so.0", "libopenblas.so", "libcblas.so.3", "libblas.so.3"):
        try:
            lib = ctypes.CDLL(name)
        except OSError:
            continue
        setter = getattr(lib, "openblas_set_num_threads", None)
        if setter is None:
            continue
        setter.argtypes = [ctypes.c_int]
        setter.restype = None
        setter(1)


_one_openblas_thread()


def require_binary(name):
    """Return an sh.Command for *name*, or skip the caller.

    Looks at EON_<NAME>_BIN first so a build directory can be pointed at
    without installing, then falls back to PATH.
    """
    import sh

    override = os.environ.get("EON_%s_BIN" % name.upper())
    if override:
        if not Path(override).is_file():
            pytest.skip("%s names %s, which does not exist"
                        % ("EON_%s_BIN" % name.upper(), override))
        return sh.Command(override)

    found = shutil.which(name)
    if found is None:
        pytest.skip(
            "%s is not on PATH and EON_%s_BIN is unset; build and install the "
            "client, or point the variable at the binary in the build "
            "directory" % (name, name.upper())
        )
    return sh.Command(found)


@pytest.fixture(scope="session")
def eonclient():
    return require_binary("eonclient")


@pytest.fixture(scope="session")
def eon():
    return require_binary("eon")
