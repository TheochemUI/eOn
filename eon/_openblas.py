"""Keep this process on one OpenBLAS thread.

OpenBLAS's pthread build runs a server thread and guards it with
``server_lock``. Forking eonclient keeps only the calling thread, so the
server thread is gone and ``server_lock`` can stay held. The next
``numpy.linalg.solve`` then waits forever in ``blas_thread_init``. One
thread never starts that server. A client this process execs starts its
own pool; the limit is not an environment variable.
"""

from __future__ import annotations

import ctypes

_SONAMES = (
    "libopenblas.so.0",
    "libopenblas.so",
    "libcblas.so.3",
    "libblas.so.3",
)

_applied = False


def use_one_openblas_thread() -> None:
    """Ask every loaded OpenBLAS for a single thread. Missing libs are skipped."""
    global _applied
    if _applied:
        return
    _applied = True
    for name in _SONAMES:
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
