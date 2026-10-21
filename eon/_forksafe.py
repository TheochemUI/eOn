"""Keep a process usable after it forks with an OpenBLAS pool running.

OpenBLAS 0.3.30 registers a pthread_atfork prepare handler that shuts the
pool down and can leave server_lock held. The next factorization in the
parent then waits inside blas_thread_init. sh launches eonclient with
os.fork, and subprocess takes the same path when it cannot vfork, so both
hit it. Drop to one thread before that handler runs, and put the previous
count back in the parent. The child keeps one thread until it execs.
"""

from __future__ import annotations

import ctypes
import os
import sys
import threading

_installed = False


def _library_names():
    if sys.platform == "win32":
        return ("libopenblas.dll", "openblas.dll")
    if sys.platform == "darwin":
        return ("libopenblas.dylib", "libopenblas.0.dylib")
    return ("libopenblas.so.0", "libopenblas.so")


def openblas_library():
    """The OpenBLAS shared library this process can call, or None."""
    for name in _library_names():
        try:
            return ctypes.CDLL(name)
        except OSError:
            continue
    return None


def install() -> None:
    """Register the fork handlers once. A missing OpenBLAS is a no-op."""
    global _installed
    if _installed or not hasattr(os, "register_at_fork"):
        _installed = True
        return
    lib = openblas_library()
    set_threads = getattr(lib, "openblas_set_num_threads", None) if lib else None
    get_threads = getattr(lib, "openblas_get_num_threads", None) if lib else None
    if set_threads is None or get_threads is None:
        _installed = True
        return
    set_threads.argtypes = [ctypes.c_int]
    set_threads.restype = None
    get_threads.restype = ctypes.c_int
    saved = threading.local()

    def before() -> None:
        # Runs before OpenBLAS's own prepare handler. Registering after the
        # library has loaded makes this the first before hook.
        count = int(get_threads())
        saved.count = count
        if count != 1:
            set_threads(1)

    def after_parent() -> None:
        count = int(getattr(saved, "count", 1))
        if count != 1:
            set_threads(count)

    def after_child() -> None:
        set_threads(1)

    os.register_at_fork(
        before=before,
        after_in_parent=after_parent,
        after_in_child=after_child,
    )
    _installed = True
