"""Forking after OpenBLAS has a pool must not stick server_lock.

sh launches eonclient with os.fork. The one-thread limit in conftest hides
that. This test turns the limit off in a fresh interpreter, imports eon so
the library handler is installed, and forks. Without the handler the
following slogdet waits in blas_thread_init.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

_CHILD = r"""
import ctypes
import os
import sys

import eon
from eon._forksafe import openblas_library

lib = openblas_library()
if lib is None:
    sys.stderr.write("no OpenBLAS\n")
    sys.exit(2)
get_threads = lib.openblas_get_num_threads
get_threads.restype = ctypes.c_int
set_threads = lib.openblas_set_num_threads
set_threads.argtypes = [ctypes.c_int]
set_threads.restype = None
set_threads(4)
if int(get_threads()) < 2:
    sys.stderr.write("OpenBLAS stayed at one thread\n")
    sys.exit(2)
import numpy as np
rng = np.random.default_rng(0)
matrix = rng.standard_normal((200, 200))
vector = rng.standard_normal(200)
np.linalg.solve(matrix, vector)
np.linalg.slogdet(matrix)
pid = os.fork()
if pid == 0:
    os.execv("/bin/true", ["true"])
os.waitpid(pid, 0)
np.linalg.solve(matrix, vector)
sign, _logdet = np.linalg.slogdet(matrix)
if sign == 0:
    sys.stderr.write("slogdet sign was 0\n")
    sys.exit(3)
if int(get_threads()) < 2:
    sys.stderr.write("parent pool was left at one thread\n")
    sys.exit(4)
"""


def test_fork_after_openblas_threads_does_not_deadlock():
    if not hasattr(os, "fork") or not os.path.exists("/bin/true"):
        pytest.skip("this platform has no os.fork and /bin/true")
    from eon._forksafe import openblas_library

    if openblas_library() is None:
        pytest.skip("OpenBLAS is not loadable, so the fork deadlock cannot be reproduced")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (str(ROOT), env.get("PYTHONPATH", "")) if part
    )
    for key in (
        "OPENBLAS_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env.pop(key, None)
    proc = subprocess.run(
        [sys.executable, "-c", _CHILD],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    if proc.returncode == 2:
        pytest.skip(proc.stderr.strip() or "OpenBLAS stayed at one thread")
    assert proc.returncode == 0, proc.stderr + proc.stdout
