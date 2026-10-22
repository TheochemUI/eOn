"""Helpers for a write that dies after a short prefix."""

import builtins
import os
from contextlib import contextmanager
from io import open as io_open


@contextmanager
def disk_full(limit=24):
    real_open = builtins.open
    real_fdopen = os.fdopen

    class _Partial:
        def __init__(self, handle):
            self._handle = handle
            self._left = limit

        def write(self, data):
            if self._left <= 0:
                raise OSError("disk full")
            chunk = data[: self._left]
            if chunk:
                self._handle.write(chunk)
                self._handle.flush()
            self._left = 0
            raise OSError("disk full")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self._handle.close()
            return False

        def close(self):
            self._handle.close()

        def flush(self):
            self._handle.flush()

        def __getattr__(self, name):
            return getattr(self._handle, name)

    def _writes(mode):
        text = mode if isinstance(mode, str) else ""
        return any(flag in text for flag in ("w", "a", "x", "+"))

    def open_hook(file, mode="r", *args, **kwargs):
        handle = real_open(file, mode, *args, **kwargs)
        if _writes(mode):
            return _Partial(handle)
        return handle

    def fdopen_hook(fd, mode="r", *args, **kwargs):
        handle = real_fdopen(fd, mode, *args, **kwargs)
        if _writes(mode):
            return _Partial(handle)
        return handle

    builtins.open = open_hook
    os.fdopen = fdopen_hook
    try:
        # pathlib.Path.write_text goes through io.open.
        import io as pyio

        pyio.open = open_hook
        yield
    finally:
        builtins.open = real_open
        os.fdopen = real_fdopen
        import io as pyio

        pyio.open = io_open

"""A stopped search-log append leaves the previous lines."""

import pytest
from types import SimpleNamespace

from eon.akmcstate import AKMCState


def test_stopped_search_append_leaves_the_log(tmp_path):
    path = tmp_path / "search dir" / "searches.txt"
    path.parent.mkdir()
    path.write_text("header\n----\n       1     random    0.10000    1.00000        10        20         0    good\n")
    state = AKMCState.__new__(AKMCState)
    state.search_result_path = str(path)
    state.config = SimpleNamespace(akmc_server_side_process_search=False)
    state.info = SimpleNamespace(set=lambda *args, **kwargs: None)
    original = path.read_text()
    result = {
        "wuid": 2,
        "type": "random",
        "results": {
            "barrier_reactant_to_product": 0.2,
            "displacement_saddle_distance": 1.0,
            "force_calls_saddle": 3,
            "force_calls_minimization": 4,
            "force_calls_prefactors": 5,
        },
    }
    with disk_full():
        with pytest.raises(OSError):
            state.append_search_result(result, "good", None)
    assert path.read_text() == original
