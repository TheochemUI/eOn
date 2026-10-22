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

"""A stopped append leaves a parsable process table."""

import pytest

from types import SimpleNamespace

from eon.akmcstate import AKMCState


class _Info:
    def get(self, section, option, default=None):
        if option == "kT":
            return 0.025
        return default

    def set(self, *args, **kwargs):
        return None


def _state(path):
    state = AKMCState.__new__(AKMCState)
    state.number = 0
    state.proctable_path = str(path)
    state.procs = None
    state.info = _Info()
    state.statelist = SimpleNamespace(kT=0.025)
    state.config = SimpleNamespace(akmc_eq_rate=0.0)
    path.write_text(
        AKMCState.processtable_header
        + AKMCState.processtable_line % (1, 0.2, 1.0e12, 2, 0.0, 1.0e12, 0.2, 1.0, 0)
    )
    return state


def test_stopped_append_leaves_the_process_table(tmp_path):
    table = tmp_path / "processtable"
    state = _state(table)
    original = table.read_text()
    with disk_full():
        with pytest.raises(OSError):
            state.append_process_table(7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0)
    assert table.read_text() == original
    state.procs = None
    state.load_process_table(force=True)
    assert 7 not in state.procs
    assert 1 in state.procs
