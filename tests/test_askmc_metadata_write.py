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

"""askmc_data.txt keeps its previous counters when a write stops."""

import pytest

from eon.askmc import ASKMC


def _ask(path):
    return ASKMC(
        kT=0.025,
        states=None,
        confidence=0.9,
        alpha=1.5,
        gamma=2.0,
        barrier_test_on=False,
        connection_test_on=False,
        sb_recycling_on=False,
        path_root=str(path),
        thermal_window=20.0,
        recycle_path=None,
    )


def test_metadata_keeps_the_previous_counters(tmp_path):
    ask = _ask(tmp_path)
    ask.save_askmc_metadata(2, 4)
    original = (tmp_path / "askmc_data.txt").read_text()
    with disk_full():
        with pytest.raises(OSError):
            ask.save_askmc_metadata(3, 5)
    assert (tmp_path / "askmc_data.txt").read_text() == original
    assert "sb_check_count = 2" in original
    assert "num_rate_changes = 4" in original
