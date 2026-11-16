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

"""Recycling metadata keeps the previous file when a write stops."""

import pytest
from types import SimpleNamespace

from eon.recycling import Recycling, SB_Recycling


def test_superbasin_metadata_keeps_the_previous_file(tmp_path):
    recycler = SB_Recycling.__new__(SB_Recycling)
    recycler.path = str(tmp_path)
    recycler.sb_state_nums = [[0, 1]]
    recycler.in_progress = True
    recycler.write_metadata()
    original = (tmp_path / "recycling_data.txt").read_text()
    recycler.sb_state_nums = [[0, 2]]
    with disk_full():
        with pytest.raises(OSError):
            recycler.write_metadata()
    assert (tmp_path / "recycling_data.txt").read_text() == original


def test_saddle_metadata_keeps_the_previous_file(tmp_path):
    recycler = Recycling.__new__(Recycling)
    recycler.metadata_path = str(tmp_path / "meta dir" / "recycling_metadata.txt")
    (tmp_path / "meta dir").mkdir()
    recycler.ref_state = SimpleNamespace(number=1)
    recycler.process_number = 0
    recycler.num_procs = 2
    recycler.moved = [0]
    recycler.unmoved = [1]
    recycler.process_atoms = [0]
    recycler.write_recycling_metadata()
    original = (tmp_path / "meta dir" / "recycling_metadata.txt").read_text()
    recycler.process_number = 4
    with disk_full():
        with pytest.raises(OSError):
            recycler.write_recycling_metadata()
    text = (tmp_path / "meta dir" / "recycling_metadata.txt").read_text()
    assert text == original
