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

"""A stopped append leaves the previous state table intact."""

import pytest

from types import SimpleNamespace

from eon.statelist import StateList


def test_stopped_append_leaves_the_state_table(tmp_path):
    root = tmp_path / "states dir"
    (root / "0").mkdir(parents=True)
    table = root / "state_table"
    table.write_text("      0         -1.00000\n")
    cfg = SimpleNamespace(
        path_states=str(root),
        comp_eps_e=1.0e-3,
        comp_eps_r=0.1,
        comp_use_identical=False,
    )
    listing = StateList(object, config=cfg)
    original = table.read_text()
    with disk_full():
        with pytest.raises(OSError):
            listing.append_state_table(-2.0)
    assert table.read_text() == original
    listing.state_table = None
    assert listing.get_num_states() == 1
