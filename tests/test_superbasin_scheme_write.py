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

"""Superbasin scheme files keep the previous text when a write stops."""

import pytest
from types import SimpleNamespace

from eon.superbasinscheme import EnergyLevel, TransitionCounting


class _State:
    def __init__(self, number, path):
        self.number = number
        self.path = str(path)


def test_transition_counts_keep_the_previous_file(tmp_path):
    start_dir = tmp_path / "state 0"
    start_dir.mkdir()
    start = _State(0, start_dir)
    end = _State(1, tmp_path / "state 1")
    data = start_dir / "sb.txt"
    data.write_text("1 1\n")
    scheme = TransitionCounting.__new__(TransitionCounting)
    scheme.config = SimpleNamespace(sb_state_file="sb.txt")
    scheme.count = {start: {end: 2}}
    with disk_full():
        with pytest.raises(OSError):
            scheme.write_data()
    assert data.read_text() == "1 1\n"


def test_energy_levels_keep_the_previous_file(tmp_path):
    state_dir = tmp_path / "state 0"
    state_dir.mkdir()
    state = _State(0, state_dir)
    data = state_dir / "levels.txt"
    data.write_text("1.0\n")
    scheme = EnergyLevel.__new__(EnergyLevel)
    scheme.config = SimpleNamespace(sb_state_file="levels.txt")
    scheme.levels = {state: 1.5}
    with disk_full():
        with pytest.raises(OSError):
            scheme.write_data()
    assert data.read_text() == "1.0\n"
