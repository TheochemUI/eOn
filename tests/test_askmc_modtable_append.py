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

"""A stopped modified-table append leaves the previous rows."""

import logging

import pytest

from eon.askmc import ASKMC, mod_processtable_header, processtable_line


class _State:
    def __init__(self, path):
        self.path = str(path)
        self.number = 0


def test_stopped_append_leaves_the_modified_table(tmp_path):
    state_dir = tmp_path / "state 0"
    state_dir.mkdir()
    table = state_dir / "askmc_processtable"
    ask = ASKMC(
        kT=0.025,
        states=None,
        confidence=0.9,
        alpha=1.5,
        gamma=2.0,
        barrier_test_on=False,
        connection_test_on=False,
        sb_recycling_on=False,
        path_root=str(tmp_path),
        thermal_window=20.0,
        recycle_path=None,
    )
    state = _State(state_dir)
    ask.append_modified_process_table(
        state, 1, 0.2, 1.0e12, 2, 0.0, 1.0e12, 0.2, 1.0, 1
    )
    original = table.read_text()
    with disk_full():
        with pytest.raises(OSError):
            ask.append_modified_process_table(
                state, 4, 0.4, 1.0e12, 3, 0.0, 1.0e12, 0.4, 0.2, 1
            )
    assert table.read_text() == original
    loaded = ask.get_modified_process_table(state)
    assert 4 not in loaded
    assert 1 in loaded


def _ask(tmp_path):
    return ASKMC(
        kT=0.025,
        states=None,
        confidence=0.9,
        alpha=1.5,
        gamma=2.0,
        barrier_test_on=False,
        connection_test_on=False,
        sb_recycling_on=False,
        path_root=str(tmp_path),
        thermal_window=20.0,
        recycle_path=None,
    )


def _torn_tail(kind):
    whole = processtable_line % (9, 1.23456, 1.0e12, 4, 0.0, 1.0e12, 0.3, 0.5, 0)
    if kind == "newline":
        return whole.rstrip("\n")
    return "9 1.23456\n"


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_tail_is_dropped(tmp_path, caplog, kind):
    state_dir = tmp_path / "state 0"
    state_dir.mkdir()
    ask = _ask(tmp_path)
    state = _State(state_dir)
    ask.append_modified_process_table(
        state, 1, 0.2, 1.0e12, 2, 0.0, 1.0e12, 0.2, 1.0, 1
    )
    table = state_dir / "askmc_processtable"
    table.write_bytes(table.read_bytes() + _torn_tail(kind).encode())
    caplog.set_level(logging.WARNING, logger="askmc")
    loaded = ask.get_modified_process_table(state)
    assert 1 in loaded
    assert 9 not in loaded
    assert "torn final line" in caplog.text


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_append_cuts_off_a_torn_tail(tmp_path, kind):
    state_dir = tmp_path / "state 0"
    state_dir.mkdir()
    ask = _ask(tmp_path)
    state = _State(state_dir)
    ask.append_modified_process_table(
        state, 1, 0.2, 1.0e12, 2, 0.0, 1.0e12, 0.2, 1.0, 1
    )
    table = state_dir / "askmc_processtable"
    table.write_bytes(table.read_bytes() + _torn_tail(kind).encode())
    ask.append_modified_process_table(
        state, 7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0
    )
    text = table.read_text()
    assert "1.23456" not in text
    assert text.startswith(mod_processtable_header)
    assert text.endswith(processtable_line % (7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0))
    loaded = ask.get_modified_process_table(state)
    assert 1 in loaded
    assert 7 in loaded
    assert 9 not in loaded
