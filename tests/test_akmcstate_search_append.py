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

import logging

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


class _Info:
    def __init__(self):
        self.store = {}

    def get(self, section, option, default=None):
        return self.store.get((section, option), default)

    def set(self, section, option, value):
        self.store[(section, option)] = value


def _result(wuid):
    return {
        "wuid": wuid,
        "type": "random",
        "results": {
            "barrier_reactant_to_product": 0.2,
            "displacement_saddle_distance": 1.0,
            "force_calls_saddle": 3,
            "force_calls_minimization": 4,
            "force_calls_prefactors": 5,
        },
    }


def _state(path):
    state = AKMCState.__new__(AKMCState)
    state.number = 0
    state.search_result_path = str(path)
    state.config = SimpleNamespace(akmc_server_side_process_search=False)
    state.info = _Info()
    return state


def _line(wuid, comment):
    return "%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (
        wuid, "random", 0.2, 1.0, 3, 4, 5, comment,
    )


def _torn_tail(kind):
    if kind == "newline":
        return _line(2, "TORNTAIL").rstrip("\n")
    return "2 TORNTAIL\n"


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_tail_is_dropped(tmp_path, caplog, kind):
    path = tmp_path / "search dir" / "searches.txt"
    path.parent.mkdir()
    good = "header\n----\n       1     random    0.10000    1.00000        10        20         0    good\n"
    path.write_text(good + _torn_tail(kind))
    state = _state(path)
    caplog.set_level(logging.WARNING, logger="state")
    assert state.get_number_of_searches() == 1
    assert "torn final line" in caplog.text


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_append_cuts_off_a_torn_tail(tmp_path, kind):
    path = tmp_path / "search dir" / "searches.txt"
    path.parent.mkdir()
    good = "header\n----\n       1     random    0.10000    1.00000        10        20         0    good\n"
    path.write_text(good + _torn_tail(kind))
    state = _state(path)
    state.append_search_result(_result(3), "kept", None)
    text = path.read_text()
    assert "TORNTAIL" not in text
    assert text.startswith("header\n----\n")
    assert text.endswith(_line(3, "kept"))
    state.info = _Info()
    assert state.get_number_of_searches() == 2
