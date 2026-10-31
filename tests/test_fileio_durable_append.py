"""A torn table tail is dropped, and a failed append is cut back."""

import builtins
import logging
import os
from contextlib import contextmanager
from io import open as io_open

import pytest

from eon.fileio import DurableTable


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
        import io as pyio

        pyio.open = open_hook
        yield
    finally:
        builtins.open = real_open
        os.fdopen = real_fdopen
        import io as pyio

        pyio.open = io_open


def _good_line():
    return "% 7d %16.5f\n" % (0, -1.0)


def _torn_tail(kind):
    whole = "% 7d %16.5f\n" % (9, 3.33333)
    if kind == "newline":
        return whole.rstrip("\n")
    return "3.33333\n"


def _table(path):
    return DurableTable(path, header_lines=0, min_fields=2, exact=True, label="state_table demo")


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_tail_is_dropped(tmp_path, caplog, kind):
    path = tmp_path / "state_table"
    path.write_text(_good_line() + _torn_tail(kind))
    caplog.set_level(logging.WARNING, logger="io")
    lines = _table(path).lines()
    assert lines == [_good_line().rstrip("\n")]
    assert "torn final line" in caplog.text
    assert "3.33333" in path.read_text()


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_append_cuts_off_a_torn_tail(tmp_path, kind):
    path = tmp_path / "state_table"
    good = _good_line()
    path.write_text(good + _torn_tail(kind))
    added = "% 7d %16.5f\n" % (1, -2.0)
    _table(path).append(added)
    text = path.read_text()
    assert "3.33333" not in text
    assert text.startswith(good)
    assert text.endswith(added)


def test_failed_append_is_cut_back(tmp_path):
    path = tmp_path / "state_table"
    path.write_text(_good_line())
    original = path.read_text()
    with disk_full():
        with pytest.raises(OSError):
            _table(path).append("% 7d %16.5f\n" % (1, -2.0))
    assert path.read_text() == original
    assert len(path.read_bytes()) == len(original.encode())


def test_failed_append_keeps_the_intact_prefix(tmp_path):
    path = tmp_path / "state_table"
    good = _good_line()
    path.write_text(good + _torn_tail("newline"))
    with disk_full():
        with pytest.raises(OSError):
            _table(path).append("% 7d %16.5f\n" % (1, -2.0))
    assert path.read_text() == good


def test_append_fsyncs(tmp_path, monkeypatch):
    path = tmp_path / "state_table"
    path.write_text(_good_line())
    seen = []
    monkeypatch.setattr(os, "fsync", lambda fd: seen.append(fd))
    added = "% 7d %16.5f\n" % (1, -2.0)
    _table(path).append(added)
    assert seen
    assert path.read_text() == _good_line() + added


def test_minimum_field_count_keeps_a_longer_line(tmp_path, caplog):
    path = tmp_path / "searches.txt"
    path.write_text("h1\nh2\n1 2 3 4 5 6 7 8 extra words\n")
    caplog.set_level(logging.WARNING, logger="io")
    lines = DurableTable(path, header_lines=2, min_fields=8, exact=False).lines()
    assert lines[-1].split()[:8] == ["1", "2", "3", "4", "5", "6", "7", "8"]
    assert "torn final line" not in caplog.text


def test_too_few_fields_is_torn(tmp_path, caplog):
    path = tmp_path / "searches.txt"
    path.write_text("h1\nh2\n1 2 3\n")
    caplog.set_level(logging.WARNING, logger="io")
    lines = DurableTable(path, header_lines=2, min_fields=8, exact=False).lines()
    assert lines == ["h1", "h2"]
    assert "torn final line" in caplog.text


def test_short_header_is_kept(tmp_path, caplog):
    path = tmp_path / "processtable"
    path.write_text("proc header")
    caplog.set_level(logging.WARNING, logger="io")
    lines = DurableTable(path, header_lines=1, min_fields=9, exact=True).lines()
    assert lines == ["proc header"]
    assert "torn final line" not in caplog.text
