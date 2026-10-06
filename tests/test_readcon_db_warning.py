"""A missing libreadcon_db.so is one warning, and a foreign kdb import is not."""

from __future__ import annotations

import logging
from pathlib import Path

from eon import concorpus, process_catalog


def test_missing_library_warns_once_at_info(monkeypatch, caplog, tmp_path: Path):
    monkeypatch.setattr(concorpus, "_unavailable", False)
    monkeypatch.setattr(concorpus, "_corpora", {})
    dlopen_error = (
        'dlopen("libreadcon_db.so", 2): cannot open shared object file'
    )

    def refuse(name, *args, **kwargs):
        raise OSError(dlopen_error)

    import ctypes

    monkeypatch.setattr(ctypes, "CDLL", refuse)

    import builtins

    real_import = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "readcon_db":
            raise ImportError("No module named 'readcon_db'", name="readcon_db")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", blocked)
    caplog.set_level(logging.INFO, logger="eon.concorpus")
    text = "hello\n"
    concorpus.mirror_con_text(tmp_path / "a.con", text)
    concorpus.mirror_con_text(tmp_path / "b.con", text)
    warnings = [
        record
        for record in caplog.records
        if record.name == "eon.concorpus" and record.levelno >= logging.WARNING
    ]
    assert len(warnings) == 1
    message = warnings[0].getMessage()
    assert dlopen_error in message
    assert warnings[0].levelno == logging.WARNING
    info_lines = [
        record.getMessage()
        for record in caplog.records
        if record.levelno >= logging.INFO
    ]
    assert any(dlopen_error in line for line in info_lines)
    debug_only = [
        record
        for record in caplog.records
        if record.levelno == logging.DEBUG and "readcon_db" in record.getMessage()
    ]
    assert debug_only == []
    assert "Python module kdb not found" not in caplog.text


def test_foreign_kdb_import_is_not_a_missing_module(monkeypatch, caplog):
    import builtins

    real_import = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "amsel":
            raise ImportError("No module named 'irc'", name="irc")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", blocked)
    caplog.set_level(logging.INFO, logger="kdb")

    class Config:
        kdb_path = "kdb"

    assert process_catalog._open_store(Config()) is None
    text = caplog.text
    assert "irc" in text
    assert "Python module kdb not found" not in text
    assert "not found" not in text.lower()
    assert "amsel is not installed" not in text
