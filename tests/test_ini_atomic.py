"""A failed write_ini leaves the previous config.ini."""

from __future__ import annotations

import configparser
from pathlib import Path

import pytest

from eon_schema.config import write_ini


def test_failed_write_leaves_existing_ini(tmp_path: Path, monkeypatch):
    dest = tmp_path / "config.ini"
    original = "[Main]\njob = point\n"
    dest.write_text(original, encoding="utf-8")

    def boom(self, fh):
        fh.write("partial")
        raise RuntimeError("disk full")

    monkeypatch.setattr(configparser.ConfigParser, "write", boom)
    with pytest.raises(RuntimeError, match="disk full"):
        write_ini(dest, {"Main": {"job": "akmc"}})
    assert dest.read_text(encoding="utf-8") == original
