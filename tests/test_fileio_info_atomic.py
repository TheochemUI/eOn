"""A failed info.txt write must leave the previous file."""

from __future__ import annotations

import configparser

import pytest

from eon import fileio as io


def test_failed_info_write_keeps_previous_text(tmp_path, monkeypatch):
    parser = configparser.ConfigParser()
    parser.add_section("Main")
    parser.set("Main", "job", "akmc")

    class Cfg:
        path_results = str(tmp_path)

    io.write_info_txt(Cfg(), parser)
    original = (tmp_path / "info.txt").read_text()
    assert "job" in original

    def boom(self, fh):
        fh.write("partial")
        raise RuntimeError("disk full")

    monkeypatch.setattr(configparser.ConfigParser, "write", boom)
    with pytest.raises(RuntimeError, match="disk full"):
        io.write_info_txt(Cfg(), parser)
    assert (tmp_path / "info.txt").read_text() == original
