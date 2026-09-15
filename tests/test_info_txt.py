"""info.txt writes through pathlib under path_results."""

import configparser

import pytest

pytest.importorskip("readcon")

from eon.fileio import info_txt_path, write_info_txt


class _Cfg:
    def __init__(self, path_results):
        self.path_results = path_results


def test_write_info_txt_creates_parent(tmp_path):
    cfg = _Cfg(tmp_path / "results")
    parser = configparser.RawConfigParser()
    parser.add_section("Simulation Information")
    parser.set("Simulation Information", "current_state", "3")
    write_info_txt(cfg, parser)
    path = info_txt_path(cfg)
    assert path.is_file()
    text = path.read_text()
    assert "current_state" in text
    assert "3" in text
