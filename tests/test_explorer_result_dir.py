"""Saving one result twice does not fail when the directory exists."""

from types import SimpleNamespace

from eon import fileio as io
from eon.explorer import ProcessSearch


def test_saving_a_result_twice_keeps_the_files(tmp_path):
    search = ProcessSearch.__new__(ProcessSearch)
    search.config = SimpleNamespace(path_incomplete=str(tmp_path / "incomplete dir"))
    result = {"name": "0_1", "saddle.con": io.StringIO("frame\n")}
    search.save_result(result)
    search.save_result(result)
    text = (tmp_path / "incomplete dir" / "0_1" / "saddle.con").read_text()
    assert text == "frame\n"
