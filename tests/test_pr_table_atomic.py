"""A killed process-table rewrite leaves the previous rows in place."""

import pytest

from eon.prstate import PRState


def test_killed_table_rewrite_keeps_the_previous_rows(tmp_path, monkeypatch):
    path = tmp_path / "processtable"
    path.write_text("OLD\n")
    state = PRState.__new__(PRState)
    state.procs = {1: {"product": 2, "product_energy": -1.5, "time": 3.0}}
    state.proctable_path = str(path)
    state.processtable_header = PRState.processtable_header
    state.processtable_line = PRState.processtable_line

    def killed(_src, _dst):
        raise OSError("killed")

    monkeypatch.setattr("eon.fileio.os.replace", killed)
    with pytest.raises(OSError, match="killed"):
        state.save_process_table()
    assert path.read_text() == "OLD\n"
