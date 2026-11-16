"""A new escape-rate product is rated from its own transition time."""

from types import SimpleNamespace

import pytest

from eon import escaperate


class _State:
    number = 0
    procs = {}

    def load_process_table(self, force=False):
        return None

    def add_process(self, result):
        return 7


def test_new_product_keeps_its_own_time(tmp_path):
    cfg = SimpleNamespace(path_states=str(tmp_path))
    table = escaperate.end_state_table_path(cfg, 0)
    escaperate.save_end_state_table(
        table,
        [{"state": 0, "views": 1, "rate": 0.1, "time": 10.0, "process_id": 3}],
    )
    got = escaperate.accumulate_end_state(
        "geom",
        _State(),
        {"results": {"transition_time_s": 2.0}},
        cfg,
    )
    assert got == 7
    rows = escaperate.load_end_state_table(table)
    assert rows[0]["time"] == pytest.approx(10.0)
    new = rows[1]
    assert new["process_id"] == 7
    assert new["views"] == 1
    assert new["time"] == pytest.approx(2.0)
    assert new["rate"] == pytest.approx(0.5)
