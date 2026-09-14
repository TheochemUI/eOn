"""Escape-rate end-state table: read from start, ConfigClass paths, xxh64 products."""
from __future__ import annotations

import os
from types import SimpleNamespace

import pytest


def test_end_state_table_path_uses_config_states():
    from eon.escaperate import end_state_table_path

    cfg = SimpleNamespace(path_states="/sim/states")
    assert end_state_table_path(cfg, 0) == os.path.join("/sim/states", "0", "end_state_table")
    assert end_state_table_path(cfg, 3) == os.path.join("/sim/states", "3", "end_state_table")


def test_load_end_state_table_reads_existing_rows(tmp_path):
    from eon.escaperate import load_end_state_table

    path = tmp_path / "end_state_table"
    path.write_text(
        "state       views         rate        time        process_id\n"
        "0             2             0.5             4.0             99\n"
        "1             1             0.25             4.0             100\n"
    )
    rows = load_end_state_table(str(path))
    assert len(rows) == 2
    assert rows[0]["state"] == 0
    assert rows[0]["views"] == 2
    assert rows[0]["process_id"] == 99
    assert rows[1]["process_id"] == 100


def test_load_end_state_table_empty_when_missing(tmp_path):
    from eon.escaperate import load_end_state_table

    assert load_end_state_table(str(tmp_path / "nope")) == []


def test_save_then_load_preserves_history(tmp_path):
    from eon.escaperate import load_end_state_table, save_end_state_table

    path = str(tmp_path / "end_state_table")
    rows = [
        {"state": 0, "views": 2, "rate": 0.5, "time": 4.0, "process_id": 99},
        {"state": 1, "views": 1, "rate": 0.25, "time": 4.0, "process_id": 100},
    ]
    save_end_state_table(path, rows)
    loaded = load_end_state_table(path)
    assert loaded[0]["views"] == 2
    assert loaded[1]["process_id"] == 100

    rows[0]["views"] = 3
    rows[0]["time"] = 6.0
    rows[0]["rate"] = 0.5
    save_end_state_table(path, rows)
    loaded = load_end_state_table(path)
    assert loaded[0]["views"] == 3
    assert loaded[0]["time"] == pytest.approx(6.0)
    assert len(loaded) == 2


def test_find_matching_process_uses_proc_product_path(tmp_path, monkeypatch):
    from eon import escaperate

    product_id = 0x123456789ABCDEF
    procdir = tmp_path / "procdata"
    procdir.mkdir()
    product_path = procdir / ("product_%d.con" % product_id)
    product_path.write_text("product-bytes")
    # Sequential names from the old loop must not be consulted.
    (procdir / "product_0.con").write_text("wrong")

    class FakeState:
        procs = {product_id: {"product": -1}}

        def load_process_table(self, force=False):
            return None

        def proc_product_path(self, pid):
            return str(procdir / ("product_%d.con" % pid))

    loaded = []

    def fake_loadcon(path):
        loaded.append(path)
        return path

    def fake_match(a, b, *args, **kwargs):
        return b == str(product_path)

    monkeypatch.setattr(escaperate.io, "loadcon", fake_loadcon)
    monkeypatch.setattr(escaperate.atoms, "match", fake_match)

    cfg = SimpleNamespace(
        comp_eps_r=0.1,
        comp_neighbor_cutoff=3.0,
        comp_check_rotation=False,
        comp_use_identical=False,
    )
    found = escaperate.find_matching_process("current", FakeState(), cfg)
    assert found == product_id
    assert loaded == [str(product_path)]


def test_accumulate_repeat_does_not_add_process(tmp_path, monkeypatch):
    from eon import escaperate

    pid = 42
    added = []

    class FakeState:
        number = 0
        procs = {pid: {}}

        def load_process_table(self, force=False):
            return None

        def proc_product_path(self, _pid):
            return str(tmp_path / "product.con")

        def add_process(self, result):
            added.append(result)
            return 99

    (tmp_path / "product.con").write_text("x")
    monkeypatch.setattr(escaperate.io, "loadcon", lambda path: "geom")
    monkeypatch.setattr(escaperate.atoms, "match", lambda *a, **k: True)

    cfg = SimpleNamespace(
        path_states=str(tmp_path),
        comp_eps_r=0.1,
        comp_neighbor_cutoff=3.0,
        comp_check_rotation=False,
        comp_use_identical=False,
    )
    table = escaperate.end_state_table_path(cfg, 0)
    escaperate.save_end_state_table(
        table,
        [{"state": 0, "views": 1, "rate": 1.0, "time": 1.0, "process_id": pid}],
    )
    result = {"results": {"transition_time_s": 1.0}}
    got = escaperate.accumulate_end_state("geom", FakeState(), result, cfg)
    assert got == pid
    assert added == []
    rows = escaperate.load_end_state_table(table)
    assert rows[0]["views"] == 2
    assert rows[0]["time"] == pytest.approx(2.0)
