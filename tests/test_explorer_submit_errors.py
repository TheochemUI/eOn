"""A submit that returns in the same clock tick is not a failure."""

import logging
from types import SimpleNamespace

import pytest

from eon.explorer import ClientMinModeExplorer, ServerMinModeExplorer
from eon import fileio as io


def _client(submit):
    explorer = ClientMinModeExplorer.__new__(ClientMinModeExplorer)
    explorer.config = SimpleNamespace(
        comm_job_buffer_size=1,
        comm_type="inprocess",
        process_search_known_saddles=False,
        saddle_method="min_mode",
        recycling_on=False,
        config_path="config.ini",
    )
    explorer.comm = SimpleNamespace(queued_search_count=lambda: 0, submit_jobs=submit)
    explorer.reactant = object()
    explorer.state = SimpleNamespace(
        number=0,
        info=SimpleNamespace(get=lambda *args, **kwargs: ""),
    )
    explorer.job_table = SimpleNamespace(add_row=lambda row: None, write=lambda: None)
    explorer.wuid = 0
    explorer.save_wuid = lambda: None
    explorer.generate_displacement = lambda: (None, None, "dynamics")
    explorer.nrecycled = 0
    return explorer


def _server(submit):
    explorer = ServerMinModeExplorer.__new__(ServerMinModeExplorer)
    explorer.config = SimpleNamespace(
        comm_job_buffer_size=1,
        comm_job_max_size=0,
        recycling_on=False,
    )
    explorer.comm = SimpleNamespace(
        queued_search_count=lambda: 0,
        in_progress_search_count=lambda: 0,
        submit_jobs=submit,
    )
    explorer.state = SimpleNamespace(number=0)
    explorer.process_searches = {
        0: SimpleNamespace(
            search_id=0,
            get_job=lambda state_number: ({"blob": 1}, "saddle_search"),
        )
    }
    explorer.wuid_to_search_id = {}
    explorer.job_info = {}
    explorer.wuid = 0
    explorer.search_id = 0
    explorer.save_wuid = lambda: None
    explorer.save = lambda: None
    explorer.nrecycled = 0
    return explorer


def test_zero_duration_submit_is_not_a_failure(monkeypatch, caplog):
    monkeypatch.setattr("eon.explorer.time", lambda: 1000.0)
    monkeypatch.setattr(io, "savecon", lambda *args, **kwargs: None)
    monkeypatch.setattr(io, "modify_config", lambda *args, **kwargs: io.StringIO("x"))
    makers = (
        _client(lambda searches, invariants: None),
        _server(lambda jobs, invariants: None),
    )
    with caplog.at_level(logging.ERROR):
        for explorer in makers:
            explorer.make_jobs()
    assert "Failed to submit" not in caplog.text


def test_submit_interrupt_propagates(monkeypatch):
    monkeypatch.setattr("eon.explorer.time", lambda: 1000.0)
    monkeypatch.setattr(io, "savecon", lambda *args, **kwargs: None)
    monkeypatch.setattr(io, "modify_config", lambda *args, **kwargs: io.StringIO("x"))

    def interrupt(jobs, invariants):
        raise KeyboardInterrupt

    for explorer in (_client(interrupt), _server(interrupt)):
        with pytest.raises(KeyboardInterrupt):
            explorer.make_jobs()
