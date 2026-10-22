"""A torn search-log line is dropped, and the next append cuts it off."""

import logging

import pytest
from types import SimpleNamespace

from eon.akmcstate import AKMCState


class _Info:
    def __init__(self):
        self.store = {}

    def get(self, section, option, default=None):
        return self.store.get((section, option), default)

    def set(self, section, option, value):
        self.store[(section, option)] = value


def _result(wuid):
    return {
        "wuid": wuid,
        "type": "random",
        "results": {
            "barrier_reactant_to_product": 0.2,
            "displacement_saddle_distance": 1.0,
            "force_calls_saddle": 3,
            "force_calls_minimization": 4,
            "force_calls_prefactors": 5,
        },
    }


def _state(path):
    state = AKMCState.__new__(AKMCState)
    state.number = 0
    state.search_result_path = str(path)
    state.config = SimpleNamespace(akmc_server_side_process_search=False)
    state.info = _Info()
    return state


def _line(wuid, comment):
    return "%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (
        wuid, "random", 0.2, 1.0, 3, 4, 5, comment,
    )


def _torn_tail(kind):
    if kind == "newline":
        return _line(2, "TORNTAIL").rstrip("\n")
    return "2 TORNTAIL\n"


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_search_line_is_replaced(tmp_path, caplog, kind):
    path = tmp_path / "search dir" / "searches.txt"
    path.parent.mkdir()
    good = "header\n----\n       1     random    0.10000    1.00000        10        20         0    good\n"
    path.write_text(good + _torn_tail(kind))
    state = _state(path)
    caplog.set_level(logging.WARNING, logger="io")
    assert state.get_number_of_searches() == 1
    assert "torn final line" in caplog.text
    state.info = _Info()
    state.append_search_result(_result(3), "kept", None)
    text = path.read_text()
    assert "TORNTAIL" not in text
    assert text.startswith("header\n----\n")
    assert text.endswith(_line(3, "kept"))
    state.info = _Info()
    assert state.get_number_of_searches() == 2
