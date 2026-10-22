"""A torn process-table row is dropped, and the next append cuts it off."""

import logging

import pytest
from types import SimpleNamespace

from eon.akmcstate import AKMCState


class _Info:
    def get(self, section, option, default=None):
        if option == "kT":
            return 0.025
        return default

    def set(self, *args, **kwargs):
        return None


def _state(path):
    state = AKMCState.__new__(AKMCState)
    state.number = 0
    state.proctable_path = str(path)
    state.procs = None
    state.info = _Info()
    state.statelist = SimpleNamespace(kT=0.025)
    state.config = SimpleNamespace(akmc_eq_rate=0.0)
    path.write_text(
        AKMCState.processtable_header
        + AKMCState.processtable_line % (1, 0.2, 1.0e12, 2, 0.0, 1.0e12, 0.2, 1.0, 0)
    )
    return state


def _torn_tail(kind):
    whole = AKMCState.processtable_line % (9, 1.23456, 1.0e12, 4, 0.0, 1.0e12, 0.3, 0.5, 0)
    if kind == "newline":
        return whole.rstrip("\n")
    return "9 1.23456\n"


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_process_row_is_replaced(tmp_path, caplog, kind):
    table = tmp_path / "processtable"
    state = _state(table)
    table.write_bytes(table.read_bytes() + _torn_tail(kind).encode())
    caplog.set_level(logging.WARNING, logger="io")
    state.load_process_table(force=True)
    assert 1 in state.procs
    assert 9 not in state.procs
    assert "torn final line" in caplog.text
    state.append_process_table(7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0)
    text = table.read_text()
    assert "1.23456" not in text
    assert text.startswith(AKMCState.processtable_header)
    assert text.endswith(
        AKMCState.processtable_line % (7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0)
    )
    state.procs = None
    state.load_process_table(force=True)
    assert 1 in state.procs
    assert 7 in state.procs
    assert 9 not in state.procs
