"""A torn modified-process row is dropped, and the next append cuts it off."""

import logging

import pytest

from eon.askmc import ASKMC, mod_processtable_header, processtable_line


class _State:
    def __init__(self, path):
        self.path = str(path)
        self.number = 0


def _ask(tmp_path):
    return ASKMC(
        kT=0.025,
        states=None,
        confidence=0.9,
        alpha=1.5,
        gamma=2.0,
        barrier_test_on=False,
        connection_test_on=False,
        sb_recycling_on=False,
        path_root=str(tmp_path),
        thermal_window=20.0,
        recycle_path=None,
    )


def _torn_tail(kind):
    whole = processtable_line % (9, 1.23456, 1.0e12, 4, 0.0, 1.0e12, 0.3, 0.5, 0)
    if kind == "newline":
        return whole.rstrip("\n")
    return "9 1.23456\n"


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_modified_row_is_replaced(tmp_path, caplog, kind):
    state_dir = tmp_path / "state 0"
    state_dir.mkdir()
    ask = _ask(tmp_path)
    state = _State(state_dir)
    ask.append_modified_process_table(
        state, 1, 0.2, 1.0e12, 2, 0.0, 1.0e12, 0.2, 1.0, 1
    )
    table = state_dir / "askmc_processtable"
    table.write_bytes(table.read_bytes() + _torn_tail(kind).encode())
    caplog.set_level(logging.WARNING, logger="io")
    loaded = ask.get_modified_process_table(state)
    assert 1 in loaded
    assert 9 not in loaded
    assert "torn final line" in caplog.text
    ask.append_modified_process_table(
        state, 7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0
    )
    text = table.read_text()
    assert "1.23456" not in text
    assert text.startswith(mod_processtable_header)
    assert text.endswith(processtable_line % (7, 0.3, 1.0e12, 3, 0.0, 1.0e12, 0.3, 0.5, 0))
    loaded = ask.get_modified_process_table(state)
    assert 1 in loaded
    assert 7 in loaded
    assert 9 not in loaded
