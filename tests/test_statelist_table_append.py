"""A torn state-table row is dropped, and the next append cuts it off."""

import logging

import pytest
from types import SimpleNamespace

from eon.statelist import StateList


def _listing(root):
    cfg = SimpleNamespace(
        path_states=str(root),
        comp_eps_e=1.0e-3,
        comp_eps_r=0.1,
        comp_use_identical=False,
    )
    return StateList(object, config=cfg)


def _torn_tail(kind):
    whole = "% 7d %16.5f\n" % (9, 3.33333)
    if kind == "newline":
        return whole.rstrip("\n")
    return "3.33333\n"


@pytest.mark.parametrize("kind", ["newline", "fields"])
def test_torn_state_row_is_replaced(tmp_path, caplog, kind):
    root = tmp_path / "states dir"
    (root / "0").mkdir(parents=True)
    table = root / "state_table"
    good = "      0         -1.00000\n"
    table.write_text(good + _torn_tail(kind))
    listing = _listing(root)
    caplog.set_level(logging.WARNING, logger="io")
    assert listing.get_num_states() == 1
    assert listing.state_table == [-1.0]
    assert "torn final line" in caplog.text
    listing.append_state_table(-2.0)
    text = table.read_text()
    assert "3.33333" not in text
    assert text.startswith(good)
    assert text.endswith("% 7d %16.5f\n" % (1, -2.0))
    listing.state_table = None
    assert listing.get_num_states() == 2
    assert listing.state_table == [-1.0, -2.0]
