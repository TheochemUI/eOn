"""A bad token in the query list does not drop the states already recorded."""

from types import SimpleNamespace

from eon.process_catalog import _mark_queried, was_queried


def _config(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    return SimpleNamespace(kdb_scratch_path=scratch)


def test_bad_token_keeps_the_recorded_state(tmp_path):
    config = _config(tmp_path)
    path = tmp_path / "scratch" / "queried"
    path.write_text("0\njunk\n")
    state = SimpleNamespace(number=0)
    assert was_queried(state, config) is True
    _mark_queried(SimpleNamespace(number=5), config)
    recorded = path.read_text().split()
    assert recorded == ["0", "5"]
