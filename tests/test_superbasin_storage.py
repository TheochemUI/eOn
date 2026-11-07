"""Merging basins archives the old basin when storage was not created."""

from types import SimpleNamespace

from eon.superbasinscheme import RateThreshold


class _State:
    def __init__(self, number):
        self.number = number
        self.path = ""


class _States:
    def connect_states(self, group):
        return None

    def get_state(self, number):
        return None


def test_absorbing_a_basin_stores_it(tmp_path):
    root = tmp_path / "sb dir"
    root.mkdir()
    states = [_State(i) for i in range(3)]
    cfg = SimpleNamespace(sb_max_size=0, comp_use_identical=True)
    scheme = RateThreshold(str(root), _States(), 0.025, 1.0, config=cfg)
    scheme.make_basin([states[0], states[1]])
    scheme.make_basin([states[0], states[2]])
    assert (root / "storage" / "1").is_file()
