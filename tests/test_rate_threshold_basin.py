"""A fast crossing merges two states into a superbasin."""

from types import SimpleNamespace

from eon.superbasinscheme import RateThreshold


class _State:
    def __init__(self, number, barrier, prefactor):
        self.number = number
        self.procs = {
            1: {
                "product": 1 if number == 0 else 0,
                "barrier": barrier,
                "prefactor": prefactor,
            }
        }

    def get_process_table(self):
        return self.procs


class _States:
    def __init__(self):
        self.connected = None

    def get_num_states(self):
        return 0

    def connect_states(self, states):
        self.connected = list(states)

    def connect_state_sets(self, start_states, end_states):
        self.sets = (set(start_states), set(end_states))


def _scheme(tmp_path, threshold):
    cfg = SimpleNamespace(
        comp_use_identical=False,
        sb_max_size=0,
        amsel_discover_decide=False,
        sb_amsel_discover_decide=False,
    )
    return RateThreshold(tmp_path / "basins", _States(), 0.025, threshold, config=cfg)


def test_slow_crossing_does_not_merge(tmp_path):
    scheme = _scheme(tmp_path, 1.0e20)
    start = _State(0, 0.5, 1.0e12)
    end = _State(1, 0.5, 1.0e12)
    scheme.register_transition(start, end)
    assert scheme.superbasins == []


def test_fast_crossing_creates_a_basin(tmp_path):
    scheme = _scheme(tmp_path, 1.0)
    start = _State(0, 0.01, 1.0e12)
    end = _State(1, 0.01, 1.0e12)
    scheme.register_transition(start, end)
    assert len(scheme.superbasins) == 1
    assert set(scheme.superbasins[0].state_numbers) == {0, 1}
    text = (tmp_path / "basins" / "1").read_text().split()
    assert text == ["0", "1"] or set(text) == {"0", "1"}
