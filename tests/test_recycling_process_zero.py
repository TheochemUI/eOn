"""A process id of zero can be the fastest exit from a superbasin."""

from types import SimpleNamespace

from eon.recycling import SB_Recycling


class _State:
    def __init__(self, number, procs=None):
        self.number = number
        self.procs = {} if procs is None else procs

    def load_process_table(self):
        return None


class _Basin:
    def __init__(self, states):
        self.states = states


class _Scheme:
    def __init__(self, basin):
        self.basin = basin

    def get_containing_superbasin(self, state):
        return self.basin


def test_process_id_zero_is_the_fastest_exit(tmp_path, monkeypatch):
    monkeypatch.setattr(
        SB_Recycling,
        "generate_corresponding_states",
        lambda self: None,
    )
    fast = _State(0, {0: {"product": 4, "rate": 10.0}})
    slow = _State(1, {5: {"product": 4, "rate": 1.0}})
    arrived = _State(4)
    recycler = SB_Recycling(
        SimpleNamespace(),
        slow,
        arrived,
        0.5,
        False,
        str(tmp_path / "recycle"),
        "mcacm",
        _Scheme(_Basin([slow, fast])),
    )
    chosen = {
        pair[0]: pair[1]
        for pair in recycler.sb_state_nums
    }
    assert chosen[fast.number] == arrived.number
    assert chosen[slow.number] is None
