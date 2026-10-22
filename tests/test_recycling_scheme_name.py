"""Superbasin recycling accepts the scheme name mcamc."""

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


def test_mcamc_records_the_exited_basin(tmp_path, monkeypatch):
    monkeypatch.setattr(
        SB_Recycling,
        "generate_corresponding_states",
        lambda self: None,
    )
    previous = _State(0, {1: {"product": 4, "rate": 3.0}})
    other = _State(2, {1: {"product": 9, "rate": 1.0}})
    current = _State(4)
    recycler = SB_Recycling(
        states=SimpleNamespace(),
        previous_state=previous,
        current_state=current,
        move_distance=0.5,
        recycle_save=False,
        path=str(tmp_path / "recycle"),
        sb_scheme="mcamc",
        superbasining=_Scheme(_Basin([previous, other])),
    )
    assert [4] == [pair[1] for pair in recycler.sb_state_nums if pair[1] == 4]
