"""The first sighting of a process counts toward Nf."""

from eon.askmc import ASKMC


class _State:
    def __init__(self, number, path, procs):
        self.number = number
        self.path = str(path)
        self.procs = procs

    def load_process_table(self):
        return None

    def get_energy(self):
        return 0.0


def test_first_visit_meets_nf(tmp_path, monkeypatch):
    proc = {
        "saddle_energy": 0.2,
        "prefactor": 1.0e12,
        "product": 1,
        "product_energy": 0.0,
        "product_prefactor": 1.0e12,
        "barrier": 0.2,
        "rate": 1.0,
        "repeats": 0,
    }
    origin_dir = tmp_path / "state 0"
    next_dir = tmp_path / "state 1"
    origin_dir.mkdir()
    next_dir.mkdir()
    origin = _State(0, origin_dir, {3: proc})
    nxt = _State(1, next_dir, {})
    states = {0: origin, 1: nxt}
    ask = ASKMC(
        kT=0.025,
        states=type("S", (), {"get_state": staticmethod(lambda n: states[n])})(),
        confidence=0.5,
        alpha=1.5,
        gamma=2.0,
        barrier_test_on=False,
        connection_test_on=False,
        sb_recycling_on=False,
        path_root=str(tmp_path),
        thermal_window=20.0,
        recycle_path=None,
    )
    assert ask.Nf == 1
    seen = []
    monkeypatch.setattr(ask, "raiseup", lambda *args, **kwargs: seen.append(args))
    ask.register_transition(origin, nxt)
    assert seen
