"""A product that cannot be stored is not a transition."""

from io import StringIO

from eon.escaperate import register_results as escape_register
from eon.parallelreplica import register_results as replica_register


class _State:
    def __init__(self):
        self.number = 0
        self.time = 5.0
        self.procs = {}

    def get_time(self):
        return self.time

    def inc_time(self, dt):
        self.time += dt

    def zero_time(self):
        self.time = 0.0

    def load_process_table(self, force=False):
        return None

    def add_process(self, result):
        return None


class _States:
    def __init__(self, state):
        self.state = state

    def get_state(self, number):
        return self.state


class _Comm:
    def __init__(self, blobs):
        self.blobs = blobs
        self.cancelled = []

    def get_results(self, path, keep):
        yield from self.blobs

    def cancel_state(self, number):
        self.cancelled.append(number)
        return 1


def _blob():
    text = "2.0 speedup\n1 transition_found\n3.0 transition_time_s\n"
    return {
        "name": "0_1",
        "results.dat": StringIO(text),
        "product.con": StringIO("product\n"),
    }


def test_rejected_replica_product_is_not_a_transition(tmp_path):
    state = _State()
    comm = _Comm([_blob()])
    config = type("C", (), {"path_jobs_in": str(tmp_path / "in")})()
    _count, transition, _total = replica_register(comm, state, _States(state), config)
    assert transition is None
    assert state.get_time() == 5.0
    assert comm.cancelled == []


def test_rejected_escape_product_is_not_a_transition(tmp_path, monkeypatch):
    monkeypatch.setattr("eon.escaperate.io.loadcon", lambda _path: "geom")
    state = _State()
    config = type("C", (), {"path_jobs_in": str(tmp_path / "in"), "path_states": str(tmp_path)})()
    _count, transition, _total = escape_register(
        _Comm([_blob()]), state, _States(state), config
    )
    assert transition is None
    assert state.get_time() == 5.0
