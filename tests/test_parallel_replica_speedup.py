"""Parallel replica reports the mean speedup, and the register returns the sum."""

from io import StringIO

import pytest

from eon.escaperate import parallelreplica, register_results


class _State:
    def __init__(self):
        self.number = 0
        self.time = 0.0

    def get_time(self):
        return self.time

    def inc_time(self, dt):
        self.time += dt

    def zero_time(self):
        self.time = 0.0


class _States:
    def __init__(self, state):
        self.state = state

    def get_state(self, number):
        assert number == 0
        return self.state


class _Comm:
    def __init__(self, blobs):
        self.blobs = blobs

    def get_results(self, path, keep):
        assert keep("0_1") is True
        for blob in self.blobs:
            yield blob


def _result(name, speedup, sim_time):
    text = f"{speedup} speedup\n0 transition_found\n{sim_time} simulation_time_s\n"
    return {"name": name, "results.dat": StringIO(text)}


def test_register_results_returns_the_speedup_sum(tmp_path):
    state = _State()
    config = type("C", (), {"path_jobs_in": str(tmp_path / "jobs_in")})()
    comm = _Comm([_result("0_1", 2.0, 4.0), _result("0_2", 4.0, 1.0)])
    count, transition, total = register_results(comm, state, _States(state), config)
    assert count == 2
    assert transition is None
    assert total == pytest.approx(6.0)
    assert state.get_time() == pytest.approx(5.0)


def test_replica_logs_the_average_speedup(tmp_path, monkeypatch, caplog):
    root = tmp_path / "root"
    root.mkdir()
    results = tmp_path / "results"
    results.mkdir()
    state = _State()
    monkeypatch.setattr("eon.escaperate.get_pr_metadata", lambda config: (0, 1.0, 3))
    monkeypatch.setattr("eon.escaperate.get_statelist", lambda config: _States(state))
    monkeypatch.setattr("eon.escaperate.communicator.get_communicator", lambda config: object())
    monkeypatch.setattr(
        "eon.escaperate.register_results",
        lambda *args, **kwargs: (2, None, 6.0),
    )
    monkeypatch.setattr("eon.escaperate.make_searches", lambda *args, **kwargs: 9)
    monkeypatch.setattr("eon.escaperate.io.write_info_txt", lambda *args, **kwargs: None)
    monkeypatch.setattr("eon.escaperate.io.save_prng_state", lambda *args, **kwargs: None)
    config = type(
        "C",
        (),
        {"path_root": str(root), "path_results": str(results)},
    )()
    caplog.set_level("INFO", logger="pr")
    parallelreplica(config)
    assert any(message == "Average speedup is 3.000000" for message in caplog.messages)
