"""A transitioning search is part of the parallel-replica speedup mean."""

from io import StringIO

import pytest

from eon.parallelreplica import register_results


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

    def add_process(self, result):
        return 4


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
        for blob in self.blobs:
            yield blob

    def cancel_state(self, number):
        return 0


def _result(name, speedup, found, **extra):
    lines = [f"{speedup} speedup", f"{found} transition_found"]
    for key, value in extra.items():
        lines.append(f"{value} {key}")
    return {"name": name, "results.dat": StringIO("\n".join(lines) + "\n")}


def test_transition_is_included_in_the_speedup_count(tmp_path):
    state = _State()
    config = type("C", (), {"path_jobs_in": str(tmp_path / "jobs_in")})()
    comm = _Comm(
        [
            _result("0_1", 2.0, 0, simulation_time_s=4.0),
            _result("0_2", 4.0, 1, transition_time_s=3.0, simulation_time_s=1.0),
        ]
    )
    count, transition, total = register_results(comm, state, _States(state), config)
    assert count == 2
    assert transition["process_id"] == 4
    assert total == pytest.approx(6.0)
    assert total / count == pytest.approx(3.0)
