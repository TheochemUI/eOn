"""dynamics.txt reopens at the next integer step and names the rate column."""

from eon.fileio import Dynamics


def test_empty_dynamics_reopens_at_step_zero(tmp_path):
    path = tmp_path / "dynamics.txt"
    Dynamics(str(path))
    again = Dynamics(str(path))
    assert again.next_step == 0
    again.append(1, 2, 3, 0.1, 0.2, 0.4, 5.0, -1.0)
    resumed = Dynamics(str(path))
    assert resumed.next_step == 1
    row = resumed.get()[0]
    assert row["barrier"] == 0.4
    assert row["rate"] == 5.0
    assert row["prefactor"] == 5.0
