"""The numpy MCAMC solver must not alter the matrices it is given."""

import numpy as np

from eon.mcamc.mcamc import np_mcamc


def test_numpy_solver_keeps_the_rate_matrices():
    q_matrix = np.array([[0.0, 2.0], [3.0, 0.0]])
    exit_matrix = np.array([[4.0], [5.0]])
    costs = np.array([6.0, 8.0])
    q_kept = q_matrix.copy()
    exit_kept = exit_matrix.copy()
    costs_kept = costs.copy()

    times, branching, _residual = np_mcamc(q_matrix, exit_matrix, costs)
    assert np.allclose(q_matrix, q_kept)
    assert np.allclose(exit_matrix, exit_kept)
    assert np.allclose(costs, costs_kept)

    times_again, branching_again, _residual_again = np_mcamc(
        q_matrix, exit_matrix, costs
    )
    assert np.allclose(times_again, times)
    assert np.allclose(branching_again, branching)
