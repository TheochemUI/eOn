"""``Displace.get_displacement`` on a Cu cell, against its per-atom definition.

A saddle search batch draws one displacement per search from the same
reactant. The draws must stay bit for bit what they were for a given
random state, and the neighbour geometry of the reactant must be built
once for the batch, not once per draw.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

pytest.importorskip("linkcell")
pytest.importorskip("minimage")

import akmc_synthetic as synth  # noqa: E402
from driver_oracles import count_events, ref_get_displacement  # noqa: E402

from eon import displace  # noqa: E402

RADIUS = 5.0
STD_DEV = 0.1


@pytest.fixture()
def disp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = synth.setup_run(tmp_path, 4)
    reactant = synth.build_cell(4)
    reactant.free[::5, :] = 0.0
    return displace.Random(
        reactant, std_dev=STD_DEV, radius=RADIUS, hole_epicenters=None, config=config
    )


def draw(fn, seed, *args):
    np.random.seed(seed)
    return fn(*args)


def assert_same_draw(got, want):
    (moved_g, direction_g), (moved_w, direction_w) = got, want
    assert np.array_equal(moved_g.r, moved_w.r)
    assert direction_g.dtype == direction_w.dtype
    assert np.array_equal(direction_g, direction_w)


@pytest.mark.parametrize("seed", [0, 1, 17])
def test_a_single_epicenter_draws_what_the_definition_draws(disp, seed):
    for epicenter in (7, 100, 255):
        assert_same_draw(
            draw(disp.get_displacement, seed, epicenter),
            draw(ref_get_displacement, seed, disp, epicenter),
        )


@pytest.mark.parametrize("kind", [list, np.array])
def test_listed_epicenters_draw_what_the_definition_draws(disp, kind):
    epicenters = kind([3, 4, 5, 6, 90, 91])
    assert_same_draw(
        draw(disp.get_displacement, 4, epicenters),
        draw(ref_get_displacement, 4, disp, epicenters),
    )


def test_draws_leave_the_random_state_where_the_definition_leaves_it(disp):
    draw(disp.get_displacement, 9, 12)
    after = np.random.random_sample()
    draw(ref_get_displacement, 9, disp, 12)
    assert np.random.random_sample() == after


def _pair_searches(disp, draws):
    disp = displace.Random(
        disp.reactant, std_dev=STD_DEV, radius=RADIUS, hole_epicenters=None,
        config=disp.config,
    )
    with count_events() as counts:
        np.random.seed(2)
        for epicenter in range(draws):
            disp.get_displacement(epicenter)
    return sum(v for k, v in counts.items() if k.endswith(".pairs_within"))


def test_the_neighbour_cells_are_built_once_for_a_batch_of_draws(disp):
    one = _pair_searches(disp, 1)
    many = _pair_searches(disp, 6)
    assert many <= one


def test_a_reactant_that_moves_between_draws_is_not_served_stale(disp):
    draw(disp.get_displacement, 0, 20)
    disp.reactant.r[20] += np.array([0.9, -0.4, 0.3])
    disp.reactant.r %= np.diag(disp.reactant.box)
    disp.neighbors_list = None
    assert_same_draw(
        draw(disp.get_displacement, 5, 20),
        draw(ref_get_displacement, 5, disp, 20),
    )
