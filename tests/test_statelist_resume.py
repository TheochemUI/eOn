"""A pass that stops inside get_product_state leaves states/<n> behind.

The next new state must not take that directory and its reactant.
"""

import numpy as np
import pytest

from eon import fileio as io
from eon.akmcstatelist import AKMCStateList
from eon.prstatelist import PRStateList
from tests.test_library_bodies import _atoms, _buf, _config, _dat, _result


class _Stop(Exception):
    pass


def _list(tmp_path):
    con = tmp_path / "reactant.con"
    if not con.is_file():
        io.savecon(str(con), _atoms(2.5))
    cfg = _config(tmp_path)
    return AKMCStateList(300.0 / 11604.5, 20.0, 40.0, initial_state=str(con), config=cfg)


def _two_processes(states):
    reactant = _atoms(2.5)
    state = states.get_state(0)
    first = state.add_process(_result(reactant, _atoms(2.7), _atoms(3.1), -0.7))
    second = state.add_process(_result(reactant, _atoms(2.9), _atoms(3.6), -0.6))
    assert first is not None and second is not None and first != second
    return first, second


def _product_of(states, process_id):
    return io.loadcon(states.get_state(0).proc_product_path(process_id))


def _pr_list(tmp_path):
    con = tmp_path / "pos.con"
    if not con.is_file():
        io.savecon(str(con), _atoms(2.5))
    return PRStateList(str(con), config=_config(tmp_path))


def _pr_process(state, x):
    fields = {
        "potential_energy_reactant": -1.0,
        "potential_energy_product": -0.8,
        "transition_time_s": 1.5,
    }
    return state.add_process(
        {
            "reactant.con": _buf(_atoms(2.5)),
            "product.con": _buf(_atoms(x)),
            "results.dat": _dat(**fields),
            "results": fields,
        }
    )


def _stop_after_the_zero_row(states, monkeypatch):
    append = states.append_state_table

    def stop(energy):
        if states.get_num_states() > 0:
            raise _Stop
        append(energy)

    monkeypatch.setattr(states, "append_state_table", stop)


def test_a_stop_before_the_process_points_at_the_state(tmp_path, monkeypatch):
    states = _list(tmp_path)
    first, second = _two_processes(states)

    def stop(*args):
        raise _Stop

    monkeypatch.setattr(states, "register_process", stop)
    with pytest.raises(_Stop):
        states.get_product_state(0, first)
    assert (tmp_path / "states" / "1").is_dir()

    resumed = _list(tmp_path)
    product = resumed.get_product_state(0, second)
    assert product.number == 1
    assert np.allclose(product.get_reactant().r, _product_of(resumed, second).r)


def test_a_stop_before_the_state_table_row(tmp_path, monkeypatch):
    states = _list(tmp_path)
    first, second = _two_processes(states)
    _stop_after_the_zero_row(states, monkeypatch)
    with pytest.raises(_Stop):
        states.get_product_state(0, first)

    resumed = _list(tmp_path)
    product = resumed.get_product_state(0, second)
    assert product.number == 2
    assert np.allclose(product.get_reactant().r, _product_of(resumed, second).r)
    assert resumed.get_num_states() == 3

    followed = resumed.get_product_state(0, first)
    assert followed.number == 1
    assert np.allclose(followed.get_reactant().r, _product_of(resumed, first).r)


def test_a_parallel_replica_stop_before_the_state_table_row(tmp_path, monkeypatch):
    states = _pr_list(tmp_path)
    first = _pr_process(states.get_state(0), 4.2)
    second = _pr_process(states.get_state(0), 7.5)
    _stop_after_the_zero_row(states, monkeypatch)
    with pytest.raises(_Stop):
        states.get_product_state(0, first)

    resumed = _pr_list(tmp_path)
    product = resumed.get_product_state(0, second)
    assert product.number == 2
    assert np.allclose(product.get_reactant().r, _atoms(7.5).r)
    assert resumed.get_num_states() == 3

    followed = resumed.get_product_state(0, first)
    assert followed.number == 1
    assert np.allclose(followed.get_reactant().r, _atoms(4.2).r)


def test_a_row_on_another_state_keeps_the_directory(tmp_path, monkeypatch):
    states = _list(tmp_path)
    first, second = _two_processes(states)
    one = states.get_state(states.get_product_state(0, first).number)
    third = one.add_process(_result(_atoms(3.1), _atoms(3.3), _atoms(5.0), -0.5))
    assert third is not None
    _stop_after_the_zero_row(states, monkeypatch)
    with pytest.raises(_Stop):
        states.get_product_state(1, third)

    # A directory reused under the old numbering still names state 0.
    info = io.ini(str(tmp_path / "states" / "2" / "info"))
    info.read()
    info.set("MetaData", "previous state", 0)

    resumed = _list(tmp_path)
    product = resumed.get_product_state(0, second)
    assert product.number == 3
    assert np.allclose(product.get_reactant().r, _product_of(resumed, second).r)

    followed = resumed.get_product_state(1, third)
    assert followed.number == 2
    assert np.allclose(followed.get_reactant().r, _atoms(5.0).r)
