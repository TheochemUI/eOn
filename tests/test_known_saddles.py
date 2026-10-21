"""Known saddles: the server list a client climb stops against."""

from io import StringIO

import numpy as np

from eon import fileio as io
from tests.test_library_bodies import _atoms, _result, _states


def _bad(reactant, saddle, termination, wuid, energy_saddle=-0.4):
    result = _result(reactant, saddle, _atoms(3.3), energy_saddle=energy_saddle)
    result["wuid"] = wuid
    result["results"]["termination_reason"] = termination
    return result


def test_save_known_saddles_writes_free_atoms_only():
    reactant = _atoms(2.5, names=("Pt", "Pt", "Pt"))
    reactant.free[0] = 0
    saddle = reactant.copy()
    saddle.r[2] = [0.0, 1.5, 0.25]
    buf = StringIO()
    io.save_known_saddles(buf, [(-0.7, saddle)], reactant)
    lines = buf.getvalue().splitlines()
    assert lines[0] == "1 2"
    assert float(lines[1]) == -0.7
    assert [list(map(float, line.split())) for line in lines[2:]] == [
        [2.5, 0.0, 0.0],
        [0.0, 1.5, 0.25],
    ]


def test_failed_saddles_are_stored_once_and_listed(tmp_path):
    _, _, state, _, proc, reactant = _states(tmp_path)
    saddle = _atoms(2.9)
    state.register_bad_saddle(_bad(reactant, saddle, 6, wuid=5))
    state.register_bad_saddle(_bad(reactant, saddle, 8, wuid=6))
    state.register_bad_saddle(_bad(reactant, _atoms(3.5), 5, wuid=7))
    index = (state.known_bad_path / "index").read_text().splitlines()
    assert index == ["5 -0.4000000000"]
    assert (state.known_bad_path / "saddle_5.con").is_file()

    state._known_bad = None
    known = state.known_saddles()
    assert [energy for energy, _ in known] == [state.procs[proc]["saddle_energy"], -0.4]
    assert np.allclose(known[1][1].r, saddle.r)


def test_known_saddle_result_is_a_repeat_or_a_failed_saddle(tmp_path):
    _, _, state, _, proc, reactant = _states(tmp_path)
    repeats = state.procs[proc]["repeats"]
    bad = state.get_bad_saddle_count()

    hit = _bad(reactant, _atoms(2.7), state.KNOWN_SADDLE, wuid=8, energy_saddle=-0.7)
    assert state.register_known_saddle(hit) is None
    state.load_process_table()
    assert state.procs[proc]["repeats"] == repeats + 1
    assert state.get_bad_saddle_count() == bad

    miss = _bad(reactant, _atoms(2.9), state.KNOWN_SADDLE, wuid=9)
    assert state.register_known_saddle(miss) is None
    assert state.get_bad_saddle_count() == bad + 1
    assert not state.known_bad_path.exists()
