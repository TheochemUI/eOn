"""In-process dynamics, a mismatched table, and a process standard output."""

import math
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.askmc import ASKMC
from eon.communicator_inprocess import _require_pyeonclient, _run_inprocess_job
from eon.structure import Structure, as_atom_free, coerce_free
from tests.test_library_bodies import _atoms, _config, _states


def _matter():
    pc = _require_pyeonclient()
    cluster = _atoms(2.5)
    params = pc.Parameters()
    params.potential = pc.PotType.LJ
    params.quiet = True
    params.write_log = False
    pot = pc.make_potential(params)
    from pyeonclient.bridge import structure_to_matter

    return pc, pot, params, structure_to_matter(cluster, pot, params)


def test_dynamics_monte_carlo_hopping_and_a_saddle(tmp_path):
    pc, pot, params, matter = _matter()
    for kind, name in (
        (pc.JobType.Dynamics, "dynamics"),
        (pc.JobType.Monte_Carlo, "monte_carlo"),
        (pc.JobType.Basin_Hopping, "basin_hopping"),
    ):
        payload = _run_inprocess_job(pc, kind, matter, pot, params, {})
        assert payload["job_type"] == name
        assert math.isfinite(payload["energy"])
    mode = np.zeros((2, 3))
    mode[1, 0] = 1.0
    direction = StringIO()
    io.save_mode(direction, mode)
    saddle = _run_inprocess_job(
        pc,
        pc.JobType.Saddle_Search,
        matter,
        pot,
        params,
        {"direction.dat": direction},
    )
    assert saddle["job_type"] == "saddle_search"
    assert math.isfinite(saddle["energy"])


def test_structure_edges_and_a_mismatched_table(tmp_path):
    assert coerce_free([], 0).shape == (0, 3)
    row = coerce_free([1.0], 1)
    assert row.shape == (1, 3)
    assert as_atom_free(np.array([0.0, 1.0])).tolist() == [False, True]
    with pytest.raises(ValueError, match="free must"):
        coerce_free([1.0, 0.0], 4)
    cluster = _atoms(2.5)
    cluster.atom_ids = np.array([1], dtype=np.uint64)
    assert cluster.ids_or_sequential().tolist() == [1, 2]
    with pytest.raises(ValueError, match="length-3"):
        cluster.append([1.0, 0.0, 0.0], [1.0, 0.0], "Pt", 195.0)

    path = tmp_path / "two.poscar"
    io.saveposcar(str(path), _atoms(2.5))
    text = path.read_text()
    path.write_text(text + "\n" + text)
    frames = io.loadposcars(str(path))
    assert len(frames) >= 1
    assert frames[0].names[0] == "Pt"

    table_path = tmp_path / "odd.tbl"
    table_path.write_text("alpha beta\n--------\n1 2\n")
    table = io.Table(str(table_path), ["alpha", "gamma"])
    with pytest.raises(io.TableException, match="mismatch"):
        len(table)


def test_stdout_is_stored_and_a_modified_row_must_exist(tmp_path):
    cfg, states, state, product, proc, _reactant = _states(tmp_path)
    from eon.state import State

    State.add_process(state, {"stdout.dat": StringIO("client log\n")})
    stored = list(Path(state.procdata_path).glob("stdout_*.dat"))
    assert stored
    assert "client log" in stored[0].read_text()

    ask = ASKMC(
        states.kT,
        states,
        0.5,
        1.5,
        2.0,
        False,
        False,
        False,
        str(tmp_path),
        20.0,
        str(tmp_path / "recycle"),
    )
    forward = dict(state.get_process(proc))
    forward["view_count"] = 1
    ask.save_modified_process_table(state, {proc: forward})
    compiled = ask.compile_process_table(state)
    assert compiled[proc]["view_count"] == 1
