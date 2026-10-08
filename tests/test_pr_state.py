"""A parallel-replica state records a product and its transition time."""

from io import StringIO
from types import SimpleNamespace

import numpy as np

from eon import fileio as io
from eon.prstate import PRState
from eon.structure import Structure


def test_add_process_returns_an_id_and_reloads(tmp_path):
    atoms = Structure(1)
    atoms.names = ["Cu"]
    atoms.mass[:] = 63.5
    atoms.box = np.diag([10.0, 10.0, 10.0])
    atoms.r[0] = [1.0, 0.0, 0.0]
    product = atoms.copy()
    product.r[0, 0] = 1.4
    reactant_con = tmp_path / "reactant.con"
    io.savecon(str(reactant_con), atoms)
    state_dir = tmp_path / "0"
    state_dir.mkdir()
    state = PRState(
        str(state_dir),
        0,
        SimpleNamespace(kT=0.025),
        reactant_path=str(reactant_con),
        config=SimpleNamespace(),
    )
    reactant_buf = StringIO()
    io.savecon(reactant_buf, atoms)
    reactant_buf.seek(0)
    product_buf = StringIO()
    io.savecon(product_buf, product)
    product_buf.seek(0)
    results = StringIO("0 termination_reason\n")
    proc_id = state.add_process(
        {
            "results": {
                "potential_energy_reactant": -2.0,
                "potential_energy_product": -1.7,
                "transition_time_s": 3.5,
            },
            "reactant.con": reactant_buf,
            "product.con": product_buf,
            "results.dat": results,
        }
    )
    assert proc_id is not None
    assert state.get_energy() == -2.0
    state.procs = None
    state.load_process_table(force=True)
    assert state.procs[proc_id]["time"] == 3.5
    assert state.procs[proc_id]["product"] == -1
    state.inc_time(1.5)
    assert state.get_time() == 1.5
    state.zero_time()
    assert float(state.get_time()) == 0.0
