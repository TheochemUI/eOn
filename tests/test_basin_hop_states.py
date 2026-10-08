"""Basin hopping records a new minimum and can draw it back."""

from io import StringIO

import numpy as np

from eon import fileio as io
from eon.basinhopping import BHStates
from eon.config import ConfigClass
from eon.structure import Structure


def _config(directory):
    path = directory / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = basin_hopping",
                "temperature = 300",
                "[Paths]",
                f"main_directory = {directory}",
                f"results = {directory}",
                f"states = {directory / 'states'}",
                "",
            ]
        )
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    cfg.bh_initial_state_pool_size = 5
    return cfg


def _min():
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    atoms.r[0] = [1.0, 2.0, 3.0]
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    return atoms, buf


def test_new_minimum_can_be_drawn(tmp_path):
    atoms, buf = _min()
    states = BHStates(_config(tmp_path))
    assert states.get_random_minimum() is None
    states.add_state({"min.con": buf}, {"minimum_energy": -1.25})
    drawn = states.get_random_minimum()
    assert drawn is not None
    loaded = io.loadcon(drawn)
    assert loaded.names == ["Pt"]
    assert np.allclose(loaded.r, atoms.r)
    again = StringIO()
    io.savecon(again, atoms)
    again.seek(0)
    states.add_state({"min.con": again}, {"minimum_energy": -1.25})
    assert states.energy_table.rows[0]["repeats"] == 1
