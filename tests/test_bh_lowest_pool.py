"""The basin-hopping pool is the lowest-energy minima."""

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
    cfg.bh_initial_state_pool_size = 1
    return cfg


def _minimum(x):
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    atoms.r[0] = [x, 2.0, 3.0]
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    return buf


def test_pool_of_one_draws_the_lowest_energy(tmp_path):
    states = BHStates(_config(tmp_path))
    states.add_state({"min.con": _minimum(1.0)}, {"minimum_energy": -1.0})
    states.add_state({"min.con": _minimum(6.0)}, {"minimum_energy": -8.0})
    drawn = io.loadcon(states.get_random_minimum())
    assert drawn.r[0, 0] == 6.0
