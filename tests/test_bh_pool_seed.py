"""The basin-hopping pool draw follows the numpy seed."""

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
                "random_seed = 1",
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
    cfg.bh_initial_state_pool_size = 4
    return cfg


def _minimum(x):
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    atoms.r[0] = [float(x), 2.0, 3.0]
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    return buf


def test_repeated_numpy_seed_repeats_the_draw(tmp_path):
    import random

    states = BHStates(_config(tmp_path))
    for x, energy in ((1.0, -1.0), (2.0, -2.0), (3.0, -3.0), (4.0, -4.0)):
        states.add_state({"min.con": _minimum(x)}, {"minimum_energy": energy})

    def draw(numpy_seed, python_seed):
        np.random.seed(numpy_seed)
        random.seed(python_seed)
        return float(io.loadcon(states.get_random_minimum()).r[0, 0])

    assert draw(1, 1) == draw(1, 99)
