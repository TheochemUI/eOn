"""A basin-hopping state directory may already exist."""

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
    return cfg


def test_existing_state_directory_receives_the_minimum(tmp_path):
    cfg = _config(tmp_path)
    states = BHStates(cfg)
    state_dir = tmp_path / "states" / "0"
    state_dir.mkdir()
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    atoms.r[0] = [1.0, 2.0, 3.0]
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    assert states.add_state({"min.con": buf}, {"minimum_energy": -1.25}) is True
    loaded = io.loadcon(str(state_dir / "minimum.con"))
    assert np.allclose(loaded.r, atoms.r)
