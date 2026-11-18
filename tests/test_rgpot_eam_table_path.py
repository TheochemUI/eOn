"""[RgpotPot] table_path has to reach the in-process EAM potential.

A config that names an eam/alloy file and backend = eam must evaluate
that file. An empty path means the ini key was dropped before the
potential was built.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def _setfl(path: Path) -> None:
    # One element, five points. The pair term is nonzero inside the cutoff.
    lines = [
        "synthetic Cu setfl",
        "one element",
        "five points",
        "1 Cu",
        "5 1.0 5 0.5 2.0",
        "29 63.546 3.615 fcc",
        "0.0 0.2 0.4 0.6 0.8",
        "0.05 0.05 0.04 0.02 0.0",
        "0.0 0.4 0.2 0.1 0.0",
    ]
    path.write_text("\n".join(lines) + "\n")


def test_rgpot_eam_table_path_is_evaluated(tmp_path: Path) -> None:
    table = tmp_path / "cu.setfl"
    _setfl(table)
    missing = tmp_path / "missing.setfl"
    ini = f"""[Potential]
potential = rgpot

[RgpotPot]
backend = eam
table_path = {table}
table_format = setfl
neighbor_skin = 0.5
"""
    params = pyec.Parameters()
    params.load_ini_text(ini)
    pot = pyec.make_potential(params)
    matter = pyec.Matter(pot, params)
    matter.resize(2)
    matter.cell = np.eye(3, dtype=np.float64) * 20.0
    matter.positions = np.array(
        [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=np.float64
    )
    matter.masses = np.array([63.546, 63.546], dtype=np.float64)
    matter.atomic_numbers = np.array([29, 29], dtype=np.int64)
    matter.fixed = np.array([0, 0], dtype=np.int64)
    matter.periodic = True
    energy = float(matter.potential_energy)
    assert math.isfinite(energy)
    assert energy != 0.0

    blank = ini.replace(str(table), str(missing))
    refused = pyec.Parameters()
    refused.load_ini_text(blank)
    with pytest.raises(ValueError, match="cannot open"):
        pyec.make_potential(refused)
