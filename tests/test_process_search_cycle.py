"""A process search walks saddle, reactant minimum, and product minimum."""

from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from eon import fileio as io
from eon.explorer import ProcessSearch
from eon.structure import Structure


def _pair(atom1):
    atoms = Structure(2)
    atoms.names = ["Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[0] = [0.0, 0.0, 0.0]
    atoms.r[1] = np.asarray(atom1, dtype=float)
    return atoms


def _con(atoms):
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    return buf


def _mode_file(atoms):
    mode = np.zeros_like(atoms.r)
    mode[1] = [0.1, 0.0, 0.0]
    buf = StringIO()
    io.save_mode(buf, mode)
    buf.seek(0)
    return mode, buf


def _result(name, job_type, atoms, energy, code=0):
    mode, mode_buf = _mode_file(atoms)
    record = {
        "name": name,
        "results": {
            "job_type": job_type,
            "termination_reason": code,
            "potential_energy": energy,
            "potential_energy_reactant": -1.0,
            "potential_energy_saddle": -0.7,
            "total_force_calls": 4,
        },
        "saddle.con": _con(atoms),
        "mode.dat": mode_buf,
        "min.con": _con(atoms),
    }
    return record, mode


def test_search_connects_reactant_to_product(tmp_path: Path):
    reactant = _pair([2.5, 0.0, 0.0])
    product = _pair([2.5, 0.6, 0.0])
    saddle = _pair([2.5, 0.3, 0.0])
    mode, _mode_buf = _mode_file(reactant)
    config_path = tmp_path / "config.ini"
    config_path.write_text("[Main]\njob = process_search\ntemperature = 300\n")
    cfg = SimpleNamespace(
        process_search_default_prefactor=1.0e12,
        path_incomplete=str(tmp_path / "incomplete"),
        config_path=str(config_path),
        process_search_minimization_offset=0.01,
        comp_eps_r=0.05,
        comp_neighbor_cutoff=3.3,
        comp_check_rotation=False,
        comp_use_identical=False,
        comp_remove_translation=True,
    )
    search = ProcessSearch(
        reactant, reactant.copy(), mode, "random", 7, 0, config=cfg
    )
    job, kind = search.get_job(0)
    assert kind == "saddle_search"
    assert "pos.con" in job
    assert "displacement.con" in job

    saddle_record, _ignored = _result("0_1", "saddle_search", saddle, -0.7)
    assert search.process_result(saddle_record) is None
    assert search.data["barrier_reactant_to_product"] == pytest.approx(0.3)

    job, kind = search.get_job(0)
    assert kind == "min1"
    min1, _ignored = _result("0_2", "minimization", reactant, -1.0)
    assert search.process_result(min1) is None

    job, kind = search.get_job(0)
    assert kind == "min2"
    min2, _ignored = _result("0_3", "minimization", product, -0.9)
    final = search.process_result(min2)
    assert final is not None
    assert final["results"]["termination_reason"] == 0
    assert final["results"]["potential_energy_product"] == -0.9
    assert final["results"]["barrier_reactant_to_product"] == pytest.approx(0.3)
    assert final["results"]["barrier_product_to_reactant"] == pytest.approx(0.2)
    assert "product.con" in final
