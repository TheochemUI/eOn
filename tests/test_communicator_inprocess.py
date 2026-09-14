"""In-process communicator: Matter end-to-end (skip without pyeonclient)."""

from __future__ import annotations

from io import StringIO

import numpy as np
import pytest

pc = pytest.importorskip("pyeonclient")


def _lj_pos_con() -> StringIO:
    """Tiny 2-atom .con via Structure + fileio."""
    from eon.structure import Structure
    from eon import fileio as fio

    s = Structure(2)
    s.r = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]])
    s.box = np.eye(3) * 20.0
    s.mass = np.array([1.0, 1.0])
    s.free = np.ones(2)
    s.names = ["H", "H"]
    buf = StringIO()
    fio.savecon(buf, s)
    buf.seek(0)
    return buf


def test_get_communicator_inprocess(tmp_path, monkeypatch):
    from eon.config import ConfigClass
    from eon import communicator as comm

    comm.reset_communicators()

    cfg = ConfigClass()
    cfg.comm_type = "inprocess"
    cfg.path_scratch = str(tmp_path / "scratch")
    cfg.comm_job_bundle_size = 1
    c = comm.get_communicator(cfg)
    assert type(c).__name__ == "LocalInProcess"


def test_inprocess_minimize_job(tmp_path):
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    config.path_scratch = str(tmp_path / "scratch")
    c = LocalInProcess(str(tmp_path / "scratch"), bundle_size=1, config=config)

    job = {"id": "t0", "pos.con": _lj_pos_con()}
    # minimal config.ini for LJ
    ini = StringIO(
        "[Main]\njob = minimization\n[Potential]\npotential = lj\n"
    )
    invariants = {"config.ini": (ini, 0o644)}
    c.submit_jobs([job], invariants)
    results = c.get_results()
    assert len(results) == 1
    r0 = results[0]
    assert "min.con" in r0
    assert "results.dat" in r0
    assert r0.get("_structure") is not None
    assert r0["_matter"].n_atoms == 2
    assert np.isfinite(r0["_energy"])


@pytest.mark.parametrize(
    "job_ini,expect",
    [
        ("[Main]\njob = point\n[Potential]\npotential = lj\n", "point"),
        ("[Main]\njob = hessian\n[Potential]\npotential = lj\n", "hessian"),
    ],
)
def test_inprocess_job_type_matrix(tmp_path, job_ini, expect):
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    config.path_scratch = str(tmp_path / "scratch")
    c = LocalInProcess(str(tmp_path / "scratch"), bundle_size=1, config=config)
    job = {"id": "t1", "pos.con": _lj_pos_con()}
    c.submit_jobs([job], {"config.ini": (StringIO(job_ini), 0o644)})
    r0 = c.get_results()[0]
    text = r0["results.dat"].getvalue()
    assert expect in text
    assert r0.get("_structure") is not None
