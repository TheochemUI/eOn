"""In-process communicator: Matter end-to-end (skip without pyeonclient)."""

from __future__ import annotations

from io import StringIO

import numpy as np
import pytest


def _lj_structure():
    """Tiny 2-atom Structure (numpy working set)."""
    from eon.structure import Structure

    s = Structure(2)
    s.r = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]])
    s.box = np.eye(3) * 20.0
    s.mass = np.array([1.0, 1.0])
    s.free = np.ones(2)
    s.names = ["H", "H"]
    return s


def _lj_conframe():
    return _lj_structure().to_conframe()


def _lj_con_text():
    """The .con text the AKMC drivers send as pos.con."""
    import eon.fileio as fio

    buf = StringIO()
    fio.savecon(buf, _lj_structure())
    buf.seek(0)
    return buf


def _assert_product_frame(record):
    frame = record["product"]
    assert type(frame).__name__ == "ConFrame"
    assert "readcon" in type(frame).__module__
    assert "saddle.con" not in record
    assert "pos.con" not in record
    from eon.structure import Structure
    import eon.fileio as fio

    back = Structure.from_conframe(frame)
    assert len(back) == 2
    assert np.isfinite(back.r).all()
    # min.con is the same geometry as .con text, for the drivers that read it.
    legacy = fio.loadcon(StringIO(record["min.con"].getvalue()))
    assert np.allclose(legacy.r, back.r)


def _need_client():
    pytest.importorskip("pyeonclient")
    pytest.importorskip("readcon")


def test_geometry_structure_conframe_or_con_text():
    pytest.importorskip("readcon")
    from eon.communicator import CommunicatorError
    from eon.communicator_inprocess import _structure_from_job

    s = _lj_structure()
    assert _structure_from_job({"structure": s}, {}) is s
    assert _structure_from_job({"pos": s}, {}) is s
    assert _structure_from_job({}, {"reactant": s}) is s
    # explorer sends job['pos.con'] and invariants['pos.con'] as .con text.
    got = _structure_from_job({"pos.con": _lj_con_text()}, {})
    assert np.allclose(got.r, s.r)
    got = _structure_from_job({}, {"pos.con": (_lj_con_text(), 0o644)})
    assert np.allclose(got.r, s.r)
    with pytest.raises(CommunicatorError, match="could not read .con text"):
        _structure_from_job({"pos.con": StringIO("Lattice")}, {})
    with pytest.raises(CommunicatorError, match="Structure, ConFrame or .con text"):
        _structure_from_job({}, {})
    with pytest.raises(CommunicatorError, match="got int"):
        _structure_from_job({"structure": 1}, {})


def test_get_communicator_inprocess(tmp_path, monkeypatch):
    _need_client()
    from eon.config import ConfigClass
    from eon import communicator as comm

    comm.reset_communicators()

    cfg = ConfigClass()
    cfg.comm_type = "inprocess"
    cfg.path_scratch = str(tmp_path / "scratch")
    cfg.comm_job_bundle_size = 1
    c = comm.get_communicator(cfg)
    assert type(c).__name__ == "LocalInProcess"


@pytest.mark.parametrize("geometry", ["structure", "conframe", "con_text"])
def test_inprocess_minimize_job(tmp_path, geometry):
    _need_client()
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    scratch = tmp_path / "scratch"
    config.path_scratch = str(scratch)
    c = LocalInProcess(str(scratch), bundle_size=1, config=config)

    geom = {
        "structure": _lj_structure,
        "conframe": _lj_conframe,
        "con_text": _lj_con_text,
    }[geometry]()
    job = {"id": "t0", "pos.con": geom}
    ini = StringIO(
        "[Main]\njob = minimization\n[Potential]\npotential = lj\n"
    )
    invariants = {"config.ini": (ini, 0o644)}
    c.submit_jobs([job], invariants)
    results = c.get_results()
    assert len(results) == 1
    r0 = results[0]
    _assert_product_frame(r0)
    assert "results.dat" in r0
    assert r0.get("_structure") is not None
    assert r0["_matter"].n_atoms == 2
    assert np.isfinite(r0["_energy"])
    assert list(scratch.rglob("*.con")) == []


def test_inprocess_rejects_unreadable_con_text(tmp_path):
    _need_client()
    from eon.communicator import CommunicatorError
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    config.path_scratch = str(tmp_path / "scratch")
    c = LocalInProcess(str(tmp_path / "scratch"), bundle_size=1, config=config)
    job = {"id": "bad", "pos.con": StringIO("not a structure\n")}
    ini = StringIO("[Main]\njob = point\n[Potential]\npotential = lj\n")
    with pytest.raises(CommunicatorError, match="could not read .con text"):
        c.submit_jobs([job], {"config.ini": (ini, 0o644)})


@pytest.mark.parametrize(
    "job_ini,expect",
    [
        ("[Main]\njob = point\n[Potential]\npotential = lj\n", "point"),
        ("[Main]\njob = hessian\n[Potential]\npotential = lj\n", "hessian"),
        (
            "[Main]\njob = process_search\n[Potential]\npotential = lj\n"
            "[Saddle Search]\nmax_iterations = 2\n",
            "process_search",
        ),
    ],
)
def test_inprocess_job_type_matrix(tmp_path, job_ini, expect):
    _need_client()
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    scratch = tmp_path / "scratch"
    config.path_scratch = str(scratch)
    c = LocalInProcess(str(scratch), bundle_size=1, config=config)
    job = {"id": "t1", "structure": _lj_structure()}
    c.submit_jobs([job], {"config.ini": (StringIO(job_ini), 0o644)})
    r0 = c.get_results()[0]
    text = r0["results.dat"].getvalue()
    assert expect in text
    _assert_product_frame(r0)
    assert r0.get("_structure") is not None
    if r0.get("saddle") is not None:
        assert type(r0["saddle"]).__name__ == "ConFrame"
    assert list(scratch.rglob("*.con")) == []


def test_inprocess_conframe_key(tmp_path):
    _need_client()
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    scratch = tmp_path / "scratch"
    config.path_scratch = str(scratch)
    c = LocalInProcess(str(scratch), bundle_size=1, config=config)
    job = {"id": "t2", "conframe": _lj_conframe()}
    ini = StringIO("[Main]\njob = point\n[Potential]\npotential = lj\n")
    c.submit_jobs([job], {"config.ini": (ini, 0o644)})
    r0 = c.get_results()[0]
    _assert_product_frame(r0)
    assert "saddle" not in r0
    assert list(scratch.rglob("*")) == []
