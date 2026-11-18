"""In-process communicator: Matter end-to-end (skip without pyeonclient)."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import numpy as np
import pytest


# .con text carries six decimals, so a coordinate read back from it sits
# within half a unit of the sixth decimal of the frame it was written from.
_CON_TEXT_ATOL = 0.5e-6


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
    from readcon import ConFrame

    frame = record["product"]
    assert isinstance(frame, ConFrame)
    # saddle.con accompanies a saddle, as the drivers expect.
    assert ("saddle.con" in record) == ("saddle" in record)
    assert "pos.con" not in record
    from eon.structure import Structure
    import eon.fileio as fio

    back = Structure.from_conframe(frame)
    assert len(back) == 2
    assert np.isfinite(back.r).all()
    # min.con is the same geometry as .con text, for the drivers that read it.
    legacy = fio.loadcon(StringIO(record["min.con"].getvalue()))
    assert np.allclose(legacy.r, back.r, rtol=0.0, atol=_CON_TEXT_ATOL)
    if "saddle" in record:
        saddle = Structure.from_conframe(record["saddle"])
        legacy_saddle = fio.loadcon(StringIO(record["saddle.con"].getvalue()))
        assert np.allclose(legacy_saddle.r, saddle.r, rtol=0.0, atol=_CON_TEXT_ATOL)


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


def test_config_reads_state_session_reuse_switch(tmp_path):
    from eon.config import ConfigClass

    ini = tmp_path / "config.ini"
    ini.write_text("[Communicator]\nreuse_state_sessions = false\n")
    cfg = ConfigClass()
    cfg.init(str(ini))
    assert cfg.comm_reuse_state_sessions is False


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


def test_inprocess_stages_potential_files_in_evaluation_workdir(tmp_path):
    _need_client()
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import ConfigClass

    potfiles = tmp_path / "potfiles"
    potfiles.mkdir()
    (potfiles / "support.table").write_text("potential data\n")
    cfg = ConfigClass()
    cfg.path_scratch = str(tmp_path / "scratch")
    cfg.path_pot = str(potfiles)
    c = LocalInProcess(cfg.path_scratch, bundle_size=1, config=cfg)
    job = {"id": "support", "structure": _lj_structure()}
    ini = StringIO("[Main]\njob = point\n[Potential]\npotential = lj\n")
    original_cwd = Path.cwd()

    c.submit_jobs([job], {"config.ini": (ini, 0o644)})

    staged = Path(cfg.path_scratch) / ".inprocess-potential" / "support.table"
    assert staged.read_text() == "potential data\n"
    assert Path.cwd() == original_cwd


def test_process_searches_reuse_one_state_potential_and_minimized_reactant(
    tmp_path, monkeypatch,
):
    _need_client()
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import config

    monkeypatch.chdir(tmp_path)
    scratch = tmp_path / "scratch"
    config.path_scratch = str(scratch)
    c = LocalInProcess(str(scratch), bundle_size=1, config=config)
    ini = StringIO(
        "[Main]\njob = process_search\nrandom_seed = 1\n"
        "[Potential]\npotential = lj\n"
        "[Process Search]\nminimize_first = true\n"
        "[Saddle Search]\nmax_iterations = 1\n"
    )
    invariants = {
        "config.ini": (ini, 0o644),
        "pos.con": (_lj_con_text(), 0o644),
    }

    c.submit_jobs([{"id": "17_1"}], invariants)
    c.get_results()
    first_potential = c._state_sessions["17"]["potential"]
    first_calls = first_potential.force_call_counter

    invariants["config.ini"] = (
        StringIO(
            "[Main]\njob = process_search\nrandom_seed = 2\n"
            "[Potential]\npotential = lj\n"
            "[Process Search]\nminimize_first = true\n"
            "[Saddle Search]\nmax_iterations = 1\n"
        ),
        0o644,
    )
    invariants["pos.con"][0].seek(0)
    c.submit_jobs([{"id": "17_2"}], invariants)
    c.get_results()
    second_potential = c._state_sessions["17"]["potential"]

    assert second_potential is first_potential
    assert c.state_session_count == 1
    assert c.state_minimization_count("17") == 1
    assert second_potential.force_call_counter > first_calls
    assert list(tmp_path.glob("*.con")) == []


def test_process_search_session_reuse_can_be_disabled_for_comparison(tmp_path):
    _need_client()
    from eon.communicator_inprocess import LocalInProcess
    from eon.config import ConfigClass

    cfg = ConfigClass()
    cfg.path_scratch = str(tmp_path / "scratch")
    cfg.comm_reuse_state_sessions = False
    c = LocalInProcess(cfg.path_scratch, bundle_size=1, config=cfg)
    ini = StringIO(
        "[Main]\njob = process_search\n"
        "[Potential]\npotential = lj\n"
        "[Process Search]\nminimize_first = true\n"
        "[Saddle Search]\nmax_iterations = 1\n"
    )
    invariants = {
        "config.ini": (ini, 0o644),
        "pos.con": (_lj_con_text(), 0o644),
    }

    c.submit_jobs([{"id": "23_1"}], invariants)
    first = c.get_results()[0]["_matter"].get_potential()
    ini.seek(0)
    invariants["pos.con"][0].seek(0)
    c.submit_jobs([{"id": "23_2"}], invariants)
    second = c.get_results()[0]["_matter"].get_potential()

    assert second is not first
    assert c.state_minimization_count("23") == 2
