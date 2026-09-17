"""In-memory ConFrame export: NEB path_frames (eOn-pxxt) and min/saddle frames (eOn-1li3)."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")
readcon = pytest.importorskip("readcon")


def _lj_params(**kwargs):
    p = pyec.Parameters()
    p.potential = pyec.PotType.LJ
    p.quiet = True
    p.write_log = False
    p.remove_net_force = True
    p.opt_max_iterations = kwargs.get("opt_max_iterations", 30)
    p.opt_converged_force = kwargs.get("opt_converged_force", 1e-2)
    for k, v in kwargs.items():
        if hasattr(p, k):
            setattr(p, k, v)
    return p


def _h2(pot, params, x2: float):
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.cell = np.eye(3, dtype=np.float64) * 20.0
    m.positions = np.array([[0.0, 0.0, 0.0], [x2, 0.0, 0.0]], dtype=np.float64)
    m.masses = np.array([1.008, 1.008], dtype=np.float64)
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.fixed = np.array([0, 0], dtype=np.int64)
    m.periodic = False
    return m


def test_path_frames_api_surface():
    assert hasattr(pyec.NudgedElasticBand, "path_frames")
    assert hasattr(pyec.NudgedElasticBand, "to_conframes")
    # pyec.NEB is a factory returning _core.NEB; path_frames lives on the class
    core_neb = pyec._core.NEB
    assert hasattr(core_neb, "path_frames")
    assert hasattr(core_neb, "to_conframes")
    assert hasattr(pyec.Matter, "movie_frames")
    assert hasattr(pyec.Matter, "take_movie_frames")
    assert hasattr(pyec.Matter, "clear_movie_frames")
    assert hasattr(pyec.MinModeSaddleSearch, "run_retain_frames")
    assert hasattr(pyec.MinModeSaddleSearch, "climb_frames")


def _short_neb(tmp_path):
    """Build LJ H2 NEB, short compute (no endpoint min). Returns (neb, params)."""
    params = _lj_params()
    params.job = pyec.JobType.Nudged_Elastic_Band
    params.neb_images = 3
    params.neb_minimize_endpoints = False
    params.neb_max_iterations = 5
    params.opt_max_iterations = 5
    params.opt_converged_force = 1.0  # cheap force target for short band
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    a = _h2(pot, params, 1.2)
    b = _h2(pot, params, 1.8)
    neb = pyec.NudgedElasticBand(a, b, params, pot)
    # Acceptance: path_frames after compute() (not only a single force update).
    status = neb.compute()
    assert status is not None
    return neb, params


def test_neb_path_frames_in_memory_no_neb_con(tmp_path, monkeypatch):
    """path_frames after compute(): n_path frames, NEB stamps, no durable neb.con."""
    monkeypatch.chdir(tmp_path)
    neb, params = _short_neb(tmp_path)

    before = set(tmp_path.iterdir())
    frames = neb.path_frames()
    after = set(tmp_path.iterdir())
    # path_frames must not create durable neb.con (temp handoff is process-private)
    durable_con = {p.name for p in after - before if p.suffix == ".con"}
    assert "neb.con" not in durable_con, f"path_frames left neb.con: {durable_con}"

    assert len(frames) == neb.n_path
    assert len(frames) == params.neb_images + 2
    for i, fr in enumerate(frames):
        assert fr.energy is not None
        assert np.isfinite(float(fr.energy))
        assert fr.frame_index == i
        assert fr.neb_bead == i
        # scalar stamps from neb_frame_metadata
        meta = fr.metadata
        assert "reaction_coordinate" in meta
        assert "relative_energy" in meta
        assert "parallel_force" in meta
        # geometry matches live path image
        img = neb.image(i)
        pos_fr = np.array([[at.x, at.y, at.z] for at in fr.atoms], dtype=np.float64)
        assert pos_fr.shape == (img.n_atoms, 3)
        assert np.allclose(pos_fr, np.asarray(img.positions), atol=1e-8)

    # alias
    frames2 = neb.to_conframes()
    assert len(frames2) == len(frames)


def test_neb_path_frames_roundtrip_vs_writepathcon(tmp_path, monkeypatch):
    """In-memory stamps match writePathCon / disk re-read via neb_write_results path."""
    monkeypatch.chdir(tmp_path)
    neb, params = _short_neb(tmp_path)

    mem = neb.path_frames()
    # writePathCon is exercised by neb_write_results → neb.con
    pyec.neb_write_results(neb, params, force_calls_neb=0)
    disk_path = tmp_path / "neb.con"
    assert disk_path.is_file()
    disk = readcon.read_con(str(disk_path))
    assert len(disk) == len(mem)
    for i, (mfr, dfr) in enumerate(zip(mem, disk)):
        assert mfr.frame_index == dfr.frame_index == i
        assert mfr.neb_bead == dfr.neb_bead == i
        assert mfr.energy == pytest.approx(dfr.energy)
        assert mfr.metadata["reaction_coordinate"] == pytest.approx(
            dfr.metadata["reaction_coordinate"]
        )
        assert mfr.metadata["relative_energy"] == pytest.approx(
            dfr.metadata["relative_energy"]
        )
        assert mfr.metadata["parallel_force"] == pytest.approx(
            dfr.metadata["parallel_force"]
        )
        pos_m = np.array([[at.x, at.y, at.z] for at in mfr.atoms])
        pos_d = np.array([[at.x, at.y, at.z] for at in dfr.atoms])
        assert np.allclose(pos_m, pos_d, atol=1e-8)


def test_matter_movie_frames_without_write_movie(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    params = _lj_params(opt_max_iterations=15, opt_converged_force=1e-3)
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = _h2(pot, params, 1.5)
    # stretch so minimize takes steps
    m.positions = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]], dtype=np.float64)

    before_con = {p.name for p in tmp_path.iterdir() if p.suffix == ".con"}
    out, ok = m.relax(inplace=True, quiet=True, write_movie=False, retain_frames=True)
    after_con = {p.name for p in tmp_path.iterdir() if p.suffix == ".con"}
    assert after_con == before_con, (
        f"retain_frames must not write movie .con files: {after_con - before_con}"
    )

    frames = m.movie_frames()
    assert len(frames) >= 1
    for i, fr in enumerate(frames):
        assert fr.frame_index == i
        assert fr.energy is not None
        assert np.isfinite(float(fr.energy))
        meta = fr.metadata
        assert "step_size" in meta
        assert "convergence" in meta

    taken = m.take_movie_frames()
    assert len(taken) == len(frames)
    assert len(m.movie_frames()) == 0


def test_matter_movie_frames_non_inplace_retain(tmp_path, monkeypatch):
    """Default relax path (inplace=False) must move retained frames to the returned Matter."""
    monkeypatch.chdir(tmp_path)
    params = _lj_params(opt_max_iterations=15, opt_converged_force=1e-3)
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = _h2(pot, params, 1.5)
    m.positions = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]], dtype=np.float64)
    # Caller unchanged; frames live on the returned working copy.
    out, ok = m.relax(inplace=False, quiet=True, write_movie=False, retain_frames=True)
    assert len(m.movie_frames()) == 0, "input Matter must not retain frames when inplace=False"
    frames = out.movie_frames()
    assert len(frames) >= 1, (
        "non-inplace relax(retain_frames=True) lost movie_frames "
        "(Matter must move ConFrames on return)"
    )
    assert frames[0].energy is not None
    assert "convergence" in frames[0].metadata
    assert "step_size" in frames[0].metadata


def test_saddle_climb_frames_without_write_movies(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    params = _lj_params()
    params.job = pyec.JobType.Saddle_Search
    params.saddle_max_iterations = 8
    params.opt_max_iterations = 8
    params.opt_converged_force = 0.5
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = _h2(pot, params, 1.4)
    # non-zero mode along bond stretch
    mode = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    e0 = float(m.potential_energy)
    ss = pyec.MinModeSaddleSearch(m, mode, e0, params, pot)

    before = set(tmp_path.iterdir())
    status = ss.run_retain_frames(max_iterations=5)
    after = set(tmp_path.iterdir())
    # climb.con / movies must not appear unless write_movies
    durable = {p.name for p in after - before if p.suffix in {".con", ".dat"}}
    assert "climb" not in durable and "climb.con" not in durable

    frames = ss.climb_frames()
    assert len(frames) >= 1, f"expected climb frames; status={status} msg={ss.status_message}"
    for i, fr in enumerate(frames):
        assert fr.frame_index == i
        assert fr.energy is not None
        meta = fr.metadata
        assert "step_size" in meta
        assert "delta_e" in meta
        assert "convergence" in meta
        assert "eigenvalue" in meta
