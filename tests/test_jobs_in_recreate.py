"""Recreating jobs_in must tolerate a directory that is already there."""

from pathlib import Path
from types import SimpleNamespace

from eon.basinhopping import register_results as hop_register
from eon.escaperate import register_results as escape_register
from eon.explorer import ClientMinModeExplorer, ServerMinModeExplorer
from eon.parallelreplica import register_results as replica_register


class _Comm:
    def get_results(self, path, keep):
        return iter(())


def _reappearing_rmtree(module_shutil, monkeypatch):
    real = module_shutil.rmtree

    def rmtree(path, *args, **kwargs):
        real(path, *args, **kwargs)
        Path(path).mkdir(parents=True)

    monkeypatch.setattr(module_shutil, "rmtree", rmtree)


def test_hop_register_tolerates_jobs_in_reappearing(tmp_path, monkeypatch):
    import eon.basinhopping as hop

    jobs = tmp_path / "jobs_in"
    jobs.mkdir()
    _reappearing_rmtree(hop.shutil, monkeypatch)
    hop_register(_Comm(), None, SimpleNamespace(path_jobs_in=str(jobs)))
    assert jobs.is_dir()


def test_replica_register_tolerates_jobs_in_reappearing(tmp_path, monkeypatch):
    import eon.parallelreplica as replica

    jobs = tmp_path / "jobs_in"
    jobs.mkdir()
    _reappearing_rmtree(replica.shutil, monkeypatch)
    replica_register(
        _Comm(), None, None, SimpleNamespace(path_jobs_in=str(jobs))
    )
    assert jobs.is_dir()


def test_escape_register_tolerates_jobs_in_reappearing(tmp_path, monkeypatch):
    import eon.escaperate as escape

    jobs = tmp_path / "jobs_in"
    jobs.mkdir()
    _reappearing_rmtree(escape.shutil, monkeypatch)
    escape_register(
        _Comm(), None, None, SimpleNamespace(path_jobs_in=str(jobs))
    )
    assert jobs.is_dir()


def test_server_explorer_tolerates_jobs_in_reappearing(tmp_path, monkeypatch):
    import eon.explorer as explorer

    jobs = tmp_path / "jobs_in"
    jobs.mkdir()
    incomplete = tmp_path / "incomplete"
    _reappearing_rmtree(explorer.shutil, monkeypatch)
    ex = ServerMinModeExplorer.__new__(ServerMinModeExplorer)
    ex.config = SimpleNamespace(
        path_jobs_in=str(jobs),
        path_incomplete=str(incomplete),
        comm_job_bundle_size=1,
    )
    ex.comm = _Comm()
    ex.save = lambda: None
    ex.register_results()
    assert jobs.is_dir()


def test_explorer_mkdir_tolerates_an_existing_jobs_in(tmp_path, monkeypatch):
    jobs = tmp_path / "jobs_in"
    jobs.mkdir()
    real_is_dir = Path.is_dir
    seen = {"n": 0}

    def is_dir(self):
        if self == jobs:
            seen["n"] += 1
            # The check before mkdir misses the directory. Later checks,
            # including the one mkdir uses to honor exist_ok, see it.
            return seen["n"] > 2
        return real_is_dir(self)

    monkeypatch.setattr(Path, "is_dir", is_dir)
    ex = ClientMinModeExplorer.__new__(ClientMinModeExplorer)
    ex.config = SimpleNamespace(
        path_jobs_in=str(jobs),
        comm_type="local",
        comm_job_bundle_size=1,
    )
    ex.comm = _Comm()
    ex.superbasin = None
    ex.job_table = SimpleNamespace(write=lambda: None)
    ex.register_results()
    assert jobs.is_dir()
