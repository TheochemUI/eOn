"""A second harvest replaces a result directory instead of nesting the job."""

from types import SimpleNamespace

from eon.communicator import Local, Script


def _plant(root, name, text):
    job = root / name
    job.mkdir(parents=True)
    (job / "results.dat").write_text(text)
    return job


def test_local_result_replaces_an_existing_directory(tmp_path):
    scratch = tmp_path / "scratch"
    results = tmp_path / "jobs_in"
    _plant(scratch, "job 1", "1 termination_reason\n")
    _plant(results, "job 1", "stale\n")
    comm = Local(
        scratch,
        "/usr/bin/true",
        1,
        1,
        config=SimpleNamespace(debug_keep_all_results=False),
    )
    harvested = list(comm.get_results(results, lambda name: True))
    assert (results / "job 1" / "results.dat").read_text() == "1 termination_reason\n"
    assert not (results / "job 1" / "job 1").exists()
    assert harvested[0]["results.dat"].getvalue() == "1 termination_reason\n"


def test_script_result_replaces_an_existing_directory(tmp_path):
    scratch = tmp_path / "scratch"
    results = tmp_path / "jobs_in"
    _plant(scratch, "0_1", "1 termination_reason\n")
    _plant(results, "0_1", "stale\n")
    comm = Script.__new__(Script)
    comm.scratchpath = str(scratch)
    comm.jobids = {101: "1"}
    comm.config = SimpleNamespace(debug_keep_all_results=False)
    comm.get_queued_jobs = lambda: []
    harvested = list(comm.get_results(results, lambda name: True))
    assert (results / "0_1" / "results.dat").read_text() == "1 termination_reason\n"
    assert not (results / "0_1" / "0_1").exists()
    assert harvested[0]["name"] == "0_1"
    assert harvested[0]["results.dat"].getvalue() == "1 termination_reason\n"
