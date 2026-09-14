"""Local communicator must not put client stderr on an undrained PIPE.

A client that writes more than one OS pipe buffer (~64 KiB on Linux) to an
undrained stderr pipe blocks in anon_pipe_write and never exits, while the
parent only polls p.poll(). Redirecting stderr to stderr.dat beside
stdout.dat removes the bound.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

# Bytes well above a typical Linux pipe buffer (64 KiB).
BIG = 1 << 20


CLIENT = (
    "import sys\n"
    "n = int(sys.argv[1])\n"
    "sys.stderr.write('x' * n)\n"
    "sys.stderr.flush()\n"
    "sys.stdout.write('done\\n')\n"
    "sys.stdout.flush()\n"
)


def _reap(p: subprocess.Popen, timeout_s: float = 5.0) -> bool:
    t0 = time.monotonic()
    while time.monotonic() - t0 < timeout_s:
        if p.poll() is not None:
            return True
        time.sleep(0.05)
    p.kill()
    try:
        p.communicate(timeout=2.0)
    except Exception:
        pass
    return False


def test_undrained_stderr_pipe_deadlocks(tmp_path: Path):
    """Documents the defect shape: undrained PIPE never exits past 64 KiB."""
    (tmp_path / "client.py").write_text(CLIENT)
    with open(tmp_path / "stdout.dat", "w") as out:
        p = subprocess.Popen(
            [sys.executable, "client.py", str(BIG)],
            cwd=tmp_path,
            stdout=out,
            stderr=subprocess.PIPE,
        )
        exited = _reap(p, timeout_s=3.0)
    assert not exited, (
        "undrained PIPE should hang once stderr exceeds the buffer"
    )


def test_stderr_to_file_completes(tmp_path: Path):
    """The fix shape: stderr to a file so nothing blocks on the pipe buffer."""
    (tmp_path / "client.py").write_text(CLIENT)
    with open(tmp_path / "stdout.dat", "w") as out, open(
        tmp_path / "stderr.dat", "w"
    ) as err:
        p = subprocess.Popen(
            [sys.executable, "client.py", str(BIG)],
            cwd=tmp_path,
            stdout=out,
            stderr=err,
        )
        exited = _reap(p, timeout_s=5.0)
    assert exited, "client did not exit with stderr redirected to a file"
    assert (tmp_path / "stderr.dat").stat().st_size == BIG


def test_local_submit_jobs_writes_stderr_file(tmp_path: Path, monkeypatch):
    """Local.submit_jobs must put client stderr on a file, not a PIPE."""
    from eon import communicator as comm
    from eon.config import ConfigClass

    client_py = tmp_path / "client.py"
    client_py.write_text(CLIENT)
    wrap = tmp_path / "client_wrap"
    wrap.write_text(
        "#!/usr/bin/env bash\n"
        f"exec {sys.executable} {client_py} {BIG}\n"
    )
    wrap.chmod(0o755)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    job = scratch / "0_0"
    job.mkdir()

    cfg = ConfigClass()
    cfg.path_scratch = str(scratch)
    cfg.comm_local_client = str(wrap)
    cfg.comm_local_ncpus = 1
    cfg.comm_job_bundle_size = 1

    local = comm.Local(
        str(scratch), str(wrap), ncpus=1, bundle_size=1, config=cfg
    )

    def fake_bundles(data, invariants):
        yield str(job)

    monkeypatch.setattr(local, "make_bundles", fake_bundles)
    t0 = time.monotonic()
    local.submit_jobs(data=[], invariants={})
    elapsed = time.monotonic() - t0

    assert elapsed < 5.0, f"submit_jobs hung for {elapsed:.1f}s"
    assert (job / "stderr.dat").is_file()
    assert (job / "stderr.dat").stat().st_size == BIG
    assert (job / "stdout.dat").is_file()
