"""Cluster (Script) and Local communicators under a fake batch system.

The fake queue is a file of job ids; submit appends, cancel removes, and
squeue prints it. That is enough to check that cancel_state cancels only
queued jobs, keeps finished results for harvest, survives a failing
cancel, and that Local cleanup reaches the client's own children.
"""
from __future__ import annotations

import os
import signal
import subprocess
import time
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from eon.communicator import Local, Script

pytestmark = pytest.mark.skipif(os.name != "posix", reason="POSIX shell scripts")


def _config(tmp_path):
    return SimpleNamespace(
        debug_keep_all_results=False,
        debug_results_path="debug",
        path_root=str(tmp_path),
        path_pot="",
    )


def _write(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(0o755)


def _fake_scripts(root: Path) -> Path:
    scripts = root / "scripts"
    scripts.mkdir()
    queue = root / "queue"
    queue.write_text("")
    _write(
        scripts / "submit_job.sh",
        f"""#!/bin/sh
n=$(cat "{root}/counter" 2>/dev/null || echo 100)
n=$((n + 1))
echo $n > "{root}/counter"
echo $n >> "{queue}"
echo $n
""",
    )
    _write(scripts / "queued_jobs.sh", f'#!/bin/sh\ncat "{queue}"\n')
    _write(
        scripts / "cancel_job.sh",
        f"""#!/bin/sh
[ "$1" = "${{FAIL_CANCEL:-}}" ] && {{ echo "cannot cancel $1" >&2; exit 1; }}
grep -v "^$1$" "{queue}" > "{queue}.new" || true
mv "{queue}.new" "{queue}"
""",
    )
    return scripts


def _script_comm(tmp_path):
    scripts = _fake_scripts(tmp_path)
    return Script(
        str(tmp_path / "scratch"),
        1,
        "eon",
        str(scripts),
        "queued_jobs.sh",
        "cancel_job.sh",
        "submit_job.sh",
        config=_config(tmp_path),
    )


def _jobs(n):
    return [{"id": f"0_{i}", "pos.con": StringIO("x\n")} for i in range(n)]


def _finish(tmp_path, jobid):
    """Take a job out of the fake queue and give it a result, as Slurm would."""
    queue = tmp_path / "queue"
    ids = [x for x in queue.read_text().split() if x != str(jobid)]
    queue.write_text("".join(i + "\n" for i in ids))


def test_cancel_state_cancels_queued_and_keeps_finished(tmp_path, monkeypatch):
    comm = _script_comm(tmp_path)
    comm.submit_jobs(_jobs(3), {})
    assert sorted(comm.jobids) == [101, 102, 103]
    assert comm.get_queue_size() == 3

    # Job 101 (eOn job 0_0) finishes and writes its result before the cancel.
    (tmp_path / "scratch" / "0_0" / "results.dat").write_text("0 termination_reason\n")
    _finish(tmp_path, 101)
    # The batch system refuses to cancel 103.
    monkeypatch.setenv("FAIL_CANCEL", "103")

    assert comm.cancel_state(0) == 1

    scratch = sorted(p.name for p in (tmp_path / "scratch").iterdir() if p.is_dir())
    # 0_1 was cancelled and removed; the finished 0_0 and the uncancelled 0_2 stay.
    assert scratch == ["0_0", "0_2"]
    assert sorted(comm.jobids) == [101, 103]

    results_dir = tmp_path / "results"
    results_dir.mkdir()
    got = list(comm.get_results(str(results_dir), lambda name: True))
    assert [r["name"] for r in got] == ["0_0"]
    assert sorted(comm.jobids) == [103]


def test_get_results_leaves_queued_jobs_in_scratch(tmp_path):
    comm = _script_comm(tmp_path)
    comm.submit_jobs(_jobs(2), {})
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    assert list(comm.get_results(str(results_dir), lambda name: True)) == []
    assert sorted(p.name for p in (tmp_path / "scratch").iterdir() if p.is_dir()) == ["0_0", "0_1"]


CLIENT = """#!/bin/sh
sleep 300 &
echo $! > grandchild.pid
wait
"""


def test_local_cleanup_kills_the_clients_children(tmp_path):
    client = tmp_path / "client.sh"
    _write(client, CLIENT)
    comm = Local(str(tmp_path / "scratch"), str(client), 1, 1, config=_config(tmp_path))
    job = tmp_path / "scratch" / "0_0"
    job.mkdir()
    p = subprocess.Popen(
        [str(client)], cwd=job, start_new_session=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    comm.joblist.append((p, str(job)))
    pidfile = job / "grandchild.pid"
    for _ in range(100):
        if pidfile.exists() and pidfile.read_text().strip():
            break
        time.sleep(0.05)
    grandchild = int(pidfile.read_text())

    comm.cleanup()
    p.wait(timeout=5)
    for _ in range(100):
        try:
            os.kill(grandchild, 0)
        except ProcessLookupError:
            break
        # A killed but unreaped child of init shows as a zombie; that is gone.
        stat = Path(f"/proc/{grandchild}/stat")
        if stat.exists() and stat.read_text().split()[2] == "Z":
            break
        time.sleep(0.05)
    else:
        os.kill(grandchild, signal.SIGKILL)
        pytest.fail("the client's child survived cleanup()")
