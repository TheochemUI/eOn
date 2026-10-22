import multiprocessing
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.locking import LockFile
from eon.structure import Structure


def test_acquirelock_writes_pid(tmp_path: Path):
    lock = LockFile(tmp_path / "lockfile")
    assert lock.acquirelock()
    assert lock.lock_path.is_file()
    assert int(lock.lock_path.read_text().strip()) == os.getpid()
    assert lock.islocked()
    second = LockFile(lock.lock_path)
    assert second.acquirelock() is False
    assert second.pid == os.getpid()
    lock.removelock()
    assert not lock.lock_path.exists()
    assert not lock.islocked()


def test_aquirelock_is_a_deprecated_alias(tmp_path: Path):
    lock = LockFile(tmp_path / "lockfile")
    with pytest.warns(DeprecationWarning, match="acquirelock"):
        assert lock.aquirelock()
    assert int(lock.lock_path.read_text().strip()) == os.getpid()
    lock.removelock()


def test_stale_lock_is_released(tmp_path: Path):
    dead = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    pid = dead.pid
    dead.kill()
    dead.wait()
    path = tmp_path / "lockfile"
    path.write_text("%i\n" % pid)
    lock = LockFile(path)
    assert not lock.islocked()
    assert not path.exists()


def test_stale_pid_is_replaced(tmp_path: Path):
    dead = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    pid = dead.pid
    dead.kill()
    dead.wait()
    path = tmp_path / "lockfile"
    path.write_text("%i\n" % pid)
    lock = LockFile(path)
    assert lock.acquirelock()
    assert path.read_text().strip() == str(os.getpid())
    lock.removelock()


def test_truncated_lock_is_replaced(tmp_path: Path):
    old = time.time() - 5
    path = tmp_path / "lockfile"
    path.write_bytes(b"12")
    os.utime(path, (old, old))
    lock = LockFile(path)
    assert lock.islocked() is False
    assert not path.exists()
    assert lock.acquirelock()
    assert path.read_text().strip() == str(os.getpid())
    lock.removelock()

    empty = tmp_path / "empty"
    empty.write_bytes(b"")
    os.utime(empty, (old, old))
    again = LockFile(empty)
    assert again.acquirelock()
    assert empty.read_text().strip() == str(os.getpid())
    again.removelock()


def test_reused_pid_without_flock_is_replaced(tmp_path: Path):
    path = tmp_path / "lockfile"
    path.write_text("%i\n" % os.getpid())
    lock = LockFile(path)
    assert lock.acquirelock()
    assert path.read_text().strip() == str(os.getpid())
    held = LockFile(path)
    assert held.acquirelock() is False
    lock.removelock()


def test_removelock_leaves_a_replaced_pid(tmp_path: Path):
    path = tmp_path / "lockfile"
    lock = LockFile(path)
    assert lock.acquirelock()
    path.write_text("424242\n")
    lock.removelock()
    assert path.read_text().strip() == "424242"


def test_removelock_missing_ok(tmp_path: Path):
    lock = LockFile(tmp_path / "absent")
    lock.removelock()
    assert not lock.lock_path.exists()


def _contend(path, begin, finished, slot, results):
    code = 0
    try:
        begin.wait(30)
        lock = LockFile(path)
        won = bool(lock.acquirelock())
        if won:
            won = lock.lock_path.read_text().strip() == str(os.getpid())
        results[slot] = 1 if won else 0
        finished.wait(30)
        if won:
            lock.removelock()
    except Exception:
        code = 1
        results[slot] = 0
    finally:
        os._exit(code)


def test_two_processes_one_owner(tmp_path: Path):
    ctx = multiprocessing.get_context("fork")
    for i in range(30):
        path = tmp_path / ("lock-%02d" % i)
        begin = ctx.Barrier(2)
        finished = ctx.Barrier(2)
        results = ctx.Array("i", [0, 0])
        procs = [
            ctx.Process(target=_contend, args=(path, begin, finished, slot, results))
            for slot in (0, 1)
        ]
        for proc in procs:
            proc.start()
        for proc in procs:
            proc.join(30)
            if proc.is_alive():
                proc.kill()
                proc.join(5)
            assert proc.exitcode == 0
        assert results[0] + results[1] == 1


def _lj_client_text():
    repo = str(Path(__file__).resolve().parents[1])
    return "\n".join(
        [
            "#!" + sys.executable,
            "import sys, time",
            "from pathlib import Path",
            "sys.path.insert(0, %r)" % repo,
            "import numpy as np",
            "from eon import fileio as io",
            "def lj(pos, epsilon=1.0, sigma=1.0):",
            "    energy = 0.0",
            "    force = np.zeros_like(pos)",
            "    n = len(pos)",
            "    for i in range(n):",
            "        for j in range(i + 1, n):",
            "            delta = pos[j] - pos[i]",
            "            dist = float(np.linalg.norm(delta))",
            "            if dist < 1e-8:",
            "                continue",
            "            inv = sigma / dist",
            "            sr6 = inv ** 6",
            "            sr12 = sr6 * sr6",
            "            energy += 4.0 * epsilon * (sr12 - sr6)",
            "            mag = 24.0 * epsilon * (sr6 - 2.0 * sr12) / dist",
            "            pair = mag * delta / dist",
            "            force[i] += pair",
            "            force[j] -= pair",
            "    return energy, force",
            "atoms = io.loadcon('pos.con')",
            "pos = np.array(atoms.r, dtype=float)",
            "energy = 0.0",
            "for _ in range(8):",
            "    energy, force = lj(pos)",
            "    pos = pos + 1e-4 * force",
            "energy, _force = lj(pos)",
            "atoms.r[:] = pos",
            "io.savecon('min.con', atoms)",
            "Path('results.dat').write_text(",
            "    '0 termination_reason\\n%.8f minimum_energy\\n8 total_force_calls\\n' % energy",
            ")",
            "time.sleep(1.5)",
            "",
        ]
    )


def _basin_hopping_run(root: Path):
    client = root.parent / "lj_client.py"
    client.write_text(_lj_client_text())
    client.chmod(0o755)
    atoms = Structure(2)
    atoms.names = ["H", "H"]
    atoms.mass[:] = 1.0
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[1] = [1.2, 0.0, 0.0]
    io.savecon(str(root / "pos.con"), atoms)
    ini = root / "config.ini"
    ini.write_text(
        "\n".join(
            [
                "[Main]",
                "job = basin_hopping",
                "temperature = 300",
                "random_seed = 1",
                "[Potential]",
                "potential = lj",
                "[Communicator]",
                "type = local",
                "client_path = %s" % client,
                "number_of_CPUs = 1",
                "num_jobs = 1",
                "jobs_per_bundle = 1",
                "[Basin Hopping]",
                "initial_state_pool_size = 0",
                "displacement = 0.05",
                "steps = 1",
                "[Paths]",
                "main_directory = %s" % root,
                "results = %s" % root,
                "states = %s" % (root / "states"),
                "[Optimizer]",
                "opt_method = cg",
                "max_iterations = 20",
                "converged_force = 0.5",
                "",
            ]
        )
    )
    return ini


def _drive_basin_hopping(ini, begin, codes, slot, err_path):
    code = 99
    try:
        os.chdir(Path(ini).parent)
        sys.stderr = open(err_path, "w", buffering=1)
        sys.stdout = open(os.devnull, "w")
        sys.argv = ["eon", str(ini)]
        begin.wait(60)
        from eon.basinhopping import main as hop_main

        try:
            hop_main()
            code = 0
        except SystemExit as exc:
            code = 0 if exc.code is None else int(exc.code)
    except Exception:
        code = 99
    finally:
        codes[slot] = code
        os._exit(0)


def test_two_basin_hopping_drivers_one_lock(tmp_path: Path):
    root = tmp_path / "lj"
    root.mkdir()
    ini = _basin_hopping_run(root)
    ctx = multiprocessing.get_context("fork")
    begin = ctx.Barrier(2)
    codes = ctx.Array("i", [-1, -1])
    procs = []
    for slot in (0, 1):
        err_path = tmp_path / ("driver-%d.err" % slot)
        proc = ctx.Process(
            target=_drive_basin_hopping,
            args=(str(ini), begin, codes, slot, str(err_path)),
        )
        procs.append(proc)
        proc.start()
    for proc in procs:
        proc.join(60)
        if proc.is_alive():
            proc.kill()
            proc.join(5)
        assert proc.exitcode == 0
    assert sorted(codes[:]) == [0, 1]
    logs = [(tmp_path / ("driver-%d.err" % slot)).read_text() for slot in (0, 1)]
    log_file = root / "bh.log"
    if log_file.is_file():
        logs.append(log_file.read_text())
    assert any("Couldn't get lock" in text for text in logs)
    assert (root / "wuid.dat").is_file()
    result = root / "jobs" / "scratch" / "0" / "results.dat"
    assert result.is_file()
    assert "minimum_energy" in result.read_text()
