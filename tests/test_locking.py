import os
import subprocess
import sys
from pathlib import Path

from eon.locking import LockFile


def test_aquirelock_writes_pid(tmp_path: Path):
    lock = LockFile(tmp_path / "lockfile")
    assert lock.aquirelock()
    assert lock.lock_path.is_file()
    assert int(lock.lock_path.read_text().strip()) > 0
    assert lock.islocked()
    assert not lock.aquirelock()
    lock.removelock()
    assert not lock.lock_path.exists()
    assert not lock.islocked()


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


def test_removelock_missing_ok(tmp_path: Path):
    lock = LockFile(tmp_path / "absent")
    lock.removelock()


def test_forked_child_exit_keeps_the_lock(tmp_path: Path):
    lock = LockFile(tmp_path / "lockfile")
    assert lock.aquirelock()
    path = lock.lock_path
    parent = os.getpid()
    pid = os.fork()
    if pid == 0:
        try:
            # The exit hook copied into the child must leave the file alone.
            lock.removelock()
            text = path.read_text().strip()
            ok = path.is_file() and text == str(parent) and text != str(os.getpid())
            if hasattr(os, "register_at_fork"):
                ok = ok and lock._disarmed
        except Exception:
            ok = False
        os._exit(0 if ok else 1)
    _child, status = os.waitpid(pid, 0)
    assert os.WIFEXITED(status)
    assert os.WEXITSTATUS(status) == 0
    assert path.is_file()
    assert int(path.read_text().strip()) == parent
    lock.removelock()
    assert not path.exists()
