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
