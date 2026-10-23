import atexit
import os
import sys
import time
import warnings
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None


# A cut-off write is indistinguishable from a driver that has created the
# file and not yet stored its pid. A file that stays incomplete for this
# long was left by a driver that died mid-write.
_PARTIAL_GRACE_S = 1.0


def _write_all(fd, payload):
    view = memoryview(payload)
    while view:
        written = os.write(fd, view)
        view = view[written:]


def _process_alive(pid):
    """True when *pid* names a running process.

    On Windows the check is OpenProcess(SYNCHRONIZE). Elsewhere it is
    ``os.kill(pid, 0)``. PermissionError means the process exists and
    this user may not signal it.
    """
    if pid is None or pid <= 0:
        return False
    if sys.platform == "win32":
        import ctypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(0x100000, False, int(pid))  # SYNCHRONIZE
        if handle:
            kernel32.CloseHandle(handle)
            return True
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def _parse_pid(payload):
    """Return the pid stored in a finished lock file, or None.

    A finished file is ``"<pid>\\n"``. Anything else is empty or cut off.
    """
    if not payload or not payload.endswith(b"\n"):
        return None
    try:
        text = payload.decode("ascii").strip()
    except UnicodeDecodeError:
        return None
    if not text.isdigit():
        return None
    pid = int(text)
    if pid <= 0:
        return None
    return pid


def _cloexec(flags):
    if hasattr(os, "O_CLOEXEC"):
        return flags | os.O_CLOEXEC
    return flags


def _flock_exclusive(fd, blocking):
    """Take an exclusive flock.

    Returns True when this call holds the lock, False when a non-blocking
    attempt found another holder, and None where fcntl is missing.
    """
    if fcntl is None:
        return None
    try:
        if blocking:
            fcntl.flock(fd, fcntl.LOCK_EX)
        else:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return False
    return True


def _flock_release(fd):
    if fcntl is None:
        return
    try:
        fcntl.flock(fd, fcntl.LOCK_UN)
    except OSError:
        pass


class LockFile:
    def __init__(self, lock_path):
        self.lock_path = Path(lock_path)
        self.pid = None
        self._fd = None
        # Pid this object wrote. removelock unlinks only when this
        # process still owns that pid. A forked child copies the value,
        # so the pid check and the disarm flag both have to refuse it.
        self._owner_pid = None
        self._disarmed = False

    def islocked(self):
        kind, payload, holder = self._classify()
        if kind == "live":
            self.pid = holder
            return True
        if kind == "pending":
            # A driver has created the file and is still writing the pid.
            self.pid = holder
            return True
        if kind == "stale":
            self._unlink_if_same(payload)
        self.pid = None
        return False

    def acquirelock(self):
        """Create the lock, or take it over when the holder is gone.

        The create is ``O_CREAT|O_EXCL|O_WRONLY``. The pid is written and
        fsynced, and on platforms with fcntl the descriptor is held under
        an exclusive flock for the life of this process. If the file
        already exists, a dead pid or a truncated file is removed and the
        create is tried once more. A live pid whose flock is not held is
        the same as a dead holder: the pid was reused.
        """
        if self._try_create():
            return True
        return self._take_stale_once()

    def aquirelock(self):
        warnings.warn(
            "aquirelock is deprecated; use acquirelock",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.acquirelock()

    def _arm_fork_guard(self):
        if not hasattr(os, "register_at_fork"):
            return

        def after_child():
            atexit.unregister(self.removelock)
            self._disarmed = True

        os.register_at_fork(after_in_child=after_child)

    def removelock(self):
        """Unlink the lock when this process still owns it.

        A forked child inherits this method, the stored pid, and the
        atexit hook. Exit in the child must not delete the parent's
        file. A missing file, or a file another driver has replaced, is
        left alone. A LockFile that never acquired the lock does not
        delete whatever happens to sit at the path.
        """
        if self._disarmed:
            return
        owner = self._owner_pid
        if owner is None or os.getpid() != owner:
            return
        try:
            text = self.lock_path.read_text()
        except FileNotFoundError:
            text = None
        if text is not None and text.strip() == str(owner):
            try:
                self.lock_path.unlink()
            except FileNotFoundError:
                pass
        self._close()
        self._owner_pid = None
        try:
            atexit.unregister(self.removelock)
        except Exception:
            pass

    def _try_create(self):
        # Arm before the file exists so a fork during the write still
        # drops the exit handler in the child.
        self._disarmed = False
        self._arm_fork_guard()
        flags = _cloexec(os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        try:
            fd = os.open(os.fspath(self.lock_path), flags, 0o644)
        except FileExistsError:
            return False
        payload = b""
        try:
            # Flock before the pid is visible so a peer that opens the
            # empty file cannot treat it as abandoned.
            held = _flock_exclusive(fd, blocking=True)
            if fcntl is not None and held is not True:
                raise OSError("could not flock %s" % self.lock_path)
            payload = ("%i\n" % os.getpid()).encode("ascii")
            _write_all(fd, payload)
            os.fsync(fd)
        except Exception:
            os.close(fd)
            self._unlink_if_same(payload if payload else b"")
            raise
        self._fd = fd
        self._owner_pid = os.getpid()
        self.pid = self._owner_pid
        atexit.register(self.removelock)
        return True

    def _take_stale_once(self):
        deadline = time.monotonic() + _PARTIAL_GRACE_S
        while True:
            kind, payload, holder = self._classify()
            if kind == "missing":
                return self._try_create()
            if kind == "live":
                self.pid = holder
                return False
            if kind == "stale" or time.monotonic() >= deadline:
                self._unlink_if_same(payload)
                if self._try_create():
                    return True
                kind, _payload, holder = self._classify()
                if kind == "live":
                    self.pid = holder
                return False
            time.sleep(0.01)

    def _classify(self):
        """Return ``(kind, payload, pid)``.

        ``kind`` is ``missing``, ``live``, ``stale``, or ``pending``.
        ``payload`` is the bytes last read, so a later unlink can refuse
        to delete a file that changed. ``pid`` is set when the bytes name one.
        """
        flags = _cloexec(os.O_RDWR)
        try:
            fd = os.open(os.fspath(self.lock_path), flags)
        except FileNotFoundError:
            return "missing", None, None
        except PermissionError:
            flags = _cloexec(os.O_RDONLY)
            try:
                fd = os.open(os.fspath(self.lock_path), flags)
            except FileNotFoundError:
                return "missing", None, None
        held = None
        try:
            held = _flock_exclusive(fd, blocking=False)
            os.lseek(fd, 0, os.SEEK_SET)
            payload = os.read(fd, 128)
            if held:
                _flock_release(fd)
        finally:
            os.close(fd)

        pid = _parse_pid(payload)
        if pid is None:
            if held is False:
                return "pending", payload, None
            try:
                age = time.time() - self.lock_path.stat().st_mtime
            except FileNotFoundError:
                return "missing", None, None
            if age >= _PARTIAL_GRACE_S:
                return "stale", payload, None
            return "pending", payload, None
        if not _process_alive(pid):
            return "stale", payload, pid
        if held is True:
            # The pid is in the process table and nobody holds the flock.
            # The original owner has exited and the id was reused.
            return "stale", payload, pid
        return "live", payload, pid

    def _unlink_if_same(self, expected):
        if expected is None:
            return False
        try:
            data = self.lock_path.read_bytes()
        except FileNotFoundError:
            return False
        if data != expected:
            return False
        try:
            self.lock_path.unlink()
        except FileNotFoundError:
            return False
        return True

    def _close(self):
        fd = self._fd
        self._fd = None
        if fd is None:
            return
        try:
            os.close(fd)
        except OSError:
            pass
