import atexit
import os
import socket
import sys
import time
import warnings
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None

try:
    import msvcrt
except ImportError:  # pragma: no cover - POSIX
    msvcrt = None


# A cut-off write is indistinguishable from a driver that has created the
# file and not yet stored its pid. A file that stays incomplete for this
# long was left by a driver that died mid-write.
_PARTIAL_GRACE_S = 1.0


class _Seen:
    def __init__(self, kind, payload, pid, flock_held):
        self.kind = kind
        self.payload = payload
        self.pid = pid
        self.flock_held = flock_held


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


def _host_token():
    raw = socket.gethostname().split(".")[0]
    raw = raw.encode("ascii", "replace").decode("ascii")
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in raw)
    return safe or "localhost"


def _owner_bytes():
    """Bytes stored in a live lock: ``"<pid> <host>\\n"``."""
    return ("%d %s\n" % (os.getpid(), _host_token())).encode("ascii")


def _parse_owner(payload):
    """Return ``(pid, host)`` from a finished lock file.

    A finished file is ``"<pid>\\n"`` or ``"<pid> <host>\\n"``. A cut-off
    write has no trailing newline. ``host`` is ``""`` when the file was
    written by an older driver that stored only the pid.
    """
    if not payload or not payload.endswith(b"\n"):
        return None, ""
    try:
        text = payload.decode("ascii").strip()
    except UnicodeDecodeError:
        return None, ""
    parts = text.split()
    if not parts or not parts[0].isdigit():
        return None, ""
    pid = int(parts[0])
    if pid <= 0:
        return None, ""
    host = parts[1] if len(parts) > 1 else ""
    return pid, host


def _parse_pid(payload):
    """Return the pid stored in a finished lock file, or None."""
    pid, _host = _parse_owner(payload)
    return pid


def _cloexec(flags):
    if hasattr(os, "O_CLOEXEC"):
        return flags | os.O_CLOEXEC
    return flags


def _have_kernel_lock():
    return fcntl is not None or msvcrt is not None


def _kernel_lock(fd, blocking):
    """Take an exclusive kernel lock on *fd*.

    Returns True when this call holds the lock, False when a non-blocking
    attempt found another holder, and None where no kernel lock exists.
    """
    if fcntl is not None:
        try:
            flags = fcntl.LOCK_EX
            if not blocking:
                flags |= fcntl.LOCK_NB
            fcntl.flock(fd, flags)
        except OSError:
            return False
        return True
    if msvcrt is not None:
        try:
            os.lseek(fd, 0, os.SEEK_SET)
            mode = msvcrt.LK_LOCK if blocking else msvcrt.LK_NBLCK
            msvcrt.locking(fd, mode, 1)
        except OSError:
            return False
        return True
    return None


def _kernel_unlock(fd):
    if fcntl is not None:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        except OSError:
            pass
        return
    if msvcrt is not None:
        try:
            os.lseek(fd, 0, os.SEEK_SET)
            msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
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
        info = self._inspect()
        if info.kind == "held" or info.kind == "pending":
            self.pid = info.pid
            return True
        if info.kind == "stale":
            self._drop_stale(info.payload)
        self.pid = None
        return False

    def acquirelock(self):
        """Create the lock, or take it over when the holder is gone.

        The create is ``O_CREAT|O_EXCL|O_RDWR``. The pid and host are
        written and fsynced. On platforms with flock or msvcrt the
        descriptor is held under an exclusive lock for the life of this
        process. If the file already exists, a dead pid or a truncated
        file is removed only after a hard link proves it is still that
        same inode, and the create is tried once more. A live pid whose
        kernel lock is not held is the same as a dead holder: the pid
        was reused. Where no kernel lock exists, the exclusive create
        is the lock.
        """
        self._disarmed = False
        self._arm_fork_guard()
        deadline = time.monotonic() + _PARTIAL_GRACE_S
        while True:
            if self._try_create():
                return True
            info = self._inspect()
            if info.kind == "missing":
                continue
            if info.kind == "held":
                self.pid = info.pid if info.pid is not None else 0
                return False
            if info.kind == "pending":
                if info.flock_held and time.monotonic() >= deadline:
                    self.pid = info.pid if info.pid is not None else 0
                    return False
                if time.monotonic() < deadline:
                    time.sleep(0.01)
                    continue
            if info.payload is None or not self._drop_stale(info.payload):
                info = self._inspect()
                if info.kind == "held":
                    self.pid = info.pid if info.pid is not None else 0
                    return False
                if info.kind == "missing":
                    continue
                return False

    def aquirelock(self):
        warnings.warn(
            "aquirelock is deprecated; use acquirelock",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.acquirelock()

    def __enter__(self):
        if not self.acquirelock():
            raise BlockingIOError("Couldn't get lock")
        return self

    def __exit__(self, exc_type, exc, tb):
        self.removelock()
        return False

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
            payload = self.lock_path.read_bytes()
        except FileNotFoundError:
            payload = None
        if payload is not None and _parse_pid(payload) == owner:
            self._unlink_our_inode()
        self._close()
        self._owner_pid = None
        try:
            atexit.unregister(self.removelock)
        except Exception:
            pass

    def _try_create(self):
        flags = _cloexec(os.O_CREAT | os.O_EXCL | os.O_RDWR)
        try:
            fd = os.open(os.fspath(self.lock_path), flags, 0o644)
        except FileExistsError:
            return False
        try:
            # msvcrt locks a byte range, so the file needs a byte first.
            # flock locks the open file and does not.
            if msvcrt is not None and fcntl is None:
                _write_all(fd, b"\n")
                os.fsync(fd)
            held = _kernel_lock(fd, blocking=True)
            if _have_kernel_lock() and held is not True:
                raise OSError("could not lock %s" % self.lock_path)
            self._write_owner(fd)
        except Exception:
            self._discard_fd(fd)
            raise
        self._remember(fd)
        return True

    def _write_owner(self, fd):
        payload = _owner_bytes()
        os.ftruncate(fd, 0)
        os.lseek(fd, 0, os.SEEK_SET)
        _write_all(fd, payload)
        os.fsync(fd)

    def _remember(self, fd):
        self._fd = fd
        self._owner_pid = os.getpid()
        self.pid = self._owner_pid
        atexit.register(self.removelock)

    def _inspect(self):
        """Return what the lock file is right now.

        ``kind`` is ``missing``, ``held``, ``stale``, or ``pending``.
        ``payload`` is the bytes last read, so a later removal can refuse
        to delete a file that changed. ``flock_held`` is True when another
        process holds the kernel lock.
        """
        flags = _cloexec(os.O_RDWR)
        try:
            fd = os.open(os.fspath(self.lock_path), flags)
        except FileNotFoundError:
            return _Seen("missing", None, None, False)
        except PermissionError:
            flags = _cloexec(os.O_RDONLY)
            try:
                fd = os.open(os.fspath(self.lock_path), flags)
            except FileNotFoundError:
                return _Seen("missing", None, None, False)
        got = None
        try:
            got = _kernel_lock(fd, blocking=False)
            os.lseek(fd, 0, os.SEEK_SET)
            payload = os.read(fd, 256)
        finally:
            if got is True:
                _kernel_unlock(fd)
            os.close(fd)

        pid, _host = _parse_owner(payload)
        if got is False:
            if pid is not None and _process_alive(pid):
                return _Seen("held", payload, pid, True)
            return _Seen("pending", payload, pid, True)
        if pid is None:
            try:
                age = time.time() - self.lock_path.stat().st_mtime
            except FileNotFoundError:
                return _Seen("missing", None, None, False)
            if age < _PARTIAL_GRACE_S:
                return _Seen("pending", payload, None, False)
            return _Seen("stale", payload, None, False)
        if not _process_alive(pid):
            return _Seen("stale", payload, pid, False)
        if got is True:
            # The pid is in the process table and nobody holds the kernel
            # lock. The original owner has exited and the id was reused.
            return _Seen("stale", payload, pid, False)
        return _Seen("held", payload, pid, False)

    def _drop_stale(self, expected):
        """Remove *expected* only while the path still names that inode.

        A peer may already have published a new lock under this path.
        Linking the path and comparing inodes refuses to delete that
        replacement. Two recoveries of the same dead file both drop it,
        then exactly one exclusive create wins.
        """
        if expected is None:
            return False
        flags = _cloexec(os.O_RDWR)
        try:
            fd = os.open(os.fspath(self.lock_path), flags)
        except FileNotFoundError:
            return False
        claim = None
        got = None
        try:
            got = _kernel_lock(fd, blocking=False)
            if got is False:
                return False
            os.lseek(fd, 0, os.SEEK_SET)
            data = os.read(fd, 256)
            if data != expected:
                return False
            st = os.fstat(fd)
            claim = self._claim_path()
            try:
                os.link(os.fspath(self.lock_path), os.fspath(claim))
            except FileNotFoundError:
                return False
            st_claim = os.stat(claim)
            if (st_claim.st_ino, st_claim.st_dev) != (st.st_ino, st.st_dev):
                return False
            try:
                os.unlink(os.fspath(self.lock_path))
            except FileNotFoundError:
                pass
            try:
                st_path = os.stat(self.lock_path)
            except FileNotFoundError:
                st_path = None
            if st_path is not None and (st_path.st_ino, st_path.st_dev) == (
                st.st_ino,
                st.st_dev,
            ):
                os.unlink(os.fspath(self.lock_path))
            os.lseek(fd, 0, os.SEEK_SET)
            if os.read(fd, 256) != expected:
                try:
                    os.link(os.fspath(claim), os.fspath(self.lock_path))
                except FileExistsError:
                    pass
                return False
            return True
        finally:
            if got is True:
                _kernel_unlock(fd)
            os.close(fd)
            if claim is not None:
                try:
                    claim.unlink()
                except FileNotFoundError:
                    pass

    def _unlink_our_inode(self):
        fd = self._fd
        if fd is None:
            try:
                self.lock_path.unlink()
            except FileNotFoundError:
                pass
            return
        try:
            st = os.fstat(fd)
        except OSError:
            return
        try:
            st_path = os.stat(self.lock_path)
        except FileNotFoundError:
            return
        if (st.st_ino, st.st_dev) != (st_path.st_ino, st_path.st_dev):
            return
        claim = self._claim_path()
        try:
            os.rename(os.fspath(self.lock_path), os.fspath(claim))
        except FileNotFoundError:
            return
        try:
            st_claim = os.stat(claim)
        except FileNotFoundError:
            return
        if (st_claim.st_ino, st_claim.st_dev) == (st.st_ino, st.st_dev):
            try:
                os.unlink(os.fspath(claim))
            except FileNotFoundError:
                pass
            return
        try:
            os.rename(os.fspath(claim), os.fspath(self.lock_path))
        except OSError:
            pass

    def _claim_path(self):
        name = ".%s.%d.%d.claim" % (
            self.lock_path.name,
            os.getpid(),
            time.time_ns(),
        )
        return self.lock_path.with_name(name)

    def _discard_fd(self, fd):
        try:
            st = os.fstat(fd)
        except OSError:
            st = None
        try:
            os.close(fd)
        except OSError:
            pass
        if st is None:
            return
        claim = self._claim_path()
        try:
            os.rename(os.fspath(self.lock_path), os.fspath(claim))
        except FileNotFoundError:
            return
        try:
            st_claim = os.stat(claim)
        except FileNotFoundError:
            return
        if (st_claim.st_ino, st_claim.st_dev) == (st.st_ino, st.st_dev):
            try:
                os.unlink(os.fspath(claim))
            except FileNotFoundError:
                pass
            return
        try:
            os.rename(os.fspath(claim), os.fspath(self.lock_path))
        except OSError:
            pass

    def _close(self):
        fd = self._fd
        self._fd = None
        if fd is None:
            return
        try:
            os.close(fd)
        except OSError:
            pass
