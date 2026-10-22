import atexit
import os
import sys
from pathlib import Path


class LockFile:
    def __init__(self, lock_path):
        self.lock_path = Path(lock_path)
        self.pid = None
        self._owner_pid = None
        self._disarmed = False

    def islocked(self):
        if not self.lock_path.is_file():
            return False
        self.pid = int(self.lock_path.read_text().strip())

        if sys.platform == "win32":
            import ctypes

            kernel32 = ctypes.windll.kernel32
            handle = kernel32.OpenProcess(0x100000, False, self.pid)  # SYNCHRONIZE
            if handle:
                kernel32.CloseHandle(handle)
                alive = True
            else:
                alive = False
        else:
            try:
                os.kill(self.pid, 0)
                alive = True
            except OSError:
                alive = False

        if not alive:
            self.removelock()
            return False
        return True

    def aquirelock(self):
        if self.islocked():
            return False

        self._owner_pid = os.getpid()
        self._disarmed = False
        # Arm before the file exists so a fork during the write still
        # drops the exit handler in the child.
        self._arm_fork_guard()
        self.lock_path.write_text("%i\n" % self._owner_pid)
        atexit.register(self.removelock)

        return True

    def _arm_fork_guard(self):
        if not hasattr(os, "register_at_fork"):
            return

        def after_child():
            atexit.unregister(self.removelock)
            self._disarmed = True

        os.register_at_fork(after_in_child=after_child)

    def removelock(self):
        # A forked child inherits this method and the pid written above.
        # Exit in the child must not delete the parent's lock file.
        if self._disarmed:
            return
        owner = self._owner_pid
        if owner is not None and os.getpid() != owner:
            return
        if owner is not None:
            try:
                text = self.lock_path.read_text()
            except FileNotFoundError:
                self._release_claim()
                return
            if text.strip() != str(owner):
                return
        self.lock_path.unlink(missing_ok=True)
        if owner is not None:
            self._release_claim()

    def _release_claim(self):
        self._owner_pid = None
        try:
            atexit.unregister(self.removelock)
        except Exception:
            pass
