import atexit
import os
import sys
from pathlib import Path


class LockFile:
    def __init__(self, lock_path):
        self.lock_path = Path(lock_path)
        self.pid = None

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

        self.lock_path.write_text("%i\n" % os.getpid())
        atexit.register(self.removelock)

        return True

    def removelock(self):
        self.lock_path.unlink(missing_ok=True)
