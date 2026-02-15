
import os
import atexit

class LockFile:
    def __init__(self, lock_path):
        self.lock_path = lock_path
        self.pid = None

    def islocked(self):
        # check to see if the lock file is there
        if os.path.isfile(self.lock_path):
            # if the lock file is there lets get the pid from it
            f = open(self.lock_path)
            self.pid = int(f.read().strip())
            f.close()

            # check and see if process is still alive
            import sys
            if sys.platform == 'win32':
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

            # clean up a stale lock
            if not alive:
                self.removelock()
                return False
            return True
        else:
            return False

    def aquirelock(self):
        if self.islocked():
            return False

        fd = open(self.lock_path, 'w')
        fd.write("%i\n" % os.getpid())
        fd.close()
        atexit.register(self.removelock)

        return True

    def removelock(self):
        os.unlink(self.lock_path)
