##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

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
