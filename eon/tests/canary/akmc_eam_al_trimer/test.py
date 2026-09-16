#!/usr/bin/env python

import os
import sys
from pathlib import Path
sys.path.insert(0, "../../")
from ndiff import ndiff

test_path = str(Path(__file__).resolve().parent)
test_name = Path(test_path).name

os.system("../../../akmc.py --reset --force --quiet")
for i in range(15):
    retval = os.system("python ../../../akmc.py --quiet")
    if retval:
        print("%s: problem running eon" % test_name)
        sys.exit(1)

same, max_rel_err, reason = ndiff("dynamics.test", "dynamics.txt", 0.01)
if same:
    print("%s: passed maximum relative error of %.3e"%(test_name,max_rel_err))
    os.system("python ../../../akmc.py --reset --force --quiet")
    sys.exit(0)
else:
    print("%s: failed %s"%(test_name,reason))
    sys.exit(1)
