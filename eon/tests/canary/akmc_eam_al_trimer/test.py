#!/usr/bin/env python

import os
import sys
sys.path.insert(0, "../../")
from ndiff import ndiff

test_path = os.path.split(os.path.realpath(__file__))[0]
test_name = os.path.basename(test_path)

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
