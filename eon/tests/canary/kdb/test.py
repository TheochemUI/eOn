#!/usr/bin/env python

import os
import sys
from pathlib import Path

sys.path.insert(0, "../../")
from ndiff import ndiff

test_path = str(Path(__file__).resolve().parent)
test_name = Path(test_path).name


def _run(command):
    retval = os.system(command)
    if retval:
        print("%s: problem running %s" % (test_name, command))
        sys.exit(1)


def _saw_kdb_search():
    for path in Path(".").rglob("search_results.txt"):
        for line in path.read_text(errors="replace").splitlines():
            parts = line.split()
            if len(parts) > 1 and parts[1] == "kdb":
                return True
    return False


_run("eon-server --reset --force --quiet")
for i in range(15):
    _run("eon-server --quiet")

if not _saw_kdb_search():
    print("%s: no search has type kdb" % test_name)
    sys.exit(1)

same, max_rel_err, reason = ndiff("dynamics.test", "dynamics.txt", 0.01)
if same:
    print("%s: passed maximum relative error of %.3e" % (test_name, max_rel_err))
    _run("eon-server --reset --force --quiet")
    sys.exit(0)
else:
    print("%s: failed %s" % (test_name, reason))
    sys.exit(1)
