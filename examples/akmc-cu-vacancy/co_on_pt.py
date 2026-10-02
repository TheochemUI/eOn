#!/usr/bin/env python3
"""Server entry for an adsorbate region with the element choice fixed.

The displacement hook runs ``python <script> <confile>`` and does not
forward flags. This wrapper adds them.
"""

import subprocess
import sys
from pathlib import Path

if len(sys.argv) < 2:
    raise SystemExit("usage: co_on_pt.py structure.con")

script = Path(__file__).with_name("adsorbate_region.py")
cmd = [
    sys.executable,
    str(script),
    sys.argv[1],
    "--adsorbate-elements",
    "C",
    "--adsorbate-elements",
    "O",
    "--cutoff",
    "4.0",
]
raise SystemExit(subprocess.call(cmd))
