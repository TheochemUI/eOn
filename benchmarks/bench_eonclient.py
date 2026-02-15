"""ASV benchmarks for eonclient binary performance."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data" / "one_pt_saddle_search"
INPUT_FILES = ["config.ini", "pos.con", "displacement.con", "direction.dat"]


class TimeSaddleSearchMorseDimer:
    """Benchmark a single saddle search with Morse potential and dimer method."""

    timeout = 120
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_eon_")
        for fname in INPUT_FILES:
            shutil.copy2(DATA_DIR / fname, self.tmpdir)

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_saddle_search_dimer(self):
        """Wall-clock time for one eonclient saddle search invocation."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )

    def peakmem_saddle_search_dimer(self):
        """Peak memory of one eonclient saddle search invocation."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )
