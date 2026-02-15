"""ASV benchmarks for eonclient binary performance."""

import shutil
import subprocess
import tempfile
from pathlib import Path

BENCH_DATA = Path(__file__).parent / "data"


def _copy_data(data_dir, tmpdir):
    """Copy all files from a benchmark data directory to tmpdir."""
    for f in data_dir.iterdir():
        if f.is_file():
            shutil.copy2(f, tmpdir)


class TimeSaddleSearchMorseDimer:
    """Benchmark a single saddle search with Morse potential and dimer method."""

    timeout = 120
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_eon_")
        _copy_data(BENCH_DATA / "one_pt_saddle_search", self.tmpdir)

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


class TimePointMorsePt:
    """Benchmark a single point energy/force evaluation (337-atom Pt slab)."""

    timeout = 30
    repeat = 10
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_eon_")
        _copy_data(BENCH_DATA / "point_morse_pt", self.tmpdir)

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_point_evaluation(self):
        """Wall-clock time for a single point energy+force calculation."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )

    def peakmem_point_evaluation(self):
        """Peak memory for a single point calculation."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )


class TimeMinimizationLJCluster:
    """Benchmark geometry optimization of a 997-atom LJ cluster."""

    timeout = 300
    repeat = 3
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_eon_")
        _copy_data(BENCH_DATA / "min_lj_cluster", self.tmpdir)

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_minimization_lbfgs(self):
        """Wall-clock time for LBFGS minimization of LJ cluster."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )

    def peakmem_minimization_lbfgs(self):
        """Peak memory for LJ cluster minimization."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )


class TimeNEBMorsePt:
    """Benchmark a NEB calculation with Morse potential (337-atom Pt slab, 5 images)."""

    timeout = 300
    repeat = 3
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_eon_")
        _copy_data(BENCH_DATA / "neb_morse_pt", self.tmpdir)

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_neb(self):
        """Wall-clock time for NEB with 5 images."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )

    def peakmem_neb(self):
        """Peak memory for NEB calculation."""
        subprocess.run(
            ["eonclient"],
            cwd=self.tmpdir,
            check=True,
            capture_output=True,
        )
