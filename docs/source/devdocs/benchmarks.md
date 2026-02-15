---
myst:
  html_meta:
    "description": "Guide for developers on running and writing ASV benchmarks for eOn, including local usage, CI integration, and adding new benchmarks."
    "keywords": "eOn benchmarks, ASV, airspeed velocity, performance testing, CI"
---

# Benchmarks

```{versionadded} 2.5
```

We use [ASV (Airspeed Velocity)](https://asv.readthedocs.io/) to track
performance of the `eonclient` binary across commits. Benchmarks live in the
`benchmarks/` directory and are configured via `asv.conf.json` at the repository
root.

## Benchmark suite

The current suite covers four workloads, each measuring wall-clock time and peak
memory:

| Class | System | Job type |
|-------|--------|----------|
| `TimeSaddleSearchMorseDimer` | 337-atom Pt slab (Morse) | Saddle search (dimer) |
| `TimePointMorsePt` | 337-atom Pt slab (Morse) | Single point evaluation |
| `TimeMinimizationLJCluster` | 997-atom LJ cluster | LBFGS minimization |
| `TimeNEBMorsePt` | 337-atom Pt slab (Morse) | NEB (5 images) |

Input data for each benchmark is stored under `benchmarks/data/<name>/` and
contains a `config.ini` plus the necessary `.con` geometry files.

## Running locally

ASV expects `eonclient` to be on `PATH`. Build and install it first:

```bash
meson setup builddir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release
meson install -C builddir
```

Then install ASV and run the benchmarks against the current working tree:

```bash
pip install asv
asv machine --yes
asv run -E "existing:$(which python)" --set-commit-hash $(git rev-parse HEAD) --quick
```

The `--quick` flag runs each benchmark once. Drop it for full statistical
sampling (controlled by each class's `repeat` attribute).

To compare two commits:

```bash
asv compare <hash1> <hash2>
```

Results are stored in `.asv/results/` and can be browsed as HTML with:

```bash
asv publish
asv preview
```

## CI integration

Every pull request targeting `main` triggers the **Benchmark PR** workflow
(`.github/workflows/ci_benchmark.yml`). It:

1. Builds and installs `eonclient` at the `main` HEAD
2. Runs the full benchmark suite against `main`
3. Builds and installs `eonclient` at the PR HEAD
4. Runs the suite again against the PR
5. Compares the two runs using
   [asv-spyglass](https://github.com/hamdanal/asv-spyglass) and posts a
   summary table as a PR comment

The comment is updated in-place on subsequent pushes to the same PR.

## Adding a new benchmark

1. Create a data directory under `benchmarks/data/<name>/` containing a
   `config.ini` and any required `.con` files. You can reuse geometry files
   from `client/tests/` or `tests/data/`.

2. Add a new class in `benchmarks/bench_eonclient.py` following the existing
   pattern:

   ```python
   class TimeMyBenchmark:
       """Short description of what this benchmarks."""

       timeout = 120
       repeat = 5
       number = 1
       warmup_time = 0

       def setup(self):
           self.tmpdir = tempfile.mkdtemp(prefix="asv_eon_")
           _copy_data(BENCH_DATA / "my_benchmark", self.tmpdir)

       def teardown(self):
           shutil.rmtree(self.tmpdir, ignore_errors=True)

       def time_my_benchmark(self):
           """Wall-clock time."""
           subprocess.run(
               ["eonclient"],
               cwd=self.tmpdir,
               check=True,
               capture_output=True,
           )

       def peakmem_my_benchmark(self):
           """Peak memory."""
           subprocess.run(
               ["eonclient"],
               cwd=self.tmpdir,
               check=True,
               capture_output=True,
           )
   ```

3. Methods prefixed with `time_` measure wall-clock seconds; `peakmem_`
   measures peak RSS in bytes. ASV discovers these by convention.

4. Adjust `timeout` and `repeat` to match the expected cost of the workload.
   Cheap benchmarks (point evaluation) can use higher `repeat`; expensive ones
   (NEB, saddle search) should use lower values.
