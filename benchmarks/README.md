# eOn ASV benchmarks

[Airspeed Velocity (ASV)](https://asv.readthedocs.io/) suite that times the
`eonclient` binary (saddle search, point energy, minimization, nudged elastic
band, and related jobs). Entry points: `bench_eonclient.py` and `data/`.

Use this suite for client performance changes, including classical pots that
go through vesin for neighbor lists. Pull requests get a matrix against
`main`; continuous history lives on the `asv-results` branch; the dashboard
is at eondocs.org/bench.

## Which benches use vesin pots

Classical pair pots (LJ / Morse / LJCluster / ZBL force pairs) go through
`rgpot` over vesin:

| Class | Fixture | Potential |
|-------|---------|-----------|
| `TimePointMorsePt` | `data/point_morse_pt` | `morse_pt` (force NL) |
| `TimeSaddleSearchMorseDimer` | `data/one_pt_saddle_search` | `morse_pt` |
| `TimeNEBMorsePt` | `data/neb_morse_pt` | `morse_pt` |
| `TimeMinimizationLJCluster` | `data/min_lj_cluster` | `ljcluster` |

Python server geometry (`eon.geometry.neighbor_list`) is covered by unit tests
(`tests/test_geometry_vesin.py`). ASV here measures client binary wall time.

## Local

```bash
pixi install
pixi run meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release
pixi run meson install -C bbdir
pixi run bash -c 'pip install asv asv-tachyon'
pixi run asv machine --yes
pixi run asv run -E "existing:$(which python)" --quick
pixi run asv publish
pixi run asv-tachyon install .asv/html
pixi run asv-tachyon serve .asv/html --open
```

## History storage

Long-lived results live on the orphan git branch `asv-results`
(`.asv/results/**` only). The ASV dashboard workflow restores that branch
before each run and force-pushes updates after publish so history survives
Actions cache eviction.

## Dashboard CI

`.github/workflows/ci_asv_dashboard.yml` runs on `main`, weekly, and
`workflow_dispatch`. It:

1. Restores `asv-results`
2. Builds `eonclient`, runs ASV, publishes with [asv-tachyon](https://github.com/HaoZeke/asv_tachyon)
3. Pushes results back to `asv-results`
4. Deploys the UI to `gh-pages` at `/bench/`
5. Optionally deploys to Netlify (`NETLIFY_AUTH_TOKEN` + `NETLIFY_SITE_ID`) for bench.eondocs.org

| Surface | URL |
|---------|-----|
| Docs | https://eondocs.org/ |
| Bench (path) | https://eondocs.org/bench/ |
| Bench (host) | https://bench.eondocs.org/ |
| Fork preview | https://haozeke.github.io/eOn/bench/ |

PR comments use [asv-perch](https://github.com/HaoZeke/asv-perch) via
`ci_benchmark.yml` and `ci_bench_commenter.yml`.
