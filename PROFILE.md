# Profile of the eOn server driver on a Cu vacancy-SIA cell

This records where the outer Python driver (`python -m eon config.ini`,
one `akmc()` cycle per invocation) spends its time before and after the
changes on this branch, and how the results were measured.

BEFORE is `253cb59a` (the base of this branch). AFTER is the tip of the
branch. Both ran from the same conda-style environment, with the same
eonclient, LAMMPS and amsel builds, on the same machine, interleaved
(base, new, base, new, ...).

## What was and was not reproduced

- amsel is 0.1.2 (the Nov 5 build from the store bundles), not bc9f0403.
  https://github.com/lode-org/amsel was not reachable from the build box.
- The in-process rgpot path is not available for LAMMPS potentials, so the
  live runs use the subprocess communicator (`Local.submit_jobs`).
- The live cell is the Cu vacancy-SIA cell of 2916 atoms with Mishin Cu01
  through LAMMPS, 16 searches per cycle. Eight cycles of it never reach a
  KMC step, so seconds per KMC step cannot be read from the live runs. It is
  measured with the synthetic driver below, which exercises the production
  `akmc`, `AKMCState`, `StateList`, superbasin and file code at state level
  with a deterministic lattice-walker standing in for the client.
- On the live cell the outer Python driver was 27% of wall before, not the 59%
  seen on the production runs. The shares below are therefore a lower bound
  for how much the driver matters on a longer, busier run.
- No hotspot was found in amsel. The amsel entry points on the profile
  (`amsel_discover_exit`, `discover_decide_for_superbasin`, the module import)
  are 0.3% to 2% of the driver and most of that is eon-side conversion.

## Live profile (py-spy, 4 runs each, 100 Hz, `--idle`)

Share of wall, inclusive. "Waiting on the client" is `Local.submit_jobs`
sleeping between polls plus `Popen.wait`.

| Where | before | after |
| --- | ---: | ---: |
| waiting on the client | 72.7% | 95.1% |
| Python driver (everything else) | 27.3% | 4.9% |
| `Displace.get_displacement` (`generate_displacement`) | 22.9% | 0.8% |
| `neighbor_list_vectors` | 21.3% | <0.1% |
| `pbc` | 6.3% | <0.1% |
| `neighbor_list` / `_pair_lists` | 3.3% / 2.9% | 0.4% / 0.4% |
| `ClientMinModeExplorer.register_results` | 0.95% | 0.73% |
| `savecon` / `loadcon` | 0.74% / 0.46% | 0.61% / 0.71% |
| `ConfigClass.__init__` (YAML) | 0.69% | <0.1% |
| `AKMCState.add_process` | 0.68% | 0.48% |
| `Structure.to_conframe` | 0.53% | 0.34% |

Before, the driver was dominated by displacement generation: for every
search it recomputed the pair search of the reactant and walked every pair in
Python (`neighbor_list_vectors`, about 0.3 s per search). After, the
neighbour table is built once per batch and the per-search work is a few
array operations.

## Synthetic driver profile (cProfile, 2916 atoms, 50 cycles, 16 searches each)

The harness part of each cycle (the stand-in client writing results) is
excluded. Driver time under cProfile: 23.5 s before, 7.1 s after.
Cumulative seconds and call counts of the driver entry points:

| Function | before | after |
| --- | ---: | ---: |
| `AKMCState.add_process` (577 calls) | 18.63 s | 5.45 s |
| `AKMCState.find_repeat` (241 calls) | 12.38 s | 2.54 s |
| `loadcon` (954 calls) | 14.89 s | 3.96 s |
| `Structure.from_conframe` (954 calls) | 13.65 s | 2.68 s |
| `akmc.kmc_step` (50 calls) | 4.38 s | 1.33 s |
| `StateList.get_product_state` (18 calls) | 3.80 s | 1.01 s |
| `State.get_reactant` (178 calls) | 2.86 s | 0.82 s |
| `atoms.match` (1010 calls) | 5.56 s | 0.56 s |
| `pbc` (2144 calls) | 5.20 s | 0.21 s |
| `atoms.remove_net_translation` (1010 calls) | 2.60 s | 0.25 s |
| `ini.get` (15699 calls; 16204 file parses before, 2395 after) | 1.65 s | 0.53 s |

## Timings

Live, one `python -m eon config.ini` invocation, 16 searches, n = 8
interleaved pairs, mean +- standard error:

| | before | after |
| --- | ---: | ---: |
| seconds per cycle | 17.39 +- 1.19 | 13.68 +- 1.14 |
| searches per hour | 3313 +- 226 | 4211 +- 349 |

The paired per-cycle saving is 3.71 +- 0.10 s, which is the 4.8 s the driver
cost before less the 0.7 s it costs now plus noise from the client.

Synthetic driver, n = 8 interleaved pairs, 50 cycles, 800 searches, 18 KMC
steps recorded, driver seconds only (no client):

| | before | after |
| --- | ---: | ---: |
| seconds per run | 11.60 +- 0.21 | 5.93 +- 0.07 |
| seconds per KMC step | 0.645 +- 0.012 | 0.329 +- 0.004 |
| seconds per cycle | 0.232 +- 0.004 | 0.1185 +- 0.0014 |
| searches per hour, driver bound | 248 000 +- 4 400 | 486 000 +- 5 700 |

The speedup of the synthetic driver is 1.96 +- 0.02 times (mean of the paired
ratios). The synthetic number is a ceiling for the searches per hour of the
driver alone: the force kernels and the client dominate the live figure.

## What changed

- `Displace.get_displacement` builds the neighbour vector table of the
  reactant once per batch (memoised on the structure, the radius and the
  box), and the void-bias force is computed on arrays with a norm and a
  cube that are bit-identical to the per-vector loop.
- `pbc` wraps an orthogonal cell in one array pass, after checking the edge
  lengths bit for bit against minimage on probe rows; other boxes keep the
  per-row path.
- The pair lists, per-center vectors and coordination numbers come from one
  sort of the pair keys.
- `Structure.from_conframe` and `to_conframe` build their arrays in bulk.
- `ini.get` re-parses an info file only when its bytes changed.
- `StateList.get_product_state` reads the state table once.
- `config.yaml` is read with libyaml's loader when PyYAML has it.

Every change keeps the output bit for bit: the tests compare against the
per-atom or per-pair definitions, and `test_akmc_driver_trajectory` replays a
recorded 150 cycle trajectory (state ids, rates, times, state table, process
tables, superbasin files) and requires the same bytes.

## Left alone

- `loadcon` is still 56% of the synthetic driver: each result is parsed
  (about 4 ms at 2916 atoms), and a new process parses its saddle twice (once
  in `find_repeat`, once to validate the result).
- `get_product_state` still opens every state's info file once per lookup.
- `find_repeat` parses and compares each close state in turn.
- `concorpus.mirror_con_path` re-inserts the same file (about 1.3 ms each).
