# 4,000-atom Cu benchmark

The input is the 10-by-10-by-10 FCC Cu V/SIA construction from amsel's
`scripts/cu-vsia/build_catalog.py`: one vacancy at half-lattice offset
`(2, 0, 0)` from a `[100]` split dumbbell. The lattice parameter is
3.615 angstrom, the initial dumbbell half-separation is 0.75 angstrom, and
the potential is Mishin `Cu01.eam.alloy`. All runs were serial with
`OMP_NUM_THREADS=1`.

## Command

```sh
CC=gcc CXX=g++ FC=gfortran meson setup build-inprocess \
  -Dwith_tests=true -Dwith_pyeonclient=true \
  -Dwith_fortran=false -Dwith_cuh2=false
meson compile -C build-inprocess
REPEATS=8 scripts/bench_rgpot_vs_socket_forces.sh
```

The script writes the generated input, exact configuration, logs, and CSV
rows under `benchmark-cu`. The worker saddle runs set
`EON_LAMMPS_WORKER=1`; the native runs use the production factory.

## Force-call latency

Each process performs two warmup calls before the eight recorded calls.
Times are milliseconds.

| repeat | forked worker | direct liblammps | native rgpot EAM |
| ---: | ---: | ---: | ---: |
| 0 | 3.880932 | 3.491696 | 7.577150 |
| 1 | 3.935077 | 3.760647 | 8.669386 |
| 2 | 3.809015 | 3.590590 | 9.621892 |
| 3 | 3.958825 | 3.771827 | 9.439038 |
| 4 | 3.926453 | 3.708053 | 9.258102 |
| 5 | 3.875352 | 3.796562 | 9.032919 |
| 6 | 4.013994 | 3.648566 | 6.340457 |
| 7 | 3.889316 | 3.855352 | 6.084263 |

| path | n | mean (ms) | SE (ms) |
| --- | ---: | ---: | ---: |
| forked worker | 8 | 3.911120500 | 0.021859689 |
| direct liblammps | 8 | 3.702911625 | 0.042320890 |
| native rgpot EAM | 8 | 8.252900875 | 0.498039845 |

The native and LAMMPS energies differ by exactly 0 eV for this frame. The
largest absolute force-component difference is
`9.62771529167e-16 eV/angstrom`, below the required `1e-8` tolerance.

## Process-search wall time

These are complete `process_search` launches at 600 K with seeds 8003
through 8008. Every run minimizes the same reactant and takes ten dimer
iterations. The worker and native runs have matching force-call counts and
matching printed reactant, saddle, and product energies.

| seed | worker wall (s) | native rgpot wall (s) | force calls |
| ---: | ---: | ---: | ---: |
| 8003 | 2.350708246231 | 1.431632041931 | 141 |
| 8004 | 2.269335508347 | 1.421497106552 | 146 |
| 8005 | 2.169852972031 | 1.273207426071 | 141 |
| 8006 | 2.224896430969 | 1.354187250137 | 141 |
| 8007 | 2.067152738571 | 1.303329944611 | 141 |
| 8008 | 2.282618999481 | 1.455107688904 | 143 |

| path | n | mean (s) | SE (s) |
| --- | ---: | ---: | ---: |
| forked worker | 6 | 2.227427483 | 0.040418542 |
| native rgpot EAM | 6 | 1.373160243 | 0.030385083 |

The independently supplied amsel AKMC baseline for the same atom count,
temperature, and seeds reports `1.700194743 s` per search with
`0.126061957 s` SE and `9.368992165 ms` wall time per force call with
`0.236220007 ms` SE. That AKMC baseline includes its surrounding session
lane, whereas the reproducible rows above isolate eOn launches.

## Reconciling 3.9, 9.37, and 38 ms per call

The three measurements include different boundaries:

- The `3.911 ms` worker microbenchmark above starts timing after two warmup
  calls. It excludes worker creation, liblammps/Open MPI initialization, and
  EAM table loading.
- The session lane's `9.369 ms` includes one approximately `1.17 s` LAMMPS
  startup per eonclient. Amortized over the profile's mean `177.5` calls,
  startup is `6.592 ms/call`. Adding the measured `3.242 ms` tight LAMMPS
  loop gives `9.834 ms/call`, consistent with that lane.
- The live throughput profile measured `6.726 s` blocked in the worker over
  `177.5` calls, or about `38 ms/call`. The production Matter path called
  `setFixedMask` for every force evaluation. The old `applySetforce` issued
  four failing `unfix` commands and four group deletions before discovering
  that all 12,000 coordinates were free. Those commands are absent from the
  force microbenchmark. After subtracting tight-loop and startup costs, this
  command replay accounts for the remaining approximately `28 ms/call` on
  the profiled build and load.

The Python server and job-directory work is a separate `0.666 s` per client
launch in the supplied profile. It is not a per-force-call IPC cost, and the
amsel state-exit bookkeeping was `0.049 s` over the complete six-step sample.
The fixed-mask hot path removes the repeated LAMMPS commands generically;
the in-process communicator removes the startup, pipe, and job-directory
boundaries.

The same command was repeated after the fixed-mask change. Eight-call means
were `3.859586125 ± 0.073382787 ms` for the warm worker,
`3.442153375 ± 0.037555698 ms` for direct liblammps, and
`6.059137875 ± 0.169135639 ms` for native rgpot EAM. Six complete
ten-iteration searches were `2.046097159 ± 0.036781452 s` through the worker
and `1.240384579 ± 0.029743410 s` through rgpot. Energy remained identical
and force-component L-infinity error remained
`9.62771529167e-16 eV/angstrom`. The variation from the earlier table is
machine load; both repeats satisfy the required `1e-8 eV` parity.

Native rgpot EAM remains the preferred default even though its isolated
force call is slower on this machine (`8.253 ms` versus `3.703 ms` for direct
liblammps in the table above). It avoids process startup and is the intended
library API. Optimization of its per-call kernel belongs in rgpot.

## End-to-end live AKMC throughput

The live comparison used the same 4,000-atom Cu V/SIA frame, Mishin table,
600 K temperature, and seeds 8003 through 8008. Each serial repetition ran
six real dimer searches through `akmc.py`, registered their real process
records, and took one amsel KMC step. Thus every mode has six independent
repetitions, 36 searches, and six KMC steps. This is a bounded one-state
throughput measurement, not a claim about time to recombination.

The three modes were:

- `old`: local communicator and `EON_LAMMPS_WORKER=1`;
- `uncached`: in-process communicator with
  `reuse_state_sessions = false`;
- `cached`: the new defaults, in-process communicator with
  `reuse_state_sessions = true`.

The current amsel trajectory runner was used only as the live driver. A
benchmark adapter changed the generated communicator type and cache switch
before `akmc.py` read the file. It removed the runner's newer
`nonlocal_distance_abort` and `nonlocal_count_abort` additions and disabled
its optional VTST Hessian pass, because neither was present in the supplied
throughput profile. These changes were identical across all three modes.
The search itself used the real native rgpot EAM or forked LAMMPS potential;
there were no mocks.

The exact common invocation was:

```sh
python3 scripts/eon/amsel_trajectory.py \
  --base-dir "$run" \
  --initial-config initial.con \
  --potfiles-dir potfiles \
  --potential lammps \
  --eonclient build-inprocess/client/eonclient \
  --akmc-py akmc_mode_adapter.py \
  --python "$(command -v python3)" \
  --temperature 600 \
  --seed "$seed" \
  --max-states 1 \
  --discovery-passes 1 \
  --jobs-per-pass 6 \
  --parallel-clients 1 \
  --candidate-atoms 3998,3999 \
  --min-mode-method dimer \
  --displace-magnitude 0.1 \
  --dimer-rotations-max 20 \
  --saddle-max-iterations 500 \
  --saddle-converged-force 0.01 \
  --no-use-vtst-correction \
  --no-kdb-superbasin \
  --no-amsel-bias-only \
  --no-close-range-recomb-neb \
  --max-consecutive-bad-passes 2 \
  --no-basin-hop-fallback \
  --stall-recovery-retries 0 \
  --discover-decide-fallback-n 1 \
  --no-resume-states
```

Runs were executed one at a time. For `old`, the environment selected
`EON_LAMMPS_WORKER=1` and left `type = local`. For the two new modes, the
adapter wrote `type = inprocess` and the stated cache switch. Wall time used
`time.monotonic_ns()` immediately around this command. Force calls are the
sum of the minimization, saddle, and prefactor columns in the registered
`search_results.txt`.

### Raw live results

| mode | seed | wall (s) | searches | KMC steps | force calls | searches/hour | KMC steps/hour |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| old | 8003 | 107.211097 | 6 | 1 | 2418 | 201.471682 | 33.578614 |
| old | 8004 | 122.973624 | 6 | 1 | 2781 | 175.647422 | 29.274570 |
| old | 8005 | 117.770702 | 6 | 1 | 2831 | 183.407245 | 30.567874 |
| old | 8006 | 94.234438 | 6 | 1 | 2007 | 229.215566 | 38.202594 |
| old | 8007 | 99.960931 | 6 | 1 | 2737 | 216.084421 | 36.014070 |
| old | 8008 | 119.045064 | 6 | 1 | 2796 | 181.443895 | 30.240649 |
| uncached | 8003 | 81.446837 | 6 | 1 | 2418 | 265.203671 | 44.200612 |
| uncached | 8004 | 114.247099 | 6 | 1 | 2787 | 189.063882 | 31.510647 |
| uncached | 8005 | 111.587771 | 6 | 1 | 2725 | 193.569598 | 32.261600 |
| uncached | 8006 | 80.078118 | 6 | 1 | 2007 | 269.736608 | 44.956101 |
| uncached | 8007 | 93.155060 | 6 | 1 | 2736 | 231.871462 | 38.645244 |
| uncached | 8008 | 120.244202 | 6 | 1 | 2907 | 179.634440 | 29.939073 |
| cached | 8003 | 73.581985 | 6 | 1 | 2038 | 293.550113 | 48.925019 |
| cached | 8004 | 91.405966 | 6 | 1 | 2407 | 236.308426 | 39.384738 |
| cached | 8005 | 89.852268 | 6 | 1 | 2345 | 240.394600 | 40.065767 |
| cached | 8006 | 62.089976 | 6 | 1 | 1627 | 347.882239 | 57.980373 |
| cached | 8007 | 72.948478 | 6 | 1 | 2356 | 296.099394 | 49.349899 |
| cached | 8008 | 91.717989 | 6 | 1 | 2527 | 235.504510 | 39.250752 |

### Mean and standard error

| mode | searches/hour | KMC steps/hour | wall/KMC step (s) | force calls/KMC step |
| --- | ---: | ---: | ---: | ---: |
| old | 197.878372 ± 8.756165 | 32.979729 ± 1.459361 | 110.199309 ± 4.716593 | 2595.000000 ± 132.656197 |
| uncached | 221.513277 ± 16.257595 | 36.918879 ± 2.709599 | 100.126514 ± 7.152954 | 2596.666667 ± 135.159823 |
| cached | 274.956547 ± 18.579169 | 45.826091 ± 3.096528 | 80.266110 ± 5.085055 | 2216.666667 ± 135.159823 |

Moving the communicator in-process, without caching, reduced mean wall time
per KMC step by `1.101x` and raised throughput by `1.119x`. Reusing the
potential and minimized reactant then added another `1.247x` wall-time gain,
`1.241x` throughput gain, and removed 380 mean force calls per step
(`14.6%`). Together, the new defaults reduced wall time by `1.373x` and
raised live search and KMC-step throughput by `1.390x`.
