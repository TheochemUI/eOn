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
