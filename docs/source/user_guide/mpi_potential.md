---
myst:
  html_meta:
    "description": "Instructions for using a modified version of VASP with eOn via the MPI Potential interface for parallel calculations."
    "keywords": "eOn MPI potential, VASP interface, parallel VASP, ab-initio"
---

# MPI Potential

```{admonition} conda-forge availability
:class: warning
**Not** included in the `conda-forge` package. Requires building from source
with `-Dwith_mpi=True`.
```

```{versionchanged} 2.15
Two changes that ship together:

1. **C bindings.** `MPIPot.cpp` and the MPI control flow in
   `ClientEON.cpp` were ported off the removed MPI C++ bindings
   (`MPI::COMM_WORLD`, `MPI::INT`, ...). The C bindings used now
   (`MPI_Send`, `MPI_Recv`, `MPI_Iprobe`, `MPI_Comm_create`, ...)
   are guaranteed by every spec-compliant MPI implementation, so
   `-Dwith_mpi=true` builds against modern conda-forge MPICH 4.x
   and OpenMPI 5.x without `mpicxx.h`.

2. **MPItrampoline by default.** The `-Dwith_mpi=true` build now
   prefers [MPItrampoline](https://github.com/eschnett/MPItrampoline)
   over a direct link to a specific MPI flavour. MPItrampoline
   provides a stable wrapper ABI -- think of it as FlexiBLAS for
   MPI -- so one eonclient binary works against any spec-compliant
   libmpi at run time. The eOn meson detection order is:

   1. `dependency('MPItrampoline', method: 'cmake')`
   2. `dependency('mpi-c')` (MPItrampoline ships a shim that
      shadows the system `mpi-c.pc`)
   3. `dependency('mpi')` (last-resort direct link to system MPI)

   `subprojects/mpitrampoline.wrap` pins upstream v5.5.1 for users
   who want meson to fetch and CMake-build MPItrampoline locally.
```

## MPItrampoline workflow

Install MPItrampoline once per system, then point it at the actual
MPI implementation you want to use at run time.

```{code-block} shell
# 1) Build MPItrampoline against the trampoline ABI (no MPI needed yet):
git clone https://github.com/eschnett/MPItrampoline
cmake -S MPItrampoline -B build-mpitrampoline \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX=$HOME/mpitrampoline
cmake --build build-mpitrampoline -j
cmake --install build-mpitrampoline

# 2) For each MPI implementation you want to wrap, build MPIwrapper:
git clone https://github.com/eschnett/MPIwrapper
cmake -S MPIwrapper -B build-mpiwrapper-openmpi \
    -DMPIEXEC_EXECUTABLE=$(which mpiexec) \
    -DCMAKE_INSTALL_PREFIX=$HOME/mpiwrappers/openmpi
cmake --build build-mpiwrapper-openmpi -j
cmake --install build-mpiwrapper-openmpi

# 3) Build eOn with MPItrampoline on the cmake path:
export CMAKE_PREFIX_PATH=$HOME/mpitrampoline:$CMAKE_PREFIX_PATH
meson setup builddir -Dwith_mpi=true
meson compile -C builddir

# 4) At run time, tell MPItrampoline which wrapper to forward to:
export MPITRAMPOLINE_LIB=$HOME/mpiwrappers/openmpi/lib/libmpiwrapper.so
export MPITRAMPOLINE_MPIEXEC=$HOME/mpiwrappers/openmpi/bin/mpiwrapper-mpiexec
mpiexec -n 4 ./eonclient   # the trampoline-shipped mpiexec
```

The same `eonclient` binary then works against MPICH, OpenMPI,
Intel MPI, Spectrum MPI, or Cray MPICH by swapping the
`MPITRAMPOLINE_LIB` env var; no eOn rebuild needed.

If MPItrampoline is unavailable, the meson build falls back to a
direct `dependency('mpi')` link against the system MPI as before.

```{note}
This is only for modified VASP at the moment..
```

## VASP

If you have access to the modified VASP source code that is compatible with eOn,
you can compile a version of VASP that will work with eOn. This can be
accomplished by compiling with `make MPMD=1`.

It is necessary to set the environment variable `eOn_NUMBER_OF_CLIENTS` to be
the number of clients in the MPI job. If the client is being run directly,
instead of being run with the server, such as with the MPI communicator, then
the environment variable `eOn_CLIENT_STANDALONE` must be set to 1.

Each group of VASP ranks will write its output to a directory named `vasp###`
where the number that follows is a zero-padded number that ranges from zero to
`eOn_NUMBER_OF_CLIENTS`. There is a script in the tools directory named
`mkvasp.py` that takes the number of clients and a path and then creates these
directories in that path.

An `INCAR` file must be prepared that has the following lines in addition to any
other settings you wish to specify:

```{code-block} bash
IBRION=3
POTIM=0
LMPMD=.TRUE.
NSW=99999999
EDIFF=1E-7
EDIFFG=-1e-6
LWAVE=.FALSE.
LCHARG=.FALSE.
```


### Example

Here is an example of a script that will run 8 VASP ranks and 1 client rank::

```{code-block} bash
#!/bin/sh
mkdir vasp000
export eOn_CLIENT_STANDALONE=1
export eOn_NUMBER_OF_CLIENTS=1
mpirun -n 48 vasp_mpmd : -n 1 client_mpi
```
