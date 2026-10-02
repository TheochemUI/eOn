---
myst:
  html_meta:
    "description": "Instructions for using a modified version of VASP with eOn via the MPI Potential interface for parallel calculations."
    "keywords": "eOn MPI potential, VASP interface, parallel VASP, ab-initio"
---

# MPI potential

```{admonition} conda-forge availability
:class: warning
The `conda-forge` package omits MPI. Build from source with `-Dwith_mpi=enabled`.
```

```{note}
This is only for modified VASP at the moment..
```

## VASP

If you have access to the modified VASP source code that is compatible with eOn,
you can compile a version of VASP that will work with eOn. This can be
accomplished by compiling with `make MPMD=1`.

That VASP binary joins the eOn client in one MPI job. A standalone client
sets `EON_CLIENT_STANDALONE`. The client count is then 1, so the run does
not read `EON_NUMBER_OF_CLIENTS`. The client loads `config_0.ini` when that
file exists, and `config.ini` otherwise. The MPI communicator leaves
`EON_CLIENT_STANDALONE` unset and sets `EON_NUMBER_OF_CLIENTS` and
`EON_SERVER_PATH`, as on the [communicator](project:communicator.md) page.
Each client still gets the VASP output directory described next.

Each group of VASP ranks will write its output to a directory named `vasp###`
where the number that follows is a zero-padded number that ranges from zero to
`EON_NUMBER_OF_CLIENTS`. There is a script in the tools directory named
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


`[Potential] potential = mpi` selects this interface. A build without
`-Dwith_mpi=enabled` has no such potential, and the client stops with
`No known potential could be constructed`. The number of `vasp_mpmd` ranks
must be a positive multiple of the client count. `[Potential] mpi_poll_period`
is the wait between probes, in seconds. The default is 0.25.

### Example

Eight VASP ranks and one client:

```{code-block} bash
#!/bin/sh
mkdir vasp000
export EON_CLIENT_STANDALONE=1
mpirun -n 8 vasp_mpmd : -n 1 eonclient
```

```{code-block} ini
[Main]
job = point

[Potential]
potential = mpi
```

`job = point` reads `pos.con` from that working directory.
