# MPI Potential

```{note}
This is only for modified VASP at the moment..
```

## VASP

If you have access to the modified VASP source code that is compatible with EON,
you can compile a version of VASP that will work with EON. This can be
accomplished by compiling with `make MPMD=1`.

It is necessary to set the environment variable `EON_NUMBER_OF_CLIENTS` to be
the number of clients in the MPI job. If the client is being run directly,
instead of being run with the server, such as with the MPI communicator, then
the environment variable `EON_CLIENT_STANDALONE` must be set to 1.

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


### Example

Here is an example of a script that will run 8 VASP ranks and 1 client rank::

```{code-block} bash
#!/bin/sh
mkdir vasp000
export EON_CLIENT_STANDALONE=1
export EON_NUMBER_OF_CLIENTS=1
mpirun -n 48 vasp_mpmd : -n 1 client_mpi
```
