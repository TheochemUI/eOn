# Communicator

`eON` has a server client architecture for running its calculations. The
simulation data is stored on the server and clients are sent jobs and return the
results. Each time `eON` is run it first checks to see if any results have come
back from clients and processes them accordingly and then submits more jobs if
needed. In `eON` there are several different ways to run jobs. One can run them
locally on the server, via MPI,  or using a job queuing system such as
[SGE](http://www.oracle.com/us/products/tools/oracle-grid-engine-075549.html).

```{note}
As of version 2.0 onwards, we recommend using dedicated workflow management tools (like [AiiDA](https://www.aiida.net/) or [Snakemake](https://snakemake.readthedocs.io/) or [Fireworks](https://materialsproject.github.io/fireworks)) instead of using `eON` to generate submission scripts.
```

## Configuration

```{code-block} ini
[Communicator]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.CommunicatorConfig
```

### Examples

An example communicator section using the local communicator with an Eon client
binary named `eonclient-custom` that either exists in the `$PATH` or in the same
directory as the configuration file and uses makes use of 8 CPUs.


```{code-block} toml
[Communicator]
type = "local"
client_path = "eonclient-custom"
number_of_cpus = 8
```

## Additional Topics

```{versionchanged} 2.0
Potentials which can be run in parallel, like those accessed through ASE (e.g. ORCA) are always run in parallel, for the others, there is little to no benefit for this additional overhead.
```

### MPI

```{warning}
Not tested on 2.0, only ever supported AKMC.
```

The MPI communicator allows for the server and client to be run as a MPI job.
The number of clients that are run and thus the number of jobs is set at runtime
by the MPI environment.

A MPI aware client must be compiled, which will be named ``eonclientmpi``
instead of ``eonclient``. It can only be used to run MPI jobs.

To run EON with MPI, two environment variables must be set. The
variable `EON_NUMBER_OF_CLIENTS` determines how many of the ranks
should become clients and `EON_SERVER_PATH` is the path to the
server Python script. In MPI mode the clients need to be started
instead of the server and one of them will become the server process.
Currently only AKMC is supported. Below is an example of running using
the MPI communicator:

```{code-block} bash
#!/bin/bash
export EON_NUMBER_OF_CLIENTS=7
export EON_SERVER_PATH=~/eon/akmc.py
mpirun -n 8 ~/eon/client/eonclientmpi
```

### Cluster

```{warning}
Not tested on 2.0
```

An example communicator section for the cluster communicator using the provided
`sge6.2` scripts and a name prefix of `al_diffusion_`:

```{code-block} toml
[Communicator]
type = "cluster"
name_prefix = "al_diffusion_"
script_path = "/home/user/eon/tools/clusters/sge6.2"
```
