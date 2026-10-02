---
myst:
  html_meta:
    "description": "Guide to eOn's server-client architecture and communication options for running parallel calculations locally, via MPI, or on a cluster."
    "keywords": "eOn communicator, parallel, MPI, cluster, server-client"
---

# Communicator

`eOn` has a server client architecture for running its calculations. The
simulation data is stored on the server and clients are sent jobs and return the
results. Each time `eOn` is run it first checks to see if any results have come
back from clients and processes them, then submits more jobs if
needed. In `eOn` there are several different ways to run jobs. One can run them
locally on the server, via MPI,  or using a job queuing system such as
[SGE](http://www.oracle.com/us/products/tools/oracle-grid-engine-075549.html).

```{note}
From 2.0 on, prefer a workflow manager over eOn-generated submit
scripts. For AiiDA the plugin is {doc}`aiida` (`pip install aiida-eon`).
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


```{code-block} ini
[Communicator]
type = "local"
client_path = "eonclient-custom"
number_of_cpus = 8
```

### In-process (`type = local_lib` / inprocess)

`LocalInProcess` runs jobs as `pyeonclient.Matter` in the server process.
There is no `eonclient` subprocess. It needs a build with
`-Dwith_pyeonclient=true`.

```{code-block} ini
[Communicator]
type = "inprocess"
```

Dispatch follows `Parameters.job`:

| Job | What runs |
|---|---|
| minimization (default) | `Matter.relax` |
| point | energy / forces only |
| process_search / saddle_search | `ProcessSearch` on the Matter |

A job carries geometry as a `Structure` or a `readcon.ConFrame`
(`structure`, `conframe`, `reactant`, `pos`, or `pos.con`). `.con` text is
refused. The result dict returns `product` and, when a saddle exists,
`saddle` as `ConFrame` objects. Those frames are not written to disk.
`_matter` stays the potential-bearing client object and `_structure` the
numpy working set. `results.dat` is still synthesized as text so the classic
explorer can parse scalars.

The same record is the dict `job_result`. Its fields are
`termination_reason` (the status integer), `termination_reason_text`
(`GOOD`, `FAIL`, or `cancelled`), `job_type`, `potential_energy`, and
`total_force_calls`.

`cancel_state` returns 1 only while `submit_jobs` is inside a batch. That
call sets a token the next job sees, and a second token that a compiled
relax, band, or saddle search polls during the current call. A call while
no batch is running returns 0. The token is cleared when `submit_jobs`
returns.

## Additional topics

```{versionchanged} 2.0
Potentials which can be run in parallel, like those accessed through ASE (e.g. ORCA) are always run in parallel, for the others, there is little to no benefit for this additional overhead.
```

### MPI

```{note}
Open MPI 5 on 3.4 ran an adaptive kinetic Monte Carlo check with one
server rank and two client ranks.
```

The MPI communicator runs the server and the clients as one MPI job. The
number of clients, and so the number of jobs in flight, comes from the MPI
launch.

Build the client with `-Dwith_mpi=enabled`; the resulting `eonclient` only runs
under MPI. Two environment variables set the layout. `EON_NUMBER_OF_CLIENTS`
is how many ranks become clients, and `EON_SERVER_PATH` is a Python script
that starts the server. Launch the clients, not the server: one extra rank
turns into the server and runs that script. With `[Communicator] type = mpi`,
adaptive kinetic Monte Carlo, parallel replica, and basin hopping wait on
that communicator. `examples/prd_mpi` and `examples/bh_mpi` use this launch.

```{code-block} python
# server.py
import eon.server
eon.server.main()
```

```{code-block} bash
#!/bin/bash
export EON_NUMBER_OF_CLIENTS=7
export EON_SERVER_PATH=$PWD/server.py
mpirun -n 8 /path/to/eonclient
```

A client whose job fails, through a potential error, a `config.ini` it cannot
load or a job type it does not know, logs the error into that job's directory
and hands the directory back without `results.dat`. The server logs
`returned no results` for it and the rank takes the next job.

### Cluster

```{warning}
Not tested on 2.0
```

An example communicator section for the cluster communicator using the provided
`sge` scripts and a name prefix of `al_diffusion_`:

```{code-block} ini
[Communicator]
type = "cluster"
name_prefix = "al_diffusion_"
script_path = "/home/user/eon/tools/clusters/sge"
```

The Slurm scripts in `tools/clusters/slurm` take their site options from
the environment of the server process. `EON_SBATCH_ARGS` sets extra
`sbatch` options (account, partition, nodes, tasks, time) and `EON_CLIENT`
the command each job runs, `eonclient` by default:

```{code-block} bash
export EON_SBATCH_ARGS="-A myaccount -p cpu -N 1 --ntasks=48 -t 01:00:00"
export EON_CLIENT=eonclient
eon-server
```

```{code-block} ini
[Communicator]
type = "cluster"
script_path = "/path/to/eon/tools/clusters/slurm"
```

A potential that needs several MPI ranks, such as CPMD behind `ext_pot`,
starts them from its wrapper inside the job (`srun -n 48 cpmd.x ...`); the
client itself stays one process. When a state reaches its confidence the
server cancels only the jobs still queued. Jobs that already finished are
harvested on the next pass and filtered by state.
