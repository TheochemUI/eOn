# User Guide

Here we collect a brief introduction to each algorithm, with curated references
for more information along with the  configuration settings as implemented
within `eON`.

Each of the sections and methods may be included in the configuration file with the appropriate section header.

## Overview

The `eON` program can run the following methods or tasks to explore the
configuration space of molecular systems and to accelerate the simulation of
their dynamics over long times. To run `eON`, a configuration file must be generated
using the options specified in the documentation. Each section header is denoted
by square brackets, and is followed by the key/value pairs. For example, to set
the *job_type* key of the *[Main]* section to the value *process_search*, your
configuration file would include the lines:

```{code-block} toml
[Main]
job_type = "process_search"
```

```{versionchanged} 2.1_TBA
The format for the configuration files were changed from INI [to TOML](https://toml.io/en/), where this causes relevant user changes it will be highlighted.
- The default file is now `config.toml` instead of `config.ini`
```

There are specific options for each method, and a set of general options which
are shared between methods. Examples of these general options include
specificion of the interatomic potential and the parameters for doing structural
optimizations and comparisons.

````{margin}
```{note}
As of version 2.0 onwards, we recommend using dedicated workflow management tools (like [AiiDA](https://www.aiida.net/) or [Snakemake](https://snakemake.readthedocs.io/) or [Fireworks](https://materialsproject.github.io/fireworks)) instead of using `eON` to generate submission scripts.
```
````

EON is designed to run in serial on one computer or in parallel using a
communicator to send jobs from a server to clients and receive the results back.
Several parallelization schemes have been implemented, including local
communication through files, cluster communication via a queuing system or via
mpi. 

Methods run in parallel are broken up by the eon server into tasks which
are run by client program. The server then compiles the information sent back
by the clients in a way that can be used by the sampling or dynamics methods.

## Configurationa Sections

```{toctree}
:maxdepth: 1
:caption: Common

main
potential
structure_comparison
optimizer
```


```{toctree}
:maxdepth: 1
:caption: Potential Energy Surface Traversal

minimization
neb
dimer
lanczos
hessian
prefactor
```

```{toctree}
:maxdepth: 1
:caption: Compound Jobs

akmc
saddle_search
basin_hopping
process_search
recycling
coarse_graining
```

```{toctree}
:maxdepth: 1
:caption: Dynamics

dynamics
parallel_replica
hyperdynamics
```

```{toctree}
:maxdepth: 1
:caption: Additional configuration

communicator
debug
paths
```

```{toctree}
:maxdepth: 1
:caption: External linkage

lammps_pot
kdb
mpi_potential
```
<!-- TODO(rg) Maybe add the whole config as well -->
