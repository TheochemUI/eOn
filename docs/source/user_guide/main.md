# General Simulation Parameters

The following contains general options which specify the calculation to be done
and general parameters which are shared between job types.

## Jobs

More details are in subsequent sections, but a brief overview of the jobs are as
grouped below.

### aKMC

Adaptive Kinetic Monte Carlo (aKMC)
: uses saddle searches to find possible reactive events and KMC to determine the state-to-state kinetics.

Within such a job, we may also use

Recycling
: to use reactive events (processes) between states to speed up aKMC.

or construct a

Kinetic database
: of reactive events so they can be reused between AKMC runs.

Or even reach for

Coarse Graining
: to identify and escape groups of states connected by low barriers.

### Dynamics

Beyond this, we may also run

Molecular dyanmics
: which are standard.

We can make this a little more interesting with

Hyperdynamics
: which uses a bias potential to accelerate transitions between states.

## Collective sampling

Parallel Replica Dynamics
: which uses a set of trajectories to accelerate the rate of escape from states.

Basin hopping
: which lowers the barrier between states to accelerate Monte Carlo sampling.

## PES traversal

We can find minima on the potential energy surface using a

Minimization
: which optimizes the geometry of a structure.

We can go close to a minimum energy path with the

Nudged elastic band
: which locates minimum energy pathways using the eponymous method, essentially a chain-of-states search strategy

with this, we can then run a

Saddle search job
: to find a nearby saddle point.

and they may do so by using

Dimer
: method to find the lowest curvature mode.

or the

Lanczos
: method to find the lowest curvature mode.

```{versionadded} 2.1_TBA
with acceleration provided by Gaussian Process Regression or neural networks.
```

## Miscellaneous

Finally we note the following overview of some other sections 

Main
: which has options not specific to a single job type.

Communicator
: which detail methods by which the code can be run in parallel.

Potential
: which configures interatomic potentials, both bundled with `eON` and interfaced

Optimizer
: which has options related to the generic step interface for optimization of atomic structures

Prefactor
: which has options to control the calculation of harmonic transition state (hTST) prefactors

Structure Comparison
: which configures the computation of similarity measures and equivalence thresholds

Paths
: which enumerate directories in which calculations will be run

Debug
: which holds options for additional output

## Configuration

```{code-block} ini
[Main]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.MainConfig
```

