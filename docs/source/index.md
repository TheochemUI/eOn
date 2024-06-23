# eON: Long Timescale Dynamics Software

The EON software package contains a set of algorithms used primarily to model
the evolution of atomic scale systems over long time scales. Standard molecular
dynamics algorithms, based upon solving Newton's equations, are limited by the
femtosecond time scale of atomic vibrations. EON simulations are designed for
rare event systems where the interesting dynamics can be described by fast
transition between stable states. In each algorithm, the residence time in the
stable states is modeled with statistical mechanics, and the important
state-to-state dynamics are modeled stochastically.

The algorithms currently implemented include parallel replica dynamics,
hyperdyamics, adaptive kinetic Monte Carlo, and basin hopping.

## Supported systems

eON can handle both molecular systems (e.g. gas-phase reactions) and extended
(surface) systems, with robust periodic image boundary support.

```{figure} fig/esys_trans.png
---
alt: Collection of systems which can be modeled
class: full-width
align: center
---
An overview of some systems modelled with eON
```

However, the systems which are best modelled using EON are those in which the
important kinetics are governed by rare events. Diffusion in solids and chemical
reactions at surface are particularly suitable when there is a clear separation
of time scales between atomic vibrations at the diffusion or catalytic events of
interest.


```{figure} fig/alripe.png
---
alt: Ripening dynamics
class: full-width
align: center
---
Al(100) ripening dynamics
```

In the example showing ripening dynamics on an `Al(100)` surface, a compact
island forms after `65720` transitions in a time scale of a `ms` at `300K`.

# Interatomic Interactions

There are a variety of empirical potentials included with `eON`. You can also
use the potentials built into the LAMMPS library. `eON` also provides an
interace to the VASP and GPAW density functional theory codes.

```{versionadded} 2.0
`eON` now supports additional potentials
* via an embeded Python interpreter, all the potentials accessible from the atomic simulation environment, ASE.
* via native Fortran-C interface, different forms of the tight binding `XTB` potentials
* via an I/O and server-client interface, potentials from the Amsterdam Modeling Suite (AMS)
```

# Parallel Interfaces

The algorithms in EON can be run in parallel using a set of communication
options including local communication, distribution over a cluster using a
queueing system, multiple process multiple data MPI jobs, and distributed
computing environments.

<!-- Is this still true -->

# Additional Topics

```{toctree}
:maxdepth: 2
:caption: Contents

team
install/index
tutorials/index
user_guide/index
devdocs/index
apidocs/index
```

## Indices and tables

- [](genindex)
- [](modindex)
- [](search)
