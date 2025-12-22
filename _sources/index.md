---
myst:
  html_meta:
    "description": "Official documentation for the EON software package, a tool for modeling long-timescale dynamics in atomic systems with methods like aKMC, NEB, and Parallel Replica Dynamics."
    "keywords": "EON, long-timescale dynamics, aKMC, NEB, Parallel Replica Dynamics"
---

# eOn: Long Timescale Dynamics Software

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

eOn can handle both molecular systems (e.g. gas-phase reactions) and extended
(surface) systems, with robust periodic image boundary support.

```{figure} fig/esys_trans.png
---
alt: Collection of systems which can be modeled
class: full-width
align: center
---
An overview of some systems modeled with eOn
```

However, the systems which are best modeled using EON are those in which the
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

## Interatomic Interactions

There are a variety of empirical potentials included with `eOn`. You can also
use the potentials built into the LAMMPS library. `eOn` also provides an
interace to the VASP and GPAW density functional theory codes.

```{versionadded} 2.0
`eOn` now supports additional potentials
* via an embeded Python interpreter, all the potentials in [ASE](https://ase-lib.org/).
* via native Fortran-C interface, different forms of the tight binding `XTB` potentials
* via native interface, the [Metatomic](https://docs.metatensor.org/metatomic/) potentials
* via an I/O and server-client interface, Amsterdam Modeling Suite (AMS) potentials
```

# Getting started

See [the installation instructions](https://eondocs.org/install/), but in a line:

```{code-block} bash
micromamba install -c conda-forge eon
# single point calculation, Lennard-Jones
eonclient -s molecule.con -p lj
# or with a config.ini and pos.con file
eonclient # reads config.ini and runs
# or for akmc, needs config.ini and pos.con
python -m eon.server
```

## Getting help

We support a variety of methods to provide assistance:

- **Github Issues** :: For bug reports and software errors, [open issues](https://github.com/TheochemUI/eOn/issues)
- **Community Forum** :: EON has a section on the [Materials Science Community Discourse](https://matsci.org/c/eon/)

## Supporting packages

Additional visualization and parsing may be found in the `rgpycrumbs` diagnostic
suite ([Home](https://rgpycrumbs.rgoswami.me/tools/eon/index.html),
[Github](https://github.com/HaoZeke/rgpycrumbs),
[PyPI](https://pypi.org/project/rgpycrumbs/)).

# User Guide

```{toctree}
:maxdepth: 2
:caption: Contents

team
install/index
tutorials/index
user_guide/index
devdocs/index
apidocs/index
releases/index
```

## Indices and tables

- [](genindex)
- [](modindex)
- [](search)
