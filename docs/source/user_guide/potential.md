---
myst:
  html_meta:
    "description": "A list and description of interatomic potentials supported by eOn, including vendored potentials and interfaces to external codes like VASP, LAMMPS, and ASE."
    "keywords": "eOn potential, interatomic potential, VASP, LAMMPS, ASE, EMT, EAM"
---

# Potential

`eOn` supports a large number of potentials, some vendored within the executable
and libraries and others via interfaces.

```{note}
Some of these require compile-time flags, detailed in the [installation instructions](project:../install/index.md).
The `conda-forge` package (`conda install -c conda-forge eon`) ships with
**Metatomic**, **XTB**, **EXT_POT**, and the vendored potentials.
LAMMPS, ASE, VASP, AMS, and MPI potentials require building from source with
the corresponding `-Dwith_*` flags.
```

## Supported Potentials

### External

VASP {cite:p}`pot-kresseEfficientIterativeSchemes1996`
: Vienna Ab-Initio Simulation Program (VASP) I/O interface. {bdg-warning}`source build`

LAMMPS {cite:p}`pot-plimptonFastParallelAlgorithms1995,pot-thompsonLAMMPSFlexibleSimulation2022`
: Library interface, detailed [documentation here](project:../user_guide/lammps_pot.md). {bdg-warning}`source build`

EXT_POT
: File-based interface to any external calculator. Detailed [documentation here](project:ext_pot.md). {bdg-success}`conda-forge`

```{versionadded} 2.0
AMS(-IO)
: Amsterdam modeling suite {cite:p}`pot-teveldeChemistryADF2001`, both I/O and library. {bdg-warning}`source build`
ASE_ORCA
: Atomic simulation environment {cite:p}`pot-larsenAtomicSimulationEnvironment2017` interface to ORCA {cite:p}`pot-neeseORCAQuantumChemistry2020`. {bdg-warning}`source build`
ASE_NWChem
: Atomic simulation environment {cite:p}`pot-larsenAtomicSimulationEnvironment2017` interface to NWChem {cite:p}`pot-apraNWChemPresentFuture2020`. {bdg-warning}`source build`
XTB
: Extended Tight binding models via native Fortran-C interfce {cite:p}`pot-bannwarthExtendedTightbindingQuantum2021`. {bdg-success}`conda-forge`
Metatomic
: Common interface to atomistic machine learning models. {bdg-success}`conda-forge`
SocketNWChem
: Socket oriented communicator for efficient integration with NWChem {cite:p}`pot-apraNWChemPresentFuture2020`. {bdg-success}`conda-forge`
```

### Vendored

CuH2
: Copper Hydride system

FeHe
: Iron-hydrides

EAM_Al
: Embedded atom method parameterized for Aluminum.

QSC {cite:p}`pot-kimuraQuantumSuttonChenManyBody1998`
: Quantum Sutton-Chen potential, for FCC metals.

EMT
: Effective medium theory, for metals.

LJ {cite:p}`pot-jonesDeterminationMolecularFields1924`
: Lennard-Jones in reduced units

Morse_Pt
: Hard sphere morse potential for Platinum

Lenosky_Si {cite:p}`pot-lenoskyHighlyOptimizedEmpirical2000`
: Lenosky potential, for silicon.

SW_SI {cite:p}`pot-stillingerComputerSimulationLocal1985`
: Stillinger-Weber potential, for silicon.

Tersoff_SI {cite:p}`pot-tersoffEmpiricalInteratomicPotential1988`
: Tersoff pair potential with angular terms, for silicon.

EDIP {cite:p}`pot-justoInteratomicPotentialSilicon1998`
: Environment-Dependent Interatomic Potential, for carbon.

TIP4P {cite:p}`pot-jorgensenComparisonSimplePotential1983`
: Point charge model for water, also for water-hydrogen and water on platinum.

SPCE {cite:p}`pot-berendsenMissingTermEffective1987`
: Extended simple point charge model for water

```{deprecated} 2.0
These potentials are missing in the SVN sources..
bopfox
: Bond order potential, for metals
```

## Configuration

```{code-block} toml
[Potential]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.PotentialConfig
```

## Potential configurations

Several potentials have additional configuration stanzas.

### Metatomic

```{eval-rst}
.. autopydantic_model:: eon.schema.Metatomic
```

### XTB

```{eval-rst}
.. autopydantic_model:: eon.schema.XTBPot
```

### ZBL

```{eval-rst}
.. autopydantic_model:: eon.schema.ZBLPot
```

### NWChem

Support for `nwchem` works best with the socket potential structure as noted in
the [reproduction details](https://github.com/theochemUI/otgpd_repro) of the
Optimal transport Gaussian Process
{cite:p}`pot-goswamiAdaptivePruningIncreased2025b`, and can lead to manyfold
increases in speed compared to file or ASE interfaces
{cite:p}`pot-goswamiEfficientExplorationChemical2025`.

```{eval-rst}
.. autopydantic_model:: eon.schema.SocketNWChemPot
```

An older ASE interface exists as well.

### ASE potentials

There are several specific ASE potentials supported,

```{eval-rst}
.. autopydantic_model:: eon.schema.ASE_NWCHEM
```

```{eval-rst}
.. autopydantic_model:: eon.schema.ASE_ORCA
```

### AMS potentials

Both a direct server model and a file based integration exist.

```{eval-rst}
.. autopydantic_model:: eon.schema.AMSConfig
```

```{eval-rst}
.. autopydantic_model:: eon.schema.AMSIOConfig
```

Along with helpers to set environment variables for these calculations.

```{eval-rst}
.. autopydantic_model:: eon.schema.AMSEnvConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: POT_
keyprefix: pot-
---
```
