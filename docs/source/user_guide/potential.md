# Potential

`eON` supports a large number of potentials, some vendored within the executable
and libraries and others via interfaces.

```{note}
Some of these require flags to be set, details are in the [installation instructions](project:../install/index.md)
```

## Supported Potentials

### External

VASP {cite:p}`pot-kresseEfficientIterativeSchemes1996`
: Vienna Ab-Initio Simulation Program (VASP) I/O interface.

LAMMPS {cite:p}`pot-plimptonFastParallelAlgorithms1995,pot-thompsonLAMMPSFlexibleSimulation2022`
: Library interface, detailed [documentation here](project:../user_guide/lammps_pot.md)

EXT_POT
: Writes box size and coordinates to the file `from_eon_to_ext` and makes a system call to `ext_pot` which must populate `from_ext_to_eon`

```{versionadded} 2.0
AMS(-IO)
: Amsterdam modeling suite {cite:p}`pot-teveldeChemistryADF2001`, both I/O and library
ASE_ORCA
: Atomic simulation environment {cite:p}`pot-larsenAtomicSimulationEnvironment2017` interface to ORCA {cite:p}`pot-neeseORCAQuantumChemistry2020`
ASE_NWChem
: Atomic simulation environment {cite:p}`pot-larsenAtomicSimulationEnvironment2017` interface to NWChem {cite:p}`pot-apraNWChemPresentFuture2020`
XTB
: Extended Tight binding models via native Fortran-C interfce {cite:p}`pot-bannwarthExtendedTightbindingQuantum2021`
Metatomic
: Common interface to atomistic machine learning models
SocketNWChem
: Socket oriented communicator for efficient integration with NWChem {cite:p}`pot-apraNWChemPresentFuture2020`
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

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: POT_
keyprefix: pot-
---
```
