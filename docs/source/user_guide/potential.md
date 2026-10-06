---
myst:
  html_meta:
    "description": "A list and description of interatomic potentials supported by eOn, including vendored potentials and interfaces to external codes like VASP, LAMMPS, and ASE."
    "keywords": "eOn potential, interatomic potential, VASP, LAMMPS, ASE, EMT, EAM"
---

# Potential

`eOn` supports many potentials, some vendored within the executable and
libraries and others via interfaces.

```{note}
Some of these require compile-time flags, detailed in the [installation instructions](project:../install/index.md).
The `conda-forge` package (`conda install -c conda-forge eon`) includes
**Metatomic**, **XTB**, **EXT_POT**, and the vendored potentials.
The client always compiles the LAMMPS interface and loads `liblammps` with `dlopen`.
The option `-Dwith_lammps` does not exist.
ASE, VASP, AMS, and MPI potentials require building from source with
the corresponding `-Dwith_*` flags.
```

## Supported potentials

### External

RGPOT
: In-process rgpot calculators. Set `potential = rgpot` and configure `[RgpotPot]`. NWChem and CPMD are backends. See {doc}`rgpot_pot`. {bdg-success}`conda-forge`

VASP {cite:p}`pot-kresseEfficientIterativeSchemes1996`
: Vienna Ab-Initio Simulation Program (VASP) I/O interface. `-Dwith_vasp` defaults false. {bdg-warning}`source build`

LAMMPS {cite:p}`pot-plimptonFastParallelAlgorithms1995,pot-thompsonLAMMPSFlexibleSimulation2022`
: Library interface, detailed [documentation here](project:../user_guide/lammps_pot.md). `in.lammps` must be in the working directory, or the potential raises an error. {bdg-success}`conda-forge`

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
: Extended Tight binding models via native Fortran-C interfce {cite:p}`pot-bannwarthExtendedTightbindingQuantum2021`. `-Dwith_xtb` defaults false. Prefer `potential = rgpot` in `[Potential]` and `backend = xtb` in `[RgpotPot]`. {bdg-success}`conda-forge`
Metatomic
: Common interface to atomistic machine learning models. `-Dwith_metatomic` defaults false. {bdg-success}`conda-forge`
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

EMT
: Effective medium theory, for metals.

LJ {cite:p}`pot-jonesDeterminationMolecularFields1924`
: Lennard-Jones in reduced units, served by `rgpot`. Neighbor pairs via
  [vesin](neighbor_lists.md); timed by ASV
  `TimeMinimizationLJCluster` (ljcluster).

LJCluster {cite:p}`pot-jonesDeterminationMolecularFields1924`
: Lennard-Jones cluster variant, served by `rgpot`.

Morse_Pt
: Hard sphere morse potential for Platinum, served by `rgpot`. Neighbor
  pairs via [vesin](neighbor_lists.md); timed by ASV `TimePointMorsePt` /
  saddle / NEB Morse fixtures.

ZBL
: Ziegler-Biersack-Littmark screened nuclear repulsion, served by `rgpot`.

```{versionadded} 3.2.1
DFTD3 / DFTD4
: Grimme DFT-D via rgpot 3.2 (`potential = dftd3` / `dftd4`,
  `[D3Pot]` / `[D4Pot]`). The option `-Dwith_dftd3` does not exist.
  The client compiles these pots when `RGPOT_HAS_DFTD3` or `RGPOT_HAS_DFTD4` is defined.

EXPR
: rgpot ExprPot. `potential = expr` with `[ExprPot]` `expression` and
  comma-separated `terms` (`0.5*lj + d3`). Terms: `lj`, `ljcluster`,
  `morse`, `zbl`, `d3`/`dftd3`, `d4`/`dftd4`, `mopac`.

MOPAC
: rgpot 3.2 MOPACPot (libmopacc). `potential = mopac`, `[MOPACPot]`
  `charge`, `spin`, `model` (4 is AM1), `engine_path`.
```

Lenosky_Si {cite:p}`pot-lenoskyHighlyOptimizedEmpirical2000`
: Lenosky potential, for silicon.

SW_SI {cite:p}`pot-stillingerComputerSimulationLocal1985`
: Stillinger-Weber potential, for silicon.

Tersoff_SI {cite:p}`pot-tersoffEmpiricalInteratomicPotential1988`
: Tersoff pair potential with angular terms, for silicon.

EDIP {cite:p}`pot-justoInteratomicPotentialSilicon1998`
: Environment-Dependent Interatomic Potential, for carbon.

TIP4P {cite:p}`pot-jorgensenComparisonSimplePotential1983`
: Point charge model for water, also for water-hydrogen and water on platinum. Source builds need `-Dwith_water=true` because `with_water` defaults false.

SPCE {cite:p}`pot-berendsenMissingTermEffective1987`
: Extended simple point charge model for water. Source builds need `-Dwith_water=true` because `with_water` defaults false.

## Configuration

A known potential name matches in any case.
`RGPOT` is stored as `rgpot`.
An exact listed spelling wins, so `SocketNWChem` stays distinct from `socketnwchem`.

```{code-block} ini
[Potential]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.PotentialConfig
```

## Potential configurations

Several potentials have additional configuration stanzas.

### RGPOT

```{eval-rst}
.. autopydantic_model:: eon.schema.RgpotPot
```

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

```{warning}
NWChem's Fortran i-PI socket driver truncates UNIX socket names to approximately
30 characters. The full socket path is ``/tmp/ipi_<unix_socket_path>``, so
``unix_socket_path`` should be kept short (under ~20 characters). For example,
``eon_nwchem`` works but ``eon_nwchem_test_socket`` is truncated,
and the connection fails with no clear error message.
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
