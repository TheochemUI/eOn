---
myst:
  html_meta:
    "description": "Glossary of method names and short forms used in the eOn documentation."
    "keywords": "eOn glossary, NEB, AKMC, minimizer, potentials"
---

# Glossary

Short forms used in this book. Each entry points at a page that uses the term.

```{glossary}
:sorted:

AKMC
  Adaptive kinetic Monte Carlo. Saddle searches supply the rare events,
  and kinetic Monte Carlo advances the state.
  See {doc}`user_guide/akmc`.

AMS
  Amsterdam Modeling Suite. The potential page documents the I/O and
  library interfaces.
  See {doc}`user_guide/potential`.

BFGS
  Broyden-Fletcher-Goldfarb-Shanno. xtsci-optimize supplies the dense
  inverse-BFGS. LBFGS is the limited-memory form.
  See {doc}`user_guide/minimization`.

CNA
  Common-neighbor labels. 0 is fcc (421), 1 is hcp (422), and 2 is
  other.
  See {doc}`releases/v3.3.0/index`.

CON
  Coordinate file extension. Jobs read frames such as `reactant.con` and
  write paths such as `neb.con`.
  See {doc}`devdocs/in-process`.

COSMO
Lab-COSMO
  Prefix on the metatomic model id `lab-cosmo/pet-mad`.
  See {doc}`user_guide/metatomic_pot`.

CPMD
  Car-Parrinello molecular dynamics. With `ranks_per_image`, one
  instanton force call spreads beads across CPMD calculator groups.
  See {doc}`user_guide/instanton`.

DFT
  Density functional theory. The introduction names VASP among these
  codes.
  See {doc}`index`.

EAM
  Embedded atom method. The vendored aluminum potential is `EAM_Al`.
  See {doc}`user_guide/potential`.

FIRE
  Fast inertial relaxation engine. `opt_method = fire` selects it, and
  the schema also offers a FIRE 2.0 variant.
  See {doc}`user_guide/optimizer`.

FPE
  Floating-point exception. v2.12 hardens the signal handler and
  suspends trapping around libtorch calls.
  See {doc}`releases/v2.12.0/index`.

GHA
  GitHub Actions. The release workflow and the Windows build run there.
  See {doc}`devdocs/release`.

GPG
  GNU Privacy Guard. Maintainers sign release tags with `git tag -s`.
  See {doc}`devdocs/release`.

GPR
  Gaussian process dimer. `min_mode_method = gprdimer` runs it when that
  checkout is present.
  See {doc}`user_guide/saddle_search`.

IDPP
  Image dependent pair potential. It builds an initial band where linear
  interpolation would overlap atoms.
  See {doc}`user_guide/neb`.

INI
  Job configuration syntax. Bracketed sections hold the keys for a run.
  See {doc}`user_guide/main`.

IRA
  Iterative Rotations and Assignments. `match_method = ira` rigid-rotates
  and permutes one NEB endpoint onto the other.
  See {doc}`user_guide/neb`.

JOTA
  Journal of Optimization Theory and Applications. The minimizer page
  cites the 1999 and 2001 volumes for the Zhang-Deng-Chen and Zhang-Xu
  secant.
  See {doc}`user_guide/minimization`.

LBFGS
  Limited-memory Broyden-Fletcher-Goldfarb-Shanno. The schema names this
  minimizer `lbfgs`.
  See {doc}`user_guide/minimization`.

MLIP
  Machine-learning potential. Metatomic loads these models.
  See {doc}`user_guide/metatomic_pot`.

MPI
  Message Passing Interface. Source builds pass `-Dwith_mpi=enabled` for
  the MPI potential.
  See {doc}`user_guide/mpi_potential`.

NEB
  Nudged elastic band. Images along a path between a reactant and a
  product relax under spring forces.
  See {doc}`user_guide/neb`.

NLCG
  Nonlinear conjugate gradient. xtsci-optimize keeps the conjugacy across
  outer steps.
  See {doc}`user_guide/minimization`.

OIDC
  OpenID Connect. The PyPI publish step authenticates with it.
  See {doc}`devdocs/release`.

ORCA
  ORCA quantum chemistry program. eOn has a direct interface.
  See {doc}`user_guide/potential`.

PBC
  Periodic boundary conditions. Molecular calculators can refuse a
  periodic cell, and Matter can wrap positions.
  See {doc}`user_guide/artn`.

PDM
  Tracks versions for the documentation build.
  See {doc}`devdocs/docbuild`.

PET-MAD
  Metatomic model published as `lab-cosmo/pet-mad`. The name is one model.
  See {doc}`user_guide/metatomic_pot`.

QTST
  Path-integral quantum transition-state theory. The instanton page
  writes the rate as PI-QTST.
  See {doc}`user_guide/instanton`.

RFO
  Banerjee/Baker rational function optimization. `[Xtsci] method = rfo`
  selects it.
  See {doc}`user_guide/minimization`.

RMS
  Root-mean-square force. `convergence_metric = rms` uses
  `||F||_2 / sqrt(3 N_free)`. The default metric is the plain Euclidean
  norm.
  See {doc}`user_guide/minimization`.

RPC
  Cap'n Proto remote procedure call. Serve mode exposes a potential
  through the rgpot protocol.
  See {doc}`user_guide/serve_mode`.

RPMD
  Ring-polymer molecular dynamics. A transmission coefficient turns the
  PI-QTST rate into this rate.
  See {doc}`user_guide/instanton`.

RTDB
  NWChem runtime database. The direct calculator skips a reset when the
  method parameters stay the same.
  See {doc}`user_guide/rgpot_integration`.

SCF
  Self-consistent field. NWChem keeps this state from one force call to
  the next.
  See {doc}`user_guide/rgpot_integration`.

SIDPP
S-IDPP
  Sequential IDPP. The path grows inward from the reactant and from the
  product.
  See {doc}`user_guide/neb`.

SIAM
  Society for Industrial and Applied Mathematics. The minimizer page
  cites SIAM J. Numer. Anal. 1986 for the Grippo-Lampariello-Lucidi
  window.
  See {doc}`user_guide/minimization`.

SIMD
  Single instruction, multiple data. Highway SIMD is an optional
  subproject.
  See {doc}`releases/v2.13.0/index`.

SOFI
  Point-group detection in `IRACompare::findSymmetry`.
  See {doc}`releases/v2.13.0/index`.

VASP
  Vienna Ab-Initio Simulation Program. eOn drives it through a file
  interface.
  See {doc}`user_guide/potential`.

XTB
  Extended tight-binding models, linked through the native Fortran
  interface.
  See {doc}`user_guide/potential`.

ZBL
  Ziegler-Biersack-Littmark screened nuclear repulsion. `sidpp_zbl` wraps
  an IDPP path with it.
  See {doc}`user_guide/potential`.

```
