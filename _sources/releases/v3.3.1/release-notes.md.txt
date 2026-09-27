# eOn 3.3.1

A patch on `v3.3.0`. See monorepo `CHANGELOG.md` section 3.3.1 for the
full list: 92 fixes, 47 changes, 44 removals and 4 additions.

## Build

The `v3.3.0` tag does not build eonclib: plugin sources used
`Parameters` with nothing but a forward declaration in scope. They now
include `Parameters.h`. Potential constructors that take `Parameters`
move to eoncbase, so plugin `.so` files do not link eonclib. The MPI
potential links against MPI.
Wheels pin nanobind to `>=2.2,<3`, and wheel repair finds
`librgpot.so*` next to `libreadcon_core`. The vendored copy goes in as
`libeon_rgpot.so.3`, so the pip `rgpot` wheel imports beside
pyeonclient in either order.

## Nudged elastic band (NEB)

- `minimize_endpoints` defaults to `false`, matching SVN and the
  pyeonclient NEB API.
- Fixed atoms carry no projected force. When the endpoints placed a
  fixed atom differently, the spring force leaked onto it and the norm
  convergence metric counted it.
- Reactant and product with different atom counts are rejected before
  the path interpolation.
- The default L-BFGS takes the two-loop step, and the NEB L-BFGS no
  longer runs the `auto_scale` finite-difference probe.
- Energy-weighted springs set `E_ref` to the higher endpoint energy.
- OCI-NEB restores a climbing image that ends a min-mode walk at
  positive curvature. `ci_mmf_restore_unhelpful` extends the restore to
  alignment rejects and force increases.
- Non-periodic images pass a zero cell into `forceBatch`.

The NEB regression tests run between two LJ13 minima joined by a
saddle, so the climbing-image case converges at `-O2` with `NDEBUG` as
well as at `-O3`.

## Basin hopping, saddle searches and dynamics

Basin hopping applies the Metropolis factor to uphill hops at a
positive temperature and rejects them at zero or below. It keeps the minimized energy as the reference
and takes the dimer direction from the minimum image. The classic dimer
keeps the orientation its torque check accepted. ARTn keeps its library
mutex locked from create to destroy and passes the cell to pARTn
untransposed.
Parallel replica keeps its bond boost across trajectory copies. Langevin
and Nose-Hoover dynamics count and move each free Cartesian component
on its own, so a partly fixed atom stays fixed in its frozen direction.

## Optimizers

`opt_method = xtsci` selects the xtsci-optimize engine, and
`[Xtsci] method` picks the solver inside it: L-BFGS, BFGS, SR1/SR2,
Newton, RFO, NLCG, Adam or PSO. One session persists across the host
loop, so L-BFGS pairs and NLCG conjugacy survive between steps.

## Runtime and Python

`eonc::Runtime` is a composition root that moves but does not copy. `ClientEON` builds one
and moves it into `Job`, and pyeonclient's `Session` owns it. The Python
server and tools use `pathlib` for path I/O throughout.

## Channels

PyPI project is `eon-akmc`. conda-forge package is `eon`, built from
the GitHub fat tarball `eon-v3.3.1.tar.xz`.
