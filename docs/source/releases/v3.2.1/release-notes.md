# eOn 3.2.1

A patch on `v3.2.0`. See monorepo `CHANGELOG.md` section 3.2.1 for the
fragment list.

## IDPP honors frozen atoms

IDPP, collective IDPP, and SIDPP no longer write the full 3N vector.
Frozen sites next to a mover stay put.

## pyeonclient optimizer

`Parameters` exposes `opt_method`, `neb_opt_method`,
`refine_opt_method`, and `refine_threshold`.

## rgpot 3.2

Wrap pin is **v3.2.0**. `potential = dftd3` / `dftd4`,
`potential = expr` with `[ExprPot]` terms (`0.5*lj + d3`), and
`potential = mopac` (libmopacc, default AM1). D3/D4 wraps stay off on
Windows, where MSVC rustc rejects those wrap flags.

## Geometry (minimage / linkcell)

`eon.geometry.pbc` uses minimage when that extra is installed.
`neighbor_list` stays vesin. `neighbor_list_linkcell` is the check
path against `linkcell.knearest`. See
[Neighbor lists](../../user_guide/neighbor_lists.md).

## ASE / Matter python surface

`from_ase` sets PBC before wrapping. Only `FixAtoms` fully freeze.
`setCell` / `setForces` dirty caches. Positions views are not writeable
in place. `to_ase` attaches SinglePoint energy/forces. NEB keeps the
GIL when the pot is not thread-safe.
