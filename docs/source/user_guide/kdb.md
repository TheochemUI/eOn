---
myst:
  html_meta:
    "description": "Guide to using the Kinetic Database (KDB) in eOn to store and reuse information about kinetic processes, speeding up aKMC simulations."
    "keywords": "eOn Kinetic Database, KDB, aKMC acceleration, process recycling, readcon-db, amsel"
---

# Kinetic database

One of the bottlenecks in an aKMC simulation is performing the saddle point
searches. The kinetic database of
{cite:t}`kdb-terrellDatabaseAtomisticReaction2012` stores each good process
as it is registered and uses it to start the next search of a matching state.

In the following figure, the hydrogen of a carboxyl group on an Au(111) surface
transfers to the other oxygen (a). In this process, the hydrogen is determined
to be the only moving atom, and the two oxygen atoms to be its neighbors. The
other atoms are stripped from the system and the resulting configurations are
stored in the database (b). If in the future the system passes through a state
with a local configuration closely resembling either minimum of (b), the kinetic
database will suggest a saddle to converge with a dimer search. The dimer starts
with this suggested configuration and mode, and if it is a good suggestion,
converges very rapidly to the saddle.

```{figure} ../fig/akmc-1.png
---
alt: Carboxyl group on an Au(111)
class: full-width
align: center
---
Snapshots of Carboxyl group on an Au(111). (a) Hydrogen of carboxyl transfers to another oxygen. (b) Other atoms are stripped and stored.
```

## Where a process is stored

The reactant, saddle, and product frames go into the run's readcon-db
corpus (`readcon.db` next to `config.ini`). The barrier in eV, the
prefactor in s^-1, the mode, and the readcon-db frame keys go into
`amsel.KdbStore`. The catalog directory is `Paths.kdb` (default
`<main_directory>/kdb/`). `KdbStore` opens that directory.

A suggestion refines from the stored saddle. The direction is the stored
mode when its absolute cosine with the reactant-to-saddle vector is at
least `kdb_mac`, and that vector otherwise. Each stored process is
offered once. The next search is a random displacement when none remain.
With `kdb_only = true` and an empty catalog, no random search is submitted.

`use_kdb = true` does not import the PyPI `kdb` package. It needs
`amsel` and `readcon-db`. A missing `amsel` logs `amsel is not installed`
and leaves the state unmarked, so a later iteration can try again.

The three match numbers are `kdb_nf`, `kdb_dc`, and `kdb_mac`:

- `kdb_nf` is the neighbor fudge, a fraction. The default is 0.2.
- `kdb_dc` is the distance cutoff in angstroms. The default is 0.3.
  A stored reactant matches the current state when every atom is within
  `kdb_dc * (1 + kdb_nf)` angstroms.
- `kdb_mac` is the minimum absolute cosine between the stored mode and
  the reactant-to-saddle displacement. The default is 0.7. At or above
  that cosine the refine uses the stored mode, flipped when the cosine
  is negative. Below it, the refine uses the reactant-to-saddle vector.

A query log line contains the `kdb_nf` value from the ini.

## Configuration

```{code-block} ini
[KDB]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.KDBConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: KDB_
keyprefix: kdb-
---
```
