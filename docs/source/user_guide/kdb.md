# Kinetic Database

One of the bottlenecks in an aKMC simulation is performing the saddle point
searches. The kinetic database has been demonstrated by
{cite:t}`kdb-terrellDatabaseAtomisticReaction2012` to ameliorate this cost by
storing information about processes as they are found and using it to predict
future saddle points.

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

## Dependencies

The kinetic database is contained in a library separate from `eOn`. It is part
of the `tsase` python module [located here](http://theory.cm.utexas.edu/tsase). 

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
