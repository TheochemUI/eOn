---
myst:
  html_meta:
    "description": "Documentation for the structure comparison settings in eOn, used to configure similarity measures and equivalence cutoffs."
    "keywords": "eOn structure comparison, point group, bond order, configuration comparison"
---

# Structure comparison

## The client decides a match

`job = structure_comparison` reads `matter1.con` and `matter2.con`. It
writes `results.dat` with `match`, `distance`, and `per_atom_norm`.

`match` follows this section. Leave `check_rotation` off and set
`indistinguishable_atoms`. Same-element atoms may then pair in any order
inside `distance_difference`. The pairing is one-to-one. A crossed pair
can match. Two atoms cannot take one site. Set `check_rotation` alone
and the client removes the centroids. It then rotates with its Kabsch
routine and compares atoms in file order. Set both flags and the client
compares sorted neighbor-distance lists.

`remove_translation` defaults to true. Leave `check_rotation` off. Leave
no atom fixed on all three Cartesian components. The client then shifts
the first structure onto the second. The shift equals the mean
minimum-image displacement. A periodic cell that drifted as a whole still
matches. `check_rotation` skips the shift. A fully fixed atom skips it
too. The adaptive kinetic Monte Carlo server uses the same test for a
repeated saddle and for a state it already stores.

`eon.atoms.rot_match` is a separate helper. It calls iterative rotations
and assignments through `pyeonclient` when that routine is present. A
missing routine selects the helper's own Kabsch test. An error from the
call selects that test too. `eon.atoms.crystal_spacegroup` reads the cell
and the atoms.
It returns an international symbol, an international number, and a Hall
number. The function does not decide whether two clusters match.

## Configuration

```{code-block} ini
[Structure Comparison]
```



```{eval-rst}
.. autopydantic_model:: eon.schema.StructureComparisonConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: SC_
keyprefix: sc-
---
```
