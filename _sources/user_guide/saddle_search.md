---
myst:
  html_meta:
    "description": "Configuration for initiating saddle searches in eOn by displacing atoms from a minimum energy state."
    "keywords": "eOn saddle search, transition state, displacement, potential energy surface"
---

# Saddle Search

A saddle search is initiated by making a local displacement of atoms from their
position at the minimum of the current state. This displacement can be done
using the different strategies indicated by the {any}`client_displace_type`
option, and the following parameters. If the user knows something about the
local environment where reactions are likely to take place in the system, this
information can be used to make saddle searches more efficient by getting them
started in the right part of configuration space.

## Displacement Strategies

### Epicenters and Weight-Based Selection

Each saddle search begins by choosing an **epicenter** — a single atom around
which the initial displacement is constructed. The epicenter is selected
probabilistically from several strategies, each controlled by a relative weight:

| Weight parameter                       | Strategy                                                       |
| -------------------------------------- | -------------------------------------------------------------- |
| `displace_random_weight`               | Pick any atom at random                                        |
| `displace_least_coordinated_weight`    | Pick an atom with the lowest coordination number               |
| `displace_not_FCC_HCP_weight`          | Pick an atom whose local environment is neither FCC nor HCP    |
| `displace_under_coordinated_weight`    | Pick an atom with coordination ≤ `displace_max_coordination`   |
| `displace_listed_atom_weight`          | Pick an atom from `displace_atom_list`                         |
| `displace_listed_type_weight`          | Pick an atom whose type is in `displace_type_list`             |

The weights do not need to sum to 1 — they are normalised internally. Setting a
weight to 0 disables that strategy. For example, setting only
`displace_listed_atom_weight = 1.0` and all other weights to 0 ensures that
every saddle search starts from an atom in the explicit list.

### Radius and Magnitude

Once the epicenter is chosen, all atoms within `displace_radius` of the
epicenter are displaced. Each displaced atom receives a random perturbation
drawn from a Gaussian distribution with standard deviation
`displace_magnitude` (in Ångströms) independently in each Cartesian direction.

### Displacing All Listed Atoms

When `displace_all_listed` is **true** and a listed-atom strategy is selected,
*every* atom in `displace_atom_list` (or `displace_type_list`) is displaced —
not just one chosen at random. Atoms within `displace_radius` of *any*
displaced atom are also included. Set `displace_radius` to 0 to restrict the
displacement strictly to the listed atoms.

### Dynamic Atom Lists via Scripts

For systems where the relevant atoms change from state to state (e.g. a
migrating vacancy), a static list is insufficient. The
`displace_atom_kmc_state_script` option lets you specify a Python script that
is executed once per new AKMC state to determine the atom list dynamically. See
the {doc}`../tutorials/displacement_scripts` tutorial for worked examples
covering vacancy diffusion and adsorbate-on-surface scenarios.

### Client-Side Displacement

The `client_displace_type` option selects how the **client** (C++ code) picks
the epicenter when the displacement is performed client-side rather than by the
server:

- `random` — uniform random atom
- `last_atom` — the last atom in the configuration
- `min_coordinated` — the atom with the fewest neighbours
- `not_fcc_or_hcp` — an atom whose local structure is neither FCC nor HCP
- `listed_atoms` — an atom from `displace_atom_list` (parsed from the INI
  config; no server displacement file needed)
- `load` — read a displacement vector from a file written by the server

## Configuration

```{code-block} ini
[Saddle Search]
```


```{versionchanged} 3.1_TBA
In TOML, this will be `[Saddle_Search]`
```


```{eval-rst}
.. autopydantic_model:: eon.schema.SaddleSearchConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: SS_
keyprefix: ss-
---
```
