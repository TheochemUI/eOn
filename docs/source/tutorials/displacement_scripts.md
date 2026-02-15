---
myst:
  html_meta:
    "description": "Tutorial on using targeted displacement scripts and atom lists for saddle searches in eOn."
    "keywords": "eOn displacement, saddle search, vacancy diffusion, adsorbate, PTM, OVITO, targeted displacement"
---

# Targeted Displacement for Saddle Searches

By default, eOn picks a random atom as the epicenter for each saddle search
displacement. This works well for small systems or when all atoms are equally
likely to participate in reactions. For larger, more structured systems —
surfaces with defects, supported catalysts, or bulk crystals with vacancies —
random selection wastes effort: most saddle searches start far from the action
and fail to find relevant transitions.

eOn provides two mechanisms for targeting displacements at the atoms that
matter:

1. **Static atom list** (`displace_atom_list`) — a fixed set of atom indices
   written directly in `config.ini`. Best when the active atoms are known in
   advance and don't change between AKMC states.

2. **Dynamic script** (`displace_atom_kmc_state_script`) — a Python script
   that is executed once per new AKMC state. The script receives the current
   geometry, analyses it, and returns the indices of atoms that should be
   displaced. Best when the active region moves (e.g. a migrating vacancy).

## Script Interface

A displacement script must satisfy the following contract:

- It receives the path to a `.con` file as its **sole positional argument**.
- It must print a **comma-separated list of 0-based atom indices** to
  **stdout** (e.g. `3, 7, 12, 45`).
- All other output (logging, warnings, progress bars) must go to **stderr**.
  Anything on stdout that is not a valid index list will cause a parse error.
- The script is run **once per new AKMC state**; the result is cached in
  `state.info` and reused for all saddle searches launched from that state.
- The path in `displace_atom_kmc_state_script` can be **relative** (resolved
  against the eOn root directory) or **absolute**.
- Scripts can use [PEP 723](https://peps.python.org/pep-0723/) inline metadata
  so they are runnable with `uvx` without a separate virtual environment.

## Example 1: Vacancy Diffusion in Cu

The `examples/akmc-cu-vacancy/` directory contains a complete AKMC setup for a
Cu crystal with a single vacancy, using OVITO's Polyhedral Template Matching
(PTM) to identify defect atoms.

### The Idea

In a nearly perfect FCC crystal, the only atoms likely to participate in
vacancy migration are those whose local environment deviates from bulk FCC —
i.e. the atoms adjacent to the vacancy. PTM classifies each atom's local
neighbourhood into crystal structure types (FCC, HCP, BCC, icosahedral, or
"other"). Atoms that are *not* FCC are the natural candidates for displacement.

### Configuration

The relevant `config.ini` block:

```{code-block} ini
:caption: examples/akmc-cu-vacancy/config.ini (excerpt)

[Saddle Search]
method = min_mode
min_mode_method = dimer
max_energy = 300
displace_listed_atom_weight = 1.0
displace_radius = 3.0
displace_magnitude = 0.01
displace_atom_kmc_state_script = ptmdisp.py
displace_all_listed = true
```

Key settings:

- `displace_listed_atom_weight = 1.0` with all other weights at their defaults
  (0) ensures that *only* listed atoms are used as epicenters.
- `displace_all_listed = true` means every atom returned by the script is
  displaced simultaneously, not just one picked at random.
- `displace_radius = 3.0` also includes neighbours within 3 Å of any listed
  atom.
- `displace_atom_kmc_state_script = ptmdisp.py` points to the script (relative
  to the eOn root directory).

### The PTM Script

`ptmdisp.py` is a [PEP 723](https://peps.python.org/pep-0723/) inline script.
Its core logic:

1. Reads the `.con` file with ASE.
2. Converts the structure to an OVITO pipeline.
3. Applies the PTM modifier to classify each atom.
4. Selects atoms that are **not** FCC.
5. Prints their indices to stdout.

You can run it standalone to inspect the output:

```{code-block} bash
uvx ptmdisp.py pos.con
# Output: 42, 43, 51, 67, ...

uvx ptmdisp.py pos.con --verbose
# Prints logging info to stderr, indices to stdout
```

### Data Flow

```
New AKMC state
  │
  ▼
eOn server runs: ptmdisp.py reactant.con
  │
  ▼
Script prints indices to stdout  →  parsed into displace_atom_list
  │
  ▼
Indices cached in state.info
  │
  ▼
DisplacementManager uses list for all saddle searches in this state
```

## Example 2: Adsorbate on a Catalyst Surface

Consider a CO molecule adsorbed on a Pt(111) surface. During an AKMC
simulation the adsorbate can diffuse, rotate, or desorb — but the relevant
atoms are always the adsorbate itself plus a handful of nearby surface atoms.
Displacing random bulk Pt atoms far from the CO would be wasteful.

### Strategy

1. Identify the adsorbate atoms (C and O) by element.
2. Expand the selection to include all Pt atoms within a cutoff radius of any
   adsorbate atom.
3. Return the combined set of indices.

### The Script

The file `examples/akmc-cu-vacancy/adsorbate_region.py` implements this
approach using only ASE (no OVITO dependency):

```{code-block} bash
uvx adsorbate_region.py pos.con --adsorbate-elements C O --cutoff 4.0
# Output: 0, 1, 23, 24, 31, ...
```

The script:

1. Reads the `.con` file with `ase.io.read`.
2. Finds atoms matching the specified adsorbate elements.
3. Computes distances from every other atom to the nearest adsorbate atom.
4. Selects all atoms within the cutoff distance.
5. Merges the adsorbate indices with the nearby-atom indices and prints the
   result.

```{tip}
You can also select adsorbate atoms by z-coordinate instead of element, which
is useful when the adsorbate and surface share the same element. Pass
`--z-above <threshold>` to select atoms above a given height.
```

### Configuration

```{code-block} ini
:caption: config.ini (excerpt)

[Saddle Search]
displace_listed_atom_weight = 1.0
displace_radius = 3.0
displace_magnitude = 0.01
displace_atom_kmc_state_script = adsorbate_region.py --adsorbate-elements C O --cutoff 4.0
displace_all_listed = true
```

## Static List Alternative

For systems where the active atoms are known ahead of time and do not change
between states, a static `displace_atom_list` is simpler — no script needed:

```{code-block} ini
:caption: config.ini (excerpt)

[Saddle Search]
displace_atom_list = 0, 1, 2
displace_listed_atom_weight = 1.0
displace_all_listed = true
displace_radius = 3.0
displace_magnitude = 0.01
```

This displaces atoms 0, 1, and 2 (plus their neighbours within 3 Å) on every
saddle search. It is suitable for small molecules, known defect sites, or any
situation where the active region is fixed.

## Client-Side `listed_atoms`

The displacement can also be handled entirely by the C++ client, without the
server writing a displacement file. Set:

```{code-block} ini
[Saddle Search]
client_displace_type = listed_atoms
displace_atom_list = 0, 1, 2
```

In this mode the client reads `displace_atom_list` directly from the INI config
and uses those atoms as epicenter candidates. This is useful for simple cases
where server-side scripting is unnecessary.
