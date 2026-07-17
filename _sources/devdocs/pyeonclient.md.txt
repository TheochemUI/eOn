# pyeonclient (nanobind client bindings)

`pyeonclient` exposes the eOn **C++ client core** to Python with
[nanobind](https://nanobind.readthedocs.io/) only (no pybind11).

| Type | Role |
|------|------|
| `Parameters` | `config.ini` / JSON parameters |
| `Potential` / `make_potential` | forcefield factory |
| `Matter` | structure, energy/forces, `.con` I/O |

## ABI

- **Stable ABI (abi3)** for normal CPython ≥ 3.12.
- **Free-threaded** (`NB_FREE_THREADED`) when the build interpreter has
  `Py_GIL_DISABLED` (CPython ≥ 3.13t). There is no stable ABI for free-threaded
  CPython yet (nanobind / PEP 803); the build picks free-threaded and skips
  abi3 in that case.

## Build

```bash
meson setup build -Dwith_pyeonclient=true
meson compile -C build
```

Requires `nanobind>=2.2` on the build Python (`pixi.toml` pins it).

## Tests

```bash
pytest tests/test_pyeonclient_matter.py -v
```

## Server bridge

```python
from eon import pyeonclient_bridge as br
# br.to_matter(structure, pot, params)
# br.from_matter(matter)
```

## Why not pybind11 embed for ASE?

pybind11 was used only so **standalone `eonclient`** could *embed* a Python
interpreter and call ASE calculators. nanobind does **not** ship a first-class
`embed.h` API; the modern polarity is inverted:

| Mode | Who owns Python | ASE / Matter |
|------|-----------------|--------------|
| **Server + pyeonclient** | Python process | Matter in-process; ASE calculators can live as Python objects (no C++ embed) |
| **Standalone eonclient** | C++ process | Prefer C++ pots, or `NbGuard` + CPython init; migrate ASE pots off pybind11 |

`client/NbGuard.h` initializes the interpreter with `Py_InitializeEx` (no
pybind11). ASE pot migration to nanobind marshalling is follow-on work; new
code should not add pybind11 surfaces.

## Matter through-and-through

```text
explorer / AKMC
    │  Structure (server)
    ▼
communicator type = local_lib | inprocess
    │  structure_to_matter
    ▼
pyeonclient.Matter.relax()   # (more Job types as bound)
    │  matter_to_structure
    ▼
results (Structure + results.dat fields)
```

File-based `local` / `cluster` remain for HPC; geometry on disk is still
readcon `.con` when those paths are used.

## PyPI

See [pyeonclient-pypi](pyeonclient-pypi.md) for the standalone package and wheel matrix.
