---
myst:
  html_meta:
    "description": "How eOn relates to rgpot: in-process engines, RPC serve mode, and optional potserv clients."
    "keywords": "eOn, rgpot, NWChemPot, potserv, Cap'n Proto, dlopen, serve mode"
---

# rgpot integration (three roles)

eOn and [rgpot](https://github.com/OmniPotentRPC/rgpot) meet in three
independent roles: different Meson options, different binaries or symbols,
and different runtime topologies.

| Role | Meson option | What runs in eOn | Wire / load | Typical use |
| --- | --- | --- | --- | --- |
| Direct in-process | always, except Windows | Potential type `RGPOT`: links rgpot NWChemPot / CPMDPot | `dlopen` of `libnwchemc.so` / `libcpmdc.so` in the eOn address space | Production NWChem/CPMD forces in one process |
| eOn as RPC server | `-Dwith_serve=true` | `eonclient --serve` implements rgpot's Potential RPC | Cap'n Proto server in eOn | ChemGP / other tools drive any eOn pot over the network |
| RPC client to potserv | Outside eOn | Code outside eOn connects to rgpot `potserv` | Cap'n Proto client to potserv, which then `dlopen`s engines | Engine debugging via potserv |

```{important}
`RGPOT` in eOn is the direct in-process role: it links NWChemPot/CPMDPot and
loads the engine `.so` in the eOn process. Potserv is a separate rgpot binary
that hosts the same engines over Cap'n Proto RPC.
```

SocketNWChem (`potential = SocketNWChem`) is a fourth path: eOn speaks the
i-PI socket protocol to an external NWChem process. That path stays warm
across POSDATA because NWChem keeps its SCF state. Direct RGPOT relies on
nwchemc skipping redundant RTDB resets when method params are unchanged so
multi-step optimize/NEB stays competitive with the socket.

## Build flags (summary)

```{code-block} bash
# Direct NWChem/CPMD via rgpot frontends (dlopen engines). No option:
# non-Windows builds always link this arm. Cap'n Proto is required.
meson setup bbdir-rgpot

# eOn exposes its potentials as an rgpot-compatible RPC *server*
meson setup bbdir-serve -Dwith_serve=true

# Serve on top of the direct arm (the direct arm is already on)
meson setup bbdir-both -Dwith_serve=true
```

| Build | Compile define | Links from rgpot subproject |
| --- | --- | --- |
| direct RGPOT (non-Windows) | `WITH_RGPOT` | `nwchempot_dep`, `cpmdpot_dep` (frontends + schema) |
| `with_serve` | `WITH_SERVE_MODE` | `ptlrpc_dep` (RPC server stack) |

Both are convenience archives inside rgpot's single versioned `librgpot`,
so neither adds a shared object next to it at run time. Packagers that
link eOn against an installed rgpot need `librgpot` alone; the QM *engines* in
the table below stay separate because eOn `dlopen`s them.

Runtime for direct mode: set `NWCHEMC_LIBRARY` / `RGPOT_NWCHEMC_ENGINE` (or
`[RgpotPot] engine_path`) to a real `libnwchemc.so` built with NWChem embed
support.

## Docs map

- Direct pot: [RgpotPot](project:rgpot_pot.md)
- RGPOT metatomic engine: [RGPOT Metatomic](project:rgpot_metatomic.md)
- Serve mode: [Serve mode](project:serve_mode.md)
- Socket NWChem: [Potentials](project:potential.md)
- External process pot pattern: [ExtPot](project:ext_pot.md)

## Molecular ORCA energies and forces

Select the in-process rgpot engine in `config.ini` using the path to
the `orca_engine` Meson target and a flat Cap’n Proto `OrcaParams` message.

```ini
[Potential]
potential = rgpot

[RgpotPot]
backend = orca
engine_path = /path/to/liborca_engine.so
params_path = /path/to/orca.params
```

Encode the message with the effective schema from the same rgpot build.
For example, the following input specifies conventional hybrid DFT for a
neutral singlet and requests a gCP correction.

```sh
printf '%s\n' '(method = "B3LYP", basis = "def2-SV(P)", keywords = "GCP(DFT/SV(P)) TightSCF NoRI", charge = 0, multiplicity = 1)' |
    capnp encode /path/to/Potentials.capnp OrcaParams > orca.params
export RGPOT_ORCA=/path/to/orca
```

The engine runs the ORCA executable for energy and gradient evaluations
of an isolated molecule. Use a separate result-cache directory for each
ORCA installation and calculation protocol. With cache-enabled builds
of both eOn and the engine, `[Lineb] result_cache` attaches the database
through the plugin’s `rgpot_engine_set_cache` entry point. The job stores
the differences in `rgpot_engine_cache_counts` for setup, solver and output
in `evaluation_counts.rgpot_results.phases` in the CON frame metadata.
A cache miss runs the provider and a cache hit returns its stored result.

For a dynamic engine, the native job result records the version of the
rgpot library selected by Meson as `rgpot_version`. An engine exporting
`rgpot_engine_version()` also supplies its own build string, stored as
`rgpot_engine_version`. LI-NEB path and replica CON frames carry these
values in `potential.linked_version` and `potential.engine_version`.
The ORCA executable is a separate program. Record its version and the
encoded electronic-structure parameters with the calculation.
