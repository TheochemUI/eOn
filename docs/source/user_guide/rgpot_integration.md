---
myst:
  html_meta:
    "description": "How eOn relates to rgpot: in-process engines, RPC serve mode, and optional potserv clients."
    "keywords": "eOn, rgpot, NWChemPot, potserv, Cap'n Proto, dlopen, serve mode"
---

# rgpot integration (three roles)

eOn and [rgpot](https://github.com/OmniPotentRPC/rgpot) meet in **three
independent roles**. They are different Meson options, different binaries or
symbols, and different runtime topologies. Do not conflate them.

```{mermaid}
flowchart LR
  subgraph direct ["Direct in-process (-Dwith_rgpot)"]
    E1[eOn RGPOT pot] --> NP[rgpot NWChemPot / CPMDPot]
    NP -->|dlopen| L1[libnwchemc.so / libcpmdc.so]
  end
  subgraph serve ["eOn as RPC server (-Dwith_serve)"]
    Client[External client e.g. ChemGP] -->|Cap'n Proto| E2[eonclient --serve]
    E2 --> Any[Any eOn Potential]
  end
  subgraph rpccli ["Optional: RPC client to potserv"]
    E3[External potserv client] -->|Cap'n Proto| PS[rgpot potserv]
    PS -->|dlopen| L2[libnwchemc.so / …]
  end
```

| Role | Meson option | What runs in eOn | Wire / load | Typical use |
| --- | --- | --- | --- | --- |
| **Direct in-process** | `-Dwith_rgpot=true` | Potential type **`RGPOT`**: links rgpot **NWChemPot** / **CPMDPot** | `dlopen` of `libnwchemc.so` / `libcpmdc.so` in the eOn address space | Production NWChem/CPMD forces without a second process |
| **eOn as RPC server** | `-Dwith_serve=true` | `eonclient --serve` implements rgpot's Potential RPC | Cap'n Proto **server** in eOn | ChemGP / other tools drive *any* eOn pot over the network |
| **RPC client → potserv** | *Not* an eOn pot type | Something *outside* eOn (or a thin harness) connects to **rgpot `potserv`** | Cap'n Proto **client** → potserv, which then `dlopen`s engines | Debugging engines via potserv; not the preferred eOn path for NWChem |

```{important}
**RGPOT in eOn is the direct in-process role.** It does **not** speak Cap'n
Proto to potserv. Potserv is an rgpot binary for hosting engines over RPC; eOn
uses the same engines by linking NWChemPot/CPMDPot and loading the `.so` itself.
```

Compare also **SocketNWChem** (`potential = SocketNWChem`): eOn speaks the i-PI
socket protocol to an external **NWChem process** (not rgpot, not potserv).
That path stays warm across POSDATA because NWChem keeps its SCF state. Direct
RGPOT relies on **nwchemc** skipping redundant RTDB resets when method params
are unchanged so multi-step optimize/NEB stays competitive with the socket.

## Build flags (summary)

```{code-block} bash
# Direct NWChem/CPMD via rgpot frontends (dlopen engines)
meson setup bbdir-rgpot -Dwith_rgpot=true

# eOn exposes its potentials as an rgpot-compatible RPC *server*
meson setup bbdir-serve -Dwith_serve=true

# Both (independent features; subproject shared)
meson setup bbdir-both -Dwith_rgpot=true -Dwith_serve=true
```

| Option | Compile define | Links from rgpot subproject |
| --- | --- | --- |
| `with_rgpot` | `WITH_RGPOT` | `nwchempot_dep`, `cpmdpot_dep` (frontends + schema) |
| `with_serve` | `WITH_SERVE_MODE` | `ptlrpc_dep` (RPC server stack) |

Runtime for direct mode: set `NWCHEMC_LIBRARY` / `RGPOT_NWCHEMC_ENGINE` (or
`[RgpotPot] engine_path`) to a real `libnwchemc.so` built with NWChem embed
support.

## Docs map

- Direct pot: [RgpotPot](project:rgpot_pot.md)
- Serve mode: [Serve mode](project:serve_mode.md)
- Socket NWChem (non-rgpot): [Potentials](project:potential.md)
- External process pot pattern: [ExtPot](project:ext_pot.md)
