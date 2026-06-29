---
myst:
  html_meta:
    "description": "eOn RGPOT potential: in-process rgpot NWChemPot/CPMDPot via dlopen (not potserv RPC)."
    "keywords": "eOn, RGPOT, rgpot, NWChemPot, CPMDPot, libnwchemc, dlopen"
---

# RgpotPot (direct in-process rgpot)

```{versionadded} TBD
```

When eOn is built with **`-Dwith_rgpot=true`**, potential type **`RGPOT`** links
[rgpot](https://github.com/OmniPotentRPC/rgpot) **NWChemPot** / **CPMDPot** and
loads `libnwchemc.so` / `libcpmdc.so` with **`dlopen` in the eOn process**.

This is **not**:

- Cap'n Proto to **potserv** (that is an *external* RPC client role), nor
- **`eonclient --serve`** (that is eOn acting as an *RPC server*; see
  [Serve mode](project:serve_mode.md)), nor
- **SocketNWChem** (i-PI socket to a standalone NWChem binary).

For the three-role overview, see [rgpot integration](project:rgpot_integration.md).

## Build

```{code-block} bash
meson setup bbdir-rgpot -Dwith_rgpot=true -Dwith_tests=true
meson compile -C bbdir-rgpot
```

Requires Cap'n Proto **headers/libs** (method params are Cap'n Proto messages
passed into the C ABI) and the `rgpot` Meson subproject
(`subprojects/rgpot.wrap`). The build pulls `nwchempot_dep` / `cpmdpot_dep`
only — not `ptlrpc_dep` (that is for serve mode).

Engines are resolved at **runtime**:

- `NWCHEMC_LIBRARY` or `RGPOT_NWCHEMC_ENGINE` (NWChem embed library), or
- `[RgpotPot] engine_path` / `engine_library` in the config.

`enginePath` in the Cap'n Proto blob is stripped before the ABI call (host-only
locator); see rgpot NWChemPot.

## Config

```{code-block} ini
[Main]
job = point
potential = RGPOT

[RgpotPot]
backend = nwchemc
basis = sto-3g
theory = scf
scf_type = rhf
charge = 0
multiplicity = 1
# engine_path = /path/to/libnwchemc.so
```

For CPMD:

```{code-block} ini
[RgpotPot]
backend = cpmdc
functional = BLYP
cutoff_ry = 70.0
```

Optional `input_block` (or env `RGPOT_NWCHEM_INPUT_BLOCK`) supplies NWChem
`inputBlocks` (e.g. explicit `dft` / `xc` stanzas). When `theory=dft` and
`scf_type` looks like an XC label (e.g. `b3lyp`), a minimal DFT block is
emitted automatically.

## vs SocketNWChem

| | SocketNWChem | RGPOT (this pot) |
| --- | --- | --- |
| Protocol | i-PI socket; eOn listens | In-process `dlopen` via rgpot frontends |
| Engine process | External NWChem | `libnwchemc.so` / `libcpmdc.so` in eOn |
| Multi-call SCF | Warm NWChem across POSDATA | nwchemc warm params cache (skip full RTDB reset when method blob unchanged) |

## Implementation notes

- `RgpotPot` (eOn `Potential`) owns an opaque `RGPotEngine` TU that includes
  **only** rgpot headers — avoids Cap'n Proto type name `Potential` colliding
  with eOn's `Potential` class.
- Forces and energies use eOn units (eV, eV/Å) after rgpot conversion from
  Hartree / Hartree·bohr⁻¹.
