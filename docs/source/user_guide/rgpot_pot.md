---
myst:
  html_meta:
    "description": "eOn RGPOT potential: in-process rgpot NWChemPot/CPMDPot via dlopen."
    "keywords": "eOn, RGPOT, rgpot, NWChemPot, CPMDPot, libnwchemc, dlopen"
---

# RgpotPot (direct in-process rgpot)

```{versionadded} 2.16.0
```

On non-Windows builds, potential type `RGPOT` links
[rgpot](https://github.com/OmniPotentRPC/rgpot) NWChemPot / CPMDPot and loads
`libnwchemc.so` / `libcpmdc.so` with `dlopen` in the eOn process.

Sibling roles live elsewhere: potserv Cap'n Proto clients, `eonclient --serve`
([Serve mode](project:serve_mode.md)), and SocketNWChem (i-PI to a standalone
NWChem binary). Overview: [rgpot integration](project:rgpot_integration.md).

## Build

```{code-block} bash
meson setup bbdir-rgpot -Dwith_tests=true
meson compile -C bbdir-rgpot
```

Requires Cap'n Proto headers and libs (method params are Cap'n Proto messages
passed into the C ABI) and the `rgpot` Meson subproject
(`subprojects/rgpot.wrap`). The build pulls `nwchempot_dep` / `cpmdpot_dep`.
Serve mode uses `ptlrpc_dep` under `-Dwith_serve`.

Engines are resolved at runtime:

- `NWCHEMC_LIBRARY` or `RGPOT_NWCHEMC_ENGINE` (NWChem embed library), or
- `[RgpotPot] engine_path` / `engine_library` in the config.

`enginePath` in the Cap'n Proto blob is stripped before the ABI call (host-only
locator); see rgpot NWChemPot.

## Config

```{code-block} ini
[Main]
job = point

[Potential]
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

For Metatomic (dlopen engine; no fat metatomic/torch link into eOn):

```{code-block} ini
[Potential]
potential = RGPOT

[RgpotPot]
backend = metatomic
model_path = /path/to/model.pt
device = cpu
# engine_path = /path/to/libmetatomic_engine.so
# or: export RGPOT_METATOMIC_ENGINE=...
```

For xTB (preferred packaging path; leaves `-Dwith_xtb=false`):

```{code-block} ini
[Potential]
potential = RGPOT

[RgpotPot]
backend = xtb
paramset = GFN2xTB
accuracy = 1.0
# engine_path = /path/to/libxtb_engine.so
# or: export RGPOT_XTB_ENGINE=...
```

Optional `input_block` (or env `RGPOT_NWCHEM_INPUT_BLOCK`) supplies NWChem
`inputBlocks` (e.g. explicit `dft` / `xc` stanzas). When `theory=dft` and
`scf_type` looks like an XC label (e.g. `b3lyp`), a minimal DFT block is
emitted automatically.

For `backend = cpmdc`, `input_block` (or env `RGPOT_CPMD_INPUT_BLOCK`) carries
CPMD `&SECTION` text. cpmdc places it ahead of the sections it generates, so
a periodic `&SYSTEM` or a full `&DFT` given here takes the place of the
isolated cold deck. `permanent_dir` sets the CPMD `FILEPATH`, where the
`RESTART` files go; `scratch_dir` is the fallback. The pseudopotential
directory comes from `CPMDC_PSEUDO_DIR` or `CPMD_PP_LIBRARY_PATH`.

```{code-block} ini
[RgpotPot]
backend = cpmdc
permanent_dir = /scratch/cpmd-restart
input_block = &SYSTEM
    ANGSTROM
    CELL VECTORS
      10.26 0.0 0.0
      0.0 10.26 0.0
      0.0 0.0 10.26
    CUTOFF
      30.0
  &END
```

### One CPMD session per NEB image

`ranks_per_image` splits the MPI world into calculator groups of that many
ranks. Each group runs its own CPMD session on its own subcommunicator, and a
NEB hands image *j* of each band update to group *j* mod *G*. Every rank then
receives every image's energy and forces, so all ranks hold the same band. A
single force call (an endpoint, a minimization, the dimer of OCI-NEB) runs on
group 0 and is shared the same way. With as many groups as images each group
keeps its own image's wavefunction from one iteration to the next.

Launch one `eonclient` per rank, with the world a multiple of
`ranks_per_image`; seven images on six ranks each is 42 ranks:

```{code-block} ini
[RgpotPot]
backend = cpmdc
ranks_per_image = 6
```

```{code-block} bash
mpirun -np 42 eonclient
```

This needs rgpot built with MPI (`-Drgpot:with_mpi=enabled`) and a libcpmdc
that exports `cpmdc_bind_calculator`, on an OpenCPMD with cpmdc's
`opencpmd_mp_comm_set.patch`.

Installed rgpot ≥ 2.5.0 is preferred via `pkg-config` (`dependency('rgpot')`);
the Meson wrap is the fallback for hermetic/dev builds.

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
