# RgpotPot (rgpot potserv client)

When eOn is built with `-Dwith_rgpot=true`, the potential type **`RGPOT`** connects
to an [rgpot](https://github.com/OmniPotentRPC/rgpot) **potserv** process over
Cap'n Proto and evaluates energy/forces via **NWChem** or **CPMD** backends
(`dlopen` of `libnwchemc` / `libcpmdc` on the server).

## Build

```bash
meson setup bbdir-rgpot -Dwith_rgpot=true -Dwith_tests=true
meson compile -C bbdir-rgpot
```

Requires Cap'n Proto and the `rgpot` Meson subproject (`subprojects/rgpot.wrap`).

## Server

```bash
# Real or fake engine on the potserv host
NWCHEMC_LIBRARY=/path/to/libnwchemc.so potserv 12345 NWChem
CPMDC_LIBRARY=/path/to/libcpmdc.so     potserv 12346 CPMD
```

## INI

```ini
[Potential]
potential = RGPOT

[RgpotPot]
host = 127.0.0.1
port = 12345
backend = NWChem
nwchem_basis = sto-3g
nwchem_theory = scf
nwchem_scf_type = rhf
# backend = CPMD
# cpmd_functional = BLYP
# cpmd_task = gradient
```

Environment overrides: `RGPOT_POTSERV_HOST`, `RGPOT_POTSERV_PORT`, `RGPOT_BACKEND`.

## vs SocketNWChem

| | SocketNWChem | RGPOT |
|--|--------------|--------|
| Protocol | i-PI socket; eOn listens | Cap'n Proto; potserv listens |
| Engine | NWChem process | potserv + libnwchemc/libcpmdc |
| CPMD | no | yes (`backend = CPMD`) |

See `scripts/rgpot_nwchem_vs_socket.sh` and `examples/rgpot_cpmd_blyp/`.
