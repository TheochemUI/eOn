# CPMD BLYP via in-process RgpotPot

`RgpotPot` runs rgpot's `CPMDPot` frontend inside `eonclient` (no potserv, no
sockets); the CPMD engine library is `dlopen`ed at runtime.

1. Build a CPMD engine: real `libcpmdc.so` from
   [OmniPotentRPC/cpmdc](https://github.com/OmniPotentRPC/cpmdc), or the fake
   engine from an rgpot build with `-Dwith_rpc=true -Dwith_tests=true`
   (`libcpmdc_fake_engine.so`).
2. Build eOn (non-Windows builds link RGPOT; Cap'n Proto is required).
3. Point the potential at the engine and run a point job with this
   `config.ini` and a `pos.con` (e.g. copy from
   `client/unit_tests/data/systems/...`):

   ```bash
   export CPMDC_LIBRARY=/path/to/libcpmdc.so   # or libcpmdc_fake_engine.so
   eonclient
   ```

Engine lookup order: `[RgpotPot] engine_path` / `engine_library`, then the
`CPMDC_LIBRARY` / `RGPOT_CPMDC_ENGINE` / `RGPOT_CPMD_ENGINE` environment
variables. `RGPOT_BACKEND=CPMD` overrides the configured backend.

For out-of-process RPC instead, use rgpot's `potserv <port> CPMD` with an RPC
client; that path is separate from this in-process example.
