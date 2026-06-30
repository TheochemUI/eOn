# CPMD BLYP via rgpot potserv

1. Build rgpot with RPC + tests (fake engine) or real `libcpmdc.so`.
2. Start server:
   ```bash
   export CPMDC_LIBRARY=/path/to/libcpmdc.so   # or libcpmdc_fake_engine.so
   potserv 12346 CPMD
   ```
3. Build eOn with `-Dwith_rgpot=true` and run a point job with this `config.ini`
   and a `pos.con` (e.g. copy from `client/unit_tests/data/systems/...`).

Env overrides: `RGPOT_POTSERV_HOST`, `RGPOT_POTSERV_PORT`, `RGPOT_BACKEND=CPMD`.
