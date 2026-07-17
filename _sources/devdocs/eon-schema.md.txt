---
myst:
  html_meta:
    "description": "eon-schema: shared L0 Cap'n Proto, L1 job config, L2 API specs, and INI helpers for eon-akmc and pyeonclient."
    "keywords": "eon-schema, Cap'n Proto, pydantic, config.ini, pyeonclient, eon-akmc"
---

# eon-schema (shared parameter package)

```{versionadded} 2.17
```

**eon-schema** is the single Python home for eOn parameter surfaces used by both
**eon-akmc** and **pyeonclient**. It is a PyPI split (`pip install eon-schema`);
the conda-forge fat tarball still ships the full monorepo (including this package
under `packages/eon-schema/` and the authoring Cap’n Proto under `schema/`).

## What it owns

| Layer | Import | Purpose |
|-------|--------|---------|
| **L0** | `eon_schema.ssot` | Cap’n Proto field graph + catalog (vendored from monorepo `schema/`) |
| **L1** | `eon_schema.config` | Full job-config pydantic models (`MainConfig`, `Metatomic`, `Config`, …) |
| **L2** | `eon_schema.api`, `eon_schema.fields` | In-process specs (`DimerSpec`, `NebSpec`, enums) |
| **INI** | `eon_schema.config.write_ini`, `hydrate_ini`, … | Author / hydrate `config.ini` without eon-akmc |

## What consumers re-export

| Package | Re-export | Notes |
|---------|-----------|--------|
| **eon-akmc** | `eon.schema` → L1 | Sphinx `autopydantic_model:: eon.schema.*` unchanged |
| **pyeonclient** | `pyeonclient.models` → L2 | Requires `pyeonclient[models]` → `eon-schema>=0.2` |

Do **not** re-author L1 under `eon/schema.py` or duplicate L2 under
`pyeonclient/models.py`. Those modules are thin shims.

## Dependencies

- **eon-akmc** hard-depends on `eon-schema>=0.2`
- **pyeonclient[models]** depends on `eon-schema>=0.2`
- **pydantic≥2** is a hard dependency of eon-schema (L1 + L2)

Downstream tooling (e.g. **rgpycrumbs**) should write and hydrate job configs
via eon-schema INI helpers instead of importing full eon-akmc
`eon.config.ConfigClass`. That cutover is tracked in the rgpycrumbs issue
tracker, not completed in this tree.

## INI helpers

```{code-block} python
from eon_schema.config import (
    MainConfig,
    PotentialConfig,
    DimerConfig,
    write_models_ini,
    hydrate_ini,
    unknown_ini_keys,
    write_ini,
)

write_models_ini(
    "config.ini",
    MainConfig(job="saddle_search"),
    PotentialConfig(potential="lj"),
    DimerConfig(dimer_improved=True),
    validate=True,
)

cfg = hydrate_ini({"Main": {"job": "minimization"}}, base_sections=["Main"])
assert not unknown_ini_keys(cfg)
write_ini("config.ini", cfg)
```

Covered L0 sections (Main, Potential, Optimizer, Structure Comparison, Process
Search) can be key-checked with `unknown_ini_keys`. Other pot sections
(SocketNWChemPot, Metatomic, …) are free-form relative to the Cap’n Proto
catalog.

## Authoring

1. **L0 fields/defaults:** edit monorepo `schema/eon_params.capnp`, run
   `python tools/params_ssot/codegen.py`, then
   `./packages/eon-schema/scripts/sync_ssot_into_package.sh`.
2. **L1 models:** edit
   `packages/eon-schema/src/eon_schema/config/models.py` (keep parity with L0
   for covered groups; `tests/test_params_ssot.py`).
3. **L2 specs:** edit `packages/eon-schema/src/eon_schema/api/`.

Package README and release checklist:
{file}`../../../packages/eon-schema/README.md`,
{file}`../../../packages/eon-schema/PUBLISHING.md`.

## Fat vs split

| Artifact | Role |
|----------|------|
| Fat `eon-v*.tar.xz` | conda-forge feedstock; full tree |
| PyPI `eon-schema` | Shared schema only; independent `0.y.z` |
| PyPI `eon-akmc` / `pyeonclient` | Depend on / re-export eon-schema |

See also [Release process](project:release.md) (fat vs split table).
