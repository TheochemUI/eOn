# eon-schema

**Shared schema package for both eon-akmc and pyeonclient.** One home for:

| Layer | Module | Role |
|-------|--------|------|
| **L0** | `eon_schema.ssot` | Cap’n Proto params field graph (vendored from monorepo `schema/`) |
| **L0** | `eon_schema.jobs` | Cap’n Proto `JobRequest` / `JobResult` / `Geometry` + `results.dat` adapters |
| **L1** | `eon_schema.config` | Full job-config pydantic models (`MainConfig`, `Metatomic`, `Config`, …) |
| **L2** | `eon_schema.api` / `fields` | In-process specs (`DimerSpec`, `NebSpec`, enums) |

Consumers:

- **eon-akmc** — depends on `eon-schema`; `eon.schema` re-exports L1 for Sphinx / legacy imports.
- **pyeonclient** — `pyeonclient[models]` pulls `eon-schema`; `pyeonclient.models` re-exports L2.

Do **not** re-author job models under `eon/` or duplicate L2 under `pyeonclient/`.

## Fat tree vs this package

| | Fat monorepo tarball | This PyPI package |
|--|----------------------|-------------------|
| Purpose | conda-forge `eon`, full C++/server builds | Shared Python schema for both clients |
| Cap’n Proto | Full monorepo including `schema/` | Vendored under `src/eon_schema/ssot/` |
| L1 models | Also present via `eon.schema` re-export | **Authoritative** under `eon_schema.config` |
| Feedstock | **Uses fat tarball only** | Not required for 0.2.x |

**L0 authoring** is monorepo `schema/eon_params.capnp`. After edits:

```bash
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh
```

**L1 authoring** is `packages/eon-schema/src/eon_schema/config/models.py`.
Keep field parity with L0 for covered sections (`tests/test_params_ssot.py`).

See **[PUBLISHING.md](PUBLISHING.md)** for fat vs split release trains.

## Install

```bash
pip install eon-schema            # L0 + L1 + L2 (pydantic is a hard dep)
# monorepo:
pip install -e packages/eon-schema
```

## Quick use

```python
# L0 params
from eon_schema.ssot import capnp_path, load_catalog

# L0 job envelope
from eon_schema.jobs import job_result_capnp_path, results_dat_to_dict

# L1 (same objects as eon.schema.*)
from eon_schema.config import MainConfig, Metatomic, Config

# L2 (same objects as pyeonclient.models.*)
from eon_schema.api import DimerSpec
from eon_schema.fields import MinModeMethod, Accelerant

print(capnp_path())
print(job_result_capnp_path())
print(MainConfig().job)
print(DimerSpec(method=MinModeMethod.improved, accelerant=Accelerant.gp).core_kwargs())
```

eon-akmc / Sphinx still use:

```python
from eon.schema import MainConfig, Metatomic  # re-export
```

pyeonclient still uses:

```python
from pyeonclient.models import DimerSpec, NebSpec  # re-export
```

## Layout

```text
src/eon_schema/
  ssot/          # L0 vendored params Cap'n Proto + catalog
  jobs/          # L0 JobRequest/JobResult schema + results.dat adapters
  config/        # L1 job-config models (models.py)
  fields/        # enums
  api/           # L2 DimerSpec, NebSpec
  _deps.py
```

## INI helpers (config write without eon-akmc)

```python
from eon_schema.config import (
    MainConfig,
    PotentialConfig,
    DimerConfig,
    write_models_ini,
    write_ini,
    hydrate_ini,
    unknown_ini_keys,
    defaults_from_catalog,
)

# L1 models → config.ini (Dimer L1 prefixes mapped to INI keys)
write_models_ini(
    "config.ini",
    MainConfig(job="saddle_search"),
    PotentialConfig(potential="lj"),
    DimerConfig(dimer_improved=True),
    extra={"Debug": {"write_movies": True}},
    validate=True,
)

# Or raw sections (rgpycrumbs-style) with L0 hydrate + key check
user = {"Main": {"job": "minimization", "temperature": 400}}
cfg = hydrate_ini(user, base_sections=["Main", "Optimizer"])
assert not unknown_ini_keys(cfg)
write_ini("config.ini", cfg)

print(defaults_from_catalog("Main")["job"])
```

Downstream: **rgpycrumbs** should depend on ``eon-schema`` for
``write_eon_config`` / seed_dimers / MLflow log_params instead of full eon-akmc.
