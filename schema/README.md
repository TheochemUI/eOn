# eOn Cap'n Proto parameter SSoT

**Authoring home:** `schema/eon_params.capnp`

Edit field names, types, defaults, and wire ordinals here first, then:

```bash
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh   # if packages/eon-schema is present
```

That regenerates:

- `schema/eon_params_catalog.json`
- `eon/_params_ssot_catalog.py`
- `client/generated/ParametersSSOT*`
- vendored copy under `packages/eon-schema/src/eon_schema/ssot/` (PyPI split)

## Release shapes (monorepo)

| Shape | Artifact | Consumer |
|-------|----------|----------|
| **Fat tree** | `eon-vX.Y.Z.tar.xz` (`git archive` of this monorepo) | conda-forge `eon-feedstock`, EasyBuild, full source builds |
| **Splits** | PyPI `eon-schema`, `pyeonclient`, `eon-akmc`, … | pip/uv focused installs |

Layout under `schema/` or `packages/` can change; the release contract is that
the **fat** tarball is a complete monorepo archive the feedstock can build, while
splits publish independently with their own versions when useful.

## Job request / result envelope

**Authoring home:** `schema/eon_job_result.capnp`

Runtime control-plane types for kill-file-IPC (not parameter authoring):

- `Geometry` — flat `positions` (3N), `box` (9), Z, frozen mask
- `JobRequest` / `JobResult` — status codes, barriers, force-call buckets, optional geometries

Python: `eon_schema.jobs` (path helper + `results.dat` adapters).

`.con` / `results.dat` remain **durable adapters**, not the in-process primary API.
