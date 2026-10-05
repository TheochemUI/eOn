# eOn Cap'n Proto parameter SSoT

**Authoring home:** `schema/eon_params.capnp`

Edit field names, types, defaults, and wire ordinals here first. A struct
tagged `# project: ini,json` is projected into the INI reader, the JSON codec,
and the default assignment. Untagged structs stay handwritten.

```bash
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh   # capnp sources into the split
```

Codegen writes:

- `schema/eon_params_catalog.json`
- `eon/_params_ssot_catalog.py`
- `include/eon/generated/ParametersSSOTDefaults.h`
- `include/eon/generated/ParametersSSOTFieldIndex.inc`
- `include/eon/generated/ParametersSSOTIni.inc`
- `include/eon/generated/ParametersSSOTJson.inc`
- `include/eon/generated/ParametersSSOTApply.inc`
- `packages/eon-schema/src/eon_schema/ssot/eon_params_catalog.json`

The sync script copies `eon_params.capnp` and `eon_job_result.capnp` into the
split package. Codegen already writes the vendored catalog JSON.

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
- `JobRequest` / `JobResult` — status codes, barriers, force-call buckets, optional geometries, optimizer provenance, append-only `landfoldArtifacts`

Python: `eon_schema.jobs` (path helper + `results.dat` adapters).
Landfold trajectories use `trajectory_manifest`: exact per-frame geometry
digests and an rgpot potential identity, not `results.dat` lines.
`landfold_consume` checks those digests and returns embedding coordinates
and, when every frame has an energy, FES samples.

`.con` / `results.dat` remain **durable adapters**, not the in-process primary API.
