# Cap'n Proto L0 SSoT (vendored for this package)

**Authoring home:** monorepo `schema/eon_params.capnp` (params) and
`schema/eon_job_result.capnp` (job request/result envelope).

This directory is a **vendored copy** shipped in the PyPI `eon-schema` wheel
so installs work without the full monorepo. Do not edit field graphs here.

After editing monorepo `schema/`:

```bash
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh
```

Fat conda-forge releases use a full monorepo `git archive` tarball (not this package alone).
