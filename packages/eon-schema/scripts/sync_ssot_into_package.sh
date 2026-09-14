#!/usr/bin/env bash
# Copy monorepo SSoT (schema/) into the eon-schema split package so the
# PyPI wheel ships a self-contained Cap'n Proto + catalog copy.
# Authoring is always schema/eon_params.capnp and schema/eon_job_result.capnp.
set -euo pipefail
PKG="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$(cd "$PKG/../.." && pwd)"
SRC="$ROOT/schema"
DEST="$PKG/src/eon_schema/ssot"
JOBS_DEST="$PKG/src/eon_schema/jobs"
cp -a "$SRC/eon_params.capnp" "$DEST/"
cp -a "$SRC/eon_params_catalog.json" "$DEST/"
mkdir -p "$JOBS_DEST"
cp -a "$SRC/eon_job_result.capnp" "$JOBS_DEST/"
cat > "$DEST/README.md" << 'NOTE'
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
NOTE
echo "synced $SRC → $DEST + $JOBS_DEST (vendored into eon-schema split)"
