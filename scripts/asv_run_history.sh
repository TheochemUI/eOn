#!/usr/bin/env bash
# Run ASV for the last N commits against an existing env (eonclient already on PATH).
set -euo pipefail
N="${1:-5}"
EXTRA="${2:-}"
PY="$(which python)"
echo "python=$PY"
echo "eonclient=$(command -v eonclient || true)"
mapfile -t SHAS < <(git rev-list --max-count="$N" HEAD)
# oldest first
for ((i=${#SHAS[@]}-1; i>=0; i--)); do
  sha="${SHAS[i]}"
  echo "=== asv run $EXTRA --set-commit-hash $sha ==="
  # shellcheck disable=SC2086
  asv run -E "existing:${PY}" --set-commit-hash "${sha}" --record-samples ${EXTRA}
done
