#!/usr/bin/env bash
# Regenerate nickel-exported workflows and fail if they differ from git.
# Entry: bash ci/gha/check-drift.sh   OR   pixi r -e cigen gha-drift
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

if ! command -v nickel >/dev/null 2>&1; then
  echo "error: nickel not on PATH (use: pixi r -e cigen gha-drift)" >&2
  exit 1
fi

bash ci/gha/gen.sh
TRACKED=(
  .github/workflows/release.yml
  .github/workflows/release-prepare.yml
  .github/workflows/towncrier.yml
)
if ! git diff --exit-code -- "${TRACKED[@]}" >/dev/null 2>&1; then
  echo "error: nickel-exported workflows out of date vs ci/gha/*.ncl" >&2
  echo "fix: edit Nickel under ci/gha/, then: pixi r -e cigen gen-gha && git add .github/workflows" >&2
  echo "--- git diff --stat (tracked nickel exports) ---" >&2
  git diff --stat -- "${TRACKED[@]}" >&2 || true
  exit 1
fi
echo "OK: nickel-exported workflows match ci/gha"
