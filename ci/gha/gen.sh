#!/usr/bin/env bash
# Regenerate release-related GHA workflows from Nickel (rgpot ci/gha style).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
command -v nickel >/dev/null 2>&1 || {
  echo "error: nickel not on PATH (install nickel-lang or pixi r -e cigen …)" >&2
  exit 1
}
nickel export --format yaml ci/gha/release.ncl -o .github/workflows/release.yml
nickel export --format yaml ci/gha/release_prepare.ncl -o .github/workflows/release-prepare.yml
nickel export --format yaml ci/gha/towncrier_check.ncl -o .github/workflows/towncrier.yml
echo "wrote:"
echo "  .github/workflows/release.yml"
echo "  .github/workflows/release-prepare.yml"
echo "  .github/workflows/towncrier.yml"
