#!/usr/bin/env bash
# Export nickel-backed release workflows (SPECS table; rgpot ci/gha/gen.sh style).
# Entry: ./ci/gha/gen.sh   OR   pixi r -e cigen gen-gha
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

mapfile -t SPECS <<'EOF'
release.ncl|.github/workflows/release.yml
release_prepare.ncl|.github/workflows/release-prepare.yml
towncrier_check.ncl|.github/workflows/towncrier.yml
EOF

if ! command -v nickel >/dev/null 2>&1; then
  echo "error: nickel not on PATH (install nickel-lang or: pixi r -e cigen gen-gha)" >&2
  exit 1
fi

for line in "${SPECS[@]}"; do
  src="${line%%|*}"
  dst="${line##*|}"
  echo "nickel export $src -> $dst"
  nickel export --format yaml "ci/gha/$src" -o "$dst"
done
echo "OK: ${#SPECS[@]} workflows exported (hand ci_*.yml not in scope)"
