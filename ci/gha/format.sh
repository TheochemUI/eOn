#!/usr/bin/env bash
# Format Nickel sources under ci/gha/ (optional; prek may run nickel-format separately).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
if ! command -v nickel >/dev/null 2>&1; then
  echo "error: nickel not on PATH (use: pixi r -e cigen nickel-fmt)" >&2
  exit 1
fi
shopt -s nullglob
files=(ci/gha/*.ncl ci/gha/lib/*.ncl)
if ((${#files[@]} == 0)); then
  echo "no .ncl files under ci/gha/"
  exit 0
fi
nickel format "${files[@]}"
echo "OK: formatted ${#files[@]} nickel files"
