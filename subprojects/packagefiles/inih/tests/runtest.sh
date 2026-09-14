#!/bin/sh

set -eu

# Strip CR so MSVC text-mode stdout (CRLF) matches LF baselines.
cd "$(dirname "${1}")"
"$2" | tr -d '\r' | diff "${1}" -
