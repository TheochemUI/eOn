#!/bin/bash -l
set -euo pipefail
if [ -n "${CPMD_MODULES:-}" ]; then
    unset LD_LIBRARY_PATH LD_PRELOAD LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH
    module purge -s
    read -r -a modules <<< "$CPMD_MODULES"
    module load "${modules[@]}"
fi
# shlex preserves quoted launcher arguments and paths containing spaces.
exec python3 -c 'import os, shlex, sys; a=shlex.split(os.environ["CPMD_MPI_LAUNCH"])+sys.argv[1:]; os.execvp(a[0], a)' "$@"
