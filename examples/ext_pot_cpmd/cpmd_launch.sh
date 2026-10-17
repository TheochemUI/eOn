#!/usr/bin/env bash
set -euo pipefail
trap 'status=$?; printf "CPMD launcher failed at line %s (exit %s)\n" "$LINENO" "$status" >&2' ERR
printf 'CPMD launcher: initializing\n' >&2
if [ -n "${CPMD_MODULES:-}" ]; then
    if ! type module >/dev/null 2>&1; then
        if [ -n "${CPMD_MODULE_INIT:-}" ]; then
            source "$CPMD_MODULE_INIT"
        else
            for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do
                if [ -r "$init" ]; then
                    source "$init"
                    break
                fi
            done
        fi
    fi
    type module >/dev/null 2>&1 || { echo 'module command unavailable: set CPMD_MODULE_INIT or use an already configured MPI with CPMD_MODULES empty' >&2; exit 2; }
    unset LD_LIBRARY_PATH LD_PRELOAD LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH
    printf 'CPMD launcher: loading requested compiler/MPI modules\n' >&2
    module purge -s
    read -r -a modules <<< "$CPMD_MODULES"
    module load "${modules[@]}"
fi
test -x "${1:?CPMD executable missing}" || { printf 'CPMD executable missing or not executable: %s\n' "$1" >&2; exit 2; }
printf 'CPMD launcher: executing %s\n' "$1" >&2
# shlex preserves quoted launcher arguments and paths containing spaces.
exec python3 -c 'import os, shlex, sys; a=shlex.split(os.environ["CPMD_MPI_LAUNCH"])+sys.argv[1:]; os.execvp(a[0], a)' "$@"
