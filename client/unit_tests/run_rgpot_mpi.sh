#!/bin/bash
# run_rgpot_mpi.sh EXE MODE
# MODE is fault, params, abort, single or uneven. Exit 77 when cpmdc or
# mpirun is absent.
# A hang in MPI_Comm_split is a failure. A Segmentation fault on the
# way out is a failure. Rank 0's text has to carry the engine error.
set -u
EXE=$1
MODE=$2
if [ -z "${EXE}" ] || [ -z "${MODE}" ]; then
  echo "usage: run_rgpot_mpi.sh EXE fault|params|abort|single|uneven" >&2
  exit 2
fi
if ! command -v mpirun >/dev/null 2>&1; then
  echo "mpirun not on PATH"
  exit 77
fi
if [ -z "${CPMDC_LIBRARY:-}${RGPOT_CPMDC_ENGINE:-}${RGPOT_CPMD_ENGINE:-}" ]; then
  echo "no libcpmdc in the environment"
  exit 77
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
case "$tmp" in
  /*) ;;
  *) tmp=$PWD/$tmp ;;
esac
export TMPDIR=$tmp
export PMIX_MCA_psec=none
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_rmaps_base_oversubscribe=1
# The compiler records the prefix libdir ahead of $ORIGIN, and that
# directory is also the install prefix. Preload the librgpot_pot.so
# next to this executable so the ranks run the build just linked.
exe_dir=$(CDPATH= cd -- "$(dirname "$EXE")" && pwd)
preload_arg=()
potso="$exe_dir/potentials/Rgpot/librgpot_pot.so"
# Preload only when the loader would pick another librgpot_pot.so. A
# preload the binary does not need corrupts the heap under some runtime
# loaders (EESSI 2026.06) before main runs.
resolved=$(ldd "$EXE" 2>/dev/null | awk '/librgpot_pot\.so/ {print $3; exit}')
if [ -f "$potso" ] && [ -n "$resolved" ] &&
  [ "$(readlink -f "$resolved")" != "$(readlink -f "$potso")" ]; then
  if [ -n "${LD_PRELOAD:-}" ]; then
    preload_arg=(-x "LD_PRELOAD=${potso}:${LD_PRELOAD}")
  else
    preload_arg=(-x "LD_PRELOAD=${potso}")
  fi
fi
cd "$tmp"
set +e
timeout 25 mpirun -np 2 --bind-to none --map-by :OVERSUBSCRIBE \
  -x TMPDIR -x PMIX_MCA_psec -x OMPI_ALLOW_RUN_AS_ROOT \
  -x OMPI_ALLOW_RUN_AS_ROOT_CONFIRM -x CPMDC_LIBRARY \
  -x RGPOT_CPMDC_ENGINE -x RGPOT_CPMD_ENGINE -x LD_LIBRARY_PATH -x PATH \
  "${preload_arg[@]}" \
  "$EXE" "$MODE" >stdout 2>stderr
rc=$?
set -e
cat stderr stdout
if [ "$rc" -eq 77 ]; then
  exit 77
fi
if [ "$rc" -eq 124 ]; then
  echo "mpirun still running after 25s"
  exit 1
fi
if grep -q "Segmentation fault" stderr stdout; then
  echo "segmentation fault on the way out"
  exit 1
fi
if [ "$rc" -gt 128 ]; then
  echo "mpirun died on signal ($rc)"
  exit 1
fi
if [ "$MODE" = "fault" ]; then
  grep -q "rank=0 fault owner=1 energy=.*engine-rank1" stderr || {
    echo "rank 0 did not print the owner, the shared energy, and the engine error"
    exit 1
  }
  e0=$(sed -n 's/.*rank=0 fault owner=1 energy=\([^ ]*\).*/\1/p' stderr | head -1)
  e1=$(sed -n 's/.*rank=1 fault owner=1 energy=\([^ ]*\).*/\1/p' stderr | head -1)
  if [ -z "$e0" ] || [ "$e0" != "$e1" ]; then
    echo "ranks did not share one energy (rank0=$e0 rank1=$e1)"
    exit 1
  fi
  [ "$rc" -ne 0 ] || {
    echo "mpirun exited 0"
    exit 1
  }
  exit 0
fi
if [ "$MODE" = "params" ]; then
  grep -q "rank=0 params-agreed .*cannot open params_path" stderr || {
    echo "rank 0 did not report the unreadable params_path"
    exit 1
  }
  [ "$rc" -ne 0 ] || {
    echo "mpirun exited 0"
    exit 1
  }
  exit 0
fi
if [ "$MODE" = "abort" ]; then
  grep -q "rank=0 abort requested" stderr || {
    echo "rank 0 did not reach the abort request"
    exit 1
  }
  [ "$rc" -ne 0 ] || {
    echo "mpirun exited 0"
    exit 1
  }
  exit 0
fi
if [ "$MODE" = "single" ] || [ "$MODE" = "uneven" ]; then
  grep -q "rank=0 $MODE done" stderr || {
    echo "rank 0 did not finish the $MODE call"
    exit 1
  }
  [ "$rc" -eq 0 ] || {
    echo "mpirun exited $rc"
    exit 1
  }
  exit 0
fi
echo "unknown mode $MODE"
exit 2
