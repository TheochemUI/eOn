#!/bin/bash
# run_rgpot_mpi.sh EXE MODE
# MODE is fault or params. Exit 77 when cpmdc or mpirun is absent.
# A hang in MPI_Comm_split is a failure. A Segmentation fault on the
# way out is a failure. Rank 0's text has to carry the engine error.
set -u
EXE=$1
MODE=$2
if [ -z "${EXE}" ] || [ -z "${MODE}" ]; then
  echo "usage: run_rgpot_mpi.sh EXE fault|params" >&2
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
cd "$tmp"
set +e
timeout 25 mpirun -np 2 --bind-to none --map-by :OVERSUBSCRIBE \
  -x TMPDIR -x PMIX_MCA_psec -x OMPI_ALLOW_RUN_AS_ROOT \
  -x OMPI_ALLOW_RUN_AS_ROOT_CONFIRM -x CPMDC_LIBRARY \
  -x RGPOT_CPMDC_ENGINE -x RGPOT_CPMD_ENGINE -x LD_LIBRARY_PATH -x PATH \
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
  grep -q "rank=0 fault .*engine-rank1" stderr || {
    echo "rank 0 did not print the engine error"
    exit 1
  }
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
echo "unknown mode $MODE"
exit 2
