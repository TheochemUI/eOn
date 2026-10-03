#!/bin/bash
# check_no_mpi_closure.sh EONCLIENT
# A plain eonclient must not map an MPI library or its transports at load:
# their constructors (UCX memory hooks, PSM2) cost 0.2 to 0.4 s per start.
# Calculator groups load them on demand. Exit 77 without ldd.
set -u
EXE=$1
if ! command -v ldd >/dev/null 2>&1; then
  echo "ldd not on PATH"
  exit 77
fi
closure=$(ldd "$EXE") || { echo "ldd failed on $EXE"; exit 1; }
hits=$(printf '%s\n' "$closure" | grep -E 'lib(mpi|mpi_cxx|open-pal|open-rte|pmix|ucp|ucs|ucm|uct|psm2|fabric)\.so' || true)
if [ -n "$hits" ]; then
  echo "eonclient maps MPI at load:"
  printf '%s\n' "$hits"
  exit 1
fi
echo "no MPI library in the load closure of $EXE"
