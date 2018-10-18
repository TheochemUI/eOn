#!/bin/bash
set -e

export EON_NUMBER_OF_CLIENTS=3
cores_per_pot=4

cores_per_node=8

eon_path=../eon
pot_path=vasp_mpmd

export EON_SERVER_PATH=${eon_path}/akmc.py
client_path=${eon_path}/client/client_mpi

client_ranks=$(expr $EON_NUMBER_OF_CLIENTS + 1)
pot_ranks=$(expr $cores_per_pot \* $EON_NUMBER_OF_CLIENTS)

command="mpirun -n $pot_ranks -wdir vasp $pot_path : -n $client_ranks $client_path"

#if no arguments
if [ -z $1 ]; then
    $command
else
    total_ranks=$(expr $pot_ranks + $client_ranks)
    if [ "$1" = "-n" ]; then
        echo $(expr $cores_per_node \* \( \( $total_ranks - 1 \) / $cores_per_node \) + $cores_per_node)
    fi

    if [ "$1" = "--dry-run" ]; then
        echo mpirun command: $command
        echo total nodes: $(expr $total_ranks / $cores_per_node)
        echo total ranks: $total_ranks
    fi
fi
