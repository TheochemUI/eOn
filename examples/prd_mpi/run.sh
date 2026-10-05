#!/bin/bash
export EON_NUMBER_OF_CLIENTS=7
export EON_SERVER_PATH="$PWD/server.py"
mpirun -n 8 eonclient
