#!/bin/bash
export EON_NUMBER_OF_CLIENTS=2
export EON_SERVER_PATH="$PWD/server.py"
mpirun -n 3 eonclient
