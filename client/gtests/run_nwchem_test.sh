#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
# Arguments passed from Meson
TEST_EXECUTABLE=$1
NWCHEM_EXECUTABLE=$2
NWCHEM_INPUT_FILE=$3
TEST_WORKDIR=$4

# Expected results for verification
EXPECTED_ENERGY="-7078.991829615012"
EXPECTED_MAX_FORCE="4.260180846246"
TOLERANCE="1e-5"

# --- Cleanup Function ---
# This function will run when the script exits, ensuring the server is killed
# even if the test fails or is interrupted.
cleanup() {
    echo "--- Cleaning up test processes ---"
    # Check if the server process ID exists and kill it
    if kill -0 "$SERVER_PID" 2>/dev/null; then
        kill "$SERVER_PID"
    fi
}
trap cleanup EXIT

# --- Test Execution ---

# 1. Navigate to the correct working directory
cd "$TEST_WORKDIR"
echo "--- Running test in: $(pwd) ---"
rm -f results.dat

# 2. Start the EON test server in the background
echo "--- Starting EON test server in background ---"
"$TEST_EXECUTABLE" &
SERVER_PID=$! # Capture the Process ID of the server

# 3. Wait for a moment to ensure the server is up and listening
#  A more robust solution would poll for the socket file, but a short sleep is usually enough.
sleep 1

# 4. Run the NWChem client in the foreground
#    It will connect to the server, do its work, and exit.
echo "--- Starting NWChem client ---"
mpirun -np 4 "$NWCHEM_EXECUTABLE" "$NWCHEM_INPUT_FILE"

# 5. Wait for the server process to finish and get its exit code
echo "--- Waiting for EON test server to exit ---"
wait "$SERVER_PID"
TEST_EXIT_CODE=$?

CALCULATED_ENERGY=$(grep 'Energy' results.dat | awk '{print $1}')
CALCULATED_MAX_FORCE=$(grep 'Max_Force' results.dat | awk '{print $1}')

# 6. Extract and compare
echo "Expected Energy: $EXPECTED_ENERGY"
echo "Calculated Energy: $CALCULATED_ENERGY"
echo "Expected Max Force: $EXPECTED_MAX_FORCE"
echo "Calculated Max Force: $CALCULATED_MAX_FORCE"

# 8. Compare floating point numbers using awk
echo "--- Comparing Energy (tolerance: $TOLERANCE) ---"
awk -v e1="$EXPECTED_ENERGY" -v e2="$CALCULATED_ENERGY" -v tol="$TOLERANCE" \
'BEGIN { diff = e1 - e2; if (diff < 0) diff = -diff; if (diff > tol) exit 1 }'

echo "--- Comparing Max Force (tolerance: $TOLERANCE) ---"
awk -v f1="$EXPECTED_MAX_FORCE" -v f2="$CALCULATED_MAX_FORCE" -v tol="$TOLERANCE" \
'BEGIN { diff = f1 - f2; if (diff < 0) diff = -diff; if (diff > tol) exit 1 }'


# 7. Exit the script with the server's exit code
echo "--- Test finished with exit code: $TEST_EXIT_CODE ---"
exit "$TEST_EXIT_CODE"
