#!/bin/bash

tests="akmc_eam_al_trimer min_morse_pt_island potentials"

echo "Running tests"
echo
failed=0
for test in $tests
do
    cd $test
    ./test.py
    if [ $? != 0 ]; then
        failed=$(($failed + 1))
    fi 
    cd ..
done
echo 
echo "$failed tests failed"
