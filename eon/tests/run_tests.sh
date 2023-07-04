#!/bin/bash

tests="akmc_eam_al_trimer min_morse_pt_island bh_lj13"

echo "Running tests"
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
exit $failed
