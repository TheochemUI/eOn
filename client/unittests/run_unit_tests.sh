#!/bin/bash

#PLEASE DO NOT RUN THIS SCRIPT DIRECTLY
#Please run this through running the target 'make check' in the previous directory

echo "Starting unit-testing"

FILE="Test.out"

if [ ! -f "$FILE" ]; then
	echo "Couldn't find main unit-testing executable: exiting"
	exit 1
fi

VERSION=`svn info http://theory.cm.utexas.edu/svn/eon | awk '/^Revision:/{print $2}'`
BUILDDATE=`date`
BUILDHOST=`hostname -d`
BUILDUSER=`whoami`

OUTPUTDIR=`date --rfc-3339="seconds" | awk -F" " '{ print $1"_"$2 }' | cut -d"-" -f1-3`

echo "Creating unit-testing detailed output directory"

mkdir $OUTPUTDIR

if [ ! -d $OUTPUTDIR ]; then
	echo "Couldn't create a directory for the run: exiting"
	exit 1
fi

echo "Entering unit-testing directiory"

pushd ./$OUTPUTDIR > /dev/null

echo "Creating log file for the run"

LOG="run_log.txt"

touch $LOG

if [ ! -f "$LOG" ]; then
	echo "Couldn't create a log for the run: cleaning up then exiting"
	popd > /dev/null
	rm -rf ./$OUTPUTDIR > /dev/null
	exit 1
fi

printf "Version: $VERSION\nBuild Date: $BUILDDATE\nBuild Host: $BUILDHOST\nBuild User: BUILDUSER\n" >> $LOG

echo "Begin running unit-tests"

echo "Running tests on functions in file HelperFunctions.cpp"

RUN_HELPERFUNCTIONS=`./../Test.out HelperFunctionTest`

if [ "$RUN_HELPERFUNCTIONS" -eq "-1" ]; then
        echo "There was an error trying to run this test..."
        OUTPUT_HELPERFUNCTIONS="ERROR: there was an error trying to test HelperFunctions.cpp"
        STATUS_HELPERFUNCTIONS=1
elif [ $((RUN_HELPERFUNCTIONS)) -eq 0 ]; then
	echo "All tests ended in success!"
	OUTPUT_HELPERFUNCTIONS="PASS: all functions in HelperFunctions.cpp ended in success"
	STATUS_HELPERFUNCTIONS=0
else
	if [ $((RUN_HELPERFUNCTIONS)) -eq 1 ]; then
		echo "$RUN_HELPERFUNCTIONS test failed..."
		OUTPUT_HELPERFUNCTIONS="FAIL: $RUN_HELPERFUNCTIONS function in HelperFunctions.cpp ended in failure"
	else
		echo "$RUN_HELPERFUNCTIONS tests failed..."
		OUTPUT_HELPERFUNCTIONS="FAIL: $RUN_HELPERFUNCTIONS functions in HelperFunctions.cpp ended in failure"
	fi
	STATUS_HELPERFUNCTIONS=1
fi

printf "\nUNIT TESTING RUN RESULTS\n---------------------------------------------------------------------\n"
printf "\nRUN RESULTS\n---------------------------------------------------------------------\n" >> $LOG
echo $OUTPUT_HELPERFUNCTIONS | tee -a $LOG
echo "---------------------------------------------------------------------" | tee -a $LOG

RUN_RESULTS=$((STATUS_HELPERFUNCTIONS))

if [ $((RUN_RESULTS)) -eq 0 ]; then
	echo "All unit-tests ended in success!" | tee -a $LOG
	echo "Look at the output files of the individual tests for more details" >> $LOG
else
	if [ $((RUN_RESULTS)) -eq 1 ]; then
		echo "$RUN_RESULTS file either erred out of testing or had some failures" | tee -a $LOG
		echo "Look at the output file of the failed test for more details" >> $LOG
	else
		echo "$RUN_RESULTS files either erred out of testing or had some failures" | tee -a $LOG
		echo "Look at the output files of the failed tests for more details" >> $LOG
	fi
fi

echo "Look at the file $LOG in the directory ./unittests/$OUTPUTDIR for more details"

popd > /dev/null
echo -e "\nexiting"
exit 0
