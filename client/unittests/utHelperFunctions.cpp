//-----------------------------------------------------------------------------------
//// eOn is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// A copy of the GNU General Public License is available at
//// http://www.gnu.org/licenses/
////-----------------------------------------------------------------------------------


#include "utHelperFunctions.h"

#include "../HelperFunctions.h"

#include <math.h>

int HelperFunctionsTest::test() {


	this->currentStatus+=HelperFunctionsTest::test_random();
/*
	this->currentStatus+=HelperFunctionsTest::test_randomDouble();

	currentStatus+=HelperFunctionsTest::test_randomInt();

	currentStatus+=HelperFunctionsTest::test_gaussRandom();

	currentStatus+=HelperFunctionsTest::test_dot();

	currentStatus+=HelperFunctionsTest::test_length();

	currentStatus+=HelperFunctionsTest::test_add();

	currentStatus+=HelperFunctionsTest::test_subtract();

	currentStatus+=HelperFunctionsTest::test_multiplyScalar();

	currentStatus+=HelperFunctionsTest::test_divideScalar();

	currentStatus+=HelperFunctionsTest::test_copyRightIntoLeft();

	currentStatus+=HelperFunctionsTest::test_nomralize();

	currentStatus+=HelperFunctionsTest::test_makeOrthogonal();

	currentStatus+=HelperFunctionsTest::test_makeProjection();

	currentStatus+=HelperFunctionsTest::test_rotationMatch();

	currentStatus+=HelperFunctionsTest::test_maxAtomMotion();

	currentStatus+=HelperFunctionsTest::test_maxAtomMotionV();

	currentStatus+=HelperFunctionsTest::test_numAtomsMoved();

	currentStatus+=HelperFunctionsTest::test_maxAtomMotionApplied();

	currentStatus+=HelperFunctionsTest::test_maxAtomMotionAppliedV();

	currentStatus+=HelperFunctionsTest::test_maxMotionApplied();

	currentStatus+=HelperFunctionsTest::test_maxMotionAppliedV();

	currentStatus+=HelperFunctionsTest::test_getTime();

	currentStatus+=HelperFunctionsTest::test_existsFile();

	currentStatus+=HelperFunctionsTest::test_getRelevantFile();

	currentStatus+=HelperFunctionsTest::test_loadMasses();

	currentStatus+=HelperFunctionsTest::test_loadMode();

	currentStatus+=HelperFunctionsTest::test_saveMode();

	currentStatus+=HelperFunctionsTest::test_split_string_int();

	currentStatus+=HelperFunctionsTest::test_identical();

	currentStatus+=HelperFunctionsTest::test_sortedR();

	currentStatus+=HelperFunctionsTest::test_pushApart();
*/
	std::string testName("HelperFunctionsTest");
	this->saveTest(testName);	

	return this->currentStatus;

}

int HelperFunctionsTest::test_random()
{
	long seed = 15;

	double random1 = helper_functions::random();

	std::string output;

	if(std::isnan(random1) || !std::isfinite(random1))
	{
		output = "FAIL: helper_functions::random() couldn't initialize a double";
		this->pushToResults(output);		
		return 1;
	
	} else if(random1 < 0 || random1 > 1) {

		output = "FAIL: helper_functions::random() gave a number out of range of (0,1)";
                this->pushToResults(output);
                return 1;

	}

	double random2 = helper_functions::random(seed);

        if(std::isnan(random2) || !std::isfinite(random2))
        {
                output = "FAIL: helper_functions::random(seed) couldn't initialize a double using a seed argument";
                this->pushToResults(output);
                return 1;

        } else if(random2 < 0 || random2 > 1) {

                output = "FAIL: helper_functions::random(seed) gave a number out of range of (0,1) using a seed argument";
                this->pushToResults(output);
                return 1;

        }


	if(random1 == random2) 
	{

		output = "FAIL: helper_functions::random(seed) does not give unique value when a seed is passed in";
		this->pushToResults(output);
		return 1;
	}

	output = "PASS: helper_functions::random() works as specified by documentation";
	this->pushToResults(output);
	return 0;

}

/*
int HelperFunctionsTest::test_randomDouble(
{
	
	return 1;

}
*/
