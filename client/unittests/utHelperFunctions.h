//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef UT_HELPER_FUNCTIONS_H
#define UT_HELPER_FUNCTIONS_H

#include "utObjectiveTest.h"

class HelperFunctionsTest : public ObjectiveTest {

public:
  HelperFunctionsTest(){};
  ~HelperFunctionsTest(){};
  int test();

private:
  int test_random();
  int test_randomDouble();
  int test_randomInt();
  int test_gaussRandom();
  int test_dot();
  int test_length();
  int test_add();
  int test_subtract();
  int test_multiplyScalar();
  int test_divideScalar();
  int test_copyRightIntoLeft();
  int test_normalize();
  int test_makeOrthogonal();
  int test_makeProjection();
  int test_rotationMatch();
  int test_maxAtomMotion();
  int test_maxAtomMotionV();
  int test_numAtomsMoved();
  int test_maxAtomMotionApplied();
  int test_maxAtomMotionAppliedV();
  int test_maxMotionApplied();
  int test_maxMotionAppliedV();
  int test_getTime();
  int test_existsFile();
  int test_getRelevantFile();
  int test_loadMasses();
  int test_loadMode();
  int test_saveMode();
  int test_split_string_int();
  int test_identical();
  int test_sortedR();
  int test_pushApart();
};

#endif
