/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/

#ifndef UT_HELPER_FUNCTIONS_H
#define UT_HELPER_FUNCTIONS_H

#include "utObjectiveTest.h"

class HelperFunctionsTest : public ObjectiveTest {

public:
  HelperFunctionsTest() {};
  ~HelperFunctionsTest() {};
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
