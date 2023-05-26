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

  this->currentStatus += HelperFunctionsTest::test_random();

  this->currentStatus += HelperFunctionsTest::test_randomDouble();

  this->currentStatus += HelperFunctionsTest::test_randomInt();

  this->currentStatus += HelperFunctionsTest::test_gaussRandom();

  this->currentStatus += HelperFunctionsTest::test_dot();

  this->currentStatus += HelperFunctionsTest::test_length();
  /*
          this->currentStatus+=HelperFunctionsTest::test_add();

          this->currentStatus+=HelperFunctionsTest::test_subtract();

          this->currentStatus+=HelperFunctionsTest::test_multiplyScalar();

          this->currentStatus+=HelperFunctionsTest::test_divideScalar();

          this->currentStatus+=HelperFunctionsTest::test_copyRightIntoLeft();

          this->currentStatus+=HelperFunctionsTest::test_nomralize();

          this->currentStatus+=HelperFunctionsTest::test_makeOrthogonal();

          this->currentStatus+=HelperFunctionsTest::test_makeProjection();

          this->currentStatus+=HelperFunctionsTest::test_rotationMatch();

          this->currentStatus+=HelperFunctionsTest::test_maxAtomMotion();

          this->currentStatus+=HelperFunctionsTest::test_maxAtomMotionV();

          this->currentStatus+=HelperFunctionsTest::test_numAtomsMoved();

          this->currentStatus+=HelperFunctionsTest::test_maxAtomMotionApplied();

          this->currentStatus+=HelperFunctionsTest::test_maxAtomMotionAppliedV();

          this->currentStatus+=HelperFunctionsTest::test_maxMotionApplied();

          this->currentStatus+=HelperFunctionsTest::test_maxMotionAppliedV();

          this->currentStatus+=HelperFunctionsTest::test_getTime();

          this->currentStatus+=HelperFunctionsTest::test_existsFile();

          this->currentStatus+=HelperFunctionsTest::test_getRelevantFile();

          this->currentStatus+=HelperFunctionsTest::test_loadMasses();

          this->currentStatus+=HelperFunctionsTest::test_loadMode();

          this->currentStatus+=HelperFunctionsTest::test_saveMode();

          this->currentStatus+=HelperFunctionsTest::test_split_string_int();

          this->currentStatus+=HelperFunctionsTest::test_identical();

          this->currentStatus+=HelperFunctionsTest::test_sortedR();

          this->currentStatus+=HelperFunctionsTest::test_pushApart();
  */
  std::string testName("HelperFunctionsTest");
  this->saveTest(testName);

  return this->currentStatus;
}

int HelperFunctionsTest::test_random() {
  long seed = 15;

  double random1 = helper_functions::random();

  std::string output;

  if (std::isnan(random1) || !std::isfinite(random1)) {
    output = "FAIL: helper_functions::random() couldn't initialize a double";
    this->pushToResults(output);
    return 1;

  } else if (random1 < 0 || random1 > 1) {

    output =
        "FAIL: helper_functions::random() gave a number out of range of (0,1)";
    this->pushToResults(output);
    return 1;
  }

  double random2 = helper_functions::random(seed);

  if (std::isnan(random2) || !std::isfinite(random2)) {
    output = "FAIL: helper_functions::random(seed) couldn't initialize a "
             "double using a seed argument";
    this->pushToResults(output);
    return 1;

  } else if (random2 < 0 || random2 > 1) {

    output = "FAIL: helper_functions::random(seed) gave a number out of range "
             "of (0,1) using a seed argument";
    this->pushToResults(output);
    return 1;
  }

  if (random1 == random2) {

    output = "FAIL: helper_functions::random(seed) does not give unique value "
             "when a seed is passed in";
    this->pushToResults(output);
    return 1;
  }

  output =
      "PASS: helper_functions::random() works as specified by documentation";
  this->pushToResults(output);
  return 0;
}

int HelperFunctionsTest::test_randomDouble() {
  int maxInt = 5;

  long maxLong = 10;

  double maxDouble = 15;

  double random1 = helper_functions::randomDouble();

  std::string output;

  if (std::isnan(random1) || !std::isfinite(random1)) {
    output =
        "FAIL: helper_functions::randomDouble() couldn't initialize a double";
    this->pushToResults(output);
    return 1;

  } else if (random1 < 0 || random1 > 1) {

    output = "FAIL: helper_functions::randomDouble() gave a number out of "
             "range of (0,1)";
    this->pushToResults(output);
    return 1;
  }

  double random2 = helper_functions::randomDouble(maxInt);

  if (std::isnan(random2) || !std::isfinite(random2)) {
    output = "FAIL: helper_functions::randomDouble(int) couldn't initialize a "
             "double";
    this->pushToResults(output);
    return 1;

  } else if (random2 < 0 || random2 > maxInt) {

    output = "FAIL: helper_functions::randomDouble(int) gave a number out of "
             "range of (0,int)";
    this->pushToResults(output);
    return 1;
  }

  double random3 = helper_functions::randomDouble(maxLong);

  if (std::isnan(random3) || !std::isfinite(random3)) {
    output = "FAIL: helper_functions::randomDouble(long) couldn't initialize a "
             "double";
    this->pushToResults(output);
    return 1;

  } else if (random3 < 0 || random3 > maxLong) {

    output = "FAIL: helper_functions::randomDouble(long) gave a number out of "
             "range of (0,long)";
    this->pushToResults(output);
    return 1;
  }

  double random4 = helper_functions::randomDouble(maxDouble);

  if (std::isnan(random4) || !std::isfinite(random4)) {
    output = "FAIL: helper_functions::randomDouble(double) couldn't initialize "
             "a double";
    this->pushToResults(output);
    return 1;

  } else if (random4 < 0 || random4 > maxDouble) {

    output = "FAIL: helper_functions::randomDouble(double) gave a number out "
             "of range of (0,double)";
    this->pushToResults(output);
    return 1;
  }

  output = "PASS: helper_functions::randomDouble() works as specified by "
           "documentation";
  this->pushToResults(output);
  return 0;
}

int HelperFunctionsTest::test_randomInt() {

  int lower = 1, upper = 4;

  long random1 = helper_functions::randomInt(lower, upper);

  std::string output;

  if (std::isnan(random1) || !std::isfinite(random1)) {
    output = "FAIL: helper_functions::randomInt(lower,upper) couldn't "
             "initialize a long";
    this->pushToResults(output);
    return 1;

  } else if (random1 < lower || random1 > upper) {

    output = "FAIL: helper_functions::randomInt(lower,upper) gave a number out "
             "of range of (lower,upper)";
    this->pushToResults(output);
    return 1;
  }

  output = "PASS: helper_functions::randomInt(lower,upper) works as specified "
           "by documentation";
  this->pushToResults(output);
  return 0;
}

int HelperFunctionsTest::test_gaussRandom() {

  double avg = 1, std = 0.1;

  double lb = avg - 4 * std, ub = avg + 4 * std;

  double random1 = helper_functions::gaussRandom(avg, std);

  std::string output;

  if (std::isnan(random1) || !std::isfinite(random1)) {
    output = "FAIL: helper_functions::gaussRandom(avg,std) couldn't initialize "
             "a double";
    this->pushToResults(output);
    return 1;

  } else if (random1 < lb || random1 > ub) {

    output = "FAIL: helper_functions::gaussRandom(avg,std) returned a double 4 "
             "standard deviations out; this is probably an error";
    this->pushToResults(output);
    return 1;
  }

  output = "PASS: helper_functions::gaussRandom(avg,std) works as specified by "
           "documentation";
  this->pushToResults(output);
  return 0;
}

int HelperFunctionsTest::test_dot() {

  double x[2] = {1, 0}, y[2] = {0, 1};

  long size1 = 2;

  double w[3] = {1, 2, 3}, z[3] = {3, 2, 1};

  long size2 = 3;

  const double *xp, *yp, *wp, *zp;

  xp = x;

  yp = y;

  wp = w;

  zp = z;

  double orthoDot = helper_functions::dot(xp, yp, size1);

  std::string output;

  if (std::isnan(orthoDot) || !std::isfinite(orthoDot)) {
    output = "FAIL: helper_functions::dot(v1, v2, size) couldn't initialize a "
             "double";
    this->pushToResults(output);
    return 1;

  } else if (orthoDot != 0) {

    output = "FAIL: helper_functions::dot(v1, v2, size) dotted two orthogonal "
             "vectors and didn't return zero";
    this->pushToResults(output);
    return 1;
  }

  double regDot1 = helper_functions::dot(wp, zp, size2);

  if (std::isnan(regDot1) || !std::isfinite(regDot1)) {
    output = "FAIL: helper_functions::dot(v1, v2, size) couldn't initialize a "
             "double";
    this->pushToResults(output);
    return 1;

  } else if (regDot1 != 10) {

    output = "FAIL: helper_functions::dot(v1, v2, size) dotted two vectors and "
             "returned the wrong value";
    this->pushToResults(output);
    return 1;
  }

  double regDot2 = helper_functions::dot(zp, wp, size2);

  if (std::isnan(regDot2) || !std::isfinite(regDot2)) {
    output = "FAIL: helper_functions::dot(v1, v2, size) couldn't initialize a "
             "double";
    this->pushToResults(output);
    return 1;

  } else if (regDot2 != 10) {

    output = "FAIL: helper_functions::dot(v1, v2, size) dotted two vectors and "
             "returned the wrong value";
    this->pushToResults(output);
    return 1;
  }

  if (regDot1 != regDot2) {

    output = "FAIL: helper_functions::dot(v1, v2, size) did not equal "
             "helper_functions::dot(v2, v1, size)";
    this->pushToResults(output);
    return 1;
  }

  output = "PASS: helper_functions::dot(v1,v2,size) works as specified by "
           "documentation";
  this->pushToResults(output);
  return 0;
}

int HelperFunctionsTest::test_length() {

  double x[2] = {3, 4}, y[4] = {1, 2, 3, 4};

  long size1 = 2, size2 = 4;

  const double *xp, *yp;

  xp = x;

  yp = y;

  double norm1 = helper_functions::length(xp, size1);

  std::string output;

  if (std::isnan(norm1) || !std::isfinite(norm1)) {
    output =
        "FAIL: helper_functions::length(v1, size) couldn't initialize a double";
    this->pushToResults(output);
    return 1;

  } else if (norm1 != 5) {

    output = "FAIL: helper_functions::length(v1, size) took the norm of a "
             "vector and returned the wrong value";
    this->pushToResults(output);
    return 1;
  }

  double normY = sqrt(30);

  double norm2 = helper_functions::length(yp, size2);

  if (std::isnan(norm2) || !std::isfinite(norm2)) {
    output =
        "FAIL: helper_functions::length(v1, size) couldn't initialize a double";
    this->pushToResults(output);
    return 1;

  } else if (norm2 != normY) {

    output = "FAIL: helper_functions::length(v1, size) took the norm of a "
             "vector and returned the wrong value";
    this->pushToResults(output);
    return 1;
  }

  output = "PASS: helper_functions::length(v1,size) works as specified by "
           "documentation";
  this->pushToResults(output);
  return 0;
}
/*
int HelperFunctionsTest::test_add()
{

        double x[3] = {1,2,3}, y[3] = {1,2,3};

        long size = 3;

        const double *xp, *yp;

        xp = x;

        yp = y;

        double result[3], *res, shouldRes[3] = {2,4,6};

        helper_functions::add(res,xp,yp,size);

        for(int i=0;i<3;i++)
        {
                result[i] = res[i];
        }

        std::string output;

        if(result != shouldRes) {

                output = "FAIL: helper_functions::add(result,xp,yp,size)
returned an incorrect result"; this->pushToResults(output); return 1;

        }

        output = "PASS: helper_functions::add(result,v1,v2,size) works as
specified by documentation"; this->pushToResults(output); return 0;


}
*/
