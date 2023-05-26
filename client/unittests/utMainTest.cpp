//-----------------------------------------------------------------------------------
//// eOn is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// A copy of the GNU General Public License is available at
//// http://www.gnu.org/licenses/
////-----------------------------------------------------------------------------------

#include "utObjectiveTest.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>

int main(int argc, char *argv[]) {
  std::string testType;
  if (argc != 2) {
    std::cout << "-1" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    testType = argv[1];
  }

  ObjectiveTest *unitTest = ObjectiveTest::getObjectiveTest(testType);
  if (unitTest == NULL) {
    std::cout << "-1" << std::endl;
    exit(EXIT_FAILURE);
  }

  int testResults = 0;
  try {
    testResults = unitTest->test();
  } catch (int error) {

    std::string errorOut("[ERROR] job exited on error ");
    std::ostringstream errorCode;
    errorCode << error;
    errorOut += errorCode.str();

    unitTest->pushToResults(errorOut);
    unitTest->saveTest(testType);
    exit(EXIT_FAILURE);
  }

  std::cout << testResults << std::endl;

  return 0;
}
