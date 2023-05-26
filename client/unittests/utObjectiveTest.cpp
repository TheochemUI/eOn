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
#include <fstream>
#include <iostream>

#include "utHelperFunctions.h"

const std::string ObjectiveTest::HELPER_FUNCTIONS = "HelperFunctionTest";

ObjectiveTest *ObjectiveTest::getObjectiveTest(std::string testType) {
  ObjectiveTest *objectiveTest = NULL;
  if (testType == ObjectiveTest::HELPER_FUNCTIONS) {
    objectiveTest = new HelperFunctionsTest();
  }

  return objectiveTest;
};

void ObjectiveTest::pushToResults(std::string update) {

  this->testResults.push_back(update);
}

void ObjectiveTest::saveTest(std::string testName) {

  std::string resultsFileName(testName);
  resultsFileName += ".txt";
  std::ofstream resultsFile(resultsFileName);

  resultsFile << std::endl;

  if (currentStatus == 0) {

    resultsFile << "PASSED: all " << testName << "'s ended in success"
                << std::endl;

  } else {

    resultsFile << "FAILED: some tests in " << testName << " ended in failure"
                << std::endl;
  }

  resultsFile << std::endl
              << "-------------------------------------------------------------"
                 "--------";

  resultsFile << std::endl
              << std::endl
              << "Detailed output: " << std::endl
              << std::endl;

  for (const auto lineOut : this->testResults)
    resultsFile << lineOut << std::endl;

  resultsFile << std::endl;

  resultsFile.close();
};
