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
