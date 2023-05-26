//-----------------------------------------------------------------------------------
//// eOn is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// A copy of the GNU General Public License is available at
//// http://www.gnu.org/licenses/
////-----------------------------------------------------------------------------------

#ifndef OBJECTIVE_TEST_H
#define OBJECTIVE_TEST_H

#include <string>
#include <vector>

class ObjectiveTest {

public:
  virtual ~ObjectiveTest(){};
  virtual int test() = 0;

  std::vector<std::string> testResults;
  int currentStatus = 0;
  void pushToResults(std::string update);
  void saveTest(std::string testType);

  static const std::string HELPER_FUNCTIONS;

  static ObjectiveTest *getObjectiveTest(std::string testType);
};

#endif
