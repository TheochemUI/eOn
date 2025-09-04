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

#ifndef OBJECTIVE_TEST_H
#define OBJECTIVE_TEST_H

#include <string>
#include <vector>

class ObjectiveTest {

public:
  virtual ~ObjectiveTest() {};
  virtual int test() = 0;

  std::vector<std::string> testResults;
  int currentStatus = 0;
  void pushToResults(std::string update);
  void saveTest(std::string testType);

  static const std::string HELPER_FUNCTIONS;

  static ObjectiveTest *getObjectiveTest(std::string testType);
};

#endif
