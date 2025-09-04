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
#ifndef NEBTEST_H
#define NEBTEST_H

#include <gtest/gtest.h>
#include <string>

#include "../MatrixHelpers.hpp"
#include "../Matter.h"
#include "../Parameters.h"

using namespace std::string_literals; // For ""s

namespace tests {
class NEBTest : public ::testing::Test {
protected:
  Parameters *params;
  Matter *m1;
  Matter *m2;
  double threshold;
  void SetUp() override;
  void TearDown() override;

public:
  NEBTest();
  virtual ~NEBTest();
};
} /* namespace tests */

#endif /* NEBTEST_H */
