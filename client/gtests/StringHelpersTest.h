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
#ifndef STRINGHELPERSTEST_H
#define STRINGHELPERSTEST_H

#include <gtest/gtest.h>
#include <string>

namespace tests {
class StringHelpersTest : public ::testing::Test {
public:
  StringHelpersTest();
  virtual ~StringHelpersTest();
  std::string number_string{"1 2.6 4 5"s};
};
} /* namespace tests */

#endif /* STRINGHELPERSTEST_H */
