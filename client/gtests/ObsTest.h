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
#ifndef OBSTEST_H
#define OBSTEST_H

#include <gtest/gtest.h>

namespace tests {
class ObsTest : public ::testing::Test {
public:
  ObsTest();
  virtual ~ObsTest();
};
} /* namespace tests */

#endif /* OBSTEST_H */
