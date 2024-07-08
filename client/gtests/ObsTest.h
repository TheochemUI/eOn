/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
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
