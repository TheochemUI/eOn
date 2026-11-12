/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
*/
#include "catch2/catch_amalgamated.hpp"
#include "eon/InstantonExternal.h"

TEST_CASE("an external instanton is off until one is installed") {
  REQUIRE(eonc::external_instanton().evaluate == nullptr);
  int calls = 0;
  eonc::ExternalInstanton hook{};
  hook.user = &calls;
  hook.evaluate = [](const double *, int, const double *, int, double,
                     double *action, void *user) {
    *action = 1.25;
    *static_cast<int *>(user) += 1;
    return 0;
  };
  eonc::set_external_instanton(hook);
  double action = 0.0;
  const auto installed = eonc::external_instanton();
  REQUIRE(installed.evaluate(nullptr, 0, nullptr, 0, 10.0, &action,
                             installed.user) == 0);
  REQUIRE(calls == 1);
  REQUIRE(action == 1.25);
  eonc::set_external_instanton({});
  REQUIRE(eonc::external_instanton().evaluate == nullptr);
}
