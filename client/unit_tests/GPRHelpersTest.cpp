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

#include "eon/GPRHelpers.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("eon_parameters_to_gpr stores check_derivatives as true or false",
          "[gpr][helpers]") {
  Parameters parameters;
  auto &opt = ParametersLoadAccess::gpr_dimer_options(parameters).opt_params;

  opt.check_derivatives = false;
  auto off = eonc::helpers::eon_parameters_to_gpr(parameters);
  REQUIRE(off.check_derivative.value == "false");
  REQUIRE(off.check_derivative.value.size() == 5);

  opt.check_derivatives = true;
  auto on = eonc::helpers::eon_parameters_to_gpr(parameters);
  REQUIRE(on.check_derivative.value == "true");
  REQUIRE(on.check_derivative.value.size() == 4);
}

} // namespace tests
