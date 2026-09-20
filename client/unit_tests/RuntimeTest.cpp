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
#include "eon/Runtime.h"
#include "catch2/catch_amalgamated.hpp"
#include "eon/PotRegistry.h"

#include <type_traits>
#include <utility>

using eonc::PotRegistry;
using eonc::Runtime;

TEST_CASE("Runtime is move-only", "[runtime]") {
  STATIC_REQUIRE(!std::is_copy_constructible_v<Runtime>);
  STATIC_REQUIRE(!std::is_copy_assignable_v<Runtime>);
  STATIC_REQUIRE(std::is_move_constructible_v<Runtime>);
  STATIC_REQUIRE(std::is_move_assignable_v<Runtime>);

  Runtime a;
  PotRegistry *owned = &a.pots();
  Runtime b = std::move(a);
  REQUIRE(&b.pots() == owned);
  Runtime c;
  c = std::move(b);
  REQUIRE(&c.pots() == owned);
}

TEST_CASE("Runtime pots() is not the process get() singleton", "[runtime]") {
  Runtime rt;
  REQUIRE(&rt.pots() != &PotRegistry::get());
  Runtime other;
  REQUIRE(&rt.pots() != &other.pots());
}
