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
#include "catch2/catch_amalgamated.hpp"
#include "client/BaseStructures.h"
#include "client/Parser.hpp"

using namespace eonc;

TEST_CASE("get_enum_toml valid input", "[get_enum_toml]") {
  const auto config = toml::table{{"Main", toml::table{{"job", "point"}}}};

  REQUIRE_NOTHROW([&]() {
    auto jtype = get_enum_toml<JobType>(config["Main"]["job"]).value();
    REQUIRE(jtype == JobType::Point);
  }());
}

TEST_CASE("get_enum_toml invalid enum value", "[get_enum_toml]") {
  const auto config = toml::table{{"Main", toml::table{{"job", "Invalid"}}}};

  REQUIRE_THROWS_WITH(get_enum_toml<JobType>(config["Main"]["job"]),
                      "Invalid does not map to a valid option");
}

// TODO(rg) :: Low prio, Inexcplicably fails on MacOS
// TEST_CASE("get_enum_toml missing node", "[get_enum_toml]") {
//   const auto config = toml::table{{"Main", toml::table{{}}}};

//   REQUIRE_THROWS_WITH(get_enum_toml<JobType>(config["Main"]["job"]),
//                       "Node not found in TOML");
// }

TEST_CASE("get_enum_toml non-string node", "[get_enum_toml]") {
  const auto config = toml::table{{"Main", toml::table{{"job", 42}}}};

  REQUIRE_THROWS_WITH(get_enum_toml<JobType>(config["Main"]["job"]),
                      "Node not found in TOML");
}

TEST_CASE("get_enum_toml case insensitive match", "[get_enum_toml]") {
  const auto config = toml::table{{"Main", toml::table{{"job", "pOint"}}}};

  REQUIRE_NOTHROW([&]() {
    auto jtype = get_enum_toml<JobType>(config["Main"]["job"]).value();
    REQUIRE(jtype == JobType::Point);
  }());
}
