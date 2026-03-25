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

#include "../StringHelpers.hpp"
#include "catch2/catch_amalgamated.hpp"
#include <string>
#include <vector>

using namespace std::string_literals;

namespace tests {

static const std::string number_string{"1 2.6 4 5"s};

TEST_CASE("IsNotNum", "[StringHelpers]") {
  std::string tester{"Coordinates of component   blah"s};
  auto split_strings = eonc::helpers::get_split_strings(tester);
  REQUIRE(split_strings.size() == 4);
  for (auto substring : split_strings) {
    REQUIRE(eonc::helpers::isNumber(substring) == false);
  }
}

TEST_CASE("IsNum", "[StringHelpers]") {
  auto split_strings = eonc::helpers::get_split_strings(number_string);
  REQUIRE(split_strings.size() == 4);
  for (auto substring : split_strings) {
    REQUIRE(eonc::helpers::isNumber(substring) == true);
  }
  REQUIRE(eonc::helpers::isNumber("2.6"s) == true);
}

TEST_CASE("SplitStrings", "[StringHelpers]") {
  auto split_strings = eonc::helpers::get_split_strings(number_string);
  REQUIRE(split_strings.size() == 4);
  CHECK(split_strings[0] == "1"s);
  CHECK(split_strings[1] == "2.6"s);
  CHECK(split_strings[2] == "4"s);
  CHECK(split_strings[3] == "5"s);
}

TEST_CASE("ValsBasic", "[StringHelpers]") {
  auto split_strings =
      eonc::helpers::get_val_from_string<double>(number_string);
  REQUIRE(split_strings.size() == 4);
  CHECK(split_strings[0] == 1.0);
  CHECK(split_strings[1] == 2.60);
  CHECK(split_strings[2] == 4.0);
  CHECK(split_strings[3] == 5.0);
  auto split_strings_size_t =
      eonc::helpers::get_val_from_string<size_t>(number_string);
  REQUIRE(split_strings_size_t.size() == 4);
  CHECK(split_strings_size_t[0] == 1);
  CHECK(split_strings_size_t[1] == 2);
  CHECK(split_strings_size_t[2] == 4);
  CHECK(split_strings_size_t[3] == 5);
}

TEST_CASE("ValsTrunc", "[StringHelpers]") {
  auto split_strings =
      eonc::helpers::get_val_from_string<double>(number_string, 2);
  REQUIRE(split_strings.size() == 2);
  CHECK(split_strings[0] == 1.0);
  CHECK(split_strings[1] == 2.60);
  auto split_strings_size_t =
      eonc::helpers::get_val_from_string<size_t>(number_string, 3);
  REQUIRE(split_strings_size_t.size() == 3);
  CHECK(split_strings_size_t[0] == 1);
  CHECK(split_strings_size_t[1] == 2);
  CHECK(split_strings_size_t[2] == 4);
}

TEST_CASE("MixedString", "[StringHelpers]") {
  std::string tester{"Coordinates of component   2"s};
  auto split_strings = eonc::helpers::get_split_strings(tester);
  REQUIRE(split_strings.size() == 4);
  CHECK(split_strings[0] == "Coordinates"s);
  CHECK(split_strings[1] == "of"s);
  CHECK(split_strings[2] == "component"s);
  CHECK(split_strings[3] == "2"s);
  auto tval = eonc::helpers::get_val_from_string<double>(tester);
  REQUIRE(tval.size() == 1);
  CHECK(tval[0] == 2);
}

} /* namespace tests */
