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
#include "ServeMode.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

static helper_functions::test::QuillTestLogger _quill_setup;

TEST_CASE("parseServeSpec single endpoint", "[serve]") {
  auto eps = parseServeSpec("lj:12345");
  REQUIRE(eps.size() == 1);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].host == "localhost");
  CHECK(eps[0].port == 12345);
}

TEST_CASE("parseServeSpec multiple endpoints", "[serve]") {
  auto eps = parseServeSpec("lj:12345,eam_al:12346");
  REQUIRE(eps.size() == 2);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].port == 12345);
  CHECK(eps[1].potential == PotType::EAM_AL);
  CHECK(eps[1].port == 12346);
}

TEST_CASE("parseServeSpec with host", "[serve]") {
  auto eps = parseServeSpec("lj:0.0.0.0:9999");
  REQUIRE(eps.size() == 1);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].host == "0.0.0.0");
  CHECK(eps[0].port == 9999);
}

TEST_CASE("parseServeSpec mixed format", "[serve]") {
  auto eps = parseServeSpec("lj:12345, eam_al:0.0.0.0:12346");
  REQUIRE(eps.size() == 2);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].host == "localhost");
  CHECK(eps[0].port == 12345);
  CHECK(eps[1].potential == PotType::EAM_AL);
  CHECK(eps[1].host == "0.0.0.0");
  CHECK(eps[1].port == 12346);
}

TEST_CASE("parseServeSpec unknown potential is skipped", "[serve]") {
  auto eps = parseServeSpec("nonexistent:12345");
  CHECK(eps.empty());
}

TEST_CASE("parseServeSpec empty spec", "[serve]") {
  auto eps = parseServeSpec("");
  CHECK(eps.empty());
}

TEST_CASE("parseServeSpec whitespace handling", "[serve]") {
  auto eps = parseServeSpec("  lj : 12345 , eam_al : 12346  ");
  // "lj " will be trimmed, ": 12345" -- the port parsing needs the colon
  // The actual parsing trims the token but the colon position is found first.
  // "  lj : 12345 " -> trimmed -> "lj : 12345"
  // first_colon at 2, pot_str="lj", rest=" 12345"
  // Since port is " 12345", stoi skips leading whitespace.
  REQUIRE(eps.size() == 2);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].port == 12345);
}

TEST_CASE("parseServeSpec case insensitive", "[serve]") {
  auto eps = parseServeSpec("LJ:12345");
  REQUIRE(eps.size() == 1);
  CHECK(eps[0].potential == PotType::LJ);
}

TEST_CASE("ServeEndpoint struct members", "[serve]") {
  ServeEndpoint ep;
  ep.potential = PotType::LJ;
  ep.host = "localhost";
  ep.port = 12345;
  CHECK(ep.potential == PotType::LJ);
  CHECK(ep.host == "localhost");
  CHECK(ep.port == 12345);
}
