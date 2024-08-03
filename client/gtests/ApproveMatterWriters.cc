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
#include "ApprovalTests.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include "client/io/ConWriter.hpp"
#include <iostream>
#include <string>

std::string readFileContent(const std::string &filename) {
  std::ifstream file(filename);
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

std::vector<eonc::Matter> getTestMatter() {
  const auto config =
      toml::table{{"Potential", toml::table{{"potential", "lj"}}}};
  auto pot_default = eonc::makePotential(config);
  auto matter = eonc::Matter(pot_default);
  std::string confile("pos.con");
  matter.con2matter(confile);
  for (size_t idx = 0; idx < 3; idx++) {
    matter.setFixed(idx, true);
  }

  return {matter};
}

TEST_CASE("VerifyMatter2Con") {
  auto testMatter = getTestMatter();
  REQUIRE(testMatter.size() > 0);
  std::string filename = "test_output.con";
  eonc::io::ConWriter conWriter;
  conWriter.write(testMatter[0], filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}
