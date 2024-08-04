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
#include "TestApprovalMain.hpp"

#include "client/Matter.h"
#include "client/io/ConWriter.hpp"
#include "client/io/ConvelWriter.hpp"
#include "client/io/TibbleWriter.hpp"
#include "client/io/XYZWriter.hpp"

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
  return {matter};
}

TEST_CASE("VerifyMatter2Con") {
  auto testMatter = getTestMatter()[0];
  std::string filename = "test_output.con";
  auto conWriter = std::make_unique<eonc::io::ConWriter>();
  conWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}

TEST_CASE("VerifyMatter2XYZ") {
  auto testMatter = getTestMatter()[0];
  std::string filename = "test_output.xyz";
  auto XYZWriter = std::make_unique<eonc::io::XYZWriter>();
  XYZWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}

TEST_CASE("VerifyMatter2Convel") {
  auto testMatter = getTestMatter()[0];
  testMatter.setVelocities(testMatter.getPositionsFree().array() + 3);
  std::string filename = "test_output.con";
  auto ConvelWriter = std::make_unique<eonc::io::ConvelWriter>();
  ConvelWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}

TEST_CASE("VerifyMatter2Tibble") {
  auto testMatter = getTestMatter()[0];
  std::string filename = "test_output.txt";
  auto TibbleWriter = std::make_unique<eonc::io::TibbleWriter>();
  TibbleWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}
