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

#include "client/io/WriteCreator.hpp"
#include "client/matter/Matter.h"
#include "client/matter/MatterCreator.hpp"

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
  auto cachelot = cachelot::cache::Cache::Create(
      eonc::cache_memory, eonc::page_size, eonc::hash_initial, true);
  auto matter = eonc::Matter(pot_default, &cachelot);
  std::string confile("pos.con");
  eonc::mat::ConFileParser cfp;
  cfp.parse(matter, confile);
  return {matter};
}

TEST_CASE("VerifyMatter2Con") {
  auto testMatter = getTestMatter()[0];
  const auto config = toml::table{{"Main", toml::table{{"write", "con"}}}};
  std::string filename = "test_output.con";
  auto conWriter = eonc::io::mkWriter(config);
  conWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}

TEST_CASE("VerifyMatter2XYZ") {
  auto testMatter = getTestMatter()[0];
  const auto config = toml::table{{"Main", toml::table{{"write", "xyz"}}}};
  std::string filename = "test_output.xyz";
  auto XYZWriter = eonc::io::mkWriter(config);
  XYZWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}

TEST_CASE("VerifyMatter2Convel") {
  auto testMatter = getTestMatter()[0];
  const auto config = toml::table{{"Main", toml::table{{"write", "convel"}}}};
  testMatter.setVelocities(testMatter.getPositionsFree().array() + 3);
  std::string filename = "test_output.con";
  auto ConvelWriter = eonc::io::mkWriter(config);
  ConvelWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}

TEST_CASE("VerifyMatter2Tibble") {
  auto testMatter = getTestMatter()[0];
  const auto config = toml::table{{"Main", toml::table{{"write", "tibble"}}}};
  std::string filename = "test_output.txt";
  auto TibbleWriter = eonc::io::mkWriter(config);
  TibbleWriter->write(testMatter, filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}
