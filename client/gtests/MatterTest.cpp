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
#include "Matter.h"
#include "Parameters.h"
#include "catch2/catch_amalgamated.hpp"
#include <memory>

using namespace Catch::Matchers;

TEST_CASE("TestCell", "[MatterTest]") {
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot_default, params);
  std::string confile("pos.con");
  m1->con2matter(confile);

  Matrix3S _cell;
  Matrix3S _cellInverse;
  // clang-format off
    _cell << // Comma initialized
        25.0, 0.0, 0.0,
        0.0, 25.0, 0.0,
        0.0, 0.0, 25.0;
  // clang-format on

  REQUIRE_THAT(m1->getCell()(0, 0), WithinAbs(_cell(0, 0), 0.01));
  REQUIRE_THAT(m1->getCell()(1, 1), WithinAbs(_cell(1, 1), 0.01));
  REQUIRE_THAT(m1->getCell()(2, 2), WithinAbs(_cell(2, 2), 0.01));
}

TEST_CASE("SetGetAtomicNrs", "[MatterTest]") {
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot_default, params);
  std::string confile("pos.con");
  m1->con2matter(confile);

  Vector<int> _atmnrs{{8, 8, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 16}};
  Vector<int> _atmnrs2{{16, 16, 12, 12, 12, 12, 2, 2, 2, 2, 2, 2, 32}};

  REQUIRE(m1->getAtomicNrs() == _atmnrs);

  for (auto &atmnr : _atmnrs) {
    atmnr *= 2;
  }
  m1->setAtomicNrs(_atmnrs);

  REQUIRE(m1->getAtomicNrs() == _atmnrs2);
}

TEST_CASE("SetPotential", "[MatterTest]") {
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot_default, params);
  std::string confile("pos.con");
  m1->con2matter(confile);

  double m1_ipot = m1->getPotentialEnergy();
  params->pot.potential = PotType::MORSE_PT;
  auto pot = helper_functions::makePotential(params->pot.potential, params);

  REQUIRE(m1->getPotential() != pot);
  m1->setPotential(pot);

  double m1_fpot = m1->getPotentialEnergy();
  REQUIRE_THAT(m1_ipot, WithinAbs(-8.9245813315, 0.01));
  REQUIRE_THAT(m1_fpot, WithinAbs(1611.8672392832, 0.01));
  REQUIRE(m1_ipot != m1_fpot);
}
