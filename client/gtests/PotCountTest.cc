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
#include "client/matter/Matter.h"
#include "client/matter/MatterCreator.hpp"
#include "client/thirdparty/toml.hpp"

using namespace eonc;

TEST_CASE("Potential instance counting and force calls", "[potential]") {
  const auto config =
      toml::table{{"Potential", toml::table{{"potential", "LJ"}}}};

  std::string confile = "pos.con";
  eonc::mat::ConFileParser cfp;

  // Create two potential instances
  auto pot1 = makePotential(config);
  auto pot2 = makePotential(config);

  // Check initial instance counts
  REQUIRE(pot1->getInstances() == 2);
  REQUIRE(pot2->getInstances() == 2);
  REQUIRE(pot1->getTotalForceCalls() == 0);
  REQUIRE(pot2->getTotalForceCalls() == 0);

  // Create Matter objects
  Matter mat1(pot1);
  cfp.parse(mat1, confile);

  Matter mat2(pot2);
  cfp.parse(mat2, confile);

  mat2.getPotentialEnergy();
  REQUIRE(pot2->getTotalForceCalls() == 1);
  REQUIRE(pot1->getTotalForceCalls() == 1);

  // Matter no longer caches, changes
  mat2.getPotentialEnergy();
  REQUIRE(pot2->getTotalForceCalls() == 2);
  REQUIRE(pot1->getTotalForceCalls() == 2);

  // Change positions and call getPotentialEnergy again
  mat2.setPositions(mat1.getPositions().array() + 3.0);
  mat2.getPotentialEnergy();
  REQUIRE(pot2->getTotalForceCalls() == 3);

  // Check final counts and instances
  REQUIRE(pot1->getInstances() == 2);
  REQUIRE(pot1->getTotalForceCalls() == 3);
}

TEST_CASE("Multiple potential instances", "[potential]") {
  const auto config =
      toml::table{{"Potential", toml::table{{"potential", "LJ"}}}};

  std::string confile = "pos.con";
  eonc::mat::ConFileParser cfp;

  // Create multiple potential instances
  auto pot1 = makePotential(config);
  auto pot2 = makePotential(config);
  auto pot3 = makePotential(config);

  REQUIRE(pot1->getInstances() == 3);
  // Scope creep!!!!! Since this is only compiled once, and TEST_CASE isn't
  // doing isolation well the earlier calls are still in the registry
  REQUIRE(pot1->getTotalForceCalls() == 3);

  Matter mat1(pot1);
  cfp.parse(mat1, confile);

  Matter mat2(pot2);
  cfp.parse(mat2, confile);

  Matter mat3(pot3);
  cfp.parse(mat3, confile);

  mat1.getPotentialEnergy();
  REQUIRE(pot1->getTotalForceCalls() == 4);

  // NOTE(rg) :: These are without cache
  mat2.getPotentialEnergy();
  REQUIRE(pot2->getTotalForceCalls() == 5);

  mat3.getPotentialEnergy();
  REQUIRE(pot3->getTotalForceCalls() == 6);
  REQUIRE(pot1->getInstances() == 3);
}
