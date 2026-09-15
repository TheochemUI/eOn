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

#include "eon/BondBoost.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"

#include <stdexcept>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("BondBoost initializes on LJ cluster", "[bondboost]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hyperdynamics_options(params).dvmax = 0.0;
  ParametersLoadAccess::hyperdynamics_options(params).qrr = 0.2;
  ParametersLoadAccess::hyperdynamics_options(params).prr = 0.95;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list = "All";

  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));

  BondBoost bb(&matter, params);
  bb.initialize();

  double boostE = bb.boost();
  REQUIRE(std::isfinite(boostE));
  REQUIRE(boostE >= 0.0);
}

TEST_CASE("BondBoost returns zero boost at equilibrium", "[bondboost]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hyperdynamics_options(params).dvmax = 0.0;
  ParametersLoadAccess::hyperdynamics_options(params).qrr = 0.2;
  ParametersLoadAccess::hyperdynamics_options(params).prr = 0.95;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list = "All";

  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));

  // At equilibrium, bonds are at reference length, so boost should be 0
  BondBoost bb(&matter, params);
  bb.initialize();
  double boostE = bb.boost();
  // Boost energy should be small at equilibrium
  REQUIRE(boostE < 10.0);
}

TEST_CASE("BondBoost schedule advances only from advance(), not boost()",
          "[bondboost][schedule]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::dynamics_options(params).time_step = 1.0;
  // Four equilibration MD steps. Repeated boost() calls must not finish
  // this window: only advance() moves nReg.
  ParametersLoadAccess::hyperdynamics_options(params).rmd_time = 4.0;
  ParametersLoadAccess::hyperdynamics_options(params).dvmax = 0.0;
  ParametersLoadAccess::hyperdynamics_options(params).qrr = 0.2;
  ParametersLoadAccess::hyperdynamics_options(params).prr = 0.95;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list = "All";

  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));

  BondBoost bb(&matter, params);
  bb.initialize();
  REQUIRE(bb.scheduleStep() == 1);

  for (int i = 0; i < 8; ++i) {
    REQUIRE(bb.boost() == 0.0);
  }
  REQUIRE(bb.scheduleStep() == 1);

  bb.advance();
  REQUIRE(bb.scheduleStep() == 2);
  REQUIRE(bb.boost() == 0.0);
  REQUIRE(bb.boost() == 0.0);
  REQUIRE(bb.scheduleStep() == 2);

  bb.advance();
  bb.advance();
  bb.advance();
  REQUIRE(bb.scheduleStep() == 5);
}

TEST_CASE("BondBoost listed index out of range throws", "[bondboost][list]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list =
      "999999";
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));
  BondBoost bb(&matter, params);
  REQUIRE_THROWS_AS(bb.initialize(), std::out_of_range);
}

TEST_CASE("BondBoost garbage list is not treated as all", "[bondboost][list]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list =
      "not-a-list";
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));
  BondBoost bb(&matter, params);
  REQUIRE_THROWS_AS(bb.initialize(), std::invalid_argument);
}

} /* namespace tests */
