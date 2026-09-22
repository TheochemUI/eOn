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
#include "eon/DynamicsSaddleSearch.h"
#include "eon/Matter.h"
#include "eon/RandomNumbers.h"

#include <memory>
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

TEST_CASE("Dynamics saddle search applies bond-boost forces",
          "[bondboost][dynamics]") {
  Parameters base;
  ParametersLoadAccess::potential_options(base).potential = PotType::LJ;
  ParametersLoadAccess::main_options(base).randomSeed = 42;
  ParametersLoadAccess::saddle_search_options(base).dynamics.temperature =
      300.0;
  const double dt = base.dynamics_options().time_step;
  REQUIRE(dt > 0.0);
  // One equilibration sample, then several boosted steps. A zero rmd_time
  // never records equilibrium lengths, so the bias force stays zero.
  ParametersLoadAccess::dynamics_options(base).steps = 6;
  ParametersLoadAccess::parallel_replica_options(base).dephase_time = 0.0;
  ParametersLoadAccess::saddle_search_options(base)
      .dynamics.state_check_interval = 1.0e6;
  ParametersLoadAccess::saddle_search_options(base).dynamics.record_interval =
      0.0;
  ParametersLoadAccess::hyperdynamics_options(base).rmd_time = dt;
  ParametersLoadAccess::hyperdynamics_options(base).dvmax = 5.0;
  ParametersLoadAccess::hyperdynamics_options(base).qrr = 0.2;
  ParametersLoadAccess::hyperdynamics_options(base).prr = 0.95;
  ParametersLoadAccess::hyperdynamics_options(base).boost_atom_list = "All";

  auto pot = eonc::helpers::makePotential(PotType::LJ, base);

  auto finalPositions = [&](const char *bias) {
    eonc::rng::random(42);
    Parameters params = base;
    ParametersLoadAccess::hyperdynamics_options(params).bias_potential = bias;
    auto matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("reactant.con"));

    DynamicsSaddleSearch search(matter, params);
    const int status = search.run();
    REQUIRE(status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT);
    // The boost object is gone. Accelerations must not call through it.
    REQUIRE_NOTHROW(matter->getBiasForces());
    REQUIRE(matter->getBiasForces().isZero(0.0));
    return AtomMatrix(matter->getPositions());
  };

  const AtomMatrix plain = finalPositions(Hyperdynamics::NONE);
  const AtomMatrix plainAgain = finalPositions(Hyperdynamics::NONE);
  const AtomMatrix boosted = finalPositions(Hyperdynamics::BOND_BOOST);

  REQUIRE(plain.allFinite());
  REQUIRE(boosted.allFinite());
  REQUIRE((plain - plainAgain).norm() < 1e-10);
  REQUIRE((plain - boosted).norm() > 1e-4);
}

} /* namespace tests */
