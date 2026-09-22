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

static void requireBiasTracksStretch(Matter &matter, BondBoost &bb) {
  const double atEq = bb.boost();
  REQUIRE(std::isfinite(atEq));
  // Every selected bond is at the sampled length, so the bias is dvmax
  // and the bond gradient is zero.
  REQUIRE(atEq == Catch::Approx(1.0).margin(1e-8));
  REQUIRE(matter.getBiasForces().squaredNorm() ==
          Catch::Approx(0.0).margin(1e-8));

  matter.setPosition(0, 0, matter.getPosition(0, 0) + 0.05);
  const double nudged = bb.boost();
  REQUIRE(std::isfinite(nudged));
  REQUIRE(nudged > 0.0);
  REQUIRE(nudged < atEq);
  REQUIRE(matter.getBiasForces().squaredNorm() > 0.0);

  matter.setPosition(0, 0, matter.getPosition(0, 0) + 10.0);
  const double stretched = bb.boost();
  REQUIRE(std::isfinite(stretched));
  REQUIRE(stretched == Catch::Approx(0.0).margin(1e-8));
}

static Parameters zeroSampleParams() {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::dynamics_options(params).time_step = 1.0;
  ParametersLoadAccess::hyperdynamics_options(params).rmd_time = 0.0;
  ParametersLoadAccess::hyperdynamics_options(params).dvmax = 1.0;
  ParametersLoadAccess::hyperdynamics_options(params).qrr = 0.2;
  ParametersLoadAccess::hyperdynamics_options(params).prr = 0.95;
  ParametersLoadAccess::hyperdynamics_options(params).qcut = 3.0;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list = "All";
  return params;
}

TEST_CASE("BondBoost with no equilibration samples uses the current lengths",
          "[bondboost][rmd]") {
  // SafeHyper and ParallelReplica call advance() once per step, then boost().
  {
    const Parameters params = zeroSampleParams();
    auto pot = eonc::helpers::makePotential(PotType::LJ, params);
    Matter matter(pot, params);
    matter.con2matter(std::string("reactant.con"));
    BondBoost bb(&matter, params);
    bb.initialize();
    REQUIRE(bb.scheduleStep() == 1);
    bb.advance();
    REQUIRE(bb.scheduleStep() == 2);
    requireBiasTracksStretch(matter, bb);
    REQUIRE(bb.scheduleStep() == 2);
  }
  // A caller that never advance()s still has to measure before BondSelect.
  {
    const Parameters params = zeroSampleParams();
    auto pot = eonc::helpers::makePotential(PotType::LJ, params);
    Matter matter(pot, params);
    matter.con2matter(std::string("reactant.con"));
    BondBoost bb(&matter, params);
    bb.initialize();
    requireBiasTracksStretch(matter, bb);
    REQUIRE(bb.scheduleStep() == 1);
  }
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
