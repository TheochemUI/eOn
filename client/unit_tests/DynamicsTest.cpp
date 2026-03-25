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

#include "Dynamics.h"
#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class DynamicsFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  Matter *matter;

  DynamicsFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    params.potential_options.potential = PotType::LJ;
    params.main_options.temperature = 300.0;
    params.main_options.randomSeed = 42;
    eonc::helpers::random(42);

    pot = eonc::helpers::makePotential(PotType::LJ, params);
    matter = new Matter(pot, params);
    matter->con2matter(std::string("reactant.con"));
  }

  ~DynamicsFixture() { delete matter; }
};

TEST_CASE_METHOD(DynamicsFixture, "Dynamics oneStep produces finite energy",
                 "[dynamics]") {
  // Reference: LJ 2-atom cluster, E_init = -39.965351
  double E_before = matter->getPotentialEnergy();
  REQUIRE(E_before == Catch::Approx(-39.965351).epsilon(1e-4));

  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();
  dyn.oneStep();

  double E_after = matter->getPotentialEnergy();
  REQUIRE(std::isfinite(E_after));
  // Atom positions must have changed
  REQUIRE(E_after != E_before);
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics 5 steps with fixed seed is deterministic",
                 "[dynamics]") {
  // Run 5 steps, record final energy
  eonc::helpers::random(42); // reset seed
  Dynamics dyn(matter, params);
  dyn.setTemperature(100.0);
  dyn.setThermalVelocity();

  for (int i = 0; i < 5; i++) {
    dyn.oneStep();
  }
  double E_run1 = matter->getPotentialEnergy();
  AtomMatrix pos_run1 = matter->getPositions();

  // Reset and run again with same seed
  matter->con2matter(std::string("reactant.con"));
  eonc::helpers::random(42);
  Dynamics dyn2(matter, params);
  dyn2.setTemperature(100.0);
  dyn2.setThermalVelocity();

  for (int i = 0; i < 5; i++) {
    dyn2.oneStep();
  }
  double E_run2 = matter->getPotentialEnergy();

  // Must be identical (deterministic with same seed)
  REQUIRE(E_run1 == Catch::Approx(E_run2).margin(1e-10));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics Andersen thermostat produces finite energy",
                 "[dynamics][andersen]") {
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  // Run several steps with Andersen thermostat
  for (int i = 0; i < 10; i++) {
    dyn.oneStep(300.0);
  }

  double E = matter->getPotentialEnergy();
  REQUIRE(std::isfinite(E));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics multiple steps accumulate kinetic energy",
                 "[dynamics]") {
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  double KE = matter->getKineticEnergy();
  REQUIRE(KE > 0.0);
  REQUIRE(std::isfinite(KE));

  dyn.oneStep();
  double KE2 = matter->getKineticEnergy();
  REQUIRE(std::isfinite(KE2));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics zero temperature produces zero velocity",
                 "[dynamics]") {
  Dynamics dyn(matter, params);
  dyn.setTemperature(0.0);
  dyn.setThermalVelocity();

  double KE = matter->getKineticEnergy();
  REQUIRE(KE == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics Nose-Hoover thermostat runs without crash",
                 "[dynamics][nose_hoover]") {
  params.thermostat_options.kind = "nose_hoover";
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  for (int i = 0; i < 5; i++) {
    dyn.oneStep(300.0);
  }

  double E = matter->getPotentialEnergy();
  REQUIRE(std::isfinite(E));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics Langevin thermostat runs without crash",
                 "[dynamics][langevin]") {
  params.thermostat_options.kind = "langevin";
  params.thermostat_options.langevin_friction = 0.01;
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  for (int i = 0; i < 5; i++) {
    dyn.oneStep(300.0);
  }

  double E = matter->getPotentialEnergy();
  REQUIRE(std::isfinite(E));
}

} /* namespace tests */
