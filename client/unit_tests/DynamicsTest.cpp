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

#include "eon/Dynamics.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"
#include "eon/Parameters.h"

#include <cmath>

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
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::main_options(params).temperature = 300.0;
    ParametersLoadAccess::main_options(params).randomSeed = 42;
    eonc::rng::random(42);

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
  eonc::rng::random(42); // reset seed
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
  eonc::rng::random(42);
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
  ParametersLoadAccess::thermostat_options(params).kind = "nose_hoover";
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  for (int i = 0; i < 5; i++) {
    dyn.oneStep(300.0);
  }

  double E = matter->getPotentialEnergy();
  REQUIRE(std::isfinite(E));
  REQUIRE(std::isfinite(matter->getKineticEnergy()));
  REQUIRE(matter->getKineticEnergy() > 0.0);
}

TEST_CASE_METHOD(DynamicsFixture, "Dynamics run with steps=0 does not move",
                 "[dynamics][steps]") {
  ParametersLoadAccess::dynamics_options(params).steps = 0;
  AtomMatrix before = matter->getPositions();
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();
  dyn.run();
  AtomMatrix after = matter->getPositions();
  REQUIRE((after - before).norm() == Catch::Approx(0.0).margin(1e-15));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics Langevin thermostat runs without crash",
                 "[dynamics][langevin]") {
  ParametersLoadAccess::thermostat_options(params).kind = "langevin";
  ParametersLoadAccess::thermostat_options(params).langevin_friction = 0.01;
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  for (int i = 0; i < 5; i++) {
    dyn.oneStep(300.0);
  }

  double E = matter->getPotentialEnergy();
  REQUIRE(std::isfinite(E));
}

TEST_CASE("Nose-Hoover targets unfixed axes of a partly fixed atom",
          "[dynamics][nose_hoover]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::main_options(params).temperature = 300.0;
  ParametersLoadAccess::thermostat_options(params).kind = Dynamics::NOSE_HOOVER;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter matter(pot, params);
  matter.resize(1);
  matter.setAtomicNr(0, 18);
  matter.setMass(0, 1.0);
  matter.setPeriodic(false);
  // Free in x and y. numberOfFreeAtoms() still counts this atom.
  matter.setFixed(0, 2, 1);
  REQUIRE(matter.getFixed(0) == 0);
  REQUIRE(matter.numberOfFreeAtoms() == 1);

  const double temperature = 300.0;
  const double kB = params.constants().kB;
  // Equipartition on the two free axes: each holds kT/2.
  const double keTarget = kB * temperature;
  const double speed = std::sqrt(keTarget);
  AtomMatrix vel(1, 3);
  vel << speed, speed, 99.0;
  matter.setVelocities(vel);
  REQUIRE(matter.getKineticEnergy() == Catch::Approx(keTarget).epsilon(1e-12));
  REQUIRE(matter.getVelocities()(0, 2) == Catch::Approx(0.0).margin(0.0));

  Dynamics dyn(&matter, DynamicsConfig::fromParams(params));
  dyn.setTemperature(temperature);
  dyn.oneStep();

  // G1 is zero when 2*KE equals n_free*kT, and a lone atom has no force,
  // so the free-axis kinetic energy does not climb toward three axes.
  REQUIRE(matter.getKineticEnergy() == Catch::Approx(keTarget).epsilon(1e-9));
  REQUIRE(matter.getPositions()(0, 2) == Catch::Approx(0.0).margin(1e-15));
  REQUIRE(matter.getVelocities()(0, 2) == Catch::Approx(0.0).margin(0.0));
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics Langevin holds a per-axis frozen coordinate",
                 "[dynamics][langevin][fixed]") {
  // Whole-atom getFixed is false; only z is frozen.
  matter->setFixed(0, 2, 1);
  matter->setPeriodic(false);
  ParametersLoadAccess::thermostat_options(params).kind = "langevin";
  ParametersLoadAccess::thermostat_options(params).langevin_friction = 0.01;
  eonc::rng::random(42);
  Dynamics dyn(matter, params);
  dyn.setTemperature(300.0);
  dyn.setThermalVelocity();

  const AtomMatrix before = matter->getPositions();
  const double frozenZ = before(0, 2);
  for (int i = 0; i < 8; i++) {
    dyn.oneStep();
  }
  const AtomMatrix after = matter->getPositions();
  REQUIRE(after(0, 2) == frozenZ);
  const double freeMove = std::abs(after(0, 0) - before(0, 0)) +
                          std::abs(after(0, 1) - before(0, 1));
  REQUIRE(freeMove > 1e-8);
  REQUIRE(matter->getVelocities()(0, 2) == 0.0);
}

} /* namespace tests */
