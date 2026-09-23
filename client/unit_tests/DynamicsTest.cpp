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
#include "eon/DynamicsSaddleSearch.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"

#include <filesystem>

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

TEST_CASE_METHOD(DynamicsFixture,
                 "DynamicsSaddleSearch skipped NEB keeps an n-atom mode",
                 "[dynamics][saddle_search]") {
  namespace fs = std::filesystem;
  // Relax is a no-op so a short trajectory cannot fall back into the
  // reactant basin before the state check. Caps here are the test's own.
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  ParametersLoadAccess::neb_options(params).max_iterations = 0;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  ParametersLoadAccess::dimer_options(params).rotations_max = 2;
  ParametersLoadAccess::dimer_options(params).rotations_min = 1;
  const double dt = params.dynamics_options().time_step;
  REQUIRE(dt > 0.0);
  ParametersLoadAccess::dynamics_options(params).steps = 40;
  ParametersLoadAccess::parallel_replica_options(params).dephase_time = 0.0;
  ParametersLoadAccess::saddle_search_options(params).dynamics.temperature =
      5000.0;
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.state_check_interval = dt;
  ParametersLoadAccess::saddle_search_options(params).dynamics.record_interval =
      dt;

  const auto tmp = fs::temp_directory_path() / "eon_dyn_skip_neb";
  fs::remove_all(tmp);
  fs::create_directories(tmp);
  struct CwdGuard {
    fs::path old;
    explicit CwdGuard(const fs::path &next) : old(fs::current_path()) {
      fs::current_path(next);
    }
    ~CwdGuard() { fs::current_path(old); }
  } guard(tmp);

  eonc::rng::random(42);
  auto shared = std::make_shared<Matter>(*matter);
  eonc::DynamicsSaddleSearch search(shared, params);
  const int status = search.run();
  REQUIRE(status !=
          eonc::MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT);
  const AtomMatrix mode = search.getEigenvector();
  REQUIRE(mode.rows() == matter->numberOfAtoms());
  REQUIRE(mode.cols() == 3);
  REQUIRE(mode.array().isFinite().all());
}

} /* namespace tests */
