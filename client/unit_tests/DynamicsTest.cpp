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
#include "eon/DynamicsSaddleSearch.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"

#include <cmath>
#include <filesystem>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Zero force, so relax stays on the stored frame at the default cap.
struct FlatPot final : Potential {
  FlatPot()
      : Potential(PotType::LJ) {}
  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    (void)positions;
    (void)atomicNrs;
    (void)box;
    *energy = 0.0;
    *variance = 0.0;
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
  }
};

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

TEST_CASE_METHOD(
    DynamicsFixture,
    "refineTransition returns the first snapshot that left the reactant",
    "[dynamics]") {
  auto pot = std::make_shared<FlatPot>();
  auto reactant = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(reactant->con2matter(std::string("reactant.con"))));

  auto saddle = std::make_shared<Matter>(*reactant);
  eonc::DynamicsSaddleSearch search(saddle, params);
  auto product = std::make_shared<Matter>(*reactant);

  auto frame = [&](bool reactantFrame) {
    auto snap = std::make_shared<Matter>(*reactant);
    if (!reactantFrame) {
      AtomMatrix pos = snap->getPositions();
      pos(0, 0) += 1.0;
      snap->setPositions(pos);
    }
    return snap;
  };

  std::vector<std::shared_ptr<Matter>> two{frame(true), frame(false)};
  REQUIRE(search.refineTransition(two, product) == 1);

  std::vector<std::shared_ptr<Matter>> four{frame(true), frame(true),
                                            frame(false), frame(false)};
  const int image = search.refineTransition(four, product);
  REQUIRE(image == 2);
  REQUIRE(four[static_cast<size_t>(image - 1)]->compare(*search.reactant));
  REQUIRE_FALSE(four[static_cast<size_t>(image)]->compare(*search.reactant));
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
  // Geometric, not a basin test: any MD step must count as a new state.
  ParametersLoadAccess::structure_comparison_options(params)
      .distance_difference = 1e-6;
  auto shared = std::make_shared<Matter>(pot, params);
  shared->con2matter(std::string("reactant.con"));
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
    explicit CwdGuard(const fs::path &next)
        : old(fs::current_path()) {
      fs::current_path(next);
    }
    ~CwdGuard() { fs::current_path(old); }
  } guard(tmp);

  eonc::rng::random(42);
  eonc::DynamicsSaddleSearch search(shared, params);
  const int status = search.run();
  REQUIRE(status !=
          eonc::MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT);
  const AtomMatrix mode = search.getEigenvector();
  REQUIRE(mode.rows() == matter->numberOfAtoms());
  REQUIRE(mode.cols() == 3);
  REQUIRE(mode.array().isFinite().all());
}

TEST_CASE_METHOD(DynamicsFixture,
                 "Dynamics saddle search survives a zero state-check interval",
                 "[dynamics][saddle_search]") {
  // state_check_interval 0 floors to 0 steps. The search must not take
  // step % 0; the check still runs, and a cold short trajectory stays put.
  ParametersLoadAccess::dynamics_options(params).time_step = 1.0;
  ParametersLoadAccess::dynamics_options(params).steps = 1;
  ParametersLoadAccess::parallel_replica_options(params).dephase_time = 0.0;
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.state_check_interval = 0.0;
  ParametersLoadAccess::saddle_search_options(params).dynamics.temperature =
      0.0;

  auto saddle = std::make_shared<Matter>(*matter);
  DynamicsSaddleSearch search(saddle, params);
  int status = search.run();
  REQUIRE(status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT);
}

TEST_CASE_METHOD(DynamicsFixture,
                 "refineTransition returns no index for an empty snapshot list",
                 "[dynamics][saddle]") {
  auto cluster = std::make_shared<Matter>(*matter);
  DynamicsSaddleSearch search(cluster, params);
  std::vector<std::shared_ptr<Matter>> empty;
  auto prod = std::make_shared<Matter>(*cluster);
  REQUIRE(search.refineTransition(empty, prod) == -1);

  std::vector<std::shared_ptr<Matter>> one{std::make_shared<Matter>(*cluster)};
  REQUIRE(search.refineTransition(one, prod) == 0);
}

TEST_CASE_METHOD(DynamicsFixture,
                 "dynamics saddle search keeps the detecting frame when no "
                 "snapshots were recorded",
                 "[dynamics][saddle]") {
  // Relax must not pull the kicked geometry back onto the reactant, or the
  // state check never enters the empty-snapshot branch.
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  ParametersLoadAccess::structure_comparison_options(params)
      .distance_difference = 1e-8;
  ParametersLoadAccess::dynamics_options(params).time_step = 1.0;
  ParametersLoadAccess::dynamics_options(params).steps = 1;
  ParametersLoadAccess::parallel_replica_options(params).dephase_time = 0.0;
  ParametersLoadAccess::saddle_search_options(params).dynamics.temperature =
      300.0;
  ParametersLoadAccess::saddle_search_options(params).dynamics.record_interval =
      0.0;
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.state_check_interval = 1.0;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  ParametersLoadAccess::neb_options(params).max_iterations = 1;
  ParametersLoadAccess::dimer_options(params).max_iterations = 1;
  ParametersLoadAccess::dimer_options(params).rotations_max = 1;

  auto cluster = std::make_shared<Matter>(pot, params);
  cluster->con2matter(std::string("reactant.con"));
  REQUIRE(cluster->numberOfAtoms() > 0);
  cluster->setMasses(Eigen::VectorXd::Constant(cluster->numberOfAtoms(), 1.0));

  const auto dir =
      std::filesystem::temp_directory_path() / "eon_dyn_saddle_empty";
  std::filesystem::create_directories(dir);
  struct CwdGuard {
    std::filesystem::path previous;
    std::filesystem::path dir;
    CwdGuard(std::filesystem::path next)
        : previous{std::filesystem::current_path()},
          dir{std::move(next)} {
      std::filesystem::current_path(dir);
    }
    ~CwdGuard() {
      std::filesystem::current_path(previous);
      std::filesystem::remove_all(dir);
    }
  } guard(dir);

  eonc::rng::random(42);
  DynamicsSaddleSearch search(cluster, params);
  const int status = search.run();
  REQUIRE(status != MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT);
  REQUIRE(std::isfinite(search.time));
  REQUIRE(search.time == Catch::Approx(1.0).margin(1e-12));
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
