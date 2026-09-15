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

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"
#include <filesystem>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class SaddleSearchFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;

  SaddleSearchFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::optimizer_options(params).method = OptType::CG;
    ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.1;
    ParametersLoadAccess::dimer_options(params).improved = true;
    ParametersLoadAccess::dimer_options(params).converged_angle = 0.01;
    ParametersLoadAccess::dimer_options(params).max_iterations = 50;
    ParametersLoadAccess::saddle_search_options(params).minmode_method =
        LowestEigenmode::MINMODE_DIMER;
    ParametersLoadAccess::saddle_search_options(params).max_iterations = 100;
    ParametersLoadAccess::saddle_search_options(params).converged_force = 0.05;
    ParametersLoadAccess::saddle_search_options(params).max_energy = 20.0;

    pot = eonc::helpers::makePotential(PotType::LJ, params);
    matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("reactant.con"));

    // Displace one atom to break symmetry
    auto pos = matter->getPositions();
    pos(0, 0) += 0.3;
    pos(0, 1) -= 0.2;
    matter->setPositions(pos);
  }
};

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch runs without crashing",
                 "[saddle_search]") {
  // Create a random initial mode
  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  int status = search.run();

  // Status should be a valid enum value (0 through 21)
  REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status <= MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST);
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch reports finite eigenvalue after run",
                 "[saddle_search]") {
  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  search.run();

  double eigenvalue = search.getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch hits max iterations with low limit",
                 "[saddle_search]") {
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-10;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  int status = search.run();

  // Should hit max iterations or some non-GOOD status
  REQUIRE(status != MinModeSaddleSearch::STATUS_GOOD);
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "write_movies writes per-iteration mode files",
                 "[saddle_search][ra6]") {
  namespace fs = std::filesystem;
  ParametersLoadAccess::debug_options(params).write_movies = true;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-10;

  const auto tmp = fs::temp_directory_path() / "eon_ra6_modes";
  fs::create_directories(tmp);
  const auto old = fs::current_path();
  fs::current_path(tmp);

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();
  MinModeSaddleSearch search(matter, mode, matter->getPotentialEnergy(), params,
                             pot);
  search.run();

  REQUIRE(fs::exists("mode_000.dat"));
  fs::current_path(old);
  fs::remove_all(tmp);
}

// Issue #20: unfeasible / unconverged climb must not report STATUS_GOOD.
// finalizeClimbStatus is the shipped policy used by MinModeSaddleSearch::run.
// If the guard body is deleted (always returns climbStatus), these fail.
TEST_CASE("finalizeClimbStatus refuses STATUS_GOOD when unconverged (#20)",
          "[saddle_search][issue_20]") {
  using S = MinModeSaddleSearch;
  REQUIRE(
      S::finalizeClimbStatus(S::STATUS_GOOD, /*objectiveConverged=*/false) ==
      S::STATUS_BAD_MAX_ITERATIONS);
  REQUIRE(S::finalizeClimbStatus(S::STATUS_GOOD, /*objectiveConverged=*/true) ==
          S::STATUS_GOOD);
  // Non-GOOD statuses are left alone regardless of convergence flag.
  REQUIRE(S::finalizeClimbStatus(S::STATUS_BAD_MAX_ITERATIONS, false) ==
          S::STATUS_BAD_MAX_ITERATIONS);
  REQUIRE(S::finalizeClimbStatus(S::STATUS_BAD_HIGH_ENERGY, false) ==
          S::STATUS_BAD_HIGH_ENERGY);
  REQUIRE(S::finalizeClimbStatus(S::STATUS_DIMER_RESTORED_BEST, false) ==
          S::STATUS_DIMER_RESTORED_BEST);
  // Removing the unconverged branch of finalizeClimbStatus makes the first
  // REQUIRE fail (would return STATUS_GOOD).
}

TEST_CASE_METHOD(
    SaddleSearchFixture,
    "MinModeSaddleSearch run never returns STATUS_GOOD if unconverged (#20)",
    "[saddle_search][issue_20]") {
  // Force a non-converged climb: tiny iteration budget + absurd force target.
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-20;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-20;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);
  int status = search.run(1);

  // Real entry path: with an unconverged climb objective the search must not
  // report success (covers finalizeClimbStatus applied inside run()).
  REQUIRE(status != MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status != MinModeSaddleSearch::STATUS_INIT);
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch forces on fixed atoms remain zero",
                 "[saddle_search][fixed_atoms]") {
  // Fix all atoms except the first
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    matter->setFixed(i, true);
  }

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  search.run();

  // After saddle search, forces on fixed atoms must still be zero
  const AtomMatrix &forces = matter->getForces();
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    REQUIRE(forces.row(i).norm() < 1e-10);
  }
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch with Lanczos eigenmode",
                 "[saddle_search][lanczos]") {
  ParametersLoadAccess::saddle_search_options(params).minmode_method =
      LowestEigenmode::MINMODE_LANCZOS;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 50;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);
  int status = search.run();

  REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status <= MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST);
  REQUIRE(std::isfinite(search.getEigenvalue()));
}

TEST_CASE_METHOD(SaddleSearchFixture, "MinModeSaddleSearch with classic Dimer",
                 "[saddle_search][classic_dimer]") {
  ParametersLoadAccess::dimer_options(params).improved = false; // classic dimer, not improved
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 50;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);
  int status = search.run();

  REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(std::isfinite(search.getEigenvalue()));
}

} /* namespace tests */
