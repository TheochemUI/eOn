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
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

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
    params.potential_options.potential = PotType::LJ;
    params.optimizer_options.method = OptType::CG;
    params.optimizer_options.converged_force = 0.01;
    params.optimizer_options.max_move = 0.1;
    params.dimer_options.improved = true;
    params.dimer_options.converged_angle = 0.01;
    params.dimer_options.max_iterations = 50;
    params.saddle_search_options.minmode_method =
        LowestEigenmode::MINMODE_DIMER;
    params.saddle_search_options.max_iterations = 100;
    params.saddle_search_options.converged_force = 0.05;
    params.saddle_search_options.max_energy = 20.0;

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
  params.saddle_search_options.max_iterations = 2;
  params.saddle_search_options.converged_force = 1e-10;

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
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_LANCZOS;
  params.saddle_search_options.max_iterations = 50;

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
  params.dimer_options.improved = false; // classic dimer, not improved
  params.saddle_search_options.max_iterations = 50;

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
