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
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class MinimizationFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;

  MinimizationFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    params.potential_options.potential = PotType::LJ;
    params.optimizer_options.method = OptType::LBFGS;
    params.optimizer_options.converged_force = 0.01;
    params.optimizer_options.max_iterations = 500;
    params.optimizer_options.max_move = 0.2;

    pot = eonc::helpers::makePotential(PotType::LJ, params);
    matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("reactant.con"));
  }
};

TEST_CASE_METHOD(MinimizationFixture,
                 "Matter::relax converges on displaced LJ cluster",
                 "[minimization]") {
  double E_init = matter->getPotentialEnergy();

  // Displace an atom
  auto pos = matter->getPositions();
  pos(0, 0) += 0.3;
  matter->setPositions(pos);

  double E_displaced = matter->getPotentialEnergy();
  REQUIRE(E_displaced > E_init);

  bool converged = matter->relax(/*quiet=*/true);
  REQUIRE(converged);

  double E_relaxed = matter->getPotentialEnergy();
  REQUIRE(E_relaxed < E_displaced);
  REQUIRE(E_relaxed == Catch::Approx(E_init).margin(0.01));
}

TEST_CASE_METHOD(MinimizationFixture,
                 "Minimization preserves fixed atom positions",
                 "[minimization][fixed_atoms]") {
  // Fix all atoms except the first
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    matter->setFixed(i, true);
  }

  // Record initial positions of fixed atoms
  AtomMatrix initPos = matter->getPositions();

  // Displace the free atom
  auto pos = matter->getPositions();
  pos(0, 0) += 0.5;
  pos(0, 1) -= 0.3;
  matter->setPositions(pos);

  matter->relax(/*quiet=*/true);

  // Fixed atom positions must not change
  AtomMatrix finalPos = matter->getPositions();
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    double drift = (finalPos.row(i) - initPos.row(i)).norm();
    REQUIRE(drift < 1e-10);
  }

  // Forces on fixed atoms must be zero
  const AtomMatrix &forces = matter->getForces();
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    REQUIRE(forces.row(i).norm() < 1e-10);
  }
}

} /* namespace tests */
