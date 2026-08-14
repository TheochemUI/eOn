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

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Guards the RgpotAdapter + parameter mapping for every migrated arm: the
// same pinned numbers as the pre-migration in-tree kernels, evaluated
// through helpers::makePotential exactly as jobs do. Kernel numerics
// themselves are pinned upstream in rgpot's own test suite.

TEST_CASE("Adapter LJ matches the pinned reference", "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965351).epsilon(1e-4));
  REQUIRE(matter->getForces().rowwise().norm().maxCoeff() ==
          Catch::Approx(0.004704).epsilon(1e-2));
}

TEST_CASE("Adapter LJCluster matches the pinned reference",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJCLUSTER;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965379).epsilon(1e-4));
}

TEST_CASE("Adapter Morse evaluates and reports shared-safe",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(params);

  REQUIRE(pot->isThreadSafe());
  REQUIRE(pot->isSharedInstanceThreadSafe());
  REQUIRE_FALSE(pot->needsPerImageInstance());

  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
  REQUIRE(matter->getForces().allFinite());
}

TEST_CASE("Adapter ZBL maps config.ini cutoffs through",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::ZBL;
  params.zbl_options.cut_inner = 2.0;
  params.zbl_options.cut_global = 2.5;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
}

// A batch of one repeated system passes even when an implementation
// evaluates the first system and broadcasts the answer, so the two systems
// here differ. LJ takes the default per-system loop rather than a native
// batch, which is the arm every migrated kernel inherits.
TEST_CASE("Adapter forceBatch agrees with per-system force calls",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(params);

  Matter m0(pot, params);
  m0.con2matter(std::string("reactant.con"));
  Matter m1(pot, params);
  m1.con2matter(std::string("reactant.con"));
  AtomMatrix displaced = m1.getPositionsCopy();
  displaced(0, 0) += 0.15;
  m1.setPositions(displaced);

  const long nAtoms = m0.numberOfAtoms();
  const auto nrs0 = m0.getAtomicNrs();
  const auto nrs1 = m1.getAtomicNrs();
  const auto cell0 = m0.getCell();
  const auto cell1 = m1.getCell();
  const AtomMatrix pos0 = m0.getPositionsCopy();
  const AtomMatrix pos1 = m1.getPositionsCopy();

  AtomMatrix singleF0 = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix singleF1 = AtomMatrix::Zero(nAtoms, 3);
  double singleU0{0}, singleU1{0}, var{0};
  pot->force(nAtoms, pos0.data(), nrs0.data(), singleF0.data(), &singleU0, &var,
             cell0.data());
  pot->force(nAtoms, pos1.data(), nrs1.data(), singleF1.data(), &singleU1, &var,
             cell1.data());
  // The displacement has to matter, or the comparison below is vacuous.
  REQUIRE(singleU0 != Catch::Approx(singleU1));

  AtomMatrix batchF0 = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix batchF1 = AtomMatrix::Zero(nAtoms, 3);
  const double *positions[2] = {pos0.data(), pos1.data()};
  const int *atomicNrs[2] = {nrs0.data(), nrs1.data()};
  double *forces[2] = {batchF0.data(), batchF1.data()};
  const double *boxes[2] = {cell0.data(), cell1.data()};
  double energies[2] = {0.0, 0.0};
  double variances[2] = {0.0, 0.0};
  pot->forceBatch(2, nAtoms, positions, atomicNrs, forces, energies, variances,
                  boxes);

  REQUIRE(energies[0] == Catch::Approx(singleU0));
  REQUIRE(energies[1] == Catch::Approx(singleU1));
  REQUIRE((batchF0 - singleF0).cwiseAbs().maxCoeff() < 1e-10);
  REQUIRE((batchF1 - singleF1).cwiseAbs().maxCoeff() < 1e-10);
}

} /* namespace tests */
