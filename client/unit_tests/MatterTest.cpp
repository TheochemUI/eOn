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
#include <memory>

using namespace Catch::Matchers;

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Helper to create an LJ Matter loaded from reactant.con (13-atom H cluster)
static std::pair<std::shared_ptr<Matter>, Parameters> makeLJCluster() {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("reactant.con"));
  return {m, params};
}

TEST_CASE("TestCell", "[MatterTest]") {
  auto [m1, params] = makeLJCluster();
  // neb_morse reactant.con has a ~102x103x103 cell
  REQUIRE(m1->getCell()(0, 0) > 100.0);
  REQUIRE(m1->getCell()(1, 1) > 100.0);
  REQUIRE(m1->getCell()(2, 2) > 100.0);
}

TEST_CASE("SetGetAtomicNrs", "[MatterTest]") {
  auto [m1, params] = makeLJCluster();
  // 13 hydrogen atoms
  REQUIRE(m1->numberOfAtoms() == 13);
  auto nrs = m1->getAtomicNrs();
  for (int i = 0; i < 13; i++) {
    REQUIRE(nrs(i) == 1);
  }

  // Double all atomic numbers
  auto nrs2 = nrs;
  for (auto &n : nrs2) {
    n *= 2;
  }
  m1->setAtomicNrs(nrs2);
  auto result = m1->getAtomicNrs();
  for (int i = 0; i < 13; i++) {
    REQUIRE(result(i) == 2);
  }
}

TEST_CASE("SetPotential changes energy", "[MatterTest]") {
  auto [m1, params] = makeLJCluster();

  double e_lj = m1->getPotentialEnergy();
  REQUIRE(std::isfinite(e_lj));
  REQUIRE(e_lj < 0.0); // LJ cluster has negative binding energy

  params.potential_options.potential = PotType::MORSE_PT;
  auto pot_morse = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  REQUIRE(m1->getPotential() != pot_morse);
  m1->setPotential(pot_morse);

  double e_morse = m1->getPotentialEnergy();
  REQUIRE(std::isfinite(e_morse));
  REQUIRE(e_lj != e_morse);
}

TEST_CASE("PBC wrapping brings displaced atom back into cell",
          "[MatterTest][pbc]") {
  auto [m1, params] = makeLJCluster();

  auto cell = m1->getCell();
  auto origPos = m1->getPositions();

  // Displace atom 0 far outside the cell
  auto pos = m1->getPositionsCopy();
  pos(0, 0) += cell(0, 0) * 3.0;
  pos(0, 1) += cell(1, 1) * 2.0;
  m1->setPositions(pos);

  auto wrappedPos = m1->getPositions();
  for (int axis = 0; axis < 3; axis++) {
    REQUIRE(wrappedPos(0, axis) >= -0.01);
    REQUIRE(wrappedPos(0, axis) < cell(axis, axis) + 0.01);
  }
}

TEST_CASE("Copy constructor preserves positions, cell, and atomic numbers",
          "[MatterTest][copy]") {
  auto [m1, params] = makeLJCluster();

  Matter m2(*m1);

  REQUIRE(m2.numberOfAtoms() == m1->numberOfAtoms());
  REQUIRE(m2.getPositions().isApprox(m1->getPositions(), 1e-12));
  REQUIRE(m2.getCell().isApprox(m1->getCell(), 1e-12));
  REQUIRE(m2.getAtomicNrs() == m1->getAtomicNrs());
}

TEST_CASE("setPositions marks forces stale", "[MatterTest][force_cache]") {
  auto [m1, params] = makeLJCluster();

  double e1 = m1->getPotentialEnergy();
  long calls1 = m1->getForceCalls();
  REQUIRE(calls1 > 0);

  auto pos = m1->getPositionsCopy();
  pos(0, 0) += 0.5;
  m1->setPositions(pos);

  double e2 = m1->getPotentialEnergy();
  long calls2 = m1->getForceCalls();
  REQUIRE(calls2 > calls1);
  REQUIRE(e1 != e2);
}

TEST_CASE("getForces zeroes fixed atoms", "[MatterTest][forces]") {
  auto [m1, params] = makeLJCluster();
  m1->setFixed(0, true);

  auto forces = m1->getForces();
  REQUIRE(forces.row(0).norm() < 1e-10);
  // Free atom forces should be non-zero
  REQUIRE(forces.row(1).norm() > 0.0);
}

TEST_CASE("getForcesFree is const and drops fixed rows",
          "[MatterTest][forces][const]") {
  auto [m1, params] = makeLJCluster();
  m1->setFixed(0, true);
  const Matter &cref = *m1;
  AtomMatrix freeF = cref.getForcesFree();
  REQUIRE(freeF.rows() == m1->numberOfFreeAtoms());
  REQUIRE(freeF.rows() == m1->numberOfAtoms() - 1);
  // All remaining rows should be non-zero for this LJ cluster
  REQUIRE(freeF.row(0).norm() > 0.0);
}

TEST_CASE("MinimumImage PBC centers fractional coords",
          "[MatterTest][pbc][convention]") {
  auto [m1, params] = makeLJCluster();
  auto cell = m1->getCell();
  auto pos = m1->getPositionsCopy();
  // Place atom 0 just past +0.5 in fractional x so Legacy keeps it in [0,1)
  // near 1.0 while MinimumImage maps toward the cell origin (negative frac).
  pos(0, 0) = cell(0, 0) * 0.9;
  pos(0, 1) = cell(1, 1) * 0.1;
  pos(0, 2) = cell(2, 2) * 0.1;
  m1->setPositions(pos);

  Matter legacy(*m1);
  legacy.setPbcConvention(eonc::PbcConvention::Legacy);
  // Force re-wrap by setPositions (applies PBC when enabled)
  legacy.setPositions(pos);
  auto legacyPos = legacy.getPositions();

  Matter mic(*m1);
  mic.setPbcConvention(eonc::PbcConvention::MinimumImage);
  mic.setPositions(pos);
  auto micPos = mic.getPositions();

  REQUIRE(legacy.getPbcConvention() == eonc::PbcConvention::Legacy);
  REQUIRE(mic.getPbcConvention() == eonc::PbcConvention::MinimumImage);
  // MinimumImage should differ from Legacy for this placement (0.9 -> -0.1 frac)
  REQUIRE((legacyPos.row(0) - micPos.row(0)).norm() > 1e-6);
}

TEST_CASE("getForcesRaw includes fixed atoms", "[MatterTest][forces]") {
  auto [m1, params] = makeLJCluster();
  m1->setFixed(0, true);

  auto rawForces = m1->getForcesRaw();
  auto maskedForces = m1->getForces();
  // Raw should have non-zero on fixed atom
  REQUIRE(rawForces.row(0).norm() > maskedForces.row(0).norm());
}

TEST_CASE("distanceTo computes correct distance", "[MatterTest]") {
  auto [m1, params] = makeLJCluster();
  Matter m2(*m1);

  // Same positions = zero distance
  REQUIRE(m1->distanceTo(m2) == Catch::Approx(0.0).margin(1e-10));

  // Displace one atom
  auto pos = m2.getPositionsCopy();
  pos(0, 0) += 1.0;
  m2.setPositions(pos);
  REQUIRE(m1->distanceTo(m2) > 0.0);
}

TEST_CASE("compare identifies identical structures", "[MatterTest]") {
  auto [m1, params] = makeLJCluster();
  Matter m2(*m1);

  REQUIRE(m1->compare(m2) == true);

  // Displace significantly
  auto pos = m2.getPositionsCopy();
  pos(0, 0) += 5.0;
  m2.setPositions(pos);
  REQUIRE(m1->compare(m2) == false);
}

TEST_CASE("getKineticEnergy returns finite value", "[MatterTest]") {
  auto [m1, params] = makeLJCluster();
  double KE = m1->getKineticEnergy();
  REQUIRE(std::isfinite(KE));
  REQUIRE(KE >= 0.0); // kinetic energy is non-negative
}

TEST_CASE("relax converges LJ cluster", "[MatterTest][relax]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.converged_force = 0.001;
  params.optimizer_options.max_iterations = 50;
  params.optimizer_options.max_move = 0.2;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot, params);
  m1->con2matter(std::string("reactant.con"));

  auto pos = m1->getPositionsCopy();
  pos(0, 0) += 0.01;
  m1->setPositions(pos);

  double E_before = m1->getPotentialEnergy();
  bool converged = m1->relax(/*quiet=*/true);
  double E_after = m1->getPotentialEnergy();

  REQUIRE(converged);
  REQUIRE(E_after <= E_before);
  // SVN reference: -39.965352 for relaxed LJ cluster
  REQUIRE(E_after == Catch::Approx(-39.965352).epsilon(1e-4));
}

} /* namespace tests */
