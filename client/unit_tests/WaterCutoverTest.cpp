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

void loadWater(Matter &matter, const AtomMatrix &positions,
               const VectorXi &numbers) {
  matter.resize(positions.rows());
  matter.setPositions(positions);
  matter.setAtomicNrs(numbers);
  Matrix3d cell = Matrix3d::Zero();
  cell(0, 0) = 40.0;
  cell(1, 1) = 40.0;
  cell(2, 2) = 40.0;
  matter.setCell(cell);
}

AtomMatrix distortedMonomer() {
  AtomMatrix pos(3, 3);
  pos << 0.80, 0.10, 0.50, -0.70, 0.05, 0.55, 0.0, 0.0, 0.0;
  return pos;
}

VectorXi waterNumbers() {
  VectorXi z(3);
  z << 1, 1, 8;
  return z;
}

TEST_CASE("TIP4P monomer matches the rgpot pin", "[water][tip4p]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::TIP4P;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  Matter matter(pot, params);
  loadWater(matter, distortedMonomer(), waterNumbers());
  REQUIRE(matter.getPotentialEnergy() ==
          Catch::Approx(0.141344038570).epsilon(1e-8));
  double maxForce = matter.getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(4.479603395144).epsilon(1e-6));
  REQUIRE(matter.getForces().colwise().sum().norm() < 1e-8);
  REQUIRE_FALSE(matter.getPotential()->isSharedInstanceThreadSafe());
}

TEST_CASE("SPC/E monomer matches the rgpot pin", "[water][spce]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::SPCE;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  Matter matter(pot, params);
  loadWater(matter, distortedMonomer(), waterNumbers());
  REQUIRE(matter.getPotentialEnergy() ==
          Catch::Approx(0.532963694950).epsilon(1e-8));
  double maxForce = matter.getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(8.458875831491).epsilon(1e-6));
  REQUIRE(matter.getForces().colwise().sum().norm() < 1e-8);
  REQUIRE_FALSE(matter.getPotential()->isSharedInstanceThreadSafe());
}

TEST_CASE("TIP4P on platinum matches the rgpot pin", "[water][tip4p_pt]") {
  AtomMatrix pos(4, 3);
  pos << 0.80, 0.10, 2.50, -0.70, 0.05, 2.55, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0;
  VectorXi z(4);
  z << 1, 1, 8, 78;
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::TIP4P_PT;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  Matter matter(pot, params);
  loadWater(matter, pos, z);
  REQUIRE(matter.getPotentialEnergy() ==
          Catch::Approx(1.642427227258).epsilon(1e-8));
  double maxForce = matter.getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(15.396689248363).epsilon(1e-6));
  REQUIRE(matter.getForces().allFinite());
  REQUIRE_FALSE(matter.getPotential()->isSharedInstanceThreadSafe());
}

} // namespace tests
