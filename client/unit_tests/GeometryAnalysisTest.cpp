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

#include "eon/GeometryAnalysis.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"

#include <cmath>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class GeomFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> m1;
  std::shared_ptr<Matter> m2;

  GeomFixture() {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    pot = eonc::helpers::makePotential(PotType::LJ, params);
    m1 = std::make_shared<Matter>(pot, params);
    m2 = std::make_shared<Matter>(pot, params);
    m1->con2matter(std::string("reactant.con"));
    m2->con2matter(std::string("reactant.con"));
  }
};

TEST_CASE_METHOD(GeomFixture, "identical structures match",
                 "[geometry][identical]") {
  REQUIRE(eonc::geometry::identical(*m1, *m2, 0.1));
}

TEST_CASE_METHOD(GeomFixture, "displaced structure is not identical",
                 "[geometry][identical]") {
  auto pos = m2->getPositions();
  pos(0, 0) += 1.0;
  m2->setPositions(pos);
  REQUIRE_FALSE(eonc::geometry::identical(*m1, *m2, 0.1));
}

TEST_CASE_METHOD(GeomFixture, "rotationExtract returns valid rotation matrix",
                 "[geometry][rotation]") {
  auto r1 = m1->getPositions();
  auto r2 = m1->getPositions();
  auto R = eonc::geometry::rotationExtract(r1, r2);
  // Identity rotation for same positions
  REQUIRE(R.isApprox(Matrix3d::Identity(), 1e-6));
}

TEST_CASE_METHOD(GeomFixture, "translationRemove centers structures",
                 "[geometry][translation]") {
  auto pos = m2->getPositions();
  pos.rowwise() += Eigen::RowVector3d(1.0, 2.0, 3.0);
  m2->setPositions(pos);

  eonc::geometry::translationRemove(*m2, *m1);
  // After removing translation, centroids should be close
  auto c1 = m1->getPositions().colwise().mean();
  auto c2 = m2->getPositions().colwise().mean();
  REQUIRE((c1 - c2).norm() < 0.1);
}

TEST_CASE_METHOD(GeomFixture, "maxAtomMotion returns zero for same positions",
                 "[geometry][motion]") {
  auto pos = m1->getPositions();
  double motion = eonc::geometry::maxAtomMotion(pos - pos);
  REQUIRE(motion == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE_METHOD(GeomFixture, "maxAtomMotionApplied limits displacement",
                 "[geometry][motion]") {
  AtomMatrix disp(m1->numberOfAtoms(), 3);
  disp.setConstant(10.0); // large displacement

  auto limited = eonc::geometry::maxAtomMotionApplied(disp, 0.1);
  double maxMotion = eonc::geometry::maxAtomMotion(limited);
  REQUIRE(maxMotion <= 0.1 + 1e-10);
}

TEST_CASE_METHOD(GeomFixture, "numAtomsMoved counts displaced atoms",
                 "[geometry][motion]") {
  AtomMatrix disp(m1->numberOfAtoms(), 3);
  disp.setZero();
  // Displace first 3 atoms
  disp(0, 0) = 1.0;
  disp(1, 1) = 1.0;
  disp(2, 2) = 1.0;

  long moved = eonc::geometry::numAtomsMoved(disp, 0.5);
  REQUIRE(moved == 3);
}

TEST_CASE_METHOD(GeomFixture, "pushApart separates overlapping atoms",
                 "[geometry][pushApart]") {
  // Move atom 1 very close to atom 0
  auto pos = m1->getPositions();
  pos(1, 0) = pos(0, 0) + 0.01;
  pos(1, 1) = pos(0, 1);
  pos(1, 2) = pos(0, 2);
  m1->setPositions(pos);

  eonc::geometry::pushApart(m1, 1.0);

  // After pushApart, distance should be >= minDistance
  auto newPos = m1->getPositions();
  REQUIRE(newPos.allFinite());
  double dist = (newPos.row(0) - newPos.row(1)).norm();
  REQUIRE(std::isfinite(dist));
  REQUIRE(dist >= 0.5); // should be pushed apart
  REQUIRE(m1->distance(0, 1) >= 0.5);
}

TEST_CASE("pushApart separates a coincident pair", "[geometry][pushApart]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->resize(2);
  m->setAtomicNr(0, 1);
  m->setAtomicNr(1, 1);
  m->setCell(Matrix3d::Identity() * 10.0);
  m->setPeriodic(false);
  AtomMatrix pos(2, 3);
  pos.setZero();
  pos.row(0) << 1.0, 2.0, 3.0;
  pos.row(1) << 1.0, 2.0, 3.0;
  m->setPositions(pos);

  eonc::geometry::pushApart(m, 1.0);

  const auto after = m->getPositions();
  REQUIRE(after.allFinite());
  REQUIRE(m->distance(0, 1) >= 1.0 - 1e-8);
}

TEST_CASE("pushApart follows the minimum image", "[geometry][pushApart]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->resize(2);
  m->setAtomicNr(0, 1);
  m->setAtomicNr(1, 1);
  m->setCell(Matrix3d::Identity() * 10.0);
  m->setPeriodic(true);
  AtomMatrix pos(2, 3);
  pos.setZero();
  // 0.1 across the periodic face, 9.9 the long way through the cell.
  pos(0, 0) = 0.05;
  pos(1, 0) = 9.95;
  m->setPositions(pos);
  REQUIRE(m->distance(0, 1) == Catch::Approx(0.1).margin(1e-8));

  eonc::geometry::pushApart(m, 1.0);

  const auto after = m->getPositions();
  REQUIRE(after.allFinite());
  REQUIRE(m->distance(0, 1) >= 1.0 - 1e-8);
  // Still on opposite sides of the face. A step along the long diagonal
  // throws both atoms into the middle of the cell.
  REQUIRE(after(0, 0) > 0.05);
  REQUIRE(after(0, 0) < 2.0);
  REQUIRE(after(1, 0) < 9.95);
  REQUIRE(after(1, 0) > 8.0);
  REQUIRE(after(0, 1) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(after(1, 1) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(after(0, 2) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(after(1, 2) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE_METHOD(GeomFixture, "sortedR matches an identical geometry",
                 "[geometry][sortedR]") {
  REQUIRE(eonc::geometry::sortedR(*m1, *m2, 1.0));
}

TEST_CASE_METHOD(GeomFixture,
                 "sortedR matches a homonuclear diatomic to itself",
                 "[geometry][sortedR]") {
  m1->resize(2);
  m2->resize(2);
  m1->setAtomicNr(0, 1);
  m1->setAtomicNr(1, 1);
  m2->setAtomicNr(0, 1);
  m2->setAtomicNr(1, 1);

  AtomMatrix pos(2, 3);
  pos.setZero();
  pos(1, 0) = 1.0;
  m1->setPositions(pos);
  m2->setPositions(pos);
  REQUIRE(eonc::geometry::sortedR(*m1, *m2, 1.0));

  pos(1, 0) = 2.5;
  m2->setPositions(pos);
  REQUIRE_FALSE(eonc::geometry::sortedR(*m1, *m2, 1.0));
}

TEST_CASE_METHOD(GeomFixture, "projectOutRotTrans removes rigid body modes",
                 "[geometry][projection]") {
  auto pos = m1->getPositions();
  long ndof = pos.rows() * 3;
  VectorXd step(ndof);
  // Pure translation step
  for (long i = 0; i < pos.rows(); i++) {
    step(3 * i) = 1.0;
    step(3 * i + 1) = 0.0;
    step(3 * i + 2) = 0.0;
  }
  eonc::geometry::projectOutRotTrans(step, pos);
  // After projection, translation component should be ~0
  double sumX = 0;
  for (long i = 0; i < pos.rows(); i++) {
    sumX += step(3 * i);
  }
  REQUIRE(std::abs(sumX) < 1e-6);
}

} /* namespace tests */
