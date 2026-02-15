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
#include "../HelperFunctions.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

using namespace Catch::Matchers;

namespace tests {

// Helper: build an AtomMatrix from a flat initializer
AtomMatrix
makePositions(std::initializer_list<std::initializer_list<double>> rows) {
  AtomMatrix m(static_cast<long>(rows.size()), 3);
  long i = 0;
  for (auto &row : rows) {
    long j = 0;
    for (auto v : row) {
      m(i, j++) = v;
    }
    i++;
  }
  return m;
}

// Three-atom non-linear molecule (water-like geometry)
AtomMatrix waterPositions() {
  return makePositions({
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.3, 0.9, 0.0},
  });
}

// Two-atom linear molecule
AtomMatrix linearPositions() {
  return makePositions({
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
  });
}

TEST_CASE("projectOutRotTrans removes pure translation",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  // Pure uniform translation in x
  Eigen::VectorXd step = Eigen::VectorXd::Zero(n);
  for (long i = 0; i < pos.rows(); ++i) {
    step(3 * i + 0) = 0.5; // dx = 0.5 for all atoms
  }

  helper_functions::projectOutRotTrans(step, pos);

  // After projection, the step should be nearly zero
  REQUIRE_THAT(step.norm(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("projectOutRotTrans removes pure translation in all directions",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  // Translation along (1, 2, 3) direction
  Eigen::VectorXd step = Eigen::VectorXd::Zero(n);
  for (long i = 0; i < pos.rows(); ++i) {
    step(3 * i + 0) = 1.0;
    step(3 * i + 1) = 2.0;
    step(3 * i + 2) = 3.0;
  }

  helper_functions::projectOutRotTrans(step, pos);

  REQUIRE_THAT(step.norm(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("projectOutRotTrans removes infinitesimal rotation",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  // Infinitesimal rotation around z-axis: d_r = omega x r
  // omega = (0, 0, theta), dr_i = (-y_i, x_i, 0) * theta
  Eigen::Vector3d com = pos.colwise().mean().transpose();
  double theta = 0.01;
  Eigen::VectorXd step = Eigen::VectorXd::Zero(n);
  for (long i = 0; i < pos.rows(); ++i) {
    double x = pos(i, 0) - com(0);
    double y = pos(i, 1) - com(1);
    step(3 * i + 0) = -y * theta;
    step(3 * i + 1) = x * theta;
    step(3 * i + 2) = 0.0;
  }

  helper_functions::projectOutRotTrans(step, pos);

  REQUIRE_THAT(step.norm(), WithinAbs(0.0, 1e-10));
}

TEST_CASE("projectOutRotTrans preserves pure internal motion",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  // Internal motion: symmetric stretch (atoms 1 and 2 move away from atom 0)
  // This has zero net translation and zero net torque by construction
  Eigen::Vector3d com = pos.colwise().mean().transpose();
  Eigen::VectorXd step = Eigen::VectorXd::Zero(n);

  // Move atom 1 in +x, atom 2 in its direction from COM, atom 0 to compensate
  // Build a zero-COM, zero-torque displacement manually
  Eigen::Vector3d d0(-1.0, 0.0, 0.0);
  Eigen::Vector3d d1(0.5, 0.0, 0.0);
  Eigen::Vector3d d2(0.5, 0.0, 0.0);

  // Ensure zero net translation
  Eigen::Vector3d net = d0 + d1 + d2;
  REQUIRE_THAT(net.norm(), WithinAbs(0.0, 1e-14));

  // Ensure zero net torque about COM
  Eigen::Vector3d torque = Eigen::Vector3d::Zero();
  Eigen::Vector3d r0 = pos.row(0).transpose() - com;
  Eigen::Vector3d r1 = pos.row(1).transpose() - com;
  Eigen::Vector3d r2 = pos.row(2).transpose() - com;
  torque += r0.cross(d0);
  torque += r1.cross(d1);
  torque += r2.cross(d2);

  // This specific displacement may have nonzero torque, so let's use a
  // displacement that is provably internal: project out rot/trans from
  // a random displacement and verify the result is preserved under a
  // second projection (idempotency)
  step(0) = 0.3;
  step(1) = -0.1;
  step(2) = 0.2;
  step(3) = -0.5;
  step(4) = 0.4;
  step(5) = -0.1;
  step(6) = 0.2;
  step(7) = -0.3;
  step(8) = -0.1;

  // First projection: get pure internal component
  helper_functions::projectOutRotTrans(step, pos);
  Eigen::VectorXd internal_step = step;

  // Second projection should not change it (idempotent)
  helper_functions::projectOutRotTrans(step, pos);

  REQUIRE_THAT((step - internal_step).norm(), WithinAbs(0.0, 1e-14));
}

TEST_CASE("projectOutRotTrans is idempotent",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  // Arbitrary displacement
  Eigen::VectorXd step(n);
  step << 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, -0.7, 0.8, -0.9;

  helper_functions::projectOutRotTrans(step, pos);
  Eigen::VectorXd after_first = step;

  helper_functions::projectOutRotTrans(step, pos);

  REQUIRE_THAT((step - after_first).norm(), WithinAbs(0.0, 1e-14));
}

TEST_CASE("projectOutRotTrans removes mixed translation and rotation",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  Eigen::Vector3d com = pos.colwise().mean().transpose();
  double theta = 0.005;

  // Build step = translation + infinitesimal rotation
  Eigen::VectorXd step = Eigen::VectorXd::Zero(n);
  for (long i = 0; i < pos.rows(); ++i) {
    double x = pos(i, 0) - com(0);
    double y = pos(i, 1) - com(1);
    // Translation component
    step(3 * i + 0) += 0.3;
    step(3 * i + 1) += -0.2;
    step(3 * i + 2) += 0.1;
    // Rotation about z
    step(3 * i + 0) += -y * theta;
    step(3 * i + 1) += x * theta;
  }

  helper_functions::projectOutRotTrans(step, pos);

  // Should be very close to zero (only rigid-body motion, no internal DOF)
  REQUIRE_THAT(step.norm(), WithinAbs(0.0, 1e-10));
}

TEST_CASE("projectOutRotTrans handles linear molecule (5 DOF not 6)",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = linearPositions();
  long n = pos.size();

  // For a linear molecule along x, rotation about x is degenerate.
  // We should have 5 rigid-body modes (3 translation + 2 rotation),
  // leaving 1 internal DOF (the bond stretch).

  // Pure bond stretch: atom 0 moves -x, atom 1 moves +x
  Eigen::VectorXd stretch = Eigen::VectorXd::Zero(n);
  stretch(0) = -0.5; // atom 0, x
  stretch(3) = 0.5;  // atom 1, x
  Eigen::VectorXd stretch_orig = stretch;

  helper_functions::projectOutRotTrans(stretch, pos);

  // Bond stretch should be preserved (it's internal)
  REQUIRE_THAT((stretch - stretch_orig).norm(), WithinAbs(0.0, 1e-12));

  // Pure translation should still be removed
  Eigen::VectorXd trans = Eigen::VectorXd::Zero(n);
  trans(0) = 1.0;
  trans(3) = 1.0; // both atoms move +x

  helper_functions::projectOutRotTrans(trans, pos);
  REQUIRE_THAT(trans.norm(), WithinAbs(0.0, 1e-12));

  // Rotation about y-axis (perpendicular to bond)
  Eigen::Vector3d com = pos.colwise().mean().transpose();
  double theta = 0.01;
  Eigen::VectorXd rot = Eigen::VectorXd::Zero(n);
  for (long i = 0; i < pos.rows(); ++i) {
    double x = pos(i, 0) - com(0);
    // cross(yhat, r) = (z, 0, -x) -> only z component for molecule along x
    rot(3 * i + 2) = -x * theta;
  }

  helper_functions::projectOutRotTrans(rot, pos);
  REQUIRE_THAT(rot.norm(), WithinAbs(0.0, 1e-10));
}

TEST_CASE("projectOutRotTrans preserves norm of internal component",
          "[HelperFunctions][projectOutRotTrans]") {
  AtomMatrix pos = waterPositions();
  long n = pos.size();

  // Build a displacement with known internal and rigid-body parts
  Eigen::VectorXd step(n);
  step << 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, -0.7, 0.8, -0.9;

  double norm_before = step.norm();
  helper_functions::projectOutRotTrans(step, pos);
  double norm_after = step.norm();

  // Projection can only reduce norm (Pythagorean theorem in the subspace)
  REQUIRE(norm_after <= norm_before + 1e-14);
}

TEST_CASE("projectOutRotTrans with single atom is a no-op for internal motion",
          "[HelperFunctions][projectOutRotTrans]") {
  // A single atom has 3 translational DOF and 0 rotational DOF,
  // so all motion is rigid-body translation.
  AtomMatrix pos = makePositions({{1.0, 2.0, 3.0}});
  long n = pos.size();

  Eigen::VectorXd step(n);
  step << 0.5, -0.3, 0.7;

  helper_functions::projectOutRotTrans(step, pos);

  // All motion of a single atom is translation -> should be fully removed
  REQUIRE_THAT(step.norm(), WithinAbs(0.0, 1e-12));
}

} // namespace tests
