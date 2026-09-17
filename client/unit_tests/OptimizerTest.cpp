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

#include "eon/Optimizer.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ConjugateGradients.h"
#include "eon/FIRE.h"
#include "eon/LBFGS.h"
#include "eon/ObjectiveFunction.h"
#include "eon/Parameters.h"
#include "eon/Quickmin.h"
#include "eon/SteepestDescent.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

/// Quadratic objective function: f(x) = 0.5 * (x[0]^2 + x[1]^2)
/// Gradient = [x[0], x[1]], minimum at origin.
class QuadraticObjectiveFunction : public ObjectiveFunction {
  VectorXd m_positions;

public:
  QuadraticObjectiveFunction(const Parameters &params)
      : ObjectiveFunction(params),
        m_positions(VectorXd::Zero(2)) {}

  double getEnergy() override { return 0.5 * m_positions.squaredNorm(); }

  VectorXd getGradient(bool /*fdstep*/ = false) override {
    // gradient of 0.5*x^2 is x, but force = -gradient
    return m_positions;
  }

  void setPositions(const VectorXd &x) override { m_positions = x; }

  VectorXd getPositions() override { return m_positions; }

  int degreesOfFreedom() override { return 2; }

  bool isConverged() override {
    return getConvergence() <
           ParametersLoadAccess::optimizer_options(params).converged_force;
  }

  double getConvergence() override { return m_positions.norm(); }

  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return a - b;
  }
};

static Parameters makeOptParams() {
  Parameters params;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-6;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 1000;
  ParametersLoadAccess::optimizer_options(params).time_step = 0.1;
  ParametersLoadAccess::optimizer_options(params).max_time_step = 1.0;
  ParametersLoadAccess::optimizer_options(params).lbfgs.memory = 20;
  ParametersLoadAccess::optimizer_options(params).lbfgs.auto_scale = true;
  ParametersLoadAccess::optimizer_options(params).lbfgs.inverse_curvature =
      0.01;
  ParametersLoadAccess::optimizer_options(params).lbfgs.angle_reset = true;
  ParametersLoadAccess::optimizer_options(params).lbfgs.distance_reset = true;
  ParametersLoadAccess::optimizer_options(params).sd.alpha = 0.1;
  ParametersLoadAccess::optimizer_options(params).sd.two_point = false;
  ParametersLoadAccess::optimizer_options(params).cg.no_overshooting = false;
  ParametersLoadAccess::optimizer_options(params).cg.knock_out_max_move = false;
  ParametersLoadAccess::optimizer_options(params).cg.line_search = false;
  ParametersLoadAccess::optimizer_options(params).cg.max_iter_before_reset = 0;
  ParametersLoadAccess::optimizer_options(params).cg.line_converged = 0.1;
  ParametersLoadAccess::optimizer_options(params).cg.line_search_max_iter = 5;
  ParametersLoadAccess::main_options(params).finiteDifference = 0.01;
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.bowl_breakout = false;
  return params;
}

TEST_CASE("FIRE throws when the time step collapses", "[optimizer][fire]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).time_step = 1e-7;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  FIRE opt(objf, params);
  REQUIRE_THROWS_AS(
      opt.step(ParametersLoadAccess::optimizer_options(params).max_move),
      std::runtime_error);
  try {
    opt.step(ParametersLoadAccess::optimizer_options(params).max_move);
    FAIL("expected throw");
  } catch (const std::runtime_error &e) {
    REQUIRE_THAT(std::string(e.what()),
                 Catch::Matchers::ContainsSubstring("m_dt is too small"));
  }
}

TEST_CASE("FIRE optimizer converges on quadratic", "[optimizer][fire]") {
  auto params = makeOptParams();
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  FIRE opt(objf, params);
  int status =
      opt.run(1000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 1e-4);
  CHECK(status == 1); // converged
}

TEST_CASE("LBFGS optimizer converges on quadratic", "[optimizer][lbfgs]") {
  auto params = makeOptParams();
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  LBFGS opt(objf, params);
  int status =
      opt.run(1000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == 1);
}

TEST_CASE("LBFGS Zhang-Xu cautious nonmonotone converges on quadratic",
          "[optimizer][lbfgs][zhangxu]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).lbfgs.secant = "zhangxu";
  ParametersLoadAccess::optimizer_options(params).lbfgs.curvature = "cautious";
  ParametersLoadAccess::optimizer_options(params).lbfgs.accept = "nonmonotone";
  ParametersLoadAccess::optimizer_options(params).lbfgs.h0 = "adaptive";
  ParametersLoadAccess::optimizer_options(params).lbfgs.extra_updates = 1;
  ParametersLoadAccess::optimizer_options(params).lbfgs.angle_reset = false;
  ParametersLoadAccess::optimizer_options(params).lbfgs.distance_reset = false;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  LBFGS opt(objf, params);
  int status =
      opt.run(1000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == 1);
}

TEST_CASE("CG optimizer converges on quadratic", "[optimizer][cg]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).converged_force =
      1e-3; // CG needs looser tol
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  ConjugateGradients opt(objf, params);
  int status =
      opt.run(5000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == 1);
}

TEST_CASE("CG with line search converges on quadratic",
          "[optimizer][cg][line_search]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-3;
  ParametersLoadAccess::optimizer_options(params).cg.line_search = true;
  ParametersLoadAccess::optimizer_options(params).cg.line_search_max_iter = 10;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  ConjugateGradients opt(objf, params);
  opt.run(5000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.1);
}

TEST_CASE("CG with no_overshooting converges on quadratic",
          "[optimizer][cg][no_overshoot]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-3;
  ParametersLoadAccess::optimizer_options(params).cg.no_overshooting = true;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  ConjugateGradients opt(objf, params);
  opt.run(5000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.5);
}

TEST_CASE("Quickmin optimizer reduces energy on quadratic",
          "[optimizer][quickmin]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-2;
  ParametersLoadAccess::optimizer_options(params).time_step = 0.01;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.5;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 1.0, 0.5;
  objf->setPositions(start);

  double E_init = objf->getEnergy();
  Quickmin opt(objf, params);
  opt.run(100, ParametersLoadAccess::optimizer_options(params).max_move);
  double E_final = objf->getEnergy();

  // Quickmin should at least reduce energy, even if it doesn't converge
  // tightly on a simple quadratic (it's designed for MD, not optimization)
  REQUIRE(E_final < E_init);
}

TEST_CASE("Quickmin zero-force step stays finite", "[optimizer][quickmin]") {
  auto params = makeOptParams();
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  objf->setPositions(VectorXd::Zero(2));
  Quickmin opt(objf, params);
  REQUIRE(opt.step(ParametersLoadAccess::optimizer_options(params).max_move) ==
          1);
  REQUIRE(objf->getPositions().norm() == Catch::Approx(0.0).margin(1e-15));
}

TEST_CASE("FIRE zero-force step stays finite", "[optimizer][fire]") {
  auto params = makeOptParams();
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  objf->setPositions(VectorXd::Zero(2));
  FIRE opt(objf, params);
  REQUIRE_NOTHROW(
      opt.step(ParametersLoadAccess::optimizer_options(params).max_move));
  REQUIRE(std::isfinite(objf->getPositions().norm()));
}

TEST_CASE("SteepestDescent optimizer converges on quadratic",
          "[optimizer][sd]") {
  auto params = makeOptParams();
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-3;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  SteepestDescent opt(objf, params);
  int status =
      opt.run(5000, ParametersLoadAccess::optimizer_options(params).max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == 1);
}

} /* namespace tests */
