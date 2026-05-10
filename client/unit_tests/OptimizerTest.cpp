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

#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "FIRE.h"
#include "LBFGS.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include "Quickmin.h"
#include "StatusTypes.h"
#include "SteepestDescent.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

using eonc::StepResult;

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
    return getConvergence() < params.optimizer_options.converged_force;
  }

  double getConvergence() override { return m_positions.norm(); }

  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return a - b;
  }
};

static Parameters makeOptParams() {
  Parameters params;
  params.optimizer_options.converged_force = 1e-6;
  params.optimizer_options.max_move = 0.2;
  params.optimizer_options.max_iterations = 1000;
  params.optimizer_options.time_step = 0.1;
  params.optimizer_options.max_time_step = 1.0;
  params.optimizer_options.lbfgs.memory = 20;
  params.optimizer_options.lbfgs.auto_scale = true;
  params.optimizer_options.lbfgs.inverse_curvature = 0.01;
  params.optimizer_options.lbfgs.angle_reset = true;
  params.optimizer_options.lbfgs.distance_reset = true;
  params.optimizer_options.sd.alpha = 0.1;
  params.optimizer_options.sd.two_point = false;
  params.optimizer_options.cg.no_overshooting = false;
  params.optimizer_options.cg.knock_out_max_move = false;
  params.optimizer_options.cg.line_search = false;
  params.optimizer_options.cg.max_iter_before_reset = 0;
  params.optimizer_options.cg.line_converged = 0.1;
  params.optimizer_options.cg.line_search_max_iter = 5;
  params.main_options.finiteDifference = 0.01;
  params.saddle_search_options.confine_positive.bowl_breakout = false;
  return params;
}

TEST_CASE("FIRE optimizer converges on quadratic", "[optimizer][fire]") {
  auto params = makeOptParams();
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  FIRE opt(objf, params);
  StepResult status = opt.run(1000, params.optimizer_options.max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 1e-4);
  CHECK(status == StepResult::Converged);
}

TEST_CASE("LBFGS optimizer converges on quadratic", "[optimizer][lbfgs]") {
  auto params = makeOptParams();
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  LBFGS opt(objf, params);
  StepResult status = opt.run(1000, params.optimizer_options.max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == StepResult::Converged);
}

TEST_CASE("CG optimizer converges on quadratic", "[optimizer][cg]") {
  auto params = makeOptParams();
  params.optimizer_options.converged_force = 1e-3; // CG needs looser tol
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  ConjugateGradients opt(objf, params);
  StepResult status = opt.run(5000, params.optimizer_options.max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == StepResult::Converged);
}

TEST_CASE("CG with line search converges on quadratic",
          "[optimizer][cg][line_search]") {
  auto params = makeOptParams();
  params.optimizer_options.converged_force = 1e-3;
  params.optimizer_options.cg.line_search = true;
  params.optimizer_options.cg.line_search_max_iter = 10;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  ConjugateGradients opt(objf, params);
  opt.run(5000, params.optimizer_options.max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.1);
}

TEST_CASE("CG with no_overshooting converges on quadratic",
          "[optimizer][cg][no_overshoot]") {
  auto params = makeOptParams();
  params.optimizer_options.converged_force = 1e-3;
  params.optimizer_options.cg.no_overshooting = true;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  ConjugateGradients opt(objf, params);
  opt.run(5000, params.optimizer_options.max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.5);
}

TEST_CASE("Quickmin optimizer reduces energy on quadratic",
          "[optimizer][quickmin]") {
  auto params = makeOptParams();
  params.optimizer_options.converged_force = 1e-2;
  params.optimizer_options.time_step = 0.01;
  params.optimizer_options.max_move = 0.5;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 1.0, 0.5;
  objf->setPositions(start);

  double E_init = objf->getEnergy();
  Quickmin opt(objf, params);
  opt.run(100, params.optimizer_options.max_move);
  double E_final = objf->getEnergy();

  // Quickmin should at least reduce energy, even if it doesn't converge
  // tightly on a simple quadratic (it's designed for MD, not optimization)
  REQUIRE(E_final < E_init);
}

TEST_CASE("SteepestDescent optimizer converges on quadratic",
          "[optimizer][sd]") {
  auto params = makeOptParams();
  params.optimizer_options.converged_force = 1e-3;
  auto objf = std::make_shared<QuadraticObjectiveFunction>(params);
  VectorXd start(2);
  start << 5.0, 3.0;
  objf->setPositions(start);

  SteepestDescent opt(objf, params);
  StepResult status = opt.run(5000, params.optimizer_options.max_move);
  auto final_pos = objf->getPositions();

  REQUIRE(final_pos.norm() < 0.01);
  CHECK(status == StepResult::Converged);
}

} /* namespace tests */
