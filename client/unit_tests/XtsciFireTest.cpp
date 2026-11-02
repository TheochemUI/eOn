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
#include "eon/ObjectiveFunction.h"
#include "eon/Parameters.h"
#include "eon/XtsciOptimizer.h"

#include <memory>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {

/// f(x) = sum_i k_i x_i^2 / 2 with distinct k_i, so the FIRE mix and
/// the uphill branch both act.
class Anisotropic final : public ObjectiveFunction {
public:
  explicit Anisotropic(const Parameters &a_params)
      : ObjectiveFunction(a_params),
        m_k{(VectorXd(3) << 1.0, 3.0, 7.0).finished()},
        m_x{VectorXd::Constant(3, 1.0)} {}
  double getEnergy() override {
    return 0.5 * m_k.cwiseProduct(m_x.cwiseAbs2()).sum();
  }
  VectorXd getGradient(bool) override { return m_k.cwiseProduct(m_x); }
  void setPositions(const VectorXd &a_x) override { m_x = a_x; }
  VectorXd getPositions() override { return m_x; }
  int degreesOfFreedom() override { return static_cast<int>(m_x.size()); }
  bool isConverged() override { return getConvergence() < 1e-12; }
  double getConvergence() override {
    return getGradient(false).cwiseAbs().maxCoeff();
  }
  VectorXd difference(const VectorXd &a_a, const VectorXd &a_b) override {
    return a_a - a_b;
  }
  const VectorXd &k() const { return m_k; }

private:
  VectorXd m_k;
  VectorXd m_x;
};

Parameters xtsciParams(const std::string &a_method, double a_timeStep) {
  Parameters params;
  auto &opts = ParametersLoadAccess::optimizer_options(params);
  opts.method = OptType::XTSCI;
  opts.xtsci.method = a_method;
  opts.converged_force = 1e-12;
  opts.max_iterations = 1000;
  opts.time_step = a_timeStep;
  return params;
}

VectorXd runWhole(const std::string &a_method, double a_timeStep,
                  size_t a_iterations, double a_maxMove) {
  const Parameters params = xtsciParams(a_method, a_timeStep);
  auto objf = std::make_shared<Anisotropic>(params);
  XtsciOptimizer optim(objf, params);
  optim.run(a_iterations, a_maxMove);
  return objf->getPositions();
}

} // namespace

// Neither max_move below clamps a step of this quadratic, so only the
// time step may change where FIRE ends up.
TEST_CASE("Optimizer::run gives FIRE the time step, not max_move",
          "[optimizer][rgmin][fire]") {
  const std::string method =
      GENERATE(std::string("fire"), std::string("fire2"));
  CAPTURE(method);
  const VectorXd a = runWhole(method, 0.1, 15, 5.0);
  const VectorXd b = runWhole(method, 0.1, 15, 50.0);
  const VectorXd c = runWhole(method, 0.2, 15, 5.0);
  CAPTURE(a.transpose(), b.transpose(), c.transpose());
  REQUIRE((a - b).cwiseAbs().maxCoeff() < 1e-12);
  REQUIRE((a - c).cwiseAbs().maxCoeff() > 1e-3);
}

} // namespace tests
