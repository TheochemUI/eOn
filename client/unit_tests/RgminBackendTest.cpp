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
#include "eon/ObjectiveFunction.h"
#include "eon/Optimizer.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#ifdef WITH_RGMIN
#include "eon/XtsciOptimizer.h"
#endif

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {

/// f(x) = |x|^2 / 2, optionally claiming its gradient is not the
/// derivative of its energy the way the band and min-mode objectives do.
class Quadratic final : public ObjectiveFunction {
public:
  Quadratic(const Parameters &a_params, bool a_conservative)
      : ObjectiveFunction(a_params),
        m_conservative{a_conservative},
        m_x{VectorXd::Constant(3, 1.0)} {}
  double getEnergy() override { return 0.5 * m_x.squaredNorm(); }
  VectorXd getGradient(bool) override { return m_x; }
  void setPositions(const VectorXd &a_x) override { m_x = a_x; }
  VectorXd getPositions() override { return m_x; }
  int degreesOfFreedom() override { return static_cast<int>(m_x.size()); }
  bool isConverged() override { return getConvergence() < 1e-6; }
  double getConvergence() override { return m_x.cwiseAbs().maxCoeff(); }
  VectorXd difference(const VectorXd &a_a, const VectorXd &a_b) override {
    return a_a - a_b;
  }
  bool gradientMatchesEnergy() const override { return m_conservative; }

private:
  bool m_conservative;
  VectorXd m_x;
};

Parameters ljParams(OptType a_method, const std::string &a_backend) {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::optimizer_options(params).method = a_method;
  ParametersLoadAccess::optimizer_options(params).backend = a_backend;
  ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 1000;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  return params;
}

} // namespace

TEST_CASE("opt_backend picks the rgmin method for each eOn optimizer",
          "[optimizer][rgmin]") {
  Parameters params;
  const Quadratic objf(params, true);
  for (OptType method : {OptType::LBFGS, OptType::FIRE, OptType::CG,
                         OptType::SD, OptType::QM, OptType::XTSCI}) {
    CAPTURE(static_cast<int>(method));
    REQUIRE_FALSE(eonc::helpers::create::rgminMethodFor(objf, method, params));
  }
  ParametersLoadAccess::optimizer_options(params).backend = "rgmin";
  using eonc::helpers::create::rgminMethodFor;
  REQUIRE(rgminMethodFor(objf, OptType::LBFGS, params) == "lbfgs");
  REQUIRE(rgminMethodFor(objf, OptType::FIRE, params) == "fire");
  REQUIRE(rgminMethodFor(objf, OptType::CG, params) == "polak_ribiere");
  REQUIRE(rgminMethodFor(objf, OptType::SD, params) == "steepest");
  REQUIRE(rgminMethodFor(objf, OptType::QM, params) == "quickmin");
  // opt_method = xtsci already names its rgmin method in [Xtsci] method.
  REQUIRE_FALSE(rgminMethodFor(objf, OptType::XTSCI, params));
}

TEST_CASE("opt_backend = rgmin leaves non-conservative objectives in eOn",
          "[optimizer][rgmin]") {
  Parameters params = ljParams(OptType::LBFGS, "rgmin");
  auto band = std::make_shared<Quadratic>(params, false);
  REQUIRE_FALSE(
      eonc::helpers::create::rgminMethodFor(*band, OptType::LBFGS, params));
  auto optim = eonc::helpers::create::mkOptim(band, OptType::LBFGS, params);
  REQUIRE(optim != nullptr);
#ifdef WITH_RGMIN
  REQUIRE(dynamic_cast<XtsciOptimizer *>(optim.get()) == nullptr);
#endif
}

#ifdef WITH_RGMIN

TEST_CASE("opt_backend = rgmin relaxes a displaced LJ cluster",
          "[optimizer][rgmin][minimization]") {
  const auto [method, token] =
      GENERATE(std::pair{OptType::LBFGS, std::string("lbfgs")},
               std::pair{OptType::FIRE, std::string("fire")},
               std::pair{OptType::CG, std::string("polak_ribiere")});
  CAPTURE(token);
  Parameters params = ljParams(method, "rgmin");
  auto lj = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto pot = std::make_shared<eonc::helpers::test::CountingPotential>(lj);
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));
  const double eMin = matter.getPotentialEnergy();
  auto pos = matter.getPositions();
  pos(0, 0) += 0.3;
  pos(3, 1) -= 0.2;
  matter.setPositions(pos);
  const double eDisplaced = matter.getPotentialEnergy();
  REQUIRE(eDisplaced > eMin);

  auto optim = eonc::helpers::create::mkOptim(
      std::make_shared<Quadratic>(params, true), method, params);
  auto *rgmin = dynamic_cast<XtsciOptimizer *>(optim.get());
  REQUIRE(rgmin != nullptr);
  REQUIRE(rgmin->method() == token);

  const size_t counterBefore = pot->forceCallCounter;
  const size_t kernelBefore = pot->kernelCalls;
  REQUIRE(matter.relax(/*quiet=*/true));
  const size_t counted = pot->forceCallCounter - counterBefore;
  const size_t kernel = pot->kernelCalls - kernelBefore;
  CAPTURE(counted, kernel);
  REQUIRE(kernel > 0);
  REQUIRE(counted == kernel);
  REQUIRE(matter.maxForce() < 0.01);
  REQUIRE(matter.getPotentialEnergy() == Catch::Approx(eMin).margin(0.01));
}

TEST_CASE("opt_method = xtsci still reads [Xtsci] method",
          "[optimizer][rgmin]") {
  Parameters params = ljParams(OptType::XTSCI, "eon");
  ParametersLoadAccess::optimizer_options(params).xtsci.method = "fire";
  auto optim = eonc::helpers::create::mkOptim(
      std::make_shared<Quadratic>(params, true), OptType::XTSCI, params);
  auto *rgmin = dynamic_cast<XtsciOptimizer *>(optim.get());
  REQUIRE(rgmin != nullptr);
  REQUIRE(rgmin->method() == "fire");
}

#else

TEST_CASE("opt_backend = rgmin without rgmin names the build option",
          "[optimizer][rgmin]") {
  Parameters params = ljParams(OptType::LBFGS, "rgmin");
  try {
    static_cast<void>(eonc::helpers::create::mkOptim(
        std::make_shared<Quadratic>(params, true), OptType::LBFGS, params));
    FAIL("mkOptim accepted opt_backend = rgmin without rgmin");
  } catch (const std::runtime_error &e) {
    REQUIRE(std::string(e.what()).find("-Dwith_rgmin") != std::string::npos);
  }
}

#endif

} // namespace tests
