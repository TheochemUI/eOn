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
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/XtsciOptimizer.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

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

/// Guenole et al. 2020 algorithm 2 with semi-implicit Euler and the
/// parameters of its Table 1, written out from the paper.
struct Guenole2020 {
  VectorXd v;
  double dt;
  double dtMax;
  double dtMin;
  double alpha{0.25};
  size_t nPos{0};
  size_t steps{0};
  size_t halfStepsInDelay{0};
  size_t halfStepsAfterDelay{0};

  Guenole2020(Eigen::Index a_dim, double a_dt0)
      : v{VectorXd::Zero(a_dim)},
        dt{a_dt0},
        dtMax{10.0 * a_dt0},
        dtMin{0.02 * a_dt0} {}

  VectorXd displacement(const VectorXd &a_force) {
    constexpr size_t nDelay = 20;
    VectorXd dx = VectorXd::Zero(v.size());
    if (a_force.dot(v) > 0.0) {
      ++nPos;
      if (nPos > nDelay) {
        dt = std::min(dt * 1.1, dtMax);
        alpha *= 0.99;
      }
    } else {
      nPos = 0;
      const bool moving = v.norm() > 0.0;
      if (steps >= nDelay) {
        if (dt * 0.5 >= dtMin) {
          dt *= 0.5;
        }
        alpha = 0.25;
        halfStepsAfterDelay += moving ? 1 : 0;
      } else {
        halfStepsInDelay += moving ? 1 : 0;
      }
      dx -= 0.5 * dt * v;
      v.setZero();
    }
    ++steps;
    v += dt * a_force;
    const double fnorm = a_force.norm();
    if (fnorm > 0.0) {
      v = (1.0 - alpha) * v + alpha * v.norm() / fnorm * a_force;
    }
    dx += dt * v;
    return dx;
  }
};

std::vector<VectorXd> stepTrajectory(const std::string &a_method,
                                     double a_timeStep, size_t a_steps) {
  const Parameters params = xtsciParams(a_method, a_timeStep);
  auto objf = std::make_shared<Anisotropic>(params);
  XtsciOptimizer optim(objf, params);
  std::vector<VectorXd> path;
  for (size_t i = 0; i < a_steps; ++i) {
    optim.step(5.0);
    path.push_back(objf->getPositions());
  }
  return path;
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

TEST_CASE("fire2 runs Guenole et al. 2020 algorithm 2",
          "[optimizer][rgmin][fire]") {
  constexpr double dt0 = 0.3;
  constexpr size_t nSteps = 60;
  const Parameters params = xtsciParams("fire2", dt0);
  Anisotropic reference_objf(params);
  Guenole2020 reference(3, dt0);
  std::vector<VectorXd> expected;
  VectorXd x = reference_objf.getPositions();
  for (size_t i = 0; i < nSteps; ++i) {
    x += reference.displacement(-reference_objf.k().cwiseProduct(x));
    expected.push_back(x);
  }
  CAPTURE(reference.halfStepsInDelay, reference.halfStepsAfterDelay);
  REQUIRE(reference.halfStepsInDelay > 0);
  REQUIRE(reference.halfStepsAfterDelay > 0);

  const auto fire2 = stepTrajectory("fire2", dt0, nSteps);
  const auto fire = stepTrajectory("fire", dt0, nSteps);
  double fire2Gap = 0.0;
  double fireGap = 0.0;
  for (size_t i = 0; i < nSteps; ++i) {
    fire2Gap =
        std::max(fire2Gap, (fire2[i] - expected[i]).cwiseAbs().maxCoeff());
    fireGap = std::max(fireGap, (fire[i] - expected[i]).cwiseAbs().maxCoeff());
  }
  CAPTURE(fire2Gap, fireGap);
  REQUIRE(fire2Gap < 1e-10);
  REQUIRE(fireGap > 1e-2);
}

TEST_CASE("Optimizer::run keeps fire2 on the session schedule",
          "[optimizer][rgmin][fire]") {
  constexpr size_t nSteps = 40;
  const VectorXd stepped = stepTrajectory("fire2", 0.3, nSteps).back();
  const VectorXd whole = runWhole("fire2", 0.3, nSteps, 5.0);
  CAPTURE(stepped.transpose(), whole.transpose());
  REQUIRE((stepped - whole).cwiseAbs().maxCoeff() < 1e-12);
}

// The two schedules relax a displaced LJ13 to the same minimum by
// different paths, so the geometries they stop at differ.
TEST_CASE("fire2 relaxes a displaced LJ cluster on its own path",
          "[optimizer][rgmin][fire][minimization]") {
  std::vector<AtomMatrix> ends;
  for (const std::string method : {"fire", "fire2"}) {
    CAPTURE(method);
    Parameters params = xtsciParams(method, 0.1);
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
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
    const size_t before = pot->forceCallCounter;
    REQUIRE(matter.relax(/*quiet=*/true));
    CAPTURE(pot->forceCallCounter - before);
    REQUIRE(matter.maxForce() < 0.01);
    REQUIRE(matter.getPotentialEnergy() == Catch::Approx(eMin).margin(0.01));
    ends.push_back(matter.getPositions());
  }
  const double gap = (ends[0] - ends[1]).cwiseAbs().maxCoeff();
  CAPTURE(gap);
  REQUIRE(gap > 1e-5);
}

} // namespace tests
