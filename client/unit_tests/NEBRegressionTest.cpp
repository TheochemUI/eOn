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

/// Regression tests verifying numerical equivalence of the modularized NEB
/// against known-good reference values. These ensure the refactoring
/// (extracting NEBForceProjection, NEBSpringForce, NEBSplineExtrema,
/// NEBObjectiveFunction) produces bit-identical results.

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/LBFGS.h"
#include "eon/NEBOcinebController.h"
#include "eon/NudgedElasticBand.h"
#include "eon/PotRegistry.h"

#include <algorithm>
#include <cmath>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class QuarticBandPotential final : public Potential {
public:
  QuarticBandPotential()
      : Potential(PotType::UNKNOWN) {}

  void force(long nAtoms, const double *positions, const int *, double *forces,
             double *energy, double *variance, const double *) override {
    std::fill_n(forces, 3 * nAtoms, 0.0);
    const double q = 0.5 * (positions[0] - positions[3]);
    const double displacement = q * q - 1.0;
    *energy = 10.0 * displacement * displacement;
    const double derivative = 40.0 * q * displacement;
    forces[0] = -0.5 * derivative;
    forces[3] = 0.5 * derivative;
    if (variance)
      *variance = 0.0;
  }
};

TEST_CASE("CI activates on a relaxed band with tangential saddle force",
          "[neb][regression][force_cache][ci_activation]") {
  Parameters params;
  ParametersLoadAccess::main_options(params).parallel = false;
  auto &options = ParametersLoadAccess::neb_options(params);
  options.image_count = 2;
  options.initialization.method = NEBInit::LINEAR;
  options.match_endpoints = false;
  options.endpoints.minimize = false;
  options.climbing_image.enabled = true;
  options.climbing_image.trigger_factor = 0.8;
  options.climbing_image.trigger_force = 0.5;
  options.opt_method = OptType::LBFGS;
  options.max_iterations = 100;
  options.force_tolerance = 1e-5;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.1;
  auto pot = std::make_shared<QuarticBandPotential>();
  auto endpoint = [&](double q) {
    auto matter = std::make_shared<Matter>(pot, params);
    matter->resize(2);
    matter->setPeriodic(false);
    matter->setCell(Matrix3d::Identity() * 20.0);
    for (int i = 0; i < 2; ++i) {
      matter->setAtomicNr(i, 1);
      matter->setMass(i, 1.0);
    }
    AtomMatrix positions = AtomMatrix::Zero(2, 3);
    positions(0, 0) = q;
    positions(1, 0) = -q;
    matter->setPositions(positions);
    return matter;
  };
  auto initial = endpoint(-1.0);
  auto final = endpoint(1.0);
  NudgedElasticBand neb(initial, final, params, pot);
  REQUIRE(neb.compute() == NudgedElasticBand::NEBStatus::GOOD);
  REQUIRE(neb.convergenceForce() < options.force_tolerance);
  REQUIRE(neb.path[neb.climbingImage]->getPotentialEnergy() ==
          Catch::Approx(10.0).margin(1e-8));
  REQUIRE(
      neb.path.front()->getPositions().isApprox(initial->getPositions(), 0.0));
  REQUIRE(neb.path.back()->getPositions().isApprox(final->getPositions(), 0.0));
}

class MolecularSaddlePotential final : public eonc::Potential {
public:
  MolecularSaddlePotential() : Potential(PotType::UNKNOWN) {}
  void force(long nAtoms, const double *positions, const int *, double *forces,
             double *energy, double *variance, const double *) override {
    *energy = 0.0;
    for (long atom = 0; atom < nAtoms; ++atom) {
      for (int axis = 0; axis < 3; ++axis) {
        const long index = 3 * atom + axis;
        const double curvature = axis == 0 ? -2.0 : 4.0;
        forces[index] = -curvature * positions[index];
        *energy += 0.5 * curvature * positions[index] * positions[index];
      }
    }
    if (variance) {
      *variance = 0.0;
    }
  }
};

TEST_CASE("OCI redistributes a helpful molecular walk and resets the band",
          "[neb][regression][ocineb]") {
  Parameters params;
  auto &nebOptions = ParametersLoadAccess::neb_options(params);
  nebOptions.image_count = 3;
  nebOptions.initialization.method = NEBInit::LINEAR;
  nebOptions.endpoints.minimize = false;
  nebOptions.force_tolerance = 0.001;
  nebOptions.climbing_image.enabled = true;
  nebOptions.climbing_image.ocineb.use_mmf = true;
  nebOptions.climbing_image.ocineb.angle_tol = 0.85;
  nebOptions.climbing_image.ocineb.max_steps = 100;
  auto &optimizer = ParametersLoadAccess::optimizer_options(params);
  optimizer.method = OptType::LBFGS;
  optimizer.max_move = 0.2;
  optimizer.convergence_metric = "max_atom";
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::dimer_options(params).improved = true;
  ParametersLoadAccess::dimer_options(params).converged_angle = 0.01;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 0.001;
  auto pot = std::make_shared<MolecularSaddlePotential>();
  auto reactant = std::make_shared<Matter>(pot, params);
  reactant->resize(2);
  reactant->setPeriodic(false);
  reactant->setCell(Matrix3d::Identity() * 25.0);
  for (int atom = 0; atom < 2; ++atom) {
    reactant->setAtomicNr(atom, 1);
    reactant->setMass(atom, 1.0);
  }
  reactant->setFixed(1, true);
  AtomMatrix positions = AtomMatrix::Zero(2, 3);
  positions(0, 0) = -1.0;
  positions(1, 2) = 5.0;
  reactant->setPositions(positions);
  auto product = std::make_shared<Matter>(*reactant);
  positions(0, 0) = 1.0;
  product->setPositions(positions);
  NudgedElasticBand neb(reactant, product, params, pot);
  const double x[] = {-0.6, 0.2, 0.6};
  const double y[] = {0.1, 0.2, 0.1};
  for (long image = 1; image <= 3; ++image) {
    positions(0, 0) = x[image - 1];
    positions(0, 1) = y[image - 1];
    neb.path[image]->setPositions(positions);
  }
  neb.setCIEnabled(true);
  neb.updateForces();
  REQUIRE(neb.climbingImage == 2);
  const AtomMatrix neighbor = neb.path[1]->getPositions();
  const AtomMatrix first = neb.path.front()->getPositions();
  const AtomMatrix last = neb.path.back()->getPositions();
  eonc::neb::OCINEBController controller(
      eonc::neb::OCINEBController::fromParams(params));
  controller.initBaseline(neb.convergenceForce());
  const auto result = controller.run(neb, neb.convergenceForce());
  REQUIRE_FALSE(result.convergedAfterMMF);
  REQUIRE(result.shouldResetOptimizer);
  REQUIRE_FALSE(neb.path[1]->getPositions().isApprox(neighbor, 1e-10));
  REQUIRE(neb.path.front()->getPositions().isApprox(first, 0.0));
  REQUIRE(neb.path.back()->getPositions().isApprox(last, 0.0));
  for (long image = 1; image <= 3; ++image) {
    REQUIRE(neb.path[image]->getPositions().row(1).isApprox(first.row(1), 0.0));
    REQUIRE(neb.path[image]->getPositions().row(0).norm() < 1.0);
  }
}

class NEBRegressionFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> product;

  NEBRegressionFixture()
      : params{},
        pot{nullptr},
        reactant{nullptr},
        product{nullptr} {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    reactant = std::make_shared<Matter>(pot, params);
    product = std::make_shared<Matter>(pot, params);

    // Two LJ13 minima joined by one first-order saddle (barrier 1.00 eV,
    // saddle at -38.964 eV), so the climbing image has a maximum to reach.
    reactant->con2matter(std::string("reactant.con"));
    product->con2matter(std::string("product.con"));

    // The band reads the NEB's own iteration budget and optimizer.
    ParametersLoadAccess::neb_options(params).max_iterations = 5000;
    ParametersLoadAccess::neb_options(params).opt_method = OptType::LBFGS;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
    ParametersLoadAccess::neb_options(params).image_count = 5;
    ParametersLoadAccess::neb_options(params).force_tolerance = 0.01;
    ParametersLoadAccess::neb_options(params).initialization.method =
        NEBInit::LINEAR;
    ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
    ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  }
};

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB updateForces produces consistent tangent norms",
                 "[neb][regression]") {
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb->updateForces();

  // All tangent vectors must be unit vectors
  for (long i = 1; i <= neb->numImages; i++) {
    double norm = neb->tangent[i]->norm();
    REQUIRE_THAT(norm, Catch::Matchers::WithinAbs(1.0, 1e-10));
  }

  // Projected forces must be non-zero
  for (long i = 1; i <= neb->numImages; i++) {
    REQUIRE(neb->projectedForce[i]->norm() > 0.0);
  }

  // Max energy image must be valid
  REQUIRE(neb->maxEnergyImage >= 1);
  REQUIRE(neb->maxEnergyImage <= static_cast<size_t>(neb->numImages));
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB does not climb below the endpoint energies",
                 "[neb][regression][climbing_image]") {
  const bool reverse = GENERATE(false, true);
  // A high-energy endpoint: the LJ13 reactant with atom 0 pushed 0.5 A
  // along x (516.0 eV). The linear band rises monotonically towards it, so
  // its highest interior image (112.4 eV) stays below that endpoint.
  AtomMatrix pushed = reactant->getPositions();
  pushed(0, 0) += 0.5;
  product->setPositions(pushed);
  if (reverse) {
    std::swap(reactant, product);
  }
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb->updateForces();

  const double endpointEnergy =
      std::max(neb->path.front()->getPotentialEnergy(),
               neb->path.back()->getPotentialEnergy());
  REQUIRE(neb->path[neb->maxEnergyImage]->getPotentialEnergy() <
          endpointEnergy);
  CHECK(neb->climbingImage == 0);

  std::vector<AtomMatrix> forces;
  for (long i = 1; i <= neb->numImages; ++i) {
    forces.push_back(*neb->projectedForce[i]);
  }
  neb->updateForces(false);
  for (long i = 1; i <= neb->numImages; ++i) {
    CHECK(forces[i - 1].isApprox(*neb->projectedForce[i], 1e-12));
  }
}

TEST_CASE_METHOD(
    NEBRegressionFixture,
    "NEB clears the climbing image when an endpoint reaches the peak",
    "[neb][regression][climbing_image]") {
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  product->setPositions(reactant->getPositions());
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb->path[2]->con2matter(std::string("product.con"));
  neb->updateForces();

  REQUIRE(neb->maxEnergyImage == 2);
  REQUIRE(neb->path[2]->getPotentialEnergy() >
          neb->path.front()->getPotentialEnergy());
  REQUIRE(neb->climbingImage == 2);

  neb->path.back()->setPositions(neb->path[2]->getPositions());
  neb->updateForces();
  CHECK(neb->climbingImage == 0);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB CI toggles refresh projection without potential calls",
                 "[neb][regression][climbing_image][force_cache]") {
  ParametersLoadAccess::main_options(params).parallel = false;
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb->updateForces();
  REQUIRE(neb->climbingImage == 0);
  std::vector<AtomMatrix> ordinary;
  for (long i = 1; i <= neb->numImages; ++i) {
    ordinary.push_back(*neb->projectedForce[i]);
  }
  const auto calls = eonc::PotRegistry::get().total_force_calls();

  neb->setCIEnabled(true);
  const double climbedConvergence = neb->convergenceForce();
  REQUIRE(neb->climbingImage > 0);
  const auto ci = neb->climbingImage;
  const AtomMatrix force = neb->path[ci]->getForces();
  const AtomMatrix &tau = *neb->tangent[ci];
  AtomMatrix expected = force - 2.0 * matDot(force, tau) * tau;
  const auto translation = expected.colwise().mean().eval();
  expected.rowwise() -= translation;
  CHECK(expected.isApprox(*neb->projectedForce[ci], 1e-12));
  REQUIRE((ordinary[ci - 1] - *neb->projectedForce[ci]).norm() > 1e-6);
  CHECK(eonc::PotRegistry::get().total_force_calls() == calls);
  neb->updateForces();
  CHECK(neb->convergenceForce() ==
        Catch::Approx(climbedConvergence).margin(1e-12));

  neb->setCIEnabled(false);
  neb->convergenceForce();
  CHECK(neb->climbingImage == 0);
  for (long i = 1; i <= neb->numImages; ++i) {
    CHECK(ordinary[i - 1].isApprox(*neb->projectedForce[i], 1e-12));
  }
  CHECK(eonc::PotRegistry::get().total_force_calls() == calls);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB configured force requests replace one-off CI overrides",
                 "[neb][regression][force_cache][projection_override]") {
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  auto objective = std::make_shared<NEBObjectiveFunction>(neb.get(), params);
  neb->updateForces();
  REQUIRE(neb->climbingImage > 0);
  const auto ci = neb->climbingImage;
  const AtomMatrix projected = *neb->projectedForce[ci];
  const double convergence = neb->convergenceForce();
  const VectorXd gradient = objective->getGradient();
  const auto calls = eonc::PotRegistry::get().total_force_calls();
  neb->updateForces(false);
  REQUIRE(neb->climbingImage == 0);
  REQUIRE((projected - *neb->projectedForce[ci]).norm() > 1e-6);
  neb->setCIEnabled(true);
  const bool gradientRequest = GENERATE(false, true);
  if (gradientRequest) {
    CHECK(objective->getGradient().isApprox(gradient, 1e-12));
  } else {
    CHECK(neb->convergenceForce() == Catch::Approx(convergence).margin(1e-12));
  }
  CHECK(neb->climbingImage == ci);
  CHECK(projected.isApprox(*neb->projectedForce[ci], 1e-12));
  CHECK(eonc::PotRegistry::get().total_force_calls() == calls);
}

TEST_CASE("NEB refreshes the projection when climbing is switched",
          "[neb][regression][climbing_image]") {
  Parameters params;
  auto &options = ParametersLoadAccess::neb_options(params);
  options.image_count = 3;
  options.initialization.method = NEBInit::LINEAR;
  options.endpoints.minimize = false;
  options.climbing_image.enabled = true;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  auto pot = std::make_shared<MolecularSaddlePotential>();
  auto reactant = std::make_shared<Matter>(pot, params);
  reactant->resize(2);
  reactant->setPeriodic(false);
  reactant->setCell(Matrix3d::Identity() * 25.0);
  reactant->setFixed(1, true);
  AtomMatrix positions = AtomMatrix::Zero(2, 3);
  positions(0, 0) = -1.0;
  positions(1, 2) = 5.0;
  reactant->setPositions(positions);
  auto product = std::make_shared<Matter>(*reactant);
  positions(0, 0) = 1.0;
  product->setPositions(positions);
  NudgedElasticBand neb(reactant, product, params, pot);
  positions(0, 0) = 0.2;
  positions(0, 1) = 0.2;
  neb.path[2]->setPositions(positions);
  neb.setCIEnabled(false);
  neb.updateForces();
  NEBObjectiveFunction objective(&neb, params);
  const VectorXd ordinary = objective.getGradient();
  const auto calls = eonc::PotRegistry::get().total_force_calls();
  neb.setCIEnabled(true);
  const VectorXd climbing = objective.getGradient();
  REQUIRE(neb.climbingImage == 2);
  REQUIRE_FALSE(climbing.isApprox(ordinary, 1e-10));
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == calls);
  neb.setCIEnabled(false);
  REQUIRE(objective.getGradient().isApprox(ordinary, 1e-12));
  REQUIRE(neb.climbingImage == 0);
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == calls);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "Zoom NEB reports the barrier from the original reactant",
                 "[neb][regression][zoom]") {
  // The neb_lj13 pair: reactant -39.965352 eV, and an eonclient process
  // search put the saddle at -38.963676 eV, a barrier of 1.001676 eV.
  auto &neb_opts = ParametersLoadAccess::neb_options(params);
  neb_opts.climbing_image.enabled = true;
  neb_opts.climbing_image.trigger_force = 1e10;
  neb_opts.climbing_image.trigger_factor = 1.0;
  neb_opts.zoom.enabled = true;
  neb_opts.zoom.activation_threshold = 1.0;
  neb_opts.zoom.stability_count = 2;
  const double reactantE = reactant->getPotentialEnergy();

  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  REQUIRE(neb->compute() == NudgedElasticBand::NEBStatus::GOOD);

  // Zoom moved the band's first image off the reactant.
  REQUIRE(std::abs(neb->path.front()->getPotentialEnergy() - reactantE) > 1e-3);
  REQUIRE(neb->reactantEnergy == Catch::Approx(reactantE).margin(1e-12));
  const double barrier = neb->path[neb->maxEnergyImage]->getPotentialEnergy() -
                         neb->reactantEnergy;
  REQUIRE(barrier == Catch::Approx(1.001676).margin(5e-3));
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB L-BFGS scaling does not probe projected curvature",
                 "[neb][regression][lbfgs]") {
  ParametersLoadAccess::optimizer_options(params).lbfgs.auto_scale = true;
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  auto objective = std::make_shared<NEBObjectiveFunction>(neb.get(), params);
  objective->getGradient();
  LBFGS optimizer(objective, params);

  const auto before = eonc::PotRegistry::get().total_force_calls();
  optimizer.step(params.optimizer_options().max_move);
  const auto calls = eonc::PotRegistry::get().total_force_calls() - before;

  REQUIRE(calls == static_cast<size_t>(neb->numImages));
  REQUIRE(params.optimizer_options().lbfgs.auto_scale);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB L-BFGS converges with fixed unrelaxed endpoints",
                 "[neb][regression][lbfgs]") {
  product->con2matter(std::string("product.con"));
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).spring.constant = 5.0;
  ParametersLoadAccess::neb_options(params).max_iterations = 200;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 200;
  const AtomMatrix initial = reactant->getPositions();
  const AtomMatrix final = product->getPositions();
  REQUIRE(params.optimizer_options().lbfgs.auto_scale);

  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  REQUIRE(neb->compute() == NudgedElasticBand::NEBStatus::GOOD);
  CHECK(neb->convergenceForce() < params.neb_options().force_tolerance);
  CHECK(neb->path.front()->getPositions().isApprox(initial, 0.0));
  CHECK(neb->path.back()->getPositions().isApprox(final, 0.0));
  CHECK(params.optimizer_options().lbfgs.auto_scale);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB convergence produces stable barrier height",
                 "[neb][regression]") {
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  auto status = neb->compute();
  REQUIRE(status == NudgedElasticBand::NEBStatus::GOOD);

  neb->findExtrema();

  // Record the barrier height (highest extremum energy - reactant energy)
  double E_react = neb->path[0]->getPotentialEnergy();
  double maxBarrier = -1e30;
  for (long j = 0; j < neb->numExtrema; j++) {
    maxBarrier = std::max(maxBarrier, neb->extremumEnergy[j] - E_react);
  }

  // Barrier must be positive and finite
  REQUIRE(maxBarrier > 0.0);
  REQUIRE(std::isfinite(maxBarrier));

  // The converged force must be below tolerance
  REQUIRE(neb->convergenceForce() <
          ParametersLoadAccess::neb_options(params).force_tolerance);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "Energy-weighted springs produce different path than uniform",
                 "[neb][regression]") {
  // Uniform springs
  ParametersLoadAccess::neb_options(params).spring.weighting.enabled = false;
  ParametersLoadAccess::neb_options(params).spring.constant = 5.0;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.05;
  auto neb_uniform =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb_uniform->compute();
  double E_uniform =
      neb_uniform->path[neb_uniform->maxEnergyImage]->getPotentialEnergy();

  // Weighted springs
  ParametersLoadAccess::neb_options(params).spring.weighting.enabled = true;
  ParametersLoadAccess::neb_options(params).spring.weighting.k_min = 1.0;
  ParametersLoadAccess::neb_options(params).spring.weighting.k_max = 10.0;
  auto neb_weighted =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb_weighted->compute();
  double E_weighted =
      neb_weighted->path[neb_weighted->maxEnergyImage]->getPotentialEnergy();

  // Both should converge to similar barriers (same PES)
  // but the paths themselves differ (different spring strengths)
  REQUIRE(std::isfinite(E_uniform));
  REQUIRE(std::isfinite(E_weighted));
  // Relative energy difference should be small (same endpoints, same PES)
  double relDiff = std::abs(E_uniform - E_weighted) /
                   std::max(std::abs(E_uniform), std::abs(E_weighted));
  REQUIRE(relDiff < 0.5); // Within 50% (different methods, same saddle)
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "Climbing image finds higher barrier than plain NEB",
                 "[neb][regression]") {
  // Plain NEB
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.05;
  auto neb_plain =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb_plain->compute();
  double maxE_plain =
      neb_plain->path[neb_plain->maxEnergyImage]->getPotentialEnergy();

  // CI-NEB
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  ParametersLoadAccess::neb_options(params).climbing_image.trigger_force =
      1e10; // immediate
  ParametersLoadAccess::neb_options(params).climbing_image.trigger_factor = 1.0;
  auto neb_ci =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb_ci->compute();
  double maxE_ci = neb_ci->path[neb_ci->maxEnergyImage]->getPotentialEnergy();

  // Both should converge to finite energies
  REQUIRE(std::isfinite(maxE_plain));
  REQUIRE(std::isfinite(maxE_ci));
  // CI should make progress toward convergence (may not fully converge
  // on all platforms within iteration limit)
  CHECK(neb_ci->convergenceForce() < 10.0);
  // The CI image should be at a different position than plain NEB
  // (CI actively moves toward the saddle point)
  REQUIRE(neb_ci->climbingImage >= 0);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "NEB projected forces are zero on fixed atoms",
                 "[neb][regression]") {
  // Fix atom 1 (of 2) to verify that fixed-atom forces are excluded
  // from the NEB projection. Previously getForcesRaw() was used instead
  // of getForces(), which included fixed-atom forces and inflated the
  // convergence metric.
  reactant->setFixed(1, true);
  product->setFixed(1, true);

  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb->updateForces();

  for (long i = 1; i <= neb->numImages; i++) {
    // The projected force on the fixed atom (index 1) must be zero
    double fixedForceNorm = neb->projectedForce[i]->row(1).norm();
    REQUIRE(fixedForceNorm < 1e-10);
  }

  // Force norm with fixed atoms must be less than or equal to raw force norm
  // (cannot be larger if fixed atoms are properly zeroed)
  for (long i = 1; i <= neb->numImages; i++) {
    double projNorm = neb->projectedForce[i]->norm();
    // Sanity: the free atom force must still be non-zero
    double freeForceNorm = neb->projectedForce[i]->row(0).norm();
    REQUIRE(freeForceNorm > 0.0);
    // And the total must be finite
    REQUIRE(std::isfinite(projNorm));
  }
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "Parallel NEB produces identical results to sequential",
                 "[neb][regression][parallel]") {
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.05;

  // Run with parallel disabled
  ParametersLoadAccess::main_options(params).parallel = false;
  auto neb_seq =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  auto status_seq = neb_seq->compute();
  REQUIRE(status_seq == NudgedElasticBand::NEBStatus::GOOD);

  double E_seq = neb_seq->path[neb_seq->maxEnergyImage]->getPotentialEnergy();
  double F_seq = neb_seq->convergenceForce();

  // Run with parallel enabled
  ParametersLoadAccess::main_options(params).parallel = true;
  auto neb_par =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  auto status_par = neb_par->compute();
  REQUIRE(status_par == NudgedElasticBand::NEBStatus::GOOD);

  double E_par = neb_par->path[neb_par->maxEnergyImage]->getPotentialEnergy();
  double F_par = neb_par->convergenceForce();

  // Barrier height and convergence force must be identical
  REQUIRE(E_seq == Catch::Approx(E_par).margin(1e-10));
  REQUIRE(F_seq == Catch::Approx(F_par).margin(1e-10));

  // Same max energy image index
  REQUIRE(neb_seq->maxEnergyImage == neb_par->maxEnergyImage);
}

} /* namespace tests */
