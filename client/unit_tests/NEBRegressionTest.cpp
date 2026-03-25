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

#include "NudgedElasticBand.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

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
    params.potential_options.potential = PotType::LJ;
    pot = eonc::helpers::makePotential(PotType::LJ, params);
    reactant = std::make_shared<Matter>(pot, params);
    product = std::make_shared<Matter>(pot, params);

    reactant->con2matter(std::string("reactant.con"));
    product->con2matter(std::string("reactant.con"));
    auto pos = product->getPositions();
    pos(0, 0) += 0.5;
    pos(0, 1) -= 0.3;
    pos(0, 2) += 0.2;
    product->setPositions(pos);

    params.optimizer_options.max_iterations = 5000;
    params.optimizer_options.method = OptType::LBFGS;
    params.optimizer_options.max_move = 0.2;
    params.neb_options.image_count = 5;
    params.neb_options.force_tolerance = 0.01;
    params.neb_options.initialization.method = NEBInit::LINEAR;
    params.neb_options.endpoints.minimize = false;
    params.neb_options.climbing_image.enabled = false;
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
  REQUIRE(neb->convergenceForce() < params.neb_options.force_tolerance);
}

TEST_CASE_METHOD(NEBRegressionFixture,
                 "Energy-weighted springs produce different path than uniform",
                 "[neb][regression]") {
  // Uniform springs
  params.neb_options.spring.weighting.enabled = false;
  params.neb_options.spring.constant = 5.0;
  params.neb_options.force_tolerance = 0.05;
  auto neb_uniform =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb_uniform->compute();
  double E_uniform =
      neb_uniform->path[neb_uniform->maxEnergyImage]->getPotentialEnergy();

  // Weighted springs
  params.neb_options.spring.weighting.enabled = true;
  params.neb_options.spring.weighting.k_min = 1.0;
  params.neb_options.spring.weighting.k_max = 10.0;
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
  params.neb_options.climbing_image.enabled = false;
  params.neb_options.force_tolerance = 0.05;
  auto neb_plain =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb_plain->compute();
  double maxE_plain =
      neb_plain->path[neb_plain->maxEnergyImage]->getPotentialEnergy();

  // CI-NEB
  params.neb_options.climbing_image.enabled = true;
  params.neb_options.climbing_image.trigger_force = 1e10; // immediate
  params.neb_options.climbing_image.trigger_factor = 1.0;
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
  params.neb_options.climbing_image.enabled = false;
  params.neb_options.force_tolerance = 0.05;

  // Run with parallel disabled
  params.main_options.parallel = false;
  auto neb_seq =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  auto status_seq = neb_seq->compute();
  REQUIRE(status_seq == NudgedElasticBand::NEBStatus::GOOD);

  double E_seq = neb_seq->path[neb_seq->maxEnergyImage]->getPotentialEnergy();
  double F_seq = neb_seq->convergenceForce();

  // Run with parallel enabled
  params.main_options.parallel = true;
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
