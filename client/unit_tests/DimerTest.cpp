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

#include "Dimer.h"
#include "Davidson.h"
#include "DimerRotationDispatch.h"
#include "EigenmodeStrategy.h"
#include "ImprovedDimer.h"
#include "LORRotation.h"
#include "Lanczos.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <algorithm>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class DimerFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;
  AtomMatrix mode;

  DimerFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    params.potential_options.potential = PotType::LJ;
    params.optimizer_options.method = OptType::CG;
    params.optimizer_options.converged_force = 0.001;
    params.dimer_options.converged_angle = 0.01;
    params.dimer_options.max_iterations = 50;
    params.saddle_search_options.minmode_method =
        LowestEigenmode::MINMODE_DIMER;

    pot = eonc::helpers::makePotential(PotType::LJ, params);
    matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("reactant.con"));

    // Displace one atom slightly to break symmetry and create a saddle region
    auto pos = matter->getPositions();
    pos(0, 0) += 0.3;
    pos(0, 1) -= 0.2;
    matter->setPositions(pos);

    // Random initial direction
    long nAtoms = matter->numberOfAtoms();
    mode = AtomMatrix::Random(nAtoms, 3);
    mode.normalize();
  }
};

// --- Classic Dimer tests ---

TEST_CASE_METHOD(DimerFixture,
                 "Classic Dimer computes eigenvalue on displaced cluster",
                 "[dimer][eigenmode]") {
  params.dimer_options.improved = false;
  auto dimer = std::make_unique<Dimer>(matter, params, pot);
  dimer->compute(matter, mode);

  double eigenvalue = dimer->getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
  REQUIRE(eigenvalue < 0.0);
}

// --- ImprovedDimer tests ---

TEST_CASE_METHOD(DimerFixture, "ImprovedDimer computes negative eigenvalue",
                 "[dimer][eigenmode][improved]") {
  params.dimer_options.improved = true;
  auto dimer = std::make_unique<ImprovedDimer>(matter, params, pot);
  dimer->compute(matter, mode);

  double eigenvalue = dimer->getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
  REQUIRE(eigenvalue < 0.0);
}

// --- Lanczos tests ---

TEST_CASE_METHOD(DimerFixture, "Lanczos computes negative eigenvalue",
                 "[lanczos][eigenmode]") {
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_LANCZOS;
  auto lanczos = std::make_unique<Lanczos>(matter, params, pot);
  lanczos->compute(matter, mode);

  double eigenvalue = lanczos->getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
  REQUIRE(eigenvalue < 0.0);
}

TEST_CASE_METHOD(DimerFixture, "Davidson computes negative eigenvalue",
                 "[davidson][eigenmode]") {
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_DAVIDSON;
  auto davidson = std::make_unique<Davidson>(matter, params, pot);
  davidson->compute(matter, mode);

  double eigenvalue = davidson->getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
  REQUIRE(eigenvalue < 0.0);
}

TEST_CASE_METHOD(DimerFixture, "Davidson and Lanczos agree on lowest mode sign",
                 "[davidson][lanczos][eigenmode]") {
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_LANCZOS;
  Lanczos lanczos(matter, params, pot);
  lanczos.compute(matter, mode);
  const double ewL = lanczos.getEigenvalue();

  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_DAVIDSON;
  Davidson davidson(matter, params, pot);
  davidson.compute(matter, mode);
  const double ewD = davidson.getEigenvalue();

  REQUIRE(std::isfinite(ewL));
  REQUIRE(std::isfinite(ewD));
  // Both should find a negative curvature direction on this saddle-ish LJ
  // setup.
  REQUIRE(ewL < 0.0);
  REQUIRE(ewD < 0.0);
  // Magnitudes within a loose factor (FD noise + method differences).
  REQUIRE(std::fabs(ewD - ewL) <
          0.5 * (std::fabs(ewL) + std::fabs(ewD) + 1e-6));
}

// --- EigenmodeStrategy variant dispatch tests ---

TEST_CASE_METHOD(DimerFixture,
                 "buildEigenmodeStrategy returns ImprovedDimer by default",
                 "[eigenmode][strategy]") {
  params.dimer_options.improved = true;
  params.saddle_search_options.minmode_method = LowestEigenmode::MINMODE_DIMER;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  REQUIRE(strategy != nullptr);
  REQUIRE(std::holds_alternative<ImprovedDimer>(*strategy));
}

TEST_CASE_METHOD(DimerFixture,
                 "buildEigenmodeStrategy returns Dimer when improved=false",
                 "[eigenmode][strategy]") {
  params.dimer_options.improved = false;
  params.saddle_search_options.minmode_method = LowestEigenmode::MINMODE_DIMER;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  REQUIRE(strategy != nullptr);
  REQUIRE(std::holds_alternative<Dimer>(*strategy));
}

TEST_CASE_METHOD(DimerFixture, "buildEigenmodeStrategy returns Lanczos variant",
                 "[eigenmode][strategy]") {
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_LANCZOS;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  REQUIRE(strategy != nullptr);
  REQUIRE(std::holds_alternative<Lanczos>(*strategy));
}

TEST_CASE_METHOD(DimerFixture,
                 "buildEigenmodeStrategy returns Davidson variant",
                 "[eigenmode][strategy][davidson]") {
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_DAVIDSON;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  REQUIRE(strategy != nullptr);
  REQUIRE(std::holds_alternative<Davidson>(*strategy));
}

// --- asImprovedDimer tests ---

TEST_CASE_METHOD(DimerFixture,
                 "asImprovedDimer returns non-null for ImprovedDimer variant",
                 "[eigenmode][strategy]") {
  params.dimer_options.improved = true;
  params.saddle_search_options.minmode_method = LowestEigenmode::MINMODE_DIMER;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  auto *ptr = eonc::asImprovedDimer(*strategy);
  REQUIRE(ptr != nullptr);
}

TEST_CASE_METHOD(DimerFixture,
                 "asImprovedDimer returns null for classic Dimer variant",
                 "[eigenmode][strategy]") {
  params.dimer_options.improved = false;
  params.saddle_search_options.minmode_method = LowestEigenmode::MINMODE_DIMER;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  auto *ptr = eonc::asImprovedDimer(*strategy);
  REQUIRE(ptr == nullptr);
}

// --- Variant dispatch compute test ---

TEST_CASE_METHOD(DimerFixture,
                 "eigenmodeCompute dispatches correctly via variant",
                 "[eigenmode][strategy]") {
  params.dimer_options.improved = true;
  params.saddle_search_options.minmode_method = LowestEigenmode::MINMODE_DIMER;

  auto strategy = eonc::buildEigenmodeStrategy(matter, params, pot);
  eonc::eigenmodeCompute(*strategy, matter, mode);

  double eigenvalue = eonc::eigenmodeGetEigenvalue(*strategy);
  REQUIRE(std::isfinite(eigenvalue));
  REQUIRE(eigenvalue < 0.0);

  AtomMatrix eigenvector = eonc::eigenmodeGetEigenvector(*strategy);
  REQUIRE(eigenvector.rows() == matter->numberOfAtoms());
  REQUIRE(eigenvector.cols() == 3);
  REQUIRE(eigenvector.norm() > 0.0);
}

// --- Fixed-atom tests ---

class DimerFixedAtomFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;
  AtomMatrix mode;

  DimerFixedAtomFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    params.potential_options.potential = PotType::LJ;
    params.optimizer_options.method = OptType::CG;
    params.optimizer_options.converged_force = 0.001;
    params.dimer_options.converged_angle = 0.01;
    params.dimer_options.max_iterations = 50;
    params.saddle_search_options.minmode_method =
        LowestEigenmode::MINMODE_DIMER;

    pot = eonc::helpers::makePotential(PotType::LJ, params);
    matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("reactant.con"));

    // Fix all atoms except the first
    for (int i = 1; i < matter->numberOfAtoms(); i++) {
      matter->setFixed(i, true);
    }

    // Displace the free atom
    auto pos = matter->getPositions();
    pos(0, 0) += 0.3;
    pos(0, 1) -= 0.2;
    matter->setPositions(pos);

    long nAtoms = matter->numberOfAtoms();
    mode = AtomMatrix::Random(nAtoms, 3);
    mode.normalize();
  }
};

TEST_CASE_METHOD(DimerFixedAtomFixture,
                 "ImprovedDimer forces on fixed atoms are zero",
                 "[dimer][fixed_atoms]") {
  params.dimer_options.improved = true;
  auto dimer = std::make_unique<ImprovedDimer>(matter, params, pot);
  dimer->compute(matter, mode);

  // After dimer computation, forces on fixed atoms must be zero
  const AtomMatrix &forces = matter->getForces();
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    REQUIRE(forces.row(i).norm() < 1e-10);
  }
  // Free atom must have non-zero force
  REQUIRE(forces.row(0).norm() > 0.0);
}

TEST_CASE_METHOD(DimerFixedAtomFixture,
                 "Dimer eigenvalue is finite with fixed atoms",
                 "[dimer][fixed_atoms]") {
  params.dimer_options.improved = true;
  auto dimer = std::make_unique<ImprovedDimer>(matter, params, pot);
  dimer->compute(matter, mode);

  double eigenvalue = dimer->getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
}

// --- LOR (Leng et al. JCP 2013) rotation backend ---

TEST_CASE_METHOD(DimerFixture, "LOR rotation finds finite lowest curvature",
                 "[dimer][lor][eigenmode]") {
  params.dimer_options.improved = true;
  params.dimer_options.rotation_backend = DimerRotationBackend::LOR;
  params.dimer_options.max_iterations = 20;
  params.dimer_options.rotations_max = 20;
  auto lor = std::make_unique<LORRotation>(matter, params, pot);
  lor->compute(matter, mode);

  double ev = lor->getEigenvalue();
  REQUIRE(std::isfinite(ev));
  REQUIRE(lor->totalForceCalls > 0);
  REQUIRE(lor->statsRotations >= 0);
  // Softest mode on displaced LJ (same fixture as classic/improved dimer).
  REQUIRE(ev < 0.0);
}

TEST_CASE_METHOD(DimerFixture,
                 "LOR curvature history is non-increasing within tolerance",
                 "[dimer][lor][eigenmode]") {
  params.dimer_options.improved = true;
  params.dimer_options.max_iterations = 20;
  params.dimer_options.rotations_max = 20;
  // Start from Lanczos softest mode so LOR operates in a negative-C basin.
  Lanczos lanczos(matter, params, pot);
  lanczos.compute(matter, mode);
  params.dimer_options.rotation_backend = DimerRotationBackend::LOR;
  LORRotation lor(matter, params, pot);
  lor.compute(matter, lanczos.getEigenvector());

  REQUIRE_FALSE(lor.curvatureHistory.empty());
  // Production appendHistory: cn <= previous + 1e-4 (float margin 1e-6).
  for (size_t i = 1; i < lor.curvatureHistory.size(); ++i) {
    REQUIRE(lor.curvatureHistory[i] <=
            lor.curvatureHistory[i - 1] + 1e-4 + 1e-6);
  }
  // Returned eigenpair matches best accepted history sample.
  const double histMin = *std::min_element(lor.curvatureHistory.begin(),
                                           lor.curvatureHistory.end());
  REQUIRE(lor.getEigenvalue() == Catch::Approx(histMin).margin(1e-4));
  REQUIRE(lor.getEigenvalue() < 0.0);
  REQUIRE(lor.getEigenvector().norm() == Catch::Approx(1.0).margin(1e-6));
}

TEST_CASE_METHOD(
    DimerFixture,
    "LOR mode agrees with classical ImprovedDimer (sign-insensitive)",
    "[dimer][lor][eigenmode]") {
  params.dimer_options.improved = true;
  params.dimer_options.max_iterations = 50;
  params.dimer_options.rotations_max = 30;

  // Shared deterministic seed (not Lanczos output — avoids self-agreement).
  AtomMatrix seed = AtomMatrix::Zero(matter->numberOfAtoms(), 3);
  seed(0, 0) = 1.0;
  seed.normalize();

  params.dimer_options.rotation_backend = DimerRotationBackend::Classical;
  ImprovedDimer classical(matter, params, pot);
  classical.compute(matter, seed);
  AtomMatrix mClass = classical.getEigenvector();

  Lanczos lanczos(matter, params, pot);
  lanczos.compute(matter, seed);
  AtomMatrix mLanc = lanczos.getEigenvector();

  params.dimer_options.rotation_backend = DimerRotationBackend::LOR;
  LORRotation lor(matter, params, pot);
  lor.compute(matter, seed);
  AtomMatrix mLor = lor.getEigenvector();

  auto absCos = [](const AtomMatrix &a, const AtomMatrix &b) {
    const double n1 = a.norm();
    const double n2 = b.norm();
    if (n1 < 1e-10 || n2 < 1e-10) {
      return 0.0;
    }
    return std::abs((a.array() * b.array()).sum() / (n1 * n2));
  };
  const double cosClass = absCos(mClass, mLor);
  const double cosLanc = absCos(mLanc, mLor);
  // Must agree with classical or Lanczos softest mode (|cos| > 0.7).
  REQUIRE(std::max(cosClass, cosLanc) > 0.7);
  REQUIRE(lor.getEigenvalue() < 0.0);
  REQUIRE(classical.getEigenvalue() < 0.0);
  REQUIRE(lanczos.getEigenvalue() < 0.0);
}

TEST_CASE_METHOD(DimerFixture,
                 "ImprovedDimer rotation_backend=lor is live path",
                 "[dimer][lor][eigenmode]") {
  params.dimer_options.improved = true;
  params.dimer_options.rotation_backend = DimerRotationBackend::LOR;
  params.dimer_options.rotations_max = 20;
  ImprovedDimer dimer(matter, params, pot);
  dimer.compute(matter, mode);
  REQUIRE(std::isfinite(dimer.getEigenvalue()));
  REQUIRE(dimer.getEigenvalue() < 0.0);
  REQUIRE(dimer.totalForceCalls > 0);
}

TEST_CASE_METHOD(DimerFixture,
                 "rotation_backend lanczos/davidson report force calls",
                 "[dimer][lor][force_calls]") {
  params.dimer_options.improved = true;
  params.dimer_options.rotations_max = 20;

  params.dimer_options.rotation_backend = DimerRotationBackend::Lanczos;
  ImprovedDimer dimL(matter, params, pot);
  dimL.compute(matter, mode);
  REQUIRE(dimL.totalForceCalls > 0);
  REQUIRE(std::isfinite(dimL.getEigenvalue()));

  params.dimer_options.rotation_backend = DimerRotationBackend::Davidson;
  ImprovedDimer dimD(matter, params, pot);
  dimD.compute(matter, mode);
  REQUIRE(dimD.totalForceCalls > 0);
  REQUIRE(std::isfinite(dimD.getEigenvalue()));

  // Direct solvers also expose totalForceCalls.
  Lanczos lan(matter, params, pot);
  lan.compute(matter, mode);
  REQUIRE(lan.totalForceCalls > 0);
  Davidson dav(matter, params, pot);
  dav.compute(matter, mode);
  REQUIRE(dav.totalForceCalls > 0);
}

TEST_CASE_METHOD(DimerFixture, "LOR residual convergence flag via dispatch",
                 "[dimer][lor][convergence]") {
  params.dimer_options.improved = true;
  params.dimer_options.rotation_backend = DimerRotationBackend::LOR;
  // Impossible residual with minimal budget → not convergedOnResidual.
  params.dimer_options.lor_residual_tol = 1e-15;
  params.dimer_options.rotations_max = 2;
  ImprovedDimer tight(matter, params, pot);
  tight.compute(matter, mode);
  REQUIRE_FALSE(tight.rotationDidConverge);

  // Loose residual → converges on residual when budget allows.
  params.dimer_options.lor_residual_tol = 10.0;
  params.dimer_options.rotations_max = 20;
  ImprovedDimer loose(matter, params, pot);
  loose.compute(matter, mode);
  REQUIRE(loose.rotationDidConverge);
  REQUIRE(loose.totalForceCalls > 0);
}

TEST_CASE_METHOD(DimerFixture, "classical rotation_backend skips alt dispatch",
                 "[dimer][lor][dispatch]") {
  params.dimer_options.rotation_backend = DimerRotationBackend::Classical;
  auto none = eonc::runAlternativeRotation(DimerRotationBackend::Classical,
                                           matter, params, pot, mode);
  REQUIRE_FALSE(none.has_value());
}

TEST_CASE_METHOD(DimerFixture, "non-improved Dimer uses LOR rotation_backend",
                 "[dimer][lor][eigenmode]") {
  params.dimer_options.improved = false;
  params.dimer_options.rotation_backend = DimerRotationBackend::LOR;
  params.dimer_options.rotations_max = 20;
  Dimer dimer(matter, params, pot);
  dimer.compute(matter, mode);
  REQUIRE(std::isfinite(dimer.getEigenvalue()));
  REQUIRE(dimer.totalForceCalls > 0);
  REQUIRE(dimer.getEigenvalue() < 0.0);
}

// Pure linear-algebra check of shipped force-translation identity for H·P3
// (Algorithm I H·P reuse). Uses a synthetic SPD H so the oracle is exact H*v,
// but the code under test is LORRotation::translateHUnitOrthoP3 (production).
TEST_CASE("LOR translateHUnitOrthoP3 matches H*P3 for linear H",
          "[dimer][lor][translation]") {
  constexpr int n = 12;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd Hmat = A.transpose() * A + Eigen::MatrixXd::Identity(n, n);

  VectorXd N = VectorXd::Random(n);
  N.normalize();
  VectorXd Theta = VectorXd::Random(n);
  Theta = Theta - N.dot(Theta) * N;
  Theta.normalize();
  VectorXd P = VectorXd::Random(n);

  const VectorXd HN = Hmat * N;
  const VectorXd HTheta = Hmat * Theta;
  const VectorXd HP = Hmat * P;

  const VectorXd P_ortho = P - N.dot(P) * N - Theta.dot(P) * Theta;
  REQUIRE(P_ortho.norm() > 1e-8);
  const VectorXd P3 = P_ortho / P_ortho.norm();
  const VectorXd HP3_exact = Hmat * P3;
  const VectorXd HP3_code =
      LORRotation::translateHUnitOrthoP3(N, Theta, P, HN, HTheta, HP);

  // Wrong ambient GS on HP would give O(1) error on this toy; linearity is ~0.
  const double err = (HP3_exact - HP3_code).norm();
  const VectorXd HP3_wrong_gs =
      HP - N.dot(HP) * N - Theta.dot(HP) * Theta; // skeptic's counterexample
  const double err_wrong = (HP3_exact - HP3_wrong_gs).norm();
  REQUIRE(err < 1e-9);
  REQUIRE(err_wrong > 0.1); // documents that GS-on-HP is not H·P3
}

} /* namespace tests */
