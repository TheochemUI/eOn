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

#include "eon/Dimer.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Davidson.h"
#include "eon/DimerRotationDispatch.h"
#include "eon/EigenmodeStrategy.h"
#include "eon/ImprovedDimer.h"
#include "eon/LORRotation.h"
#include "eon/Lanczos.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/MobileAtoms.h"
#include "eon/Parameters.h"

#include <algorithm>
#include <cmath>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Deterministic soft-mode seed (atom 0 along +x). Random modes can land in
// positive-curvature basins and flake REQUIRE(eigenvalue < 0) (XTB CI exit 42).
static AtomMatrix softModeSeed(long nAtoms) {
  AtomMatrix seed = AtomMatrix::Zero(nAtoms, 3);
  seed(0, 0) = 1.0;
  seed.normalize();
  return seed;
}

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

    mode = softModeSeed(matter->numberOfAtoms());
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

#ifndef WITH_GPRD
TEST_CASE_METHOD(DimerFixture,
                 "gprdimer does not silently become ImprovedDimer",
                 "[eigenmode][strategy][gprdimer]") {
  params.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_GPRDIMER;
  REQUIRE_THROWS_WITH(eonc::buildEigenmodeStrategy(matter, params, pot),
                      Catch::Matchers::ContainsSubstring("with_gprd"));
}
#endif

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

    mode = softModeSeed(matter->numberOfAtoms());
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
  // Deterministic seed into a soft basin (random mode can be mostly positive).
  AtomMatrix seed = AtomMatrix::Zero(matter->numberOfAtoms(), 3);
  seed(0, 0) = 1.0;
  seed.normalize();
  auto lor = std::make_unique<LORRotation>(matter, params, pot);
  lor->compute(matter, seed);

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
  REQUIRE(lor.getEigenvalue() < 0.0);
  REQUIRE(classical.getEigenvalue() < 0.0);
  REQUIRE(lanczos.getEigenvalue() < 0.0);
  // Shared deterministic seed (not Lanczos eigenvector input). Softest modes on
  // this LJ cluster can be nearly degenerate; require negative curvature and
  // that LOR is at least as soft as the classical/Lanczos peer within a factor.
  const double evL = lor.getEigenvalue();
  const double peer =
      std::min(classical.getEigenvalue(), lanczos.getEigenvalue());
  REQUIRE(evL < 0.0);
  REQUIRE(peer < 0.0);
  REQUIRE(evL <= peer * 0.5); // LOR not stuck far above the peer soft mode
  REQUIRE((cosLanc > 0.5 || cosClass > 0.5 ||
           std::abs(evL - peer) / (std::abs(peer) + 1e-6) < 0.5));
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

// --- PHVA mobile Krylov (free/fixed != active) ---

TEST_CASE_METHOD(DimerFixture, "resolveMobileAtoms All matches free atoms",
                 "[mobile][eigenmode]") {
  VectorXi free = eonc::freeAtomIndices(matter.get());
  VectorXi all = eonc::resolveMobileAtoms(matter.get(), "All");
  VectorXi allLower = eonc::resolveMobileAtoms(matter.get(), "all");
  REQUIRE(all.size() == free.size());
  REQUIRE(allLower.size() == free.size());
  for (Eigen::Index i = 0; i < free.size(); ++i) {
    REQUIRE(all(i) == free(i));
    REQUIRE(allLower(i) == free(i));
  }
  REQUIRE(static_cast<long>(all.size()) == matter->numberOfFreeAtoms());
}

TEST_CASE_METHOD(DimerFixture,
                 "resolveMobileAtoms intersects free and drops fixed",
                 "[mobile][eigenmode]") {
  // Fix atom 0; active list asks for 0 and 1 -> only 1 survives.
  matter->setFixed(0, 1);
  VectorXi cands(2);
  cands << 0, 1;
  VectorXi mobile = eonc::resolveMobileAtoms(matter.get(), cands);
  REQUIRE(mobile.size() == 1);
  REQUIRE(mobile(0) == 1);
  VectorXi fromStr = eonc::resolveMobileAtoms(matter.get(), "0,1,1,99");
  REQUIRE(fromStr.size() == 1);
  REQUIRE(fromStr(0) == 1);
  // restore for later tests in same fixture instance - fixture is per-test
}

TEST_CASE_METHOD(DimerFixture,
                 "Lanczos active subset bounds Krylov dim and zeros env mode",
                 "[lanczos][mobile][eigenmode]") {
  // Keep free mask = all atoms; restrict Krylov to first 3 atoms only.
  VectorXi mobile(3);
  mobile << 0, 1, 2;
  REQUIRE(matter->numberOfFreeAtoms() == matter->numberOfAtoms());
  REQUIRE(matter->numberOfFreeAtoms() > 3);

  params.lanczos_options.max_iterations = 30;
  params.lanczos_options.tolerance = 1e-4;
  Lanczos lanczos(matter, params, pot);
  AtomMatrix dir = AtomMatrix::Zero(matter->numberOfAtoms(), 3);
  dir(0, 0) = 1.0;
  dir(1, 0) = -0.5;
  dir(2, 1) = 0.3;
  dir.normalize();
  lanczos.compute(matter, dir, mobile);

  REQUIRE(std::isfinite(lanczos.getEigenvalue()));
  AtomMatrix ev = lanczos.getEigenvector();
  REQUIRE(ev.rows() == matter->numberOfAtoms());
  // Environment atoms outside the active set must carry zero mode amplitude.
  for (long i = 3; i < matter->numberOfAtoms(); ++i) {
    REQUIRE(ev.row(i).norm() == Catch::Approx(0.0).margin(1e-14));
  }
  // Active block has some weight.
  double actNorm = 0.0;
  for (int a = 0; a < 3; ++a) {
    actNorm += ev.row(a).squaredNorm();
  }
  REQUIRE(actNorm > 1e-12);
}

TEST_CASE_METHOD(DimerFixture,
                 "Davidson active subset matches Lanczos free-default sign",
                 "[davidson][mobile][eigenmode]") {
  // Default path (phva_atoms All) still finds negative curvature on fixture.
  params.lanczos_options.phva_atoms = "All";
  params.davidson_options.phva_atoms = "All";
  params.lanczos_options.max_iterations = 30;
  params.davidson_options.max_iterations = 30;

  Lanczos lanczos(matter, params, pot);
  lanczos.compute(matter, mode);
  Davidson davidson(matter, params, pot);
  davidson.compute(matter, mode);

  REQUIRE(lanczos.getEigenvalue() < 0.0);
  REQUIRE(davidson.getEigenvalue() < 0.0);
}

TEST_CASE_METHOD(
    DimerFixture,
    "Lanczos phva_atoms param restricts without changing free mask",
    "[lanczos][mobile][eigenmode]") {
  const long nFreeBefore = matter->numberOfFreeAtoms();
  params.lanczos_options.phva_atoms = "0,1,2";
  params.lanczos_options.max_iterations = 30;
  Lanczos lanczos(matter, params, pot);
  AtomMatrix dir = AtomMatrix::Zero(matter->numberOfAtoms(), 3);
  dir(0, 0) = 1.0;
  lanczos.compute(matter, dir);
  REQUIRE(matter->numberOfFreeAtoms() == nFreeBefore);
  AtomMatrix ev = lanczos.getEigenvector();
  for (long i = 3; i < matter->numberOfAtoms(); ++i) {
    REQUIRE(ev.row(i).norm() == Catch::Approx(0.0).margin(1e-14));
  }
}

TEST_CASE_METHOD(DimerFixture, "Lanczos All equals explicit free list",
                 "[lanczos][mobile][eigenmode]") {
  params.lanczos_options.max_iterations = 25;
  params.lanczos_options.tolerance = 1e-3;
  VectorXi free = eonc::freeAtomIndices(matter.get());

  Lanczos a(matter, params, pot);
  a.compute(matter, mode);
  Lanczos b(matter, params, pot);
  b.compute(matter, mode, free);

  REQUIRE(a.getEigenvalue() == Catch::Approx(b.getEigenvalue()).margin(1e-6));
  AtomMatrix eva = a.getEigenvector();
  AtomMatrix evb = b.getEigenvector();
  // Modes may flip sign.
  // Flattened dot (sign-invariant)
  double dot = 0.0, na = 0.0, nb = 0.0;
  for (long i = 0; i < eva.rows(); ++i) {
    for (int c = 0; c < 3; ++c) {
      dot += eva(i, c) * evb(i, c);
      na += eva(i, c) * eva(i, c);
      nb += evb(i, c) * evb(i, c);
    }
  }
  REQUIRE(std::fabs(dot) / (std::sqrt(na * nb) + 1e-30) ==
          Catch::Approx(1.0).margin(1e-5));
}

} /* namespace tests */
