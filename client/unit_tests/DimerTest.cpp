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

#include "Davidson.h"
#include "Dimer.h"
#include "EigenmodeStrategy.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

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

TEST_CASE_METHOD(DimerFixture,
                 "Davidson and Lanczos agree on lowest mode sign",
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
  // Both should find a negative curvature direction on this saddle-ish LJ setup.
  REQUIRE(ewL < 0.0);
  REQUIRE(ewD < 0.0);
  // Magnitudes within a loose factor (FD noise + method differences).
  REQUIRE(std::fabs(ewD - ewL) < 0.5 * (std::fabs(ewL) + std::fabs(ewD) + 1e-6));
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

} /* namespace tests */
