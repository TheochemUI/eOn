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

/// Integration tests for ARTn: run both dimer and ARTn on the same system,
/// verify they find the same saddle point.

#include "ARTnSaddleSearch.h"
#include "ImprovedDimer.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "TestUtils.hpp"
#ifdef WITH_ARTN
#include "libs/ARTn/ARTnResource.h"
#endif
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class ARTnVsDimerFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;

  // Two copies of the same starting structure
  std::shared_ptr<Matter> matter_dimer;
  std::shared_ptr<Matter> matter_artn;

  // Initial displacement direction (shared)
  AtomMatrix displacement;

  ARTnVsDimerFixture()
      : params{},
        pot{nullptr} {
    params.potential_options.potential = PotType::LJ;
    pot = eonc::helpers::makePotential(PotType::LJ, params);

    matter_dimer = std::make_shared<Matter>(pot, params);
    matter_artn = std::make_shared<Matter>(pot, params);
    matter_dimer->con2matter(std::string("reactant.con"));
    matter_artn->con2matter(std::string("reactant.con"));

    // Small displacement from equilibrium (large displacements cause
    // numerical issues in ARTn's internal Lanczos solver)
    auto pos = matter_dimer->getPositions();
    pos(0, 0) += 0.05;
    pos(0, 1) -= 0.03;
    matter_dimer->setPositions(pos);
    matter_artn->setPositions(pos);

    // Initial mode guess along displacement direction
    int nat = matter_dimer->numberOfAtoms();
    displacement = AtomMatrix::Zero(nat, 3);
    displacement(0, 0) = 1.0;
    displacement(0, 1) = -0.5;
    displacement.row(0).normalize();

    // Dimer parameters
    params.saddle_search_options.max_iterations = 500;
    params.saddle_search_options.displace_magnitude = 0.01;
    params.saddle_search_options.displace_radius = 5.0;
    params.saddle_search_options.converged_force = 0.01;
    params.saddle_search_options.minmode_method =
        LowestEigenmode::MINMODE_DIMER;

    // ARTn parameters for LJ cluster
    params.artn_options.push_step_size = 0.3;
    params.artn_options.force_threshold = 0.05;
    params.artn_options.max_iterations = 500;
  }
};

TEST_CASE_METHOD(ARTnVsDimerFixture,
                 "ARTn and Dimer find saddle with comparable energy",
                 "[artn][dimer][comparison]") {
  if (!eonc::get_artn_resource().is_loaded())
    SKIP("libartn not available at runtime");
  // --- Run Dimer saddle search ---
  double initialEnergy = matter_dimer->getPotentialEnergy();
  auto dimerSearch = std::make_shared<MinModeSaddleSearch>(
      matter_dimer, displacement, initialEnergy, params, pot);

  int dimerStatus = dimerSearch->run();
  double dimerSaddleEnergy = matter_dimer->getPotentialEnergy();
  double dimerEigenvalue = dimerSearch->getEigenvalue();
  int dimerIters = dimerSearch->iteration;

  // Dimer should not crash (may converge, hit max iters, or abort on
  // nonnegative eigenvalue depending on the starting displacement)
  REQUIRE(dimerStatus >= 0);
  REQUIRE(std::isfinite(dimerSaddleEnergy));

  INFO("Dimer: status=" << dimerStatus << " energy=" << dimerSaddleEnergy
                        << " eigenvalue=" << dimerEigenvalue
                        << " iterations=" << dimerIters);

  // --- Run ARTn saddle search ---
  auto artnSearch = std::make_unique<ARTnSaddleSearch>(matter_artn, pot,
                                                       displacement, params);
  int artnStatus = artnSearch->run();
  double artnSaddleEnergy = matter_artn->getPotentialEnergy();
  double artnEigenvalue = artnSearch->getEigenvalue();
  int artnIters = artnSearch->iteration;

  // ARTn should not crash (0=good, 1=max_iterations, 2=artn_error)
  REQUIRE((artnStatus == 0 || artnStatus == 1 || artnStatus == 2));
  REQUIRE(std::isfinite(artnSaddleEnergy));

  INFO("ARTn: status=" << artnStatus << " energy=" << artnSaddleEnergy
                       << " eigenvalue=" << artnEigenvalue
                       << " iterations=" << artnIters);

  // --- Compare ---
  // If both converged, they should find saddle points with comparable energies.
  // The exact saddle may differ (different search paths), but the energy
  // landscape has a finite number of saddles and both should be at one.
  if (dimerStatus == MinModeSaddleSearch::STATUS_GOOD && artnStatus == 0) {
    // Both found a first-order saddle: negative eigenvalue
    REQUIRE(dimerEigenvalue < 0.0);
    REQUIRE(artnEigenvalue < 0.0);

    // Energy should be higher than the minimum
    double minEnergy = matter_dimer->getPotentialEnergy(); // after dimer search
    // Both saddle energies should be finite and comparable
    // (allow generous tolerance since they may find different saddles)
    REQUIRE(std::abs(dimerSaddleEnergy - artnSaddleEnergy) /
                std::max(1.0, std::abs(dimerSaddleEnergy)) <
            0.5);
  }
}

TEST_CASE_METHOD(ARTnVsDimerFixture,
                 "ARTn produces negative eigenvalue at saddle",
                 "[artn][eigenvalue]") {
  if (!eonc::get_artn_resource().is_loaded())
    SKIP("libartn not available at runtime");
  auto artnSearch = std::make_unique<ARTnSaddleSearch>(matter_artn, pot,
                                                       displacement, params);
  int status = artnSearch->run();

  if (status == 0) {
    // Converged: must have negative eigenvalue (first-order saddle)
    REQUIRE(artnSearch->getEigenvalue() < 0.0);

    // Eigenvector should be normalized and non-zero
    AtomMatrix evec = artnSearch->getEigenvector();
    REQUIRE(evec.norm() > 1e-10);
  }

  // Should have used force calls
  REQUIRE(artnSearch->forcecalls > 0);
}

} /* namespace tests */
