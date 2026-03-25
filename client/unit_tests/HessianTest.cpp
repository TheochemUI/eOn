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

#include "Hessian.h"
#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Hessian on LJ cluster is symmetric", "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->numberOfAtoms() == 13);

  // Use a small subset (3 atoms) to keep the test fast
  VectorXi subAtoms(3);
  subAtoms << 0, 1, 2;

  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), subAtoms);

  REQUIRE(H.rows() == 9);
  REQUIRE(H.cols() == 9);

  for (long i = 0; i < H.rows(); i++) {
    for (long j = 0; j < i; j++) {
      REQUIRE_THAT(H(i, j), Catch::Matchers::WithinAbs(H(j, i), 1e-6));
    }
  }
}

TEST_CASE("Hessian getFreqs returns finite eigenvalues", "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  VectorXi subAtoms(3);
  subAtoms << 0, 1, 2;

  Hessian hess(params, matter.get());
  VectorXd freqs = hess.getFreqs(matter.get(), subAtoms);

  REQUIRE(freqs.size() == 9);
  for (long i = 0; i < freqs.size(); i++) {
    REQUIRE(std::isfinite(freqs(i)));
  }
}

TEST_CASE("Hessian on Pt frozen layers system handles mixed fixed/free",
          "[hessian][morse_pt]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("../Pt_Heptamer_FrozenLayers/pos.con"));

  // Build moved atom list (free atoms only)
  long nAtoms = matter->numberOfAtoms();
  std::vector<long> freeList;
  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i))
      freeList.push_back(i);
  }
  REQUIRE(freeList.size() > 0);
  REQUIRE(freeList.size() < static_cast<size_t>(nAtoms));

  VectorXi moved(freeList.size());
  for (size_t i = 0; i < freeList.size(); i++) {
    moved[i] = freeList[i];
  }

  Hessian hess(params, matter.get());
  VectorXd freqs = hess.getFreqs(matter.get(), moved);

  // SVN crashes here (VectorXi size mismatch). Our fix handles it.
  long expectedSize = static_cast<long>(freeList.size()) * 3;
  REQUIRE(freqs.size() == expectedSize);
  for (long i = 0; i < freqs.size(); i++) {
    REQUIRE(std::isfinite(freqs(i)));
  }
}

} /* namespace tests */
