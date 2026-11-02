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
#include "eon/MinModeSaddleSearch.h"
#include "eon/Potential.h"
#include "eon/XtsciIrc.h"

#include <memory>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("rgsaddle IRC from the LJ13 saddle relaxes into the NEB minima",
          "[saddle_search][rgsaddle][irc]") {
  if (!eonc::ircAvailable()) {
    SKIP("rgsaddle older than C ABI minor 14 has no IRC");
  }
  Parameters params;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-4;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 2000;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-4;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));

  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  auto saddle = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(reactant->con2matter(std::string("reactant.con"))));
  REQUIRE(eonc::io::io_ok(product->con2matter(std::string("product.con"))));
  REQUIRE(eonc::io::io_ok(saddle->con2matter(std::string("saddle.con"))));

  // A nudged elastic band gives a loose saddle at its climbing image. Like
  // the saddle search job, the test tightens it with a min-mode search,
  // which also supplies the mode.
  const AtomMatrix guess =
      saddle->pbc(product->getPositions() - reactant->getPositions());
  eonc::MinModeSaddleSearch search(saddle, guess.normalized(),
                                   reactant->getPotentialEnergy(), params, pot);
  REQUIRE(search.run() == eonc::MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(search.getEigenvalue() < 0.0);

  const auto path = eonc::followIrc(*saddle, search.getEigenvector(), params);
  REQUIRE_FALSE(path.minus.positions.empty());
  REQUIRE_FALSE(path.plus.positions.empty());
  const double saddleEnergy = saddle->getPotentialEnergy();
  for (const auto *branch : {&path.minus, &path.plus}) {
    for (size_t i = 0; i < branch->energies.size(); ++i) {
      CAPTURE(i, branch->energies[i], saddleEnergy);
      REQUIRE(branch->energies[i] < saddleEnergy);
      if (i > 0) {
        REQUIRE(branch->energies[i] <= branch->energies[i - 1] + 1e-8);
      }
    }
  }

  REQUIRE(reactant->relax(true));
  REQUIRE(product->relax(true));
  Matter ends[2] = {Matter(*saddle), Matter(*saddle)};
  ends[0].setPositions(path.minus.positions.back());
  ends[1].setPositions(path.plus.positions.back());
  REQUIRE(ends[0].relax(true));
  REQUIRE(ends[1].relax(true));
  // A single session for both branches started the second one at the
  // 1e-4 step floor, and it ran out of steps far from the minimum.
  for (int k = 0; k < 2; ++k) {
    const eonc::IrcBranch &branch = k == 0 ? path.minus : path.plus;
    CAPTURE(k, branch.steps, branch.reachedMinimum, branch.energies.back(),
            ends[k].getPotentialEnergy());
    REQUIRE(branch.reachedMinimum);
    REQUIRE(branch.energies.back() - ends[k].getPotentialEnergy() < 1e-3);
  }

  const bool minusIsReactant = ends[0].compare(*reactant);
  const bool minusIsProduct = ends[0].compare(*product);
  CAPTURE(ends[0].getPotentialEnergy(), ends[1].getPotentialEnergy(),
          reactant->getPotentialEnergy(), product->getPotentialEnergy());
  REQUIRE(minusIsReactant != minusIsProduct);
  if (minusIsReactant) {
    REQUIRE(ends[1].compare(*product));
  } else {
    REQUIRE(ends[1].compare(*reactant));
  }
}

} // namespace tests
