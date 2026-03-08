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

#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class ASEPotTest {
public:
  ASEPotTest()
      : params{},
        matter{nullptr},
        pot{nullptr},
        threshold{1e-6} {
    params.potential_options.potential = PotType::ASE_POT;
    auto script = std::filesystem::canonical("ase_lj.py").string();
    params.potential_options.extPotPath = script;

    pot = eonc::helpers::makePotential(params.potential_options.potential,
                                          params);
    matter = std::make_shared<Matter>(pot, params);

    const std::string confile("pos.con");
    const bool file_read_ok = matter->con2matter(confile);
    REQUIRE(file_read_ok);
  }

  ~ASEPotTest() = default;

protected:
  Parameters params;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<Potential> pot;
  double threshold;
};

TEST_CASE_METHOD(ASEPotTest, "ASE LJ energy and forces match reference",
                 "[PotTest][ASE]") {
  // Reference: two Al atoms at (0,0,0) and (1.5,0,0) in 10x10x10 box
  // LJ(epsilon=1, sigma=1, rc=10, smooth=False) via ASE
  const double expected_energy = -0.320339200131812;
  AtomMatrix expected_forces(2, 3);
  expected_forces.row(0) << 1.158021344587014, 0.0, 0.0;
  expected_forces.row(1) << -1.158021344587014, 0.0, 0.0;

  double calculated_energy = 0.0;
  AtomMatrix calculated_forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);

  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), calculated_forces.data(),
             &calculated_energy, nullptr, matter->getCell().data());

  REQUIRE_THAT(calculated_energy, WithinAbs(expected_energy, threshold));

  auto matEq =
      std::bind(eonc::helpers::eigenEquality<AtomMatrix>, _1, _2, threshold);
  REQUIRE(matEq(calculated_forces, expected_forces));
}

} // namespace tests
