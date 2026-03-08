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
#include "catch2/catch_amalgamated.hpp"
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

class ExtPotTest {
public:
  ExtPotTest()
      : params{},
        matter{nullptr},
        pot_ext{nullptr},
        threshold{1e-6} {
    params.potential_options.potential = PotType::EXT_POT;
    // Use absolute path so the test works regardless of cwd
    auto ext_pot_script = std::filesystem::canonical("ext_pot").string();
    params.potential_options.extPotPath = ext_pot_script;

    pot_ext = eonc::helpers::makePotential(
        params.potential_options.potential, params);
    matter = std::make_shared<Matter>(pot_ext, params);

    const std::string confile("pos.con");
    const bool file_read_ok = matter->con2matter(confile);
    REQUIRE(file_read_ok);
  }

  ~ExtPotTest() = default;

protected:
  Parameters params;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<Potential> pot_ext;
  double threshold;
};

TEST_CASE_METHOD(ExtPotTest, "ExtPot harmonic spring energy and forces",
                 "[PotTest][ExtPot]") {
  // Two Al atoms at (0,0,0) and (1.5,0,0)
  // Harmonic spring: k=1.0, r0=1.0
  // r = 1.5, E = 0.5 * 1.0 * (1.5 - 1.0)^2 = 0.125
  // F on atom 0 = k*(r-r0)/r * (dx,dy,dz) = (0.5/1.5)*1.5 in x = +0.5
  //   actually: fmag = k*(r-r0)/r = 1.0*0.5/1.5 = 1/3
  //   F0 = fmag * (1.5, 0, 0) = (0.5, 0, 0)
  //   F1 = -F0 = (-0.5, 0, 0)
  const double expected_energy = 0.125;
  AtomMatrix expected_forces(2, 3);
  expected_forces.row(0) << 0.5, 0.0, 0.0;
  expected_forces.row(1) << -0.5, 0.0, 0.0;

  double calculated_energy = 0.0;
  AtomMatrix calculated_forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);

  pot_ext->force(matter->numberOfAtoms(), matter->getPositions().data(),
                 matter->getAtomicNrs().data(), calculated_forces.data(),
                 &calculated_energy, nullptr, matter->getCell().data());

  REQUIRE_THAT(calculated_energy, WithinAbs(expected_energy, threshold));

  auto matEq =
      std::bind(eonc::helpers::eigenEquality<AtomMatrix>, _1, _2, threshold);
  REQUIRE(matEq(calculated_forces, expected_forces));
}

} // namespace tests
