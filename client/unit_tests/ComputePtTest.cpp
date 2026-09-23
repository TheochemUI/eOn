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

#include "catch2/catch_amalgamated.hpp"
#include "eon/potentials/Water/potential_base.hpp"

#include <array>

namespace tests {

namespace {

class PtLjProbe : public forcefields::PotentialBase {
public:
  void eval(int nAtoms, double positions[], double forces[], double &energy,
            double const periods[], bool const fixed[]) {
    computePt(nAtoms, positions, forces, energy, periods, fixed);
  }
};

struct PairResult {
  double energy{};
  std::array<double, 6> forces{};
};

PairResult evalPair(bool fixed0, bool fixed1) {
  PtLjProbe pot;
  double positions[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  double periods[3] = {40.0, 40.0, 40.0};
  bool fixed[2] = {fixed0, fixed1};
  PairResult out;
  pot.eval(2, positions, out.forces.data(), out.energy, periods, fixed);
  return out;
}

} // namespace

TEST_CASE("computePt keeps an LJ pair when one atom is fixed", "[water][pt]") {
  const PairResult free = evalPair(false, false);
  REQUIRE(free.energy != 0.0);
  REQUIRE(free.forces[0] != 0.0);
  REQUIRE(free.forces[0] == -free.forces[3]);
  REQUIRE(free.forces[1] == 0.0);
  REQUIRE(free.forces[2] == 0.0);

  const PairResult fix0 = evalPair(true, false);
  const PairResult fix1 = evalPair(false, true);
  REQUIRE(fix0.energy == free.energy);
  REQUIRE(fix1.energy == free.energy);
  REQUIRE(fix0.forces == free.forces);
  REQUIRE(fix1.forces == free.forces);

  const PairResult both = evalPair(true, true);
  REQUIRE(both.energy == 0.0);
  for (double component : both.forces) {
    REQUIRE(component == 0.0);
  }
}

} // namespace tests
