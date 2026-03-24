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

/// Micro-benchmarks for potential force evaluation.
/// Uses Catch2 BENCHMARK macro for reliable timing.

#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

/// Load a Matter from the test data directory (meson test sets CWD).
static std::shared_ptr<Matter>
loadMatter(PotType potType, const std::string &conFile = "reactant.con") {
  auto params = std::make_shared<Parameters>();
  params->potential_options.potential = potType;
  auto pot = eonc::helpers::makePotential(*params);
  auto matter = std::make_shared<Matter>(pot, *params);
  matter->con2matter(conFile);
  return matter;
}

TEST_CASE("Benchmark: Morse Pt force evaluation", "[.benchmark]") {
  auto matter = loadMatter(PotType::MORSE_PT, "reactant.con");

  BENCHMARK("Morse force (337 atoms)") {
    matter->setPositions(matter->getPositions());
    return matter->getPotentialEnergy();
  };
}

TEST_CASE("Benchmark: LJ force evaluation", "[.benchmark]") {
  auto matter = loadMatter(PotType::LJ, "reactant.con");

  BENCHMARK("LJ force (337 atoms)") {
    matter->setPositions(matter->getPositions());
    return matter->getPotentialEnergy();
  };
}

TEST_CASE("Benchmark: Matter PBC application", "[.benchmark]") {
  auto matter = loadMatter(PotType::MORSE_PT, "reactant.con");
  AtomMatrix diff = matter->getPositions();
  diff.array() += 0.1;

  BENCHMARK("PBC (337 atoms)") { return matter->pbc(diff); };
}

TEST_CASE("Benchmark: Matter kinetic energy", "[.benchmark]") {
  auto matter = loadMatter(PotType::MORSE_PT, "reactant.con");
  AtomMatrix vel = AtomMatrix::Random(matter->numberOfAtoms(), 3) * 0.01;
  matter->setVelocities(vel);

  BENCHMARK("getKineticEnergy (337 atoms)") {
    return matter->getKineticEnergy();
  };
}

} // namespace tests
