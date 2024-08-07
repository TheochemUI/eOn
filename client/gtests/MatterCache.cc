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
#include "client/matter/Matter.h"
#include <memory>

using namespace Catch::Matchers;
using namespace eonc;

using namespace std::chrono;

// Mock potential class for testing
class MockPotential : public Potential<MockPotential> {
public:
  MockPotential() {}

  void forceImpl(const ForceInput &params, ForceOut *efvd) override {
    // std::this_thread::sleep_for(
    //     std::chrono::milliseconds(50)); // Simulate heavy computation
    efvd->energy = params.pos[0];
  }
};

TEST_CASE("Matter caching", "[Matter]") {
  // auto potential = std::make_shared<MockPotential>();
  // Matter matter(potential);
  const auto config =
      toml::table{{"Potential", toml::table{{"potential", "lj"}}}};
  auto pot_default = eonc::makePotential(config);
  auto matter = eonc::Matter(pot_default);

  SECTION("Initial cache miss and subsequent cache hit") {
    matter.resize(10); // Resize to 10 atoms for testing
    matter.positions = AtomMatrix::Random(10, 3);

    auto start = high_resolution_clock::now();
    double energy1 = matter.getPotentialEnergy();
    auto end = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(end - start).count();

    start = high_resolution_clock::now();
    double energy2 = matter.getPotentialEnergy();
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(end2 - start).count();

    REQUIRE(energy1 == energy2);
    // Cache hit should be faster
    REQUIRE(duration2 < duration1);
  }

  SECTION("Cache invalidation on position change") {
    matter.resize(10);
    matter.positions = AtomMatrix::Random(10, 3);
    AtomMatrix positions2 = AtomMatrix::Random(10, 3);

    auto start = high_resolution_clock::now();
    double energy1 = matter.getPotentialEnergy();
    auto end = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(end - start).count();

    matter.setPositions(positions2);

    start = high_resolution_clock::now();
    double energy2 = matter.getPotentialEnergy();
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(end2 - start).count();
    // Energy should be different for different positions
    REQUIRE(energy1 == energy2);
    // Cache miss should take longer
    REQUIRE(duration2 > duration1);
  }
}
