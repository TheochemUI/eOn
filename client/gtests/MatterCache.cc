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
#include "client/matter/MatterCreator.hpp"
#include <chrono>
#include <memory>

using namespace Catch::Matchers;
using namespace eonc;

using namespace std::chrono;

TEST_CASE("Matter caching", "[Matter]") {
  const auto config =
      toml::table{{"Potential", toml::table{{"potential", "lj"}}}};
  auto pot_default = eonc::makePotential(config);
  auto CACHELOT_EONCTEST = cachelot::cache::Cache::Create(
      eonc::cache_memory, eonc::page_size, eonc::hash_initial, true);
  auto matter = eonc::Matter(pot_default, &CACHELOT_EONCTEST);
  std::string confile("pos.con");
  eonc::mat::ConFileParser cfp;
  cfp.parse(matter, confile);

  // Initial cache miss
  auto start = high_resolution_clock::now();
  auto energy1 = matter.getPotentialEnergy();
  auto frcs = matter.getForces();
  auto end = high_resolution_clock::now();
  auto base_call = duration_cast<nanoseconds>(end - start).count();

  SECTION("Cache hit") {
    start = high_resolution_clock::now();
    double energy2 = matter.getPotentialEnergy();
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<nanoseconds>(end2 - start).count();

    REQUIRE(energy1 == energy2);
    // Cache hit should be faster
    REQUIRE(duration2 < base_call);
  }

  SECTION("Cache invalidation on position change") {
    matter.setPositions(matter.positions.array() * 2);

    start = high_resolution_clock::now();
    matter.getPotentialEnergy();
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<nanoseconds>(end2 - start).count();
    // Cache miss should take longer
    // TODO(rg) :: Needs to be tested correctly..
    // REQUIRE(duration2 > base_call);
  }

  SECTION("Cache hit on older elements") {
    cfp.parse(matter, confile);

    auto start = high_resolution_clock::now();
    matter.getPotentialEnergy();
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<nanoseconds>(end2 - start).count();

    // Cache hit on older keys
    REQUIRE(duration2 < base_call);

    matter.setPositions(matter.positions.array() * 2);
    start = high_resolution_clock::now();
    matter.getPotentialEnergy();
    end2 = high_resolution_clock::now();
    duration2 = duration_cast<nanoseconds>(end2 - start).count();
    REQUIRE(duration2 < base_call);
  }
}
