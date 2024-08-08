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
#include "client/Eigen.h"
#include "client/Potential.h"
#include "client/potentials/LJ/LJ.h"
#include "client/potentials/Morse/Morse.h"

using namespace Catch::Matchers;
using namespace eonc;

TEST_CASE("Cache key uniqueness and consistency", "[cache]") {
  const auto config_lj =
      toml::table{{"Potential", toml::table{{"potential", "lj"}}}};
  const auto config_morse_pt =
      toml::table{{"Potential", toml::table{{"potential", "morse_pt"}}}};

  // Create potential instances
  auto pot_lj =
      dynamic_cast<eonc::Potential<LJ> *>(makePotential(config_lj).get());
  auto pot_morse_pt = dynamic_cast<eonc::Potential<Morse> *>(
      makePotential(config_morse_pt).get());

  Eigen::MatrixXd positions(5, 3);
  positions << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  Vector<size_t> atomic_numbers(5);
  atomic_numbers << 1, 2, 3, 4, 5;

  size_t hash_lj = pot_lj->getHash(positions, atomic_numbers);
  size_t hash_morse_pt = pot_morse_pt->getHash(positions, atomic_numbers);

  SECTION("Unique keys for different potential types") {
    REQUIRE(hash_lj != hash_morse_pt);
  }

  // Modify positions slightly
  Eigen::MatrixXd positions_modified = positions;
  positions_modified(0, 0) += 1;

  // Compute hash for modified positions for LJ potential
  size_t hash_lj_modified = pot_lj->getHash(positions_modified, atomic_numbers);

  SECTION("Consistent keys for same potential type with same inputs") {
    // Recompute hash for original positions
    size_t hash_lj_recomputed = pot_lj->getHash(positions, atomic_numbers);
    REQUIRE(hash_lj == hash_lj_recomputed);
  }

  SECTION("Different keys for modified inputs") {
    REQUIRE(hash_lj != hash_lj_modified);
  }

  // Create another instance of the same potential type
  auto pot_lj_2 =
      dynamic_cast<eonc::Potential<LJ> *>(makePotential(config_lj).get());

  SECTION("Consistent keys across different instances of same potential type") {
    size_t hash_lj_2 = pot_lj_2->getHash(positions, atomic_numbers);
    REQUIRE(hash_lj == hash_lj_2);
  }
  SECTION("Determinism when returning to older positions") {
    // Compute hash for modified positions for LJ potential
    size_t hash_lj_modified =
        pot_lj->getHash(positions_modified, atomic_numbers);

    // Return to original positions and compute hash again
    size_t hash_lj_recomputed = pot_lj->getHash(positions, atomic_numbers);
    REQUIRE(hash_lj == hash_lj_recomputed);

    // Compute hash for modified positions again
    size_t hash_lj_modified_again =
        pot_lj->getHash(positions_modified, atomic_numbers);
    REQUIRE(hash_lj_modified == hash_lj_modified_again);
  }
}
