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

/// Ported from legacy unittests/utHelperFunctions.cpp to Catch2.

#include "HelperFunctions.h"
#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <cmath>
#include <cstdio>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("HelperFunctions: random() returns value in [0,1)", "[helpers]") {
  double r = eonc::helpers::random();
  REQUIRE_FALSE(std::isnan(r));
  REQUIRE(std::isfinite(r));
  REQUIRE(r >= 0.0);
  REQUIRE(r < 1.0);
}

TEST_CASE("HelperFunctions: random(seed) returns value in [0,1)", "[helpers]") {
  double r = eonc::helpers::random(42);
  REQUIRE_FALSE(std::isnan(r));
  REQUIRE(std::isfinite(r));
  REQUIRE(r >= 0.0);
  REQUIRE(r < 1.0);
}

TEST_CASE("HelperFunctions: randomDouble() returns value in [0,1)",
          "[helpers]") {
  double r = eonc::helpers::randomDouble();
  REQUIRE_FALSE(std::isnan(r));
  REQUIRE(std::isfinite(r));
  REQUIRE(r >= 0.0);
  REQUIRE(r < 1.0);
}

TEST_CASE("HelperFunctions: loadOrSynthesizeDisplacement from mode (#189/#79)",
          "[helpers][displacement]") {
  // Standalone saddle_search needs displacement without AKMC (#189).
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter initial(pot, params);
  REQUIRE(eonc::io::io_ok(initial.con2matter(std::string("reactant.con"))));
  const long nAtoms = initial.numberOfAtoms();
  REQUIRE(nAtoms > 0);

  const std::string modePath = "test_mode_for_synth.dat";
  {
    FILE *f = fopen(modePath.c_str(), "w");
    REQUIRE(f != nullptr);
    for (long i = 0; i < nAtoms; ++i) {
      if (initial.getFixed(i)) {
        fprintf(f, "0 0 0\n");
      } else {
        fprintf(f, "1 0 0\n");
      }
    }
    fclose(f);
  }
  Matter target(pot, params);
  const double scale = 0.1;
  AtomMatrix before = initial.getPositionsCopy();
  REQUIRE(eonc::helpers::loadOrSynthesizeDisplacement(
      target, initial, "missing_displacement.con", modePath, scale));
  // At least one free atom should move from the synthesized mode
  REQUIRE((target.getPositions() - before).norm() > 1e-6);
  std::remove(modePath.c_str());
}

TEST_CASE("HelperFunctions: randomDouble(max) respects upper bound",
          "[helpers]") {
  double r = eonc::helpers::randomDouble(5.0);
  REQUIRE(std::isfinite(r));
  REQUIRE(r >= 0.0);
  REQUIRE(r <= 5.0);
}

TEST_CASE("HelperFunctions: randomInt(lo, hi) respects bounds", "[helpers]") {
  for (int trial = 0; trial < 100; trial++) {
    long r = eonc::helpers::randomInt(1, 4);
    REQUIRE(r >= 1);
    REQUIRE(r <= 4);
  }
}

TEST_CASE("HelperFunctions: gaussRandom() produces finite values",
          "[helpers]") {
  double avg = 1.0, sd = 0.1;
  double r = eonc::helpers::gaussRandom(avg, sd);
  REQUIRE(std::isfinite(r));
  // Within 6 sigma (extremely unlikely to fail)
  REQUIRE(r > avg - 6.0 * sd);
  REQUIRE(r < avg + 6.0 * sd);
}

TEST_CASE("HelperFunctions: split_string_int parses CSV", "[helpers]") {
  auto result = eonc::helpers::split_string_int("1,2,3", ",");
  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == 1);
  REQUIRE(result[1] == 2);
  REQUIRE(result[2] == 3);
}

TEST_CASE("HelperFunctions: split_string_int empty string", "[helpers]") {
  auto result = eonc::helpers::split_string_int("", ",");
  REQUIRE(result.empty());
}

TEST_CASE("HelperFunctions: maxAtomMotionV", "[helpers]") {
  Eigen::VectorXd v(6);
  v << 1.0, 0.0, 0.0, 0.0, 3.0, 4.0;
  double maxMotion = eonc::helpers::maxAtomMotionV(v);
  REQUIRE(maxMotion == Catch::Approx(5.0));
}

} // namespace tests
