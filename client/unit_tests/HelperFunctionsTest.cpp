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

#include "eon/HelperFunctions.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Eigen.h"
#include "eon/EpiCenters.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/RandomNumbers.h"

#include <atomic>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("getRelevantFile accepts a name with no extension", "[helpers]") {
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() / "eon_relevant_file";
  fs::create_directories(dir);
  const auto old = fs::current_path();
  fs::current_path(dir);
  REQUIRE(eonc::helpers::getRelevantFile("config") == "config");
  {
    std::ofstream{dir / "config_cp"};
  }
  REQUIRE(eonc::helpers::getRelevantFile("config") == "config_cp");
  fs::current_path(old);
  fs::remove_all(dir);
}

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

  const auto tmp = std::filesystem::temp_directory_path() /
                   "eon_mode_for_synth.dat";
  {
    FILE *f = fopen(tmp.c_str(), "w");
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
      target, initial, "missing_displacement.con", tmp.string(), scale));
  // At least one free atom should move from the synthesized mode
  REQUIRE((target.getPositions() - before).norm() > 1e-6);
  std::filesystem::remove(tmp);
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

TEST_CASE("ran2 streams are independent across threads", "[helpers][rng]") {
  std::atomic<int> bad{0};
  auto worker = [&](long seed) {
    eonc::rng::random(seed);
    for (int i = 0; i < 2000; ++i) {
      const double r = eonc::rng::random();
      if (!(r > 0.0 && r < 1.0)) {
        ++bad;
      }
    }
  };
  std::thread a(worker, 11);
  std::thread b(worker, 17);
  a.join();
  b.join();
  REQUIRE(bad.load() == 0);
}

TEST_CASE("SaddleSearchJob listed_atoms moves a free atom",
          "[helpers][job][listed_atoms]") {
  // SaddleSearchJob / ProcessSearchJob call applyClientDisplacement
  // when client_displace_type = listed_atoms.
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  params.saddle_search_options.displace_type =
      std::string(eonc::EpiCenters::DISP_LISTED_ATOMS);
  params.saddle_search_options.displace_atom_list = {0};
  params.saddle_search_options.displace_radius = 0.0;
  params.saddle_search_options.displace_magnitude = 0.2;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter initial(pot, params);
  initial.resize(2);
  AtomMatrix pos(2, 3);
  pos << 0.0, 0.0, 0.0, 5.0, 0.0, 0.0;
  initial.setPositions(pos);
  initial.setFixed(1, 1);
  Matrix3d cell = Matrix3d::Identity() * 20.0;
  initial.setCell(cell);

  Matter target(pot, params);
  AtomMatrix mode;
  REQUIRE(
      eonc::helpers::applyClientDisplacement(target, initial, params, &mode));
  const AtomMatrix delta = target.getPositions() - initial.getPositions();
  REQUIRE(delta.row(0).norm() > 0.0);
  REQUIRE(delta.row(1).norm() == Catch::Approx(0.0));
  REQUIRE_FALSE(initial.getFixed(0));
}

TEST_CASE("applyClientDisplacement load type is a no-op",
          "[helpers][listed_atoms]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  params.saddle_search_options.displace_type =
      std::string(eonc::EpiCenters::DISP_LOAD);
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter initial(pot, params);
  initial.resize(1);
  AtomMatrix pos(1, 3);
  pos << 0.0, 0.0, 0.0;
  initial.setPositions(pos);
  Matter target(pot, params);
  REQUIRE_FALSE(
      eonc::helpers::applyClientDisplacement(target, initial, params, nullptr));
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

TEST_CASE("HelperFunctions: convergenceMetricLabel known spellings",
          "[helpers][convergence]") {
  REQUIRE(eonc::helpers::convergenceMetricLabel("norm") == "||Force||");
  REQUIRE(eonc::helpers::convergenceMetricLabel("max_atom") ==
          "Max atom force");
  REQUIRE(eonc::helpers::convergenceMetricLabel("max_component") ==
          "Max force comp");
  REQUIRE_FALSE(eonc::helpers::convergenceMetricLabel("typo").has_value());
}

TEST_CASE("HelperFunctions: requireKnownConvergenceMetric throws on typo",
          "[helpers][convergence]") {
  REQUIRE_NOTHROW(
      eonc::helpers::requireKnownConvergenceMetric("norm", "[test]"));
  REQUIRE_THROWS_AS(
      eonc::helpers::requireKnownConvergenceMetric("nope", "[test]"),
      std::invalid_argument);
  try {
    eonc::helpers::requireKnownConvergenceMetric("nope", "[test]");
    FAIL("expected throw");
  } catch (const std::invalid_argument &e) {
    REQUIRE_THAT(std::string(e.what()), Catch::Matchers::ContainsSubstring(
                                            "unknown convergence_metric"));
    REQUIRE_THAT(std::string(e.what()),
                 Catch::Matchers::ContainsSubstring("nope"));
  }
}

} // namespace tests
