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

#include "ParametersJSON.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Parameters round-trips through JSON", "[params][json]") {
  Parameters p1;
  p1.potential_options.potential = PotType::MORSE_PT;
  p1.main_options.job = JobType::Minimization;
  p1.main_options.temperature = 500.0;
  p1.main_options.randomSeed = 12345;
  p1.optimizer_options.method = OptType::LBFGS;
  p1.optimizer_options.converged_force = 0.005;
  p1.optimizer_options.max_iterations = 999;
  p1.neb_options.image_count = 7;
  p1.neb_options.force_tolerance = 0.02;

  auto j = eonc::config::to_json(p1);

  Parameters p2;
  eonc::config::from_json(j, p2);

  REQUIRE(p2.potential_options.potential == PotType::MORSE_PT);
  REQUIRE(p2.main_options.job == JobType::Minimization);
  REQUIRE(p2.main_options.temperature == Catch::Approx(500.0));
  REQUIRE(p2.main_options.randomSeed == 12345);
  REQUIRE(p2.optimizer_options.method == OptType::LBFGS);
  REQUIRE(p2.optimizer_options.converged_force == Catch::Approx(0.005));
  REQUIRE(p2.optimizer_options.max_iterations == 999);
  REQUIRE(p2.neb_options.image_count == 7);
}

TEST_CASE("JSON to_json produces valid JSON string", "[params][json]") {
  Parameters p;
  p.potential_options.potential = PotType::LJ;
  p.main_options.job = JobType::Point;

  auto j = eonc::config::to_json(p);
  std::string s = j.dump(2);

  REQUIRE(s.find("potential") != std::string::npos);
  REQUIRE(s.find("job") != std::string::npos);
  REQUIRE(s.size() > 100);
}

TEST_CASE("JSON from_json handles missing keys gracefully", "[params][json]") {
  nlohmann::json j = nlohmann::json::object();
  j["Main"]["job"] = "point";
  j["Potential"]["potential"] = "lj";

  Parameters p;
  eonc::config::from_json(j, p);

  REQUIRE(p.main_options.job == JobType::Point);
  REQUIRE(p.potential_options.potential == PotType::LJ);
}

TEST_CASE("JSON round-trip preserves saddle search options", "[params][json]") {
  Parameters p1;
  p1.saddle_search_options.max_energy = 15.0;
  p1.saddle_search_options.max_iterations = 42;
  p1.saddle_search_options.displace_radius = 4.5;

  auto j = eonc::config::to_json(p1);
  Parameters p2;
  eonc::config::from_json(j, p2);

  REQUIRE(p2.saddle_search_options.max_energy == Catch::Approx(15.0));
  REQUIRE(p2.saddle_search_options.max_iterations == 42);
  REQUIRE(p2.saddle_search_options.displace_radius == Catch::Approx(4.5));
}

TEST_CASE("JSON to_json includes dynamics section", "[params][json]") {
  Parameters p1;
  p1.dynamics_options.time = 500.0;

  auto j = eonc::config::to_json(p1);
  std::string s = j.dump();
  REQUIRE(s.find("Dynamics") != std::string::npos);
}

TEST_CASE("JSON round-trip preserves debug compatibility flags",
          "[params][json]") {
  Parameters p1;
  p1.debug_options.write_movies = true;
  p1.debug_options.write_movies_interval = 4;
  p1.debug_options.write_deprecated_outs = true;

  auto j = eonc::config::to_json(p1);
  Parameters p2;
  eonc::config::from_json(j, p2);

  REQUIRE(j["Debug"]["write_deprecated_outs"] == true);
  REQUIRE(p2.debug_options.write_movies == true);
  REQUIRE(p2.debug_options.write_movies_interval == 4);
  REQUIRE(p2.debug_options.write_deprecated_outs == true);
}

TEST_CASE("Parameters INI load handles all sections", "[params][ini]") {
  auto tmppath = std::filesystem::temp_directory_path() / "_test_full.ini";
  std::string tmpfile = tmppath.string();
  {
    std::ofstream f(tmpfile);
    f << R"(
[Main]
job = saddle_search
temperature = 500
random_seed = 12345
finite_difference = 0.005

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.005
max_iterations = 999
max_move = 0.15
time_step = 0.05
max_time_step = 0.5

[Saddle Search]
max_energy = 15.0
max_iterations = 200
displace_radius = 4.0
displace_magnitude = 0.02
min_mode_method = dimer
converged_force = 0.01

[Dimer]
improved = true
converged_angle = 0.005
max_iterations = 100
separation = 0.005

[Dynamics]
time_step = 0.5
time = 200.0
thermostat = langevin
andersen_alpha = 0.5
andersen_collision_steps = 20

[Nudged Elastic Band]
images = 7
spring = 3.0
max_iterations = 300
converged_force = 0.005
climbing_image = true

[Parallel Replica]
dephase_time = 50.0
state_check_interval = 30.0
record_interval = 10.0
refine_transition = true
corr_time = 50.0

[Hessian]
min_displacement = 0.3
within_radius = 4.0

[Structure Comparison]
distance_difference = 0.2
energy_difference = 0.05
neighbor_cutoff = 4.0
)";
  }

  Parameters p;
  p.load(tmpfile);

  REQUIRE(p.main_options.job == JobType::Saddle_Search);
  REQUIRE(p.main_options.temperature == Catch::Approx(500.0));
  REQUIRE(p.main_options.randomSeed == 12345);
  REQUIRE(p.potential_options.potential == PotType::MORSE_PT);
  REQUIRE(p.optimizer_options.method == OptType::LBFGS);
  REQUIRE(p.optimizer_options.converged_force == Catch::Approx(0.005));
  REQUIRE(p.optimizer_options.max_iterations == 999);
  REQUIRE(p.saddle_search_options.max_energy == Catch::Approx(15.0));
  REQUIRE(p.neb_options.image_count == 7);
  REQUIRE(p.neb_options.spring.constant == Catch::Approx(3.0));
  REQUIRE(p.structure_comparison_options.distance_difference ==
          Catch::Approx(0.2));

  std::filesystem::remove(tmpfile);
}

} /* namespace tests */
