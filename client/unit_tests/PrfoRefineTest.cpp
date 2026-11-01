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

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Job.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Runtime.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace fs = std::filesystem;

namespace {

std::string slurp(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// Runs one job in a copy of a systems/ directory and keeps the outputs.
struct JobDir {
  fs::path dir;
  std::map<std::string, std::string> results;

  JobDir(const std::string &system, const std::string &config) {
    static int counter = 0;
    dir = fs::temp_directory_path() /
          ("eon_prfo_refine_" + std::to_string(counter++));
    fs::remove_all(dir);
    fs::create_directories(dir);
    const char *root = std::getenv("EON_TEST_SYSTEMS_DIR");
    const fs::path src = root ? fs::path(root) / system
                              : fs::current_path().parent_path() / system;
    REQUIRE(fs::is_directory(src));
    for (const auto &entry : fs::directory_iterator(src)) {
      fs::copy_file(entry.path(), dir / entry.path().filename());
    }
    std::ofstream(dir / "config.ini") << config;
    const fs::path back = fs::current_path();
    fs::current_path(dir);
    auto params = std::make_unique<Parameters>();
    params->load("config.ini");
    eonc::helpers::makeJob(std::move(params), eonc::Runtime{})->run();
    fs::current_path(back);
    std::istringstream in(slurp(dir / "results.dat"));
    std::string line;
    while (std::getline(in, line)) {
      const auto pos = line.find(' ');
      if (pos != std::string::npos) {
        results[line.substr(pos + 1)] = line.substr(0, pos);
      }
    }
  }

  // Max-atom force and energy of a written geometry, on the Lennard-Jones (LJ)
  // surface.
  std::pair<double, double> forceAndEnergy(const std::string &con) const {
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    auto pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    Matter m(pot, params);
    REQUIRE(eonc::io::io_ok(m.con2matter((dir / con).string())));
    return {m.maxForce(), m.getPotentialEnergy()};
  }
};

std::string dimerConfig(const std::string &extra) {
  return "[Main]\njob = saddle_search\nrandom_seed = 1\nparallel = false\n"
         "[Potential]\npotential = lj\n"
         "[Optimizer]\nopt_method = lbfgs\nconverged_force = 0.01\n"
         "max_move = 0.1\nmax_iterations = 1000\n"
         "[Saddle Search]\nmin_mode_method = dimer\n"
         "client_displace_type = load\nconverged_force = 0.01\n"
         "max_iterations = 1000\n" +
         extra +
         "[Dimer]\nopt_method = cg\nimproved = true\nconverged_angle = 0.1\n"
         "rotations_max = 20\n";
}

std::string nebConfig(const std::string &extra) {
  return "[Main]\njob = nudged_elastic_band\nparallel = false\n"
         "[Potential]\npotential = lj\n"
         "[Optimizer]\nopt_method = lbfgs\nconverged_force = 0.01\n"
         "max_move = 0.1\nmax_iterations = 1000\n"
         "[Nudged Elastic Band]\nimages = 7\nclimbing_image_method = true\n"
         "minimize_endpoints = false\nconverged_force = 0.01\n"
         "[Saddle Search]\n" +
         extra;
}

const std::string kRefine = "prfo_refine = true\nprfo_refine_force = 1e-5\n";

} // namespace

TEST_CASE("prfo_refine tightens the LJ13 dimer saddle",
          "[saddle][rgsaddle][prfo_refine]") {
  const JobDir plain("lj13_saddle", dimerConfig(""));
  const JobDir refined("lj13_saddle", dimerConfig(kRefine));
  REQUIRE(plain.results.at("termination_reason") == "0");
  REQUIRE(refined.results.at("termination_reason") == "0");
  const auto [f0, e0] = plain.forceAndEnergy("saddle.con");
  const auto [f1, e1] = refined.forceAndEnergy("saddle.con");
  CAPTURE(f0, f1, e0, e1);
  REQUIRE(f1 <= 1e-5);
  REQUIRE(f1 < f0);
  REQUIRE(e1 == Catch::Approx(e0).margin(1e-3));
  REQUIRE(std::stod(refined.results.at("final_eigenvalue")) < 0.0);
  REQUIRE(std::stol(refined.results.at("total_force_calls")) >
          std::stol(plain.results.at("total_force_calls")));
}

TEST_CASE("prfo_refine refines the LJ13 climbing image",
          "[neb][rgsaddle][prfo_refine]") {
  const JobDir plain("neb_lj13", nebConfig(""));
  const JobDir refined("neb_lj13", nebConfig(kRefine));
  REQUIRE(refined.results.at("termination_reason") == "0");
  REQUIRE(plain.results.count("force_calls_prfo_refine") == 0);
  REQUIRE(std::stod(refined.results.at("prfo_refine_converged")) == 1.0);
  const auto [fBand, eBand] = refined.forceAndEnergy("sp.con");
  const auto [f1, e1] = refined.forceAndEnergy("sp_refined.con");
  CAPTURE(fBand, f1, eBand, e1);
  REQUIRE(f1 <= 1e-5);
  REQUIRE(f1 < fBand);
  REQUIRE(e1 == Catch::Approx(eBand).margin(1e-2));
  REQUIRE(std::stod(refined.results.at("prfo_refine_max_force")) ==
          Catch::Approx(f1).margin(1e-9));
  // The band itself is unchanged.
  REQUIRE(slurp(refined.dir / "neb.con") == slurp(plain.dir / "neb.con"));
  REQUIRE(std::stod(refined.results.at("total_force_calls")) ==
          std::stod(plain.results.at("total_force_calls")) +
              std::stod(refined.results.at("force_calls_prfo_refine")));
}

TEST_CASE("prfo_refine = false leaves the default search unchanged",
          "[saddle][neb][prfo_refine]") {
  const std::string off = "prfo_refine = false\nprfo_refine_force = 1e-7\n";
  const JobDir dimerDefault("lj13_saddle", dimerConfig(""));
  const JobDir dimerOff("lj13_saddle", dimerConfig(off));
  for (const char *file : {"results.dat", "saddle.con", "mode.dat"}) {
    CAPTURE(file);
    REQUIRE(slurp(dimerOff.dir / file) == slurp(dimerDefault.dir / file));
  }
  const JobDir nebDefault("neb_lj13", nebConfig(""));
  const JobDir nebOff("neb_lj13", nebConfig(off));
  REQUIRE(nebOff.results == nebDefault.results);
  REQUIRE(slurp(nebOff.dir / "neb.con") == slurp(nebDefault.dir / "neb.con"));
  REQUIRE_FALSE(fs::exists(nebOff.dir / "sp_refined.con"));
}

} // namespace tests
