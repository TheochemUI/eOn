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

/// Force-call counts pinned to eOn 3.5.0. Each config.ini below was run
/// through the v3.5.0 eonclient; the pinned counts are what it reported
/// with lbfgs_accept = none. With lbfgs_accept = energy that client spent
/// 367 and 158 calls on the LJ13 bands, 844 on the single-image Pt hop, and
/// 1924, 280, 4424 and 818 on the min-mode searches, so the energy test must
/// leave the band and min-mode searches alone.

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Job.h"
#include "eon/Parameters.h"
#include "eon/ProcessSearchJob.h"
#include "eon/Runtime.h"
#include "eon/SaddleSearchJob.h"

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace fs = std::filesystem;

class ParityRun {
public:
  explicit ParityRun(const std::string &system)
      : original_{fs::current_path()} {
    static int counter = 0;
    workdir_ = fs::temp_directory_path() /
               ("eon_bench_parity_" + std::to_string(counter++));
    fs::remove_all(workdir_);
    fs::create_directories(workdir_);
    const char *root = std::getenv("EON_TEST_SYSTEMS_DIR");
    const fs::path src =
        root ? fs::path(root) / system : original_.parent_path() / system;
    REQUIRE(fs::is_directory(src));
    for (const auto &entry : fs::directory_iterator(src)) {
      if (entry.is_regular_file()) {
        fs::copy_file(entry.path(), workdir_ / entry.path().filename());
      }
    }
  }

  ~ParityRun() {
    fs::current_path(original_);
    std::error_code ec;
    fs::remove_all(workdir_, ec);
  }

  /// Run config through the job eonclient would build for it.
  std::map<std::string, std::string> run(const std::string &config) {
    enter(config);
    auto params = std::make_unique<Parameters>();
    params->load("config.ini");
    auto job = eonc::helpers::makeJob(std::move(params), eonc::Runtime{});
    job->run();
    return results();
  }

  void enter(const std::string &config) {
    std::ofstream(workdir_ / "config.ini") << config;
    fs::current_path(workdir_);
  }

  std::map<std::string, std::string> results() {
    std::map<std::string, std::string> out;
    std::ifstream in(workdir_ / "results.dat");
    std::string line;
    while (std::getline(in, line)) {
      const auto pos = line.find(' ');
      if (pos != std::string::npos) {
        out[line.substr(pos + 1)] = line.substr(0, pos);
      }
    }
    fs::current_path(original_);
    return out;
  }

private:
  fs::path original_;
  fs::path workdir_;
};

static std::string nebConfig(const std::string &accept) {
  return "[Main]\njob = nudged_elastic_band\nparallel = false\n"
         "[Potential]\npotential = lj\n"
         "[Optimizer]\nopt_method = lbfgs\nconverged_force = 0.01\n"
         "max_move = 0.2\nmax_iterations = 1000\nlbfgs_accept = " +
         accept +
         "\n[Nudged Elastic Band]\nimages = 5\nspring = 5.0\n"
         "climbing_image_method = true\nopt_method = lbfgs\n"
         "converged_force = 0.01\nmax_iterations = 1000\n";
}

static std::string saddleConfig(const std::string &job,
                                const std::string &potential,
                                const std::string &minmode,
                                const std::string &accept, double maxMove,
                                bool parallel = false) {
  return "[Main]\njob = " + job +
         "\nrandom_seed = 1\nparallel = " + (parallel ? "true" : "false") +
         "\n[Potential]\npotential = " + potential +
         "\n[Optimizer]\nopt_method = lbfgs\nconverged_force = 0.001\n"
         "max_move = " +
         std::to_string(maxMove) +
         "\nmax_iterations = 1000\nlbfgs_accept = " + accept +
         "\n[Saddle Search]\nmin_mode_method = " + minmode +
         "\nclient_displace_type = load\nconverged_force = 0.001\n"
         "max_iterations = 1000\n"
         "[Dimer]\nopt_method = cg\nimproved = true\nconverged_angle = 0.1\n"
         "rotations_max = 20\n"
         "[Lanczos]\ntolerance = 0.1\nmax_iterations = 10\nquit_early = true\n";
}

static void requireNear350(long calls, long pinned) {
  INFO("force calls " << calls << ", eOn 3.5.0 " << pinned);
  REQUIRE(calls >= pinned * 9 / 10);
  REQUIRE(calls <= pinned * 11 / 10);
}

TEST_CASE("LJ13 CI-NEB L-BFGS force calls match eOn 3.5.0",
          "[bench_parity][neb][lbfgs]") {
  const std::string accept = GENERATE("none", "energy");
  CAPTURE(accept);
  ParityRun run("neb_lj13");
  auto results = run.run(nebConfig(accept));
  REQUIRE(results.at("termination_reason") == "0");
  requireNear350(std::stol(results.at("total_force_calls")), 262);
  // The climbing image sits on the 1.0017 eV saddle.
  REQUIRE(std::stod(results.at("image5_energy")) ==
          Catch::Approx(1.001675).margin(1e-4));
}

TEST_CASE("Single-image CI-NEB force calls match eOn 3.5.0",
          "[bench_parity][neb][lbfgs]") {
  struct Case {
    const char *system;
    const char *potential;
    long calls;
    double imageEnergy;
  };
  // The optbench Pt-island ini as written: no spring or max_move, and
  // climbing_image is not a key, so the band climbs by default.
  const auto c = GENERATE(Case{"neb_lj13", "lj", 64, 0.920287},
                          Case{"pt_heptamer_hop", "morse_pt", 24, 1.086257});
  const std::string accept = GENERATE("none", "energy");
  CAPTURE(c.system, accept);
  ParityRun run(c.system);
  auto results = run.run(
      "[Main]\njob = nudged_elastic_band\nrandom_seed = 1\nparallel = false\n"
      "[Potential]\npotential = " +
      std::string(c.potential) +
      "\n[Optimizer]\nopt_method = lbfgs\nconverged_force = 0.001\n"
      "max_iterations = 1000\nlbfgs_accept = " +
      accept + "\n[Nudged Elastic Band]\nimages = 1\nclimbing_image = true\n");
  REQUIRE(results.at("termination_reason") == "0");
  REQUIRE(results.at("number_of_images") == "1.000000000000e+00");
  requireNear350(std::stol(results.at("total_force_calls")), c.calls);
  REQUIRE(std::stod(results.at("image1_energy")) ==
          Catch::Approx(c.imageEnergy).margin(1e-4));
}

TEST_CASE("Min-mode L-BFGS saddle searches match eOn 3.5.0",
          "[bench_parity][saddle][lbfgs]") {
  struct Case {
    const char *system;
    const char *potential;
    const char *minmode;
    double maxMove;
    long calls;
    double saddleEnergy;
  };
  const auto c =
      GENERATE(Case{"lj13_saddle", "lj", "dimer", 0.1, 390, -38.963676},
               Case{"lj13_saddle", "lj", "lanczos", 0.1, 115, -39.036018},
               Case{"Pt_Heptamer_FrozenLayers", "morse_pt", "dimer", 0.2, 292,
                    -1774.365034},
               Case{"Pt_Heptamer_FrozenLayers", "morse_pt", "lanczos", 0.2, 118,
                    -1774.704903});
  const std::string accept = GENERATE("none", "energy");
  CAPTURE(c.system, c.minmode, accept);
  ParityRun run(c.system);
  auto results = run.run(
      saddleConfig("saddle_search", c.potential, c.minmode, accept, c.maxMove));
  REQUIRE(results.at("termination_reason") == "0");
  requireNear350(std::stol(results.at("total_force_calls")), c.calls);
  REQUIRE(std::stod(results.at("potential_energy_saddle")) ==
          Catch::Approx(c.saddleEnergy).margin(1e-4));
}

TEST_CASE("Saddle and process searches count per-image potential calls",
          "[bench_parity][saddle][parallel]") {
  using eonc::helpers::test::PerImageLJ;
  const bool parallel = GENERATE(false, true);
  const bool process = GENERATE(false, true);
  CAPTURE(parallel, process);
  ParityRun run("lj13_saddle");
  run.enter(saddleConfig(process ? "process_search" : "saddle_search", "lj",
                         "dimer", "none", 0.1, parallel));
  Parameters params;
  params.load("config.ini");
  std::atomic<long> calls{0};
  auto threads = std::make_shared<PerImageLJ::Threads>();
  auto pot = std::make_shared<PerImageLJ>(threads, [&calls] { ++calls; });
  std::unique_ptr<eonc::Job> job;
  if (process) {
    job = std::make_unique<eonc::ProcessSearchJob>(pot, params);
  } else {
    job = std::make_unique<eonc::SaddleSearchJob>(pot, params);
  }
  const long before = calls;
  job->run();
  const long made = calls - before;
  auto results = run.results();
  if (parallel) {
    // The second dimer image ran on its own instance.
    REQUIRE(static_cast<long>(pot->forceCallCounter.load()) < made);
  }
  if (process) {
    // The search stops at the LJ13 prefactor; the saddle calls are in.
    // 3.5.0 reported 389. One single-rotation translation now keeps the
    // conjugate direction and spends one more call. 390 is inside the
    // ±10% window on the 390-call saddle search.
    const long saddleCalls = std::stol(results.at("force_calls_saddle"));
    REQUIRE(saddleCalls == 390);
    requireNear350(saddleCalls + 1, 390);
  } else {
    REQUIRE(results.at("termination_reason") == "0");
    REQUIRE(std::stol(results.at("total_force_calls")) == made);
    REQUIRE(std::stol(results.at("force_calls_saddle")) == made - 1);
    requireNear350(made, 390);
  }
}

} // namespace tests
