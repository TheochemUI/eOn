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
#include "eon/Parameters.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

namespace tests {

#define EON_REQUIRE_TEST_DATA(src)                                             \
  do {                                                                         \
    if (!copyTestData(src)) {                                                  \
      SKIP("test system data not found for '"                                  \
           << (src) << "' (set EON_TEST_SYSTEMS_DIR or run via meson test)");  \
    }                                                                          \
  } while (0)

static eonc::helpers::test::QuillTestLogger _quill_setup;

static std::map<std::string, std::string>
parseResultsDat(const std::string &path) {
  std::map<std::string, std::string> result;
  std::ifstream f(path);
  std::string line;
  while (std::getline(f, line)) {
    auto pos = line.find(' ');
    if (pos != std::string::npos) {
      result[line.substr(pos + 1)] = line.substr(0, pos);
    }
  }
  return result;
}

class JobIntegrationFixture {
protected:
  std::filesystem::path workdir;
  std::filesystem::path originalDir;
  std::unique_ptr<Parameters> params;
  size_t forceCalls_{0};

  JobIntegrationFixture() : originalDir{std::filesystem::current_path()} {
    static int counter = 0;
    workdir = std::filesystem::temp_directory_path() /
              ("eon_test_jobcov_" + std::to_string(counter++));
    std::filesystem::create_directories(workdir);
  }

  ~JobIntegrationFixture() {
    std::filesystem::current_path(originalDir);
    std::filesystem::remove_all(workdir);
  }

  void writeConfig(const std::string &content) {
    std::ofstream f(workdir / "config.ini");
    f << content;
  }

  static std::filesystem::path systemsDir() {
    if (const char *e = std::getenv("EON_TEST_SYSTEMS_DIR")) {
      return std::filesystem::path(e);
    }
    return std::filesystem::current_path().parent_path();
  }

  bool copyTestData(const std::string &srcDir) {
    namespace fs = std::filesystem;
    fs::path src;
    if (srcDir.empty() || srcDir == ".") {
      src = fs::current_path();
    } else {
      std::string name = srcDir;
      if (name.rfind("../", 0) == 0) {
        name = name.substr(3);
      }
      src = systemsDir() / name;
    }
    if (!fs::is_directory(src)) {
      return false;
    }
    for (auto &entry : fs::directory_iterator(src)) {
      std::error_code ec;
      if (!entry.is_regular_file(ec) || ec) {
        continue;
      }
      fs::copy_file(entry.path(), workdir / entry.path().filename(),
                    fs::copy_options::overwrite_existing, ec);
    }
    return true;
  }

  std::map<std::string, std::string> runJob() {
    std::filesystem::current_path(workdir);
    params = std::make_unique<Parameters>();
    params->load("config.ini");
    auto job = eonc::helpers::makeJob(std::move(params));
    job->run();
    std::filesystem::current_path(originalDir);
    return parseResultsDat((workdir / "results.dat").string());
  }
};

TEST_CASE_METHOD(JobIntegrationFixture,
                 "BasinHoppingJob linear gaussian displacement writes uniques",
                 "[job][basin_hopping][displacement]") {
  EON_REQUIRE_TEST_DATA(".");
  writeConfig(R"(
[Main]
job = basin_hopping
random_seed = 7

[Potential]
potential = lj

[Basin Hopping]
steps = 1
temperature = 300.0
displacement = 0.3
displacement_algorithm = linear
displacement_distribution = gaussian
push_apart_distance = 0.4
write_unique = true

[Optimizer]
opt_method = lbfgs
converged_force = 0.05
max_iterations = 20
)");
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  auto results = runJob();
  REQUIRE(results.count("minimum_energy") > 0);
  REQUIRE(std::isfinite(std::stod(results["minimum_energy"])));
  REQUIRE(std::filesystem::exists(workdir / "min.con"));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "BasinHoppingJob quadratic displacement stays finite",
                 "[job][basin_hopping][displacement]") {
  EON_REQUIRE_TEST_DATA(".");
  writeConfig(R"(
[Main]
job = basin_hopping
random_seed = 11

[Potential]
potential = lj

[Basin Hopping]
steps = 1
temperature = 200.0
displacement = 0.25
displacement_algorithm = quadratic
displacement_distribution = uniform
push_apart_distance = 0.4

[Optimizer]
opt_method = lbfgs
converged_force = 0.05
max_iterations = 20
)");
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  auto results = runJob();
  REQUIRE(std::isfinite(std::stod(results["minimum_energy"])));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ReplicaExchangeJob linear ladder with one replica",
                 "[job][replica_exchange][integration]") {
  EON_REQUIRE_TEST_DATA(".");
  writeConfig(R"(
[Main]
job = replica_exchange
temperature = 300
random_seed = 3
parallel = false

[Potential]
potential = lj

[Dynamics]
time_step = 1.0
time = 2.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0

[Replica Exchange]
replicas = 1
temperature_distribution = linear
temperature_low = 100.0
temperature_high = 400.0
sampling_time = 2.0
exchange_period = 2.0
)");
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  auto results = runJob();
  REQUIRE(results.count("force_calls_sampling") > 0);
  REQUIRE(std::stoi(results["force_calls_sampling"]) > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ReplicaExchangeJob two-replica linear ladder",
                 "[job][replica_exchange][integration]") {
  EON_REQUIRE_TEST_DATA(".");
  writeConfig(R"(
[Main]
job = replica_exchange
temperature = 300
random_seed = 5
parallel = false

[Potential]
potential = lj

[Dynamics]
time_step = 1.0
time = 2.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0

[Replica Exchange]
replicas = 2
temperature_distribution = linear
temperature_low = 100.0
temperature_high = 500.0
sampling_time = 2.0
exchange_period = 1.0
)");
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  auto results = runJob();
  REQUIRE(std::stoi(results["force_calls_sampling"]) > 0);
}

} /* namespace tests */
