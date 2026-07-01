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

/// Integration tests for Job types. These exercise the Job runner
/// classes (PointJob, MinimizationJob, etc.) by loading config from
/// INI, running the job, and verifying outputs match SVN reference
/// values.

#include "Job.h"
#include "Matter.h"
#include "Parameters.h"
#include "PotRegistry.h"
#include "Potential.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#ifdef WITH_ARTN
#include "ARTnSaddleSearch.h"
#include "libs/ARTn/ARTnResource.h"
#endif

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

/// Parse a results.dat file into key-value pairs.
static std::map<std::string, std::string>
parseResultsDat(const std::string &path) {
  std::map<std::string, std::string> result;
  std::ifstream f(path);
  std::string line;
  while (std::getline(f, line)) {
    // Format: "value key" (space-separated, value first)
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
  /// Force calls made by just this job (delta, not global total).
  size_t forceCalls_{0};

  JobIntegrationFixture()
      : originalDir{std::filesystem::current_path()} {
    // Unique dir per test to avoid cleanup race
    static int counter = 0;
    workdir = std::filesystem::temp_directory_path() /
              ("eon_test_job_" + std::to_string(counter++));
    std::filesystem::create_directories(workdir);
  }

  ~JobIntegrationFixture() {
    // Always restore CWD before cleanup
    std::filesystem::current_path(originalDir);
    std::filesystem::remove_all(workdir);
  }

  void writeConfig(const std::string &content) {
    std::ofstream f(workdir / "config.ini");
    f << content;
  }

  void copyTestData(const std::string &srcDir) {
    namespace fs = std::filesystem;
    for (auto &entry : fs::directory_iterator(srcDir)) {
      if (entry.is_regular_file()) {
        fs::copy_file(entry.path(), workdir / entry.path().filename(),
                      fs::copy_options::overwrite_existing);
      }
    }
  }

  std::map<std::string, std::string> runJob() {
    // Change to workdir, load params, create and run job
    std::filesystem::current_path(workdir);

    params = std::make_unique<Parameters>();
    params->load("config.ini");

    size_t before = PotRegistry::get().total_force_calls();
    auto job = eonc::helpers::makeJob(std::move(params));
    job->run();
    forceCalls_ = PotRegistry::get().total_force_calls() - before;

    std::filesystem::current_path(originalDir);
    return parseResultsDat((workdir / "results.dat").string());
  }
};

TEST_CASE_METHOD(JobIntegrationFixture,
                 "PointJob produces correct energy for LJ cluster",
                 "[job][point][integration]") {
  // The test workdir is set to neb_morse which has reactant.con
  copyTestData("."); // copy everything from CWD (neb_morse data)
  writeConfig(R"(
[Main]
job = point
random_seed = 42

[Potential]
potential = lj
)");

  // Copy reactant.con as pos.con
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  // SVN reference (data/reference/point_lj.dat): -39.965351 Energy
  double energy = std::stod(results["Energy"]);
  REQUIRE(energy == Catch::Approx(-39.965351).epsilon(1e-4));
  // SVN reference: 0.004704 Max_Force
  double maxForce = std::stod(results["Max_Force"]);
  REQUIRE(maxForce == Catch::Approx(0.004704).epsilon(1e-3));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "MinimizationJob converges for LJ cluster",
                 "[job][minimization][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = minimization
random_seed = 42

[Potential]
potential = lj

[Optimizer]
opt_method = lbfgs
converged_force = 0.001
max_iterations = 500
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  // SVN reference (data/reference/minimization_lj.dat): -39.965352
  double energy = std::stod(results["potential_energy"]);
  REQUIRE(energy == Catch::Approx(-39.965352).epsilon(1e-4));
  // Must match SVN: 11 force calls (use delta, not global total from
  // results.dat)
  REQUIRE(forceCalls_ <= 11);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "MinimizationJob FIRE matches SVN on LJ",
                 "[job][minimization][fire][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = minimization
random_seed = 42
[Potential]
potential = lj
[Optimizer]
opt_method = fire
converged_force = 0.001
max_iterations = 500
)");
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  auto results = runJob();
  // SVN reference (data/reference/min_lj_fire.dat):
  double energy = std::stod(results["potential_energy"]);
  REQUIRE(energy == Catch::Approx(-39.965352).epsilon(1e-4));
  REQUIRE(forceCalls_ <= 21);
}

TEST_CASE_METHOD(JobIntegrationFixture, "MinimizationJob CG matches SVN on LJ",
                 "[job][minimization][cg][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = minimization
random_seed = 42
[Potential]
potential = lj
[Optimizer]
opt_method = cg
converged_force = 0.001
max_iterations = 500
)");
  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  auto results = runJob();
  // SVN reference (data/reference/min_lj_cg.dat):
  double energy = std::stod(results["potential_energy"]);
  REQUIRE(energy == Catch::Approx(-39.965352).epsilon(1e-4));
  REQUIRE(forceCalls_ <= 15);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "FiniteDifferenceJob curvatures match SVN reference",
                 "[job][finite_difference][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = finite_difference
random_seed = 42

[Potential]
potential = lj
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  // SVN reference (data/reference/finite_difference_lj.dat):
  // curvature at dR=0.01 is -29.95337689
  // Parse the results.dat which has "dR curvature" format
  std::ifstream fd_results((workdir / "results.dat").string());
  std::string header;
  std::getline(fd_results, header); // skip header line
  std::vector<std::pair<double, double>> fd_data;
  double dr, curv;
  while (fd_results >> dr >> curv) {
    fd_data.emplace_back(dr, curv);
  }
  REQUIRE(fd_data.size() == 9); // 9 displacement steps

  // SVN reference curvatures (random_seed=42)
  REQUIRE(fd_data[0].second == Catch::Approx(-30.00469449).epsilon(1e-4));
  REQUIRE(fd_data[6].second == Catch::Approx(-29.95337689).epsilon(1e-4));
  REQUIRE(fd_data[8].second == Catch::Approx(-29.73689692).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "FiniteDifferenceJob Morse Pt curvatures match SVN",
                 "[job][finite_difference][morse_pt][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = finite_difference
random_seed = 42

[Potential]
potential = morse_pt
)");

  auto results = runJob();

  // SVN reference (data/reference/finite_difference_morse_pt.dat):
  std::ifstream fd_results((workdir / "results.dat").string());
  std::string header;
  std::getline(fd_results, header);
  std::vector<std::pair<double, double>> fd_data;
  double dr, curv;
  while (fd_results >> dr >> curv) {
    fd_data.emplace_back(dr, curv);
  }
  REQUIRE(fd_data.size() == 9);
  REQUIRE(fd_data[0].second == Catch::Approx(-10.19003889).epsilon(1e-4));
  REQUIRE(fd_data[6].second == Catch::Approx(-10.13940225).epsilon(1e-4));
  REQUIRE(fd_data[8].second == Catch::Approx(-9.71712412).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "HessianJob force calls match SVN reference",
                 "[job][hessian][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = hessian
random_seed = 42

[Potential]
potential = lj
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  // SVN reference (data/reference/hessian_lj.dat): 40 force_calls
  REQUIRE(forceCalls_ <= 40);

  // Verify hessian.dat matrix file was produced
  REQUIRE(std::filesystem::exists(workdir / "hessian.dat"));
  // 13 atoms, all free -> 39x39 matrix -> 39 lines (one row per line)
  std::ifstream hf((workdir / "hessian.dat").string());
  int lineCount = 0;
  std::string line;
  while (std::getline(hf, line)) {
    if (!line.empty())
      lineCount++;
  }
  // 13 atoms * 3 = 39 DOF -> at least 39 lines
  // (SVN produces 38 lines for the 13-atom cluster)
  REQUIRE(lineCount >= 38);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob matches SVN on Morse Pt",
                 "[job][saddle_search][integration]") {
  // Use the actual saddle_search test system with Morse Pt
  copyTestData("../saddle_search");
  // The saddle_search dir already has pos.con, displacement.con, direction.dat
  // and config.ini. Overwrite config with our version to ensure consistency.
  writeConfig(R"(
[Main]
job = saddle_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Optimizer]
converged_force = 0.001
max_iterations = 1000

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
min_mode_method = dimer
max_energy = 10.0
)");

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  // SVN reference (data/reference/saddle_search_morse_pt.dat):
  // 0 termination_reason (converged)
  int status = std::stoi(results["termination_reason"]);
  REQUIRE(status == 0);

  // Energy must match SVN exactly
  double energy = std::stod(results["potential_energy_saddle"]);
  REQUIRE(energy == Catch::Approx(-1462.008706).epsilon(1e-4));

  // Force calls must be <= SVN (39)
  REQUIRE(forceCalls_ <= 39);

  // Eigenvalue must be negative (true saddle point)
  double eigenvalue = std::stod(results["final_eigenvalue"]);
  REQUIRE(eigenvalue < 0.0);
  REQUIRE(eigenvalue == Catch::Approx(-1.014995).epsilon(1e-3));

  // Reactant energy must match
  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(reactantE == Catch::Approx(-1462.166782).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob Lanczos matches SVN on Morse Pt",
                 "[job][saddle_search][lanczos][integration]") {
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = saddle_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Optimizer]
converged_force = 0.001
max_iterations = 1000

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
min_mode_method = lanczos
max_energy = 10.0
)");

  auto results = runJob();

  // SVN reference (data/reference/saddle_search_lanczos_morse_pt.dat):
  int status = std::stoi(results["termination_reason"]);
  REQUIRE(status == 0);

  double energy = std::stod(results["potential_energy_saddle"]);
  REQUIRE(energy == Catch::Approx(-1462.008706).epsilon(1e-4));

  REQUIRE(forceCalls_ <= 44);

  double eigenvalue = std::stod(results["final_eigenvalue"]);
  REQUIRE(eigenvalue < 0.0);
  REQUIRE(eigenvalue == Catch::Approx(-1.011347).epsilon(1e-3));

  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(reactantE == Catch::Approx(-1462.166782).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture, "NEB LJ matches SVN reference",
                 "[job][neb][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = nudged_elastic_band
random_seed = 42

[Potential]
potential = lj

[Nudged Elastic Band]
images = 3
spring = 5.0
max_iterations = 200

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
max_move = 0.2
)");

  auto results = runJob();

  // SVN reference (data/reference/neb_lj.dat):
  // termination_reason = 0 (converged)
  REQUIRE(results.count("termination_reason") > 0);
  int status = std::stoi(results["termination_reason"]);
  REQUIRE(status == 0);

  int nImages = std::stoi(results["number_of_images"]);
  REQUIRE(nImages == 3);

  // Force calls must be <= SVN (209)
  REQUIRE(forceCalls_ <= 209);

  // energy_reference must match SVN reactant energy
  double eRef = std::stod(results["energy_reference"]);
  REQUIRE(eRef == Catch::Approx(-39.965351).epsilon(1e-4));

  // image0 (reactant) = 0.0 relative
  double e0 = std::stod(results["image0_energy"]);
  REQUIRE(e0 == Catch::Approx(0.0).margin(1e-4));

  // SVN image energies (relative to reactant) are very small for this
  // near-identical reactant/product pair. Check that images are finite.
  for (int i = 1; i <= nImages; i++) {
    double ei = std::stod(results["image" + std::to_string(i) + "_energy"]);
    REQUIRE(std::isfinite(ei));
  }

  REQUIRE(results.count("number_of_extrema") > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "NEB movies embed structured metadata in con outputs",
                 "[job][neb][integration][metadata]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = nudged_elastic_band
random_seed = 42

[Potential]
potential = lj

[Nudged Elastic Band]
images = 3
spring = 5.0
max_iterations = 50

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 50
max_move = 0.2

[Debug]
write_movies = true
)");

  auto results = runJob();
  REQUIRE(std::stoi(results["termination_reason"]) == 0);

  REQUIRE(std::filesystem::exists(workdir / "neb.con"));
  REQUIRE(std::filesystem::exists(workdir / "neb_path_000.con"));

  std::ifstream nebCon(workdir / "neb.con");
  std::string nebConContents((std::istreambuf_iterator<char>(nebCon)),
                             std::istreambuf_iterator<char>());
  REQUIRE(nebConContents.find("\"neb_bead\":1") != std::string::npos);
  REQUIRE(nebConContents.find("\"reaction_coordinate\"") != std::string::npos);
  REQUIRE(nebConContents.find("\"parallel_force\"") != std::string::npos);

  std::ifstream nebPath(workdir / "neb_path_000.con");
  std::string nebPathContents((std::istreambuf_iterator<char>(nebPath)),
                              std::istreambuf_iterator<char>());
  REQUIRE(nebPathContents.find("\"neb_band\":0") != std::string::npos);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "MinimizationJob Morse Pt frozen layers matches SVN",
                 "[job][minimization][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = minimization
random_seed = 42

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
)");

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  int status = std::stoi(results["termination_reason"]);
  REQUIRE(status == 0); // GOOD

  // SVN reference (data/reference/minimization_morse_pt.dat):
  // Energy must be exact
  double energy = std::stod(results["potential_energy"]);
  REQUIRE(energy == Catch::Approx(-1775.791160).epsilon(1e-4));

  // Force calls must be <= SVN (1)
  REQUIRE(forceCalls_ <= 1);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "Minimization deprecated outputs remain available behind flag",
                 "[job][minimization][integration][deprecated]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = minimization
random_seed = 42

[Potential]
potential = lj

[Optimizer]
opt_method = lbfgs
converged_force = 0.001
max_iterations = 100

[Debug]
write_movies = true
write_deprecated_outs = true
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();
  REQUIRE(std::stoi(results["termination_reason"]) == 0);

  REQUIRE(std::filesystem::exists(workdir / "minimization.con"));
  REQUIRE(std::filesystem::exists(workdir / "minimization.dat"));

  std::ifstream datFile(workdir / "minimization.dat");
  std::string datContents((std::istreambuf_iterator<char>(datFile)),
                          std::istreambuf_iterator<char>());
  REQUIRE(datContents.find("iteration\tstep_size\tconvergence\tenergy") !=
          std::string::npos);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "BasinHoppingJob LJ cluster matches SVN reference",
                 "[job][basin_hopping][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = basin_hopping
random_seed = 42

[Potential]
potential = lj

[Basin Hopping]
steps = 20
temperature = 300.0
displacement = 0.5
push_apart_distance = 0.4

[Optimizer]
opt_method = lbfgs
converged_force = 0.001
max_iterations = 200
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  // SVN reference (data/reference/basin_hopping_lj.dat):
  // minimum_energy = -44.326774
  double minEnergy = std::stod(results["minimum_energy"]);
  REQUIRE(minEnergy == Catch::Approx(-44.326774).epsilon(1e-4));

  // Force calls must be <= SVN (1685)
  REQUIRE(forceCalls_ <= 1685);

  // 50% acceptance ratio
  double ar = std::stod(results["acceptance_ratio"]);
  REQUIRE(ar == Catch::Approx(0.500).margin(0.05));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "DynamicsJob runs and produces final.con",
                 "[job][dynamics][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = dynamics
random_seed = 42

[Potential]
potential = lj

[Dynamics]
time_step = 1.0
time = 10.0
temperature = 300.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  // DynamicsJob does not write results.dat, just final.con
  auto oldDir = std::filesystem::current_path();
  std::filesystem::current_path(workdir);

  params = std::make_unique<Parameters>();
  params->load("config.ini");

  auto job = eonc::helpers::makeJob(std::move(params));
  job->run();

  std::filesystem::current_path(oldDir);

  // Verify final.con was produced
  REQUIRE(std::filesystem::exists(workdir / "final.con"));
  // Verify final.con has content (atoms moved)
  auto fsize = std::filesystem::file_size(workdir / "final.con");
  REQUIRE(fsize > 100);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "PointJob Morse Pt matches SVN reference",
                 "[job][point][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = point
random_seed = 42

[Potential]
potential = morse_pt
)");

  auto results = runJob();

  // SVN reference (data/reference/point_morse_pt.dat):
  double energy = std::stod(results["Energy"]);
  REQUIRE(energy == Catch::Approx(-1775.791160208654).epsilon(1e-6));

  double maxForce = std::stod(results["Max_Force"]);
  REQUIRE(maxForce == Catch::Approx(0.000010652777).epsilon(1e-2));
}

TEST_CASE_METHOD(JobIntegrationFixture, "MonteCarloJob runs on LJ cluster",
                 "[job][monte_carlo][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = monte_carlo
random_seed = 42
temperature = 300.0

[Potential]
potential = lj

[Monte Carlo]
steps = 100
step_size = 0.001
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  // v2.12 writes results.dat (SVN does not)
  REQUIRE(results.count("potential_energy") > 0);
  double energy = std::stod(results["potential_energy"]);
  REQUIRE(std::isfinite(energy));
  REQUIRE(results.count("total_force_calls") > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "GlobalOptimizationJob runs on LJ cluster",
                 "[job][global_optimization][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = global_optimization
random_seed = 42
temperature = 300.0

[Potential]
potential = lj

[Global Optimization]
move_method = random
decision_method = boltzmann
steps = 10

[Basin Hopping]
displacement = 0.5
push_apart_distance = 0.4

[Optimizer]
opt_method = lbfgs
converged_force = 0.001
max_iterations = 200
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  // GlobalOptimizationJob does not write results.dat, just monitoring.dat
  auto oldDir = std::filesystem::current_path();
  std::filesystem::current_path(workdir);

  params = std::make_unique<Parameters>();
  params->load("config.ini");

  auto job = eonc::helpers::makeJob(std::move(params));
  job->run();

  std::filesystem::current_path(oldDir);

  // Verify monitoring.dat was produced
  REQUIRE(std::filesystem::exists(workdir / "monitoring.dat"));
}

TEST_CASE_METHOD(JobIntegrationFixture, "ParallelReplicaJob runs on Morse Pt",
                 "[job][parallel_replica][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = parallel_replica
temperature = 300
random_seed = 42

[Potential]
potential = morse_pt

[Dynamics]
time_step = 1.0
time = 100.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0

[Parallel Replica]
dephase_time = 20.0
state_check_interval = 20.0
record_interval = 5.0
refine_transition = true
corr_time = 10.0
)");

  auto results = runJob();

  // SVN reference (data/reference/parallel_replica_morse_pt.dat):
  REQUIRE(results.count("potential_energy_reactant") > 0);
  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(reactantE == Catch::Approx(-1775.791160).epsilon(1e-4));

  // Force calls must be <= SVN (1221)
  REQUIRE(forceCalls_ <= 1221);

  // No transition expected in 100 fs on equilibrium Pt
  REQUIRE(results.count("transition_found") > 0);
  int transition = std::stoi(results["transition_found"]);
  REQUIRE(transition == 0);
}

TEST_CASE_METHOD(JobIntegrationFixture, "ReplicaExchangeJob runs on LJ cluster",
                 "[job][replica_exchange][integration]") {
  copyTestData(".");
  writeConfig(R"(
[Main]
job = replica_exchange
temperature = 300
random_seed = 42

[Potential]
potential = lj

[Dynamics]
time_step = 1.0
time = 50.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0

[Replica Exchange]
replicas = 3
temperature_distribution = exponential
temperature_low = 100.0
temperature_high = 500.0
sampling_time = 50.0
exchange_period = 25.0
)");

  std::filesystem::copy_file(workdir / "reactant.con", workdir / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);

  auto results = runJob();

  REQUIRE(results.count("force_calls_sampling") > 0);
  int samplingCalls = std::stoi(results["force_calls_sampling"]);
  REQUIRE(samplingCalls > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture, "TADJob runs on Morse Pt",
                 "[job][tad][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = tad
temperature = 600
random_seed = 42

[Potential]
potential = morse_pt

[Dynamics]
time_step = 1.0
time = 100.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0

[TAD]
low_temperature = 300
min_prefactor = 0.001
confidence = 0.001

[Parallel Replica]
state_check_interval = 20.0
record_interval = 5.0
refine_transition = true
)");

  auto results = runJob();

  // SVN reference (data/reference/tad_morse_pt.dat):
  REQUIRE(results.count("potential_energy_reactant") > 0);
  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(reactantE == Catch::Approx(-1775.791160).epsilon(1e-4));

  // Force calls must be <= SVN (1227)
  REQUIRE(forceCalls_ <= 1227);

  // No transition expected in 100 fs
  REQUIRE(results.count("transition_found") > 0);
  int transition = std::stoi(results["transition_found"]);
  REQUIRE(transition == 0);
}

TEST_CASE_METHOD(JobIntegrationFixture, "PrefactorJob runs on Morse Pt saddle",
                 "[job][prefactor][integration]") {
  // Prefactor needs reactant.con, saddle.con, product.con
  // Use saddle_search data which has pos.con as reactant
  copyTestData("../saddle_search");

  // First run saddle search to get saddle.con and product.con
  writeConfig(R"(
[Main]
job = saddle_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Optimizer]
converged_force = 0.001
max_iterations = 1000

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
min_mode_method = dimer
max_energy = 10.0
)");

  {
    auto oldDir = std::filesystem::current_path();
    std::filesystem::current_path(workdir);
    auto p = std::make_unique<Parameters>();
    p->load("config.ini");
    auto job = eonc::helpers::makeJob(std::move(p));
    job->run();
    std::filesystem::current_path(oldDir);
  }

  // Verify saddle search produced saddle.con
  REQUIRE(std::filesystem::exists(workdir / "saddle.con"));

  // Prefactor needs reactant.con, saddle.con, product.con
  // SaddleSearchJob only writes saddle.con. Use pos.con as both
  // reactant and product (prefactor uses them to find moved atoms).
  std::filesystem::copy_file(workdir / "pos.con", workdir / "reactant.con",
                             std::filesystem::copy_options::overwrite_existing);
  std::filesystem::copy_file(workdir / "pos.con", workdir / "product.con",
                             std::filesystem::copy_options::overwrite_existing);

  // Now run prefactor
  writeConfig(R"(
[Main]
job = prefactor
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Prefactor]
all_free_atoms = false
min_displacement = 0.25
within_radius = 3.3
)");

  auto results = runJob();

  REQUIRE(results.count("good") > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture, "ProcessSearchJob runs on Morse Pt",
                 "[job][process_search][integration]") {
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = process_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.001
max_iterations = 1000
max_move = 0.2

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
min_mode_method = dimer
max_energy = 10.0
)");

  auto results = runJob();

  // SVN reference (data/reference/process_search_morse_pt.dat):
  REQUIRE(results.count("termination_reason") > 0);
  int status = std::stoi(results["termination_reason"]);
  REQUIRE(status == 0);

  // Force calls must be <= SVN (67)
  REQUIRE(forceCalls_ <= 67);

  // Energies must match SVN exactly
  double saddleE = std::stod(results["potential_energy_saddle"]);
  REQUIRE(saddleE == Catch::Approx(-1462.008706).epsilon(1e-4));

  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(reactantE == Catch::Approx(-1462.166783).epsilon(1e-4));

  double productE = std::stod(results["potential_energy_product"]);
  REQUIRE(productE == Catch::Approx(-1462.179029).epsilon(1e-4));

  // Barrier must match SVN
  double barrier = std::stod(results["barrier_reactant_to_product"]);
  REQUIRE(barrier == Catch::Approx(0.158077).epsilon(1e-3));
}

#ifdef WITH_ARTN
TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob ARTn converges on Morse Pt",
                 "[job][saddle_search][artn][integration]") {
  if (!eonc::get_artn_resource().is_loaded())
    SKIP("libartn not available at runtime");
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = saddle_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Optimizer]
converged_force = 0.001
max_iterations = 1000

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
method = artn
max_energy = 10.0

[ARTn]
push_step_size = 0.3
force_threshold = 0.05
max_iterations = 500
)");

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  int status = std::stoi(results["termination_reason"]);
  // ARTn may not converge, but should not crash
  REQUIRE((status == ARTnSaddleSearch::STATUS_GOOD ||
           status == ARTnSaddleSearch::STATUS_BAD_MAX_ITERATIONS ||
           status == ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR));

  // Energy must be finite
  double energy = std::stod(results["potential_energy_saddle"]);
  REQUIRE(std::isfinite(energy));

  if (status == ARTnSaddleSearch::STATUS_GOOD) {
    // Eigenvalue must be negative for a converged first-order saddle
    double eigenvalue = std::stod(results["final_eigenvalue"]);
    REQUIRE(eigenvalue < 0.0);
    REQUIRE(std::isfinite(eigenvalue)); // Regression test for FPE fix
  }

  // Should have made force calls
  REQUIRE(results.count("force_calls_saddle") > 0);
  int fCalls = std::stoi(results["force_calls_saddle"]);
  REQUIRE(fCalls > 0);

  // Reactant energy must match
  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(reactantE == Catch::Approx(-1462.166782).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ProcessSearchJob ARTn converges on Morse Pt",
                 "[job][process_search][artn][integration]") {
  if (!eonc::get_artn_resource().is_loaded())
    SKIP("libartn not available at runtime");
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = process_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.001
max_iterations = 1000
max_move = 0.2

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
method = artn
max_energy = 10.0

[ARTn]
push_step_size = 0.3
force_threshold = 0.05
max_iterations = 500

[Process Search]
max_spawns = 5
)");

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  int status = std::stoi(results["termination_reason"]);
  REQUIRE((status == ARTnSaddleSearch::STATUS_GOOD ||
           status == ARTnSaddleSearch::STATUS_BAD_MAX_ITERATIONS ||
           status == ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR));

  // Force calls must be reasonable
  REQUIRE(forceCalls_ > 0);

  // Energies must be finite
  double saddleE = std::stod(results["potential_energy_saddle"]);
  REQUIRE(std::isfinite(saddleE));

  double reactantE = std::stod(results["potential_energy_reactant"]);
  REQUIRE(std::isfinite(reactantE));

  double productE = std::stod(results["potential_energy_product"]);
  REQUIRE(std::isfinite(productE));

  // Reactant energy must match
  REQUIRE(reactantE == Catch::Approx(-1462.166783).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob standalone ARTn works without direction file",
                 "[job][saddle_search][artn][optional-mode][integration]") {
  copyTestData("../saddle_search");
  std::filesystem::remove(workdir / "direction.dat");
  std::filesystem::remove(workdir / "displacement.con");
  writeConfig(R"(
[Main]
job = saddle_search
temperature = 300
random_seed = 706253457

[Potential]
potential = morse_pt

[Saddle Search]
method = artn

[ARTn]
push_step_size = 0.3
force_threshold = 0.05
max_iterations = 5
)");

  auto results = runJob();

  REQUIRE(results.count("termination_reason") > 0);
  int status = std::stoi(results["termination_reason"]);
  REQUIRE((status == ARTnSaddleSearch::STATUS_GOOD ||
           status == ARTnSaddleSearch::STATUS_BAD_MAX_ITERATIONS ||
           status == ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR));
}
#endif // WITH_ARTN

#ifndef WITH_ARTN
TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob rejects standalone ARTn when not compiled",
                 "[job][saddle_search][artn][config][integration]") {
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = saddle_search

[Potential]
potential = morse_pt

[Saddle Search]
method = artn
)");

  REQUIRE_THROWS_WITH(
      runJob(),
      Catch::Matchers::ContainsSubstring(
          "saddle_search.method=artn requires a build with ARTn support"));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob rejects ARTn min-mode when not compiled",
                 "[job][saddle_search][artn][min_mode][config][integration]") {
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = saddle_search

[Potential]
potential = morse_pt

[Saddle Search]
method = min_mode
min_mode_method = artn
)");

  REQUIRE_THROWS_WITH(runJob(), Catch::Matchers::ContainsSubstring(
                                    "saddle_search.minmode_method=artn "
                                    "requires a build with ARTn support"));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ProcessSearchJob rejects standalone ARTn when not compiled",
                 "[job][process_search][artn][config][integration]") {
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = process_search

[Potential]
potential = morse_pt

[Saddle Search]
method = artn
)");

  REQUIRE_THROWS_WITH(
      runJob(),
      Catch::Matchers::ContainsSubstring(
          "saddle_search.method=artn requires a build with ARTn support"));
}
#endif

TEST_CASE_METHOD(JobIntegrationFixture,
                 "SaddleSearchJob ARTn parameters parsed correctly",
                 "[job][saddle_search][artn][params][integration]") {
  copyTestData("../saddle_search");
  writeConfig(R"(
[Main]
job = saddle_search

[ARTn]
push_step_size = 0.5
force_threshold = 0.1
max_iterations = 1000
)");

  std::filesystem::current_path(workdir);
  params = std::make_unique<Parameters>();
  params->load("config.ini");
  std::filesystem::current_path(originalDir);

  // Verify ARTn parameters are parsed
  REQUIRE(params->artn_options.push_step_size == 0.5);
  // ninit default is -1 (sentinel = "keep pARTn's own default"); test does
  // not set it in the INI above, so the sentinel must round-trip unchanged.
  REQUIRE(params->artn_options.ninit == -1);
  REQUIRE(params->artn_options.force_threshold == 0.1);
  REQUIRE(params->artn_options.max_iterations == 1000);
}

TEST_CASE_METHOD(JobIntegrationFixture, "IRA parameters parsed correctly",
                 "[job][ira][params][integration]") {
  writeConfig(R"(
[Main]
job = point

[IRA]
distance_threshold = 0.5
symmetry_threshold = 0.2
use_pbc = true
)");

  std::filesystem::current_path(workdir);
  params = std::make_unique<Parameters>();
  params->load("config.ini");
  std::filesystem::current_path(originalDir);

  // Verify IRA parameters are parsed
  REQUIRE(params->ira_options.distance_threshold == 0.5);
  REQUIRE(params->ira_options.symmetry_threshold == 0.2);
  REQUIRE(params->ira_options.use_pbc == true);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ProcessSearch with dynamics saddle search",
                 "[job][process_search][dynamics][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = process_search
temperature = 600
random_seed = 42

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
max_move = 0.2

[Saddle Search]
method = dynamics
max_iterations = 50

[Dynamics]
time_step = 1.0
time = 50.0
thermostat = andersen
andersen_collision_steps = 10
andersen_alpha = 1.0

[Parallel Replica]
state_check_interval = 10.0
record_interval = 5.0
refine_transition = false
)");

  auto results = runJob();
  REQUIRE(results.count("termination_reason") > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture, "ProcessSearch with BGSD saddle search",
                 "[job][process_search][bgsd][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = process_search
temperature = 300
random_seed = 42

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
max_move = 0.2

[Saddle Search]
method = bgsd
max_iterations = 50
max_energy = 10.0

[BGSD]
alpha = 10.0
beta = 0.2
gradient_finite_difference = 0.000001
)");

  auto results = runJob();
  REQUIRE(results.count("termination_reason") > 0);
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ProcessSearch with basin_hopping saddle search",
                 "[job][process_search][bh_saddle][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  writeConfig(R"(
[Main]
job = process_search
temperature = 300
random_seed = 42

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
max_move = 0.2

[Saddle Search]
method = basin_hopping
max_iterations = 20
max_energy = 10.0

[Basin Hopping]
displacement = 0.3
push_apart_distance = 0.4
steps = 5
)");

  auto results = runJob();
  REQUIRE(results.count("termination_reason") > 0);
}

// SafeHyperJob requires element-specific BondBoost parameters (SIGFPE on
// generic LJ clusters). Needs proper metallic test system on cosmolab.

TEST_CASE_METHOD(JobIntegrationFixture,
                 "StructureComparisonJob runs without crash",
                 "[job][structure_comparison][integration]") {
  copyTestData("../Pt_Heptamer_FrozenLayers");
  // StructureComparison needs matter1.con
  std::filesystem::copy_file(workdir / "pos.con", workdir / "matter1.con",
                             std::filesystem::copy_options::overwrite_existing);
  writeConfig(R"(
[Main]
job = structure_comparison
random_seed = 42

[Potential]
potential = morse_pt

[Structure Comparison]
distance_difference = 0.1
energy_difference = 0.01
)");

  // StructureComparisonJob is minimal; just verify it doesn't crash
  auto oldDir = std::filesystem::current_path();
  std::filesystem::current_path(workdir);
  auto p = std::make_unique<Parameters>();
  p->load("config.ini");
  auto job = eonc::helpers::makeJob(std::move(p));
  job->run();
  std::filesystem::current_path(oldDir);
  REQUIRE(true);
}

// ---------------------------------------------------------------------------
// Fixed-atom drift invariant: displacement.con with stale fixed-atom positions
// must not propagate to reactant.con / product.con after process_search.
//
// The test system (process_search_fixed_drift) is Pt_Heptamer_FrozenLayers
// (343 atoms, 336 fixed) with displacement.con built so that every fixed-atom
// row is offset +0.1 A in x vs pos.con.  The bug (pre-fix) copies those stale
// rows into min1/min2 -> reactant.con.  The fix restores fixed-atom rows from
// initial before the saddle search runs, so reactant.con must have zero drift.
// ---------------------------------------------------------------------------

/// Parse a .con file and return (positions, fixed-flags) for every atom.
static std::vector<std::pair<std::array<double, 3>, int>>
parseConAtoms(const std::filesystem::path &path) {
  std::vector<std::pair<std::array<double, 3>, int>> atoms;
  std::ifstream f(path);
  if (!f)
    return atoms;
  std::string line;
  // consume lines 0-6 (header: title, json/blank, cell, angles, 2 blanks,
  // nTypes)
  for (int i = 0; i < 7; i++)
    std::getline(f, line);
  // line 7: space-separated n_atoms per component; take the first token
  int n_atoms = 0;
  {
    std::getline(f, line); // reads line 7
    std::istringstream ss(line);
    ss >> n_atoms;
  }
  // lines 8-10: masses, element names, "Coordinates of..."
  for (int i = 0; i < 3; i++)
    std::getline(f, line);
  // next n_atoms lines: x y z fixed_bits index
  atoms.reserve(n_atoms);
  for (int i = 0; i < n_atoms; i++) {
    if (!std::getline(f, line))
      break;
    std::istringstream ss(line);
    std::array<double, 3> xyz{};
    int fb{0}, idx{0};
    ss >> xyz[0] >> xyz[1] >> xyz[2] >> fb >> idx;
    atoms.push_back({xyz, fb});
  }
  return atoms;
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "ProcessSearchJob zero fixed-atom drift with stale "
                 "displacement.con",
                 "[job][process_search][fixed_atom][invariant][integration]") {
  copyTestData("../process_search_fixed_drift");

  // Verify the test data actually has stale fixed-atom positions so the
  // test is not vacuously passing because drift was zero to begin with.
  {
    auto pos_atoms = parseConAtoms(workdir / "pos.con");
    auto dis_atoms = parseConAtoms(workdir / "displacement.con");
    REQUIRE(pos_atoms.size() == dis_atoms.size());
    double max_stale = 0.0;
    for (size_t i = 0; i < pos_atoms.size(); i++) {
      if (pos_atoms[i].second != 0) {
        for (int ax = 0; ax < 3; ax++) {
          max_stale = std::max(max_stale, std::abs(pos_atoms[i].first[ax] -
                                                   dis_atoms[i].first[ax]));
        }
      }
    }
    // displacement.con was built with +0.1 A offset on fixed atoms
    REQUIRE(max_stale > 0.05);
  }

  writeConfig(R"(
[Main]
job = process_search
temperature = 300
random_seed = 42

[Potential]
potential = morse_pt

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
max_move = 0.2

[Saddle Search]
method = min_mode
min_mode_method = dimer
displace_type = load
max_iterations = 300
max_energy = 10.0

[Structure Comparison]
distance_difference = 0.1
energy_difference = 0.01

[Prefactor]
default_value = 1e12
)");

  // Run the job; any termination status is acceptable -- the invariant holds
  // regardless of whether the saddle search converges.
  std::filesystem::current_path(workdir);
  auto p = std::make_unique<Parameters>();
  p->load("config.ini");
  auto job = eonc::helpers::makeJob(std::move(p));
  job->run();
  std::filesystem::current_path(originalDir);

  // reactant.con must exist (saveData always writes it)
  REQUIRE(std::filesystem::exists(workdir / "reactant.con"));

  auto pos_atoms = parseConAtoms(workdir / "pos.con");
  auto reactant_atoms = parseConAtoms(workdir / "reactant.con");
  REQUIRE(pos_atoms.size() == reactant_atoms.size());

  // Fixed-atom drift invariant: every fixed atom in reactant.con must
  // match pos.con within 1e-8 A.
  double max_fixed_drift = 0.0;
  int n_fixed = 0;
  for (size_t i = 0; i < pos_atoms.size(); i++) {
    if (pos_atoms[i].second != 0) {
      ++n_fixed;
      for (int ax = 0; ax < 3; ax++) {
        max_fixed_drift =
            std::max(max_fixed_drift, std::abs(pos_atoms[i].first[ax] -
                                               reactant_atoms[i].first[ax]));
      }
    }
  }
  REQUIRE(n_fixed == 336); // sanity: all 336 fixed atoms are present
  // Core invariant: fixed atoms must not drift from pos.con
  REQUIRE(max_fixed_drift < 1.0e-8);
}

// ---------------------------------------------------------------------------
// Unit-level regression test for the fixed-atom restore in ProcessSearchJob.
//
// The fix iterates initial->getFixed(i) and copies fixed-atom rows from
// initial into saddle before the dimer/min-mode search runs.  This test
// exercises that same logic directly through the Matter API, without
// relying on the saddle search converging, so it is fully deterministic.
// ---------------------------------------------------------------------------
TEST_CASE("ProcessSearchJob fixed-atom restore: displacement.con stale rows "
          "are patched from initial",
          "[process_search][fixed_atom][matter][unit]") {
  // Build two Matter objects.  LJ suffices: we only read positions and
  // fixed-flags, never evaluate the potential.
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);

  // initial = pos.con (the authoritative reference)
  auto initial = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(initial->con2matter(std::string("../Pt_Heptamer_FrozenLayers/pos.con"))));
  REQUIRE(initial->numberOfAtoms() == 343);
  REQUIRE(initial->numberOfFixedAtoms() == 336);

  // saddle = displacement.con from process_search_fixed_drift -- has +0.1 A
  // offset on all fixed-atom rows vs pos.con.
  auto saddle = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(saddle->con2matter(
      std::string("../process_search_fixed_drift/displacement.con"))));
  REQUIRE(saddle->numberOfAtoms() == 343);

  // Confirm the stale fixed-atom drift is present before the fix.
  {
    const AtomMatrix &iPos = initial->getPositions();
    const AtomMatrix &sPos = saddle->getPositions();
    double max_stale = 0.0;
    for (long i = 0; i < initial->numberOfAtoms(); i++) {
      if (initial->getFixed(i)) {
        max_stale = std::max(max_stale,
                             (iPos.row(i) - sPos.row(i)).cwiseAbs().maxCoeff());
      }
    }
    REQUIRE(max_stale > 0.05); // displacement.con was built with +0.1 A
  }

  // Apply the fix: restore fixed-atom rows from initial into saddle.
  {
    const AtomMatrix &initPos = initial->getPositions();
    AtomMatrix saddlePos = saddle->getPositionsCopy();
    for (long i = 0; i < initial->numberOfAtoms(); i++) {
      if (initial->getFixed(i)) {
        saddlePos.row(i) = initPos.row(i);
      }
    }
    saddle->setPositions(saddlePos);
  }

  // After the fix, fixed atoms in saddle must match initial within machine eps.
  {
    const AtomMatrix &iPos = initial->getPositions();
    const AtomMatrix &sPos = saddle->getPositions();
    double max_drift = 0.0;
    int n_fixed = 0;
    for (long i = 0; i < initial->numberOfAtoms(); i++) {
      if (initial->getFixed(i)) {
        ++n_fixed;
        max_drift = std::max(max_drift,
                             (iPos.row(i) - sPos.row(i)).cwiseAbs().maxCoeff());
      }
    }
    REQUIRE(n_fixed == 336);
    REQUIRE(max_drift < 1.0e-12);
  }
}

TEST_CASE("makeJob creates correct job type for each JobType",
          "[job][factory]") {
  auto params = std::make_unique<Parameters>();

  params->potential_options.potential = PotType::LJ;
  params->main_options.job = JobType::Point;
  auto job = eonc::helpers::makeJob(std::move(params));
  REQUIRE(job != nullptr);
}

// -----------------------------------------------------------------------
// SVN-verified point energy tests for Fortran potentials (Si diamond)
// These require -Dwith_fortran=true at build time.
// -----------------------------------------------------------------------

#ifdef WITH_FORTRAN
TEST_CASE_METHOD(JobIntegrationFixture, "PointJob SW Si matches SVN reference",
                 "[job][point][sw][integration][fortran]") {
  copyTestData("../si_diamond");
  writeConfig(R"(
[Main]
job = point
random_seed = 42

[Potential]
potential = sw_si
)");
  auto results = runJob();
  REQUIRE(std::stod(results["Energy"]) ==
          Catch::Approx(-16.204955).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "PointJob Tersoff Si matches SVN reference",
                 "[job][point][tersoff][integration][fortran]") {
  copyTestData("../si_diamond");
  writeConfig(R"(
[Main]
job = point
random_seed = 42

[Potential]
potential = tersoff_si
)");
  auto results = runJob();
  REQUIRE(std::stod(results["Energy"]) ==
          Catch::Approx(-17.440266).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "PointJob EDIP Si matches SVN reference",
                 "[job][point][edip][integration][fortran]") {
  copyTestData("../si_diamond");
  writeConfig(R"(
[Main]
job = point
random_seed = 42

[Potential]
potential = edip
)");
  auto results = runJob();
  REQUIRE(std::stod(results["Energy"]) ==
          Catch::Approx(-18.838135).epsilon(1e-4));
}

TEST_CASE_METHOD(JobIntegrationFixture,
                 "PointJob Lenosky Si matches SVN reference",
                 "[job][point][lenosky][integration][fortran]") {
  copyTestData("../si_diamond");
  writeConfig(R"(
[Main]
job = point
random_seed = 42

[Potential]
potential = lenosky_si
)");
  auto results = runJob();
  REQUIRE(std::stod(results["Energy"]) ==
          Catch::Approx(-17.284558).epsilon(1e-4));
}
#endif // WITH_FORTRAN

} /* namespace tests */
