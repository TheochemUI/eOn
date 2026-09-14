#include "eon/ParametersSSOT.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>

namespace tests {
static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Parameters constructor applies Cap'n Proto SSoT defaults",
          "[params][ssot]") {
  Parameters p;
  // From schema/eon_params.capnp MainOptions / PotentialOptions
  REQUIRE(p.main_options.temperature == Catch::Approx(300.0));
  REQUIRE(p.main_options.finiteDifference == Catch::Approx(0.01));
  REQUIRE(p.main_options.job == JobType::Process_Search);
  REQUIRE(p.potential_options.potential == PotType::LJ);
  REQUIRE(p.optimizer_options.method == OptType::CG);
  REQUIRE(p.optimizer_options.max_iterations == 1000);
  REQUIRE(p.process_search_options.minimize_first == true);
  REQUIRE(p.structure_comparison_options.distance_difference ==
          Catch::Approx(0.1));
}

TEST_CASE("ssot_has_field knows covered catalog keys", "[params][ssot]") {
  REQUIRE(eonc::config::ssot_has_field("Main", "temperature"));
  REQUIRE(eonc::config::ssot_has_field("Potential", "potential"));
  REQUIRE_FALSE(eonc::config::ssot_has_field("Main", "not_a_real_field"));
  REQUIRE_FALSE(eonc::config::ssot_has_field("Dimer", "rotation_angle"));
}

TEST_CASE("Parameters::load INI still overrides SSoT defaults",
          "[params][ssot][ini]") {
  namespace fs = std::filesystem;
  auto dir = fs::temp_directory_path() / "eon_ssot_ini_XXXXXX";
  // unique path
  dir = fs::temp_directory_path() /
        ("eon_ssot_ini_" + std::to_string(std::rand()));
  fs::create_directories(dir);
  auto ini = dir / "config.ini";
  {
    std::ofstream out(ini);
    out << "[Main]\njob = minimization\ntemperature = 450.0\n"
           "random_seed = 7\n"
           "[Potential]\npotential = eam_al\n"
           "[Optimizer]\nopt_method = lbfgs\nmax_iterations = 42\n";
  }
  Parameters p;
  REQUIRE(p.load(ini.string()) == 0);
  REQUIRE(p.main_options.job == JobType::Minimization);
  REQUIRE(p.main_options.temperature == Catch::Approx(450.0));
  REQUIRE(p.main_options.randomSeed == 7);
  REQUIRE(p.potential_options.potential == PotType::EAM_AL);
  REQUIRE(p.optimizer_options.method == OptType::LBFGS);
  REQUIRE(p.optimizer_options.max_iterations == 42);
  fs::remove_all(dir);
}
} // namespace tests
