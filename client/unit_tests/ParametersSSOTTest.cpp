#include "eon/ParametersSSOT.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <utility>

namespace tests {
static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Parameters constructor applies Cap'n Proto SSoT defaults",
          "[params][ssot]") {
  Parameters p;
  // From schema/eon_params.capnp MainOptions / PotentialOptions
  REQUIRE(p.main_options().temperature == Catch::Approx(300.0));
  REQUIRE(p.main_options().finiteDifference == Catch::Approx(0.01));
  REQUIRE(p.main_options().job == JobType::Process_Search);
  REQUIRE(p.potential_options().potential == PotType::LJ);
  REQUIRE(p.optimizer_options().method == OptType::CG);
  REQUIRE(p.optimizer_options().max_iterations == 1000);
  REQUIRE(p.process_search_options().minimize_first == true);
  REQUIRE(p.structure_comparison_options().distance_difference ==
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
  REQUIRE(p.main_options().job == JobType::Minimization);
  REQUIRE(p.main_options().temperature == Catch::Approx(450.0));
  REQUIRE(p.main_options().randomSeed == 7);
  REQUIRE(p.potential_options().potential == PotType::EAM_AL);
  REQUIRE(p.optimizer_options().method == OptType::LBFGS);
  REQUIRE(p.optimizer_options().max_iterations == 42);
  REQUIRE(p.last_load_source() == ini.string());
  REQUIRE(p.last_load_error() == 0);
  Parameters copied = p;
  REQUIRE(copied.last_load_source() == p.last_load_source());
  REQUIRE(copied.last_load_error() == 0);
  REQUIRE(copied.main_options().temperature == Catch::Approx(450.0));
  Parameters assigned;
  assigned = p;
  REQUIRE(assigned.last_load_source() == p.last_load_source());
  fs::remove_all(dir);
}

TEST_CASE("Parameters::load accepts AMS DFTB and FORCEFIELD engines",
          "[params][ams]") {
  namespace fs = std::filesystem;
  auto dir = fs::temp_directory_path() /
             ("eon_ams_engine_" + std::to_string(std::rand()));
  fs::create_directories(dir);
  const auto ini = dir / "config.ini";
  const auto load_body = [&](const char *body) {
    {
      std::ofstream out(ini);
      out << body;
    }
    Parameters loaded;
    const int rc = loaded.load(ini.string());
    return std::make_pair(rc, std::move(loaded));
  };

  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = DFTB\n"
                                        "resources = DFTB\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.potential_options().potential == PotType::AMS);
    REQUIRE(loaded.ams_options().engine == "DFTB");
    REQUIRE(loaded.ams_options().resources == "DFTB");
    REQUIRE(loaded.ams_options().forcefield.empty());
    REQUIRE(loaded.ams_options().model.empty());
    REQUIRE(loaded.ams_options().xc.empty());
    REQUIRE(loaded.last_load_error() == 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = forcefield\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.ams_options().engine == "forcefield");
    REQUIRE(loaded.last_load_error() == 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams_io\n"
                                        "[AMS_IO]\n"
                                        "engine = FORCEFIELD\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.potential_options().potential == PotType::AMS_IO);
    REQUIRE(loaded.ams_options().engine == "FORCEFIELD");
    REQUIRE(loaded.last_load_error() == 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = DFTB\n");
    REQUIRE(rc != 0);
    REQUIRE(loaded.last_load_error() != 0);
  }
  {
    const auto rejected = load_body("[Potential]\n"
                                    "potential = ams\n"
                                    "[AMS]\n"
                                    "engine = MOPAC\n");
    REQUIRE(rejected.first != 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = MOPAC\n"
                                        "model = PM3\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.ams_options().model == "PM3");
  }
  {
    const auto rejected = load_body("[Potential]\n"
                                    "potential = ams\n"
                                    "[AMS]\n"
                                    "engine = MOPAC\n"
                                    "forcefield = ff\n"
                                    "model = PM3\n"
                                    "xc = LDA\n");
    REQUIRE(rejected.first != 0);
  }
  fs::remove_all(dir);
}

TEST_CASE("Parameters load-state Impl records a missing file",
          "[params][pimpl]") {
  Parameters p;
  REQUIRE(p.last_load_source().empty());
  REQUIRE(p.last_load_error() == 0);
  REQUIRE(p.load("no-such-eon-config.ini") != 0);
  REQUIRE(p.last_load_source() == "no-such-eon-config.ini");
  REQUIRE(p.last_load_error() != 0);
  Parameters moved = std::move(p);
  REQUIRE(moved.last_load_source() == "no-such-eon-config.ini");
  REQUIRE(p.last_load_source().empty());
}
} // namespace tests
