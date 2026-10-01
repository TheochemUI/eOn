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

#include "eon/ParametersJSON.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Parameters round-trips through JSON", "[params][json]") {
  Parameters p1;
  ParametersLoadAccess::potential_options(p1).potential = PotType::MORSE_PT;
  ParametersLoadAccess::main_options(p1).job = JobType::Minimization;
  ParametersLoadAccess::main_options(p1).temperature = 500.0;
  ParametersLoadAccess::main_options(p1).randomSeed = 12345;
  ParametersLoadAccess::optimizer_options(p1).method = OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(p1).converged_force = 0.005;
  ParametersLoadAccess::optimizer_options(p1).max_iterations = 999;
  ParametersLoadAccess::neb_options(p1).image_count = 7;
  ParametersLoadAccess::neb_options(p1).force_tolerance = 0.02;

  auto j = eonc::config::to_json(p1);

  Parameters p2;
  eonc::config::from_json(j, p2);

  REQUIRE(p2.potential_options().potential == PotType::MORSE_PT);
  REQUIRE(p2.main_options().job == JobType::Minimization);
  REQUIRE(p2.main_options().temperature == Catch::Approx(500.0));
  REQUIRE(p2.main_options().randomSeed == 12345);
  REQUIRE(p2.optimizer_options().method == OptType::LBFGS);
  REQUIRE(p2.optimizer_options().converged_force == Catch::Approx(0.005));
  REQUIRE(p2.optimizer_options().max_iterations == 999);
  REQUIRE(p2.neb_options().image_count == 7);
}

TEST_CASE("JSON to_json produces valid JSON string", "[params][json]") {
  Parameters p;
  ParametersLoadAccess::potential_options(p).potential = PotType::LJ;
  ParametersLoadAccess::main_options(p).job = JobType::Point;

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

  REQUIRE(p.main_options().job == JobType::Point);
  REQUIRE(p.potential_options().potential == PotType::LJ);
}

TEST_CASE("JSON round-trip preserves saddle search options", "[params][json]") {
  Parameters p1;
  ParametersLoadAccess::saddle_search_options(p1).max_energy = 15.0;
  ParametersLoadAccess::saddle_search_options(p1).max_iterations = 42;
  ParametersLoadAccess::saddle_search_options(p1).displace_radius = 4.5;

  auto j = eonc::config::to_json(p1);
  Parameters p2;
  eonc::config::from_json(j, p2);

  REQUIRE(p2.saddle_search_options().max_energy == Catch::Approx(15.0));
  REQUIRE(p2.saddle_search_options().max_iterations == 42);
  REQUIRE(p2.saddle_search_options().displace_radius == Catch::Approx(4.5));
}

TEST_CASE("JSON to_json includes dynamics section", "[params][json]") {
  Parameters p1;
  ParametersLoadAccess::dynamics_options(p1).time = 500.0;

  auto j = eonc::config::to_json(p1);
  std::string s = j.dump();
  REQUIRE(s.find("Dynamics") != std::string::npos);
}

TEST_CASE("JSON round-trip preserves dimer rotation_backend",
          "[params][json][lor]") {
  Parameters p1;
  ParametersLoadAccess::dimer_options(p1).rotation_backend =
      DimerRotationBackend::LOR;

  auto j = eonc::config::to_json(p1);
  REQUIRE(j["Dimer"].contains("rotation_backend"));

  Parameters p2;
  eonc::config::from_json(j, p2);
  REQUIRE(p2.dimer_options().rotation_backend == DimerRotationBackend::LOR);

  // Case-insensitive load (serve / Python schema use lowercase tokens)
  nlohmann::json j_lor = {{"Dimer", {{"rotation_backend", "lor"}}}};
  Parameters p3;
  eonc::config::from_json(j_lor, p3);
  REQUIRE(p3.dimer_options().rotation_backend == DimerRotationBackend::LOR);
}

TEST_CASE("JSON round-trip preserves debug compatibility flags",
          "[params][json]") {
  Parameters p1;
  ParametersLoadAccess::debug_options(p1).write_movies = true;
  ParametersLoadAccess::debug_options(p1).write_movies_interval = 4;
  ParametersLoadAccess::debug_options(p1).write_deprecated_outs = true;

  auto j = eonc::config::to_json(p1);
  Parameters p2;
  eonc::config::from_json(j, p2);

  REQUIRE(j["Debug"]["write_deprecated_outs"] == true);
  REQUIRE(p2.debug_options().write_movies == true);
  REQUIRE(p2.debug_options().write_movies_interval == 4);
  REQUIRE(p2.debug_options().write_deprecated_outs == true);
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

  REQUIRE(p.main_options().job == JobType::Saddle_Search);
  REQUIRE(p.main_options().temperature == Catch::Approx(500.0));
  REQUIRE(p.main_options().randomSeed == 12345);
  REQUIRE(p.potential_options().potential == PotType::MORSE_PT);
  REQUIRE(p.optimizer_options().method == OptType::LBFGS);
  REQUIRE(p.optimizer_options().converged_force == Catch::Approx(0.005));
  REQUIRE(p.optimizer_options().max_iterations == 999);
  REQUIRE(p.saddle_search_options().max_energy == Catch::Approx(15.0));
  REQUIRE(p.neb_options().image_count == 7);
  REQUIRE(p.neb_options().spring.constant == Catch::Approx(3.0));
  REQUIRE(p.structure_comparison_options().distance_difference ==
          Catch::Approx(0.2));

  std::filesystem::remove(tmpfile);
}

TEST_CASE("Parameters load APIs accept string_view", "[params][ini]") {
  constexpr std::string_view ini = "[Main]\njob = point\n";
  Parameters p;
  REQUIRE(p.load_ini_text(ini) == 0);
  REQUIRE(p.main_options().job == JobType::Point);

  constexpr std::string_view js = R"({"Main":{"job":"minimization"}})";
  REQUIRE(p.load_json(js) == 0);
  REQUIRE(p.main_options().job == JobType::Minimization);
}

TEST_CASE("JSON from_json rejects unknown convergence_metric",
          "[params][json]") {
  nlohmann::json j = nlohmann::json::object();
  j["Optimizer"]["convergence_metric"] = "typo";
  Parameters p;
  REQUIRE_THROWS_AS(eonc::config::from_json(j, p), std::invalid_argument);
  REQUIRE(eonc::config::load_json(j.dump(), p) != 0);
}

TEST_CASE("JSON from_json accepts known convergence_metric", "[params][json]") {
  nlohmann::json j = nlohmann::json::object();
  j["Optimizer"]["convergence_metric"] = "Max_Atom";
  Parameters p;
  eonc::config::from_json(j, p);
  REQUIRE(p.optimizer_options().convergence_metric == "max_atom");
  REQUIRE(p.optimizer_options().convergence_metric_label == "Max atom force");
}

TEST_CASE("Parameters.load rejects unknown enumerated strings",
          "[params][ini]") {
  const auto tmpdir = std::filesystem::temp_directory_path();
  {
    const auto tmpfile = tmpdir / "eon_bad_metric.ini";
    {
      std::ofstream ofs(tmpfile);
      ofs << "[Optimizer]\nconvergence_metric = typo\n";
    }
    Parameters p;
    REQUIRE(p.load(tmpfile.string()) != 0);
    std::filesystem::remove(tmpfile);
  }
  {
    const auto tmpfile = tmpdir / "eon_bad_disp_algo.ini";
    {
      std::ofstream ofs(tmpfile);
      ofs << "[Basin Hopping]\ndisplacement_algorithm = spiral\n";
    }
    Parameters p;
    REQUIRE(p.load(tmpfile.string()) != 0);
    std::filesystem::remove(tmpfile);
  }
  {
    const auto tmpfile = tmpdir / "eon_bad_disp_dist.ini";
    {
      std::ofstream ofs(tmpfile);
      ofs << "[Basin Hopping]\ndisplacement_distribution = cauchy\n";
    }
    Parameters p;
    REQUIRE(p.load(tmpfile.string()) != 0);
    std::filesystem::remove(tmpfile);
  }
  {
    const auto tmpfile = tmpdir / "eon_bad_decision.ini";
    {
      std::ofstream ofs(tmpfile);
      ofs << "[Global Optimization]\ndecision_method = coin_flip\n";
    }
    Parameters p;
    REQUIRE(p.load(tmpfile.string()) != 0);
    std::filesystem::remove(tmpfile);
  }
}

TEST_CASE("JSON reads path-integral keys from Dynamics and Thermostat",
          "[params][json]") {
  nlohmann::json dynamics = {
      {"Dynamics",
       {{"path_beads", 16},
        {"path_springs", "Eco"},
        {"path_pile_tau", 50.0},
        {"path_seed", 3},
        {"thermostat", "Pile"}}},
  };
  Parameters fromDynamics;
  eonc::config::from_json(dynamics, fromDynamics);
  REQUIRE(fromDynamics.thermostat_options().path_beads == 16);
  REQUIRE(fromDynamics.thermostat_options().path_springs == "eco");
  REQUIRE(fromDynamics.thermostat_options().path_pile_tau ==
          Catch::Approx(50.0 / fromDynamics.constants().timeUnit));
  REQUIRE(fromDynamics.thermostat_options().path_seed == 3);
  REQUIRE(fromDynamics.thermostat_options().kind == "pile");

  nlohmann::json thermostatWins = {
      {"Dynamics", {{"path_beads", 16}, {"path_springs", "trotter"}}},
      {"Thermostat", {{"path_beads", 4}, {"path_springs", "ECO"}}},
  };
  Parameters overridden;
  eonc::config::from_json(thermostatWins, overridden);
  REQUIRE(overridden.thermostat_options().path_beads == 4);
  REQUIRE(overridden.thermostat_options().path_springs == "eco");

  nlohmann::json badSprings = {{"Dynamics", {{"path_springs", "foo"}}}};
  Parameters bad;
  REQUIRE_THROWS_AS(eonc::config::from_json(badSprings, bad),
                    std::invalid_argument);

  nlohmann::json instanton = {
      {"Instanton",
       {{"springs", "Trotter"},
        {"beads", 32},
        {"mode", "rate"},
        {"initial_hessians", "finite_difference"},
        {"temperatures", "10, 20"}}},
  };
  Parameters fromInstanton;
  eonc::config::from_json(instanton, fromInstanton);
  REQUIRE(fromInstanton.instanton_options().springs == "trotter");
  REQUIRE(fromInstanton.instanton_options().initial_hessians ==
          "finite_difference");
  REQUIRE(fromInstanton.instanton_options().beads == 32);
  REQUIRE(fromInstanton.instanton_options().mode == "rate");
  REQUIRE(fromInstanton.instanton_options().temperatures.size() == 2);
  REQUIRE(fromInstanton.instanton_options().temperatures[0] ==
          Catch::Approx(10.0));
  REQUIRE(fromInstanton.instanton_options().temperatures[1] ==
          Catch::Approx(20.0));

  nlohmann::json listed = {
      {"Instanton", {{"temperatures", {5.0, 15.0}}}},
  };
  Parameters fromList;
  eonc::config::from_json(listed, fromList);
  REQUIRE(fromList.instanton_options().temperatures.size() == 2);
  REQUIRE(fromList.instanton_options().temperatures[0] == Catch::Approx(5.0));
  REQUIRE(fromList.instanton_options().temperatures[1] == Catch::Approx(15.0));

  Parameters written;
  ParametersLoadAccess::thermostat_options(written).path_beads = 12;
  ParametersLoadAccess::instanton_options(written).springs = "trotter";
  auto roundTrip = eonc::config::to_json(written);
  Parameters loaded;
  eonc::config::from_json(roundTrip, loaded);
  REQUIRE(loaded.thermostat_options().path_beads == 12);
  REQUIRE(loaded.instanton_options().springs == "trotter");
}

TEST_CASE("JSON round-trips every Instanton and Hessian key",
          "[params][json]") {
  Parameters written;
  auto &o = ParametersLoadAccess::instanton_options(written);
  o.mode = "rate";
  o.reactant_filename = "r.con";
  o.product_filename = "p.con";
  o.initial_path = "band.con";
  o.beads = 48;
  o.beta_hbar_omega = 12.5;
  o.max_iterations = 77;
  o.force_tolerance = 2e-4;
  o.hessian_stride = 3;
  o.saddle_filename = "s.con";
  o.temperature = 150.0;
  o.temperatures = {200.0, 100.0};
  o.half_ring = false;
  o.initial_hessians = "finite_difference";
  o.energy_shift = -1.25;
  o.bead_ladder = true;
  o.hessian_final = "interpolated";
  o.springs = "trotter";
  auto &h = ParametersLoadAccess::hessian_options(written);
  h.phva_atoms = "0,1,2";
  h.zero_freq_value = 2e-3;
  h.fd_scheme = "central";
  h.resume = true;
  h.checkpoint_path = "hessian.ckpt";
  h.write_modes = false;

  Parameters loaded;
  eonc::config::from_json(eonc::config::to_json(written), loaded);

  const auto &l = loaded.instanton_options();
  REQUIRE(l.mode == "rate");
  REQUIRE(l.reactant_filename == "r.con");
  REQUIRE(l.product_filename == "p.con");
  REQUIRE(l.initial_path == "band.con");
  REQUIRE(l.beads == 48);
  REQUIRE(l.beta_hbar_omega == Catch::Approx(12.5));
  REQUIRE(l.max_iterations == 77);
  REQUIRE(l.force_tolerance == Catch::Approx(2e-4));
  REQUIRE(l.hessian_stride == 3);
  REQUIRE(l.saddle_filename == "s.con");
  REQUIRE(l.temperature == Catch::Approx(150.0));
  REQUIRE(l.temperatures == std::vector<double>{200.0, 100.0});
  REQUIRE_FALSE(l.half_ring);
  REQUIRE(l.initial_hessians == "finite_difference");
  REQUIRE(l.energy_shift == Catch::Approx(-1.25));
  REQUIRE(l.bead_ladder);
  REQUIRE(l.hessian_final == "interpolated");
  REQUIRE(l.springs == "trotter");
  const auto &lh = loaded.hessian_options();
  REQUIRE(lh.phva_atoms == "0,1,2");
  REQUIRE(lh.zero_freq_value == Catch::Approx(2e-3));
  REQUIRE(lh.fd_scheme == "central");
  REQUIRE(lh.resume);
  REQUIRE(lh.checkpoint_path == "hessian.ckpt");
  REQUIRE_FALSE(lh.write_modes);
}

TEST_CASE("JSON derives the thermostat times the way the ini does",
          "[params][json]") {
  // No path_pile_tau or andersen_collision_period: both loaders convert
  // the femtosecond defaults to internal time.
  Parameters fromIni;
  REQUIRE(fromIni.load_ini_text("[Dynamics]\nthermostat = pile\n") == 0);
  Parameters fromJson;
  REQUIRE(fromJson.load_json(R"({"Dynamics":{"thermostat":"pile"}})") == 0);
  REQUIRE(fromIni.thermostat_options().path_pile_tau > 0.0);
  REQUIRE(fromJson.thermostat_options().path_pile_tau ==
          Catch::Approx(fromIni.thermostat_options().path_pile_tau));
  REQUIRE(fromJson.thermostat_options().andersen_tcol ==
          Catch::Approx(fromIni.thermostat_options().andersen_tcol));

  // A written file loads back with the same internal values.
  Parameters roundTrip;
  eonc::config::from_json(eonc::config::to_json(fromIni), roundTrip);
  REQUIRE(roundTrip.thermostat_options().path_pile_tau ==
          Catch::Approx(fromIni.thermostat_options().path_pile_tau));
  REQUIRE(roundTrip.thermostat_options().andersen_tcol ==
          Catch::Approx(fromIni.thermostat_options().andersen_tcol));
}

TEST_CASE("JSON round-trips RgpotPot canonical keys", "[params][json]") {
  Parameters written;
  auto &o = ParametersLoadAccess::rgpot_options(written);
  o.backend = "cpmdc";
  o.functional = "PBE";
  o.cutoff_ry = 55.5;
  o.charge = 4;
  o.params_path = "/data/si3n4.bin";
  o.ranks_per_image = 6;

  auto j = eonc::config::to_json(written);
  REQUIRE(j["RgpotPot"].contains("cutoff_ry"));
  REQUIRE_FALSE(j["RgpotPot"].contains("cutOffRy"));

  Parameters loaded;
  eonc::config::from_json(j, loaded);
  REQUIRE(loaded.rgpot_options().backend == "cpmdc");
  REQUIRE(loaded.rgpot_options().functional == "PBE");
  REQUIRE(loaded.rgpot_options().cutoff_ry == Catch::Approx(55.5));
  REQUIRE(loaded.rgpot_options().charge == 4);
  REQUIRE(loaded.rgpot_options().params_path == "/data/si3n4.bin");
  REQUIRE(loaded.rgpot_options().ranks_per_image == 6);
}

TEST_CASE("JSON reads RgpotPot cutOffRy and cpmd_functional aliases",
          "[params][json]") {
  nlohmann::json j = {
      {"RgpotPot", {{"cutOffRy", 18.0}, {"cpmd_functional", "PBE"}}}};
  Parameters p;
  eonc::config::from_json(j, p);
  REQUIRE(p.rgpot_options().cutoff_ry == Catch::Approx(18.0));
  REQUIRE(p.rgpot_options().functional == "PBE");
}

TEST_CASE("JSON RgpotPot omits xtb_charge resets it from charge",
          "[params][json]") {
  Parameters reset;
  ParametersLoadAccess::rgpot_options(reset).charge = 5;
  ParametersLoadAccess::rgpot_options(reset).xtb_charge = 9.0;
  nlohmann::json missing = {{"RgpotPot", {{"charge", 5}}}};
  eonc::config::from_json(missing, reset);
  REQUIRE(reset.rgpot_options().xtb_charge == Catch::Approx(5.0));

  Parameters kept;
  ParametersLoadAccess::rgpot_options(kept).charge = 5;
  ParametersLoadAccess::rgpot_options(kept).xtb_charge = 9.0;
  nlohmann::json present = {
      {"RgpotPot", {{"charge", 5}, {"xtb_charge", 2.5}}}};
  eonc::config::from_json(present, kept);
  REQUIRE(kept.rgpot_options().xtb_charge == Catch::Approx(2.5));
}

} /* namespace tests */
