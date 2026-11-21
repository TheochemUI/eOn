/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/LinebJob.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/LinebSeed.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/Tunneling.h"

#ifdef WITH_RGSADDLE_LINEB
#include <rgsaddle.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("linebAvailable matches the linked rgsaddle lineb C API",
          "[lineb][rgsaddle]") {
#ifdef WITH_RGSADDLE_LINEB
  REQUIRE(eonc::linebAvailable());
#else
  REQUIRE_FALSE(eonc::linebAvailable());
#endif
}

TEST_CASE("lineb config default fills paper spring and endpoint scale",
          "[lineb][rgsaddle]") {
#ifndef WITH_RGSADDLE_LINEB
  SKIP("rgsaddle without rgsaddle_lineb_solve");
#else
  rgsaddle_lineb_config_t cfg{};
  REQUIRE(rgsaddle_lineb_config_default(1e-3, 2e-3, &cfg) == RGSADDLE_OK);
  REQUIRE(cfg.opt_tol == Catch::Approx(1e-3));
  REQUIRE(cfg.perp_tol == Catch::Approx(2e-3));
  REQUIRE(cfg.spring == Catch::Approx(0.1));
  REQUIRE(cfg.endpoint_scale == Catch::Approx(0.2));
  REQUIRE(cfg.time_step == Catch::Approx(0.1));
  REQUIRE(cfg.max_move == Catch::Approx(0.2));
  REQUIRE(cfg.turning_width == Catch::Approx(1e-3));
  REQUIRE(cfg.quadrature_pairs == 4);
#endif
}

TEST_CASE("LinebJob refuses a build without the lineb C API", "[lineb]") {
  if (eonc::linebAvailable()) {
    SKIP("lineb C API is present");
  }
  Parameters params;
  ParametersLoadAccess::main_options(params).job = JobType::Lineb;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::potential_options(params).potential = PotType::MORSE_PT;
  ParametersLoadAccess::lineb_options(params).initial_path = "neb.con";
  ParametersLoadAccess::lineb_options(params).energy = 0.1;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::MORSE_PT, params));
  eonc::LinebJob job(pot, params);
  REQUIRE_THROWS_WITH(
      job.run(), Catch::Matchers::ContainsSubstring("rgsaddle_lineb_solve"));
}

TEST_CASE("LinebJob runs a fixed-energy seeded Morse band when available",
          "[lineb][rgsaddle][integration]") {
  if (!eonc::linebAvailable()) {
    SKIP("rgsaddle without rgsaddle_lineb_solve");
  }

  Parameters params;
  ParametersLoadAccess::main_options(params).job = JobType::Lineb;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::potential_options(params).potential = PotType::MORSE_PT;

  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::MORSE_PT, params));
  auto &lineb = ParametersLoadAccess::lineb_options(params);
  lineb.opt_tol = 1e-2;
  lineb.perp_tol = 1e-2;
  lineb.max_steps = 200;
  lineb.beads = 0;
  lineb.molecular = false;
  lineb.hbar = eonc::tunneling::kHbar;
  lineb.boltzmann = eonc::tunneling::kBoltzmann;

  SECTION("a one-sided cluster displacement cannot produce a tunneling path") {
    auto reactant = std::make_shared<Matter>(pot, params);
    auto product = std::make_shared<Matter>(pot, params);
    REQUIRE(eonc::io::io_ok(reactant->con2matter(std::string("reactant.con"))));
    REQUIRE(eonc::io::io_ok(product->con2matter(std::string("product.con"))));

    // Interpolate the original cluster displacement under the Pt Morse
    // potential.
    std::vector<std::shared_ptr<Matter>> seed;
    seed.push_back(std::make_shared<Matter>(*reactant));
    auto mid = std::make_shared<Matter>(*reactant);
    const AtomMatrix r0 = reactant->getPositions();
    const AtomMatrix r1 = product->getPositions();
    mid->setPositions(0.5 * (r0 + r1));
    seed.push_back(mid);
    seed.push_back(std::make_shared<Matter>(*product));

    const double e0 = seed[0]->getPotentialEnergy();
    const double e1 = seed[1]->getPotentialEnergy();
    const double e2 = seed[2]->getPotentialEnergy();
    lineb.energy = e0 + 0.25 * (e1 - e0);
    REQUIRE(e0 < lineb.energy);
    REQUIRE(e1 > lineb.energy);
    REQUIRE(e2 > e1);
    std::vector<eonc::io::ConFrameMetadata> metas(seed.size());
    for (size_t i = 0; i < seed.size(); ++i) {
      metas[i].energy = seed[i]->getPotentialEnergy();
    }
    REQUIRE(
        eonc::io::io_ok(eonc::io::writeNebPath("lineb_seed.con", seed, metas)));

    lineb.initial_path = "lineb_seed.con";
    eonc::LinebJob job(pot, params);
    REQUIRE_THROWS_WITH(
        job.run(),
        Catch::Matchers::ContainsSubstring("product-side crossing") ||
            Catch::Matchers::ContainsSubstring("STALLED"));
  }

  SECTION("both sides of the Morse barrier reach the target energy") {
    const char *systems = std::getenv("EON_TEST_SYSTEMS_DIR");
    REQUIRE(systems != nullptr);
    const auto example = std::filesystem::path(systems) /
                         "../../../../examples/lineb-morse/seed.con";
    lineb.initial_path = example.string();
    auto seed = eonc::lineb::loadSeedBand(params, pot, lineb.initial_path);
    REQUIRE(seed.size() == 3);
    const double e0 = seed.front()->getPotentialEnergy();
    const double e1 = seed[1]->getPotentialEnergy();
    lineb.energy = e0 + 0.25 * (e1 - e0);
    REQUIRE(e0 < lineb.energy);
    REQUIRE(e1 > lineb.energy);
    REQUIRE(seed.back()->getPotentialEnergy() < lineb.energy);
    eonc::LinebJob job(pot, params);
    const auto files = job.run();
    REQUIRE_FALSE(files.empty());

    std::ifstream results("results.dat");
    REQUIRE(results.good());
    std::string blob((std::istreambuf_iterator<char>(results)),
                     std::istreambuf_iterator<char>());
    REQUIRE(blob.find("action") != std::string::npos);
    REQUIRE(blob.find("log_rate_over_reactant") != std::string::npos);
    REQUIRE(std::filesystem::exists("lineb.con"));
    const auto frame = readcon::read_first_frame("lineb.con");
    const auto metadata = nlohmann::json::parse(frame.metadata_json());
    REQUIRE(metadata.at("lineb").at("spectrum").is_null());
    REQUIRE(metadata.at("lineb").at("replica_clock") == "quadrature");
#ifdef RGSADDLE_HAS_LINEB_RATE_DIAGNOSTICS
    lineb.beads = 8;
    eonc::LinebJob rateJob(pot, params);
    const auto callsBeforeRateJob = pot->forceCallCounter.load();
    const auto rateFiles = rateJob.run();
    REQUIRE_FALSE(rateFiles.empty());
    const auto rateFrame = readcon::read_first_frame("lineb.con");
    const auto rateMetadata = nlohmann::json::parse(rateFrame.metadata_json());
    const auto &rate = rateMetadata.at("lineb");
    REQUIRE(rate.at("potential_calls").at("total").get<size_t>() ==
            pot->forceCallCounter.load() - callsBeforeRateJob);
    // One mobile atom contributes three coordinates per replica. The
    // imaginary-time translation removes one mode from the full ring.
    REQUIRE(rate.at("spectrum").at("modes") == 3 * lineb.beads - 1);
    REQUIRE(rate.at("spectrum").at("negative_modes").is_number_integer());
    REQUIRE(rate.at("log_rate_over_reactant").is_number());
#ifdef RGSADDLE_HAS_LINEB_ORBIT
    REQUIRE(std::find(rateFiles.begin(), rateFiles.end(), "lineb-orbit.con") !=
            rateFiles.end());
    const auto frames = readcon::read_all_frames("lineb-orbit.con");
    REQUIRE(frames.size() == static_cast<size_t>(lineb.beads));
    std::vector<std::array<double, 3>> coordinates, gradients;
    double mass = 0;
    for (const auto &bead : frames) {
      const auto info = nlohmann::json::parse(bead.metadata_json());
      REQUIRE(info.at("temperature_K") == rate.at("temperature"));
      REQUIRE(info.at("tau_over_beta").get<double>() ==
              Catch::Approx(static_cast<double>(coordinates.size()) /
                            frames.size()));
      REQUIRE(info.at("evaluation_counts") ==
              rateMetadata.at("evaluation_counts"));
      const auto &observed = info.at("lineb_orbit");
      REQUIRE(observed.at("mobile_atom_ids").size() == 1);
      const auto id = observed.at("mobile_atom_ids").at(0).get<uint64_t>();
      const auto atoms = bead.atoms();
      const auto atom =
          std::find_if(atoms.begin(), atoms.end(),
                       [id](const auto &a) { return a.atom_id == id; });
      REQUIRE(atom != atoms.end());
      coordinates.push_back({atom->x, atom->y, atom->z});
      gradients.push_back(observed.at("gradients_eV_per_angstrom")
                              .at(0)
                              .get<std::array<double, 3>>());
      mass = atom->mass;
    }
    const auto n = frames.size();
    const double spring =
        mass * std::pow(n / (2 * rate.at("half_period").get<double>()), 2);
    double maxAbs = 0, sumSquares = 0;
    for (size_t b = 0; b < n; ++b) {
      for (size_t c = 0; c < 3; ++c) {
        const double residual =
            gradients[b][c] +
            spring * ((coordinates[b][c] - coordinates[(b + n - 1) % n][c]) +
                      (coordinates[b][c] - coordinates[(b + 1) % n][c]));
        maxAbs = std::max(maxAbs, std::abs(residual));
        sumSquares += residual * residual;
      }
    }
    const auto &summary = rate.at("stationarity");
    REQUIRE(summary.at("max_abs_eV_per_angstrom").get<double>() ==
            Catch::Approx(maxAbs).margin(1e-10));
    REQUIRE(summary.at("rms_eV_per_angstrom").get<double>() ==
            Catch::Approx(std::sqrt(sumSquares / (3 * n))).margin(1e-10));
#ifdef WITH_RGSADDLE_LINEB_TURNING
    lineb.interpolated_force_clock = true;
    eonc::LinebJob clockJob(pot, params);
    clockJob.run();
    const auto clockFrame = readcon::read_first_frame("lineb-orbit.con");
    const auto clockMetadata =
        nlohmann::json::parse(clockFrame.metadata_json());
    const auto &clock = clockMetadata.at("lineb");
    REQUIRE(clock.at("replica_clock") == "interpolated_force_verlet");
    REQUIRE(clock.at("surface_calls").at("dynamics").get<size_t>() ==
            seed.size() + 1);
    const auto &turn = clock.at("turning");
    const auto launch = turn.at("launch_image").get<size_t>();
    REQUIRE((launch == 0 || launch == seed.size() - 1));
    REQUIRE(turn.at("return_energy_minus_target_eV").is_number());
    REQUIRE(turn.at("max_cartesian_displacement_angstrom").is_number());
#ifdef WITH_RGSADDLE_LINEB_RATE
    lineb.interpolated_hessian = true;
    lineb.image_hessian_step = 0.005;
    eonc::LinebJob imageRateJob(pot, params);
    imageRateJob.run();
    const auto imageFrame = readcon::read_first_frame("lineb-orbit.con");
    const auto imageMetadata =
        nlohmann::json::parse(imageFrame.metadata_json());
    const auto &imageRate = imageMetadata.at("lineb");
    REQUIRE(imageRate.at("hessian") == "interpolated_images");
    REQUIRE(imageRate.at("image_hessian_step_angstrom").get<double>() ==
            Catch::Approx(lineb.image_hessian_step));
    const auto &cost = imageRate.at("rate_calls");
    REQUIRE(cost.at("curvature").get<size_t>() == 2 * 3 * seed.size());
    REQUIRE(cost.at("bead").get<long>() == lineb.beads);
    REQUIRE(imageRate.at("surface_calls").at("rate").get<size_t>() ==
            cost.at("curvature").get<size_t>() + cost.at("bead").get<size_t>());
#endif
#endif
#endif
#endif
  }
}

} // namespace tests
