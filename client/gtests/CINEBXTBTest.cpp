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

#include "NudgedElasticBand.h"
#include "catch2/catch_amalgamated.hpp"
#include <spdlog/sinks/null_sink.h>
#include <spdlog/spdlog.h>

namespace tests {

// Set up a null logger so NEB internals don't crash on spdlog::get("combi").
struct LoggerSetup {
  LoggerSetup() {
    if (!spdlog::get("combi")) {
      auto sink = std::make_shared<spdlog::sinks::null_sink_mt>();
      auto logger = std::make_shared<spdlog::logger>("combi", sink);
      spdlog::register_logger(logger);
      spdlog::set_default_logger(logger);
    }
  }
};
static LoggerSetup _logger_setup;

// Regression test: CI-NEB with XTB on a small (9-atom) molecule.
// Reproduces a bug where removing EIGEN_DEFAULT_TO_ROW_MAJOR silently changed
// the storage order of bare MatrixXd types, corrupting force projections and
// causing the NEB to diverge from the first step (issue introduced in 6e8461c3).
TEST_CASE("CI-NEB XTB regression", "[neb][xtb]") {
  Parameters params;
  params.potential_options.potential = PotType::XTB;
  params.xtb_options.paramset = "GFN2xTB";
  params.xtb_options.acc = 1.0;
  params.xtb_options.elec_temperature = 300.0;
  params.xtb_options.maxiter = 250;

  params.neb_options.image_count = 10;
  params.neb_options.spring.weighting.enabled = true;
  params.neb_options.spring.weighting.k_min = 0.972;
  params.neb_options.spring.weighting.k_max = 9.72;
  params.neb_options.spring.weighting.trigger = 0.5;
  params.neb_options.initialization.method = NEBInit::LINEAR;
  params.neb_options.endpoints.minimize = false;
  params.neb_options.climbing_image.enabled = true;
  params.neb_options.climbing_image.converged_only = true;
  params.neb_options.climbing_image.trigger_force = 0.5;
  params.neb_options.climbing_image.trigger_factor = 0.8;
  params.neb_options.force_tolerance = 0.0514221;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.max_iterations = 100;
  params.optimizer_options.max_move = 0.1;

  auto pot =
      helper_functions::makePotential(params.potential_options.potential, params);
  auto initial = std::make_shared<Matter>(pot, params);
  auto final_state = std::make_shared<Matter>(pot, params);

  std::string reactFile("reactant.con");
  std::string prodFile("product.con");
  initial->con2matter(reactFile);
  final_state->con2matter(prodFile);

  auto neb =
      std::make_unique<NudgedElasticBand>(initial, final_state, params, pot);
  auto status = neb->compute();

  REQUIRE(static_cast<int>(status) ==
          static_cast<int>(NudgedElasticBand::NEBStatus::GOOD));

  neb->findExtrema();
  REQUIRE(neb->numExtrema >= 1);
}

} /* namespace tests */
