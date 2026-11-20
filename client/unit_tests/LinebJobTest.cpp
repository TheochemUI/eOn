/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/LinebJob.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/Tunneling.h"

#ifdef WITH_RGSADDLE_LINEB
#include <rgsaddle.h>
#endif

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
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
  REQUIRE_THROWS_WITH(job.run(),
                      Catch::Matchers::ContainsSubstring("rgsaddle_lineb_solve"));
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
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(reactant->con2matter(std::string("reactant.con"))));
  REQUIRE(eonc::io::io_ok(product->con2matter(std::string("product.con"))));

  // Three-image linear seed between the minima, energy below the mid-point.
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
  const double barrier = std::max(e1 - e0, e1 - e2);
  REQUIRE(barrier > 0.0);
  const double energy = e0 + 0.25 * (e1 - e0);

  std::vector<eonc::io::ConFrameMetadata> metas(seed.size());
  for (size_t i = 0; i < seed.size(); ++i) {
    metas[i].energy = seed[i]->getPotentialEnergy();
  }
  REQUIRE(eonc::io::io_ok(eonc::io::writeNebPath("lineb_seed.con", seed, metas)));

  auto &lineb = ParametersLoadAccess::lineb_options(params);
  lineb.initial_path = "lineb_seed.con";
  lineb.energy = energy;
  lineb.opt_tol = 1e-2;
  lineb.perp_tol = 1e-2;
  lineb.max_steps = 200;
  lineb.beads = 0;
  lineb.molecular = false;
  lineb.hbar = eonc::tunneling::kHbar;
  lineb.boltzmann = eonc::tunneling::kBoltzmann;

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
}

} // namespace tests
