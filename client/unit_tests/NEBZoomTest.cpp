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

#include "eon/NEBZoom.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"

#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>

namespace tests {
static eonc::helpers::test::QuillTestLogger _quill_setup;

using Zoom = eonc::neb_options_t::zoom_options_t;

TEST_CASE("Zoom window manual offset is centered on the climbing image",
          "[neb][zoom]") {
  Zoom cfg;
  cfg.mode = Zoom::Mode::Manual;
  cfg.offset = 2;
  const std::vector<double> energy(11, 0.0);
  const auto window = eonc::neb::zoom::selectWindow(energy, 5, cfg);
  REQUIRE(window.valid);
  REQUIRE(window.lo == 3);
  REQUIRE(window.hi == 7);
}

TEST_CASE("Zoom window auto keeps images above the barrier fraction",
          "[neb][zoom]") {
  Zoom cfg;
  cfg.mode = Zoom::Mode::Auto;
  cfg.alpha = 0.5;
  cfg.offset = 1;
  const std::vector<double> energy{0.0, 0.2, 0.4, 1.2, 2.0, 1.2, 0.4, 0.2, 0.0};
  const auto window = eonc::neb::zoom::selectWindow(energy, 4, cfg);
  REQUIRE(window.valid);
  REQUIRE(window.lo == 3);
  REQUIRE(window.hi == 5);
}

TEST_CASE(
    "Zoom window auto falls back to offset when only the climber qualifies",
    "[neb][zoom]") {
  Zoom cfg;
  cfg.mode = Zoom::Mode::Auto;
  cfg.alpha = 0.5;
  cfg.offset = 1;
  const std::vector<double> energy{0.0, 0.1, 0.2, 0.3, 5.0, 0.3, 0.2, 0.1, 0.0};
  const auto window = eonc::neb::zoom::selectWindow(energy, 4, cfg);
  REQUIRE(window.valid);
  REQUIRE(window.lo == 3);
  REQUIRE(window.hi == 5);
}

TEST_CASE("Zoom rejects a climbing image on the endpoint", "[neb][zoom]") {
  Zoom cfg;
  cfg.mode = Zoom::Mode::Manual;
  cfg.offset = 2;
  const std::vector<double> energy{0.0, 1.0, 0.0};
  const auto window = eonc::neb::zoom::selectWindow(energy, 0, cfg);
  REQUIRE_FALSE(window.valid);
}

static std::vector<std::shared_ptr<eonc::Matter>>
linePath(const std::shared_ptr<eonc::Potential> &pot, const Parameters &params,
         int count) {
  std::vector<std::shared_ptr<eonc::Matter>> path;
  path.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    auto image = std::make_shared<eonc::Matter>(pot, params);
    image->resize(1);
    image->setPeriodic(false);
    AtomMatrix pos(1, 3);
    pos.setZero();
    pos(0, 0) = static_cast<double>(i);
    image->setPositions(pos);
    path.push_back(std::move(image));
  }
  return path;
}

TEST_CASE("Zoom redistribution packs the band onto the window", "[neb][zoom]") {
  Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::LJ;
  auto pot = eonc::helpers::makePotential(eonc::PotType::LJ, params);
  auto path = linePath(pot, params, 11);

  eonc::neb::zoom::Window window;
  window.lo = 3;
  window.hi = 7;
  window.valid = true;
  REQUIRE(eonc::neb::zoom::redistributePath(path, window,
                                            Zoom::Interpolation::Cubic));

  REQUIRE(path.front()->getPositions()(0, 0) ==
          Catch::Approx(3.0).margin(1e-8));
  REQUIRE(path.back()->getPositions()(0, 0) == Catch::Approx(7.0).margin(1e-8));
  for (std::size_t i = 1; i < path.size(); ++i) {
    const double x = path[i]->getPositions()(0, 0);
    const double prev = path[i - 1]->getPositions()(0, 0);
    REQUIRE(x == Catch::Approx(prev + 0.4).margin(1e-6));
  }
}

TEST_CASE("Zoom linear and cubic agree on a straight window", "[neb][zoom]") {
  Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::LJ;
  auto pot = eonc::helpers::makePotential(eonc::PotType::LJ, params);
  auto cubic = linePath(pot, params, 11);
  auto linear = linePath(pot, params, 11);
  eonc::neb::zoom::Window window;
  window.lo = 3;
  window.hi = 7;
  window.valid = true;
  REQUIRE(eonc::neb::zoom::redistributePath(cubic, window,
                                            Zoom::Interpolation::Cubic));
  REQUIRE(eonc::neb::zoom::redistributePath(linear, window,
                                            Zoom::Interpolation::Linear));
  for (std::size_t i = 0; i < cubic.size(); ++i) {
    REQUIRE(cubic[i]->getPositions()(0, 0) ==
            Catch::Approx(linear[i]->getPositions()(0, 0)).margin(1e-6));
  }
}

TEST_CASE("Zoom ini keys load into neb options", "[neb][zoom]") {
  const auto file = std::filesystem::temp_directory_path() / "eon-zoom-neb.ini";
  {
    std::ofstream out(file);
    out << "[Nudged Elastic Band]\n"
        << "zoom_neb = true\n"
        << "zoom_alpha = 0.25\n"
        << "zoom_offset = 3\n"
        << "zoom_mode = manual\n"
        << "zoom_after = 0.2\n"
        << "zoom_interpolation = linear\n"
        << "zoom_ci_stability = 4\n"
        << "zoom_max_iterations = 40\n";
  }
  Parameters params;
  REQUIRE(params.load(file.string()) == 0);
  const auto &zoom = params.neb_options().zoom;
  REQUIRE(zoom.enabled);
  REQUIRE(zoom.alpha == Catch::Approx(0.25));
  REQUIRE(zoom.offset == 3);
  REQUIRE(zoom.mode == Zoom::Mode::Manual);
  REQUIRE(zoom.activation_threshold == Catch::Approx(0.2));
  REQUIRE(zoom.interpolation == Zoom::Interpolation::Linear);
  REQUIRE(zoom.stability_count == 4);
  REQUIRE(zoom.max_iterations == 40);
  std::filesystem::remove(file);
}

} // namespace tests
