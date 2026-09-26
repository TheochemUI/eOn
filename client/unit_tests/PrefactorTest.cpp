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
#include "eon/Matter.h"
#include "eon/Prefactor.h"

#include <filesystem>
#include <fstream>
#include <sstream>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

static std::pair<std::shared_ptr<Matter>, Parameters> makeLJCluster() {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("reactant.con"));
  return {m, params};
}

TEST_CASE("getPrefactors rejects a null Matter", "[prefactor]") {
  Parameters params;
  double pref1 = 0.0;
  double pref2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, nullptr, nullptr, nullptr,
                                         pref1, pref2) == -1);
}

TEST_CASE("cutoff filter with no displacement reports no moved atoms",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  ParametersLoadAccess::prefactor_options(params).filter_scheme =
      eonc::Prefactor::FILTER_CUTOFF;
  ParametersLoadAccess::prefactor_options(params).min_displacement = 100.0;
  double pref1 = 0.0;
  double pref2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2) == -1);
}

TEST_CASE("movedAtoms includes the displaced atom and a neighbor",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  AtomMatrix pos = saddle->getPositions();
  pos.row(1) += Eigen::RowVector3d(0.6, 0.0, 0.0);
  saddle->setPositions(pos);
  ParametersLoadAccess::prefactor_options(params).min_displacement = 0.2;
  ParametersLoadAccess::prefactor_options(params).within_radius = 3.3;
  VectorXi moved = eonc::Prefactor::movedAtoms(params, min1.get(), saddle.get(),
                                               min2.get());
  REQUIRE(moved.size() > 1);
  bool sawDisplaced = false;
  for (int i = 0; i < moved.size(); ++i) {
    if (moved[i] == 1) {
      sawDisplaced = true;
    }
  }
  REQUIRE(sawDisplaced);
}

TEST_CASE("movedAtomsPct pulls in neighbors inside the radius", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  AtomMatrix pos = saddle->getPositions();
  pos.row(2) += Eigen::RowVector3d(0.0, 0.8, 0.0);
  saddle->setPositions(pos);
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.5;
  ParametersLoadAccess::prefactor_options(params).within_radius = 4.0;
  VectorXi moved = eonc::Prefactor::movedAtomsPct(params, min1.get(),
                                                  saddle.get(), min2.get());
  REQUIRE(moved.size() > 1);
}

TEST_CASE("allFreeAtoms drops fixed rows", "[prefactor]") {
  auto [matter, params] = makeLJCluster();
  (void)params;
  matter->setFixed(0, true);
  VectorXi free = eonc::Prefactor::allFreeAtoms(matter.get());
  REQUIRE(free.size() == matter->numberOfAtoms() - 1);
  for (int i = 0; i < free.size(); ++i) {
    REQUIRE(free[i] != 0);
  }
}

TEST_CASE("logFreqs appends a wrapped frequency table", "[prefactor]") {
  namespace fs = std::filesystem;
  const auto original = fs::current_path();
  const auto dir = fs::temp_directory_path() / "eon_prefactor_freqs";
  fs::create_directories(dir);
  fs::current_path(dir);
  VectorXd freqs(6);
  freqs << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  eonc::Prefactor::logFreqs(freqs, "minimum 1");
  fs::current_path(original);
  std::ifstream in(dir / "freqs.dat");
  std::stringstream buf;
  buf << in.rdbuf();
  const std::string text = buf.str();
  fs::remove_all(dir);
  REQUIRE(text.find("minimum 1") != std::string::npos);
  REQUIRE(text.find("1.000000") != std::string::npos);
  REQUIRE(text.find("6.000000") != std::string::npos);
}

} /* namespace tests */
