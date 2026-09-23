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
#include "eon/BasinHoppingSaddleSearch.h"
#include "eon/ConFileIO.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"
#include "eon/RandomNumbers.h"
#include "eon/fpe_handler.h"
#include <cfenv>
#include <cmath>
#include <filesystem>
#include <system_error>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class SaddleSearchFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;

  SaddleSearchFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::optimizer_options(params).method = OptType::CG;
    ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.1;
    ParametersLoadAccess::dimer_options(params).improved = true;
    ParametersLoadAccess::dimer_options(params).converged_angle = 0.01;
    ParametersLoadAccess::dimer_options(params).max_iterations = 50;
    ParametersLoadAccess::saddle_search_options(params).minmode_method =
        LowestEigenmode::MINMODE_DIMER;
    ParametersLoadAccess::saddle_search_options(params).max_iterations = 100;
    ParametersLoadAccess::saddle_search_options(params).converged_force = 0.05;
    ParametersLoadAccess::saddle_search_options(params).max_energy = 20.0;

    pot = eonc::helpers::makePotential(PotType::LJ, params);
    matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("reactant.con"));

    // Displace one atom to break symmetry
    auto pos = matter->getPositions();
    pos(0, 0) += 0.3;
    pos(0, 1) -= 0.2;
    matter->setPositions(pos);
  }
};

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch runs without crashing",
                 "[saddle_search]") {
  // Create a random initial mode
  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  int status = search.run();

  // Status should be a valid enum value (0 through 21)
  REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status <= MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST);
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch reports finite eigenvalue after run",
                 "[saddle_search]") {
  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  search.run();

  double eigenvalue = search.getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch hits max iterations with low limit",
                 "[saddle_search]") {
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-10;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  int status = search.run();

  // Should hit max iterations or some non-GOOD status
  REQUIRE(status != MinModeSaddleSearch::STATUS_GOOD);
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "write_movies writes per-iteration mode files",
                 "[saddle_search][ra6]") {
  namespace fs = std::filesystem;
  ParametersLoadAccess::debug_options(params).write_movies = true;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-10;

  const auto tmp = fs::temp_directory_path() / "eon_ra6_modes";
  fs::create_directories(tmp);
  const auto old = fs::current_path();
  fs::current_path(tmp);

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();
  MinModeSaddleSearch search(matter, mode, matter->getPotentialEnergy(), params,
                             pot);
  search.run();

  REQUIRE(fs::exists("mode_000.dat"));
  fs::current_path(old);
  fs::remove_all(tmp);
}

// Issue #20: unfeasible / unconverged climb must not report STATUS_GOOD.
// finalizeClimbStatus is the shipped policy used by MinModeSaddleSearch::run.
// If the guard body is deleted (always returns climbStatus), these fail.
TEST_CASE("finalizeClimbStatus refuses STATUS_GOOD when unconverged (#20)",
          "[saddle_search][issue_20]") {
  using S = MinModeSaddleSearch;
  REQUIRE(
      S::finalizeClimbStatus(S::STATUS_GOOD, /*objectiveConverged=*/false) ==
      S::STATUS_BAD_MAX_ITERATIONS);
  REQUIRE(S::finalizeClimbStatus(S::STATUS_GOOD, /*objectiveConverged=*/true) ==
          S::STATUS_GOOD);
  // Non-GOOD statuses are left alone regardless of convergence flag.
  REQUIRE(S::finalizeClimbStatus(S::STATUS_BAD_MAX_ITERATIONS, false) ==
          S::STATUS_BAD_MAX_ITERATIONS);
  REQUIRE(S::finalizeClimbStatus(S::STATUS_BAD_HIGH_ENERGY, false) ==
          S::STATUS_BAD_HIGH_ENERGY);
  REQUIRE(S::finalizeClimbStatus(S::STATUS_DIMER_RESTORED_BEST, false) ==
          S::STATUS_DIMER_RESTORED_BEST);
  // Removing the unconverged branch of finalizeClimbStatus makes the first
  // REQUIRE fail (would return STATUS_GOOD).
}

TEST_CASE_METHOD(
    SaddleSearchFixture,
    "MinModeSaddleSearch run never returns STATUS_GOOD if unconverged (#20)",
    "[saddle_search][issue_20]") {
  // Force a non-converged climb: tiny iteration budget + absurd force target.
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-20;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-20;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);
  int status = search.run(1);

  // Real entry path: with an unconverged climb objective the search must not
  // report success (covers finalizeClimbStatus applied inside run()).
  REQUIRE(status != MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status != MinModeSaddleSearch::STATUS_INIT);
}

// Basin hopping used to ignore MinModeSaddleSearch::run and return success.
// One climb step cannot meet the force tolerance, so the hop must report
// that failure on both the return value and getStatus().
TEST_CASE_METHOD(SaddleSearchFixture,
                 "BasinHoppingSaddleSearch reports a failed climb",
                 "[saddle_search][basin_hopping]") {
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-20;

  namespace fs = std::filesystem;
  const auto tmp = fs::temp_directory_path() / "eon_bh_climb_status";
  fs::create_directories(tmp);
  struct CwdGuard {
    fs::path old;
    explicit CwdGuard(const fs::path &next)
        : old(fs::current_path()) {
      fs::current_path(next);
    }
    ~CwdGuard() { fs::current_path(old); }
  };

  int status = MinModeSaddleSearch::STATUS_INIT;
  {
    CwdGuard guard(tmp);
    auto hop = std::make_shared<Matter>(*matter);
    BasinHoppingSaddleSearch search(matter, hop, pot, params);
    status = search.run();
    REQUIRE(search.getStatus() == status);
  }
  fs::remove_all(tmp);

  REQUIRE(status != MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status != MinModeSaddleSearch::STATUS_INIT);
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "Basin hop rejection is not labeled Initialized",
                 "[saddle_search][basin_hopping]") {
  BasinHoppingSaddleSearch search(matter, matter, pot, params);
  // Shared MinMode table still calls code 1 Initialized. This search does not.
  REQUIRE(std::string(MinModeSaddleSearch::statusMessage(1)) == "Initialized");
  REQUIRE(std::string(search.describeStatus(1)) == "Basin hop rejected");
  REQUIRE(std::string(search.describeStatus(0)) == "Success");
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch forces on fixed atoms remain zero",
                 "[saddle_search][fixed_atoms]") {
  // Fix all atoms except the first
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    matter->setFixed(i, true);
  }

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);

  search.run();

  // After saddle search, forces on fixed atoms must still be zero
  const AtomMatrix &forces = matter->getForces();
  for (int i = 1; i < matter->numberOfAtoms(); i++) {
    REQUIRE(forces.row(i).norm() < 1e-10);
  }
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "MinModeSaddleSearch with Lanczos eigenmode",
                 "[saddle_search][lanczos]") {
  ParametersLoadAccess::saddle_search_options(params).minmode_method =
      LowestEigenmode::MINMODE_LANCZOS;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 50;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);
  int status = search.run();

  REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(status <= MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST);
  REQUIRE(std::isfinite(search.getEigenvalue()));
}

TEST_CASE("basin hop Metropolis divides only after the uphill test",
          "[saddle_search][basin_hopping]") {
  const double kB = 8.6173324e-5;
  const double de = 0.05;
  const double temperature = 300.0;
  using Search = eonc::BasinHoppingSaddleSearch;
  REQUIRE(Search::metropolisProbability(-1.0, kB, 0.0) == 1.0);
  REQUIRE(Search::metropolisProbability(0.0, kB, 0.0) == 1.0);
  REQUIRE(Search::metropolisProbability(1.0, kB, 0.0) == 0.0);
  REQUIRE(Search::metropolisProbability(1.0, kB, -25.0) == 0.0);
  REQUIRE(Search::metropolisProbability(1.0, 0.0, 300.0) == 0.0);
  REQUIRE(Search::metropolisProbability(de, kB, temperature) ==
          Catch::Approx(std::exp(-de / (kB * temperature))));

#if !defined(_WIN32) && !(defined(__APPLE__) && defined(__aarch64__))
  eonc::enableFPE();
  REQUIRE(Search::metropolisProbability(1.0, kB, 0.0) == 0.0);
  REQUIRE(Search::metropolisProbability(-2.0, kB, 0.0) == 1.0);
  REQUIRE(Search::metropolisProbability(1.0, kB, -10.0) == 0.0);
  eonc::disableFPE();
  feclearexcept(FE_ALL_EXCEPT);
#endif
}

TEST_CASE_METHOD(
    SaddleSearchFixture,
    "BasinHoppingSaddleSearch rejects uphill hops at non-positive temperature",
    "[saddle_search][basin_hopping]") {
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  ParametersLoadAccess::neb_options(params).max_iterations = 1;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;

  auto displaced = std::make_shared<Matter>(pot, params);
  *displaced = *matter;
  auto pos = displaced->getPositions();
  // Park atom 0 on atom 1. The contact is uphill for LJ, and max_iterations
  // is 0 so relax does not walk it back downhill before the accept test.
  pos.row(0) = pos.row(1);
  pos(0, 0) += 0.05;
  displaced->setPositions(pos);
  REQUIRE(displaced->getPotentialEnergy() > matter->getPotentialEnergy());

  struct RestoreCwd {
    std::filesystem::path old;
    std::filesystem::path tmp;
    explicit RestoreCwd(std::filesystem::path next)
        : old(std::filesystem::current_path()),
          tmp(std::move(next)) {
      std::filesystem::create_directories(tmp);
      std::filesystem::current_path(tmp);
    }
    ~RestoreCwd() {
      std::error_code ec;
      std::filesystem::current_path(old, ec);
      std::filesystem::remove_all(tmp, ec);
    }
  } cwd(std::filesystem::temp_directory_path() / "eon_bh_metro");

  ParametersLoadAccess::main_options(params).temperature = -100.0;
  {
    eonc::BasinHoppingSaddleSearch search(matter, displaced, pot, params);
    REQUIRE(search.run() == 1);
  }

#if !defined(_WIN32) && !(defined(__APPLE__) && defined(__aarch64__))
  eonc::enableFPE();
#endif
  ParametersLoadAccess::main_options(params).temperature = 0.0;
  eonc::BasinHoppingSaddleSearch search(matter, displaced, pot, params);
  const int status = search.run();
#if !defined(_WIN32) && !(defined(__APPLE__) && defined(__aarch64__))
  eonc::disableFPE();
  feclearexcept(FE_ALL_EXCEPT);
#endif
  REQUIRE(status == 1);
}

TEST_CASE_METHOD(SaddleSearchFixture, "MinModeSaddleSearch with classic Dimer",
                 "[saddle_search][classic_dimer]") {
  ParametersLoadAccess::dimer_options(params).improved =
      false; // classic dimer, not improved
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 50;

  long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();

  double reactantEnergy = matter->getPotentialEnergy();
  MinModeSaddleSearch search(matter, mode, reactantEnergy, params, pot);
  int status = search.run();

  REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(std::isfinite(search.getEigenvalue()));
}

TEST_CASE_METHOD(SaddleSearchFixture,
                 "basin hopping keeps the last interior NEB image",
                 "[saddle_search][basin_hopping]") {
  const long nImages = 5;
  std::vector<std::shared_ptr<Matter>> path;
  path.reserve(static_cast<size_t>(nImages + 2));
  for (long i = 0; i <= nImages + 1; i++) {
    auto image = std::make_shared<Matter>(*matter);
    image->setComputedPotential(static_cast<double>(i), 0.0);
    path.push_back(std::move(image));
  }
  // Endpoints are higher. Only interiors 1..numImages are eligible, and the
  // last of those is the maximum. An exclusive upper bound would return 4.
  path[0]->setComputedPotential(100.0, 0.0);
  path[static_cast<size_t>(nImages + 1)]->setComputedPotential(100.0, 0.0);
  REQUIRE(BasinHoppingSaddleSearch::highestEnergyInteriorImage(path, nImages) ==
          static_cast<int>(nImages));

  path[1]->setComputedPotential(3.0, 0.0);
  path[2]->setComputedPotential(9.0, 0.0);
  path[3]->setComputedPotential(4.0, 0.0);
  path[4]->setComputedPotential(8.0, 0.0);
  path[5]->setComputedPotential(2.0, 0.0);
  REQUIRE(BasinHoppingSaddleSearch::highestEnergyInteriorImage(path, nImages) ==
          2);

  // image_count 1: the only interior bead, not the reactant or the product.
  path[1]->setComputedPotential(1.0, 0.0);
  path[2]->setComputedPotential(50.0, 0.0);
  REQUIRE(BasinHoppingSaddleSearch::highestEnergyInteriorImage(path, 1) == 1);
  REQUIRE(BasinHoppingSaddleSearch::highestEnergyInteriorImage(path, 0) == 0);
}

namespace {

struct CwdGuard {
  std::filesystem::path old;
  explicit CwdGuard(const std::filesystem::path &next)
      : old(std::filesystem::current_path()) {
    std::filesystem::current_path(next);
  }
  ~CwdGuard() { std::filesystem::current_path(old); }
  CwdGuard(const CwdGuard &) = delete;
  CwdGuard &operator=(const CwdGuard &) = delete;
};

} // namespace

TEST_CASE_METHOD(SaddleSearchFixture,
                 "basin hopping writes the last interior image",
                 "[saddle_search][basin_hopping]") {
  namespace fs = std::filesystem;
  auto reactant = std::make_shared<Matter>(pot, params);
  reactant->con2matter(std::string("reactant.con"));

  // Force tolerance and the default image count stay put. The iteration
  // caps only bound this band check; they are not a looser convergence target.
  ParametersLoadAccess::main_options(params).temperature = 1.0e20;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 20;
  ParametersLoadAccess::neb_options(params).max_iterations = 20;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 20;
  ParametersLoadAccess::dimer_options(params).max_iterations = 20;
  eonc::rng::random(42);

  const auto tmp = fs::temp_directory_path() / "eon_bh_last_image";
  fs::remove_all(tmp);
  fs::create_directories(tmp);

  auto runBand = [&](long imageCount) {
    ParametersLoadAccess::neb_options(params).image_count = imageCount;
    auto displaced = std::make_shared<Matter>(*matter);
    CwdGuard guard(tmp);
    std::error_code ec;
    fs::remove("neb_initial_band.con", ec);
    BasinHoppingSaddleSearch search(reactant, displaced, pot, params);
    const int status = search.run();
    const auto frames = readcon::read_all_frames("neb_initial_band.con");
    // Reactant plus every interior bead, including path[numImages].
    REQUIRE(frames.size() == static_cast<size_t>(imageCount) + 1);
    REQUIRE(status >= MinModeSaddleSearch::STATUS_GOOD);
    REQUIRE(status <= MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST);
    REQUIRE(std::isfinite(search.getEigenvalue()));
  };

  runBand(5);
  runBand(1);
}

TEST_CASE("basin hopping dimer direction uses the minimum image",
          "[saddle_search][basin_hopping][pbc]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  Matter image(pot, params);
  image.resize(2);
  image.setAtomicNr(0, 1);
  image.setAtomicNr(1, 1);
  image.setCell(Matrix3d::Identity() * 10.0);
  image.setPeriodic(true);

  AtomMatrix prev(2, 3);
  AtomMatrix next(2, 3);
  prev.setZero();
  next.setZero();
  // Atom 0 has crossed the cell: 0.2 and 9.8. The short step is -0.4.
  prev(0, 0) = 0.2;
  prev(0, 1) = 1.0;
  prev(0, 2) = 2.0;
  next(0, 0) = 9.8;
  next(0, 1) = 1.0;
  next(0, 2) = 2.0;
  // Atom 1 stays inside the cell. The step must not be rewritten.
  prev(1, 0) = 1.0;
  next(1, 0) = 1.5;

  AtomMatrix direction =
      BasinHoppingSaddleSearch::initialDimerDirection(image, prev, next);
  REQUIRE(direction(0, 0) == Catch::Approx(-0.2).margin(1e-12));
  REQUIRE(direction(0, 1) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(direction(0, 2) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(direction(1, 0) == Catch::Approx(0.25).margin(1e-12));
  REQUIRE(direction(1, 1) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(direction(1, 2) == Catch::Approx(0.0).margin(1e-12));

  image.setPeriodic(false);
  AtomMatrix raw =
      BasinHoppingSaddleSearch::initialDimerDirection(image, prev, next);
  REQUIRE(raw(0, 0) == Catch::Approx(4.8).margin(1e-12));
  REQUIRE(raw(1, 0) == Catch::Approx(0.25).margin(1e-12));
}

} /* namespace tests */
