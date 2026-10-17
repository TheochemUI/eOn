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
#include "eon/Dimer.h"
#include "eon/ImprovedDimer.h"
#include "eon/LORRotation.h"
#include "eon/LowestEigenmode.h"
#include "eon/Matter.h"
#include "eon/NudgedElasticBand.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/ProcessSearchJob.h"
#include "eon/RandomNumbers.h"
#include "eon/Runtime.h"
#include "eon/SaddleSearchJob.h"

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace tests {

static eonc::helpers::test::QuillTestLogger quillSetup;

// LJ wrapper. force() counts on an atomic the registry cannot see.
// needsPerImageInstance makes dimer, NEB, and process search build a
// second instance when Main.parallel is left at its default.
class CountingLj final : public eonc::Potential {
public:
  CountingLj(std::unique_ptr<eonc::Potential> inner,
             std::shared_ptr<std::atomic<size_t>> calls,
             eonc::Parameters params)
      : eonc::Potential(eonc::PotType::LJ),
        inner_(std::move(inner)),
        calls_(std::move(calls)),
        params_(std::move(params)) {}

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return true;
  }
  [[nodiscard]] std::shared_ptr<eonc::Potential> clonePotential() const override {
    auto inner = eonc::helpers::makePotential(eonc::PotType::LJ, params_);
    return std::make_shared<CountingLj>(std::move(inner), calls_, params_);
  }

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    calls_->fetch_add(1, std::memory_order_relaxed);
    inner_->force(nAtoms, positions, atomicNrs, forces, energy, variance, box);
  }

private:
  std::unique_ptr<eonc::Potential> inner_;
  std::shared_ptr<std::atomic<size_t>> calls_;
  eonc::Parameters params_;
};

struct Counted {
  eonc::Runtime rt;
  std::shared_ptr<std::atomic<size_t>> calls;
  std::shared_ptr<CountingLj> pot;
  eonc::Parameters params;
};

struct CwdGuard {
  std::filesystem::path previous;
  std::filesystem::path dir;
  bool active{false};

  explicit CwdGuard(std::string name) {
    previous = std::filesystem::current_path();
    dir = std::filesystem::temp_directory_path() / std::move(name);
    std::filesystem::create_directories(dir);
    std::filesystem::current_path(dir);
    active = true;
  }
  ~CwdGuard() {
    std::error_code ec;
    if (active) {
      std::filesystem::current_path(previous, ec);
    }
    if (!dir.empty()) {
      std::filesystem::remove_all(dir, ec);
    }
  }
  CwdGuard(const CwdGuard &) = delete;
  CwdGuard &operator=(const CwdGuard &) = delete;
};

eonc::Parameters makeParams(bool parallelOff) {
  eonc::Parameters params;
  auto &main = eonc::ParametersLoadAccess::main_options(params);
  main.randomSeed = 1;
  if (parallelOff) {
    main.parallel = false;
  }
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::LJ;
  auto &opt = eonc::ParametersLoadAccess::optimizer_options(params);
  opt.max_iterations = 1;
  opt.converged_force = 1.0e6;
  auto &saddle = eonc::ParametersLoadAccess::saddle_search_options(params);
  saddle.method = "min_mode";
  saddle.minmode_method = eonc::LowestEigenmode::MINMODE_DIMER;
  saddle.max_iterations = 1;
  saddle.converged_force = 0.0;
  saddle.displace_type = "last_atom";
  saddle.displace_magnitude = 0.1;
  auto &dimer = eonc::ParametersLoadAccess::dimer_options(params);
  dimer.improved = true;
  dimer.rotations_max = 1;
  dimer.rotations_min = 1;
  dimer.max_iterations = 1;
  eonc::ParametersLoadAccess::prefactor_options(params).default_value = 1.0;
  eonc::ParametersLoadAccess::process_search_options(params).minimize_first =
      false;
  auto &neb = eonc::ParametersLoadAccess::neb_options(params);
  neb.image_count = 3;
  neb.max_iterations = 1;
  neb.climbing_image.enabled = false;
  neb.endpoints.minimize = false;
  neb.initialization.method = eonc::NEBInit::LINEAR;
  return params;
}

Counted beginRun(bool parallelOff) {
  eonc::PotRegistry::get().reset();
  eonc::rng::random(1);
  Counted run;
  run.calls = std::make_shared<std::atomic<size_t>>(0);
  run.params = makeParams(parallelOff);
  {
    eonc::PotentialConstructionScope scope(run.rt.pots());
    auto inner =
        eonc::helpers::makePotential(eonc::PotType::LJ, run.params);
    run.pot =
        std::make_shared<CountingLj>(std::move(inner), run.calls, run.params);
  }
  return run;
}

std::shared_ptr<eonc::Matter> twoAtom(const std::shared_ptr<eonc::Potential> &pot,
                                      const eonc::Parameters &params,
                                      double shift) {
  auto matter = std::make_shared<eonc::Matter>(pot, params);
  matter->resize(2);
  matter->setAtomicNr(0, 1);
  matter->setAtomicNr(1, 1);
  matter->setMass(0, 1.0);
  matter->setMass(1, 1.0);
  matter->setCell(Matrix3d::Identity() * 20.0);
  matter->setPeriodic(false);
  AtomMatrix pos = AtomMatrix::Zero(2, 3);
  pos(1, 0) = 1.5 + shift;
  matter->setPositions(pos);
  return matter;
}

AtomMatrix axisMode() {
  AtomMatrix mode = AtomMatrix::Zero(2, 3);
  mode(0, 0) = 1.0;
  return mode;
}

std::string workName(const char *label) {
  static int seq = 0;
  return std::string("eon-force-count-") + label + "-" + std::to_string(++seq) +
         "-" +
         std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count());
}

void expectRegistry(const Counted &run, bool parallel) {
  const auto reported = run.rt.pots().total_force_calls();
  const auto independent = run.calls->load();
  REQUIRE(independent > 0);
  REQUIRE(reported == independent);
  REQUIRE(run.pot->registryForceCalls() == independent);
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == 0);
  if (parallel) {
    REQUIRE(run.pot->forceCallCounter.load() < independent);
  } else {
    REQUIRE(run.pot->forceCallCounter.load() == independent);
  }
}

std::map<std::string, std::string> readResults(const std::filesystem::path &path) {
  std::map<std::string, std::string> out;
  std::ifstream in(path);
  REQUIRE(in.good());
  std::string line;
  while (std::getline(in, line)) {
    const auto pos = line.find(' ');
    if (pos == std::string::npos) {
      continue;
    }
    out.emplace(line.substr(pos + 1), line.substr(0, pos));
  }
  return out;
}

void runImproved(bool parallelOff) {
  CwdGuard cwd(workName("idimer"));
  auto run = beginRun(parallelOff);
  auto matter = twoAtom(run.pot, run.params, 0.0);
  eonc::ImprovedDimer dimer(matter, run.params, run.pot);
  dimer.compute(matter, axisMode());
  expectRegistry(run, !parallelOff);
}

void runClassical(bool parallelOff) {
  CwdGuard cwd(workName("dimer"));
  auto run = beginRun(parallelOff);
  eonc::ParametersLoadAccess::dimer_options(run.params).improved = false;
  eonc::ParametersLoadAccess::dimer_options(run.params).rotation_backend =
      eonc::DimerRotationBackend::Classical;
  auto matter = twoAtom(run.pot, run.params, 0.0);
  eonc::Dimer dimer(matter, run.params, run.pot);
  dimer.compute(matter, axisMode());
  expectRegistry(run, !parallelOff);
}

void runLor(bool parallelOff) {
  CwdGuard cwd(workName("lor"));
  auto run = beginRun(parallelOff);
  eonc::ParametersLoadAccess::dimer_options(run.params).rotation_backend =
      eonc::DimerRotationBackend::LOR;
  auto matter = twoAtom(run.pot, run.params, 0.0);
  eonc::LORRotation lor(matter, run.params, run.pot);
  lor.compute(matter, axisMode());
  expectRegistry(run, !parallelOff);
}

void runNeb(bool parallelOff) {
  CwdGuard cwd(workName("neb"));
  auto run = beginRun(parallelOff);
  auto reactant = twoAtom(run.pot, run.params, 0.0);
  auto product = twoAtom(run.pot, run.params, 0.4);
  eonc::NudgedElasticBand neb(reactant, product, run.params, run.pot);
  static_cast<void>(neb.compute());
  expectRegistry(run, !parallelOff);
}

void runSaddle(bool parallelOff) {
  CwdGuard cwd(workName("saddle"));
  auto run = beginRun(parallelOff);
  auto seed = twoAtom(run.pot, run.params, 0.0);
  eonc::SaddleSearchJob job(run.pot, run.params);
  static_cast<void>(job.runFromMatter(seed));
  expectRegistry(run, !parallelOff);
  const auto results = readResults("results.dat");
  REQUIRE(results.count("total_force_calls") == 1);
  const auto written = std::stoull(results.at("total_force_calls"));
  const auto independent = run.calls->load();
  REQUIRE(written > 0);
  REQUIRE(independent >= written);
  REQUIRE(independent - written <= 2);
}

void runProcess(bool parallelOff) {
  CwdGuard cwd(workName("process"));
  auto run = beginRun(parallelOff);
  auto seed = twoAtom(run.pot, run.params, 0.0);
  eonc::ProcessSearchJob job(run.pot, run.params);
  static_cast<void>(job.runFromMatter(seed));
  expectRegistry(run, !parallelOff);
  const auto results = readResults("results.dat");
  REQUIRE(results.count("total_force_calls") == 1);
  const auto written = std::stoull(results.at("total_force_calls"));
  const auto summed = std::stoull(results.at("force_calls_minimization")) +
                      std::stoull(results.at("force_calls_saddle")) +
                      std::stoull(results.at("force_calls_prefactors"));
  REQUIRE(written == summed);
  REQUIRE(run.rt.pots().total_force_calls() >= written);
}

} // namespace tests

TEST_CASE("Main.parallel default counts the second potential",
          "[force][parallel]") {
  eonc::Parameters defaults;
  REQUIRE(defaults.main_options().parallel);
  SECTION("improved dimer") { tests::runImproved(false); }
  SECTION("classical dimer") { tests::runClassical(false); }
  SECTION("LOR") { tests::runLor(false); }
  SECTION("NEB") { tests::runNeb(false); }
  SECTION("saddle search") { tests::runSaddle(false); }
  SECTION("process search") { tests::runProcess(false); }
}

TEST_CASE("parallel false counts the shared potential", "[force][parallel]") {
  SECTION("improved dimer") { tests::runImproved(true); }
  SECTION("classical dimer") { tests::runClassical(true); }
  SECTION("LOR") { tests::runLor(true); }
  SECTION("NEB") { tests::runNeb(true); }
  SECTION("saddle search") { tests::runSaddle(true); }
  SECTION("process search") { tests::runProcess(true); }
}
