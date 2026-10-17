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
#include "eon/ARTnSaddleSearch.h"
#ifdef WITH_GPRD
#include "eon/AtomicGPDimer.h"
#endif
#include "eon/BasinHoppingJob.h"
#include "eon/Dimer.h"
#include "eon/FiniteDifferenceJob.h"
#include "eon/ForceCallTimer.h"
#include "eon/HessianJob.h"
#include "eon/ImprovedDimer.h"
#include "eon/LORRotation.h"
#include "eon/LowestEigenmode.h"
#include "eon/Matter.h"
#include "eon/MonteCarloJob.h"
#include "eon/NudgedElasticBand.h"
#include "eon/ParallelReplicaJob.h"
#include "eon/Parameters.h"
#include "eon/PointJob.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/PrefactorJob.h"
#include "eon/ProcessSearchJob.h"
#include "eon/RandomNumbers.h"
#include "eon/Runtime.h"
#include "eon/SaddleSearchJob.h"
#include "eon/libs/ARTn/ARTnResource.h"
#include "relax/HostSurfacePotential.h"

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
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
             eonc::Parameters params, bool batch = false)
      : eonc::Potential(eonc::PotType::LJ),
        inner_(std::move(inner)),
        calls_(std::move(calls)),
        params_(std::move(params)),
        batch_(batch) {}

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return true;
  }
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return batch_;
  }
  [[nodiscard]] std::shared_ptr<eonc::Potential> clonePotential() const override {
    auto inner = eonc::helpers::makePotential(eonc::PotType::LJ, params_);
    return std::make_shared<CountingLj>(std::move(inner), calls_, params_,
                                        batch_);
  }

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    calls_->fetch_add(1, std::memory_order_relaxed);
    inner_->force(nAtoms, positions, atomicNrs, forces, energy, variance, box);
  }

  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override {
    if (nSystems > 0) {
      calls_->fetch_add(static_cast<size_t>(nSystems),
                        std::memory_order_relaxed);
    }
    inner_->forceBatch(nSystems, nAtoms, positions, atomicNrs, forces, energies,
                       variances, boxes);
  }

private:
  std::unique_ptr<eonc::Potential> inner_;
  std::shared_ptr<std::atomic<size_t>> calls_;
  eonc::Parameters params_;
  bool batch_;
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

Counted beginRun(bool parallelOff, bool batch = false) {
  eonc::PotRegistry::get().reset();
  eonc::rng::random(1);
  Counted run;
  run.calls = std::make_shared<std::atomic<size_t>>(0);
  run.params = makeParams(parallelOff);
  {
    eonc::PotentialConstructionScope scope(run.rt.pots());
    auto inner =
        eonc::helpers::makePotential(eonc::PotType::LJ, run.params);
    run.pot = std::make_shared<CountingLj>(std::move(inner), run.calls,
                                           run.params, batch);
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

void expectRegistry(const char *name, const Counted &run, bool parallel) {
  CAPTURE(name);
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

void expectReported(const char *name, const Counted &run,
                    std::uint64_t reported) {
  CAPTURE(name);
  const auto independent = run.calls->load();
  const auto registry = run.pot->registryForceCalls();
  REQUIRE(independent > 0);
  REQUIRE(registry == independent);
  REQUIRE(reported == independent);
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == 0);
}

void saveCon(const std::shared_ptr<eonc::Matter> &matter, const char *name) {
  REQUIRE(eonc::io::io_ok(matter->matter2con(name)));
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
  expectRegistry(parallelOff ? "improved dimer parallel=false"
                             : "improved dimer parallel=default",
                 run, !parallelOff);
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
  expectRegistry(parallelOff ? "classical dimer parallel=false"
                             : "classical dimer parallel=default",
                 run, !parallelOff);
}

void runLor(bool parallelOff) {
  CwdGuard cwd(workName("lor"));
  auto run = beginRun(parallelOff);
  eonc::ParametersLoadAccess::dimer_options(run.params).rotation_backend =
      eonc::DimerRotationBackend::LOR;
  auto matter = twoAtom(run.pot, run.params, 0.0);
  eonc::LORRotation lor(matter, run.params, run.pot);
  lor.compute(matter, axisMode());
  expectRegistry(parallelOff ? "LOR parallel=false" : "LOR parallel=default",
                 run, !parallelOff);
}

void runNeb(bool parallelOff) {
  CwdGuard cwd(workName("neb"));
  auto run = beginRun(parallelOff);
  auto reactant = twoAtom(run.pot, run.params, 0.0);
  auto product = twoAtom(run.pot, run.params, 0.4);
  eonc::NudgedElasticBand neb(reactant, product, run.params, run.pot);
  static_cast<void>(neb.compute());
  expectRegistry(parallelOff ? "NEB parallel=false" : "NEB parallel=default",
                 run, !parallelOff);
}

void runSaddle(bool parallelOff) {
  CwdGuard cwd(workName("saddle"));
  auto run = beginRun(parallelOff);
  auto seed = twoAtom(run.pot, run.params, 0.0);
  eonc::SaddleSearchJob job(run.pot, run.params);
  static_cast<void>(job.runFromMatter(seed));
  expectRegistry(parallelOff ? "saddle search parallel=false"
                             : "saddle search parallel=default",
                 run, !parallelOff);
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
  expectRegistry(parallelOff ? "process search parallel=false"
                             : "process search parallel=default",
                 run, !parallelOff);
  const auto results = readResults("results.dat");
  REQUIRE(results.count("total_force_calls") == 1);
  const auto written = std::stoull(results.at("total_force_calls"));
  const auto summed = std::stoull(results.at("force_calls_minimization")) +
                      std::stoull(results.at("force_calls_saddle")) +
                      std::stoull(results.at("force_calls_prefactors"));
  REQUIRE(written == summed);
  REQUIRE(run.rt.pots().total_force_calls() >= written);
}

void runPoint() {
  CwdGuard cwd(workName("point"));
  auto run = beginRun(false);
  saveCon(twoAtom(run.pot, run.params, 0.0), "pos.con");
  eonc::PointJob job(run.pot, run.params);
  static_cast<void>(job.run());
  const auto results = readResults("results.dat");
  expectReported("point", run, std::stoull(results.at("total_force_calls")));
}

void runPrefactor() {
  CwdGuard cwd(workName("prefactor"));
  auto run = beginRun(false, true);
  saveCon(twoAtom(run.pot, run.params, 0.0), "reactant.con");
  saveCon(twoAtom(run.pot, run.params, 0.6), "saddle.con");
  saveCon(twoAtom(run.pot, run.params, 1.2), "product.con");
  eonc::ParametersLoadAccess::hessian_options(run.params).write_modes = false;
  eonc::PrefactorJob job(run.pot, run.params);
  static_cast<void>(job.run());
  const auto results = readResults("results.dat");
  expectReported("prefactor", run,
                 std::stoull(results.at("total_force_calls")));
}

void runBasin() {
  CwdGuard cwd(workName("basin"));
  auto run = beginRun(false);
  eonc::ParametersLoadAccess::basin_hopping_options(run.params).steps = 1;
  eonc::ParametersLoadAccess::basin_hopping_options(run.params)
      .quenching_steps = 0;
  eonc::ParametersLoadAccess::basin_hopping_options(run.params)
      .adjust_displacement = false;
  eonc::ParametersLoadAccess::basin_hopping_options(run.params)
      .swap_probability = 0.0;
  saveCon(twoAtom(run.pot, run.params, 0.0), "pos.con");
  eonc::BasinHoppingJob job(run.pot, run.params);
  static_cast<void>(job.run());
  const auto results = readResults("results.dat");
  expectReported("basin hopping", run,
                 std::stoull(results.at("total_force_calls")));
}

void runFiniteDifference() {
  CwdGuard cwd(workName("fd"));
  auto run = beginRun(false);
  saveCon(twoAtom(run.pot, run.params, 0.0), "pos.con");
  eonc::FiniteDifferenceJob job(run.pot, run.params);
  static_cast<void>(job.run());
  const auto results = readResults("results.dat");
  expectReported("finite difference", run,
                 std::stoull(results.at("total_force_calls")));
}

void runHessian() {
  CwdGuard cwd(workName("hessian"));
  auto run = beginRun(false, true);
  eonc::ParametersLoadAccess::hessian_options(run.params).write_modes = false;
  saveCon(twoAtom(run.pot, run.params, 0.0), "pos.con");
  eonc::HessianJob job(run.pot, run.params);
  static_cast<void>(job.run());
  const auto results = readResults("results.dat");
  expectReported("hessian", run, std::stoull(results.at("total_force_calls")));
}

void runMonteCarlo() {
  CwdGuard cwd(workName("mc"));
  auto run = beginRun(false);
  eonc::ParametersLoadAccess::monte_carlo_options(run.params).steps = 1;
  saveCon(twoAtom(run.pot, run.params, 0.0), "pos.con");
  eonc::MonteCarloJob job(run.pot, run.params);
  static_cast<void>(job.run());
  const auto results = readResults("results.dat");
  expectReported("monte carlo", run,
                 std::stoull(results.at("total_force_calls")));
}

void runParallelReplica() {
  CwdGuard cwd(workName("preplica"));
  auto run = beginRun(false);
  auto &dyn = eonc::ParametersLoadAccess::dynamics_options(run.params);
  dyn.time_step = 0.01;
  dyn.steps = 1;
  auto &replica =
      eonc::ParametersLoadAccess::parallel_replica_options(run.params);
  replica.dephase_time = 0.01;
  replica.dephase_loop_max = 1;
  replica.refine_transition = false;
  replica.corr_time = 0.0;
  replica.state_check_interval = 0.01;
  replica.record_interval = 0.01;
  auto seed = twoAtom(run.pot, run.params, 0.0);
  eonc::ParallelReplicaJob job(run.pot, run.params);
  static_cast<void>(job.runFromMatter(seed));
  const auto results = readResults("results.dat");
  expectReported("parallel replica", run,
                 std::stoull(results.at("total_force_calls")));
}

void runBatch() {
  CwdGuard cwd(workName("batch"));
  auto run = beginRun(false);
  auto a = twoAtom(run.pot, run.params, 0.0);
  auto b = twoAtom(run.pot, run.params, 0.3);
  const auto nrsA = a->getAtomicNrs();
  const auto nrsB = b->getAtomicNrs();
  const AtomMatrix &posA = a->getPositions();
  const AtomMatrix &posB = b->getPositions();
  const Matrix3d box = Matrix3d::Zero();
  AtomMatrix forceA = AtomMatrix::Zero(2, 3);
  AtomMatrix forceB = AtomMatrix::Zero(2, 3);
  const double *positions[2] = {posA.data(), posB.data()};
  const int *atomicNrs[2] = {nrsA.data(), nrsB.data()};
  double *forces[2] = {forceA.data(), forceB.data()};
  const double *boxes[2] = {box.data(), box.data()};
  double energies[2] = {};
  double variances[2] = {};
  run.pot->forceBatch(2, 2, positions, atomicNrs, forces, energies, variances,
                      boxes);
  expectReported("batch", run, run.pot->registryForceCalls());
}

int zeroSurface(void *, eon_relax_surface_request_t *req) {
  if (req == nullptr) {
    return -1;
  }
  const long n = static_cast<long>(req->n_images);
  const long atoms = static_cast<long>(req->n_atoms);
  if (req->energies) {
    for (long s = 0; s < n; ++s) {
      req->energies[s] = 0.0;
    }
  }
  if (req->forces) {
    for (long k = 0; k < 3 * atoms * n; ++k) {
      req->forces[k] = 0.0;
    }
  }
  if (req->variances) {
    for (long s = 0; s < n; ++s) {
      req->variances[s] = 0.0;
    }
  }
  return 0;
}

void runHostBatch() {
  eonc::PotRegistry::get().reset();
  eonc::Runtime rt;
  eonc::PotentialConstructionScope scope(rt.pots());
  eonc::HostSurfacePotential pot(zeroSurface, nullptr, 0);
  double pos0[3]{0.1, 0.0, 0.0};
  double pos1[3]{0.2, 0.0, 0.0};
  double frc0[3]{};
  double frc1[3]{};
  double box0[9]{10, 0, 0, 0, 10, 0, 0, 0, 10};
  double box1[9]{10, 0, 0, 0, 10, 0, 0, 0, 10};
  int nrs0[1]{1};
  int nrs1[1]{1};
  const double *positions[2] = {pos0, pos1};
  const int *atomicNrs[2] = {nrs0, nrs1};
  double *forces[2] = {frc0, frc1};
  const double *boxes[2] = {box0, box1};
  double energies[2]{};
  double variances[2]{};
  pot.forceBatch(2, 1, positions, atomicNrs, forces, energies, variances,
                 boxes);
  REQUIRE(pot.forceCallCounter.load() == 2);
  REQUIRE(pot.registryForceCalls() == 2);
  REQUIRE(rt.pots().total_force_calls() == 2);
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == 0);
}

void runTimer() {
  CwdGuard cwd(workName("timer"));
  auto run = beginRun(false);
  auto matter = twoAtom(run.pot, run.params, 0.0);
  size_t bucket = 0;
  {
    eonc::ForceCallTimer timer(bucket);
    static_cast<void>(matter->getForces());
  }
  REQUIRE(run.calls->load() == 1);
  REQUIRE(bucket == run.calls->load());
  REQUIRE(run.rt.pots().total_force_calls() == run.calls->load());
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == 0);
}

int artnCreate() { return 0; }
void artnDestroy() {}
int artnSetParam(const char *, int, const int *, const void *) { return 0; }
void artnSetup(const int, bool *cerr) {
  if (cerr) {
    *cerr = false;
  }
}
int artnOut(void **cval, int stored) {
  if (cval == nullptr) {
    return 1;
  }
  auto *value = static_cast<int *>(std::malloc(sizeof(int)));
  if (value == nullptr) {
    return 1;
  }
  *value = stored;
  *cval = value;
  return 0;
}
int artnGetParam(const char *, void **cval) { return artnOut(cval, 0); }
int artnGetRun(const char *, void **cval) { return artnOut(cval, 0); }
int gArtnSteps = 0;
void artnStep(const int, const double, double *const, int const *,
              double *const, const double *, const int *, double *,
              bool *lconv) {
  gArtnSteps += 1;
  if (lconv) {
    *lconv = gArtnSteps >= 1;
  }
}
int artnGetData(const char *name, void **cval) {
  if (name == nullptr || cval == nullptr) {
    return 1;
  }
  if (std::strcmp(name, "has_error") == 0) {
    auto *flag = static_cast<bool *>(std::malloc(sizeof(bool)));
    if (!flag) {
      return 1;
    }
    *flag = false;
    *cval = flag;
    return 0;
  }
  if (std::strcmp(name, "has_sad") == 0) {
    auto *flag = static_cast<bool *>(std::malloc(sizeof(bool)));
    if (!flag) {
      return 1;
    }
    *flag = true;
    *cval = flag;
    return 0;
  }
  if (std::strcmp(name, "tau_sad") == 0 || std::strcmp(name, "eigen_sad") == 0) {
    const std::size_t n = 6;
    auto *buf = static_cast<double *>(std::malloc(n * sizeof(double)));
    if (!buf) {
      return 1;
    }
    for (std::size_t i = 0; i < n; ++i) {
      buf[i] = 0.0;
    }
    if (std::strcmp(name, "eigen_sad") == 0) {
      buf[0] = 1.0;
    }
    *cval = buf;
    return 0;
  }
  if (std::strcmp(name, "eigval_sad") == 0) {
    auto *eig = static_cast<double *>(std::malloc(sizeof(double)));
    if (!eig) {
      return 1;
    }
    *eig = -1.0;
    *cval = eig;
    return 0;
  }
  return 1;
}

class MockARTn final : public eonc::IARTnResource {
public:
  void require_loaded() override {}
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }
  [[nodiscard]] artn_create_fn get_create_fn() const override {
    return &artnCreate;
  }
  [[nodiscard]] setup_artn_fn get_setup_fn() const override {
    return &artnSetup;
  }
  [[nodiscard]] artn_fn get_artn_fn() const override { return nullptr; }
  [[nodiscard]] artn_destroy_fn get_destroy_fn() const override {
    return &artnDestroy;
  }
  [[nodiscard]] set_param_fn get_set_param_fn() const override {
    return &artnSetParam;
  }
  [[nodiscard]] get_param_fn get_get_param_fn() const override {
    return &artnGetParam;
  }
  [[nodiscard]] get_runparam_fn get_get_runparam_fn() const override {
    return &artnGetRun;
  }
  [[nodiscard]] get_data_fn get_get_data_fn() const override {
    return &artnGetData;
  }
  [[nodiscard]] print_caller_fn get_print_caller_fn() const override {
    return nullptr;
  }
  [[nodiscard]] artn_step_fn get_artn_step_fn() const override {
    return &artnStep;
  }
  [[nodiscard]] get_error_fn get_get_error_fn() const override {
    return nullptr;
  }
};

void runArtn() {
  CwdGuard cwd(workName("artn"));
  gArtnSteps = 0;
  auto run = beginRun(false);
  eonc::ParametersLoadAccess::saddle_search_options(run.params).method = "artn";
  eonc::ParametersLoadAccess::artn_options(run.params).max_iterations = 3;
  auto matter = twoAtom(run.pot, run.params, 0.0);
  const auto before = run.calls->load();
  eonc::ARTnSaddleSearch search(matter, run.pot, axisMode(), run.params);
  MockARTn library;
  REQUIRE(search.run(library) == eonc::ARTnSaddleSearch::STATUS_GOOD);
  const auto during = run.calls->load() - before;
  REQUIRE(during > 0);
  REQUIRE(static_cast<size_t>(search.getForceCalls()) == during);
  REQUIRE(run.pot->registryForceCalls() == run.calls->load());
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == 0);
}

#ifdef WITH_GPRD
void runGpr() {
  CwdGuard cwd(workName("gpr"));
  auto run = beginRun(false);
  eonc::ParametersLoadAccess::saddle_search_options(run.params).minmode_method =
      eonc::LowestEigenmode::MINMODE_GPRDIMER;
  auto &gpr = eonc::ParametersLoadAccess::gpr_dimer_options(run.params);
  gpr.max_outer_iterations = 1;
  gpr.max_inner_iterations = 1;
  gpr.init_rotations_max = 1;
  gpr.opt_params.max_iterations = 20;
  auto matter = twoAtom(run.pot, run.params, 0.0);
  const auto before = run.calls->load();
  eonc::AtomicGPDimer dimer(matter, run.params, run.pot);
  dimer.compute(matter, axisMode());
  const auto during = run.calls->load() - before;
  REQUIRE(during > 0);
  REQUIRE(static_cast<size_t>(dimer.totalForceCalls) == during);
  REQUIRE(run.pot->registryForceCalls() == run.calls->load());
  REQUIRE(eonc::PotRegistry::get().total_force_calls() == 0);
}
#endif

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

TEST_CASE("every job total matches the counting wrapper", "[force][registry]") {
  SECTION("point") { tests::runPoint(); }
  SECTION("prefactor") { tests::runPrefactor(); }
  SECTION("basin hopping") { tests::runBasin(); }
  SECTION("finite difference") { tests::runFiniteDifference(); }
  SECTION("hessian") { tests::runHessian(); }
  SECTION("monte carlo") { tests::runMonteCarlo(); }
  SECTION("parallel replica") { tests::runParallelReplica(); }
  SECTION("batch") { tests::runBatch(); }
  SECTION("host surface batch") { tests::runHostBatch(); }
  SECTION("force call timer") { tests::runTimer(); }
  SECTION("ARTn") { tests::runArtn(); }
#ifdef WITH_GPRD
  SECTION("GPR dimer") { tests::runGpr(); }
#endif
}
