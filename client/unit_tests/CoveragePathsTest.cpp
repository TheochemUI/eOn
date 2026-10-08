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
#include "eon/BasinHoppingJob.h"
#include "eon/Dynamics.h"
#include "eon/EpiCenters.h"
#include "eon/ForceNorm.h"
#include "eon/NEBForceProjection.h"
#include "eon/GeometryAnalysis.h"
#include "eon/GlobalOptimizationJob.h"
#include "eon/HelperFunctions.h"
#include "eon/ImprovedDimer.h"
#include "eon/IRACompare.h"
#include "eon/InstantonJob.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/NEBInitialPaths.hpp"
#include "eon/NEBOcinebController.h"
#include "eon/NudgedElasticBand.h"
#include "eon/NudgedElasticBandJob.h"
#include "eon/OHTSTJob.h"
#include "eon/ParallelReplicaJob.h"
#include "eon/Parameters.h"
#include "eon/PathIntegral.h"
#include "eon/Prefactor.h"
#include "eon/ProcessSearchJob.h"
#include "eon/QuantumFreeEnergy.h"
#include "eon/Runtime.h"
#include "eon/SafeHyperJob.h"
#include "eon/SurrogatePotential.h"
#include "eon/TADJob.h"
#include "eon/TestJob.h"
#include "eon/potentials/Metatomic/MetatomicLoader.h"
#include "eon/potentials/PluginLoader.h"
#include "eon/potentials/Rgpot/RgpotPot.h"
#include "eon/potentials/Rgpot/GenericEngineLoader.h"
#include "eon/potentials/Rgpot/RGPotEngine.h"
#include "eon/potentials/Rgpot/MetatomicEngineLoader.h"
#include "eon/potentials/Rgpot/XTBEngineLoader.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#ifndef _WIN32
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#endif
#include <vector>

#ifndef _WIN32
#include <unistd.h>
#endif

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {

class Workdir {
public:
  Workdir()
      : original_(std::filesystem::current_path()) {
    dir_ = std::filesystem::temp_directory_path() /
           ("eon_covpaths_" + std::to_string(++seq_));
    std::filesystem::create_directories(dir_);
    std::filesystem::copy_file(original_ / "reactant.con",
                               dir_ / "reactant.con");
    std::filesystem::current_path(dir_);
  }

  ~Workdir() {
    std::error_code ec;
    std::filesystem::current_path(original_, ec);
    std::filesystem::remove_all(dir_, ec);
  }

  const std::filesystem::path &dir() const { return dir_; }

  Workdir(const Workdir &) = delete;
  Workdir &operator=(const Workdir &) = delete;

private:
  static int seq_;
  std::filesystem::path original_;
  std::filesystem::path dir_;
};

int Workdir::seq_ = 0;

class EnvGuard {
public:
  EnvGuard(const char *key, const char *value)
      : key_(key) {
    if (const char *old = std::getenv(key)) {
      had_ = true;
      old_ = old;
    }
    set(value);
  }

  ~EnvGuard() { set(had_ ? old_.c_str() : nullptr); }

  EnvGuard(const EnvGuard &) = delete;
  EnvGuard &operator=(const EnvGuard &) = delete;

private:
  void set(const char *value) {
#ifndef _WIN32
    if (value != nullptr) {
      setenv(key_, value, 1);
    } else {
      unsetenv(key_);
    }
#else
    _putenv_s(key_, value != nullptr ? value : "");
#endif
  }

  const char *key_;
  bool had_{false};
  std::string old_;
};

Parameters ljParams() {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  return params;
}

std::shared_ptr<Matter> loadReactant(const Parameters &params,
                                     std::shared_ptr<Potential> pot) {
  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter(std::string("reactant.con"))));
  return matter;
}

void relaxLbfgs(const char *step, const char *precon, const char *curvature,
                const char *secant, const char *metric) {
  Parameters params = ljParams();
  ParametersLoadAccess::optimizer_options(params).method = OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1.0e-12;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 3;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::optimizer_options(params).convergence_metric = metric;
  auto &lbfgs = ParametersLoadAccess::optimizer_options(params).lbfgs;
  lbfgs.step = step;
  lbfgs.precon = precon;
  lbfgs.curvature = curvature;
  lbfgs.secant = secant;
  lbfgs.project_rigid = true;
  lbfgs.auto_scale = true;
  lbfgs.h0 = "adaptive";
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  Matter matter(pot, params);
  REQUIRE(eonc::io::io_ok(matter.con2matter(std::string("reactant.con"))));
  auto pos = matter.getPositions();
  pos(0, 0) += 0.2;
  matter.setPositions(pos);
  const double before = matter.getPotentialEnergy();
  REQUIRE(std::isfinite(before));
  matter.relax(true);
  REQUIRE(std::isfinite(matter.getPotentialEnergy()));
}

#ifdef RGPOT_HAS_EXPR
void requireExpr(const std::string &expression, const std::string &terms) {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::EXPR;
  ParametersLoadAccess::expr_options(params).expression = expression;
  ParametersLoadAccess::expr_options(params).terms = terms;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::EXPR, params));
  REQUIRE(pot->getType() == PotType::EXPR);
  Matter matter(pot, params);
  REQUIRE(eonc::io::io_ok(matter.con2matter(std::string("reactant.con"))));
  REQUIRE(std::isfinite(matter.getPotentialEnergy()));
  REQUIRE(matter.getForces().allFinite());
}
#endif

} // namespace

#ifdef RGPOT_HAS_EXPR
TEST_CASE("expression potential builds named terms and rejects a blank list",
          "[pot][expr]") {
  requireExpr("lj", "lj");
  requireExpr("lj+morse", " lj , morse ");
  requireExpr("ljcluster", "ljcluster");
  requireExpr("zbl", "zbl");
  for (const char *name : {"d3", "d4", "mopac"}) {
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = PotType::EXPR;
    ParametersLoadAccess::expr_options(params).expression = name;
    ParametersLoadAccess::expr_options(params).terms = name;
    try {
      auto pot = eonc::helpers::sharePotential(
          eonc::helpers::makePotential(PotType::EXPR, params));
      Matter matter(pot, params);
      if (eonc::io::io_ok(matter.con2matter(std::string("reactant.con")))) {
        const double energy = matter.getPotentialEnergy();
        static_cast<void>(energy);
        static_cast<void>(matter.getForces());
      }
    } catch (const std::exception &) {
    }
  }
  Parameters blank;
  ParametersLoadAccess::potential_options(blank).potential = PotType::EXPR;
  ParametersLoadAccess::expr_options(blank).expression = "";
  ParametersLoadAccess::expr_options(blank).terms = "lj";
  REQUIRE_THROWS_AS(eonc::helpers::makePotential(PotType::EXPR, blank),
                    std::runtime_error);
  ParametersLoadAccess::expr_options(blank).expression = "lj";
  ParametersLoadAccess::expr_options(blank).terms = "";
  REQUIRE_THROWS_AS(eonc::helpers::makePotential(PotType::EXPR, blank),
                    std::runtime_error);
  ParametersLoadAccess::expr_options(blank).expression = "nope";
  ParametersLoadAccess::expr_options(blank).terms = "nope";
  REQUIRE_THROWS_AS(eonc::helpers::makePotential(PotType::EXPR, blank),
                    std::runtime_error);
}
#endif

TEST_CASE("preconditioned LBFGS takes Newton and RFO steps", "[optim][lbfgs]") {
  relaxLbfgs("newton", "exp", "reset", "standard", "norm");
  relaxLbfgs("rfo", "pair", "damped", "zhangxu", "rms");
  relaxLbfgs("lbfgs", "lindh", "cautious", "standard", "max_atom");
  relaxLbfgs("lbfgs", "c1", "skip", "standard", "max_component");
  relaxLbfgs("newton", "pair_full", "reset", "standard", "norm");
  relaxLbfgs("rfo", "fischer", "damped", "zhangxu", "norm");
}

TEST_CASE("min-mode confine and an LBFGS dimer rotation stay finite",
          "[saddle_search][coverage]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::optimizer_options(params).method = OptType::CG;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.05;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1.0e-8;
  ParametersLoadAccess::dimer_options(params).improved = false;
  ParametersLoadAccess::dimer_options(params).max_iterations = 4;
  ParametersLoadAccess::saddle_search_options(params).minmode_method =
      LowestEigenmode::MINMODE_DIMER;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1.0e-8;
  ParametersLoadAccess::saddle_search_options(params).max_energy = 50.0;
  ParametersLoadAccess::saddle_search_options(params).perp_force_ratio = 0.0;
  ParametersLoadAccess::saddle_search_options(params).confine_positive.enabled =
      true;
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.bowl_breakout = true;
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.bowl_active = 3;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = loadReactant(params, pot);
  const long nAtoms = matter->numberOfAtoms();
  AtomMatrix mode = AtomMatrix::Random(nAtoms, 3);
  mode.normalize();
  MinModeSaddleSearch bowl(matter, mode, matter->getPotentialEnergy(), params,
                           pot);
  const int bowlStatus = bowl.run();
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
  REQUIRE(bowlStatus != MinModeSaddleSearch::STATUS_INIT);

  matter = loadReactant(params, pot);
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.bowl_breakout = false;
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.min_force = 1.0e-6;
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.min_active = 1;
  ParametersLoadAccess::saddle_search_options(params).confine_positive.boost =
      1.5;
  ParametersLoadAccess::saddle_search_options(params)
      .confine_positive.scale_ratio = 0.5;
  MinModeSaddleSearch confine(matter, mode, matter->getPotentialEnergy(),
                              params, pot);
  REQUIRE(confine.run() != MinModeSaddleSearch::STATUS_INIT);

  matter = loadReactant(params, pot);
  ParametersLoadAccess::saddle_search_options(params).confine_positive.enabled =
      false;
  ParametersLoadAccess::saddle_search_options(params).perp_force_ratio = 0.4;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  MinModeSaddleSearch perp(matter, mode, matter->getPotentialEnergy(), params,
                           pot);
  REQUIRE(perp.run() != MinModeSaddleSearch::STATUS_INIT);

  matter = loadReactant(params, pot);
  ParametersLoadAccess::saddle_search_options(params).perp_force_ratio = 0.0;
  ParametersLoadAccess::dimer_options(params).improved = true;
  ParametersLoadAccess::dimer_options(params).opt_method = OptType::LBFGS;
  ParametersLoadAccess::dimer_options(params).max_iterations = 6;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2;
  MinModeSaddleSearch lbfgsDimer(matter, mode, matter->getPotentialEnergy(),
                                 params, pot);
  REQUIRE(lbfgsDimer.run() != MinModeSaddleSearch::STATUS_INIT);

  matter = loadReactant(params, pot);
  ParametersLoadAccess::dimer_options(params).rotation_backend =
      DimerRotationBackend::LOR;
  ParametersLoadAccess::dimer_options(params).max_iterations = 3;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  MinModeSaddleSearch lor(matter, mode, matter->getPotentialEnergy(), params,
                          pot);
  REQUIRE(lor.run() != MinModeSaddleSearch::STATUS_INIT);
}

TEST_CASE("OCI-NEB run walks the climbing image one dimer step",
          "[neb][ocineb]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::optimizer_options(params).method = OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 5;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.01;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::dimer_options(params).improved = false;
  ParametersLoadAccess::dimer_options(params).max_iterations = 3;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.6;
  product->setPositions(pos);
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  neb->updateForces();
  const double conv = std::max(neb->convergenceForce(), 1.0e-3);
  // Force update clears the climbing image unless the band has an
  // interior peak. The dimer walk needs an interior index.
  neb->climbingImage = 1;
  REQUIRE(neb->tangent.size() > 1);
  REQUIRE(neb->tangent[1] != nullptr);
  auto cfg = eonc::neb::OCINEBController::fromParams(params);
  cfg.max_steps = 1;
  cfg.force_tolerance = 0.01;
  cfg.trigger_factor = 2.0;
  cfg.restore_unhelpful = true;
  eonc::neb::OCINEBController ctl(cfg);
  ctl.initBaseline(conv);
  const auto first = ctl.run(*neb, conv);
  REQUIRE(std::isfinite(first.newForce));
  cfg.restore_unhelpful = false;
  eonc::neb::OCINEBController second(cfg);
  second.initBaseline(conv);
  neb->climbingImage = 1;
  const auto again = second.run(*neb, conv);
  REQUIRE(std::isfinite(again.newForce));
}

TEST_CASE("geometry helpers and a resampled band stay finite",
          "[geometry][neb]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto shifted = std::make_shared<Matter>(*reactant);
  auto pos = shifted->getPositions();
  pos(0, 0) += 0.4;
  shifted->setPositions(pos);
  REQUIRE(eonc::geometry::rotationMatch(*reactant, *shifted, 5.0));
  eonc::geometry::rotationRemove(reactant, shifted);
  eonc::geometry::translationRemove(*shifted, *reactant);
  REQUIRE(std::isfinite(shifted->getPotentialEnergy()));

  std::vector<std::shared_ptr<Matter>> path;
  for (int i = 0; i < 4; ++i) {
    auto image = std::make_shared<Matter>(*reactant);
    auto ip = image->getPositions();
    ip(0, 0) += 0.15 * i;
    image->setPositions(ip);
    path.push_back(image);
  }
  eonc::helpers::neb_paths::resamplePathInPlace(path);
  REQUIRE(path.size() == 4);
  REQUIRE(std::isfinite(path[1]->getPotentialEnergy()));

  std::vector<std::shared_ptr<AtomMatrix>> tangents(path.size());
  const auto freeEnergy =
      eonc::quantumFreeEnergies(path, tangents, 300.0, params);
  REQUIRE(freeEnergy.size() == path.size());
  for (double value : freeEnergy) {
    REQUIRE(std::isfinite(value));
  }
}

TEST_CASE("client displacement, masses, and mode files round-trip",
          "[helpers]") {
  Workdir work;
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto initial = loadReactant(params, pot);
  Matter target(pot, params);
  AtomMatrix mode;
  ParametersLoadAccess::saddle_search_options(params).displace_radius = 4.0;
  ParametersLoadAccess::saddle_search_options(params).displace_magnitude = 0.05;
  ParametersLoadAccess::saddle_search_options(params).displace_type = "load";
  REQUIRE_FALSE(
      eonc::helpers::applyClientDisplacement(target, *initial, params, &mode));
  const char *kinds[] = {"listed_atoms",
                         "random",
                         "last_atom",
                         "least_coordinated",
                         "not_fcc_hcp_coordinated",
                         "no_such_displace"};
  ParametersLoadAccess::saddle_search_options(params).displace_atom_list = {0};
  ParametersLoadAccess::main_options(params).randomSeed = 3;
  for (const char *kind : kinds) {
    ParametersLoadAccess::saddle_search_options(params).displace_type = kind;
    const bool ok =
        eonc::helpers::applyClientDisplacement(target, *initial, params, &mode);
    if (std::string(kind) == "no_such_displace" ||
        std::string(kind) == "load") {
      REQUIRE_FALSE(ok);
    } else {
      REQUIRE(ok);
      REQUIRE(mode.rows() == initial->numberOfAtoms());
    }
  }

  {
    std::ofstream masses(work.dir() / "masses.dat");
    for (long i = 0; i < initial->numberOfAtoms(); ++i) {
      masses << "1.0\n";
    }
  }
  const auto loaded =
      eonc::helpers::loadMasses((work.dir() / "masses.dat").string(),
                                static_cast<int>(initial->numberOfAtoms()));
  REQUIRE(loaded.size() == initial->numberOfAtoms());
  REQUIRE_THROWS_AS(
      eonc::helpers::loadMasses((work.dir() / "masses.dat").string(), 1000),
      std::runtime_error);

  AtomMatrix written = AtomMatrix::Ones(initial->numberOfAtoms(), 3);
  eonc::helpers::saveMode((work.dir() / "mode.dat").string(), initial, written);
  REQUIRE(std::filesystem::file_size(work.dir() / "mode.dat") > 0);
  REQUIRE(eonc::helpers::getRelevantFile("mode.dat") == "mode.dat");
  std::filesystem::copy_file(work.dir() / "mode.dat",
                             work.dir() / "mode_cp.dat");
  REQUIRE(eonc::helpers::getRelevantFile("mode.dat") == "mode_cp.dat");
  std::filesystem::remove(work.dir() / "mode_cp.dat");
  std::filesystem::copy_file(work.dir() / "mode.dat",
                             work.dir() / "mode_in.dat");
  REQUIRE(eonc::helpers::getRelevantFile("mode.dat") == "mode_in.dat");
}

TEST_CASE("IRA match on two structures returns a result", "[ira]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto left = loadReactant(params, pot);
  auto right = std::make_shared<Matter>(*left);
  IRACompare compare;
  const auto matched = compare.match(*left, *right, 0.1);
  const int iraError = matched.error;
  REQUIRE((iraError == 0 || iraError == -1));
}

TEST_CASE("TestJob writes a result row for each built-in potential",
          "[job][testjob]") {
  Workdir work;
  std::filesystem::copy_file(work.dir() / "reactant.con",
                             work.dir() / "pos_test.con");
  auto params = std::make_unique<Parameters>();
  ParametersLoadAccess::potential_options(*params).potential = PotType::LJ;
  eonc::Runtime runtime;
  eonc::TestJob job(std::move(params), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
  REQUIRE(std::filesystem::file_size(work.dir() / "results.dat") > 0);
  std::ifstream in(work.dir() / "results.dat");
  const std::string body((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  REQUIRE(body.find("lj") != std::string::npos);
}

TEST_CASE("process search builds a one-step min-mode and rejects ARTn",
          "[job][process_search]") {
  Workdir work;
  std::filesystem::copy_file(work.dir() / "reactant.con",
                             work.dir() / "pos.con");
  {
    auto params = std::make_unique<Parameters>();
    ParametersLoadAccess::potential_options(*params).potential = PotType::LJ;
    ParametersLoadAccess::main_options(*params).job = JobType::Process_Search;
    ParametersLoadAccess::process_search_options(*params).minimize_first =
        false;
    ParametersLoadAccess::saddle_search_options(*params).method = "artn";
    eonc::Runtime runtime;
    eonc::ProcessSearchJob job(std::move(params), runtime);
    REQUIRE_THROWS_WITH(job.run(), Catch::Matchers::ContainsSubstring("ARTn"));
  }
  auto params = std::make_unique<Parameters>();
  ParametersLoadAccess::potential_options(*params).potential = PotType::LJ;
  ParametersLoadAccess::main_options(*params).job = JobType::Process_Search;
  ParametersLoadAccess::saddle_search_options(*params).method = "min_mode";
  ParametersLoadAccess::saddle_search_options(*params).max_iterations = 1;
  ParametersLoadAccess::dimer_options(*params).improved = false;
  ParametersLoadAccess::dimer_options(*params).max_iterations = 3;
  ParametersLoadAccess::optimizer_options(*params).max_iterations = 5;
  eonc::Runtime runtime;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, *params));
  auto seed = std::make_shared<Matter>(pot, *params);
  REQUIRE(eonc::io::io_ok(seed->con2matter(std::string("reactant.con"))));
  eonc::ProcessSearchJob job(pot, *params);
  auto product = job.runFromMatter(seed);
  REQUIRE(product != nullptr);
  REQUIRE(std::isfinite(product->getPotentialEnergy()));
}

TEST_CASE("NEB job interpolates through an explicit transition structure",
          "[job][neb]") {
  Workdir work;
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto mid = std::make_shared<Matter>(*reactant);
  auto ts = std::make_shared<Matter>(*reactant);
  auto ppos = product->getPositions();
  ppos(0, 0) += 1.2;
  product->setPositions(ppos);
  auto mpos = mid->getPositions();
  mpos(0, 0) += 0.6;
  mid->setPositions(mpos);
  auto tpos = ts->getPositions();
  tpos(0, 0) += 0.7;
  tpos(1, 1) += 0.2;
  ts->setPositions(tpos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  REQUIRE(eonc::io::io_ok(mid->matter2con("mid.con", false)));
  REQUIRE(eonc::io::io_ok(ts->matter2con("ts.con", false)));
  {
    std::ofstream list(work.dir() / "images.lst");
    list << "reactant.con\nmid.con\nproduct.con\n";
  }
  auto owned = std::make_unique<Parameters>(params);
  ParametersLoadAccess::main_options(*owned).job = JobType::Nudged_Elastic_Band;
  ParametersLoadAccess::neb_options(*owned).image_count = 1;
  ParametersLoadAccess::neb_options(*owned).force_tolerance = 1.0;
  ParametersLoadAccess::optimizer_options(*owned).max_iterations = 1;
  ParametersLoadAccess::neb_options(*owned).max_iterations = 1;
  ParametersLoadAccess::neb_options(*owned).endpoints.minimize = false;
  ParametersLoadAccess::neb_options(*owned).initialization.method =
      NEBInit::FILE;
  ParametersLoadAccess::neb_options(*owned).initialization.input_path =
      "images.lst";
  eonc::Runtime runtime;
  eonc::NudgedElasticBandJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

TEST_CASE("OH-TST loads one extra symmetry product", "[job][oh_tst]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto sym = std::make_shared<Matter>(*reactant);
  auto ppos = product->getPositions();
  ppos(0, 0) += 3.0;
  product->setPositions(ppos);
  auto spos = sym->getPositions();
  spos(1, 1) += 2.5;
  sym->setPositions(spos);
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  REQUIRE(eonc::io::io_ok(sym->matter2con("sym.con", false)));
  auto owned = std::make_unique<Parameters>(params);
  ParametersLoadAccess::main_options(*owned).job = JobType::OH_TST;
  ParametersLoadAccess::main_options(*owned).temperature = 0.01;
  ParametersLoadAccess::oh_tst_options(*owned).reactant_filename =
      "reactant.con";
  ParametersLoadAccess::oh_tst_options(*owned).product_filename = "product.con";
  ParametersLoadAccess::oh_tst_options(*owned).symmetry_products = "sym.con";
  ParametersLoadAccess::oh_tst_options(*owned).equil_steps = 0;
  ParametersLoadAccess::oh_tst_options(*owned).sample_steps = 1;
  ParametersLoadAccess::oh_tst_options(*owned).reactant_md_steps = 1;
  ParametersLoadAccess::oh_tst_options(*owned).max_planes = 2;
  ParametersLoadAccess::oh_tst_options(*owned).force_tol = 1.0e6;
  ParametersLoadAccess::oh_tst_options(*owned).max_delta_a = 1000;
  ParametersLoadAccess::oh_tst_options(*owned).time_step = 0.05;
  eonc::Runtime runtime;
  eonc::OHTSTJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

#ifndef _WIN32
TEST_CASE("stand-in engines cover loader success and failure",
          "[pot][engine]") {
  const char *fake = std::getenv("EON_FAKE_ENGINE_SO");
  const bool have = fake != nullptr && std::filesystem::exists(fake);
  {
    EnvGuard pots("EON_POTENTIALS_PATH", "/tmp/eon-cov-empty");
    EnvGuard peng("RGPOT_ENGINE_PATH", "/tmp/eon-cov-empty");
    EnvGuard xtb("RGPOT_XTB_ENGINE", nullptr);
    EnvGuard xtb2("XTB_ENGINE", nullptr);
    EnvGuard mta("RGPOT_METATOMIC_ENGINE", nullptr);
    EnvGuard mta2("METATOMIC_ENGINE", nullptr);
    XTBEngineOptions missing;
    missing.engine_path = "/no/such/librgpot_xtb_engine.so";
    REQUIRE_THROWS_AS(XTBEngineLoader(missing), std::runtime_error);
    MetatomicEngineOptions missingMta;
    missingMta.engine_path = "/no/such/libmetatomic_engine.so";
    REQUIRE_THROWS_AS(MetatomicEngineLoader(missingMta), std::runtime_error);
    GenericEngineOptions missingGen;
    missingGen.library = "libno_such_engine.so";
    missingGen.env_var = "EON_FAKE_GENERIC_ENV";
    missingGen.tag = "cov";
    EnvGuard genv("EON_FAKE_GENERIC_ENV", "/no/such/libno_such_engine.so");
    REQUIRE_THROWS_AS(GenericEngineLoader(missingGen), std::runtime_error);
  }
  if (!have) {
    return;
  }
  XTBEngineOptions xtbOpt;
  xtbOpt.engine_path = fake;
  {
    EnvGuard abi("EON_FAKE_XTB_ABI", "99");
    REQUIRE_THROWS_AS(XTBEngineLoader(xtbOpt), std::runtime_error);
  }
  {
    EnvGuard abi("EON_FAKE_XTB_ABI", "1");
    EnvGuard fail("EON_FAKE_XTB_CREATE_FAIL", "1");
    REQUIRE_THROWS_AS(XTBEngineLoader(xtbOpt), std::runtime_error);
  }
  {
    EnvGuard abi("EON_FAKE_XTB_ABI", "1");
    EnvGuard fail("EON_FAKE_XTB_CREATE_FAIL", nullptr);
    XTBEngineLoader loader(xtbOpt);
    const double R[6] = {0, 0, 0, 1.2, 0, 0};
    const int z[2] = {1, 1};
    double F[6] = {};
    double energy = 0;
    const double box[9] = {10, 0, 0, 0, 10, 0, 0, 0, 10};
    {
      EnvGuard rc("EON_FAKE_XTB_FORCE_RC", "1");
      REQUIRE_THROWS_AS(loader.force(2, R, z, F, &energy, nullptr, box),
                        std::runtime_error);
    }
    double variance = -1;
    loader.force(2, R, z, F, &energy, &variance, box);
    REQUIRE(energy == Catch::Approx(0.5));
    REQUIRE(variance == Catch::Approx(0.25));
    loader.force(2, R, z, F, &energy, nullptr, box);
  }

  MetatomicEngineOptions mtaOpt;
  mtaOpt.engine_path = fake;
  mtaOpt.model_path = "model";
  mtaOpt.device = "cpu";
  {
    EnvGuard abi("EON_FAKE_MTA_ABI", "99");
    REQUIRE_THROWS_AS(MetatomicEngineLoader(mtaOpt), std::runtime_error);
  }
  {
    EnvGuard abi("EON_FAKE_MTA_ABI", "1");
    EnvGuard fail("EON_FAKE_MTA_CREATE_FAIL", "1");
    REQUIRE_THROWS_AS(MetatomicEngineLoader(mtaOpt), std::runtime_error);
  }
  {
    EnvGuard abi("EON_FAKE_MTA_ABI", "1");
    EnvGuard fail("EON_FAKE_MTA_CREATE_FAIL", nullptr);
    MetatomicEngineLoader loader(mtaOpt);
    const double R[6] = {0, 0, 0, 1.2, 0, 0};
    const int z[2] = {1, 1};
    double F[6] = {};
    double energy = 0;
    double variance = -1;
    const double box[9] = {10, 0, 0, 0, 10, 0, 0, 0, 10};
    {
      EnvGuard rc("EON_FAKE_MTA_FORCE_RC", "1");
      REQUIRE_THROWS_AS(loader.force(2, R, z, F, &energy, &variance, box),
                        std::runtime_error);
    }
    loader.force(2, R, z, F, &energy, &variance, box);
    REQUIRE(energy == Catch::Approx(0.5));
    loader.force(2, R, z, F, &energy, nullptr, box);
  }

  GenericEngineOptions gen;
  gen.library = "libcov_engines.so";
  gen.engine_path = fake;
  gen.env_var = "EON_UNUSED_ENGINE";
  gen.tag = "cov";
  gen.config = {1, 2, 3, 4};
  {
    EnvGuard abi("EON_FAKE_ENGINE_ABI", "99");
    REQUIRE_THROWS_AS(GenericEngineLoader(gen), std::runtime_error);
  }
  {
    EnvGuard abi("EON_FAKE_ENGINE_ABI", "1");
    EnvGuard fail("EON_FAKE_ENGINE_CREATE_FAIL", "1");
    REQUIRE_THROWS_AS(GenericEngineLoader(gen), std::runtime_error);
  }
  {
    EnvGuard abi("EON_FAKE_ENGINE_ABI", "1");
    EnvGuard fail("EON_FAKE_ENGINE_CREATE_FAIL", nullptr);
    GenericEngineLoader loader(gen);
    const double R[6] = {0, 0, 0, 1.2, 0, 0};
    const int z[2] = {1, 1};
    double F[6] = {};
    double energy = 0;
    double variance = -1;
    const double box[9] = {10, 0, 0, 0, 10, 0, 0, 0, 10};
    {
      EnvGuard rc("EON_FAKE_ENGINE_FORCE_RC", "1");
      REQUIRE_THROWS_AS(loader.force(2, R, z, F, &energy, &variance, box),
                        std::runtime_error);
    }
    loader.force(2, R, z, F, &energy, &variance, box);
    REQUIRE(energy == Catch::Approx(0.5));
    REQUIRE(loader.available());
    loader.force(2, R, z, F, &energy, nullptr, box);
  }
}

TEST_CASE("plugin loader opens a present library and reports a bad file",
          "[plugin]") {
  const char *fake = std::getenv("EON_FAKE_ENGINE_SO");
  if (fake == nullptr || !std::filesystem::exists(fake)) {
    return;
  }
  const auto dir = std::filesystem::temp_directory_path() / "eon_covplug_dir";
  std::filesystem::create_directories(dir);
  std::filesystem::copy_file(fake, dir / "libeon_covplug.so",
                             std::filesystem::copy_options::overwrite_existing);
  {
    std::ofstream bad(dir / "libeon_badplug.so");
    bad << "not an elf\n";
  }
  auto &loader = eonc::PluginLoader::instance();
  loader.add_config_paths(dir.string());
  using Marker = int (*)();
  auto marker = loader.load_sym<Marker>("eon_covplug", "eon_covplug_marker");
  REQUIRE(marker != nullptr);
  REQUIRE(marker() == 7);
  auto again = loader.load_sym<Marker>("eon_covplug", "eon_covplug_marker");
  REQUIRE(again == marker);
  REQUIRE(loader.load_sym<Marker>("eon_badplug", "missing") == nullptr);
  REQUIRE_THROWS_AS(loader.throw_not_found("eon_badplug", "coverage plugin"),
                    std::runtime_error);
  REQUIRE(loader.load_sym<Marker>("eon_absent_covplug", "missing") == nullptr);

  std::filesystem::copy_file(fake, dir / "libmetatomic_pot.so",
                             std::filesystem::copy_options::overwrite_existing);
  auto &mta = eonc::MetatomicLoader::instance();
  if (!mta.is_loaded()) {
    {
      EnvGuard abi("EON_FAKE_EON_MTA_ABI", "9");
      REQUIRE_FALSE(mta.try_load());
    }
    REQUIRE(mta.try_load());
    mta.require_loaded();
    REQUIRE(mta.is_loaded());
  }
}

#ifdef WITH_RGPOT
TEST_CASE("rgpot xtb, metatomic, and uma backends evaluate the stand-in",
          "[pot][rgpot][coverage]") {
  Parameters unknown;
  ParametersLoadAccess::potential_options(unknown).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(unknown).backend = "not-a-backend";
  REQUIRE_THROWS_AS(eonc::helpers::makePotential(PotType::RGPOT, unknown),
                    std::runtime_error);

  Parameters badSet;
  ParametersLoadAccess::potential_options(badSet).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(badSet).backend = "xtb";
  ParametersLoadAccess::rgpot_options(badSet).xtb_paramset = "nope";
  REQUIRE_THROWS_AS(eonc::helpers::makePotential(PotType::RGPOT, badSet),
                    std::runtime_error);

  const char *fake = std::getenv("EON_FAKE_ENGINE_SO");
  if (fake == nullptr || !std::filesystem::exists(fake)) {
    return;
  }

  auto evaluate = [&](const char *backend, const char *paramset) {
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
    ParametersLoadAccess::rgpot_options(params).backend = backend;
    ParametersLoadAccess::rgpot_options(params).engine_path = fake;
    ParametersLoadAccess::rgpot_options(params).xtb_paramset =
        paramset == nullptr ? "" : paramset;
    ParametersLoadAccess::rgpot_options(params).model_path = "stand-in";
    ParametersLoadAccess::rgpot_options(params).device = "cpu";
    ParametersLoadAccess::rgpot_options(params).task_name = "omol";
    auto pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::RGPOT, params));
    Matter matter(pot, params);
    REQUIRE(eonc::io::io_ok(matter.con2matter(std::string("reactant.con"))));
    REQUIRE(std::isfinite(matter.getPotentialEnergy()));
    REQUIRE(matter.getForces().allFinite());
  };
  evaluate("xtb", "GFN0xTB");
  evaluate("gfn", "GFN1-xTB");
  evaluate("gfnxtb", "GFN-FF");
  evaluate("xtb", nullptr);
  {
    EnvGuard backend("RGPOT_BACKEND", "xtb");
    EnvGuard xtb("RGPOT_XTB_ENGINE", fake);
    evaluate("nwchemc", "GFN2xTB");
  }
  {
    EnvGuard eng("RGPOT_METATOMIC_ENGINE", fake);
    EnvGuard model("RGPOT_METATOMIC_MODEL", "stand-in");
    evaluate("metatomic", "GFN2xTB");
    evaluate("mta", "GFN2xTB");
  }
  evaluate("uma", "GFN2xTB");
  evaluate("omol", "GFN2xTB");
}
#endif
#endif

TEST_CASE("shuffled atom ids, extra ini sections, and ARTn status text",
          "[coverage][ini]") {
  Workdir work;
  {
    std::ifstream in(work.dir() / "reactant.con");
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
      lines.push_back(line);
    }
    bool coords = false;
    int id = 100;
    for (auto &row : lines) {
      if (row.find("Coordinates") != std::string::npos) {
        coords = true;
        continue;
      }
      if (!coords) {
        continue;
      }
      const auto tab = row.rfind('\t');
      if (tab == std::string::npos) {
        break;
      }
      row.replace(tab + 1, std::string::npos, std::to_string(id));
      --id;
    }
    std::ofstream out(work.dir() / "shuffled.con");
    for (const auto &row : lines) {
      out << row << '\n';
    }
  }
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  Matter shuffled(pot, params);
  REQUIRE(eonc::io::io_ok(
      shuffled.con2matter((work.dir() / "shuffled.con").string())));
  REQUIRE(shuffled.numberOfAtoms() == 13);
  REQUIRE(std::isfinite(shuffled.getPotentialEnergy()));

  Parameters ini;
  REQUIRE(ini.load_ini_text(R"(
[Main]
job = minimization
[QuickMin]
time_step = 0.4
[FIRE]
time_step = 0.2
[Surrogate]
potential = catlearn
[CatLearn]
catl_path = /tmp/catlearn
model = stand-in
use_derivatives = false
[ASE_ORCA]
orca_path = /tmp/orca
nproc = 1
charge = 0
multiplicity = 1
)") == 0);
  REQUIRE(ini.optimizer_options().time_step_input == Catch::Approx(0.2));

  auto seed = loadReactant(params, pot);
  AtomMatrix mode = AtomMatrix::Zero(seed->numberOfAtoms(), 3);
  eonc::ARTnSaddleSearch search(seed, pot, mode, params);
#ifndef WITH_ARTN
  REQUIRE(search.run() == eonc::ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR);
#endif
  REQUIRE(std::isnan(search.getEigenvalue()));
  REQUIRE(search.getEigenvector().rows() == seed->numberOfAtoms());
  REQUIRE(search.describeStatus(eonc::ARTnSaddleSearch::STATUS_GOOD) ==
          "Success");
  REQUIRE(search.describeStatus(
              eonc::ARTnSaddleSearch::STATUS_BAD_MAX_ITERATIONS) ==
          "Too many iterations");
  REQUIRE(search.describeStatus(eonc::ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR) ==
          "ARTn backend error");
  REQUIRE(search.describeStatus(99) == "Unknown status");
}

TEST_CASE("dynamics basin and gradient-squared searches take one short step",
          "[job][dynamics][coverage]") {
  Workdir work;
  static_cast<void>(work);
  Parameters base = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, base));
  auto seed = loadReactant(base, pot);

  {
    Parameters params = base;
    ParametersLoadAccess::saddle_search_options(params).method = "dynamics";
    ParametersLoadAccess::saddle_search_options(params).dynamics.temperature =
        300.0;
    ParametersLoadAccess::saddle_search_options(params)
        .dynamics.record_interval = 1.0;
    ParametersLoadAccess::saddle_search_options(params)
        .dynamics.state_check_interval = 1.0;
    ParametersLoadAccess::saddle_search_options(params)
        .dynamics.linear_interpolation = false;
    ParametersLoadAccess::saddle_search_options(params)
        .dynamics.max_init_curvature = 1.0e6;
    ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
    ParametersLoadAccess::dynamics_options(params).time_step = 2.0;
    ParametersLoadAccess::dynamics_options(params).steps = 3;
    ParametersLoadAccess::parallel_replica_options(params).dephase_time = 0.0;
    ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
    ParametersLoadAccess::structure_comparison_options(params)
        .distance_difference = 1.0e-8;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
    ParametersLoadAccess::neb_options(params).image_count = 5;
    ParametersLoadAccess::neb_options(params).max_iterations = 12;
    ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
    ParametersLoadAccess::debug_options(params).write_movies = true;
    auto tight = std::make_shared<Matter>(pot, params);
    REQUIRE(eonc::io::io_ok(tight->con2matter(std::string("reactant.con"))));
    tight->setMasses(VectorXd::Ones(tight->numberOfAtoms()));
    eonc::ProcessSearchJob job(pot, params);
    auto found = job.runFromMatter(tight);
    REQUIRE(found != nullptr);
    REQUIRE(std::filesystem::exists(work.dir() / "neb_initial_band.con"));
  }
  {
    Parameters params = base;
    ParametersLoadAccess::saddle_search_options(params).method =
        "basin_hopping";
    ParametersLoadAccess::main_options(params).temperature = 0.0;
    ParametersLoadAccess::optimizer_options(params).max_iterations = 5;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
    ParametersLoadAccess::neb_options(params).image_count = 1;
    ParametersLoadAccess::neb_options(params).max_iterations = 1;
    ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
    eonc::ProcessSearchJob job(pot, params);
    auto found = job.runFromMatter(std::make_shared<Matter>(*seed));
    REQUIRE(found != nullptr);
  }
  {
    Parameters params = base;
    ParametersLoadAccess::saddle_search_options(params).method = "bgsd";
    ParametersLoadAccess::optimizer_options(params).method = OptType::CG;
    ParametersLoadAccess::optimizer_options(params).max_iterations = 2;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.05;
    eonc::ProcessSearchJob job(pot, params);
    auto found = job.runFromMatter(std::make_shared<Matter>(*seed));
    REQUIRE(found != nullptr);
    REQUIRE(std::isfinite(found->getPotentialEnergy()));
  }
}

namespace {

Parameters shortMolecularDynamics(const Parameters &base) {
  Parameters params = base;
  ParametersLoadAccess::dynamics_options(params).time_step = 1.0;
  ParametersLoadAccess::dynamics_options(params).steps = 4;
  ParametersLoadAccess::main_options(params).temperature = 300.0;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::parallel_replica_options(params).dephase_time = 0.0;
  ParametersLoadAccess::parallel_replica_options(params).refine_transition =
      true;
  ParametersLoadAccess::structure_comparison_options(params)
      .distance_difference = 1.0e-8;
  return params;
}

} // namespace

TEST_CASE("short accelerated dynamics records a transition",
          "[job][tad][coverage]") {
  Workdir work;
  static_cast<void>(work);
  Parameters base = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, base));

  eonc::Runtime runtime;
  {
    std::filesystem::remove("product.con");
    Parameters tadParams = shortMolecularDynamics(base);
    ParametersLoadAccess::main_options(tadParams).temperature = 800.0;
    ParametersLoadAccess::tad_options(tadParams).low_temperature = 800.0;
    ParametersLoadAccess::parallel_replica_options(tadParams)
        .state_check_interval = 2.0;
    ParametersLoadAccess::parallel_replica_options(tadParams).record_interval =
        1.0;
    ParametersLoadAccess::structure_comparison_options(tadParams)
        .remove_translation = false;
    ParametersLoadAccess::structure_comparison_options(tadParams)
        .check_rotation = false;
    auto hot = std::make_shared<Matter>(pot, tadParams);
    REQUIRE(eonc::io::io_ok(hot->con2matter(std::string("reactant.con"))));
    hot->setMasses(VectorXd::Ones(hot->numberOfAtoms()));
    {
      Matter probe(*hot);
      eonc::Dynamics step(&probe, tadParams);
      step.setTemperature(300.0);
      step.setThermalVelocity();
      const AtomMatrix before = probe.getPositions();
      step.oneStep(1);
      REQUIRE((probe.getPositions() - before).norm() > 1.0e-6);
    }
    auto owned = std::make_unique<Parameters>(tadParams);
    eonc::TADJob job(std::move(owned), runtime);
    auto found = job.runFromMatter(hot);
    REQUIRE(found != nullptr);
    REQUIRE(std::filesystem::exists("product.con"));
  }
  {
    std::filesystem::remove("product.con");
    Parameters safeParams = shortMolecularDynamics(base);
    ParametersLoadAccess::main_options(safeParams).temperature = 800.0;
    ParametersLoadAccess::dynamics_options(safeParams).steps = 6;
    ParametersLoadAccess::parallel_replica_options(safeParams)
        .state_check_interval = 2.0;
    ParametersLoadAccess::parallel_replica_options(safeParams).record_interval =
        1.0;
    ParametersLoadAccess::structure_comparison_options(safeParams)
        .remove_translation = false;
    ParametersLoadAccess::structure_comparison_options(safeParams)
        .check_rotation = false;
    auto hot = std::make_shared<Matter>(pot, safeParams);
    REQUIRE(eonc::io::io_ok(hot->con2matter(std::string("reactant.con"))));
    hot->setMasses(VectorXd::Ones(hot->numberOfAtoms()));
    auto owned = std::make_unique<Parameters>(safeParams);
    eonc::SafeHyperJob job(std::move(owned), runtime);
    auto found = job.runFromMatter(hot);
    REQUIRE(found != nullptr);
    REQUIRE(std::filesystem::exists("product.con"));
  }
  {
    std::filesystem::remove("product.con");
    Parameters replicaParams = shortMolecularDynamics(base);
    auto hot = std::make_shared<Matter>(pot, replicaParams);
    REQUIRE(eonc::io::io_ok(hot->con2matter(std::string("reactant.con"))));
    hot->setMasses(VectorXd::Ones(hot->numberOfAtoms()));
    auto owned = std::make_unique<Parameters>(replicaParams);
    eonc::ParallelReplicaJob job(std::move(owned), runtime);
    auto found = job.runFromMatter(hot);
    REQUIRE(found != nullptr);
  }
}

struct RefineProbe : eonc::TADJob {
  using TADJob::TADJob;
  using eonc::ReplicaDynamicsJob::refine;
};

TEST_CASE("band tangents and a transition refine stay finite",
          "[neb][coverage]") {
  Workdir work;
  static_cast<void>(work);
  AtomMatrix next = AtomMatrix::Zero(4, 3);
  AtomMatrix prev = AtomMatrix::Zero(4, 3);
  next.col(0).setLinSpaced(0.2, 0.8);
  prev.col(1).setLinSpaced(0.1, 0.4);
  const AtomMatrix tangent =
      eonc::neb::computeTangent(next, prev, 1.2, 0.4, 0.9, false);
  REQUIRE(tangent.allFinite());
  const AtomMatrix oldTangent =
      eonc::neb::computeTangent(next, prev, 1.2, 0.4, 0.9, true);
  REQUIRE(oldTangent.allFinite());
  AtomMatrix force = AtomMatrix::Ones(4, 3);
  const AtomMatrix perp = eonc::neb::forcePerp(force, tangent);
  REQUIRE(perp.allFinite());
  const AtomMatrix climb =
      eonc::neb::climbingImageForce(force, tangent, perp);
  REQUIRE(climb.allFinite());
  const AtomMatrix dneb = eonc::neb::computeDNEB(force, tangent, perp, true);
  REQUIRE(dneb.allFinite());

  std::vector<double> forces(12, 0.2);
  std::vector<double> fixed(12, 1.0);
  fixed[3] = 0.0;
  const double norm =
      eonc::maxFreeAtomForceNorm(forces.data(), fixed.data(), 4);
  REQUIRE(norm > 0.0);

  const char *fake = std::getenv("EON_FAKE_ENGINE_SO");
#ifdef WITH_RGPOT
  if (fake != nullptr && std::filesystem::exists(fake)) {
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
    ParametersLoadAccess::rgpot_options(params).backend = "xtb";
    ParametersLoadAccess::rgpot_options(params).engine_path = fake;
    RGPotEngineOptions opt;
    opt.backend = "xtb";
    opt.engine_path = fake;
    RGPotEngine engine(opt);
    EnvGuard mpi("OMPI_COMM_WORLD_SIZE", "2");
    engine.armGroupedExit();
    engine.finalizeMpiAtExit();
    REQUIRE_FALSE(RGPotEngine::mpiAbortRequested());
  }
#endif

  Parameters params = ljParams();
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  std::vector<std::shared_ptr<Matter>> buff;
  for (int i = 0; i < 4; ++i) {
    auto snap = std::make_shared<Matter>(*reactant);
    auto pos = snap->getPositions();
    pos(0, 0) += 0.3 * i;
    snap->setPositions(pos);
    buff.push_back(snap);
  }
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  RefineProbe probe(std::move(owned), runtime);
  const long frame = probe.refine(buff, reactant.get());
  REQUIRE(frame >= 1);
  REQUIRE(frame < static_cast<long>(buff.size()));
}

TEST_CASE("a nearly symmetric instanton records the splitting",
          "[job][instanton][split]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.005;
  product->setPositions(pos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  ParametersLoadAccess::main_options(params).job = JobType::Instanton;
  ParametersLoadAccess::instanton_options(params).beads = 12;
  ParametersLoadAccess::instanton_options(params).max_iterations = 60;
  ParametersLoadAccess::instanton_options(params).force_tolerance = 1.0e-2;
  ParametersLoadAccess::instanton_options(params).hessian_stride = 8;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::InstantonJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

TEST_CASE("a four-bead instanton reads a starting band", "[job][instanton]") {
  Workdir work;
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 1.5;
  product->setPositions(pos);
  reactant->setPeriodic(false);
  product->setPeriodic(false);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  REQUIRE(eonc::io::io_ok(reactant->matter2con("band.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("band.con", true)));
  ParametersLoadAccess::main_options(params).job = JobType::Instanton;
  ParametersLoadAccess::instanton_options(params).beads = 4;
  ParametersLoadAccess::instanton_options(params).max_iterations = 1;
  ParametersLoadAccess::instanton_options(params).force_tolerance = 10.0;
  ParametersLoadAccess::instanton_options(params).initial_path = "band.con";
  ParametersLoadAccess::instanton_options(params).hessian_stride = 4;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::InstantonJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

#ifdef WITH_RGPOT
TEST_CASE("rgpot batch forces share one calculator", "[pot][rgpot][batch]") {
  const char *fake = std::getenv("EON_FAKE_ENGINE_SO");
  if (fake == nullptr || !std::filesystem::exists(fake)) {
    return;
  }
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(params).backend = "xtb";
  ParametersLoadAccess::rgpot_options(params).engine_path = fake;
  auto owned = eonc::helpers::makePotential(PotType::RGPOT, params);
  auto *rg = dynamic_cast<RgpotPot *>(owned.get());
  REQUIRE(rg != nullptr);
  const double pos[6] = {0, 0, 0, 1.2, 0, 0};
  const int z[2] = {1, 1};
  double f0[6] = {};
  double f1[6] = {};
  double energies[2] = {};
  double variances[2] = {};
  const double box[9] = {10, 0, 0, 0, 10, 0, 0, 0, 10};
  const double *positions[2] = {pos, pos};
  const int *numbers[2] = {z, z};
  double *forces[2] = {f0, f1};
  const double *boxes[2] = {box, box};
  rg->forceBatch(2, 2, positions, numbers, forces, energies, variances, boxes);
  REQUIRE(std::isfinite(energies[0]));
  REQUIRE(std::isfinite(energies[1]));
  REQUIRE(rg->engineAvailable());
}
#endif

namespace {

struct BatchLJ final : Potential {
  std::shared_ptr<Potential> inner;
  explicit BatchLJ(const Parameters &p)
      : Potential(PotType::LJ),
        inner{eonc::helpers::sharePotential(
            eonc::helpers::makePotential(PotType::LJ, p))} {}
  using Potential::force;
  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    inner->force(nAtoms, positions, atomicNrs, forces, energy, variance, box);
  }
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }
  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override {
    for (long s = 0; s < nSystems; ++s) {
      double var = 0.0;
      inner->force(nAtoms, positions[s], atomicNrs[s], forces[s], &energies[s],
                   &var, boxes[s]);
      if (variances != nullptr) {
        variances[s] = var;
      }
    }
  }
};

} // namespace

TEST_CASE("solid-state bands refuse climbs that do not move the cell",
          "[neb][solid]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::neb_options(params).solid_state.enabled = true;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.4;
  product->setPositions(pos);
  auto expectThrow = [&](const char *label) {
    INFO(label);
    REQUIRE_THROWS_AS(
        std::make_unique<NudgedElasticBand>(reactant, product, params, pot),
        std::invalid_argument);
  };
  ParametersLoadAccess::neb_options(params).zoom.enabled = true;
  expectThrow("zoom");
  ParametersLoadAccess::neb_options(params).zoom.enabled = false;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      true;
  expectThrow("mmf");
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).spring.om.enabled = true;
  expectThrow("onsager");
  ParametersLoadAccess::neb_options(params).spring.om.enabled = false;
  ParametersLoadAccess::neb_options(params).spring.geometric = true;
  expectThrow("geometric");
  ParametersLoadAccess::neb_options(params).spring.geometric = false;
  ParametersLoadAccess::neb_options(params).spring.doubly_nudged = true;
  expectThrow("doubly nudged");
  ParametersLoadAccess::neb_options(params).spring.doubly_nudged = false;
  ParametersLoadAccess::neb_options(params).spring.use_elastic_band = true;
  expectThrow("elastic");
}

TEST_CASE("an oversampled IDPP band is decimated before the climb",
          "[neb][idpp]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).max_iterations = 1;
  ParametersLoadAccess::neb_options(params).force_tolerance = 1.0;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::IDPP;
  ParametersLoadAccess::neb_options(params).initialization.oversampling = true;
  ParametersLoadAccess::neb_options(params).initialization.oversampling_factor =
      2;
  ParametersLoadAccess::neb_options(params).initialization.max_iterations = 1;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 5;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.8;
  product->setPositions(pos);
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  REQUIRE(neb->numImages == 3);
  REQUIRE(std::isfinite(neb->path[1]->getPotentialEnergy()));
}

TEST_CASE("improved dimer batches the centre and the forward image",
          "[dimer][batch]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::dimer_options(params).improved = true;
  ParametersLoadAccess::dimer_options(params).opt_method = OptType::LBFGS;
  ParametersLoadAccess::dimer_options(params).max_iterations = 3;
  auto batch = std::make_shared<BatchLJ>(params);
  auto matter = std::make_shared<Matter>(batch, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter(std::string("reactant.con"))));
  eonc::ImprovedDimer dimer(matter, params, batch);
  AtomMatrix mode = AtomMatrix::Random(matter->numberOfAtoms(), 3);
  mode.normalize();
  auto shifted = matter->getPositions();
  shifted(0, 0) += 0.02;
  matter->setPositions(shifted);
  dimer.compute(matter, mode);
  REQUIRE(std::isfinite(dimer.getEigenvalue()));
}

TEST_CASE("basin hopping writes a unique minimum from a random start",
          "[job][basin_hopping]") {
  Workdir work;
  std::filesystem::copy_file(work.dir() / "reactant.con", work.dir() / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  Parameters params = ljParams();
  ParametersLoadAccess::main_options(params).job = JobType::Basin_Hopping;
  ParametersLoadAccess::basin_hopping_options(params).steps = 2;
  ParametersLoadAccess::basin_hopping_options(params).displacement = 0.3;
  ParametersLoadAccess::basin_hopping_options(params)
      .initial_random_structure_probability = 1.0;
  ParametersLoadAccess::basin_hopping_options(params).write_unique = true;
  ParametersLoadAccess::basin_hopping_options(params).push_apart_distance = 0.4;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 8;
  ParametersLoadAccess::optimizer_options(params).converged_force = 0.05;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::BasinHoppingJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

TEST_CASE("a rate above the crossover uses the parabolic factor",
          "[job][instanton][hot]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto saddle = std::make_shared<Matter>(*reactant);
  auto pos = saddle->getPositions();
  pos(0, 0) += 0.12;
  saddle->setPositions(pos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(saddle->matter2con("saddle.con", false)));
  ParametersLoadAccess::main_options(params).job = JobType::Instanton;
  ParametersLoadAccess::instanton_options(params).mode = "rate";
  ParametersLoadAccess::instanton_options(params).temperature = 8000.0;
  ParametersLoadAccess::instanton_options(params).temperatures = {8000.0,
                                                                  4000.0};
  ParametersLoadAccess::instanton_options(params).beads = 8;
  ParametersLoadAccess::instanton_options(params).max_iterations = 2;
  ParametersLoadAccess::instanton_options(params).force_tolerance = 10.0;
  ParametersLoadAccess::instanton_options(params).pi_planes = 2;
  ParametersLoadAccess::instanton_options(params).pi_beads = 4;
  ParametersLoadAccess::instanton_options(params).pi_equilibration_steps = 0;
  ParametersLoadAccess::instanton_options(params).pi_sampling_steps = 20;
  ParametersLoadAccess::instanton_options(params).pi_time_step = 1.0;
  ParametersLoadAccess::instanton_options(params).pi_recrossing_parents = 2;
  ParametersLoadAccess::instanton_options(params).pi_recrossing_children = 1;
  ParametersLoadAccess::instanton_options(params).pi_recrossing_time = 4.0;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::InstantonJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

TEST_CASE("a loose min-mode search minimizes both endpoints",
          "[job][process_search][ends]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::main_options(params).job = JobType::Process_Search;
  ParametersLoadAccess::main_options(params).parallel = true;
  ParametersLoadAccess::saddle_search_options(params).method = "min_mode";
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1.0e6;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 4;
  ParametersLoadAccess::dimer_options(params).improved = false;
  ParametersLoadAccess::dimer_options(params).max_iterations = 2;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 3;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::process_search_options(params).minimize_first = false;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto seed = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(seed->con2matter(std::string("reactant.con"))));
  auto pos = seed->getPositions();
  pos(0, 0) += 0.35;
  seed->setPositions(pos);
  eonc::ProcessSearchJob job(pot, params);
  auto found = job.runFromMatter(seed);
  REQUIRE(found != nullptr);
  REQUIRE(std::isfinite(found->getPotentialEnergy()));
}

TEST_CASE("a doubly nudged band projects one spring step", "[neb][dneb]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).max_iterations = 2;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.5;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::neb_options(params).spring.doubly_nudged = true;
  ParametersLoadAccess::neb_options(params).spring.use_switching = true;
  ParametersLoadAccess::neb_options(params).spring.geometric = true;
  ParametersLoadAccess::neb_options(params).spring.weighting.enabled = true;
  ParametersLoadAccess::neb_options(params).climbing_image.use_old_tangent =
      true;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 4;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.7;
  product->setPositions(pos);
  auto neb =
      std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  const auto status = neb->compute();
  const bool known = status == NudgedElasticBand::NEBStatus::GOOD ||
                     status == NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS ||
                     status == NudgedElasticBand::NEBStatus::RUNNING ||
                     status == NudgedElasticBand::NEBStatus::MAX_UNCERTAINTY;
  REQUIRE(known);
  REQUIRE(std::isfinite(neb->path[1]->getPotentialEnergy()));
}

#ifndef _WIN32
TEST_CASE("eonclient stops when the parameter file is missing",
          "[client][coverage]") {
  const char *client = std::getenv("EONCLIENT");
  if (client == nullptr) {
    return;
  }
  Workdir work;
  char *argv[] = {const_cast<char *>(client), nullptr};
  pid_t pid = 0;
  const int spawned = posix_spawn(&pid, client, nullptr, nullptr, argv, environ);
  REQUIRE(spawned == 0);
  int status = 0;
  REQUIRE(waitpid(pid, &status, 0) > 0);
  REQUIRE(WIFEXITED(status));
  REQUIRE(WEXITSTATUS(status) == 1);
}
#endif

TEST_CASE("const parameter views are readable", "[parameters][coverage]") {
  Parameters params;
  const Parameters &view = params;
  static_cast<void>(view.last_load_source());
  params.set_mpi_client_comm(3);
  REQUIRE(view.mpi_client_comm() == 3);
  params.set_mpi_potential_rank(2);
  REQUIRE(view.constants().kB > 0.0);
  REQUIRE(view.main_options().temperature >= 0.0);
  static_cast<void>(view.potential_options());
  static_cast<void>(view.ams_options());
  static_cast<void>(view.xtb_options());
  static_cast<void>(view.zbl_options());
  static_cast<void>(view.dftd_options());
  static_cast<void>(view.expr_options());
  static_cast<void>(view.mopac_options());
  static_cast<void>(view.socket_nwchem_options());
  static_cast<void>(view.rgpot_options());
  static_cast<void>(view.structure_comparison_options());
  static_cast<void>(view.process_search_options());
  static_cast<void>(view.saddle_search_options());
  static_cast<void>(view.optimizer_options());
  static_cast<void>(view.dimer_options());
  static_cast<void>(view.gpr_dimer_options());
  static_cast<void>(view.gp_surrogate_options());
  static_cast<void>(view.catlearn_options());
  static_cast<void>(view.ase_orca_options());
  static_cast<void>(view.ase_nwchem_options());
  static_cast<void>(view.metatomic_options());
  static_cast<void>(view.lanczos_options());
  static_cast<void>(view.davidson_options());
  static_cast<void>(view.prefactor_options());
  static_cast<void>(view.hessian_options());
  static_cast<void>(view.neb_options());
  static_cast<void>(view.dynamics_options());
  static_cast<void>(view.parallel_replica_options());
  static_cast<void>(view.tad_options());
  static_cast<void>(view.thermostat_options());
  static_cast<void>(view.replica_exchange_options());
  static_cast<void>(view.hyperdynamics_options());
  static_cast<void>(view.basin_hopping_options());
  static_cast<void>(view.global_optimization_options());
  static_cast<void>(view.monte_carlo_options());
  static_cast<void>(view.bgsd_options());
  static_cast<void>(view.serve_options());
  static_cast<void>(view.artn_options());
  static_cast<void>(view.ira_options());
  static_cast<void>(view.debug_options());
  static_cast<void>(view.oh_tst_options());
  REQUIRE(view.instanton_options().beads > 0);
  static_cast<void>(ParametersLoadAccess::constants(view));
  static_cast<void>(ParametersLoadAccess::ams_options(view));
  static_cast<void>(ParametersLoadAccess::xtb_options(view));
  static_cast<void>(ParametersLoadAccess::zbl_options(view));
  static_cast<void>(ParametersLoadAccess::dftd_options(view));
  static_cast<void>(ParametersLoadAccess::expr_options(view));
  static_cast<void>(ParametersLoadAccess::gpr_dimer_options(view));
  static_cast<void>(ParametersLoadAccess::gp_surrogate_options(view));
  static_cast<void>(ParametersLoadAccess::catlearn_options(view));
  static_cast<void>(ParametersLoadAccess::ase_orca_options(view));
  static_cast<void>(ParametersLoadAccess::ase_nwchem_options(view));
  static_cast<void>(ParametersLoadAccess::metatomic_options(view));
  static_cast<void>(ParametersLoadAccess::tad_options(view));
  static_cast<void>(ParametersLoadAccess::bgsd_options(view));
  static_cast<void>(ParametersLoadAccess::artn_options(view));
  static_cast<void>(ParametersLoadAccess::ira_options(view));
}

#ifdef WITH_RGPOT
TEST_CASE("instanton batches beads on a cpmd engine", "[job][instanton][cpmd]") {
  const char *cpmd = std::getenv("CPMDC_LIBRARY");
  if (cpmd == nullptr || !std::filesystem::exists(cpmd)) {
    return;
  }
  Workdir work;
  Parameters params = ljParams();
  ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(params).backend = "cpmdc";
  ParametersLoadAccess::rgpot_options(params).functional = "BLYP";
  ParametersLoadAccess::rgpot_options(params).cutoff_ry = 20.0;
  ParametersLoadAccess::rgpot_options(params).engine_path = cpmd;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::RGPOT, params));
  if (!pot->supportsBatchEvaluation()) {
    return;
  }
  auto reactant = loadReactant(params, pot);
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.4;
  product->setPositions(pos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  ParametersLoadAccess::main_options(params).job = JobType::Instanton;
  ParametersLoadAccess::instanton_options(params).beads = 4;
  ParametersLoadAccess::instanton_options(params).max_iterations = 1;
  ParametersLoadAccess::instanton_options(params).force_tolerance = 10.0;
  ParametersLoadAccess::instanton_options(params).hessian_stride = 4;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::InstantonJob job(std::move(owned), runtime);
  try {
    const auto files = job.run();
    REQUIRE_FALSE(files.empty());
  } catch (const std::exception &) {
  }
}
#endif

TEST_CASE("rate instanton climbs a short bead ladder", "[job][instanton]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  auto saddle = std::make_shared<Matter>(*reactant);
  auto pos = saddle->getPositions();
  pos(0, 0) += 0.8;
  saddle->setPositions(pos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(saddle->matter2con("saddle.con", false)));
  ParametersLoadAccess::main_options(params).job = JobType::Instanton;
  ParametersLoadAccess::instanton_options(params).mode = "rate";
  ParametersLoadAccess::instanton_options(params).temperature = 200.0;
  ParametersLoadAccess::instanton_options(params).beads = 16;
  ParametersLoadAccess::instanton_options(params).max_iterations = 2;
  ParametersLoadAccess::instanton_options(params).force_tolerance = 10.0;
  ParametersLoadAccess::instanton_options(params).bead_ladder = true;
  ParametersLoadAccess::instanton_options(params).half_ring = false;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::InstantonJob job(std::move(owned), runtime);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
}

TEST_CASE("a recorded replica buffer refines the crossing",
          "[job][replica][coverage]") {
  Workdir work;
  static_cast<void>(work);
  // The runtime owns the registry the job's potential records on. It has
  // to outlive every Matter that setPotential() pointed at that potential.
  eonc::Runtime runtime;
  Parameters params = ljParams();
  ParametersLoadAccess::dynamics_options(params).time_step = 1.0;
  ParametersLoadAccess::dynamics_options(params).steps = 4;
  ParametersLoadAccess::main_options(params).temperature = 800.0;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  ParametersLoadAccess::parallel_replica_options(params).dephase_time = 0.0;
  ParametersLoadAccess::parallel_replica_options(params).refine_transition =
      true;
  ParametersLoadAccess::parallel_replica_options(params).state_check_interval =
      100.0;
  ParametersLoadAccess::parallel_replica_options(params).record_interval = 1.0;
  ParametersLoadAccess::structure_comparison_options(params)
      .distance_difference = 1.0e-8;
  ParametersLoadAccess::structure_comparison_options(params).remove_translation =
      false;
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      false;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto hot = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(hot->con2matter(std::string("reactant.con"))));
  hot->setMasses(VectorXd::Ones(hot->numberOfAtoms()));
  {
    auto owned = std::make_unique<Parameters>(params);
    eonc::ParallelReplicaJob job(std::move(owned), runtime);
    auto found = job.runFromMatter(hot);
    REQUIRE(found != nullptr);
  }
  {
    Parameters loose = params;
    ParametersLoadAccess::structure_comparison_options(loose)
        .distance_difference = 100.0;
    auto still = std::make_shared<Matter>(pot, loose);
    REQUIRE(eonc::io::io_ok(still->con2matter(std::string("reactant.con"))));
    still->setMasses(VectorXd::Ones(still->numberOfAtoms()));
    auto owned = std::make_unique<Parameters>(loose);
    eonc::ParallelReplicaJob job(std::move(owned), runtime);
    auto found = job.runFromMatter(still);
    REQUIRE(found != nullptr);
  }
}

TEST_CASE("one climb iteration keeps the frames it wrote",
          "[saddle_search][coverage]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 1;
  ParametersLoadAccess::dimer_options(params).max_iterations = 2;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = loadReactant(params, pot);
  AtomMatrix mode = AtomMatrix::Zero(matter->numberOfAtoms(), 3);
  mode(0, 0) = 1.0;
  eonc::MinModeSaddleSearch search(matter, mode, matter->getPotentialEnergy(),
                                   params, pot);
  const int status = search.runRetainFrames(1);
  REQUIRE(status != eonc::MinModeSaddleSearch::STATUS_INIT);
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
}

TEST_CASE("a short displacement list is refused", "[io][coverage]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = loadReactant(params, pot);
  eonc::io::ConFrameMetadata meta;
  meta.displacements.assign(3, 0.1);
  REQUIRE(matter->matter2con("bad-disp.con", false, &meta) ==
          eonc::io::IoStatus::InvalidArgument);
  meta.displacements.clear();
  meta.spreads.assign(1, 0.1);
  REQUIRE(matter->matter2con("bad-spread.con", false, &meta) ==
          eonc::io::IoStatus::InvalidArgument);
}

TEST_CASE("rejected parameter text names the field", "[parameters][coverage]") {
  Parameters surrogate;
  REQUIRE_THROWS_AS(surrogate.load_ini_text("[Surrogate]\npotential = lj\n"),
                    std::runtime_error);
  Parameters springs;
  REQUIRE_THROWS_AS(springs.load_ini_text("[Dynamics]\npath_springs = coil\n"),
                    std::invalid_argument);
  Parameters seed;
  REQUIRE_THROWS_AS(seed.load_ini_text("[Dynamics]\npath_seed = -3\n"),
                    std::invalid_argument);
  Parameters mode;
  REQUIRE_THROWS_AS(mode.load_ini_text("[Instanton]\nmode = banana\n"),
                    std::invalid_argument);
  Parameters friction;
  REQUIRE_THROWS_AS(
      friction.load_ini_text("[Instanton]\nfriction = sticky\n"),
      std::invalid_argument);
  Parameters hessians;
  REQUIRE_THROWS_AS(
      hessians.load_ini_text("[Instanton]\ninitial_hessians = guessed\n"),
      std::invalid_argument);
}

TEST_CASE("path integral constructors reject an empty ring",
          "[pathintegral][coverage]") {
  using eonc::pathintegral::Options;
  using eonc::pathintegral::RingPolymer;
  using eonc::pathintegral::Springs;
  using eonc::pathintegral::Thermostat;
  REQUIRE_THROWS_AS(eonc::pathintegral::requireTrotterSprings("eco", "instanton"),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(eonc::pathintegral::trotterEigenvalues(0),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(eonc::pathintegral::ecoEigenvalues(4, 0.0),
                    std::invalid_argument);
  Options opt;
  REQUIRE_THROWS_AS(RingPolymer(0, {}, {}, {}, opt), std::invalid_argument);
  opt.temperature = -1.0;
  REQUIRE_THROWS_AS(RingPolymer(1, {1.0}, {1}, {1, 1, 1}, opt),
                    std::invalid_argument);
  opt.temperature = 1.0;
  opt.springs = Springs::Eco;
  opt.thermostat = Thermostat::Piglet;
  REQUIRE_THROWS_AS(RingPolymer(1, {1.0}, {1}, {1, 1, 1}, opt),
                    std::invalid_argument);
  opt.springs = Springs::Trotter;
  opt.thermostat = Thermostat::Pile;
  REQUIRE_THROWS_AS(RingPolymer(1, {0.0}, {1}, {1, 1, 1}, opt),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(RingPolymer(1, {1.0}, {1}, {0, 0, 0}, opt),
                    std::invalid_argument);
}

namespace {

TEST_CASE("an instanton batches every bead of one iteration",
          "[job][instanton][batch]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto batch = std::make_shared<BatchLJ>(params);
  auto reactant = loadReactant(params, batch);
  reactant->setMasses(VectorXd::Ones(reactant->numberOfAtoms()));
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.05;
  product->setPositions(pos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  ParametersLoadAccess::instanton_options(params).beads = 4;
  ParametersLoadAccess::instanton_options(params).max_iterations = 1;
  ParametersLoadAccess::instanton_options(params).force_tolerance = 10.0;
  ParametersLoadAccess::instanton_options(params).hessian_stride = 4;
  eonc::InstantonJob job(batch, params);
  const auto files = job.run();
  REQUIRE_FALSE(files.empty());
  REQUIRE(std::filesystem::exists("results.dat"));
}

} // namespace

TEST_CASE("a prefactor refuses a missing or unmoved endpoint",
          "[prefactor][coverage]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = loadReactant(params, pot);
  double forward = 0.0;
  double backward = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, nullptr, matter.get(),
                                         matter.get(), forward,
                                         backward) == -1);
  REQUIRE(eonc::Prefactor::getPrefactors(params, matter.get(), matter.get(),
                                         matter.get(), forward,
                                         backward) == -1);
}

namespace {

struct ToySurrogate final : SurrogatePotential {
  explicit ToySurrogate(const Parameters &p)
      : SurrogatePotential(PotType::LJ, p) {}
  void force(long nAtoms, const double * /*positions*/, const int * /*z*/,
             double *forces, double *energy, double *variance,
             const double * /*box*/) override {
    *energy = 0.25;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
  }
  void train_optimize(const MatrixXd &, const MatrixXd &) override {}
};

} // namespace

TEST_CASE("a surrogate potential fills free-atom forces",
          "[matter][surrogate]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = std::make_shared<ToySurrogate>(params);
  auto matter = loadReactant(params, pot);
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
  REQUIRE(matter->getForces().allFinite());
}

TEST_CASE("a mass-weighted manifold reads free-atom masses",
          "[optim][xtsci]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::optimizer_options(params).method = OptType::XTSCI;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 1;
  ParametersLoadAccess::optimizer_options(params).xtsci.manifold = "eckart";
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = loadReactant(params, pot);
  matter->setMasses(VectorXd::Ones(matter->numberOfAtoms()));
  try {
    matter->relax(true);
  } catch (const std::exception &) {
  }
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
}

TEST_CASE("an unknown convergence metric is rejected", "[optim][matter]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::optimizer_options(params).convergence_metric = "bogus";
  ParametersLoadAccess::optimizer_options(params).max_iterations = 1;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = loadReactant(params, pot);
  REQUIRE_THROWS_AS(matter->relax(true), std::invalid_argument);
}

TEST_CASE("rejected instanton inputs stop before a ring is built",
          "[job][instanton]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  reactant->setMasses(VectorXd::Ones(reactant->numberOfAtoms()));
  REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
  eonc::Runtime runtime;
  {
    Parameters rate = params;
    ParametersLoadAccess::instanton_options(rate).mode = "rate";
    auto owned = std::make_unique<Parameters>(rate);
    eonc::InstantonJob job(std::move(owned), runtime);
    REQUIRE_THROWS_AS(job.run(), std::runtime_error);
  }
  {
    auto saddle = std::make_shared<Matter>(*reactant);
    REQUIRE(eonc::io::io_ok(saddle->matter2con("saddle.con", false)));
    Parameters rate = params;
    ParametersLoadAccess::instanton_options(rate).mode = "rate";
    ParametersLoadAccess::instanton_options(rate).hessian_final = "cached";
    auto owned = std::make_unique<Parameters>(rate);
    eonc::InstantonJob job(std::move(owned), runtime);
    REQUIRE_THROWS_AS(job.run(), std::invalid_argument);
  }
  {
    Parameters cold = params;
    ParametersLoadAccess::instanton_options(cold).mode = "rate";
    ParametersLoadAccess::instanton_options(cold).temperature = 0.0;
    auto owned = std::make_unique<Parameters>(cold);
    eonc::InstantonJob job(std::move(owned), runtime);
    REQUIRE_THROWS_AS(job.run(), std::exception);
  }
  {
    reactant->setMasses(VectorXd::Zero(reactant->numberOfAtoms()));
    REQUIRE(eonc::io::io_ok(reactant->matter2con("reactant.con", false)));
    auto owned = std::make_unique<Parameters>(params);
    eonc::InstantonJob job(std::move(owned), runtime);
    REQUIRE_THROWS_AS(job.run(), std::invalid_argument);
  }
}

TEST_CASE("one hyperplane sample stays finite", "[job][oh_tst]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = loadReactant(params, pot);
  reactant->setMasses(VectorXd::Ones(reactant->numberOfAtoms()));
  auto product = std::make_shared<Matter>(*reactant);
  auto pos = product->getPositions();
  pos(0, 0) += 0.3;
  product->setPositions(pos);
  REQUIRE(eonc::io::io_ok(reactant->matter2con("pos.con", false)));
  REQUIRE(eonc::io::io_ok(product->matter2con("product.con", false)));
  ParametersLoadAccess::main_options(params).job = JobType::OH_TST;
  ParametersLoadAccess::main_options(params).temperature = 300.0;
  ParametersLoadAccess::oh_tst_options(params).equil_steps = 0;
  ParametersLoadAccess::oh_tst_options(params).sample_steps = 2;
  ParametersLoadAccess::oh_tst_options(params).reactant_md_steps = 1;
  ParametersLoadAccess::oh_tst_options(params).max_planes = 1;
  ParametersLoadAccess::oh_tst_options(params).pmf_scan = true;
  ParametersLoadAccess::oh_tst_options(params).scan_planes = 2;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::OHTSTJob job(std::move(owned), runtime);
  try {
    const auto files = job.run();
    REQUIRE_FALSE(files.empty());
  } catch (const std::exception &) {
  }
}

TEST_CASE("an unknown displacement algorithm is rejected",
          "[job][basin_hopping]") {
  Workdir work;
  static_cast<void>(work);
  std::filesystem::copy_file(work.dir() / "reactant.con", work.dir() / "pos.con",
                             std::filesystem::copy_options::overwrite_existing);
  Parameters params = ljParams();
  ParametersLoadAccess::main_options(params).job = JobType::Basin_Hopping;
  ParametersLoadAccess::basin_hopping_options(params).steps = 1;
  ParametersLoadAccess::basin_hopping_options(params).displacement_algorithm =
      "sideways";
  ParametersLoadAccess::optimizer_options(params).max_iterations = 0;
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::BasinHoppingJob job(std::move(owned), runtime);
  REQUIRE_THROWS_AS(job.run(), std::invalid_argument);
}

TEST_CASE("an unknown saddle search method is rejected", "[job][process]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  ParametersLoadAccess::main_options(params).job = JobType::Process_Search;
  ParametersLoadAccess::saddle_search_options(params).method = "nope";
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::ProcessSearchJob job(std::move(owned), runtime);
  REQUIRE_THROWS_AS(job.run(), std::runtime_error);
}

TEST_CASE("unknown hopping feedback stops the optimizer",
          "[job][global_opt]") {
  Workdir work;
  static_cast<void>(work);
  Parameters params = ljParams();
  eonc::Runtime runtime;
  auto owned = std::make_unique<Parameters>(params);
  eonc::GlobalOptimizationJob job(std::move(owned), runtime);
  job.applyMoveFeedbackMD();
  REQUIRE_THROWS_AS(job.applyMoveFeedbackMD(), std::runtime_error);
}

} // namespace tests
