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
#include "eon/ProcessSearchJob.h"
#include "eon/PotCapabilities.h"
#ifdef WITH_ARTN
#include "eon/ARTnSaddleSearch.h"
#endif
#include "eon/BasinHoppingSaddleSearch.h"
#include "eon/BiasedGradientSquaredDescent.h"
#include "eon/DynamicsSaddleSearch.h"
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Optimizer.h"
#include "eon/Prefactor.h"
#include <exception>
#include <filesystem>
#include <thread>

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "eon/ConFileIO.h"
#include "eon/EonLogger.h"
#include "eon/JobResult.h"

namespace eonc {

std::vector<std::string> ProcessSearchJob::run() {
  std::string reactantFilename = eonc::helpers::getRelevantFile("pos.con");
  std::string displacementFilename("displacement.con");
  std::string modeFilename("direction.dat");
  size_t fctmp{0};
  initial = std::make_shared<Matter>(pot, params);
  if (params.saddle_search_options().method == "min_mode" ||
      params.saddle_search_options().method == "basin_hopping" ||
      params.saddle_search_options().method == "bgsd") {
    displacement = std::make_shared<Matter>(pot, params);
  } else if (params.saddle_search_options().method == "dynamics") {
    displacement = nullptr;
  }
  saddle = std::make_shared<Matter>(pot, params);
  // Give min2 its own potential for parallel endpoint minimization.
  // A clone keeps this job's potential; makePotential rebuilds from the
  // configuration and is the fallback for backends that cannot clone.
  std::shared_ptr<Potential> min2Pot = pot;
  if (pot->needsPerImageInstance() && params.main_options().parallel) {
    auto cloned = pot->clonePotential();
    min2Pot = cloned ? cloned
                     : eonc::helpers::sharePotential(
                           eonc::helpers::makePotential(params));
  }
  min1 = std::make_shared<Matter>(pot, params);
  min2 = std::make_shared<Matter>(min2Pot, params);

  if (!eonc::io::io_ok(initial->con2matter(reactantFilename))) {
    EONC_LOG_CRITICAL("Failed to load {}", reactantFilename);
    throw std::runtime_error("failed to load " + reactantFilename);
  }

  if (params.process_search_options().minimize_first) {
    QUILL_LOG_DEBUG(log, "Minimizing initial structure\n");
    fctmp = initial->getPotentialCalls();
    initial->relax();
    fCallsMin += initial->getPotentialCalls() - fctmp;
    QUILL_LOG_DEBUG(log, "Initial minimization took {} fcalls",
                    initial->getPotentialCalls() - fctmp);
  }

  barriersValues[0] = barriersValues[1] = 0;
  prefactorsValues[0] = prefactorsValues[1] = 0;

  AtomMatrix mode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
  if (params.saddle_search_options().method == "min_mode" ||
      params.saddle_search_options().method == "basin_hopping" ||
      params.saddle_search_options().method == "bgsd") {
    if (params.saddle_search_options().displace_type ==
        eonc::EpiCenters::DISP_LOAD) {
      // Load displacement.con, or synthesize from pos.con + direction.dat
      // (#79).
      if (!eonc::helpers::loadOrSynthesizeDisplacement(
              *saddle, *initial, displacementFilename, modeFilename,
              params.saddle_search_options().displace_magnitude)) {
        EONC_LOG_CRITICAL("Failed to load {} (and no usable {})",
                          displacementFilename, modeFilename);
        throw std::runtime_error("failed to load " + displacementFilename);
      }
      *min1 = *min2 = *initial;
    } else if (eonc::helpers::applyClientDisplacement(*saddle, *initial, params,
                                                      &mode)) {
      *min1 = *min2 = *initial;
    } else {
      *saddle = *min1 = *min2 = *initial;
    }
    if (displacement) {
      *displacement = *saddle;
    }
  } else {
    // ARTn and dynamics start from the initial minimum
    *saddle = *min1 = *min2 = *initial;
  }
  min2->setPotential(min2Pot);

  const bool useARTnAsMinMode =
      params.saddle_search_options().method == "min_mode" &&
      params.saddle_search_options().minmode_method == "artn";

  if (params.saddle_search_options().method == "min_mode") {
    if (params.saddle_search_options().displace_type ==
            eonc::EpiCenters::DISP_LOAD &&
        std::filesystem::exists(modeFilename)) {
      mode = eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
    }
#ifdef WITH_ARTN
    // ARTn as a min-mode drop-in: eOn displaces, seeds the mode, ARTn
    // takes over from the displaced structure.
    if (useARTnAsMinMode) {
      saddleSearch =
          std::make_unique<ARTnSaddleSearch>(saddle, pot, mode, params);
    } else
#endif
    {
      saddleSearch = std::make_unique<MinModeSaddleSearch>(
          saddle, mode, initial->getPotentialEnergy(), params, pot);
    }
#ifdef WITH_ARTN
  } else if (params.saddle_search_options().method == "artn") {
    // ARTn handles its own push from the minimum, eigenmode estimation,
    // and perpendicular relaxation internally.
    AtomMatrix artnMode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
    if (params.saddle_search_options().displace_type ==
            eonc::EpiCenters::DISP_LOAD &&
        std::filesystem::exists(modeFilename)) {
      artnMode =
          eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
    }
    saddleSearch =
        std::make_unique<ARTnSaddleSearch>(saddle, pot, artnMode, params);
#endif
  } else if (params.saddle_search_options().method == "basin_hopping") {
    saddleSearch =
        std::make_unique<BasinHoppingSaddleSearch>(min1, saddle, pot, params);
  } else if (params.saddle_search_options().method == "dynamics") {
    saddleSearch = std::make_unique<DynamicsSaddleSearch>(saddle, params);
  } else if (params.saddle_search_options().method == "bgsd") {
    saddleSearch = std::make_unique<BiasedGradientSquaredDescent>(
        saddle, initial->getPotentialEnergy(), params);
  }

#ifndef WITH_ARTN
  // Post-dispatch guard for both ARTn entry points so users without a
  // WITH_ARTN build get a clean per-case error instead of a silent
  // fall-through. Two distinct messages so downstream tooling and the
  // integration tests can match on the specific entry point.
  if (params.saddle_search_options().method == "artn") {
    throw std::runtime_error(
        "saddle_search.method=artn requires a build with ARTn support "
        "(reconfigure with -Dwith_artn=true)");
  }
  if (useARTnAsMinMode) {
    throw std::runtime_error(
        "saddle_search.minmode_method=artn requires a build with ARTn "
        "support (reconfigure with -Dwith_artn=true)");
  }
#endif

  if (!saddleSearch) {
    throw std::runtime_error("unknown saddle_search.method");
  }

  (void)runPrepared();
  return returnFiles;
}

std::shared_ptr<Matter>
ProcessSearchJob::runFromMatter(std::shared_ptr<Matter> seed) {
  if (!seed) {
    throw std::runtime_error("ProcessSearchJob::runFromMatter: null Matter");
  }
  initial = seed;
  initial->setPotential(pot);
  // A clone keeps this job's potential; makePotential rebuilds from the
  // configuration and is the fallback for backends that cannot clone.
  std::shared_ptr<Potential> min2Pot = pot;
  if (pot->needsPerImageInstance() && params.main_options().parallel) {
    auto cloned = pot->clonePotential();
    min2Pot = cloned ? cloned
                     : eonc::helpers::sharePotential(
                           eonc::helpers::makePotential(params));
  }
  displacement = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);
  min1 = std::make_shared<Matter>(pot, params);
  min2 = std::make_shared<Matter>(min2Pot, params);
  AtomMatrix mode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
  if (!eonc::helpers::applyClientDisplacement(*saddle, *initial, params,
                                              &mode)) {
    *saddle = *initial;
  }
  *displacement = *saddle;
  *min1 = *min2 = *initial;
  min2->setPotential(min2Pot);
  if (params.saddle_search_options().method == "min_mode") {
    saddleSearch = std::make_unique<MinModeSaddleSearch>(
        saddle, mode, initial->getPotentialEnergy(), params, pot);
  } else if (params.saddle_search_options().method == "basin_hopping") {
    saddleSearch =
        std::make_unique<BasinHoppingSaddleSearch>(min1, saddle, pot, params);
  } else if (params.saddle_search_options().method == "dynamics") {
    saddleSearch = std::make_unique<DynamicsSaddleSearch>(saddle, params);
  } else if (params.saddle_search_options().method == "bgsd") {
    saddleSearch = std::make_unique<BiasedGradientSquaredDescent>(
        saddle, initial->getPotentialEnergy(), params);
  } else {
    throw std::runtime_error(
        "ProcessSearchJob::runFromMatter: unsupported saddle_search.method");
  }
  return runPrepared();
}

std::shared_ptr<Matter> ProcessSearchJob::runPrepared() {
  if (!saddleSearch) {
    throw std::runtime_error("unknown saddle_search.method");
  }

  int status = doProcessSearch();

  printEndState(status);
  saveData(status);

  return min2 ? min2 : saddle;
}

int ProcessSearchJob::doProcessSearch() {
  Matter matterTemp(pot, params);
  long status;
  size_t fctmp{0};

  fctmp = pot->forceCallCounter;
  status = saddleSearch->run();
  if (params.saddle_search_options().method == "min_mode" &&
      params.saddle_search_options().minmode_method ==
          LowestEigenmode::MINMODE_GPRDIMER) {
    fCallsSaddle += saddleSearch->getForceCalls();
  } else if (params.saddle_search_options().method == "artn") {
    fCallsSaddle += saddleSearch->getForceCalls();
  } else {
    fCallsSaddle += pot->forceCallCounter - fctmp;
  }
  EONC_LOG_DEBUG("Got {} calls in the saddle search, with previous {}",
                 fCallsSaddle, fctmp);

  if (status != MinModeSaddleSearch::STATUS_GOOD) {
    return status;
  }

  AtomMatrix posSaddle = saddle->getPositions();
  AtomMatrix displacedPos;

  // Matter's copy assignment copies the potential; keep each endpoint's
  // own instance so the two minimizations can run at the same time.
  const auto min1Pot = min1->getPotential();
  const auto min2Pot = min2->getPotential();
  *min1 = *saddle;
  min1->setPotential(min1Pot);

  displacedPos =
      posSaddle - saddleSearch->getEigenvector() *
                      params.process_search_options().minimization_offset;
  min1->setPositions(displacedPos);

  *min2 = *saddle;
  min2->setPotential(min2Pot);
  displacedPos =
      posSaddle + saddleSearch->getEigenvector() *
                      params.process_search_options().minimization_offset;
  min2->setPositions(displacedPos);

  // Minimize both endpoints concurrently when the shared potential instance is
  // safe to call from multiple threads, or when each endpoint owns a separate
  // potential instance.
  QUILL_LOG_DEBUG(log, "Starting Minimization 1 & 2");
  bool converged1{false}, converged2{false};
  long fc1_before = min1->getPotentialCalls();
  long fc2_before = min2->getPotentialCalls();

  // Two threads may share an instance only when it is thread safe; a
  // per-image potential needs the endpoints to hold distinct instances.
  bool canParallel = eonc::potAllowsSharedInstance(*pot) ||
                     (pot->needsPerImageInstance() &&
                      min1->getPotential().get() != min2->getPotential().get());
  if (params.main_options().parallel && canParallel) {
    // An exception may not leave a thread function (std::terminate), and a
    // joinable std::thread may not be destroyed: t1 hands its error back and
    // the caller joins before rethrowing either side's.
    std::exception_ptr t1Error;
    std::thread t1([&] {
      try {
        converged1 = min1->relax(false, params.debug_options().write_movies,
                                 false, "min1");
      } catch (...) {
        t1Error = std::current_exception();
      }
    });
    try {
      converged2 = min2->relax(false, params.debug_options().write_movies,
                               false, "min2");
    } catch (...) {
      t1.join();
      throw;
    }
    t1.join();
    if (t1Error)
      std::rethrow_exception(t1Error);
  } else {
    converged1 =
        min1->relax(false, params.debug_options().write_movies, false, "min1");
    converged2 =
        min2->relax(false, params.debug_options().write_movies, false, "min2");
  }

  if (min1->getPotential().get() == min2->getPotential().get()) {
    fCallsMin += min1->getPotentialCalls() - fc1_before;
  } else {
    fCallsMin += (min1->getPotentialCalls() - fc1_before) +
                 (min2->getPotentialCalls() - fc2_before);
  }
  QUILL_LOG_DEBUG(log, "Min1: {} fcalls, Min2: {} fcalls",
                  min1->getPotentialCalls() - fc1_before,
                  min2->getPotentialCalls() - fc2_before);

  if (!converged1 || !converged2) {
    return MinModeSaddleSearch::STATUS_BAD_MINIMA;
  }

  auto sameAs = [](const Matter &a, const Matter &b) {
    Matter probe(a);
    return probe.compare(b);
  };

  if (!sameAs(*initial, *min1) && sameAs(*initial, *min2)) {
    matterTemp = *min1;
    *min1 = *min2;
    *min2 = matterTemp;
  }

  if (!sameAs(*initial, *min1)) {
    // Report how far off the endpoint landed. Whether the minimisation
    // stopped just outside the state-identity tolerance or relaxed into a
    // different state entirely calls for opposite fixes, and the status
    // alone does not distinguish them.
    const double tol =
        params.structure_comparison_options().distance_difference;
    auto countMoved = [&](const Matter &m) {
      long moved = 0;
      for (long i = 0; i < initial->numberOfAtoms(); ++i) {
        if (initial
                ->pbc(initial->getPositions().row(i) - m.getPositions().row(i))
                .norm() > tol) {
          ++moved;
        }
      }
      return moved;
    };
    QUILL_LOG_INFO(log,
                   "initial != min1: {} of {} atoms past the {} A tolerance "
                   "for min1 ({} for min2); largest separation {} A",
                   countMoved(*min1), initial->numberOfAtoms(), tol,
                   countMoved(*min2), initial->perAtomNorm(*min1));
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  if (sameAs(*initial, *min2)) {
    QUILL_LOG_DEBUG(log, "both minima are the initial state");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  if (!params.process_search_options().minimize_first) {
    min1 = initial;
  }

  barriersValues[0] = saddle->getPotentialEnergy() - min1->getPotentialEnergy();
  barriersValues[1] = saddle->getPotentialEnergy() - min2->getPotentialEnergy();

  if ((params.saddle_search_options().max_energy < barriersValues[0]) ||
      (params.saddle_search_options().max_energy < barriersValues[1])) {
    return MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER;
  }

  if (barriersValues[0] < 0.0 || barriersValues[1] < 0.0) {
    return MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER;
  }

  if (!params.prefactor_options().default_value) {
    fctmp = min1->getPotentialCalls();
    int prefStatus;
    double pref1, pref2;
    prefStatus = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2);
    if (prefStatus == -1) {
      EONC_LOG_ERROR("Prefactor: bad calculation");
      return MinModeSaddleSearch::STATUS_FAILED_PREFACTOR;
    }
    fCallsPrefactors += min1->getPotentialCalls() - fctmp;

    if ((pref1 > params.prefactor_options().max_value) ||
        (pref1 < params.prefactor_options().min_value)) {
      EONC_LOG_ERROR("Bad reactant-to-saddle prefactor: {}", pref1);
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    if ((pref2 > params.prefactor_options().max_value) ||
        (pref2 < params.prefactor_options().min_value)) {
      EONC_LOG_ERROR("Bad product-to-saddle prefactor: {}", pref2);
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    prefactorsValues[0] = pref1;
    prefactorsValues[1] = pref2;

  } else {
    prefactorsValues[0] = params.prefactor_options().default_value;
    prefactorsValues[1] = params.prefactor_options().default_value;
  }
  return MinModeSaddleSearch::STATUS_GOOD;
}

void ProcessSearchJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  const bool dynamics = params.saddle_search_options().method == "dynamics";
  double simTime = 0.0;
  double mdTemp = 0.0;
  if (dynamics) {
    auto ds = dynamic_cast<DynamicsSaddleSearch &>(*saddleSearch);
    simTime = ds.time * params.constants().timeUnit;
    mdTemp = params.saddle_search_options().dynamics.temperature;
  }
  const double disp = params.saddle_search_options().method == "min_mode"
                          ? displacement->perAtomNorm(*saddle)
                          : 0.0;
  auto env = JobResultEnvelope::fromProcessSearch(
      status, std::string(saddleSearch->describeStatus(status)),
      params.potential_options().potential, params.main_options().randomSeed,
      fCallsMin, fCallsSaddle, fCallsPrefactors, saddle->getPotentialEnergy(),
      min1->getPotentialEnergy(), min2->getPotentialEnergy(), barriersValues[0],
      barriersValues[1], disp, prefactorsValues[0], prefactorsValues[1],
      dynamics, simTime, mdTemp);
  env.reactant_frame = eonc::io::matterToConFrame(*min1);
  env.saddle_frame = eonc::io::matterToConFrame(*saddle);
  env.product_frame = eonc::io::matterToConFrame(*min2);
  env.writeResultsDat(resultsFilename);

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  if (!eonc::io::io_ok(min1->matter2con(reactantFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", reactantFilename);
  }

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  eonc::helpers::saveMode(modeFilename, saddle, saddleSearch->getEigenvector());

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  if (!eonc::io::io_ok(saddle->matter2con(saddleFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", saddleFilename);
  }

  std::string productFilename("product.con");
  returnFiles.push_back(productFilename);
  if (!eonc::io::io_ok(min2->matter2con(productFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", productFilename);
  }
}

void ProcessSearchJob::printEndState(int status) {
  auto msg = saddleSearch->describeStatus(status);
  if (status == MinModeSaddleSearch::STATUS_GOOD) {
    QUILL_LOG_DEBUG(log, "[Saddle Search] {}", msg);
  } else {
    QUILL_LOG_ERROR(log, "[Saddle Search] {}", msg);
  }
}

} // namespace eonc
