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
#include "ProcessSearchJob.h"
#ifdef WITH_ARTN
#include "ARTnSaddleSearch.h"
#endif
#include "BasinHoppingSaddleSearch.h"
#include "BiasedGradientSquaredDescent.h"
#include "DynamicsSaddleSearch.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "MinModeSaddleSearch.h"
#include "Optimizer.h"
#include "Prefactor.h"
#include <thread>

#include <format>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "EonLogger.h"

std::vector<std::string> ProcessSearchJob::run() {
  std::string reactantFilename("pos.con");
  std::string displacementFilename("displacement.con");
  std::string modeFilename("direction.dat");
  size_t fctmp{0};
  initial = std::make_shared<Matter>(pot, params);
  if (params.saddle_search_options.method == "min_mode" ||
      params.saddle_search_options.method == "basin_hopping" ||
      params.saddle_search_options.method == "bgsd") {
    displacement = std::make_shared<Matter>(pot, params);
  } else if (params.saddle_search_options.method == "dynamics") {
    displacement = nullptr;
  }
  saddle = std::make_shared<Matter>(pot, params);
  // Give min2 its own potential for parallel endpoint minimization
  auto min2Pot = (pot->needsPerImageInstance() && params.main_options.parallel)
                     ? eonc::helpers::makePotential(params)
                     : pot;
  min1 = std::make_shared<Matter>(pot, params);
  min2 = std::make_shared<Matter>(min2Pot, params);

  if (!initial->con2matter(reactantFilename)) {
    EONC_LOG_CRITICAL("Failed to load {}", reactantFilename);
    exit(1);
  }

  if (params.process_search_options.minimize_first) {
    QUILL_LOG_DEBUG(log, "Minimizing initial structure\n");
    fctmp = initial->getPotentialCalls();
    initial->relax();
    fCallsMin += initial->getPotentialCalls() - fctmp;
    QUILL_LOG_DEBUG(log, "Initial minimization took {} fcalls",
                    initial->getPotentialCalls() - fctmp);
  }

  barriersValues[0] = barriersValues[1] = 0;
  prefactorsValues[0] = prefactorsValues[1] = 0;

  if (params.saddle_search_options.method == "min_mode" ||
      params.saddle_search_options.method == "basin_hopping" ||
      params.saddle_search_options.method == "bgsd") {
    if (params.saddle_search_options.displace_type ==
        eonc::EpiCenters::DISP_LOAD) {
      if (!saddle->con2matter(displacementFilename)) {
        EONC_LOG_CRITICAL("Failed to load {}", displacementFilename);
        exit(1);
      }
      *min1 = *min2 = *initial;
    } else {
      *saddle = *min1 = *min2 = *initial;
    }
  } else {
    // ARTn and dynamics start from the initial minimum
    *saddle = *min1 = *min2 = *initial;
  }

  AtomMatrix mode;
  const bool useARTnAsMinMode =
      params.saddle_search_options.method == "min_mode" &&
      params.saddle_search_options.minmode_method == "artn";

  if (params.saddle_search_options.method == "min_mode") {
    if (params.saddle_search_options.displace_type ==
        eonc::EpiCenters::DISP_LOAD) {
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
  } else if (params.saddle_search_options.method == "artn") {
    // ARTn handles its own push from the minimum, eigenmode estimation,
    // and perpendicular relaxation internally.
    AtomMatrix artnMode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
    if (params.saddle_search_options.displace_type ==
        eonc::EpiCenters::DISP_LOAD) {
      artnMode =
          eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
    }
    saddleSearch =
        std::make_unique<ARTnSaddleSearch>(saddle, pot, artnMode, params);
#endif
  } else if (params.saddle_search_options.method == "basin_hopping") {
    saddleSearch =
        std::make_unique<BasinHoppingSaddleSearch>(min1, saddle, pot, params);
  } else if (params.saddle_search_options.method == "dynamics") {
    saddleSearch = std::make_unique<DynamicsSaddleSearch>(saddle, params);
  } else if (params.saddle_search_options.method == "bgsd") {
    saddleSearch = std::make_unique<BiasedGradientSquaredDescent>(
        saddle, initial->getPotentialEnergy(), params);
  }

#ifndef WITH_ARTN
  // Post-dispatch guard for both ARTn entry points so users without a
  // WITH_ARTN build get a clean per-case error instead of a silent
  // fall-through. Two distinct messages so downstream tooling and the
  // integration tests can match on the specific entry point.
  if (params.saddle_search_options.method == "artn") {
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

  int status = doProcessSearch();

  printEndState(status);
  saveData(status);

  return returnFiles;
}

int ProcessSearchJob::doProcessSearch() {
  Matter matterTemp(pot, params);
  long status;
  size_t fctmp{0};

  fctmp = pot->forceCallCounter;
  status = saddleSearch->run();
  if (params.saddle_search_options.method == "min_mode" &&
      params.saddle_search_options.minmode_method ==
          LowestEigenmode::MINMODE_GPRDIMER) {
    fCallsSaddle += saddleSearch->getForceCalls();
  } else if (params.saddle_search_options.method == "artn") {
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

  *min1 = *saddle;

  displacedPos =
      posSaddle - saddleSearch->getEigenvector() *
                      params.process_search_options.minimization_offset;
  min1->setPositions(displacedPos);

  *min2 = *saddle;
  displacedPos =
      posSaddle + saddleSearch->getEigenvector() *
                      params.process_search_options.minimization_offset;
  min2->setPositions(displacedPos);

  // Minimize both endpoints concurrently (independent Matter objects)
  QUILL_LOG_DEBUG(log, "Starting Minimization 1 & 2");
  bool converged1{false}, converged2{false};
  long fc1_before = min1->getPotentialCalls();
  long fc2_before = min2->getPotentialCalls();

  bool canParallel = pot->isThreadSafe() || pot->needsPerImageInstance();
  if (params.main_options.parallel && canParallel) {
    std::thread t1([&] {
      converged1 =
          min1->relax(false, params.debug_options.write_movies, false, "min1");
    });
    converged2 =
        min2->relax(false, params.debug_options.write_movies, false, "min2");
    t1.join();
  } else {
    converged1 =
        min1->relax(false, params.debug_options.write_movies, false, "min1");
    converged2 =
        min2->relax(false, params.debug_options.write_movies, false, "min2");
  }

  fCallsMin += (min1->getPotentialCalls() - fc1_before) +
               (min2->getPotentialCalls() - fc2_before);
  QUILL_LOG_DEBUG(log, "Min1: {} fcalls, Min2: {} fcalls",
                  min1->getPotentialCalls() - fc1_before,
                  min2->getPotentialCalls() - fc2_before);

  if (!converged1 || !converged2) {
    return MinModeSaddleSearch::STATUS_BAD_MINIMA;
  }

  if (!(initial->compare(*min1)) && initial->compare(*min2)) {
    matterTemp = *min1;
    *min1 = *min2;
    *min2 = matterTemp;
  }

  if (!initial->compare(*min1)) {
    QUILL_LOG_DEBUG(log, "initial != min1");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  if (initial->compare(*min2)) {
    QUILL_LOG_DEBUG(log, "both minima are the initial state");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  if (!params.process_search_options.minimize_first) {
    min1 = initial;
  }

  barriersValues[0] = saddle->getPotentialEnergy() - min1->getPotentialEnergy();
  barriersValues[1] = saddle->getPotentialEnergy() - min2->getPotentialEnergy();

  if ((params.saddle_search_options.max_energy < barriersValues[0]) ||
      (params.saddle_search_options.max_energy < barriersValues[1])) {
    return MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER;
  }

  if (barriersValues[0] < 0.0 || barriersValues[1] < 0.0) {
    return MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER;
  }

  if (!params.prefactor_options.default_value) {
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

    if ((pref1 > params.prefactor_options.max_value) ||
        (pref1 < params.prefactor_options.min_value)) {
      EONC_LOG_ERROR("Bad reactant-to-saddle prefactor: {}", pref1);
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    if ((pref2 > params.prefactor_options.max_value) ||
        (pref2 < params.prefactor_options.min_value)) {
      EONC_LOG_ERROR("Bad product-to-saddle prefactor: {}", pref2);
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    prefactorsValues[0] = pref1;
    prefactorsValues[1] = pref2;

  } else {
    prefactorsValues[0] = params.prefactor_options.default_value;
    prefactorsValues[1] = params.prefactor_options.default_value;
  }
  return MinModeSaddleSearch::STATUS_GOOD;
}

void ProcessSearchJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  std::ofstream out(resultsFilename, std::ios::binary);
  if (out) {
    out << std::format("{} termination_reason\n", status);
    out << std::format("{} termination_reason_text\n",
                       saddleSearch->describeStatus(status));
    out << std::format("{} random_seed\n", params.main_options.randomSeed);
    out << std::format(
        "{} potential_type\n",
        magic_enum::enum_name<PotType>(params.potential_options.potential));
    out << std::format("{} total_force_calls\n",
                       fCallsMin + fCallsSaddle + fCallsPrefactors);
    out << std::format("{} force_calls_minimization\n", fCallsMin);
    out << std::format("{} force_calls_saddle\n", fCallsSaddle);
    out << std::format("{:.12e} potential_energy_saddle\n",
                       saddle->getPotentialEnergy());
    out << std::format("{:.12e} potential_energy_reactant\n",
                       min1->getPotentialEnergy());
    out << std::format("{:.12e} potential_energy_product\n",
                       min2->getPotentialEnergy());
    out << std::format("{:.12e} barrier_reactant_to_product\n",
                       barriersValues[0]);
    out << std::format("{:.12e} barrier_product_to_reactant\n",
                       barriersValues[1]);
    if (params.saddle_search_options.method == "min_mode") {
      out << std::format("{:.12e} displacement_saddle_distance\n",
                         displacement->perAtomNorm(*saddle));
    } else {
      out << std::format("{:.12e} displacement_saddle_distance\n", 0.0);
    }
    if (params.saddle_search_options.method == "dynamics") {
      auto ds = dynamic_cast<DynamicsSaddleSearch &>(*saddleSearch);
      out << std::format("{:.12e} simulation_time\n",
                         ds.time * params.constants.timeUnit);
      out << std::format("{:.12e} md_temperature\n",
                         params.saddle_search_options.dynamics.temperature);
    }
    out << std::format("{} force_calls_prefactors\n", fCallsPrefactors);
    out << std::format("{:.12e} prefactor_reactant_to_product\n",
                       prefactorsValues[0]);
    out << std::format("{:.12e} prefactor_product_to_reactant\n",
                       prefactorsValues[1]);
  }

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  min1->matter2con(reactantFilename);

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  eonc::helpers::saveMode(modeFilename, saddle, saddleSearch->getEigenvector());

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  saddle->matter2con(saddleFilename);

  std::string productFilename("product.con");
  returnFiles.push_back(productFilename);
  min2->matter2con(productFilename);
}

void ProcessSearchJob::printEndState(int status) {
  auto msg = saddleSearch->describeStatus(status);
  if (status == MinModeSaddleSearch::STATUS_GOOD) {
    QUILL_LOG_DEBUG(log, "[Saddle Search] {}", msg);
  } else {
    QUILL_LOG_ERROR(log, "[Saddle Search] {}", msg);
  }
}
