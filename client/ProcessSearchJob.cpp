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
#include "BasinHoppingSaddleSearch.h"
#include "BiasedGradientSquaredDescent.h"
#include "DynamicsSaddleSearch.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "MinModeSaddleSearch.h"
#include "Optimizer.h"
#include "Prefactor.h"
#include <memory>
#include <spdlog/spdlog.h>

std::vector<std::string> ProcessSearchJob::run(void) {
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");
  size_t fctmp{0}; // force call temporary
  initial = std::make_shared<Matter>(pot, params);
  if (params.saddle_search_options.method == "min_mode" ||
      params.saddle_search_options.method == "basin_hopping" ||
      params.saddle_search_options.method == "bgsd") {
    displacement = std::make_shared<Matter>(pot, params);
  } else if (params.saddle_search_options.method == "dynamics") {
    displacement = NULL;
  }
  saddle = std::make_shared<Matter>(pot, params);
  min1 = std::make_shared<Matter>(pot, params);
  min2 = std::make_shared<Matter>(pot, params);

  if (!initial->con2matter(reactantFilename)) {
    printf("Stop\n");
    exit(1);
  }

  if (params.process_search_options.minimize_first) {
    SPDLOG_LOGGER_DEBUG(log, "Minimizing initial structure\n");
    fctmp = initial->getPotentialCalls();
    initial->relax();
    fCallsMin += initial->getPotentialCalls() - fctmp;
    SPDLOG_LOGGER_DEBUG(log, "Initial minimization took {} fcalls",
                        initial->getPotentialCalls() - fctmp);
  }

  barriersValues[0] = barriersValues[1] = 0;
  prefactorsValues[0] = prefactorsValues[1] = 0;

  if (params.saddle_search_options.method == "min_mode" ||
      params.saddle_search_options.method == "basin_hopping" ||
      params.saddle_search_options.method == "bgsd") {
    if (params.saddle_search_options.displace_type == EpiCenters::DISP_LOAD) {
      // displacement was passed from the server
      if (!saddle->con2matter(displacementFilename)) {
        printf("Stop\n");
        exit(1);
      }
      *min1 = *min2 = *initial;
    } else {
      // displacement and mode will be made on the client
      // in SaddleSearch->initialize(...)
      *saddle = *min1 = *min2 = *initial;
    }
  } else {
    *saddle = *min1 = *min2 = *initial;
  }

  AtomMatrix mode;

  if (params.saddle_search_options.method == "min_mode") {
    if (params.saddle_search_options.displace_type == EpiCenters::DISP_LOAD) {
      // mode was passed from the server
      mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
    }
    saddleSearch = std::make_unique<MinModeSaddleSearch>(
        saddle, mode, initial->getPotentialEnergy(), params, pot);
  } else if (params.saddle_search_options.method == "basin_hopping") {
    saddleSearch =
        std::make_unique<BasinHoppingSaddleSearch>(min1, saddle, pot, params);
  } else if (params.saddle_search_options.method == "dynamics") {
    saddleSearch = std::make_unique<DynamicsSaddleSearch>(saddle, params);
  } else if (params.saddle_search_options.method == "bgsd") {
    saddleSearch = std::make_unique<BiasedGradientSquaredDescent>(
        saddle, initial->getPotentialEnergy(), params);
  }

  int status = doProcessSearch();

  printEndState(status);
  saveData(status);

  // might have been forced to be equal if the structure passed to the client
  // when determining barrier and the prefactor

  return returnFiles;
}

int ProcessSearchJob::doProcessSearch(void) {
  Matter matterTemp(pot, params);
  long status;
  size_t fctmp{0};

  fctmp = pot->forceCallCounter;
  status = saddleSearch->run();
  fCallsSaddle += pot->forceCallCounter - fctmp;
  SPDLOG_DEBUG("Got {} calls in the saddle search, with previous {}",
               fCallsSaddle, fctmp);

  if (status != MinModeSaddleSearch::STATUS_GOOD) {
    return status;
  }

  // relax from the saddle point

  AtomMatrix posSaddle = saddle->getPositions();
  AtomMatrix displacedPos;

  *min1 = *saddle;

  displacedPos =
      posSaddle - saddleSearch->getEigenvector() *
                      params.process_search_options.minimization_offset;
  min1->setPositions(displacedPos);

  SPDLOG_LOGGER_DEBUG(log, "Starting Minimization 1");
  fctmp = min1->getPotentialCalls();
  bool converged =
      min1->relax(false, params.debug_options.write_movies, false, "min1");
  fCallsMin += min1->getPotentialCalls() - fctmp;
  SPDLOG_LOGGER_DEBUG(log, "Min1 minimization took {} fcalls",
                      min1->getPotentialCalls() - fctmp);

  if (!converged) {
    return MinModeSaddleSearch::STATUS_BAD_MINIMA;
  }

  *min2 = *saddle;
  displacedPos =
      posSaddle + saddleSearch->getEigenvector() *
                      params.process_search_options.minimization_offset;
  min2->setPositions(displacedPos);

  SPDLOG_LOGGER_DEBUG(log, "Starting Minimization 2");
  fctmp = min2->getPotentialCalls();
  converged =
      min2->relax(false, params.debug_options.write_movies, false, "min2");
  fCallsMin += min2->getPotentialCalls() - fctmp;
  SPDLOG_LOGGER_DEBUG(log, "Min2 minimization took {} fcalls",
                      min2->getPotentialCalls() - fctmp);

  if (!converged) {
    return MinModeSaddleSearch::STATUS_BAD_MINIMA;
  }

  // if min2 corresponds to initial state, swap min1 && min2
  if (!(initial->compare(*min1)) && initial->compare(*min2)) {
    matterTemp = *min1;
    *min1 = *min2;
    *min2 = matterTemp;
  }

  if ((initial->compare(*min1)) == false) {
    SPDLOG_LOGGER_DEBUG(log, "initial != min1");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  if (initial->compare(*min2)) {
    // both minima are the initial state
    SPDLOG_LOGGER_DEBUG(log, "both minima are the initial state");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  // use the structure passed to the client when determining
  // the barrier and prefactor for the forward process
  if (!params.process_search_options.minimize_first) {
    min1 = initial;
  }

  // Calculate the barriers
  barriersValues[0] = saddle->getPotentialEnergy() - min1->getPotentialEnergy();
  barriersValues[1] = saddle->getPotentialEnergy() - min2->getPotentialEnergy();

  if ((params.saddle_search_options.max_energy < barriersValues[0]) ||
      (params.saddle_search_options.max_energy < barriersValues[1])) {
    return MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER;
  }

  if (barriersValues[0] < 0.0 || barriersValues[1] < 0.0) {
    return MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER;
  }

  // calculate the prefactor
  if (!params.prefactor_options.default_value) {
    fctmp = min1->getPotentialCalls();
    int prefStatus;
    double pref1, pref2;
    // XXX: no get() calls
    prefStatus = Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                          min2.get(), pref1, pref2);
    if (prefStatus == -1) {
      printf("Prefactor: bad calculation\n");
      return MinModeSaddleSearch::STATUS_FAILED_PREFACTOR;
    }
    fCallsPrefactors += min1->getPotentialCalls() - fctmp;

    /* Check that the prefactors are in the correct range */
    if ((pref1 > params.prefactor_options.max_value) ||
        (pref1 < params.prefactor_options.min_value)) {
      cout << "Bad reactant-to-saddle prefactor: " << pref1 << endl;
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    if ((pref2 > params.prefactor_options.max_value) ||
        (pref2 < params.prefactor_options.min_value)) {
      cout << "Bad product-to-saddle prefactor: " << pref2 << endl;
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    prefactorsValues[0] = pref1;
    prefactorsValues[1] = pref2;

  } else {
    // use the default prefactor value specified
    prefactorsValues[0] = params.prefactor_options.default_value;
    prefactorsValues[1] = params.prefactor_options.default_value;
  }
  return MinModeSaddleSearch::STATUS_GOOD;
}

void ProcessSearchJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  {
    auto out = fmt::output_file(resultsFilename);
    out.print("{} termination_reason\n", status);
    out.print("{} random_seed\n", params.main_options.randomSeed);
    out.print("{} potential_type\n",
              std::string{magic_enum::enum_name<PotType>(
                  params.potential_options.potential)});
    out.print("{} total_force_calls\n",
              fCallsMin + fCallsSaddle + fCallsPrefactors);
    out.print("{} force_calls_minimization\n", fCallsMin);
    out.print("{} force_calls_saddle\n", fCallsSaddle);
    out.print("{:.12e} potential_energy_saddle\n",
              saddle->getPotentialEnergy());
    out.print("{:.12e} potential_energy_reactant\n",
              min1->getPotentialEnergy());
    out.print("{:.12e} potential_energy_product\n",
              min2->getPotentialEnergy());
    out.print("{:.12e} barrier_reactant_to_product\n", barriersValues[0]);
    out.print("{:.12e} barrier_product_to_reactant\n", barriersValues[1]);
    if (params.saddle_search_options.method == "min_mode") {
      out.print("{:.12e} displacement_saddle_distance\n",
                displacement->perAtomNorm(*saddle));
    } else {
      out.print("{:.12e} displacement_saddle_distance\n", 0.0);
    }
    if (params.saddle_search_options.method == "dynamics") {
      auto ds = dynamic_cast<DynamicsSaddleSearch &>(*saddleSearch);
      out.print("{:.12e} simulation_time\n",
                ds.time * params.constants.timeUnit);
      out.print("{:.12e} md_temperature\n",
                params.saddle_search_options.dynamics.temperature);
    }
    out.print("{} force_calls_prefactors\n", fCallsPrefactors);
    out.print("{:.12e} prefactor_reactant_to_product\n", prefactorsValues[0]);
    out.print("{:.12e} prefactor_product_to_reactant\n", prefactorsValues[1]);
  }

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  min1->matter2con(reactantFilename);

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  helper_functions::saveMode(modeFilename, saddle, saddleSearch->getEigenvector());

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  saddle->matter2con(saddleFilename);

  std::string productFilename("product.con");
  returnFiles.push_back(productFilename);
  min2->matter2con(productFilename);

  return;
}

void ProcessSearchJob::printEndState(int status) {
  SPDLOG_LOGGER_DEBUG(log, "[Saddle Search] Final status: ");

  if (status == MinModeSaddleSearch::STATUS_GOOD)
    SPDLOG_LOGGER_DEBUG(log, "Success");
  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
    SPDLOG_LOGGER_ERROR(log,
                        "Initial displacement unable to reach convex region");
  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
    SPDLOG_LOGGER_ERROR(log, "Barrier too high");
  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
    SPDLOG_LOGGER_ERROR(log, "Too many iterations in concave region");
  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
    SPDLOG_LOGGER_ERROR(log, "Too many iterations");
  else if (status == MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
    SPDLOG_LOGGER_ERROR(log, "Saddle is not connected to initial state");
  else if (status == MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
    SPDLOG_LOGGER_ERROR(log, "Prefactors not within window");
  else if (status == MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
    SPDLOG_LOGGER_ERROR(log, "Hessian calculation failed");
  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
    SPDLOG_LOGGER_ERROR(log, "Energy barrier not within window");
  else if (status == MinModeSaddleSearch::STATUS_BAD_MINIMA)
    SPDLOG_LOGGER_ERROR(log, "Minimizations from saddle did not converge");
  else if (status == MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
    SPDLOG_LOGGER_CRITICAL(log, "Nonnegative initial mode, aborting");
  else if (status == MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
    SPDLOG_LOGGER_ERROR(log, "Negative barrier detected");
  else if (status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
    SPDLOG_LOGGER_ERROR(log, "No reaction found during MD trajectory");
  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
    SPDLOG_LOGGER_ERROR(
        log, "Converged to stationary point with zero negative modes");
  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
    SPDLOG_LOGGER_ERROR(log,
                        "No forward barrier was found along minimized band");
  else if (status == MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
    SPDLOG_LOGGER_CRITICAL(log, "Zero mode abort.");
  else if (status == MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
    SPDLOG_LOGGER_ERROR(log, "Optimizer error.");
  else
    SPDLOG_LOGGER_ERROR(log, "Unknown status: %i!", status);
  return;
}
