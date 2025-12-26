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
#include "NudgedElasticBandJob.h"
#include "ConjugateGradients.h"
#include "NEBInitialPaths.hpp"
#include "Potential.h"

using namespace std;

std::vector<std::string> NudgedElasticBandJob::run(void) {
  NudgedElasticBand::NEBStatus status;
  int f1;

  string reactantFilename = helper_functions::getRelevantFile("reactant.con");
  string productFilename = helper_functions::getRelevantFile("product.con");

  string transitionStateFilename = helper_functions::getRelevantFile("ts.con");
  bool tsInterpolate = false;
  auto transitionState = std::make_shared<Matter>(pot, params);
  FILE *fhTransitionState = fopen("ts.con", "r");
  if (fhTransitionState != nullptr) {
    tsInterpolate = true;
    fclose(fhTransitionState);
    transitionState->con2matter(transitionStateFilename);
  } else {
    transitionState = nullptr;
  }

  auto initial = std::make_shared<Matter>(pot, params);
  auto final_state = std::make_shared<Matter>(pot, params);

  initial->con2matter(reactantFilename);
  final_state->con2matter(productFilename);

  // Endpoint minimization logic:
  // - If params->neb_options.endpoints.minimize is false: never minimize
  // endpoints.
  // - If params->neb_options.endpoints.minimize is true and
  // params->neb_options.initialization.input_path is empty: minimize endpoints.
  // - If params->nebMinimEP is true and params->nebIpath is NOT empty:
  //     -> minimize endpoints only if params->nebMinimEPIpath is true.
  // Log what decision was made so users can see behavior.
  bool shouldMinimizeEndpoints = false;
  if (!params->neb_options.endpoints.minimize) {
    SPDLOG_LOGGER_DEBUG(
        m_log, "minimize_endpoints == false: not minimizing endpoints.");
    shouldMinimizeEndpoints = false;
  } else {
    if (params->neb_options.initialization.input_path.empty()) {
      SPDLOG_LOGGER_DEBUG(m_log, "minimize_endpoints == true and nebIpath "
                                 "empty: minimizing endpoints.");
      shouldMinimizeEndpoints = true;
    } else {
      // nebIpath provided: only minimize if neb_options.endpoints.use_path_file
      // explicitly allowed.
      if (params->neb_options.endpoints.use_path_file) {
        SPDLOG_LOGGER_DEBUG(
            m_log,
            "minimize_endpoints == true and nebIpath provided, but "
            "minimize_endpoints_for_ipath == true: minimizing endpoints.");
        shouldMinimizeEndpoints = true;
      } else {
        SPDLOG_LOGGER_DEBUG(
            m_log,
            "minimize_endpoints == true but nebIpath provided and "
            "minimize_endpoints_for_ipath == false: not minimizing endpoints.");
        shouldMinimizeEndpoints = false;
      }
    }
  }

  if (shouldMinimizeEndpoints) {
    SPDLOG_LOGGER_DEBUG(m_log, "Minimizing reactant");
    // TODO(rg): Maybe when we have even more parameters, false can be set by
    // the user too..
    initial->relax(false, params->writeMovies, params->checkpoint, "react_neb",
                   "react_neb");
    // TODO(rg): How do we report the total E/F now? Currently this is just the
    // total total, people might want "per-stage" totals (but they can also get
    // them from the log.)
    SPDLOG_LOGGER_DEBUG(m_log, "Minimized reactant in ");
    SPDLOG_LOGGER_DEBUG(m_log, "Minimizing product");
    final_state->relax(false, params->writeMovies, params->checkpoint,
                       "prod_neb", "prod_neb");
  }

  auto neb =
      std::make_unique<NudgedElasticBand>(initial, final_state, params, pot);

  if (tsInterpolate) {
    AtomMatrix reactantToTS = transitionState->pbc(
        transitionState->getPositions() - initial->getPositions());
    AtomMatrix TSToProduct = transitionState->pbc(
        final_state->getPositions() - transitionState->getPositions());
    for (int image = 1; image <= neb->numImages; image++) {
      int mid = neb->numImages / 2 + 1;
      if (image < mid) {
        double frac = ((double)image) / ((double)mid);
        neb->path[image]->setPositions(initial->getPositions() +
                                       frac * reactantToTS);
      } else if (image > mid) {
        double frac =
            (double)(image - mid) / (double)(neb->numImages - mid + 1);
        neb->path[image]->setPositions(transitionState->getPositions() +
                                       frac * TSToProduct);
      } else if (image == mid) {
        neb->path[image]->setPositions(transitionState->getPositions());
      }
    }
  }

  // f1 = Potential::fcalls;
  status = neb->compute();
  // fCallsNEB += Potential::fcalls - f1;

  if (status == NudgedElasticBand::NEBStatus::GOOD) {
    neb->printImageData();
    neb->findExtrema();
  }

  printEndState(status);
  saveData(status, neb.get());

  return returnFiles;
}

void NudgedElasticBandJob::saveData(NudgedElasticBand::NEBStatus status,
                                    NudgedElasticBand *neb) {
  FILE *fileResults, *fileNEB;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");
  if (!fileResults) {
    SPDLOG_LOGGER_ERROR(m_log, "Failed to open {} for writing",
                        resultsFilename);
    return;
  }

  fprintf(fileResults, "%d termination_reason\n", static_cast<int>(status));
  fprintf(
      fileResults, "%s potential_type\n",
      std::string{magic_enum::enum_name<PotType>(params->potential)}.c_str());
  // fprintf(fileResults, "%ld total_force_calls\n", Potential::fcalls);
  // fprintf(fileResults, "%ld force_calls_neb\n", fCallsNEB);
  fprintf(fileResults, "%f energy_reference\n",
          neb->path[0]->getPotentialEnergy());
  fprintf(fileResults, "%li number_of_images\n", neb->numImages);

  for (long i = 0; i <= neb->numImages + 1; i++) {
    fprintf(fileResults, "%f image%li_energy\n",
            neb->path[i]->getPotentialEnergy() -
                neb->path[0]->getPotentialEnergy(),
            i);
    fprintf(fileResults, "%f image%li_force\n",
            neb->path[i]->getForces().norm(), i);

    // Only interior images have a meaningful projected force
    double proj_norm =
        (i >= 1 && i <= neb->numImages) ? neb->projectedForce[i]->norm() : 0.0;
    fprintf(fileResults, "%f image%li_projected_force\n", proj_norm, i);
  }

  fprintf(fileResults, "%li number_of_extrema\n", neb->numExtrema);
  for (long i = 0; i < neb->numExtrema; i++) {
    fprintf(fileResults, "%f extremum%li_position\n", neb->extremumPosition[i],
            i);
    fprintf(fileResults, "%f extremum%li_energy\n", neb->extremumEnergy[i], i);
  }
  fclose(fileResults);

  // Save the Full NEB Path
  std::string nebFilename("neb.con");
  returnFiles.push_back(nebFilename);
  fileNEB = fopen(nebFilename.c_str(), "wb");
  for (long i = 0; i <= neb->numImages + 1; i++) {
    neb->path[i]->matter2con(fileNEB);
  }
  fclose(fileNEB);

  // Save Discrete Saddle Point (Highest Energy Image)
  std::string spFilename("sp.con");
  neb->path[neb->maxEnergyImage]->matter2con(spFilename);
  returnFiles.push_back(spFilename);

  // Setup Dimer Configurations for Each Spline Peak
  if (params->neb_options.mmf_peaks.enabled && neb->numExtrema > 0) {
    int peakCount = 0;
    for (long i = 0; i < neb->numExtrema; i++) {
      // Filter 1: Only look at maxima (negative curvature)
      // Filter 2: Energy threshold (e.g., peak must be > 0.05 eV above
      // reactant)
      double relativeEnergy =
          neb->extremumEnergy[i] - neb->path[0]->getPotentialEnergy();

      if (neb->extremumCurvature[i] < 0 &&
          relativeEnergy > params->neb_options.mmf_peaks.tolerance) {
        double posFraction = neb->extremumPosition[i];
        int leftIdx = static_cast<int>(std::floor(posFraction));
        double f = posFraction - leftIdx;

        if (leftIdx < 0 || leftIdx >= neb->numImages + 1)
          continue;

        // 1. Write Interpolated Position (.con)
        Matter peakPos = helper_functions::neb_paths::interpolateImage(
            *neb->path[leftIdx], *neb->path[leftIdx + 1], f);
        std::string peakPosFile = fmt::format("peak{:02d}_pos.con", peakCount);
        peakPos.matter2con(peakPosFile);
        returnFiles.push_back(peakPosFile);

        // 2. Write Interpolated Tangent as standard mode.dat
        AtomMatrix peakMode = (1.0 - f) * (*neb->tangent[leftIdx]) +
                              f * (*neb->tangent[leftIdx + 1]);
        peakMode.normalize();

        std::string peakModeFile =
            fmt::format("peak{:02d}_mode.dat", peakCount);
        FILE *fMode = fopen(peakModeFile.c_str(), "w");
        if (fMode) {
          for (int row = 0; row < peakMode.rows(); ++row) {
            fprintf(fMode, "%12.6f %12.6f %12.6f\n", peakMode(row, 0),
                    peakMode(row, 1), peakMode(row, 2));
          }
          fclose(fMode);
          returnFiles.push_back(peakModeFile);
        }

        SPDLOG_LOGGER_INFO(
            m_log,
            "Generated MMF peak {:02d} at position {:.3f} (Energy: {:.3f} eV)",
            peakCount, posFraction, relativeEnergy);
        peakCount++;
      }
    }
  }

  returnFiles.push_back("neb.dat");
  neb->printImageData(true, std::numeric_limits<size_t>::max());
}

void NudgedElasticBandJob::printEndState(NudgedElasticBand::NEBStatus status) {
  SPDLOG_LOGGER_DEBUG(m_log, "Final state: ");
  if (status == NudgedElasticBand::NEBStatus::GOOD)
    SPDLOG_LOGGER_DEBUG(m_log, "Nudged elastic band, successful.");
  else if (status == NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS)
    SPDLOG_LOGGER_DEBUG(m_log, "Nudged elastic band, too many iterations.");
  else
    SPDLOG_LOGGER_WARN(m_log, "Unknown status: {}!", static_cast<int>(status));
  return;
}
