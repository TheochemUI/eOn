
#include "SaddleSearchJob.h"
#include "EpiCenters.h"
#include "Potential.h"

#include <stdio.h>
#include <string>

using namespace std;

std::vector<std::string> SaddleSearchJob::run(void) {
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");

  if (params->checkpoint) {
    FILE *disp, *mode;
    disp = fopen("displacement_cp.con", "r");
    mode = fopen("mode_cp.dat", "r");
    if (disp != NULL && mode != NULL) {
      displacementFilename = "displacement_cp.con";
      modeFilename = "mode_cp.dat";
      SPDLOG_LOGGER_DEBUG(log, "Resuming from checkpoint");
    } else {
      SPDLOG_LOGGER_DEBUG(log, "No checkpoint files found");
    }
  }

  initial = std::make_shared<Matter>(pot, params);
  displacement = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);

  initial->con2matter(reactantFilename);

  if (params->saddleDisplaceType == EpiCenters::DISP_LOAD) {
    // displacement was passed from the server
    saddle->con2matter(displacementFilename);
  } else {
    // displacement and mode will be made on the client
    // in saddleSearch->initialize(...)
    *saddle = *initial;
  }
  AtomMatrix mode;
  if (params->saddleDisplaceType == EpiCenters::DISP_LOAD) {
    // mode was passed from the server
    mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
  }

  saddleSearch = std::make_unique<MinModeSaddleSearch>(
      saddle, mode, initial->getPotentialEnergy(), params, pot);

  int status;
  status = doSaddleSearch();
  printEndState(status);
  saveData(status);

  return returnFiles;
}

int SaddleSearchJob::doSaddleSearch() {
  Matter matterTemp(pot, params);
  long status;
  int f1{0};
  f1 = this->pot->forceCallCounter;
  try {
    status = saddleSearch->run();
  } catch (int e) {
    if (e == 100) {
      status = MinModeSaddleSearch::STATUS_POTENTIAL_FAILED;
    } else {
      printf("unknown exception: %i\n", e);
      throw e;
    }
  }

  if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_GPRDIMER) {
    fCallsSaddle = saddleSearch->forcecalls;
  } else {
    fCallsSaddle += this->pot->forceCallCounter - f1;
  }

  return status;
}

void SaddleSearchJob::saveData(int status) {
  FILE *fileResults, *fileSaddle, *fileMode;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");
  /// XXX: min_fcalls isn't quite right it should get them from
  //      the minimizer. But right now the minimizers are in
  //      the SaddleSearch object. They will be taken out eventually.

  fprintf(fileResults, "%d termination_reason\n", status);
  fprintf(fileResults, "saddle_search job_type\n");
  fprintf(fileResults, "%ld random_seed\n", params->randomSeed);
  fprintf(
      fileResults, "%s potential_type\n",
      std::string{magic_enum::enum_name<PotType>(params->potential)}.c_str());
  // XXX(rg): Adding one because the "total" number of calls should include the
  // a call to getPotentialEnergy() for the saddle
  fprintf(fileResults, "%li total_force_calls\n", this->pot->forceCallCounter + 1);
  fprintf(fileResults, "%d force_calls_saddle\n", fCallsSaddle);
  fprintf(fileResults, "%i iterations\n", saddleSearch->iteration);
  if (status != MinModeSaddleSearch::STATUS_POTENTIAL_FAILED) {
    fprintf(fileResults, "%f potential_energy_saddle\n",
            saddle->getPotentialEnergy());
    fprintf(fileResults, "%f final_eigenvalue\n",
            saddleSearch->getEigenvalue());
  }
  fprintf(fileResults, "%f potential_energy_reactant\n",
          initial->getPotentialEnergy());
  fclose(fileResults);

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  fileMode = fopen(modeFilename.c_str(), "wb");
  helper_functions::saveMode(fileMode, saddle, saddleSearch->getEigenvector());
  fclose(fileMode);

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  fileSaddle = fopen(saddleFilename.c_str(), "wb");
  saddle->matter2con(fileSaddle);
  fclose(fileSaddle);
}

void SaddleSearchJob::printEndState(int status) {
  SPDLOG_LOGGER_DEBUG(log, "[Saddle Search] Final status: ");

  if (status == MinModeSaddleSearch::STATUS_GOOD)
    SPDLOG_LOGGER_DEBUG(log, "Success");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
    SPDLOG_LOGGER_DEBUG(log,
                        "Initial displacement unable to reach convex region");

  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
    SPDLOG_LOGGER_DEBUG(log, "Barrier too high");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
    SPDLOG_LOGGER_DEBUG(log, "Too many iterations in concave region");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
    SPDLOG_LOGGER_DEBUG(log, "Too many iterations");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
    SPDLOG_LOGGER_DEBUG(log, "Saddle is not connected to initial state");

  else if (status == MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
    SPDLOG_LOGGER_DEBUG(log, "Prefactors not within window");

  else if (status == MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
    SPDLOG_LOGGER_DEBUG(log, "Hessian calculation failed");

  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
    SPDLOG_LOGGER_DEBUG(log, "Energy barrier not within window");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MINIMA)
    SPDLOG_LOGGER_DEBUG(log, "Minimizations from saddle did not converge");

  else if (status == MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
    SPDLOG_LOGGER_DEBUG(log, "Nonnegative initial mode, aborting");

  else if (status == MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
    SPDLOG_LOGGER_DEBUG(log, "Negative barrier detected");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
    SPDLOG_LOGGER_DEBUG(log, "No reaction found during MD trajectory");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
    SPDLOG_LOGGER_DEBUG(
        log, "Converged to stationary point with zero negative modes");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
    SPDLOG_LOGGER_DEBUG(log,
                        "No forward barrier was found along minimized band");

  else if (status == MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
    SPDLOG_LOGGER_DEBUG(log, "Zero mode abort.");

  else if (status == MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
    SPDLOG_LOGGER_DEBUG(log, "Optimizer error.");

  else
    SPDLOG_LOGGER_DEBUG(log, "Unknown status: {}!", status);

  return;
}
