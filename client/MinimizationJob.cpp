#include "MinimizationJob.h"
#include "BaseStructures.h"
#include "HelperFunctions.h"
#include "Log.h"
#include "Matter.h"
#include "Optimizer.h"

std::vector<std::string> MinimizationJob::run(void) {
  string posInFilename("pos.con");
  string posOutFilename("min.con");

  if (params->checkpoint) {
    FILE *pos;
    pos = fopen("pos_cp.con", "r");
    if (pos != NULL) {
      posInFilename = "pos_cp.con";
      log("[Minimization] Resuming from checkpoint\n");
    } else {
      log("[Minimization] No checkpoint files found\n");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  Matter *pos = new Matter(params);
  pos->con2matter(posInFilename);

  printf("\nBeginning minimization of %s\n", posInFilename.c_str());

  bool converged;
  try {
    converged = pos->relax(false, params->writeMovies, params->checkpoint,
                           "minimization", "pos");
    if (converged) {
      status = RunStatus::GOOD;
      printf("Minimization converged within tolerence\n");
    } else {
      status = RunStatus::MAX_ITERATIONS;
      printf("Minimization did not converge to tolerence!\n"
             "Maybe try to increase max_iterations?\n");
    }
  } catch (int e) {
    if (e == 100) {
      status = RunStatus::POTENTIAL_FAILED;
    } else {
      throw e;
    }
  }

  printf("Saving result to %s\n", posOutFilename.c_str());
  pos->matter2con(posOutFilename);
  if (status != RunStatus::POTENTIAL_FAILED) {
    printf("Final Energy: %f\n", pos->getPotentialEnergy());
  }

  FILE *fileResults;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");

  fprintf(fileResults, "%s termination_reason\n",
          (helper_functions::getRunStatusName(status)).c_str());
  fprintf(fileResults, "minimization job_type\n");
  fprintf(fileResults, "%s potential_type\n",
          helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  if (status != RunStatus::POTENTIAL_FAILED) {
    fprintf(fileResults, "%f potential_energy\n", pos->getPotentialEnergy());
  }
  fclose(fileResults);

  return returnFiles;
}
